#include "yamlmodelwrapper.h"
#include "model/modelinfo.h"

YAMLModelBinary::YAMLModelBinary(ModelInfoFromYAMLFile& info,
                                 bool make_copy, const char *model_name, 
                                 std::string model_params, StateFreqType freq, 
                                 std::string freq_params, PhyloTree *tree, 
                                 PhyloTree* report_to_tree)
        : super(info, make_copy, tree, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelCodon::YAMLModelCodon(ModelInfoFromYAMLFile& info,
                               bool        make_copy,    const char*   model_name, 
                               std::string model_params, StateFreqType freq, 
                               std::string freq_params,  PhyloTree*    tree, 
                               PhyloTree*  report_to_tree)
        : super(info, make_copy, tree, report_to_tree) {
    setReversible(info.isReversible());
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNA::YAMLModelDNA(ModelInfoFromYAMLFile& info,
                           bool        make_copy,    const char*   model_name, 
                           std::string model_params, StateFreqType freq, 
                           std::string freq_params,  PhyloTree*    tree, 
                           PhyloTree*  report_to_tree)
        : super(info, make_copy, tree, report_to_tree) {
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNAError::YAMLModelDNAError(ModelInfoFromYAMLFile& info,
                                     bool        make_copy,    const char* model_name, 
                                     std::string model_params, StateFreqType freq, 
                                     std::string freq_params,  PhyloTree*    tree, 
                                     PhyloTree* report_to_tree)                    
        : super(info, make_copy, tree, report_to_tree) {
    init(model_name, model_params, freq,
            freq_params, report_to_tree);        
    setRateMatrixFromModel();
}

bool YAMLModelDNAError::getVariables(double *variables) {
    bool changed = super::getVariables(variables);
    if (changed && !fix_epsilon) {
        epsilon = model_info->getVariableValue("epsilon");
    }
    return changed;
}

YAMLModelMorphology::YAMLModelMorphology(ModelInfoFromYAMLFile& info,
                                         bool        make_copy,    const char* model_name, 
                                         std::string model_params, StateFreqType freq, 
                                         std::string freq_params,  PhyloTree*    tree, 
                                         PhyloTree*  report_to_tree)
        : super(info, make_copy, tree, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelProtein::YAMLModelProtein(ModelInfoFromYAMLFile& info, 
                                   bool        make_copy,    const char*   model_name, 
                                   std::string model_params, StateFreqType freq, 
                                   std::string freq_params,  ModelsBlock*  block,
                                   PhyloTree*  tree,         PhyloTree*    report_to_tree)
        : super(info, make_copy, tree, report_to_tree) {
    setModelsBlock(block);
    setNumberOfStates(20);
    setReversible(info.isReversible());
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setNumberOfStates(20); //Why again?
    setRateMatrixFromModel();
}

YAMLModelMixture::YAMLModelMixture(ModelInfoFromYAMLFile& info,
                                   bool          make_copy, const char*  model_name, 
                                   StateFreqType freq,      ModelsBlock* models_block,
                                   PhyloTree*    tree,      PhyloTree*   report_to_tree)
    : super(info, make_copy, tree, report_to_tree) {
    //as per init() and initMixture()
    ASSERT(info.isMixtureModel());

    StateFreqType old_freq = info.getFrequencyType();
    StateFreqType new_freq = old_freq;
    if (new_freq == StateFreqType::FREQ_UNKNOWN) {
        new_freq = freq;
    }
    if (new_freq == StateFreqType::FREQ_UNKNOWN) {
        new_freq = StateFreqType::FREQ_USER_DEFINED;
    }
    info.setFrequencyType(new_freq);

    bool optimize_weights = false;

    const std::string no_parameters; //parameters are passed to the mixture.
    full_name = "MIX";
    full_name += OPEN_BRACKET;
    const char* separator = "";
    DoubleVector weights;
    for (auto child: model_info->getMixedModels()) {
        auto model = ModelListFromYAMLFile::getModelByReference
                    (*child, tree, info.getFrequencyType(),
                    models_block, no_parameters, 
                    report_to_tree);
        models.push_back(model);
        weights.push_back(child->getModelWeight());
        if (!optimize_weights) {
            optimize_weights = !child->isModelWeightFixed();
        }
        //what about weights?!
        full_name += separator;
        full_name += child->getName();
        separator = ",";
    }
    full_name += CLOSE_BRACKET;

    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, 
                  "optimize_weights=" << optimize_weights);

    checkProportionsAndWeights(weights);
    setOptimizationSteps(optimize_weights);
    checkModelReversibility();
    decomposeRateMatrix();

    phylo_tree = tree;

    setRateMatrixFromModel();
}

bool YAMLModelMixture::isMixtureModel() {
    return true;
}

void YAMLModelMixture::setRateMatrixFromModel() {
    //
    //Called when variables have been changed, for at least some of
    //the child models associated with this mixture, during initialization.
    //Each child model has its own copy of the ModelInfoFromYAMLFile
    //(variables and all), which might be out of date; the variables
    //updated in the mixture models in model_info.getMixedModels()
    //need to be copied across to the copies in the child models
    //(and the rate matrices of the child models need to be recalculated).
    //
    //See also afterVariablesChanged() which copies in the other direction
    //during optimization.
    //
    for (ModelMarkov* model : models) {
        model->setRateMatrixFromModel();
    }
}

void YAMLModelMixture::afterVariablesChanged() {
    //
    //Called when variables have been changed, for at least some of
    //the child models associated with this mixture, during optimization.
    //
}

void YAMLModelMixture::afterWeightsChanged() {
    int nmix = getNMixtures();
    if (1<nmix) {
        int i = 0;
        model_info->updateModelVariablesByType(prop, getNMixtures(), true,
                                               ModelParameterType::WEIGHT, i);
        return;
    }
}

YAMLRateFree::YAMLRateFree(PhyloTree* tree, PhyloTree* report_to_tree,
                           const ModelInfoFromYAMLFile& info)
        : super(info, report_to_tree) {
    std::string algorithm = info.getOptimizationAlgorithm();
    optimize_alg = algorithm.empty() ? optimize_alg : algorithm;
    gamma_shape  = 1;

    //num_rate_cats, gamma_shape,
    //freerate_params, !fused_mix_rate,
    //params.optimize_alg_freerate, tree
}

void YAMLRateFree::updateRateClassFromModelVariables() {
    int rate_count = model_info.getNumberOfRateCategories();
    int prop_count = model_info.getNumberOfProportions();
    int rate_ix    = 0;
    int prop_ix    = 0;
    model_info.readModelVariablesByType(rates, rate_count, true,
                                        ModelParameterType::RATE, rate_ix);
    model_info.readModelVariablesByType(prop,  prop_count, true,
                                        ModelParameterType::PROPORTION, prop_ix);
}

void YAMLRateFree::sortUpdatedRates() {
    super::sortUpdatedRates();
    int rate_count = model_info.getNumberOfRateCategories();
    int prop_count = model_info.getNumberOfProportions();
    int rate_ix    = 0;
    int prop_ix    = 0;
    model_info.updateModelVariablesByType(rates, rate_count, true,
                                          ModelParameterType::RATE, rate_ix);
    model_info.updateModelVariablesByType(prop,  prop_count, true,
                                          ModelParameterType::PROPORTION, prop_ix);
}