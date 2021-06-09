#include "yamlmodelwrapper.h"

YAMLModelBinary::YAMLModelBinary(const char *model_name, std::string model_params,
                StateFreqType freq, std::string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelCodon::YAMLModelCodon(const char *model_name, std::string model_params,
                StateFreqType freq, std::string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    setReversible(info.isReversible());
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNA::YAMLModelDNA(const char *model_name, string model_params,
                StateFreqType freq, string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNAError::YAMLModelDNAError(const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree *tree, PhyloTree* report_to_tree,
                    const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    init(model_name, model_params, freq,
            freq_params, report_to_tree);        
    setRateMatrixFromModel();
}

bool YAMLModelDNAError::getVariables(double *variables) {
    bool changed = super::getVariables(variables);
    if (changed && !fix_epsilon) {
        epsilon = model_info.getVariableValue("epsilon");
    }
    return changed;
}

YAMLModelMorphology::YAMLModelMorphology(const char *model_name, std::string model_params,
                StateFreqType freq, std::string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelProtein::YAMLModelProtein(ModelsBlock* block,
                    const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree *tree, PhyloTree* report_to_tree,
                    const ModelInfoFromYAMLFile& info)
        : super(info, tree, report_to_tree) {
    setModelsBlock(block);
    setNumberOfStates(20);
    setReversible(info.isReversible());
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setNumberOfStates(20); //Why again?
    setRateMatrixFromModel();
}

YAMLModelMixture::YAMLModelMixture(ModelInfoFromYAMLFile& info,
                                   PhyloTree *tree,
                                   ModelsBlock* models_block,
                                   PhyloTree* report_to_tree)
    : super(info, tree, report_to_tree) {
    //as per init() and initMixture()
    ASSERT(info.isMixtureModel());

    StateFreqType freq = info.getFrequencyType();
    if (freq == StateFreqType::FREQ_UNKNOWN) {
        freq = StateFreqType::FREQ_USER_DEFINED;
    }

    if (freq == StateFreqType::FREQ_MIXTURE) {
        //Todo: Do what?!        
    }

    bool optimize_weights = false;

    const std::string no_parameters; //parameters are passed to the mixture.
    full_name = "MIX";
    full_name += OPEN_BRACKET;
    const char* separator = "";
    DoubleVector weights;
    for (auto child_name_and_info: model_info.getMixedModels()) {
        ModelInfoFromYAMLFile& child = child_name_and_info.second;
        ModelInfoFromYAMLFile* child_model_info = nullptr;
        auto model = ModelListFromYAMLFile::getModelByReference
                    (child, tree, info.getFrequencyType(),
                    models_block, no_parameters, 
                    child_model_info, report_to_tree);
        mixed_model_infos.push_back(child_model_info);
        child_model_info->setParentModel(&model_info);
        models.push_back(model);
        weights.push_back(child.getModelWeight());
        if (!optimize_weights) {
            optimize_weights = !child.isModelWeightFixed();
        }
        //what about weights?!
        full_name += separator;
        full_name += child.getName();
        separator = ",";
    }
    full_name += CLOSE_BRACKET;

    setRateMatrixFromModel();

    //Do more of what ModelMixture::initMixture() does.

    checkProportionsAndWeights(weights);
    setOptimizationSteps(optimize_weights);
    checkModelReversibility();
    decomposeRateMatrix();

    phylo_tree = tree;
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
    int model_num = 0;
    MapOfModels& map_of_models = model_info.getMixedModels();
    ASSERT(map_of_models.size() == mixed_model_infos.size());
    ASSERT(map_of_models.size() == models.size());
    for (auto child_name_and_info: map_of_models) {
        ModelInfoFromYAMLFile& child = child_name_and_info.second;
        ModelInfoFromYAMLFile* copy  = mixed_model_infos[model_num];
        copy->copyVariablesFrom(&child);
        ModelMarkov* model = models[model_num];
        model->setRateMatrixFromModel();
        ++model_num;
    }
}

void YAMLModelMixture::afterVariablesChanged() {
    //
    //Called when variables have been changed, for at least some of
    //the child models associated with this mixture, during optimization.
    //Each child model has its own copy of the ModelInfoFromYAMLFile
    //(variables and all), which will have been updated, and the updates
    //need to be copied back to what is in model_info's map of mixed models.
    //
    int model_num = 0;
    MapOfModels& map_of_models = model_info.getMixedModels();
    ASSERT(map_of_models.size() == mixed_model_infos.size());
    ASSERT(map_of_models.size() == models.size());
    for (auto child_name_and_info: map_of_models) {
        ModelInfoFromYAMLFile& child = child_name_and_info.second;   //original
        ModelInfoFromYAMLFile* copy  = mixed_model_infos[model_num]; //changed
        child.copyVariablesFrom(copy);  
        ++model_num;
    }
}

void YAMLModelMixture::afterWeightsChanged() {
    int nmix = getNMixtures();
    if (1<nmix) {
        int i = 0;
        model_info.updateModelVariablesByType(prop, getNMixtures(), true,
                                              ModelParameterType::WEIGHT, i);
        return;
    }
}

YAMLRateFree::YAMLRateFree(PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
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
