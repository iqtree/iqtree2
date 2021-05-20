#include "yamlmodelwrapper.h"

YAMLModelBinary::YAMLModelBinary(const char *model_name, std::string model_params,
                StateFreqType freq, std::string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelCodon::YAMLModelCodon(const char *model_name, std::string model_params,
                StateFreqType freq, std::string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
    setReversible(info.isReversible());
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNA::YAMLModelDNA(const char *model_name, string model_params,
                StateFreqType freq, string freq_params,
                PhyloTree *tree, PhyloTree* report_to_tree,
                const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelDNAError::YAMLModelDNAError(const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree *tree, PhyloTree* report_to_tree,
                    const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
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
                const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
    init(model_name, model_params, freq, freq_params, report_to_tree);
    setRateMatrixFromModel();
}

YAMLModelProtein::YAMLModelProtein(ModelsBlock* block,
                    const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree *tree, PhyloTree* report_to_tree,
                    const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
    setModelsBlock(block);
    setNumberOfStates(20);
    setReversible(info.isReversible());
    init(model_name, model_params, freq,
            freq_params, report_to_tree);
    setNumberOfStates(20); //Why again?
    setRateMatrixFromModel();
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
