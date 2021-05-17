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