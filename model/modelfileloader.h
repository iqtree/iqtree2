//
//modelfileloader.h
//Created by James Barbetti on 01-Apr-2021.
//

#ifndef modelfileloader_h
#define modelfileloader_h

#include <string> //for std::string
#include <yaml-cpp/yaml.h> //for YAML::Node
#include "modelinfo.h" //for ModelInfoFromYAMLFile

class ModelFileLoader {
public:
    const char*       file_path;
    const std::string model_name;
   
    ModelFileLoader(const char* path);
    std::string stringScalar(const YAML::Node& node,
                             const char* key,
                             const char* default_value);
    bool booleanScalar(const YAML::Node& node,
                       const char* key,
                       const bool  default_value);
    int  integerScalar(const YAML::Node& node,
                       const char* key,
                       const int   default_value);
    
    void   complainIfSo (bool check_me, std::string error_message);
    void   complainIfNot(bool check_me, std::string error_message);
    double toDouble     (const YAML::Node& i, double default_val);
    
    ModelParameterRange parseRange(const YAML::Node& node, const char* key,
                                   const ModelParameterRange& default_value);
    
    void parseYAMLModelParameters(const YAML::Node& params,
                                  ModelInfoFromYAMLFile& info,
                                  PhyloTree* report_to_tree);
    void parseModelParameter (const YAML::Node& param,
                              std::string name,
                              ModelInfoFromYAMLFile& info,
                              PhyloTree* report_to_tree);
    void parseMatrixParameter(const YAML::Node& param,
                              std::string name,
                              ModelInfoFromYAMLFile& info,
                              PhyloTree* report_to_tree);

    YAMLFileParameter
         addDummyFrequencyParameterTo(ModelInfoFromYAMLFile& info,
                                      PhyloTree* report_to_tree);
    void parseYAMLModelConstraints(const YAML::Node& params,
                                   ModelInfoFromYAMLFile& info,
                                   PhyloTree* report_to_tree);
    void parseYAMLMixtureModels(const YAML::Node& mixture_models,
                                ModelInfoFromYAMLFile& info,
                                ModelListFromYAMLFile& list,
                                PhyloTree* report_to_tree);
    
    template <class S>
    void dumpRateMatrixTo(const ModelInfoFromYAMLFile& info, S &out) {
        std::stringstream dump;
        for (auto r : info.rate_matrix_expressions) {
            const char* separator = "";
            for (auto c: r) {
                dump << separator << c;
                separator = " : ";
            }
            dump << "\n";
        }
        out << "Rate matrix for " << info.model_name
            << " is...\n" << dump.str();
    }

    void parseRateMatrix(const YAML::Node& rate_matrix,
                         ModelInfoFromYAMLFile& info,
                         PhyloTree* report_to_tree);
        
    void parseYAMLSubstitutionModel(const YAML::Node& substitution_model,
                                    const std::string& name_of_model,
                                    ModelInfoFromYAMLFile& info,
                                    ModelListFromYAMLFile& list,
                                    ModelInfoFromYAMLFile* parent_model,
                                    PhyloTree* report_to_tree);
};

#endif /* modelfileloader_h */
