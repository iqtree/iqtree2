//
//modelfileloader.h
//Created by James Barbetti on 01-Apr-2021.
//

#ifndef modelfileloader_h
#define modelfileloader_h

#include <string> //for std::string
#include <yaml-cpp/yaml.h> //for YAML::Node
#include "modelinfofromyamlfile.h" //for ModelInfoFromYAMLFile class

namespace ModelExpression {
    class Assignment;
};

class ModelFileLoader {
protected:    
    const char*       file_path;
    const std::string model_name;
    void  handleInheritance        (ModelInfoFromYAMLFile& info, 
                                    ModelListFromYAMLFile& list,
                                    PhyloTree*             report_to_tree);
    void  setModelStateFrequency   (const YAML::Node&      substitution_model, 
                                    ModelInfoFromYAMLFile& info,
                                    PhyloTree*             report_to_tree);
    void setParameterSubscriptRange(ModelInfoFromYAMLFile& info,
                                    YAMLFileParameter&     p);
    void setParameterType          (const YAML::Node&      param,
                                    bool                   overriding,
                                    ModelInfoFromYAMLFile& info,
                                    YAMLFileParameter&     p);
    void setParameterTolerance     (const YAML::Node&      param,
                                    bool                   overriding,
                                    ModelInfoFromYAMLFile& info,
                                    YAMLFileParameter&     p);
    void setParameterValue         (const YAML::Node&      param,
                                    bool                   overriding,
                                    ModelInfoFromYAMLFile& info,
                                    YAMLFileParameter&     p);
    bool isAParameterOverride      (ModelInfoFromYAMLFile& info,
                                    YAMLFileParameter& p);
    void parseYAMLSubModels        (const YAML::Node& substitution_model,
                                    ModelInfoFromYAMLFile& info,
                                    ModelListFromYAMLFile& list,
                                    PhyloTree* report_to_tree);
    void parseYAMLModelStringProperties(const YAML::Node& substitution_model,
                                        ModelInfoFromYAMLFile& info,
                                        PhyloTree* report_to_tree);
    void parseYAMLModelWeightAndScale(const YAML::Node& substitution_model,
                                      ModelInfoFromYAMLFile& info,
                                      PhyloTree* report_to_tree);

public:
   
    ModelFileLoader(const char* path);
    std::string stringScalar(const YAML::Node& node,
                             const char* key,
                             const char* default_value);
    bool   booleanScalar(const YAML::Node& node,
                         const char*       key,
                         const bool        default_value);
    int    integerScalar(const YAML::Node& node,
                         const char*       key,
                         const int         default_value);
    double doubleScalar (const YAML::Node& node,
                         const char*       key,
                         const double      default_value);
    
    void   complainIfSo (bool check_me, std::string error_message);
    void   complainIfNot(bool check_me, std::string error_message);
    double toDouble     (const YAML::Node& i, double default_val,
                         ModelInfoFromYAMLFile& info, std::string context);
    
    ModelParameterRange parseRange(const YAML::Node& node, const char* key,
                                   const ModelParameterRange& default_value,
                                   ModelInfoFromYAMLFile& info);
    
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
    double setConstraint(ModelExpression::Assignment* a,
                         ModelInfoFromYAMLFile&       info,
                         const std::string& constraint_string,
                         PhyloTree*         report_to_tree);
    
    void parseYAMLMixtureModels(const YAML::Node& mixture_models,
                                ModelInfoFromYAMLFile& info,
                                ModelListFromYAMLFile& list,
                                PhyloTree* report_to_tree);
    void dumpMatrixTo(const char* name, ModelInfoFromYAMLFile& info,
                      const StringMatrix& matrix,
                      int rank,
                      const std::string& formula, std::stringstream &out);



    void parseRateMatrix(const YAML::Node& rate_matrix,
                         ModelInfoFromYAMLFile& info,
                         PhyloTree* report_to_tree);
        
    void parseYAMLModel(const YAML::Node& substitution_model,
                        const std::string& name_of_model,
                        ModelInfoFromYAMLFile& info,
                        ModelListFromYAMLFile& list,
                        ModelInfoFromYAMLFile* parent_model,
                        PhyloTree* report_to_tree);

};

#endif /* modelfileloader_h */
