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

    bool doesStringEndInNumber      (const std::string& input,
                                     std::string& stem, 
                                     std::string& numeric_suffix) const;
    bool parseModelNameAndParameters(const std::string& input,
                                     std::string& model_name, 
                                     std::string& parameters) const;         
    void setModelStateFrequency     (const YAML::Node&      substitution_model, 
                                     ModelInfoFromYAMLFile& info,
                                     LoggingTarget*         logging_target);
    void setParameterSubscriptRange (ModelInfoFromYAMLFile& info,
                                     YAMLFileParameter&     p);
    void setParameterType           (const YAML::Node&      param,
                                     bool                   overriding,
                                     ModelInfoFromYAMLFile& info,
                                     YAMLFileParameter&     p);
    void setParameterTolerance      (const YAML::Node&      param,
                                     bool                   overriding,
                                     ModelInfoFromYAMLFile& info,
                                     YAMLFileParameter&     p);
    void setParameterValue          (const YAML::Node&      param,
                                     bool                   overriding,
                                     ModelInfoFromYAMLFile& info,
                                     YAMLFileParameter&     p);
    bool isAParameterOverride       (const ModelInfoFromYAMLFile& info,
                                     YAMLFileParameter&           p);
    void parseYAMLSubModels         (Params&                params,
                                     const YAML::Node&      substitution_model,
                                     ModelInfoFromYAMLFile& info,
                                     ModelListFromYAMLFile& list,
                                     LoggingTarget*         logging_target);
    void parseYAMLModelInheritance  (Params&                params,
                                     const YAML::Node&      substitution_model,
                                     ModelInfoFromYAMLFile& info,
                                     const ModelListFromYAMLFile& list,
                                     LoggingTarget*         logging_target);
    void parseYAMLModelStringProperties(const YAML::Node& substitution_model,
                                        ModelInfoFromYAMLFile& info,
                                        LoggingTarget*       logging_target);
    void parseYAMLModelWeightAndScale(const YAML::Node&      substitution_model,
                                      ModelInfoFromYAMLFile& info,
                                      LoggingTarget*         logging_target);
    void inheritOneModel            (ModelInfoFromYAMLFile& info,
                                     bool           must_be_rate_models,
                                     const ModelInfoFromYAMLFile* ancestor,
                                     LoggingTarget* logging_target,
                                     bool&          have_first_parent);

public:
    explicit ModelFileLoader(const char* path);
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
                                  LoggingTarget* logging_target);
    void parseModelParameter (const YAML::Node& param,
                              const std::string& name,
                              ModelInfoFromYAMLFile& info,
                              LoggingTarget* logging_target);
    void parseMatrixParameter(const YAML::Node& param,
                              std::string name,
                              ModelInfoFromYAMLFile& info,
                              LoggingTarget* logging_target);

    YAMLFileParameter
         addDummyFrequencyParameterTo(ModelInfoFromYAMLFile& info,
                                      LoggingTarget* logging_target);
    void parseYAMLModelConstraints(const YAML::Node& params,
                                   ModelInfoFromYAMLFile& info,
                                   LoggingTarget* logging_target);
    double setConstraint(ModelExpression::Assignment* a,
                         ModelInfoFromYAMLFile&       info,
                         const std::string& constraint_string,
                         LoggingTarget*     logging_target);
    
    void parseYAMLMixtureModels(Params& params,
                                const YAML::Node& mixture_models,
                                ModelInfoFromYAMLFile& info,
                                ModelListFromYAMLFile& list,
                                LoggingTarget* logging_target);
    void dumpMatrixTo(const char* name, ModelInfoFromYAMLFile& info,
                      const StringMatrix& matrix,
                      int rank,
                      const std::string& formula, std::stringstream &out);

    void parseRateMatrix(const YAML::Node& rate_matrix,
                         ModelInfoFromYAMLFile& info,
                         LoggingTarget* logging_target);
        
    void parseYAMLModel(Params& params,
                        const YAML::Node& substitution_model,
                        const std::string& name_of_model,
                        ModelInfoFromYAMLFile& info,
                        ModelListFromYAMLFile& list,
                        ModelInfoFromYAMLFile* parent_model,
                        LoggingTarget* logging_target);

    void handleInheritance(Params&                      params,
                           ModelInfoFromYAMLFile&       info, 
                           const ModelListFromYAMLFile& list,
                           std::string                  inheritance_list,
                           bool                         must_be_rate_models,
                           LoggingTarget*               logging_target);

};

#endif /* modelfileloader_h */
