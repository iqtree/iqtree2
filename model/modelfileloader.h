//
//  modelfileloader.h
//  alignment
//
//  Created by James Barbetti on 1/4/21.
//

#ifndef modelfileloader_h
#define modelfileloader_h

#include <string> //for std::string
#include <yaml-cpp/yaml.h> //for YAML::Node
#include "modelinfo.h" //for ModelInfoFromYAMLFile

class ModelFileLoader {
public:
    const char* file_path;
    const std::string model_name;
   
    ModelFileLoader(const char* path);
    std::string stringScalar(const YAML::Node& node,
                             const char* key);
    bool booleanScalar(const YAML::Node& node,
                       const char* key);
    void complainIfNot(bool check_me,
                       std::string error_message);
    double toDouble(const YAML::Node& i, double default_val);
    
    ParameterRange parseRange(const YAML::Node& node, const char* key);
    
    void parseYAMLModelParameters(const YAML::Node& params,
                                  ModelInfoFromYAMLFile& info);
    
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
    
    /* Example of a rate matrix
     
     rateMatrix:
     - [      , r1*r2, r2*f3, r3*f4 ]
     - [ r1*f1,      , r4*f3, r5*f4 ]
     - [ r2*f1, r4*f2,      , f4    ]
     - [ r3*f1, r5*f2, f3   ,       ]

     */
    
    void parseRateMatrix(const YAML::Node& rate_matrix,
                         ModelInfoFromYAMLFile& info);
        
    void parseRateMatrixExpressions(ModelInfoFromYAMLFile& info);
    
    void parseYAMLSubstitutionModel(const YAML::Node& substitution_model,
                                    const std::string& name_of_model,
                                    ModelInfoFromYAMLFile& info);
};

#endif /* modelfileloader_hpp */
