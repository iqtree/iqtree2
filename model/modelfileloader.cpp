//
// modelfileloader.cpp
// Created by James Barbetti on 01-Apr-2021.
//

#include "modelfileloader.h"
#include <utils/stringfunctions.h>
   
ModelFileLoader::ModelFileLoader(const char* path): file_path(path) {
}
        
std::string ModelFileLoader::stringScalar(const YAML::Node& node,
                                const char* key) {
    auto cite = node[key];
    return cite ? cite.Scalar() : "";
}
    
bool ModelFileLoader::booleanScalar(const YAML::Node& node,
                          const char* key) {
    std::string s = string_to_lower(stringScalar(node, key));
    return s == "true" || s == "yes" || s == "t" || s == "y" || s == "1";
}
    
void ModelFileLoader::complainIfNot(bool check_me,
                          std::string error_message) {
    if (!check_me) {
        outError(error_message);
    }
}
    
double ModelFileLoader::toDouble(const YAML::Node& i, double default_val) {
    if (!i.IsScalar()) {
        return default_val;
    }
    std::string double_string = i.Scalar();
    return convert_double_nothrow(double_string.c_str(), default_val);
}
    
ModelParameterRange ModelFileLoader::parseRange(const YAML::Node& node, const char* key) {
    ModelParameterRange range;
    auto r = node[key];
    if (!r) {
        return range;
    }
    //Todo: what if range is a string?
    if (r.IsSequence()) {
        int ix = 0;
        for (auto i : r) {
            if (ix==0) {
                range.first  = toDouble(i, 0);
            } else if (ix==1) {
                range.second = toDouble(i, 0);
            } else {
                //Throw: Range ought to be a two-value sequence
                //       Three values... one too many.
                //(for now, 3rd and later items in sequence
                // are just ignored).
            }
            ++ix;
        }
        range.is_set = ( 1 < ix );
    }
    return range;
}
    
void ModelFileLoader::parseYAMLModelParameters(const YAML::Node& params,
                              ModelInfoFromYAMLFile& info) {
    //
    //Assumes: params is a sequence of parameter declarations
    //
    for (auto param: params) {
        YAMLFileParameter p;
        p.name       = stringScalar(param, "name");
        auto bracket = p.name.find('(');
        if ( bracket != std::string::npos ) {
            p.is_subscripted = true;
            const char* range = p.name.c_str() + bracket + 1;
            int ix;
            p.minimum_subscript = convert_int(range, ix);
            if (strncmp(range + ix, "..", 2)==0) {
                range = range + ix + 2;
                p.maximum_subscript = convert_int(range, ix);
            } else {
                p.maximum_subscript = p.minimum_subscript;
                p.minimum_subscript   = 1;
            }
            if (strncmp(range + ix, ")", 1)!=0) {
                throw "Subscript range does not end with right parenthesis";
            }
            p.name = p.name.substr(0, bracket);
        }
        p.type_name = string_to_lower(stringScalar(param, "type"));
        if (p.type_name=="rate") {
            p.type = ModelParameterType::RATE;
        } else if (p.type_name=="frequency") {
            p.type = ModelParameterType::FREQUENCY;
        } else if (p.type_name=="weight") {
            p.type = ModelParameterType::WEIGHT;
        } else {
            p.type = ModelParameterType::OTHER;
        }
        auto count = p.maximum_subscript - p.minimum_subscript + 1;
        ASSERT(0<count);
        double dv   = 0.0; //default initial value
        if (p.type==ModelParameterType::FREQUENCY) {
            dv = 1.0 / (double)count;
                        //Todo: should be 1.0 divided by number of states
                        //determined from the data type (info.data_type_name ?)
                        //Or 1 divided by the number of parameters.
        } else if (p.type==ModelParameterType::RATE) {
            dv = 1.0;
        } else if (p.type==ModelParameterType::WEIGHT) {
            dv = 1.0 / (double)count;    //Todo: Should be 1.0 divided by # of parameters
        }
        std::string value_string = stringScalar(param, "initValue");
        p.range                  = parseRange  (param, "range");
        p.value                  = convert_double_nothrow(value_string.c_str(), dv);
        std::cout << "Parsed parameter " << p.name
                  << " of type " << p.type_name
                  << ", with range " << p.range.first
                  << " to " << p.range.second
                  << ", and initial value " << p.value << std::endl;
        info.addParameter(p);
    }
}
    
void ModelFileLoader::parseRateMatrix(const YAML::Node& rate_matrix,
                     ModelInfoFromYAMLFile& info) {
    //Assumes rate_matrix is a sequence (of rows)
    for (auto row : rate_matrix) {
        ++info.rate_matrix_rank;
        std::stringstream s;
        s << "Row " << info.rate_matrix_rank << " of rate matrix "
          << " for model " << info.model_name
          << " in " << info.model_file_path;
        std::string context = s.str();
        complainIfNot(row.IsSequence(),
                      context + " is not a sequence" );
        StrVector expr_row;
        for (auto col : row) {
            if (col.IsNull()) {
                expr_row.emplace_back("");
                //Gets whatever hasn't been assigned in this row
            }
            else if (!col.IsScalar()) {
                std::stringstream s2;
                s2 << "Column " << (expr_row.size()+1)
                   << " of " << context << " is not a scalar";
                outError(s2.str());
            } else {
                expr_row.emplace_back(col.Scalar());
            }
        }
        info.rate_matrix_expressions.emplace_back(expr_row);
    }
    dumpRateMatrixTo(info, cout);
}
        
void ModelFileLoader::parseRateMatrixExpressions(ModelInfoFromYAMLFile& info) {
    //Todo:
    //This is to execute after the rate matrix expressions are known.
    //It is to convert rate matrix expression strings into
    //a vector of expression trees, to check they are valid,
    //and then discard the expression trees.
}
    
void ModelFileLoader::parseYAMLSubstitutionModel(const YAML::Node& substitution_model,
                                const std::string& name_of_model,
                                ModelInfoFromYAMLFile& info) {
    info.model_file_path = file_path;
    info.model_name      = name_of_model;
    info.citation        = stringScalar(substitution_model,  "citation");
    info.reversible      = booleanScalar(substitution_model, "reversible");
    //
    //Todo: read off the datatype (if there is one).
    //
    info.data_type_name  = stringScalar(substitution_model,  "datatype");
    
    //
    //Todo: extract other information from the subsstitution model.
    //      Such as parameters and rate matrices and so forth
    //
    auto params = substitution_model["parameters"];
    if (params) {
        complainIfNot(params.IsSequence(),
                      "Parameters of model " + model_name +
                      " in file " + file_path + " not a sequence");
        parseYAMLModelParameters(params, info);
    }
    
    //
    //Todo: It should be okay for a model *not* to specify a rateMatrix,
    //      so long as it subclasses another model that *does* specify one.
    //
    auto rateMatrix = substitution_model["rateMatrix"];
    complainIfNot(rateMatrix, "Model " + model_name +
                  " in file " + file_path +
                  " does not specify a rateMatrix" );
    parseRateMatrix(rateMatrix, info);
    
    //
    //Todo: How do you specify paramters
    //
    auto stateFrequency = substitution_model["stateFrequency"];
    if (stateFrequency) {
        //Check that dimension of the specified parameter is the
        //same as the rank of the rate matrix (it must be!).
    } else {
        //If we have parameters with a type of frequency, we're all good.
    }
}
