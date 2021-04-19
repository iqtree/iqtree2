//
// modelfileloader.cpp
// Created by James Barbetti on 01-Apr-2021.
//

#include "modelfileloader.h"
#include "modelexpression.h" //for ModelExpression::ModelException
#include <utils/stringfunctions.h>
   
ModelFileLoader::ModelFileLoader(const char* path): file_path(path) {
}
        
std::string ModelFileLoader::stringScalar(const YAML::Node& node,
                                          const char* key,
                                          const char* default_value) {
    auto   scalar_node = node[key];
    return scalar_node ? scalar_node.Scalar() : default_value;
}
    
bool ModelFileLoader::booleanScalar(const YAML::Node& node,
                                    const char* key,
                                    const bool default_value) {
    std::string s = string_to_lower(stringScalar(node, key, ""));
    if (s.empty()) {
        return default_value;
    }
    return s == "true" || s == "yes" || s == "t" || s == "y" || s == "1";
}

int ModelFileLoader::integerScalar(const YAML::Node& node,
                                   const char* key,
                                   const int default_value) {
    std::string s = stringScalar(node, key, "");
    if (s.empty()) {
        return default_value;
    }
    if (s[0]<'0' || '9'<s[0]) {
        return default_value;
    }
    return atoi(s.c_str());
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
    
ModelParameterRange ModelFileLoader::parseRange(const YAML::Node& node,
                                                const char* key,
                                                const ModelParameterRange& default_value ) {
    auto r = node[key];
    if (!r) {
        return default_value;
    }
    ModelParameterRange range;
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
    for (const YAML::Node& param: params) {
        const YAML::Node& name_node = param["name"];
        if (name_node) {
            if (name_node.IsScalar()) {
                parseModelParameter(param, name_node.Scalar(), info);
            }
            else if (name_node.IsSequence()) {
                for (auto current_name: name_node) {
                    if (current_name.IsScalar()) {
                        parseModelParameter(param, current_name.Scalar(), info);
                    }
                }
            }
            else {
                outError("Model parameter must have a name");
            }
        }
    }
}

void ModelFileLoader::parseModelParameter(const YAML::Node& param,
                                          std::string name,
                                          ModelInfoFromYAMLFile& info) {
    YAMLFileParameter p;
    p.name       = name;
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
            const char* msg = "Subscript range does not end with right parenthesis";
            throw ModelExpression::ModelException(msg);
        }
        p.name = p.name.substr(0, bracket);
    }
    
    bool overriding = false;
    for (const YAMLFileParameter& oldp: info.parameters) {
        if (oldp.name == p.name) {
            complainIfNot(oldp.is_subscripted == p.is_subscripted,
                          "Canot redefine subscripted parameter"
                          " as unsubscripted (or vice versa)");
            complainIfNot(oldp.minimum_subscript == p.minimum_subscript,
                          "Cannot redefine parameter subscript range");
            complainIfNot(oldp.maximum_subscript == p.maximum_subscript,
                          "Cannot redefine parameter subscript range");
            p = oldp;
            overriding = true;
            break;
        }
    }
    
    p.type_name = string_to_lower(stringScalar(param, "type", p.type_name.c_str()));
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
    //Todo: What if name was a list, and initValue is also a list?!
    std::string value_string = stringScalar(param, "initValue", "");
    p.range                  = parseRange  (param, "range", p.range);
    if (value_string!="") {
        p.value = convert_double_nothrow(value_string.c_str(), dv);
    } else if (!overriding) {
        p.value = dv;
    }
    p.description = stringScalar(param, "description", p.description.c_str());
    std::cout << "Parsed parameter " << p.name
              << " of type " << p.type_name
              << ", with range " << p.range.first
              << " to " << p.range.second
              << ", and initial value " << p.value << std::endl;
    info.addParameter(p);
}

void ModelFileLoader::parseYAMLModelConstraints(const YAML::Node& constraints,
                                                ModelInfoFromYAMLFile& info) {
    for (const YAML::Node& constraint: constraints) {
        //constraints are assignments of the form: name = value
        //and are equivalent to parameter name/initialValue pairs
        //Todo: For now, I don't want to support (x,y) = (1,2).
        //
        std::stringstream complaint;
        if (!constraint.IsScalar()) {
            complaint << "Constraint setting"
                      << " for model " << info.model_name
                      << " was not a scalar.";
            outError(complaint.str());
        }
        std::string constraint_string = constraint.Scalar();
        ModelExpression::InterpretedExpression interpreter(info, constraint_string);
        ModelExpression::Expression* x = interpreter.expression();
        
        if (!x->isAssignment()) {
            complaint << "Constraint setting for model " << info.model_name
            << " was not an asignment: " << constraint_string;
            outError(complaint.str());
        }
        ModelExpression::Assignment* a = dynamic_cast<ModelExpression::Assignment*>(x);
        
        if (!a->getTarget()->isVariable()) {
            delete x;
            complaint << "Constraint setting for model " << info.model_name
                      << " did not assign a variable: " << constraint_string;
            outError(complaint.str());
        }
        ModelExpression::Variable* v = a->getTargetVariable();
        double setting = a->getExpression()->evaluate();
        ModelVariable& mv = info.assign(v->getName(), setting);
        mv.markAsFixed();
        std::cout << "Assigned " << v->getName()
                  << " := " << setting << std::endl;
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
                                                 ModelInfoFromYAMLFile& info,
                                                 ModelListFromYAMLFile& list) {
    
    std::string parent_model = stringScalar(substitution_model, "frommodel", "");
    if (parent_model != "") {
        if (list.hasModel(parent_model)) {
            info = list.getModel(parent_model);
        } else {
            std::stringstream complaint;
            complaint << "Model " << name_of_model << " specifies frommodel "
                      << parent_model << ", but that model was not found.";
            outError(complaint.str());
        }
    }
    
    info.model_file_path = file_path;
    info.model_name      = name_of_model;
    info.citation        = stringScalar(substitution_model,  "citation",   info.citation.c_str());
    info.DOI             = stringScalar(substitution_model,  "doi",        info.DOI.c_str());
    info.reversible      = booleanScalar(substitution_model, "reversible", info.reversible);
    info.data_type_name  = stringScalar(substitution_model,  "datatype",   info.data_type_name.c_str());
    //Note: doco currently says this will be called "forData".
    //
    //Todo: read off the in-lined datatype (if there is one).
    //
    
    info.num_states      = integerScalar(substitution_model, "numStates", 0);
    if (info.num_states==0) {
        info.num_states = 4;
    }
    
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
    
    auto constraints = substitution_model["constraints"];
    if (constraints) {
        complainIfNot(constraints.IsSequence(),
                      "Constraints for model " + model_name +
                      " in file " + file_path + " not a sequence");
        parseYAMLModelConstraints(constraints, info);
    }
    
    auto rateMatrix = substitution_model["rateMatrix"];
    if (info.rate_matrix_expressions.empty()) {
        complainIfNot(rateMatrix, "Model " + model_name +
                      " in file " + file_path +
                      " does not specify a rateMatrix" );
    }
    
    //If this model subclasses another it doesn't have to specify
    //a rate matrix (if it doesn't it inherits from its parent model).
    if (rateMatrix) {
        parseRateMatrix(rateMatrix, info);
    }
    
    auto stateFrequency = substitution_model["stateFrequency"];
    if (stateFrequency) {
        //
        //Check that dimension of the specified parameter is the
        //same as the rank of the rate matrix (it must be!).
        //
        std::string freq     = stateFrequency.IsScalar() ? stateFrequency.Scalar() : "";
        std::string low_freq = string_to_lower(freq);
        if (low_freq=="estimate") {
            info.frequency_type = StateFreqType::FREQ_ESTIMATE;
        } else if (low_freq=="empirical") {
            info.frequency_type = StateFreqType::FREQ_EMPIRICAL;
        } else if (low_freq=="uniform") {
            info.frequency_type = StateFreqType::FREQ_EQUAL;
        } else if (info.isFrequencyParameter(low_freq)) {
            info.frequency_type = StateFreqType::FREQ_USER_DEFINED;
        }
    } else {
        //If we have parameters with a type of frequency, we're all good.
        //If we don't, then what?   We might have inherited from a parent
        //model, too.  That'd be okay.
    }
}
