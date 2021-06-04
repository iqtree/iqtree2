//
// modelfileloader.cpp
// Created by James Barbetti on 01-Apr-2021.
//

#include "modelfileloader.h"
#include "modelexpression.h" //for ModelExpression::ModelException
#include <utils/stringfunctions.h>
#include <tree/phylotree.h> //for TREE_LOG_LINE macro
   
typedef ModelExpression::InterpretedExpression Interpreter;

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

double ModelFileLoader::doubleScalar(const YAML::Node& node,
                                  const char* key, 
                                  const double default_value) {
    std::string s = stringScalar(node, key, "");
    if (s.empty()) {
        return default_value;
    }
    char ch = s[0];
    if (ch<'0' || '9'<ch) {
        if (ch!='-' && ch!='+' && ch!='.') {
            return default_value;
        }
    }
    double value = convert_double_nothrow(s.c_str(), default_value);
    return value;
}

void ModelFileLoader::complainIfSo(bool        check_me,
                                   std::string error_message) {
    if (check_me) {
        outError(error_message);
    }
}

void ModelFileLoader::complainIfNot(bool check_me,
                          std::string error_message) {
    if (!check_me) {
        outError(error_message);
    }
}
    
double ModelFileLoader::toDouble(const YAML::Node& i, double default_val,
                                 ModelInfoFromYAMLFile& info, std::string context) {
    if (!i.IsScalar()) {
        return default_val;
    }
    std::string expr = i.Scalar();
    return info.evaluateExpression(expr, context);
}
    
ModelParameterRange ModelFileLoader::parseRange(const YAML::Node& node,
                                                const char* key,
                                                const ModelParameterRange& default_value,
                                                ModelInfoFromYAMLFile& info) {
    auto r = node[key];
    if (!r) {
        return default_value;
    }
    ModelParameterRange range;
    if (r.IsSequence()) {
        std::string context = std::string(key) + " parameter";
        int ix = 0;
        for (auto i : r) {
            if (ix==0) {
                range.first  = toDouble(i, 0, info, "lower bound for " + context);
            } else if (ix==1) {
                range.second = toDouble(i, 0, info, "upper bound for " + context);
            } else {
                throw ModelExpression::ModelException
                      ("Range may only have two bounds (lower, upper)");
            }
            ++ix;
        }
        range.is_set = ( 1 <= ix );
        if ( ix == 1 ) {
            range.second = range.first;
        } else if (range.second < range.first) {
            std::stringstream complaint;
            complaint << "Range has lower bound (" << range.first << ")"
                      << " greater than its upper bound (" << range.second << ")";
            throw ModelExpression::ModelException(complaint.str());
        }
    } else {
        
    }
    return range;
}
    
void ModelFileLoader::parseYAMLModelParameters(const YAML::Node& params,
                                               ModelInfoFromYAMLFile& info,
                                               PhyloTree* report_to_tree) {
    //
    //Assumes: params is a sequence of parameter declarations
    //
    for (const YAML::Node& param: params) {
        const YAML::Node& name_node = param["name"];
        if (name_node) {
            if (name_node.IsScalar()) {
                parseModelParameter(param, name_node.Scalar(), info,
                                    report_to_tree);
            }
            else if (name_node.IsSequence()) {
                for (auto current_name: name_node) {
                    if (current_name.IsScalar()) {
                        std::string name = current_name.Scalar();
                        parseModelParameter(param, name, info,
                                            report_to_tree);
                    }
                }
            }
            else {
                outError("Model parameter must have a name");
            }
        }
    }
}

void ModelFileLoader::setParameterSubscriptRange(ModelInfoFromYAMLFile& info,
                                                 YAMLFileParameter& p) {
    const char* parsing_what = "subscript expression";
    try {
        Interpreter interpreter(info, p.subscript_expression );
        auto x = interpreter.expression();
        if (x->isRange()) {
            auto r = dynamic_cast<ModelExpression::RangeOperator*>(x);
            parsing_what = "minimum subscript";
            p.minimum_subscript = r->getIntegerMinimum();
            parsing_what = "maximum subscript";
            p.maximum_subscript = r->getIntegerMaximum();
        } else {
            p.minimum_subscript = 1;
            p.maximum_subscript = x->evaluateAsInteger();
        }
        if (p.maximum_subscript<p.minimum_subscript) {
            std::stringstream msg;
            msg << "Illegal subscript expression"
                << " for " << info.getName()
                << " parameter " << p.name << ":\n"
                << " minimum subscript (" << p.minimum_subscript << ")"
                << " was greater than maximum (" << p.maximum_subscript << ")";
        }
    }
    catch (ModelExpression::ModelException x) {
        std::stringstream msg;
        msg << "Error parsing " << parsing_what
            << " for " << info.getName()
            << " parameter " << p.name << ":\n"
            << x.getMessage();
        outError(msg.str());
    }
}

void ModelFileLoader::setParameterType(const YAML::Node&      param, 
                                       bool                   overriding,
                                       ModelInfoFromYAMLFile& info,
                                       YAMLFileParameter&     p) {
    if (p.type_name=="calculated rate") {
        p.type      = ModelParameterType::RATE;
    } else if (p.type_name=="rate") {
        p.type = ModelParameterType::RATE;
    } else if (p.type_name=="frequency") {
        p.type = ModelParameterType::FREQUENCY;
    } else if (p.type_name=="proportion") {
        p.type = ModelParameterType::PROPORTION;
    } else if (p.type_name=="shape") {
        p.type = ModelParameterType::SHAPE;
    } else if (p.type_name=="weight") {
        p.type = ModelParameterType::WEIGHT;
    } else {
        p.type_name = "other";
        p.type = ModelParameterType::OTHER;
    }
}

void ModelFileLoader::setParameterTolerance(const YAML::Node&      param, 
                                            bool                   overriding,
                                            ModelInfoFromYAMLFile& info,
                                            YAMLFileParameter&     p) {
    if (!overriding) {
        switch (p.type) {
            case ModelParameterType::FREQUENCY:
                p.tolerance = 1e-4;
                break;
            case ModelParameterType::PROPORTION:
                p.tolerance = 1e-4;
                break;
            case ModelParameterType::RATE:
                p.tolerance = 1e-6; //though FreeRate 
                break;
            case ModelParameterType::SHAPE:
                p.tolerance = 0.001;
                break;
            case ModelParameterType::WEIGHT:
                p.tolerance = 1e-4;
                break;
            case ModelParameterType::OTHER:
                p.tolerance = 0;
                break;
        }
    }
    p.tolerance = doubleScalar(param, "tolerance", p.tolerance);
}

void ModelFileLoader::setParameterValue(const YAML::Node&      param, 
                                        bool                   overriding,
                                        ModelInfoFromYAMLFile& info,
                                        YAMLFileParameter&     p) {
    auto count = p.maximum_subscript - p.minimum_subscript + 1;
    ASSERT(0<count);
    double dv   = 0.0; //default initial value
    if (p.type==ModelParameterType::FREQUENCY) {
        dv = 1.0 / static_cast<double>(count);
    } else if (p.type==ModelParameterType::PROPORTION) {
        dv = 1.0;
    } else if (p.type==ModelParameterType::RATE) {
        dv = 1.0;
    } else if (p.type==ModelParameterType::SHAPE) {
        dv = 1.0;
    } else if (p.type==ModelParameterType::WEIGHT) {
        dv = 1.0 / static_cast<double>(count);    
    }
    //Todo: What if name was a list, and initValue is also a list?!
    p.init_expression = stringScalar(param, "initValue", "");
    p.range            = parseRange  (param, "range", p.range, info);
    if (p.init_expression!="") {
        int dummy_subscript = p.is_subscripted ? p.minimum_subscript : -1;
        info.forceAssign("subscript", dummy_subscript);
        p.value = info.evaluateExpression(p.init_expression,
                                          "initial value of " + p.name +
                                          " parameter" );
    } else if (!overriding) {
        if (p.range.is_set) {
            //Initial value has to be inside the legal range.
            if (dv<p.range.first) {
                dv = p.range.first;
            } if (p.range.second<dv) {
                dv = p.range.second;
            }
        }
        p.value = dv;
    }
}

void ModelFileLoader::parseModelParameter(const YAML::Node& param,
                                          std::string name,
                                          ModelInfoFromYAMLFile& info,
                                          PhyloTree* report_to_tree) {
    YAMLFileParameter p;
    p.name       = name;
    auto name_len = p.name.length();
    auto bracket  = p.name.find('(');
    if ( bracket != std::string::npos ) {
        auto expr_len = name_len-bracket;
        p.is_subscripted = true;
        p.subscript_expression = p.name.substr(bracket, expr_len );
        p.name = p.name.substr(0, bracket);
        setParameterSubscriptRange(info, p);
    } else {
        p.is_subscripted    = false;
        p.minimum_subscript = 0;
        p.maximum_subscript = 1;
    }
    
    p.type_name = string_to_lower(stringScalar(param, "type", p.type_name.c_str()));
    if (p.type_name=="matrix") {
        complainIfSo(p.is_subscripted,
                     "Matrix subscripts are implied by the matrix value itself, "
                     " but " + p.name + " parameter of model " + info.getName() +
                     " was explicitly subscripted (which is not supported).");
        auto value   = param["value"];
        auto formula = param["formula"];
        auto rank    = param["rank"];
        complainIfNot(value || (formula && rank), p.name +
                      " matrix parameter's value must be defined"
                      " in model " + info.getName() + ".");
        parseMatrixParameter(param, name, info, report_to_tree);
        return;
    }
    
    bool overriding = isAParameterOverride(info, p);

    setParameterType     (param, overriding, info, p);
    setParameterTolerance(param, overriding, info, p);
    setParameterValue    (param, overriding, info, p);
    
    p.description = stringScalar(param, "description", p.description.c_str());

    p.logParameterState("Parsed", report_to_tree);

    info.addParameter(p);
}

bool ModelFileLoader::isAParameterOverride(ModelInfoFromYAMLFile& info,
                                           YAMLFileParameter& p) {
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
            return true;
            break;
        }
    }
    return false;
}


void ModelFileLoader::parseMatrixParameter(const YAML::Node& param,
                                           std::string name,
                                           ModelInfoFromYAMLFile& info,
                                           PhyloTree* report_to_tree) {
    //Assumed: parameter has a "value" entry, and it is a sequence of sequences
    
    auto         value        = param["value"];
    int          column_count = 0;
    StringMatrix expressions;
    
    auto         formula_node = param["formula"];
    std::string  formula;

    auto         rank_node    = param["rank"];
    int          rank         = 0;

    if (value) {
        complainIfNot(value.IsSequence(),
                      "value of " + name +
                      " matrix of model " + info.getName() +
                      " was not a matrix");
        for (auto row : value) {
            ++rank;
            std::stringstream s;
            s << "Row " << rank << " of " + name + " matrix "
              << " for model " << info.model_name
              << " in " << info.model_file_path;
            std::string context = s.str();
            complainIfNot(row.IsSequence(),
                          context + " is not a sequence" );
            StrVector expr_row;
            for (auto col : row) {
                if (col.IsNull()) {
                    expr_row.emplace_back("");
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
            expressions.emplace_back(expr_row);
            column_count = ( expr_row.size() <= column_count )
                         ? column_count : static_cast<int>(expr_row.size());
        }
        if (column_count<expressions.size()) {
            column_count = static_cast<int>(expressions.size());
        }
    }
    else if (!formula_node || !rank_node) {
        outError(name +
                 " matrix of model " + info.getName() +
                 " had no value, and lacked either a rank or a formula");
    }
    if (rank_node) {
        complainIfNot(rank_node.IsScalar(), "rank of " + name +
                      " matrix of model " + info.getName() +
                      " was not a scalar");
        std::string rank_str = rank_node.Scalar();
        
        Interpreter interpreter(info, rank_str);
        double rank_dbl = interpreter.evaluate();
        rank            = (int)floor(rank_dbl);
        complainIfNot(0<rank, "rank of " + name + " matrix of model "
                      + info.getName() + " was invalid (" + rank_str + ")");
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "Rank of " << info.getName() << "." << name <<
                      " was " << rank_str << " ... or " << rank);
    }
    if (formula_node) {
        complainIfNot(formula_node.IsScalar(), "formula of " + name +
                      " matrix of model " + info.getName() +
                      " was not a scalar");
        formula = formula_node.Scalar();
    }
    
    std::string lower_name = string_to_lower(name);
    if (lower_name=="ratematrix") {
        expressions.makeSquare(true);
        info.rate_matrix_rank           = rank;
        info.rate_matrix_expressions    = expressions;
        info.rate_matrix_formula        = formula;
    }
    else if (lower_name=="tiplikelihood") {
        expressions.makeRectangular(column_count);
        info.tip_likelihood_rank        = rank;
        info.tip_likelihood_expressions = expressions;
        info.tip_likelihood_formula     = formula;
    }
    else {
        outError(name + " matrix parameter not recognized"
                 " in " + info.getName() + " model");
    }
    if (YAMLMatrixVerbosity <= verbose_mode) {
        std::stringstream matrix_stream;
        dumpMatrixTo(lower_name.c_str(), info, expressions, rank,
                     formula, matrix_stream);
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, matrix_stream.str());
    }
}

YAMLFileParameter
    ModelFileLoader::addDummyFrequencyParameterTo(ModelInfoFromYAMLFile& info,
                                                  PhyloTree* report_to_tree) {
    YAMLFileParameter   p;
    p.name              = "freq";
    p.is_subscripted    = true;
    p.minimum_subscript = 1;
    p.maximum_subscript = info.getNumStates();
    p.type_name         = "frequency";
    p.type              = ModelParameterType::FREQUENCY;
    auto num_states     = info.getNumStates();
    ASSERT(0 < num_states);
    p.value            = 1 / static_cast<double>(num_states);
    info.addParameter(p);
    return p;
}

void ModelFileLoader::parseYAMLMixtureModels(const YAML::Node& mixture_models,
                                             ModelInfoFromYAMLFile& info,
                                             ModelListFromYAMLFile& list,
                                             PhyloTree* report_to_tree) {
    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, "Processing mixtures" );
    info.mixed_models = new MapOfModels();
    for (const YAML::Node& model: mixture_models) {
        std::string child_model_name = stringScalar(model, "substitutionmodel", "");
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, "Processing mixture model" );
        ModelInfoFromYAMLFile child_info;
        parseYAMLModel(model, child_model_name, child_info,
                       list, &info, report_to_tree);
        (*info.mixed_models)[info.getName()] = child_info;
    }
}

void ModelFileLoader::parseYAMLModelConstraints(const YAML::Node& constraints,
                                                ModelInfoFromYAMLFile& info,
                                                PhyloTree* report_to_tree) {
    int constraint_num = 1;
    for (const YAML::Node& constraint: constraints) {
        //
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
        try { 
            ModelExpression::InterpretedExpression
                interpreter(info, constraint_string);
            ModelExpression::Expression* x = interpreter.expression();
            if (!x->isAssignment()) {
                complaint << "Constraint setting for model " << info.model_name
                          << " was not an asignment: " << constraint_string;
                outError(complaint.str());
            }
            ModelExpression::Assignment* a =
                dynamic_cast<ModelExpression::Assignment*>(x);
            setConstraint(a, info, constraint_string, report_to_tree);
        }
        catch (ModelExpression::ModelException x) {
            std::stringstream complaint;
            complaint << "An error occurred parsing constraint " << constraint_num
                      << ": " << x.getMessage();
            outError(complaint.str());
        }
        ++constraint_num;
    }
}

double ModelFileLoader::setConstraint(ModelExpression::Assignment* a,
                                      ModelInfoFromYAMLFile& info,
                                      const std::string& constraint_string,
                                      PhyloTree* report_to_tree) {
    ModelExpression::Expression* assigned = a->getTarget();
    if (!assigned->isVariable()) {
        std::stringstream complaint;
        complaint << "Constraint setting for model " << info.model_name
                  << " did not assign a variable: "  << constraint_string;
        outError(complaint.str());
    }
    ModelExpression::Variable*   v = a->getTargetVariable();
    ModelExpression::Expression* x = a->getExpression();
    double setting;
    if (x->isAssignment()) {
        auto a2 = dynamic_cast<ModelExpression::Assignment*>(x);
        setting = setConstraint(a2, info, constraint_string,
                                report_to_tree);
    } else {
        setting = a->getExpression()->evaluate();
    }
    ModelVariable& mv = info.assign(v->getName(), setting);
    mv.markAsFixed();
    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                  "Assigned " << v->getName()
                  << " := " << setting);

    //And now for the horrible bit. Variable v might have been referenced
    //in a range expression, for a parameter.
    for (auto it = info.parameters.begin(); it!=info.parameters.end(); ++it) {
        YAMLFileParameter& param = *it;
        if (param.is_subscripted && !param.subscript_expression.empty()) {
            Interpreter interpreter(info, param.subscript_expression);
            std::pair<int,int> range;
            bool rangeChanged = false;
            try {
                if (interpreter.evaluateIntegerRange(range)) {
                    if (range.first != param.minimum_subscript ||
                        range.second != param.maximum_subscript) {
                        rangeChanged = true;
                    }
                }
            }
            catch (ModelExpression::ModelException x) {
                throw ModelExpression::ModelException
                      ("Error evaluating subscript range for " + param.name + 
                       ": " + x.getMessage() );
            }
            if (rangeChanged) {
                info.updateParameterSubscriptRange(param, range.first, 
                                                   range.second, report_to_tree);
            }
        }
    }

    return setting;
}

void ModelFileLoader::parseRateMatrix(const YAML::Node& rate_matrix,
                                      ModelInfoFromYAMLFile& info,
                                      PhyloTree* report_to_tree) {
    //Assumes rate_matrix is a sequence (of rows)
    size_t column_count = 0;
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
        column_count = ( expr_row.size() < column_count )
                     ? column_count : expr_row.size();
    }
    size_t row_count = info.rate_matrix_expressions.size();
    if ( row_count != column_count) {
        std::stringstream s2;
        s2 << "Rate matrix "
           << " for model " << info.model_name
           << " in " << info.model_file_path
           << " was not square: it had " << row_count << " rows"
           << " and " << column_count << " columns.";
        outWarning(s2.str());
    }
    info.rate_matrix_expressions.makeSquare(true);
    
    //
    //Todo: Are off-diagonal entries in the matrix allowed
    //to be blank?
    //
    if (YAMLMatrixVerbosity <= verbose_mode) {
        std::stringstream matrix_stream;
        dumpMatrixTo("rate", info, info.rate_matrix_expressions,
                     info.rate_matrix_rank, info.rate_matrix_formula,
                     matrix_stream);
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, matrix_stream.str());
    }
}

void ModelFileLoader::dumpMatrixTo(const char* name, ModelInfoFromYAMLFile& info,
                                   const StringMatrix& matrix, int rank,
                                   const std::string& formula, std::stringstream &out) {
    ModelVariable& row_var    = info.forceAssign("row",    0);
    ModelVariable& column_var = info.forceAssign("column", 0);

    std::string       with_formula;
    std::stringstream dump;
    if (!matrix.empty()) {
        //If there's a matrix of expressions, dump the expressions
        int row = 0;
        for (auto r : matrix) {
            row_var.setValue(row+1);
            const char* separator = "";
            int col = 0;
            for (auto c: r) {
                column_var.setValue(col+1);
                dump << separator << c;
                separator = " : ";
                ++col;
            }
            dump << "\n";
            ++row;
        }
    }
    else {
        //If there is a formula, dump the formula, and its
        //current value, for each entry in the matrix
        with_formula = " (with formula " + formula + ")";
        bool broken = false;
        for (int row=0; row<rank && !broken; ++row) {
            row_var.setValue(row+1);
            const char* separator = "";
            for (int col=0; col<rank && !broken; ++col) {
                column_var.setValue(col+1);
                try {
                    Interpreter interpreter(info, formula);
                    double value = interpreter.evaluate();
                    dump << separator << value;
                }
                catch (ModelExpression::ModelException x) {
                    dump << separator << " ERROR(" << x.getMessage() << ")"
                         << " for entry [" << (row+1) << "," << (col+1) << "]";
                    broken = true;
                }
                separator = " : ";
            }
            dump << "\n";
        }
    }
    out << name << " matrix for " << info.model_name
        << with_formula
        << " is...\n" << dump.str();
}

void ModelFileLoader::handleInheritance(ModelInfoFromYAMLFile& info,
                                        ModelListFromYAMLFile& list,
                                        PhyloTree* report_to_tree) {
    bool have_first_parent = false;
    for (string ancestral_model : split_string(info.superclass_model_name, "+")) {
        if (list.hasModel(ancestral_model)) {
            const ModelInfoFromYAMLFile& ancestor =
                    list.getModel(ancestral_model);
            if (info.is_rate_model && !ancestor.is_rate_model) {
                std::stringstream complaint;
                complaint << "Rate model " << info.model_name
                            << " cannot be based on" 
                            << "a substitution model " 
                            << ancestral_model << ".";
                outError(complaint.str());
            }
            if (!have_first_parent) {
                std::string save_name = info.model_name;
                info = list.getModel(ancestral_model);
                info.model_name = save_name;
                TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                            "Model " << info.model_name
                            << " is based on model " << ancestral_model);
                have_first_parent = true;
            } else {
                info.inheritModel(ancestor);
                TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                            "Model " << info.model_name
                            << " is also based on" 
                            << " model " << ancestral_model);
            }
        } else {
            std::stringstream complaint;
            complaint << "Model " << info.model_name 
                        << " specifies a superclass model of "
                        << info.superclass_model_name << ","
                        << " but that model was not found.";
            if (info.is_rate_model) {
                complaint << "\nRecognized rate models are: ";
                complaint << list.getListOfRateModelNames() << ".";
            } else {
                complaint << "\nRecognized substitution models are: ";
                complaint << list.getListOfSubstitutionModelNames() << ".";
            }
            outError(complaint.str());
        }
    }
}

void ModelFileLoader::setModelStateFrequency(const YAML::Node& substitution_model,
                                             ModelInfoFromYAMLFile& info,
                                             PhyloTree* report_to_tree) {
    auto stateFrequency = substitution_model["stateFrequency"];
    std::string low_freq;
    if (stateFrequency) {
        complainIfSo(info.is_rate_model,
                     "Cannot specify state frequencies"
                     " for rate model " + model_name +
                      " in file " + file_path);
        //
        //Check that dimension of the specified parameter is the
        //same as the rank of the rate matrix (it must be!).
        //
        std::string freq     = stateFrequency.IsScalar() ? stateFrequency.Scalar() : "";
        low_freq = string_to_lower(freq);
    }

    if (stateFrequency && stateFrequency.IsSequence()) {
        info.frequency_type = StateFreqType::FREQ_USER_DEFINED;
        YAMLFileParameter freq_param = addDummyFrequencyParameterTo(info, report_to_tree);
        int subscript = freq_param.minimum_subscript;
        for (auto f: stateFrequency) {
            complainIfNot(f.IsScalar(), "Model " + model_name +
                            " in file " + file_path +
                            " has unrecognized frequency ");
            complainIfNot(subscript<=freq_param.maximum_subscript,
                            "Too many frequencies specified for "
                            "Model " + model_name +
                            " in file " + file_path);
            ModelExpression::InterpretedExpression x(info, f.Scalar());
            auto   var_name  = freq_param.getSubscriptedVariableName(subscript);
            double var_value = x.evaluate();
            info.assign(var_name, var_value);
            TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                            "Assigned frequency: " << var_name
                            << " := " << var_value  );
            ++subscript;
        }
    }
    else if (low_freq.empty()) {   
        if (info.frequency_type != StateFreqType::FREQ_UNKNOWN) {
            //We inherited a frequency type from a parent model. Use that.
        }
        else if (info.hasFrequencyParameters(info.num_states)) {
            //If we have parameters with a type of frequency, we're all good.
            info.frequency_type = StateFreqType::FREQ_USER_DEFINED;
        } else {
            //If we don't, then what?  Use estimate as a last resort.
            info.frequency_type = StateFreqType::FREQ_ESTIMATE;
        }
        return;
    } else if (low_freq=="estimate") {
        info.frequency_type = StateFreqType::FREQ_ESTIMATE;
    } else if (low_freq=="empirical") {
        info.frequency_type = StateFreqType::FREQ_EMPIRICAL;
    } else if (low_freq=="uniform") {
        info.frequency_type = StateFreqType::FREQ_EQUAL;
    } else if (info.isFrequencyParameter(low_freq)) {
        info.frequency_type = StateFreqType::FREQ_USER_DEFINED;
    } 
    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                    "After setting frequency type"
                    " of " << info.model_name <<
                    " to " << low_freq << "...");
    info.logVariablesTo(*report_to_tree);
}


void ModelFileLoader::parseYAMLModel(const YAML::Node& substitution_model,
                                     const std::string& name_of_model,
                                     ModelInfoFromYAMLFile& info,
                                     ModelListFromYAMLFile& list,
                                     ModelInfoFromYAMLFile* parent_model,
                                     PhyloTree* report_to_tree) {
    info.parent_model          = parent_model;
    info.is_modifier_model     = false;
    info.superclass_model_name = stringScalar(substitution_model, "frommodel", "");
    info.rate_distribution     = stringScalar(substitution_model, "ratedistribution", "");
    info.model_name            = name_of_model;
                            
    if (!info.superclass_model_name.empty()) {
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity, 
                      name_of_model << " parent is " 
                      << info.superclass_model_name << "." );
        if (string_to_upper(info.superclass_model_name)=="ANY") {
            info.superclass_model_name.clear();
            info.is_modifier_model = true;
        }
        else {
            if (info.model_name.empty()) {            
                info.model_name = info.superclass_model_name;
            }
            handleInheritance(info, list, report_to_tree);
        }
    }
    info.model_file_path = file_path;
    info.citation        = stringScalar (substitution_model, "citation",    info.citation.c_str());
    info.DOI             = stringScalar (substitution_model, "doi",         info.DOI.c_str());
    info.url             = stringScalar (substitution_model, "doi",         info.url.c_str());
    info.reversible      = booleanScalar(substitution_model, "reversible",  info.reversible);
    info.opt_algorithm   = stringScalar (substitution_model, "optimization", info.opt_algorithm.c_str());
    info.data_type_name  = stringScalar (substitution_model, "datatype",    info.data_type_name.c_str());
    //Note: doco currently says this will be called "forData".
        
    if (!info.is_rate_model) {
        int num_states_requested = integerScalar(substitution_model,
                                                "numStates", info.num_states);
        info.setNumberOfStatesAndSequenceType(num_states_requested, report_to_tree);
    } else {
        info.forceAssign("categories", 1);   
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
        parseYAMLModelParameters(params, info, report_to_tree);
    }

    parseYAMLSubModels(substitution_model, info, list, report_to_tree);
    
    auto constraints = substitution_model["constraints"];
    if (constraints) {
        complainIfNot(constraints.IsSequence(),
                      "Constraints for model " + model_name +
                      " in file " + file_path + " not a sequence");
        parseYAMLModelConstraints(constraints, info, report_to_tree);
    }
    
    //
    //Note: if rateMatrix was read as a parameter,
    //      rate_matrix_expressions will aready have been set.
    //
    auto rateMatrix = substitution_model["rateMatrix"];
    if (info.rate_matrix_expressions.empty() && 
        info.rate_matrix_formula.empty() && !info.isMixtureModel() &&
        !info.is_modifier_model && info.linked_models==nullptr && 
        !info.is_rate_model) {
        complainIfNot(rateMatrix, "Model " + info.model_name +
                      " in file " + file_path +
                      " does not specify a rateMatrix" );
    }
    
    //If this model subclasses another it doesn't have to specify
    //a rate matrix (if it doesn't it inherits from its parent model).
    if (rateMatrix) {
        complainIfSo(info.is_rate_model,
                     "Cannot specify rate matrix"
                     " for rate model " + model_name +
                      " in file " + file_path);
        parseRateMatrix(rateMatrix, info, report_to_tree);
    }
    
    setModelStateFrequency(substitution_model, info, report_to_tree);

    const char* recognized_string_property_names[] = {
        "errormodel", //One of "+E", "+EA", "+EC", "+EG", "+ET"
                      //recognized by ModelDNAError.
    };

    for (const char* prop_name : recognized_string_property_names ) {
        auto prop_node = substitution_model[prop_name];
        if (prop_node) {
            if (prop_node.IsScalar()) {
                std::string prop_value = prop_node.Scalar();
                info.string_properties[prop_name] = prop_value;
                TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                              "string property " << prop_name <<
                              " set to " << prop_value);
            }
            //Todo: what about lists?
        }
    }

    parseYAMLModelWeightAndScale(substitution_model, info, report_to_tree );


}

void ModelFileLoader::parseYAMLSubModels(const YAML::Node& substitution_model,
                                         ModelInfoFromYAMLFile& info,
                                         ModelListFromYAMLFile& list,
                                         PhyloTree* report_to_tree) {
    //Mixtures have to be handled before constraints, as constraints
    //that are setting parameters in mixed models... would otherwise
    //not be resolved correctly.
    auto mixtures = substitution_model["model_mixture"];
    if (mixtures) {
        complainIfSo(info.is_rate_model,
                     "mixed rate models not supported."
                     " cannot parse model " + model_name +
                      " in file " + file_path);
        complainIfNot(mixtures.IsSequence(),
                      "model_mixture for model " + model_name +
                      " in file " + file_path + " not a sequence");
        parseYAMLMixtureModels(mixtures, info, list, report_to_tree);
    }

    auto linked = substitution_model["linked_models"];
    if (linked) {
        complainIfSo(info.is_rate_model,
                     "linked substitution models"
                     " are not supported for rate models."
                     " cannot parse model " + model_name +
                      " in file " + file_path);
        complainIfNot(linked.IsSequence(),
                      "linked_models for model " + model_name +
                      " in file " + file_path + " not a sequence");
        //Todo: parseYAMLLinkedModels(linked, info, list, report_to_tree);
    }
}

void ModelFileLoader::parseYAMLModelWeightAndScale
        (const YAML::Node& substitution_model,
        ModelInfoFromYAMLFile& info,
        PhyloTree* report_to_tree) {
    auto weight = substitution_model["weight"];
    if (weight) {
        complainIfNot(info.parent_model!=nullptr,
                    "Model " + model_name +
                    " in file " + file_path +
                    " is not part of a mixture model" );
        complainIfNot(weight.IsScalar(),
                    "Weight for model " + model_name +
                    " in mixture model " + info.parent_model->getName() +
                    " in file " + file_path +
                    " was not a scalar" );
        info.weight_formula = weight.Scalar();
        info.model_weight = info.getModelWeight();
        bool is_fixed = info.isModelWeightFixed();
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      (is_fixed ? "Fixed" : "Variable") <<
                      " weight of " << info.parent_model->getName() <<
                      "." << model_name << " was " << info.weight_formula <<
                      "=" << info.model_weight);
    }

    auto scale  = substitution_model["scale"];
    if (scale) {
        //Todo: can you legitimately set the scale
        //      for a substitution model that is not part
        //      of a mixture model.
        //
        complainIfNot(info.parent_model!=nullptr,
                      "Model " + model_name +
                      " in file " + file_path +
                      " is not part of a mixture model" );
        //Todo: set the scale.
    }

    if (info.parent_model!=nullptr) {
        if (!scale) {
            //Default the scale
        }
        complainIfNot(weight,
                      "No weight specified"
                      " for model " + model_name +
                      " in mixture " + info.parent_model->getName() +
                      " in file " + file_path);
    }

}