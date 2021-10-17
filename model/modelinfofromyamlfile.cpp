//
// modelinfofromyamlfile.cpp
// 
/***************************************************************************
 *   Created by James Barbetti on 30-Apr-2021                              *
 *   james_barbetti@yahoo.com                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "modelinfofromyamlfile.h"
#include "modelexpression.h"
#include <utils/stringfunctions.h> //for string_to_lower, startsWith, endsWith
#include <tree/phylotree.h>        //for TREE_LOG_LINE macro
#include "yamlmodelwrapper.h"      //for YAMLRateFree and friends.

//YAML Logging Levels
VerboseMode YAMLModelVerbosity     = VerboseMode::VB_MIN;
VerboseMode YAMLFrequencyVerbosity = VerboseMode::VB_MAX;
VerboseMode YAMLMatrixVerbosity    = VerboseMode::VB_MAX;
VerboseMode YAMLParsingVerbosity   = VerboseMode::VB_MAX;
VerboseMode YAMLRateVerbosity      = VerboseMode::VB_MAX;
VerboseMode YAMLVariableVerbosity  = VerboseMode::VB_MAX;

YAMLFileParameter::YAMLFileParameter()
    : is_subscripted(false), minimum_subscript(0), maximum_subscript(0)
    , type(ModelParameterType::OTHER), tolerance(0.001), value(0.0)  {
}

std::string YAMLFileParameter::getSubscriptedVariableName(int subscript) const {
    if (!is_subscripted) {
        return name;
    }
    std::stringstream subscripted_name;
    subscripted_name << name << "(" << subscript << ")";
    return subscripted_name.str();
}

bool YAMLFileParameter::isMatchFor(const std::string& match_name) const {
    return string_to_lower(name) == match_name;
}

bool YAMLFileParameter::isMatchFor(const std::string& match_name,
    ModelParameterType match_type) const {
    return (type == match_type &&
        string_to_lower(name) == match_name);
}

void YAMLFileParameter::logParameterState(const char* verb, 
                                          LoggingTarget* logging_target) const {
    std::stringstream msg;
    msg << verb << " ";
    if (is_subscripted) {
        msg  << "subscripted parameter " << name
             << "(" << minimum_subscript << " .. " << maximum_subscript << ")";
    } else {
        msg  << "parameter " << name;
    }

    msg  << " of type " << type_name
         << ", with range " << range.first
         << " to " << range.second
         << ", and initial value " << value;

    if (0<tolerance) {
        msg << ", and tolerance " << tolerance;
    }
    TREE_LOG_LINE(*logging_target, YAMLVariableVerbosity, msg.str());
}

ModelVariable::ModelVariable() : type(ModelParameterType::OTHER)
, value(0), is_fixed(false) {
}

ModelVariable::ModelVariable(ModelParameterType t,
    const ModelParameterRange& r,
    double v)
    : range(r), type(t), value(v), is_fixed(false) {
}

void ModelVariable::setTypeName(const std::string& type_name) {
    type = modelParameterTypeFromString(type_name);
}

void ModelVariable::setMinimum(double min_value) {
    range.first = min_value;
}

void ModelVariable::setMaximum(double max_value) {
    range.second = max_value;
}

bool ModelVariable::constrainValueToRange() {
    if (value < range.first ) {
        value = range.first;
        return true;
    } else if (range.second < value) {
        value = range.second;
        return true;
    }
    return false;
}

void ModelVariable::setValue(double v) {
    value = v;
}

void ModelVariable::markAsFixed() {
    is_fixed = true;
}

ModelParameterType ModelVariable::getType() const {
    return type;
}

ModelParameterRange ModelVariable::getRange() const {
    return range;
}

std::string ModelVariable::getTypeName() const {
    return modelParameterTypeToString(type);
}

double ModelVariable::getValue() const {
    return value;
}

bool ModelVariable::isFixed() const {
    return is_fixed;
}

void StringMatrix::makeRectangular(size_t column_count) {
    size_t row_count = size();
    for (int row_num = 0; row_num < row_count; ++row_num) {
        StrVector& row = at(row_num);
        if (row.size() != column_count) {
            row.resize(row_count, "");
        }
    }
}

void StringMatrix::makeSquare(bool reflect) {
    size_t row_count = size();
    size_t col_count = row_count;
    for (StrVector& row : *this) {
        if (col_count < row.size()) {
            col_count = row.size();
        }
    }
    for (size_t row_num = 0;
        row_num < size(); ++row_num) {
        StrVector& row = at(row_num);
        int old_col_count = static_cast<int>(row.size());
        if (old_col_count < col_count) {
            row.resize(col_count, "");
            if (reflect) {
                for (int col_num = old_col_count;
                    col_num < col_count; ++col_num) {
                    if (col_num == row_num) {
                        continue;
                    }
                    if (row_num < at(col_num).size()) {
                        row[col_num] = at(col_num)[row_num];
                    }
                }
            }
        }
    }
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile()
    : is_modifier_model(false), is_rate_model(false)
    , sequence_type(SeqType::SEQ_UNKNOWN), num_states(0)
    , reversible(false), rate_matrix_rank(0)
    , tip_likelihood_rank(0)
    , frequency_type(StateFreqType::FREQ_UNKNOWN)
    , parent_model(nullptr), linked_models(nullptr)
    , mixed_models(nullptr), weight_formula("1"), model_weight(1.0)
    , subtree_models(nullptr)
    , specified_rate_model_info(nullptr) {
}

void ModelInfoFromYAMLFile::copyMixedAndLinkedModels
        (const ModelInfoFromYAMLFile& rhs) {
    delete mixed_models;
    mixed_models = nullptr;
    if (rhs.mixed_models != nullptr) {
        mixed_models = new MapOfModels(*rhs.mixed_models);
    }

    delete linked_models;
    linked_models = nullptr;
    if (rhs.linked_models != nullptr) {
        linked_models = new MapOfModels(*rhs.linked_models);
    }

    delete subtree_models;
    subtree_models = nullptr;
    if (rhs.subtree_models != nullptr) {
        subtree_models = new MapOfModels(*rhs.subtree_models);
    }
    clade_names          = rhs.clade_names;
    distinct_clade_names = rhs.distinct_clade_names;

    delete specified_rate_model_info;
    specified_rate_model_info = nullptr;
    if (rhs.specified_rate_model_info != nullptr) {
        specified_rate_model_info = new ModelInfoFromYAMLFile
                                        (*rhs.specified_rate_model_info);
    }

}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const ModelInfoFromYAMLFile& rhs)
    : model_name(rhs.model_name), model_file_path(rhs.model_file_path)
    , superclass_model_name(rhs.superclass_model_name)
    , is_modifier_model(rhs.is_modifier_model), is_rate_model(rhs.is_rate_model)
    , rate_distribution(rhs.rate_distribution)
    , citation(rhs.citation), DOI(rhs.DOI), url(rhs.url)
    , data_type_name(rhs.data_type_name), sequence_type(rhs.sequence_type)
    , num_states(rhs.num_states), reversible(rhs.reversible)
    , rate_matrix_rank(rhs.rate_matrix_rank)
    , rate_matrix_expressions(rhs.rate_matrix_expressions)
    , rate_matrix_formula(rhs.rate_matrix_formula)
    , tip_likelihood_rank(rhs.tip_likelihood_rank)
    , tip_likelihood_expressions(rhs.tip_likelihood_expressions)
    , tip_likelihood_formula(rhs.tip_likelihood_formula)
    , parameters(rhs.parameters), frequency_type(rhs.frequency_type)
    , variables(rhs.variables), parent_model(rhs.parent_model)
    , linked_models(nullptr), mixed_models(nullptr)
    , weight_formula(rhs.weight_formula), model_weight(rhs.model_weight)
    , subtree_models(nullptr), clade_names(rhs.clade_names)
    , distinct_clade_names(rhs.distinct_clade_names)
    , specified_rate_model_info(nullptr) {
    copyMixedAndLinkedModels(rhs);
}

ModelInfoFromYAMLFile& ModelInfoFromYAMLFile::operator=(const ModelInfoFromYAMLFile& rhs) {
    if (&rhs != this) {
        model_name                 = rhs.model_name;
        model_file_path            = rhs.model_file_path;
        superclass_model_name      = rhs.superclass_model_name;
        is_modifier_model          = rhs.is_modifier_model;
        is_rate_model              = rhs.is_rate_model;
        rate_distribution          = rhs.rate_distribution;
        citation                   = rhs.citation;
        DOI                        = rhs.DOI; 
        url                        = rhs.url;
        data_type_name             = rhs.data_type_name;
        sequence_type              = rhs.sequence_type;
        num_states                 = rhs.num_states;
        reversible                 = rhs.reversible;
        rate_matrix_rank           = rhs.rate_matrix_rank;
        rate_matrix_expressions    = rhs.rate_matrix_expressions;
        rate_matrix_formula        = rhs.rate_matrix_formula;
        tip_likelihood_rank        = rhs.tip_likelihood_rank;
        tip_likelihood_expressions = rhs.tip_likelihood_expressions;
        tip_likelihood_formula     = rhs.tip_likelihood_formula;
        parameters                 = rhs.parameters; 
        frequency_type             = rhs.frequency_type;
        variables                  = rhs.variables;
        weight_formula             = rhs.weight_formula;
        model_weight               = rhs.model_weight;
        //Note: parent_model is NOT set.
        copyMixedAndLinkedModels(rhs);
    }
    return *this;
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const std::string& path)
    : model_file_path(path), is_modifier_model(false), is_rate_model(false)
    , sequence_type(SeqType::SEQ_UNKNOWN)
    , num_states(0), reversible(true)
    , rate_matrix_rank(0), tip_likelihood_rank(0)
    , frequency_type(StateFreqType::FREQ_UNKNOWN)
    , parent_model(nullptr), linked_models(nullptr)
    , mixed_models(nullptr), model_weight(1.0)
    , subtree_models(nullptr), specified_rate_model_info(nullptr) {
}

ModelInfoFromYAMLFile::~ModelInfoFromYAMLFile() {
    delete mixed_models;
    mixed_models = nullptr;
    delete linked_models;
    linked_models = nullptr;
    delete subtree_models;
    subtree_models = nullptr;
    delete specified_rate_model_info;
    specified_rate_model_info = nullptr;
}

void ModelInfoFromYAMLFile::specifyRateModel(const ModelInfoFromYAMLFile& ancestor,
                                             LoggingTarget* logging_target) {
    if (specified_rate_model_info==nullptr) {
        specified_rate_model_info = new ModelInfoFromYAMLFile(ancestor);
    } else {
        specified_rate_model_info->inheritModel(ancestor, logging_target);
    }
}

bool ModelInfoFromYAMLFile::isGammaModel() const {
    for (auto p : parameters) {
        if (p.type == ModelParameterType::SHAPE) {
            //Todo: decide. Should it have to have name 
            //that contains the substring, "gamma"?
            return true;
        }
    }
    return false;
}

bool ModelInfoFromYAMLFile::isInvariantModel() const {
    for (auto p : parameters) {
        if (p.type == ModelParameterType::INVARIANT_PROPORTION) {
            return true;
        }
    }
    return false;
}

bool ModelInfoFromYAMLFile::isKategoryModel()  const {
    return contains(rate_distribution, "kategory");
}
    
bool ModelInfoFromYAMLFile::isMixtureModel() const {
    return mixed_models != nullptr;
}

bool ModelInfoFromYAMLFile::isDivergentModel() const {
    return subtree_models != nullptr;
}

bool ModelInfoFromYAMLFile::isModelFinder() const {
    return false;
}

bool ModelInfoFromYAMLFile::isModelFinderOnly() const {
    return false;
}

bool ModelInfoFromYAMLFile::isReversible() const {
    return reversible;
}

int  ModelInfoFromYAMLFile::getKategoryRateCount(int rate_count, int min_count) const {
    rate_count = getNumberOfRateCategories();
    if (rate_count<min_count) {
         outError("Wrong number of rate categories");
    }
    return rate_count;
}

void ModelInfoFromYAMLFile::updateName(const std::string& name) {
    model_name = name;
}

std::string ModelInfoFromYAMLFile::getQualifiedName() const {
    auto qualified_name = model_name;
    for (auto ancestor=parent_model; ancestor!=nullptr;
         ancestor=ancestor->parent_model) {
        qualified_name = ancestor->model_name + "." + qualified_name;
    }
    return qualified_name;
}

std::string ModelInfoFromYAMLFile::getLongName() const {
    return model_name + " from YAML model file " +
        model_file_path;
}

bool ModelInfoFromYAMLFile::hasDot(const char* name) const {
    return strstr(name, ".") != nullptr;
}

void ModelInfoFromYAMLFile::breakAtDot(const char* name,
    std::string& sub_model_name,
    const char*& remainder) const {
    auto dotpos = strstr(name, ".");
    if (dotpos == nullptr) {
        sub_model_name.clear();
        remainder = name;
    }
    else {
        sub_model_name = std::string(name, dotpos - name);
        remainder = name + 1;
    }
}

MapOfModels::const_iterator
ModelInfoFromYAMLFile::findMixedModel(const std::string& name) const {
    auto it = mixed_models->find(name);
    if (it == mixed_models->end()) {
        std::stringstream complaint;
        complaint << "Could not evaluate mixed model name " << name
                  << " for model " << getLongName();
        throw ModelExpression::ModelException(complaint.str());
    }
    return it;
}

MapOfModels::const_iterator
    ModelInfoFromYAMLFile::findSubtreeModel
        (const std::string& name) const {
    auto it = subtree_models->find(name);
    if (it == subtree_models->end()) {
        std::stringstream complaint;
        complaint << "Could not evaluate subtree model name " << name
                  << " for model " << getLongName();
        throw ModelExpression::ModelException(complaint.str());
    }
    return it;
}

void ModelInfoFromYAMLFile::setNumberOfStatesAndSequenceType
        (int requested_num_states,
         LoggingTarget* logging_target) {
    if (requested_num_states != 0) {
        num_states = requested_num_states;
    }
    if (num_states == 0) {
        num_states = 4;
    }
    if (!data_type_name.empty()) {
        TREE_LOG_LINE(*logging_target, YAMLParsingVerbosity,
                      "Data Type Name is " << data_type_name);
        auto seq_type_requested = getSeqType(data_type_name.c_str());
        if (seq_type_requested != SeqType::SEQ_UNKNOWN) {
            sequence_type = seq_type_requested;
            if (sequence_type == SeqType::SEQ_CODON) {
                num_states = 61;
            } else {
                num_states = getNumStatesForSeqType(sequence_type, num_states);
            }
        }
    }
    if (sequence_type == SeqType::SEQ_UNKNOWN) {
        switch (num_states) {
        case 2:   sequence_type = SeqType::SEQ_BINARY;  break;
        case 4:   sequence_type = SeqType::SEQ_DNA;     break;
        case 20:  sequence_type = SeqType::SEQ_PROTEIN; break;
        case 61:  sequence_type = SeqType::SEQ_CODON;   break;
        default:  /* still, no idea */
            outWarning("Could not determine sequence type"
                " for model " + model_name);
            break;
        }
    }
    TREE_LOG_LINE(*logging_target, YAMLParsingVerbosity,
                  "Number of states is " << num_states);

    forceAssign("num_states", num_states);
    forceAssign("numStates", num_states);
}

double ModelInfoFromYAMLFile::evaluateExpression(const std::string& expr,
                                                 const std::string& context) {
    const char* verb = "parsing";
    try {
        ModelExpression::InterpretedExpression interpreter(*this, expr);
        verb = "evaluating";
        return interpreter.evaluate();
    }
    catch (ModelExpression::ModelException& x) {
        std::stringstream msg;
        msg << "Error " << verb << " " << context
            << " for " << model_name << ":\n"
            << x.getMessage();
        outError(msg.str());
    }
    return 0;
}

const YAMLFileParameter* ModelInfoFromYAMLFile::findParameter
(const char* name) const {
    size_t      param_count = parameters.size();
    std::string lower_name = string_to_lower(name);
    for (size_t i = 0; i < param_count; ++i) {
        if (parameters[i].isMatchFor(lower_name)) {
            return &parameters[i];
        }
    }
    if (parent_model!=nullptr) {
        return parent_model->findParameter(name);    
    }
    return nullptr;
}

const YAMLFileParameter* ModelInfoFromYAMLFile::findParameter
                         (const std::string& name) const {
    return findParameter(name.c_str());
}

const YAMLFileParameter* ModelInfoFromYAMLFile::findParameter
                         (const char* name, ModelParameterType type) const {
    size_t      param_count = parameters.size();
    std::string lower_name = string_to_lower(name);
    for (size_t i = 0; i < param_count; ++i) {
        if (parameters[i].isMatchFor(lower_name, type)) {
            return &parameters[i];
        }
    }
    if (parent_model!=nullptr) {
        return parent_model->findParameter(name, type);    
    }
    return nullptr;
}

const YAMLFileParameter* ModelInfoFromYAMLFile::findParameter
                         (const std::string& name, ModelParameterType type) const {
    return findParameter(name.c_str(), type);
}

void  ModelInfoFromYAMLFile::moveParameterToBack
      (const char* name, ModelParameterType type) {
    size_t      param_count = parameters.size();
    std::string lower_name = string_to_lower(name);
    size_t      i;
    for (i = 0; i < param_count; ++i) {
        if (parameters[i].isMatchFor(lower_name, type)) {
            break;
        }
    }
    if (i == param_count) {
        return; //Didn't find the parameter to be moved
    }
    YAMLFileParameter lift = parameters[i];
    for (size_t j = i + 1; j < param_count; ++j) {
        if (parameters[j].type == type) {
            parameters[i] = parameters[j];
            i = j;
        }
    }
    parameters[i] = lift;
}

const ModelInfoFromYAMLFile::Variables& 
    ModelInfoFromYAMLFile::getVariables() const {
    return variables;
}

const ModelVariable* 
    ModelInfoFromYAMLFile::getVariableByName
        (const char* name) const {
    //Todo: look up linked_models too
    if (hasDot(name) && mixed_models != nullptr ) {
        std::string sub_model_name;
        const char* var_name = nullptr;
        breakAtDot(name, sub_model_name, var_name);
        auto it = findMixedModel(sub_model_name);
        return (*it)->getVariableByName(var_name);
    }
    auto it = variables.find(name);
    if (variables.find(name) != variables.end()) {
        return &(it->second);
    }
    if (parent_model==nullptr) {
        return nullptr;
    }
    return parent_model->getVariableByName(name);
}

const ModelVariable* ModelInfoFromYAMLFile::getVariableByName(const std::string& name) const {
    return getVariableByName(name.c_str());
}

bool ModelInfoFromYAMLFile::hasVariable(const char* name) const {
    return getVariableByName(name) != nullptr;
}

bool ModelInfoFromYAMLFile::hasVariable(const std::string& name) const {
    return hasVariable(name.c_str());
}

double ModelInfoFromYAMLFile::getVariableValue(const char* name) const {
    const ModelVariable* pvar = getVariableByName(name);
    if (pvar==nullptr) {
        return 0.0;
    }
    return pvar->getValue();
}

double ModelInfoFromYAMLFile::getVariableValue(const std::string& name) const {
    return getVariableValue(name.c_str());
}

bool ModelInfoFromYAMLFile::isVariableFixed(const std::string& name) const {
    const ModelVariable* pvar = getVariableByName(name);
    if (pvar==nullptr) {
        return false;
    }
    return pvar->isFixed();
}

void ModelInfoFromYAMLFile::addParameter(const YAMLFileParameter& p,
                                         LoggingTarget* logging_target) {
    bool replaced = false;
    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        if (it->name == p.name) {
            *it = p;
            replaced = true;
            break;
        }
    }
    //Problem here: What if this is a child model in a mixture model
    //and there is a parameter with this name, in the parent model?
    if (!replaced) {
        parameters.emplace_back(p);
    }
    if (p.is_subscripted) {
        for (int i = p.minimum_subscript; i <= p.maximum_subscript; ++i) {
            forceAssign("subscript", i);
            //Doesn't do forceAssign outside the loop and assign the 
            //ModelVariable& that comes back from it, because the code
            //here might *add* a variable, and depending on the container
            //type of Variables... the ModelVariable& could go wrong
            //(and that'd be *nasty* insidious coupling).

            setSubscriptedVariable(p, i, logging_target);
        }
    }
    else {
        variables[p.name] = ModelVariable(p.type, p.range, p.value);
    }
}

void ModelInfoFromYAMLFile::setSubscriptedVariable
        (const YAMLFileParameter& p, int i, 
         LoggingTarget* logging_target) {
    typedef ModelExpression::InterpretedExpression Interpreter;
    typedef ModelExpression::ModelException        Exception;
    std::string var_name = p.getSubscriptedVariableName(i);
    double      v        = p.value;
    if (!p.init_expression.empty()) {
        try {
            Interpreter ix(*this, p.init_expression);
            v = ix.evaluate();
        }
        catch (Exception& x) {
            std::stringstream complaint;
            complaint << "Error initializing " << p.name 
                        << "(" << i << ")";
            throw Exception(complaint.str());
        }
    }
    if (i<11) {
        TREE_LOG_LINE(*logging_target, YAMLVariableVerbosity, 
                    "Setting " << getName() << "." 
                    << var_name << "=" << p.init_expression 
                    << (p.init_expression.empty() ? "" : "=")
                    << v );
    }
    variables[var_name] = ModelVariable(p.type, p.range, v);
}

bool ModelInfoFromYAMLFile::removeSubscriptedVariable
        (const YAMLFileParameter& param, int i,
         LoggingTarget* logging_target) {
    std::string dead_var = param.getSubscriptedVariableName(i);
    auto it = variables.find(dead_var);
    if (it != variables.end()) {
        variables.erase(it);
        return true;
    }
    return false;
}

void ModelInfoFromYAMLFile::updateParameterSubscriptRange(YAMLFileParameter& p, 
                                                          int min_subscript,
                                                          int max_subscript,
                                                          LoggingTarget* logging_target) {
    typedef ModelExpression::InterpretedExpression Interpreter;
    typedef ModelExpression::ModelException        Exception;
    if (max_subscript<min_subscript) {
        std::stringstream complaint;
        complaint << "New subscript range for " << p.name
                  << " (" << min_subscript << ".." << max_subscript <<")"
                  << " is invalid.";
        throw Exception(complaint.str());
    }
    for (int i=p.minimum_subscript; i<=p.maximum_subscript; ++i) {
        if (i<min_subscript || max_subscript<i) {
            std::string dead_var_name = p.getSubscriptedVariableName(i);
            auto it = variables.find(dead_var_name);
            if (it!=variables.end()) {
                variables.erase(it);
            }
        }
    }
    for (int i=min_subscript; i<=max_subscript; ++i) {
        forceAssign("subscript", i);
        //Doesn't do forceAssign outside the loop and assign the 
        //ModelVariable& that comes back from it, because the code
        //here might *add* a variable, and depending on the container
        //type of Variables... the ModelVariable& could go wrong
        //(and that'd be *nasty* insidious coupling).

        std::string new_var_name = p.getSubscriptedVariableName(i);
        const char* verb = "Added";
        auto it = variables.find(new_var_name);
        if (it!=variables.end()) {
            if (it->second.isFixed()) {
                continue;
            }
            verb = "Updated";
        }
        double v = 0;
        if (!p.init_expression.empty()) {
            try {
                Interpreter ix(*this, p.init_expression);
                v = ix.evaluate();
            }
            catch (Exception& x) {
                std::stringstream complaint;
                complaint << "Error initializing " << p.name
                            << "(" << i <<")=(" << p.init_expression << ")"
                            << ": " << x.getMessage();
                throw Exception(complaint.str());
            }
        }
        variables[new_var_name] = ModelVariable(p.type, p.range, v);
        TREE_LOG_LINE(*logging_target, YAMLVariableVerbosity,
                      verb << " subscripted variable " << new_var_name << "=" << v);
    }
}

bool  ModelInfoFromYAMLFile::hasFrequencyParameters(int min_variable_count) const {
    for (auto p : parameters) {
        if (p.type == ModelParameterType::FREQUENCY) {
            if (p.is_subscripted) {
                min_variable_count -= (p.maximum_subscript - p.minimum_subscript + 1);
            } else {
                --min_variable_count;
            }
            if (min_variable_count<=0) {
                return true;
            }
        }
    }
    return false;
}

bool ModelInfoFromYAMLFile::isFrequencyParameter(const std::string& param_name) const {
    for (auto p : parameters) {
        if (string_to_lower(p.name) == string_to_lower(param_name)) {
            return p.type == ModelParameterType::FREQUENCY;
        }
    }
    return false;
}

StateFreqType ModelInfoFromYAMLFile::getFrequencyType() const {
    return frequency_type;
}

void ModelInfoFromYAMLFile::setFrequencyType(StateFreqType new_type) {
    frequency_type = new_type;
} 

double ModelInfoFromYAMLFile::getModelWeight() const {
    return model_weight;
}

double ModelInfoFromYAMLFile::getModelWeight() {
    ASSERT(!weight_formula.empty());
    try {
        ModelExpression::InterpretedExpression expr(*this, weight_formula);
        model_weight = expr.evaluate();
    }
    catch (ModelExpression::ModelException& x) {
        std::stringstream complaint;
        complaint << "Error determining weight of model " << getName()
                  << ": " << x.getMessage();
        throw ModelExpression::ModelException(complaint.str());
    }
    return model_weight;
}

bool ModelInfoFromYAMLFile::isModelWeightFixed() {
    ModelExpression::InterpretedExpression expr(*this, weight_formula);
    return expr.expression()->isFixed();
}

void ModelInfoFromYAMLFile::setBounds(int param_count, 
                                      const std::vector<ModelParameterType>& types,
                                      double* lower_bound, double* upper_bound, 
                                      bool* bound_check,
                                      PhyloTree* report_to_tree) const {
    int i = 1; //Parameter numbering starts at 1, see ModelMarkov
    for (auto param_type : types) {
        bool is_freq = (param_type == ModelParameterType::FREQUENCY);
        for (auto p : parameters) {
            if (p.type == param_type) {
                for (int sub = p.minimum_subscript;
                    sub <= p.maximum_subscript; ++sub) {
                    std::string var_name = p.getSubscriptedVariableName(sub);
                    if (i <= param_count) { 
                        lower_bound[i] = p.range.first;
                        upper_bound[i] = p.range.second;
                        bound_check[i] = false;
                        TREE_LOG_LINE(*report_to_tree, YAMLVariableVerbosity,
                                    "Set bound " << i 
                                    << " : " << var_name 
                                    << " [" << lower_bound[i]
                                    << ".." << upper_bound[i] << "]");
                    } else if (is_freq) {
                        ASSERT(i <= param_count + 1);
                        TREE_LOG_LINE(*report_to_tree, YAMLVariableVerbosity,
                                    "Did *not* set (last freq?) bound " << i 
                                    << " : " << var_name );
                    } else {
                        ASSERT(i <= param_count);
                    }
                    ++i;
                }
            }
        }
    }
}

void ModelInfoFromYAMLFile::updateModelVariables(const double* updated_values,
                                                 int first_freq_index,
                                                 int last_param_index,
                                                 PhyloTree* report_to_tree) {
    int i = 1; //Rate parameter numbering starts at 1, see ModelMarkov
    ModelParameterType supported_types[] = {
        ModelParameterType::PROPORTION,
        ModelParameterType::INVARIANT_PROPORTION, 
        ModelParameterType::RATE, 
        ModelParameterType::FREQUENCY 
    };
    //PROPORTION must be before RATE (e.g. RateFree)
    //INVARIANT_PROPORTION must be be after (because it is, e.g. in RateFreeInvar)
    //FREQUENCY  must be after RATE  (e.g. ModelMarkov)
    //Todo: Where do weight parameters go?    
    //      Especially in mixture models
    for (auto param_type : supported_types) {
        if (param_type == ModelParameterType::FREQUENCY) {
            i = first_freq_index;
        }
        updateModelVariablesByType(updated_values, last_param_index, 
                                   false, param_type, i, report_to_tree);
    }
}

bool ModelInfoFromYAMLFile::updateModelVariablesByType(const double* updated_values,
                                                       int param_count,
                                                       bool even_fixed_ones,                                                
                                                       ModelParameterType param_type,
                                                       int &i,
                                                       PhyloTree* report_to_tree) {
    bool anyChanges                  = false;
    int  skipped_frequency_variables = 0;
    for (auto p : parameters) {
        if (p.type != param_type) {
            continue;
        }
        for (int sub = p.minimum_subscript;
            sub <= p.maximum_subscript; ++sub) {
            std::string var_name = p.getSubscriptedVariableName(sub);
            ModelVariable& var   = this->variables[var_name];
            if (var.isFixed() && !even_fixed_ones) {
                continue;
            }
            if (param_count<i) {
                if (param_type == ModelParameterType::FREQUENCY) {
                    //The last frequency variable won't be set this way
                    ++skipped_frequency_variables;
                    if (skipped_frequency_variables==1) {
                        continue;
                    }
                }
                std::stringstream complaint;
                complaint << "Internal logic error: Cannot use variable "
                            << " with 1-based-index " << i << " to assign "
                            << " model variable " << var_name
                            << " of " << getLongName()
                            << " (too many model variables found).";
                outError(complaint.str());
            }
            double value     = updated_values[i];
            double old_value = var.getValue();
            double delta     = value - old_value;
            if (p.tolerance<=fabs(delta)) {
                var.setValue(value);
                anyChanges = true;
                TREE_LOG_LINE(*report_to_tree, YAMLVariableVerbosity,
                            "Updated " << var_name << "=" << value
                            << " (delta " << delta << ") "
                            << " from parameter " << i 
                            << " of " << param_count);
            } else {
                TREE_LOG_LINE(*report_to_tree, YAMLVariableVerbosity,
                            "Did not update " << var_name << "=" << old_value
                            << " (to " << value << ") "
                            << " from parameter " << i 
                            << " of " << param_count 
                            << " (delta was " << delta << ")");
            }
            ++i;
        }
    }
    return anyChanges;
}

void ModelInfoFromYAMLFile::readModelVariablesByType
        (double* write_them_here, int param_count,
         bool even_fixed_ones, ModelParameterType param_type,
         int &i, PhyloTree* report_to_tree) const {                                                
    for (auto p : parameters) {
        if (p.type != param_type) {
            continue;
        }
        for (int sub = p.minimum_subscript;
            sub <= p.maximum_subscript; ++sub) {
            std::string var_name = p.getSubscriptedVariableName(sub);
            auto it = variables.find(var_name);
            if (it==variables.end()) {
                std::stringstream complaint;
                complaint << "Internal logic error: Cannot set variable "
                            << " with 1-based-index " << i << " from "
                            << " model variable " << var_name
                            << " of " << getLongName()
                            << " (model does not define variable).";
                outError(complaint.str());
                return;
            }
            const ModelVariable& var = it->second;
            if (var.isFixed() && !even_fixed_ones) {
                continue;
            }
            if (param_count<i) {
                std::stringstream complaint;
                complaint << "Internal logic error: Cannot set variable "
                            << " with 1-based-index " << i << " from "
                            << " model variable " << var_name
                            << " of " << getLongName()
                            << " (too many model variables found).";
                outError(complaint.str());
            }
            double v = var.getValue();
            write_them_here[i] = v;
            TREE_LOG_LINE(*report_to_tree, YAMLVariableVerbosity,
                          "Wrote from " << var_name << "=" << v
                          << " (subscript " << sub << ")"
                          << " to parameter " << i 
                          << " of " << param_count);
            ++i;
        }
    }
}

void ModelInfoFromYAMLFile::logVariablesTo
        (LoggingTarget* logging_target) const {
    if (verbose_mode < YAMLVariableVerbosity) {
        return;
    }
    std::stringstream var_list;
    const char* sep = "Variables: ";
    for (auto itv : variables) {
        var_list << sep << itv.first << "=" << itv.second.getValue();
        sep = ", ";
    }
    std::string list = var_list.str();
    if (list.find("nan") != std::string::npos) {
        list += " ...?";
    }
    TREE_LOG_LINE(*logging_target, YAMLModelVerbosity, list);
}

ModelVariable& ModelInfoFromYAMLFile::assign(const std::string& var_name,
    double value_to_set) {
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        //Todo: look up linked_models too
        bool dotted = hasDot(var_name.c_str());
        if (dotted) {
            std::string sub_model_name;
            const char* sub_model_var_name = nullptr;
            breakAtDot(var_name.c_str(), sub_model_name,
                       sub_model_var_name);
        
            if (subtree_models != nullptr) {
                auto it_model = findSubtreeModel(sub_model_name);
                return (*it_model)->assign(std::string(sub_model_var_name),
                                    value_to_set);
            }
            if (mixed_models != nullptr) {
                auto it_model = findMixedModel(sub_model_name);
                return (*it_model)->assign(std::string(sub_model_var_name),
                                    value_to_set);
            }
        }
        std::stringstream complaint;
        complaint << "Could not assign"
            << " to unrecognized variable " << var_name
            << " of model " << model_name << ".";
        outError(complaint.str());
    }
    it->second.setValue(value_to_set);
    return it->second;
}

ModelVariable& ModelInfoFromYAMLFile::forceAssign(const std::string& var_name,
                                                  double value_to_set) {
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        //Todo: look up linked_models too
        bool dotted = hasDot(var_name.c_str());
        if (dotted) {
            std::string sub_model_name;
            const char* sub_model_var_name = nullptr;
            breakAtDot(var_name.c_str(), sub_model_name,
                sub_model_var_name);
            if (subtree_models != nullptr) {
                auto it_model = findSubtreeModel(sub_model_name);
                return (*it_model)->forceAssign(std::string(sub_model_var_name),
                                                value_to_set);
            }
            else if (mixed_models != nullptr) {
                auto it_model = findMixedModel(sub_model_name);
                return (*it_model)->forceAssign(std::string(sub_model_var_name),
                                                value_to_set);
            }
        }
        ModelVariable var(ModelParameterType::OTHER,
                          ModelParameterRange(), value_to_set);
        variables[var_name] = var;
        it = variables.find(var_name);
    }
    it->second.setValue(value_to_set);
    return it->second;
}

ModelVariable& ModelInfoFromYAMLFile::forceAssign(const std::string& var_name,
                                                  int value_to_set) {
    return forceAssign(var_name, static_cast<double>(value_to_set));
}

ModelVariable& ModelInfoFromYAMLFile::forceAssign(const char* var_name,
                                                  int value_to_set) {
    return forceAssign(std::string(var_name), static_cast<double>(value_to_set));
}

const StrVector& ModelInfoFromYAMLFile::getVariableNamesByPosition() const {
    variable_names.clear();
    ModelParameterType supported_types[] = {
        ModelParameterType::WEIGHT, 
        ModelParameterType::SHAPE, 
        ModelParameterType::PROPORTION, 
        ModelParameterType::RATE, 
        ModelParameterType::FREQUENCY,
        ModelParameterType::OTHER,
    };
    //FREQUENCY must be after RATE.
    //Todo: Where do weight parameters go?    Esepecially in mixture
    //      models
    for (auto param_type : supported_types) {
        for (auto p : parameters) {
            if (p.type == param_type) {
                for (int sub = p.minimum_subscript;
                    sub <= p.maximum_subscript; ++sub) {
                    std::string var_name = p.getSubscriptedVariableName(sub);
                    variable_names.emplace_back(var_name);
                }
            }
        }
    }
    return variable_names;
}

std::string ModelInfoFromYAMLFile::getVariableNameByPosition(int position) const {
    if (position<0 || variable_names.size()<=position) {
        std::stringstream complaint;
        complaint << "Cannot set variable with (1-based) position"
                  << " " << (position+1) << ". Only variables with positions"
                  << " 1 through " << (variable_names.size()+1) 
                  << " can be supplied.";
        const char* prefix = "\nVariable names are: ";
        for (size_t pos=0; pos<variable_names.size(); ++pos) {
            complaint << prefix << variable_names[pos];
            prefix = ", ";
        }
        complaint << ".";
        throw ModelExpression::ModelException(complaint.str());
    } else {
        return variable_names[position];
    }
}

ModelVariable& ModelInfoFromYAMLFile::assignByPosition(size_t position,
    double value_to_set) {
    if (variable_names.size() < variables.size()) {
        getVariableNamesByPosition();
    }
    if (variable_names.size() <= position) {
        std::stringstream complaint;
        complaint << "Could not assign parameter " << (position + 1)
            << " as there are only " << variable_names.size()
            << " parameters in model " << model_name << ".";
        outError(complaint.str());
    }
    return assign(variable_names[position], value_to_set);
}

bool ModelInfoFromYAMLFile::assignLastFrequency(double value) {
    //Ee-uw.  It's sad that this is needed!
    for (auto pit = parameters.rbegin();
        pit != parameters.rend(); ++pit) {
        if (pit->type == ModelParameterType::FREQUENCY) {
            for (int sub = pit->maximum_subscript;
                pit->minimum_subscript <= sub; --sub) {
                std::string var_name = pit->getSubscriptedVariableName(sub);
                if (!this->variables[var_name].isFixed()) {
                    this->variables[var_name].setValue(value);
                    return true;
                }
                //For the last frequency parameter, if there is one,
                //i will be equal to ( param_count + 1 ).
            }
        }
    }
    return false;
}

const YAMLFileParameter& ModelInfoFromYAMLFile::getInvariantProportionParameter() const {
    for (const YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::INVARIANT_PROPORTION) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has invariant proportion parameter");
    return parameters.front();
} 

YAMLFileParameter& ModelInfoFromYAMLFile::getInvariantProportionParameter()  {
    for ( YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::INVARIANT_PROPORTION) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has invariant proportion parameter");
    return parameters.front();
} 

const YAMLFileParameter& ModelInfoFromYAMLFile::getProportionParameter() const {
    for (const YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::PROPORTION) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has no proportion parameter");
    return parameters.front();
} 

YAMLFileParameter& ModelInfoFromYAMLFile::getProportionParameter()  {
    for ( YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::PROPORTION) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has no proportion parameter");
    return parameters.front();
} 

const ModelVariable* ModelInfoFromYAMLFile::getInvariantProportionVariable() const {
    for (const YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::INVARIANT_PROPORTION) {
            auto i = parameter.minimum_subscript;
            std::string var_name = parameter.getSubscriptedVariableName(i);
            return getVariableByName(var_name);
        }
    }
    return nullptr;
}

const YAMLFileParameter& ModelInfoFromYAMLFile::getRateParameter() const {
    for (const YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::RATE) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has no rate parameter");
    return parameters.front();
} 

YAMLFileParameter& ModelInfoFromYAMLFile::getRateParameter()  {
    for ( YAMLFileParameter& parameter: parameters) {
        if (parameter.type ==  ModelParameterType::RATE) {
            return parameter;
        }
    }
    ASSERT(0 && "model info has no rate parameter");
    return parameters.front();
} 

std::string ModelInfoFromYAMLFile::getStringProperty(const char* name,
    const char* default_value) const {
    std::string low_name = string_to_lower(name);
    auto        it       = string_properties.find(low_name);
    return it == string_properties.end() ? default_value : it->second;
}

bool ModelInfoFromYAMLFile::hasStringProperty
        (const char* name, std::string& value) const {
    std::string low_name = string_to_lower(name);
    auto it = string_properties.find(low_name);
    if (it==string_properties.end()) {
        return false;
    }
    value = it->second;
    return true;
}

int ModelInfoFromYAMLFile::getNumStates() const {
    //Todo: decide from the data type (or at least from data_type_name!)
    return 4; //but, for now, hardcoded!
}

int  ModelInfoFromYAMLFile::getTipLikelihoodMatrixRank() const {
    return tip_likelihood_rank;
}

void ModelInfoFromYAMLFile::computeTipLikelihoodsForState
        (int state, int num_states, double* likelihoods) {
    std::stringstream complaint;
    int tip_states = getTipLikelihoodMatrixRank();
    if (state < 0) {
        complaint << "Cannot calculate tip likelihoods for state " << state << ".";
    }
    else if (num_states <= state && tip_likelihood_formula.empty()) {
        complaint << "Cannot calculate tip likelihoods for state " << state
            << " as there are only " << num_states << " states.";
    }
    else if (tip_states <= state && tip_likelihood_formula.empty()) {
        complaint << "Cannot calculate tip likelihoods for state " << state
            << " as tip likelihoods were provided"
            << " only for " << tip_states << " states.";
    }
    if (!complaint.str().empty()) {
        outError(complaint.str());
    }
    StrVector dummyRow;
    const StrVector& expr_row = tip_likelihood_expressions.empty()
        ? dummyRow : tip_likelihood_expressions[state];

    typedef ModelExpression::InterpretedExpression Interpreter;

    forceAssign("row", state+1);
    ModelVariable& column_var = forceAssign("column", 0);
    for (int column = 0; column < num_states; ++column) {
        column_var.setValue(column+1);
        std::string expr_string;
        if (column < expr_row.size()) {
            expr_string = expr_row[column];
        }
        if (expr_string.empty() && !tip_likelihood_formula.empty()) {
            expr_string = tip_likelihood_formula;
        }
        if (expr_string.empty()) {
            likelihoods[column] = (column == state) ? 1.0 : 0.0;
        }
        else {
            try {
                Interpreter interpreter(*this, expr_string);
                likelihoods[column] = interpreter.evaluate();
            }
            catch (ModelExpression::ModelException& x) {
                std::stringstream msg;
                msg << "Error parsing expression"
                    << " for tip likelihood matrix entry"
                    << " for (0-based) row " << state << ","
                    << " and (0-based) column " << column << ":\n"
                    << x.getMessage();
                outError(msg.str());
            }
        }
    }
}

int ModelInfoFromYAMLFile::getRateMatrixRank() const {
    return rate_matrix_rank;
}

const std::string& ModelInfoFromYAMLFile::getRateMatrixFormula() const {
    return rate_matrix_formula;
}

RateHeterogeneity* ModelInfoFromYAMLFile::getSpecifiedRateModel(PhyloTree* tree) {
    ASSERT (!is_rate_model);    //rate models don't *have* rate models, 
                                //they are rate models
    if (!hasSpecifiedRateModel()) {
        return nullptr;
    }
    return specified_rate_model_info->getRateHeterogeneity(tree);
}

RateHeterogeneity* ModelInfoFromYAMLFile::getRateHeterogeneity(PhyloTree* tree) {
    ASSERT(is_rate_model);

    if (false) {
        return new YAMLRateMeyerDiscrete(tree, tree, *this);
        //return new YAMLRateMeyerHaeseler(tree, tree, *this);
    }
    if (isKategoryModel()) {
        return new YAMLRateKategory(tree, tree, *this);
    }
    bool isInvar = isInvariantModel();
    if (hasRateHeterotachy()) {
        if (isInvar) {
            return new YAMLRateHeterotachyInvar(tree, tree, *this);
        } else {
            return new YAMLRateHeterotachy(tree, tree, *this);
        }
    }
    //Note:gamma rate model is treated as a *kind* of free rate model
    if (isInvar) {
        if (getNumberOfRateCategories()==0) {
            return new YAMLRateInvar(tree, tree, *this);
        } else {
            return new YAMLRateFreeInvar(tree, tree, *this);
        }
    }
    else {
        return new YAMLRateFree(tree, tree, *this);
    }
}

std::string ModelInfoFromYAMLFile::getParameterList(ModelParameterType param_type) const {
    std::stringstream list;
    appendParameterList(param_type, list);
    return list.str();
}

bool ModelInfoFromYAMLFile::appendParameterList(ModelParameterType param_type,
    std::stringstream& list) const {
    const char* separator       = "";
    bool        anything_listed = false;
    for (auto p : parameters) {
        if (p.type == param_type) {
            if (p.is_subscripted) {
                for (int sub = p.minimum_subscript;
                    sub <= p.maximum_subscript; ++sub) {
                    std::string var_name = p.getSubscriptedVariableName(sub);
                    auto it = variables.find(var_name);
                    if (it != variables.end()) {
                        const ModelVariable& v = it->second;
                        list << separator << var_name << "=" << v.getValue();
                        if (v.isFixed()) {
                            list << "(*)";
                        }
                        separator = ", ";
                        anything_listed = true;
                    }
                    else {
                        std::stringstream complaint;
                        complaint << "Variable " << var_name << " not found ";
                        outError(complaint.str());
                    }
                }
            }
            else {
                const std::string& var_name = p.name;
                auto it = variables.find(var_name);
                if (it != variables.end()) {
                    const ModelVariable& v = it->second;
                    list << separator << p.name << "=" << v.getValue();
                    if (v.isFixed()) {
                        list << "(*)";
                    }
                    separator = ", ";
                    anything_listed = true;
                }
                else {
                    std::stringstream complaint;
                    complaint << "Variable " << var_name << " not found ";
                    outError(complaint.str());
                }
            }
        }
    }

    if (this->subtree_models != nullptr) {
        for (ModelInfoFromYAMLFile* model : *subtree_models) {
            std::stringstream sub_list;
            sub_list << separator << model->getName();
            sub_list << "={";
            if (model->appendParameterList(param_type, sub_list)) {
                list << sub_list.str() << "}";
                separator = ", ";
                anything_listed = true;
            }
        }
    }

    //
    //Todo: Decide, is this right?  e.g. Rate parameters
    //      for a JC+GTR model might read "GTR={r(1)=4, r(2)=3, ...}"
    //
    if (this->mixed_models != nullptr) {
        for (ModelInfoFromYAMLFile* model : *mixed_models) {
            std::stringstream sub_list;
            sub_list << separator << model->getName();
            sub_list << "={";
            if (model->appendParameterList(param_type, sub_list)) {
                list << sub_list.str() << "}";
                separator = ", ";
                anything_listed = true;
            }
        }
    }
    return anything_listed;
}

const std::string& ModelInfoFromYAMLFile::getRateMatrixExpression
(int row, int col) const {
    if (rate_matrix_expressions.empty()) {
        return rate_matrix_formula;
    }
    ASSERT(0 <= row);
    ASSERT(row < rate_matrix_expressions.size());
    ASSERT(0 <= col);
    ASSERT(col < rate_matrix_expressions.size());
    const StrVector& matrix_row = rate_matrix_expressions[row];
    if (col < matrix_row.size()) {
        return matrix_row[col];
    }
    //
    //possibly a triangular matrix?!
    //Todo: Don't the stationary frequencies have
    //to be taken into account, here?
    //
    std::swap(row, col);
    const StrVector& other_matrix_row = rate_matrix_expressions[row];
    ASSERT(col < other_matrix_row.size());
    return other_matrix_row[col];
}

const std::string& ModelInfoFromYAMLFile::getName() const {
    return model_name;
}

void ModelInfoFromYAMLFile::appendTo(const std::string& append_me,
                                     const char* with_sep,
                                     std::string& to_me) {
    if (append_me.empty()) {
        return;
    }
    if (to_me.empty()) {
        to_me = append_me;
        return;
    }
    to_me += with_sep;
    to_me += append_me;
}

bool ModelInfoFromYAMLFile::checkIntConsistent(const std::string& value_source,
                                               const char* int_name,
                                               int new_value,
                                               int &old_value,
                                               std::stringstream& complaint) {
    if (old_value==0) {
        //Wasn't set, take new value
        old_value = new_value;
        return true;
    }
    if (new_value==0) {
        //Wasn't modified, take old value
        return true;
    }
    if (new_value==old_value) {
        //Old value and new value agreed
        return true;
    }
    complaint << "Cannot have " << int_name << " of both " << old_value
              << " and " << new_value << " (as per " << value_source << "). "; 
    return false;
}

void ModelInfoFromYAMLFile::inheritModelRateDistributions
        (const ModelInfoFromYAMLFile& mummy) {
    if (mummy.rate_distribution.empty()) {
        return;
    }
    //Comma-Separated List set union
    if (rate_distribution.empty()) {
        rate_distribution = mummy.rate_distribution;
        return;
    } 
    std::set<std::string> dist_set;
    for (auto child_distribution : split_string(rate_distribution,",") ) {
        dist_set.insert(child_distribution);
    }
    for (auto mummy_distribution : split_string(mummy.rate_distribution,",") ) {
        if (dist_set.insert(mummy_distribution).second) {
            rate_distribution += ",";
            rate_distribution += mummy_distribution;
        }
    }
}

void ModelInfoFromYAMLFile::inheritModel(const ModelInfoFromYAMLFile& mummy,
                                         LoggingTarget* logging_target) {
    if (this==&mummy) {
        //Ouch.  You're not allowed to go R+R+R... Or similar.
        //It'd probably *work* but it's not a case that programmers
        //should have to worry about.
        TREE_LOG_LINE(*logging_target, VerboseMode::VB_MIN,
                      "Model " << getName() << " cannot inherit from itself.");
        return;
    }
    inheritModelRateDistributions(mummy);

    if (is_rate_model) {
        //Todo: Recalculate, how many categories there are!
        //Or something.  I'm not sure what needs to happen here.
        std::set<std::string> rate_vars;
        mummy.addNamesOfVariablesOfTypeToSet(ModelParameterType::RATE, rate_vars);
        addNamesOfVariablesOfTypeToSet(ModelParameterType::RATE, rate_vars);
        int rate_cats = static_cast<int>(rate_vars.size());
        forceAssign("categories", rate_cats );
        model_name += "+";
        model_name += mummy.getName();
    }

    std::stringstream complaint;
    inheritModelProperties(mummy, complaint);
    inheritModelParameters(mummy, complaint, logging_target);
    inheritModelVariables (mummy, complaint);
    inheritModelMatrices  (mummy, complaint);

    ModelExpression::ModelException::throwIfNonBlank(complaint);
}

void ModelInfoFromYAMLFile::addNamesOfVariablesOfTypeToSet
        (ModelParameterType type,
         std::set<std::string>& var_names) const {
    for (auto it = variables.begin();
            it!=variables.end(); ++it) {
        if (it->second.getType() == type) {
            var_names.insert(it->first);
        }
    }
}


void ModelInfoFromYAMLFile::inheritModelProperties(const ModelInfoFromYAMLFile& mummy,
                                                   std::stringstream& complaint) {
    appendTo(mummy.citation, ", ", citation);
    appendTo(mummy.DOI,      ", ", DOI);
    appendTo(mummy.url,      ", ", url);
    if (data_type_name=="") {
        data_type_name=mummy.data_type_name;
        sequence_type=mummy.sequence_type;
    } 
    else if (mummy.data_type_name!="") {
        if (mummy.data_type_name != data_type_name) {
            complaint << "Cannot have data type of " << data_type_name << " and "
                      << mummy.data_type_name 
                      << " (as per " << mummy.model_name << "). ";
        }
    }
    //Todo: what about reversible?!
    checkIntConsistent(mummy.model_name, "states", 
                       mummy.num_states, num_states, complaint);

    checkIntConsistent(mummy.model_name, "rate matrix rank", 
                       mummy.rate_matrix_rank, rate_matrix_rank, complaint);
    if (!mummy.rate_matrix_expressions.empty()) {
        rate_matrix_expressions = mummy.rate_matrix_expressions;
    }
    if (!mummy.rate_matrix_formula.empty()) {
        rate_matrix_formula = mummy.rate_matrix_formula;
    }

    checkIntConsistent(mummy.model_name, "tip likelihood matrix rank", 
                       mummy.tip_likelihood_rank, tip_likelihood_rank,
                       complaint);
    if (!mummy.tip_likelihood_expressions.empty()) {
        tip_likelihood_expressions = mummy.tip_likelihood_expressions;
    }
    if (!mummy.tip_likelihood_formula.empty()) {
        tip_likelihood_formula = mummy.tip_likelihood_formula;
    }
    if (mummy.frequency_type != StateFreqType::FREQ_UNKNOWN) {
        frequency_type = mummy.frequency_type;
    }
}

void ModelInfoFromYAMLFile::inheritModelParameters(const ModelInfoFromYAMLFile& mummy,
                                                   std::stringstream& complaint,
                                                   LoggingTarget* logging_target) {
    for (const YAMLFileParameter& mummy_param : mummy.parameters) {
        //Ee-uw.  Add parameter
        const YAMLFileParameter* child_param = findParameter(mummy_param.name);
        if (child_param == nullptr) {
            addParameter(mummy_param, logging_target);
            continue;
        }
        if (child_param->type != mummy_param.type) {
            complaint << "Cannot override type of parameter " << mummy_param.name 
                      << " to " << child_param->type_name 
                      << " from " << mummy_param.type_name
                      << " (as per " << mummy.getName() << "). ";
            continue;
        }
        if (child_param->is_subscripted != mummy_param.is_subscripted) {
            const char* child_state = child_param->is_subscripted
                                    ? "subscripted" : "un-subscripted";
            const char* mummy_state = mummy_param.is_subscripted
                                    ? "subscripted" : "un-subscripted";
            complaint << "Cannot override parameter " << mummy_param.name 
                      << " to " << child_state << " from " << mummy_state
                      << " (as per " << mummy.getName() << "). ";
            continue;
        } 
        else if (child_param->is_subscripted) {
            if (child_param->minimum_subscript != mummy_param.minimum_subscript ||
                child_param->maximum_subscript != mummy_param.maximum_subscript) {
                //Todo: What if the subscript *expressions* are the SAME, but
                //      the subscripts are different?
                //      Evaluating mummy_param.subscript_expression
                complaint << "Cannot change subscript range"
                          << " of parameter " << mummy_param.name
                          << " to " << child_param->minimum_subscript
                          << " .." << child_param->maximum_subscript
                          << " from " << mummy_param.minimum_subscript
                          << ".." << mummy_param.maximum_subscript
                          << " (as per " << mummy.getName() << "). ";
                continue;
            }
            //Todo: What if mummy *constrained* some of the expressions?
            //Todo: What if mummy *defaulted* some of the expressions?
        }
    }
}

void ModelInfoFromYAMLFile::inheritModelVariables(const ModelInfoFromYAMLFile& mummy,
                                                  std::stringstream& complaint) {
    for (auto mapping : mummy.variables) {
        std::string var_name = mapping.first;
        const ModelVariable& mummy_var = mapping.second;
        if (hasVariable(mapping.first)) {
            ModelVariable& child_var = variables[var_name];
            if (mummy_var.isFixed() && !child_var.isFixed()) {
                //Todo: Decide how to handle this:
                //      mummy constrained the variable, child didn't.
                //      Probably need to do *something*.  Warning, maybe?
            }
        } else {
            //Only add variables we don't already have
            variables[var_name] = mapping.second;
        }
    }
    //Todo: What about mixed_models?
    for (auto prop_mapping: mummy.string_properties) {
        string_properties[prop_mapping.first] = prop_mapping.second;
    }
    getVariableNamesByPosition();
}

void ModelInfoFromYAMLFile::inheritModelMatrices(const ModelInfoFromYAMLFile& mummy,
                                                 std::stringstream& complaint) {
    if (0<mummy.rate_matrix_rank) {
        rate_matrix_rank           = mummy.rate_matrix_rank;
        rate_matrix_expressions    = mummy.rate_matrix_expressions;
    }
    if (!mummy.rate_matrix_formula.empty()) {
        rate_matrix_formula        = mummy.rate_matrix_formula;
    }
    if (0<mummy.tip_likelihood_rank) {
       tip_likelihood_rank        = mummy.tip_likelihood_rank;
        tip_likelihood_expressions = mummy.tip_likelihood_expressions;
    }
    if (!mummy.tip_likelihood_formula.empty()) {
        tip_likelihood_formula     = mummy.tip_likelihood_formula;
    }
}

const std::string& ModelInfoFromYAMLFile::getOptimizationAlgorithm() const {
    return opt_algorithm;
}

int ModelInfoFromYAMLFile::getNumberOfRateCategories() const {
    int count = 0;
    for (auto v: variables) {
        if (v.second.getType() == ModelParameterType::RATE) {
            ++count;
        }
    }
    return count;
}

int ModelInfoFromYAMLFile::getNumberOfVariableRates() const {
    int count = 0;
    for (auto v: variables) {
        if (v.second.getType() == ModelParameterType::RATE) {
            if (!v.second.isFixed()) {
                ++count;
            }
        }
    }
    return count;
}

int ModelInfoFromYAMLFile::getNumberOfVariableShapes() const {
    int count = 0;
    for (auto v: variables) {
        if (v.second.getType() == ModelParameterType::SHAPE) {
            if (!v.second.isFixed()) {
                ++count;
            }
        }
    }
    return count;
}

int ModelInfoFromYAMLFile::getNumberOfProportions() const {
    int count = 0;
    for (auto v: variables) {
        auto t = v.second.getType();
        if (t == ModelParameterType::PROPORTION) {
            ++count;
        } else if (t == ModelParameterType::INVARIANT_PROPORTION) {
            ++count;
        }
    }
    return count;
}

int ModelInfoFromYAMLFile::getNumberOfInvariantProportions() const {
    int count = 0;
    for (auto v: variables) {
        auto t = v.second.getType();
        if (t == ModelParameterType::INVARIANT_PROPORTION) {
            ++count;
        }
    }
    return count;
}

int ModelInfoFromYAMLFile::getNumberOfVariableProportions() const {
    int count = 0;
    for (auto v: variables) {
        auto t = v.second.getType();
        if ( t == ModelParameterType::PROPORTION ||
             t == ModelParameterType::INVARIANT_PROPORTION ) {
            if (!v.second.isFixed()) {
                ++count;
            }
        }
    }
    return count;
}


bool ModelInfoFromYAMLFile::acceptParameterList(Params& params,
                                                std::string parameter_list,
                                                LoggingTarget* logging_target) {
    typedef ModelExpression::InterpretedExpression Interpreter;
    typedef ModelExpression::Expression            Expression;
    typedef ModelExpression::Assignment            Assignment;
    typedef ModelExpression::Variable              Variable;
    typedef ModelExpression::ModelException        Exception;
                                                
    trimString(parameter_list);
    if (startsWith(parameter_list, "{") &&
        endsWith(parameter_list, "}")) {
        auto len = parameter_list.length();
        parameter_list = parameter_list.substr(1, len-2);
    }
    size_t    param_list_length = parameter_list.length();
    size_t    i                 = 0;
    std::vector<Interpreter*> expr_list;
    while (i<param_list_length) {
        size_t      j     = findEndOfParameter(parameter_list, 
                                               param_list_length, i);
        std::string param = parameter_list.substr(i, j-i);
        expr_list.push_back(new Interpreter(*this, param));
        i = j + 1;
    }
    bool fix      = !params.optimize_from_given_params;
    getVariableNamesByPosition();
    try {
        int  position = 0;
        for (i=0; i<expr_list.size(); ++i) {
            Interpreter* ix = expr_list[i];
            Expression*  ex = ix->expression();
            if (ex->isAssignment()) {
                Assignment*    a        = dynamic_cast<Assignment*>(ex);
                Variable*      xv       = a->getTargetVariable();
                string         var_name = xv->getName();
                Expression*    x        = a->getExpression();

                assign( var_name, x, fix, "by name", logging_target);
            } else {
                string         var_name = getVariableNameByPosition(position);

                assign( var_name, ex, fix, "by position", logging_target);
                ++position;
            }
            delete ix;
            expr_list[i] = 0;
        }   
    }
    catch (Exception& problem) {
        std::stringstream complaint;
        complaint << "An error occurred parsing parameter list"
                  << " ... " << parameter_list << " ...:\n"
                  << problem.getMessage();
        outError(complaint.str());
    }
    changeSubscriptRangesOfParameters(logging_target);
    return !expr_list.empty();
}

size_t ModelInfoFromYAMLFile::findEndOfParameter(const std::string& parameter_list,
                                                 size_t param_list_length,
                                                 size_t i) const {
    int    bracket_depth = 0;
    size_t j             = i;
    for (;j<param_list_length; ++j) {
        char ch = parameter_list[j];
        if (ch==',') {
            if (0==bracket_depth) {
                break;
            }
        }
        else if (ch=='(') {
            ++bracket_depth;
        }
        else if (ch==')') {
            --bracket_depth;
        }
    }
    return j;
}

void ModelInfoFromYAMLFile::changeSubscriptRangesOfParameters
        (LoggingTarget* logging_target) {
    typedef ModelExpression::InterpretedExpression Interpreter;
    typedef ModelExpression::RangeOperator         Range;                                 

    for (YAMLFileParameter& p : parameters) {
        if (p.is_subscripted && !p.subscript_expression.empty()) {
            Interpreter expr(*this, p.subscript_expression);
            if (expr.expression()->isRange()) {
                Range* op = dynamic_cast<Range*>(expr.expression());
                changeParameterSubscriptRange(op->getIntegerMinimum(), 
                                              op->getIntegerMaximum(),p,
                                              logging_target);
            }
        }
    }
}

void ModelInfoFromYAMLFile::changeParameterSubscriptRange
        (int new_min, int new_max, YAMLFileParameter& param,
         LoggingTarget* logging_target) {
    if (new_max < new_min) {
        std::stringstream complaint;
        complaint << "An error occurred, setting the subscript "
                  << " range for parameter " << param.name
                  << " of model " << getName() << ": lower "
                  << " bound (" << new_min <<") exceeded "
                  << " the upper bound (" << new_max << ").";
        outError(complaint.str());
        return;
    }
    int old_min = param.minimum_subscript;
    int old_max = param.maximum_subscript;
    param.minimum_subscript = new_min;
    param.maximum_subscript = new_max;
    for (int i=new_min; i<=new_max; ++i) {
        if (i<old_min || old_max<i) {
            setSubscriptedVariable(param, i, logging_target);
        } 
    }
    for (int i=old_min; i<=old_max; ++i) {
        if (i<new_min || new_max<i) {
            removeSubscriptedVariable(param, i, logging_target);
        }
    }
}

MapOfModels& ModelInfoFromYAMLFile::getMixedModels() {
    ASSERT(mixed_models!=nullptr);
    return *mixed_models;
}

const MapOfModels& ModelInfoFromYAMLFile::getMixedModels() const {
    ASSERT(mixed_models!=nullptr);
    return *mixed_models;
}

MapOfModels& ModelInfoFromYAMLFile::getSubtreeModels() {
    ASSERT(subtree_models!=nullptr);
    return *subtree_models;
}

const MapOfModels& ModelInfoFromYAMLFile::getSubtreeModels() const {
    ASSERT(subtree_models!=nullptr);
    return *subtree_models;
}

void ModelInfoFromYAMLFile::copyVariablesFrom(const ModelInfoFromYAMLFile* original) {
    for (auto it = original->variables.begin();
         it != original->variables.end(); ++it) {
        std::string    var_name  = it->first;
        const ModelVariable& var_value = it->second;
        variables[var_name]     = var_value;
    }
}

void ModelInfoFromYAMLFile::setParentModel(ModelInfoFromYAMLFile* parent) {
    parent_model = parent;
}

ModelVariable& ModelInfoFromYAMLFile::assign(const std::string& var_name,  
                                   ModelExpression::Expression *x,
                                   bool fix, const char* how,
                                   LoggingTarget* logging_target) {
    typedef ModelExpression::RangeOperator         Range;
    typedef ModelExpression::ModelException        Exception;                                    
    if (x->isRange()) {
        Range*         range     = dynamic_cast<Range*>(x);
        ModelVariable& mv        = variables[var_name];
        double         min_value = range->getMinimum();
        double         max_value = range->getMaximum();
        if (max_value<min_value) {
            std::stringstream complaint;
            complaint << "Range supplied for " << var_name 
                        << " is invalid.  Requested minimum " << min_value
                        << " is greater than requested maximum " << max_value;
            throw Exception(complaint.str());

        } 
        mv.setMinimum( min_value );
        mv.setMaximum( max_value );
        mv.constrainValueToRange();
        if (min_value==max_value) {
            mv.markAsFixed(); //Todo: permanently, though?
        }
        TREE_LOG_LINE(*logging_target, YAMLVariableVerbosity,
                        "Set " << var_name 
                        << " range to " << min_value 
                        << ".." << max_value 
                        << " " << how << "." );
        return mv;
    }  
    else {
        double         setting  = x->evaluate();
        ModelVariable& mv       = assign(var_name, setting);
        const char* verb = "Set ";
        if (fix && ! x->isEstimate()) {
            mv.markAsFixed();
            verb = "Set and fixed ";
        }
        TREE_LOG_LINE(*logging_target, YAMLVariableVerbosity,
                        verb << var_name 
                        << " to " << setting 
                        << " " << how << ".");
        return mv;
    }
}

void ModelInfoFromYAMLFile::writeInfo(const char* caption,
                                      ModelParameterType param_type,
                                      std::ostream& out ) const {
    auto params = getParameterList(param_type);
    if (!params.empty()) {
        out << caption << ": " << params << std::endl;
    }
}

void ModelInfoFromYAMLFile::saveToCheckpoint(Checkpoint* checkpoint) const {
    //Todo: what about variables of linked models?
    int variable_count = static_cast<int>(variables.size());
    CKP_SAVE(variable_count);
    checkpoint->startList(variable_count);
    for (auto it=variables.begin(); it!=variables.end(); ++it) {
        checkpoint->addListElement();
        checkpoint->put("name",  it->first);
        checkpoint->put("type",  it->second.getTypeName());
        checkpoint->put("value", it->second.getValue());
        checkpoint->put("fixed", it->second.isFixed());
    }
    checkpoint->endList();
}

void ModelInfoFromYAMLFile::restoreFromCheckpoint(Checkpoint* checkpoint) {
    //Todo: what about variables of linked models?
    int variable_count = 0;
    CKP_RESTORE(variable_count);
    checkpoint->startList(variable_count);
    for (int var_num=0; var_num<variable_count; ++var_num) {
        checkpoint->addListElement();
        std::string name;
        std::string type;
        double      value;
        bool        fixed;
        CKP_RESTORE(name);
        CKP_RESTORE(type);
        CKP_RESTORE(value);
        CKP_RESTORE(fixed);
        ModelVariable& var = forceAssign(name, value);
        var.setTypeName(type);
        if (fixed) {
            var.markAsFixed();
        }
    }
    checkpoint->endList();
}

ASCType ModelInfoFromYAMLFile::extractASCType(std::string& leftover_name) const {
    ASCType type;
    if (checkAscertainmentBiasCorrection(false, type)) {
        return type;
    } else {
        return ASCType::ASC_NONE;
    }
}

bool ModelInfoFromYAMLFile::checkAscertainmentBiasCorrection
        (bool warnIfInvalid, ASCType &type) const {
    std::string asc_name;
    if (!hasStringProperty(PROPERTY_NAME_ASC, asc_name)) {
        return false;
    }
    asc_name = string_to_lower(asc_name);
    bool unrecognized = false;
    bool is_holder    = false;
    if (contains(asc_name, "holder")) {
        is_holder     = true;
    }
    else {
        unrecognized |= !contains(asc_name, "lewis");
    }    
    bool is_variant = false;
    if (contains(asc_name, "variant")) {
        is_variant    = true;
    }
    else {
        unrecognized |= !contains(asc_name, "informative");
    }
    type = is_holder
         ? (is_variant ? ASCType::ASC_VARIANT_MISSING 
                       : ASCType::ASC_INFORMATIVE_MISSING )
         : (is_variant ? ASCType::ASC_VARIANT
                       : ASCType::ASC_INFORMATIVE );
    if (unrecognized && warnIfInvalid) {
        outWarning("Ascertainment bias correction"
                " (" + asc_name + ") for model " + getName() +
                " not recognized.");
    }
    return !unrecognized;
}

bool ModelInfoFromYAMLFile::hasRateHeterotachy() const {
    return contains(rate_distribution,"heterotachy");
}

bool ModelInfoFromYAMLFile::hasSpecifiedRateModel() const {
    return specified_rate_model_info!=nullptr;
}

bool ModelInfoFromYAMLFile::hasAscertainmentBiasCorrection() const {
    ASCType type;
    return checkAscertainmentBiasCorrection(false, type);
}

bool isYAMLRateHeterotachyModel
        (Params& params, const std::string& model_name, 
         int& num_mixlen) {
    if (model_name.empty()) {
        return false; //no model
    }
    auto yaml_path = params.yaml_model_file;
    if (yaml_path.empty()) {
        //std::cout << "XX NO YAML file" << std::endl;
        return false; //not using YAML model file
    }
    LoggingTarget dummy_logging_target; //No tree yet!
    ModelListFromYAMLFile yaml_list;
    yaml_list.loadFromFile(params, yaml_path.c_str(), 
                           &dummy_logging_target);

    std::string subst_model_name = split_string(model_name, "+")[0];

    if (!yaml_list.isSubstitutionModelNameRecognized(subst_model_name)) {
        //std::cout << "XX Name not recognized: " << subst_model_name << std::endl;
        return false; //not a YAML model
    }
    return yaml_list.isRateHeterotachyRequired(params, model_name, num_mixlen,
                                            &dummy_logging_target);
}

void ModelInfoFromYAMLFile::addCladeName
        (const std::string& clade_name) {
    std::string lower_name = string_to_lower(clade_name);
    if (distinct_clade_names.insert(lower_name).second) {
        clade_names.push_back(lower_name);
    } else {
        std::stringstream complaint;
        complaint << "Could not add the same clade, " << clade_name
                  << ", twice, to model " << getQualifiedName()
                  << ".";
        throw ModelExpression::ModelException(complaint.str());
    }
}

const StrVector& ModelInfoFromYAMLFile::getCladeNames() const {
    //std::cout << "Clades for " << getQualifiedName()
    //          << " are " << clade_names.join(",") << std::endl;
    return clade_names;
}