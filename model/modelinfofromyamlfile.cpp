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
#include <tree/phylotree.h> //for TREE_LOG_LINE macro

VerboseMode YAMLModelVerbosity     = VerboseMode::VB_MIN;
VerboseMode YAMLVariableVerbosity  = VerboseMode::VB_MAX;
VerboseMode YAMLFrequencyVerbosity = VerboseMode::VB_MAX;
VerboseMode YAMLMatrixVerbosity    = VerboseMode::VB_MAX;

YAMLFileParameter::YAMLFileParameter()
    : is_subscripted(false), minimum_subscript(0), maximum_subscript(0)
    , type(ModelParameterType::OTHER), value(0.0) {
}

std::string YAMLFileParameter::getSubscriptedVariableName(int subscript) const {
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


ModelVariable::ModelVariable() : type(ModelParameterType::OTHER)
, value(0), is_fixed(false) {
}

ModelVariable::ModelVariable(ModelParameterType t,
    const ModelParameterRange& r,
    double v)
    : range(r), type(t), value(v), is_fixed(false) {
}

void ModelVariable::setValue(double v) {
    value = v;
}

void ModelVariable::markAsFixed() {
    is_fixed = true;
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
    : is_modifier_model(false)
    , sequence_type(SeqType::SEQ_UNKNOWN), num_states(0)
    , reversible(false), rate_matrix_rank(0)
    , tip_likelihood_rank(0)
    , frequency_type(StateFreqType::FREQ_UNKNOWN)
    , mixed_models(nullptr), linked_models(nullptr) {
}

void ModelInfoFromYAMLFile::copyMixedAndLinkedModels(const ModelInfoFromYAMLFile& rhs) {
    delete mixed_models;
    if (rhs.mixed_models != nullptr) {
        mixed_models = new MapOfModels(*rhs.mixed_models);
    }
    delete linked_models;
    if (rhs.linked_models != nullptr) {
        linked_models = new MapOfModels(*rhs.linked_models);
    }
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const ModelInfoFromYAMLFile& rhs)
    : model_name(rhs.model_name), model_file_path(rhs.model_file_path)
    , parent_model_name(rhs.parent_model_name)
    , is_modifier_model(rhs.is_modifier_model)
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
    , variables(rhs.variables)
    , mixed_models(nullptr), linked_models(nullptr) {
    copyMixedAndLinkedModels(rhs);
}

ModelInfoFromYAMLFile& ModelInfoFromYAMLFile::operator=(const ModelInfoFromYAMLFile& rhs) {
    if (&rhs != this) {
        model_name                 = rhs.model_name;
        model_file_path            = rhs.model_file_path;
        parent_model_name          = rhs.parent_model_name;
        is_modifier_model          = rhs.is_modifier_model;
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
        copyMixedAndLinkedModels(rhs);
    }
    return *this;
}

ModelInfoFromYAMLFile::ModelInfoFromYAMLFile(const std::string& path)
    : model_file_path(path), is_modifier_model(false)
    , sequence_type(SeqType::SEQ_UNKNOWN)
    , num_states(0), reversible(true)
    , rate_matrix_rank(0), tip_likelihood_rank(0)
    , frequency_type(StateFreqType::FREQ_UNKNOWN)
    , mixed_models(nullptr), linked_models(nullptr) {
}

ModelInfoFromYAMLFile::~ModelInfoFromYAMLFile() {
    delete mixed_models;
    mixed_models = nullptr;
    delete linked_models;
    linked_models = nullptr;
}

bool ModelInfoFromYAMLFile::isMixtureModel() const {
    return mixed_models != nullptr;
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

void ModelInfoFromYAMLFile::updateName(const std::string& name) {
    model_name = name;
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

ModelInfoFromYAMLFile::MapOfModels::const_iterator
ModelInfoFromYAMLFile::findMixedModel(const std::string& name) const {
    auto it = mixed_models->find(name);
    if (it == mixed_models->end()) {
        std::stringstream complaint;
        complaint << "Could not evaluate mixedl model name " << name
            << " for model " << getLongName();
        throw ModelExpression::ModelException(complaint.str());
    }
    return it;
}

ModelInfoFromYAMLFile::MapOfModels::iterator
ModelInfoFromYAMLFile::findMixedModel(const std::string& name) {
    auto it = mixed_models->find(name);
    if (it == mixed_models->end()) {
        std::stringstream complaint;
        complaint << "Could not evaluate mixed model name " << name
            << " for model " << getLongName();
        throw ModelExpression::ModelException(complaint.str());
    }
    return it;
}

void ModelInfoFromYAMLFile::setNumberOfStatesAndSequenceType(int requested_num_states,
                                                             PhyloTree* report_to_tree) {
    if (requested_num_states != 0) {
        num_states = requested_num_states;
    }
    if (num_states == 0) {
        num_states = 4;
    }
    if (!data_type_name.empty()) {
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "Data Type Name is" << data_type_name);
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
    TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                  "Number of states is " << num_states);

    forceAssign("num_states", num_states);
    forceAssign("numStates", num_states);
}

double ModelInfoFromYAMLFile::evaluateExpression(std::string& expr,
    std::string context) {
    const char* verb = "parsing";
    try {
        ModelExpression::InterpretedExpression interpreter(*this, expr);
        verb = "evaluating";
        return interpreter.evaluate();
    }
    catch (ModelExpression::ModelException x) {
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


bool ModelInfoFromYAMLFile::hasVariable(const char* name) const {
    //Todo: look up linked_models too
    if (hasDot(name) && mixed_models != nullptr ) {
        std::string sub_model_name;
        const char* var_name = nullptr;
        breakAtDot(name, sub_model_name, var_name);
        auto it = findMixedModel(sub_model_name);
        return it->second.hasVariable(var_name);
    }
    return variables.find(name) != variables.end();
}

bool ModelInfoFromYAMLFile::hasVariable(const std::string& name) const {
    return hasVariable(name.c_str());
}

double ModelInfoFromYAMLFile::getVariableValue(const char* name) const {
    auto found = variables.find(name);
    if (found == variables.end()) {
        //Todo: look up linked_models too
        if (hasDot(name) && mixed_models != nullptr) {
            std::string sub_model_name;
            const char* var_name = nullptr;
            breakAtDot(name, sub_model_name, var_name);
            auto it = findMixedModel(sub_model_name);
            return it->second.getVariableValue(var_name);
        }
        return 0.0;
    }
    return found->second.getValue();
}

double ModelInfoFromYAMLFile::getVariableValue(const std::string& name) const {
    return getVariableValue(name.c_str());
}

bool ModelInfoFromYAMLFile::isVariableFixed(const std::string& name) const {
    auto found = variables.find(name);
    if (found == variables.end()) {
        //Todo: look up linked_models too
        if (hasDot(name.c_str()) && mixed_models != nullptr) {
            std::string sub_model_name;
            const char* var_name = nullptr;
            breakAtDot(name.c_str(), sub_model_name, var_name);
            auto it = findMixedModel(sub_model_name);
            return it->second.isVariableFixed(var_name);
        }
        return false;
    }
    return found->second.isFixed();
}

void ModelInfoFromYAMLFile::addParameter(const YAMLFileParameter& p) {
    bool replaced = false;
    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        if (it->name == p.name) {
            *it = p;
            replaced = true;
            break;
        }
    }
    if (!replaced) {
        parameters.emplace_back(p);
    }
    if (p.is_subscripted) {
        for (int i = p.minimum_subscript; i <= p.maximum_subscript; ++i) {
            std::string var_name = p.getSubscriptedVariableName(i);
            variables[var_name] = ModelVariable(p.type, p.range, p.value);
        }
    }
    else {
        variables[p.name] = ModelVariable(p.type, p.range, p.value);
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

void ModelInfoFromYAMLFile::setBounds(int param_count, double* lower_bound,
    double* upper_bound, bool* bound_check) const {
    int i = 1; //Rate parameter numbering starts at 1, see ModelMarkov
    for (auto p : parameters) {
        if (p.type == ModelParameterType::RATE) {
            for (int sub = p.minimum_subscript;
                sub <= p.maximum_subscript; ++sub) {
                ASSERT(i <= param_count);
                lower_bound[i] = p.range.first;
                upper_bound[i] = p.range.second;
                bound_check[i] = false;
                ++i;
            }
        }
    }
}

void ModelInfoFromYAMLFile::updateVariables(const double* updated_values,
    int first_freq_index,
    int param_count) {
    int i = 1; //Rate parameter numbering starts at 1, see ModelMarkov
    ModelParameterType supported_types[] = {
        ModelParameterType::RATE, ModelParameterType::FREQUENCY };
    //FREQUENCY must be after RATE.
    //Todo: Where do weight parameters go?    Esepecially in mixture
    //      models
    for (auto param_type : supported_types) {
        if (param_type == ModelParameterType::FREQUENCY) {
            i = first_freq_index;
        }
        for (auto p : parameters) {
            if (p.type == param_type) {
                for (int sub = p.minimum_subscript;
                    sub <= p.maximum_subscript; ++sub) {
                    if (i <= param_count) {
                        std::string var_name = p.getSubscriptedVariableName(sub);
                        double      var_value = updated_values[i];
                        if (!this->variables[var_name].isFixed()) {
                            this->variables[var_name].setValue(var_value);
                        }
                        ++i;
                    }
                    //For the last frequency parameter, if there is one,
                    //i will be equal to ( param_count + 1 ).
                }
            }
        }
    }
}

void ModelInfoFromYAMLFile::logVariablesTo(PhyloTree& report_to_tree) const {
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
    TREE_LOG_LINE(report_to_tree, YAMLModelVerbosity, list);
}

ModelVariable& ModelInfoFromYAMLFile::assign(const std::string& var_name,
    double value_to_set) {
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        //Todo: look up linked_models too
        if (hasDot(var_name.c_str()) && mixed_models != nullptr) {
            std::string sub_model_name;
            const char* sub_model_var_name = nullptr;
            breakAtDot(var_name.c_str(), sub_model_name,
                sub_model_var_name);
            auto it = findMixedModel(sub_model_name);
            return it->second.assign(std::string(sub_model_var_name),
                                     value_to_set);
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
        if (hasDot(var_name.c_str()) && mixed_models != nullptr) {
            std::string sub_model_name;
            const char* sub_model_var_name = nullptr;
            breakAtDot(var_name.c_str(), sub_model_name,
                sub_model_var_name);
            auto it = findMixedModel(sub_model_name);
            return it->second.forceAssign(std::string(sub_model_var_name),
                                          value_to_set);
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
    return forceAssign(var_name, (double)value_to_set);
}

ModelVariable& ModelInfoFromYAMLFile::forceAssign(const char* var_name,
                                                  int value_to_set) {
    return forceAssign(std::string(var_name), (double)value_to_set);
}

const StrVector& ModelInfoFromYAMLFile::getVariableNamesByPosition() const {
    variable_names.clear();
    ModelParameterType supported_types[] = {
        ModelParameterType::RATE, ModelParameterType::FREQUENCY };
    //FREQUENCY must be after RATE.
    //Todo: Where do weight parameters go?    Esepecially in mixture
    //      models
    for (auto param_type : supported_types) {
        for (auto p : parameters) {
            if (p.type == param_type) {
                for (int sub = p.minimum_subscript;
                    sub <= p.maximum_subscript; ++sub) {
                    std::string var_name = p.getSubscriptedVariableName(sub);
                }
            }
        }
    }
    return variable_names;
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

std::string ModelInfoFromYAMLFile::getStringProperty(const char* name,
    const char* default_value) const {
    auto it = string_properties.find(name);
    return it == string_properties.end() ? default_value : it->second;
}

int ModelInfoFromYAMLFile::getNumStates() const {
    //Todo: decide from the data type (or at least from data_type_name!)
    return 4; //but, for now, hardcoded!
}

int  ModelInfoFromYAMLFile::getTipLikelihoodMatrixRank() const {
    return tip_likelihood_rank;
}

void ModelInfoFromYAMLFile::computeTipLikelihoodsForState(int state, int num_states,
    double* likelihoods) {
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
    StrVector& expr_row = tip_likelihood_expressions.empty()
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
            catch (ModelExpression::ModelException x) {
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

std::string ModelInfoFromYAMLFile::getParameterList(ModelParameterType param_type) const {
    std::stringstream list;
    appendParameterList(param_type, list);
    return list.str();
}

void ModelInfoFromYAMLFile::appendParameterList(ModelParameterType param_type,
    std::stringstream& list) const {
    const char* separator = "";
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
                }
                else {
                    std::stringstream complaint;
                    complaint << "Variable " << var_name << " not found ";
                    outError(complaint.str());
                }
            }
        }
    }

    //
    //Todo: Decide, is this right?  e.g. Rate parameters
    //      for a JC+GTR model might read "GTR={r(1)=4, r(2)=3, ...}"
    //
    if (this->mixed_models != nullptr) {
        for (auto it = mixed_models->begin(); it != mixed_models->end(); ++it) {
            list << separator << it->first;
            list << "={";
            it->second.appendParameterList(param_type, list);
            list << "}";
            separator = ", ";
        }
    }
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

void ModelInfoFromYAMLFile::inheritModel(const ModelInfoFromYAMLFile& mummy) {
    std::stringstream complaint;
    appendTo(mummy.citation, ", ", citation);
    appendTo(mummy.DOI,      ", ", DOI);
    appendTo(mummy.url,      ", ", url);
    if (data_type_name=="") {
        data_type_name=mummy.data_type_name;
        sequence_type=mummy.sequence_type;
    } else if (mummy.data_type_name!="") {
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
    for (const YAMLFileParameter& param : mummy.parameters) {
        //Ee-uw.  Add parameter
        const YAMLFileParameter* old_param = findParameter(param.name);
        if (old_param == nullptr) {
            addParameter(param);
            continue;
        }
        if (old_param->type != param.type) {
            complaint << "Cannot change type of parameter " << param.name 
                      << " from " << old_param->type_name 
                      << " to " << param.type_name
                      << " (as per " << mummy.getName() << "). ";
            continue;
        }
        if (old_param->is_subscripted != param.is_subscripted) {
            const char* old_state = old_param->is_subscripted
                                  ? "subscripted" : "un-subscripted";
            const char* new_state = param.is_subscripted
                                  ? "subscripted" : "un-subscripted";
            complaint << "Cannot change parameter " << param.name 
                      << " from " << old_state << " to " << new_state
                      << " (as per " << mummy.getName() << "). ";
            continue;
        } 
        else if (old_param->is_subscripted) {
            if (old_param->minimum_subscript != param.minimum_subscript ||
                old_param->maximum_subscript != param.maximum_subscript) {
                complaint << "Cannot change subscript range"
                          << " of parameter " << param.name
                          << " from " << old_param->minimum_subscript
                          << " .." << old_param->maximum_subscript
                          << " to " << param.minimum_subscript
                          << ".." << param.maximum_subscript
                          << " (as per " << mummy.getName() << "). ";
                continue;
            }            
        }
    }
    if (mummy.frequency_type != StateFreqType::FREQ_UNKNOWN) {
        frequency_type = mummy.frequency_type;
    }
    for (auto mapping : mummy.variables) {
        if (hasVariable(mapping.first)) {
            if (variables[mapping.first].isFixed()) {
                continue;
            }
        }
        variables[mapping.first] = mapping.second;
    }
    //Todo: What about mixed_models?
    for (auto prop_mapping: mummy.string_properties) {
        string_properties[prop_mapping.first] = prop_mapping.second;
    }
    getVariableNamesByPosition();
    std::string problem = complaint.str();
    if (!problem.empty()) {
        throw ModelExpression::ModelException(problem);
    }
}
