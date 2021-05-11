//
// modellistfromyamlfile.cpp
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
#include "modelsubst.h"            //for OPEN_BRACKET and CLOSE_BRACKET
#include "modelfileloader.h"
#include "modelexpression.h"       //for InterpretedExpression

#include "modeldna.h"              //for ModelDNA
#include "modeldnaerror.h"         //for ModelDNAError
#include "modelprotein.h"          //for ModelProtein
#include "modelcodon.h"            //for ModelCodon
#include "modelbin.h"              //for ModelBIN
#include "modelmorphology.h"       //for ModelMorphology

#include <utils/my_assert.h>       //for ASSERT macro
#include <utils/stringfunctions.h> //for convert_int
#include <utils/tools.h>           //for outError

void ModelListFromYAMLFile::loadFromFile (const char* file_path,
                                          PhyloTree* report_to_tree) {
    try {
        YAML::Node yaml_model_list = YAML::LoadFile(file_path);
        ModelFileLoader loader(file_path);
        if (!yaml_model_list.IsSequence()) {
            throw YAML::Exception(yaml_model_list.Mark(),
                                  "list '[...]' expected");
        }
        for (auto node : yaml_model_list) {
            if (!(node["substitutionmodel"])) {
                continue;
            }
            std::string yaml_model_name = node["substitutionmodel"].Scalar();
            TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                          "Parsing YAML model " << yaml_model_name);
            ModelInfoFromYAMLFile &y = models_found[yaml_model_name]
                                     = ModelInfoFromYAMLFile();
            loader.parseYAMLSubstitutionModel(node, yaml_model_name, y, *this,
                                              nullptr, report_to_tree);
        }
    }
    catch (YAML::ParserException e) {
        outError(e.what());
    }
    catch (YAML::Exception e) {
        outError(e.what());
    }
    catch (ModelExpression::ModelException x) {
        outError(x.getMessage());
    }
}

bool ModelListFromYAMLFile::isModelNameRecognized (const char* model_name) {
    size_t i = 0;
    while (model_name[i]!='\0' && model_name[i]!='{') {
        ++i;
    }
    auto found      = models_found.find(std::string(model_name, i));
    bool recognized = found != models_found.end();
    return recognized;
}

template <class S> class YAMLModelWrapper: public S {
protected:
    ModelInfoFromYAMLFile model_info;
    PhyloTree*            report_tree;
public:
    typedef S super;
    using   S::freq_type;
    using   S::num_params;
    using   S::num_states;
    using   S::rates;
    using   S::state_freq;
    
    using   S::getNDim;
    using   S::getNumberOfRates;
    using   S::setRateMatrix;
    
    YAMLModelWrapper(const ModelInfoFromYAMLFile& info,
                     PhyloTree* report_to_tree)
        : super(report_to_tree, report_to_tree)
        , model_info(info), report_tree(report_to_tree) {
    }
    
    void acceptParameterList(std::string parameter_list) {
        //parameter_list is passed by value so it can be modified
        trimString(parameter_list);
        if (startsWith(parameter_list, "{") &&
            endsWith(parameter_list, "}")) {
            auto len = parameter_list.length();
            parameter_list = parameter_list.substr(1, len-2);
        }
        size_t param_list_length = parameter_list.length();
        size_t i                 = 0;
        int    bracket_depth     = 0;
        typedef ModelExpression::InterpretedExpression Expr;
        std::vector<Expr> expr_list;
        while (i<param_list_length) {
            size_t j = i;
            for (;j<param_list_length &&
                 (parameter_list[j]!=',' || 0<bracket_depth); ++j) {
                char ch = parameter_list[j];
                if (ch=='(') {
                    ++bracket_depth;
                }
                else if (ch==')') {
                    --bracket_depth;
                }
            }
            std::string param = parameter_list.substr(i, j-i);
            expr_list.emplace_back(model_info, param);
            i = j + 1;
        }
        bool fix = !report_tree->params->optimize_from_given_params;
        size_t position = 0;
        for (Expr& ix : expr_list) {
            ModelExpression::Expression* x = ix.expression();
            if (x->isAssignment()) {
                typedef ModelExpression::Assignment A;
                typedef ModelExpression::Variable V;
                A*             a        = dynamic_cast<A*>(x);
                V*             xv       = a->getTargetVariable();
                string         var_name = xv->getName();
                double         setting  = a->getExpression()->evaluate();
                ModelVariable& mv       = model_info.assign(var_name, setting);
                if (fix) {
                    mv.markAsFixed();
                }
                TREE_LOG_LINE(*report_tree, YAMLModelVerbosity,
                              "Set " << xv->getName() << " to " << setting 
                              << " by name." );
            } else {
                double         setting  = x->evaluate();
                ModelVariable& mv       = model_info.assignByPosition(position, setting);
                string         var_name = model_info.getVariableNamesByPosition()[position];
                //No need for an index check for var_name's look-up, because
                //assignByPosition would already have bombed if it were out of bounds.
                if (fix) {
                    mv.markAsFixed();
                }
                TREE_LOG_LINE(*report_tree, YAMLModelVerbosity,
                              "Set " << var_name << " to " << setting 
                              << " by position." );
                ++position;
            }
        }
        if (!expr_list.empty()) {
            setRateMatrixFromModel();
        }
    }
    
    virtual void setBounds(double *lower_bound, double *upper_bound,
                            bool *bound_check) {
        //
        int ndim = getNDim();
        for (int i = 1; i <= ndim; ++i) {
            lower_bound[i] = MIN_RATE;
            upper_bound[i] = MAX_RATE;
            bound_check[i] = false;
        }
        model_info.setBounds(ndim, lower_bound,
                             upper_bound, bound_check);
    }
    
    virtual bool getVariables(double *variables) {
        bool changed = false;
        if (num_params > 0) {
            int num_all = getNumberOfRates();
            for (int i = 0; i < num_all; i++) {
                if (rates[i] != variables[i] ) {
                    TREE_LOG_LINE(*report_tree, VerboseMode::VB_MAX,
                                  "  estimated rates[" << i << "] changing"
                                  " from " << rates[i] << " to " << variables[i]);
                    rates[i] = variables[i];
                    changed  = true;
                }
            }
        }
        int ndim = getNDim();
        int first_freq_index = (ndim-num_states+2);
        if (freq_type == StateFreqType::FREQ_ESTIMATE) {
            auto read_freq = variables+first_freq_index;
            for (int i=0; i<num_states-1; ++i) {
                if (state_freq[i]!=read_freq[i]) {
                    TREE_LOG_LINE(*report_tree, VerboseMode::VB_MAX,
                                  "  estimated freqs[" << i << "] changing"
                                  << " from " << state_freq[i]
                                  << " to " << read_freq[i]);
                    state_freq[i] = read_freq[i];
                    changed       = true;
                }
            }
            //Set the last frequency to the residual
            //(one minus the sum of the others)
            if (scaleStateFreq()) {
                changed = true;
                model_info.assignLastFrequency(state_freq[num_states-1]);
            }
        } else {
            changed |= freqsFromParams(state_freq, variables+num_params+1,
                                       freq_type);
        }
        TREE_LOG_LINE(*report_tree, VerboseMode::VB_MAX, "");
        if (changed) {
            model_info.updateVariables(variables, first_freq_index, ndim);
            model_info.logVariablesTo(*report_tree);
            setRateMatrixFromModel();
        }
        return changed;
    }
    
    virtual bool scaleStateFreq() {
        // make the frequencies sum to 1.0
        bool   changed = false;
        double sum     = 0.0;
        for (int i = 0; i < num_states-1; ++i) {
            sum += state_freq[i];
        }
        if (1.0<sum) {
            sum     += state_freq[num_states-1];
            changed  = true;
            for (int i = 0; i < num_states; ++i) {
                state_freq[i] /= sum;
            }
        } else {
            //Set last state frequency to 1.0 minus
            //the sum of the others
            double residual = 1.0 - sum;
            if (state_freq[num_states-1] != residual) {
                state_freq[num_states-1] = residual;
                changed = true;
            }
        }
        return changed;
    }

    virtual void setVariables(double *variables) {
        if (num_params > 0) {
            for (int i = 0; i < num_params; ++i) {
                variables[i] = rates[i];
            }
        }
        if (freq_type == StateFreqType::FREQ_ESTIMATE) {
            // 2015-09-07: relax the sum of state_freq to be 1.0,
            // this will be done at the end of optimization
            int ndim = getNDim();
            memcpy(variables+(ndim-num_states+2), state_freq,
                   (num_states-1)*sizeof(double));
        } else {
            paramsFromFreqs(variables+num_params+1,
                            state_freq, freq_type);
        }
    }
    
    void setRateMatrixFromModel() {
        auto rank = model_info.getRateMatrixRank();
        ASSERT( rank == num_states );
        
        DoubleVector      rates;
        const char*       separator = "";
        std::stringstream trace;
        trace << "Rate Matrix: { ";
        
        model_info.forceAssign("num_states", num_states);
        ModelVariable& row_var    = model_info.forceAssign("row",    0);
        ModelVariable& column_var = model_info.forceAssign("column", 0);
        
        for (int row = 0; row < rank; ++row) {
            row_var.setValue((double)row+1);
            for (int col = 0; col < rank; ++col) {
                column_var.setValue((double)col+1);
                if (col != row) {
                    std::string expr_string =
                        model_info.getRateMatrixExpression(row,col);
                    typedef ModelExpression::InterpretedExpression Interpreter;
                    try {
                        Interpreter interpreter(model_info, expr_string);
                        double entry = interpreter.evaluate();
                        rates.push_back(entry);
                        trace << separator << entry;
                    }
                    catch (ModelExpression::ModelException x) {
                        std::stringstream msg;
                        msg << "Error parsing expression"
                            << " for " << model_info.getName()
                            << " rate matrix entry"
                            << " for row "    << (row + 1) << ","
                            << " and column " << (col + 1) << ": "
                            << x.getMessage();
                        outError(msg.str());
                    }
                } else {
                    trace << separator << "-";
                }
                separator = ", ";
            }
        }
        trace << " }";
        TREE_LOG_LINE(*report_tree, VerboseMode::VB_MAX, trace.str());
        setRateMatrix(rates.data());
    }
    
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk) {
        int state_num = static_cast<int>(state);
        if ( state_num < model_info.getTipLikelihoodMatrixRank()) {
            model_info.computeTipLikelihoodsForState(state, num_states, state_lk);
        } else if (state_num < num_states) {
            // single state
            memset(state_lk, 0, num_states*sizeof(double));
            state_lk[state] = 1.0;
        } else {
            // unknown state
            for (int i = 0; i < num_states; i++) {
                state_lk[i] = 1.0;
            }
        }
    }
    
    virtual void writeInfo(ostream &out) {
        if (YAMLVariableVerbosity <= verbose_mode ) {
            auto rates = model_info.getParameterList(ModelParameterType::RATE);
            if (!rates.empty()) {
                out << "Rate parameters: " << rates << std::endl;
            }
        }
        if (YAMLFrequencyVerbosity <= verbose_mode) {
            auto freqs = model_info.getParameterList(ModelParameterType::FREQUENCY);
            if (!freqs.empty()) {
                out << "State frequencies: " << freqs << std::endl;
            }
        }
    }
};
    
class YAMLModelDNA: public YAMLModelWrapper<ModelDNA> {
public:
    typedef YAMLModelWrapper<ModelDNA> super;
    YAMLModelDNA(const char *model_name, string model_params,
                 StateFreqType freq, string freq_params,
                 PhyloTree *tree, PhyloTree* report_to_tree,
                 const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
        init(model_name, model_params, freq,
             freq_params, report_to_tree);
        setRateMatrixFromModel();
    }
};

class YAMLModelDNAError: public YAMLModelWrapper<ModelDNAError> {
public:
    typedef YAMLModelWrapper<ModelDNAError> super;
    YAMLModelDNAError(const char *model_name, string model_params,
                      StateFreqType freq, string freq_params,
                      PhyloTree *tree, PhyloTree* report_to_tree,
                      const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
        init(model_name, model_params, freq,
             freq_params, report_to_tree);        
        setRateMatrixFromModel();
    }
    bool getVariables(double *variables) {
        bool changed = super::getVariables(variables);
        if (changed && !fix_epsilon) {
            epsilon = model_info.getVariableValue("epsilon");
        }
        return changed;
    }
};

class YAMLModelProtein: public YAMLModelWrapper<ModelProtein> {
public:
    typedef YAMLModelWrapper<ModelProtein> super;
    YAMLModelProtein(ModelsBlock* block,
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
};

class YAMLModelBinary: public YAMLModelWrapper<ModelBIN> {
public:
    typedef YAMLModelWrapper<ModelBIN> super;
    YAMLModelBinary(const char *model_name, std::string model_params,
                 StateFreqType freq, std::string freq_params,
                 PhyloTree *tree, PhyloTree* report_to_tree,
                 const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
        init(model_name, model_params, freq, freq_params, report_to_tree);
        setRateMatrixFromModel();
    }
};

class YAMLModelMorphology: public YAMLModelWrapper<ModelMorphology> {
public:
    typedef YAMLModelWrapper<ModelMorphology> super;
    YAMLModelMorphology(const char *model_name, std::string model_params,
                 StateFreqType freq, std::string freq_params,
                 PhyloTree *tree, PhyloTree* report_to_tree,
                 const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
        init(model_name, model_params, freq, freq_params, report_to_tree);
        setRateMatrixFromModel();
    }
};

class YAMLModelCodon: public YAMLModelWrapper<ModelCodon> {
public:
    typedef YAMLModelWrapper<ModelCodon> super;
    YAMLModelCodon(const char *model_name, std::string model_params,
                 StateFreqType freq, std::string freq_params,
                 PhyloTree *tree, PhyloTree* report_to_tree,
                 const ModelInfoFromYAMLFile& info): super(info, report_to_tree) {
        setReversible(info.isReversible());
        init(model_name, model_params, freq, freq_params, report_to_tree);
        setRateMatrixFromModel();
    }
};

bool ModelListFromYAMLFile::hasModel(const std::string& model_name) const {
    return models_found.find(model_name) != models_found.end();
}

StrVector ModelListFromYAMLFile::getModelNames() const {
    StrVector answer;
    for (auto it = models_found.begin(); it != models_found.end(); ++it) {
        answer.push_back(it->first);
    }
    return answer;
}

const ModelInfoFromYAMLFile&
    ModelListFromYAMLFile::getModel(const std::string& model_name) const {
    auto it = models_found.find(model_name);
    ASSERT(it != models_found.end());
    return it->second;
}

namespace {
    void extractModelNameAndParameters(const char* model_plus_params,
                                       std::string& name,
                                       std::string& params) {
        size_t i = 0;
        while (model_plus_params[i]!='\0' && model_plus_params[i]!='{') {
            ++i;
        }
        name   = std::string(model_plus_params,i);
        params = std::string(model_plus_params+i);
        //includes the { at front and the } at back.
    }
};

ModelMarkov* ModelListFromYAMLFile::getModelByName(const char* model_name,   PhyloTree *tree,
                                                   const char* model_params, StateFreqType freq_type,
                                                   const char* freq_params,  ModelsBlock* models_block,
                                                   PhyloTree* report_to_tree) {
    std::string name;
    std::string parameter_list;
    extractModelNameAndParameters(model_name, name, parameter_list);
    ModelInfoFromYAMLFile& model_info = models_found[model_name];
    if (0<strlen(model_params) || 0<strlen(freq_params)) {
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "Model Params: " << model_params
                      << " Freq Params: " << freq_params);
        if (parameter_list.empty()) {
            parameter_list = model_params;
        }
    }
    if (freq_type == StateFreqType::FREQ_UNKNOWN) {
        freq_type = model_info.frequency_type;
    }
    switch (model_info.sequence_type) {
        case SeqType::SEQ_BINARY:
            return getBinaryModel(model_info, parameter_list,
                                 freq_type, tree, report_to_tree);
        case SeqType::SEQ_CODON:
            return getCodonModel(model_info, parameter_list,
                                 freq_type, tree, report_to_tree);

        case SeqType::SEQ_DNA:
            return getDNAModel(model_info, parameter_list,
                               freq_type, tree, report_to_tree);

        case SeqType::SEQ_MORPH:
            return getMorphologicalModel(model_info, parameter_list,
                                         freq_type, tree, report_to_tree);

        case SeqType::SEQ_PROTEIN:
            return getProteinModel(model_info, parameter_list,
                                   freq_type, tree, models_block,
                                   report_to_tree);

        default:
            outError("YAML model uses unsupported sequence type");
            return nullptr;
    };
}

ModelMarkov* ModelListFromYAMLFile::getBinaryModel(ModelInfoFromYAMLFile& model_info,
                                                   const std::string& parameter_list,
                                                   StateFreqType freq_type,
                                                   PhyloTree* tree,
                                                   PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_BINARY);
    YAMLModelBinary* model;
    model = new YAMLModelBinary("", dummy_rate_params, freq_type,
                                dummy_freq_params, tree, report_to_tree,
                                model_info);
    model->acceptParameterList(parameter_list);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getCodonModel(ModelInfoFromYAMLFile& model_info,
                                                  const std::string& parameter_list,
                                                  StateFreqType freq_type,
                                                  PhyloTree* tree,
                                                  PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_CODON);
    YAMLModelCodon* model;

    model = new YAMLModelCodon("", dummy_rate_params, freq_type,
                               dummy_freq_params, tree, report_to_tree,
                               model_info);
    model->acceptParameterList(parameter_list);
    return model;
}
    
ModelMarkov* ModelListFromYAMLFile::getDNAModel(ModelInfoFromYAMLFile& model_info,
                                                const std::string& parameter_list,
                                                StateFreqType freq_type,
                                                PhyloTree* tree,
                                                PhyloTree* report_to_tree) {
    ModelMarkov* model = nullptr;
    const YAMLFileParameter* p =
        model_info.findParameter("epsilon",
                                 ModelParameterType::RATE);
    if (p!=nullptr) {
        if (p->is_subscripted) {
            outError("epsilon parameter for DNA+error model"
                     "may not be subscripted");
        }
        double epsilon = p->value;
        bool epsilon_is_fixed = model_info.isVariableFixed("epsilon");
        model_info.moveParameterToBack("epsilon",
                                       ModelParameterType::RATE);
        YAMLModelDNAError* emodel;
        emodel = new YAMLModelDNAError("", dummy_rate_params, freq_type,
                                      dummy_freq_params, tree,
                                      report_to_tree, model_info);
        std::string error_model = model_info.getStringProperty("errormodel", "+E");
        emodel->setEpsilon(epsilon, epsilon_is_fixed, error_model);
        TREE_LOG_LINE(*report_to_tree, YAMLModelVerbosity,
                      "epsilon is " << epsilon
                      << ", fixed is " << epsilon_is_fixed
                      << ", and errormodel is " << error_model);
        emodel->acceptParameterList(parameter_list);
        model = emodel;
    }
    else {
        YAMLModelDNA* dmodel;
        dmodel = new YAMLModelDNA("", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree,
                                  report_to_tree, model_info);
        dmodel->acceptParameterList(parameter_list);
        model = dmodel;
    }
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getMorphologicalModel(ModelInfoFromYAMLFile& model_info,
                                                          const std::string& parameter_list,
                                                          StateFreqType freq_type,
                                                          PhyloTree* tree,
                                                          PhyloTree* report_to_tree) {
    YAMLModelMorphology* model;
    model = new YAMLModelMorphology("", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree, report_to_tree,
                                  model_info);
    model->acceptParameterList(parameter_list);
    return model;
}

ModelMarkov* ModelListFromYAMLFile::getProteinModel(ModelInfoFromYAMLFile& model_info,
                                                    const std::string& parameter_list,
                                                    StateFreqType freq_type,
                                                    PhyloTree* tree,
                                                    ModelsBlock* models_block,
                                                    PhyloTree* report_to_tree) {
    insistOnAlignmentSequenceType(tree->aln, SeqType::SEQ_PROTEIN);
    YAMLModelProtein* model;
    model = new YAMLModelProtein(models_block, "", dummy_rate_params, freq_type,
                                  dummy_freq_params, tree, report_to_tree,
                                  model_info);
    model->acceptParameterList(parameter_list);
    return model;
}

void ModelListFromYAMLFile::insistOnAlignmentSequenceType(const Alignment* alignment, 
                                                          SeqType desired_type) const {
    if (alignment->seq_type != desired_type) {
        std::stringstream complaint;
        complaint << "Cannot use " << getSeqTypeName(desired_type) << " model"
                  << " with a " << getSeqTypeName(alignment->seq_type) << "alignment."
                  << "\nPlease re-run, with"
                  << " a " << getSeqTypeName(alignment->seq_type) << " model,"
                  << " or passing -seqtype " << getSeqTypeShortName(desired_type, false)
                  << " on the command-line.";
        outError(complaint.str());
    }
}