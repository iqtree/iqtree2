//
//  modeldnaerror.cpp
//  model
//
//  Created by Minh Bui on 26/6/20.
//

#include "modeldnaerror.h"
#include <utils/stringfunctions.h> //for convert_double


// Bound for sequencing error probability (epsilon)
#define MIN_EPSILON 0.0001
#define MAX_EPSILON 0.5

ModelDNAError::ModelDNAError(PhyloTree* tree, PhyloTree* report_to_tree)
: ModelDNA(tree, report_to_tree), seqerr_name("+E")
{
    epsilon     = 0.05;
    fix_epsilon = false;
}

ModelDNAError::ModelDNAError(const char* model_name, const string& model_params,
                             StateFreqType freq, const string& freq_params, 
                             const string& seqerr, PhyloTree* tree, 
                             PhyloTree* report_to_tree)
    : ModelDNA(model_name, model_params, freq, freq_params, tree, report_to_tree)
    , seqerr_name(seqerr)
{
    epsilon = 0.05;
    fix_epsilon = false;
    // now parse the epsilon parameter
    string::size_type pos = seqerr.find(OPEN_BRACKET);
    if (pos != string::npos) {
        auto end_pos = seqerr.find(CLOSE_BRACKET);
        if (end_pos == string::npos) {
            outError("Missing closing bracket in " + seqerr);
        }
        epsilon = convert_double(seqerr.substr(pos+1, end_pos-pos-1).c_str());
        if (epsilon < 0.0 || epsilon > 1.0) {
            outError("Sequencing error probability " +
                     convertDoubleToString(epsilon) +
                     " is not between 0 and 1");
        }
        if (!Params::getInstance().optimize_from_given_params) {
            fix_epsilon = true;
        }
        seqerr_name = seqerr.substr(0, pos);
    }
}

void ModelDNAError::setEpsilon(double e, bool fixed, const std::string& seqerr) {
    epsilon = e;
    fix_epsilon = fixed;
    seqerr_name = seqerr;
}


void ModelDNAError::startCheckpoint() {
    checkpoint->startStruct("ModelDNAError");
}

void ModelDNAError::saveCheckpoint() {
    startCheckpoint();
    if (!fix_epsilon) {
        CKP_SAVE(epsilon);
    }
    endCheckpoint();
    ModelDNA::saveCheckpoint();
}

void ModelDNAError::restoreCheckpoint() {
    ModelDNA::restoreCheckpoint();
    startCheckpoint();
    if (!fix_epsilon) {
        CKP_RESTORE(epsilon);
    }
    endCheckpoint();
}

string ModelDNAError::getName() const {
    string retname = ModelDNA::getName();
    retname += seqerr_name;
    return retname;
}

std::string ModelDNAError::getNameParams() const {
    string retname = ModelDNA::getNameParams();
    retname += seqerr_name + "{" + convertDoubleToString(epsilon) + "}";
    return retname;
}

void ModelDNAError::writeInfo(ostream &out) {
    ModelDNA::writeInfo(out);
    auto prec = out.precision(6);
    out << "Sequencing error probability: " << epsilon << endl;
    out.precision(prec);
}

int ModelDNAError::getNDim() const {
    if (fix_epsilon) {
        return ModelDNA::getNDim();
    }
    else {
        return ModelDNA::getNDim() + 1;
    }
}

int  ModelDNAError::getErrorNucleotideState() const {
    if (seqerr_name == "+EA") {
        return 0;
    } else if (seqerr_name == "+EC") {
        return 1;
    } else if (seqerr_name == "+EG") {
        return 2;
    } else if (seqerr_name == "+ET") {
        return 3;
    } else if (seqerr_name == "+E")  {
        return -1;
    }
    outError("Unknown sequencing error model " + seqerr_name);
    return -1;
}

void ModelDNAError::computeTipLikelihood(PML::StateType state, 
                                         double *state_lk) const {
    if (epsilon == 0.0) {
        ModelDNA::ModelSubst::computeTipLikelihood(state, state_lk);
        return;
    }
    int b = getErrorNucleotideState();
    
    // true for observed states, false for unobserved state
    bool observed[4] = {false, false, false, false};
    int num_observed = 0; // number of observed states
    if (state < 4) {
        // single state
        observed[state] = true;
        num_observed = 1;
    } else if (state < 18) {
        // ambiguous (polymorphic) state
        int cstate = state-num_states+1;
        for (int i = 0; i < num_states; ++i) {
            if ((cstate) & (1 << i)) {
                observed[i] = true;
                num_observed++;
            }
        }
    } else {
        // unknown state
        for (int i = 0; i < num_states; i++) {
            observed[i] = true;
        }
        num_observed = num_states;
    }

    double observed_lk; // likelihood of observed state
    double unobserved_lk; // likelihood of unobserved state
    if (b >= 0) {
        // nucleotide-specific error towards nucleotide b (Nicola de Maio)
        observed_lk = observed[b] ? 1.0 : 1.0-epsilon;
        unobserved_lk = observed[b] ? epsilon : 0.0;
    } else {
        // uniform error model (Felsenstein 2004)
        observed_lk = 1.0 - (4-num_observed)*epsilon/3.0;
        unobserved_lk = num_observed*epsilon/3.0;
    }
    for (int i = 0; i < num_states; i++) {
        state_lk[i] = observed[i] ? observed_lk : unobserved_lk;
    }
}

void ModelDNAError::setBounds(double *lower_bound, double *upper_bound,
                              bool *bound_check) {
    ModelDNA::setBounds(lower_bound, upper_bound, bound_check);
    if (!fix_epsilon) {
        int id = ModelDNA::getNDim()+1;
        lower_bound[id] = MIN_EPSILON;
        upper_bound[id] = MAX_EPSILON;
        bound_check[id] = false;
    }
}

bool ModelDNAError::getVariables(const double *variables) {
    bool changed = ModelDNA::getVariables(variables);
    if (!fix_epsilon) {
        int id = ModelDNA::getNDim()+1;
        changed |= (epsilon != variables[id]);
        epsilon = variables[id];
    }
    return changed;
}

void ModelDNAError::setVariables(double *variables) {
    ModelDNA::setVariables(variables);
    if (!fix_epsilon) {
        int id = ModelDNA::getNDim()+1;
        variables[id] = epsilon;
    }
}
