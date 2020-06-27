//
//  modeldnaerror.cpp
//  model
//
//  Created by Minh Bui on 26/6/20.
//

#include "modeldnaerror.h"

#define MIN_EPSILON 0.0001
#define MAX_EPSILON 0.5

ModelDNAError::ModelDNAError(PhyloTree *tree)
: ModelDNA(tree)
{
    epsilon = 0.05;
    fix_epsilon = false;
}

ModelDNAError::ModelDNAError(const char *model_name, string model_params,
                             StateFreqType freq, string freq_params, string seqerr, PhyloTree *tree)
: ModelDNA(model_name, model_params, freq, freq_params, tree)
{
    epsilon = 0.05;
    fix_epsilon = false;
    seqerr_name = seqerr;
}

void ModelDNAError::startCheckpoint() {
    checkpoint->startStruct("ModelDNAError");
}

void ModelDNAError::saveCheckpoint() {
    startCheckpoint();
    if (!fix_epsilon)
        CKP_SAVE(epsilon);
    endCheckpoint();
    ModelDNA::saveCheckpoint();
}

void ModelDNAError::restoreCheckpoint() {
    ModelDNA::restoreCheckpoint();
    startCheckpoint();
    if (!fix_epsilon)
        CKP_RESTORE(epsilon);
    endCheckpoint();
}

string ModelDNAError::getName() {
    string retname = ModelDNA::getName();
    retname += seqerr_name;
    return retname;
}

string ModelDNAError::getNameParams() {
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

int ModelDNAError::getNDim() {
    if (fix_epsilon)
        return ModelDNA::getNDim();
    else
        return ModelDNA::getNDim() + 1;
}

void ModelDNAError::computeTipLikelihood(PML::StateType state, double *state_lk) {
    if (epsilon == 0.0)
        return ModelDNA::ModelSubst::computeTipLikelihood(state, state_lk);

    int i;
    if (state < num_states) {
        // single state
        for (i = 0; i < num_states; i++)
            state_lk[i] = epsilon/3.0;
        state_lk[state] = 1.0 - epsilon;
    } else if (state < 18) {
        // ambiguous (polymorphic) state
        int cstate = state-num_states+1;
        int obs_states = 0; // number of observed states
        for (i = 0; i < num_states; i++)
            if ((cstate) & (1 << i))
                obs_states++;
        
        double obs_lk = 1.0 - (4-obs_states)*epsilon/3.0;
        double unobs_lk = obs_states*epsilon/3.0;
        for (i = 0; i < num_states; i++) {
            if ((cstate) & (1 << i))
                state_lk[i] = obs_lk;
            else
                state_lk[i] = unobs_lk;
        }
    } else {
        // unknown state
        for (int i = 0; i < num_states; i++)
            state_lk[i] = 1.0;
    }
}

void ModelDNAError::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelDNA::setBounds(lower_bound, upper_bound, bound_check);
    if (!fix_epsilon) {
        int id = ModelDNA::getNDim()+1;
        lower_bound[id] = MIN_EPSILON;
        upper_bound[id] = MAX_EPSILON;
        bound_check[id] = false;
    }
}

bool ModelDNAError::getVariables(double *variables) {
    bool changed = ModelDNA::getVariables(variables);
    if (!fix_epsilon) {
        int id = ModelDNA::getNDim()+1;
        changed |= (epsilon != variables[id]);
        epsilon = variables[getNDim()];
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
