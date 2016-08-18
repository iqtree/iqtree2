//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"

ModelPoMoMixture::ModelPoMoMixture(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params, string pomo_rate_str)
	:  
        ModelGTR(tree),
        ModelPoMo(model_name, model_params, freq_type, freq_params, tree, is_reversible, pomo_params),
        ModelMixture(tree)
    
{
    // initialize pomo_mixture
}


ModelPoMoMixture::~ModelPoMoMixture() {
}

void ModelPoMoMixture::saveCheckpoint() {
    ModelPoMo::saveCheckpoint();
}

void ModelPoMoMixture::restoreCheckpoint() {
    ModelPoMo::restoreCheckpoint();
}


int ModelPoMoMixture::getNDim() {
    return ModelPoMo::getNDim();
}


int ModelPoMoMixture::getNDimFreq() {
    return ModelPoMo::getNDimFreq();
}


double ModelPoMoMixture::targetFunk(double x[]) {
    return ModelPoMo::targetFunk(x);
}



void ModelPoMoMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelPoMo::setBounds(lower_bound, upper_bound, bound_check);
}


void ModelPoMoMixture::writeInfo(ostream &out) {
    ModelPoMo::writeInfo(out);
}


void ModelPoMoMixture::setVariables(double *variables) {
    ModelPoMo::setVariables(variables);
}


bool ModelPoMoMixture::getVariables(double *variables) {
    return ModelPoMo::getVariables(variables);
}


