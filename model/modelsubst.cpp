//
// C++ Implementation: substmodel
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "modelsubst.h"
#include "utils/tools.h"

ModelSubst::ModelSubst(int nstates) : Optimization(), CheckpointFactory()
{
	num_states = nstates;
	name = "JC";
	full_name = "JC (Juke and Cantor, 1969)";
	state_freq = new double[num_states];
	for (int i = 0; i < num_states; i++)
		state_freq[i] = 1.0 / num_states;
	freq_type = FREQ_EQUAL;
    fixed_parameters = false;
//    linked_model = NULL;
}

void ModelSubst::startCheckpoint() {
    checkpoint->startStruct("ModelSubst");
}

void ModelSubst::saveCheckpoint() {
    startCheckpoint();
//    CKP_SAVE(num_states);
//    CKP_SAVE(name);
//    CKP_SAVE(full_name);
//    CKP_SAVE(freq_type);
    if (freq_type == FREQ_ESTIMATE && !fixed_parameters)
        CKP_ARRAY_SAVE(num_states, state_freq);
    endCheckpoint();
    CheckpointFactory::saveCheckpoint();
}

void ModelSubst::restoreCheckpoint() {
    CheckpointFactory::restoreCheckpoint();
    startCheckpoint();
//    CKP_RESTORE(num_states);
//    CKP_RESTORE(name);
//    CKP_RESTORE(full_name);
//    int freq_type = this->freq_type;
//    CKP_RESTORE(freq_type);
//    this->freq_type = (StateFreqType)freq_type;
    if (freq_type == FREQ_ESTIMATE && !fixed_parameters)
        CKP_ARRAY_RESTORE(num_states, state_freq);
    endCheckpoint();

    decomposeRateMatrix();
}

// here the simplest Juke-Cantor model is implemented, valid for all kind of data (DNA, AA,...)
void ModelSubst::computeTransMatrix(double time, double *trans_matrix, int mixture, int selected_row) {
	double non_diagonal = (1.0 - exp(-time*num_states/(num_states - 1))) / num_states;
	double diagonal = 1.0 - non_diagonal * (num_states - 1);
	int nstates_sqr = num_states * num_states;

	for (int i = 0; i < nstates_sqr; i++)
		if (i % (num_states+1) == 0) 
			trans_matrix[i] = diagonal; 
		else 
			trans_matrix[i] = non_diagonal;
}


double ModelSubst::computeTrans(double time, int state1, int state2) {
	double expt = exp(-time * num_states / (num_states-1));
	if (state1 != state2) {
		return (1.0 - expt) / num_states;
	}
	return (1.0 + (num_states-1)*expt) / num_states;

/*	double non_diagonal = (1.0 - exp(-time*num_states/(num_states - 1))) / num_states;
	if (state1 != state2)
		return non_diagonal;
	return 1.0 - non_diagonal * (num_states - 1);*/
}

double ModelSubst::computeTrans(double time, int model_id, int state1, int state2) {
	return computeTrans(time, state1, state2);
}

double ModelSubst::computeTrans(double time, int state1, int state2, double &derv1, double &derv2) {
	double coef = -double(num_states) / (num_states-1);
	double expt = exp(time * coef);
	if (state1 != state2) {
		derv1 = expt / (num_states-1);
		derv2 = derv1 * coef;
		return (1.0 - expt) / num_states;
	}

	derv1 = -expt;
	derv2 = derv1 * coef;
	return (1.0 + (num_states-1)*expt) / num_states;
}

double ModelSubst::computeTrans(double time, int model_id, int state1, int state2, double &derv1, double &derv2) {
	return computeTrans(time, state1, state2, derv1, derv2);
}

void ModelSubst::getRateMatrix(double *rate_mat) {
	int nrate = getNumRateEntries();
	for (int i = 0; i < nrate; i++)
		rate_mat[i] = 1.0;
}

void ModelSubst::getQMatrix(double *q_mat, int mixture) {
	int i, j, k;
	for (i = 0, k = 0; i < num_states; i++)
		for (j = 0; j < num_states; j++, k++)
			if (i == j) q_mat[k] = -1.0; else q_mat[k] = 1.0/3;
}

void ModelSubst::getStateFrequency(double *state_freq, int mixture) {
	double freq = 1.0 / num_states;
	for (int i = 0; i < num_states; i++)
		state_freq[i] = freq;
}

void ModelSubst::setStateFrequency(double *state_freq) {
    memcpy(this->state_freq, state_freq, sizeof(double)*num_states);
}

void ModelSubst::computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2, int mixture)
{
	double expf = exp(-time*num_states/(num_states - 1));
	double non_diag = (1.0 - expf) / num_states;
	double diag = 1.0 - non_diag * (num_states - 1);
	double derv1_non_diag = expf / (num_states-1);
	double derv1_diag = -expf;
	double derv2_non_diag = -derv1_non_diag*num_states/(num_states-1);
	double derv2_diag = -derv1_diag*num_states/(num_states-1);

	int nstates_sqr = num_states * num_states;
	int i;
	for (i = 0; i < nstates_sqr; i++)
		if (i % (num_states+1) == 0) { 
			trans_matrix[i] = diag;
			trans_derv1[i] = derv1_diag;
			trans_derv2[i] = derv2_diag;
		} else { 
			trans_matrix[i] = non_diag;
			trans_derv1[i] = derv1_non_diag;
			trans_derv2[i] = derv2_non_diag;
		}

	// DEBUG
	/*int j;
	if (verbose_mode == VB_DEBUG) {
		cout.precision(4);
		cout << "time = " << time << endl;
		for (i = 0; i < num_states; i++, cout << endl) {
			for (j = 0; j < num_states; j++) {
				cout.width(8);
				cout << right << trans_matrix[i*num_states+j] << " ";
			}
			cout << "| ";
			for (j = 0; j < num_states; j++) {
				cout << right << trans_derv1[i*num_states+j] << " ";
				cout.width(8);
			}
			cout << "| ";
			for (j = 0; j < num_states; j++) {
				cout.width(8);
				cout << right << trans_derv2[i*num_states+j] << " ";
			}
		}
		cout.precision(10);
	}*/

}

void ModelSubst::multiplyWithInvEigenvector(double *state_lk) {
    int nmixtures = getNMixtures();
    double *inv_eigenvectors = getInverseEigenvectors();
    double saved_state_lk[num_states];
    memcpy(saved_state_lk, state_lk, sizeof(double)*num_states);
    memset(state_lk, 0, sizeof(double)*num_states*nmixtures);
    for (int m = 0; m < nmixtures; m++) {
        double *inv_evec = &inv_eigenvectors[m*num_states*num_states];
        double *this_state_lk = &state_lk[m*num_states];
        for (int i = 0; i < num_states; i++)
            for (int j = 0; j < num_states; j++, inv_evec++)
              this_state_lk[i] += (*inv_evec) * saved_state_lk[j];
    }
}

void ModelSubst::computeTipLikelihood(PML::StateType state, double *state_lk) {
    if (state < num_states) {
        // single state
        memset(state_lk, 0, num_states*sizeof(double));
        state_lk[state] = 1.0;
    } else {
        // unknown state
        for (int i = 0; i < num_states; i++)
            state_lk[i] = 1.0;
    }
}

double *ModelSubst::newTransMatrix() {
	return new double[num_states * num_states];
}

ModelSubst::~ModelSubst()
{
    // mem space pointing to target model and thus avoid double free here
//    if (linked_model && linked_model != this)
//        return;

    if (state_freq) delete [] state_freq;
}



