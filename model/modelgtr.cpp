/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
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
#include "modelgtr.h"
#include <stdlib.h>
#include <string.h>

//const double MIN_FREQ_RATIO = MIN_FREQUENCY;
//const double MAX_FREQ_RATIO = 1.0/MIN_FREQUENCY;

ModelGTR::ModelGTR(PhyloTree *tree, bool count_rates)
 : ModelSubst(tree->aln->num_states), EigenDecomposition()
{
    half_matrix = true;
	int i;
	int nrate = getNumRateEntries();
//	int ncoeff = num_states*num_states*num_states;
	
	highest_freq_state = num_states-1;
	name = "GTR";
	full_name = "GTR (Tavare, 1986)";
	phylo_tree = tree;
	
	rates = new double[nrate];
	memset(rates, 0, sizeof(double) * nrate);

	freq_type = FREQ_UNKNOWN;
	
	eigenvalues = aligned_alloc<double>(num_states);

	eigenvectors = aligned_alloc<double>(num_states*num_states);
//	for (i = 0; i < num_states; i++)
//		eigenvectors[i] = new double[num_states];

	inv_eigenvectors = aligned_alloc<double>(num_states*num_states);
//	for (i = 0; i < num_states; i++)
//		inv_eigenvectors[i] = new double[num_states];
		
//	eigen_coeff = aligned_alloc<double>(ncoeff);

//	if (count_rates) 
//		computeEmpiricalRate();
//	else
		for (i=0; i < nrate; i++) rates[i] = 1.0;
	//eigen_coeff_derv1 = new double[ncoeff];
	//eigen_coeff_derv2 = new double[ncoeff];
	num_params = getNumRateEntries() - 1;
}

void ModelGTR::saveCheckpoint() {
    checkpoint->startStruct("ModelGTR");
    checkpoint->endStruct();
    ModelSubst::saveCheckpoint();
}

void ModelGTR::restoreCheckpoint() {
    ModelSubst::restoreCheckpoint();
    checkpoint->startStruct("ModelGTR");
    checkpoint->endStruct();
}


void ModelGTR::setTree(PhyloTree *tree) {
	phylo_tree = tree;
}

string ModelGTR::getName() {
	if (getFreqType() == FREQ_EMPIRICAL)
		return name + "+F";
	else if (getFreqType() == FREQ_CODON_1x4)
		return name += "+F1X4";
	else if (getFreqType() == FREQ_CODON_3x4)
		return name + "+F3X4";
	else if (getFreqType() == FREQ_CODON_3x4C)
		return name + "+F3X4C";
	else if (getFreqType() == FREQ_ESTIMATE && phylo_tree->aln->seq_type != SEQ_DNA)
		return name + "+FO";
	else if (getFreqType() == FREQ_EQUAL && phylo_tree->aln->seq_type != SEQ_DNA)
		return name + "+FQ";
    else
        return name;
}

string ModelGTR::getNameParams() {

	ostringstream retname;
	retname << name;
//	if (num_states != 4) retname << num_states;
	retname << '{';
	int nrates = getNumRateEntries();
	for (int i = 0; i < nrates; i++) {
		if (i>0) retname << ',';
		retname << rates[i];
	}
	retname << '}';
    getNameParamsFreq(retname);
    return retname.str();    
}
    
void ModelGTR::getNameParamsFreq(ostream &retname) {
	if (getFreqType() == FREQ_EMPIRICAL || (getFreqType() == FREQ_USER_DEFINED && phylo_tree->aln->seq_type == SEQ_DNA)) {
		retname << "+F";
        retname << "{" << state_freq[0];
        for (int i = 1; i < num_states; i++)
            retname << "," << state_freq[i];
        retname << "}";
	} else if (getFreqType() == FREQ_CODON_1x4)
		retname << "+F1X4";
	else if (getFreqType() == FREQ_CODON_3x4)
		retname << "+F3X4";
	else if (getFreqType() == FREQ_CODON_3x4C)
		name += "+F3X4C";
	else if (getFreqType() == FREQ_ESTIMATE) {
		retname << "+FO";
        retname << "{" << state_freq[0];
        for (int i = 1; i < num_states; i++)
            retname << "," << state_freq[i];
        retname << "}";
    } else if (getFreqType() == FREQ_EQUAL && phylo_tree->aln->seq_type != SEQ_DNA)
		retname << "+FQ";
}

void ModelGTR::init(StateFreqType type) {
	//if (type == FREQ_UNKNOWN) return;
	int i;
	freq_type = type;
	assert(freq_type != FREQ_UNKNOWN);
	switch (freq_type) {
	case FREQ_EQUAL:
		if (phylo_tree->aln->seq_type == SEQ_CODON) {
			int nscodon = phylo_tree->aln->getNumNonstopCodons();
            double freq_codon = (1.0-(num_states-nscodon)*MIN_FREQUENCY)/(nscodon);
			for (i = 0; i < num_states; i++)
				if (phylo_tree->aln->isStopCodon(i))
					state_freq[i] = MIN_FREQUENCY;
				else
					state_freq[i] = freq_codon;
		} else {
            double freq_state = 1.0/num_states;
			for (i = 0; i < num_states; i++)
				state_freq[i] = freq_state;
		}
		break;	
	case FREQ_ESTIMATE:
	case FREQ_EMPIRICAL:
		if (phylo_tree->aln->seq_type == SEQ_CODON) {
			double ntfreq[12];
			phylo_tree->aln->computeCodonFreq(freq_type, state_freq, ntfreq);
//			phylo_tree->aln->computeCodonFreq(state_freq);
		} else
			phylo_tree->aln->computeStateFreq(state_freq);
		for (i = 0; i < num_states; i++)
			if (state_freq[i] > state_freq[highest_freq_state])
				highest_freq_state = i;
		break;
	case FREQ_USER_DEFINED:
		if (state_freq[0] == 0.0) outError("State frequencies not specified");
		break;
	default: break;
	}
	decomposeRateMatrix();
	if (verbose_mode >= VB_MAX)
		writeInfo(cout);

}

void ModelGTR::writeInfo(ostream &out) {
	if (num_states == 4) {
		out << "Rate parameters:";
		//out.precision(3);
		//out << fixed;
		out << "  A-C: " << rates[0];
		out << "  A-G: " << rates[1];
		out << "  A-T: " << rates[2];
		out << "  C-G: " << rates[3];
		out << "  C-T: " << rates[4];
		out << "  G-T: " << rates[5];
		out << endl;
		//if (freq_type != FREQ_ESTIMATE) return;
		out << "Base frequencies: ";
		out << "  A: " << state_freq[0];
		out << "  C: " << state_freq[1];
		out << "  G: " << state_freq[2];
		out << "  T: " << state_freq[3];
		out << endl;
	}
//	if (verbose_mode >= VB_DEBUG) {
//		int i, j;
//		out.precision(6);
//		out << "eigenvalues: " << endl;
//		for (i = 0; i < num_states; i++) out << " " << eigenvalues[i];
//		out << endl << "eigenvectors: " << endl;
//		for (i = 0; i < num_states; i++)  {
//			for (j = 0; j < num_states; j++)
//				out << " " << eigenvectors[i*num_states+j];
//			out << endl;
//		}
//		out << endl << "inv_eigenvectors: " << endl;
//		for (i = 0; i < num_states; i++)  {
//			for (j = 0; j < num_states; j++)
//				out << " " << inv_eigenvectors[i*num_states+j];
//			out << endl;
//		}
//	}
	//out.unsetf(ios::fixed);
}

void ModelGTR::computeTransMatrix(double time, double *trans_matrix) {
	/* compute P(t) */
	double evol_time = time / total_num_subst;
	double *exptime = new double[num_states];
	int i, j, k;

	for (i = 0; i < num_states; i++)
		exptime[i] = exp(evol_time * eigenvalues[i]);

	int row_offset;
	for (i = 0, row_offset = 0; i < num_states; i++, row_offset+=num_states) {
		double *trans_row = trans_matrix + row_offset;
		for (j = i+1; j < num_states; j ++) { 
			// compute upper triangle entries
			double *trans_entry = trans_row + j;
//			double *coeff_entry = eigen_coeff + ((row_offset+j)*num_states);
			*trans_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				*trans_entry += eigenvectors[i*num_states+k] * inv_eigenvectors[k*num_states+j] * exptime[k];
			}
			if (*trans_entry < 0.0) {
				*trans_entry = 0.0;
			}
			// update lower triangle entries
			trans_matrix[j*num_states+i] = (state_freq[i]/state_freq[j]) * (*trans_entry);
		}
		trans_row[i] = 0.0; // initialize diagonal entry
		// taking the sum of row
		double sum = 0.0;
		for (j = 0; j < num_states; j++)
			sum += trans_row[j];
		trans_row[i] = 1.0 - sum; // update diagonal entry
	}
	delete [] exptime;
}

void ModelGTR::computeTransMatrixFreq(double time, double* trans_matrix)
{
	computeTransMatrix(time, trans_matrix);
	for (int state1 = 0; state1 < num_states; state1++) {
		double *trans_mat_state = trans_matrix + (state1 * num_states);
		for (int state2 = 0; state2 < num_states; state2++)
			trans_mat_state[state2] *= state_freq[state1];
	}
}

double ModelGTR::computeTrans(double time, int state1, int state2) {
	double evol_time = time / total_num_subst;
	int i;

//	double *coeff_entry = eigen_coeff + ((state1*num_states+state2)*num_states);
	double trans_prob = 0.0;
	for (i = 0; i < num_states; i++) {
		trans_prob += eigenvectors[state1*num_states+i] * inv_eigenvectors[i*num_states+state2] * exp(evol_time * eigenvalues[i]);
	}
	return trans_prob;
}

double ModelGTR::computeTrans(double time, int state1, int state2, double &derv1, double &derv2) {
	double evol_time = time / total_num_subst;
	int i;

//	double *coeff_entry = eigen_coeff + ((state1*num_states+state2)*num_states);
	double trans_prob = 0.0;
	derv1 = derv2 = 0.0;
	for (i = 0; i < num_states; i++) {
		double trans = eigenvectors[state1*num_states+i] * inv_eigenvectors[i*num_states+state2] * exp(evol_time * eigenvalues[i]);
		double trans2 = trans * eigenvalues[i];
		trans_prob += trans;
		derv1 += trans2;
		derv2 += trans2 * eigenvalues[i];
	}
	return trans_prob;
}


void ModelGTR::computeTransDerv(double time, double *trans_matrix, 
	double *trans_derv1, double *trans_derv2) 
{
	/* compute P(t) */

	double evol_time = time / total_num_subst;
	double *exptime = new double[num_states];
	int i, j, k;

	for (i = 0; i < num_states; i++)
		exptime[i] = exp(evol_time * eigenvalues[i]);

	for (i = 0; i < num_states; i ++) {
		for (j = 0; j < num_states; j ++) {
			int offset = (i*num_states+j);
			double *trans_entry = trans_matrix + offset;
			double *derv1_entry = trans_derv1 + offset;
			double *derv2_entry = trans_derv2 + offset;

//			int coeff_offset = offset*num_states;
//			double *coeff_entry       = eigen_coeff + coeff_offset;
			*trans_entry = 0.0;
			*derv1_entry = 0.0;
			*derv2_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				double trans = eigenvectors[i*num_states+k] * inv_eigenvectors[k*num_states+j] * exptime[k];
				double trans2 = trans * eigenvalues[k];
				*trans_entry += trans;
				*derv1_entry += trans2;
				*derv2_entry += trans2 * eigenvalues[k];
			}
			if (*trans_entry < 0.0) {
				*trans_entry = 0.0;
			}
		}
	}
	delete [] exptime;
}

void ModelGTR::computeTransDervFreq(double time, double rate_val, double* trans_matrix, double* trans_derv1, double* trans_derv2)
{
	int nstates = num_states;
	double rate_sqr = rate_val*rate_val;
	computeTransDerv(time * rate_val, trans_matrix, trans_derv1, trans_derv2);
	for (int state1 = 0; state1 < nstates; state1++) {
		double *trans_mat_state = trans_matrix + (state1 * nstates);
		double *trans_derv1_state = trans_derv1 + (state1 * nstates);
		double *trans_derv2_state = trans_derv2 + (state1 * nstates);
		for (int state2 = 0; state2 < nstates; state2++) {
			trans_mat_state[state2] *= state_freq[state1];
			trans_derv1_state[state2] *= (state_freq[state1] * rate_val);
			trans_derv2_state[state2] *= (state_freq[state1] * rate_sqr);
		}
	}
}


void ModelGTR::getRateMatrix(double *rate_mat) {
	int nrate = getNumRateEntries();
	memcpy(rate_mat, rates, nrate * sizeof(double));
}

void ModelGTR::setRateMatrix(double* rate_mat)
{
	int nrate = getNumRateEntries();
	memcpy(rates, rate_mat, nrate * sizeof(double));
}

void ModelGTR::getStateFrequency(double *freq) {
	assert(state_freq);
	assert(freq_type != FREQ_UNKNOWN);
	memcpy(freq, state_freq, sizeof(double) * num_states);
    // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
    double sum = 0.0;
    int i;
    for (i = 0; i < num_states; i++) sum += freq[i];
    sum = 1.0/sum;
    for (i = 0; i < num_states; i++) freq[i] *= sum;
}

void ModelGTR::setStateFrequency(double* freq)
{
	assert(state_freq);
	memcpy(state_freq, freq, sizeof(double) * num_states);
}

void ModelGTR::getQMatrix(double *q_mat) {
	double **rate_matrix = (double**) new double[num_states];
	int i, j, k = 0;

	for (i = 0; i < num_states; i++)
		rate_matrix[i] = new double[num_states];

	for (i = 0, k = 0; i < num_states; i++) {
		rate_matrix[i][i] = 0.0;
		for (j = i+1; j < num_states; j++, k++) {
			rate_matrix[i][j] = (state_freq[i] <= ZERO_FREQ || state_freq[j] <= ZERO_FREQ) ? 0 : rates[k];
			rate_matrix[j][i] = rate_matrix[i][j];
		}
	}

	computeRateMatrix(rate_matrix, state_freq, num_states);
	for (i = 0; i < num_states; i++)
		memmove(q_mat + (i*num_states), rate_matrix[i], num_states * sizeof(double));

	for (i = num_states-1; i >= 0; i--)
		delete [] rate_matrix[i];
	delete [] rate_matrix;

}

int ModelGTR::getNDim() { 
	assert(freq_type != FREQ_UNKNOWN);
	int ndim = num_params;
	if (freq_type == FREQ_ESTIMATE) 
		ndim += num_states-1;
	return ndim;
}

int ModelGTR::getNDimFreq() { 
	if (freq_type == FREQ_EMPIRICAL) 
        return num_states-1;
	else if (freq_type == FREQ_CODON_1x4) 
        return 3;
	else if (freq_type == FREQ_CODON_3x4 || freq_type == FREQ_CODON_3x4C) 
        return 9;
    
    return 0;
}

void ModelGTR::scaleStateFreq(bool sum_one) {
	int i;
	if (sum_one) {
		// make the frequencies sum to 1
		double sum = 0.0;
		for (i = 0; i < num_states; i++) sum += state_freq[i];
		for (i = 0; i < num_states; i++) state_freq[i] /= sum;		
	} else {
		// make the last frequency equal to 0.1
		if (state_freq[num_states-1] == 0.1) return;
		assert(state_freq[num_states-1] > 1.1e-6);
		for (i = 0; i < num_states; i++) 
			state_freq[i] /= state_freq[num_states-1]*10.0;
	}
}

void ModelGTR::setVariables(double *variables) {
	int nrate = getNDim();
	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0)
		memcpy(variables+1, rates, nrate*sizeof(double));
	if (freq_type == FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
		int ndim = getNDim();
		memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
        
//		int i, j;
//		for (i = 0, j = 1; i < num_states; i++)
//			if (i != highest_freq_state) {
//				variables[nrate+j] = state_freq[i] / state_freq[highest_freq_state];
//				j++;
//			}
		//scaleStateFreq(false);
//		memcpy(variables+nrate+1, state_freq, (num_states-1)*sizeof(double));
		//scaleStateFreq(true);
	}
}

bool ModelGTR::getVariables(double *variables) {
	int nrate = getNDim();
	int i;
	bool changed = false;
	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0) {
		for (i = 0; i < nrate; i++)
			changed |= (rates[i] != variables[i+1]);
		memcpy(rates, variables+1, nrate * sizeof(double));
	}

	if (freq_type == FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
		int ndim = getNDim();
		for (i = 0; i < num_states-1; i++)
			changed |= (state_freq[i] != variables[i+ndim-num_states+2]);
		memcpy(state_freq, variables+(ndim-num_states+2), (num_states-1)*sizeof(double));

//		memcpy(state_freq, variables+nrate+1, (num_states-1)*sizeof(double));
		//state_freq[num_states-1] = 0.1;
		//scaleStateFreq(true);

//		double sum = 0.0;
//		for (int i = 0; i < num_states-1; i++)
//			sum += state_freq[i];
//		state_freq[num_states-1] = 1.0 - sum;
//		double sum = 1.0;
//		int i, j;
//		for (i = 1; i < num_states; i++)
//			sum += variables[nrate+i];
//		for (i = 0, j = 1; i < num_states; i++)
//			if (i != highest_freq_state) {
//				state_freq[i] = variables[nrate+j] / sum;
//				j++;
//			}
//		state_freq[highest_freq_state] = 1.0/sum;
	}
	return changed;
}

double ModelGTR::targetFunk(double x[]) {
	bool changed = getVariables(x);
	if (state_freq[num_states-1] < 0) return 1.0e+12;
	if (changed) {
		decomposeRateMatrix();
		assert(phylo_tree);
		phylo_tree->clearAllPartialLH();
	}
	return -phylo_tree->computeLikelihood();
}

bool ModelGTR::isUnstableParameters() {
	int nrates = getNumRateEntries();
	int i;
    // NOTE: zero rates are not consider unstable anymore
	for (i = 0; i < nrates; i++)
		if (/*rates[i] < MIN_RATE+TOL_RATE || */rates[i] > MAX_RATE-TOL_RATE)
			return true;
	for (i = 0; i < num_states; i++)
		if (state_freq[i] < MIN_RATE+TOL_RATE)
			return true;
	return false;
}

void ModelGTR::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i, ndim = getNDim();

	for (i = 1; i <= ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = MIN_RATE;
		upper_bound[i] = MAX_RATE;
		bound_check[i] = false;
	}

	if (freq_type == FREQ_ESTIMATE) {
		for (i = ndim-num_states+2; i <= ndim; i++) {
//            lower_bound[i] = MIN_FREQUENCY/state_freq[highest_freq_state];
//			upper_bound[i] = state_freq[highest_freq_state]/MIN_FREQUENCY;
            lower_bound[i]  = MIN_FREQUENCY;
//            upper_bound[i] = 100.0;
            upper_bound[i] = 1.0;
            bound_check[i] = false;
        }
	}
}

double ModelGTR::optimizeParameters(double gradient_epsilon) {
	int ndim = getNDim();
	
	// return if nothing to be optimized
	if (ndim == 0) return 0.0;
    
	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters..." << endl;

	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	double score;

    for (int i = 0; i < num_states; i++)
        if (state_freq[i] > state_freq[highest_freq_state])
            highest_freq_state = i;

	// by BFGS algorithm
	setVariables(variables);
	setBounds(lower_bound, upper_bound, bound_check);
	//packData(variables, lower_bound, upper_bound, bound_check);
//    if (phylo_tree->params->optimize_alg.find("BFGS-B") == string::npos)
        score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));
//    else
//        score = -L_BFGS_B(ndim, variables+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_RATE));

	bool changed = getVariables(variables);
    // BQM 2015-09-07: normalize state_freq
	if (freq_type == FREQ_ESTIMATE) { 
        scaleStateFreq(true);
        changed = true;
    }
    if (changed) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
	
	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



void ModelGTR::decomposeRateMatrix(){
	int i, j, k = 0;

	if (num_params == -1) {
		// manual compute eigenvalues/vectors for F81-style model
		eigenvalues[0] = 0.0;
		double mu = 0.0;
		for (i = 0; i < num_states; i++)
			mu += state_freq[i]*state_freq[i];
		mu = total_num_subst/(1.0 - mu);

		// compute eigenvalues
		for (i = 1; i < num_states; i++)
			eigenvalues[i] = -mu;

//		double *f = new double[num_states];
//		for (i = 0; i < num_states; i++) f[i] = sqrt(state_freq[i]);
		// compute eigenvectors
		memset(eigenvectors, 0, num_states*num_states*sizeof(double));
		memset(inv_eigenvectors, 0, num_states*num_states*sizeof(double));
		eigenvectors[0] = 1.0;
		for (i = 1; i < num_states; i++)
			eigenvectors[i] = -1.0;
//			eigenvectors[i] = f[i]/f[num_states-1];
		for (i = 1; i < num_states; i++) {
			eigenvectors[i*num_states] = 1.0;
			eigenvectors[i*num_states+i] = state_freq[0]/state_freq[i];
		}

		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++)
				inv_eigenvectors[i*num_states+j] = state_freq[j]*eigenvectors[j*num_states+i];
		writeInfo(cout);
		// sanity check
		double *q = new double[num_states*num_states];
		getQMatrix(q);
		double zero;
		for (j = 0; j < num_states; j++) {
			for (i = 0, zero = 0.0; i < num_states; i++) {
				for (k = 0; k < num_states; k++) zero += q[i*num_states+k] * eigenvectors[k*num_states+j];
				zero -= eigenvalues[j] * eigenvectors[i*num_states+j];
				if (fabs(zero) > 1.0e-5) {
					cout << "\nERROR: Eigenvector doesn't satisfy eigenvalue equation! (gap=" << fabs(zero) << ")" << endl;
					abort();
				}
			}
		}
		delete [] q;
	} else {
		double **rate_matrix = new double*[num_states];

		for (i = 0; i < num_states; i++)
			rate_matrix[i] = new double[num_states];

        if (half_matrix) {
            for (i = 0, k = 0; i < num_states; i++) {
                rate_matrix[i][i] = 0.0;
                for (j = i+1; j < num_states; j++, k++) {
                    rate_matrix[i][j] = (state_freq[i] <= ZERO_FREQ || state_freq[j] <= ZERO_FREQ) ? 0 : rates[k];
                    rate_matrix[j][i] = rate_matrix[i][j];
                }
            }
        } else {
            // full matrix
            for (i = 0; i < num_states; i++) {
                memcpy(rate_matrix[i], &rates[i*num_states], num_states*sizeof(double));
                rate_matrix[i][i] = 0.0;
            }
//            IntVector codonid;
//            codonid.reserve(num_states);
//            int baseid[] = {3,1,0,2};
//            for (i=0; i<4; i++)
//                for (j=0; j<4; j++)
//                    for (k=0; k<4; k++)
//                        codonid.push_back(baseid[i]*16+baseid[j]*4+baseid[k]);
//            cout.precision(4);
//            cout << "rate_matrix=" << endl;
//            for (i = 0; i < num_states; i++) {
//                for (j = 0; j < num_states; j++)
//                    cout << " " << rate_matrix[codonid[i]][codonid[j]];
//                cout << endl;
//            }
//            cout << "state_freq=";
//            for (i = 0; i < num_states; i++)
//                cout << " " << state_freq[codonid[i]];
//            cout << endl;
        }
		/* eigensystem of 1 PAM rate matrix */
		eigensystem_sym(rate_matrix, state_freq, eigenvalues, eigenvectors, inv_eigenvectors, num_states);
		//eigensystem(rate_matrix, state_freq, eigenvalues, eigenvectors, inv_eigenvectors, num_states);
		for (i = num_states-1; i >= 0; i--)
			delete [] rate_matrix[i];
		delete [] rate_matrix;
	}
//	for (i = 0; i < num_states; i++)
//		for (j = 0; j < num_states; j++) {
//			int offset = (i*num_states+j)*num_states;
//			double sum = 0.0;
//			for (k = 0; k < num_states; k++) {
//				eigen_coeff[offset+k] = eigenvectors[i*num_states+k] * inv_eigenvectors[k*num_states+j];
//				sum += eigen_coeff[offset+k];
//				//eigen_coeff_derv1[offset+k] = eigen_coeff[offset+k] * eigenvalues[k];
//				//eigen_coeff_derv2[offset+k] = eigen_coeff_derv1[offset+k] * eigenvalues[k];
//			}
//			if (i == j) {
//				if (fabs(sum-1.0) > 1e-6) {
//					cout << "sum = " << sum << endl;
//					assert(0);
//				}
//			}
//			else assert(fabs(sum) < 1e-6);
//		}
//

} 

void ModelGTR::readRates(istream &in) throw(const char*, string) {
	int nrates = getNumRateEntries();
	string str;
	in >> str;
	if (str == "equalrate") {
		for (int i = 0; i < nrates; i++)
			rates[i] = 1.0;
	} else {
		try {
			rates[0] = convert_double(str.c_str());
		} catch (string &str) {
			outError(str);
		}
		if (rates[0] < 0.0)
			throw "Negative rates not allowed";
		for (int i = 1; i < nrates; i++) {
			if (!(in >> rates[i]))
				throw "Rate entries could not be read";
			if (rates[i] < 0.0)
				throw "Negative rates not allowed";
		}
	}
}

void ModelGTR::readRates(string str) throw(const char*) {
	int nrates = getNumRateEntries();
	int end_pos = 0;
	cout << __func__ << " " << str << endl;
	if (str.find("equalrate") != string::npos) {
		for (int i = 0; i < nrates; i++)
			rates[i] = 1.0;
	} else for (int i = 0; i < nrates; i++) {
		int new_end_pos;
		try {
			rates[i] = convert_double(str.substr(end_pos).c_str(), new_end_pos);
		} catch (string &str) {
			outError(str);
		}
		end_pos += new_end_pos;
		if (rates[i] <= 0.0)
			outError("Non-positive rates found");
		if (i == nrates-1 && end_pos < str.length())
			outError("String too long ", str);
		if (i < nrates-1 && end_pos >= str.length())
			outError("Unexpected end of string ", str);
		if (end_pos < str.length() && str[end_pos] != ',')
			outError("Comma to separate rates not found in ", str);
		end_pos++;
	}
	num_params = 0;

}

void ModelGTR::readStateFreq(istream &in) throw(const char*) {
	int i;
	for (i = 0; i < num_states; i++) {
		if (!(in >> state_freq[i])) 
			throw "State frequencies could not be read";
		if (state_freq[i] < 0.0)
			throw "Negative state frequencies found";
	}
	double sum = 0.0;
	for (i = 0; i < num_states; i++) sum += state_freq[i];
	if (fabs(sum-1.0) > 1e-2)
		throw "State frequencies do not sum up to 1.0";
}

void ModelGTR::readStateFreq(string str) throw(const char*) {
	int i;
	int end_pos = 0;
	for (i = 0; i < num_states; i++) {
		int new_end_pos;
		state_freq[i] = convert_double(str.substr(end_pos).c_str(), new_end_pos);
		end_pos += new_end_pos;
		//cout << i << " " << state_freq[i] << endl;
		if (state_freq[i] < 0.0 || state_freq[i] > 1)
			outError("State frequency must be in [0,1] in ", str);
		if (i == num_states-1 && end_pos < str.length())
			outError("Unexpected end of string ", str);
		if (end_pos < str.length() && str[end_pos] != ',' && str[end_pos] != ' ')
			outError("Comma/Space to separate state frequencies not found in ", str);
		end_pos++;
	}
	double sum = 0.0;
	for (i = 0; i < num_states; i++) sum += state_freq[i];
	if (fabs(sum-1.0) > 1e-2)
		outError("State frequencies do not sum up to 1.0 in ", str);
}

void ModelGTR::readParameters(const char *file_name) { 
	try {
		ifstream in(file_name);
		if (in.fail()) {
			outError("Invalid model name ", file_name);
        }
		cout << "Reading model parameters from file " << file_name << endl;
		readRates(in);
		readStateFreq(in);
		in.close();
	}
	catch (const char *str) {
		outError(str);
	} 
	num_params = 0;
	writeInfo(cout);
}


ModelGTR::~ModelGTR() {
	freeMem();
}

void ModelGTR::freeMem()
{
//	int i;
	//delete eigen_coeff_derv2;
	//delete eigen_coeff_derv1;
//	aligned_free(eigen_coeff);

//	for (i = num_states-1; i>=0; i--)
//		delete [] inv_eigenvectors[i];
	aligned_free(inv_eigenvectors);
//	for (i = num_states-1; i>=0; i--)
//		delete [] eigenvectors[i];
	aligned_free(eigenvectors);

	aligned_free(eigenvalues);

	if (rates) delete [] rates;
}

//double *ModelGTR::getEigenCoeff() const
//{
//    return eigen_coeff;
//}
//
double *ModelGTR::getEigenvalues() const
{
    return eigenvalues;
}

double *ModelGTR::getEigenvectors() const
{
    return eigenvectors;
}

double* ModelGTR::getInverseEigenvectors() const {
	return inv_eigenvectors;
}

//void ModelGTR::setEigenCoeff(double *eigenCoeff)
//{
//    eigen_coeff = eigenCoeff;
//}

void ModelGTR::setEigenvalues(double *eigenvalues)
{
    this->eigenvalues = eigenvalues;
}

void ModelGTR::setEigenvectors(double *eigenvectors)
{
    this->eigenvectors = eigenvectors;
}

