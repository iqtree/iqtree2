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
#include "gtrmodel.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

GTRModel::GTRModel(PhyloTree *tree)
 : SubstModel(tree->aln->num_states), EigenDecomposition(), Optimization()
{
	int i, j;
	int nrate = num_states*(num_states-1)/2;
	int ncoeff = num_states*num_states*num_states;
	
	name = "GTR";
	phylo_tree = tree;
	
	rates = new double[nrate];
	memset(rates, 0, sizeof(double) * nrate);
	state_freq = new double[num_states];
	freq_type = FREQ_UNKNOWN;
	
	eigenvalues = new double[num_states];

	eigenvectors = (double**) new double[num_states];
	for (i = 0; i < num_states; i++)
		eigenvectors[i] = new double[num_states];

	inv_eigenvectors = (double**) new double[num_states];
	for (i = 0; i < num_states; i++)
		inv_eigenvectors[i] = new double[num_states];
		

	eigen_coeff = new double[ncoeff];


	phylo_tree->aln->computeEmpiricalRate(rates);
	//eigen_coeff_derv1 = new double[ncoeff];
	//eigen_coeff_derv2 = new double[ncoeff];
			
}

void GTRModel::setTree(PhyloTree *tree) {
	phylo_tree = tree;
}


void GTRModel::init(StateFreqType type) {
	int i;
	freq_type = type;
	assert(freq_type != FREQ_UNKNOWN);
	switch (freq_type) {
	case FREQ_EQUAL:
	case FREQ_ESTIMATE:
		for (i = 0; i < num_states; i++)
			state_freq[i] = 1.0/num_states;
		break;	
	case FREQ_EMPIRICAL:
		phylo_tree->aln->computeStateFreq(state_freq);
		break;
	}

	decomposeRateMatrix();
}

void GTRModel::writeInfo(ostream &out) {
	if (num_states != 4) return;
	cout << endl;
	cout << "1. Transversion rate from A to C = " << rates[0] << endl;
	cout << "2. Transition   rate from A to G = " << rates[1] << endl;
	cout << "3. Transversion rate from A to T = " << rates[2] << endl;
	cout << "4. Transversion rate from C to G = " << rates[3] << endl;
	cout << "5. Transition   rate from C to T = " << rates[4] << endl;
	cout << "6. Transversion rate from G to T = " << rates[5] << endl;
	if (freq_type != FREQ_ESTIMATE) return;
	cout << "   Base frequency A = " << state_freq[0] << endl;
	cout << "   Base frequency C = " << state_freq[1] << endl;
	cout << "   Base frequency G = " << state_freq[2] << endl;
	cout << "   Base frequency T = " << state_freq[3] << endl;
}


void GTRModel::computeTransMatrix(double time, double *trans_matrix) {
	/* compute P(t) */

	double evol_time = time / total_num_subst;
	double exptime[num_states];
	int i, j, k;

	for (i = 0; i < num_states; i++)
		exptime[i] = exp(evol_time * eigenvalues[i]);

	for (i = 0; i < num_states; i ++) {
		for (j = 0; j < num_states; j ++) {
			double *trans_entry = trans_matrix + (i*num_states+j);
			double *coeff_entry = eigen_coeff + ((i*num_states+j)*num_states);
			*trans_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				*trans_entry += coeff_entry[k] * exptime[k];
			}
			if (*trans_entry < 0.0) {
				*trans_entry = 0.0;
			}
		}
	}
}

void GTRModel::computeTransDerv(double time, double *trans_matrix, 
	double *trans_derv1, double *trans_derv2) 
{
	/* compute P(t) */

	double evol_time = time / total_num_subst;
	double exptime[num_states];
	int i, j, k;

	for (i = 0; i < num_states; i++)
		exptime[i] = exp(evol_time * eigenvalues[i]);

	for (i = 0; i < num_states; i ++) {
		for (j = 0; j < num_states; j ++) {
			int offset = (i*num_states+j);
			double *trans_entry = trans_matrix + offset;
			double *derv1_entry = trans_derv1 + offset;
			double *derv2_entry = trans_derv2 + offset;

			int coeff_offset = offset*num_states;
			double *coeff_entry       = eigen_coeff + coeff_offset;
			//double *coeff_derv1_entry = eigen_coeff_derv1 + coeff_offset;
			//double *coeff_derv2_entry = eigen_coeff_derv2 + coeff_offset;
			*trans_entry = 0.0;
			*derv1_entry = 0.0;
			*derv2_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				double trans = coeff_entry[k] * exptime[k];
				double trans2 = trans * eigenvalues[k];
				*trans_entry += trans;
				*derv1_entry += trans2;
				*derv2_entry += trans2 * eigenvalues[k];
				//*derv1_entry += coeff_derv1_entry[k] * exptime[k];
				//*derv2_entry += coeff_derv2_entry[k] * exptime[k];
			}
			if (*trans_entry < 0.0) {
				*trans_entry = 0.0;
			}
		}
	}
	/*if (verbose_mode == VB_DEBUG) {
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


void GTRModel::getStateFrequency(double *freq) {
	assert(state_freq);
	assert(freq_type != FREQ_UNKNOWN);
	memcpy(freq, state_freq, sizeof(double) * num_states);
}


int GTRModel::getNDim() { 
	assert(freq_type != FREQ_UNKNOWN);
	int ndim = num_states*(num_states-1)/2 - 1; 
	if (freq_type == FREQ_ESTIMATE) 
		ndim += num_states-1;
	return ndim;
}

void GTRModel::setVariables(double *variables) {
	int ndim = getNDim();
	if (freq_type != FREQ_ESTIMATE)
		memcpy(variables+1, rates, ndim*sizeof(double));
	else {
		int nrate = ndim-num_states+1;
		memcpy(variables+1, rates, nrate*sizeof(double));
		memcpy(variables+nrate+1, state_freq, (num_states-1)*sizeof(double));
	}
}

void GTRModel::getVariables(double *variables) {
	int ndim = getNDim();
	if (freq_type != FREQ_ESTIMATE)
		memcpy(rates, variables+1, ndim * sizeof(double));
	else {
		int nrate = ndim-num_states+1;
		double sum = 0.0;
		memcpy(rates, variables+1, nrate*sizeof(double));
		memcpy(state_freq, variables+nrate+1, (num_states-1)*sizeof(double));
		for (int i = 0; i < num_states-1; i++) 
			sum += state_freq[i];
		state_freq[num_states-1] = 1.0 - sum;
	}
}

double GTRModel::targetFunk(double x[]) {
	getVariables(x);
	if (state_freq[num_states-1] <= 0.0) return 1.0e+9;
	decomposeRateMatrix();
	assert(phylo_tree);
	phylo_tree->clearAllPartialLh();
	return -phylo_tree->computeLikelihood();
}


double GTRModel::optimizeParameters() {
	int ndim = getNDim();
	
	// return if nothing to be optimized
	if (ndim == 0) return 0.0;

	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters..." << endl;


	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	int i;
	double score;
	
	// by BFGS algorithm
	setVariables(variables);
	for (i = 1; i <= ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = 1e-4;
		upper_bound[i] = 100.0;
		bound_check[i] = false;
	}
	if (freq_type == FREQ_ESTIMATE) {
		for (i = ndim-num_states+2; i <= ndim; i++) 
			upper_bound[i] = 1.0;
	}
	//packData(variables, lower_bound, upper_bound, bound_check);
	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, 1e-6);

	getVariables(variables);
	decomposeRateMatrix();
	phylo_tree->clearAllPartialLh();
	
	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



void GTRModel::decomposeRateMatrix() {
	double *rate_matrix[num_states];
	int i, j, k = 0;

	for (i = 0; i < num_states; i++)
		rate_matrix[i] = new double[num_states];

	for (i = 0, k = 0; i < num_states; i++) {
		rate_matrix[i][i] = 0.0;
		for (j = i+1; j < num_states; j++, k++) {
			rate_matrix[i][j] = rates[k];
			rate_matrix[j][i] = rates[k];
		}
	}
	/* eigensystem of 1 PAM rate matrix */
	eigensystem_sym(rate_matrix, state_freq, eigenvalues, eigenvectors, inv_eigenvectors, num_states); 

	

	for (i = 0; i < num_states; i++)
		for (j = 0; j < num_states; j++) {
			int offset = (i*num_states+j)*num_states;
			double sum = 0.0;
			for (k = 0; k < num_states; k++) {
				eigen_coeff[offset+k] = eigenvectors[i][k] * inv_eigenvectors[k][j];
				sum += eigen_coeff[offset+k];
				//eigen_coeff_derv1[offset+k] = eigen_coeff[offset+k] * eigenvalues[k];
				//eigen_coeff_derv2[offset+k] = eigen_coeff_derv1[offset+k] * eigenvalues[k];
			}
			if (i == j)
				assert(fabs(sum-1.0) < 1e-6);
			else assert(fabs(sum) < 1e-6);
		}


	for (i = num_states-1; i >= 0; i--)
		delete [] rate_matrix[i];
} 

GTRModel::~GTRModel()
{
	int i, j;
	//delete eigen_coeff_derv2;
	//delete eigen_coeff_derv1;
	delete [] eigen_coeff;

	for (i = num_states-1; i>=0; i--)
		delete [] inv_eigenvectors[i];
	delete [] inv_eigenvectors;
	for (i = num_states-1; i>=0; i--)
		delete [] eigenvectors[i];
	delete [] eigenvectors;

	delete [] eigenvalues;

	if (state_freq) delete [] state_freq;
	if (rates) delete [] rates;
}


