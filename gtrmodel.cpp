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


GTRModel::GTRModel(PhyloTree *tree, bool count_rates)
 : ModelSubst(tree->aln->num_states), EigenDecomposition()
{
	int i;
	int nrate = getNumRateEntries();
	int ncoeff = num_states*num_states*num_states;
	
	name = "GTR";
	full_name = "GTR (Tavare, 1986)";
	phylo_tree = tree;
	
	rates = new double[nrate];
	memset(rates, 0, sizeof(double) * nrate);

	freq_type = FREQ_UNKNOWN;
	
	eigenvalues = new double[num_states];

	eigenvectors = (double**) new double[num_states];
	for (i = 0; i < num_states; i++)
		eigenvectors[i] = new double[num_states];

	inv_eigenvectors = (double**) new double[num_states];
	for (i = 0; i < num_states; i++)
		inv_eigenvectors[i] = new double[num_states];
		
	eigen_coeff = new double[ncoeff];

	if (count_rates) 
		phylo_tree->aln->computeEmpiricalRate(rates);
	else
		for (i=0; i < nrate; i++) rates[i] = 1.0;
	//eigen_coeff_derv1 = new double[ncoeff];
	//eigen_coeff_derv2 = new double[ncoeff];
	num_params = getNumRateEntries() - 1;
}

void GTRModel::setTree(PhyloTree *tree) {
	phylo_tree = tree;
}

string GTRModel::getNameParams() {
	ostringstream retname;
	retname << "GTR";
	if (num_states != 4) retname << num_states;
	retname << '{';
	int nrates = getNumRateEntries();
	for (int i = 0; i < nrates; i++) {
		if (i>0) retname << ',';
		retname << rates[i];
	}
	retname << '}';
	return retname.str();
}

void GTRModel::init(StateFreqType type) {
	//if (type == FREQ_UNKNOWN) return;
	int i;
	freq_type = type;
	assert(freq_type != FREQ_UNKNOWN);
	switch (freq_type) {
	case FREQ_EQUAL:
		for (i = 0; i < num_states; i++)
			state_freq[i] = 1.0/num_states;
		break;	
	case FREQ_ESTIMATE:
	case FREQ_EMPIRICAL:
		phylo_tree->aln->computeStateFreq(state_freq);
		break;
	case FREQ_USER_DEFINED:
		if (state_freq[0] == 0.0) outError("State frequencies not specified");
		break;
	default: break;
	}
	decomposeRateMatrix();
}

void GTRModel::writeInfo(ostream &out) {
	if (num_states != 4) return;
	out << "Rate parameters:" << endl;
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
	out << "Base frequencies: " << endl;
	out << "  A: " << state_freq[0];
	out << "  C: " << state_freq[1];
	out << "  G: " << state_freq[2];
	out << "  T: " << state_freq[3];
	out << endl;
	//out.unsetf(ios::fixed);
}

void GTRModel::computeTransMatrix(double time, double *trans_matrix) {
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
			double *coeff_entry = eigen_coeff + ((row_offset+j)*num_states);
			*trans_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				*trans_entry += coeff_entry[k] * exptime[k];
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

void GTRModel::computeTransMatrixFreq(double time, double* trans_matrix)
{
	computeTransMatrix(time, trans_matrix);
	for (int state1 = 0; state1 < num_states; state1++) {
		double *trans_mat_state = trans_matrix + (state1 * num_states);
		for (int state2 = 0; state2 < num_states; state2++)
			trans_mat_state[state2] *= state_freq[state1];
	}
}

double GTRModel::computeTrans(double time, int state1, int state2) {
	double evol_time = time / total_num_subst;
	int i;

	double *coeff_entry = eigen_coeff + ((state1*num_states+state2)*num_states);
	double trans_prob = 0.0;
	for (i = 0; i < num_states; i++) {
		trans_prob += coeff_entry[i] * exp(evol_time * eigenvalues[i]);
	}
	return trans_prob;
}

double GTRModel::computeTrans(double time, int state1, int state2, double &derv1, double &derv2) {
	double evol_time = time / total_num_subst;
	int i;

	double *coeff_entry = eigen_coeff + ((state1*num_states+state2)*num_states);
	double trans_prob = 0.0;
	derv1 = derv2 = 0.0;
	for (i = 0; i < num_states; i++) {
		double trans = coeff_entry[i] * exp(evol_time * eigenvalues[i]);
		double trans2 = trans * eigenvalues[i];
		trans_prob += trans;
		derv1 += trans2;
		derv2 += trans2 * eigenvalues[i];
	}
	return trans_prob;
}


void GTRModel::computeTransDerv(double time, double *trans_matrix, 
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

			int coeff_offset = offset*num_states;
			double *coeff_entry       = eigen_coeff + coeff_offset;
			*trans_entry = 0.0;
			*derv1_entry = 0.0;
			*derv2_entry = 0.0;
			for (k = 0; k < num_states; k ++) {
				double trans = coeff_entry[k] * exptime[k];
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

void GTRModel::computeTransDervFreq(double time, double rate_val, double* trans_matrix, double* trans_derv1, double* trans_derv2)
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


void GTRModel::getRateMatrix(double *rate_mat) {
	int nrate = getNumRateEntries();
	memcpy(rate_mat, rates, nrate * sizeof(double));
}

void GTRModel::setRateMatrix(double* rate_mat)
{
	int nrate = getNumRateEntries();
	memcpy(rates, rate_mat, nrate * sizeof(double));
}

void GTRModel::getStateFrequency(double *freq) {
	assert(state_freq);
	assert(freq_type != FREQ_UNKNOWN);
	memcpy(freq, state_freq, sizeof(double) * num_states);
}

void GTRModel::setStateFrequency(double* freq)
{
	assert(state_freq);
	memcpy(state_freq, freq, sizeof(double) * num_states);
}

void GTRModel::getQMatrix(double *q_mat) {
	double **rate_matrix = (double**) new double[num_states];
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

	computeRateMatrix(rate_matrix, state_freq, num_states);
	for (i = 0; i < num_states; i++)
		memmove(q_mat + (i*num_states), rate_matrix[i], num_states * sizeof(double));

	for (i = num_states-1; i >= 0; i--)
		delete [] rate_matrix[i];
	delete [] rate_matrix;

}

int GTRModel::getNDim() { 
	assert(freq_type != FREQ_UNKNOWN);
	int ndim = num_params;
	if (freq_type == FREQ_ESTIMATE) 
		ndim += num_states-1;
	return ndim;
}


void GTRModel::scaleStateFreq(bool sum_one) {
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

void GTRModel::setVariables(double *variables) {
	int nrate = getNDim();
	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0)
		memcpy(variables+1, rates, nrate*sizeof(double));
	if (freq_type == FREQ_ESTIMATE) {
		//scaleStateFreq(false);
		memcpy(variables+nrate+1, state_freq, (num_states-1)*sizeof(double));
		//scaleStateFreq(true);
	}
}

void GTRModel::getVariables(double *variables) {
	int nrate = getNDim();
	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0)
		memcpy(rates, variables+1, nrate * sizeof(double));

	if (freq_type == FREQ_ESTIMATE) {
		memcpy(state_freq, variables+nrate+1, (num_states-1)*sizeof(double));
		//state_freq[num_states-1] = 0.1;
		//scaleStateFreq(true);

		double sum = 0.0;
		for (int i = 0; i < num_states-1; i++) 
			sum += state_freq[i];
		state_freq[num_states-1] = 1.0 - sum;
	}
}

double GTRModel::targetFunk(double x[]) {
	getVariables(x);
	if (state_freq[num_states-1] < 1e-4) return 1.0e+12;
	decomposeRateMatrix();
	assert(phylo_tree);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}


double GTRModel::optimizeParameters(double epsilon) {
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
	int i;
	double score;

	// by BFGS algorithm
	setVariables(variables);
	for (i = 1; i <= ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = MIN_RATE;
		upper_bound[i] = MAX_RATE;
		bound_check[i] = false;
	}

	if (freq_type == FREQ_ESTIMATE) {
		for (i = ndim-num_states+2; i <= ndim; i++) 
			upper_bound[i] = 1.0;
	}
	//packData(variables, lower_bound, upper_bound, bound_check);
	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(epsilon, TOL_RATE));

	getVariables(variables);
	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(true);
	decomposeRateMatrix();
	phylo_tree->clearAllPartialLH();
	
	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



void GTRModel::decomposeRateMatrix(){
	double **rate_matrix = (double**) new double[num_states];
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
	//eigensystem(rate_matrix, state_freq, eigenvalues, eigenvectors, inv_eigenvectors, num_states);  

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
	delete [] rate_matrix;
} 

void GTRModel::readRates(istream &in) throw(const char*) {
	int nrates = getNumRateEntries();
	for (int i = 0; i < nrates; i++) {
		if (!(in >> rates[i]))
			throw "Rate entries could not be read";
		if (rates[i] < 0.0)
			throw "Negative rates found";
	}
}

void GTRModel::readRates(string str) throw(const char*) {
	int nrates = getNumRateEntries();
	int end_pos = 0;
	cout << __func__ << " " << str << endl;
	for (int i = 0; i < nrates; i++) {
		int new_end_pos;
		try {
			rates[i] = convert_double(str.substr(end_pos).c_str(), new_end_pos);
		} catch (string str) {
			outError(str);
		}
		end_pos += new_end_pos;
		if (rates[i] <= 0.0)
			outError("Negative rates found");
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

void GTRModel::readStateFreq(istream &in) throw(const char*) {
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

void GTRModel::readStateFreq(string str) throw(const char*) {
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
		if (end_pos < str.length() && str[end_pos] != ',')
			outError("Comma to separate state frequencies not found in ", str);
		end_pos++;
	}
	double sum = 0.0;
	for (i = 0; i < num_states; i++) sum += state_freq[i];
	if (fabs(sum-1.0) > 1e-2)
		outError("State frequencies do not sum up to 1.0 in ", str);
}

void GTRModel::readParameters(const char *file_name) { 
	try {
		cout << "Reading model parameters from file " << file_name << endl;
		ifstream in(file_name);
		readRates(in);
		readStateFreq(in);
	}
	catch (const char *str) {
		outError(str);
	} 
	num_params = 0;
	writeInfo(cout);
}


GTRModel::~GTRModel() {
	freeMem();
}

void GTRModel::freeMem()
{
	int i;
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

	if (rates) delete [] rates;
}

double *GTRModel::getEigenCoeff() const
{
    return eigen_coeff;
}

double *GTRModel::getEigenvalues() const
{
    return eigenvalues;
}

double **GTRModel::getEigenvectors() const
{
    return eigenvectors;
}

double** GTRModel::getInverseEigenvectors() const {
	return inv_eigenvectors;
}

void GTRModel::setEigenCoeff(double *eigenCoeff)
{
    eigen_coeff = eigenCoeff;
}

void GTRModel::setEigenvalues(double *eigenvalues)
{
    this->eigenvalues = eigenvalues;
}

void GTRModel::setEigenvectors(double **eigenvectors)
{
    this->eigenvectors = eigenvectors;
}

