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
#include "substmodel.h"
#include "tools.h"

SubstModel::SubstModel(int nstates)
{
	num_states = nstates;
	name = "JC";
	full_name = "JC (Juke and Cantor, 1969)";
}

// here the simplest Juke-Cantor model is implemented, valid for all kind of data (DNA, AA,...)
void SubstModel::computeTransMatrix(double time, double *trans_matrix) {
	double non_diagonal = (1.0 - exp(-time*num_states/(num_states - 1))) / num_states;
	double diagonal = 1.0 - non_diagonal * (num_states - 1);
	int nstates_sqr = num_states * num_states;

	for (int i = 0; i < nstates_sqr; i++)
		if (i % (num_states+1) == 0) 
			trans_matrix[i] = diagonal; 
		else 
			trans_matrix[i] = non_diagonal;
}

double SubstModel::computeTrans(double time, int state1, int state2) {
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

double SubstModel::computeTrans(double time, int state1, int state2, double &derv1, double &derv2) {
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

void SubstModel::getRateMatrix(double *rate_mat) {
	int nrate = getNumRateEntries();
	for (int i = 0; i < nrate; i++)
		rate_mat[i] = 1.0;
}

void SubstModel::getStateFrequency(double *state_freq) {
	double freq = 1.0 / num_states;
	for (int i = 0; i < num_states; i++)
		state_freq[i] = freq;
}

void SubstModel::computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2)
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


double *SubstModel::newTransMatrix() {
	return new double[num_states * num_states];
}

SubstModel::~SubstModel()
{
}



