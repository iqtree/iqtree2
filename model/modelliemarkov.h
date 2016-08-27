/*
 * modelliemarkov.h
 *
 *  Created on: 24/05/2016
 *      Author: Michael Woodhams
 */

#ifndef MODELLIEMARKOV_H_
#define MODELLIEMARKOV_H_

#include "modelnonrev.h"

class ModelLieMarkov: public ModelNonRev {
public:
        ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, bool count_rates = false);
        virtual ~ModelLieMarkov();

	/**
		this function is served for model testing
		@param model_name name of the model
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

	static bool validModelName(string model_name);
	void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
	
	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

    /** decompose rate matrix using closed formula derived by Michael Woodhams */
    void decomposeRateMatrixClosedForm();

    /** decompose rate matrix using Eigen library */
    void decomposeRateMatrixEigen3lib();

	/**
		compute the transition probability matrix.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix);

protected:
	const double **basis;
	int symmetry; // RY->0, WS->1, MK->2
	int model_num; // 0->1.1, etc to 36->12.12
	virtual void setRates();

	static void parseModelName(string model_name, int* model_num, int* symmetry);
	/*
	const static double ***BASES;
	const static int *MODEL_PARAMS;
	const static string *SYMMETRY;
	const static string *MODEL_NAMES;
	const static int NUM_RATES;
	const static int NUM_LM_MODELS;
	*/
};

#endif /* MODELLIEMARKOV_H_ */
