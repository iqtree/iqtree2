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
  ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, bool count_rates);
	static bool validModelName(string model_name);
	void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
private:
	const double **basis;
	int symmetry; // RY->0, WS->1, MK->2
	int model_num; // 0->1.1, etc to 36->12.12
	void setRates();

	static void parseModelName(string model_name, int* model_num, int* symmetry);
	const static double ***BASES;
	const static int *MODEL_PARAMS;
	const static string *SYMMETRY;
	const static string *MODEL_NAMES;
	const static int NUM_RATES;
	const static int NUM_LM_MODELS;
};

#endif /* MODELLIEMARKOV_H_ */
