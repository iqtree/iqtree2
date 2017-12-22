/*
 * modelunrest.cpp
 *
 *  Created on: 24/05/2016
 *      Author: mdw2
 */

#include "modelunrest.h"

ModelUnrest::ModelUnrest(PhyloTree *tree, string model_params)
	: ModelMarkov(tree, false)
{
	num_params = getNumRateEntries() - 1;
	model_parameters = new double [num_params];
	for (int i=0; i< num_params; i++) model_parameters[i] = 1;
	setRates();
	if (model_params != "") {
		cout << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
		// TODO: parse model_params into model_parameters, then call setRates().
	}
    name = "UNREST";
    full_name = "Unrestricted model (non-reversible)";
    ModelMarkov::init(FREQ_ESTIMATE);
}

/* static */ bool ModelUnrest::validModelName(string model_name) {
	return (model_name == "UNREST");
}

void ModelUnrest::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i, ndim = getNDim();

	for (i = 1; i <= ndim; i++) {
		lower_bound[i] = 0.01;
		upper_bound[i] = 100.0;
		bound_check[i] = false;
	}
}

/*
 * Set rates from model_parameters
 */
void ModelUnrest::setRates() {
	// For UNREST, parameters are simply the off-diagonal rate matrix entries
	// (except [4,3] = rates[11], which is constrained to be 1)
	memcpy(rates, model_parameters, num_params*sizeof(double));
	rates[num_params]=1;
}
