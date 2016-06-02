/*
 * modelunrest.cpp
 *
 *  Created on: 24/05/2016
 *      Author: mdw2
 */

#include "modelunrest.h"

ModelUnrest::ModelUnrest(PhyloTree *tree, string model_params, bool count_rates)
	: ModelNonRev(tree)
{
	num_params = getNumRateEntries() - 1;
	model_parameters = new double [num_params];
	for (int i=0; i<= num_params; i++) model_parameters[i] = 1;
	this->setRates();
	/*
	 * I'm not sure how to correctly handle count_rates, so for now I'm just
	 * avoiding the problem. Actual IQTree programmers can fix this.
	 * Whatever happens should leave model_parameters[] and rates[]
	 * consistent with each other.
	 */
	if (count_rates)
		cerr << "WARNING: count_rates=TRUE not implemented in ModelUnrest constructor -- ignored" << endl;
		/* phylo_tree->aln->computeEmpiricalRateNonRev(rates); */
	if (model_params != "") {
		cerr << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
		// TODO: parse model_params into model_parameters, then call setRates().
	}
    name = "UNREST";
    full_name = "Unrestricted model (non-reversible)";
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
