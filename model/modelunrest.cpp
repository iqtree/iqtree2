/*
 * modelunrest.cpp
 *
 *  Created on: 24/05/2016
 *      Author: mdw2
 */

#include "modelunrest.h"
#include <stdlib.h>
#include <string.h>

ModelUnrest::ModelUnrest(PhyloTree *tree, string model_params)
	: ModelMarkov(tree, false)
{
	num_params = getNumRateEntries() - 1;
	//model_parameters = new double [num_params];
	for (int i=0; i< num_params; i++) rates[i] = 1.0;
	//setRates();
	if (model_params != "") {
		int end_pos = 0;
		cout << __func__ << " " << model_params << endl;
		for (int i = 0; i < 12; i++) {
			int new_end_pos;
			try {
				rates[i] = convert_double(model_params.substr(end_pos).c_str(), new_end_pos);
			} catch (string &model_params) {
				outError(model_params);
			}
		
			end_pos += new_end_pos;
			if (rates[i] <= 0.0)
				outError("Non-positive rates found");
			if (i == 11 && end_pos < model_params.length())
				outError("String too long ", model_params);
			if (i < 11 && end_pos >= model_params.length())
				outError("Unexpected end of string ", model_params);
			if (end_pos < model_params.length() && model_params[end_pos] != ',')
				outError("Comma to separate rates not found in ", model_params);
			end_pos++;
		}
		num_params = 0;
		writeInfo(cout);
	}
    name = "UNREST";
    full_name = "Unrestricted model (non-reversible)";
    ModelMarkov::init(FREQ_ESTIMATE);
}

void ModelUnrest::writeInfo(ostream &out) {
		out << "UNREST rate values:";
		out << "  A-C: " << rates[0];
		out << "  A-G: " << rates[1];
		out << "  A-T: " << rates[2];
		out << "  C-A: " << rates[3];
		out << "  C-G: " << rates[4];
		out << "  C-T: " << rates[5];
		out << "  G-A: " << rates[6];
		out << "  G-C: " << rates[7];
		out << "  G-T: " << rates[8];
		out << "  T-A: " << rates[9];
		out << "  T-C: " << rates[10];
		out << "  T-G: " << rates[11];
		out << endl;
}

/* static */ bool ModelUnrest::validModelName(string model_name) {
	return (model_name == "UNREST");
}

void ModelUnrest::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i, ndim = getNDim();

	for (i = 1; i <= ndim; i++) {
		lower_bound[i] = MIN_RATE;
		upper_bound[i] = MAX_RATE;
		bound_check[i] = false;
	}
}



/*
void ModelUnrest::setRates() {
	// For UNREST, parameters are simply the off-diagonal rate matrix entries
	// (except [4,3] = rates[11], which is constrained to be 1)
	memcpy(rates, model_parameters, num_params*sizeof(double));
	rates[num_params]=1;
}
*/

void ModelUnrest::setStateFrequency(double* freq) {
    // DOES NOTHING
}

void ModelUnrest::startCheckpoint() {
    checkpoint->startStruct("ModelUnrest");
}

void ModelUnrest::saveCheckpoint() {
    startCheckpoint();
    if (!fixed_parameters)
        CKP_ARRAY_SAVE(getNumRateEntries(), rates);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelUnrest::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();
    startCheckpoint();
    if (!fixed_parameters)
        CKP_ARRAY_RESTORE(getNumRateEntries(), rates);
    endCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}
