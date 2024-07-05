/*
 * modelunrest.cpp
 *
 *  Created on: 24/05/2016
 *      Author: mdw2
 */

#include "modelunrest.h"
#include <stdlib.h>
#include <string.h>

ModelUnrest::ModelUnrest(PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params)
	: ModelMarkov(tree, false)
{
    num_params = getNumRateEntries() - 1;
    //ModelMarkov::setReversible in the ModelMarkov
    //constructor sets all the rates to 0.0.  But...
    for (int i=0; i <= num_params; i++) {
        rates[i] = 1.0;
    }
	
	if (model_params != "") {
		//cout << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
		// TODO: parse model_params into model_parameters, then call setRates().
        // detect the seperator
        char separator = ',';
        if (model_params.find('/') != std::string::npos)
            separator = '/';
        
        // parse input into vector
        DoubleVector tmp_rates;
        convert_double_vec_with_distributions(model_params.c_str(), tmp_rates, false, separator);
        
        // validate the number of params (11 or 12)
        if (tmp_rates.size() != num_params && tmp_rates.size() != num_params + 1)
            outError("Model UNREST requires "+convertIntToString(num_params)+" parameters. Please check and try again!");
        
        // set rates from input params
        for (int i = 0; i < tmp_rates.size(); i++) {
            rates[i] = tmp_rates[i];
            
            // check to fix parameters
            fixed_parameters = !Params::getInstance().optimize_from_given_params;
        }
        
        // if the user supplies 11 params -> set the last rate at 1.0
        if (tmp_rates.size() == num_params)
            setRates();
	}
    name = "UNREST";
    full_name = "Unrestricted model (non-reversible)";
    
    // parse state_freqs if specified
    if (freq_params != "")
        outWarning("In the UNREST model, state frequencies should be embedded into the substitution rates. Thus, AliSim skips the user-specified state frequencies.");
    
    ModelMarkov::init(FREQ_ESTIMATE);
    
    // change the state freq type to user defined if users specify the model parameters
    if (model_params != "")
        this->freq_type = FREQ_USER_DEFINED;
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

void ModelUnrest::setRates() {
    // For UNREST, parameters are simply the off-diagonal rate matrix entries
    // (except [4,3] = rates[11], which is constrained to be 1)
    rates[num_params] = 1;
    return;
}

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
