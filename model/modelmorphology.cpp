/*
 * modelmorphology.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: minh
 */

#include "modelmorphology.h"

ModelMorphology::ModelMorphology(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree)
: ModelGTR(tree, false)
{
	init(model_name, model_params, freq, freq_params);
}

void ModelMorphology::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	name = model_name;
	full_name = model_name;
	freq = FREQ_EQUAL;
	if (name == "MK") {
		// all were initialized
	} else if (name == "ORDERED") {
		int k = 0;
		// only allow for substitution from state i to state i+1 and back.
		for (int i = 0; i < num_states-1; i++) {
			rates[k++] = 1.0;
			for (int j = i+2; j < num_states; j++, k++)
				rates[k] = 0.0;
		}
	} else {
		// if name does not match, read the user-defined model
		readParameters(model_name);
        num_params = 0;
	}
	ModelGTR::init(freq);
}

void ModelMorphology::readRates(istream &in) throw(const char*, string) {
	int nrates = getNumRateEntries();
	int row = 1, col = 0;
	// since states for protein is stored in lower-triangle, special treatment is needed
	for (int i = 0; i < nrates; i++, col++) {
		if (col == row) {
			row++; col = 0;
		}
		// switch col and row
		int id = col*(2*num_states-col-1)/2 + (row-col-1);
		if (id >= nrates) {
			cout << row << " " << col << endl;
		}
		assert(id < nrates && id >= 0); // make sure that the conversion is correct
		if (!(in >> rates[id]))
			throw name+string(": Rate entries could not be read");
		if (rates[id] < 0.0)
			throw "Negative rates found";
	}
}


ModelMorphology::~ModelMorphology() {
}

