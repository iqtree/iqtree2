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
		outError("Unknown morphological model name");
	}
	ModelGTR::init(freq);
}

ModelMorphology::~ModelMorphology() {
}

