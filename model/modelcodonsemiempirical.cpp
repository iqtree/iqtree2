/*
 * modelcodonsemiempirical.cpp
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#include "modelcodonsemiempirical.h"

ModelCodonSemiEmpirical::ModelCodonSemiEmpirical(const char *model_name, string model_params,
		StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates) :
		ModelCodon(tree, count_rates)
{
	init(model_name, model_params, freq, freq_params);
}


ModelCodonSemiEmpirical::~ModelCodonSemiEmpirical() {
}


void ModelCodonSemiEmpirical::init(const char *model_name, string model_params, StateFreqType freq, string freq_params) {
	name = full_name = model_name;
	size_t pos = name.find('+');
	ASSERT(pos != string::npos);
	if (name.substr(0,3) == "ECM") {
		ModelCodonEmpirical::init(name.substr(0,pos), "", FREQ_USER_DEFINED, "");
		ModelCodonParametric::init(name.substr(pos), model_params, freq, freq_params);
	}

}
