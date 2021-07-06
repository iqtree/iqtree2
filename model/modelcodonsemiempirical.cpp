/*
 * modelcodonsemiempirical.cpp
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#include "modelcodonsemiempirical.h"

ModelCodonSemiEmpirical::ModelCodonSemiEmpirical
		(const char *model_name, string model_params,
		 StateFreqType freq, string freq_params, 
		 PhyloTree *tree, bool count_rates) :
		ModelCodonEmpirical(model_name, model_params, freq, 
		                   freq_params, tree, count_rates),
		parametric(model_name, model_params, freq, 
		           freq_params, tree, count_rates)
{
	init(model_name, model_params, freq, freq_params, tree);
}


ModelCodonSemiEmpirical::~ModelCodonSemiEmpirical() {
}


void ModelCodonSemiEmpirical::init(const char *model_name, string model_params, 
                                   StateFreqType freq, string freq_params,
								   PhyloTree* report_to_tree) {
	parametric.name = full_name = model_name;
	size_t pos = name.find('+');
	ASSERT(pos != string::npos);
	if (name.substr(0,3) == "ECM") {
		init(name.substr(0,pos).c_str(), "", StateFreqType::FREQ_USER_DEFINED, 
		     "", report_to_tree);
		parametric.init(name.substr(pos).c_str(), model_params, freq, 
		                 freq_params, report_to_tree);
	}

}
