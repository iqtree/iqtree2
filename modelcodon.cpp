/*
 * modelcodon.cpp
 *
 *  Created on: May 24, 2013
 *      Author: minh
 */

#include "modelcodon.h"

ModelCodon::ModelCodon(const char *model_name, string model_params,
		StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates) : GTRModel(tree, count_rates)
{
	init(model_name, model_params, freq, freq_params);
}


void ModelCodon::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	assert(num_states == 61); // make sure that you create model for codon
	StateFreqType def_freq = FREQ_UNKNOWN;
	name = model_name;
	full_name = model_name;
	string name_upper = model_name;
	for (string::iterator it = name_upper.begin(); it != name_upper.end(); it++)
		(*it) = toupper(*it);
	if (name_upper == "JC61") {
		name = "JC61";
		def_freq = FREQ_EQUAL;
		full_name = "JC (Juke and Cantor, 1969)";
	} else {
		//cout << "User-specified model "<< model_name << endl;
		readParameters(model_name);
			//name += " (user-defined)";
	}

	if (freq_params != "") {
		readStateFreq(freq_params);
	}
	if (model_params != "") {
		readRates(model_params);
	}

	if (freq == FREQ_UNKNOWN ||  def_freq == FREQ_EQUAL) freq = def_freq;
	GTRModel::init(freq);
}

ModelCodon::~ModelCodon() {
}

