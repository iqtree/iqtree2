/*
 * modelcodonparametric.cpp
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#include "modelcodonparametric.h"

ModelCodonParametric::ModelCodonParametric(const char *model_name, string model_params,
		StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates) :
		ModelCodon(tree, count_rates)
{
	init(model_name, model_params, freq, freq_params);
}


ModelCodonParametric::~ModelCodonParametric() {
}

void ModelCodonParametric::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	StateFreqType def_freq = FREQ_UNKNOWN;
	name = full_name = model_name;
	string name_upper = model_name;
	for (string::iterator it = name_upper.begin(); it != name_upper.end(); it++)
		(*it) = toupper(*it);
	if (name_upper == "JCC") {
		name = "JCC";
		def_freq = FREQ_EQUAL;
		full_name = "JC-like codon model";
	} else if (name_upper == "MG") {
		initMG94();
	} else if (name_upper == "GY") {
		initGY94();
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
	ModelCodon::init(freq);
}


void ModelCodonParametric::initMG94() {
	/* Muse-Gaut 1994 model with 1 parameters: omega */
	int i,j;
	IntVector group;
	for (i = 0; i < num_states-1; i++) {
		for (j = i+1; j < num_states; j++) {
			if (isMultipleSubst(i, j))
				group.push_back(0); // multiple substitution
			else if (isSynonymous(i, j))
				group.push_back(1); // synonymous substitution
			else
				group.push_back(2); // non-synonymous substitution
		}
	}
	setRateGroup(group);
	// set zero rate for multiple substitution and 1 for synonymous substitution
	setRateGroupConstraint("x0=0,x1=1");
}


void ModelCodonParametric::initGY94() {
	/* Yang-Nielsen 1998 model (also known as Goldman-Yang 1994) with 2 parameters: omega and kappa */
	int i,j;
	IntVector group;
	for (i = 0; i < num_states-1; i++) {
		for (j = i+1; j < num_states; j++) {
			if (isMultipleSubst(i, j))
				group.push_back(0); // multiple substitution
			else if (isSynonymous(i, j)) {
				if (isTransversion(i, j))
					group.push_back(1); // synonymous transversion
				else
					group.push_back(2); // synonymous transition
			} else {
				if (isTransversion(i, j))
					group.push_back(3); // non-synonymous transversion
				else
					group.push_back(4); // non-synonymous transition
			}
		}
	}
	setRateGroup(group);
	// set zero rate for multiple substitution
	// 1 for synonymous transversion
	// and kappa*omega for non-synonymous transition
	setRateGroupConstraint("x0=0,x1=1,x4=x2*x3");
}

void ModelCodonParametric::writeInfo(ostream &out) {
	double *variables = new double[getNDim()+1];
	setVariables(variables);
	if (name == "MG") {
		out << "Nonsynonymous/synonymous ratio (omega): " << variables[1] << endl;
	} else if (name == "GY") {
		out << "Transition/transversion ratio (kappa): " << variables[1] << endl;
		out << "Nonsynonymous/synonymous ratio (omega): " << variables[2] << endl;
	}
	delete [] variables;
}

