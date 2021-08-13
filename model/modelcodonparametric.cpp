/*
 * modelcodonparametric.cpp
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#include "modelcodonparametric.h"

ModelCodonParametric::ModelCodonParametric
		(const char *model_name, const string& model_params,
		 StateFreqType freq, const string& freq_params, 
		 PhyloTree *tree, bool count_rates) : ModelCodon(tree, tree)
{
	init(model_name, model_params, freq, freq_params, tree);
}

ModelCodonParametric::~ModelCodonParametric() {
}

void ModelCodonParametric::init(const char*   model_name, std::string model_params, 
                                StateFreqType freq,       std::string freq_params,
								PhyloTree*    report_to_tree)
{
	StateFreqType def_freq = StateFreqType::FREQ_UNKNOWN;
	name = full_name = model_name;
	string name_upper = string_to_upper(model_name);
	if (name_upper == "JCC") {
		name = "JCC";
		def_freq = StateFreqType::FREQ_EQUAL;
		full_name = "JC-like codon model";
	} else if (name_upper == "MG") {
		initMG94();
	} else if (name_upper == "GY") {
		initGY94();
	} else {
		//cout << "User-specified model "<< model_name << endl;
		readParameters(model_name, true, report_to_tree);
			//name += " (user-defined)";
	}
	if (freq_params != "") {
		readStateFreq(freq_params, report_to_tree);
	}
	if (model_params != "") {
		readRates(model_params);
	}
	if (freq == StateFreqType::FREQ_UNKNOWN ||  
		def_freq == StateFreqType::FREQ_EQUAL) {
		freq = def_freq;
	}
	super::init(model_name, model_params, freq, 
	            freq_params, report_to_tree);
}

void ModelCodonParametric::initMG94() {
	/* Muse-Gaut 1994 model with 1 parameters: omega */
	IntVector groups;
	for (int i = 0; i < num_states-1; i++) {
		for (int j = i+1; j < num_states; j++) {
			auto attr = rate_attr[i*num_states+j];
			int  group;
			if (( attr & CA_MULTI_NT) != 0 ) {
				group = 0; // multiple substitution
			}
			else if (( attr & CA_SYNONYMOUS) != 0) {
				group = 1 ; // synonymous substitution
			}
			else {
				group = 2 ; // non-synonymous substitution
			}
			groups.push_back(group);
		}
	}
	setRateGroup(groups);
	// set zero rate for multiple substitution and 1 for synonymous substitution
	setRateGroupConstraint("x0=0,x1=1");
}


void ModelCodonParametric::initGY94() {
	/* Yang-Nielsen 1998 model (also known as Goldman-Yang 1994) with 2 parameters: omega and kappa */
	IntVector groups;
	for (int i = 0; i < num_states-1; i++) {
		for (int j = i+1; j < num_states; j++) {
			auto attr = rate_attr[i*num_states+j];
			int  group;
			if (( attr & CA_MULTI_NT) != 0 ) {
				group = 0; // multiple substitution
			} else if (attr & CA_SYNONYMOUS) {
				if ((attr & CA_TRANSVERSION) != 0) {
					group = 1; // synonymous transversion
				} else {
					group = 2; // synonymous transition
				}
			} else {
				if ((attr & CA_TRANSVERSION) != 0) {
					group = 3; // non-synonymous transversion
				} else {
					group = 4; // non-synonymous transition
				}
			}
			groups.push_back(group);
		}
	}
	setRateGroup(groups);
	// set zero rate for multiple substitution
	// 1 for synonymous transversion
	// and kappa*omega for non-synonymous transition
	setRateGroupConstraint("x0=0,x1=1,x4=x2*x3");
}

void ModelCodonParametric::setRateGroup(IntVector& upper_triangle_entry_to_rate) {
	throw "setRateGroup is not implemented";
}

void ModelCodonParametric::setRateGroupConstraint(const char* constraint) {
	throw "setRateGroupConstraint is not implemented";
}

void ModelCodonParametric::writeInfo(ostream &out) {
	DoubleVector variables(getNDim()+1, 0);
	setVariables(variables.data());
	if (name == "MG") {
		out << "Nonsynonymous/synonymous ratio (omega): " 
		    << variables[1] << endl;
	} else if (name == "GY") {
		out << "Transition/transversion ratio (kappa): " 
		    << variables[1] << endl;
		out << "Nonsynonymous/synonymous ratio (omega): " 
		    << variables[2] << endl;
	}
}

