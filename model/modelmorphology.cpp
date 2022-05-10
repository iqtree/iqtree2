/*
 * modelmorphology.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: minh
 */

#include "modelmorphology.h"

ModelMorphology::ModelMorphology(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree)
: ModelMarkov(tree)
{
	init(model_name, model_params, freq, freq_params);
}

void ModelMorphology::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	name = model_name;
	full_name = model_name;
	if (name == "MK") {
		// all were initialized
        num_params = 0;
	} else if (name == "ORDERED") {
		int i, j, k = 0;
		// only allow for substitution from state i to state i+1 and back.
		for (i = 0; i < num_states-1; i++) {
			rates[k++] = 1.0;
			for (j = i+2; j < num_states; j++, k++)
				rates[k] = 0.0;
		}
        num_params = 0;
    } else if (name == "GTR" || name == "GTRX") {
        outWarning("GTRX multistate model will estimate " + convertIntToString(getNumRateEntries()-1) + " substitution rates that might be overfitting!");
        outWarning("Please only use GTRX with very large data and always test for model fit!");
        name = "GTRX";
	} else {
		// if name does not match, read the user-defined model
		readParameters(model_name);
        num_params = 0;
        freq = FREQ_USER_DEFINED;
	}
    
    // parse user-specified state frequencies (if any)
    if (freq_params != "")
    {
        freq_type = FREQ_USER_DEFINED;
        readStateFreq(freq_params);
    }
    
	ModelMarkov::init(freq);
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
        
        string tmp_value;
        in >> tmp_value;
        if (tmp_value.length() == 0)
            throw name+string(": Rate entries could not be read");
        rates[id] = convert_double_with_distribution(tmp_value.c_str());
        
		if (rates[id] < 0.0)
			throw "Negative rates found";
	}
}

int ModelMorphology::getNDim() {
    int ndim = num_params;
    if (freq_type == FREQ_ESTIMATE)
        ndim += num_states-1;
    return ndim;
}

ModelMorphology::~ModelMorphology() {
}

void ModelMorphology::startCheckpoint() {
    checkpoint->startStruct("ModelMorph");
}

void ModelMorphology::saveCheckpoint() {
    startCheckpoint();
    if (num_params > 0)
        CKP_ARRAY_SAVE(getNumRateEntries(), rates);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelMorphology::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();
    startCheckpoint();
    if (num_params > 0)
        CKP_ARRAY_RESTORE(getNumRateEntries(), rates);
    endCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

string ModelMorphology::getNameParams(bool show_fixed_params) {
    if (num_params == 0) return name;
    ostringstream retname;
    retname << name << '{';
    int nrates = getNumRateEntries();
    for (int i = 0; i < nrates; i++) {
        if (i>0) retname << ',';
        retname << rates[i];
    }
    retname << '}';
    getNameParamsFreq(retname);
    return retname.str();
}

void ModelMorphology::writeParameters(ostream &out) {
    int i;
    if (freq_type == FREQ_ESTIMATE) {
        for (i = 0; i < num_states; i++)
            out << "\t" << state_freq[i];
    }
    if (num_params == 0) return;
    int nrateout = getNumRateEntries() - 1;
    for (i = 0; i < nrateout; i++)
        out << "\t" << rates[i];
}

void ModelMorphology::writeInfo(ostream &out) {
    if (num_params > 0) {
        out << "Rate parameters:";
        int nrate = getNumRateEntries();
        for (int i = 0; i < nrate; i++)
            out << " " << rates[i];
        out << endl;
    }
    if (freq_type != FREQ_EQUAL) {
        out << "State frequencies:";
        for (int i = 0; i < num_states; i++)
            out << " " << state_freq[i];
        out << endl;
    }
}
