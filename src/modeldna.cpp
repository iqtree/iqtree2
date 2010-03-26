/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "modeldna.h"

ModelDNA::ModelDNA(const char *model_name, StateFreqType freq, PhyloTree *tree)
: GTRModel(tree)
{
	init(model_name, freq);
}


void ModelDNA::init(const char *model_name, StateFreqType freq)
{
	assert(num_states == 4); // make sure that you create model for DNA
	StateFreqType def_freq = FREQ_UNKNOWN;
	if (strcmp(model_name, "JC") == 0 || strcmp(model_name, "JC69") == 0) {
		setRateType("000000");
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "F81") == 0) {
		setRateType("000000");
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "K2P") == 0 || strcmp(model_name, "K80") == 0) {
		setRateType("010010");
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "HKY") == 0 || strcmp(model_name, "HKY85") == 0) {
		setRateType("010010");
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "K3P") == 0 || strcmp(model_name, "K81") == 0) {
		setRateType("012210");
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "K81uf") == 0) {
		setRateType("012210");
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "TN") == 0 || strcmp(model_name, "TN93") == 0) {
		setRateType("010020");
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "TNef") == 0) {
		setRateType("010020");
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "TIM") == 0) {
		setRateType("012230");		
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "TIMef") == 0) {
		setRateType("012230");		
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "TVM") == 0) {
		setRateType("412310");		
		def_freq = FREQ_ESTIMATE;
	} else if (strcmp(model_name, "TVMef") == 0) {
		setRateType("412310");		
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "SYM") == 0) {
		setRateType("123450");
		def_freq = FREQ_EQUAL;
	} else if (strcmp(model_name, "GTR") == 0) {
		setRateType("123450");
		def_freq = FREQ_ESTIMATE;
	} else {
		string str = "Invalid model name: ";
		str += model_name;
		outError(str);
	}
	
	name = model_name;
	if (freq == FREQ_UNKNOWN || def_freq == FREQ_EQUAL) freq = def_freq;
	GTRModel::init(freq);
	//decomposeRateMatrix();
}


void ModelDNA::setRateType(const char *rate_str) {
	char first_type = 127;
	char last_type = 0;
	char t = first_type;
	int num_ch = strlen(rate_str);
	int i, j;

	assert(num_ch == num_states*(num_states-1)/2);
	for (i = 0; i < num_ch; i++) {
		if (rate_str[i] > last_type) last_type = rate_str[i];
		if (rate_str[i] < first_type) first_type = rate_str[i];
	}
	assert(first_type == rate_str[num_ch-1]);

	num_params = last_type - first_type;
	param_spec = "";
	for (i = 0; i < num_ch; i++) {
		param_spec.push_back(rate_str[i]-first_type);
	}
	assert(param_spec.length() == num_ch);
	double *avg_rates = new double[num_params+1];
	int *num_rates = new int[num_params+1];
	memset(avg_rates, 0, sizeof(double) * (num_params+1));
	memset(num_rates, 0, sizeof(int) * (num_params+1));
	for (i = 0; i < param_spec.size(); i++) {
		avg_rates[param_spec[i]] += rates[i];
		num_rates[param_spec[i]]++;
	}
	for (i = 0; i <= num_params; i++)
		avg_rates[i] /= num_rates[i];
	for (i = 0; i < param_spec.size(); i++) {
		rates[i] = avg_rates[param_spec[i]] / avg_rates[0];
	}
	if (verbose_mode >= VB_DEBUG) {
		cout << "Initialized rates: ";
		for (i = 0; i < param_spec.size(); i++) 
			cout << rates[i] << " ";
		cout << endl;
	}
	delete [] num_rates;
	delete [] avg_rates;
}


int ModelDNA::getNDim() {
	assert(freq_type != FREQ_UNKNOWN);
	int ndim = num_params; 
	if (freq_type == FREQ_ESTIMATE) 
		ndim += num_states-1;
	return ndim;
}

void ModelDNA::writeParameters(ostream &out) {
	int i;
	if (freq_type == FREQ_ESTIMATE) {
		for (i = 0; i < num_states; i++)
			out << "\t" << state_freq[i];
	}
	if (num_params == 0) return;
	if (num_params <= 1)
		out << "\t" << rates[1];
	else {
		int nrateout = num_states*(num_states-1)/2 - 1;
		for (i = 0; i < nrateout; i++)
			out << "\t" << rates[i];
	}
}


void ModelDNA::getVariables(double *variables) {
	int i;
	int num_all = param_spec.length();
	for (i = 0; i < num_all; i++)
		if (param_spec[i] > 0)
			rates[i] = variables[param_spec[i]];
	if (freq_type == FREQ_ESTIMATE) {
		int ndim = getNDim();
		memcpy(state_freq, variables+(ndim-num_states+2), (num_states-1)*sizeof(double));
		double sum = 0;
		for (i = 0; i < num_states-1; i++) 
			sum += state_freq[i];
		state_freq[num_states-1] = 1.0 - sum;
	}
}

void ModelDNA::setVariables(double *variables) {
	int num_all = param_spec.length();
	for (int i = 0; i < num_all; i++)
		if (param_spec[i] > 0)
			variables[param_spec[i]] = rates[i];
	if (freq_type == FREQ_ESTIMATE) {
		int ndim = getNDim();
		memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
	}
}
