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

ModelDNA::ModelDNA(const char *model_name, StateFreqType freq, PhyloTree *tree, bool count_rates)
: GTRModel(tree, count_rates)
{
	init(model_name, freq);
}


void ModelDNA::init(const char *model_name, StateFreqType freq)
{
	assert(num_states == 4); // make sure that you create model for DNA
	StateFreqType def_freq = FREQ_UNKNOWN;
	name = model_name;
	full_name = model_name;
	string name_upper = model_name;
	for (string::iterator it = name_upper.begin(); it != name_upper.end(); it++)
		(*it) = toupper(*it);
	if (name_upper == "JC" || name_upper == "JC69") {
		name = "JC";
		setRateType("000000");
		def_freq = FREQ_EQUAL;
		full_name = "JC (Juke and Cantor, 1969)";
	} else if (name_upper == "F81") {
		name = "F81";
		setRateType("000000");
		def_freq = FREQ_ESTIMATE;
		full_name = "F81 (Felsenstein, 1981)";
	} else if (name_upper == "K2P" || name_upper == "K80") {
		name = "K2P";
		setRateType("010010");
		def_freq = FREQ_EQUAL;
		full_name = "K2P (Kimura, 1980)";
	} else if (name_upper == "HKY" || name_upper == "HKY85") {
		name = "HKY";
		setRateType("010010");
		def_freq = FREQ_ESTIMATE;
		full_name = "HKY (Hasegawa, Kishino and Yano, 1985)";
	} else if (name_upper == "K3P" || name_upper == "K81") {
		name = "K3P";
		setRateType("012210");
		def_freq = FREQ_EQUAL;
		full_name = "K3P (Kimura, 1981)";
	} else if (name_upper == "K81UF" || name_upper == "K81U") {
		name = "K3Pu";
		setRateType("012210");
		def_freq = FREQ_ESTIMATE;
		full_name = "K3P unequal frequencies (Kimura, 1981)";
	} else if (name_upper == "TN" || name_upper == "TRN" || name_upper == "TN93") {
		name = "TN";
		setRateType("010020");
		def_freq = FREQ_ESTIMATE;
		full_name = "TN (Tamura and Nei, 1993)";
	} else if (name_upper == "TNEF" || name_upper == "TRNEF" || name_upper == "TNE" || name_upper == "TRNE") {
		name = "TNe";
		setRateType("010020");
		def_freq = FREQ_EQUAL;
		full_name = "TN equal frequencies (Tamura and Nei, 1993)";
	} else if (name_upper == "TPM2") {
		name = "TPM2";
		setRateType("121020");
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM2 ()";
	} else if (name_upper == "TPM2U" || name_upper == "TPM2UF") {
		name = "TPM2u";
		setRateType("121020");
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM2 unequal frequencies ()";
	} else if (name_upper == "TPM3") {
		name = "TPM3";
		setRateType("120120");
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM3 ()";
	} else if (name_upper == "TPM3U" || name_upper == "TPM3UF") {
		name = "TPM3u";
		setRateType("120120");
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM3 unequal frequencies ()";
	} else if (name_upper == "TIM" || name_upper == "TIM1") {
		name = "TIM";
		setRateType("012230");		
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM ()";
	} else if (name_upper == "TIMEF" || name_upper == "TIME" || name_upper == "TIM1EF" || name_upper == "TIM1E") {
		name = "TIMe";
		setRateType("012230");		
		def_freq = FREQ_EQUAL;
		full_name = "TIM equal frequencies";
	} else if (name_upper == "TIM2") {
		name = "TIM2";
		setRateType("121030");
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM2 ()";
	} else if (name_upper == "TIM2EF" || name_upper == "TIM2E") {
		name = "TIM2e";
		setRateType("121030");
		def_freq = FREQ_EQUAL;
		full_name = "TIM2 equal frequencies";
	} else if (name_upper == "TIM3") {
		name = "TIM3";
		setRateType("120130");
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM3 ()";
	} else if (name_upper == "TIM3EF" || name_upper == "TIM3E") {
		name = "TIM3e";
		setRateType("120130");
		def_freq = FREQ_EQUAL;
		full_name = "TIM3 equal frequencies";
	} else if (name_upper == "TVM") {
		name = "TVM";
		setRateType("412310");		
		def_freq = FREQ_ESTIMATE;
		full_name = "TVM";
	} else if (name_upper == "TVMEF" || name_upper == "TVME") {
		name = "TVMe";
		setRateType("412310");		
		def_freq = FREQ_EQUAL;
		full_name = "TVM equal frequencies";
	} else if (name_upper == "SYM") {
		name = "SYM";
		setRateType("123450");
		def_freq = FREQ_EQUAL;
		full_name = "SYM (Zharkihk, 1994)";
	} else if (name_upper == "GTR" || name_upper == "REV") {
		name = "GTR";
		setRateType("123450");
		def_freq = FREQ_ESTIMATE;
		full_name = "GTR (Tavare, 1986)";
	} else {
		//cout << "User-specified model "<< model_name << endl;
		if (setRateType(model_name))
			def_freq = FREQ_ESTIMATE;
		else {
			readParameters(model_name);
			//name += " (user-defined)";
		}
	}
	
	if (freq == FREQ_UNKNOWN ||  def_freq == FREQ_EQUAL) freq = def_freq;
	GTRModel::init(freq);
}


bool ModelDNA::setRateType(const char *rate_str) {
	char first_type = 127;
	char last_type = 0;
	//char t = first_type;
	int num_ch = strlen(rate_str);
	int i;

	if (num_ch != getNumRateEntries()) {
		//outError("Model specification has wrong length!");
		return false;
	}
	if (rate_str[num_ch-1] != '0') {
		//outError("Model specification must end with '0'");
		return false;
	}
	for (i = 0; i < num_ch; i++) {
		if (rate_str[i] > last_type) last_type = rate_str[i];
		if (rate_str[i] < first_type) first_type = rate_str[i];
	}
	if (first_type != rate_str[num_ch-1]) {
		//outError("Model specification must contain digits!");
		return false;
	}

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
		avg_rates[(int)param_spec[i]] += rates[i];
		num_rates[(int)param_spec[i]]++;
	}
	for (i = 0; i <= num_params; i++)
		avg_rates[i] /= num_rates[i];
	for (i = 0; i < param_spec.size(); i++) {
		rates[i] = avg_rates[(int)param_spec[i]] / avg_rates[0];
	}
	if (verbose_mode >= VB_DEBUG) {
		cout << "Initialized rates: ";
		for (i = 0; i < param_spec.size(); i++) 
			cout << rates[i] << " ";
		cout << endl;
	}
	delete [] num_rates;
	delete [] avg_rates;
	return true;
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
		int nrateout = getNumRateEntries() - 1;
		for (i = 0; i < nrateout; i++)
			out << "\t" << rates[i];
	}
}


void ModelDNA::getVariables(double *variables) {
	int i;
	int num_all = param_spec.length();
	for (i = 0; i < num_all; i++)
		if (param_spec[i] > 0)
			rates[i] = variables[(int)param_spec[i]];
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
			variables[(int)param_spec[i]] = rates[i];
	if (freq_type == FREQ_ESTIMATE) {
		int ndim = getNDim();
		memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
	}
}
