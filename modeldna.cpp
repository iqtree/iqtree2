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

ModelDNA::ModelDNA(PhyloTree *tree, bool count_rates)
: GTRModel(tree, count_rates)
{
}

ModelDNA::ModelDNA(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates)
: GTRModel(tree, count_rates)
{
	init(model_name, model_params, freq, freq_params);
}

string getDNAModelInfo(string model_name, string &full_name, string &rate_type, StateFreqType &def_freq) {
	string name_upper = model_name;
	for (string::iterator it = name_upper.begin(); it != name_upper.end(); it++)
		(*it) = toupper(*it);
	string name = model_name;
	full_name = name;
	rate_type = "";
	def_freq = FREQ_UNKNOWN;
	if (name_upper == "JC" || name_upper == "JC69") {
		name = "JC";
		rate_type = "000000";
		def_freq = FREQ_EQUAL;
		full_name = "JC (Juke and Cantor, 1969)";
	} else if (name_upper == "F81") {
		name = "F81";
		rate_type = "000000";
		def_freq = FREQ_ESTIMATE;
		full_name = "F81 (Felsenstein, 1981)";
	} else if (name_upper == "K2P" || name_upper == "K80") {
		name = "K2P";
		rate_type = "010010";
		def_freq = FREQ_EQUAL;
		full_name = "K2P (Kimura, 1980)";
	} else if (name_upper == "HKY" || name_upper == "HKY85") {
		name = "HKY";
		rate_type = "010010";
		def_freq = FREQ_ESTIMATE;
		full_name = "HKY (Hasegawa, Kishino and Yano, 1985)";
	} else if (name_upper == "K3P" || name_upper == "K81" || name_upper=="TPM1") {
		name = "K3P";
		rate_type = "012210";
		def_freq = FREQ_EQUAL;
		full_name = "K3P (Kimura, 1981)";
	} else if (name_upper == "K81UF" || name_upper == "K81U" || name_upper == "K3PU" ||
			name_upper == "K3PUF" || name_upper=="TPM1UF" || name_upper=="TPM1U") {
		name = "K3Pu";
		rate_type = "012210";
		def_freq = FREQ_ESTIMATE;
		full_name = "K3P unequal frequencies (Kimura, 1981)";
	} else if (name_upper == "TN" || name_upper == "TRN" || name_upper == "TN93") {
		name = "TN";
		rate_type = "010020";
		def_freq = FREQ_ESTIMATE;
		full_name = "TN (Tamura and Nei, 1993)";
	} else if (name_upper == "TNEF" || name_upper == "TRNEF" || name_upper == "TNE" || name_upper == "TRNE") {
		name = "TNe";
		rate_type = "010020";
		def_freq = FREQ_EQUAL;
		full_name = "TN equal frequencies (Tamura and Nei, 1993)";
	} else if (name_upper == "TPM2") {
		name = "TPM2";
		rate_type = "121020";
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM2 ()";
	} else if (name_upper == "TPM2U" || name_upper == "TPM2UF") {
		name = "TPM2u";
		rate_type = "121020";
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM2 unequal frequencies ()";
	} else if (name_upper == "TPM3") {
		name = "TPM3";
		rate_type = "120120";
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM3 ()";
	} else if (name_upper == "TPM3U" || name_upper == "TPM3UF") {
		name = "TPM3u";
		rate_type = "120120";
		def_freq = FREQ_ESTIMATE;
		full_name = "TPM3 unequal frequencies ()";
	} else if (name_upper == "TIM" || name_upper == "TIM1") {
		name = "TIM";
		rate_type = "012230";
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM ()";
	} else if (name_upper == "TIMEF" || name_upper == "TIME" || name_upper == "TIM1EF" || name_upper == "TIM1E") {
		name = "TIMe";
		rate_type = "012230";
		def_freq = FREQ_EQUAL;
		full_name = "TIM equal frequencies";
	} else if (name_upper == "TIM2") {
		name = "TIM2";
		rate_type = "121030";
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM2 ()";
	} else if (name_upper == "TIM2EF" || name_upper == "TIM2E") {
		name = "TIM2e";
		rate_type = "121030";
		def_freq = FREQ_EQUAL;
		full_name = "TIM2 equal frequencies";
	} else if (name_upper == "TIM3") {
		name = "TIM3";
		rate_type = "120130";
		def_freq = FREQ_ESTIMATE;
		full_name = "TIM3 ()";
	} else if (name_upper == "TIM3EF" || name_upper == "TIM3E") {
		name = "TIM3e";
		rate_type = "120130";
		def_freq = FREQ_EQUAL;
		full_name = "TIM3 equal frequencies";
	} else if (name_upper == "TVM") {
		name = "TVM";
		rate_type = "412310";
		def_freq = FREQ_ESTIMATE;
		full_name = "TVM";
	} else if (name_upper == "TVMEF" || name_upper == "TVME") {
		name = "TVMe";
		rate_type = "412310";
		def_freq = FREQ_EQUAL;
		full_name = "TVM equal frequencies";
	} else if (name_upper == "SYM") {
		name = "SYM";
		rate_type = "123450";
		def_freq = FREQ_EQUAL;
		full_name = "SYM (Zharkihk, 1994)";
	} else if (name_upper == "GTR" || name_upper == "REV") {
		name = "GTR";
		rate_type = "123450";
		def_freq = FREQ_ESTIMATE;
		full_name = "GTR (Tavare, 1986)";
	} else {
		name = "";
		rate_type = "";
		full_name = "";
	}
	return name;
}


void ModelDNA::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	assert(num_states == 4); // make sure that you create model for DNA
	StateFreqType def_freq = FREQ_UNKNOWN;
	string rate_type;
	name = getDNAModelInfo((string)model_name, full_name, rate_type, def_freq);

	if (name != "") {
		setRateType(rate_type.c_str());
	} else {
		//cout << "User-specified model "<< model_name << endl;
		if (setRateType(model_name))
			def_freq = FREQ_ESTIMATE;
		else {
			readParameters(model_name);
			//name += " (user-defined)";
		}
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


void ModelDNA::readRates(string str) throw(const char*) {
	int nrates = *max_element(param_spec.begin(), param_spec.end());
	int end_pos = 0;
	int i, j;
	for (j = 0; j < param_spec.length(); j++)
		rates[j] = 1.0;
	num_params = 0;
	for (i = 0; i < nrates && end_pos < str.length(); i++) {
		int new_end_pos;
		double rate = 0;
		if (str[end_pos] == '?') {
			param_fixed[i+1] = false;
			end_pos++;
			rate = i + 0.4;
			num_params++;
		} else {
			param_fixed[i+1] = true;
			try {
				rate = convert_double(str.substr(end_pos).c_str(), new_end_pos);
			} catch (string str) {
				outError(str);
			}
			end_pos += new_end_pos;
		}
		if (rate < 0.0)
			outError("Negative rates found");
		if (i == nrates-1 && end_pos < str.length())
			outError("String too long ", str);
		if (i < nrates-1 && end_pos >= str.length())
			outError("Unexpected end of string ", str);
		if (end_pos < str.length() && str[end_pos] != ',')
			outError("Comma to separate rates not found in ", str);
		end_pos++;
		for (j = 0; j < param_spec.length(); j++)
			if (param_spec[j] == i+1)
				rates[j] = rate;
	}
}


string ModelDNA::getNameParams() {
	if (num_params == 0) return name;
	ostringstream retname;
	retname << name << '{';
	int nrates = getNumRateEntries();
	int k = 0;
	for (int i = 0; i < nrates; i++) {
		if (param_spec[i] > k) {
			if (k>0) retname << ',';
			retname << rates[i];
			k++;
		}
	}
	retname << '}';
	return retname.str();
}

bool ModelDNA::setRateType(const char *rate_str) {
	//char first_type = 127;
	//char last_type = 0;
	//char t = first_type;
	int num_ch = strlen(rate_str);
	int i;

	if (num_ch != getNumRateEntries()) {
		//outError("Model specification has wrong length!");
		return false;
	}
	// only accept string of digits
	for (i = 0; i < num_ch; i++)
		if (!isdigit(rate_str[i])) return false;
	/*
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
	}*/

	map<char,char> param_k;
	num_params = 0;
	param_spec = "";
	// last entry get ID of 0 for easy management
	param_k[rate_str[num_ch-1]] = 0;
	for (i = 0; i < num_ch; i++) {
		if (param_k.find(rate_str[i]) == param_k.end()) {
			num_params++;
			param_k[rate_str[i]] = (char)num_params;
			param_spec.push_back(num_params);
		} else {
			param_spec.push_back(param_k[rate_str[i]]);
		}
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
	param_fixed.resize(num_params+1, false);
	param_fixed[0] = true; // fix the last entry
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
	if (num_params > 0) {
		int num_all = param_spec.length();
		if (verbose_mode >= VB_MAX) {
			for (i = 1; i <= num_params; i++)
				cout << "  estimated variables[" << i << "] = " << variables[i] << endl;
		}
		for (i = 0; i < num_all; i++)
			if (!param_fixed[param_spec[i]]) {
				rates[i] = variables[(int)param_spec[i]];
			}
	}
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
	if (num_params > 0) {
		int num_all = param_spec.length();
		for (int i = 0; i < num_all; i++)
			if (!param_fixed[param_spec[i]])
				variables[(int)param_spec[i]] = rates[i];
	}
	if (freq_type == FREQ_ESTIMATE) {
		int ndim = getNDim();
		memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
	}
}
