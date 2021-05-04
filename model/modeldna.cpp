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
#include "modelliemarkov.h"
#include <utils/stringfunctions.h> //for string_to_upper, convert_double

ModelDNA::ModelDNA(PhyloTree *tree, PhyloTree* report_to_tree)
: ModelMarkov(tree, report_to_tree)
{
}

ModelDNA::ModelDNA(const char* model_name, std::string model_params,
                   StateFreqType freq,     std::string freq_params,
                   PhyloTree*    tree,     PhyloTree*  report_to_tree)
    : ModelMarkov(tree, report_to_tree)
{
  init(model_name, model_params, freq, freq_params, report_to_tree);
}

namespace {
    const struct dna_model_alias {
        public: const char* from, *to;
    } dna_model_aliases [] =
    {
        { "JC69",   "JC" },    { "K80",    "K2P"   }, { "HKY85",  "HKY"   },
        { "K81",    "K3P" },   { "TPM1",   "K3P"   },
        { "K81UF",  "K3PU" },  { "K81U",   "K3PU"  }, { "K3PUF",  "K3PU"  },
        { "TPM1UF", "K3PU"  }, { "TPM1U",  "K3PU"  },
        { "TRN",    "TN" },    { "TN93",   "TN"    },
        { "TNEF",   "TNE" },   { "TRNEF",  "TNE"   }, { "TRNE",   "TNE"   },
        { "TPM2UF", "TPM2U" }, { "TPM3UF", "TPM3U" }, { "TIM1",   "TIM"   },
        { "TIMEF",  "TIME" },  { "TIM1EF", "TIME"  }, { "TIM1E",  "TIME"  },
        { "TIM2EF", "TIM2E" }, { "TIM3EF", "TIM3E" },
        { "TVMEF",  "TVME" },  { "REV",    "GTR"   },
    };
    const struct freq_lookup  {
    public: const char* name;
            const char* rate_type;
            StateFreqType def_freq;
            const char* full_name;
    } dna_model_lookups [] =
    {
        { "JC",    "000000", StateFreqType::FREQ_EQUAL,    "JC (Juke and Cantor, 1969)"},
        { "F81",   "000000", StateFreqType::FREQ_ESTIMATE, "F81 (Felsenstein, 1981)"},
        { "K2P",   "010010", StateFreqType::FREQ_EQUAL,    "K2P (Kimura, 1980)"},
        { "HKY",   "010010", StateFreqType::FREQ_ESTIMATE, "HKY (Hasegawa, Kishino and Yano, 1985)"},
        { "K3P",   "012210", StateFreqType::FREQ_EQUAL,    "K3P (Kimura, 1981)" },
        { "K3Pu",  "012210", StateFreqType::FREQ_ESTIMATE, "K3P unequal frequencies (Kimura, 1981)"},
        { "TN",    "010020", StateFreqType::FREQ_ESTIMATE, "TN (Tamura and Nei, 1993)"},
        { "TNe",   "010020", StateFreqType::FREQ_EQUAL,    "TN equal frequencies (Tamura and Nei, 1993)"},
        { "TPM2",  "121020", StateFreqType::FREQ_ESTIMATE, "TPM2 ()"},
        { "TPM2u", "121020", StateFreqType::FREQ_ESTIMATE, "TPM2 unequal frequencies ()"},
        { "TPM3",  "120120", StateFreqType::FREQ_ESTIMATE, "TPM3 ()"},
        { "TPM3u", "120120", StateFreqType::FREQ_ESTIMATE, "TPM3 unequal frequencies ()"},
        { "TIM",   "012230", StateFreqType::FREQ_ESTIMATE, "TIM ()"},
        { "TIMe",  "012230", StateFreqType::FREQ_EQUAL,    "TIM equal frequencies"},
        { "TIM2",  "121030", StateFreqType::FREQ_ESTIMATE, "TIM2 ()"},
        { "TIM2e", "121030", StateFreqType::FREQ_EQUAL,    "TIM2 equal frequencies"},
        { "TIM3",  "120130", StateFreqType::FREQ_ESTIMATE, "TIM3 ()"},
        { "TIM3e", "120130", StateFreqType::FREQ_EQUAL,    "TIM3 equal frequencies"},
        { "TVM",   "412310", StateFreqType::FREQ_ESTIMATE, "TVM"},
        { "TVMe",  "412310", StateFreqType::FREQ_EQUAL,    "TVM equal frequencies"},
        { "SYM",   "123450", StateFreqType::FREQ_EQUAL,    "SYM (Zharkihk, 1994)"},
        { "GTR",   "123450", StateFreqType::FREQ_ESTIMATE, "GTR (Tavare, 1986)"},
    };
};

string getDNAModelInfo(string model_name, string &full_name,
                       string &rate_type, StateFreqType &def_freq) {
    string name_upper = string_to_upper(model_name);
    string name       = model_name;
    full_name         = name;
    rate_type         = "";
    def_freq          = StateFreqType::FREQ_UNKNOWN;
    
    std::string search = name_upper;
    for (int i=0; i<sizeof(dna_model_aliases)/sizeof(dna_model_aliases[0]); ++i) {
        if (dna_model_aliases[i].from == search) {
            search = dna_model_aliases[i].to;
            break;
        }
    }
    name      = "";
    rate_type = "";
    full_name = "";
    for (int i=0; i<sizeof(dna_model_lookups)/sizeof(dna_model_lookups[0]); ++i) {
        std::string this_one = string_to_upper(dna_model_lookups[i].name);
        if (this_one == search) {
            name      = dna_model_lookups[i].name;
            rate_type = dna_model_lookups[i].rate_type;
            def_freq  = dna_model_lookups[i].def_freq;
            full_name = dna_model_lookups[i].full_name;
            break;
        }
    }
    return name;
}

void ModelDNA::init(const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree* report_to_tree)
{
    ASSERT(num_states == 4); // make sure that you create model for DNA
    StateFreqType def_freq = StateFreqType::FREQ_UNKNOWN;
    string rate_type;
    // First try: the time reversible models
    name = getDNAModelInfo((string)model_name, full_name, rate_type, def_freq);
    if (name == "") {
        // Second try: Lie Markov models. (Note, we're still missing UNREST
        // model. 12.12 is equivalent, but user may not realize that.)
        int model_num, symmetry; // returned by getLieMarkovModelInfo, but not used here
        ModelLieMarkov::getLieMarkovModelInfo((string)model_name, name, full_name,
                                              model_num, symmetry, def_freq);
    }
    if (name != "") {
        setRateType(rate_type.c_str());
    } else {
        //cout << "User-specified model "<< model_name << endl;
        if (setRateType(model_name)) {
            // model was six digits (e.g. 010010 for K2P/HKY)
            name = model_name;
            full_name = "Time reversible ("+name+")";
        } else if (strlen(model_name)!=0) {
            readParameters(model_name, true, report_to_tree);
            name = full_name = model_name;
            freq = StateFreqType::FREQ_USER_DEFINED;
            //name += " (user-defined)";
        }
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
    ModelMarkov::init(freq, report_to_tree);
    //model_parameters = new double [getNDim()+1]; 
    // see setVariables for explaination of +1
    //setVariables(model_parameters);
}

void ModelDNA::startCheckpoint() {
    checkpoint->startStruct("ModelDNA");
}

void ModelDNA::saveCheckpoint() {
    // construct model_parameters from rates and base freqs. 
    // This is one-indexed, so parameters are in model_parameters[1]
    // up to model_parameters[num_params]
    //    setVariables(model_parameters);
    startCheckpoint();
    if (!fixed_parameters) {
        CKP_ARRAY_SAVE(6, rates);
    }
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelDNA::restoreCheckpoint() {
    // curiously, this seems to be the only plase ModelDNA uses model_parameters.
    ModelMarkov::restoreCheckpoint();
    startCheckpoint();
    if (!fixed_parameters) {
        CKP_ARRAY_RESTORE(6, rates);
    }
    endCheckpoint();
    //getVariables(model_parameters);       // updates rates and state_freq
    string rate_spec = param_spec;
    for (auto i = rate_spec.begin(); i != rate_spec.end(); ++i) {
        *i = *i + '0';
    }
    if (!rate_spec.empty()) {
        if (!setRateType(rate_spec)) {
            ASSERT(0 && "Cannot set rate_spec");
        }
    }
    decomposeRateMatrix();
    if (phylo_tree) {
        phylo_tree->clearAllPartialLH();
    }
}

int ModelDNA::getNumberOfRates() const {
    return static_cast<int>(param_spec.length());
}

void ModelDNA::readRates(string str) THROW_SPEC(const char*) {
    int nrates = *max_element(param_spec.begin(), param_spec.end());
    int end_pos = 0;
    int i, j;
    for (j = 0; j < param_spec.length(); ++j) {
        rates[j] = 1.0;
    }
    num_params = 0;
    for (i = 0; i <= nrates && end_pos < str.length(); ++i) {
        int new_end_pos;
        double rate = 0;
        int id = (i < nrates) ? i+1 : 0;
        if (str[end_pos] == '?') {
            param_fixed[id] = false;
            end_pos++;
            rate = 1.0;
            num_params++;
        } else {
            if (Params::getInstance().optimize_rate_matrix) {
                num_params++;
                param_fixed[id] = false;
            } else
                if (Params::getInstance().optimize_from_given_params) {
                    num_params++;
                    param_fixed[id] = false;
                } else {
                    param_fixed[id] = true;
                }
            try {
                rate = convert_double(str.substr(end_pos).c_str(), new_end_pos);
            } catch (string str) {
                outError(str);
            }
            end_pos += new_end_pos;
        }
        if (rate < 0.0) {
            outError("Negative rates found");
        }
        if (i == nrates && end_pos < str.length()) {
            outError("More than " + convertIntToString(nrates) + 
            " rate parameters specified in " + str);
        }
        if (i < nrates-1 && end_pos >= str.length()) {
            outError("Unexpected end of string ", str);
        }
        if (end_pos < str.length() && str[end_pos] != ',') {
            outError("Comma to separate rates not found in ", str);
        }
        end_pos++;
        for (j = 0; j < param_spec.length(); ++j) {
            if (param_spec[j] == id) {
                rates[j] = rate;
            }
        }
    }
}

string ModelDNA::getNameParams() {
    if (num_params == 0) {
        return name;
    }
    ostringstream retname;
    retname << name;
    if (!fixed_parameters) {
        retname << '{';
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
    }
    getNameParamsFreq(retname);
    return retname.str();
}

bool ModelDNA::setRateType(string rate_str) {
	//char first_type = 127;
	//char last_type = 0;
	//char t = first_type;
	int num_ch = static_cast<int>(rate_str.length());
	int i;

	if (num_ch != getNumRateEntries()) {
		//outError("Model specification has wrong length!");
		return false;
	}
	// only accept string of digits
    for (i = 0; i < num_ch; ++i) {
        if (!isdigit(rate_str[i])) {
            return false;
        }
    }
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
    for (i = 0; i < num_ch; ++i) {
        if (param_k.find(rate_str[i]) == param_k.end()) {
            num_params++;
            param_k[rate_str[i]] = (char)num_params;
            param_spec.push_back(num_params);
        } else {
            param_spec.push_back(param_k[rate_str[i]]);
        }
    }

	ASSERT(param_spec.length() == num_ch);
	double *avg_rates = new double[num_params+1];
	int *num_rates = new int[num_params+1];
	memset(avg_rates, 0, sizeof(double) * (num_params+1));
	memset(num_rates, 0, sizeof(int) * (num_params+1));
	for (i = 0; i < param_spec.size(); ++i) {
		avg_rates[(int)param_spec[i]] += rates[i];
		num_rates[(int)param_spec[i]]++;
	}
    for (i = 0; i <= num_params; ++i) {
        avg_rates[i] /= num_rates[i];
    }
	for (i = 0; i < param_spec.size(); ++i) {
        if (avg_rates[0] > 0.0) {
            rates[i] = avg_rates[(int)param_spec[i]] / avg_rates[0];
        }
        else {
            rates[i] = avg_rates[(int)param_spec[i]];
        }
	}
	if (verbose_mode >= VerboseMode::VB_DEBUG) {
		cout << "Initialized rates: ";
        for (i = 0; i < param_spec.size(); ++i) {
			cout << rates[i] << " ";
        }
		cout << endl;
	}
    if (param_fixed.size() == num_params + 1) {
        num_params = 0;
        for (auto p : param_fixed) {
            if (!p) {
                num_params++;
            }
        }
    } else {
        param_fixed.resize(num_params+1, false);
        param_fixed[0] = true; // fix the last entry
    }
    delete [] num_rates;
    delete [] avg_rates;
    return true;
}

int ModelDNA::getNDim() {
    if (fixed_parameters)
    {
        return 0;
    }
    ASSERT(freq_type != StateFreqType::FREQ_UNKNOWN);
    // possible TO-DO: cache nFreqParams(freq_type) to avoid repeat calls.
    // return (num_params+nFreqParams(freq_type));

    // if (linked_model && linked_model != this) {
    //     return 0;
    // }
    
    int ndim = num_params;
    if (freq_type == StateFreqType::FREQ_ESTIMATE) {
        ndim += num_states - 1;
    } else {
        ndim += nFreqParams(freq_type);
    }
    return ndim;
}

void ModelDNA::writeParameters(ostream& out) {
    int i;
    if (freq_type == StateFreqType::FREQ_ESTIMATE) {
        for (i = 0; i < num_states; ++i) {
            out << "\t" << state_freq[i];
        }
    }
    if (num_params == 0) {
        return;
    }
    if (num_params <= 1) {
        out << "\t" << rates[1];
    }
    else {
        int nrateout = getNumRateEntries() - 1;
        for (i = 0; i < nrateout; ++i)
            out << "\t" << rates[i];
    }
}

void ModelDNA::computeTipLikelihood(PML::StateType state, double *state_lk) {
    if (static_cast<int>(state) < num_states || state >= 18) {
        ModelSubst::computeTipLikelihood(state, state_lk);
        return;
    }
    // special treatment for ambiguous (polymorphic) state
    memset(state_lk, 0, num_states*sizeof(double));
    int cstate = state-num_states+1;
    for (int i = 0; i < num_states; i++) {
        if ((cstate) & (1 << i)) {
            state_lk[i] = 1.0;
        }
    }
}

/*
 * getVariables *changes* the state of the model, setting from *variables
 * Returns true if the model state has changed, false if not.
 */
bool ModelDNA::getVariables(double *variables) {
    bool changed = false;
    if (num_params > 0) {
        int num_all = static_cast<int>(param_spec.length());
        if (verbose_mode >= VerboseMode::VB_MAX) {
            for (int i = 1; i <= num_params; i++) {
                cout << "  estimated variables[" << i << "] = "
                     << variables[i] << endl;
            }
        }
        for (int i = 0; i < num_all; i++) {
            if (!param_fixed[param_spec[i]]) {
                changed |= (rates[i] != variables[(int)param_spec[i]]);
                            rates[i]  = variables[(int)param_spec[i]];
            }
        }
    }
    if (freq_type == StateFreqType::FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1,
        // this will be done at the end of optimization
        int ndim = getNDim();
        changed |= memcmpcpy(state_freq, variables+(ndim-num_states+2),
                             (num_states-1)*sizeof(double));
        //double sum = 0;
        //for (int i = 0; i < num_states-1; i++)
        //    sum += state_freq[i];
        //state_freq[num_states-1] = 1.0 - sum;
    } else {
        // BQM: for special DNA freq stuffs from MDW
        changed |= freqsFromParams(state_freq,variables+num_params+1,freq_type);
    }
    return changed;

        // BUG FIX 2015.08.28
//        int nrate = getNDim();
//        if (freq_type == StateFreqType::FREQ_ESTIMATE) nrate -= (num_states-1);
//              double sum = 1.0;
//              int i, j;
//              for (i = 1; i < num_states; i++)
//                      sum += variables[nrate+i];
//              for (i = 0, j = 1; i < num_states; i++)
//                      if (i != highest_freq_state) {
//                              state_freq[i] = variables[nrate+j] / sum;
//                              j++;
//                      }
//              state_freq[highest_freq_state] = 1.0/sum;
}

/*
 * setVariables *reads* the state of the model and writes into "variables"
 * Model does not change state. *variables should have length getNDim()+1
 * If param_spec is (e.g.) 012210 (e.g. K3P model) then in general
 * we'd have rates 0 and 5 (A<->C and G<->T) written to variables[0],
 * rates 1 and 4 to variables[1] and rates 2 and 3 to variables[2].
 * However one of these (typically 0) is 'fixed' (param_fixed)
 * to always have value 1, and this doesn't get written. 
 * num_parameters in this case will be two, for two free rates parameters.
 * Base frequency parameters get written after the rate parameters, so
 * K3P+FO model (i.e. fits base frequencies with no constraints) would
 * use variables[3] to variables[5] (3 values) to store base freq info. 
 */
void ModelDNA::setVariables(double *variables) {
    if (num_params > 0) {
        int num_all = static_cast<int>(param_spec.length());
        for (int i = 0; i < num_all; i++) {
            if (!param_fixed[param_spec[i]]) {
                variables[(int)param_spec[i]] = rates[i];
            }
        }
    }
    // and copy parameters for base frequencies

    if (freq_type == StateFreqType::FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1,
        // this will be done at the end of optimization
        int ndim = getNDim();
        memcpy(variables+(ndim-num_states+2), state_freq,
               (num_states-1)*sizeof(double));
    } else {
        paramsFromFreqs(variables+num_params+1, state_freq, freq_type);
    }

        // BUG FIX 2015.08.28
//        int nrate = getNDim();
//        if (freq_type == StateFreqType::FREQ_ESTIMATE) nrate -= (num_states-1);
//              int i, j;
//              for (i = 0, j = 1; i < num_states; i++)
//                      if (i != highest_freq_state) {
//                              variables[nrate+j] = state_freq[i] / state_freq[highest_freq_state];
//                              j++;
//                      }
}
