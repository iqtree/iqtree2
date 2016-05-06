/*
 * phylotesting.cpp
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "phylotree.h"
#include "iqtree.h"
#include "phylosupertree.h"
#include "phylotesting.h"

#include "model/modelgtr.h"
#include "model/modeldna.h"
#include "myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
#include "model/ratefree.h"
//#include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "model/modelmorphology.h"
#include "model/modelmixture.h"
#include "timeutil.h"

#include "phyloanalysis.h"
#include "gsl/mygsl.h"
#include "vectorclass/vectorclass.h"


/******* Binary model set ******/
const char* bin_model_names[] = { "JC2", "GTR2" };


/******* Morphological model set ******/
const char* morph_model_names[] = {"MK", "ORDERED"};


/******* DNA model set ******/
const char* dna_model_names[] = { "JC", "F81", "K80", "HKY", "TNe",
		"TN", "K81", "K81u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIMe", "TIM",
		"TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR" };

/* DNA models supported by PhyML/PartitionFinder */
const char* dna_model_names_old[] ={"JC", "F81", "K80", "HKY", "TNe",
 	 	 "TN", "K81", "K81u", "TIMe", "TIM", "TVMe", "TVM", "SYM", "GTR"};

/* DNA model supported by RAxML */
const char* dna_model_names_rax[] ={"GTR"};

/* DNA model supported by MrBayes */
const char *dna_model_names_mrbayes[] = {"JC", "F81", "K80", "HKY", "SYM", "GTR"};


/****** Protein model set ******/
const char* aa_model_names[] = { "Dayhoff", "mtMAM", "JTT", "WAG",
		"cpREV", "mtREV", "rtREV", "mtART", "mtZOA", "VT", "LG", "DCMut", "PMB",
		"HIVb", "HIVw", "JTTDCMut", "FLU", "Blosum62" };
        
/* Protein models supported by PhyML/PartitionFinder */
const char *aa_model_names_phyml[] = { "Dayhoff", "mtMAM", "JTT", "WAG",
		"cpREV", "mtREV", "rtREV", "mtART", "VT", "LG", "DCMut",
		"HIVb", "HIVw", "Blosum62" };

/* Protein models supported by RAxML */
const char *aa_model_names_rax[] = { "Dayhoff", "mtMAM", "JTT", "WAG",
		"cpREV", "mtREV", "rtREV", "mtART", "mtZOA", "PMB", "HIVb", "HIVw", "JTTDCMut", "FLU", "VT", "LG", "DCMut", "Blosum62" };

const char* aa_model_names_mrbayes[] = {"Poisson", "Dayhoff", "mtMAM", "JTT", "WAG",
		"cpREV", "mtREV", "rtREV", "VT", "Blosum62" };

const char *aa_model_names_nuclear[] = {"WAG", "Dayhoff","JTT", "LG", "VT", "DCMut", "PMB", "JTTDCMut", "Blosum62"};

const char *aa_model_names_mitochondrial[] = {"mtREV", "mtMAM", "mtART", "mtZOA"};

const char *aa_model_names_chloroplast[] = {"cpREV"};

const char *aa_model_names_viral[] = {"HIVb", "HIVw", "FLU", "rtREV"};

const char* aa_freq_names[] = {"", "+F"};


/****** Codon models ******/
//const char *codon_model_names[] = {"GY", "MG", "MGK", "KOSI07", "SCHN05","KOSI07_GY1KTV","SCHN05_GY1KTV"};
//short int std_genetic_code[]    = {   0,    0,     0,        1,        1,              1,              1};
const char *codon_model_names[] = {"MG", "MGK", "GY", "KOSI07", "SCHN05"};
short int std_genetic_code[]    = {   0,    0,     0,        1,        1};

const char *codon_freq_names[] = {"", "+F1X4", "+F3X4", "+F"};

const double TOL_LIKELIHOOD_MODELTEST = 0.01;
const double TOL_GRADIENT_MODELTEST   = 0.001;

/**
 * copy from cvec to strvec
 */
void copyCString(const char **cvec, int n, StrVector &strvec, bool touppercase = false) {
	strvec.resize(n);
	for (int i = 0; i < n; i++) {
		strvec[i] = cvec[i];
        if (touppercase)
            std::transform(strvec[i].begin(), strvec[i].end(), strvec[i].begin(), ::toupper);
    }
}

int getSeqType(const char *model_name, SeqType &seq_type) {
    bool empirical_model = false;
    int i;
    string model_str = model_name;
    std::transform(model_str.begin(), model_str.end(), model_str.begin(), ::toupper);
    StrVector model_list;

    seq_type = SEQ_UNKNOWN;
    
    copyCString(bin_model_names, sizeof(bin_model_names)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_BINARY;
            break;
        }
    copyCString(morph_model_names, sizeof(morph_model_names)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_MORPH;
            break;
        }
    copyCString(dna_model_names, sizeof(dna_model_names)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_DNA;
            break;
        }
    copyCString(aa_model_names, sizeof(aa_model_names)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_PROTEIN;
            empirical_model = true;
            break;
        }
    copyCString(codon_model_names, sizeof(codon_model_names)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str.substr(0,model_list[i].length()) == model_list[i]) {
            seq_type = SEQ_CODON;
            if (std_genetic_code[i]) empirical_model = true;
            break;
        }
        
    return (empirical_model) ? 2 : 1;
}

string getSeqType(string model_name) {
    SeqType seq_type;
    getSeqType(model_name.c_str(), seq_type);
    switch (seq_type) {
    case SEQ_BINARY: return "BIN"; break;
    case SEQ_MORPH: return "MORPH"; break;
    case SEQ_DNA: return "DNA"; break;
    case SEQ_PROTEIN: return "AA"; break;
    case SEQ_CODON: return "CODON"; break;
    default: break;
    }
    return "";
}

void computeInformationScores(double tree_lh, int df, int ssize, double &AIC, double &AICc, double &BIC) {
	AIC = -2 * tree_lh + 2 * df;
	AICc = AIC + 2.0 * df * (df + 1) / max(ssize - df - 1, 1);
	BIC = -2 * tree_lh + df * log(ssize);
}

double computeInformationScore(double tree_lh, int df, int ssize, ModelTestCriterion mtc) {
	double AIC, AICc, BIC;
	computeInformationScores(tree_lh, df, ssize, AIC, AICc, BIC);
	if (mtc == MTC_AIC)
		return AIC;
	if (mtc == MTC_AICC)
		return AICc;
	if (mtc == MTC_BIC)
		return BIC;
	return 0.0;
}

string criterionName(ModelTestCriterion mtc) {
	if (mtc == MTC_AIC)
		return "AIC";
	if (mtc == MTC_AICC)
		return "AICc";
	if (mtc == MTC_BIC)
		return "BIC";
	return "";
}

void printSiteLh(const char*filename, PhyloTree *tree, double *ptn_lh,
		bool append, const char *linename) {
	int i;
	double *pattern_lh;
	if (!ptn_lh) {
		pattern_lh = new double[tree->getAlnNPattern()];
		tree->computePatternLikelihood(pattern_lh);
	} else
		pattern_lh = ptn_lh;

	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		if (append) {
			out.open(filename, ios::out | ios::app);
		} else {
			out.open(filename);
			out << 1 << " " << tree->getAlnNSite() << endl;
		}
		IntVector pattern_index;
		tree->aln->getSitePatternIndex(pattern_index);
		if (!linename)
			out << "Site_Lh   ";
		else {
			out.width(10);
			out << left << linename;
		}
		for (i = 0; i < tree->getAlnNSite(); i++)
			out << " " << pattern_lh[pattern_index[i]];
		out << endl;
		out.close();
		if (!append)
			cout << "Site log-likelihoods printed to " << filename << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, filename);
	}

	if (!ptn_lh)
		delete[] pattern_lh;
}

void printSiteLhCategory(const char*filename, PhyloTree *tree, SiteLoglType wsl) {

    if (wsl == WSL_NONE || wsl == WSL_SITE)
        return;
    // error checking
    if (!tree->getModel()->isMixture()) {
        if (wsl != WSL_RATECAT) {
            outWarning("Switch now to '-wslr' as it is the only option for non-mixture model");
            wsl = WSL_RATECAT;
        }
    } else {
        // mixture model
        if (wsl == WSL_MIXTURE_RATECAT && tree->getModelFactory()->fused_mix_rate) {
            outWarning("-wslmr is not suitable for fused mixture model, switch now to -wslm");
            wsl = WSL_MIXTURE;
        }
    }
	int ncat = tree->getNumLhCat(wsl);
	double *pattern_lh, *pattern_lh_cat;
	int i;
	pattern_lh = new double[tree->getAlnNPattern()];
	pattern_lh_cat = new double[tree->getAlnNPattern()*ncat];
	tree->computePatternLikelihood(pattern_lh, NULL, pattern_lh_cat, wsl);

    
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(filename);
		out << "Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)" << endl;
        if (wsl == WSL_RATECAT) {
            out << "P(D|M,rr[i]) is the probability of site D given the model M and the relative rate" << endl;
            out << "of evolution rr[i], where i is the class of rate to be considered." << endl;
            out << "We have P(D|M) = \\sum_i P(i) x P(D|M,rr[i])." << endl << endl;
            out << "Site   logP(D|M)       ";
            for (i = 0; i < ncat; i++)
                out << "log{P(" << i+1 << ")xP(D|M,rr[" << i+1 << "]=" << tree->getRate()->getRate(i)<< ")} ";
        } else if (wsl == WSL_MIXTURE) {
            out << "P(D|M[i]) is the probability of site D given the model M[i]," << endl;
            out << "where i is the mixture class to be considered." << endl;
            out << "We have P(D|M) = \\sum_i P(i) x P(D|M[i])." << endl << endl;
            out << "Site   logP(D|M)       ";
            for (i = 0; i < ncat; i++)
                out << "log{P(" << i+1 << ")xP(D|M[" << i+1 << "])} ";
        } else {
            // WSL_MIXTURE_RATECAT
            out << "P(D|M[i],rr[j]) is the probability of site D given the model M[i] and the relative rate" << endl;
            out << "of evolution rr[j], where i and j are the mixture class and rate class, respectively." << endl;
            out << "We have P(D|M) = \\sum_i \\sum_j P(i) x P(j) x P(D|M[i],rr[j])." << endl << endl;
            out << "Site   logP(D|M)       ";
            for (i = 0; i < tree->getModel()->getNMixtures(); i++)
                for (int j = 0; j < tree->getRate()->getNRate(); j++) {
                    out << "log{P(" << i+1 << ")xP(" << j+1 << ")xP(D|M[" << i+1 << "],rr[" << j+1 << "]=" << tree->getRate()->getRate(j) << ")} ";
                }
        }
		out << endl;
		IntVector pattern_index;
		tree->aln->getSitePatternIndex(pattern_index);
		for (i = 0; i < tree->getAlnNSite(); i++) {
			out.width(6);
			out << left << i+1 << " ";
			out.width(15);
			out << pattern_lh[pattern_index[i]] << " ";
			for (int j = 0; j < ncat; j++) {
				out.width(15);
				out << pattern_lh_cat[pattern_index[i]*ncat+j] << " ";
			}
			out << endl;
		}
		out.close();
		cout << "Site log-likelihoods per category printed to " << filename << endl;
        if (!tree->isSuperTree()) {
            cout << "Log-likelihood of constant sites: " << endl;
            double const_prob = 0.0;
            for (i = 0; i < tree->aln->getNPattern(); i++)
                if (tree->aln->at(i).is_const) {
                    Pattern pat = tree->aln->at(i);
                    for (Pattern::iterator it = pat.begin(); it != pat.end(); it++)
                        cout << tree->aln->convertStateBackStr(*it);
                    cout << ": " << pattern_lh[i] << endl;
                    const_prob += exp(pattern_lh[i]);
                }
            cout << "Probability of const sites: " << const_prob << endl;
        }
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, filename);
	}

	delete[] pattern_lh_cat;
	delete[] pattern_lh;

}

void printSiteStateFreq(const char*filename, PhyloTree *tree) {

    int i, j, nsites = tree->getAlnNSite(), nstates = tree->aln->num_states;
    double *ptn_state_freq = new double[tree->getAlnNPattern() * nstates];
    
    tree->computePatternStateFreq(ptn_state_freq);

	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(filename);
		IntVector pattern_index;
		tree->aln->getSitePatternIndex(pattern_index);
		for (i = 0; i < nsites; i++) {
			out.width(6);
			out << left << i+1 << " ";
            double *state_freq = &ptn_state_freq[pattern_index[i]*nstates];
			for (j = 0; j < nstates; j++) {
				out.width(15);
				out << state_freq[j] << " ";
			}
			out << endl;
		}
		out.close();
		cout << "Site state frequency vectors printed to " << filename << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, filename);
	}
    delete [] ptn_state_freq;
}

bool checkModelFile(ifstream &in, bool is_partitioned, vector<ModelInfo> &infos) {
	if (!in.is_open()) return false;
	in.exceptions(ios::badbit);
	string str;
	if (is_partitioned) {
		in >> str;
		if (str != "Charset")
			return false;
	}
	in >> str;
	if (str != "Model")
		return false;
	in >> str;
	if (str != "df")
		return false;
	in >> str;
	if (str != "LnL")
		return false;
	in >> str;
	if (str != "TreeLen") {
        outWarning(".model file was produced from a previous version of IQ-TREE");
		return false;
    }
	getline(in, str);
	while (!in.eof()) {
		in >> str;
		if (in.eof())
			break;
		ModelInfo info;
		if (is_partitioned) {
			info.set_name = str;
			in >> str;
		}
		info.name = str;
		in >> info.df >> info.logl >> info.tree_len;
		getline(in, str);
        info.tree = "";
        if (*str.rbegin() == ';') {
            size_t pos = str.rfind('\t');
            if (pos != string::npos)
                info.tree = str.substr(pos+1);
//            else 
//                outWarning(".model file was produced from a previous version of IQ-TREE");
        }
		infos.push_back(info);
		//cout << str << " " << df << " " << logl << endl;
	}
	in.clear();
	return true;
}

bool checkModelFile(string model_file, bool is_partitioned, vector<ModelInfo> &infos) {
	if (!fileExists(model_file))
		return false;
	//cout << model_file << " exists, checking this file" << endl;
	ifstream in;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(model_file.c_str());
		if (!checkModelFile(in, is_partitioned, infos))
			throw false;
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (bool ret) {
		in.close();
		return ret;
	} catch (ios::failure) {
		outError("Cannot read file ", model_file);
	}
	return true;
}

/**
 * get the list of model
 * @param models (OUT) vectors of model names
 * @return maximum number of rate categories
 */
int getModelList(Params &params, Alignment *aln, StrVector &models, bool separate_rate = false) {
	StrVector model_names;
    StrVector freq_names;
	SeqType seq_type = aln->seq_type;
    
	const char *rate_options[]    = {  "", "+I", "+ASC", "+G", "+I+G", "+ASC+G", "+R", "+ASC+R"};
	bool test_options_default[]   = {true, true,  false, true,   true,    false,false,    false};
	bool test_options_morph[]     = {true,false,   true, true,  false,     true,false,    false};    
	bool test_options_noASC_I[]   = {true,false,  false, true,  false,    false,false,    false};    
	bool test_options_asc[]       ={false,false,   true,false,  false,     true,false,    false};
	bool test_options_new[]       = {true, true,  false, true,   true,    false, true,    false};
	bool test_options_morph_new[] = {true,false,   true, true,  false,     true, true,     true};
	bool test_options_noASC_I_new[] = {true,false,  false, true,  false,    false, true,    false};
	bool test_options_asc_new[]   ={false,false,   true,false,  false,     true,false,     true};
    bool *test_options = test_options_default;
//	bool test_options_codon[] =  {true,false,  false,false,  false,    false};
	const int noptions = sizeof(rate_options) / sizeof(char*);
	int i, j;
    
	if (seq_type == SEQ_BINARY) {
		copyCString(bin_model_names, sizeof(bin_model_names) / sizeof(char*), model_names);
	} else if (seq_type == SEQ_MORPH) {
		copyCString(morph_model_names, sizeof(morph_model_names) / sizeof(char*), model_names);
	} else if (seq_type == SEQ_DNA) {
		if (params.model_set == NULL) {
			copyCString(dna_model_names, sizeof(dna_model_names) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "partitionfinder") == 0 || strcmp(params.model_set, "phyml") == 0) {
			copyCString(dna_model_names_old, sizeof(dna_model_names_old) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "raxml") == 0) {
			copyCString(dna_model_names_rax, sizeof(dna_model_names_rax) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "mrbayes") == 0) {
			copyCString(dna_model_names_mrbayes, sizeof(dna_model_names_mrbayes) / sizeof(char*), model_names);
		} else {
			convert_string_vec(params.model_set, model_names);
		}
	} else if (seq_type == SEQ_PROTEIN) {
		if (params.model_set == NULL) {
			copyCString(aa_model_names, sizeof(aa_model_names) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "partitionfinder") == 0 || strcmp(params.model_set, "phyml") == 0) {
			copyCString(aa_model_names_phyml, sizeof(aa_model_names_phyml) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "raxml") == 0) {
			copyCString(aa_model_names_rax, sizeof(aa_model_names_rax) / sizeof(char*), model_names);
		} else if (strcmp(params.model_set, "mrbayes") == 0) {
			copyCString(aa_model_names_mrbayes, sizeof(aa_model_names_mrbayes) / sizeof(char*), model_names);
		} else {
			convert_string_vec(params.model_set, model_names);
		}
        copyCString(aa_freq_names, sizeof(aa_freq_names)/sizeof(char*), freq_names);
        
        if (params.model_subset) {
            StrVector submodel_names;
            if (strncmp(params.model_subset, "nuclear", 3) == 0) {
                copyCString(aa_model_names_nuclear, sizeof(aa_model_names_nuclear) / sizeof(char*), submodel_names);
            } else if (strncmp(params.model_subset, "mitochondrial", 3) == 0) {
                copyCString(aa_model_names_mitochondrial, sizeof(aa_model_names_mitochondrial) / sizeof(char*), submodel_names);
            } else if (strncmp(params.model_subset, "chloroplast", 3) == 0) {
                copyCString(aa_model_names_chloroplast, sizeof(aa_model_names_chloroplast) / sizeof(char*), submodel_names);
            } else if (strncmp(params.model_subset, "viral",3) == 0) {
                copyCString(aa_model_names_viral, sizeof(aa_model_names_viral) / sizeof(char*), submodel_names);
            } else {
                outError("Wrong -msub option");
            }
            for (i = 0; i < model_names.size(); i++) {
                bool appear = false;
                for (j = 0; j < submodel_names.size(); j++) 
                    if (model_names[i] == submodel_names[j]) {
                        appear = true;
                        break;
                    }
                if (!appear) {
                    model_names.erase(model_names.begin()+i);
                    i--;
                }
            }
        }

	} else if (seq_type == SEQ_CODON) {
		if (params.model_set == NULL) {
			if (aln->isStandardGeneticCode())
				copyCString(codon_model_names, sizeof(codon_model_names) / sizeof(char*), model_names);
			else {
                i = sizeof(codon_model_names) / sizeof(char*);
                for (j = 0; j < i; j++)
                    if (!std_genetic_code[j])
                        model_names.push_back(codon_model_names[j]);
//				copyCString(codon_model_names, sizeof(codon_model_names) / sizeof(char*) - 1, model_names);
            }
		} else
			convert_string_vec(params.model_set, model_names);
        copyCString(codon_freq_names, sizeof(codon_freq_names) / sizeof(char*), freq_names);
	}
    
	if (model_names.empty()) 
        return 1;
    
    if (params.state_freq_set)
        convert_string_vec(params.state_freq_set, freq_names);
    for (j = 0; j < freq_names.size(); j++) {
        std::transform(freq_names[j].begin(), freq_names[j].end(), freq_names[j].begin(), ::toupper);
//        for (i = 0; i < freq_names.size(); i++)
//            cout << " " << freq_names[i];
//        cout << endl;
        if (freq_names[j] != "" && freq_names[j][0] != '+')
            freq_names[j] = "+" + freq_names[j];
    }
    
    if (freq_names.size() > 0) {
        StrVector orig_model_names = model_names;
        model_names.clear();
        for (j = 0; j < orig_model_names.size(); j++) {
            if (aln->seq_type == SEQ_CODON) {
                SeqType seq_type;
                int model_type = getSeqType(orig_model_names[j].c_str(), seq_type);
                for (i = 0; i < freq_names.size(); i++) {
                    // disallow MG+F
                    if (freq_names[i] == "+F" && orig_model_names[j].find("MG") != string::npos)
                        continue;
                    if (freq_names[i] != "" || (model_type == 2 && orig_model_names[j].find("MG") == string::npos)) 
                        // empirical model also allow ""
                        model_names.push_back(orig_model_names[j] + freq_names[i]);
                }
            } else {
                for (i = 0; i < freq_names.size(); i++)
                    model_names.push_back(orig_model_names[j] + freq_names[i]);
            }
        }
    }

    bool with_new = params.model_name.find("NEW") != string::npos;
    bool with_asc = params.model_name.find("ASC") != string::npos;

//	if (seq_type == SEQ_CODON) {
//		for (i = 0; i < noptions; i++)
//			test_options[i] = test_options_codon[i];
//	} else 
    if (aln->frac_const_sites == 0.0) {
        // morphological or SNP data: activate +ASC
        if (with_new) {
            if (with_asc)
                test_options = test_options_asc_new;
            else if (seq_type == SEQ_PROTEIN)
                test_options = test_options_noASC_I_new;
            else
                test_options = test_options_morph_new;
        } else if (with_asc)
            test_options = test_options_asc;
        else if (seq_type == SEQ_PROTEIN)
            test_options = test_options_noASC_I;
        else
            test_options = test_options_morph;
	} else {
        // normal data, use +I instead
        if (with_new) {
            // change +I+G to +R
            if (with_asc)
                test_options = test_options_asc_new;
            else
                test_options = test_options_new;
        } else if (with_asc) {
            test_options = test_options_asc;
        } else
            test_options = test_options_default;
    }
    

    StrVector ratehet;
    int max_cats = params.num_rate_cats;

	if (params.ratehet_set) {
		// take the rate_options from user-specified models
		convert_string_vec(params.ratehet_set, ratehet);
		if (!ratehet.empty() && ratehet[0] == "default") {
			ratehet.erase(ratehet.begin());
			StrVector ratedef;
			for (j = 0; j < noptions; j++)
				if (test_options[j])
					ratedef.push_back(rate_options[j]);
			ratehet.insert(ratehet.begin(), ratedef.begin(), ratedef.end());
		}
        for (j = 0; j < ratehet.size(); j++) {
            if (ratehet[j] != "" && ratehet[j][0] != '+')
                ratehet[j] = "+" + ratehet[j];
            if (ratehet[j] == "+E") // for equal rate model 
                ratehet[j] = "";
        }
    } else {
        for (j = 0; j < noptions; j++)
            if (test_options[j])
                ratehet.push_back(rate_options[j]);
        
    }
    
    size_t pos;

    for (j = 0; j < ratehet.size(); j++)
        if ( (pos = ratehet[j].find("+R")) != string::npos && (pos >= ratehet[j].length()-2 || !isdigit(ratehet[j][pos+2]) )) {
            string str = ratehet[j];
            ratehet[j].insert(pos+2, convertIntToString(params.min_rate_cats));
            max_cats = max(max_cats, params.max_rate_cats);
            for (int k = params.min_rate_cats+1; k <= params.max_rate_cats; k++) {
                ratehet.insert(ratehet.begin()+j+k-params.min_rate_cats, str.substr(0, pos+2) + convertIntToString(k) + str.substr(pos+2));
            }
//            break;
        }

    if (separate_rate) {
        for (i = 0; i < model_names.size(); i++) 
            models.push_back(model_names[i]);
        for (j = 0; j < ratehet.size(); j++)
            if (ratehet[j] != "")
                models.push_back(ratehet[j]);
    } else {
        for (i = 0; i < model_names.size(); i++)
            for (j = 0; j < ratehet.size(); j++) {
                models.push_back(model_names[i] + ratehet[j]);
            }
    }
    if (params.model_extra_set) {
        StrVector extra_model_names;
        convert_string_vec(params.model_extra_set, extra_model_names);        
        models.insert(models.end(), extra_model_names.begin(), extra_model_names.end());
    }
    return max_cats;
}

/*
bool checkPartitionModel(Params &params, PhyloSuperTree *in_tree, vector<ModelInfo> &model_info) {
	return true;

	PhyloSuperTree::iterator it;
	int i, all_models = 0;
	for (it = in_tree->begin(), i = 0; it != in_tree->end(); it++, i++) {
		int count = 0;
		for (vector<ModelInfo>::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
			if (mit->set_name == in_tree->part_info[i].name)
				count++;
		int nstates = (*it)->aln->num_states;
		int num_models;
		getModelList(params, nstates, num_models);
		if (count != num_models * 4) {
			return false;
		}
		all_models += count;
	}
	return true;
}
*/
void replaceModelInfo(vector<ModelInfo> &model_info, vector<ModelInfo> &new_info) {
	vector<ModelInfo>::iterator first_info = model_info.end(), last_info = model_info.end();
	vector<ModelInfo>::iterator mit;
	// scan through models for this partition, assuming the information occurs consecutively
	for (mit = model_info.begin(); mit != model_info.end(); mit++)
		if (mit->set_name == new_info.front().set_name) {
			if (first_info == model_info.end()) first_info = mit;
		} else if (first_info != model_info.end()) {
			last_info = mit;
			break;
		}
	if (new_info.size() == (last_info - first_info)) {
		// replace sub vector
		for (mit = first_info; mit != last_info; mit++)
			*mit = new_info[mit - first_info];
	} else {
		if (first_info != model_info.end()) {
			model_info.erase(first_info, last_info);
		}
		model_info.insert(model_info.end(), new_info.begin(), new_info.end());
	}
}

void extractModelInfo(string set_name, vector<ModelInfo> &model_info, vector<ModelInfo> &part_model_info) {
	for (vector<ModelInfo>::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
		if (mit->set_name == set_name)
			part_model_info.push_back(*mit);
		else if (part_model_info.size() > 0)
			break;
}

void mergePartitions(PhyloSuperTree* super_tree, vector<IntVector> &gene_sets, StrVector &model_names) {
	cout << "Merging into " << gene_sets.size() << " partitions..." << endl;
	vector<IntVector>::iterator it;
	SuperAlignment *super_aln = (SuperAlignment*)super_tree->aln;
	vector<PartitionInfo> part_info;
	vector<PhyloTree*> tree_vec;
	for (it = gene_sets.begin(); it != gene_sets.end(); it++) {
		PartitionInfo info;
		info.name = "";
		info.position_spec = "";
		info.aln_file = "";
		info.sequence_type = "";
		info.model_name = model_names[it-gene_sets.begin()];
        info.part_rate = 1.0; // BIG FIX: make -spp works with -m TESTMERGE now!
        info.evalNNIs = 0;
		for (IntVector::iterator i = it->begin(); i != it->end(); i++) {
			if (i != it->begin()) {
				info.name += "+";
				info.position_spec += ", ";
			}
			info.name += super_tree->part_info[*i].name;
			info.position_spec += super_tree->part_info[*i].position_spec;
			if (!super_tree->part_info[*i].aln_file.empty()) {
                if (info.aln_file.empty())
                    info.aln_file = super_tree->part_info[*i].aln_file;
                else if (info.aln_file != super_tree->part_info[*i].aln_file) {
                    info.aln_file = "__NA__";
                }
			}
			if (!super_tree->part_info[*i].sequence_type.empty()) {
                if (info.sequence_type.empty())
                    info.sequence_type = super_tree->part_info[*i].sequence_type;
                else if (info.sequence_type != super_tree->part_info[*i].sequence_type) {
                    info.sequence_type = "__NA__";
                }
			}
		}
		info.cur_ptnlh = NULL;
		info.nniMoves[0].ptnlh = NULL;
		info.nniMoves[1].ptnlh = NULL;
		part_info.push_back(info);
		Alignment *aln = super_aln->concatenateAlignments(*it);
		PhyloTree *tree = super_tree->extractSubtree(*it);
		tree->setAlignment(aln);
		tree_vec.push_back(tree);
	}

	for (PhyloSuperTree::reverse_iterator tit = super_tree->rbegin(); tit != super_tree->rend(); tit++)
		delete (*tit);
	super_tree->clear();
	super_tree->insert(super_tree->end(), tree_vec.begin(), tree_vec.end());
	super_tree->part_info = part_info;

	delete super_tree->aln;
	super_tree->aln = new SuperAlignment(super_tree);
}

void printModelFile(ostream &fmodel, Params &params, PhyloTree *tree, ModelInfo &info, string &set_name) {
	string sitelh_file = params.out_prefix;
	sitelh_file += ".sitelh";
	SeqType seq_type = tree->aln->seq_type;
	if (tree->isSuperTree())
		seq_type = ((PhyloSuperTree*)tree)->front()->aln->seq_type;

    fmodel.precision(4);
    fmodel << fixed;
    if (set_name != "")
        fmodel << set_name << "\t";
    fmodel << info.name << "\t" << info.df << "\t" << info.logl << "\t" << info.tree_len;
    if (seq_type == SEQ_DNA) {
        int nrates = tree->getModel()->getNumRateEntries();
        double *rate_mat = new double[nrates];
        tree->getModel()->getRateMatrix(rate_mat);
        for (int rate = 0; rate < nrates; rate++)
            fmodel << "\t" << rate_mat[rate];
        delete [] rate_mat;
    }
    if (seq_type == SEQ_DNA || seq_type == SEQ_BINARY) {
        int nstates = (seq_type == SEQ_DNA) ? 4 : 2;
        double *freqs = new double[nstates];
        tree->getModel()->getStateFrequency(freqs);
        for (int freq = 0; freq < nstates; freq++)
            fmodel << "\t" << freqs[freq];
        delete [] freqs;
    }
    double alpha = tree->getRate()->getGammaShape();
    fmodel << "\t";
    if (alpha > 0) fmodel << alpha; else fmodel << "NA";
    fmodel << "\t";
    double pinvar = tree->getRate()->getPInvar();
    if (pinvar > 0) fmodel << pinvar; else fmodel << "NA";
    fmodel << "\t";
//    tree->printTree(fmodel);
    fmodel << info.tree;
    fmodel << endl;
    fmodel.precision(4);
    const char *model_name = (params.print_site_lh) ? info.name.c_str() : NULL;
    if (params.print_site_lh)
        printSiteLh(sitelh_file.c_str(), tree, NULL, true, model_name);
    if (params.model_test_and_tree) {
        delete tree;
        tree = NULL;
    }
}

/**
 * select models for all partitions
 * @param model_info (IN/OUT) all model information
 * @return total number of parameters
 */
void testPartitionModel(Params &params, PhyloSuperTree* in_tree, vector<ModelInfo> &model_info, ostream &fmodel, ModelsBlock *models_block ) {
//    params.print_partition_info = true;
//    params.print_conaln = true;
	int i = 0;
//	PhyloSuperTree::iterator it;
	DoubleVector lhvec; // log-likelihood for each partition
	DoubleVector dfvec; // number of parameters for each partition
    DoubleVector lenvec; // tree length for each partition
	double lhsum = 0.0;
	int dfsum = 0;
	int ssize = in_tree->getAlnNSite();
	int num_model = 0;
    int total_num_model = in_tree->size();
	if (params.model_name.find("LINK") != string::npos || params.model_name.find("MERGE") != string::npos) {
        double p = params.partfinder_rcluster/100.0;
        total_num_model += round(in_tree->size()*(in_tree->size()-1)*p/2);
        for (i = in_tree->size()-2; i > 0; i--)
            total_num_model += max(round(i*p), 1.0);
    }
    
    double start_time = getRealTime();
    
	cout << "Selecting individual models for " << in_tree->size() << " charsets using " << criterionName(params.model_test_criterion) << "..." << endl;
	//cout << " No. AIC         AICc        BIC         Charset" << endl;
	cout << " No. Model        Score       Charset" << endl;

	lhvec.resize(in_tree->size());
	dfvec.resize(in_tree->size());
	lenvec.resize(in_tree->size());

    double *dist = new double[in_tree->size()*(in_tree->size()-1)/2];
    int *distID = new int[in_tree->size()*(in_tree->size()-1)/2];
    
    // sort partition by computational cost for OpenMP effciency
	for (i = 0; i < in_tree->size(); i++) {
        distID[i] = i;
        Alignment *this_aln = in_tree->at(i)->aln;
        // computation cost is proportional to #sequences, #patterns, and #states
        dist[i] = -((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
    }
    
    if (params.num_threads > 1) 
    {
        quicksort(dist, 0, in_tree->size()-1, distID);
        if (verbose_mode >= VB_MED) {
            for (i = 0; i < in_tree->size(); i++) {
                cout << i+1 << "\t" << in_tree->part_info[distID[i]].name << endl;
            }
        }
    }

#ifdef _OPENMP
//        for (i = 0; i < in_tree->size(); i++)
//            cout << distID[i]+1 << "\t" << in_tree->part_info[distID[i]].name << "\t" << -dist[i] << endl;
#pragma omp parallel for private(i) schedule(dynamic) reduction(+: lhsum, dfsum)
#endif
	for (int j = 0; j < in_tree->size(); j++) {
        i = distID[j];
        PhyloTree *this_tree = in_tree->at(i);
		// scan through models for this partition, assuming the information occurs consecutively
		vector<ModelInfo> part_model_info;
		extractModelInfo(in_tree->part_info[i].name, model_info, part_model_info);
        stringstream this_fmodel;
		// do the computation
//#ifdef _OPENMP
		string model = testModel(params, this_tree, part_model_info, this_fmodel, models_block, in_tree->part_info[i].name);
//#else
//		string model = testModel(params, this_tree, part_model_info, fmodel, in_tree->part_info[i].name);
//#endif
		double score = computeInformationScore(part_model_info[0].logl,part_model_info[0].df,
				this_tree->getAlnNSite(),params.model_test_criterion);
		in_tree->part_info[i].model_name = model;
		lhsum += (lhvec[i] = part_model_info[0].logl);
		dfsum += (dfvec[i] = part_model_info[0].df);
        lenvec[i] = part_model_info[0].tree_len;

#ifdef _OPENMP
#pragma omp critical
#endif
        {
//#ifdef _OPENMP
            fmodel << this_fmodel.str();
//#endif
            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            cout.width(12);
            cout << left << model << " ";
            cout.width(11);
            cout << score << " " << in_tree->part_info[i].name;
            if (num_model >= 10) {
                double remain_time = (total_num_model-num_model)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " (" 
                    << convert_time(remain_time) << " left)";
            }
            cout << endl;
            replaceModelInfo(model_info, part_model_info);
        }
    }

	if (params.model_name.find("LINK") == string::npos && params.model_name.find("MERGE") == string::npos) {
		in_tree->printBestPartition((string(params.out_prefix) + ".best_scheme.nex").c_str());
		in_tree->printBestPartitionRaxml((string(params.out_prefix) + ".best_scheme").c_str());
        delete [] distID;
        delete [] dist;
		return;
	}

	/* following implements the greedy algorithm of Lanfear et al. (2012) */
//	int part1, part2;
	double inf_score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
	cout << "Full partition model " << criterionName(params.model_test_criterion) << " score: " << inf_score << " (lh=" << lhsum << "  df=" << dfsum << ")" << endl;
	SuperAlignment *super_aln = ((SuperAlignment*)in_tree->aln);
	vector<IntVector> gene_sets;
	gene_sets.resize(in_tree->size());
	StrVector model_names;
	model_names.resize(in_tree->size());
	StrVector greedy_model_trees;
	greedy_model_trees.resize(in_tree->size());
	for (i = 0; i < gene_sets.size(); i++) {
		gene_sets[i].push_back(i);
		model_names[i] = in_tree->part_info[i].model_name;
		greedy_model_trees[i] = in_tree->part_info[i].name;
	}
	cout << "Merging models to increase model fit (about " << total_num_model << " total partition schemes)..." << endl;
	int prev_part = -1;
	while (gene_sets.size() >= 2) {
		// stepwise merging charsets
		double new_score = DBL_MAX;
		double opt_lh = 0.0;
		int opt_df = 0;
        double opt_treelen = 0.0;
		int opt_part1 = 0, opt_part2 = 1;
		IntVector opt_merged_set;
		string opt_set_name = "";
		string opt_model_name = "";
        int num_pairs = 0;
        // 2015-06-24: begin rcluster algorithm
        // compute distance between gene_sets
		for (int part1 = 0; part1 < gene_sets.size()-1; part1++)
			for (int part2 = part1+1; part2 < gene_sets.size(); part2++)
			if (super_aln->partitions[gene_sets[part1][0]]->seq_type == super_aln->partitions[gene_sets[part2][0]]->seq_type)
            {
				// only merge partitions of the same data type
                dist[num_pairs] = fabs(lenvec[part1] - lenvec[part2]);
                distID[num_pairs] = (part1 << 16) | part2;
                num_pairs++;
            }
        if (num_pairs > 0 && params.partfinder_rcluster < 100) {
            // sort distance
            quicksort(dist, 0, num_pairs-1, distID);
            num_pairs = (int)round(num_pairs * (params.partfinder_rcluster/100.0));
            if (num_pairs <= 0) num_pairs = 1;
        }
        // sort partition by computational cost for OpenMP effciency
        for (i = 0; i < num_pairs; i++) {
            // computation cost is proportional to #sequences, #patterns, and #states
            Alignment *this_aln = in_tree->at(distID[i] >> 16)->aln;
            dist[i] = -((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
            this_aln = in_tree->at(distID[i] & ((1<<16)-1))->aln;
            dist[i] -= ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
        }
        if (params.num_threads > 1 && num_pairs >= 1)
            quicksort(dist, 0, num_pairs-1, distID);

#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(dynamic)
#endif
        for (int pair = 0; pair < num_pairs; pair++) {
            int part1 = distID[pair] >> 16;
            int part2 = distID[pair] & ((1<<16)-1);
            assert(part1 != part2);
            IntVector merged_set;
            merged_set.insert(merged_set.end(), gene_sets[part1].begin(), gene_sets[part1].end());
            merged_set.insert(merged_set.end(), gene_sets[part2].begin(), gene_sets[part2].end());
            string set_name = "";
            for (i = 0; i < merged_set.size(); i++) {
                if (i > 0)
                    set_name += "+";
                set_name += in_tree->part_info[merged_set[i]].name;
            }
            string model = "";
            double logl = 0.0;
            int df = 0;
            double treelen = 0.0;
            bool done_before = false;
            if (prev_part >= 0 && part1 != prev_part && part2 != prev_part) {
                // if pairs previously examined, reuse the information
                for (vector<ModelInfo>::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
                    if (mit->set_name == set_name) {
                        model = mit->name;
                        logl = mit->logl;
                        df = mit->df;
                        treelen = mit->tree_len;
                        done_before = true;
                        break;
                    }
            }
            vector<ModelInfo> part_model_info;
            stringstream this_fmodel;
            if (!done_before) {
                Alignment *aln = super_aln->concatenateAlignments(merged_set);
                PhyloTree *tree = in_tree->extractSubtree(merged_set);
                tree->setAlignment(aln);
                extractModelInfo(set_name, model_info, part_model_info);
//#ifdef _OPENMP
                model = testModel(params, tree, part_model_info, this_fmodel, models_block, set_name);
//#else
//                model = testModel(params, tree, part_model_info, fmodel, set_name);
//#endif
                logl = part_model_info[0].logl;
                df = part_model_info[0].df;
                treelen = part_model_info[0].tree_len;
                delete tree;
                delete aln;
            }
            double lhnew = lhsum - lhvec[part1] - lhvec[part2] + logl;
            int dfnew = dfsum - dfvec[part1] - dfvec[part2] + df;
            double score = computeInformationScore(lhnew, dfnew, ssize, params.model_test_criterion);
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				if (!done_before) {
//#ifdef _OPENMP
                    fmodel << this_fmodel.str();
//#endif
					replaceModelInfo(model_info, part_model_info);
                    num_model++;
					cout.width(4);
					cout << right << num_model << " ";
					cout.width(12);
					cout << left << model << " ";
					cout.width(11);
					cout << score << " " << set_name;
                    if (num_model >= 10) {
                        double remain_time = max(total_num_model-num_model, 0)*(getRealTime()-start_time)/num_model;
                        cout << "\t" << convert_time(getRealTime()-start_time) << " (" 
                            << convert_time(remain_time) << " left)";
                    }
                    cout << endl;
				}
				if (score < new_score) {
					new_score = score;
					opt_part1 = part1;
					opt_part2 = part2;
					opt_lh = logl;
					opt_df = df;
                    opt_treelen = treelen;
					opt_merged_set = merged_set;
					opt_set_name = set_name;
					opt_model_name = model;
				}
			}

        }
		if (new_score >= inf_score) break;
		inf_score = new_score;

		lhsum = lhsum - lhvec[opt_part1] - lhvec[opt_part2] + opt_lh;
		dfsum = dfsum - dfvec[opt_part1] - dfvec[opt_part2] + opt_df;
		cout << "Merging " << opt_set_name << " with " << criterionName(params.model_test_criterion) << " score: " << new_score << " (lh=" << lhsum << "  df=" << dfsum << ")" << endl;
		// change entry opt_part1 to merged one
		gene_sets[opt_part1] = opt_merged_set;
		lhvec[opt_part1] = opt_lh;
		dfvec[opt_part1] = opt_df;
        lenvec[opt_part1] = opt_treelen;
		model_names[opt_part1] = opt_model_name;
		greedy_model_trees[opt_part1] = "(" + greedy_model_trees[opt_part1] + "," + greedy_model_trees[opt_part2] + ")" +
				convertIntToString(in_tree->size()-gene_sets.size()+1) + ":" + convertDoubleToString(inf_score);
		prev_part = opt_part1;

		// delete entry opt_part2
		lhvec.erase(lhvec.begin() + opt_part2);
		dfvec.erase(dfvec.begin() + opt_part2);
		lenvec.erase(lenvec.begin() + opt_part2);
		gene_sets.erase(gene_sets.begin() + opt_part2);
		model_names.erase(model_names.begin() + opt_part2);
		greedy_model_trees.erase(greedy_model_trees.begin() + opt_part2);
	}

	string final_model_tree;
	if (greedy_model_trees.size() == 1)
		final_model_tree = greedy_model_trees[0];
	else {
		final_model_tree = "(";
		for (i = 0; i < greedy_model_trees.size(); i++) {
			if (i>0)
				final_model_tree += ",";
			final_model_tree += greedy_model_trees[i];
		}
		final_model_tree += ")";
	}

	cout << "BEST-FIT PARTITION MODEL: " << endl;
	cout << "  charpartition " << criterionName(params.model_test_criterion) << " = ";
	for (i = 0; i < gene_sets.size(); i++) {
		if (i > 0)
			cout << ", ";
		cout << model_names[i] << ":";
		for (int j = 0; j < gene_sets[i].size(); j++) {
			cout << " " << in_tree->part_info[gene_sets[i][j]].name;
		}
	}
	cout << ";" << endl;
	cout << "Agglomerative model selection: " << final_model_tree << endl;
    
    delete [] distID;
    delete [] dist;
	mergePartitions(in_tree, gene_sets, model_names);
	in_tree->printBestPartition((string(params.out_prefix) + ".best_scheme.nex").c_str());
	in_tree->printBestPartitionRaxml((string(params.out_prefix) + ".best_scheme").c_str());
}

bool isMixtureModel(ModelsBlock *models_block, string &model_str) {
    size_t mix_pos;
    for (mix_pos = 0; mix_pos < model_str.length(); mix_pos++) {
        size_t next_mix_pos = model_str.find_first_of("+*", mix_pos);
        string sub_model_str = model_str.substr(mix_pos, next_mix_pos-mix_pos);
        if (models_block->findMixModel(sub_model_str))
            return true;
        if (next_mix_pos == string::npos)
            break;
        mix_pos = next_mix_pos;
    }
    return false;
}

string testModel(Params &params, PhyloTree* in_tree, vector<ModelInfo> &model_info, ostream &fmodel, ModelsBlock *models_block,
    string set_name, bool print_mem_usage) 
{
	SeqType seq_type = in_tree->aln->seq_type;
	if (in_tree->isSuperTree())
		seq_type = ((PhyloSuperTree*)in_tree)->front()->aln->seq_type;
	if (seq_type == SEQ_UNKNOWN)
		outError("Unknown data for model testing.");
	string fmodel_str = params.out_prefix;
	fmodel_str += ".model";
	string sitelh_file = params.out_prefix;
	sitelh_file += ".sitelh";
	in_tree->params = &params;
	StrVector model_names;
	int max_cats = getModelList(params, in_tree->aln, model_names, params.model_test_separate_rate);
	int model;

    if (print_mem_usage) {
        uint64_t mem_size = in_tree->getMemoryRequired(max_cats);
        cout << "NOTE: MODEL SELECTION REQUIRES " << (mem_size / 1024) / 1024
                << " MB MEMORY!" << endl;
        if (mem_size >= getMemorySize()) {
            outError("Memory required exceeds your computer RAM size!");
        }
#ifdef BINARY32
        if (mem_size >= 2000000000) {
            outError("Memory required exceeds 2GB limit of 32-bit executable");
        }
#endif
    }


	string best_model = "";
	/* first check the model file */

	if (in_tree->isSuperTree()) {
		// select model for each partition
		PhyloSuperTree *stree = (PhyloSuperTree*)in_tree;
		testPartitionModel(params, stree, model_info, fmodel, models_block);
//        stree->linkTrees();
        stree->mapTrees();
		string res_models = "";
		for (vector<PartitionInfo>::iterator it = stree->part_info.begin(); it != stree->part_info.end(); it++) {
			if (it != stree->part_info.begin()) res_models += ",";
			res_models += (*it).model_name;
		}
		return res_models;
	}

	in_tree->optimize_by_newton = params.optimize_by_newton;
	in_tree->setLikelihoodKernel(params.SSE);

//    int num_rate_classes = 3 + params.max_rate_cats;

	RateHeterogeneity ** rate_class = new RateHeterogeneity*[4];
	rate_class[0] = new RateHeterogeneity();
	rate_class[1] = new RateInvar(-1, NULL);
	rate_class[2] = new RateGamma(params.num_rate_cats, params.gamma_shape, params.gamma_median, NULL);
	rate_class[3] = new RateGammaInvar(params.num_rate_cats, params.gamma_shape, params.gamma_median, -1, params.optimize_model_rate_joint, NULL);
    
    RateFree ** rate_class_free = new RateFree*[params.max_rate_cats-1];
    
    for (model = 0; model < params.max_rate_cats-1; model++)
        rate_class_free[model] = new RateFree(model+2, params.gamma_shape, "", false, params.optimize_alg, NULL);
        
	ModelGTR *subst_model = NULL;
	if (seq_type == SEQ_BINARY)
		subst_model = new ModelBIN("JC2", "", FREQ_UNKNOWN, "", in_tree);
	else if (seq_type == SEQ_DNA)
		subst_model = new ModelDNA("JC", "", FREQ_UNKNOWN, "", in_tree);
	else if (seq_type == SEQ_PROTEIN)
		subst_model = new ModelProtein("WAG", "", FREQ_UNKNOWN, "", in_tree);
	else if (seq_type == SEQ_MORPH)
		subst_model = new ModelMorphology("MK", "", FREQ_UNKNOWN, "", in_tree);
	else if (seq_type == SEQ_CODON)
		subst_model = new ModelCodon("GY", "", FREQ_UNKNOWN, "", in_tree, false);

	assert(subst_model);

	ModelFactory *model_fac = new ModelFactory();
	model_fac->joint_optimize = params.optimize_model_rate_joint;

	int ssize = in_tree->aln->getNSite(); // sample size
	if (params.model_test_sample_size)
		ssize = params.model_test_sample_size;
	if (set_name == "") {
		cout << "Testing " << model_names.size() << " "
			<< ((seq_type == SEQ_BINARY) ? "binary" : ((seq_type == SEQ_DNA) ? "DNA" :
				((seq_type == SEQ_PROTEIN) ? "protein": ((seq_type == SEQ_CODON) ? "codon": "morphological"))))
			<< " models (sample size: " << ssize << ") ..." << endl;
        if (params.model_test_and_tree == 0)
            cout << " No. Model         -LnL         df  AIC          AICc         BIC" << endl;
	}
	if (params.print_site_lh) {
		ofstream sitelh_out(sitelh_file.c_str());
		if (!sitelh_out.is_open())
			outError("Cannot write to file ", sitelh_file);
		sitelh_out << model_names.size() << " " << in_tree->getAlnNSite() << endl;
		sitelh_out.close();
	}
	vector<ModelInfo>::iterator it;
	for (it = model_info.begin(); it != model_info.end(); it++) {
		it->AIC_score = DBL_MAX;
		it->AICc_score = DBL_MAX;
		it->BIC_score = DBL_MAX;
	}

	uint64_t RAM_requirement = 0;
    int model_aic = -1, model_aicc = -1, model_bic = -1;
    string prev_tree_string = "";
    int prev_model_id = -1;
    int skip_model = 0;

	for (model = 0; model < model_names.size(); model++) {
		//cout << model_names[model] << endl;
        if (model_names[model][0] == '+') {
            // now switching to test rate heterogeneity
            if (best_model == "")
                switch (params.model_test_criterion) {
                case MTC_AIC:
                    best_model = model_info[model_aic].name;
                    break;
                case MTC_AICC:
                    best_model = model_info[model_aicc].name;
                    break;
                case MTC_BIC:
                    best_model = model_info[model_bic].name;
                    break;
                default: assert(0);
                }
            model_names[model] = best_model + model_names[model];
        }
		PhyloTree *tree = in_tree;
        ModelFactory *this_model_fac = NULL;
        bool mixture_model = false;
        int ncat = 0;
        string orig_name = params.model_name;
        
        if (isMixtureModel(models_block, model_names[model])) {
            // mixture model
            try {
                mixture_model = true;
                params.model_name = model_names[model];
                this_model_fac = new ModelFactory(params, tree, models_block);
                tree->setModelFactory(this_model_fac);
                tree->setModel(this_model_fac->model);
                tree->setRate(this_model_fac->site_rate);
                tree->deleteAllPartialLh();
                tree->initializeAllPartialLh();
                RAM_requirement = max(RAM_requirement, tree->getMemoryRequired());
            } catch (string &str) {
                outError("Invalid -madd model " + model_names[model] + ": " + str);
            }
        } else {
            // kernel might be changed if mixture model was tested
            in_tree->setLikelihoodKernel(params.SSE);
            // normal model
            if (model_names[model].find("+ASC") != string::npos) {
                model_fac->unobserved_ptns = in_tree->aln->getUnobservedConstPatterns();
                if (model_fac->unobserved_ptns.size() == 0) {
                    cout.width(3);
                    cout << right << model+1 << "  ";
                    cout.width(13);
                    cout << left << model_names[model] << " ";                
                    cout << "Skipped since +ASC is not applicable" << endl;
                    continue;
                }
                tree->aln->buildSeqStates(true);
                if (model_fac->unobserved_ptns.size() < tree->aln->getNumNonstopCodons())
                    outError("Invalid use of +ASC because constant patterns are observed in the alignment");
            } else {
                model_fac->unobserved_ptns = "";
                tree->aln->buildSeqStates(false);
            }
            // initialize tree
            // initialize model
            subst_model->setTree(tree);
            StateFreqType freq_type = FREQ_UNKNOWN;
            if (model_names[model].find("+F1X4") != string::npos)
                freq_type = FREQ_CODON_1x4;
            else if (model_names[model].find("+F3X4C") != string::npos)
                freq_type = FREQ_CODON_3x4C;
            else if (model_names[model].find("+F3X4") != string::npos)
                freq_type = FREQ_CODON_3x4;
            else if (model_names[model].find("+FQ") != string::npos)
                freq_type = FREQ_EQUAL;
            else if (model_names[model].find("+FO") != string::npos)
                freq_type = FREQ_ESTIMATE;
            else if (model_names[model].find("+FU") != string::npos)
                freq_type = FREQ_USER_DEFINED;
            else if (model_names[model].find("+F") != string::npos)
                freq_type = FREQ_EMPIRICAL;
                
            subst_model->init(model_names[model].substr(0, model_names[model].find('+')).c_str(), "", freq_type, "");
            tree->params = &params;

            tree->setModel(subst_model);
            // initialize rate
            size_t pos;
            if ((pos = model_names[model].find("+R")) != string::npos) {
                ncat = params.num_rate_cats;
                if (model_names[model].length() > pos+2 && isdigit(model_names[model][pos+2])) {
                    ncat = convert_int(model_names[model].c_str() + pos+2);
    //                tree->getRate()->setNCategory(ncat);
                }
                if (ncat <= 1) outError("Number of rate categories for " + model_names[model] + " is <= 1");
                if (ncat > params.max_rate_cats)
                    outError("Number of rate categories for " + model_names[model] + " exceeds " + convertIntToString(params.max_rate_cats));
                tree->setRate(rate_class_free[ncat-2]);
            } else if (model_names[model].find("+I") != string::npos && (pos = model_names[model].find("+G")) != string::npos) {
                tree->setRate(rate_class[3]);
                if (model_names[model].length() > pos+2 && isdigit(model_names[model][pos+2])) {
                    int ncat = convert_int(model_names[model].c_str() + pos+2);
                    if (ncat < 1) outError("Wrong number of category for +G in " + model_names[model]);
                    tree->getRate()->setNCategory(ncat);
                }
            } else if ((pos = model_names[model].find("+G")) != string::npos) {
                tree->setRate(rate_class[2]);
                if (model_names[model].length() > pos+2 && isdigit(model_names[model][pos+2])) {
                    ncat = convert_int(model_names[model].c_str() + pos+2);
                    if (ncat < 1) outError("Wrong number of category for +G in " + model_names[model]);
                    tree->getRate()->setNCategory(ncat);
                }
            } else if (model_names[model].find("+I") != string::npos)
                tree->setRate(rate_class[1]);
            else
                tree->setRate(rate_class[0]);

            tree->getRate()->setTree(tree);

            // initialize model factory
            model_fac->model = subst_model;
            model_fac->site_rate = tree->getRate();
            tree->setModelFactory(model_fac);
        }
        
        tree->clearAllPartialLH();


		// optimize model parameters
		ModelInfo info;        
		info.set_name = set_name;
		info.df = tree->getModelFactory()->getNParameters();
        if (mixture_model)
            info.name = model_names[model];
        else
            info.name = tree->getModelName();
		int model_id = -1;
        if (skip_model) {
            assert(prev_model_id>=0);
            size_t pos_r = info.name.find("+R");
            if (pos_r == string::npos || info.name.substr(0, pos_r) != model_info[prev_model_id].name.substr(0, pos_r))
                skip_model = 0;
        }
		for (int i = 0; i < model_info.size(); i++)
			if (info.name == model_info[i].name) {
				model_id = i;
				if (info.df != model_info[i].df)
					outError("Inconsistent model file " + fmodel_str + ", please rerun using -mredo option");
				break;
			}
		if (model_id >= 0) {
			info.logl = model_info[model_id].logl;
            info.tree_len = model_info[model_id].tree_len;
            info.tree = model_info[model_id].tree;
            prev_tree_string = model_info[model_id].tree;
        } else if (skip_model) {
            assert(prev_model_id >= 0);
            if (prev_model_id >= 0) {
            info.logl = model_info[prev_model_id].logl;
            info.tree_len = model_info[prev_model_id].tree_len;
//            info.tree = model_info[prev_model_id].tree;
//            prev_tree_string = model_info[prev_model_id].tree;
//            cout << "Skipped " << info.name << endl;
            }
		} else {
            if (params.model_test_and_tree) {
                string original_model = params.model_name;
                params.model_name = model_names[model];
                char *orig_user_tree = params.user_file;
                string new_user_tree = (string)params.out_prefix+".treefile";
                if (params.model_test_and_tree == 1 && model>0 && fileExists(new_user_tree)) {
                    params.user_file = (char*)new_user_tree.c_str();
                }
                if (in_tree->isSuperTree()) {
                    outError("-mtree option is not supported for partition model");
                }
                IQTree *iqtree = new IQTree(in_tree->aln);
                // set checkpoint
                iqtree->setCheckpoint(in_tree->getCheckpoint());
                iqtree->num_precision = in_tree->num_precision;
                
                cout << endl << "===> Testing model " << model+1 << ": " << params.model_name << endl;
                runTreeReconstruction(params, original_model, *iqtree, model_info);
                info.logl = iqtree->computeLikelihood();
                info.tree_len = iqtree->treeLength();
//                info.tree = iqtree->getTreeString();
                params.model_name = original_model;
                params.user_file = orig_user_tree;
                tree = iqtree;

                // clear all checkpointed information
                Checkpoint *newCheckpoint = new Checkpoint;
                tree->getCheckpoint()->getSubCheckpoint(newCheckpoint, "iqtree");
                tree->getCheckpoint()->clear();
                tree->getCheckpoint()->insert(newCheckpoint->begin(), newCheckpoint->end());
                tree->getCheckpoint()->putBool("finished", false);
                tree->getCheckpoint()->dump(true);
                delete newCheckpoint;

            } else {
                if (tree->getMemoryRequired() > RAM_requirement) {
                    tree->deleteAllPartialLh();
                    RAM_requirement = tree->getMemoryRequired();
                }
                tree->initializeAllPartialLh();
                if (prev_tree_string != "") {
                    tree->readTreeString(prev_tree_string);
                }
                prev_tree_string = "";
                if (model_fac->unobserved_ptns.size() > 0 && tree->aln->seq_type == SEQ_PROTEIN) {
                    // treatment for +ASC for protein data
                    tree->fixNegativeBranch(true);
                    tree->clearAllPartialLH();
                }
                info.logl = tree->getModelFactory()->optimizeParameters(false, false, TOL_LIKELIHOOD_MODELTEST, TOL_GRADIENT_MODELTEST);
                info.tree_len = tree->treeLength();
                if (prev_model_id >= 0) {
                    // check stop criterion for +R
                    size_t prev_pos_r = model_info[prev_model_id].name.find("+R");
                    size_t pos_r = info.name.find("+R");
                    if ( prev_pos_r != string::npos &&  pos_r != string::npos && 
                        model_info[prev_model_id].name.substr(0,prev_pos_r) == info.name.substr(0, pos_r) && 
                        info.logl < model_info[prev_model_id].logl) 
                    {
                        if (verbose_mode >= VB_MED)
                            cout << "reoptimizing from previous parameters of +R...." << endl;
                        assert(ncat >= 3);
                        rate_class_free[ncat-2]->setRateAndProp(rate_class_free[ncat-3]);
                        info.logl = tree->getModelFactory()->optimizeParameters(false, false, TOL_LIKELIHOOD_MODELTEST, TOL_GRADIENT_MODELTEST);
                        info.tree_len = tree->treeLength();                        
                    }
                }
//                info.tree = tree->getTreeString();
            }
			// print information to .model file
            info.tree = tree->getTreeString();
            printModelFile(fmodel, params, tree, info, set_name);
		}
		computeInformationScores(info.logl, info.df, ssize, info.AIC_score, info.AICc_score, info.BIC_score);
        if (prev_model_id >= 0) {
            // check stop criterion for +R
            size_t prev_pos_r = model_info[prev_model_id].name.find("+R");
            size_t pos_r = info.name.find("+R");
            if ( prev_pos_r != string::npos &&  pos_r != string::npos && 
            model_info[prev_model_id].name.substr(0,prev_pos_r) == info.name.substr(0, pos_r)) {
                switch (params.model_test_criterion) {
                case MTC_ALL:
                    if (info.AIC_score > model_info[prev_model_id].AIC_score && info.AICc_score > model_info[prev_model_id].AICc_score &&
                        info.BIC_score > model_info[prev_model_id].BIC_score) {
                        // skip remaining model
                        skip_model++;
                    }
                    break;
                case MTC_AIC:
                    if (info.AIC_score > model_info[prev_model_id].AIC_score) {
                        // skip remaining model
                        skip_model++;
                    }
                    break;
                case MTC_AICC:
                    if (info.AICc_score > model_info[prev_model_id].AICc_score) {
                        // skip remaining model
                        skip_model++;
                    }
                    break;
                case MTC_BIC:
                    if (info.BIC_score > model_info[prev_model_id].BIC_score) {
                        // skip remaining model
                        skip_model++;
                    }
                    break;
                }
            }
        }
        if (skip_model > 1)
            info.AIC_score = DBL_MAX;
        
		if (model_id >= 0) {
			model_info[model_id] = info;
		} else {
			model_info.push_back(info);
            model_id = model_info.size()-1;
		}
		if (model_aic < 0 || model_info[model_id].AIC_score < model_info[model_aic].AIC_score)
			model_aic = model_id;
		if (model_aicc < 0 || model_info[model_id].AICc_score < model_info[model_aicc].AICc_score)
			model_aicc = model_id;
		if (model_bic < 0 || model_info[model_id].BIC_score < model_info[model_bic].BIC_score)
			model_bic = model_id;
        
        if (mixture_model) {
            delete this_model_fac->model;
            delete this_model_fac->site_rate;
            delete this_model_fac;
            this_model_fac = NULL;
            params.model_name = orig_name;
        }
        
        in_tree->setModel(NULL);
        in_tree->setModelFactory(NULL);
        in_tree->setRate(NULL);

        prev_model_id = model_id;

		if (set_name != "") continue;

		cout.width(3);
		cout << right << model+1 << "  ";
		cout.width(13);
		cout << left << info.name << " ";
        
        if (skip_model > 1) {
            cout << "Skipped " << endl;
            continue;
        }
        
		cout.precision(3);
		cout << fixed;
		cout.width(12);
		cout << -info.logl << " ";
		cout.width(3);
		cout << info.df << " ";
		cout.width(12);
		cout << info.AIC_score << " ";
		cout.width(12);
		cout << info.AICc_score << " " << info.BIC_score;
		cout << endl;


	}

    if (model_bic < 0) 
        outError("No models were examined! Please check messages above");

	//cout.unsetf(ios::fixed);
	/*
	for (it = model_info.begin(); it != model_info.end(); it++)
		computeInformationScores(it->logl, it->df, ssize, it->AIC_score, it->AICc_score, it->BIC_score);
*/
//	for (it = model_info.begin(), model = 0; it != model_info.end(); it++, model++) {
//		if ((*it).AIC_score < model_info[model_aic].AIC_score)
//			model_aic = model;
//		if ((*it).AICc_score < model_info[model_aicc].AICc_score)
//			model_aicc = model;
//		if ((*it).BIC_score < model_info[model_bic].BIC_score)
//			model_bic = model;
//	}
	if (set_name == "") {
		cout << "Akaike Information Criterion:           " << model_info[model_aic].name << endl;
		cout << "Corrected Akaike Information Criterion: " << model_info[model_aicc].name << endl;
		cout << "Bayesian Information Criterion:         " << model_info[model_bic].name << endl;
	} else {
		/*
		cout.width(11);
		cout << left << model_info[model_aic].name << " ";
		cout.width(11);
		cout << left << model_info[model_aicc].name << " ";
		cout.width(11);
		cout << left << model_info[model_bic].name << " ";
		cout << set_name;
		*/
	}

	/* computing model weights */
	double AIC_sum = 0.0, AICc_sum = 0.0, BIC_sum = 0.0;
	for (it = model_info.begin(); it != model_info.end(); it++) {
		it->AIC_weight = exp(-0.5*(it->AIC_score-model_info[model_aic].AIC_score));
		it->AICc_weight = exp(-0.5*(it->AICc_score-model_info[model_aicc].AICc_score));
		it->BIC_weight = exp(-0.5*(it->BIC_score-model_info[model_bic].BIC_score));
		it->AIC_conf = false;
		it->AICc_conf = false;
		it->BIC_conf = false;
		AIC_sum += it->AIC_weight;
		AICc_sum += it->AICc_weight;
		BIC_sum += it->BIC_weight;
	}
	for (it = model_info.begin(); it != model_info.end(); it++) {
		it->AIC_weight /= AIC_sum;
		it->AICc_weight /= AICc_sum;
		it->BIC_weight /= BIC_sum;
	}

	int *model_rank = new int[model_info.size()];
	double *scores = new double[model_info.size()];

	/* compute confidence set for BIC */
    AIC_sum = 0.0;
    AICc_sum = 0.0;
    BIC_sum = 0.0;
	for (model = 0; model < model_info.size(); model++)
		scores[model] = model_info[model].BIC_score;
	sort_index(scores, scores+model_info.size(), model_rank);
	for (model = 0; model < model_info.size(); model++) {
		model_info[model_rank[model]].BIC_conf = true;
		BIC_sum += model_info[model_rank[model]].BIC_weight;
		if (BIC_sum > 0.95) break;
	}
	/* compute confidence set for AIC */
	for (model = 0; model < model_info.size(); model++)
		scores[model] = model_info[model].AIC_score;
	sort_index(scores, scores+model_info.size(), model_rank);
	for (model = 0; model < model_info.size(); model++) {
		model_info[model_rank[model]].AIC_conf = true;
		AIC_sum += model_info[model_rank[model]].AIC_weight;
		if (AIC_sum > 0.95) break;
	}

	/* compute confidence set for AICc */
	for (model = 0; model < model_info.size(); model++)
		scores[model] = model_info[model].AICc_score;
	sort_index(scores, scores+model_info.size(), model_rank);
	for (model = 0; model < model_info.size(); model++) {
		model_info[model_rank[model]].AICc_conf = true;
		AICc_sum += model_info[model_rank[model]].AICc_weight;
		if (AICc_sum > 0.95) break;
	}

    string best_tree; // BQM 2015-07-21: With Lars find best model
	/* sort models by their scores */
	switch (params.model_test_criterion) {
	case MTC_AIC:
		for (model = 0; model < model_info.size(); model++)
			scores[model] = model_info[model].AIC_score;
		best_model = model_info[model_aic].name;
        best_tree = model_info[model_aic].tree;
		break;
	case MTC_AICC:
		for (model = 0; model < model_info.size(); model++)
			scores[model] = model_info[model].AICc_score;
		best_model = model_info[model_aicc].name;
        best_tree = model_info[model_aicc].tree;
		break;
	case MTC_BIC:
		for (model = 0; model < model_info.size(); model++)
			scores[model] = model_info[model].BIC_score;
		best_model = model_info[model_bic].name;
        best_tree = model_info[model_bic].tree;
		break;
    default: assert(0);
	}
	sort_index(scores, scores + model_info.size(), model_rank);

	vector<ModelInfo> sorted_info;
	for (model = 0; model < model_info.size(); model++)
		sorted_info.push_back(model_info[model_rank[model]]);
	model_info = sorted_info;

	delete [] model_rank;
	delete [] scores;

	delete model_fac;
	delete subst_model;
    int rate_type;
	for (rate_type = 3; rate_type >= 0; rate_type--) {
		delete rate_class[rate_type];
    }
    delete [] rate_class;
    
	for (rate_type = params.max_rate_cats-2; rate_type >= 0; rate_type--) {
		delete rate_class_free[rate_type];
    }
    delete [] rate_class_free;
    
//	delete tree_hetero;
//	delete tree_homo;
	in_tree->deleteAllPartialLh();
    
    // BQM 2015-07-21 with Lars: load the best_tree
//	if (params.model_test_and_tree)
		in_tree->readTreeString(best_tree);

    
	if (set_name == "") {
		cout << "Best-fit model: " << best_model << " chosen according to " << 
            ((params.model_test_criterion == MTC_BIC) ? "BIC" :
			((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc")) << endl;
	}
	if (params.print_site_lh)
		cout << "Site log-likelihoods per model printed to " << sitelh_file << endl;
	return best_model;
}

int countDistinctTrees(const char *filename, bool rooted, IQTree *tree, IntVector &distinct_ids, bool exclude_duplicate) {
	StringIntMap treels;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		// remove the failbit
		in.exceptions(ios::badbit);
		int tree_id;
		for (tree_id = 0; !in.eof(); tree_id++) {
			if (exclude_duplicate) {
				tree->freeNode();
				tree->readTree(in, rooted);
				tree->setAlignment(tree->aln);
				tree->setRootNode((char*)tree->aln->getSeqName(0).c_str());
				StringIntMap::iterator it = treels.end();
				ostringstream ostr;
				tree->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
				it = treels.find(ostr.str());
				if (it != treels.end()) { // already in treels
					distinct_ids.push_back(it->second);
				} else {
					distinct_ids.push_back(-1);
					treels[ostr.str()] = tree_id;
				}
			} else {
				// ignore tree
				char ch;
				do {
					in >> ch;
				} while (!in.eof() && ch != ';');
				distinct_ids.push_back(-1);
			}
			char ch;
			in.exceptions(ios::goodbit);
			(in) >> ch;
			if (in.eof()) break;
			in.unget();
			in.exceptions(ios::failbit | ios::badbit);

		}
		in.close();
	} catch (ios::failure) {
		outError("Cannot read file ", filename);
	}
	if (exclude_duplicate)
		return treels.size();
	else
		return distinct_ids.size();
}

//const double TOL_RELL_SCORE = 0.01;

/*
    Problem: solve the following linear system equation:
    a_1*x + b_1*y = c_1
    a_2*x + b_2*y = c_2
    ....
    a_n*x + b_n*y = c_n
    
becomes minimizing weighted least square:

    sum_k { w_k*[ c_k - (a_k*x + b_k*y) ]^2 }


the solution is:

    x = [(sum_k w_k*b_k*c_k)*(sum_k w_k*a_k*b_k) - (sum_k w_k*a_k*c_k)(sum_k w_k*b_k^2)] / 
        [ (sum_k w_k*a_k*b_k)^2 - (sum_k w_k*a_k^2)*(sum_k w_k*b_k^2) ]
    
    y = [(sum_k w_k*a_k*c_k)*(sum_k w_k*a_k*b_k) - (sum_k w_k*b_k*c_k)(sum_k w_k*a_k^2)] / 
        [ (sum_k w_k*a_k*b_k)^2 - (sum_k w_k*a_k^2)*(sum_k w*k*b_k^2) ]
    
    @param n number of data points
    @param w weight vector of length n
    @param a a value vector of length n
    @param b b value vector of length n
    @param c c value vector of length n
    @param[out] x x-value
    @param[out] y y-value
    @return least square value
*/
void doWeightedLeastSquare(int n, double *w, double *a, double *b, double *c, double &x, double &y) {
    int k;
    double BC = 0.0, AB = 0.0, AC = 0.0, A2 = 0.0, B2 = 0.0;
    double denom;
    for (k = 0; k < n; k++) {
        double wa = w[k]*a[k];
        double wb = w[k]*b[k];
        AB += wa*b[k];
        BC += wb*c[k];
        AC += wa*c[k];
        A2 += wa*a[k];
        B2 += wb*b[k];
    }
    denom = 1.0/(AB*AB - A2*B2);
    x = (BC*AB - AC*B2) * denom;
    y = (AC*AB - BC*A2) * denom;
}

/**
    MLE estimates for AU test
*/
class OptimizationAUTest : public Optimization {

public:

    OptimizationAUTest(double d, double c, int nscales, double *bp, double *rr, double *rr_inv) {
        this->d = d;
        this->c = c;
        this->bp = bp;
        this->rr = rr;
        this->rr_inv = rr_inv;
        this->nscales = nscales;
        
    }

	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return 2; }


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) {
        d = x[1];
        c = x[2];
        double res = 0.0;
        for (int k = 0; k < nscales; k++) {
            double cdf = gsl_cdf_ugaussian_P(d*rr[k] + c*rr_inv[k]);
            res += bp[k] * log(1.0 - cdf) + (1.0-bp[k])*log(cdf);
        }
        return res;
    }

    void optimizeDC() {
        double x[3], lower[3], upper[3];
        bool bound_check[3];
        x[1] = d;
        x[2] = c;
        lower[1] = lower[2] = 1e-4;
        upper[1] = upper[2] = 100.0;
        bound_check[1] = bound_check[2] = false;
        minimizeMultiDimen(x, 2, lower, upper, bound_check, 1e-4);
        d = x[1];
        c = x[2];
    }

    double d, c;
    int nscales;
    double *bp;
    double *rr;
    double *rr_inv;
};

/**
    @param tree_lhs RELL score matrix of size #trees x #replicates
*/
void performAUTest(Params &params, PhyloTree *tree, double *pattern_lhs, vector<TreeInfo> &info) {
    
    /* STEP 1: specify scale factors */
    int nscales = 10;
    double r[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    double rr[] = {sqrt(0.5), sqrt(0.6), sqrt(0.7), sqrt(0.8), sqrt(0.9), 1.0, 
        sqrt(1.1), sqrt(1.2), sqrt(1.3), sqrt(1.4)};
    double rr_inv[] = {sqrt(1/0.5), sqrt(1/0.6), sqrt(1/0.7), sqrt(1/0.8), sqrt(1/0.9), 1.0, 
        sqrt(1/1.1), sqrt(1/1.2), sqrt(1/1.3), sqrt(1/1.4)};
        
    /* STEP 2: compute bootstrap proportion */
    int ntrees = info.size();
    size_t nboot = params.topotest_replicates;
    double nboot_inv = 1.0 / nboot;
    
    int nptn = tree->getAlnNPattern();
    int maxnptn = get_safe_upper_limit(nptn);
    
    double *bp = new double[ntrees*nscales];
    memset(bp, 0, sizeof(double)*ntrees*nscales);
    
    int k, tid, ptn;
#ifdef _OPENMP
    #pragma omp parallel private(k, tid, ptn)
    {
    int *rstream;
    init_random(params.ran_seed + omp_get_thread_num(), false, &rstream);
#else
    int *rstream = randstream;
#endif
    size_t boot;
    int *boot_sample = aligned_alloc<int>(maxnptn);
    memset(boot_sample, 0, maxnptn*sizeof(int));
    
    double *boot_sample_dbl = aligned_alloc<double>(maxnptn);
    
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (k = 0; k < nscales; k++) {
        string str = "SCALE=" + convertDoubleToString(r[k]);    
		for (boot = 0; boot < nboot; boot++) {
			tree->aln->createBootstrapAlignment(boot_sample, str.c_str(), rstream);
            for (ptn = 0; ptn < maxnptn; ptn++)
                boot_sample_dbl[ptn] = boot_sample[ptn];
            double max_lh = -1e20;
            int max_tid = -1;
            for (tid = 0; tid < ntrees; tid++) {
                double *pattern_lh = pattern_lhs + (tid*maxnptn);
                double tree_lh;
#ifdef BINARY32
                tree_lh = tree->dotProductSIMD<double, Vec2d, 2>(pattern_lh, boot_sample_dbl, nptn);
#else
                if (instruction_set >= 7)
                    tree_lh = tree->dotProductSIMD<double, Vec4d, 4>(pattern_lh, boot_sample_dbl, nptn);
                else
                    tree_lh = tree->dotProductSIMD<double, Vec2d, 2>(pattern_lh, boot_sample_dbl, nptn);
#endif
                if (tree_lh > max_lh) {
                    max_lh = tree_lh;
                    max_tid = tid;
                } 
            }
            bp[k*ntrees+max_tid] += nboot_inv;
        }
    }

    aligned_free(boot_sample_dbl);
    aligned_free(boot_sample);

#ifdef _OPENMP
    finish_random(rstream);
    }
#endif

    if (verbose_mode >= VB_MED) {
        cout << "scale";
        for (k = 0; k < nscales; k++)
            cout << "\t" << r[k];
        cout << endl;
        for (tid = 0; tid < ntrees; tid++) {
            cout << tid;
            for (k = 0; k < nscales; k++) {
                cout << "\t" << bp[tid+k*ntrees];
            }
            cout << endl;
        }
    }
    
    /* STEP 3: weighted least square fit */
    
    double *cc = new double[nscales];
    double *w = new double[nscales];
    double *this_bp = new double[nscales];
    cout << "TreeID\tAU\tRSS\td_WLS\tc_WLS\td_MLE\tc_MLE" << endl;
    for (tid = 0; tid < ntrees; tid++) {
        for (k = 0; k < nscales; k++) {
            this_bp[k] = bp[tid + k*ntrees];
            double bp_val = min(max(bp[tid + k*ntrees], nboot_inv),1.0-nboot_inv);
            double bp_cdf = gsl_cdf_ugaussian_Pinv(bp_val);
            double bp_pdf = gsl_ran_ugaussian_pdf(bp_cdf);
            cc[k] = gsl_cdf_ugaussian_Pinv(1.0 - bp_val);
            w[k] = bp_pdf*bp_pdf*nboot / (bp_val*(1.0-bp_val));
        }
        double c, d; // c, d in original paper
        // first obtain d and c by weighted least square
        doWeightedLeastSquare(nscales, w, rr, rr_inv, cc, d, c);
        
        // second, perform MLE estimate of d and c
        OptimizationAUTest mle(d, c, nscales, this_bp, rr, rr_inv);
        mle.optimizeDC();
        
        // compute sum of squared difference
        double rss = 0.0;
        for (k = 0; k < nscales; k++) {
            double diff = cc[k] - (rr[k]*d + rr_inv[k]*c);
            rss += w[k] * diff * diff;
        }
        
        double pchi2 = computePValueChiSquare(rss, nscales-2);
        /* STEP 4: compute p-value according to Eq. 11 */
        info[tid].au_pvalue = 1.0 - gsl_cdf_ugaussian_P(mle.d-mle.c);
        cout << tid+1 << "\t" << info[tid].au_pvalue << "\t" << rss << "\t" << d << "\t" << c << "\t" << mle.d << "\t" << mle.c;
        
        // warning if p-value of chi-square < 0.01 (rss too high)
        if (pchi2 < 0.01) 
            cout << " !!!";
        cout << endl;
    }
    
    delete [] this_bp;
    delete [] w;
    delete [] cc;
    delete [] bp;
}


void evaluateTrees(Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids)
{
	if (!params.treeset_file)
		return;
	cout << endl;
	//MTreeSet trees(params.treeset_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
	cout << "Reading trees in " << params.treeset_file << " ..." << endl;
	int ntrees = countDistinctTrees(params.treeset_file, params.is_rooted, tree, distinct_ids, params.distinct_trees);
	if (ntrees < distinct_ids.size()) {
		cout << "WARNING: " << distinct_ids.size() << " trees detected but only " << ntrees << " distinct trees will be evaluated" << endl;
	} else {
		cout << ntrees << (params.distinct_trees ? " distinct" : "") << " trees detected" << endl;
	}
	if (ntrees == 0) return;
	ifstream in(params.treeset_file);

	//if (trees.size() == 1) return;
	//string tree_file = params.treeset_file;
	string tree_file = params.out_prefix;
	tree_file += ".trees";
	ofstream treeout;
	//if (!params.fixed_branch_length) {
		treeout.open(tree_file.c_str());
	//}
	string score_file = params.out_prefix;
	score_file += ".treelh";
	ofstream scoreout;
	if (params.print_tree_lh)
		scoreout.open(score_file.c_str());
	string site_lh_file = params.out_prefix;
	site_lh_file += ".sitelh";
	if (params.print_site_lh) {
		ofstream site_lh_out(site_lh_file.c_str());
		site_lh_out << ntrees << " " << tree->getAlnNSite() << endl;
		site_lh_out.close();
	}

	double time_start = getRealTime();

	int *boot_samples = NULL;
	int boot;
	//double *saved_tree_lhs = NULL;
	double *tree_lhs = NULL; // RELL score matrix of size #trees x #replicates
	double *pattern_lh = NULL;
	double *pattern_lhs = NULL;
	double *orig_tree_lh = NULL; // Original tree log-likelihoods
	double *max_lh = NULL;
	double *lhdiff_weights = NULL;
	int nptn = tree->getAlnNPattern();
    int maxnptn = get_safe_upper_limit(nptn);
    
	if (params.topotest_replicates && ntrees > 1) {
		size_t mem_size = (size_t)params.topotest_replicates*nptn*sizeof(int) +
				ntrees*params.topotest_replicates*sizeof(double) +
				(nptn + ntrees*3 + params.topotest_replicates*2)*sizeof(double) +
				ntrees*sizeof(TreeInfo) +
				params.do_weighted_test*(ntrees * nptn * sizeof(double) + ntrees*ntrees*sizeof(double));
		cout << "Note: " << ((double)mem_size/1024)/1024 << " MB of RAM required!" << endl;
		if (mem_size > getMemorySize()-100000)
			outWarning("The required memory does not fit in RAM!");
		cout << "Creating " << params.topotest_replicates << " bootstrap replicates..." << endl;
		if (!(boot_samples = new int [params.topotest_replicates*nptn]))
			outError(ERR_NO_MEMORY);
		for (boot = 0; boot < params.topotest_replicates; boot++)
			tree->aln->createBootstrapAlignment(boot_samples + (boot*nptn), params.bootstrap_spec);
		//if (!(saved_tree_lhs = new double [ntrees * params.topotest_replicates]))
		//	outError(ERR_NO_MEMORY);
		if (!(tree_lhs = new double [ntrees * params.topotest_replicates]))
			outError(ERR_NO_MEMORY);
		if (params.do_weighted_test || params.do_au_test) {
			if (!(lhdiff_weights = new double [ntrees * ntrees]))
				outError(ERR_NO_MEMORY);
            pattern_lhs = aligned_alloc<double>(ntrees*maxnptn);
//			if (!(pattern_lhs = new double[ntrees* nptn]))
//				outError(ERR_NO_MEMORY);
		}
        pattern_lh = aligned_alloc<double>(maxnptn);
//		if (!(pattern_lh = new double[nptn]))
//			outError(ERR_NO_MEMORY);
		if (!(orig_tree_lh = new double[ntrees]))
			outError(ERR_NO_MEMORY);
		if (!(max_lh = new double[params.topotest_replicates]))
			outError(ERR_NO_MEMORY);
	}
	int tree_index, tid, tid2;
	info.resize(ntrees);
	//for (MTreeSet::iterator it = trees.begin(); it != trees.end(); it++, tree_index++) {
	for (tree_index = 0, tid = 0; tree_index < distinct_ids.size(); tree_index++) {

		cout << "Tree " << tree_index + 1;
		if (distinct_ids[tree_index] >= 0) {
			cout << " / identical to tree " << distinct_ids[tree_index]+1 << endl;
			// ignore tree
			char ch;
			do {
				in >> ch;
			} while (!in.eof() && ch != ';');
			continue;
		}
		tree->freeNode();
		tree->readTree(in, params.is_rooted);
		tree->setAlignment(tree->aln);
        tree->setRootNode(params.root);
		if (tree->isSuperTree())
			((PhyloSuperTree*) tree)->mapTrees();
//		if ((tree->sse == LK_EIGEN || tree->sse == LK_EIGEN_SSE) && !tree->isBifurcating()) {
//			cout << "NOTE: Changing to old kernel as user tree is multifurcating" << endl;
//			if (tree->sse == LK_EIGEN)
//				tree->changeLikelihoodKernel(LK_NORMAL);
//			else
//				tree->changeLikelihoodKernel(LK_SSE);
//		}

		tree->initializeAllPartialLh();
		tree->fixNegativeBranch(false);
		if (!params.fixed_branch_length) {
			tree->setCurScore(tree->optimizeAllBranches(100, 0.001));
		} else {
			tree->setCurScore(tree->computeLikelihood());
		}
		treeout << "[ tree " << tree_index+1 << " lh=" << tree->getCurScore() << " ]";
		tree->printTree(treeout);
		treeout << endl;
		if (params.print_tree_lh)
			scoreout << tree->getCurScore() << endl;

		cout << " / LogL: " << tree->getCurScore() << endl;

		if (pattern_lh) {
			double curScore = tree->getCurScore();
            memset(pattern_lh, 0, maxnptn*sizeof(double));
			tree->computePatternLikelihood(pattern_lh, &curScore);
			if (params.do_weighted_test || params.do_au_test)
				memcpy(pattern_lhs + tid*maxnptn, pattern_lh, maxnptn*sizeof(double));
		}
		if (params.print_site_lh) {
			string tree_name = "Tree" + convertIntToString(tree_index+1);
			printSiteLh(site_lh_file.c_str(), tree, pattern_lh, true, tree_name.c_str());
		}
		info[tid].logl = tree->getCurScore();

		if (!params.topotest_replicates || ntrees <= 1) {
			tid++;
			continue;
		}
		// now compute RELL scores
		orig_tree_lh[tid] = tree->getCurScore();
		double *tree_lhs_offset = tree_lhs + (tid*params.topotest_replicates);
		for (boot = 0; boot < params.topotest_replicates; boot++) {
			double lh = 0.0;
			int *this_boot_sample = boot_samples + (boot*nptn);
			for (int ptn = 0; ptn < nptn; ptn++)
				lh += pattern_lh[ptn] * this_boot_sample[ptn];
			tree_lhs_offset[boot] = lh;
		}
		tid++;
	}

	assert(tid == ntrees);

	if (params.topotest_replicates && ntrees > 1) {
		double *tree_probs = new double[ntrees];
		memset(tree_probs, 0, ntrees*sizeof(double));
		int *tree_ranks = new int[ntrees];

		/* perform RELL BP method */
		cout << "Performing RELL-BP test..." << endl;
		int *maxtid = new int[params.topotest_replicates];
		double *maxL = new double[params.topotest_replicates];
		int *maxcount = new int[params.topotest_replicates];
		memset(maxtid, 0, params.topotest_replicates*sizeof(int));
		memcpy(maxL, tree_lhs, params.topotest_replicates*sizeof(double));
		for (boot = 0; boot < params.topotest_replicates; boot++)
			maxcount[boot] = 1;
		for (tid = 1; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				if (tree_lhs_offset[boot] > maxL[boot] + params.ufboot_epsilon) {
					maxL[boot] = tree_lhs_offset[boot];
					maxtid[boot] = tid;
					maxcount[boot] = 1;
				} else if (tree_lhs_offset[boot] > maxL[boot] - params.ufboot_epsilon &&
						random_double() <= 1.0/(maxcount[boot]+1)) {
					maxL[boot] = max(maxL[boot],tree_lhs_offset[boot]);
					maxtid[boot] = tid;
					maxcount[boot]++;
				}
		}
		for (boot = 0; boot < params.topotest_replicates; boot++)
			tree_probs[maxtid[boot]] += 1.0;
		for (tid = 0; tid < ntrees; tid++) {
			tree_probs[tid] /= params.topotest_replicates;
			info[tid].rell_confident = false;
			info[tid].rell_bp = tree_probs[tid];
		}
		sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
		double prob_sum = 0.0;
		// obtain the confidence set
		for (tid = ntrees-1; tid >= 0; tid--) {
			info[tree_ranks[tid]].rell_confident = true;
			prob_sum += tree_probs[tree_ranks[tid]];
			if (prob_sum > 0.95) break;
		}

		// sanity check
		for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
			prob_sum += tree_probs[tid];
		if (fabs(prob_sum-1.0) > 0.01)
			outError("Internal error: Wrong ", __func__);

		delete [] maxcount;
		delete [] maxL;
		delete [] maxtid;

		/* now do the SH test */
		cout << "Performing KH and SH test..." << endl;
		// SH centering step
		for (boot = 0; boot < params.topotest_replicates; boot++)
			max_lh[boot] = -DBL_MAX;
		double *avg_lh = new double[ntrees];
		for (tid = 0; tid < ntrees; tid++) {
			avg_lh[tid] = 0.0;
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				avg_lh[tid] += tree_lhs_offset[boot];
			avg_lh[tid] /= params.topotest_replicates;
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot] - avg_lh[tid]);
			}
		}

		double orig_max_lh = orig_tree_lh[0];
		int orig_max_id = 0;
		double orig_2ndmax_lh = -DBL_MAX;
		int orig_2ndmax_id = -1;
		// find the max tree ID
		for (tid = 1; tid < ntrees; tid++)
			if (orig_max_lh < orig_tree_lh[tid]) {
				orig_max_lh = orig_tree_lh[tid];
				orig_max_id = tid;
			}
		// find the 2nd max tree ID
		for (tid = 0; tid < ntrees; tid++)
			if (tid != orig_max_id && orig_2ndmax_lh < orig_tree_lh[tid]) {
				orig_2ndmax_lh = orig_tree_lh[tid];
				orig_2ndmax_id = tid;
			}


		// SH compute p-value
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			// SH compute original deviation from max_lh
			info[tid].kh_pvalue = 0.0;
			info[tid].sh_pvalue = 0.0;
			int max_id = (tid != orig_max_id) ? orig_max_id : orig_2ndmax_id;
			double orig_diff = orig_tree_lh[max_id] - orig_tree_lh[tid] - avg_lh[tid];
			double *max_kh = tree_lhs + (max_id * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				if (max_lh[boot] - tree_lhs_offset[boot] > orig_diff)
					info[tid].sh_pvalue += 1.0;
				//double max_kh_here = max(max_kh[boot]-avg_lh[max_id], tree_lhs_offset[boot]-avg_lh[tid]);
				double max_kh_here = (max_kh[boot]-avg_lh[max_id]);
				if (max_kh_here - tree_lhs_offset[boot] > orig_diff)
					info[tid].kh_pvalue += 1.0;
			}
			info[tid].sh_pvalue /= params.topotest_replicates;
			info[tid].kh_pvalue /= params.topotest_replicates;
		}

		if (params.do_weighted_test) {

			cout << "Computing pairwise logl difference variance ..." << endl;
			/* computing lhdiff_weights as 1/sqrt(lhdiff_variance) */
			for (tid = 0; tid < ntrees; tid++) {
				double *pattern_lh1 = pattern_lhs + (tid * maxnptn);
				lhdiff_weights[tid*ntrees+tid] = 0.0;
				for (tid2 = tid+1; tid2 < ntrees; tid2++) {
					double lhdiff_variance = tree->computeLogLDiffVariance(pattern_lh1, pattern_lhs + (tid2*maxnptn));
					lhdiff_weights[tid*ntrees+tid2] = 1.0/sqrt(lhdiff_variance);
					lhdiff_weights[tid2*ntrees+tid] = lhdiff_weights[tid*ntrees+tid2];
				}
			}

			// Weighted KH and SH test
			cout << "Performing WKH and WSH test..." << endl;
			for (tid = 0; tid < ntrees; tid++) {
				double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
				info[tid].wkh_pvalue = 0.0;
				info[tid].wsh_pvalue = 0.0;
				double worig_diff = -DBL_MAX;
				int max_id = -1;
				for (tid2 = 0; tid2 < ntrees; tid2++)
					if (tid2 != tid) {
						double wdiff = (orig_tree_lh[tid2] - orig_tree_lh[tid])*lhdiff_weights[tid*ntrees+tid2];
						if (wdiff > worig_diff) {
							worig_diff = wdiff;
							max_id = tid2;
						}
					}
				for (boot = 0; boot < params.topotest_replicates; boot++) {
					double wmax_diff = -DBL_MAX;
					for (tid2 = 0; tid2 < ntrees; tid2++)
						if (tid2 != tid)
							wmax_diff = max(wmax_diff,
									(tree_lhs[tid2*params.topotest_replicates+boot] - avg_lh[tid2] -
									tree_lhs_offset[boot] + avg_lh[tid]) * lhdiff_weights[tid*ntrees+tid2]);
					if (wmax_diff > worig_diff)
						info[tid].wsh_pvalue += 1.0;
					wmax_diff = (tree_lhs[max_id*params.topotest_replicates+boot] - avg_lh[max_id] -
							tree_lhs_offset[boot] + avg_lh[tid]);
					if (wmax_diff >  orig_tree_lh[max_id] - orig_tree_lh[tid])
						info[tid].wkh_pvalue += 1.0;
				}
				info[tid].wsh_pvalue /= params.topotest_replicates;
				info[tid].wkh_pvalue /= params.topotest_replicates;
			}
		}
        
        delete [] avg_lh;
        
		/* now to ELW - Expected Likelihood Weight method */
		cout << "Performing ELW test..." << endl;

		for (boot = 0; boot < params.topotest_replicates; boot++)
			max_lh[boot] = -DBL_MAX;
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot]);
		}
		double *sumL = new double[params.topotest_replicates];
		memset(sumL, 0, sizeof(double) * params.topotest_replicates);
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				tree_lhs_offset[boot] = exp(tree_lhs_offset[boot] - max_lh[boot]);
				sumL[boot] += tree_lhs_offset[boot];
			}
		}
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			tree_probs[tid] = 0.0;
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				tree_probs[tid] += (tree_lhs_offset[boot] / sumL[boot]);
			}
			tree_probs[tid] /= params.topotest_replicates;
			info[tid].elw_confident = false;
			info[tid].elw_value = tree_probs[tid];
		}

		sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
		prob_sum = 0.0;
		// obtain the confidence set
		for (tid = ntrees-1; tid >= 0; tid--) {
			info[tree_ranks[tid]].elw_confident = true;
			prob_sum += tree_probs[tree_ranks[tid]];
			if (prob_sum > 0.95) break;
		}

		// sanity check
		for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
			prob_sum += tree_probs[tid];
		if (fabs(prob_sum-1.0) > 0.01)
			outError("Internal error: Wrong ", __func__);
		delete [] sumL;

        if (params.do_au_test) {
            cout << "Performing approximately unbiased (AU) test..." << endl;
            performAUTest(params, tree, pattern_lhs, info);
        }

		delete [] tree_ranks;
		delete [] tree_probs;

	}
	if (max_lh)
		delete [] max_lh;
	if (orig_tree_lh)
		delete [] orig_tree_lh;
	if (pattern_lh)
        aligned_free(pattern_lh);
	if (pattern_lhs)
        aligned_free(pattern_lhs);
	if (lhdiff_weights)
		delete [] lhdiff_weights;
	if (tree_lhs)
		delete [] tree_lhs;
	//if (saved_tree_lhs)
	//	delete [] saved_tree_lhs;
	if (boot_samples)
		delete [] boot_samples;

	if (params.print_tree_lh) {
		scoreout.close();
	}

	treeout.close();
	in.close();

	cout << "Time for evaluating all trees: " << getRealTime() - time_start << " sec." << endl;

}


void evaluateTrees(Params &params, IQTree *tree) {
	vector<TreeInfo> info;
	IntVector distinct_ids;
	evaluateTrees(params, tree, info, distinct_ids);
}



