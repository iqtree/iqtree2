/*
 * phylotesting.cpp
 * implementation of ModelFinder and PartitionFinder
 *  Created on: Aug 23, 2013
 *      Author: minh
 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include <numeric>
#include "tree/phylotree.h"
#include "tree/iqtree.h"
#include "tree/phylosupertree.h"
#include "tree/phylotreemixlen.h"
#include "phylotesting.h"

#include "model/modelmarkov.h"
#include "model/modeldna.h"
#include "nclextra/myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
#include "model/ratefree.h"
#include "model/ratefreeinvar.h"
//#include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "model/modelmorphology.h"
#include "model/modelmixture.h"
#include "model/modelliemarkov.h"
#include "model/modelpomo.h"
#include "utils/timeutil.h"
#include "model/modelfactorymixlen.h"
#include "tree/phylosupertreeplen.h"
#include "tree/phylosupertreeunlinked.h"

#include "phyloanalysis.h"
#include "gsl/mygsl.h"
#include "utils/MPIHelper.h"
//#include "vectorclass/vectorclass.h"


/******* Binary model set ******/
const char* bin_model_names[] = {"GTR2", "JC2"};


/******* Morphological model set ******/
// 2018-08-20: don't test ORDERED model due to lots of numerical issues
//const char* morph_model_names[] = {"MK", "ORDERED"};
const char* morph_model_names[] = {"MK"};


/******* DNA model set ******/
const char* dna_model_names[] = {"GTR", "SYM", "TVM",  "TVMe", "TIM3",
        "TIM3e", "TIM2", "TIM2e", "TIM", "TIMe", "TPM3u", "TPM3",
        "TPM2u",  "TPM2",  "K81u", "K81", "TN", "TNe",  "HKY",  "K80", "F81", "JC"};

/* DNA models supported by PhyML/PartitionFinder */
const char* dna_model_names_old[] ={"GTR",  "SYM", "TVM", "TVMe", "TIM", "TIMe",
         "K81u", "K81", "TN", "TNe", "HKY", "K80", "F81", "JC"};

/* DNA model supported by RAxML */
const char* dna_model_names_rax[] ={"GTR"};

/* DNA model supported by MrBayes */
const char *dna_model_names_mrbayes[] = {"GTR", "SYM", "HKY", "K80", "F81", "JC"};

/* DNA model supported by BEAST1 */
const char *dna_model_names_beast1[] = {"GTR", "TN", "HKY"};

/* DNA model supported by BEAST2 */
const char *dna_model_names_beast2[] = {"GTR", "TN", "HKY", "JC"};

/* DNA model supported by ModelOMatic */
const char *dna_model_names_modelomatic[] = {"GTR", "HKY", "K80", "F81", "JC"};

//const char* dna_freq_names[] = {"+FO"};

// Lie-Markov models without an RY, WS or MK prefix
const char *dna_model_names_lie_markov_fullsym[] = 
    {"1.1", "3.3a", "4.4a", "6.7a", "9.20b", "12.12"};
// Lie-Markov models with RY symmetry/distinguished pairing
const char *dna_model_names_lie_markov_ry[] = {
          "RY2.2b",  "RY3.3b", "RY3.3c",  "RY3.4",   "RY4.4b", 
          "RY4.5a",  "RY4.5b", "RY5.6a",  "RY5.6b",  "RY5.7a",
          "RY5.7b",  "RY5.7c", "RY5.11a", "RY5.11b", "RY5.11c", 
          "RY5.16",  "RY6.6",  "RY6.7b",  "RY6.8a",  "RY6.8b",  
          "RY6.17a", "RY6.17b","RY8.8",   "RY8.10a", "RY8.10b", 
          "RY8.16",  "RY8.17", "RY8.18",  "RY9.20a", "RY10.12",
	  "RY10.34"
};
// Lie-Markov models with WS symmetry/distinguished pairing
const char *dna_model_names_lie_markov_ws[] = {
          "WS2.2b",  "WS3.3b", "WS3.3c",  "WS3.4",   "WS4.4b", 
          "WS4.5a",  "WS4.5b", "WS5.6a",  "WS5.6b",  "WS5.7a",
          "WS5.7b",  "WS5.7c", "WS5.11a", "WS5.11b", "WS5.11c", 
          "WS5.16",  "WS6.6",  "WS6.7b",  "WS6.8a",  "WS6.8b",  
          "WS6.17a", "WS6.17b","WS8.8",   "WS8.10a", "WS8.10b", 
          "WS8.16",  "WS8.17", "WS8.18",  "WS9.20a", "WS10.12",
	  "WS10.34"
};
// Lie-Markov models with MK symmetry/distinguished pairing
const char *dna_model_names_lie_markov_mk[] = {
          "MK2.2b",  "MK3.3b", "MK3.3c",  "MK3.4",   "MK4.4b", 
          "MK4.5a",  "MK4.5b", "MK5.6a",  "MK5.6b",  "MK5.7a",
          "MK5.7b",  "MK5.7c", "MK5.11a", "MK5.11b", "MK5.11c", 
          "MK5.16",  "MK6.6",  "MK6.7b",  "MK6.8a",  "MK6.8b",  
          "MK6.17a", "MK6.17b","MK8.8",   "MK8.10a", "MK8.10b", 
          "MK8.16",  "MK8.17", "MK8.18",  "MK9.20a", "MK10.12",
	  "MK10.34"
};
// Lie-Markov models which are strand symmetric
const char *dna_model_names_lie_markov_strsym[] = {
          "1.1",    "WS2.2b", "3.3a",   "WS3.3b", "WS3.3c", "WS3.4",   
          "WS4.4b", "WS4.5a", "WS4.5b", "WS5.6a", "WS6.6"
};


/****** Protein model set ******/
const char* aa_model_names[] = {"LG", "WAG", "JTT", "Q.pfam", "Q.bird", "Q.mammal", "Q.insect", "Q.plant", "Q.yeast", "JTTDCMut", "DCMut", "VT", "PMB", "Blosum62", "Dayhoff",
        "mtREV", "mtART", "mtZOA", "mtMet" , "mtVer" , "mtInv", "mtMAM", "FLAVI",
		"HIVb", "HIVw", "FLU", "rtREV", "cpREV"};

/****** Protein mixture model set ******/
const char* aa_mixture_model_names[] = {"C10", "C20", "C30", "C40", "C50", "C60", "EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X", "CF4"};

/* Protein models supported by PhyML/PartitionFinder */
const char *aa_model_names_phyml[] = {"LG", "WAG", "JTT", "DCMut", "VT", "Blosum62", "Dayhoff",
		"mtREV", "mtART", "mtMAM",
		"HIVb", "HIVw", "rtREV", "cpREV"};

/* Protein models supported by RAxML */
const char *aa_model_names_rax[] = {"LG", "WAG", "JTT", "JTTDCMut", "DCMut", "VT", "PMB", "Blosum62", "Dayhoff",
        "mtREV", "mtART", "mtZOA", "mtMAM",
		"HIVb", "HIVw", "FLU", "rtREV", "cpREV"};

const char* aa_model_names_mrbayes[] = {"WAG", "JTT", "VT", "Blosum62", "Dayhoff",
        "mtREV", "mtMAM",
		"rtREV", "cpREV"};

const char* aa_model_names_beast1[] = {"LG", "WAG", "JTT", "Blosum62", "Dayhoff", "mtREV", "cpREV", "FLU"};

const char* aa_model_names_beast2[] = {"LG", "WAG", "JTT", "DCMut", "VT", "Blosum62", "Dayhoff",
    "mtREV", "mtART", "mtMAM", "HIVb", "HIVw", "FLU", "rtREV", "cpREV"};

const char* aa_model_names_modelomatic[] = {"LG", "WAG", "JTT", "VT", "Blosum62", "Dayhoff",
        "mtART", "mtMAM", "mtREV",
        "HIVb", "HIVw", "rtREV", "cpREV"};

const char *aa_model_names_nuclear[] = {"LG", "WAG", "JTT", "Q.pfam", "JTTDCMut","DCMut", "VT", "PMB", "Blosum62", "Dayhoff"};

const char *aa_model_names_mitochondrial[] = {"mtREV", "mtART", "mtZOA", "mtMet" , "mtVer" , "mtInv", "mtMAM"};

const char *aa_model_names_chloroplast[] = {"cpREV"};

const char *aa_model_names_viral[] = {"HIVb", "HIVw", "FLU", "rtREV", "FLAVI"};

const char* aa_freq_names[] = {"", "+F"};


/****** Codon models ******/
//const char *codon_model_names[] = {"GY", "MG", "MGK", "KOSI07", "SCHN05","KOSI07_GY1KTV","SCHN05_GY1KTV"};
//short int std_genetic_code[]    = {   0,    0,     0,        1,        1,              1,              1};
const char *codon_model_names[] = { "GY", "MGK", "MG", "KOSI07", "SCHN05"};
short int std_genetic_code[]    = {   0,    0,     0,        1,        1};
const char *codon_model_names_modelomatic[] = {"GY"};
short int std_genetic_code_modelomatic[]    = {   0};

const char *codon_freq_names[] = {"+F3X4", "+F1X4", "+F", ""};

//const double TOL_LIKELIHOOD_MODELTEST = 0.1;
const double TOL_GRADIENT_MODELTEST   = 0.0001;

extern double RunKMeans1D(int n, int k, double *points, int *weights, double *centers, int *assignments);


string getSeqTypeName(SeqType seq_type) {
    switch (seq_type) {
        case SEQ_BINARY: return "binary";
        case SEQ_DNA: return "DNA";
        case SEQ_PROTEIN: return "protein";
        case SEQ_CODON: return "codon";
        case SEQ_MORPH: return "morphological";
        case SEQ_POMO: return "PoMo";
        case SEQ_UNKNOWN: return "unknown";
        case SEQ_MULTISTATE: return "MultiState";
    }
}

string getUsualModelSubst(SeqType seq_type) {
    switch (seq_type) {
        case SEQ_DNA: return dna_model_names[0];
        case SEQ_PROTEIN: return aa_model_names[0];
        case SEQ_CODON: return string(codon_model_names[0]) + codon_freq_names[0];
        case SEQ_BINARY: return bin_model_names[0];
        case SEQ_MORPH: return morph_model_names[0];
        case SEQ_POMO: return string(dna_model_names[0]) + "+P";
        default: ASSERT(0 && "Unprocessed seq_type"); return "";
    }
}

void getRateHet(SeqType seq_type, string model_name, double frac_invariant_sites,
                string rate_set, StrVector &ratehet);

size_t CandidateModel::getUsualModel(Alignment *aln) {
    size_t aln_len = 0;
    if (aln->isSuperAlignment()) {
        SuperAlignment *super_aln = (SuperAlignment*)aln;
        for (auto it = super_aln->partitions.begin(); it != super_aln->partitions.end(); it++) {
            CandidateModel usual_model(*it);
            if (!subst_name.empty())
                subst_name += ',';
            subst_name += usual_model.subst_name;
            if (!rate_name.empty())
                rate_name += ',';
            rate_name += usual_model.rate_name;
            aln_len += (*it)->getNSite();
        }
    } else {
        subst_name = getUsualModelSubst(aln->seq_type);
        StrVector ratehet;
        getRateHet(aln->seq_type, Params::getInstance().model_name, aln->frac_invariant_sites, "1", ratehet);
        ASSERT(!ratehet.empty());
        rate_name = ratehet[0];
        aln_len = aln->getNSite();
    }
    orig_subst_name = subst_name;
    orig_rate_name = rate_name;
    return aln_len;
}

void CandidateModel::computeICScores(size_t sample_size) {
    computeInformationScores(logl, df, sample_size, AIC_score, AICc_score, BIC_score);
}

void CandidateModel::computeICScores() {
    size_t sample_size = aln->getNSite();
    if (aln->isSuperAlignment()) {
        sample_size = 0;
        SuperAlignment *super_aln = (SuperAlignment*)aln;
        for (auto a : super_aln->partitions)
            sample_size += a->getNSite();
    }
    if (hasFlag(MF_SAMPLE_SIZE_TRIPLE))
        sample_size /= 3;
    computeInformationScores(logl, df, sample_size, AIC_score, AICc_score, BIC_score);
}

double CandidateModel::computeICScore(size_t sample_size) {
    return computeInformationScore(logl, df, sample_size, Params::getInstance().model_test_criterion);
}

double CandidateModel::getScore(ModelTestCriterion mtc) {
    switch (mtc) {
        case MTC_AIC:
            return AIC_score;
        case MTC_AICC:
            return AICc_score;
        case MTC_BIC:
            return BIC_score;
        case MTC_ALL:
            ASSERT(0 && "Unhandled case");
            return 0.0;
    }
}

double CandidateModel::getScore() {
    return getScore(Params::getInstance().model_test_criterion);
}

int CandidateModelSet::getBestModelID(ModelTestCriterion mtc) {
    double best_score = DBL_MAX;
    int best_model = -1;
    for (int model = 0; model < size(); model++)
        if (at(model).hasFlag(MF_DONE) && best_score > at(model).getScore(mtc)) {
            best_score = at(model).getScore(mtc);
            best_model = model;
        }
    return best_model;
}

bool ModelCheckpoint::getBestModel(string &best_model) {
    return getString("best_model_" + criterionName(Params::getInstance().model_test_criterion), best_model);
}

bool ModelCheckpoint::getBestModelList(string &best_model_list) {
    return getString("best_model_list_" + criterionName(Params::getInstance().model_test_criterion), best_model_list);
}

void ModelCheckpoint::putBestModelList(string &best_model_list) {
    return put("best_model_list_" + criterionName(Params::getInstance().model_test_criterion), best_model_list);
}

bool  ModelCheckpoint::getBestTree(string &best_tree) {
    return getString("best_tree_" + criterionName(Params::getInstance().model_test_criterion), best_tree);
}

bool ModelCheckpoint::getOrderedModels(PhyloTree *tree, CandidateModelSet &ordered_models) {
    double best_score_AIC, best_score_AICc, best_score_BIC;
    if (tree->isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*)tree;
        ordered_models.clear();
        for (int part = 0; part != stree->size(); part++) {
            startStruct(stree->at(part)->aln->name);
            CandidateModel info;
            if (!getBestModel(info.subst_name)) return false;
            info.restoreCheckpoint(this);
            info.computeICScores(stree->at(part)->getAlnNSite());
            endStruct();
            ordered_models.push_back(info);
        }
        return true;
    } else {
        CKP_RESTORE2(this, best_score_AIC);
        CKP_RESTORE2(this, best_score_AICc);
        CKP_RESTORE2(this, best_score_BIC);
        double sum_AIC = 0, sum_AICc = 0, sum_BIC = 0;
        string str;
        bool ret = getBestModelList(str);
        if (!ret) return false;
        istringstream istr(str);
        string model;
        ordered_models.clear();
        while (istr >> model) {
            CandidateModel info;
            info.subst_name = model;
            info.restoreCheckpoint(this);
            info.computeICScores(tree->getAlnNSite());
            sum_AIC  += info.AIC_weight = exp(-0.5*(info.AIC_score-best_score_AIC));
            sum_AICc += info.AICc_weight = exp(-0.5*(info.AICc_score-best_score_AICc));
            sum_BIC  += info.BIC_weight = exp(-0.5*(info.BIC_score-best_score_BIC));
            ordered_models.push_back(info);
        }
        sum_AIC = 1.0/sum_AIC;
        sum_AICc = 1.0/sum_AICc;
        sum_BIC = 1.0/sum_BIC;
        for (auto it = ordered_models.begin(); it != ordered_models.end(); it++) {
            it->AIC_weight *= sum_AIC;
            it->AICc_weight *= sum_AICc;
            it->BIC_weight *= sum_BIC;
            it->AIC_conf = it->AIC_weight > 0.05;
            it->AICc_conf = it->AICc_weight > 0.05;
            it->BIC_conf = it->BIC_weight > 0.05;
        }
        return true;
    }
}

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

/**
 * append from cvec to strvec
 */
void appendCString(const char **cvec, int n, StrVector &strvec, bool touppercase = false) {
        strvec.reserve(strvec.size()+n);
	for (int i = 0; i < n; i++) {
	    strvec.push_back(cvec[i]);
            if (touppercase)
	      std::transform(strvec.back().begin(), strvec.back().end(), strvec.back().begin(), ::toupper);
        }
}


int detectSeqType(const char *model_name, SeqType &seq_type) {
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
    copyCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_DNA;
            break;
        }
    copyCString(dna_model_names_lie_markov_ry, sizeof(dna_model_names_lie_markov_ry)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_DNA;
            break;
        }
    copyCString(dna_model_names_lie_markov_ws, sizeof(dna_model_names_lie_markov_ws)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_DNA;
            break;
        }
    copyCString(dna_model_names_lie_markov_mk, sizeof(dna_model_names_lie_markov_mk)/sizeof(char*), model_list, true);
    for (i = 0; i < model_list.size(); i++)
        if (model_str == model_list[i]) {
            seq_type = SEQ_DNA;
            break;
        }
    copyCString(dna_model_names_lie_markov_strsym, sizeof(dna_model_names_lie_markov_strsym)/sizeof(char*), model_list, true);
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
    copyCString(aa_mixture_model_names, sizeof(aa_mixture_model_names)/sizeof(char*), model_list, true);
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

string convertSeqTypeToSeqTypeName(SeqType seq_type)
{
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

string detectSeqTypeName(string model_name) {
    SeqType seq_type;
    detectSeqType(model_name.c_str(), seq_type);
    return convertSeqTypeToSeqTypeName(seq_type);
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


/**
 * select models for all partitions
 * @param[in,out] model_info (IN/OUT) all model information
 * @return total number of parameters
 */
void testPartitionModel(Params &params, PhyloSuperTree* in_tree, ModelCheckpoint &model_info,
                        ModelsBlock *models_block, int num_threads);


/**
 compute log-adapter function according to Whelan et al. 2015
 @param orig_aln original codon alignment
 @param newaln AA alignment
 @param[out] adjusted_df adjusted degree of freedom factor
 @return adjusted log-likelihood factor
 */
double computeAdapter(Alignment *orig_aln, Alignment *newaln, int &adjusted_df) {
    int aa, codon;

    // count codon occurences
    unsigned int codon_counts[orig_aln->num_states];
    orig_aln->computeAbsoluteStateFreq(codon_counts);
    
    // compute AA frequency
//    double aa_freq[newaln->num_states];
//    newaln->computeStateFreq(aa_freq);
    
    // compute codon frequency
    double codon_freq[orig_aln->num_states];
    //orig_aln->computeStateFreq(codon_freq);
    
    double sum = 0.0;
    for (codon = 0; codon < orig_aln->num_states; codon++)
        sum += codon_counts[codon];
    sum = 1.0/sum;
    for (codon = 0; codon < orig_aln->num_states; codon++)
        codon_freq[codon] = sum*codon_counts[codon];
    
    // new rescale codon_freq s.t. codons coding for the same AA
    // have f summing up to the frequency of this AA
    for (aa = 0; aa < newaln->num_states; aa++) {
        double sum = 0;
        for (codon = 0; codon < orig_aln->num_states; codon++)
            if (newaln->convertState(orig_aln->genetic_code[(int)orig_aln->codon_table[codon]]) == aa)
                sum += codon_freq[codon];
        sum = 1.0/sum;
        for (codon = 0; codon < orig_aln->num_states; codon++)
            if (newaln->convertState(orig_aln->genetic_code[(int)orig_aln->codon_table[codon]]) == aa)
                codon_freq[codon] *= sum;
    }
    
    // now compute adapter function
    double adapter = 0.0;
    adjusted_df = 0;
    vector<bool> has_AA;
    has_AA.resize(newaln->num_states, false);
    
    for (codon = 0; codon < orig_aln->num_states; codon++) {
        if (codon_counts[codon] == 0)
            continue;
        has_AA[newaln->convertState(orig_aln->genetic_code[(int)orig_aln->codon_table[codon]])] = true;
        adapter += codon_counts[codon]*log(codon_freq[codon]);
        adjusted_df++;
    }
    for (aa = 0; aa < has_AA.size(); aa++)
        if (has_AA[aa])
            adjusted_df--;
    return adapter;
}

/**
 compute fast ML tree by stepwise addition MP + ML-NNI
 @return the tree string
 */
string computeFastMLTree(Params &params, Alignment *aln,
                       ModelCheckpoint &model_info, ModelsBlock *models_block,
                       int &num_threads, int brlen_type, string dist_file) {
    //string model_name;
    CandidateModel usual_model(aln);
    StrVector subst_names;
    StrVector rate_names;
    convert_string_vec(usual_model.subst_name.c_str(), subst_names);
    convert_string_vec(usual_model.rate_name.c_str(), rate_names);
    ASSERT(subst_names.size() == rate_names.size());
    //set<string> model_set;

    string concat_tree;
    
    IQTree *iqtree = NULL;
    
    StrVector saved_model_names;
    
    if (aln->isSuperAlignment()) {
        SuperAlignment *saln = (SuperAlignment*)aln;
        if (params.partition_type == TOPO_UNLINKED)
            iqtree = new PhyloSuperTreeUnlinked(saln);
        else if (params.partition_type == BRLEN_OPTIMIZE)
            iqtree = new PhyloSuperTree(saln);
        else
            iqtree = new PhyloSuperTreePlen(saln, brlen_type);
        for (int part = 0; part != subst_names.size(); part++) {
            saved_model_names.push_back(saln->partitions[part]->model_name);
            saln->partitions[part]->model_name = subst_names[part] + rate_names[part];
        }
    } else if (posRateHeterotachy(rate_names[0]) != string::npos) {
        iqtree = new PhyloTreeMixlen(aln, 0);
    } else {
        iqtree = new IQTree(aln);
    }
    if ((params.start_tree == STT_PLL_PARSIMONY || params.start_tree == STT_RANDOM_TREE || params.pll) && !iqtree->isInitializedPLL()) {
        /* Initialized all data structure for PLL*/
        iqtree->initializePLL(params);
    }
    iqtree->setParams(&params);
    iqtree->setLikelihoodKernel(params.SSE);
    iqtree->optimize_by_newton = params.optimize_by_newton;
    iqtree->setNumThreads(num_threads);
    iqtree->setCheckpoint(&model_info);

    iqtree->dist_file = dist_file;
    iqtree->computeInitialTree(params.SSE);
    iqtree->restoreCheckpoint();
    
    //ASSERT(iqtree->root);
    iqtree->initializeModel(params, usual_model.getName(), models_block);
    if (!iqtree->getModel()->isMixture() || aln->seq_type == SEQ_POMO) {
        usual_model.subst_name = iqtree->getSubstName();
        usual_model.rate_name = iqtree->getRateName();
    }

    iqtree->getModelFactory()->restoreCheckpoint();
    iqtree->ensureNumberOfThreadsIsSet(&params);
    iqtree->initializeAllPartialLh();
    double saved_modelEps = params.modelEps;
    params.modelEps = params.modelfinder_eps;
    string initTree;
    
    double start_time = getRealTime();
    
    cout << "Perform fast likelihood tree search using " << subst_names[0]+rate_names[0] << " model..." << endl;
    
    if (iqtree->getCheckpoint()->getBool("finishedFastMLTree")) {
        // model optimization already done: ignore this step
        iqtree->setCurScore(iqtree->computeLikelihood());
        initTree = iqtree->getTreeString();
        cout << "CHECKPOINT: Tree restored, LogL: " << iqtree->getCurScore() << endl;
    } else {
        bool saved_opt_gammai = params.opt_gammai;
        // disable thorough I+G optimization
        params.opt_gammai = false;
        initTree = iqtree->optimizeModelParameters(false, params.modelEps*50.0);
        if (iqtree->isMixlen())
            initTree = ((ModelFactoryMixlen*)iqtree->getModelFactory())->sortClassesByTreeLength();

        // do quick NNI search
        if (params.start_tree != STT_USER_TREE) {
            cout << "Perform nearest neighbor interchange..." << endl;
            iqtree->doNNISearch(true);
            initTree = iqtree->getTreeString();
        }
        params.opt_gammai = saved_opt_gammai;

        iqtree->saveCheckpoint();
        iqtree->getModelFactory()->saveCheckpoint();
        iqtree->getCheckpoint()->putBool("finishedFastMLTree", true);
        iqtree->getCheckpoint()->dump();
        //        cout << "initTree: " << initTree << endl;
        cout << "Time for fast ML tree search: " << getRealTime() - start_time << " seconds" << endl;
        cout << endl;
    }

    // restore model epsilon
    params.modelEps = saved_modelEps;
    
    // save information to the checkpoint for later retrieval
    if (iqtree->isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*)iqtree;
        int part = 0;
        for (auto it = stree->begin(); it != stree->end(); it++, part++) {
            model_info.startStruct((*it)->aln->name);
            (*it)->saveCheckpoint();
            (*it)->getModelFactory()->saveCheckpoint();
            model_info.endStruct();
        }
        SuperAlignment *saln = (SuperAlignment*)aln;
        // restore model_names
        for (int i = 0; i < saln->partitions.size(); i++)
            saln->partitions[i]->model_name = saved_model_names[i];
    } else {
        iqtree->saveCheckpoint();
        iqtree->getModelFactory()->saveCheckpoint();
    }


    delete iqtree;
    return initTree;
}

/**
 Transfer parameters from ModelFinder into the a checkpoint to speed up later stage
 */
void transferModelFinderParameters(IQTree *iqtree, Checkpoint *target) {
    
    Checkpoint *source = iqtree->getCheckpoint();
    
    // transfer the substitution model and site-rate parameters
    if (iqtree->isSuperTree()) {
        DoubleVector tree_lens;
        string struct_name;
        if (iqtree->params->partition_type == BRLEN_SCALE || iqtree->params->partition_type == BRLEN_FIX)
            struct_name = "PartitionModelPlen";
        else
            struct_name = "PartitionModel";
        target->startStruct(struct_name);
        SuperAlignment *super_aln = (SuperAlignment*)iqtree->aln;
        for (auto aln : super_aln->partitions) {
            source->transferSubCheckpoint(target, aln->name + CKP_SEP + "Model");
            source->transferSubCheckpoint(target, aln->name + CKP_SEP + "Rate");

            // transfer partition rates
            if (iqtree->params->partition_type == BRLEN_SCALE) {
                source->startStruct(aln->name);
                CandidateModel info;
                info.subst_name = aln->model_name;
                if (info.restoreCheckpoint(source))
                    tree_lens.push_back(info.tree_len);
                else
                    ASSERT(0 && "Could not restore tree_len");
                source->endStruct();
            }
        }

        if (iqtree->params->partition_type == BRLEN_SCALE) {
            // now normalize the rates
            PhyloSuperTree *tree = (PhyloSuperTree*)iqtree;
            double sum = 0.0;
            size_t nsite = 0;
            int i;
            for (i = 0; i < tree->size(); i++) {
                sum += tree_lens[i] * tree->at(i)->aln->getNSite();
                if (tree->at(i)->aln->seq_type == SEQ_CODON && tree->rescale_codon_brlen)
                    nsite += 3*tree->at(i)->aln->getNSite();
                else
                    nsite += tree->at(i)->aln->getNSite();
            }

            sum /= nsite;
            iqtree->restoreCheckpoint();
            iqtree->scaleLength(sum/iqtree->treeLength());
            iqtree->saveCheckpoint();
            sum = 1.0/sum;
            for (i = 0; i < tree->size(); i++)
                tree_lens[i] *= sum;
            target->putVector("part_rates", tree_lens);            
        }
        target->endStruct();
    } else {
        source->transferSubCheckpoint(target, "Model");
        source->transferSubCheckpoint(target, "Rate");
    }
    
    // transfer tree
    source->transferSubCheckpoint(target, "PhyloTree");
}

void runModelFinder(Params &params, IQTree &iqtree, ModelCheckpoint &model_info)
{
    if (params.model_name.find("+T") != string::npos) {
        // tree mixture
        return;
    }
    
    //    iqtree.setCurScore(-DBL_MAX);
    bool test_only = (params.model_name.find("ONLY") != string::npos) ||
        (params.model_name.substr(0,2) == "MF" && params.model_name.substr(0,3) != "MFP");
    
    bool empty_model_found = params.model_name.empty() && !iqtree.isSuperTree();
    
    if (params.model_name.empty() && iqtree.isSuperTree()) {
        // check whether any partition has empty model_name
        PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
        for (auto i = stree->begin(); i != stree->end(); i++)
            if ((*i)->aln->model_name.empty()) {
                empty_model_found = true;
                break;
            }
    }
    
    if (params.model_joint)
        empty_model_found = false;

    // Model already specifed, nothing to do here
    if (!empty_model_found && params.model_name.substr(0, 4) != "TEST" && params.model_name.substr(0, 2) != "MF")
        return;
    if (MPIHelper::getInstance().getNumProcesses() > 1)
        outError("Please use only 1 MPI process! We are currently working on the MPI parallelization of model selection.");
    // TODO: check if necessary
    //        if (iqtree.isSuperTree())
    //            ((PhyloSuperTree*) &iqtree)->mapTrees();
    double cpu_time = getCPUTime();
    double real_time = getRealTime();
    model_info.setFileName((string)params.out_prefix + ".model.gz");
    model_info.setDumpInterval(params.checkpoint_dump_interval);
    
    bool ok_model_file = false;
    if (!params.model_test_again) {
        ok_model_file = model_info.load();
    }
    
    cout << endl;
    
    ok_model_file &= model_info.size() > 0;
    if (ok_model_file)
        cout << "NOTE: Restoring information from model checkpoint file " << model_info.getFileName() << endl;
    
    
    Checkpoint *orig_checkpoint = iqtree.getCheckpoint();
    iqtree.setCheckpoint(&model_info);
    iqtree.restoreCheckpoint();
    
    int partition_type;
    if (CKP_RESTORE2((&model_info), partition_type)) {
        if (partition_type != params.partition_type)
            outError("Mismatch partition type between checkpoint and partition file command option\nRerun with -mredo to ignore .model.gz checkpoint file");
    } else {
        partition_type = params.partition_type;
        CKP_SAVE2((&model_info), partition_type);
    }
    
    ModelsBlock *models_block = readModelsDefinition(params);
    
    // compute initial tree
    if (params.modelfinder_ml_tree) {
        // 2019-09-10: Now perform NNI on the initial tree
        string tree_str = computeFastMLTree(params, iqtree.aln, model_info,
            models_block, params.num_threads, params.partition_type, iqtree.dist_file);
        iqtree.restoreCheckpoint();
    } else {
        iqtree.computeInitialTree(params.SSE);
    
        if (iqtree.isSuperTree()) {
            PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
            int part = 0;
            for (auto it = stree->begin(); it != stree->end(); it++, part++) {
                model_info.startStruct((*it)->aln->name);
                (*it)->saveCheckpoint();
                model_info.endStruct();
            }
        } else {
            iqtree.saveCheckpoint();
        }
    }
    
    // also save initial tree to the original .ckp.gz checkpoint
    //        string initTree = iqtree.getTreeString();
    //        CKP_SAVE(initTree);
    //        iqtree.saveCheckpoint();
    //        checkpoint->dump(true);
    
    CandidateModelSet candidate_models;
    int max_cats = candidate_models.generate(params, iqtree.aln, params.model_test_separate_rate, false);
    
    uint64_t mem_size = iqtree.getMemoryRequiredThreaded(max_cats);
    cout << "NOTE: ModelFinder requires " << (mem_size / 1024) / 1024 << " MB RAM!" << endl;
    if (mem_size >= getMemorySize()) {
        outError("Memory required exceeds your computer RAM size!");
    }
#ifdef BINARY32
    if (mem_size >= 2000000000) {
        outError("Memory required exceeds 2GB limit of 32-bit executable");
    }
#endif
    

    if (iqtree.isSuperTree()) {
        // partition model selection
        PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
        testPartitionModel(params, stree, model_info, models_block, params.num_threads);
        stree->mapTrees();
        string res_models = "";
        for (auto it = stree->begin(); it != stree->end(); it++) {
            if (it != stree->begin()) res_models += ",";
            res_models += (*it)->aln->model_name;
        }
        iqtree.aln->model_name = res_models;
    } else {
        // single model selection
        CandidateModel best_model;
        if (params.openmp_by_model)
            best_model = CandidateModelSet().evaluateAll(params, &iqtree,
                model_info, models_block, params.num_threads, BRLEN_OPTIMIZE);
        else
            best_model = CandidateModelSet().test(params, &iqtree,
                model_info, models_block, params.num_threads, BRLEN_OPTIMIZE);
        iqtree.aln->model_name = best_model.getName();
        
        Checkpoint *checkpoint = &model_info;
        string best_model_AIC, best_model_AICc, best_model_BIC;
        CKP_RESTORE(best_model_AIC);
        CKP_RESTORE(best_model_AICc);
        CKP_RESTORE(best_model_BIC);
        cout << "Akaike Information Criterion:           " << best_model_AIC << endl;
        cout << "Corrected Akaike Information Criterion: " << best_model_AICc << endl;
        cout << "Bayesian Information Criterion:         " << best_model_BIC << endl;
        cout << "Best-fit model: " << iqtree.aln->model_name << " chosen according to "
            << criterionName(params.model_test_criterion) << endl;
    }

    delete models_block;
    
    // force to dump all checkpointing information
    model_info.dump(true);
    
    // transfer models parameters
    transferModelFinderParameters(&iqtree, orig_checkpoint);
    iqtree.setCheckpoint(orig_checkpoint);
    
    params.startCPUTime = cpu_time;
    params.start_real_time = real_time;
    cpu_time = getCPUTime() - cpu_time;
    real_time = getRealTime() - real_time;
    cout << endl;
    cout << "All model information printed to " << model_info.getFileName() << endl;
    cout << "CPU time for ModelFinder: " << cpu_time << " seconds (" << convert_time(cpu_time) << ")" << endl;
    cout << "Wall-clock time for ModelFinder: " << real_time << " seconds (" << convert_time(real_time) << ")" << endl;
    
    //        alignment = iqtree.aln;
    if (test_only) {
        params.min_iterations = 0;
    }
}

/**
 * get the list of substitution models
 */
void getModelSubst(SeqType seq_type, bool standard_code, string model_name,
                   string model_set, char *model_subset, StrVector &model_names) {
    int i, j;
    
    if (model_set == "1") {
        model_names.push_back(getUsualModelSubst(seq_type));
        return;
    }
    
    if (iEquals(model_set, "ALL") || iEquals(model_set, "AUTO"))
        model_set = "";
    
    if (seq_type == SEQ_BINARY) {
        if (model_set.empty()) {
            copyCString(bin_model_names, sizeof(bin_model_names) / sizeof(char*), model_names);
        } else if (model_set[0] == '+') {
            // append model_set into existing models
            convert_string_vec(model_set.c_str()+1, model_names);
            appendCString(bin_model_names, sizeof(bin_model_names) / sizeof(char*), model_names);
        } else {
            convert_string_vec(model_set.c_str(), model_names);
        }
    } else if (seq_type == SEQ_MORPH) {
        if (model_set.empty()) {
            copyCString(morph_model_names, sizeof(morph_model_names) / sizeof(char*), model_names);
        } else if (model_set[0] == '+') {
            // append model_set into existing models
            convert_string_vec(model_set.c_str()+1, model_names);
            appendCString(morph_model_names, sizeof(morph_model_names) / sizeof(char*), model_names);
        } else {
            convert_string_vec(model_set.c_str(), model_names);
        }
    } else if (seq_type == SEQ_DNA || seq_type == SEQ_POMO) {
        if (model_set.empty()) {
            copyCString(dna_model_names, sizeof(dna_model_names) / sizeof(char*), model_names);
            //            copyCString(dna_freq_names, sizeof(dna_freq_names)/sizeof(char*), freq_names);
        } else if (model_set == "partitionfinder" || model_set== "phyml") {
            copyCString(dna_model_names_old, sizeof(dna_model_names_old) / sizeof(char*), model_names);
            //            copyCString(dna_freq_names, sizeof(dna_freq_names)/sizeof(char*), freq_names);
        } else if (model_set == "raxml") {
            copyCString(dna_model_names_rax, sizeof(dna_model_names_rax) / sizeof(char*), model_names);
            //            copyCString(dna_freq_names, sizeof(dna_freq_names)/sizeof(char*), freq_names);
        } else if (model_set == "mrbayes") {
            copyCString(dna_model_names_mrbayes, sizeof(dna_model_names_mrbayes) / sizeof(char*), model_names);
            //            copyCString(dna_freq_names, sizeof(dna_freq_names)/sizeof(char*), freq_names);
        } else if (model_set == "beast1") {
            copyCString(dna_model_names_beast1, sizeof(dna_model_names_beast1) / sizeof(char*), model_names);
        } else if (model_set == "beast2") {
            copyCString(dna_model_names_beast2, sizeof(dna_model_names_beast2) / sizeof(char*), model_names);
        } else if (model_set == "modelomatic") {
            copyCString(dna_model_names_modelomatic, sizeof(dna_model_names_modelomatic) / sizeof(char*), model_names);
        } else if (model_set == "liemarkov") {
            copyCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ry, sizeof(dna_model_names_lie_markov_ry) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ws, sizeof(dna_model_names_lie_markov_ws) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_mk, sizeof(dna_model_names_lie_markov_mk) / sizeof(char*), model_names);
        } else if (model_set == "liemarkovry") {
            copyCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ry, sizeof(dna_model_names_lie_markov_ry) / sizeof(char*), model_names);
        } else if (model_set == "liemarkovws") {
            copyCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ws, sizeof(dna_model_names_lie_markov_ws) / sizeof(char*), model_names);
        } else if (model_set == "liemarkovmk") {
            copyCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_mk, sizeof(dna_model_names_lie_markov_mk) / sizeof(char*), model_names);
        } else if (model_set == "strandsymmetric") {
            copyCString(dna_model_names_lie_markov_strsym, sizeof(dna_model_names_lie_markov_strsym) / sizeof(char*), model_names);
            // IMPORTANT NOTE: If you add any more -mset names for sets of Lie Markov models,
            // you also need to change getPrototypeModel function.
        } else if (model_set[0] == '+') {
            // append model_set into existing models
            convert_string_vec(model_set.c_str()+1, model_names);
            appendCString(dna_model_names, sizeof(dna_model_names) / sizeof(char*), model_names);
        } else {
            convert_string_vec(model_set.c_str(), model_names);
        }
        
        if (model_name.find("+LMRY") != string::npos) {
            appendCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ry, sizeof(dna_model_names_lie_markov_ry) / sizeof(char*), model_names);
        } else if (model_name.find("+LMWS") != string::npos) {
            appendCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ws, sizeof(dna_model_names_lie_markov_ws) / sizeof(char*), model_names);
        } else if (model_name.find("+LMMK") != string::npos) {
            appendCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_mk, sizeof(dna_model_names_lie_markov_mk) / sizeof(char*), model_names);
        } else if (model_name.find("+LMSS") != string::npos) {
            appendCString(dna_model_names_lie_markov_strsym, sizeof(dna_model_names_lie_markov_strsym) / sizeof(char*), model_names);
        } else if (model_name.find("+LM") != string::npos) {
            appendCString(dna_model_names_lie_markov_fullsym, sizeof(dna_model_names_lie_markov_fullsym) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ry, sizeof(dna_model_names_lie_markov_ry) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_ws, sizeof(dna_model_names_lie_markov_ws) / sizeof(char*), model_names);
            appendCString(dna_model_names_lie_markov_mk, sizeof(dna_model_names_lie_markov_mk) / sizeof(char*), model_names);
        }
    } else if (seq_type == SEQ_PROTEIN) {
        if (model_set.empty()) {
            copyCString(aa_model_names, sizeof(aa_model_names) / sizeof(char*), model_names);
        } else if (model_set == "partitionfinder" || model_set == "phyml") {
            copyCString(aa_model_names_phyml, sizeof(aa_model_names_phyml) / sizeof(char*), model_names);
        } else if (model_set == "raxml") {
            copyCString(aa_model_names_rax, sizeof(aa_model_names_rax) / sizeof(char*), model_names);
        } else if (model_set == "mrbayes") {
            copyCString(aa_model_names_mrbayes, sizeof(aa_model_names_mrbayes) / sizeof(char*), model_names);
        } else if (model_set == "beast1") {
            copyCString(aa_model_names_beast1, sizeof(aa_model_names_beast1) / sizeof(char*), model_names);
        } else if (model_set == "beast2") {
            copyCString(aa_model_names_beast2, sizeof(aa_model_names_beast2) / sizeof(char*), model_names);
        } else if (model_set == "modelomatic") {
            copyCString(aa_model_names_modelomatic, sizeof(aa_model_names_modelomatic) / sizeof(char*), model_names);
        } else if (model_set[0] == '+') {
            // append model_set into existing models
            convert_string_vec(model_set.c_str()+1, model_names);
            appendCString(aa_model_names, sizeof(aa_model_names) / sizeof(char*), model_names);
        } else {
            convert_string_vec(model_set.c_str(), model_names);
        }
        
        if (model_subset) {
            StrVector submodel_names;
            if (strncmp(model_subset, "nuclear", 3) == 0) {
                copyCString(aa_model_names_nuclear, sizeof(aa_model_names_nuclear) / sizeof(char*), submodel_names);
            } else if (strncmp(model_subset, "mitochondrial", 3) == 0) {
                copyCString(aa_model_names_mitochondrial, sizeof(aa_model_names_mitochondrial) / sizeof(char*), submodel_names);
            } else if (strncmp(model_subset, "chloroplast", 3) == 0) {
                copyCString(aa_model_names_chloroplast, sizeof(aa_model_names_chloroplast) / sizeof(char*), submodel_names);
            } else if (strncmp(model_subset, "viral",3) == 0) {
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
        if (model_set.empty()) {
            if (standard_code)
                copyCString(codon_model_names, sizeof(codon_model_names) / sizeof(char*), model_names);
            else {
                i = sizeof(codon_model_names) / sizeof(char*);
                for (j = 0; j < i; j++)
                    if (!std_genetic_code[j])
                        model_names.push_back(codon_model_names[j]);
                //                copyCString(codon_model_names, sizeof(codon_model_names) / sizeof(char*) - 1, model_names);
            }
        } else if (model_set == "modelomatic") {
            copyCString(codon_model_names_modelomatic, sizeof(codon_model_names_modelomatic) / sizeof(char*), model_names);
        } else if (model_set[0] == '+') {
            // append model_set into existing models
            convert_string_vec(model_set.c_str()+1, model_names);
            if (standard_code)
                appendCString(codon_model_names, sizeof(codon_model_names) / sizeof(char*), model_names);
            else {
                i = sizeof(codon_model_names) / sizeof(char*);
                for (j = 0; j < i; j++)
                    if (!std_genetic_code[j])
                        model_names.push_back(codon_model_names[j]);
            }
        } else
            convert_string_vec(model_set.c_str(), model_names);
    }
}

void getStateFreqs(SeqType seq_type, char *state_freq_set, StrVector &freq_names) {
    int j;
    
    switch (seq_type) {
        case SEQ_PROTEIN:
            copyCString(aa_freq_names, sizeof(aa_freq_names)/sizeof(char*), freq_names);
            break;
        case SEQ_CODON:
            copyCString(codon_freq_names, sizeof(codon_freq_names) / sizeof(char*), freq_names);
            break;
        default:
            break;
    }
    if (state_freq_set)
        convert_string_vec(state_freq_set, freq_names);
    for (j = 0; j < freq_names.size(); j++) {
        std::transform(freq_names[j].begin(), freq_names[j].end(), freq_names[j].begin(), ::toupper);
        if (freq_names[j] != "" && freq_names[j][0] != '+')
            freq_names[j] = "+" + freq_names[j];
    }
}

/**
 get list of rate heterogeneity
 */
void getRateHet(SeqType seq_type, string model_name, double frac_invariant_sites,
                string rate_set, StrVector &ratehet) {
    const char *rate_options[]    = {  "", "+I", "+ASC", "+G", "+I+G", "+ASC+G", "+R", "+ASC+R", "+I+R"};
    bool test_options_default[]   = {true,   true, false,  true,  true,   false, false,  false, false};
    bool test_options_fast[]      = {false, false, false, false,  true,   false, false,  false, false};
    bool test_options_morph[]     = {true,  false,  true,  true, false,    true, false,  false, false};
    bool test_options_morph_fast[]= {false, false, false, false, false,    true, false,  false, false};
    bool test_options_noASC_I[]   = {true,  false, false,  true, false,   false, false,  false, false};
    bool test_options_noASC_I_fast[]={false,false, false,  true, false,   false, false,  false, false};
    bool test_options_asc[]       ={false,  false,  true, false, false,    true, false,  false, false};
    bool test_options_new[]       = {true,   true, false,  true,  true,   false,  true,  false, true};
    bool test_options_morph_new[] = {true,  false,  true,  true, false,    true,  true,   true, false};
    bool test_options_noASC_I_new[]= {true, false, false,  true, false,   false,  true,  false, false};
    bool test_options_asc_new[]   = {false, false,  true, false, false,    true, false,   true, false};
    bool test_options_pomo[]      = {true,  false, false,  true, false,   false, false,  false, false};
    bool test_options_norate[]    = {true,  false, false, false, false,   false, false,  false, false};
    bool *test_options = test_options_default;
    //    bool test_options_codon[] =  {true,false,  false,false,  false,    false};
    const int noptions = sizeof(rate_options) / sizeof(char*);
    int i, j;
    
    bool with_new = (model_name.find("NEW") != string::npos || model_name.substr(0,2) == "MF" || model_name.empty());
    bool with_asc = model_name.find("ASC") != string::npos;

    if (seq_type == SEQ_POMO) {
        for (i = 0; i < noptions; i++)
            test_options[i] = test_options_pomo[i];
    }
    // If not PoMo, go on with normal treatment.
    else if (frac_invariant_sites == 0.0) {
        // morphological or SNP data: activate +ASC
        if (with_new && rate_set != "1") {
            if (with_asc)
                test_options = test_options_asc_new;
            else if (seq_type == SEQ_DNA || seq_type == SEQ_BINARY || seq_type == SEQ_MORPH)
                test_options = test_options_morph_new;
            else
                test_options = test_options_noASC_I_new;
        } else if (with_asc)
            test_options = test_options_asc;
        else if (seq_type == SEQ_DNA || seq_type == SEQ_BINARY || seq_type == SEQ_MORPH) {
            if (rate_set == "1")
                test_options = test_options_morph_fast;
            else
                test_options = test_options_morph;
        } else {
            if (rate_set == "1")
                test_options = test_options_noASC_I_fast;
            else
                test_options = test_options_noASC_I;
        }
    } else if (frac_invariant_sites >= 1.0) {
        // 2018-06-12: alignment with only invariant sites, no rate variation added
        test_options = test_options_norate;
    } else {
        // normal data, use +I instead
        if (with_new && rate_set != "1") {
            // change +I+G to +R
            if (with_asc)
                test_options = test_options_asc_new;
            else
                test_options = test_options_new;
        } else if (with_asc) {
            test_options = test_options_asc;
        } else if (rate_set == "1")
            test_options = test_options_fast;
        else
            test_options = test_options_default;
        if (frac_invariant_sites == 0.0) {
            // deactivate +I
            for (j = 0; j < noptions; j++)
                if (strstr(rate_options[j], "+I"))
                    test_options[j] = false;
        }
    }
    if (!rate_set.empty() && rate_set != "1" && !iEquals(rate_set, "ALL") && !iEquals(rate_set, "AUTO")) {
        // take the rate_options from user-specified models
        convert_string_vec(rate_set.c_str(), ratehet);
        if (!ratehet.empty() && iEquals(ratehet[0], "ALL")) {
            ratehet.erase(ratehet.begin());
            StrVector ratedef;
            for (j = 0; j < noptions; j++)
                if (test_options[j])
                    ratedef.push_back(rate_options[j]);
            ratehet.insert(ratehet.begin(), ratedef.begin(), ratedef.end());
        }
        for (j = 0; j < ratehet.size(); j++) {
            if (ratehet[j] != "" && ratehet[j][0] != '+' && ratehet[j][0] != '*')
                ratehet[j] = "+" + ratehet[j];
            if (ratehet[j] == "+E") // for equal rate model
                ratehet[j] = "";
        }
    } else {
        for (j = 0; j < noptions; j++)
            if (test_options[j])
                ratehet.push_back(rate_options[j]);
        
    }
}

int CandidateModelSet::generate(Params &params, Alignment *aln, bool separate_rate, bool merge_phase) {
	StrVector model_names;
    StrVector freq_names;
	SeqType seq_type = aln->seq_type;
    
	int i, j;
    string model_set;
    
    if (merge_phase) {
        model_set = params.merge_models;
    } else
        model_set = params.model_set;

    bool auto_model = iEquals(model_set, "AUTO");
    
    getModelSubst(seq_type, aln->isStandardGeneticCode(), params.model_name,
                  model_set, params.model_subset, model_names);

	if (model_names.empty()) 
        return 1;
    
    getStateFreqs(seq_type, params.state_freq_set, freq_names);
    
    // combine model_names with freq_names
    if (freq_names.size() > 0) {
        StrVector orig_model_names = model_names;
        model_names.clear();
        for (j = 0; j < orig_model_names.size(); j++) {
            if (aln->seq_type == SEQ_CODON) {
                SeqType seq_type;
                int model_type = detectSeqType(orig_model_names[j].c_str(), seq_type);
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


    

    StrVector ratehet;
    int max_cats = params.num_rate_cats;
    string ratehet_set;
    if (merge_phase) {
        ratehet_set = params.merge_rates;
    } else
        ratehet_set = params.ratehet_set;

    //bool auto_rate = iEquals(ratehet_set, "AUTO");
    
    getRateHet(seq_type, params.model_name, aln->frac_invariant_sites, ratehet_set, ratehet);
    
    // add number of rate cateogories for special rate models
    const char *rates[] = {"+R", "*R", "+H", "*H"};

    for (i = 0; i < sizeof(rates)/sizeof(char*); i++)
        if (params.model_name.find(rates[i]) != string::npos)
            ratehet.push_back(rates[i]);

    size_t pos;

    vector<int> flags;
    flags.resize(ratehet.size(), 0);
    
    for (i = 0; i < sizeof(rates)/sizeof(char*); i++)
    for (j = 0; j < ratehet.size(); j++)
        if ((pos = ratehet[j].find(rates[i])) != string::npos &&
            (pos >= ratehet[j].length()-2 || !isdigit(ratehet[j][pos+2]) ))
        {
            string str = ratehet[j];
            ratehet[j].insert(pos+2, convertIntToString(params.min_rate_cats));
            max_cats = max(max_cats, params.max_rate_cats);
            for (int k = params.min_rate_cats+1; k <= params.max_rate_cats; k++) {
                int ins_pos = j+k-params.min_rate_cats;
                ratehet.insert(ratehet.begin() + ins_pos, str.substr(0, pos+2) + convertIntToString(k) + str.substr(pos+2));
                flags.insert(flags.begin() + ins_pos, MF_WAITING);
            }
        }

    ASSERT(ratehet.size() == flags.size());
    
    string pomo_suffix = (seq_type == SEQ_POMO) ? "+P" : "";
    // TODO DS: should we allow virtual population size?

    // combine substitution models with rate heterogeneity
    if (separate_rate) {
        for (i = 0; i < model_names.size(); i++) 
            push_back(CandidateModel(model_names[i], ratehet[0] + pomo_suffix, aln));
        for (j = 0; j < ratehet.size(); j++)
            if (ratehet[j] != "")
                push_back(CandidateModel("", ratehet[j] + pomo_suffix, aln));
    } else {
        if (auto_model) {
            // all rate heterogeneity for the first model
            for (j = 0; j < ratehet.size(); j++)
                push_back(CandidateModel(model_names[0], ratehet[j] + pomo_suffix, aln, flags[j]));
            // now all models the first RHAS
            for (i = 1; i < model_names.size(); i++)
                push_back(CandidateModel(model_names[i], ratehet[0] + pomo_suffix, aln, flags[0]));
            // all remaining models
            for (i = 1; i < model_names.size(); i++)
                for (j = 1; j < ratehet.size(); j++) {
                    push_back(CandidateModel(model_names[i], ratehet[j] + pomo_suffix, aln, flags[j]));
                }
        } else {
            // testing all models
            for (i = 0; i < model_names.size(); i++)
                for (j = 0; j < ratehet.size(); j++) {
                    push_back(CandidateModel(model_names[i], ratehet[j] + pomo_suffix, aln, flags[j]));
                }
        }
    }
    if (params.model_extra_set) {
        StrVector extra_model_names;
        convert_string_vec(params.model_extra_set, extra_model_names);
        for (auto s : extra_model_names)
            push_back(CandidateModel(s, "", aln));
    }
    return max_cats;
}

void replaceModelInfo(string &set_name, ModelCheckpoint &model_info, ModelCheckpoint &new_info) {
    for (auto it = new_info.begin(); it != new_info.end(); it++) {
        model_info.put(set_name + CKP_SEP + it->first, it->second);
    }
}

void extractModelInfo(string &orig_set_name, ModelCheckpoint &model_info, ModelCheckpoint &part_model_info) {
    string set_name = orig_set_name + CKP_SEP;
    int len = set_name.length();
    for (auto it = model_info.lower_bound(set_name); it != model_info.end() && it->first.substr(0, len) == set_name; it++) {
        part_model_info.put(it->first.substr(len), it->second);
    }
}

string getSubsetName(PhyloSuperTree *super_tree, set<int> &subset) {
    string set_name;
    for (auto it = subset.begin(); it != subset.end(); it++) {
        if (it != subset.begin())
            set_name += "+";
        set_name += super_tree->at(*it)->aln->name;
    }
    return set_name;
}

int getSubsetAlnLength(PhyloSuperTree *super_tree, set<int> &subset) {
    int len = 0;
    for (auto i : subset) {
        len += super_tree->at(i)->aln->getNSite();
    }
    return len;
}

/**
 * transfer model parameters from two subsets to the target subsets
 */
void transferModelParameters(PhyloSuperTree *super_tree, ModelCheckpoint &model_info, ModelCheckpoint &part_model_info,
                             set<int> &gene_set1, set<int> &gene_set2)
{
    set<int> merged_set;
    merged_set.insert(gene_set1.begin(), gene_set1.end());
    merged_set.insert(gene_set2.begin(), gene_set2.end());
    string set_name = getSubsetName(super_tree, merged_set);
    string set1_name = getSubsetName(super_tree, gene_set1);
    string set2_name = getSubsetName(super_tree, gene_set2);
    double weight1 = getSubsetAlnLength(super_tree, gene_set1);
    double weight2 = getSubsetAlnLength(super_tree, gene_set2);
    double weight_sum = weight1 + weight2;
    weight1 = weight1/weight_sum;
    weight2 = weight2/weight_sum;
    enum MeanComp {GEOM_MEAN, ARIT_MEAN};
    enum ValType {VAL_SINGLE, VAL_VECTOR};
    vector<tuple<ValType, MeanComp,string> > info_strings = {
        make_tuple(VAL_SINGLE, ARIT_MEAN, (string)"RateGamma" + CKP_SEP + "gamma_shape"),
        make_tuple(VAL_SINGLE, ARIT_MEAN, (string)"RateGammaInvar" + CKP_SEP + "gamma_shape"),
        make_tuple(VAL_SINGLE, ARIT_MEAN, (string)"RateGammaInvar" + CKP_SEP + "p_invar"),
        make_tuple(VAL_SINGLE, ARIT_MEAN, (string)"RateInvar" + CKP_SEP + "p_invar")
        //make_tuple(VAL_VECTOR, GEOM_MEAN, (string)"ModelDNA" + CKP_SEP + "rates")
    };
    for (auto info : info_strings) {
        switch (std::get<0>(info)) {
            case VAL_SINGLE: {
                double value1, value2, value;
                bool ok1 = model_info.get(set1_name + CKP_SEP + std::get<2>(info), value1);
                bool ok2 = model_info.get(set2_name + CKP_SEP + std::get<2>(info), value2);
                if (!ok1 || !ok2)
                    continue;
                if (part_model_info.get(std::get<2>(info), value))
                    continue; // value already exist
                switch (std::get<1>(info)) {
                    case ARIT_MEAN:
                        value = weight1*value1 + weight2*value2;
                        break;
                    case GEOM_MEAN:
                        value = sqrt(value1*value2);
                        break;
                }
                part_model_info.put(std::get<2>(info), value);
                break;
            }
            case VAL_VECTOR: {
                DoubleVector value1, value2, value;
                bool ok1 = model_info.getVector(set1_name + CKP_SEP + std::get<2>(info), value1);
                bool ok2 = model_info.getVector(set2_name + CKP_SEP + std::get<2>(info), value2);
                if (!ok1 || !ok2)
                    continue;
                ASSERT(value1.size() == value2.size());
                if (part_model_info.getVector(std::get<2>(info), value))
                    continue; // value already exist
                value.reserve(value1.size());
                for (int i = 0; i < value1.size(); i++)
                switch (std::get<1>(info)) {
                    case ARIT_MEAN:
                        value.push_back(weight1*value1[i] + weight2*value2[i]);
                        break;
                    case GEOM_MEAN:
                        value.push_back(sqrt(value1[i]*value2[i]));
                        break;
                }
                part_model_info.putVector(std::get<2>(info), value);
                break;
            }
        }
    }
}

void mergePartitions(PhyloSuperTree* super_tree, vector<set<int> > &gene_sets, StrVector &model_names) {
	cout << "Merging into " << gene_sets.size() << " partitions..." << endl;
	vector<set<int> >::iterator it;
	SuperAlignment *super_aln = (SuperAlignment*)super_tree->aln;
	vector<PartitionInfo> part_info;
	vector<PhyloTree*> tree_vec;
    SuperAlignment *new_super_aln = new SuperAlignment();
	for (it = gene_sets.begin(); it != gene_sets.end(); it++) {
        Alignment *aln = super_aln->concatenateAlignments(*it);
		PartitionInfo info;
		aln->model_name = model_names[it-gene_sets.begin()];
        info.part_rate = 1.0; // BIG FIX: make -spp works with -m TESTMERGE now!
        info.evalNNIs = 0;
		for (set<int>::iterator i = it->begin(); i != it->end(); i++) {
			if (i != it->begin()) {
				aln->name += "+";
                if (!super_aln->partitions[*i]->position_spec.empty())
                    aln->position_spec += ", ";
			}
			aln->name += super_aln->partitions[*i]->name;
			aln->position_spec += super_aln->partitions[*i]->position_spec;
			if (!super_aln->partitions[*i]->aln_file.empty()) {
                if (aln->aln_file.empty())
                    aln->aln_file = super_aln->partitions[*i]->aln_file;
                else if (aln->aln_file != super_aln->partitions[*i]->aln_file) {
                    aln->aln_file = aln->aln_file + ',' + super_aln->partitions[*i]->aln_file;
                }
			}
			if (!super_aln->partitions[*i]->sequence_type.empty()) {
                if (aln->sequence_type.empty())
                    aln->sequence_type = super_aln->partitions[*i]->sequence_type;
                else if (aln->sequence_type != super_aln->partitions[*i]->sequence_type) {
                    aln->sequence_type = "__NA__";
                }
			}
		}
		info.cur_ptnlh = NULL;
		info.nniMoves[0].ptnlh = NULL;
		info.nniMoves[1].ptnlh = NULL;
		part_info.push_back(info);
		PhyloTree *tree = super_tree->extractSubtree(*it);
        tree->setParams(super_tree->params);
		tree->setAlignment(aln);
		tree_vec.push_back(tree);
        new_super_aln->partitions.push_back(aln);
	}

    // BUG FIX 2016-11-29: when merging partitions with -m TESTMERGE, sequence order is changed
    // get the taxa names from existing tree
    StrVector seq_names;
    if (super_tree->root) {
        super_tree->getTaxaName(seq_names);
    }
    new_super_aln->init(&seq_names);
    
	for (PhyloSuperTree::reverse_iterator tit = super_tree->rbegin(); tit != super_tree->rend(); tit++)
		delete (*tit);
	super_tree->clear();
	super_tree->insert(super_tree->end(), tree_vec.begin(), tree_vec.end());
	super_tree->part_info = part_info;

	delete super_tree->aln;
//    super_tree->aln = new SuperAlignment(super_tree);
    super_tree->setAlignment(new_super_aln);
}

/**
 called when some partition is changed
 */
void fixPartitions(PhyloSuperTree* super_tree) {
    SuperAlignment *super_aln = (SuperAlignment*)super_tree->aln;
    int part;
    bool aln_changed = false;
    for (part = 0; part < super_tree->size(); part++)
        if (super_aln->partitions[part] != super_tree->at(part)->aln) {
            aln_changed = true;
            super_aln->partitions[part] = super_tree->at(part)->aln;
        }
    if (!aln_changed)
        return;
    super_aln->buildPattern();
    super_aln->orderPatternByNumChars(PAT_VARIANT);
    super_tree->deleteAllPartialLh();
}

string CandidateModel::evaluate(Params &params,
    ModelCheckpoint &in_model_info, ModelCheckpoint &out_model_info,
    ModelsBlock *models_block,
    int &num_threads, int brlen_type)
{
    //string model_name = name;
    Alignment *in_aln = aln;
    IQTree *iqtree = NULL;
    if (in_aln->isSuperAlignment()) {
        SuperAlignment *saln = (SuperAlignment*)in_aln;
        if (params.partition_type == BRLEN_OPTIMIZE)
            iqtree = new PhyloSuperTree(saln);
        else
            iqtree = new PhyloSuperTreePlen(saln, brlen_type);
        StrVector subst_names;
        StrVector rate_names;
        convert_string_vec(subst_name.c_str(), subst_names);
        convert_string_vec(rate_name.c_str(), rate_names);
        ASSERT(subst_names.size() == rate_names.size());
        for (int part = 0; part != subst_names.size(); part++)
            saln->partitions[part]->model_name = subst_names[part]+rate_names[part];
    } else if (posRateHeterotachy(getName()) != string::npos)
        iqtree = new PhyloTreeMixlen(in_aln, 0);
    else
        iqtree = new IQTree(in_aln);
    iqtree->setParams(&params);
    iqtree->setLikelihoodKernel(params.SSE);
    iqtree->optimize_by_newton = params.optimize_by_newton;
    iqtree->setNumThreads(num_threads);

    iqtree->setCheckpoint(&in_model_info);
#ifdef _OPENMP
#pragma omp critical
#endif
    iqtree->restoreCheckpoint();
    ASSERT(iqtree->root);
    iqtree->initializeModel(params, getName(), models_block);
    if (!iqtree->getModel()->isMixture() || in_aln->seq_type == SEQ_POMO) {
        subst_name = iqtree->getSubstName();
        rate_name = iqtree->getRateName();
    }


    if (restoreCheckpoint(&in_model_info)) {
        delete iqtree;
        return "";
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    iqtree->getModelFactory()->restoreCheckpoint();
    
    // now switch to the output checkpoint
    iqtree->getModelFactory()->setCheckpoint(&out_model_info);
    iqtree->setCheckpoint(&out_model_info);
    
    double new_logl;
    
    if (params.model_test_and_tree) {
        //--- PERFORM FULL TREE SEARCH PER MODEL ----//

        // BQM 2017-03-29: disable bootstrap
        int orig_num_bootstrap_samples = params.num_bootstrap_samples;
        int orig_gbo_replicates = params.gbo_replicates;
        params.num_bootstrap_samples = 0;
        params.gbo_replicates = 0;
        STOP_CONDITION orig_stop_condition = params.stop_condition;
        if (params.stop_condition == SC_BOOTSTRAP_CORRELATION)
            params.stop_condition = SC_UNSUCCESS_ITERATION;

        iqtree->aln->model_name = getName();
        
        cout << endl << "===> Testing model " << getName() << endl;

        if (iqtree->root) {
            // start from previous tree
            string initTree = iqtree->getTreeString();
            iqtree->getCheckpoint()->put("initTree", initTree);
            iqtree->saveCheckpoint();
        }

        iqtree->ensureNumberOfThreadsIsSet(nullptr);

        runTreeReconstruction(params, iqtree);
        new_logl = iqtree->computeLikelihood();
        tree_len = iqtree->treeLength();
        tree = iqtree->getTreeString();

        // restore original parameters
        // 2017-03-29: restore bootstrap replicates
        params.num_bootstrap_samples = orig_num_bootstrap_samples;
        params.gbo_replicates = orig_gbo_replicates;
        params.stop_condition = orig_stop_condition;

        int count = iqtree->getCheckpoint()->eraseKeyPrefix("finished");
        cout << count << " finished checkpoint entries erased" << endl;
        iqtree->getCheckpoint()->eraseKeyPrefix("CandidateSet");

    } else {
        //--- FIX TREE TOPOLOGY AND ESTIMATE MODEL PARAMETERS ----//

        if (verbose_mode >= VB_MED)
            cout << "Optimizing model " << getName() << endl;

        iqtree->ensureNumberOfThreadsIsSet(nullptr);
        iqtree->initializeAllPartialLh();

        for (int step = 0; step < 2; step++) {
            new_logl = iqtree->getModelFactory()->optimizeParameters(brlen_type, false,
                params.modelfinder_eps, TOL_GRADIENT_MODELTEST);
            tree_len = iqtree->treeLength();
            iqtree->getModelFactory()->saveCheckpoint();
            iqtree->saveCheckpoint();

            // check if logl(+R[k]) is worse than logl(+R[k-1])
            CandidateModel prev_info;
            if (!prev_info.restoreCheckpointRminus1(&in_model_info, this)) break;
            if (prev_info.logl < new_logl + params.modelfinder_eps) break;
            if (step == 0) {
                iqtree->getRate()->initFromCatMinusOne();
            } else if (new_logl < prev_info.logl - params.modelfinder_eps*10.0) {
                outWarning("Log-likelihood " + convertDoubleToString(new_logl) + " of " +
                           getName() + " worse than " + prev_info.getName() + " " + convertDoubleToString(prev_info.logl));
            }
        }

    }

    // sum in case of adjusted df and logl already stored
    df += iqtree->getModelFactory()->getNParameters(brlen_type);
    logl += new_logl;
    string tree_string = iqtree->getTreeString();

#ifdef _OPENMP
#pragma omp critical
    {
#endif
    saveCheckpoint(&in_model_info);
#ifdef _OPENMP
    }
#endif

    delete iqtree;
    return tree_string;
}

/** model information by merging two partitions */
struct ModelPair {
    /** score after merging */
    double score;
    /** ID of partition 1 */
    int part1;
    /** ID of partition 2 */
    int part2;
    /** log-likelihood */
    double logl;
    /** degree of freedom */
    int df;
    /** tree length */
    double tree_len;
    /** IDs of merged partitions */
    set<int> merged_set;
    /** set name */
    string set_name;
    /* best model name */
    string model_name;
};

class ModelPairSet : public multimap<double, ModelPair> {

public:

    /** insert a partition pair */
    void insertPair(ModelPair &pair) {
        insert(value_type(pair.score, pair));
    }

    /** 
        find the maximum compatible partition pairs
        @param num max number of pairs to return
    */
    void getCompatiblePairs(int num, ModelPairSet &res) {
        set<int> part_ids;

        for (auto it = begin(); it != end() && res.size() < num; it++) {

            // check for compatibility
            vector<int> overlap;
            set_intersection(part_ids.begin(), part_ids.end(),
                it->second.merged_set.begin(), it->second.merged_set.end(),
                std::back_inserter(overlap));

            if (!overlap.empty()) continue;

            // take the union
            part_ids.insert(it->second.merged_set.begin(), it->second.merged_set.end());

            // put the compatible pair to the set
            res.insertPair(it->second);
        }
    }

};

string CandidateModel::evaluateConcatenation(Params &params, SuperAlignment *super_aln,
    ModelCheckpoint &model_info, ModelsBlock *models_block, int num_threads)
{
    aln = super_aln->concatenateAlignments();
    size_t ssize = getUsualModel(aln);

    string concat_tree;

    cout << "Testing " << getName() << " on supermatrix..." << endl;
    concat_tree = evaluate(params, model_info, model_info,
        models_block, num_threads, BRLEN_OPTIMIZE);

    computeICScores(ssize);

    delete aln;
    aln = NULL;
    return concat_tree;
}

/**
 * k-means clustering of partitions using partition-specific tree length
 * @return score (AIC/BIC/etc.) of the clustering
 * @param[out] gene_sets
 * @param[out[ model_names
 */
double doKmeansClustering(Params &params, PhyloSuperTree *in_tree,
    int ncluster, DoubleVector &lenvec,
    ModelCheckpoint &model_info, ModelsBlock *models_block,
    int num_threads,
    vector<set<int> > &gene_sets, StrVector &model_names)
{
    
    cout << "k-means merging into " << ncluster << " partitions..." << endl;
    
    ASSERT(lenvec.size() == in_tree->size());
    int npart = in_tree->size();
    IntVector weights;
    weights.resize(npart, 1);
    int *clusters = new int[npart];
    double *centers = new double[ncluster];
    RunKMeans1D(npart, ncluster, lenvec.data(), weights.data(), centers, clusters);
    
    SuperAlignment *super_aln = ((SuperAlignment*)in_tree->aln);

    double lhsum = 0.0;
    int dfsum = 0;
    if (params.partition_type == BRLEN_FIX || params.partition_type == BRLEN_SCALE) {
        dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
        if (params.partition_type == BRLEN_SCALE)
            dfsum -= 1;
    }

    for (int cluster = 0; cluster < ncluster; cluster++) {
        string set_name;
        set<int> merged_set;
        for (int i = 0; i < in_tree->size(); i++)
            if (clusters[i] == cluster) {
                if (!set_name.empty())
                    set_name += "+";
                set_name += in_tree->at(i)->aln->name;
                merged_set.insert(i);
            }
        gene_sets.push_back(merged_set);
        CandidateModel best_model;
        bool done_before = false;
        {
            // if pairs previously examined, reuse the information
            model_info.startStruct(set_name);
            if (model_info.getBestModel(best_model.subst_name)) {
                best_model.restoreCheckpoint(&model_info);
                done_before = true;
            }
            model_info.endStruct();
        }
        ModelCheckpoint part_model_info;
        if (!done_before) {
            Alignment *aln = super_aln->concatenateAlignments(merged_set);
            PhyloTree *tree = in_tree->extractSubtree(merged_set);
            tree->setAlignment(aln);
            extractModelInfo(set_name, model_info, part_model_info);
            tree->num_precision = in_tree->num_precision;
            tree->setParams(&params);
            tree->sse = params.SSE;
            tree->optimize_by_newton = params.optimize_by_newton;
            tree->setNumThreads(params.model_test_and_tree ? num_threads : 1);
            /*if (params.model_test_and_tree) {
             tree->setCheckpoint(new Checkpoint());
             tree->saveCheckpoint();
             } else*/
            {
            tree->setCheckpoint(&part_model_info);
            // trick to restore checkpoint
            tree->restoreCheckpoint();
            tree->saveCheckpoint();
            }
            best_model = CandidateModelSet().test(params, tree, part_model_info, models_block,
                params.model_test_and_tree ? num_threads : 1, params.partition_type,
                set_name, "", true);
            best_model.restoreCheckpoint(&part_model_info);
            model_names.push_back(best_model.getName());
            delete tree;
            delete aln;
        }
        lhsum += best_model.logl;
        dfsum += best_model.df;
        {
            if (!done_before) {
                replaceModelInfo(set_name, model_info, part_model_info);
                model_info.dump();
                cout.width(4);
                cout << right << cluster+1 << " ";
                cout.width(12);
                cout << left << best_model.getName() << " ";
                cout.width(11);
                cout << best_model.logl << " " << set_name;
                cout << endl;
            }
        }
    }
    
    size_t ssize = in_tree->getAlnNSite();
    double score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
    cout << "k-means score for " << ncluster << " partitions: " << score << " (LnL: " << lhsum << "  " << "df: " << dfsum << ")" << endl;

    delete [] centers;
    delete [] clusters;
    return score;
}

class SubsetPair : public pair<int,int> {
public:
    // distance between two partition pairs
    double distance;
};

bool comparePairs(const SubsetPair &a, const SubsetPair &b) {
    return a.distance < b.distance;
}

bool comparePartition(const pair<int,double> &a, const pair<int, double> &b) {
    return a.second > b.second;
}

/**
 find k-closest partition pairs for rcluster algorithm
 */
void findClosestPairs(SuperAlignment *super_aln, DoubleVector &lenvec, vector<set<int> > &gene_sets,
                      double log_transform, vector<SubsetPair> &closest_pairs) {
    
    for (int part1 = 0; part1 < lenvec.size()-1; part1++)
        for (int part2 = part1+1; part2 < lenvec.size(); part2++)
            if (super_aln->partitions[*gene_sets[part1].begin()]->seq_type == super_aln->partitions[*gene_sets[part2].begin()]->seq_type &&
                super_aln->partitions[*gene_sets[part1].begin()]->genetic_code == super_aln->partitions[*gene_sets[part2].begin()]->genetic_code) {
                // only merge partitions of the same data type
                SubsetPair pair;
                pair.first = part1;
                pair.second = part2;
                if (log_transform)
                    pair.distance = fabs(log(lenvec[part1]) - log(lenvec[part2]));
                else
                    pair.distance = fabs(lenvec[part1] - lenvec[part2]);
                closest_pairs.push_back(pair);
            }
    if (!closest_pairs.empty() && Params::getInstance().partfinder_rcluster < 100) {
        // sort distance
        std::sort(closest_pairs.begin(), closest_pairs.end(), comparePairs);
        size_t num_pairs = round(closest_pairs.size() * (Params::getInstance().partfinder_rcluster/100.0));
        num_pairs = min(num_pairs, Params::getInstance().partfinder_rcluster_max);
        if (num_pairs <= 0) num_pairs = 1;
        closest_pairs.erase(closest_pairs.begin() + num_pairs, closest_pairs.end());
    }
}

/**
 merge vector src into dest, eliminating duplicates
 */
void mergePairs(vector<SubsetPair> &dest, vector<SubsetPair> &src) {
    unordered_set<string> dest_set;
    for (SubsetPair s: dest)
        dest_set.insert(convertIntToString(s.first) + "-" + convertIntToString(s.second));
    for (SubsetPair s: src)
        if (dest_set.find(convertIntToString(s.first) + "-" + convertIntToString(s.second)) == dest_set.end())
            dest.push_back(s);
}

/**
 * select models for all partitions
 * @param[in,out] model_info (IN/OUT) all model information
 * @return total number of parameters
 */
void testPartitionModel(Params &params, PhyloSuperTree* in_tree, ModelCheckpoint &model_info,
    ModelsBlock *models_block, int num_threads)
{
//    params.print_partition_info = true;
//    params.print_conaln = true;
	int i = 0;
//	PhyloSuperTree::iterator it;
	DoubleVector lhvec; // log-likelihood for each partition
	DoubleVector dfvec; // number of parameters for each partition
    DoubleVector lenvec; // tree length for each partition
	double lhsum = 0.0;
	int dfsum = 0;
    if (params.partition_type == BRLEN_FIX || params.partition_type == BRLEN_SCALE) {
        dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
        if (params.partition_type == BRLEN_SCALE)
            dfsum -= 1;
    }
	size_t  ssize = in_tree->getAlnNSite();
	int64_t num_model = 0;
    int64_t total_num_model = in_tree->size();

    // 2017-06-07: -rcluster-max for max absolute number of pairs
    if (params.partfinder_rcluster_max == 0) {
        params.partfinder_rcluster_max = max((size_t)1000, 10 * in_tree->size());
    }

	if (params.partition_merge != MERGE_NONE) {
        double p = params.partfinder_rcluster/100.0;
        size_t num_pairs = round(in_tree->size()*(in_tree->size()-1)*p/2);
        if (p < 1.0)
            num_pairs = min(num_pairs, params.partfinder_rcluster_max);
        total_num_model += num_pairs;
        for (i = in_tree->size()-2; i > 0; i--)
            total_num_model += max(round(i*p), 1.0);
    }


#ifdef _OPENMP
    if (num_threads <= 0) {
        // partition selection scales well with many cores
        num_threads = min((int64_t)countPhysicalCPUCores(), total_num_model);
        num_threads = min(num_threads, params.num_threads_max);
        omp_set_num_threads(num_threads);
        cout << "NUMBER OF THREADS FOR PARTITION FINDING: " << num_threads << endl;
    }
#endif

    double start_time = getRealTime();

	SuperAlignment *super_aln = ((SuperAlignment*)in_tree->aln);

	cout << "Selecting individual models for " << in_tree->size() << " charsets using " << criterionName(params.model_test_criterion) << "..." << endl;
	//cout << " No. AIC         AICc        BIC         Charset" << endl;
	cout << " No. Model        Score       TreeLen     Charset" << endl;

	lhvec.resize(in_tree->size());
	dfvec.resize(in_tree->size());
	lenvec.resize(in_tree->size());

    // sort partition by computational cost for OpenMP effciency
    vector<pair<int,double> > partitionID;
    
	for (i = 0; i < in_tree->size(); i++) {
        Alignment *this_aln = in_tree->at(i)->aln;
        // computation cost is proportional to #sequences, #patterns, and #states
        partitionID.push_back({i, ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states});
    }
    if (num_threads > 1) {
        std::sort(partitionID.begin(), partitionID.end(), comparePartition);
    }
    bool parallel_over_partitions = false;
    int brlen_type = params.partition_type;
    if (brlen_type == TOPO_UNLINKED) {
        brlen_type = BRLEN_OPTIMIZE;
    }
    bool test_merge = (params.partition_merge != MERGE_NONE) && params.partition_type != TOPO_UNLINKED && (in_tree->size() > 1);
    
#ifdef _OPENMP
    parallel_over_partitions = !params.model_test_and_tree && (in_tree->size() >= num_threads);
#pragma omp parallel for private(i) schedule(dynamic) reduction(+: lhsum, dfsum) if(parallel_over_partitions)
#endif
	for (int j = 0; j < in_tree->size(); j++) {
        i = partitionID[j].first;
        PhyloTree *this_tree = in_tree->at(i);
		// scan through models for this partition, assuming the information occurs consecutively
		ModelCheckpoint part_model_info;
		extractModelInfo(this_tree->aln->name, model_info, part_model_info);
		// do the computation
        string part_model_name;
        if (params.model_name.empty())
            part_model_name = this_tree->aln->model_name;
        CandidateModel best_model;
		best_model = CandidateModelSet().test(params, this_tree, part_model_info, models_block,
            (parallel_over_partitions ? 1 : num_threads), brlen_type, this_tree->aln->name, part_model_name, test_merge);

        bool check = (best_model.restoreCheckpoint(&part_model_info));
        ASSERT(check);

		double score = best_model.computeICScore(this_tree->getAlnNSite());
		this_tree->aln->model_name = best_model.getName();
		lhsum += (lhvec[i] = best_model.logl);
		dfsum += (dfvec[i] = best_model.df);
        lenvec[i] = best_model.tree_len;

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            cout.width(12);
            cout << left << best_model.getName() << " ";
            cout.width(11);
            cout << score << " ";
            cout.width(11);
            cout << best_model.tree_len << " ";
            cout << this_tree->aln->name;
            if (num_model >= 10) {
                double remain_time = (total_num_model-num_model)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " (" 
                    << convert_time(remain_time) << " left)";
            }
            cout << endl;
            replaceModelInfo(this_tree->aln->name, model_info, part_model_info);
            model_info.dump();
        }
    }

    // in case ModelOMatic change the alignment
    fixPartitions(in_tree);
    
	double inf_score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
	cout << "Full partition model " << criterionName(params.model_test_criterion)
         << " score: " << inf_score << " (LnL: " << lhsum << "  df:" << dfsum << ")" << endl;

	if (!test_merge) {
		super_aln->printBestPartition((string(params.out_prefix) + ".best_scheme.nex").c_str());
		super_aln->printBestPartitionRaxml((string(params.out_prefix) + ".best_scheme").c_str());
        model_info.dump();
		return;
	}

    vector<set<int> > gene_sets;
    StrVector model_names;
    StrVector greedy_model_trees;

    gene_sets.resize(in_tree->size());
    model_names.resize(in_tree->size());
    greedy_model_trees.resize(in_tree->size());
    for (i = 0; i < gene_sets.size(); i++) {
        gene_sets[i].insert(i);
        model_names[i] = in_tree->at(i)->aln->model_name;
        greedy_model_trees[i] = in_tree->at(i)->aln->name;
    }

    if (params.partition_merge == MERGE_KMEANS) {
        // kmeans cluster based on parition tree length
        double cur_score = inf_score;
        for (int ncluster = in_tree->size()-1; ncluster >= 1; ncluster--) {
            vector<set<int> > this_gene_sets;
            StrVector this_model_names;
            //double sum = in_tree->size()/std::accumulate(lenvec.begin(), lenvec.end(), 0.0);
            double score = doKmeansClustering(params, in_tree, ncluster, lenvec, model_info,
                models_block, num_threads, this_gene_sets, this_model_names);
            if (score < cur_score) {
                cout << "Better score found: " << score << endl;
                cur_score = score;
                gene_sets = this_gene_sets;
                model_names = this_model_names;
            } else {
                //break;
            }
        }
    } else {
        cout << "Merging models to increase model fit (about " << total_num_model << " total partition schemes)..." << endl;
    }

    /* following implements the greedy algorithm of Lanfear et al. (2012) */
	while (params.partition_merge != MERGE_KMEANS && gene_sets.size() >= 2) {
		// stepwise merging charsets

        // list of all better pairs of partitions than current partitioning scheme
        ModelPairSet better_pairs;

        // 2015-06-24: begin rcluster algorithm
        // compute distance between gene_sets
        ASSERT(gene_sets.size() == lenvec.size());
        // find closest partition pairs
        vector<SubsetPair> closest_pairs;
        findClosestPairs(super_aln, lenvec, gene_sets, false, closest_pairs);
        if (params.partfinder_log_rate) {
            // additional consider pairs by log-rate
            vector<SubsetPair> log_closest_pairs;
            findClosestPairs(super_aln, lenvec, gene_sets, true, log_closest_pairs);
            mergePairs(closest_pairs, log_closest_pairs);
        }
        // sort partition by computational cost for OpenMP effciency
        for (i = 0; i < closest_pairs.size(); i++) {
            // computation cost is proportional to #sequences, #patterns, and #states
            Alignment *this_aln = in_tree->at(closest_pairs[i].first)->aln;
            closest_pairs[i].distance = -((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
            this_aln = in_tree->at(closest_pairs[i].second)->aln;
            closest_pairs[i].distance -= ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
        }
        if (num_threads > 1) {
            std::sort(closest_pairs.begin(), closest_pairs.end(), comparePairs);
        }
        size_t num_pairs = closest_pairs.size();
        
#ifdef _OPENMP
#pragma omp parallel for private(i) schedule(dynamic) if(!params.model_test_and_tree)
#endif
        for (size_t pair = 0; pair < num_pairs; pair++) {
            // information of current partitions pair
            ModelPair cur_pair;
            cur_pair.part1 = closest_pairs[pair].first;
            cur_pair.part2 = closest_pairs[pair].second;
            ASSERT(cur_pair.part1 < cur_pair.part2);
            cur_pair.merged_set.insert(gene_sets[cur_pair.part1].begin(), gene_sets[cur_pair.part1].end());
            cur_pair.merged_set.insert(gene_sets[cur_pair.part2].begin(), gene_sets[cur_pair.part2].end());
            cur_pair.set_name = getSubsetName(in_tree, cur_pair.merged_set);
            double weight1 = getSubsetAlnLength(in_tree, gene_sets[cur_pair.part1]);
            double weight2 = getSubsetAlnLength(in_tree, gene_sets[cur_pair.part2]);
            double sum = 1.0 / (weight1 + weight2);
            weight1 *= sum;
            weight2 *= sum;
            CandidateModel best_model;
            bool done_before = false;
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                // if pairs previously examined, reuse the information
                model_info.startStruct(cur_pair.set_name);
                if (model_info.getBestModel(best_model.subst_name)) {
                    best_model.restoreCheckpoint(&model_info);
                    done_before = true;
                }
                model_info.endStruct();
            }
            ModelCheckpoint part_model_info;
            double cur_tree_len = 0.0;
            if (!done_before) {
                Alignment *aln = super_aln->concatenateAlignments(cur_pair.merged_set);
                PhyloTree *tree = in_tree->extractSubtree(cur_pair.merged_set);
                //tree->scaleLength((weight1*lenvec[cur_pair.part1] + weight2*lenvec[cur_pair.part2])/tree->treeLength());
                tree->scaleLength(sqrt(lenvec[cur_pair.part1]*lenvec[cur_pair.part2])/tree->treeLength());
                cur_tree_len = tree->treeLength();
                tree->setAlignment(aln);
                extractModelInfo(cur_pair.set_name, model_info, part_model_info);
                transferModelParameters(in_tree, model_info, part_model_info, gene_sets[cur_pair.part1], gene_sets[cur_pair.part2]);
                tree->num_precision = in_tree->num_precision;
                tree->setParams(&params);
                tree->sse = params.SSE;
                tree->optimize_by_newton = params.optimize_by_newton;
                tree->setNumThreads(params.model_test_and_tree ? num_threads : 1);
                {
                    tree->setCheckpoint(&part_model_info);
                    // trick to restore checkpoint
                    tree->restoreCheckpoint();
                    tree->saveCheckpoint();
                }
                best_model = CandidateModelSet().test(params, tree, part_model_info, models_block,
                    params.model_test_and_tree ? num_threads : 1, params.partition_type, cur_pair.set_name, "", true);
                best_model.restoreCheckpoint(&part_model_info);
                delete tree;
                delete aln;
            }
            cur_pair.logl = best_model.logl;
            cur_pair.df = best_model.df;
            cur_pair.model_name = best_model.getName();
            cur_pair.tree_len = best_model.tree_len;
            double lhnew = lhsum - lhvec[cur_pair.part1] - lhvec[cur_pair.part2] + best_model.logl;
            int dfnew = dfsum - dfvec[cur_pair.part1] - dfvec[cur_pair.part2] + best_model.df;
            cur_pair.score = computeInformationScore(lhnew, dfnew, ssize, params.model_test_criterion);
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				if (!done_before) {
					replaceModelInfo(cur_pair.set_name, model_info, part_model_info);
                    model_info.dump();
                    num_model++;
					cout.width(4);
					cout << right << num_model << " ";
					cout.width(12);
					cout << left << best_model.getName() << " ";
					cout.width(11);
                    cout << cur_pair.score << " ";
                    cout.width(11);
                    cout << cur_pair.tree_len << " " << cur_pair.set_name;
                    if (num_model >= 10) {
                        double remain_time = max(total_num_model-num_model, (int64_t)0)*(getRealTime()-start_time)/num_model;
                        cout << "\t" << convert_time(getRealTime()-start_time) << " (" 
                            << convert_time(remain_time) << " left)";
                    }
                    cout << endl;
				}
                if (cur_pair.score < inf_score)
                    better_pairs.insertPair(cur_pair);
			}

        }
		if (better_pairs.empty()) break;
        ModelPairSet compatible_pairs;

        int num_comp_pairs = params.partition_merge == MERGE_RCLUSTERF ? gene_sets.size()/2 : 1;
        better_pairs.getCompatiblePairs(num_comp_pairs, compatible_pairs);
        if (compatible_pairs.size() > 1)
            cout << compatible_pairs.size() << " compatible better partition pairs found" << endl;

        // 2017-12-21: simultaneously merging better pairs
        for (auto it_pair = compatible_pairs.begin(); it_pair != compatible_pairs.end(); it_pair++) {
            ModelPair opt_pair = it_pair->second;

            lhsum = lhsum - lhvec[opt_pair.part1] - lhvec[opt_pair.part2] + opt_pair.logl;
            dfsum = dfsum - dfvec[opt_pair.part1] - dfvec[opt_pair.part2] + opt_pair.df;
            inf_score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
            ASSERT(inf_score <= opt_pair.score + 0.1);

            cout << "Merging " << opt_pair.set_name << " with " << criterionName(params.model_test_criterion)
                 << " score: " << inf_score << " (LnL: " << lhsum << "  df: " << dfsum << ")" << endl;
            // change entry opt_part1 to merged one
            gene_sets[opt_pair.part1] = opt_pair.merged_set;
            lhvec[opt_pair.part1] = opt_pair.logl;
            dfvec[opt_pair.part1] = opt_pair.df;
            lenvec[opt_pair.part1] = opt_pair.tree_len;
            model_names[opt_pair.part1] = opt_pair.model_name;
            greedy_model_trees[opt_pair.part1] = "(" + greedy_model_trees[opt_pair.part1] + "," +
                greedy_model_trees[opt_pair.part2] + ")" +
                convertIntToString(in_tree->size()-gene_sets.size()+1) + ":" +
                convertDoubleToString(inf_score);

            // delete entry opt_part2
            lhvec.erase(lhvec.begin() + opt_pair.part2);
            dfvec.erase(dfvec.begin() + opt_pair.part2);
            lenvec.erase(lenvec.begin() + opt_pair.part2);
            gene_sets.erase(gene_sets.begin() + opt_pair.part2);
            model_names.erase(model_names.begin() + opt_pair.part2);
            greedy_model_trees.erase(greedy_model_trees.begin() + opt_pair.part2);

            // decrease part ID for all pairs beyond opt_pair.part2
            auto next_pair = it_pair;
            for (next_pair++; next_pair != compatible_pairs.end(); next_pair++) {
                if (next_pair->second.part1 > opt_pair.part2)
                    next_pair->second.part1--;
                if (next_pair->second.part2 > opt_pair.part2)
                    next_pair->second.part2--;
            }
        }
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

	cout << "Agglomerative model selection: " << final_model_tree << endl;
    
    if (gene_sets.size() < in_tree->size())
        mergePartitions(in_tree, gene_sets, model_names);

    if (!iEquals(params.merge_models, "all")) {
        // test all candidate models again
        lhsum = 0.0;
        dfsum = 0;
        if (params.partition_type == BRLEN_FIX || params.partition_type == BRLEN_SCALE) {
            dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
            if (params.partition_type == BRLEN_SCALE)
                dfsum -= 1;
        }

        // sort partition by computational cost for OpenMP effciency
        partitionID.clear();
        for (i = 0; i < in_tree->size(); i++) {
            Alignment *this_aln = in_tree->at(i)->aln;
            // computation cost is proportional to #sequences, #patterns, and #states
            partitionID.push_back({i, ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states});
        }
        
        if (num_threads > 1) {
            std::sort(partitionID.begin(), partitionID.end(), comparePartition);
        }

    #ifdef _OPENMP
        parallel_over_partitions = !params.model_test_and_tree && (in_tree->size() >= num_threads);
        #pragma omp parallel for private(i) schedule(dynamic) reduction(+: lhsum, dfsum) if(parallel_over_partitions)
    #endif
        for (int j = 0; j < in_tree->size(); j++) {
            i = partitionID[j].first;
            PhyloTree *this_tree = in_tree->at(i);
            // scan through models for this partition, assuming the information occurs consecutively
            ModelCheckpoint part_model_info;
            extractModelInfo(this_tree->aln->name, model_info, part_model_info);
            // do the computation
            string part_model_name;
            if (params.model_name.empty())
                part_model_name = this_tree->aln->model_name;
            CandidateModel best_model;
            best_model = CandidateModelSet().test(params, this_tree, part_model_info, models_block,
                (parallel_over_partitions ? 1 : num_threads), brlen_type,
                this_tree->aln->name, part_model_name, false);
            
            bool check = (best_model.restoreCheckpoint(&part_model_info));
            ASSERT(check);
            
            double score = best_model.computeICScore(this_tree->getAlnNSite());
            this_tree->aln->model_name = best_model.getName();
            lhsum += (lhvec[i] = best_model.logl);
            dfsum += (dfvec[i] = best_model.df);
            lenvec[i] = best_model.tree_len;
            
    #ifdef _OPENMP
    #pragma omp critical
    #endif
            {
            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            cout.width(12);
            cout << left << best_model.getName() << " ";
            cout.width(11);
            cout << score << " " << this_tree->aln->name;
            if (num_model >= 10) {
                double remain_time = (total_num_model-num_model)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " ("
                << convert_time(remain_time) << " left)";
            }
            cout << endl;
            replaceModelInfo(this_tree->aln->name, model_info, part_model_info);
            model_info.dump();
            }
        }
    }

    inf_score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
    cout << "Best partition model " << criterionName(params.model_test_criterion) << " score: " << inf_score << " (LnL: " << lhsum << "  df:" << dfsum << ")" << endl;

    ((SuperAlignment*)in_tree->aln)->printBestPartition((string(params.out_prefix) + ".best_scheme.nex").c_str());
	((SuperAlignment*)in_tree->aln)->printBestPartitionRaxml((string(params.out_prefix) + ".best_scheme").c_str());
    model_info.dump();
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

void CandidateModelSet::filterRates(int finished_model) {
    if (Params::getInstance().score_diff_thres < 0)
        return;
    double best_score = DBL_MAX;
    ASSERT(finished_model >= 0);
    int model;
    for (model = 0; model <= finished_model; model++)
        if (at(model).subst_name == at(0).subst_name) {
            if (!at(model).hasFlag(MF_DONE + MF_IGNORED))
                return; // only works if all models done
            best_score = min(best_score, at(model).getScore());
        }
    
    double ok_score = best_score + Params::getInstance().score_diff_thres;
    set<string> ok_rates;
    for (model = 0; model <= finished_model; model++)
        if (at(model).getScore() <= ok_score) {
            string rate_name = at(model).orig_rate_name;
            ok_rates.insert(rate_name);
        }
    for (model = finished_model+1; model < size(); model++)
        if (ok_rates.find(at(model).orig_rate_name) == ok_rates.end())
            at(model).setFlag(MF_IGNORED);
}

void CandidateModelSet::filterSubst(int finished_model) {
    if (Params::getInstance().score_diff_thres < 0)
        return;
    double best_score = DBL_MAX;
    ASSERT(finished_model >= 0);
    int model;
    for (model = 0; model <= finished_model; model++)
        if (at(model).rate_name == at(0).rate_name)
            best_score = min(best_score, at(model).getScore());
    
    double ok_score = best_score + Params::getInstance().score_diff_thres;
    set<string> ok_model;
    for (model = 0; model <= finished_model; model++) {
        if (at(model).rate_name != at(0).rate_name)
            continue;
        if (at(model).getScore() <= ok_score) {
            string subst_name = at(model).orig_subst_name;
            ok_model.insert(subst_name);
        } else
            at(model).setFlag(MF_IGNORED);
    }
    for (model = finished_model+1; model < size(); model++)
        if (ok_model.find(at(model).orig_subst_name) == ok_model.end())
            at(model).setFlag(MF_IGNORED);
}

CandidateModel CandidateModelSet::test(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
    ModelsBlock *models_block, int num_threads, int brlen_type,
    string set_name, string in_model_name, bool merge_phase)
{
    ModelCheckpoint *checkpoint = &model_info;

	in_tree->params = &params;
    
    // for ModelOMatic
    Alignment *prot_aln = NULL;
    Alignment *dna_aln = NULL;
    bool do_modelomatic = params.modelomatic && in_tree->aln->seq_type == SEQ_CODON;

    if (in_model_name.empty()) {
        generate(params, in_tree->aln, params.model_test_separate_rate, merge_phase);
        if (do_modelomatic) {
            // generate models for protein
            // adapter coefficient according to Whelan et al. 2015
            prot_aln = in_tree->aln->convertCodonToAA();
            int adjusted_df;
            double adjusted_logl = computeAdapter(in_tree->aln, prot_aln, adjusted_df);
            if (set_name.empty())
                cout << "Adjusted LnL: " << adjusted_logl << "  df: " << adjusted_df << endl;
            size_t start = size();
            generate(params, prot_aln, params.model_test_separate_rate, merge_phase);
            size_t i;
            for (i = start; i < size(); i++) {
                at(i).logl = adjusted_logl;
                at(i).df = adjusted_df;
            }

            // generate models for DNA
            dna_aln = in_tree->aln->convertCodonToDNA();
            start = size();
            generate(params, dna_aln, params.model_test_separate_rate, merge_phase);
            for (i = start; i < size(); i++) {
                at(i).setFlag(MF_SAMPLE_SIZE_TRIPLE);
            }
        }
    } else {
        push_back(CandidateModel(in_model_name, "", in_tree->aln));
    }

    DoubleVector model_scores;
    int model;
	int best_model = -1;
    Alignment *best_aln = in_tree->aln;

	int ssize = in_tree->aln->getNSite(); // sample size
    //if (adjust)
    //    ssize = adjust->sample_size;
	if (params.model_test_sample_size)
		ssize = params.model_test_sample_size;
	if (set_name == "") {
        cout << "ModelFinder will test up to " << size() << " ";
        if (do_modelomatic)
            cout << "codon/AA/DNA";
        else
            cout << getSeqTypeName(in_tree->aln->seq_type);
        cout << " models (sample size: " << ssize << ") ..." << endl;
        if (params.model_test_and_tree == 0)
            cout << " No. Model         -LnL         df  AIC          AICc         BIC" << endl;
	}

    //	uint64_t RAM_requirement = 0;
    int best_model_AIC = -1, best_model_AICc = -1, best_model_BIC = -1;
    double best_score_AIC = DBL_MAX, best_score_AICc = DBL_MAX, best_score_BIC = DBL_MAX;
    string best_tree_AIC, best_tree_AICc, best_tree_BIC;

//    CKP_RESTORE(best_score_AIC);
//    CKP_RESTORE(best_score_AICc);
//    CKP_RESTORE(best_score_BIC);
//    CKP_RESTORE(best_model_AIC);
//    CKP_RESTORE(best_model_AICc);
//    CKP_RESTORE(best_model_BIC);

    CKP_RESTORE(best_tree_AIC);
    CKP_RESTORE(best_tree_AICc);
    CKP_RESTORE(best_tree_BIC);

    // detect rate hetegeneity automatically or not
    bool auto_rate = merge_phase ? iEquals(params.merge_rates, "AUTO") : iEquals(params.ratehet_set, "AUTO");
    bool auto_subst = merge_phase ? iEquals(params.merge_models, "AUTO") : iEquals(params.model_set, "AUTO");
    int rate_block = size();
    if (auto_rate) {
        for (rate_block = 0; rate_block < size(); rate_block++)
            if (rate_block+1 < size() && at(rate_block+1).subst_name != at(rate_block).subst_name)
                break;
    }
    
    int subst_block = size();
    if (auto_subst) {
        for (subst_block = size()-1; subst_block >= 0; subst_block--)
            if (at(subst_block).rate_name == at(0).rate_name)
                break;
    }
    
    
    //------------- MAIN FOR LOOP GOING THROUGH ALL MODELS TO BE TESTED ---------//

	for (model = 0; model < size(); model++) {
        if (model == rate_block+1)
            filterRates(rate_block); // auto filter rate models
        if (model == subst_block+1)
            filterSubst(subst_block); // auto filter substitution model
        if (at(model).hasFlag(MF_IGNORED)) {
            model_scores.push_back(DBL_MAX);
            continue;
        }
		//cout << model_names[model] << endl;
        if (at(model).subst_name == "") {
            // now switching to test rate heterogeneity
            if (best_model == -1)
                switch (params.model_test_criterion) {
                case MTC_AIC:
                    best_model = best_model_AIC;
                    break;
                case MTC_AICC:
                    best_model = best_model_AICc;
                    break;
                case MTC_BIC:
                    best_model = best_model_BIC;
                    break;
                default: ASSERT(0);
                }
            at(model).subst_name = at(best_model).subst_name;
        }

		// optimize model parameters
        string orig_model_name = at(model).getName();
        // keep separate output model_info to only update model_info if better model found
        ModelCheckpoint out_model_info;
		//CandidateModel info;
		//info.set_name = set_name;
        at(model).set_name = set_name;
        string tree_string;

        /***** main call to estimate model parameters ******/
        tree_string = at(model).evaluate(params,
            model_info, out_model_info, models_block, num_threads, brlen_type);

        at(model).computeICScores(ssize);
        at(model).setFlag(MF_DONE);

        CandidateModel prev_info;

        bool skip_model = false;

        if (prev_info.restoreCheckpointRminus1(checkpoint, &at(model))) {
            // check stop criterion for +R
            prev_info.computeICScores(ssize);
            switch (params.model_test_criterion) {
            case MTC_ALL:
                if (at(model).AIC_score > prev_info.AIC_score &&
                    at(model).AICc_score > prev_info.AICc_score &&
                    at(model).BIC_score > prev_info.BIC_score) {
                    // skip remaining model
                    skip_model = true;
                }
                break;
            case MTC_AIC:
                if (at(model).AIC_score > prev_info.AIC_score) {
                    // skip remaining model
                    skip_model = true;
                }
                break;
            case MTC_AICC:
                if (at(model).AICc_score > prev_info.AICc_score) {
                    // skip remaining model
                    skip_model = true;
                }
                break;
            case MTC_BIC:
                if (at(model).BIC_score > prev_info.BIC_score) {
                    // skip remaining model
                    skip_model = true;
                }
                break;
            }
        }

		if (at(model).AIC_score < best_score_AIC) {
            best_model_AIC = model;
            best_score_AIC = at(model).AIC_score;
            if (!tree_string.empty())
                best_tree_AIC = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_AIC) {
                model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
            }
        }
		if (at(model).AICc_score < best_score_AICc) {
            best_model_AICc = model;
            best_score_AICc = at(model).AICc_score;
            if (!tree_string.empty())
                best_tree_AICc = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_AICC) {
                model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
            }
        }

		if (at(model).BIC_score < best_score_BIC) {
			best_model_BIC = model;
            best_score_BIC = at(model).BIC_score;
            if (!tree_string.empty())
                best_tree_BIC = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_BIC) {
                model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
            }
        }

        switch (params.model_test_criterion) {
            case MTC_AIC: model_scores.push_back(at(model).AIC_score); break;
            case MTC_AICC: model_scores.push_back(at(model).AICc_score); break;
            default: model_scores.push_back(at(model).BIC_score); break;
        }

        CKP_SAVE(best_tree_AIC);
        CKP_SAVE(best_tree_AICc);
        CKP_SAVE(best_tree_BIC);
        checkpoint->dump();


		if (set_name == "") {
            cout.width(3);
            cout << right << model+1 << "  ";
            cout.width(13);
            cout << left << at(model).getName() << " ";
            
            cout.precision(3);
            cout << fixed;
            cout.width(12);
            cout << -at(model).logl << " ";
            cout.width(3);
            cout << at(model).df << " ";
            cout.width(12);
            cout << at(model).AIC_score << " ";
            cout.width(12);
            cout << at(model).AICc_score << " " << at(model).BIC_score;
            cout << endl;
        }

        if (skip_model) {
            // skip over all +R model of higher categories
            const char *rates[] = {"+R", "*R", "+H", "*H"};
            size_t posR;
            for (int i = 0; i < sizeof(rates)/sizeof(char*); i++)
                if ((posR = orig_model_name.find(rates[i])) != string::npos)
                    break;
            string first_part = orig_model_name.substr(0, posR+2);
            for (int next = model+1; next < size() && at(next).getName().substr(0, posR+2) == first_part; next++) {
                at(next).setFlag(MF_IGNORED);
            }
        }
	}

    ASSERT(model_scores.size() == size());

    if (best_model_BIC == -1) {
        outError("No models were examined! Please check messages above");
    }
	int *model_rank = new int[model_scores.size()];

//    string best_tree; // BQM 2015-07-21: With Lars find best model
	/* sort models by their scores */
	switch (params.model_test_criterion) {
	case MTC_AIC:
		best_model = best_model_AIC;
		break;
	case MTC_AICC:
		best_model = best_model_AICc;
		break;
	case MTC_BIC:
		best_model = best_model_BIC;
		break;
    default: ASSERT(0);
	}
	sort_index(model_scores.data(), model_scores.data() + model_scores.size(), model_rank);

    string model_list;
	for (model = 0; model < model_scores.size(); model++) {
        if (model_scores[model_rank[model]] == DBL_MAX)
            break;
        if (model > 0)
            model_list += " ";
		model_list += at(model_rank[model]).getName();
    }

    model_info.putBestModelList(model_list);
    
    model_info.put("best_model_AIC", at(best_model_AIC).getName());
    model_info.put("best_model_AICc", at(best_model_AICc).getName());
    model_info.put("best_model_BIC", at(best_model_BIC).getName());

    CKP_SAVE(best_score_AIC);
    CKP_SAVE(best_score_AICc);
    CKP_SAVE(best_score_BIC);

    checkpoint->dump();

	delete [] model_rank;
    
    // update alignment if best data type changed
    if (best_aln != in_tree->aln) {
        delete in_tree->aln;
        in_tree->aln = best_aln;
        if (best_aln == prot_aln)
            prot_aln = NULL;
        else
            dna_aln = NULL;
    }

    if (dna_aln)
        delete dna_aln;
    if (prot_aln)
        delete prot_aln;

//	in_tree->deleteAllPartialLh();

    string best_tree;
    model_info.getBestTree(best_tree);

    // BQM 2015-07-21 with Lars: load the best_tree
//	if (params.model_test_and_tree)
		in_tree->readTreeString(best_tree);

    
	return at(best_model);
}

int64_t CandidateModelSet::getNextModel() {
    int64_t next_model;
#pragma omp critical
    {
    if (size() == 0)
        next_model = -1;
    else if (current_model == -1)
        next_model = 0;
    else {
        for (next_model = current_model+1; next_model != current_model; next_model++) {
            if (next_model == size())
                next_model = 0;
            if (!at(next_model).hasFlag(MF_IGNORED + MF_WAITING + MF_RUNNING)) {
                break;
            }
        }
    }
    }
    if (next_model != current_model) {
        current_model = next_model;
        at(next_model).setFlag(MF_RUNNING);
        return next_model;
    } else
        return -1;
}

CandidateModel CandidateModelSet::evaluateAll(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
                                    ModelsBlock *models_block, int num_threads, int brlen_type,
                                    string in_model_name, bool merge_phase, bool write_info)
{
    //ModelCheckpoint *checkpoint = &model_info;

    in_tree->params = &params;
    
    Alignment *prot_aln = NULL;
    Alignment *dna_aln = NULL;
    bool do_modelomatic = params.modelomatic && in_tree->aln->seq_type == SEQ_CODON;
    
    
    
    if (in_model_name.empty()) {
        generate(params, in_tree->aln, params.model_test_separate_rate, merge_phase);
        if (do_modelomatic) {
            // generate models for protein
            // adapter coefficient according to Whelan et al. 2015
            prot_aln = in_tree->aln->convertCodonToAA();
            int adjusted_df;
            double adjusted_logl = computeAdapter(in_tree->aln, prot_aln, adjusted_df);
            if (write_info)
                cout << "Adjusted LnL: " << adjusted_logl << "  df: " << adjusted_df << endl;
            size_t start = size();
            generate(params, prot_aln, params.model_test_separate_rate, merge_phase);
            size_t i;
            for (i = start; i < size(); i++) {
                at(i).logl = adjusted_logl;
                at(i).df = adjusted_df;
            }
            
            // generate models for DNA
            dna_aln = in_tree->aln->convertCodonToDNA();
            start = size();
            generate(params, dna_aln, params.model_test_separate_rate, merge_phase);
            for (i = start; i < size(); i++) {
                at(i).setFlag(MF_SAMPLE_SIZE_TRIPLE);
            }
        }
    } else {
        push_back(CandidateModel(in_model_name, "", in_tree->aln));
    }

    if (write_info) {
        cout << "ModelFinder will test " << size() << " ";
        if (do_modelomatic)
            cout << "codon/AA/DNA";
        else
            cout << getSeqTypeName(in_tree->aln->seq_type);
        cout << " models (sample size: " << in_tree->aln->getNSite() << ") ..." << endl;
        cout << " No. Model         -LnL         df  AIC          AICc         BIC" << endl;
    }

    double best_score = DBL_MAX;

    // detect rate hetegeneity automatically or not
    bool auto_rate = merge_phase ? iEquals(params.merge_rates, "AUTO") : iEquals(params.ratehet_set, "AUTO");
    bool auto_subst = merge_phase ? iEquals(params.merge_models, "AUTO") : iEquals(params.model_set, "AUTO");
    int rate_block = size();
    if (auto_rate) {
        for (rate_block = 0; rate_block < size(); rate_block++)
            if (rate_block+1 < size() && at(rate_block+1).subst_name != at(rate_block).subst_name)
                break;
    }
    
    int subst_block = size();
    if (auto_subst) {
        for (subst_block = size()-1; subst_block >= 0; subst_block--)
            if (at(subst_block).rate_name == at(0).rate_name)
                break;
    }

    int64_t num_models = size();
#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads)
#endif
    {
    int64_t model;
    do {
        model = getNextModel();
        if (model == -1)
            break;

        // optimize model parameters
        string orig_model_name = at(model).getName();
        // keep separate output model_info to only update model_info if better model found
        ModelCheckpoint out_model_info;
        at(model).set_name = at(model).aln->name;
        string tree_string;
        
        // main call to estimate model parameters
        tree_string = at(model).evaluate(params, model_info, out_model_info,
                                         models_block, num_threads, brlen_type);
        at(model).computeICScores();
        at(model).setFlag(MF_DONE);
        
        int lower_model = getLowerKModel(model);
        if (lower_model >= 0 && at(lower_model).getScore() < at(model).getScore()) {
            // ignore all +R_k model with higher category
            for (int higher_model = model; higher_model != -1;
                higher_model = getHigherKModel(higher_model)) {
                at(higher_model).setFlag(MF_IGNORED);
            }
            
        }
#ifdef _OPENMP
#pragma omp critical
        {
#endif
        if (best_score > at(model).getScore()) {
            best_score = at(model).getScore();
            if (!tree_string.empty()) {
                //model_info.put("best_tree_" + criterionName(params.model_test_criterion), tree_string);
            }
            // only update model_info with better model
            model_info.putSubCheckpoint(&out_model_info, "");
        }
        model_info.dump();
        if (write_info) {
            cout.width(3);
            cout << right << model+1 << "  ";
            cout.width(13);
            cout << left << at(model).getName() << " ";
            
            cout.precision(3);
            cout << fixed;
            cout.width(12);
            cout << -at(model).logl << " ";
            cout.width(3);
            cout << at(model).df << " ";
            cout.width(12);
            cout << at(model).AIC_score << " ";
            cout.width(12);
            cout << at(model).AICc_score << " " << at(model).BIC_score;
            cout << endl;

        }
        if (model >= rate_block)
            filterRates(model); // auto filter rate models
        if (model >= subst_block)
            filterSubst(model); // auto filter substitution model
#ifdef _OPENMP
        }
#endif
    } while (model != -1);
    }
    
    // store the best model
    ModelTestCriterion criteria[] = {MTC_AIC, MTC_AICC, MTC_BIC};
    for (auto mtc : criteria) {
        int best_model = getBestModelID(mtc);
        model_info.put("best_score_" + criterionName(mtc), at(best_model).getScore(mtc));
        model_info.put("best_model_" + criterionName(mtc), at(best_model).getName());
    }
    
    
    /* sort models by their scores */
    multimap<double,int> model_sorted;
    for (int64_t model = 0; model < num_models; model++)
        if (at(model).hasFlag(MF_DONE)) {
            model_sorted.insert(multimap<double,int>::value_type(at(model).getScore(), model));
        }
    string model_list;
    for (auto it = model_sorted.begin(); it != model_sorted.end(); it++) {
        if (it != model_sorted.begin())
            model_list += " ";
        model_list += at(it->second).getName();
    }
    
    model_info.putBestModelList(model_list);
    model_info.dump();

    // update alignment if best data type changed
    int best_model = getBestModelID(params.model_test_criterion);
    if (at(best_model).aln != in_tree->aln) {
        delete in_tree->aln;
        in_tree->aln = at(best_model).aln;
        if (in_tree->aln == prot_aln)
            prot_aln = NULL;
        else
            dna_aln = NULL;
    }
    
    if (dna_aln)
        delete dna_aln;
    if (prot_aln)
        delete prot_aln;

    return at(best_model);
}


