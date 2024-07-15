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

#if defined(_NN) || defined(_OLD_NN)
#include "nn/neuralnetwork.h"
#endif

// *********for MPI communication*********
// the following defines the communication used for partition process
// note that only SYN_COMM communication is used for merging process
// #define ONESIDE_COMM
#define SYN_COMM

// for one-side communication, how often perform synchronization between the master and the workers
#define TIME_SYN 10 // in seconds

/******* Binary model set ******/
const char* bin_model_names[] = {"GTR2", "JC2"};


/******* Morphological model set ******/
// 2018-08-20: don't test ORDERED model due to lots of numerical issues
//const char* morph_model_names[] = {"MK", "ORDERED"};
const char* morph_model_names[] = {"MK"};


/******* DNA model set ******/
const char* dna_model_names[] = {"JC", "F81", "K80", "HKY", "TNe", "TN",
                                 "K81", "K81u", "TPM2", "TPM2u", "TPM3", "TPM3u",
                                 "TIMe", "TIM", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR"};

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
    return "unknown";
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

/**
 * restrict the number of threads should be used for a single partition
 */
/*
int numThresSinglePart(int nptn, SeqType seqType, int numThres) {
    if (numThres == 0)
        return numThres;
    int thres = numThres;
    switch(seqType) {
        case SEQ_BINARY:
            if (nptn <= 1000)
                thres = 1;
            else if (nptn <= 10000)
                thres = 8;
            break;
        case SEQ_CODON:
            if (nptn*3 <= 10)
                thres = 1;
            else if (nptn*3 <= 50)
                thres = 4;
            else if (nptn*3 <= 1000)
                thres = 8;
            else if (nptn*3 <= 10000)
                thres = 16;
            break;
        case SEQ_PROTEIN:
            if (nptn <= 100)
                thres = 1;
            else if (nptn <= 1000)
                thres = 8;
            else if (nptn <= 10000)
                thres = 16;
            break;
        default: // DNA or others
            if (nptn <= 1000)
                thres = 1;
            else if (nptn <= 10000)
                thres = 4;
            else if (nptn <= 50000)
                thres = 16;
            break;
    }
    if (thres > numThres)
        thres = numThres;
    return thres;
}
*/

/**
 * restrict the number of threads should be used for fast tree estimation
 * and for the final step after modelFinder
 */
int numThresFastTree(int nPart, int nptn, SeqType seqType, int numThres) {
    if (nptn == 1 || numThres == 0)
        return numThres;
    int thres = nPart;
    switch(seqType) {
        case SEQ_PROTEIN:
            if (nptn >= 8000 && nPart < 8)
                thres = 8;
            break;
        case SEQ_CODON:
            if (nptn*3 >= 8000 && nPart < 8)
                thres = 8;
            break;
        default:
            break;
    }
    if (thres > numThres)
        thres = numThres;
    return thres;
}

size_t CandidateModel::getUsualModel(Alignment *aln) {
    size_t aln_len = 0;
    if (aln->isSuperAlignment()) {
        SuperAlignment *super_aln = (SuperAlignment*)aln;
        for (auto it = super_aln->partitions.begin(); it != super_aln->partitions.end(); it++) {
            CandidateModel usual_model(*it);
            if (!subst_name.empty() || !rate_name.empty()) {
                subst_name += ',';
                rate_name += ',';
            }
            subst_name += usual_model.subst_name;
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
    return 0.0;
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

    if (params.constraint_tree_file) {
        iqtree->constraintTree.readConstraint(params.constraint_tree_file, aln->getSeqNames());
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

            // For MPI, we compared between the iqtree objects from all processes and select the optimal one
#ifdef _IQTREE_MPI
            int worker, winner;
            struct {
                double like;
                int rank;
            } in, out;
            if (MPIHelper::getInstance().getNumProcesses() > 1) {
                // find out which process has the maximum likelihood
                in.like = iqtree->getCurScore();
                in.rank = MPIHelper::getInstance().getProcessID();
                MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
                winner = out.rank;
                if (MPIHelper::getInstance().isMaster())
                    cout << "Optimal Fast-NNI log-likelihood: " << out.like << " from process " << winner << endl;

                // the winner sends the checkpoint to others
                if (MPIHelper::getInstance().getProcessID() == winner) {
                    iqtree->saveCheckpoint();
                    iqtree->getModelFactory()->saveCheckpoint();
                    for (worker=0; worker<MPIHelper::getInstance().getNumProcesses(); worker++) {
                        if (worker==winner)
                            continue;
                        MPIHelper::getInstance().sendCheckpoint(iqtree->getCheckpoint(), worker);
                    }
                } else {
                    // receive the checkpoint from the winner
                    // and update the iqtree object
                    MPIHelper::getInstance().recvCheckpoint(iqtree->getCheckpoint());
                    iqtree->restoreCheckpoint();
                    iqtree->getModelFactory()->restoreCheckpoint();
                    iqtree->initializeAllPartialLh();
                    iqtree->setCurScore(iqtree->computeLikelihood());
                }
            }
#endif
            
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

void runModelFinder(Params &params, IQTree &iqtree, ModelCheckpoint &model_info, string &best_subst_name, string &best_rate_name, map<string, vector<string> > nest_network, bool under_mix_finder)
{
    // if it is an alignment with partitions and the number of threads is more than the number of alignments,
    // then set the number of threads = the number of alignments (for most of the cases)
    bool autoThread = (params.num_threads == 0);
    int orig_nthreads = params.num_threads;
    int updated_nthreads = params.num_threads;
    if (!autoThread && iqtree.isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
        updated_nthreads = numThresFastTree(stree->size(), stree->at(0)->aln->getNPattern(), stree->at(0)->aln->seq_type, orig_nthreads);
        params.num_threads = updated_nthreads;
    }
    if (updated_nthreads != orig_nthreads) {
        cout << "The number of threads is changed to: " << updated_nthreads << endl;
    }
    
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
    // if (MPIHelper::getInstance().getNumProcesses() > 1)
    //    outError("Please use only 1 MPI process! We are currently working on the MPI parallelization of model selection.");
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
    
    // after loading, workers are not allowed to write checkpoint anymore
    if (MPIHelper::getInstance().isWorker())
        model_info.setFileName("");
    
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
    // cout << "params.modelfinder_ml_tree = " << params.modelfinder_ml_tree << endl << flush;
    if (!params.use_nn_model && params.modelfinder_ml_tree) {
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
    
    if (!autoThread)
        params.num_threads = orig_nthreads;
    
    if (!params.use_nn_model) {
        // also save initial tree to the original .ckp.gz checkpoint
        //        string initTree = iqtree.getTreeString();
        //        CKP_SAVE(initTree);
        //        iqtree.saveCheckpoint();
        //        checkpoint->dump(true);
        
        CandidateModelSet candidate_models;
        int max_cats = candidate_models.generate(params, iqtree.aln, params.model_test_separate_rate, false);
        
        // If the option -m MF1 is used, consider ALL candidates
        if (params.model_name == "MF1") {
            params.score_diff_thres = -1.0;
            cout << "ModelFinder 1 is activated" << endl;
        }
        
        uint64_t mem_size = iqtree.getMemoryRequiredThreaded(max_cats);
        cout << "NOTE: ModelFinder requires " << (mem_size / 1024) / 1024 << " MB RAM!" << endl;
        if (mem_size >= getMemorySize()) {
            outError("Memory required exceeds your computer RAM size!");
        }
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
        CandidateModelSet model_set;
        model_set.nest_network = nest_network;
        model_set.under_mix_finder = under_mix_finder;
        Checkpoint *checkpoint = &model_info;
        // neural network model selection (added by TD)
#if defined(_NN) || defined(_OLD_NN)
        if (params.use_nn_model) {
            cout << "We are using the neural network to select the model of sequence evolution because "
            "option --use-nn-model is set to " << params.use_nn_model << endl;
            Alignment *alignment = (iqtree.aln->removeAndFillUpGappySites())->replaceAmbiguousChars();
            NeuralNetwork nn(alignment);
            iqtree.aln->model_name = nn.doModelInference();
            best_subst_name = iqtree.aln->model_name;
            best_rate_name = "";
            double alpha = nn.doAlphaInference();
            if (alpha >= 0) { // +G
                iqtree.aln->model_name += "+G{" + to_string(alpha) + "}";
                best_rate_name = "G"; // to confirm, G or +G
            } 
            string best_model_NN;
            CKP_RESTORE(best_model_NN);
            delete alignment;
            
            cout << "Best-fit model: " << iqtree.aln->model_name << " chosen according to neural network" << endl;
        } else {
#endif
            if (params.openmp_by_model)
                best_model = model_set.evaluateAll(params, &iqtree,
                                                             model_info, models_block, params.num_threads,
                                                             BRLEN_OPTIMIZE);
            else
                best_model = model_set.test(params, &iqtree,
                                                      model_info, models_block, params.num_threads, BRLEN_OPTIMIZE);
            iqtree.aln->model_name = best_model.getName();
            best_subst_name = best_model.subst_name;
            best_rate_name = best_model.rate_name;
            // Checkpoint *checkpoint = &model_info;
            string best_model_AIC, best_model_AICc, best_model_BIC;
            CKP_RESTORE(best_model_AIC);
            CKP_RESTORE(best_model_AICc);
            CKP_RESTORE(best_model_BIC);
            cout << "Akaike Information Criterion:           " << best_model_AIC << endl;
            cout << "Corrected Akaike Information Criterion: " << best_model_AICc << endl;
            cout << "Bayesian Information Criterion:         " << best_model_BIC << endl;
            
            cout << "Best-fit model: " << iqtree.aln->model_name << " chosen according to "
            << criterionName(params.model_test_criterion) << endl;
#if defined(_NN) || defined(_OLD_NN)
        }
#endif
    }

    // remove key "OptModel" from the checkpoint file, which is only used for initialising models from the nested models.
    iqtree.getCheckpoint()->eraseKeyPrefix("OptModel");

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

    if (!autoThread) {
        if (iqtree.isSuperTree()) {
            // change the number of threads to the number of partitions for alignment with partitions
            PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
            params.num_threads = stree->size();
            cout << "The number of threads is changed to: " << params.num_threads << endl;
        } else {
            params.num_threads = updated_nthreads;
        }
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
            model_names = reorderModelNames(model_names);
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

    // put "FO" to the last
    vector<string>::iterator itr = remove(freq_names.begin(), freq_names.end(), "+FO");
    if (itr != freq_names.end()) {
        freq_names.erase(itr, freq_names.end());
        freq_names.push_back("+FO");
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
    // if (!iqtree->getModel()->isMixture() || in_aln->seq_type == SEQ_POMO) {
        subst_name = iqtree->getSubstName();
        rate_name = iqtree->getRateName();
    // }


    if (restoreCheckpoint(&in_model_info)) {
        delete iqtree;
        return "";
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    iqtree->getModelFactory()->restoreCheckpoint();
    
    bool rate_restored = iqtree->getRate()->hasCheckpoint();
    
    // now switch to the output checkpoint
    iqtree->getModelFactory()->setCheckpoint(&out_model_info);
    iqtree->setCheckpoint(&out_model_info);

    double new_logl;

    if (syncChkPoint != nullptr)
        iqtree->getModelFactory()->syncChkPoint = this->syncChkPoint;
    
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
        

        if (init_first_mix) {

            // now switch to the input checkpoint
            iqtree->getModelFactory()->setCheckpoint(&in_model_info);
            iqtree->setCheckpoint(&in_model_info);

            // get the model mixture object
            ModelMixture* modelmix = dynamic_cast<ModelMixture*> (iqtree->getModelFactory()->model);
            ASSERT(modelmix);
            double init_weight = 1.0 / modelmix->getNMixtures();
            
            // obtain the likelihood value from the (k-1)-class mixture model
            string criteria_str = criterionName(params.model_test_criterion);
            string best_model = in_model_info["best_model_" + criteria_str];
            string best_model_logl_df = in_model_info[best_model];
            stringstream ss (best_model_logl_df);
            double pre_logl;
            ss >> pre_logl;
            
            for (int step = 0; step < 10; step++) {
                
                // initialize the parameters from the (k-1)-class mixture model
                modelmix->initFromClassMinusOne(init_weight);
                
                new_logl = iqtree->getModelFactory()->optimizeParameters(brlen_type, false,
                                                                         params.modelfinder_eps, TOL_GRADIENT_MODELTEST);
                
                // check if new logl is worse than logl from the (k-1)-class mixture model
                if (pre_logl < new_logl + params.modelfinder_eps) break;
                init_weight *= 0.5;
                cout << getName() << " reinitialized from the previous (k-1)-class mixture model with initial weight: " << init_weight << endl;
            }
            tree_len = iqtree->treeLength();

            // now switch to the output checkpoint
            iqtree->getModelFactory()->setCheckpoint(&out_model_info);
            iqtree->setCheckpoint(&out_model_info);

            iqtree->getModelFactory()->saveCheckpoint();
            iqtree->saveCheckpoint();
            if (new_logl < pre_logl - params.modelfinder_eps*10.0) {
                outWarning("Log-likelihood " + convertDoubleToString(new_logl) + " of " +
                           getName() + " worse than the previous (k-1)-class mixture model " + convertDoubleToString(pre_logl));
            }
        } else {
            CandidateModel prev_info;
            bool prev_rate_present = prev_info.restoreCheckpointRminus1(&in_model_info, this);

            if (!prev_rate_present){
                iqtree->getModelFactory()->setCheckpoint(&in_model_info);
                iqtree->getModelFactory()->initFromNestedModel(nest_network);

                new_logl = iqtree->getModelFactory()->optimizeParameters(brlen_type, false,
                                                                         params.modelfinder_eps, TOL_GRADIENT_MODELTEST);
                tree_len = iqtree->treeLength();

                // now switch to the output checkpoint
                iqtree->getModelFactory()->setCheckpoint(&out_model_info);
                iqtree->setCheckpoint(&out_model_info);

                iqtree->getModelFactory()->saveCheckpoint();
                iqtree->saveCheckpoint();

            } else {
                // try to initialise +R[k+1] from +R[k] if not restored from checkpoint
                double weight_rescale = 1.0;
                if (!rate_restored) {
                    iqtree->getRate()->initFromCatMinusOne(in_model_info, weight_rescale);
                    if (verbose_mode >= VB_MED)
                        cout << iqtree->getRate()->name << " initialized from " << prev_info.rate_name << endl;
                }
                for (int step = 0; step < 5; step++) {
                    new_logl = iqtree->getModelFactory()->optimizeParameters(brlen_type, false,
                                                                             params.modelfinder_eps,
                                                                             TOL_GRADIENT_MODELTEST);
                    tree_len = iqtree->treeLength();
                    iqtree->getModelFactory()->saveCheckpoint();
                    iqtree->saveCheckpoint();

                    // check if logl(+R[k]) is worse than logl(+R[k-1])
                    // if (!prev_rate_present) break;
                    if (prev_info.logl < new_logl + params.modelfinder_eps) break;
                    weight_rescale *= 0.5;
                    iqtree->getRate()->initFromCatMinusOne(in_model_info, weight_rescale);
                    cout << iqtree->getRate()->name << " reinitialized from " << prev_info.rate_name
                         << " with factor " << weight_rescale << endl;
                }
                if (prev_rate_present && new_logl < prev_info.logl - params.modelfinder_eps * 10.0) {
                    outWarning("Log-likelihood " + convertDoubleToString(new_logl) + " of " +
                               getName() + " worse than " + prev_info.getName() + " " +
                               convertDoubleToString(prev_info.logl));
                }
            }
        }
    }
    // sum in case of adjusted df and logl already stored
    df += iqtree->getModelFactory()->getNParameters(brlen_type);
    logl += new_logl;
    string tree_string = iqtree->getTreeString();

    if (verbose_mode >= VB_MED) {
        cout << "[optimized] " << iqtree->getModelFactory()->model->getNameParams(false) << endl;
    }

    if (syncChkPoint != nullptr)
        iqtree->getModelFactory()->syncChkPoint = nullptr;

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

/*
bool comparePartition(const pair<int,double> &a, const pair<int, double> &b) {
    return a.second > b.second;
}
*/

bool compareJob(const pair<int,double> &a, const pair<int, double> &b) {
    return (a.second == b.second)?(a.first < b.first):(a.second > b.second);
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

struct workloadcmp {
  bool operator() (const pair<int,double>& a, const pair<int,double>& b) const
    {return (a.second==b.second)?a.first<b.first:a.second<b.second;}
};

void replaceModelInfo(ModelCheckpoint* model_info, ModelCheckpoint &new_info) {
    for (auto it = new_info.begin(); it != new_info.end(); it++) {
        model_info->put(it->first, it->second);
    }
}

/**
 * select models for all partitions
 * @param[in,out] model_info (IN/OUT) all model information
 * @return total number of parameters
 */
/*
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
    int step = 0;
    double pre_inf_score;
    string blkStr = string(80,' ');
    if (params.partition_type == BRLEN_FIX || params.partition_type == BRLEN_SCALE) {
        dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
        if (params.partition_type == BRLEN_SCALE)
            dfsum -= 1;
    }
	size_t  ssize = in_tree->getAlnNSite();
	int64_t num_model = 0;
    int64_t total_num_model = in_tree->size();
    
    // get the name of the algorithm
    string part_algo = "";
    if (params.partition_merge == MERGE_GREEDY)
        part_algo = "Greedy Algorithm";
    else if (params.partition_merge == MERGE_RCLUSTER)
        part_algo = "Relaxed Algorithm";
    else if (params.partition_merge == MERGE_RCLUSTERF)
        part_algo = "Fast Relaxed Algorithm";
    else if (params.partition_merge == MERGE_KMEANS)
        part_algo = "Kmean Algorithm";

    // for greedy algorithm
    if (params.partition_merge == MERGE_GREEDY) {
        params.partfinder_rcluster_max = in_tree->size() * (in_tree->size()-1) / 2;
        params.partfinder_log_rate = false;
        params.partfinder_rcluster = 100.0;
    }
    
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
	// cout << " No. Model        Score       TreeLen     Charset" << endl;

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
//            cout.width(4);
//            cout << right << num_model << " ";
//            cout.width(12);
//            cout << left << best_model.getName() << " ";
//            cout.width(11);
//            cout << score << " ";
//            cout.width(11);
//            cout << best_model.tree_len << " ";
//            cout << this_tree->aln->name;
//            if (num_model >= 10) {
//                double remain_time = (total_num_model-num_model)*(getRealTime()-start_time)/num_model;
//                double finish_percent = (double) num_model * 100.0 / total_num_model;
//                cout << "Finished subset " << num_model << "/" << total_num_model << "\t" << finish_percent << " percent done";
//                cout << "\t" << convert_time(getRealTime()-start_time) << " ("
//                    << convert_time(remain_time) << " left)\r";
//                cout << flush;
//            }
//            cout << endl;
            replaceModelInfo(this_tree->aln->name, model_info, part_model_info);
            model_info.dump();
        }
    }

    // in case ModelOMatic change the alignment
    fixPartitions(in_tree);
    
	double inf_score = computeInformationScore(lhsum, dfsum, ssize, params.model_test_criterion);
	cout << "Full partition model " << criterionName(params.model_test_criterion)
         << " score: " << inf_score << " (LnL: " << lhsum << "  df:" << dfsum << ")" << endl;

    pre_inf_score = inf_score;

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

    // following implements the greedy algorithm of Lanfear et al. (2012)
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
        size_t compute_pairs = 0;

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
                    compute_pairs++;
//					cout.width(4);
//					cout << right << num_model << " ";
//					cout.width(12);
//					cout << left << best_model.getName() << " ";
//					cout.width(11);
//                    cout << cur_pair.score << " ";
//                    cout.width(11);
//                    cout << cur_pair.tree_len << " " << cur_pair.set_name;
                    if (num_model >= 10) {
                        double remain_time = max(total_num_model-num_model, (int64_t)0)*(getRealTime()-start_time)/num_model;
                        if (remain_time < 0.0)
                            remain_time = 0.0;
                        double finish_percent = (double) pair * 100.0 / num_pairs;
                        cout << " Finished subset " << pair << "/" << num_pairs << "     " << finish_percent << "  percent done";
                        cout << "     " << convert_time(getRealTime()-start_time) << " ("
                            << convert_time(remain_time) << " left)     \r" << flush;
                    }
//                    cout << endl;

				}
                if (cur_pair.score < inf_score)
                    better_pairs.insertPair(cur_pair);
			}

        }

        // clear the message previous on this line
        cout << blkStr << "\r" << flush;

        if (better_pairs.size() > 0) {
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
                
                // cout << "Merging " << opt_pair.set_name << " with " << criterionName(params.model_test_criterion)
                //     << " score: " << inf_score << " (LnL: " << lhsum << "  df: " << dfsum << ")" << endl;
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
        
        cout << "ModelFinder2\t";
        if (part_algo.length() > 0)
            cout << part_algo << "\t";
        cout << "Step " << ++step << "\t" << compute_pairs << " Subsets\t" << criterionName(params.model_test_criterion) << " " << inf_score;
        cout << "\tdeltaBIC " << inf_score - pre_inf_score;
        cout << endl;
        pre_inf_score = inf_score;

        if (better_pairs.empty()) break;
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

	// cout << "Agglomerative model selection: " << final_model_tree << endl;
    
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

        cout << endl;
        cout << "No. Model        Score       Charset" << endl;
        int partition_id = 0;

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
            cout << right << ++partition_id << " ";
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
*/

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
    string set_name, string in_model_name, bool merge_phase,
    bool generate_candidates, bool skip_all_when_drop)
{

    ModelCheckpoint *checkpoint = &model_info;

	    in_tree->params = &params;
    
    // for ModelOMatic
    Alignment *prot_aln = NULL;
    Alignment *dna_aln = NULL;
    bool do_modelomatic = params.modelomatic && in_tree->aln->seq_type == SEQ_CODON;
    if (generate_candidates) {
        if (in_model_name.empty()) {
#if defined(_NN) || defined(_OLD_NN)
            if (params.use_nn_model && in_tree->aln->seq_type == SEQ_DNA) {
                cout << "Using NN" << endl;
                // todo: to work with multi-threading: pass along the random number streams to the rngs in the stochastic functions
                // determine substitution model using neural network
                Alignment *alignment = (in_tree->aln->removeAndFillUpGappySites())->replaceAmbiguousChars(); // todo: here
                NeuralNetwork nn(alignment);
                string model_name = nn.doModelInference(); // todo: here
                string rate_name = "";
                double alpha = nn.doAlphaInference(); // todo: here
                if (alpha >= 0) { // +G
                    rate_name = "+G{" + to_string(alpha) + "}";
                }
                string best_model_NN;
                CKP_RESTORE(best_model_NN);
                delete alignment;
                push_back(CandidateModel(model_name, rate_name, in_tree->aln));
            } else {
#endif
                // generate all models the normal way
                generate(params, in_tree->aln, params.model_test_separate_rate, merge_phase);
#if defined(_NN) || defined(_OLD_NN)
            }
            if (do_modelomatic) {
                ASSERT(!params.use_nn_model);
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
#endif
        } else {
            push_back(CandidateModel(in_model_name, "", in_tree->aln));
        }
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
        cout << " models (sample size: " << ssize << " epsilon: " << params.modelfinder_eps << ") ..." << endl;
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
        at(model).nest_network = nest_network;
        /***** main call to estimate model parameters ******/
        at(model).syncChkPoint = this->syncChkPoint;
        tree_string = at(model).evaluate(params,
            model_info, out_model_info, models_block, num_threads, brlen_type);
        at(model).syncChkPoint = nullptr;
        at(model).computeICScores(ssize);
        at(model).setFlag(MF_DONE);

        // for testing
        CandidateModel prev_info;

        bool skip_model = false;
        bool skip_all_models = false;
        
        bool check_condition = prev_info.restoreCheckpointRminus1(checkpoint, &at(model));

        if (check_condition) {
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

        if (skip_all_when_drop && model>0) {
            switch (params.model_test_criterion) {
            case MTC_ALL:
                if (at(model).AIC_score > at(model-1).AIC_score &&
                    at(model).AICc_score > at(model-1).AICc_score &&
                    at(model).BIC_score > at(model-1).BIC_score) {
                    // skip all remaining models
                    skip_all_models = true;
                }
                break;
            case MTC_AIC:
                if (at(model).AIC_score > at(model-1).AIC_score) {
                    // skip all remaining models
                    skip_all_models = true;
                }
                break;
            case MTC_AICC:
                if (at(model).AICc_score > at(model-1).AICc_score) {
                    // skip all remaining models
                    skip_all_models = true;
                }
                break;
            case MTC_BIC:
                if (at(model).BIC_score > at(model-1).BIC_score) {
                    // skip all remaining models
                    skip_all_models = true;
                }
                break;
            }
        }

        // BQM 2024-06-22: save checkpoint for starting values of next model
        model_info.putSubCheckpoint(&out_model_info, "");
        
        bool is_better_model = false;

		if (at(model).AIC_score < best_score_AIC) {
            best_model_AIC = model;
            best_score_AIC = at(model).AIC_score;
            if (!tree_string.empty())
                best_tree_AIC = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_AIC) {
                //model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
                is_better_model = true;
            }
        }
		if (at(model).AICc_score < best_score_AICc) {
            best_model_AICc = model;
            best_score_AICc = at(model).AICc_score;
            if (!tree_string.empty())
                best_tree_AICc = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_AICC) {
                //model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
                is_better_model = true;
            }
        }

		if (at(model).BIC_score < best_score_BIC) {
			best_model_BIC = model;
            best_score_BIC = at(model).BIC_score;
            if (!tree_string.empty())
                best_tree_BIC = tree_string;
            // only update model_info with better model
            if (params.model_test_criterion == MTC_BIC) {
                //model_info.putSubCheckpoint(&out_model_info, "");
                best_aln = at(model).aln;
                is_better_model = true;
            }
        }

        model_info.startStruct("OptModel");
        model_info.putSubCheckpoint(&out_model_info, at(model).getName());
        model_info.endStruct();
        
        if (under_mix_finder && is_better_model) {
            model_info.putSubCheckpoint(&out_model_info, "BestOfTheKClass");
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

        if (skip_all_models) {
            // skip over all the remaining models
            for (int next = model+1; next < size(); next++) {
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


/**
 * select models for all partitions
 * @param[in,out] model_info (IN/OUT) all model information
 * @return total number of parameters
 */
void testPartitionModel(Params &params, PhyloSuperTree* in_tree, ModelCheckpoint &model_info,
    ModelsBlock *models_block, int num_threads)
{
    PartitionFinder partitionFinder(&params, in_tree, &model_info, models_block, num_threads);
    partitionFinder.test_PartitionModel();
}

struct jobcomp {
    // for sorting the jobs
    bool operator() (const int& i1, const int& i2) const {return i1<i2;}
};

/* Constructor
 */
PartitionFinder::PartitionFinder(Params *inparams, PhyloSuperTree* intree, ModelCheckpoint *modelinfo,
                                 ModelsBlock *modelsblock, int numthreads) {
 
    params = inparams;
    in_tree = intree;
    model_info = modelinfo;
    models_block = modelsblock;
    num_threads = numthreads;
    num_processes = MPIHelper::getInstance().getNumProcesses();
}

/*
 * Show the the other worker's result of best model for the merge
 */
void PartitionFinder::showMergeResult(ModelCheckpoint& part_model_info, double tree_len, const string& model_name, string& set_name, bool done_before, int tag) {

    double remain_time;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        if (!done_before) {
            replaceModelInfo(set_name, *model_info, part_model_info);
            model_info->dump();

            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            if (tag != -1)
                cout << tag << " ";
            cout.width(12);
            cout << left << model_name << " ";
            // cout.width(11);
            // cout << score << " ";
            cout.width(11);
            cout << tree_len << " " << set_name;
            if (num_model >= 10) {
                remain_time = max(total_num_model-num_model, (int64_t)0)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " ("
                     << convert_time(remain_time) << " left)";
            }
            cout << endl;
        }
        
        // update the number of jobs done
        jobdone++;
    }
}

/*
 * Show the the other worker's result of best model for the merge
 */
void PartitionFinder::showMergeResults(ModelCheckpoint& part_model_info, vector<double>& tree_len, vector<string>& model_name, vector<string>& set_name, vector<int>& tag, int tot_jobs_done) {

    double remain_time;
    int i;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        replaceModelInfo(model_info, part_model_info);
        model_info->dump();
        
        for (i=0; i<tree_len.size(); i++) {
            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            if (tag[i] != -1)
                cout << tag[i] << " ";
            cout.width(12);
            cout << left << model_name[i] << " ";
            // cout.width(11);
            // cout << score[i] << " ";
            cout.width(11);
            cout << tree_len[i] << " " << set_name[i];
            if (num_model >= 10) {
                remain_time = max(total_num_model-num_model, (int64_t)0)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " (" << convert_time(remain_time) << " left)";
            }
            cout << endl;
        }

        // update the number of jobs done
        jobdone += tot_jobs_done;
    }

}

#ifdef _IQTREE_MPI

/**
 * Process the computation of the best model for a merge with MPI
 * Find the best model for merging two partitions in job_id
 * @param job_id ID of closest_pairs array
 * @param nthreads number of threads available for this job
 * @param need_next_jobID whether it is needed to get the next tree ID

 * @return next job ID if need_next_treeID and (MASTER or IS_ASYN_COMM = 0), otherwise -1
 */
void PartitionFinder::getBestModelForOneMergeMPI(MergeJob* job, int nthreads, bool need_next_jobID, SyncChkPoint& syncChkPt, double& run_time, double& wait_time) {
    
    CandidateModel best_model;
    ModelPair cur_pair;
    ModelCheckpoint part_model_info;
    CandidateModelSet candModelSet;
    double weight1, weight2, sum, cur_tree_len;
    bool under_mpi, done_before, check;
    string key;
    double lhnew;
    int dfnew;
    int next_job_id = -1;
    bool noChkMessage = false;
    int job_type = 2; // compute the best model for the merge
    double t_begin;
    double t_wait_begin;
    
    t_begin = getRealTime();
    wait_time = 0;
    part_model_info.clear();

    // information of current partitions pair
    cur_pair.part1 = job->id1;
    cur_pair.part2 = job->id2;
    ASSERT(cur_pair.part1 < cur_pair.part2);
    cur_pair.merged_set.insert(job->geneset1.begin(), job->geneset1.end());
    cur_pair.merged_set.insert(job->geneset2.begin(), job->geneset2.end());
    cur_pair.set_name = getSubsetName(in_tree, cur_pair.merged_set);
    weight1 = getSubsetAlnLength(in_tree, job->geneset1);
    weight2 = getSubsetAlnLength(in_tree, job->geneset2);
    sum = 1.0 / (weight1 + weight2);
    weight1 *= sum;
    weight2 *= sum;
    done_before = false;

    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        // if pairs previously examined, reuse the information
        model_info->startStruct(cur_pair.set_name);
        if (model_info->getBestModel(best_model.subst_name)) {
            best_model.restoreCheckpoint(model_info);
            done_before = true;
        }
        model_info->endStruct();
    }

    cur_tree_len = 0.0;
    if (!done_before) {
        Alignment *aln = super_aln->concatenateAlignments(cur_pair.merged_set);
        PhyloTree *tree = in_tree->extractSubtree(cur_pair.merged_set);
        //tree->scaleLength((weight1*lenvec[cur_pair.part1] + weight2*lenvec[cur_pair.part2])/tree->treeLength());
        tree->scaleLength(sqrt(job->treelen1*job->treelen2)/tree->treeLength());
        cur_tree_len = tree->treeLength();
        tree->setAlignment(aln);
#ifdef _OPENMP
#pragma omp critical
#endif
{
        extractModelInfo(cur_pair.set_name, *model_info, part_model_info);
        transferModelParameters(in_tree, *model_info, part_model_info, job->geneset1, job->geneset2);
}
        tree->num_precision = in_tree->num_precision;
        tree->setParams(params);
        tree->sse = params->SSE;
        tree->optimize_by_newton = params->optimize_by_newton;
        tree->setNumThreads(params->model_test_and_tree ? num_threads : 1);

        tree->setCheckpoint(&part_model_info);
        // trick to restore checkpoint
        tree->restoreCheckpoint();
        tree->saveCheckpoint();

        candModelSet.syncChkPoint = &(syncChkPt);
        best_model = candModelSet.test(*params, tree, part_model_info, models_block,
                                       nthreads, params->partition_type, cur_pair.set_name, "", true);
        candModelSet.syncChkPoint = nullptr;
        check = (best_model.restoreCheckpoint(&part_model_info));
        ASSERT(check);
        delete tree;
        delete aln;
    }

    cur_pair.logl = best_model.logl;
    cur_pair.df = best_model.df;
    cur_pair.model_name = best_model.getName();
    cur_pair.tree_len = best_model.tree_len;

    if (MPIHelper::getInstance().isMaster()) {
        // for Master
        showMergeResult(part_model_info, cur_pair.tree_len, cur_pair.model_name, cur_pair.set_name, done_before, syncChkPt.mytag);
        if (need_next_jobID) {
            t_wait_begin = getRealTime();
            // next_job_id = syncChkPt.getNextJobID();
            syncChkPt.getNextMergeJob(job);
            wait_time += (getRealTime() - t_wait_begin);
        }
        // collect the answers from workers
        syncChkPt.masterSyncOtherChkpts(true);
    } else {
        // for Worker
        key = "pf_tree_len"; part_model_info.put(key, cur_pair.tree_len);
        key = "pf_model_name"; part_model_info.put(key, cur_pair.model_name);
        key = "pf_set_name"; part_model_info.put(key, cur_pair.set_name);
        key = "pf_done_before"; part_model_info.putBool(key, done_before);

        // send the part_model_info to master if time is long enough
        t_wait_begin = getRealTime();
        next_job_id = syncChkPt.sendChkptToMaster(part_model_info, need_next_jobID, job_type, job);
        wait_time += (getRealTime() - t_wait_begin);
    }
    run_time = (getRealTime() - t_begin) - wait_time;
}

#endif


/*
 * Show the result of best model for the partition
 */
void PartitionFinder::showPartitionResult(ModelCheckpoint& part_model_info, int tree_id, double tree_len, const string& model_name, double score, int tag) {

    PhyloTree *this_tree = in_tree->at(tree_id);
    double remain_time;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        num_model++;
        cout.width(4);
        cout << right << num_model << " ";
        if (tag != -1)
            cout << tag << " ";
        cout.width(12);
        cout << left << model_name << " ";
        cout.width(11);
        cout << score << " ";
        cout.width(11);
        cout << tree_len << " ";
        cout << this_tree->aln->name;
        if (num_model >= 10) {
            remain_time = (double)(total_num_model-num_model)*(getRealTime()-start_time)/num_model;
            cout << "\t" << convert_time(getRealTime()-start_time) << " ("
                 << convert_time(remain_time) << " left)";
        }
        cout << endl;
        replaceModelInfo(this_tree->aln->name, *model_info, part_model_info);
        model_info->dump();

        // update the number of jobs done
        jobdone++;
    }
}

/*
 * Show a set of best-model results for the partition
 */
void PartitionFinder::showPartitionResults(ModelCheckpoint& part_model_info, vector<int>& tree_id, vector<double>& tree_len, vector<string>& model_name, vector<double>& score, vector<int>& tag) {

    PhyloTree *this_tree;
    double remain_time;
    int i;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        replaceModelInfo(model_info, part_model_info);
        model_info->dump();

        for (i=0; i<tree_id.size(); i++) {
            this_tree = in_tree->at(tree_id[i]);
            num_model++;
            cout.width(4);
            cout << right << num_model << " ";
            if (tag[i] != -1)
                cout << tag[i] << " ";
            cout.width(12);
            cout << left << model_name[i] << " ";
            cout.width(11);
            cout << score[i] << " ";
            cout.width(11);
            cout << tree_len[i] << " ";
            cout << this_tree->aln->name;
            if (num_model >= 10) {
                remain_time = (double)(total_num_model-num_model)*(getRealTime()-start_time)/num_model;
                cout << "\t" << convert_time(getRealTime()-start_time) << " ("
                     << convert_time(remain_time) << " left)";
            }
            cout << endl;
        }

        // update the number of jobs done
        jobdone += tree_id.size();
    }
}

#ifdef _IQTREE_MPI

/**
 * Process the computation of the best model for a single partition with MPI
 *
 * nthreads : number of threads available for this job
 * need_next_treeID : whether it is needed to get the next tree ID
 *
 * if need_next_treeID and (MASTER or IS_ASYN_COMM = 0)
 *    return the next Job ID from master
 * else
 *    return -1
 */
int PartitionFinder::computeBestModelforOnePartitionMPI(int tree_id, int nthreads, bool need_next_treeID, SyncChkPoint& syncChkPt, double& run_time, double& wait_time) {
    
    CandidateModel best_model;
    PhyloTree *this_tree;
    ModelCheckpoint part_model_info;
    CandidateModelSet candModelSet;
    string part_model_name, key;
    bool under_mpi, check;
    double score, remain_time;
    ostringstream ss;
    int next_tree_id = -1;
    bool noChkMessage = false;
    int job_type = 1; // compute the best model for partition
    double t_begin;
    double t_wait_begin;
    
    t_begin = getRealTime();
    wait_time = 0;
    this_tree = in_tree->at(tree_id);

#ifdef _OPENMP
#pragma omp critical
#endif
{
    extractModelInfo(this_tree->aln->name, *model_info, part_model_info);
}

    // do the computation
    if (params->model_name.empty())
        part_model_name = this_tree->aln->model_name;
    
    candModelSet.syncChkPoint = &(syncChkPt);
    
    best_model = candModelSet.test(*params, this_tree, part_model_info, models_block,
        nthreads, brlen_type, this_tree->aln->name, part_model_name, test_merge);
    
    candModelSet.syncChkPoint = nullptr;
    
    check = (best_model.restoreCheckpoint(&part_model_info));
    ASSERT(check);
    score = best_model.computeICScore(this_tree->getAlnNSite());
    this_tree->aln->model_name = best_model.getName();

    if (MPIHelper::getInstance().isMaster()) {
        // for Master

        showPartitionResult(part_model_info, tree_id, best_model.tree_len, best_model.getName(), score, syncChkPt.mytag);
        if (need_next_treeID) {
            t_wait_begin = getRealTime();
            next_tree_id = syncChkPt.getNextJobID();
            wait_time += (getRealTime() - t_wait_begin);
        }

    } else {
        
        // for Worker -- SYN communication
        #ifdef SYN_COMM
        
            key = "pf_tree_id"; part_model_info.put(key, tree_id);
            key = "pf_tree_len"; part_model_info.put(key, best_model.tree_len);
            key = "pf_model_name"; part_model_info.put(key, best_model.getName());
            key = "pf_score"; part_model_info.put(key, score);
        
            // send the part_model_info to master if time is long enough
            t_wait_begin = getRealTime();
            next_tree_id = syncChkPt.sendChkptToMaster(part_model_info, need_next_treeID, job_type);
            wait_time += (getRealTime() - t_wait_begin);

        #endif // SYN_COMM

        // for Worker -- ONESIDE communication
        #ifdef ONESIDE_COMM

            // consolidate part_model_info into the process_model_info
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                replaceModelInfo(this_tree->aln->name, process_model_info, part_model_info);
                tree_id_vec.push_back(tree_id);
                tree_len_vec.push_back(best_model.tree_len);
                model_name_vec.push_back(best_model.getName());
                score_vec.push_back(score);
                tag_vec.push_back(syncChkPt.mytag);
            }
            
            // send the process_model_info to master if time is long enough
            t_wait_begin = getRealTime();
            next_tree_id = syncChkPt.sendChkptToMaster(process_model_info, need_next_treeID, job_type);
            wait_time += (getRealTime() - t_wait_begin);
        
        #endif // ONESIDE_COMM
    }
    run_time = (getRealTime() - t_begin) - wait_time;
    
    return next_tree_id;
}

#endif // _IQTREE_MPI

// retreive the answers from checkpoint (for merging)
// and remove those jobs from the array jobIDs
void PartitionFinder::retreiveAnsFrChkpt(vector<pair<int,double> >& jobs) {

    CandidateModel best_model;
    vector<char> to_delete;

    // for merging partitions
    for (int j = 0; j < jobs.size(); j++) {
        ModelPair cur_pair;
        double lhnew;
        int dfnew;
        int job_type = 2; // compute the best model for the merge
        double t_begin;
        double t_wait_begin;
        int pair = jobs[j].first;

        // information of current partitions pair
        cur_pair.part1 = closest_pairs[pair].first;
        cur_pair.part2 = closest_pairs[pair].second;
        ASSERT(cur_pair.part1 < cur_pair.part2);
        cur_pair.merged_set.insert(gene_sets[cur_pair.part1].begin(), gene_sets[cur_pair.part1].end());
        cur_pair.merged_set.insert(gene_sets[cur_pair.part2].begin(), gene_sets[cur_pair.part2].end());
        cur_pair.set_name = getSubsetName(in_tree, cur_pair.merged_set);
        
        // check whether the pair was previously examined, reuse the information
        model_info->startStruct(cur_pair.set_name);
        if (model_info->getBestModel(best_model.subst_name)) {
            best_model.restoreCheckpoint(model_info);
            cur_pair.logl = best_model.logl;
            cur_pair.df = best_model.df;
            cur_pair.model_name = best_model.getName();
            cur_pair.tree_len = best_model.tree_len;
            lhnew = lhsum - lhvec[cur_pair.part1] - lhvec[cur_pair.part2] + best_model.logl;
            dfnew = dfsum - dfvec[cur_pair.part1] - dfvec[cur_pair.part2] + best_model.df;
            cur_pair.score = computeInformationScore(lhnew, dfnew, ssize, params->model_test_criterion);
            to_delete.push_back(1);
        } else {
            to_delete.push_back(0);
        }
        model_info->endStruct();
    }
    
    // remove the finished jobs from the list
    int k = 0;
    for (int j = 0; j < jobs.size(); j++) {
        if (!to_delete[j]) {
            if (j > k) {
                jobs[k] = jobs[j];
            }
            k++;
        }
    }
    if (k == 0)
        jobs.clear();
    else if (k < jobs.size())
        jobs.resize(k);
}

#ifdef _IQTREE_MPI

/**
 * compute and process the best model for partitions (for MPI)
 */
void PartitionFinder::getBestModelforPartitionsMPI(int nthreads, vector<int> &jobs, double* run_time, double* wait_time, double* fstep_time, int* partNum,  double& cpu_time, double& wall_time) {

    if (jobs.empty())
        return;

    bool parallel_job = false;
    
    // reset the arrays
    memset(run_time, 0, sizeof(double)*nthreads);
    memset(wait_time, 0, sizeof(double)*nthreads);
    memset(fstep_time, 0, sizeof(double)*nthreads);
    memset(partNum, 0, sizeof(int)*nthreads);

    // wall time and cpu time
    cpu_time = getCPUTime();
    wall_time = getRealTime();

#ifdef _OPENMP
    parallel_job = (jobs.size() > 1);
#pragma omp parallel for schedule(dynamic) if(parallel_job)
#endif
    for (int i=0; i<jobs.size(); i++) {
        SyncChkPoint syncChkPt(this, i);
        int next_job_id = jobs[i];
        bool need_next_jobID = true;
        double sub_run_time;
        double sub_wait_time;
        while (next_job_id != -1) {
            next_job_id = computeBestModelforOnePartitionMPI(next_job_id, (parallel_job ? 1 : nthreads), need_next_jobID, syncChkPt, sub_run_time, sub_wait_time);
            run_time[i] += sub_run_time;
            wait_time[i] += sub_wait_time;
            partNum[i]++;
        }
    }
    // Master needs to collect all the answers from workers
    if (MPIHelper::getInstance().isMaster()) {
        SyncChkPoint syncChkPt(this, 0);
        double t_start = getRealTime();
        while (jobdone < tot_job_num) {
            syncChkPt.masterSyncOtherChkpts(false);
        }
        double sub_fstep_time = (getRealTime() - t_start);
        for (size_t i = 0; i < jobs.size(); i++) {
            fstep_time[i] += sub_fstep_time;
        }
    }

#ifdef ONESIDE_COMM
    // worker sends the final process_model_info to master
    if (MPIHelper::getInstance().isWorker()) {
        double t_start = getRealTime();
        bool need_nextJobID = false;
        bool forceToSyn = true;
        int job_type = 1; // partition
        if (tree_id_vec.size() > 0) {
            SyncChkPoint syncChkPt(this, 0);
            syncChkPt.sendChkptToMaster(process_model_info, need_nextJobID, job_type, nullptr, forceToSyn);
        }
        double sub_fstep_time = (getRealTime() - t_start);
        for (size_t i = 0; i < jobs.size(); i++) {
            fstep_time[i] += sub_fstep_time;
        }
    }
#endif // ONESIDE_COMM
    cpu_time = getCPUTime() - cpu_time;
    wall_time = getRealTime() - wall_time;
}

/**
 * compute and process the best model for merges (for MPI)
 */
void PartitionFinder::getBestModelforMergesMPI(int nthreads, vector<MergeJob* >& jobs, double* run_time, double* wait_time, double* fstep_time, int* partNum, double& cpu_time, double& wall_time)  {

    if (jobs.empty())
        return;

    bool parallel_job = false;

    // reset the arrays
    memset(run_time, 0, sizeof(double)*nthreads);
    memset(wait_time, 0, sizeof(double)*nthreads);
    memset(fstep_time, 0, sizeof(double)*nthreads);
    memset(partNum, 0, sizeof(int)*nthreads);

    // wall time and cpu time
    cpu_time = getCPUTime();
    wall_time = getRealTime();

#ifdef _OPENMP
    parallel_job = (jobs.size() > 1);
#pragma omp parallel for schedule(dynamic) if(parallel_job)
#endif
    for (int i=0; i<jobs.size(); i++) {
        SyncChkPoint syncChkPt(this, i);
        MergeJob* curr_job = jobs[i];
        bool need_next_jobID = false;
        double sub_run_time;
        double sub_wait_time;
        need_next_jobID = true;
        while (!curr_job->isEmpty()) {
            getBestModelForOneMergeMPI(curr_job, (parallel_job ? 1 : nthreads), need_next_jobID, syncChkPt, sub_run_time, sub_wait_time);
            run_time[i] += sub_run_time;
            wait_time[i] += sub_wait_time;
            partNum[i]++;
        }
    }
    // master needs to wait and collect all the answers from workers
    if (MPIHelper::getInstance().isMaster()) {
        SyncChkPoint syncChkPt(this, 0);
        double t_start = getRealTime();
        while (jobdone < tot_job_num) {
            syncChkPt.masterSyncOtherChkpts(false);
        }
        double sub_fstep_time = (getRealTime() - t_start);
        for (size_t i = 0; i < jobs.size(); i++) {
            fstep_time[i] += sub_fstep_time;
        }
    }
    cpu_time = getCPUTime() - cpu_time;
    wall_time = getRealTime() - wall_time;
}

#endif // _IQTREE_MPI

/**
 * compute and process the best model for partitions (without MPI)
 * nthreads : the number of threads available for these jobs
 */
void PartitionFinder::getBestModelforPartitionsNoMPI(int nthreads, vector<pair<int,double> >& jobs) {

    if (jobs.empty())
        return;

    bool parallel_job = false;

#ifdef _OPENMP
    // parallel_job = ((!params->model_test_and_tree) && nthreads > 1 && jobs.size() > nthreads);
    parallel_job = ((!params->model_test_and_tree) && nthreads > 1);
#pragma omp parallel for schedule(dynamic) reduction(+: lhsum, dfsum) if (parallel_job)
#endif
    for (int j = 0; j < jobs.size(); j++) {
        int tree_id = jobs[j].first;
        PhyloTree *this_tree = in_tree->at(tree_id);
        // scan through models for this partition, assuming the information occurs consecutively
        ModelCheckpoint part_model_info;

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            extractModelInfo(this_tree->aln->name, *model_info, part_model_info);
        }

        // do the computation
        string part_model_name;
        if (params->model_name.empty())
            part_model_name = this_tree->aln->model_name;
        CandidateModel best_model;

        best_model = CandidateModelSet().test(*params, this_tree, part_model_info, models_block,
                                              (parallel_job ? 1 : nthreads), brlen_type, this_tree->aln->name, part_model_name, test_merge);

        bool check = (best_model.restoreCheckpoint(&part_model_info));
        ASSERT(check);

        double score = best_model.computeICScore(this_tree->getAlnNSite());
        this_tree->aln->model_name = best_model.getName();
        lhsum += (lhvec[tree_id] = best_model.logl);
        dfsum += (dfvec[tree_id] = best_model.df);
        lenvec[tree_id] = best_model.tree_len;

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
            replaceModelInfo(this_tree->aln->name, *model_info, part_model_info);
            model_info->dump();
        }
    }
}

/**
 * compute and process the best model for merges (without MPI)
 * nthreads : the number of threads available for these jobs
 */
void PartitionFinder::getBestModelforMergesNoMPI(int nthreads, vector<pair<int,double> >& jobs) {
    if (jobs.empty())
        return;

    bool parallel_job = false;

#ifdef _OPENMP
    // parallel_job = ((!params->model_test_and_tree) && nthreads > 1 && jobs.size() > nthreads);
    parallel_job = ((!params->model_test_and_tree) && nthreads > 1);
#pragma omp parallel for schedule(dynamic) if (parallel_job)
#endif
    for (int j = 0; j < jobs.size(); j++) {
        // information of current partitions pair
        int pair = jobs[j].first;
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
            model_info->startStruct(cur_pair.set_name);
            if (model_info->getBestModel(best_model.subst_name)) {
                best_model.restoreCheckpoint(model_info);
                done_before = true;
            }
            model_info->endStruct();
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

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                extractModelInfo(cur_pair.set_name, *model_info, part_model_info);
                transferModelParameters(in_tree, *model_info, part_model_info, gene_sets[cur_pair.part1], gene_sets[cur_pair.part2]);
            }

            tree->num_precision = in_tree->num_precision;
            tree->setParams(params);
            tree->sse = params->SSE;
            tree->optimize_by_newton = params->optimize_by_newton;
            tree->setNumThreads(params->model_test_and_tree ? num_threads : 1);
            {
                tree->setCheckpoint(&part_model_info);
                // trick to restore checkpoint
                tree->restoreCheckpoint();
                tree->saveCheckpoint();
            }
            best_model = CandidateModelSet().test(*params, tree, part_model_info, models_block,
                                                  parallel_job ? 1 : nthreads, params->partition_type, cur_pair.set_name, "", true);
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
        cur_pair.score = computeInformationScore(lhnew, dfnew, ssize, params->model_test_criterion);
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if (!done_before) {
                replaceModelInfo(cur_pair.set_name, *model_info, part_model_info);
                model_info->dump();
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
}

/**
 * compute the best model
 * job_type = 1 : for all partitions
 * job_type = 2 : for all merges
 */
void PartitionFinder::getBestModel(int job_type) {

    vector<pair<int,double> > jobIDs;
    vector<int> currPartJobs; // for partition jobs
#ifdef _IQTREE_MPI
    vector<MergeJob* > currMergeJobs; // for merge jobs
#endif
    vector<int>* jobs;
    vector<int> closest_p_vector;
    bool run_parallel;
    int i,w;

    if (job_type == 1) {
        // for partitions
        // sort partition by computational cost for OpenMP effciency
        for (i = 0; i < in_tree->size(); i++) {
            Alignment *this_aln = in_tree->at(i)->aln;
            // computation cost is proportional to #sequences, #patterns, and #states
            jobIDs.push_back({i, ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states});
        }
    } else if (num_processes == 1 || MPIHelper::getInstance().isMaster()) {
        // for merges
        ASSERT(gene_sets.size() == lenvec.size());
        better_pairs.clear();
        // find closest partition pairs
        closest_pairs.clear();
        findClosestPairs(super_aln, lenvec, gene_sets, false, closest_pairs);
        if (params->partfinder_log_rate) {
            // additional consider pairs by log-rate
            vector<SubsetPair> log_closest_pairs;
            findClosestPairs(super_aln, lenvec, gene_sets, true, log_closest_pairs);
            mergePairs(closest_pairs, log_closest_pairs);
        }
        
        // sort partition by computational cost for OpenMP/MPI effciency
        for (i = 0; i < closest_pairs.size(); i++) {
            // computation cost is proportional to #sequences, #patterns, and #states
            Alignment *this_aln = in_tree->at(closest_pairs[i].first)->aln;
            closest_pairs[i].distance = -((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
            this_aln = in_tree->at(closest_pairs[i].second)->aln;
            closest_pairs[i].distance -= ((double)this_aln->getNSeq())*this_aln->getNPattern()*this_aln->num_states;
            jobIDs.push_back({i, -closest_pairs[i].distance});
        }
    }

    if (!params->model_test_and_tree && jobIDs.size() > 1) {
        if ((num_threads > 1 && num_processes == 1)|| (num_processes > 1 && MPIHelper::getInstance().isMaster())) {

            // for merging partitions, retreive the answers from checkpoint
            // and remove those jobs from the array jobIDs
            if (job_type == 2)
                retreiveAnsFrChkpt(jobIDs);

            // sort the jobs
            std::sort(jobIDs.begin(), jobIDs.end(), compareJob);
        }
    }
    tot_job_num = jobIDs.size();
    jobdone = 0;

#ifdef _IQTREE_MPI
    if (num_processes == 1 || params->model_test_and_tree) {
        // not using MPI
#endif

        if (job_type == 1) {
            getBestModelforPartitionsNoMPI(num_threads, jobIDs);
        } else {
            getBestModelforMergesNoMPI(num_threads, jobIDs);
        }

#ifdef _IQTREE_MPI
    } else {

        int num_job_array;
        
        // assign the initial jobs to processors
        if (job_type == 1)
            num_job_array = partjobAssignment(jobIDs, currPartJobs);
        else
            num_job_array = mergejobAssignment(jobIDs, currMergeJobs);

        // initialize the value of base
        base = MPIHelper::getInstance().getProcessID() * num_threads;
        
        // initialize the syn time
        last_syn_time = getRealTime();
        
        // initialize the vectors
        tree_id_vec.clear();
        tree_len_vec.clear();
        model_name_vec.clear();
        score_vec.clear();
        set_name_vec.clear();
        tag_vec.clear();
        tot_jobs_done = 0;
        
        // for timing
        double* run_time = new double[num_threads];
        double* wait_time = new double[num_threads];
        double* fstep_time = new double[num_threads];
        // to check the number of partitions each processor handle
        int* num_part = new int[num_threads];

        double cpu_time;
        double wall_time;
        // initialize the checkpoint for the whole processor
        process_model_info.clear();
        
        // compute the best model
        if (job_type == 1) {
            getBestModelforPartitionsMPI(num_threads, currPartJobs, run_time, wait_time, fstep_time, num_part, cpu_time, wall_time);
        } else {
            getBestModelforMergesMPI(num_threads, currMergeJobs, run_time, wait_time, fstep_time, num_part, cpu_time, wall_time);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        // gather all timing information
        int* num_job_arrays = new int[num_processes];
        int* num_part_arrays = new int[num_processes * num_threads];
        double* run_time_arrays = new double[num_processes * num_threads];
        double* wait_time_arrays = new double[num_processes * num_threads];
        double* fstep_time_arrays = new double[num_processes * num_threads];
        double* cpu_time_array = new double[num_processes];
        double* wall_time_array = new double[num_processes];

        MPI_Gather(&wall_time, 1, MPI_DOUBLE, wall_time_array, 1, MPI_DOUBLE, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(&cpu_time, 1, MPI_DOUBLE, cpu_time_array, 1, MPI_DOUBLE, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(&num_job_array, 1, MPI_INT, num_job_arrays, 1, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(num_part, num_threads, MPI_INT, num_part_arrays, num_threads, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(run_time, num_threads, MPI_DOUBLE, run_time_arrays, num_threads, MPI_DOUBLE, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(wait_time, num_threads, MPI_DOUBLE, wait_time_arrays, num_threads, MPI_DOUBLE, PROC_MASTER, MPI_COMM_WORLD);
        MPI_Gather(fstep_time, num_threads, MPI_DOUBLE, fstep_time_arrays, num_threads, MPI_DOUBLE, PROC_MASTER, MPI_COMM_WORLD);

        // show all timing information
        if (MPIHelper::getInstance().isMaster() && tot_job_num > 0) {
            cout << endl;
            cout << "\tproc\tthres\trun time\twait time\tfinal-step time\tnumber-parts" << endl;
            for (int w=0; w<num_processes; w++) {
                for (int t=0; t<num_job_arrays[w]; t++) {
                    cout << "\t" << w << "\t" << t << "\t" << run_time_arrays[w*num_threads+t] << "\t" << wait_time_arrays[w*num_threads+t] << "\t" << fstep_time_arrays[w*num_threads+t] << "\t" << num_part_arrays[w*num_threads+t]<< endl;
                }
            }
            cout << endl;
            cout << "\tproc\tcpu_time\twall_time"<<endl;
            for (int w=0; w<num_processes; w++) {
                cout << "\t" << w << "\t" <<  cpu_time_array[w] << "\t" << wall_time_array[w] << endl;
            }
            cout << endl;
        }
        
        if (job_type == 1) {
            // distribute the checkpoints from Master to Workers
            // for merging, the checkpoints will be distributed to workers after the merging finishes.
            double time_start = getRealTime();
            MPIHelper::getInstance().broadcastCheckpoint(model_info);
            if (MPIHelper::getInstance().isMaster()) {
                cout << "\tTime used for distributing the checkpoints from Master to Workers: " << getRealTime() - time_start << endl;
                cout << endl;
            }
        }
        
        // consolidate the results
        if (job_type == 1) {
            consolidPartitionResults();
        } else if (MPIHelper::getInstance().isMaster()) {
            // for merging, only Master needs consolidation.
            // the workers will do so after the merging finishes.
            consolidMergeResults();
        }
    
        if (job_type == 1) {
            currPartJobs.clear();
        } else if (job_type == 2) {
            // clear the memory for currMergeJobs
            for (i = 0; i < currMergeJobs.size(); i++) {
                delete (currMergeJobs[i]);
            }
            currMergeJobs.clear();
        }
        
        // clear the memory for timing
        delete[] run_time;
        delete[] wait_time;
        delete[] fstep_time;
        delete[] num_part;
        delete[] num_job_arrays;
        delete[] num_part_arrays;
        delete[] run_time_arrays;
        delete[] wait_time_arrays;
        delete[] fstep_time_arrays;
        delete[] cpu_time_array;
        delete[] wall_time_array;


    }
#endif

}

#ifdef _IQTREE_MPI

/*
 * Consolidate the partition results (for MPI)
 */
void PartitionFinder::consolidPartitionResults() {
    int i;

    for (i = 0; i < in_tree->size(); i++) {
        PhyloTree *this_tree = in_tree->at(i);

        string bestModel_key = this_tree->aln->name + CKP_SEP + "best_model_" + criterionName(params->model_test_criterion);
        string bestModel;
        string bestScore_key = this_tree->aln->name + CKP_SEP + "best_score_" + criterionName(params->model_test_criterion);
        double bestScore;

        ASSERT(model_info->getString(bestModel_key, bestModel));
        ASSERT(model_info->get(bestScore_key, bestScore));

        string info_key = this_tree->aln->name + CKP_SEP + bestModel;
        string info;
        double logL;
        int df;
        double treeLen;

        ASSERT(model_info->getString(info_key, info));
        size_t pos1 = info.find_first_of(" ");
        ASSERT (pos1 != string::npos && pos1 > 0);
        size_t pos2 = info.find_first_of(" ", pos1+1);
        ASSERT (pos2 != string::npos && pos2 > pos1+1);
        logL = atof(info.substr(0,pos1).c_str());
        df = atoi(info.substr(pos1+1,pos2-pos1-1).c_str());
        treeLen = atof(info.substr(pos2+1).c_str());

        this_tree->aln->model_name = bestModel;
        lhsum += (lhvec[i] = logL);
        dfsum += (dfvec[i] = df);
        lenvec[i] = treeLen;
    }
}

/*
 * Consolidate the merge results (for MPI)
 */
void PartitionFinder::consolidMergeResults() {
    
    better_pairs.clear();
    for (size_t pair = 0; pair < closest_pairs.size(); pair++) {
        // information of current partitions pair
        ModelPair cur_pair;
        cur_pair.part1 = closest_pairs[pair].first;
        cur_pair.part2 = closest_pairs[pair].second;
        ASSERT(cur_pair.part1 < cur_pair.part2);
        cur_pair.merged_set.clear();
        cur_pair.merged_set.insert(gene_sets[cur_pair.part1].begin(), gene_sets[cur_pair.part1].end());
        cur_pair.merged_set.insert(gene_sets[cur_pair.part2].begin(), gene_sets[cur_pair.part2].end());
        cur_pair.set_name = getSubsetName(in_tree, cur_pair.merged_set);
        double weight1 = getSubsetAlnLength(in_tree, gene_sets[cur_pair.part1]);
        double weight2 = getSubsetAlnLength(in_tree, gene_sets[cur_pair.part2]);
        double sum = 1.0 / (weight1 + weight2);
        weight1 *= sum;
        weight2 *= sum;
        CandidateModel best_model;
        
        model_info->startStruct(cur_pair.set_name);
        ASSERT(model_info->getBestModel(best_model.subst_name));
        best_model.restoreCheckpoint(model_info);
        model_info->endStruct();
        
        cur_pair.logl = best_model.logl;
        cur_pair.df = best_model.df;
        cur_pair.model_name = best_model.getName();
        cur_pair.tree_len = best_model.tree_len;
        
        double lhnew = lhsum - lhvec[cur_pair.part1] - lhvec[cur_pair.part2] + best_model.logl;
        int dfnew = dfsum - dfvec[cur_pair.part1] - dfvec[cur_pair.part2] + best_model.df;
        cur_pair.score = computeInformationScore(lhnew, dfnew, ssize, params->model_test_criterion);
        
        if (cur_pair.score < inf_score) {
            better_pairs.insertPair(cur_pair);
        }
    }
}

#endif

void PartitionFinder::test_PartitionModel() {

    int i, job_type;

    lhsum = 0.0;
    dfsum = 0;
    if (params->partition_type == BRLEN_FIX || params->partition_type == BRLEN_SCALE) {
        dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
        if (params->partition_type == BRLEN_SCALE)
            dfsum -= 1;
    }
    ssize = in_tree->getAlnNSite();
    num_model = 0;
    total_num_model = in_tree->size();

#ifdef _IQTREE_MPI
    // initialize the shared memory space
    initialMPIShareMemory();
#endif

    // get the name of the algorithm
    string part_algo = "";
    if (params->partition_merge == MERGE_GREEDY)
        part_algo = "Greedy Algorithm";
    else if (params->partition_merge == MERGE_RCLUSTER)
        part_algo = "Relaxed Algorithm";
    else if (params->partition_merge == MERGE_RCLUSTERF)
        part_algo = "Fast Relaxed Algorithm";
    else if (params->partition_merge == MERGE_KMEANS)
        part_algo = "Kmean Algorithm";

    // for greedy algorithm
    if (params->partition_merge == MERGE_GREEDY) {
        params->partfinder_rcluster_max = in_tree->size() * (in_tree->size()-1) / 2;
        params->partfinder_log_rate = false;
        params->partfinder_rcluster = 100.0;
    }
    
    // 2017-06-07: -rcluster-max for max absolute number of pairs
    if (params->partfinder_rcluster_max == 0) {
        // params->partfinder_rcluster_max = max((size_t)1000, 10 * in_tree->size());
        params->partfinder_rcluster_max = 10 * in_tree->size();
    }

    // show the parameters for partition finder
    cout << endl;
    cout << "PartitionFinder's parameters:" << endl;
    cout << part_algo << endl;
    cout << "Percentage: " << params->partfinder_rcluster << endl;
    cout << "Maximum pairs: " << params->partfinder_rcluster_max << endl;
    cout << endl;

    if (params->partition_merge != MERGE_NONE) {
        double p = params->partfinder_rcluster/100.0;
        size_t num_pairs = round(in_tree->size()*(in_tree->size()-1)*p/2);
        if (p < 1.0)
            num_pairs = min(num_pairs, params->partfinder_rcluster_max);
        total_num_model += num_pairs;
        for (i = in_tree->size()-2; i > 0; i--)
            total_num_model += max(round(i*p), 1.0);
    }

    if (num_threads <= 0) {
        // partition selection scales well with many cores
        num_threads = min((int64_t)countPhysicalCPUCores(), total_num_model);
        num_threads = min(num_threads, params->num_threads_max);
        omp_set_num_threads(num_threads);
        cout << "NUMBER OF THREADS FOR PARTITION FINDING: " << num_threads << endl;
    }

    start_time = getRealTime();

    super_aln = ((SuperAlignment*)in_tree->aln);

    cout << "Selecting individual models for " << in_tree->size() << " charsets using " << criterionName(params->model_test_criterion) << "..." << endl;
    //cout << " No. AIC         AICc        BIC         Charset" << endl;
    cout << " No. Model        Score       TreeLen     Charset" << endl;

    lhvec.resize(in_tree->size());
    dfvec.resize(in_tree->size());
    lenvec.resize(in_tree->size());

    test_merge = (params->partition_merge != MERGE_NONE) && params->partition_type != TOPO_UNLINKED && (in_tree->size() > 1);

    brlen_type = params->partition_type;
    if (brlen_type == TOPO_UNLINKED) {
        brlen_type = BRLEN_OPTIMIZE;
    }

    // compute the best model for all partitions
    job_type = 1; // for all partitions
    getBestModel(job_type);

    // in case ModelOMatic change the alignment
    fixPartitions(in_tree);

    inf_score = computeInformationScore(lhsum, dfsum, ssize, params->model_test_criterion);
    cout << "Full partition model " << criterionName(params->model_test_criterion)
         << " score: " << inf_score << " (LnL: " << lhsum << "  df:" << dfsum << ")" << endl;

    if (!test_merge) {
        super_aln->printBestPartition((string(params->out_prefix) + ".best_scheme.nex").c_str());
        super_aln->printBestPartitionRaxml((string(params->out_prefix) + ".best_scheme").c_str());
        model_info->dump();
        return;
    }

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

    if (params->partition_merge == MERGE_KMEANS) {
        // kmeans cluster based on parition tree length
        double cur_score = inf_score;
        for (int ncluster = in_tree->size()-1; ncluster >= 1; ncluster--) {
            vector<set<int> > this_gene_sets;
            StrVector this_model_names;
            //double sum = in_tree->size()/std::accumulate(lenvec.begin(), lenvec.end(), 0.0);
            double score = doKmeansClustering(*params, in_tree, ncluster, lenvec, *model_info,
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
    bool perform_merge = (params->partition_merge != MERGE_KMEANS && gene_sets.size() >= 2);
    while (params->partition_merge != MERGE_KMEANS && gene_sets.size() >= 2) {
        // stepwise merging charsets

        // get the closest partition pairs, and
        // compute the best model for each pair
        job_type = 2; // for all merges
        getBestModel(job_type);
        
        bool is_pairs_empty = better_pairs.empty();
        
#ifdef _IQTREE_MPI
        MPI_Bcast(&is_pairs_empty, 1, MPI_CXX_BOOL,PROC_MASTER, MPI_COMM_WORLD);
#endif

        if (is_pairs_empty) break;
        
        if (MPIHelper::getInstance().isMaster()) {

            ModelPairSet compatible_pairs;
            
            int num_comp_pairs = params->partition_merge == MERGE_RCLUSTERF ? gene_sets.size()/2 : 1;
            better_pairs.getCompatiblePairs(num_comp_pairs, compatible_pairs);
            if (compatible_pairs.size() > 1)
                cout << compatible_pairs.size() << " compatible better partition pairs found" << endl;
            
            // 2017-12-21: simultaneously merging better pairs
            for (auto it_pair = compatible_pairs.begin(); it_pair != compatible_pairs.end(); it_pair++) {
                ModelPair opt_pair = it_pair->second;
                
                lhsum = lhsum - lhvec[opt_pair.part1] - lhvec[opt_pair.part2] + opt_pair.logl;
                dfsum = dfsum - dfvec[opt_pair.part1] - dfvec[opt_pair.part2] + opt_pair.df;
                inf_score = computeInformationScore(lhsum, dfsum, ssize, params->model_test_criterion);
                ASSERT(inf_score <= opt_pair.score + 0.1);
                
                cout << "Merging " << opt_pair.set_name << " with " << criterionName(params->model_test_criterion)
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
    }

#ifdef _IQTREE_MPI
    if (perform_merge) {
        // distribute the checkpoints from Master to Workers
        double time_start = getRealTime();
        MPIHelper::getInstance().broadcastCheckpoint(model_info);
        if (MPIHelper::getInstance().isMaster()) {
            cout << "\tTime used for distributing the checkpoints from Master to Workers: " << getRealTime() - time_start << endl;
            cout << endl;
        }
    }
#endif
    
    if (MPIHelper::getInstance().isMaster()) {
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
    }

#ifdef _IQTREE_MPI
    if (perform_merge && num_processes > 1) {
        SyncChkPoint syncChkPoint(this,0);
        // broadcast gene_sets from Master to Workers
        syncChkPoint.broadcastVecSetInt(gene_sets);
        // broadcast model_names from Master to Workers
        syncChkPoint.broadcastVecStr(model_names);
    }
#endif

    if (gene_sets.size() < in_tree->size())
        mergePartitions(in_tree, gene_sets, model_names);
    
    if (!iEquals(params->merge_models, "all")) {
        // test all candidate models again
        lhsum = 0.0;
        dfsum = 0;
        if (params->partition_type == BRLEN_FIX || params->partition_type == BRLEN_SCALE) {
            dfsum = in_tree->getNBranchParameters(BRLEN_OPTIMIZE);
            if (params->partition_type == BRLEN_SCALE)
                dfsum -= 1;
        }
        // compute the best model for all partitions
        // but this time "test_merge = false"
        test_merge = false;
        job_type = 1; // for all partitions
        getBestModel(job_type);
    }

    inf_score = computeInformationScore(lhsum, dfsum, ssize, params->model_test_criterion);
    cout << "Best partition model " << criterionName(params->model_test_criterion) << " score: " << inf_score << " (LnL: " << lhsum << "  df:" << dfsum << ")" << endl;

    ((SuperAlignment*)in_tree->aln)->printBestPartition((string(params->out_prefix) + ".best_scheme.nex").c_str());
    ((SuperAlignment*)in_tree->aln)->printBestPartitionRaxml((string(params->out_prefix) + ".best_scheme").c_str());
    model_info->dump();

#ifdef _IQTREE_MPI
    // free the MPI share memory
    freeMPIShareMemory();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    cout << "Finish the procedure test_PartitionModel()" << endl;
    
#endif
}

// -----------------------------------------------------------------------------
// The following functions and the classes SyncChkPoint and MergeJob are for MPI
// -----------------------------------------------------------------------------

#ifdef _IQTREE_MPI


/*
 *  initialize the shared memory space to be accessed by the other processors
 */
void PartitionFinder::initialMPIShareMemory() {
#ifdef ONESIDE_COMM
    if (MPIHelper::getInstance().getProcessID()==PROC_MASTER) {
        val_ptr = (int*) malloc(sizeof(int));
        MPI_Win_create(val_ptr, sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    } else {
        val_ptr = NULL;
        MPI_Win_create(val_ptr, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    }
#else
    win = NULL;
    val_ptr = NULL;
#endif // ONESIDE_COMM
}

/*
 *  free the shared memory space
 */
void PartitionFinder::freeMPIShareMemory() {
#ifdef ONESIDE_COMM
    MPI_Win_free(&win);
#endif
    if (val_ptr != nullptr) {
        delete[] val_ptr;
        val_ptr = nullptr;
    }
}

/*
 * assign initial partition jobs to processors
 * input: a set of partition jobs ordered by the estimated computational costs
 * output: the set of jobs in currJobs, number of jobs assigned <= number of threads
 */
int PartitionFinder::partjobAssignment(vector<pair<int,double> > &job_ids, vector<int> &currJobs) {

    int num_job_assigned = 0;
    int n = num_processes * num_threads;
    int* assignJobs = new int[n];
    int remain_job_num;
    nextjob = 0;
    remain_job_list.clear();
    if (MPIHelper::getInstance().isMaster()) {
        int k = 0;
        while (k < n && k < job_ids.size()) {
            assignJobs[k] = job_ids[k].first;
            k++;
        }
        if (k < n) {
            for (int i = k; i < n; i++)
                assignJobs[i] = -1;
        }
        while (k < job_ids.size()) {
            remain_job_list.push_back(job_ids[k].first);
            k++;
        }
        remain_job_num = remain_job_list.size();
    }
    currJobs.resize(num_threads);
    MPI_Scatter(assignJobs, num_threads, MPI_INT, &currJobs[0], num_threads, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&remain_job_num, 1, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
    if (MPIHelper::getInstance().isWorker())
        remain_job_list.resize(remain_job_num);
    MPI_Bcast(&remain_job_list[0], remain_job_num, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
    // resize the currJobs if there are some empty jobs
    int i = 0;
    while (i < currJobs.size() && currJobs[i] != -1)
        i++;
    if (i == 0)
        currJobs.clear();
    else if (i < currJobs.size())
        currJobs.resize(i);
    // clear the memory
    delete[] assignJobs;
    return currJobs.size();
}

/*
 * assign initial merge jobs to processors
 * input: a set of merge jobs ordered by the estimated computational costs
 * output: number of items in currJobs
 */
/*
int PartitionFinder::mergejobAssignment(vector<pair<int,double> > &job_ids, vector<MergeJob*>&currJobs) {
    int i,j,k;
    int w,id1,id2,accum_len,len;

    nextjob = 0;
    // clear any existing jobs
    for (k=0; k<currJobs.size(); k++) {
        delete(currJobs[k]);
    }
    currJobs.clear();
    // clear any remaining jobs
    for (k=0; k<remain_mergejobs.size(); k++) {
        delete(remain_mergejobs[k]);
    }
    remain_mergejobs.clear();
    
    int n = num_processes * num_threads;
    int* scounts = new int[num_processes];
    int* displs = new int[num_processes];
    int* alljoblens = NULL;
    int* joblens = new int[num_threads];
    int pid;
    char* sendbuf = NULL;
    char* recvbuf = NULL;
    int recvlen;
    if (MPIHelper::getInstance().isMaster()) {
        // assign one job to every thread
        i = 0;
        accum_len = 0;
        string data_str = "";
        alljoblens = new int[n];
        for (j=0; j<num_processes; j++) {
            len = 0;
            for (k=0; k<num_threads; k++) {
                string str;
                if (i < job_ids.size()) {
                    w = job_ids[i].first;
                    id1 = closest_pairs[w].first;
                    id2 = closest_pairs[w].second;
                    MergeJob mjob = MergeJob(id1,id2,gene_sets[id1],gene_sets[id2],lenvec[id1],lenvec[id2]);
                    mjob.toString(str);
                } else {
                    // empty job
                    MergeJob mjob;
                    mjob.toString(str);
                }
                len += str.length();
                alljoblens[i] = str.length();
                data_str.append(str);
                i++;
            }
            displs[j] = accum_len;
            scounts[j] = len;
            accum_len += len;
        }
        sendbuf = new char[accum_len];
        strcpy(sendbuf, data_str.c_str());

        // place all the unassigned jobs to the array remain_mergejobs
        while (i<job_ids.size()) {
            w = job_ids[i].first;
            id1 = closest_pairs[w].first;
            id2 = closest_pairs[w].second;
            remain_mergejobs.push_back(new MergeJob(id1,id2,gene_sets[id1],gene_sets[id2],lenvec[id1],lenvec[id2]));
            i++;
        }
    }
    MPI_Scatter(alljoblens, num_threads, MPI_INT, joblens, num_threads, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
    recvlen = 0;
    for (k = 0; k < num_threads; k++)
        recvlen += joblens[k];
    recvbuf = new char[recvlen];
    MPI_Scatterv(sendbuf, scounts, displs, MPI_CHAR, recvbuf, recvlen, MPI_CHAR, PROC_MASTER, MPI_COMM_WORLD);

    // add the merge jobs to currJobs
    k = 0;
    for (j = 0; j < num_threads; j++) {
        string str = string(recvbuf, k, joblens[j]);
        MergeJob* mjob = new MergeJob();
        mjob->loadFrString(str);
        currJobs.push_back(mjob);
        k += joblens[j];
    }
    
    // release the memory
    delete[] scounts;
    delete[] displs;
    delete[] joblens;
    delete[] recvbuf;
    if (sendbuf != NULL)
        delete[] sendbuf;
    if (alljoblens != NULL)
        delete[] alljoblens;
    return currJobs.size();
}
*/

/*
 * assign initial merge jobs to processors
 * input: a set of merge jobs ordered by the estimated computational costs
 * output: number of items in currJobs
 */
int PartitionFinder::mergejobAssignment(vector<pair<int,double> > &job_ids, vector<MergeJob*>&currJobs) {

    int i,j,k;
    int w,id1,id2;

    nextjob = 0;
    // clear any existing jobs
    for (k=0; k<currJobs.size(); k++) {
        delete(currJobs[k]);
    }
    currJobs.clear();
    // clear any remaining jobs
    for (k=0; k<remain_mergejobs.size(); k++) {
        delete(remain_mergejobs[k]);
    }
    remain_mergejobs.clear();

    if (MPIHelper::getInstance().isMaster()) {
        // MASTER: assign one job to its every thread
        i=0;
        while (i<num_threads && i<job_ids.size()) {
            w = job_ids[i].first;
            id1 = closest_pairs[w].first;
            id2 = closest_pairs[w].second;
            currJobs.push_back(new MergeJob(id1,id2,gene_sets[id1],gene_sets[id2],lenvec[id1],lenvec[id2]));
            i++;
        }
        // MASTER: send one job to every thread of the processors
        SyncChkPoint syncChkPoint(this, 0);
        for (j=1; j<num_processes; j++) {
            for (k=0; k<num_threads; k++) {
                int n = 0;
                int tag = j * num_threads + k;
                if (i < job_ids.size()) {
                    w = job_ids[i].first;
                    id1 = closest_pairs[w].first;
                    id2 = closest_pairs[w].second;
                    MergeJob mergejob(id1,id2,gene_sets[id1],gene_sets[id2],lenvec[id1],lenvec[id2]);
                    syncChkPoint.sendMergeJobToWorker(mergejob,j,tag);
                    i++;
                } else {
                    MergeJob mergejob;
                    syncChkPoint.sendMergeJobToWorker(mergejob,j,tag); // send the empty job to the worker
                }
            }
        }
        // place all the unassigned jobs to the array remain_mergejobs
        while (i<job_ids.size()) {
            w = job_ids[i].first;
            id1 = closest_pairs[w].first;
            id2 = closest_pairs[w].second;
            remain_mergejobs.push_back(new MergeJob(id1,id2,gene_sets[id1],gene_sets[id2],lenvec[id1],lenvec[id2]));
            i++;
        }
    } else {

        // WORKER: receive jobs from the master
        int n;
        int src, tag;
        SyncChkPoint syncChkPoint(this, 0);
        int proc_id = MPIHelper::getInstance().getProcessID();
        for (i=0; i<num_threads; i++) {
            MergeJob* mergeJob = new MergeJob();
            tag = proc_id * num_threads + i;
            syncChkPoint.recMergeJobFrMaster(*mergeJob, tag);
            if (mergeJob->id1 == -1) {
                // empty job
                delete(mergeJob);
            } else {
                currJobs.push_back(mergeJob);
            }
        }
    }
    return currJobs.size();
}

/*  constructor
 */
SyncChkPoint::SyncChkPoint(PartitionFinder* pf, int thres_id) {

    pfinder = pf;
    mytag = thres_id + pf->base;
}

/*
 * FOR MASTER
 * Show the other worker's result of best model
 */
void SyncChkPoint::showResult(ModelCheckpoint& part_model_info, int work_tag) {
    string key, data_num;
    int job_type;

    key = "pf_data_num";
    ASSERT(part_model_info.get(key, data_num));
    key = "pf_job_type";
    ASSERT(part_model_info.get(key, job_type));

    if (data_num == "single") {
        double tree_len,score;
        string model_name, set_name;
        int tree_id;
        bool done_before;

        key = "pf_tree_len";
        ASSERT(part_model_info.get(key, tree_len));
        key = "pf_model_name";
        ASSERT(part_model_info.get(key, model_name));

        if (job_type == 1) {
            // partition
            key = "pf_tree_id";
            ASSERT(part_model_info.get(key, tree_id));
            key = "pf_score";
            ASSERT(part_model_info.get(key, score));
            pfinder->showPartitionResult(part_model_info, tree_id, tree_len, model_name, score, work_tag);
        } else {
            // merge
            key = "pf_set_name";
            ASSERT(part_model_info.get(key, set_name));
            key = "pf_done_before";
            ASSERT(part_model_info.getBool(key, done_before));
            pfinder->showMergeResult(part_model_info, tree_len, model_name, set_name, done_before, work_tag);
        }

    } else {
        vector<int> tree_id_vec;
        vector<double> tree_len_vec;
        vector<string> model_name_vec;
        vector<double> score_vec;
        vector<int> tag_vec;
        vector<string> set_name_vec;
        int tot_jobsdone;
        key = "pf_tree_len";
        ASSERT(part_model_info.getVector(key, tree_len_vec));
        key = "pf_model_name";
        ASSERT(part_model_info.getVector(key, model_name_vec));
        key = "pf_tag";
        ASSERT(part_model_info.getVector(key, tag_vec));

        if (job_type == 1) {
            // partition
            key = "pf_tree_id";
            ASSERT(part_model_info.getVector(key, tree_id_vec));
            key = "pf_score";
            ASSERT(part_model_info.getVector(key, score_vec));
            pfinder->showPartitionResults(part_model_info, tree_id_vec, tree_len_vec, model_name_vec, score_vec, tag_vec);
        } else {
            // merge
            key = "pf_set_name";
            ASSERT(part_model_info.getVector(key, set_name_vec));
            key = "pf_tot_jobs_done";
            ASSERT(part_model_info.get(key, tot_jobsdone));
            pfinder->showMergeResults(part_model_info, tree_len_vec, model_name_vec, set_name_vec, tag_vec, tot_jobsdone);
        }
    }
}

/*
 * FOR MASTER - synchronize the checkpoints from the other processors
 * Receive checkpoint from worker and send the next Job ID to workers
 * increase the value of next_job and job_done by 1
 * update the master's checkpoint: model_info
 */
void SyncChkPoint::masterSyncOtherChkpts(bool chk_gotMessage) {

    if (MPIHelper::getInstance().isWorker() || mytag > 0)
        return;

    ModelCheckpoint proc_model_info;
    string key;
    int worker, next_jobID, job_type, work_tag, tree_id;
    bool job_finished, need_nextJobID, proceed, thread_finished, is_old_result;
    map<int,int>::iterator itr;

    next_jobID = -1;
    job_finished = false;
    is_old_result = false;

    if (chk_gotMessage) {
        // only proceed if there is a message
        while (gotMessage(work_tag, worker)) {

            // receive checkpoint from the WORKER
            recvAnyCheckpoint(&proc_model_info, worker, work_tag);
            
            key = "pf_job_type";
            ASSERT(proc_model_info.get(key, job_type));
            if (job_type == 1) {
                // for partition job
#ifdef SYN_COMM
                key = "need_nextJobID";
                need_nextJobID = proc_model_info.getBool(key);
                if (need_nextJobID) {
                    // get the next Job ID
                    next_jobID = getNextJobID();
                    // send the next job ID to the WORKER
                    MPI_Send(&next_jobID, 1, MPI_INT, worker, work_tag, MPI_COMM_WORLD);
                }
#endif  // SYN_COMM
            } else {
                // for merge job
                key = "need_nextJobID";
                need_nextJobID = proc_model_info.getBool(key);
                if (need_nextJobID) {
                    // send the next mergejob to the WORKER
                    MergeJob mergeJob;
                    getNextMergeJob(&mergeJob);
                    sendMergeJobToWorker(mergeJob, worker, work_tag);
                }
            }

            showResult(proc_model_info, work_tag);

            proc_model_info.clear();

        }
    } else {

        // receive checkpoint from any worker with any tag
        recvAnyCheckpoint(&proc_model_info, worker, work_tag);

        key = "pf_job_type";
        ASSERT(proc_model_info.get(key, job_type));
        if (job_type == 1) {
            // for partition job
#ifdef SYN_COMM
            key = "need_nextJobID";
            need_nextJobID = proc_model_info.getBool(key);
            if (need_nextJobID) {
                // get the next Job ID
                next_jobID = getNextJobID();
                // send the next job ID to the WORKER
                MPI_Send(&next_jobID, 1, MPI_INT, worker, work_tag, MPI_COMM_WORLD);
            }
#endif // SYN_COMM
        } else {
            // for merge job
            key = "need_nextJobID";
            need_nextJobID = proc_model_info.getBool(key);
            if (need_nextJobID) {
                // send the next mergejob to the WORKER
                MergeJob mergeJob;
                getNextMergeJob(&mergeJob);
                sendMergeJobToWorker(mergeJob, worker, work_tag);
            }
        }

        showResult(proc_model_info, work_tag);
        
        proc_model_info.clear();
    }
}

/*
 * FOR WORKER
 * send checkpoint to master
 *
 * return the next Job ID from master if necessary
 */
int SyncChkPoint::sendChkptToMaster(ModelCheckpoint &model_info, bool need_nextJobID, int job_type, MergeJob* mergeJob, bool forceToSyn) {

    if (MPIHelper::getInstance().getNumProcesses() == 1 || MPIHelper::getInstance().isMaster()) {
        return -1;
    }

    int next_jobID = -1;
    string key;
    bool syn_comm = true;

#ifdef ONESIDE_COMM
    if (job_type == 1) {
        syn_comm = false;
    }
#endif // ONESIDE_COMM

    if (syn_comm) {
        // workers: send checkpoint to MASTER synchronously
        key = "need_nextJobID"; model_info.putBool(key, need_nextJobID);
        key = "pf_job_type"; model_info.put(key, job_type);
        key = "pf_data_num"; model_info.put(key, "single");
        
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            sendCheckpoint(&model_info, PROC_MASTER, mytag);
            if (need_nextJobID) {
                if (job_type == 1) {
                    // receive the next job ID from MASTER synchronously
                    MPI_Recv(&next_jobID, 1, MPI_INT, PROC_MASTER, mytag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    ASSERT(mergeJob != nullptr);
                    // receive the next merge job from MASTER synchronously
                    recMergeJobFrMaster(*mergeJob, mytag);
                }
            }
        }
        model_info.clear();
        
    } else {
        // using ONESIDE communication
        // and send the checkpoint to master when the time is long enough

        if (need_nextJobID) {
            next_jobID = getNextJobID();
        }
        
        string str = "";
        
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            // send the checkpoint to master when the time is long enough
            if (forceToSyn || (getRealTime() - pfinder->last_syn_time > TIME_SYN && (pfinder->tot_jobs_done > 0 || pfinder->tree_id_vec.size() > 0))) {
                // prepare to do synchronization
                pfinder->last_syn_time = getRealTime();
                key = "pf_tree_id"; model_info.putVector(key, pfinder->tree_id_vec);
                key = "pf_tree_len"; model_info.putVector(key, pfinder->tree_len_vec);
                key = "pf_model_name"; model_info.putVector(key, pfinder->model_name_vec);
                key = "pf_score"; model_info.putVector(key, pfinder->score_vec);
                key = "pf_tag"; model_info.putVector(key, pfinder->tag_vec);
                key = "pf_job_type"; model_info.put(key, job_type);
                key = "pf_data_num"; model_info.put(key, "multiple");
                key = "pf_tot_jobs_done"; model_info.put(key, pfinder->tot_jobs_done);
                key = "pf_set_name"; model_info.putVector(key, pfinder->set_name_vec);
                
                stringstream ss;
                model_info.dump(ss);
                str = ss.str();
                
                // do synchronization with Master
                MPIHelper::getInstance().sendString(str, PROC_MASTER, mytag);
                
                // clear all vectors
                pfinder->tree_id_vec.clear();
                pfinder->tree_len_vec.clear();
                pfinder->model_name_vec.clear();
                pfinder->score_vec.clear();
                pfinder->tag_vec.clear();
                pfinder->tot_jobs_done = 0;
                
                // clear the checkpoint for this process
                model_info.clear();
            }
        }
    }

    return next_jobID;
}

/*
 * receive an integer from the master (for synchronous communication)
 */
/*
int SyncChkPoint::recvInt(int tag) {
    int mesg = -1;
#ifdef _IQTREE_MPI
    MPI_Recv(&mesg, 1, MPI_INT, PROC_MASTER, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
    return mesg;
}
*/

/*
 * get the next Job ID
 */
int SyncChkPoint::getNextJobID() {
    int one = 1, indx = -1, nxtJobID = -1;

#ifdef SYN_COMM
    if (MPIHelper::getInstance().isMaster()) {
            // get the next Job ID
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                if (pfinder->nextjob < pfinder->remain_job_list.size()) {
                    nxtJobID = pfinder->remain_job_list[pfinder->nextjob];
                    pfinder->nextjob++;
                }
            }
        }
#endif // SYN_COMM

#ifdef ONESIDE_COMM

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, PROC_MASTER, 0, pfinder->win);
        MPI_Fetch_and_op(&one, &indx, MPI_INT, PROC_MASTER, 0, MPI_SUM, pfinder->win);
        MPI_Win_unlock(PROC_MASTER, pfinder->win);
        if (indx >= 0 && indx < pfinder->remain_job_list.size()) {
            nxtJobID = pfinder->remain_job_list[indx];
        }
    }

#endif // ONESIDE_COMM

    return nxtJobID;
}

/*
 * get the next MergeJob
 */
void SyncChkPoint::getNextMergeJob(MergeJob* mergejob) {

    if (MPIHelper::getInstance().isMaster()) {
        // get the next Job ID
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            if (pfinder->nextjob < pfinder->remain_mergejobs.size()) {
                mergejob->copyFrom(pfinder->remain_mergejobs[pfinder->nextjob]);
                pfinder->nextjob++;
            } else {
                mergejob->setEmpty();
            }
        }
    }
}


void SyncChkPoint::sendCheckpoint(Checkpoint *ckp, int dest, int tag) {
    stringstream ss;
    ckp->dump(ss);
    string str = ss.str();
    MPIHelper::getInstance().sendString(str, dest, tag);
}

void SyncChkPoint::recvAnyCheckpoint(Checkpoint *ckp, int& src, int& tag) {
    string str;
    recvAnyString(str, src, tag);
    stringstream ss(str);
    ckp->load(ss);
}

void SyncChkPoint::recvAnyString(string &str, int& src, int& tag) {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int msgCount;
    MPI_Get_count(&status, MPI_CHAR, &msgCount);
    // receive the message
    char *recvBuffer = new char[msgCount];
    MPI_Recv(recvBuffer, msgCount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    str = recvBuffer;
    src = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    delete [] recvBuffer;
}

/*
 * Check for incoming messages
 * if there is a message, collect the tag value and the source
 */
bool SyncChkPoint::gotMessage(int& tag, int& source) {
    if (MPIHelper::getInstance().getNumProcesses() == 1 || mytag > 0)
        return false;
    int flag = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (flag) {
        tag = status.MPI_TAG;
        source = status.MPI_SOURCE;
        return true;
    } else
        return false;
}

void SyncChkPoint::sendMergeJobToWorker(MergeJob& mergeJob, int dest, int tag) {
    string str;
    mergeJob.toString(str);
    MPIHelper::getInstance().sendString(str, dest, tag);
}

void SyncChkPoint::recMergeJobFrMaster(MergeJob& mergeJob, int tag) {
    string str;
    int src = 0;
    MPIHelper::getInstance().recvString(str, src, tag);
    mergeJob.loadFrString(str);
}

int* SyncChkPoint::toIntArr(vector<set<int> >& gene_sets, int& buffsize) {
    set<int>::iterator itr;
    int i,k;

    buffsize = gene_sets.size();
    for (int i=0; i<gene_sets.size(); i++)
        buffsize += gene_sets[i].size();

    int* buff = new int[buffsize];
    
    k = 0;
    for (i=0; i<gene_sets.size() && k<buffsize; i++) {
        for (itr=gene_sets[i].begin(); itr!=gene_sets[i].end() && k<buffsize; itr++) {
            buff[k] = *itr;
            k++;
        }
        if (k < buffsize) {
            buff[k] = -1;
            k++;
        }
    }
    return buff;
}

void SyncChkPoint::loadFrIntArr(vector<set<int> >& gene_sets, int* buff, int buffsize) {
    int j = 0;
    int k = 0; // size of gene_sets
    bool set_end = true;
    gene_sets.clear();
    while (j < buffsize) {
        if (buff[j] == -1) {
            // end of the current set
            set_end = true;
        } else {
            if (set_end) {
                // create a new set
                set_end = false;
                k++;
                gene_sets.resize(k);
                gene_sets[k-1].clear();
            }
            gene_sets[k-1].insert(buff[j]);
        }
        j++;
    }
}

char* SyncChkPoint::toCharArr(vector<string>& model_names, int& buffsize) {
    string buff_str = "";
    char* buff = NULL;
    int i;
    for (i = 0; i < model_names.size(); i++) {
        buff_str.append(model_names[i]);
        buff_str.append(" ");
    }
    buffsize = buff_str.length() + 1;
    if (buffsize > 0) {
        buff = new char[buffsize];
        strcpy(buff,buff_str.c_str());
        buff[buff_str.length()] = '\0';
    }
    return buff;
}

void SyncChkPoint::loadFrCharArr(vector<string>& model_names, char* buff) {
    model_names.clear();
    if (buff == NULL)
        return;
    string buff_str = string(buff);
    int start_pos = 0;
    for (int j=0; j<buff_str.length(); j++) {
        if (buff_str[j] == ' ') {
            if (start_pos < j)
                model_names.push_back(buff_str.substr(start_pos, j - start_pos));
            else
                model_names.push_back("");
            start_pos = j+1;
        }
    }
}

void SyncChkPoint::broadcastVecSetInt(vector<set<int> >& gene_sets) {
    // broadcast vector<set<int> > object to all workers
    set<int>::iterator itr;
    int buffsize;
    int* buff = NULL;

    // broadcast the buffsize to workers
    if (MPIHelper::getInstance().isMaster())
        buff = toIntArr(gene_sets, buffsize);

    MPI_Bcast(&buffsize, 1, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
    
    if (buffsize > 0) {
        if (MPIHelper::getInstance().isWorker())
            buff = new int[buffsize];
        
        // broadcast buff to workers
        MPI_Bcast(buff, buffsize, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);
        
        // for workers, rebuild the gene_sets
        if (MPIHelper::getInstance().isWorker())
            loadFrIntArr(gene_sets, buff, buffsize);
    }
    
    if (buff != NULL)
        delete[] buff;
}

void SyncChkPoint::broadcastVecStr(vector<string>& model_names) {
    int buffsize;
    char* buff = NULL;

    // for Master, build the long string
    if (MPIHelper::getInstance().isMaster()) {
        buff = toCharArr(model_names, buffsize);
    }
    
    // broadcast buffsize to workers
    MPI_Bcast(&buffsize, 1, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);

    if (buffsize > 0) {
        if (MPIHelper::getInstance().isWorker()) {
            buff = new char[buffsize];
        }
        // broadcast buff to workers
        MPI_Bcast(buff, buffsize, MPI_CHAR, PROC_MASTER, MPI_COMM_WORLD);
        // for workers, rebuid the model_names
        if (MPIHelper::getInstance().isWorker()) {
            loadFrCharArr(model_names, buff);
        }
    }
    
    if (buff != NULL)
        delete[] buff;
}

MergeJob::MergeJob() {
    setEmpty();
}

MergeJob::MergeJob(int id_1, int id_2, set<int>& geneset_1, set<int>& geneset_2, double treelen_1, double treelen_2) {
    id1 = id_1;
    id2 = id_2;
    geneset1.clear();
    geneset1.insert(geneset_1.begin(), geneset_1.end());
    geneset2.clear();
    geneset2.insert(geneset_2.begin(), geneset_2.end());
    treelen1 = treelen_1;
    treelen2 = treelen_2;
}

bool MergeJob::isEmpty() {
    return (id1 == -1);
}

void MergeJob::copyFrom(MergeJob* anotherMergeJob) {
    this->id1 = anotherMergeJob->id1;
    this->id2 = anotherMergeJob->id2;
    this->geneset1.clear();
    this->geneset1.insert(anotherMergeJob->geneset1.begin(), anotherMergeJob->geneset1.end());
    this->geneset2.clear();
    this->geneset2.insert(anotherMergeJob->geneset2.begin(), anotherMergeJob->geneset2.end());
    this->treelen1 = anotherMergeJob->treelen1;
    this->treelen2 = anotherMergeJob->treelen2;
}

void MergeJob::setEmpty() {
    id1 = id2 = -1;
    geneset1.clear();
    geneset2.clear();
    treelen1 = 0.0;
    treelen2 = 0.0;
}

void MergeJob::toString(string& str) {
    set<int>::iterator itr;
    str.clear();
    str.append(convertIntToString(id1) + ";");
    str.append(convertIntToString(id2) + ";");
    str.append(convertDoubleToString(treelen1)+";");
    str.append(convertDoubleToString(treelen2)+";");
    for (itr=geneset1.begin(); itr!=geneset1.end(); itr++) {
        if (itr!=geneset1.begin())
            str.append(",");
        str.append(convertIntToString(*itr));
    }
    str.append(";");
    for (itr=geneset2.begin(); itr!=geneset2.end(); itr++) {
        if (itr!=geneset2.begin())
            str.append(",");
        str.append(convertIntToString(*itr));
    }
    str.append(";");
}

void MergeJob::loadFrString(string& str) {
    
    int start_pos = 0;
    int pos = 0;
    
    // reset all variables
    setEmpty();
    
    // read id1
    while (pos < str.length() && str[pos] != ';')
        pos++;
    if (start_pos < pos && start_pos < str.length())
        id1 = atoi(str.substr(start_pos, pos - start_pos).c_str());
    pos++;
    start_pos = pos;
    
    // read id2
    while (pos < str.length() && str[pos] != ';')
        pos++;
    if (start_pos < pos && start_pos < str.length())
        id2 = atoi(str.substr(start_pos, pos - start_pos).c_str());
    pos++;
    start_pos = pos;
    
    // read treelen1
    while (pos < str.length() && str[pos] != ';')
        pos++;
    if (start_pos < pos && start_pos < str.length())
        treelen1 = atof(str.substr(start_pos, pos - start_pos).c_str());
    pos++;
    start_pos = pos;

    // read treelen2
    while (pos < str.length() && str[pos] != ';')
        pos++;
    if (start_pos < pos && start_pos < str.length())
        treelen2 = atof(str.substr(start_pos, pos - start_pos).c_str());
    pos++;
    start_pos = pos;
    
    // read geneset1
    while (pos < str.length() && str[pos] != ';') {
        if (str[pos] == ',') {
            if (start_pos < pos && start_pos < str.length())
                geneset1.insert(atoi(str.substr(start_pos, pos - start_pos).c_str()));
            start_pos = pos+1;
        }
        pos++;
    }
    if (start_pos < pos && start_pos < str.length())
        geneset1.insert(atoi(str.substr(start_pos, pos - start_pos).c_str()));
    pos++;
    start_pos = pos;

    // read geneset2
    while (pos < str.length() && str[pos] != ';') {
        if (str[pos] == ',') {
            if (start_pos < pos && start_pos < str.length())
                geneset2.insert(atoi(str.substr(start_pos, pos - start_pos).c_str()));
            start_pos = pos+1;
        }
        pos++;
    }
    if (start_pos < pos && start_pos < str.length())
        geneset2.insert(atoi(str.substr(start_pos, pos - start_pos).c_str()));
}

#endif // _IQTREE_MPI

// to check how many classes from the model string
int getClassNum(string model_str) {
    // the number of commas inside the model string + 1
    size_t pos = 0;
    int k = 0;
    pos = model_str.find_first_of(',',pos);
    while (pos != string::npos) {
        k++;
        pos++;
        pos = model_str.find_first_of(',',pos);
    }
    return k+1;
}

// get the k-th class model
// model_str should not contain the RHAS model
string classKModel(string model_str, int k) {
    int n = getClassNum(model_str);
    if (k >= n)
        return "";
    if (n == 1) {
        return model_str;
    }
    
    int j = 0;
    size_t pos = 0;
    size_t pos_fr;
    pos = model_str.find_first_of('{');
    while (j < k) {
        pos++;
        pos = model_str.find_first_of(',',pos);
        j++;
    }
    pos++;
    pos_fr = pos;
    pos = model_str.find_first_of("},",pos);
    return model_str.substr(pos_fr,pos-pos_fr);
}

// assign new substitution to the k-th class
// return false if there are less than k classes in the model
bool changeModel(string model_str, string& new_model_str, string new_subst, int k) {
    int n = getClassNum(model_str);
    if (k >= n)
        return false;
    if (n == 1) {
        new_model_str = new_subst;
        return true;
    }
    
    int j = 0;
    size_t pos = 0;
    string left_part, right_part;
    pos = model_str.find_first_of('{');
    while (j < k) {
        pos++;
        pos = model_str.find_first_of(',',pos);
        j++;
    }
    left_part = model_str.substr(0, pos+1);
    pos++;
    pos = model_str.find_first_of("},",pos);
    right_part = model_str.substr(pos);
    new_model_str = left_part + new_subst + right_part;
    return true;
}

// assign a new substitution into a q-mixture model
void addModel(string model_str, string& new_model_str, string new_subst) {
    size_t pos;
    int n;
    n = getClassNum(model_str);
    if (n == 1) {
        new_model_str = "MIX{" + model_str + "," + new_subst + "}";
    } else {
        pos = model_str.find_last_of('}');
        new_model_str = model_str.substr(0, pos) + "," + new_subst + model_str.substr(pos);
    }
}

// This function is similar to runModelFinder, but it is designed for optimisation of Q-Mixture model
// action: 1 - estimate the RHAS model
//         2 - estimate the number of classes in a mixture model
//         3 - estimate the k-th substitution matrix
//         4 - estimate an additional substitution matrix
CandidateModel runModelSelection(Params &params, IQTree &iqtree, ModelCheckpoint &model_info, int action, bool do_init_tree, string model_str, string& best_subst_name, string& best_rate_name, map<string, vector<string> > nest_network, int class_k = 0)
{
    double cpu_time;
    double real_time;
    bool ok_model_file;
    int partition_type;
    Checkpoint *orig_checkpoint;
    ModelsBlock *models_block;
    CandidateModelSet candidate_models;
    CandidateModel best_model;
    string multi_class_str;
    string single_class_str;
    int max_cats;
    string set_name = "";
    string in_model_name = "";
    bool merge_phase = false;
    bool generate_candidates;
    bool skip_all_when_drop;
    string orig_model_set;
    string orig_ratehet_set;
    vector<string> model_names;
    vector<string> ratehet;
    vector<string> freq_names;
    int i,j;
    
    // timing
    cpu_time = getCPUTime();
    real_time = getRealTime();
    
    // handling checkpoint file
    model_info.setFileName((string)params.out_prefix + ".model.gz");
    model_info.setDumpInterval(params.checkpoint_dump_interval);
    ok_model_file = false;
    if (!params.model_test_again) {
        ok_model_file = model_info.load();
    }
    cout << endl;
    ok_model_file &= model_info.size() > 0;
    if (ok_model_file)
        cout << "NOTE: Restoring information from model checkpoint file " << model_info.getFileName() << endl;
    orig_checkpoint = iqtree.getCheckpoint();
    iqtree.setCheckpoint(&model_info);
    iqtree.restoreCheckpoint();
    if (CKP_RESTORE2((&model_info), partition_type)) {
        if (partition_type != params.partition_type)
            outError("Mismatch partition type between checkpoint and partition file command option\nRerun with -mredo to ignore .model.gz checkpoint file");
    } else {
        partition_type = params.partition_type;
        CKP_SAVE2((&model_info), partition_type);
    }
    
    models_block = readModelsDefinition(params);
    
    if (do_init_tree) {
        // compute initial tree
        iqtree.computeInitialTree(params.SSE);
        iqtree.saveCheckpoint();
    }
    
    max_cats = getClassNum(model_str) * params.max_rate_cats;
    
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
    
    // generate candidate models
    // setting the params
    orig_ratehet_set = params.ratehet_set;
    orig_model_set = params.model_set;
    
    // params.model_extra_set = NULL;
    // params.model_subset = NULL;
    // params.state_freq_set = NULL;
    generate_candidates = false;
    candidate_models.nest_network = nest_network;

    if (action == 1) {
        params.model_set = model_str;
        getRateHet(iqtree.aln->seq_type, params.model_name, iqtree.aln->frac_invariant_sites, params.ratehet_set, ratehet);

        // add number of rate cateogories for special rate models
        const char *rates[] = {"+R", "*R", "+H", "*H"};

        size_t pos;

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
                }
            }
        
        for (i=0; i<ratehet.size(); i++) {
            candidate_models.push_back(CandidateModel(model_str, ratehet[i], iqtree.aln, 0));
        }

        skip_all_when_drop = false;
    } else if (action == 2) {
        params.ratehet_set = iqtree.getModelFactory()->site_rate->name;
        // generate candidate models for the possible mixture models
        multi_class_str = "";
        single_class_str = model_str;
        for (i = 1; i <= params.max_mix_cats; i++) {
            if (!multi_class_str.empty())
                multi_class_str.append(",");
            multi_class_str.append(single_class_str);
            if (i >= params.min_mix_cats) {
                if (i > 1)
                    model_str = "MIX{" + multi_class_str + "}";
                else
                    model_str = multi_class_str;
                candidate_models.push_back(CandidateModel(model_str, iqtree.getModelFactory()->site_rate->name, iqtree.aln, 0));
            }
        }
        skip_all_when_drop = true;
    } else if (action == 3) {
        char init_state_freq_set[] = "FO";
        if (!params.state_freq_set) {
            params.state_freq_set = init_state_freq_set;
        }
        params.ratehet_set = iqtree.getModelFactory()->site_rate->name;
        getModelSubst(iqtree.aln->seq_type, iqtree.aln->isStandardGeneticCode(), params.model_name,
                      params.model_set, params.model_subset, model_names);
        
        if (model_names.empty())
            return best_model;
        
        getStateFreqs(iqtree.aln->seq_type, params.state_freq_set, freq_names);
        
        // combine model_names with freq_names
        if (freq_names.size() > 0) {
            StrVector orig_model_names = model_names;
            model_names.clear();
            for (j = 0; j < orig_model_names.size(); j++) {
                if (iqtree.aln->seq_type == SEQ_CODON) {
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

        for (i=0; i<model_names.size(); i++) {
            string new_model_str;
            if (changeModel(model_str, new_model_str, model_names[i], class_k))
                candidate_models.push_back(CandidateModel(new_model_str, iqtree.getModelFactory()->site_rate->name, iqtree.aln, 0));
        }

        skip_all_when_drop = false;
    } else {
        params.ratehet_set = best_rate_name;
        
        getModelSubst(iqtree.aln->seq_type, iqtree.aln->isStandardGeneticCode(), params.model_name,
                      params.model_set, params.model_subset, model_names);
        
        if (model_names.empty())
            return best_model; // empty
        
        getStateFreqs(iqtree.aln->seq_type, params.state_freq_set, freq_names);
        
        // combine model_names with freq_names
        if (freq_names.size() > 0) {
            StrVector orig_model_names = model_names;
            model_names.clear();
            for (j = 0; j < orig_model_names.size(); j++) {
                if (iqtree.aln->seq_type == SEQ_CODON) {
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

        for (i=0; i<model_names.size(); i++) {
            string new_model_str;
            addModel(model_str, new_model_str, model_names[i]);
            candidate_models.push_back(CandidateModel(new_model_str, best_rate_name, iqtree.aln, 0));
        }

        skip_all_when_drop = false;

        if (candidate_models.size() > 0) {
            candidate_models.at(0).init_first_mix = true;
        }
    }
    // model selection
    candidate_models.under_mix_finder = true;
    best_model = candidate_models.test(params, &iqtree, model_info, models_block, params.num_threads, BRLEN_OPTIMIZE,
                                       set_name, in_model_name, merge_phase, generate_candidates, skip_all_when_drop);
    
    iqtree.aln->model_name = best_model.getName();
    best_subst_name = best_model.subst_name;
    best_rate_name = best_model.rate_name;
    
    Checkpoint *checkpoint = &model_info;
    string best_model_AIC, best_model_AICc, best_model_BIC;
    CKP_RESTORE(best_model_AIC);
    CKP_RESTORE(best_model_AICc);
    CKP_RESTORE(best_model_BIC);
    //cout << "Akaike Information Criterion:           " << best_model_AIC << endl;
    //cout << "Corrected Akaike Information Criterion: " << best_model_AICc << endl;
    //cout << "Bayesian Information Criterion:         " << best_model_BIC << endl;
    cout << "Best-fit model: " << iqtree.aln->model_name << " chosen according to "
        << criterionName(params.model_test_criterion) << endl;

    // remove key "OptModel" from the checkpoint file, which is only used for initialising models from the nested models.
    iqtree.getCheckpoint()->eraseKeyPrefix("OptModel");

    delete models_block;
    
    // force to dump all checkpointing information
    model_info.dump(true);
    
    // transfer models parameters
    transferModelFinderParameters(&iqtree, orig_checkpoint);
    iqtree.setCheckpoint(orig_checkpoint);

    params.model_set = orig_model_set;
    params.ratehet_set = orig_ratehet_set;
    // params.startCPUTime = cpu_time;
    // params.start_real_time = real_time;
    cpu_time = getCPUTime() - cpu_time;
    real_time = getRealTime() - real_time;
    cout << endl;
    cout << "All model information printed to " << model_info.getFileName() << endl;
    cout << "CPU time for ModelFinder: " << cpu_time << " seconds (" << convert_time(cpu_time) << ")" << endl;
    cout << "Wall-clock time for ModelFinder: " << real_time << " seconds (" << convert_time(real_time) << ")" << endl;
    
    return best_model;
}

// Optimisation of Q-Mixture model, including estimation of best number of classes in the mixture
// Method updated
void optimiseQMixModel_method_update(Params &params, IQTree* &iqtree, ModelCheckpoint &model_info, string& model_str) {

    bool do_init_tree;
    string best_subst_name;
    string best_rate_name;
    int action, best_class_num, i;
    set<string> skip_models;
    string model_str1, model_i;
    ModelsBlock *models_block;
    CandidateModel best_model;
    string best_model_AIC, best_model_AICc, best_model_BIC;
    double best_score_AIC, best_score_AICc, best_score_BIC;
    Checkpoint *checkpoint;
    int ssize;
    int curr_df;
    double curr_loglike;
    double curr_score;
    bool better_model;
    double LR, df_diff, pvalue;
    string criteria_str;

    char init_state_freq_set[] = "FO";
    if (!params.state_freq_set) {
        params.state_freq_set = init_state_freq_set;
    }

    models_block = readModelsDefinition(params);
    ssize = iqtree->getAlnNSite();
    criteria_str = criterionName(params.model_test_criterion);

    // Step 0: (reorder candidate DNA models when -mset is used) build the nest-relationship network
    SeqType seq_type = iqtree->aln->seq_type;
    map<string, vector<string> > nest_network;
    if (seq_type == SEQ_DNA) {
        StrVector model_names, freq_names;
        getModelSubst(iqtree->aln->seq_type, iqtree->aln->isStandardGeneticCode(), params.model_name,
                      params.model_set, params.model_subset, model_names);
        getStateFreqs(iqtree->aln->seq_type, params.state_freq_set, freq_names);

        nest_network = generateNestNetwork(model_names, freq_names);
    } else {
        nest_network = {};
    }


    // Step 1: run ModelFinder
    params.model_name = "";
    bool under_mix_finder = true;
    runModelFinder(params, *iqtree, model_info, best_subst_name, best_rate_name, nest_network, under_mix_finder);

    // (cancel) Step 2: do tree search for this single-class model
    // runTreeReconstruction(params, iqtree);
    // curr_df = iqtree->getModelFactory()->getNParameters(BRLEN_OPTIMIZE);
    // curr_loglike = iqtree->getCurScore();
    // curr_score = computeInformationScore(curr_loglike, curr_df, ssize, params.model_test_criterion);
    string best_model_logl_df = model_info[best_subst_name+best_rate_name];
    stringstream ss (best_model_logl_df);
    ss >> curr_loglike >> curr_df;
    string best_score = model_info["best_score_" + criteria_str];
    curr_score = convert_double(best_score.c_str());

    cout << endl << "Model: " << best_subst_name << best_rate_name << "; df: " << curr_df << "; loglike: " << curr_loglike << "; " << criteria_str << " score: " << curr_score << endl;
    
    // Step 3: keep adding a new class until no further improvement
    if (params.opt_qmix_criteria == 1) {
        cout << endl << "Keep adding an additional class until the p-value from the likelihood ratio test > " << params.opt_qmix_pthres << endl;
    } else {
        cout << endl << "Keep adding an additional class until there is no better " << criteria_str <<  " value" << endl;
    }
    action = 4;
    do_init_tree = false;
    model_str = best_subst_name;
    do {
        best_model = runModelSelection(params, *iqtree, model_info, action, do_init_tree, model_str, best_subst_name, best_rate_name, nest_network);
        cout << endl << "Model: " << best_subst_name << best_rate_name << "; df: " << best_model.df << "; loglike: " << best_model.logl << "; " << criteria_str << " score: " << best_model.getScore() << ";";
        if (params.opt_qmix_criteria == 1) {
            LR = 2.0 * (best_model.logl - curr_loglike);
            df_diff = best_model.df - curr_df;
            pvalue = computePValueChiSquare(LR, df_diff);
            better_model = (pvalue <= params.opt_qmix_pthres);
            cout << " pvalue: " << pvalue << "; ";
        } else {
            // compare the models based on IC score
            better_model = (best_model.getScore() < curr_score);
        }
        cout << endl;
        if (better_model) {
            curr_df = best_model.df;
            curr_loglike = best_model.logl;
            curr_score = best_model.getScore();
            model_str = best_subst_name;
        }
    } while (better_model && getClassNum(best_subst_name)+1 <= params.max_mix_cats);
    
    best_subst_name = model_str;
    
    if (params.opt_rhas_again) {
        // Step 4: estimate the RHAS model again
        action = 1; // estimating the RHAS model
        do_init_tree = false;
        model_str = best_subst_name;
        best_model = runModelSelection(params, *iqtree, model_info, action, do_init_tree, model_str, best_subst_name, best_rate_name, nest_network);
        curr_df = best_model.df;
        curr_loglike = best_model.logl;
        curr_score = best_model.getScore();
    }

    model_str = best_subst_name+best_rate_name;
}

// Optimisation of Q-Mixture model, including estimation of best number of classes in the mixture
void optimiseQMixModel(Params &params, IQTree* &iqtree, ModelCheckpoint &model_info) {

    IQTree* new_iqtree;
    string model_str;

    if (params.model_name.substr(0,6) != "MIX+MF")
        return;
    
    bool test_only = (params.model_name == "MIX+MF");
    params.model_name = "";
    
    if (MPIHelper::getInstance().getNumProcesses() > 1)
        outError("Error! The option -m '" + params.model_name + "' does not support MPI parallelization");
    
    if (iqtree->isSuperTree())
        outError("Error! The option -m '" + params.model_name + "' cannot work on data set with partitions");
    
    if (iqtree->aln->seq_type != SEQ_DNA)
        outError("Error! The option -m '" + params.model_name + "' can only work on DNA data set");

    cout << "--------------------------------------------------------------------" << endl;
    cout << "|                Optimizing Q-mixture model                        |" << endl;
    cout << "--------------------------------------------------------------------" << endl;

    // disable the bootstrapping
    int orig_gbo_replicates = params.gbo_replicates;
    ConsensusType orig_consensus_type = params.consensus_type;
    STOP_CONDITION orig_stop_condition = params.stop_condition;
    params.gbo_replicates = 0;
    params.consensus_type = CT_NONE;
    params.stop_condition = SC_UNSUCCESS_ITERATION;

    optimiseQMixModel_method_update(params, iqtree, model_info, model_str);
    
    // restore the original values
    params.gbo_replicates = orig_gbo_replicates;
    params.consensus_type = orig_consensus_type;
    params.stop_condition = orig_stop_condition;

    cout << "-------------------------------------------------------" << endl;
    cout << "  Best-fit Q-Mixture model: " << model_str << endl;
    cout << "-------------------------------------------------------" << endl;

    params.model_name = model_str;
    iqtree->aln->model_name = model_str;

    // create a new IQTree object for this mixture model
    // allocate heterotachy tree if neccessary
    int pos = posRateHeterotachy(iqtree->aln->model_name);
    if (params.num_mixlen > 1) {
        new_iqtree = new PhyloTreeMixlen(iqtree->aln, params.num_mixlen);
    } else if (pos != string::npos) {
        new_iqtree = new PhyloTreeMixlen(iqtree->aln, 0);
    } else {
        new_iqtree = new IQTree(iqtree->aln);
    }
    new_iqtree->setCheckpoint(iqtree->getCheckpoint());
    new_iqtree->setParams(&params);
    delete(iqtree);
    iqtree = new_iqtree;
    
    if (test_only) {
        params.min_iterations = 0;
    }
}

/****************************************************/
/*    Q MATRICES NESTING CHECK                      */
/****************************************************/

int findModelIndex(const string& model) {
    int i;
    for (i = 0; i < sizeof(dna_model_names) / sizeof(dna_model_names[0]); i++) {
        if (strcmp(dna_model_names[i], model.c_str()) == 0) {
            return i;
        }
    }
    return -1;
}

StrVector reorderModelNames(StrVector model_names) {
    int i;
    struct model_index {
        string name;
        int index;
        model_index(string a_name, int an_index):name(a_name),index(an_index){}
    };
    struct sort_by_ind {
        bool operator() (const model_index& a, const model_index& b) const {return a.index < b.index;}
    };
    vector<model_index> mi;

    for (i = 0; i < model_names.size(); i++) {
        if (findModelIndex(model_names[i]) == -1) {
            outError("Incorrect model name: " + model_names[i] + " is input");
        }
        mi.push_back(model_index(model_names[i], findModelIndex(model_names[i])));
    }
    sort(mi.begin(), mi.end(), sort_by_ind());
    for (i = 0; i < model_names.size(); i++) {
        model_names[i] = mi[i].name;
    }
    return model_names;
}

bool isRateTypeNested(string rate_type1, string rate_type2) {
    if (rate_type1.length() != 6) {
        outError("Incorrect DNA model rate type code: " + rate_type1);
    }
    if (rate_type2.length() != 6) {
        outError("Incorrect DNA model rate type code: " + rate_type2);
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++ ){
            if (rate_type1[i] == rate_type1[j] && rate_type2[i] != rate_type2[j]){
                return false;
            }
        }
    }
    return true;
}

map<string, vector<string> > generateNestNetwork(StrVector model_names, StrVector freq_names) {
    int i, j, k;
    bool covered;
    string full_name1, full_name2, rate_type1, rate_type2;
    StateFreqType freq1, freq2;
    vector<string> mfname_and_mname, nested_models, nested_models_all;
    map<string, vector<string> > nest_network, nest_network_all;
    vector<StrVector> model_freq_names;

    model_freq_names = {};
    for (i = 0; i < model_names.size(); i++) {
        string new_model_name = getDNAModelInfo(model_names[i], full_name1, rate_type1, freq1);
        if (model_names[i] != new_model_name)
            model_names[i] = new_model_name;

        if (freq1 == FREQ_EQUAL) {
            mfname_and_mname = {};
            mfname_and_mname.push_back(model_names[i]);
            mfname_and_mname.push_back(model_names[i]);
            mfname_and_mname.push_back("+FQ");
            model_freq_names.push_back(mfname_and_mname);
        } else if (freq1 == FREQ_ESTIMATE) {
            for (j = 0; j < freq_names.size(); j++) {
                mfname_and_mname = {};
                if (freq_names[j] == "+FQ") {
                    mfname_and_mname.push_back(model_names[i]);
                } else {
                    mfname_and_mname.push_back(model_names[i] + freq_names[j]);
                }
                mfname_and_mname.push_back(model_names[i]);
                mfname_and_mname.push_back(freq_names[j]);
                model_freq_names.push_back(mfname_and_mname);
            }
        }
    }

    nest_network[model_freq_names[0][0]] = {};
    nest_network_all[model_freq_names[0][0]] = {};
    for (i = 1; i < model_freq_names.size(); i++) {
        nested_models = {};
        nested_models_all = {};
        for (j = nest_network.size()-1; j >= 0; j--) {
            string result1 = getDNAModelInfo(model_freq_names[i][1], full_name1, rate_type1, freq1);
            string result2 = getDNAModelInfo(model_freq_names[j][1], full_name2, rate_type2, freq2);
            if (model_freq_names[i][2] != "+FO" && model_freq_names[j][2] != model_freq_names[i][2]) {
                continue;
            } else {
                if (isRateTypeNested(rate_type1, rate_type2)) {
                    nested_models_all.push_back(model_freq_names[j][0]);
                    covered = false;
                    if (!nested_models.empty()) {
                        for (k = 0; k < nested_models.size(); k++) {
                            vector<string> nested_nested_models = nest_network_all[nested_models[k]];

                            auto it = find(nested_nested_models.begin(), nested_nested_models.end(), model_freq_names[j][0]);
                            if (it != nested_nested_models.end()) {
                                covered = true;
                            }
                        }
                    }
                    if (!covered) {
                        nested_models.push_back(model_freq_names[j][0]);
                    }
                }
            }
        }
        nest_network[model_freq_names[i][0]] = nested_models;
        nest_network_all[model_freq_names[i][0]] = nested_models_all;
    }
    return nest_network;
}