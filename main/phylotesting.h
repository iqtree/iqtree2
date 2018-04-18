/*
 * phylotesting.h
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */

#ifndef PHYLOTESTING_H_
#define PHYLOTESTING_H_

#include "utils/tools.h"
#include "utils/checkpoint.h"

class PhyloTree;
class IQTree;


class ModelInfo {
public:
	string set_name; // subset name
	string name; // model name
	double logl; // tree log likelihood
	int df;      // #parameters
    double tree_len; // tree length, added 2015-06-24 for rcluster algorithm
    string tree; // added 2015-04-28: tree string
	double AIC_score, AICc_score, BIC_score;    // scores
	double AIC_weight, AICc_weight, BIC_weight; // weights
	bool AIC_conf, AICc_conf, BIC_conf;         // in confidence set?

    /**
        compute information criterion scores (AIC, AICc, BIC)
    */
    void computeICScores(size_t sample_size);

    /**
        compute information criterion scores (AIC, AICc, BIC)
    */
    double computeICScore(size_t sample_size);

    /**
        save model into checkpoint
    */
    void saveCheckpoint(Checkpoint *ckp) {
        stringstream ostr;
        ostr.precision(10);
        ostr << logl << " " << df << " " << tree_len;
        if (!tree.empty())
            ostr << " " << tree;
        ckp->put(name, ostr.str());
    }

    /**
        restore model from checkpoint
    */
    bool restoreCheckpoint(Checkpoint *ckp) {
        string val;
        if (ckp->getString(name, val)) {
            stringstream str(val);
            str >> logl >> df >> tree_len;
            return true;
        }
        return false;
    }

    /**
        restore model from checkpoint
    */
    bool restoreCheckpointRminus1(Checkpoint *ckp, string &model_name) {
        size_t posR;
        const char *rates[] = {"+R", "*R", "+H", "*H"};
        for (int i = 0; i < sizeof(rates)/sizeof(char*); i++) {
            if ((posR = model_name.find(rates[i])) != string::npos) {
                int cat = convert_int(model_name.substr(posR+2).c_str());
                name = model_name.substr(0, posR+2) + convertIntToString(cat-1);
                return restoreCheckpoint(ckp);
            }
        }
        return false;
    }

};

//typedef vector<ModelInfo> ModelCheckpoint;

class ModelCheckpoint : public Checkpoint {

public:

    /*
        get the best model
        @param[out] best_model name of the best model
        @return TRUE if best model found, FALSE otherwise (unfinished job)
    */
    bool getBestModel(string &best_model);

    /*
        get the ordered model list according to AIC, AICc or BIC
        @param tree associated tree
        @param[out] ordered_models list of models ordered by specified criterion
        @return TRUE if ordered_models found, FALSE otherwise (unfinished job)
    */
    bool getOrderedModels(PhyloTree *tree, vector<ModelInfo> &ordered_models);

    /*
        get the best tree
        @param[out] best_tree NEWICK string of the best tree
        @return TRUE if best tree found, FALSE otherwise (unfinished job)
    */
    bool getBestTree(string &best_tree);

};

struct TreeInfo {
	double logl; // log likelihood
	double se; // standard error of deltaL (logl difference to max), or square root of variance
	double rell_bp; // bootstrap proportion by RELL method
	bool rell_confident; // confidence set for RELL-BP
	double sh_pvalue; // p-value by Shimodaira-Hasegawa test
	double wsh_pvalue; // p-value by weighted Shimodaira-Hasegawa test
	double kh_pvalue; // p-value by Kishino-Hasegawa test
	double wkh_pvalue; // p-value by weighted Kishino-Hasegawa test
	double elw_value; // ELW - expected likelihood weights test
	bool elw_confident; // to represent confidence set of ELW test
    double au_pvalue; // p-value by approximately unbiased (AU) test
};


/**
 * computing AIC, AICc, and BIC scores
 */
void computeInformationScores(double tree_lh, int df, int ssize, double &AIC, double &AICc, double &BIC);

double computeInformationScore(double tree_lh, int df, int ssize, ModelTestCriterion mtc);

string criterionName(ModelTestCriterion mtc);

/**
 * check if the model file contains correct information
 * @param model_file model file names
 * @param model_name (OUT) vector of model names
 * @param lh_scores (OUT) vector of tree log-likelihoods
 * @param df_vec (OUT) vector of degrees of freedom (or K)
 * @return TRUE if success, FALSE failed.
 */

bool checkModelFile(string model_file, bool is_partitioned, ModelCheckpoint &infos);

/**
 perform ModelFinder to find the best-fit model
 @param params program parameters
 @param iqtree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
 */
void runModelFinder(Params &params, IQTree &iqtree, ModelCheckpoint &model_info);

/**
 testing the best-fit model
 return in params.freq_type and params.rate_type
 @param set_name for partitioned analysis
 @param in_tree phylogenetic tree
 @param model_info (IN/OUT) information for all models considered
 @param set_name for partition model selection
 @param print_mem_usage true to print RAM memory used (default: false) 
 @return name of best-fit-model
 */
//string testModel(Params &params, PhyloTree* in_tree, ModelCheckpoint &model_info,
//		ModelsBlock *models_block, int num_threads, int brlen_type,
//        string set_name = "", bool print_mem_usage = false, string in_model_name = "");

/**
 * print site log likelihoods to a fileExists
 * @param filename output file name
 * @param tree phylogenetic tree
 * @param ptn_lh pattern log-likelihoods, will be computed if NULL
 * @param append TRUE to append to existing file, FALSE otherwise
 * @param linename name of the line, default "Site_Lh" if NULL
 */
void printSiteLh(const char*filename, PhyloTree *tree, double *ptn_lh = NULL,
		bool append = false, const char *linename = NULL);

/**
 * print partition log likelihoods to a file
 * @param filename output file name
 * @param tree phylogenetic tree
 * @param ptn_lh pattern log-likelihoods, will be computed if NULL
 * @param append TRUE to append to existing file, FALSE otherwise
 * @param linename name of the line, default "Site_Lh" if NULL
 */
void printPartitionLh(const char*filename, PhyloTree *tree, double *ptn_lh = NULL,
		bool append = false, const char *linename = NULL);

/**
 * print site log likelihoods per category to a file
 * @param filename output file name
 * @param tree phylogenetic tree
 */
void printSiteLhCategory(const char*filename, PhyloTree *tree, SiteLoglType wsl);

/**
 * print site posterior probabilities per rate/mixture category to a file
 * @param filename output file name
 * @param tree phylogenetic tree
 */
void printSiteProbCategory(const char*filename, PhyloTree *tree, SiteLoglType wsl);

/**
 * print site state frequency vectors (for Huaichun)
 * @param filename output file name
 * @param tree phylogenetic tree
*/
void printSiteStateFreq(const char*filename, PhyloTree *tree, double *state_freqs = NULL);

/**
 * print site state frequency vectors (for Huaichun)
 * @param filename output file name
 * @param aln alignment
*/
void printSiteStateFreq(const char* filename, Alignment *aln);

/**
    print ancestral sequences
    @param filename output file name
    @param tree phylogenetic tree
    @param ast either AST_MARGINAL or AST_JOINT
*/
void printAncestralSequences(const char*filename, PhyloTree *tree, AncestralSeqType ast);

/**
 * Evaluate user-trees with possibility of tree topology tests
 * @param params program parameters
 * @param tree current tree
 * @param info (OUT) output information
 * @param distinct_ids IDs of distinct trees
 */
void evaluateTrees(Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids);

void evaluateTrees(Params &params, IQTree *tree);

/**
    get sequence type for a model name
    @param model_name model name string
    @param seq_type (OUT) sequence type, SEQ_UNKNOWN if is not determined
    @return 1 for parametric model, 2 for empirical model
*/
int detectSeqType(const char *model_name, SeqType &seq_type);

string detectSeqTypeName(string model_name);


#endif /* PHYLOTESTING_H_ */
