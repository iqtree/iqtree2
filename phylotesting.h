/*
 * phylotesting.h
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */

#ifndef PHYLOTESTING_H_
#define PHYLOTESTING_H_

#include "tools.h"

class PhyloTree;
class IQTree;


struct ModelInfo {
	string set_name; // subset name
	string name; // model name
	double logl; // tree log likelihood
	int df;      // #parameters
	double AIC_score, AICc_score, BIC_score;    // scores
	double AIC_weight, AICc_weight, BIC_weight; // weights
	bool AIC_conf, AICc_conf, BIC_conf;         // in confidence set?
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
};


/**
 * computing AIC, AICc, and BIC scores
 */
void computeInformationScores(double tree_lh, int df, int ssize, double &AIC, double &AICc, double &BIC);

/**
 testing the best-fit model
 return in params.freq_type and params.rate_type
 @param set_name for partitioned analysis
 @return name of best-fit-model
 */
string testModel(Params &params, PhyloTree* in_tree, vector<ModelInfo> &model_info, string set_name = "");

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
 * print site log likelihoods per category to a file
 * @param filename output file name
 * @param tree phylogenetic tree
 */
void printSiteLhCategory(const char*filename, PhyloTree *tree);

/**
 * Evaluate user-trees with possibility of tree topology tests
 * @param params program parameters
 * @param tree current tree
 * @param info (OUT) output information
 * @param distinct_ids IDs of distinct trees
 */
void evaluateTrees(Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids);

void evaluateTrees(Params &params, IQTree *tree);

#endif /* PHYLOTESTING_H_ */
