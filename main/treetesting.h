/*
 * treetesting.h
 *
 *  Created on: Sep 21, 2019
 *      Author: minh
 */

#ifndef TREETESTING_H_
#define TREETESTING_H_

#include "utils/tools.h"
#include "alignment/alignment.h"

class PhyloTree;
class IQTree;

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
void evaluateTrees(istream &in, Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids);

void evaluateTrees(string treeset_file, Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids);


void printTreeTestResults(vector<TreeInfo> &info, IntVector &distinct_ids, IntVector &branch_ids, string out_file);

#endif // TREETESTING_H_
