/*
 * upperbounds.h
 *
 *  Created on: Aug 13, 2014
 *      Author: olga
 */

#ifndef UPPERBOUNDS_H_
#define UPPERBOUNDS_H_

/**
	main function to carry out Upper Bounds analysis
*/
#include "tree/iqtree.h"
#include "tree/mexttree.h"
#include "alignment/alignment.h"
#include "tree/phylotree.h"

class PhyloTree;
class IQTree;

void UpperBounds(Params* params, Alignment* alignment, IQTree* tree);

void printUB();

void printTreeUB(MTree *tree);

/**
 * extracting subtree spanned by corresponding taxa together with subalignment
 *
 * ids  - vector of taxa ids, the subtree is spanned by these taxa
 * type - specifies the copying procedure:
 * 		= 0, normal (PhyloTree) tree->copyTree;
 * 			 two branches, incident to one of the end nodes of the branch spliting two subtrees, are collapsed
 * 		= 1, remember the end of the branch, do not collapse two branches
 * 			 this is used for the subtree with "artificial root",
 * 			 this will be a leaf without any nucleotide fixed at any site
 */
PhyloTree* extractSubtreeUB(IntVector &ids, MTree* tree, Params *params, int sw = 0);

// Slightly changed functions from usual MTree.copy, to allow for non collapsed branches: type=1
void copyTreeUB(MTree *tree, MTree *treeCopy, string &taxa_set);
Node* copyTreeUBnode(MTree *tree, MTree *treeCopy, string &taxa_set, double &len, Node *node = NULL, Node *dad = NULL);

// With artificial root internal node might get ID < taxaNUM. Just to make consistent with any other part of program,
// reindex all taxa from 0 to taxaNUM-1, and then internal nodes will get ids from taxaNUM to nodesNUM.
void reindexTaxonIDs(MTree *tree);


/*
 * Generate random YH tree, which is considered as a subtree of parent tree
 * input tree - is a subtree,
 * output - log-likelihood of generated subtree
 */
MTree* generateRandomYH_UB(Params &params, PhyloTree *tree);

/*
 * This function generates a random tree that has A|B split.
 * One can also specify the length of corresponding branch.
 * t(A|B) = brLen
 */
double RandomTreeAB(PhyloTree* treeORGN, PhyloTree* treeAorgn, PhyloTree* treeBorgn, IntVector &taxaA, IntVector &taxaB, Params *params, double brLen = 0.0);

/*
 * This function computes the product of logLhs of two subtrees corresponding to a given split A|B on tree.
 */
double UpperBoundAB(IntVector &taxaA, IntVector &taxaB, PhyloTree* tree, Params *params);

/*
 * Auxiliary function for RandomTreeAB, chooses randomly node and inserts new branch with randomLen
 */
void extendingTree(MTree* tree, Params* params);

NodeVector getBranchABid(double brLen, PhyloTree* tree);

/*
 * Applying UBs to NNI search
 */
NNIMove getBestNNIForBranUB(PhyloNode *node1, PhyloNode *node2, PhyloTree *tree);
double logC(double t, PhyloTree* tree);

/**
 * Tests on fractions ai/(ai+bi) and bi/(ai+bi)
 * (fractions of sums for matching and non-matching pairs of nucleotides on the ends of branch)
 */

void sumFraction(PhyloNode *node1, PhyloNode *node2, PhyloTree *tree);

#endif /* UPPERBOUNDS_H_ */
