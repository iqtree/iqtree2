#ifndef NNISEARCH_H
#define NNISEARCH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "pll/pll.h"

const int TOPO_ONLY = 0;
const int NO_BRAN_OPT = 1;
const int ONE_BRAN_OPT = 2;
const int FIVE_BRAN_OPT = 4;

#define IQTREE_NEWZPERCYCLE		10

/* This is the info you need to copy the vector*/
typedef struct
{
  int node_number;
  int num_partitions;
  size_t *partition_sizes;
  double **lh_values;
  int **expVector;
} LH_VECTOR;


typedef struct {
	//pllInstance* tr;
	nodeptr p;
	int nniType;
    double z0[PLL_NUM_BRANCHES]; // p
    double z1[PLL_NUM_BRANCHES]; // p->next
    double z2[PLL_NUM_BRANCHES]; // p->next->next
    double z3[PLL_NUM_BRANCHES]; // q->next
    double z4[PLL_NUM_BRANCHES]; // q->next->next
	double likelihood;
} pllNNIMove;

LH_VECTOR backup_likelihood_pointers(pllInstance *tr, partitionList *pr, nodeptr p);

static int cmp_nni(const void* nni1, const void* nni2);

int compareDouble(const void * a, const void * b);

pllNNIMove *getNonConflictNNIList(pllInstance* tr);

void _update(pllInstance *tr, partitionList *pr, nodeptr p);

/**
 * TODO: read tree from string in memory
 */
 int treeReadLenString (const char *buffer, pllInstance *tr, pll_boolean readBranches, pll_boolean readNodeLabels, pll_boolean topologyOnly);

#define MAX_NUM_DELTA 10000

typedef struct {
	double delta[MAX_NUM_DELTA];
	int num_delta;
	double delta_min;
	int doNNICut;
} NNICUT;


/*
 *  Evaluate NNI moves for the current internal branch
 *  @param tr the current tree data structure
 *  @param pr partition data structure
 *  @param p the node representing the current branch
 *  @param nniList array containing evaluated NNI moves
 *  @param numBran number of branches that have been evaluated
 *  @param numPosNNI
 *  @param curLH the curren log-likelihood of the tree
 *  @return 1 if a positive NNI is found, 0 otherwise
 */
int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p,  pllNNIMove* nniList, int searchType, int* numBran, double curLH);

/**
 * Perturb the best tree
 *
 * Given the best tree, apply some NNIs to escape local optimum
 * @param tr
 * @param pr
 * @param nnis list of all NNI to apply
 * @param numNNI size of the array nnis
 * @return
 */
double perturbTree(pllInstance *tr, partitionList *pr, pllNNIMove *nnis, int numNNI);

/**
 * 	do 1 round of fastNNI
 *  return new tree log-likelihood if found improving NNI otherwise -1.0
 *
 *  @param[in] tr the tree data structure
 *  @param[in] pr partition data structure
 *  @param[out] nni_count number of NNI that have been applied
 *  @param[out] deltaNNI average improvement made by one NNI
 */
double doNNISearch(pllInstance* tr, partitionList *pr, topol* curTree, pllNNIMove* nniList, int searchType, int* nni_count, double* deltaNNI);

/**
 *  perturb the current tree by randomly carrying some negative NNI moves
 *  @param[in] tr the tree
 *  @param[in] nniList list of all possible NNIs
 */
double pertub(pllInstance* tr, pllNNIMove* nniList);

//void optimizeOneBranches(pllInstance* tr, nodeptr p, int numNRStep);

/**
 *  @brief Do 1 NNI move.
 *  @param[in] tr: the tree data structure
 *  @param[in] pr partition data structure
 *  @param[in] swap: represents one of the 2 NNI moves. Could be either 0 or 1
 *  @param[in] evalType: NO_NR, WITH_ONE_NR, WITH_FIVE_NR
 */
double doOneNNI(pllInstance * tr, partitionList *pr, nodeptr p, int swap, int evalType, double curLH);

/**
 *  Go through all 2(n-3) internal branches of the tree and
 *  evaluate all possible NNI moves
 */
void evalAllNNI(pllInstance* tr);

/*
 * 	@brief evaluate all NNIs within the subtree specified by node p
 * 	populates the list containing all possible NNI moves
 *
 * 	@param[in] tr: the tree data structure
 * 	@param[in] pr partition data structure
 * 	@param[in] p node pointer that specify the subtree
 * 	@param[out] nniList list of evaluated NNI moves
 * 	@param[out] numBran number of internal branches that have been visited
 *  @param[out] numPosNNI number of positive NNI found
 */
void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, pllNNIMove* nniList, int searchType, int* numBran, int* numPosNNI, double curLH);

/*
 *  @brief return the array which can be used to store evaluated NNI moves
 *
 *  @param[in] tr: the tree data structure
 */
pllNNIMove *getNNIList(pllInstance* tr);

/*
 *  Save the likelihood vector of p and q to the 2 pointer p_lhsave and
 *  q_lhsave.
 *  Should I use memcpy or just copy the pointer ?
 */
//void saveLHVector(nodeptr p, nodeptr q, double* p_lhsave, double* q_lhsave);




#ifdef __cplusplus
}
#endif

#endif

