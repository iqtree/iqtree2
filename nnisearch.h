#ifndef NNISEARCH_H
#define NNISEARCH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

const int TOPO_ONLY = 0;
const int NO_BRAN_OPT = 1;
const int ONE_BRAN_OPT = 2;
const int FIVE_BRAN_OPT = 4;


/*--------------------------------------------------------------*/
/* portable version for fmemopenevalType */
/*--------------------------------------------------------------*/

#if defined __APPLE__ || defined __MACH__
struct fmem {
    size_t pos;
    size_t size;
    char *buffer;
};
typedef struct fmem fmem_t;

/* simple, but portable version of fmemopen for OS X / BSD */
FILE * fmemopen(void *buf, size_t size, const char *mode);

#endif /* APPLE */


/**
 * TODO: read tree from string in memory
 */
 int treeReadLenString (const char *buffer, tree *tr, pl_boolean readBranches, pl_boolean readNodeLabels, pl_boolean topologyOnly);

#define MAX_NUM_DELTA 10000

typedef struct {
	double delta[MAX_NUM_DELTA];
	int num_delta;
	double delta_min;
	int doNNICut;
} NNICUT;

typedef struct {
	tree* tr;
	nodeptr p;
	int nniType;
    double z0[NUM_BRANCHES]; // p
    double z1[NUM_BRANCHES]; // p->next
    double z2[NUM_BRANCHES]; // p->next->next
    double z3[NUM_BRANCHES]; // q->next
    double z4[NUM_BRANCHES]; // q->next->next
	double likelihood;
} nniMove;



/*
 *  Evaluate NNI moves for the current internal branch
 *  @param tr the current tree data structure
 *  @param p the node representing the current branch
 *  @param nniList array containing evaluated NNI moves
 *  @param numBran number of branches that have been evaluated
 *  @param numPosNNI
 *  @param curLH the curren log-likelihood of the tree
 *  @return 1 if a positive NNI is found, 0 otherwise
 */
int evalNNIForBran(tree* tr, nodeptr p,  nniMove* nniList, int* numBran, double curLH, NNICUT* nnicut);

/**
 * 	do 1 round of fastNNI
 *  return new tree log-likelihood if found improving NNI otherwise -1.0
 *
 *  @param tr: the tree data structure
 *  @param nni_count: pointer to the number of NNI that has been apply (OUT parameter)
 *  @param deltaNNI: pointer to the average improvement made by one NNI (OUT parameters)
 */
double doNNISearch(tree* tr, int* nni_count, double* deltaNNI, NNICUT* nnicut, int numSmooth);

/**
 *  perturb the current tree by randomly carrying some negative NNI moves
 *  @param[in] tr the tree
 *  @param[in] nniList list of all possible NNIs
 */
double pertub(tree* tr, nniMove* nniList);

void optimizeOneBranches(tree* tr, nodeptr p, int numNRStep);

/**
 *  @brief Do 1 NNI move.
 *  @param[in] tr: the tree data structure
 *  @param[in] swap: represents one of the 2 NNI moves. Could be either 0 or 1
 *  @param[in] evalType: NO_NR, WITH_ONE_NR, WITH_FIVE_NR
 */
double doOneNNI(tree * tr, nodeptr p, int swap, int evalType);

/**
 *  Go through all 2(n-3) internal branches of the tree and
 *  evaluate all possible NNI moves
 */
void evalAllNNI(tree* tr);

/*
 * 	@brief evaluate all NNIs within the subtree specified by node p
 * 	populates the list containing all possible NNI moves
 * 	@param[in] p node pointer that specify the subtree
 * 	@param[out] nniList list of evaluated NNI moves
 * 	@param[out] numBran number of internal branches that have been visited
 *  @param[out] numPosNNI number of positive NNI found
 */
void evalNNIForSubtree(tree* tr, nodeptr p, nniMove* nniList, int* numBran, int* numPosNNI, double curLH, NNICUT* nnicut);

/*
 *  @brief return the array which can be used to store evaluated NNI moves
 */
nniMove *getNNIList(tree* tr);

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

