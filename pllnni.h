/***************************************************************************
 *   Copyright (C) 2014 by                                            *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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

#ifndef NNISEARCH_H
#define NNISEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include "tools.h"
#include <string>
#include <sstream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
//#include <unordered_set>
extern "C" {
#include "pll/pllInternal.h"
}

typedef struct {
	nodeptr p;
	int nniType;
	char* idString;
    double z0[PLL_NUM_BRANCHES]; // p
    double z1[PLL_NUM_BRANCHES]; // p->next
    double z2[PLL_NUM_BRANCHES]; // p->next->next
    double z3[PLL_NUM_BRANCHES]; // q->next
    double z4[PLL_NUM_BRANCHES]; // q->next->next
	double likelihood;
	double loglDelta;
	double negLoglDelta;
} pllNNIMove;

typedef struct {
    // FOR GENERAL TREE SEARCH
	vector<pllNNIMove> posNNIList; // positive NNIs
	unordered_set<string> aBranches; // Set of branches that are affected by the previous NNIs
	double curLogl; // Current tree log-likelihood
	int curIter; // Current iteration number
	double curPerStrength; // Current perturbation strength

	// FOR NNI SEARCH
	NNI_Type nni_type;
	int numAppliedNNIs; // total number of applied NNIs sofar
	int curNumAppliedNNIs; // number of applied NNIs at the current step
	int curNumNNISteps;
	int maxNNISteps;
	set<double> deltaLogl; // logl differences between nni1 and nni5
} SearchInfo;


/* This is the info you need to copy the vector*/
typedef struct
{
  int node_number;
  int num_partitions;
  size_t *partition_sizes;
  double **lh_values;
  int **expVector;
} LH_VECTOR;

inline bool comparePLLNNIMove(const pllNNIMove &a, const pllNNIMove &b)
{
    return a.likelihood < b.likelihood;
}

void countDistinctTrees(pllInstance* pllInst, partitionList *pllPartitions);

//static int cmp_nni(const void* nni1, const void* nni2);

int compareDouble(const void * a, const void * b);

pllNNIMove *getNonConflictNNIList(pllInstance* tr);

void _update(pllInstance *tr, partitionList *pr, nodeptr p);

#define MAX_NUM_DELTA 10000

typedef struct {
	double delta[MAX_NUM_DELTA];
	int num_delta;
	double delta_min;
	int doNNICut;
} NNICUT;

/**
 * get all the nodes affected by the NNI
 * @param p
 * @return
 */
set<int> getAffectedNodes(pllInstance* tr, nodeptr p);

string getBranString(nodeptr p);

bool containsAffectedNodes(nodeptr p, SearchInfo &searchinfo);

void updateBranchLengthForNNI(pllInstance* tr, partitionList *pr, pllNNIMove &nni);

void pllEvalAllNNIs(pllInstance *tr, partitionList *pr, SearchInfo &searchinfo);

double pllDoRandomNNIs(pllInstance *tr, partitionList *pr, int numNNI);

/**
 *  Compute the possible best logl improvement of NNI5 over NNI1
 *  @param serachinfo contains the logl improvement collected in previous iterations
 *  @return a logl delta that is larger than 95% of the collected values
 */
double estBestLoglImp(SearchInfo* searchinfo);

/**
 *  Evaluate NNI moves for the current internal branch
 *  @param tr the current tree data structure
 *  @param pr partition data structure
 *  @param p the node representing the current branch
 *  @return number of positive NNIs found
 */
int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo);

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
double pllPerturbTree(pllInstance *tr, partitionList *pr, vector<pllNNIMove> &nnis);

/**
 * 	do 1 round of fastNNI
 *  return new tree log-likelihood if found improving NNI otherwise -1.0
 *
 *  @param[in] tr the tree data structure
 *  @param[in] pr partition data structure
 *  @param[out] nniList list containing information about the 2(n-3) evaluated NNIs
 *  @param[in/out] tabuNNIs tabu list
 *  @param[out] nni_count number of NNI that have been applied
 *  @param[out] deltaNNI average improvement made by one NNI
 */
double pllDoNNISearch(pllInstance* tr, partitionList *pr, SearchInfo &searchinfo);

void pllUpdateTabuList(pllInstance *tr, SearchInfo &searchinfo);

void pllSaveQuartetForSubTree(pllInstance* tr, nodeptr p, SearchInfo &searchinfo);


/**
 *  @brief Do 1 NNI move.
 *  @param[in] tr: the tree data structure
 *  @param[in] pr partition data structure
 *  @param[in] swap: represents one of the 2 NNI moves. Could be either 0 or 1
 *  @param[in] NNI_Type
 */
double doOneNNI(pllInstance * tr, partitionList *pr, nodeptr p, int swap, NNI_Type nni_type, SearchInfo *searchinfo = NULL);

void pllGetAllInBran(pllInstance *tr, vector<nodeptr> &branlist);

void pllGetAllInBranForSubtree(pllInstance *tr, nodeptr p, vector<nodeptr> &branlist);


string convertQuartet2String(nodeptr p);
/**
 *  Go through all 2(n-3) internal branches of the tree and
 *  evaluate all possible NNI moves
 */
void evalAllNNI(pllInstance* tr);

/**
 * 	@brief evaluate all NNIs within the subtree specified by node p
 * 	populates the list containing all possible NNI moves
 *
 * 	@param[in] tr: the tree data structure
 * 	@param[in] pr partition data structure
 * 	@param[in] p node pointer that specify the subtree
 */
void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo);

/*
 *  @brief return the array which can be used to store evaluated NNI moves
 *
 *  @param[in] tr: the tree data structure
 */
pllNNIMove *getNNIList(pllInstance* tr);



/*
 * ****************************************************************************
 * pllUFBoot area
 * ****************************************************************************
 */

/**
 * DTH:
 * pllUFBootData struct
 * This one keeps all info necessary to run UFBoot in PLL mode
 */
typedef struct{
    int max_candidate_trees;
    int treels_size;
    int save_all_trees;
//    pllBoolean save_all_br_lens;
    double logl_cutoff;
    int duplication_counter;
    int n_patterns;
    struct pllHashTable * treels;
    unsigned int candidate_trees_count; /* counter of trees in pllHashTable */
    double * treels_logl; // maintain size == treels_size
//    char ** treels_newick; // maintain size == treels_size
    double ** treels_ptnlh; // maintain size == treels_size
    int ** boot_samples;
    double * boot_logl;
    int * boot_counts;
    StrVector boot_trees;
} pllUFBootData;

/**
 * DTH:
 * The PLL version of saveCurrentTree function
 * @param tr: the tree (a pointer to a pllInstance)
 * @param pr: pointer to a partitionList (this one keeps tons of tree info)
 * @param p: root?
 */
void pllSaveCurrentTree(pllInstance* tr, partitionList *pr, nodeptr p);

/**
 * DTH:
 * Extract the array of site log likelihood to be kept in ptnlh
 * And update *cur_log
 * @param tr: the tree (pointer to an pllInstance)
 * @param ptnlh: to-be-kept array of site log likelihood
 * @param cur_logl: pointer to current tree log likelihood
 */
void pllComputePatternLikelihood(pllInstance* tr, double * ptnlh, double * cur_logl);

/**
 * DTH:
 * Resize some of the arrays in UFBootData if they're full
 * Along with update treels_size (to track the size of these arrays)
 */
void pllResizeUFBootData();

/**
 * DTH:
 * (Based on function Tree2StringREC of PLL)
 * Print out the tree topology with IQTree taxa ID (starts at 0) instead of PLL taxa ID (starts at 1)
 * @param All are the same as in PLL's
 */
char *pllTree2StringREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean printBranchLengths, pllBoolean printNames,
		pllBoolean printLikelihood, pllBoolean rellTree, pllBoolean finalPrint, int perGene, pllBoolean branchLabelSupport, pllBoolean printSHSupport);

#endif

