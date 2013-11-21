#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pll/pll.h"
#define GLOBAL_VARIABLES_DEFINITION
#include "nnisearch.h"
#include "strmap.h"

#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__
#include <sys/resource.h>
#endif
/**
 verbose mode, determine how verbose should the screen be printed.
 */
enum VerboseMode {
	VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG
};

/**
 verbose level on the screen
 */
extern int verbose_mode;

/* program options */
double TOL_LIKELIHOOD_PHYLOLIB;
int numSmoothTree;
int nni0;
int nni5;
/* program options */

//pllNNIMove* nniList;



static int cmp_nni(const void* nni1, const void* nni2) {
	pllNNIMove* myNNI1 = (pllNNIMove*) nni1;
	pllNNIMove* myNNI2 = (pllNNIMove*) nni2;
	if (myNNI1->likelihood > myNNI2->likelihood)
		return 1;
	else if (myNNI1->likelihood < myNNI2->likelihood)
		return -1;
	else
		return 0;
}

int compareDouble(const void * a, const void * b) {
	if (*(double*) a > *(double*) b)
		return 1;
	else if (*(double*) a < *(double*) b)
		return -1;
	else
		return 0;
}

pllNNIMove *getNNIList(pllInstance* tr) {
	static pllNNIMove* nniList;
	if (nniList == NULL) {
		nniList = (pllNNIMove*) malloc(2 * (tr->mxtips - 3) * sizeof(pllNNIMove));
		assert(nniList != NULL);
	}
	return nniList;
}

pllNNIMove *getNonConflictNNIList(pllInstance* tr) {
	static pllNNIMove* nonConfNNIList;
	if (nonConfNNIList == NULL) {
		nonConfNNIList = (pllNNIMove*) malloc((tr->mxtips - 3) * sizeof(pllNNIMove));
		assert(nonConfNNIList != NULL);
	}
	return nonConfNNIList;
}

double perturbTree(pllInstance *tr, partitionList *pr, pllNNIMove *nnis, int numNNI) {
	int numBranches = pr->numberOfPartitions;
	int i;
	//printf("Perturbing %d NNIs \n", numNNI);
	for (i = 0; i < numNNI; i++) {
		/* First, do the topological change */
		//printf("Do pertubing NNI (%d - %d) with logl = %10.4f \n", nnis[i].p->number, nnis[i].p->back->number, nnis[i].likelihood);
		doOneNNI(tr, pr, nnis[i].p, nnis[i].nniType, TOPO_ONLY);
		/* Then apply the new branch lengths */
		int j;
		for (j = 0; j < numBranches; j++) {
			nnis[i].p->z[j] = nnis[i].z0[j];
			nnis[i].p->back->z[j] = nnis[i].z0[j];
			nnis[i].p->next->z[j] = nnis[i].z1[j];
			nnis[i].p->next->back->z[j] = nnis[i].z1[j];
			nnis[i].p->next->next->z[j] = nnis[i].z2[j];
			nnis[i].p->next->next->back->z[j] = nnis[i].z2[j];
			nnis[i].p->back->next->z[j] = nnis[i].z3[j];
			nnis[i].p->back->next->back->z[j] = nnis[i].z3[j];
			nnis[i].p->back->next->next->z[j] = nnis[i].z4[j];
			nnis[i].p->back->next->next->back->z[j] = nnis[i].z4[j];
		}
	}
	pllEvaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	pllTreeEvaluate(tr, pr, 1);
	return tr->likelihood;
}

void quicksort_nni(pllNNIMove* arr,int left, int right) {
    int i = left, j = right;
    pllNNIMove tmp, pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i].likelihood < pivot.likelihood)
            i++;
        while (pivot.likelihood < arr[j].likelihood)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quicksort_nni(arr, left, j);
    if (i < right)
        quicksort_nni(arr, i, right);
}

//TODO: Workaround for memory leak problem when calling setupTopol within doNNISearch
topol *_setupTopol(pllInstance* tr) {
	static topol* tree;
	if (tree == NULL)
		tree = setupTopol(tr->mxtips);
	return tree;
}

double doNNISearch(pllInstance* tr, partitionList *pr, int searchType, pllNNIMove *outNNIList, int* nni_count, double* deltaNNI) {
	double initLH = tr->likelihood;
	double finalLH = initLH;
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	/* data structure to store the initial tree topology + branch length */
	topol* curTree = _setupTopol(tr);
	saveTree(tr, curTree, numBranches);

	/* Initialize the NNI list that holds information about 2n-6 NNI moves */
	pllNNIMove nniList[2 * (tr->mxtips - 3)];

	/* Now fill up the NNI list */
	int numBran = 0; // number of visited internal branches during NNI evaluation
	int numPosNNI = 0; // number of positive NNI branhces found
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, nniList, searchType, &numBran, &numPosNNI, initLH);
		q = q->next;
		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr,  q->back, PLL_TRUE);
		} else {
			pllNewviewGeneric(tr, pr,  q->back, PLL_FALSE);
		}
	}

	int totalNNIs = numBran;
	/* Make sure all 2n-6 NNIs were evalauted */
	assert(totalNNIs == (2 * (tr->mxtips - 3)));
	/* Sort the NNI list ascendingly according to the log-likelihood */
	for (int i = 0; i < totalNNIs; i++) {
		outNNIList[i] = nniList[i];
	}
	if (numPosNNI == 0) {
		*nni_count = numPosNNI;
	} else  {
		qsort(nniList, totalNNIs, sizeof(pllNNIMove), cmp_nni);
		/* Generate a list of independent positive NNI */
		pllNNIMove nnis[tr->mxtips - 3];
		/* The best NNI is the first to come to the list */
		nnis[0] = nniList[totalNNIs - 1];
		/* Subsequently add positive NNIs that are non-conflicting with the previous ones */
		int numInNNI = 1; // size of the existing non-conflicting NNIs in the list;
		int k;
		for (k = totalNNIs - 2; k > totalNNIs - numPosNNI - 1; k--) {
			int conflict = PLL_FALSE;
			/* Go through all the existing non-conflicting NNIs to check whether the next NNI will conflict with one of them */
			int j;
			for (j = 0; j < numInNNI; j++) {
				if (nniList[k].p->number == nnis[j].p->number
						|| nniList[k].p->number == nnis[j].p->back->number) {
					conflict = PLL_TRUE;
					break;
				}
				if (nniList[k].p->back->number == nnis[j].p->number
						|| nniList[k].p->back->number
								== nnis[j].p->back->number) {
					conflict = PLL_TRUE;
					break;
				}
			}
			if (conflict) {
				continue;
			} else {
				nnis[numInNNI] = nniList[k];
				numInNNI++;
			}
		}

		/* Applying all independent NNI moves */
		int numNNI = numInNNI;
		int MAXROLLBACK = 2;
		int step;
		for (step = 1; step <= MAXROLLBACK; step++) {
			int i;
			for (i = 0; i < numNNI; i++) {
				/* do the topological change */
				doOneNNI(tr, pr, nnis[i].p, nnis[i].nniType, TOPO_ONLY);

//				printf(" Do NNI on branch: (%d-%d-%d) <-> (%d-%d-%d) / logl :%10.6f / length:%10.6f \n",
//						nnis[i].p->next->back->number, nnis[i].p->number,
//						nnis[i].p->next->next->back->number,
//						nnis[i].p->back->next->back->number,
//						nnis[i].p->back->number,
//						nnis[i].p->back->next->next->back->number,
//						nnis[i].likelihood, nnis[i].z0[0]);
				/*  apply branch lengths */
				int j;
				for (j = 0; j < numBranches; j++) {
					nnis[i].p->z[j] = nnis[i].z0[j];
					nnis[i].p->back->z[j] = nnis[i].z0[j];
					nnis[i].p->next->z[j] = nnis[i].z1[j];
					nnis[i].p->next->back->z[j] = nnis[i].z1[j];
					nnis[i].p->next->next->z[j] = nnis[i].z2[j];
					nnis[i].p->next->next->back->z[j] = nnis[i].z2[j];
					nnis[i].p->back->next->z[j] = nnis[i].z3[j];
					nnis[i].p->back->next->back->z[j] = nnis[i].z3[j];
					nnis[i].p->back->next->next->z[j] = nnis[i].z4[j];
					nnis[i].p->back->next->next->back->z[j] = nnis[i].z4[j];
				}
				/* update partial likelihood */
				if (numBranches > 1 && !tr->useRecom) {
					pllNewviewGeneric(tr, pr, nnis[i].p, PLL_TRUE);
					pllNewviewGeneric(tr, pr, nnis[i].p->back, PLL_TRUE);
				} else {
					pllNewviewGeneric(tr, pr, nnis[i].p, PLL_FALSE);
					pllNewviewGeneric(tr, pr, nnis[i].p->back, PLL_FALSE);
				}
			}
			pllTreeEvaluate(tr, pr, 1);
			/* new tree likelihood should not be smaller the likelihood of the computed best NNI */
			if (tr->likelihood < nnis[0].likelihood) {
				if (numNNI == 1) {
					printf(
							"ERROR: new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
							tr->likelihood, nnis[0].likelihood);
					exit(1);
				}
				if (!restoreTree(curTree, tr, pr)) {
					printf("ERROR: failed to roll back tree \n");
					exit(1);
				}
				/* Only apply the best NNI after the tree has been rolled back */
				numNNI = 1;
			} else {
				if (tr->likelihood - initLH < 0.01) {
					if (!restoreTree(curTree, tr, pr)) {
						printf("ERROR: failed to roll back tree \n");
						exit(1);
					}
					for (int i = 0; i < totalNNIs; i++) {
						outNNIList[i] = nniList[i];
					}
					numNNI = 0;
				}
				break;
			}
		}
		*nni_count = numNNI;
		*deltaNNI = (tr->likelihood - initLH) / numNNI;
		finalLH = tr->likelihood;
		//printf("numNNI = %d / logl = %10.4f \n", numNNI, tr->likelihood);
	}
	return finalLH;
}

/** @brief Optimize the length of a specific branch with variable number of Newton Raphson iterations

 Optimize the length of the branch connecting \a p and \a p->back
 for each partition (\a tr->numBranches) in library instance \a tr.

 @param tr
 The library instance

 @param pr
 Partition list

 @param p
 Endpoints of branch to be optimized

 @param newzpercycle Maximal number of Newton Raphson iterations

 */
void _update(pllInstance *tr, partitionList *pr, nodeptr p) {
	nodeptr q;
	int i;
	double z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	q = p->back;

	for (i = 0; i < numBranches; i++)
		z0[i] = q->z[i];

	if (numBranches > 1)
		makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_TRUE);
	else
		makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_FALSE);

	for (i = 0; i < numBranches; i++) {
		if (!tr->partitionConverged[i]) {
			if (PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ) {
				tr->partitionSmoothed[i] = PLL_FALSE;
			}

			p->z[i] = q->z[i] = z[i];
		}
	}
}

double doOneNNI(pllInstance *tr, partitionList *pr, nodeptr p, int swap, int evalType) {
	nodeptr q;
	nodeptr tmp;
	q = p->back;
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	}

	if (evalType == TOPO_ONLY) {
		return 0.0;
	} else if (evalType == ONE_BRAN_OPT) {
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);

		} else {
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		_update(tr, pr, p);
		pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
	} else if (evalType == NO_BRAN_OPT) {
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);

		} else {
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
	} else { // 5 branches optimization
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);
		} else {
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		nodeptr r; // temporary node poiter
		r = p->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		r = p->next->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
		_update(tr, pr, p);
		// optimize 2 branches at node q
		r = q->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		r = q->next->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		pllEvaluateGeneric(tr, pr, r, PLL_FALSE, PLL_FALSE);
	}
	return tr->likelihood;

}

int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, pllNNIMove* nniList, int searchType, int* numBran, double curLH) {
	nodeptr q = p->back;
	assert(!isTip(p->number, tr->mxtips));
	assert(!isTip(q->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	int betterNNI = 0;
	int i;
	pllNNIMove nni0; // dummy NNI to store backup information
	nni0.p = p;
	nni0.nniType = 0;
	nni0.likelihood = curLH;
	for (i = 0; i < numBranches; i++) {
		nni0.z0[i] = p->z[i];
		nni0.z1[i] = p->next->z[i];
		nni0.z2[i] = p->next->next->z[i];
		nni0.z3[i] = q->next->z[i];
		nni0.z4[i] = q->next->next->z[i];
	}

	/* do an NNI move of type 1 */
	double lh1 = doOneNNI(tr, pr, p, 1, searchType);
	pllNNIMove nni1;
	nni1.p = p;
	nni1.nniType = 1;
	// Store the optimized branch lengths
	for (i = 0; i < numBranches; i++) {
		nni1.z0[i] = p->z[i];
		nni1.z1[i] = p->next->z[i];
		nni1.z2[i] = p->next->next->z[i];
		nni1.z3[i] = q->next->z[i];
		nni1.z4[i] = q->next->next->z[i];
	}
	nni1.likelihood = lh1;

	nniList[*numBran] = nni1;

	if (nni1.likelihood > curLH) {
		betterNNI = 1;
	}

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 1, TOPO_ONLY);
	/* Restore the old branch length */
	for (i = 0; i < numBranches; i++) {
		p->z[i] = nni0.z0[i];
		q->z[i] = nni0.z0[i];
		p->next->z[i] = nni0.z1[i];
		p->next->back->z[i] = nni0.z1[i];
		p->next->next->z[i] = nni0.z2[i];
		p->next->next->back->z[i] = nni0.z2[i];
		q->next->z[i] = nni0.z3[i];
		q->next->back->z[i] = nni0.z3[i];
		q->next->next->z[i] = nni0.z4[i];
		q->next->next->back->z[i] = nni0.z4[i];
	}

	/* do an NNI move of type 2 */
	double lh2 = doOneNNI(tr, pr, p, 2, searchType);
	// Create the nniMove struct to store this move
	pllNNIMove nni2;
	nni2.p = p;
	nni2.nniType = 2;
	// Store the optimized and unoptimized central branch length
	for (i = 0; i < numBranches; i++) {
		nni2.z0[i] = p->z[i];
		nni2.z1[i] = p->next->z[i];
		nni2.z2[i] = p->next->next->z[i];
		nni2.z3[i] = q->next->z[i];
		nni2.z4[i] = q->next->next->z[i];
	}
	nni2.likelihood = lh2;

	nniList[*numBran + 1] = nni2;

	if (nni2.likelihood > curLH) {
		betterNNI = 1;
	}

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 2, TOPO_ONLY);
    /* Restore the old branch length */
    for (i = 0; i < numBranches; i++) {
      p->z[i] = nni0.z0[i];
      q->z[i] = nni0.z0[i];
      p->next->z[i] = nni0.z1[i];
      p->next->back->z[i] = nni0.z1[i];
      p->next->next->z[i] = nni0.z2[i];
      p->next->next->back->z[i] = nni0.z2[i];
      q->next->z[i] = nni0.z3[i];
      q->next->back->z[i] = nni0.z3[i];
      q->next->next->z[i] = nni0.z4[i];
      q->next->next->back->z[i] = nni0.z4[i];
  }
	return betterNNI;
}

/*int test_nni(int argc, char * argv[])
 {
 struct rusage usage;
 struct timeval start, end;
 struct timeval start_tmp, end_tmp;
 getrusage(RUSAGE_SELF, &usage);
 start = usage.ru_utime;
 tree        * tr;
 if (argc < 2)
 {
 fprintf (stderr, "syntax: %s [binary-alignment-file] [prefix] \n", argv[0]);
 return (1);
 }
 tr = (tree *)malloc(sizeof(tree));

 read the binary input, setup tree, initialize model with alignment
 read_msa(tr,argv[1]);

 tr->randomNumberSeed = 665;

 Create random tree
 //makeRandomTree(tr);
 //printf("RANDOM TREE: Number of taxa: %d\n", tr->mxtips);
 //printf("RANDOM TREE: Number of partitions: %d\n", tr->NumberOfModels);

 printf("Creating parsimony tree ...\n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 makeParsimonyTree(tr);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Parsimony tree created in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);

 //printf("PARSIMONY TREE: Number of taxa: %d\n", tr->mxtips);
 //printf("PARSIMONY TREE: Number of partitions: %d\n", tr->NumberOfModels);

 compute the LH of the full tree
 //printf ("Virtual root: %d\n", tr->start->number);

 int printBranchLengths=TRUE;
 Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, TRUE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 //fprintf(stderr, "%s\n", tr->tree_string);
 char parTree[100];
 strcpy (parTree, argv[1]);
 strcat(parTree,".parsimony_tree");
 FILE *parTreeFile = fopen(parTree, "w");
 fprintf(parTreeFile, "%s\n", tr->tree_string);
 fclose(parTreeFile);

 Model optimization
 printf("Optimizing model parameters and branch lengths ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 modOpt(tr, 0.1);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Model parameters and branch lengths optimized in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);
 evaluateGeneric(tr, tr->start, FALSE);
 printf("Likelihood of the parsimony tree: %f\n", tr->likelihood);

 8 rounds of branch length optimization
 //smoothTree(tr, 32);
 //evaluateGeneric(tr, tr->start, TRUE);
 //printf("Likelihood after branch length optimization: %f\n", tr->likelihood);
 //printBranchLengths=TRUE;
 //Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 //fprintf(stderr, "%s\n", tr->tree_string);


 int foundNNI=TRUE;
 int nniRound=1;
 printf("Starting the first NNI search ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 double curLH = tr->likelihood;
 while (TRUE) {
 double newLH = doNNISearch(tr);
 if (newLH == 0.0) {
 break;
 } else {
 printf("Likelihood after doing NNI round %d = %f\n", nniRound, newLH);
 //if (newLH - curLH < 0.001)
 //break;
 curLH = newLH;
 nniRound++;
 }

 }
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("First NNI search finished in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);

 printf("Optimizing model parameters and branch lengths ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 modOpt(tr, 0.1);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Model parameters and branch lengths optimized in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);
 evaluateGeneric(tr, tr->start, FALSE);
 printf("Likelihood after model optimization: %f\n", tr->likelihood);
 Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, TRUE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 char resultFile[100];
 strcpy (resultFile, argv[1]);
 strcat(resultFile,".result");
 FILE *treefile = fopen(resultFile, "w");
 fprintf(treefile, "%s\n", tr->tree_string);
 getrusage(RUSAGE_SELF, &usage);
 end = usage.ru_utime;
 printf("CPU time = %f seconds \n", (double) end.tv_sec + end.tv_usec/1.0e6 - start.tv_sec - start.tv_usec/1.0e6);
 printf("Result tree written to %s \n", resultFile);
 //printTopology(tr, TRUE);
 fclose(treefile);
 return (0);
 }*/

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, pllNNIMove* nniList, int searchType, int* numBran, int* numPosNNI,
		double curLH) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		int betterNNI = evalNNIForBran(tr, pr, p, nniList, searchType, numBran, curLH);
		if (betterNNI) {
			*numPosNNI = *numPosNNI + 1;
		}
		*numBran = *numBran + 2;
		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

		nodeptr q = p->next;
		while (q != p) {
			if (numBranches > 1 && !tr->useRecom) {
				pllNewviewGeneric(tr, pr, p->back, PLL_TRUE);
			} else {
				pllNewviewGeneric(tr, pr, p->back, PLL_FALSE);
			}
			evalNNIForSubtree(tr, pr, q->back, nniList, searchType, numBran, numPosNNI, curLH);
			if (numBranches > 1 && !tr->useRecom) {
				pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			} else {
				pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			}
			q = q->next;
		}

	}
}
