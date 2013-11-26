#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pll/pll.h"
#define GLOBAL_VARIABLES_DEFINITION
#include "nnisearch.h"

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

double doNNISearch(pllInstance* tr, partitionList *pr, SearchInfo &searchinfo) {
	double initLH = tr->likelihood;
	double finalLH = initLH;
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	/* data structure to store the initial tree topology + branch length */
	topol* curTree = _setupTopol(tr);
	saveTree(tr, curTree, numBranches);

	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, searchinfo);
		q = q->next;
		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr,  q->back, PLL_TRUE);
		} else {
			pllNewviewGeneric(tr, pr,  q->back, PLL_FALSE);
		}
	}

	if (searchinfo.curNumPosNNIs != 0) {
		/* Sort the NNI list ascendingly according to the log-likelihood */
		sort(searchinfo.nniList.begin(), searchinfo.nniList.end(), comparePLLNNIMove);
		/* list of independent positive NNI */
		vector<pllNNIMove> selectedNNIs;
		/* The best NNI is the first to come to the list */
		selectedNNIs.push_back(searchinfo.nniList.back());
		/* Subsequently add positive NNIs that are non-conflicting with the previous ones */
		int totalNNIs = searchinfo.nniList.size();
		int numPosNNI = searchinfo.curNumPosNNIs;
		for (int k = totalNNIs - 2; k > totalNNIs - numPosNNI - 1; k--) {
			int conflict = PLL_FALSE;
			/* Go through all the existing non-conflicting NNIs to check whether the next NNI will conflict with one of them */
			int j;
			for (j = 0; j < selectedNNIs.size(); j++) {
				if (searchinfo.nniList[k].p->number == selectedNNIs[j].p->number
						|| searchinfo.nniList[k].p->number == selectedNNIs[j].p->back->number) {
					conflict = PLL_TRUE;
					break;
				}
				if (searchinfo.nniList[k].p->back->number == selectedNNIs[j].p->number
						|| searchinfo.nniList[k].p->back->number
								== selectedNNIs[j].p->back->number) {
					conflict = PLL_TRUE;
					break;
				}
			}
			if (conflict) {
				continue;
			} else {
				selectedNNIs.push_back(searchinfo.nniList[k]);
			}
		}

		/* Applying all independent NNI moves */
		searchinfo.curNumAppliedNNIs = selectedNNIs.size();
		int MAXROLLBACK = 2;
		int step;
		for (step = 1; step <= MAXROLLBACK; step++) {
			int i;
			for (i = 0; i < searchinfo.curNumAppliedNNIs; i++) {
				/* do the topological change */
				doOneNNI(tr, pr, selectedNNIs[i].p, selectedNNIs[i].nniType, TOPO_ONLY);
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
					selectedNNIs[i].p->z[j] = selectedNNIs[i].z0[j];
					selectedNNIs[i].p->back->z[j] = selectedNNIs[i].z0[j];
					selectedNNIs[i].p->next->z[j] = selectedNNIs[i].z1[j];
					selectedNNIs[i].p->next->back->z[j] = selectedNNIs[i].z1[j];
					selectedNNIs[i].p->next->next->z[j] = selectedNNIs[i].z2[j];
					selectedNNIs[i].p->next->next->back->z[j] = selectedNNIs[i].z2[j];
					selectedNNIs[i].p->back->next->z[j] = selectedNNIs[i].z3[j];
					selectedNNIs[i].p->back->next->back->z[j] = selectedNNIs[i].z3[j];
					selectedNNIs[i].p->back->next->next->z[j] = selectedNNIs[i].z4[j];
					selectedNNIs[i].p->back->next->next->back->z[j] = selectedNNIs[i].z4[j];
				}
				/* update partial likelihood */
				if (numBranches > 1 && !tr->useRecom) {
					pllNewviewGeneric(tr, pr, selectedNNIs[i].p, PLL_TRUE);
					pllNewviewGeneric(tr, pr, selectedNNIs[i].p->back, PLL_TRUE);
				} else {
					pllNewviewGeneric(tr, pr, selectedNNIs[i].p, PLL_FALSE);
					pllNewviewGeneric(tr, pr, selectedNNIs[i].p->back, PLL_FALSE);
				}
			}
			pllTreeEvaluate(tr, pr, 1);
			/* new tree likelihood should not be smaller the likelihood of the computed best NNI */
			if (tr->likelihood < selectedNNIs[0].likelihood) {
				if (searchinfo.curNumAppliedNNIs == 1) {
					printf(
							"ERROR: new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
							tr->likelihood, selectedNNIs[0].likelihood);
					exit(1);
				}
				if (!restoreTree(curTree, tr, pr)) {
					printf("ERROR: failed to roll back tree \n");
					exit(1);
				}
				/* Only apply the best NNI after the tree has been rolled back */
				searchinfo.curNumAppliedNNIs  = 1;
			} else {
				if (tr->likelihood - initLH < 0.1) {
					if (!restoreTree(curTree, tr, pr)) {
						printf("ERROR: failed to roll back tree \n");
						exit(1);
					}
					searchinfo.curNumAppliedNNIs = 0;
				}
				break;
			}
		}
		finalLH = tr->likelihood;
	} else {
		searchinfo.curNumAppliedNNIs = 0;
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
	assert(swap == 1 || swap == 2);
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

string convertQuartet2String(nodeptr p) {
	nodeptr q = p->back;
	int pNr = p->number;
	int qNr = q->number;
	int pNei1Nr = p->next->back->number;
	int pNei2Nr = p->next->next->back->number;
	int qNei1Nr = q->next->back->number;
	int qNei2Nr = q->next->next->back->number;
	stringstream middle;
	stringstream left;
	stringstream right;
	stringstream res;
	if (pNr < qNr) {
		middle << pNr << "-" << qNr;
	} else {
		middle << qNr << "-" << pNr;
	}
	if (pNei1Nr < pNei2Nr) {
		left << pNei1Nr << "-" << pNei2Nr;
	} else {
		left << pNei2Nr << "-" << pNei1Nr;
	}
	if (qNei1Nr < qNei2Nr) {
		left << qNei1Nr << "-" << qNei2Nr;
	} else {
		left << qNei2Nr << "-" << qNei1Nr;
	}
	res << left << middle << right;
	return res.str();
}

int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	nodeptr q = p->back;
	assert(!isTip(p->number, tr->mxtips));
	assert(!isTip(q->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	int numPosNNI = 0;
	int i;
	pllNNIMove nni0; // dummy NNI to store backup information
	nni0.p = p;
	nni0.nniType = 0;
	nni0.likelihood = searchinfo.curLogl;
	for (i = 0; i < numBranches; i++) {
		nni0.z0[i] = p->z[i];
		nni0.z1[i] = p->next->z[i];
		nni0.z2[i] = p->next->next->z[i];
		nni0.z3[i] = q->next->z[i];
		nni0.z4[i] = q->next->next->z[i];
	}

	/* do an NNI move of type 1 */
	double lh1 = doOneNNI(tr, pr, p, 1, searchinfo.evalType);
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

	searchinfo.nniList.push_back(nni1);
	//nniList[*numBran] = nni1;

	if (nni1.likelihood > searchinfo.curLogl) {
		numPosNNI++;
		searchinfo.curNumPosNNIs++;
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
	double lh2 = doOneNNI(tr, pr, p, 2, searchinfo.evalType);
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

	//nniList[*numBran + 1] = nni2;
	searchinfo.nniList.push_back(nni2);
	if (nni2.likelihood > searchinfo.curLogl) {
		numPosNNI++;
		searchinfo.curNumPosNNIs++;
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
	return numPosNNI;
}

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		evalNNIForBran(tr, pr, p, searchinfo);
		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
		nodeptr q = p->next;
		while (q != p) {
			if (numBranches > 1 && !tr->useRecom) {
				pllNewviewGeneric(tr, pr, p->back, PLL_TRUE);
			} else {
				pllNewviewGeneric(tr, pr, p->back, PLL_FALSE);
			}
			evalNNIForSubtree(tr, pr, q->back, searchinfo);
			if (numBranches > 1 && !tr->useRecom) {
				pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			} else {
				pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			}
			q = q->next;
		}

	}
}
