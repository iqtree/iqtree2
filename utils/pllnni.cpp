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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <assert.h>

#define GLOBAL_VARIABLES_DEFINITION

#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
#include <sys/resource.h>
#endif

#include "tree/phylotree.h"
#include "pllnni.h"
#include "alignment/alignment.h"

/* program options */
int nni5;
extern VerboseMode verbose_mode;
int NNI_MAX_NR_STEP = 10;

/* program options */
extern Params *globalParams;
extern Alignment *globalAlignment;

/**
 * map from newick tree string to frequencies that a tree is revisited during tree search
 */
StringIntMap pllTreeCounter;


/*
 * ****************************************************************************
 * pllUFBoot area
 * ****************************************************************************
 */

pllUFBootData * pllUFBootDataPtr = NULL;


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
		ASSERT(nniList != NULL);
	}
	return nniList;
}

pllNNIMove *getNonConflictNNIList(pllInstance* tr) {
	static pllNNIMove* nonConfNNIList;
	if (nonConfNNIList == NULL) {
		nonConfNNIList = (pllNNIMove*) malloc((tr->mxtips - 3) * sizeof(pllNNIMove));
		ASSERT(nonConfNNIList != NULL);
	}
	return nonConfNNIList;
}

double pllDoRandomNNIs(pllInstance *tr, partitionList *pr, int numNNI) {
	int numInBrans = tr->mxtips - 3;
	int numNNIinStep = (int) numInBrans / 5;

	// devided in multiple round, each round collect 1/5 of numNNI
	int cnt1 = 0;
	unordered_set<int> selectedNodes;
	vector<nodeptr> selectedBrans;
	vector<nodeptr> brans;
	do {
		int cnt2 = 0;
		selectedNodes.clear();
		selectedBrans.clear();
		brans.clear();
		pllGetAllInBran(tr, brans);
		ASSERT(brans.size() == numInBrans);
		while(cnt2 < numNNIinStep && cnt2 < numNNI) {
			int branIndex = random_int(numInBrans);
			if (selectedNodes.find(brans[branIndex]->number) == selectedNodes.end() &&
					selectedNodes.find(brans[branIndex]->back->number) == selectedNodes.end()) {
				selectedNodes.insert(brans[branIndex]->number);
				selectedNodes.insert(brans[branIndex]->back->number);
				selectedBrans.push_back(brans[branIndex]);
				cnt2++;
			}
		}
		for (vector<nodeptr>::iterator it = selectedBrans.begin(); it != selectedBrans.end(); ++it) {
			int nniType = random_int(2);
			doOneNNI(tr, pr, (*it), nniType, TOPO_ONLY);
		}
		cnt1 += selectedBrans.size();
		if (numNNI - cnt1 < numNNIinStep) {
			numNNIinStep = numNNI - cnt1;
		}
	} while (cnt1 < numNNI);
	pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	pllOptimizeBranchLengths(tr, pr, 1);
	return tr->likelihood;
}

void pllGetAllInBran(pllInstance *tr, vector<nodeptr> &branlist) {
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		pllGetAllInBranForSubtree(tr, q->back, branlist);
		q = q->next;
	}
}

void pllGetAllInBranForSubtree(pllInstance *tr, nodeptr p, vector<nodeptr> &branlist) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		branlist.push_back(p);
		nodeptr q = p->next;
		while (q != p) {
			pllGetAllInBranForSubtree(tr, q->back, branlist);
			q = q->next;
		}
	}
}

double pllPerturbTree(pllInstance *tr, partitionList *pr, vector<pllNNIMove> &nnis) {
	//printf("Perturbing %d NNIs \n", numNNI);
	for (vector<pllNNIMove>::iterator it = nnis.begin(); it != nnis.end(); it++) {
		printf("Do pertubing NNI (%d - %d) with logl = %10.4f \n", (*it).p->number, (*it).p->back->number, (*it).likelihood);
		doOneNNI(tr, pr, (*it).p, (*it).nniType, TOPO_ONLY);
		updateBranchLengthForNNI(tr, pr, (*it));

	}
	pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	pllOptimizeBranchLengths(tr, pr, 1);
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

vector<string> getAffectedBranches(pllInstance* tr, nodeptr p) {
	vector<string> res;
	res.push_back(getBranString(p));
	nodeptr q = p->back;
	nodeptr p_nei = p->next;
	nodeptr q_nei = q->next;
	while (p_nei != p) {
		res.push_back(getBranString(p_nei));
		if (!isTip(p_nei->back->number, tr->mxtips)) {
			res.push_back(getBranString(p_nei->back->next));
			res.push_back(getBranString(p_nei->back->next->next));
		}
		p_nei = p_nei->next;
	}
	while (q_nei != q) {
		res.push_back(getBranString(q_nei));
		if (!isTip(q_nei->back->number, tr->mxtips)) {
			res.push_back(getBranString(q_nei->back->next));
			res.push_back(getBranString(q_nei->back->next->next));
		}
		q_nei = q_nei->next;
	}
	return res;
}

string getBranString(nodeptr p) {
	stringstream branString;
	nodeptr q = p->back;
	if (p->number < q->number) {
		branString << p->number << "-" << q->number;
	} else {
		branString << q->number << "-" << p->number;
	}
	return branString.str();
}

set<int> getAffectedNodes(pllInstance* tr, nodeptr p) {
	set<int> nodeSet;
	nodeptr q = p->back;
	nodeptr p_nei = p->next;
	nodeptr q_nei = q->next;
	nodeSet.insert(p->number);
	nodeSet.insert(q->number);
	while (p_nei != p) {
		nodeptr nei = p_nei->back;
		if (isTip(nei->number, tr->mxtips)) {
			nodeSet.insert(nei->number);
		} else {
			nodeSet.insert(nei->number);
			nodeSet.insert(nei->next->back->number);
			nodeSet.insert(nei->next->next->back->number);
		}
		p_nei = p_nei->next;
	}
	while (q_nei != q) {
		nodeptr nei = q_nei->back;
		if (isTip(nei->number, tr->mxtips)) {
			nodeSet.insert(nei->number);
		} else {
			nodeSet.insert(nei->number);
			nodeSet.insert(nei->next->back->number);
			nodeSet.insert(nei->next->next->back->number);
		}
		q_nei = q_nei->next;
	}
	return nodeSet;
}

void pllEvalAllNNIs(pllInstance *tr, partitionList *pr, SearchInfo &searchinfo) {
    /* DTH: mimic IQTREE::optimizeNNI 's first call to IQTREE::saveCurrentTree */
    if((globalParams->online_bootstrap == PLL_TRUE) &&
            (globalParams->gbo_replicates > 0)){
        tr->fastScaling = PLL_FALSE;
        pllEvaluateLikelihood(tr, pr, tr->start, PLL_FALSE, PLL_TRUE);
        pllSaveCurrentTree(tr, pr, tr->start);
    }

	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, searchinfo);
		q = q->next;
	}
}

/*
void pllSaveQuartetForSubTree(pllInstance *tr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
			evalNNIForBran(tr, pr, p, searchinfo);
		nodeptr q = p->next;
		while (q != p) {
			pllSaveQuartetForSubTree(tr, q->back, searchinfo);
			q = q->next;
		}
	}
}

void pllSaveAllQuartet(pllInstance *tr, SearchInfo &searchinfo) {
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		pllSaveQuartetForSubTree(tr, q->back, searchinfo);
	}
}
*/

double pllDoNNISearch(pllInstance* tr, partitionList *pr, SearchInfo &searchinfo) {
	double initLH = tr->likelihood;
	double finalLH = initLH;
	vector<pllNNIMove> selectedNNIs;
	unordered_set<int> selectedNodes;
    int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	/* data structure to store the initial tree topology + branch length */
	topol* curTree = _setupTopol(tr);
	saveTree(tr, curTree, numBranches);

	/* evaluate NNIs */
	pllEvalAllNNIs(tr, pr, searchinfo);

	if (globalParams->speednni) {
		searchinfo.aBranches.clear();
	}

	/* apply non-conflicting positive NNIs */
	if (searchinfo.posNNIList.size() != 0) {
		sort(searchinfo.posNNIList.begin(), searchinfo.posNNIList.end(), comparePLLNNIMove);
        if (verbose_mode >= VB_DEBUG) {
        	cout << "curScore: "  << searchinfo.curLogl << endl;
            for (int i = 0; i < searchinfo.posNNIList.size(); i++) {
                cout << "Logl of positive NNI " << i << " : " << searchinfo.posNNIList[i].likelihood << endl;
            }
        }
		for (vector<pllNNIMove>::reverse_iterator rit = searchinfo.posNNIList.rbegin(); rit != searchinfo.posNNIList.rend(); ++rit) {
			if (selectedNodes.find((*rit).p->number) == selectedNodes.end() && selectedNodes.find((*rit).p->back->number) == selectedNodes.end()) {
				selectedNNIs.push_back((*rit));
				selectedNodes.insert((*rit).p->number);
				selectedNodes.insert((*rit).p->back->number);
			} else {
				continue;
			}
		}

		/* Applying all independent NNI moves */
		searchinfo.curNumAppliedNNIs = selectedNNIs.size();
		for (vector<pllNNIMove>::iterator it = selectedNNIs.begin(); it != selectedNNIs.end(); it++) {
			/* do the topological change */
			doOneNNI(tr, pr, (*it).p, (*it).nniType, TOPO_ONLY);
			if (globalParams->speednni) {
				vector<string> aBranches = getAffectedBranches(tr, (*it).p);
				searchinfo.aBranches.insert(aBranches.begin(), aBranches.end());
			}
			updateBranchLengthForNNI(tr, pr, (*it));
            if (numBranches > 1 && !tr->useRecom) {
                pllUpdatePartials(tr, pr, (*it).p, PLL_TRUE);
                pllUpdatePartials(tr, pr, (*it).p->back, PLL_TRUE);
            } else {
                pllUpdatePartials(tr, pr, (*it).p, PLL_FALSE);
                pllUpdatePartials(tr, pr, (*it).p->back, PLL_FALSE);
            }
		}
		if (selectedNNIs.size() != 0) {
			//pllEvaluateLikelihood(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
			pllOptimizeBranchLengths(tr, pr, 1);
			if (globalParams->count_trees) {
	            countDistinctTrees(tr, pr);
			}
			int numNNI = selectedNNIs.size();
			/* new tree likelihood should not be smaller the likelihood of the computed best NNI */
			while (tr->likelihood < selectedNNIs.back().likelihood) {
				if (numNNI == 1) {
					printf("ERROR: new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
							tr->likelihood, selectedNNIs[0].likelihood);
					ASSERT(0);
				} else {
					cout << "Best logl: " << selectedNNIs.back().likelihood << " / " << "NNI step " << searchinfo.curNumNNISteps<< " / Applying " << numNNI << " NNIs give logl: " << tr->likelihood << " (worse than best)";
					cout << " / Roll back tree ... " << endl;
			        //restoreTL(rl, tr, 0, pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
				    if (!restoreTree(curTree, tr, pr)) {
				        printf("ERROR: failed to roll back tree \n");
				        ASSERT(0);
				    }
				    // If tree log-likelihood decreases only apply the best NNI
					numNNI = 1;
					int count = numNNI;
					for (vector<pllNNIMove>::reverse_iterator rit = selectedNNIs.rbegin(); rit != selectedNNIs.rend(); ++rit) {
						doOneNNI(tr, pr, (*rit).p, (*rit).nniType, TOPO_ONLY);
						updateBranchLengthForNNI(tr, pr, (*rit));
			            if (numBranches > 1 && !tr->useRecom) {
			                pllUpdatePartials(tr, pr, (*rit).p, PLL_TRUE);
			                pllUpdatePartials(tr, pr, (*rit).p->back, PLL_TRUE);
			            } else {
			                pllUpdatePartials(tr, pr, (*rit).p, PLL_FALSE);
			                pllUpdatePartials(tr, pr, (*rit).p->back, PLL_FALSE);
			            }
						count--;
						if (count == 0) {
							break;
						}
					}
		            //pllEvaluateLikelihood(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
					pllOptimizeBranchLengths(tr, pr, 1);
					//cout << "Number of NNIs reduced to " << numNNI << ": " << tr->likelihood << endl;

					/* Only apply the best NNI after the tree has been rolled back */
					searchinfo.curNumAppliedNNIs = numNNI;
				}
			}
			if (tr->likelihood - initLH < 0.1) {
				searchinfo.curNumAppliedNNIs = 0;
			}
			finalLH = tr->likelihood;
		}
	} else {
		searchinfo.curNumAppliedNNIs = 0;
	}
	//freeTopol(curTree);
	return finalLH;
}


void updateBranchLengthForNNI(pllInstance* tr, partitionList *pr, pllNNIMove &nni) {
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	/*  apply branch lengths */
	for (int i = 0; i < numBranches; i++) {
		nni.p->z[i] = nni.z0[i];
		nni.p->back->z[i] = nni.z0[i];
		nni.p->next->z[i] = nni.z1[i];
		nni.p->next->back->z[i] = nni.z1[i];
		nni.p->next->next->z[i] = nni.z2[i];
		nni.p->next->next->back->z[i] = nni.z2[i];
		nni.p->back->next->z[i] = nni.z3[i];
		nni.p->back->next->back->z[i] = nni.z3[i];
		nni.p->back->next->next->z[i] = nni.z4[i];
		nni.p->back->next->next->back->z[i] = nni.z4[i];
	}
	/* update partial likelihood */
//	if (numBranches > 1 && !tr->useRecom) {
//		pllNewviewGeneric(tr, pr, nni.p, PLL_TRUE);
//		pllNewviewGeneric(tr, pr, nni.p->back, PLL_TRUE);
//	} else {
//		pllNewviewGeneric(tr, pr, nni.p, PLL_FALSE);
//		pllNewviewGeneric(tr, pr, nni.p->back, PLL_FALSE);
//	}
}

void pllOptimizeOneBranch(pllInstance *tr, partitionList *pr, nodeptr p) {
    nodeptr  q;
    int i;
    double   z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
    int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

    #ifdef _DEBUG_UPDATE
      double
        startLH;

      pllEvaluateLikelihood (tr, p);

      startLH = tr->likelihood;
    #endif

    q = p->back;

    for(i = 0; i < numBranches; i++)
      z0[i] = q->z[i];

    if(numBranches > 1)
      makenewzGeneric(tr, pr, p, q, z0, PLL_ITERATIONS, z, PLL_TRUE);
    else
      makenewzGeneric(tr, pr, p, q, z0, PLL_ITERATIONS, z, PLL_FALSE);

    for(i = 0; i < numBranches; i++)
    {
      if(!tr->partitionConverged[i])
      {
        if(PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ)
        {
          tr->partitionSmoothed[i] = PLL_FALSE;
        }

        p->z[i] = q->z[i] = z[i];
      }
    }

    #ifdef _DEBUG_UPDATE
      pllEvaluateLikelihood (tr, p);

      if(tr->likelihood <= startLH)
        {
          if(fabs(tr->likelihood - startLH) > 0.01)
      {
        printf("%f %f\n", startLH, tr->likelihood);
        ASSERT(0);
      }
        }
    #endif
}

double doOneNNI(pllInstance *tr, partitionList *pr, nodeptr p, int swap, NNI_Type nni_type, SearchInfo *searchinfo) {
	ASSERT(swap == 0 || swap == 1);
	nodeptr q;
	nodeptr tmp;
	q = p->back;
	ASSERT(!isTip(q->number, tr->mxtips));
	ASSERT(!isTip(p->number, tr->mxtips));
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

	if (nni_type == TOPO_ONLY) {
		return 0.0;
	}

    if (numBranches > 1 && !tr->useRecom) {
        pllUpdatePartials(tr, pr, p, PLL_TRUE);
        pllUpdatePartials(tr, pr, q, PLL_TRUE);
    } else {
        pllUpdatePartials(tr, pr, p, PLL_FALSE);
        pllUpdatePartials(tr, pr, q, PLL_FALSE);
    }
    // Optimize the central branch
    pllOptimizeOneBranch(tr, pr, p);
    if((globalParams->online_bootstrap == PLL_TRUE) && (globalParams->gbo_replicates > 0)){
        tr->fastScaling = PLL_FALSE;
        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_TRUE); // DTH: modified the last arg
        pllSaveCurrentTree(tr, pr, p);
    }else{
        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_FALSE);
    }
    // if after optimizing the central branch we already obtain better logl
    // then there is no need for optimizing other branches
    if (tr->likelihood > searchinfo->curLogl) {
        return tr->likelihood;
    }
    // Optimize 4 other branches
    if (nni_type == NNI5) {
        if (numBranches > 1 && !tr->useRecom) {
            pllUpdatePartials(tr, pr, q, PLL_TRUE);
        } else {
            pllUpdatePartials(tr, pr, q, PLL_FALSE);
        }
        nodeptr r;
        r = p->next;
        if (numBranches > 1 && !tr->useRecom) {
            pllUpdatePartials(tr, pr, r, PLL_TRUE);
        } else {
            pllUpdatePartials(tr, pr, r, PLL_FALSE);
        }
        pllOptimizeOneBranch(tr, pr, r);
//        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_FALSE);
//        if (tr->likelihood > searchinfo->curLogl) {
//            return tr->likelihood;
//        }
        r = p->next->next;
        if (numBranches > 1 && !tr->useRecom)
            pllUpdatePartials(tr, pr, r, PLL_TRUE);
        else
            pllUpdatePartials(tr, pr, r, PLL_FALSE);
        pllOptimizeOneBranch(tr, pr, r);
//        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_FALSE);
//        if (tr->likelihood > searchinfo->curLogl) {
//            return tr->likelihood;
//        }
        if (numBranches > 1 && !tr->useRecom)
            pllUpdatePartials(tr, pr, p, PLL_TRUE);
        else
            pllUpdatePartials(tr, pr, p, PLL_FALSE);
        pllOptimizeOneBranch(tr, pr, p);
//        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_FALSE);
//        if (tr->likelihood > searchinfo->curLogl) {
//            return tr->likelihood;
//        }
        // optimize 2 branches at node q
        r = q->next;
        if (numBranches > 1 && !tr->useRecom)
            pllUpdatePartials(tr, pr, r, PLL_TRUE);
        else
            pllUpdatePartials(tr, pr, r, PLL_FALSE);
        pllOptimizeOneBranch(tr, pr, r);
//        pllEvaluateLikelihood(tr, pr, p, PLL_FALSE, PLL_FALSE);
//        if (tr->likelihood > searchinfo->curLogl) {
//            return tr->likelihood;
//        }
        r = q->next->next;
        if (numBranches > 1 && !tr->useRecom)
            pllUpdatePartials(tr, pr, r, PLL_TRUE);
        else
            pllUpdatePartials(tr, pr, r, PLL_FALSE);
        pllOptimizeOneBranch(tr, pr, r);
        if((globalParams->online_bootstrap == PLL_TRUE) &&
                        (globalParams->gbo_replicates > 0)){
            tr->fastScaling = PLL_FALSE;
            pllEvaluateLikelihood(tr, pr, r, PLL_FALSE, PLL_TRUE); // DTH: modified the last arg
            pllSaveCurrentTree(tr, pr, r);
        }else{
            pllEvaluateLikelihood(tr, pr, r, PLL_FALSE, PLL_FALSE);
        }
    }
	return tr->likelihood;
}

double estBestLoglImp(SearchInfo* searchinfo) {
    double res = 0.0;
    int index = floor(searchinfo->deltaLogl.size() * 5 / 100);
    set<double>::reverse_iterator ri;
    for (ri = searchinfo->deltaLogl.rbegin(); ri != searchinfo->deltaLogl.rend(); ++ri) {
        //cout << (*ri) << " ";
        --index;
        if (index == 0) {
            res = (*ri);
            break;
        }
    }
    //cout << res << endl;
    //cout << searchinfo->deltaLogl.size() << endl;
    return res;
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
		middle << "-" << pNr << "-" << qNr << "-";
	} else {
		middle << "-" << qNr << "-" << pNr << "-";
	}
	if (pNei1Nr < pNei2Nr) {
		left << pNei1Nr << "-" << pNei2Nr;
	} else {
		left << pNei2Nr << "-" << pNei1Nr;
	}
	if (qNei1Nr < qNei2Nr) {
		right << qNei1Nr << "-" << qNei2Nr;
	} else {
		right << qNei2Nr << "-" << qNei1Nr;
	}
	res << left.str() << middle.str() << right.str();
	return res.str();
}

void countDistinctTrees(pllInstance* pllInst, partitionList *pllPartitions) {
    pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_FALSE,
            PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	PhyloTree mtree;
	mtree.rooted = false;
	mtree.aln = globalAlignment;
	mtree.readTreeString(string(pllInst->tree_string));
//    mtree.root = mtree.findNodeName(globalAlignment->getSeqName(0));
    mtree.setRootNode(mtree.params->root);
	ostringstream ostr;
	mtree.printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	string tree_str = ostr.str();
	if (pllTreeCounter.find(tree_str) == pllTreeCounter.end()) {
		// not found in hash_map
	    pllTreeCounter[tree_str] = 1;
	} else {
		// found in hash_map
	    pllTreeCounter[tree_str]++;
	}
}

int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	nodeptr q = p->back;
	ASSERT(!isTip(p->number, tr->mxtips));
	ASSERT(!isTip(q->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	int numPosNNI = 0;
	int i;
    pllNNIMove nniA; // dummy NNI to store backup information
    nniA.p = p;
    nniA.nniType = 0;
    nniA.likelihood = searchinfo.curLogl;
    for (i = 0; i < PLL_NUM_BRANCHES; i++) {
        nniA.z0[i] = p->z[i];
        nniA.z1[i] = p->next->z[i];
        nniA.z2[i] = p->next->next->z[i];
        nniA.z3[i] = q->next->z[i];
        nniA.z4[i] = q->next->next->z[i];
    }

	pllNNIMove bestNNI;

	/* do an NNI move of type 1 */
	double lh1 = doOneNNI(tr, pr, p, 0, searchinfo.nni_type, &searchinfo);
    if (globalParams->count_trees) {
        countDistinctTrees(tr, pr);
    }
    pllNNIMove nniB;
    nniB.p = p;
    nniB.nniType = 0;
    // Store the optimized branch lengths
    for (i = 0; i < PLL_NUM_BRANCHES; i++) {
        nniB.z0[i] = p->z[i];
        nniB.z1[i] = p->next->z[i];
        nniB.z2[i] = p->next->next->z[i];
        nniB.z3[i] = q->next->z[i];
        nniB.z4[i] = q->next->next->z[i];
    }
    nniB.likelihood = lh1;
    nniB.loglDelta = lh1 - nniB.likelihood;
    nniB.negLoglDelta = -nniB.loglDelta;

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 0, TOPO_ONLY);
	/* Restore the old branch length */
    for (i = 0; i < PLL_NUM_BRANCHES; i++) {
        p->z[i] = nniA.z0[i];
        q->z[i] = nniA.z0[i];
        p->next->z[i] = nniA.z1[i];
        p->next->back->z[i] = nniA.z1[i];
        p->next->next->z[i] = nniA.z2[i];
        p->next->next->back->z[i] = nniA.z2[i];
        q->next->z[i] = nniA.z3[i];
        q->next->back->z[i] = nniA.z3[i];
        q->next->next->z[i] = nniA.z4[i];
        q->next->next->back->z[i] = nniA.z4[i];
    }

	/* do an NNI move of type 2 */
	double lh2 = doOneNNI(tr, pr, p, 1, searchinfo.nni_type, &searchinfo);
	if (globalParams->count_trees)
		countDistinctTrees(tr, pr);

    // Create the nniMove struct to store this move
    pllNNIMove nniC;
    nniC.p = p;
    nniC.nniType = 1;
    // Store the optimized central branch length
    for (i = 0; i < PLL_NUM_BRANCHES; i++) {
        nniC.z0[i] = p->z[i];
        nniC.z1[i] = p->next->z[i];
        nniC.z2[i] = p->next->next->z[i];
        nniC.z3[i] = q->next->z[i];
        nniC.z4[i] = q->next->next->z[i];
    }
    nniC.likelihood = lh2;
    nniC.loglDelta = lh2 - nniA.likelihood;
    nniC.negLoglDelta = -nniC.loglDelta;
    
    if (nniC.likelihood > nniB.likelihood) {
        bestNNI = nniC;
    } else {
        bestNNI = nniB;
    }

	if (bestNNI.likelihood > searchinfo.curLogl + 1e-6) {
		numPosNNI++;
		searchinfo.posNNIList.push_back(bestNNI);
	}

    /* Restore previous NNI move */
    doOneNNI(tr, pr, p, 1, TOPO_ONLY);
    /* Restore the old branch length */
    for (i = 0; i < PLL_NUM_BRANCHES; i++) {
        p->z[i] = nniA.z0[i];
        q->z[i] = nniA.z0[i];
        p->next->z[i] = nniA.z1[i];
        p->next->back->z[i] = nniA.z1[i];
        p->next->next->z[i] = nniA.z2[i];
        p->next->next->back->z[i] = nniA.z2[i];
        q->next->z[i] = nniA.z3[i];
        q->next->back->z[i] = nniA.z3[i];
        q->next->next->z[i] = nniA.z4[i];
        q->next->next->back->z[i] = nniA.z4[i];
    }

	// Re-compute the partial likelihood vectors
	// TODO: One could save these instead of recomputation
    if (numBranches > 1 && !tr->useRecom) {
        pllUpdatePartials(tr, pr, p, PLL_TRUE);
        pllUpdatePartials(tr, pr, p->back, PLL_TRUE);
    } else {
        pllUpdatePartials(tr, pr, p, PLL_FALSE);
        pllUpdatePartials(tr, pr, p->back, PLL_FALSE);
    }

	return numPosNNI;
}

//void recomputePartial(pllInstance tr, partitionList pr, nodeptr p) {
//    if (numBranches > 1 && !tr->useRecom) {
//        pllUpdatePartials(tr, pr, p, PLL_TRUE);
//        pllUpdatePartials(tr, pr, p->back, PLL_TRUE);
//    } else {
//        pllUpdatePartials(tr, pr, p, PLL_FALSE);
//        pllUpdatePartials(tr, pr, p->back, PLL_FALSE);
//    }
//}

bool isAffectedBranch(nodeptr p, SearchInfo &searchinfo) {
	string branString = getBranString(p);
	if (searchinfo.aBranches.find(branString) != searchinfo.aBranches.end()) {
		return true;
	} else {
		return false;
	}
}

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		if (globalParams->speednni && searchinfo.curNumNNISteps != 1) {
			if (isAffectedBranch(p, searchinfo)) {
				evalNNIForBran(tr, pr, p, searchinfo);
			}
		} else {
			evalNNIForBran(tr, pr, p, searchinfo);
		}
		nodeptr q = p->next;
		while (q != p) {
			evalNNIForSubtree(tr, pr, q->back, searchinfo);
			q = q->next;
		}
	}
}

/**
* DTH:
* The PLL version of saveCurrentTree function
* @param tr: the tree (a pointer to a pllInstance)
* @param pr: pointer to a partitionList (this one keeps tons of tree info)
* @param p: root?
*/
void pllSaveCurrentTree(pllInstance* tr, partitionList *pr, nodeptr p){
    double cur_logl = tr->likelihood;

    struct pllHashItem * item_ptr = (struct pllHashItem *) malloc(sizeof(struct pllHashItem));
    item_ptr->data = (int *) malloc(sizeof(int));
    item_ptr->next = NULL;
    item_ptr->str = NULL;

    unsigned int tree_index = -1;
    char * tree_str = NULL;
    pllTree2StringREC(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE,
            PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_TRUE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
    tree_str = (char *) malloc (strlen(tr->tree_string) + 1);
#ifdef CLANG_UNDER_VS
    strcpy_s(tree_str, strlen(tr->tree_string) + 1, tr->tree_string);
#else
    strcpy(tree_str, tr->tree_string);
#endif
    pllBoolean is_stored = PLL_FALSE;

    if(is_stored){ // if found the tree_str
        pllUFBootDataPtr->duplication_counter++;
        tree_index = *((int *)item_ptr->data);
        if (cur_logl <= pllUFBootDataPtr->treels_logl[tree_index] + 1e-4) {
            if (cur_logl < pllUFBootDataPtr->treels_logl[tree_index] - 5.0)
                if (verbose_mode >= VB_MED)
                    printf("Current lh %f is much worse than expected %f\n",
                            cur_logl, pllUFBootDataPtr->treels_logl[tree_index]);
            return;
        }
        if (verbose_mode >= VB_MAX)
            printf("Updated logl %f to %f\n", pllUFBootDataPtr->treels_logl[tree_index], cur_logl);
        pllUFBootDataPtr->treels_logl[tree_index] = cur_logl;

        if (pllUFBootDataPtr->boot_samples == NULL) {
            (pllUFBootDataPtr->treels_ptnlh)[tree_index] =
                    (double *) malloc(pllUFBootDataPtr->n_patterns * sizeof(double));
            pllComputePatternLikelihood(tr, (pllUFBootDataPtr->treels_ptnlh)[tree_index], &cur_logl);
            return;
        }
        if (verbose_mode >= VB_MAX)
            printf("Update treels_logl[%d] := %f\n", tree_index, cur_logl);

    } else {
		if (pllUFBootDataPtr->logl_cutoff != 0.0 && cur_logl <= pllUFBootDataPtr->logl_cutoff + 1e-4){
			free(tree_str);       //James B. 23-Jul-2020 (memory leak)
			free(item_ptr->data); //James B. 23-Jul-2020 (memory leak)
			free(item_ptr);       //James B. 23-Jul-2020 (memory leak)
			return;
		}
		if (pllUFBootDataPtr->treels_size == pllUFBootDataPtr->candidate_trees_count) {
			pllResizeUFBootData();
		}
		tree_index = pllUFBootDataPtr->candidate_trees_count;
		pllUFBootDataPtr->candidate_trees_count++;
		pllUFBootDataPtr->treels_logl[tree_index] = cur_logl;
		if (verbose_mode >= VB_MAX)
			printf("Add    treels_logl[%d] := %f\n", tree_index, cur_logl);
   }

    //if (write_intermediate_trees)
    //        printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

    double *pattern_lh = (double *) malloc(pllUFBootDataPtr->n_patterns * sizeof(double));
	if (!pattern_lh) {
		outError("Not enough dynamic memory!");
	}
    pllComputePatternLikelihood(tr, pattern_lh, &cur_logl);

    if (pllUFBootDataPtr->boot_samples == NULL) {
        // for runGuidedBootstrap
        pllUFBootDataPtr->treels_ptnlh[tree_index] = pattern_lh;
		free(tree_str); //James B. 23-Jul-2020 (memory leak)
		free(item_ptr->data); //James B. 23-Jul-2020 (memory leak)
		free(item_ptr); //James B. 23-Jul-2020 (memory leak)
	}
	else {
		// online bootstrap
		int nptn = pllUFBootDataPtr->n_patterns;
		int updated = 0;
		int nsamples = globalParams->gbo_replicates;
		for (int sample = 0; sample < nsamples; sample++) {
			double rell = 0.0;
			for (int ptn = 0; ptn < nptn; ptn++)
				rell += pattern_lh[ptn] * pllUFBootDataPtr->boot_samples[sample][ptn];

			if (rell > pllUFBootDataPtr->boot_logl[sample] + globalParams->ufboot_epsilon ||
				(rell > pllUFBootDataPtr->boot_logl[sample] - globalParams->ufboot_epsilon &&
					random_double() <= 1.0 / (pllUFBootDataPtr->boot_counts[sample] + 1))) {
						{
							is_stored = pllHashSearch(pllUFBootDataPtr->treels, tree_str, &(item_ptr->data));
							if (is_stored)
								tree_index = *((int*)item_ptr->data);
							else {
								*((int*)item_ptr->data) = tree_index = pllUFBootDataPtr->candidate_trees_count - 1;
								item_ptr->str = tree_str;
								pllHashAdd(pllUFBootDataPtr->treels, pllHashString(tree_str, pllUFBootDataPtr->treels->size), tree_str, item_ptr->data);
							}
						}
						if (rell <= pllUFBootDataPtr->boot_logl[sample] +
							globalParams->ufboot_epsilon) {
							pllUFBootDataPtr->boot_counts[sample]++;
						}
						else {
							pllUFBootDataPtr->boot_counts[sample] = 1;
						}
						if (rell > pllUFBootDataPtr->boot_logl[sample])
							pllUFBootDataPtr->boot_logl[sample] = rell;
						pllUFBootDataPtr->boot_trees[sample] = tree_index;
						updated++;
			}
		}
	}

    if (pllUFBootDataPtr->boot_samples){
        free(pattern_lh);
        pllUFBootDataPtr->treels_ptnlh[tree_index] = NULL;
    }
}

/**
* DTH:
* Extract the array of site log likelihood to be kept in ptnlh
* And update *cur_log
* @param tr: the tree (pointer to an pllInstance)
* @param ptnlh: to-be-kept array of site log likelihood
* @param cur_logl: pointer to current tree log likelihood
*/
void pllComputePatternLikelihood(pllInstance* tr, double * ptnlh, double * cur_logl){
    int i;
    double tree_logl = 0;
    for(i = 0; i < pllUFBootDataPtr->n_patterns; i++){
        ptnlh[i] = tr->lhs[i];
        tree_logl += tr->lhs[i] * tr->aliaswgt[i];
    }
    *cur_logl = tree_logl;
}

/**
* DTH:
* Resize some of the arrays in UFBootData if they're full
* Along with update treels_size (to track the size of these arrays)
*/
void pllResizeUFBootData(){
    int count = pllUFBootDataPtr->candidate_trees_count;
    pllUFBootDataPtr->treels_size = 2 * count;

    double * rtreels_logl =
            (double *) malloc(2 * count * (sizeof(double)));
    if(!rtreels_logl) outError("Not enough dynamic memory!");
    memcpy(rtreels_logl, pllUFBootDataPtr->treels_logl, count * sizeof(double));
    free(pllUFBootDataPtr->treels_logl);
    pllUFBootDataPtr->treels_logl = rtreels_logl;

    double ** rtreels_ptnlh =
        (double **) malloc(2 * count * (sizeof(double *)));
    if(!rtreels_ptnlh) outError("Not enough dynamic memory!");
    memset(rtreels_ptnlh, 0, 2 * count * sizeof(double *));
    memcpy(rtreels_ptnlh, pllUFBootDataPtr->treels_ptnlh, count * sizeof(double *));
    free(pllUFBootDataPtr->treels_ptnlh);
    pllUFBootDataPtr->treels_ptnlh = rtreels_ptnlh;
}


/**
* DTH:
* (Based on function Tree2StringREC of PLL)
* Print out the tree topology with IQTree taxa ID (starts at 0) instead of PLL taxa ID (starts at 1)
* @param All are the same as in PLL's
* 2014.4.23: REPLACE getBranchLength(tr, pr, perGene, p) BY pllGetBranchLength(tr, p, pr->numberOfPartitions)
* BECAUSE OF LIB NEW DECLARATION: pllGetBranchLength (pllInstance *tr, nodeptr p, int partition_id);
*/
char *pllTree2StringREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean  printBranchLengths, pllBoolean  printNames,
		pllBoolean  printLikelihood, pllBoolean  rellTree, pllBoolean  finalPrint, int perGene, pllBoolean  branchLabelSupport, pllBoolean  printSHSupport)
{
	char * result = treestr; // DTH: added this var to be able to remove the '\n' at the end
	char  *nameptr;

	if(isTip(p->number, tr->mxtips)){
		if(printNames){
			nameptr = tr->nameList[p->number];
			sprintf(treestr, "%s", nameptr);
		}else
			sprintf(treestr, "%d", p->number - 1);

		while (*treestr) treestr++;
	}else{
		*treestr++ = '(';
		treestr = pllTree2StringREC(treestr, tr, pr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
				finalPrint, perGene, branchLabelSupport, printSHSupport);
		*treestr++ = ',';
		treestr = pllTree2StringREC(treestr, tr, pr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
				finalPrint, perGene, branchLabelSupport, printSHSupport);
		if(p == tr->start->back){
			*treestr++ = ',';
			treestr = pllTree2StringREC(treestr, tr, pr, p->back, printBranchLengths, printNames, printLikelihood, rellTree,
				finalPrint, perGene, branchLabelSupport, printSHSupport);
		}
		*treestr++ = ')';
	}

	if(p == tr->start->back){
		if(printBranchLengths && !rellTree)
			sprintf(treestr, ":0.0;\n");
		else
			sprintf(treestr, ";\n");
	}else{
		if(rellTree || branchLabelSupport || printSHSupport){
			if(( !isTip(p->number, tr->mxtips)) &&
					( !isTip(p->back->number, tr->mxtips)))
			{
				ASSERT(p->bInf != (branchInfo *)NULL);
				if(rellTree)
					sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
				if(branchLabelSupport)
					sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
				if(printSHSupport)
					sprintf(treestr, ":%8.20f[%d]", pllGetBranchLength(tr, p, pr->numberOfPartitions), p->bInf->support);
			}else{
				if(rellTree || branchLabelSupport)
					sprintf(treestr, ":%8.20f", p->z[0]);
				if(printSHSupport)
					sprintf(treestr, ":%8.20f", pllGetBranchLength(tr, p, pr->numberOfPartitions));
			}
		}else{
			if(printBranchLengths)
				sprintf(treestr, ":%8.20f", pllGetBranchLength(tr, p, pr->numberOfPartitions));
			else
				sprintf(treestr, "%s", "\0");
		}
	}

	if(result[strlen(result) - 1] == '\n') result[strlen(result) - 1] = '\0';
	while (*treestr) treestr++;
	return  treestr;
}


