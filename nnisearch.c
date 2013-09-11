
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylolib/axml.h"
#define GLOBAL_VARIABLES_DEFINITION
#include "phylolib/globalVariables.h"
#include "nnisearch.h"
#include "phylolib.h"
#include <math.h>
#include "fmemopen.h"
//#include "tools.h"

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


double TOL_LIKELIHOOD_PHYLOLIB;
int numSmoothTree;
int fast_eval;
int fivebran;

int treeReadLenString(const char *buffer, tree *tr, pl_boolean readBranches,
        pl_boolean readNodeLabels, pl_boolean topologyOnly) {
    FILE *stream;

    stream = fmemopen((char*) buffer, strlen(buffer), "r");
    int lcount = treeReadLen(stream, tr, readBranches, readNodeLabels, topologyOnly);
    fclose(stream);
    return lcount;
}

int cmp_nni(const void* nni1, const void* nni2) {
    nniMove* myNNI1 = (nniMove*) nni1;
    nniMove* myNNI2 = (nniMove*) nni2;
    if (myNNI1->likelihood > myNNI2->likelihood)
    	return 1;
    else if (myNNI1->likelihood < myNNI2->likelihood)
    	return -1;
    else
    	return 0;
}

int compareDouble(const void * a, const void * b) {
    if (*(double*) a > *(double*) b) return 1;
    else if (*(double*) a < *(double*) b) return -1;
    else return 0;
}

double doNNISearch(tree* tr, int* nni_count, double* deltaNNI, NNICUT* nnicut, int numSmooth) {
    double curScore = tr->likelihood;

    /* Initialize the NNI list */
    nniMove* nniList = (nniMove*) malloc((tr->ntips - 3) * sizeof (nniMove));
    int i;

    /* fill up the NNI list */
    nodeptr p = tr->start->back;
    nodeptr q = p->next;
    int cnt = 0; // number of visited internal branches during NNI evaluation
    int cnt_nni = 0; // number of positive NNI found

    while (q != p) {
        evalNNIForSubtree(tr, q->back, nniList, &cnt, &cnt_nni, curScore, nnicut);
        q = q->next;
    }

    if (cnt_nni == 0)
        return 0.0;

    nniMove* impNNIList = (nniMove*) malloc(cnt_nni * sizeof (nniMove));
    int j = 0;
    for (i = 0; i < tr->ntips - 3; i++) {
        if (nniList[i].deltaLH > 0.0) {
            impNNIList[j] = nniList[i];
            j++;
        }
    }

    // sort impNNIList
    qsort(impNNIList, cnt_nni, sizeof (nniMove), cmp_nni);

    if (verbose_mode >= VB_DEBUG)
    {
    	int i;
    	for (i = cnt_nni-1; i >= 0; i--) {
    		printf("Log-likelihood of positive NNI %d : %10.6f \n", cnt_nni - i - 1, impNNIList[i].likelihood);
    	}
    }

    // creating a list of non-conflicting positive NNI
    nniMove* nonConfNNIList = (nniMove*) malloc(cnt_nni * sizeof (nniMove));

    // the best NNI will always be taken
    nonConfNNIList[0] = impNNIList[cnt_nni - 1];

    // Filter out conflicting NNI
    int numNonConflictNNI = 1; // size of the non-conflicting NNI list;
    int k;
    for (k = cnt_nni - 2; k >= 0; k--) {
        int conflict = FALSE;
        int j;
        for (j = 0; j < numNonConflictNNI; j++) {
            if (impNNIList[k].p->number == nonConfNNIList[j].p->number || impNNIList[k].p->number == nonConfNNIList[j].p->back->number) {
                conflict = TRUE;
                break;
            }
            if (impNNIList[k].p->back->number == nonConfNNIList[j].p->number || impNNIList[k].p->back->number == nonConfNNIList[j].p->back->number) {
                conflict = TRUE;
                break;
            }
        }
        if (conflict) {
            continue;
        } else {
            nonConfNNIList[numNonConflictNNI] = impNNIList[k];
            numNonConflictNNI++;
        }
    }
    if (verbose_mode >= VB_DEBUG)
    {
    	int i;
    	for (i = 0; i < numNonConflictNNI; i++) {
    		printf("Log-likelihood of non-conflicting NNI %d : %10.6f \n", i, nonConfNNIList[i].likelihood);
    		printf("Log-likelihood before NNI : %10.6f \n", nonConfNNIList[i].likelihood - nonConfNNIList[i].deltaLH);
    	}
    }

    // Applying non-conflicting NNI moves
    double delta = 1.0; // portion of NNI moves to apply
    int notImproved;
    int numNNI2Apply;
    do {
        notImproved = FALSE;
        numNNI2Apply = ceil(numNonConflictNNI * delta);
        for (i = 0; i < numNNI2Apply; i++) {
            // Just do the topological change
            doOneNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType, TOPO_ONLY, curScore);
            // Apply the store branch length
            int j;
            for (j = 0; j < tr->numBranches; j++) {
                nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
                nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
                nonConfNNIList[i].p->next->z[j] = nonConfNNIList[i].z1[j];
                nonConfNNIList[i].p->next->back->z[j] = nonConfNNIList[i].z1[j];
                nonConfNNIList[i].p->next->next->z[j] = nonConfNNIList[i].z2[j];
                nonConfNNIList[i].p->next->next->back->z[j] = nonConfNNIList[i].z2[j];
                nonConfNNIList[i].p->back->next->z[j] = nonConfNNIList[i].z3[j];
                nonConfNNIList[i].p->back->next->back->z[j] = nonConfNNIList[i].z3[j];
                nonConfNNIList[i].p->back->next->next->z[j] = nonConfNNIList[i].z4[j];
                nonConfNNIList[i].p->back->next->next->back->z[j] = nonConfNNIList[i].z4[j];
            }
            newviewGeneric(tr, nonConfNNIList[i].p, FALSE);
            newviewGeneric(tr, nonConfNNIList[i].p->back, FALSE);
        }

        // Re-optimize all branches
        treeEvaluate(tr, numSmooth);
        //treeEvaluate(tr, 32);
        if (tr->likelihood < curScore - TOL_LIKELIHOOD_PHYLOLIB) {
            if (numNNI2Apply == 1) {
                printf("This mighbt be a BUG: Tree gets worse when 1 NNI is applied?\n");
                printf("Tree supposed to get LH greater than or equal %30.20f\n", nonConfNNIList[0].likelihood);
//                printf("Rolling back the tree\n");
//                for (i = numNNI2Apply - 1; i >= 0; i--) {
//                    doOneNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType, TOPO_ONLY, curScore);
//                    // Restore the branch length
//                    int j;
//                    for (j = 0; j < tr->numBranches; j++) {
//                        nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
//                        nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
//                    }
//                    newviewGeneric(tr, nonConfNNIList[i].p, FALSE);
//                    newviewGeneric(tr, nonConfNNIList[i].p->back, FALSE);
//                }
                return tr->likelihood;
            }
            //printf("Rolling back the tree\n");
            for (i = numNNI2Apply - 1; i >= 0; i--) {
                //printf("Swaping back node (%d) -- (%d) \n", nonConfNNIList[i].p->number, nonConfNNIList[i].p->back->number);
                doOneNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType, TOPO_ONLY, curScore);
                // Restore the branch length
                int j;
                for (j = 0; j < tr->numBranches; j++) {
                    nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
                    nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
                }
                newviewGeneric(tr, nonConfNNIList[i].p, FALSE);
                newviewGeneric(tr, nonConfNNIList[i].p->back, FALSE);
            }
            notImproved = TRUE;
            delta = delta * 0.5;
        }
    } while (notImproved);

    free(nniList);
    free(impNNIList);
    free(nonConfNNIList);
    *nni_count = numNNI2Apply;
    *deltaNNI = (tr->likelihood - curScore) / numNNI2Apply;
    return tr->likelihood;
}

void optimizeFiveBranches(tree* tr, nodeptr p) {

}

void optimizeOneBranches(tree* tr, nodeptr p, int numNRStep) {
    nodeptr  q;
    int i;
    double   z[NUM_BRANCHES], z0[NUM_BRANCHES];

    q = p->back;

    for(i = 0; i < tr->numBranches; i++)
      z0[i] = q->z[i];

    if(tr->numBranches > 1)
      makenewzGeneric(tr, p, q, z0, numNRStep, z, TRUE);
    else
      makenewzGeneric(tr, p, q, z0, numNRStep, z, FALSE);

    for(i = 0; i < tr->numBranches; i++)
    {
      if(!tr->partitionConverged[i])
      {
        if(ABS(z[i] - z0[i]) > deltaz)
        {
          tr->partitionSmoothed[i] = FALSE;
        }

        p->z[i] = q->z[i] = z[i];
      }
    }
}

double doOneNNI(tree * tr, nodeptr p, int swap, int evalType, double curLH) {
	nodeptr q;
	nodeptr tmp;

	q = p->back;
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));

	//printTopology(tr, TRUE);

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
	}

	assert(pNum == p->number);
	assert(qNum == q->number);

	if (evalType == TOPO_ONLY) {
		return 0.0;
	} else if (evalType == ONE_BRAN_OPT) {
		newviewGeneric(tr, q, FALSE);
		newviewGeneric(tr, p, FALSE);
		update(tr, p);
		evaluateGeneric(tr, p, FALSE);
		//printf("log-likelihood of bran 1: %f \n", tr->likelihood);
	} else {
		newviewGeneric(tr, q, FALSE);
		newviewGeneric(tr, p, FALSE);
		update(tr, p);
		nodeptr r; // temporary node poiter
		// optimize 2 branches at node p
		r = p->next;
		newviewGeneric(tr, r, FALSE);
		update(tr, r);
		r = p->next->next;
		newviewGeneric(tr, r, FALSE);
		update(tr, r);
		// optimize 2 branches at node q
		r = q->next;
		newviewGeneric(tr, r, FALSE);
		update(tr, r);
		r = q->next->next;
		newviewGeneric(tr, r, FALSE);
		update(tr, r);
		evaluateGeneric(tr, r, FALSE);
	}
	return tr->likelihood;

}

nniMove getBestNNIForBran(tree* tr, nodeptr p, double curLH, NNICUT* nnicut) {
    nodeptr q = p->back;
    assert(!isTip(p->number, tr->mxtips));
    assert(!isTip(q->number, tr->mxtips));
    int i;
    nniMove nni0; // nni0 means no NNI move is done
    nni0.p = p;
    nni0.nniType = 0;
    nni0.deltaLH = 0.0;
    nni0.likelihood = curLH;
    //evaluateGeneric(tr, p, TRUE);
    //printf("Current log-likelihood: %f\n", tr->likelihood);
    for (i = 0; i < tr->numBranches; i++) {
        nni0.z0[i] = p->z[i];
        nni0.z1[i] = p->next->z[i];
        nni0.z2[i] = p->next->next->z[i];
        nni0.z3[i] = q->next->z[i];
        nni0.z4[i] = q->next->next->z[i];
    }

    double lh0 = curLH;
//    double multiLH = 0.0;
//    if (nnicut->doNNICut) {
//        // compute likelihood of the multifurcating tree
//        for (i = 0; i < tr->numBranches; i++) {
//            p->z[i] = 1.0;
//            p->back->z[i] = 1.0;
//        }
//        // now compute the log-likelihood
//        // This is actually time consuming!!!
//        newviewGeneric(tr, p, FALSE);
//        newviewGeneric(tr, p->back, FALSE);
//        evaluateGeneric(tr, p, FALSE);
//        multiLH = tr->likelihood;
//        //printf("curLH - multiLH : %f - %f \n", curLH, multiLH);
//        for (i = 0; i < tr->numBranches; i++) {
//            p->z[i] = z0[i];
//            p->back->z[i] = z0[i];
//        }
//        // If the log-likelihood of the zero branch configuration is some delta log-likelihood smaller than
//        // the log-likelihood of the current tree then it is very likely that there no better NNI tree
//        if (curLH - multiLH > nnicut->delta_min)
//            return nni0;
//    }

    /* TODO Save the likelihood vector at node p and q */
    //saveLHVector(p, q, p_lhsave, q_lhsave);
    /* Save the scaling factor */

    // Now try to do an NNI move of type 1
    double lh1;
    if (!fast_eval && !fivebran) {
        lh1 = doOneNNI(tr, p, 1, ONE_BRAN_OPT, curLH);
    } else if (fivebran) {
    	lh1 = doOneNNI(tr, p, 1, FIVE_BRAN_OPT, curLH);
    } else {
        lh1 = doOneNNI(tr, p, 1, NO_BRAN_OPT, curLH);
    }
    nniMove nni1;
    nni1.p = p;
    nni1.nniType = 1;
    // Store the optimized and unoptimized central branch length
    for (i = 0; i < tr->numBranches; i++) {
        nni1.z0[i] = p->z[i];
        nni1.z1[i] = p->next->z[i];
        nni1.z2[i] = p->next->next->z[i];
        nni1.z3[i] = q->next->z[i];
        nni1.z4[i] = q->next->next->z[i];
    }
    nni1.likelihood = lh1;
    nni1.deltaLH = lh1 - lh0;

    /* Restore previous NNI move */
    doOneNNI(tr, p, 1, TOPO_ONLY, curLH);
    /* Restore the old branch length */
    for (i = 0; i < tr->numBranches; i++) {
        p->z[i] = nni0.z0[i];
        q->z[i] =  nni0.z0[i];
        p->next->z[i] =  nni0.z1[i];
        p->next->back->z[i] =  nni0.z1[i];
        p->next->next->z[i] =  nni0.z2[i];
        p->next->next->back->z[i] =  nni0.z2[i];
        q->next->z[i] =  nni0.z3[i];
        q->next->back->z[i] =  nni0.z3[i];
        q->next->next->z[i] =  nni0.z4[i];
        q->next->next->back->z[i] =  nni0.z4[i];
    }

    /* Try to do an NNI move of type 2 */
    double lh2;
    if (!fast_eval && !fivebran) {
        lh2 = doOneNNI(tr, p, 2, ONE_BRAN_OPT, curLH);
    } else if (fivebran) {
    	lh2 = doOneNNI(tr, p, 2, FIVE_BRAN_OPT, curLH);
    } else {
    	lh2 = doOneNNI(tr, p, 2, NO_BRAN_OPT, curLH);
    }
    // Create the nniMove struct to store this move
    nniMove nni2;
    nni2.p = p;
    nni2.nniType = 2;
    // Store the optimized and unoptimized central branch length
    for (i = 0; i < tr->numBranches; i++) {
        nni2.z0[i] = p->z[i];
        nni2.z1[i] = p->next->z[i];
        nni2.z2[i] = p->next->next->z[i];
        nni2.z3[i] = q->next->z[i];
        nni2.z4[i] = q->next->next->z[i];
    }
    nni2.likelihood = lh2;
    nni2.deltaLH = lh2 - lh0;

    /* Restore previous NNI move */
    doOneNNI(tr, p, 2, TOPO_ONLY, curLH);
    /* Restore the old branch length */
    for (i = 0; i < tr->numBranches; i++) {
        p->z[i] = nni0.z0[i];
        q->z[i] =  nni0.z0[i];
        p->next->z[i] =  nni0.z1[i];
        p->next->back->z[i] =  nni0.z1[i];
        p->next->next->z[i] =  nni0.z2[i];
        p->next->next->back->z[i] =  nni0.z2[i];
        q->next->z[i] =  nni0.z3[i];
        q->next->back->z[i] =  nni0.z3[i];
        q->next->next->z[i] =  nni0.z4[i];
        q->next->next->back->z[i] =  nni0.z4[i];
    }
    newviewGeneric(tr, p, FALSE);
    newviewGeneric(tr, p->back, FALSE);


//    if (nnicut->doNNICut && (nnicut->num_delta) < MAX_NUM_DELTA) {
//        double lh_array[3];
//        lh_array[0] = curLH;
//        lh_array[1] = lh1;
//        lh_array[2] = lh2;
//
//        qsort(lh_array, 3, sizeof (double), compareDouble);
//
//        double deltaLH = lh_array[1] - multiLH;
//        if (deltaLH > nnicut->delta_min) {
//            printf("%f %f %f %f\n", multiLH, curLH, lh1, lh2);
//        }
//        nnicut->delta[nnicut->num_delta] = deltaLH;
//        (nnicut->num_delta)++;
//    }

    if (nni1.deltaLH > TOL_LIKELIHOOD_PHYLOLIB && nni1.deltaLH > nni2.deltaLH) {
        return nni1;
    } else if (nni1.deltaLH > TOL_LIKELIHOOD_PHYLOLIB && nni1.deltaLH < nni2.deltaLH) {
//        if (verbose_mode >= VB_DEBUG) {
//            printf("Both NNIs are positive. NNI 1 has log-lh = %10.6f and NNI 2 has log-lh = %10.6f \n", nni1.likelihood, nni2.likelihood);
//            if ( ABS(nni1.likelihood - nni2.likelihood) < 0.1 ) {
//                float p = (float) rand() / RAND_MAX;
//                printf("Probability = %10.6f \n", p);
//                if (  p < 0.5 ) {
//                    return nni1;
//                } else {
//                    return nni2;
//                }
//            }
//        }
        return nni2;
    } else if (nni1.deltaLH < TOL_LIKELIHOOD_PHYLOLIB && nni2.deltaLH > TOL_LIKELIHOOD_PHYLOLIB) {
        return nni2;
    } else {
        return nni0;
    }

    /* Restore the likelihood vector */
    //restoreLHVector(p,q, p_lhsave, q_lhsave);

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


void evalNNIForSubtree(tree* tr, nodeptr p, nniMove* nniList, int* cnt, int* cnt_nni, double curLH, NNICUT* nnicut) {
    if (!isTip(p->number, tr->mxtips)) {
        nniList[*cnt] = getBestNNIForBran(tr, p, curLH, nnicut);
        if (nniList[*cnt].deltaLH > 0.0) {
            *cnt_nni = *cnt_nni + 1;
        }
        *cnt = *cnt + 1;
        nodeptr q = p->next;
        while (q != p) {
            evalNNIForSubtree(tr, q->back, nniList, cnt, cnt_nni, curLH, nnicut);
            q = q->next;
        }
    }
}
