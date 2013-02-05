#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#include "nnisearch.h"

double doNNISearch(tree* tr) {
	double curScore = tr->likelihood;

	/* Initialize the NNI list */
	nniMove* nniList = (nniMove*) malloc((tr->ntips - 3) * sizeof(nniMove));
	int i;

	/* fill up the NNI list */
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	int cnt = 0; // number of visited internal branches during NNI evaluation
	int cnt_nni = 0; // number of positive NNI found

	while ( q != p ) {
		evalNNIForSubtree(tr, q->back, nniList, &cnt, &cnt_nni, curScore);
		q = q->next;
	}

#ifdef DEBUG_MAX
	printf("Number of internal branches traversed = %d \n", cnt);
	printf("Number of positive NNI found = %d \n", cnt_nni);
#endif

	if (cnt_nni == 0)
		return 0.0;

	nniMove* impNNIList = (nniMove*) malloc( cnt_nni * sizeof(nniMove) );
	int j = 0;
	for (i = 0; i < tr->ntips - 3; i++) {
		if (nniList[i].deltaLH > 0.0) {
			impNNIList[j] = nniList[i];
			j++;
		}
	}

#ifdef DEBUG_MAX
	puts("Before sorting");
	for (i = 0; i < cnt_nni; i++) {
		printf("NNI %d has delta = %f \n", i, impNNIList[i].deltaLH);
	}
	puts("After sorting");
#endif

	// sort impNNIList
	qsort(impNNIList, cnt_nni, sizeof(nniMove), cmp_nni);

#ifdef DEBUG_MAX
	for (i = 0; i < cnt_nni; i++) {
		printf("NNI %d has delta = %f \n", i, impNNIList[i].deltaLH);
	}
#endif

	// creating a list of non-conflicting positive NNI
	nniMove* nonConfNNIList = (nniMove*) calloc( cnt_nni, sizeof(nniMove) );

	// the best NNI will always be taken
	nonConfNNIList[0] = impNNIList[cnt_nni-1];

	// Filter out conflicting NNI
	int numNonConflictNNI = 1; // size of the non-conflicting NNI list;
	int k;
	for (k = cnt_nni-2; k >= 0; k--) {
		int conflict = FALSE;
		int j;
		for (j = 0; j < numNonConflictNNI; j++) {
			if (impNNIList[k].p->number == nonConfNNIList[j].p->number || impNNIList[k].p->number == nonConfNNIList[j].p->back->number ) {
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

	// Applying non-conflicting NNI moves
	double delta = 1.0; // portion of NNI moves to apply
	int notImproved;
	do {
		notImproved = FALSE;
		//printf("numNonConflictNNI = %d \n", numNonConflictNNI);
		int numNNI2Apply = ceil(numNonConflictNNI * delta);
		//printf("numNNI2Apply = %d \n", numNNI2Apply);
		for (i = 0; i < numNNI2Apply; i++) {
#ifdef DEBUG_MAX
			printf("Doing an improving NNI %d - %d with deltaLH = %f \n", nonConfNNIList[i].p->number, nonConfNNIList[i].p->back->number, nonConfNNIList[i].deltaLH);
#endif
			// Just do the topological change
			doOneNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType, FALSE);
			newviewGeneric(tr, nonConfNNIList[i].p, FALSE);
			newviewGeneric(tr, nonConfNNIList[i].p->back, FALSE);
			// Apply the store branch length
			int j;
			for (j = 0; j < tr->numBranches; j++) {
				nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z[j];
				nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z[j];
			}
		}
		// Re-optimize all branches
		smoothTree(tr, 2);
		//modOpt(tr, 0.1);
		evaluateGeneric(tr, tr->start, FALSE);
		if (tr->likelihood < curScore) {
			printf("Tree likelihood gets worse after applying NNI\n");
			printf("curScore = %30.20f\n", curScore);
			printf("newScore = %30.20f\n", tr->likelihood);
			if (numNNI2Apply == 1) {
				printf("What the FUCK?\n");
				exit(1);
			}
			printf("Rolling back the tree\n");
			for (i = 0; i < numNNI2Apply; i++) {
				doOneNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType,
						FALSE);
				// Restore the branch length
				int j;
				for (j = 0; j < tr->numBranches; j++) {
					nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
					nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
				}
			}
			evaluateGeneric(tr, tr->start, FALSE);
			printf("Tree likelihood after rolling back = %f \n", tr->likelihood);
			notImproved = TRUE;
			delta = delta * 0.5;
		}
	} while (notImproved);

	free(nniList);
	free(impNNIList);
	free(nonConfNNIList);

	return tr->likelihood;
}

double doOneNNI(tree * tr, nodeptr p, int swap, int optBran) {
	nodeptr q;
	nodeptr tmp;

	q = p->back;
	//printTopology(tr, TRUE);
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));

	int pNum = p->number;
	int qNum = q->number;

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
		//hookupDefault(p->next, q->next->back, tr->numBranches);
		//hookupDefault(q->next, tmp, tr->numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, tr->numBranches);
		hookup(q->next, tmp, tmp->z, tr->numBranches);
		//hookup(p->next->next, q->next->back, q->next->z, tr->numBranches);
		//hookup(q->next, tmp, tmp->z, tr->numBranches);
	}

	assert(pNum == p->number);
	assert(qNum == q->number);

	if (optBran) {
		newviewGeneric(tr, p, FALSE);
		newviewGeneric(tr, q, FALSE);
		update(tr, p);
		//printf("New branch length %f \n", getBranchLength(tr, 0, p) );
		evaluateGeneric(tr,p, FALSE);
		return tr->likelihood;
	} else {
		//newviewGeneric(tr, p, FALSE);
		//newviewGeneric(tr, q, FALSE);
		return -1.0;
	}
}

nniMove getBestNNIForBran(tree* tr, nodeptr p, double curLH) {
	nodeptr q = p->back;
	assert ( ! isTip(p->number, tr->mxtips) );
	assert ( ! isTip(q->number, tr->mxtips) );

#ifdef DEBUG_MAX
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
#endif

	/* Backup the current branch length */
	double z0[NUM_BRANCHES];
	int i;
	for(i = 0; i < tr->numBranches; i++) {
		z0[i] = p->z[i];
	}
#ifdef DEBUG_MAX
	double lhOld = tr->likelihood;
	printf("lhOld: %f \n", lhOld);
#endif

	//update(tr, p);
	//evaluateGeneric(tr, p, FALSE);
	//printf("Current tree LH = %f \n", tr->likelihood);
	//double lh0 = tr->likelihood;
	double lh0 = curLH;
	//printf("zNew: %f \n", getBranchLength(tr, 0, p));

#ifdef DEBUG_MAX
	printf("lh0: %f \n", lh0);
#endif

	nniMove nni0; // nni0 means no NNI move is done
	nni0.p = p;
	nni0.nniType = 0;
	nni0.deltaLH = 0;
	for (i = 0; i < tr->numBranches; i++) {
		nni0.z[i] = p->z[i];
	}

	/* TODO Save the likelihood vector at node p and q */
	//saveLHVector(p, q, p_lhsave, q_lhsave);
	/* Save the scaling factor */

	// Now try to do an NNI move of type 1
	double lh1 = doOneNNI(tr, p, 1, TRUE);
	nniMove nni1;
	nni1.p = p;
	nni1.nniType = 1;
	// Store the optimized und unoptimized central branch length
	for (i = 0; i < tr->numBranches; i++) {
		nni1.z[i] = p->z[i];
		nni1.z0[i] = z0[i];
	}
	nni1.likelihood = lh1;
	nni1.deltaLH = lh1 - lh0;

#ifdef DEBUG_MAX
	printf("Delta likelihood of the 1.NNI move: %f\n", nni1.deltaLH);
	//printTopology(tr, TRUE);
#endif

	/* Restore previous NNI move */
	doOneNNI(tr, p, 1, FALSE);
	/* Restore the old branch length */
	for(i = 0; i < tr->numBranches; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}

#ifdef DEBUG_MAX
	printf("Restore topology\n");
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 1: %f\n", tr->likelihood);
#endif

	/* Try to do an NNI move of type 2 */
	double lh2 = doOneNNI(tr, p, 2, TRUE);
	// Create the nniMove struct to store this move
	nniMove nni2;
	nni2.p = p;
	nni2.nniType = 2;
	// Store the optimized and unoptimized central branch length
	for (i = 0; i < tr->numBranches; i++) {
		nni2.z[i] = p->z[i];
		nni2.z0[i] = z0[i];
	}
	nni2.likelihood = lh2;
	nni2.deltaLH = lh2 - lh0;

#ifdef DEBUG_MAX
	printf("Delta likelihood of the 2.NNI move: %f\n", nni2.deltaLH);
	//printTopology(tr, TRUE);
#endif

	/* Restore previous NNI move */
	doOneNNI(tr, p, 2, FALSE);
	newviewGeneric(tr, p, FALSE);
	newviewGeneric(tr, p->back, FALSE);
	/* Restore the old branch length */
	for(i = 0; i < tr->numBranches; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}

#ifdef DEBUG_MAX
	printf("Restore topology\n");
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	evaluateGeneric(tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 2: %f\n", tr->likelihood);
#endif

	if (nni1.deltaLH > 0 && nni1.deltaLH >= nni2.deltaLH) {
		return nni1;
	} else if (nni1.deltaLH > 0 && nni1.deltaLH < nni2.deltaLH ){
		return nni2;
	} else if (nni1.deltaLH < 0 && nni2.deltaLH > 0) {
		return nni2;
	} else {
		return nni0;
	}


	/******************** NNI part *******************************/

	/* Restore the likelihood vector */
	//restoreLHVector(p,q, p_lhsave, q_lhsave);

}

int main(int argc, char * argv[])
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

	/* read the binary input, setup tree, initialize model with alignment */
	read_msa(tr,argv[1]);

	tr->randomNumberSeed = 665;

	/* Create random tree */
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

	/* compute the LH of the full tree */
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

	/* Model optimization */
	printf("Optimizing model parameters and branch lengths ... \n");
	getrusage(RUSAGE_SELF, &usage);
	start_tmp = usage.ru_utime;
	modOpt(tr, 0.1);
	getrusage(RUSAGE_SELF, &usage);
	end_tmp = usage.ru_utime;
	printf("Model parameters and branch lengths optimized in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);
	evaluateGeneric(tr, tr->start, FALSE);
	printf("Likelihood of the parsimony tree: %f\n", tr->likelihood);

	/* 8 rounds of branch length optimization */
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
}

void evalNNIForSubtree(tree* tr, nodeptr p, nniMove* nniList, int* cnt, int* cnt_nni, double curLH) {
	if ( ! isTip(p->number, tr->mxtips) ) {
		//newviewGeneric(tr, p, FALSE);
		//newviewGeneric(tr, p->back, FALSE);
		nniList[*cnt] = getBestNNIForBran(tr, p, curLH);
		if (nniList[*cnt].deltaLH != 0.0) {
			*cnt_nni = *cnt_nni + 1;
		}
		*cnt = *cnt + 1;
		nodeptr q = p->next;
		while ( q != p ) {
			evalNNIForSubtree(tr, q->back, nniList, cnt, cnt_nni, curLH);
			q = q->next;
		}
	}
}




