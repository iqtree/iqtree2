/*
 * phylotreepars.cpp
 *
 * Fast implementation of parsimony kernel
 *
 *  Created on: May 18, 2015
 *      Author: minh
 */

#include "phylotree.h"

/***********************************************************/
/****** optimized version of parsimony kernel **************/
/***********************************************************/

void PhyloTree::computePartialParsimonyFast(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (dad_branch->partial_lh_computed & 2)
        return;
    Node *node = dad_branch->node;
    int ptn;
    int nptn = aln->size();
    int nstates = aln->num_states;
    int pars_size = getBitsBlockSize();

    if (node->isLeaf() && dad) {
        // external node
        int leafid = node->id;
        memset(dad_branch->partial_pars, 0, pars_size*sizeof(UINT));

        for (ptn = 0; ptn < nptn; ptn++) {
            UINT *p = dad_branch->partial_pars+(ptn/UINT_BITS);
        	int state = aln->at(ptn)[leafid];
        	UINT bit1 = (1 << (ptn%UINT_BITS));
        	StateBitset state_app;
        	aln->getAppearance(state, state_app);
        	for (int i = 0; i < nstates; i++)
        		if (state_app[i])
        			p[i] |= bit1;
        }
    } else {
        // internal node
    }
}


int PhyloTree::computeParsimonyBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(node_branch, node);
    int tree_pars = 0;

    return tree_pars;
}




