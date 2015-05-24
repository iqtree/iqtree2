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
    int nstates = aln->num_states;
    int site;

    if (node->isLeaf() && dad) {
        // external node
        int leafid = node->id;
        int pars_size = getBitsBlockSize();
        memset(dad_branch->partial_pars, 0, pars_size*sizeof(UINT));
        int ptn;
        int nptn = aln->size();
        for (ptn = 0, site = 0; ptn < nptn; ptn++) {
            if (!aln->at(ptn).is_informative)
                continue;
        	int state = aln->at(ptn)[leafid];
        	StateBitset state_app;
        	aln->getAppearance(state, state_app);
            int freq = aln->at(ptn).frequency;
            // duplicate entries corresponding to pattern frequency
            for (int j = 0; j < freq; j++, site++) {
                UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates+1);
                UINT bit1 = (1 << (site%UINT_BITS));
                for (int i = 0; i < nstates; i++)
                    if (state_app[i])
                        p[i] |= bit1;
            }
        }
        int max_sites = ((site+UINT_BITS-1)/UINT_BITS)*UINT_BITS;
        // add dummy states
        for (; site < max_sites; site++) {
            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates+1);
            UINT bit1 = (1 << (site%UINT_BITS));
            p[0] |= bit1;
        }
    } else {
        // internal node
        UINT *u = new UINT[nstates];
        assert(node->degree() == 3); // it works only for strictly bifurcating tree
        PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor* pit = (PhyloNeighbor*) (*it);
            if ((*it)->node->name != ROOT_NAME && (pit->partial_lh_computed & 2) == 0) {
                computePartialParsimonyFast(pit, (PhyloNode*) node);
            }
            if (!left) left = pit; else right = pit;
        }
        int score = left->partial_pars[0] + right->partial_pars[0];
        int nsites = aln->num_informative_sites;
        for (site = 0; site<nsites; site+=UINT_BITS) {
            int i;
            size_t offset = ((site/UINT_BITS)*nstates+1);
            UINT *x = left->partial_pars+offset;
            UINT *y = right->partial_pars+offset;
            UINT *z = dad_branch->partial_pars+offset;
            UINT w = 0;
            for (i = 0; i < nstates; i++) {
                u[i] |= x[i] * y[i];
                w |= u[i];
            }
            w = ~w;
            for (i = 0; i < nstates; i++) {
                z[i] = u[i] | (w & (x[i] | y[i]));
            }
            score += __builtin_popcount(w);
        }
        dad_branch->partial_pars[0] = score;
        delete [] u;
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
    int site;
    int nsites = aln->num_informative_sites;
    int nstates = aln->num_states;

    int sum = dad_branch->partial_pars[0] + node_branch->partial_pars[0];
    int score = 0;
    UINT *u = new UINT[nstates];
    
    for (site = 0; site < nsites; site+=UINT_BITS) {
        int i;
        size_t offset = ((site/UINT_BITS)*nstates+1);
        UINT *x = dad_branch->partial_pars+offset;
        UINT *y = node_branch->partial_pars+offset;
        UINT w = 0;
        for (i = 0; i < nstates; i++) {
            u[i] |= x[i] * y[i];
            w |= u[i];
        }
        w = ~w;
        score += __builtin_popcount(w);
        
    }
    if (branch_subst)
        *branch_subst = score;
    score += sum;
    
    delete [] u;
    return score;
}




