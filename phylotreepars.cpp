/*
 * phylotreepars.cpp
 *
 * Fast implementation of parsimony kernel
 *
 *  Created on: May 18, 2015
 *      Author: minh
 */

#include "phylotree.h"
#include "vectorclass/vectorclass.h"

/***********************************************************/
/****** optimized version of parsimony kernel **************/
/***********************************************************/

void PhyloTree::computePartialParsimonyFast(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (dad_branch->partial_lh_computed & 2)
        return;
    Node *node = dad_branch->node;
    int nstates = aln->num_states;
    int site;

    dad_branch->partial_lh_computed |= 2;

    if (node->isLeaf() && dad) {
        // external node
        if (aln->ordered_pattern.empty())
            aln->orderPatternByNumChars();
        int leafid = node->id;
        int pars_size = getBitsBlockSize();
        memset(dad_branch->partial_pars, 0, pars_size*sizeof(UINT));
//        int ptn;
//        int nptn = aln->size();
    	int ambi_aa[] = {2, 3, 5, 6, 9, 10}; // {4+8, 32+64, 512+1024};
        int max_sites = ((aln->num_informative_sites+UINT_BITS-1)/UINT_BITS)*UINT_BITS;
        Alignment::iterator pat;
    	switch (aln->seq_type) {
    	case SEQ_DNA:
//            nptn = aln->ordered_pattern.size();
            for (pat = aln->ordered_pattern.begin(), site = 0; pat != aln->ordered_pattern.end(); pat++) {
//                Pattern *pat = &aln->ordered_pattern[ptn];
//                if (!pat->is_informative)
//                    continue;
            	int state = pat->at(leafid);
                int freq = pat->frequency;
                if (state < 4) {
                    for (int j = 0; j < freq; j++, site++) {
                        dad_branch->partial_pars[(site/UINT_BITS)*4+state] |= (1 << (site % UINT_BITS));
                    }
                } else if (state == aln->STATE_UNKNOWN) {
                    for (int j = 0; j < freq; j++, site++) {
                        UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*4);
                        UINT bit1 = (1 << (site%UINT_BITS));
                        p[0] |= bit1;
                        p[1] |= bit1;
                        p[2] |= bit1;
                        p[3] |= bit1;
                    }
                } else {
                	state -= 3;
                    for (int j = 0; j < freq; j++, site++) {
                        UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*4);
                        UINT bit1 = (1 << (site%UINT_BITS));
                        for (int i = 0; i < 4; i++)
                            if (state & (1<<i))
                                p[i] |= bit1;
                    }
                }
            }
            assert(site == aln->num_informative_sites);
            // add dummy states
            if (site < max_sites)
            	dad_branch->partial_pars[(site/UINT_BITS)*4] |= ~((1<<(site%UINT_BITS)) - 1);
//            for (; site < max_sites; site++) {
//                dad_branch->partial_pars[(site/UINT_BITS)*4] |= (1 << (site%UINT_BITS));
//            }
    		break;
    	case SEQ_PROTEIN:
            for (pat = aln->ordered_pattern.begin(), site = 0; pat != aln->ordered_pattern.end(); pat++) {
//                if (!aln->at(ptn).is_informative)
//                    continue;
            	int state = pat->at(leafid);
                int freq = pat->frequency;
                if (state < 20) {
                    for (int j = 0; j < freq; j++, site++) {
                        dad_branch->partial_pars[(site/UINT_BITS)*20+state] |= (1 << (site % UINT_BITS));
                    }
                } else if (state == aln->STATE_UNKNOWN) {
                    for (int j = 0; j < freq; j++, site++) {
                        UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*20);
                        UINT bit1 = (1 << (site%UINT_BITS));
                        for (int i = 0; i < 20; i++)
                                p[i] |= bit1;
                    }
                } else {
                	assert(state < 23);
            		state = (state-20)*2;
                    for (int j = 0; j < freq; j++, site++) {
                        UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*20);
                        UINT bit1 = (1 << (site%UINT_BITS));
                        p[ambi_aa[state]] |= bit1;
                        p[ambi_aa[state+1]] |= bit1;
                    }
                }
            }
            assert(site == aln->num_informative_sites);
            // add dummy states
            if (site < max_sites)
            	dad_branch->partial_pars[(site/UINT_BITS)*20] |= ~((1<<(site%UINT_BITS)) - 1);
//            for (; site < max_sites; site++) {
//                dad_branch->partial_pars[(site/UINT_BITS)*20] |= (1 << (site%UINT_BITS));
//            }
    		break;
    	default:
//            for (ptn = 0, site = 0; ptn < nptn; ptn++) {
            for (pat = aln->ordered_pattern.begin(), site = 0; pat != aln->ordered_pattern.end(); pat++) {
//                if (!aln->at(ptn).is_informative)
//                    continue;
            	int state = pat->at(leafid);
                int freq = pat->frequency;
                if (state < nstates) {
                    for (int j = 0; j < freq; j++, site++) {
                        dad_branch->partial_pars[(site/UINT_BITS)*nstates+state] |= (1 << (site % UINT_BITS));
                    }
                } else if (state == aln->STATE_UNKNOWN) {
                    for (int j = 0; j < freq; j++, site++) {
                        UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                        UINT bit1 = (1 << (site%UINT_BITS));
                        for (int i = 0; i < nstates; i++)
                                p[i] |= bit1;
                    }
                } else {
                	assert(0);
                }
            }
            assert(site == aln->num_informative_sites);
            // add dummy states
            if (site < max_sites)
            	dad_branch->partial_pars[(site/UINT_BITS)*nstates] |= ~((1<<(site%UINT_BITS)) - 1);
//            for (; site < max_sites; site++) {
//                dad_branch->partial_pars[(site/UINT_BITS)*nstates] |= (1 << (site%UINT_BITS));
//            }
    		break;
    	}

    } else {
        // internal node
        assert(node->degree() == 3); // it works only for strictly bifurcating tree
        PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNeighbor* pit = (PhyloNeighbor*) (*it);
            if ((*it)->node->name != ROOT_NAME && (pit->partial_lh_computed & 2) == 0) {
                computePartialParsimonyFast(pit, (PhyloNode*) node);
            }
            if (!left) left = pit; else right = pit;
        }
//        UINT score = left->partial_pars[0] + right->partial_pars[0];
        UINT score = 0;
        int nsites = aln->num_informative_sites;
        UINT *x = left->partial_pars;
        UINT *y = right->partial_pars;
        UINT *z = dad_branch->partial_pars;

        switch (nstates) {
        case 4:
			for (site = 0; site<nsites; site+=UINT_BITS) {
				UINT w;
				z[0] = x[0] & y[0];
				z[1] = x[1] & y[1];
				z[2] = x[2] & y[2];
				z[3] = x[3] & y[3];
				w = z[0] | z[1] | z[2] | z[3];
				w = ~w;
				score += vml_popcnt(w);
				z[0] |= w & (x[0] | y[0]);
				z[1] |= w & (x[1] | y[1]);
				z[2] |= w & (x[2] | y[2]);
				z[3] |= w & (x[3] | y[3]);
				x += 4;
				y += 4;
				z += 4;
			}
			break;
        default:
			for (site = 0; site<nsites; site+=UINT_BITS) {
				int i;
				UINT w = 0;
				for (i = 0; i < nstates; i++) {
					z[i] = x[i] & y[i];
					w |= z[i];
				}
				w = ~w;
				score += vml_popcnt(w);
				for (i = 0; i < nstates; i++) {
					z[i] |= w & (x[i] | y[i]);
				}
				x += nstates;
				y += nstates;
				z += nstates;
			}
			break;
        }
        *z = score + *x + *y;
//        dad_branch->partial_pars[0] = score;
    }
}


int PhyloTree::computeParsimonyBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(node_branch, node);
    int site;
    int nsites = (aln->num_informative_sites + UINT_BITS-1) / UINT_BITS;
    int nstates = aln->num_states;

    int scoreid = ((aln->num_informative_sites+UINT_BITS-1)/UINT_BITS)*nstates;
    UINT sum_end_node = (dad_branch->partial_pars[scoreid] + node_branch->partial_pars[scoreid]);
    UINT score = sum_end_node;

    UINT lower_bound = best_pars_score;
    if (branch_subst) lower_bound = INT_MAX;
    UINT *x = dad_branch->partial_pars;
    UINT *y = node_branch->partial_pars;
    switch (nstates) {
    case 4:
		for (site = 1; site <= nsites; site++) {
			UINT w = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
			w = ~w;
			score += vml_popcnt(w);
            if (score >= lower_bound)
                break;
			x += 4;
			y += 4;
		}
		break;
    default:
		for (site = 0; site < nsites; site+=UINT_BITS) {
			int i;
			UINT w = x[0] & y[0];
			for (i = 1; i < nstates; i++) {
				w |= x[i] & y[i];
			}
			w = ~w;
			score += vml_popcnt(w);
            if (score >= lower_bound)
                break;
			x += nstates;
			y += nstates;

		}
		break;
    }
    if (branch_subst)
        *branch_subst = score - sum_end_node;
//    score += sum_end_node;
    return score;
}




