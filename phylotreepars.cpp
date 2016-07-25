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
#include "phylosupertree.h"

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
        nsites = (nsites+UINT_BITS-1)/UINT_BITS;

        switch (nstates) {
        case 4:
            #ifdef _OPENMP
            #pragma omp parallel for private (site) reduction(+: score) if(nsites>200)
            #endif
			for (site = 0; site<nsites; site++) {
				UINT w;
                size_t offset = 4*site;
                UINT *x = left->partial_pars + offset;
                UINT *y = right->partial_pars + offset;
                UINT *z = dad_branch->partial_pars + offset;
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
			}
			break;
        default:
            #ifdef _OPENMP
            #pragma omp parallel for private (site) reduction(+: score) if(nsites > 800/nstates)
            #endif
			for (site = 0; site<nsites; site++) {
				int i;
				UINT w = 0;
                size_t offset = nstates*site;
                UINT *x = left->partial_pars + offset;
                UINT *y = right->partial_pars + offset;
                UINT *z = dad_branch->partial_pars + offset;
                
				for (i = 0; i < nstates; i++) {
					z[i] = x[i] & y[i];
					w |= z[i];
				}
				w = ~w;
				score += vml_popcnt(w);
				for (i = 0; i < nstates; i++) {
					z[i] |= w & (x[i] | y[i]);
				}
			}
			break;
        }
        dad_branch->partial_pars[nstates*nsites] = score + left->partial_pars[nstates*nsites] + right->partial_pars[nstates*nsites];
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
    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for private (site) reduction(+: score) if(nsites>200)
        #endif
		for (site = 0; site < nsites; site++) {
            size_t offset = 4*site;
            UINT *x = dad_branch->partial_pars + offset;
            UINT *y = node_branch->partial_pars + offset;
			UINT w = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
			w = ~w;
			score += vml_popcnt(w);
            #ifndef _OPENMP
            if (score >= lower_bound)
                break;
            #endif
		}
		break;
    default:
        #ifdef _OPENMP
        #pragma omp parallel for private (site) reduction(+: score) if(nsites > 800/nstates)
        #endif
		for (site = 0; site < nsites; site++) {
            size_t offset = nstates*site;
            UINT *x = dad_branch->partial_pars + offset;
            UINT *y = node_branch->partial_pars + offset;
			int i;
			UINT w = x[0] & y[0];
			for (i = 1; i < nstates; i++) {
				w |= x[i] & y[i];
			}
			w = ~w;
			score += vml_popcnt(w);
            #ifndef _OPENMP
            if (score >= lower_bound)
                break;
            #endif
		}
		break;
    }
    if (branch_subst)
        *branch_subst = score - sum_end_node;
//    score += sum_end_node;
    return score;
}

void PhyloTree::computeAllPartialPars(PhyloNode *node, PhyloNode *dad) {
	if (!node) node = (PhyloNode*)root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		if ((((PhyloNeighbor*)*it)->partial_lh_computed & 1) == 0)
			computePartialParsimony((PhyloNeighbor*)*it, node);
		PhyloNeighbor *rev = (PhyloNeighbor*) (*it)->node->findNeighbor(node);
		if ((rev->partial_lh_computed & 1) == 0)
			computePartialParsimony(rev, (PhyloNode*)(*it)->node);
		computeAllPartialPars((PhyloNode*)(*it)->node, node);
	}
}


/****************************************************************************
 Stepwise addition (greedy) by maximum parsimony
 ****************************************************************************/

// random generator function:
//ptrdiff_t myrandom(ptrdiff_t i) {
//    return random_int(i);
//}

// pointer object to it:
//ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

int PhyloTree::computeParsimonyTree(const char *out_prefix, Alignment *alignment) {
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    freeNode();

    root = newNode(size);

    IntVector taxon_order;
    taxon_order.resize(size);
    for (int i = 0; i < size; i++)
        taxon_order[i] = i;
    // randomize the addition order
    my_random_shuffle(taxon_order.begin(), taxon_order.end());

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree" << endl;
        Node *new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
        root->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(root, -1.0);
    }
    root = findNodeID(taxon_order[0]);
    initializeAllPartialPars();
    size_t index = 6;
    size_t pars_block_size = getBitsBlockSize();

    if (isSuperTree())
        ((PhyloSuperTree*)this)->mapTrees();
    
    UINT *tmp_partial_pars;
    tmp_partial_pars = newBitsBlock();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree";
        NodeVector nodes1, nodes2;
        getBranches(nodes1, nodes2);
        PhyloNode *target_node = NULL;
        PhyloNode *target_dad = NULL;
        best_pars_score = INT_MAX;
        // allocate a new taxon and a new adjacent internal node
        
        UINT *new_taxon_partial_pars = central_partial_pars + ((index++) * pars_block_size);
        
        PhyloNode *new_taxon = (PhyloNode*)newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
        PhyloNode *added_node = (PhyloNode*)newNode(size+leafNum-2);
        added_node->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(added_node, -1.0);
        ((PhyloNeighbor*) added_node->findNeighbor(new_taxon))->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        ((PhyloNeighbor*) new_taxon->findNeighbor(added_node))->partial_pars = tmp_partial_pars;

        // preserve two neighbors
        added_node->addNeighbor((Node*) 1, -1.0);
        added_node->addNeighbor((Node*) 2, -1.0);

        for (int nodeid = 0; nodeid < nodes1.size(); nodeid++) {
        
            // check for compatible with constraint tree
            if (!constraintTree.empty()) {
                StrVector taxset1, taxset2;
                getUnorderedTaxaName(taxset1, nodes1[nodeid], nodes2[nodeid]);
                getUnorderedTaxaName(taxset2, nodes2[nodeid], nodes1[nodeid]);
                taxset1.push_back(new_taxon->name);
                if (!constraintTree.isCompatible(taxset1, taxset2))
                    continue;
                taxset1.pop_back();
                taxset2.push_back(new_taxon->name);
                if (!constraintTree.isCompatible(taxset1, taxset2))
                    continue;                
            }

            int score = addTaxonMPFast(new_taxon, added_node, nodes1[nodeid], nodes2[nodeid]);
            if (score < best_pars_score) {
                best_pars_score = score;
                target_node = (PhyloNode*)nodes1[nodeid];
                target_dad = (PhyloNode*)nodes2[nodeid];
                memcpy(new_taxon_partial_pars, tmp_partial_pars, pars_block_size*sizeof(UINT));
            }
        }
        
        if (verbose_mode >= VB_MAX)
            cout << ", score = " << best_pars_score << endl;
        // now insert the new node in the middle of the branch node-dad
        target_node->updateNeighbor(target_dad, added_node, -1.0);
        target_dad->updateNeighbor(target_node, added_node, -1.0);
        added_node->updateNeighbor((Node*) 1, target_node, -1.0);
        added_node->updateNeighbor((Node*) 2, target_dad, -1.0);
        ((PhyloNeighbor*) added_node->findNeighbor(target_node))->partial_pars =
            ((PhyloNeighbor*) target_dad->findNeighbor(added_node))->partial_pars;
        ((PhyloNeighbor*) added_node->findNeighbor(target_dad))->partial_pars =
            ((PhyloNeighbor*) target_node->findNeighbor(added_node))->partial_pars;
            
        ((PhyloNeighbor*) added_node->findNeighbor(target_node))->partial_lh_computed = 
            ((PhyloNeighbor*) target_dad->findNeighbor(added_node))->partial_lh_computed;
        ((PhyloNeighbor*) added_node->findNeighbor(target_dad))->partial_lh_computed = 
            ((PhyloNeighbor*) target_node->findNeighbor(added_node))->partial_lh_computed;
        
        ((PhyloNeighbor*) new_taxon->findNeighbor(added_node))->partial_lh_computed |= 2;
        ((PhyloNeighbor*) new_taxon->findNeighbor(added_node))->partial_pars = new_taxon_partial_pars;

        ((PhyloNeighbor*)target_dad->findNeighbor(added_node))->clearPartialLh();
        ((PhyloNeighbor*)target_dad->findNeighbor(added_node))->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        ((PhyloNeighbor*)target_node->findNeighbor(added_node))->clearPartialLh();
        ((PhyloNeighbor*)target_node->findNeighbor(added_node))->partial_pars = central_partial_pars + ((index++) * pars_block_size);

        target_dad->clearReversePartialLh(added_node);
        target_node->clearReversePartialLh(added_node);

    }

    aligned_free(tmp_partial_pars);
    
    assert(index == 4*leafNum-6);

    nodeNum = 2 * leafNum - 2;
    initializeTree();

    setAlignment(alignment);
//    initializeAllPartialPars();
//    clearAllPartialLH();
    fixNegativeBranch(true);
    if (out_prefix) {
		string file_name = out_prefix;
		file_name += ".parstree";
		printTree(file_name.c_str(), WT_NEWLINE);
    }
    return best_pars_score;
}

int PhyloTree::addTaxonMPFast(Node *added_taxon, Node* added_node, Node* node, Node* dad) {
    Neighbor *dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*) 1, node, len / 2.0);
    added_node->updateNeighbor((Node*) 2, dad, len / 2.0);
    ((PhyloNeighbor*) added_node->findNeighbor(node))->partial_pars =
        ((PhyloNeighbor*) dad->findNeighbor(added_node))->partial_pars;
    ((PhyloNeighbor*) added_node->findNeighbor(dad))->partial_pars =
        ((PhyloNeighbor*) node->findNeighbor(added_node))->partial_pars;
    ((PhyloNeighbor*) added_node->findNeighbor(node))->partial_lh_computed = 
        ((PhyloNeighbor*) dad->findNeighbor(added_node))->partial_lh_computed;
    ((PhyloNeighbor*) added_node->findNeighbor(dad))->partial_lh_computed = 
        ((PhyloNeighbor*) node->findNeighbor(added_node))->partial_lh_computed;
    // compute the likelihood
    ((PhyloNeighbor*) added_taxon->findNeighbor(added_node))->clearPartialLh();
    int score = computeParsimonyBranch((PhyloNeighbor*) added_node->neighbors[0], (PhyloNode*) added_node);
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*) 1, len);
    added_node->updateNeighbor(dad, (Node*) 2, len);

    // set partial_pars to COMPUTED
    ((PhyloNeighbor*)node->findNeighbor(dad))->partial_lh_computed |= 2;
    ((PhyloNeighbor*)dad->findNeighbor(node))->partial_lh_computed |= 2;

    // now tranverse the tree downwards

//    FOR_NEIGHBOR_IT(node, dad, it){
//        addTaxonMPFast(added_node, target_node, target_dad, target_partial_pars, (*it)->node, node);
//    }
    return score;

}
