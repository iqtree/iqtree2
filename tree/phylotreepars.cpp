/*
 * phylotreepars.cpp
 *
 * Fast implementation of parsimony kernel
 *
 *  Created on: May 18, 2015
 *      Author: minh
 */

#include "phylotree.h"
//#include "vectorclass/vectorclass.h"
#include "phylosupertree.h"

#if defined (__GNUC__) || defined(__clang__)
#define vml_popcnt __builtin_popcount
#else
// taken from vectorclass library
static inline uint32_t vml_popcnt (uint32_t a) {	
    // popcnt instruction not available
    uint32_t b = a - ((a >> 1) & 0x55555555);
    uint32_t c = (b & 0x33333333) + ((b >> 2) & 0x33333333);
    uint32_t d = (c + (c >> 4)) & 0x0F0F0F0F;
    uint32_t e = d * 0x01010101;
    return   e >> 24;
}
#endif

/***********************************************************/
/****** optimized version of parsimony kernel **************/
/***********************************************************/

void PhyloTree::computePartialParsimonyFast(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (dad_branch->partial_lh_computed & 2)
        return;
    Node *node = dad_branch->node;
    int nstates = aln->getMaxNumStates();
    int site = 0;

    dad_branch->partial_lh_computed |= 2;

    vector<Alignment*> *partitions = NULL;
    if (aln->isSuperAlignment())
        partitions = &((SuperAlignment*)aln)->partitions;
    else {
        partitions = new vector<Alignment*>;
        partitions->push_back(aln);
    }

    if (node->isLeaf() && dad) {
        // external node
        int leafid = node->id;
        memset(dad_branch->partial_pars, 0, getBitsBlockSize()*sizeof(UINT));
        int max_sites = ((aln->num_parsimony_sites+UINT_BITS-1)/UINT_BITS)*UINT_BITS;
        int ambi_aa[] = {2, 3, 5, 6, 9, 10}; // {4+8, 32+64, 512+1024};
//        if (aln->ordered_pattern.empty())
//            aln->orderPatternByNumChars();
        ASSERT(!aln->ordered_pattern.empty());
        int start_pos = 0;
        for (vector<Alignment*>::iterator alnit = partitions->begin(); alnit != partitions->end(); alnit++) {
            int end_pos = start_pos + (*alnit)->ordered_pattern.size();
            switch ((*alnit)->seq_type) {
            case SEQ_DNA:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (state < 4) {
                        for (int j = 0; j < freq; j++, site++) {
                            dad_branch->partial_pars[(site/UINT_BITS)*nstates+state] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                            UINT bit1 = (1 << (site%UINT_BITS));
                            p[0] |= bit1;
                            p[1] |= bit1;
                            p[2] |= bit1;
                            p[3] |= bit1;
                        }
                    } else {
                        state -= 3;
                        ASSERT(state < 15);
                        for (int j = 0; j < freq; j++, site++) {
                            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                            UINT bit1 = (1 << (site%UINT_BITS));
                            for (int i = 0; i < 4; i++)
                                if (state & (1<<i))
                                    p[i] |= bit1;
                        }
                    }
                }
                //assert(site == aln->num_informative_sites);
                // add dummy states
                //if (site < max_sites)
                //    dad_branch->partial_pars[(site/UINT_BITS)*4] |= ~((1<<(site%UINT_BITS)) - 1);
                break;
            case SEQ_PROTEIN:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (state < 20) {
                        for (int j = 0; j < freq; j++, site++) {
                            dad_branch->partial_pars[(site/UINT_BITS)*nstates+state] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                            UINT bit1 = (1 << (site%UINT_BITS));
                            for (int i = 0; i < 20; i++)
                                    p[i] |= bit1;
                        }
                    } else {
                        ASSERT(state < 23);
                        state = (state-20)*2;
                        for (int j = 0; j < freq; j++, site++) {
                            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                            UINT bit1 = (1 << (site%UINT_BITS));
                            p[ambi_aa[state]] |= bit1;
                            p[ambi_aa[state+1]] |= bit1;
                        }
                    }
                }
                //assert(site == aln->num_informative_sites);
                // add dummy states
                //if (site < max_sites)
                //    dad_branch->partial_pars[(site/UINT_BITS)*20] |= ~((1<<(site%UINT_BITS)) - 1);
                break;
            default:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (aln->seq_type == SEQ_POMO && state >= (*alnit)->num_states && state < (*alnit)->STATE_UNKNOWN) {
                        state = (*alnit)->convertPomoState(state);
                    }
                     if (state < (*alnit)->num_states) {
                        for (int j = 0; j < freq; j++, site++) {
                            dad_branch->partial_pars[(site/UINT_BITS)*nstates+state] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            UINT *p = dad_branch->partial_pars+((site/UINT_BITS)*nstates);
                            UINT bit1 = (1 << (site%UINT_BITS));
                            for (int i = 0; i < (*alnit)->num_states; i++)
                                    p[i] |= bit1;
                        }
                    } else {
                        ASSERT(0);
                    }
                }
                break;
            } // end of switch
            
            start_pos = end_pos;
        } // FOR LOOP

        ASSERT(site == aln->num_parsimony_sites);
        // add dummy states
        if (site < max_sites)
            dad_branch->partial_pars[(site/UINT_BITS)*nstates] |= ~((1<<(site%UINT_BITS)) - 1);
    } else {
        // internal node
        ASSERT(node->degree() == 3); // it works only for strictly bifurcating tree
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
        int nsites = aln->num_parsimony_sites;
        nsites = (nsites+UINT_BITS-1)/UINT_BITS;

        switch (nstates) {
        case 4:
            #ifdef _OPENMP
            #pragma omp parallel for private (site) reduction(+: score) if(nsites>200)
            #endif
			for (site = 0; site<nsites; site++) {
				UINT w;
                size_t offset = nstates*site;
                UINT *x = left->partial_pars + offset;
                UINT *y = right->partial_pars + offset;
                UINT *z = dad_branch->partial_pars + offset;
				z[0] = x[0] & y[0];
				z[1] = x[1] & y[1];
				z[2] = x[2] & y[2];
				z[3] = x[3] & y[3];
				w = z[0] | z[1] | z[2] | z[3];
				w = ~w;
				score += __builtin_popcount(w);
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
    
    if (!aln->isSuperAlignment())
        delete partitions;
}


int PhyloTree::computeParsimonyBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    ASSERT(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFast(node_branch, node);
    int site;
    int nsites = (aln->num_parsimony_sites + UINT_BITS-1) / UINT_BITS;
    int nstates = aln->getMaxNumStates();

    int scoreid = ((aln->num_parsimony_sites+UINT_BITS-1)/UINT_BITS)*nstates;
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
//            #ifndef _OPENMP
//            if (score >= lower_bound)
//                break;
//            #endif
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
//            #ifndef _OPENMP
//            if (score >= lower_bound)
//                break;
//            #endif
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

    IntVector taxon_order;
    taxon_order.reserve(size);

    if (constraintTree.empty()) {
        freeNode();
        taxon_order.resize(size);
        for (int i = 0; i < size; i++)
            taxon_order[i] = i;
        // randomize the addition order
        my_random_shuffle(taxon_order.begin(), taxon_order.end());

        root = newNode(size);

        // create initial tree with 3 taxa
        for (leafNum = 0; leafNum < 3; leafNum++) {
            if (verbose_mode >= VB_MAX)
                cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree" << endl;
            Node *new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
            root->addNeighbor(new_taxon, -1.0);
            new_taxon->addNeighbor(root, -1.0);
        }
    } else {
        // first copy the constraint tree
        MTree::copyTree(&constraintTree);
        
        // convert to birfucating tree if needed
        extractBifurcatingSubTree();
        ASSERT(isBifurcating());
        
        // assign proper taxon IDs
        NodeVector nodes;
        NodeVector::iterator it;
        getTaxa(nodes);
        leafNum = nodes.size();
        vector<int> pushed;
        pushed.resize(size, 0);
        for (it = nodes.begin(); it != nodes.end(); it++) {
            (*it)->id = aln->getSeqID((*it)->name);
            ASSERT((*it)->id >= 0);
            taxon_order.push_back((*it)->id);
            pushed[(*it)->id] = 1;
        }

        // start with constraint tree
        int i;
        for (i = 0; i < size; i++)
            if (!pushed[i] && constraintTree.hasTaxon(aln->getSeqName(i))) {
                taxon_order.push_back(i);
                pushed[i] = 1;
            }
        ASSERT(taxon_order.size() == constraintTree.leafNum);
        for (int i = 0; i < size; i++)
            if (!pushed[i]) {
                taxon_order.push_back(i);
            }
        // randomize the addition order
        my_random_shuffle(taxon_order.begin()+leafNum, taxon_order.begin()+constraintTree.leafNum);
        my_random_shuffle(taxon_order.begin()+constraintTree.leafNum, taxon_order.end());

    }
    root = findNodeID(taxon_order[0]);
    initializeAllPartialPars();
    size_t index = (2*leafNum-3)*2;
    size_t pars_block_size = getBitsBlockSize();

    UINT *tmp_partial_pars;
    tmp_partial_pars = newBitsBlock();

    best_pars_score = 0;
    if (leafNum == size) {
        outWarning("Constraint tree has all taxa and is bifurcating, which strictly enforces final tree!");
    }
    
    // stepwise adding the next taxon for the remaining taxa
    for (; leafNum < size; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Adding " << aln->getSeqName(taxon_order[leafNum]) << " to the tree..." << endl;
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
    
    ASSERT(index == 4*leafNum-6);

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
//    if (isSuperTree())
//        ((PhyloSuperTree*)this)->mapTrees();
    
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
    if (leafNum < constraintTree.leafNum) {
        // still during addition of taxa from constraint tree
        if (!constraintTree.isCompatible(this))
            score = INT_MAX;
    }
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
