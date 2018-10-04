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

double PhyloTree::JukesCantorCorrection(double dist, double alpha) {
    double z = (double) aln->num_states / (aln->num_states - 1);
    double x = 1.0 - (z * dist);
    if (x > 0) {
        if (alpha <= 0.0) {
            dist = -log(x) / z;
        } else {
            //if (verbose_mode >= VB_MAX) cout << "alpha: " << alpha << endl;
            dist = alpha * (pow(x, -1.0/alpha) - 1) / z;
        }
    }
    // Branch lengths under PoMo are #events, which is ~N^2 * #substitutions
    if (aln->seq_type == SEQ_POMO)
        dist *= aln->virtual_pop_size * aln->virtual_pop_size;
    if (dist < Params::getInstance().min_branch_length)
        dist = Params::getInstance().min_branch_length;
    return dist;
}

int PhyloTree::setParsimonyBranchLengths() {
    
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    clearAllPartialLH();
    
    int sum_score = 0;
    double persite = 1.0/getAlnNSite();
    double alpha = (site_rate) ? site_rate->getGammaShape() : 1.0;
//    int pars_score;
    //int i, state;

    PhyloNeighbor *dad_branch = (PhyloNeighbor*)nodes1[0]->findNeighbor(nodes2[0]);
    PhyloNeighbor *node_branch = (PhyloNeighbor*)nodes2[0]->findNeighbor(nodes1[0]);
    PhyloNode *dad =  (PhyloNode*) nodes1[0];
    PhyloNode *node =  (PhyloNode*) nodes2[0];

    // determine state of the root
    int branch_subst = 0;
    int pars_score = computeParsimonyBranchFast(dad_branch, dad, &branch_subst);
    int site, real_site;
    int nsites = (aln->num_parsimony_sites + UINT_BITS-1) / UINT_BITS;
    int nstates = aln->getMaxNumStates();

    vector<vector<StateType> > sequences;
    sequences.resize(nodeNum, vector<StateType>(aln->num_parsimony_sites, aln->STATE_UNKNOWN));
    vector<bool> done;
    done.resize(nodeNum, false);
    done[node->id] = done[dad->id] = true;

    int subst = 0;
    
    for (site = 0, real_site = 0; site < nsites; site++) {
        size_t offset = nstates*site;
        UINT *x = dad_branch->partial_pars + offset;
        UINT *y = node_branch->partial_pars + offset;
        UINT w = x[0] & y[0];
        int state;
        for (state = 1; state < nstates; state++) {
            w |= x[state] & y[state];
        }
        UINT bit = 1;
        for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites; s++, bit = bit << 1, real_site++)
        if (w & bit) {
            // intersection is non-empty
            for (state = 0; state < nstates; state++)
                if ((x[state] & bit) && (y[state] & bit)) {
                    // assign the first state in the intersection
                    sequences[node->id][real_site] = sequences[dad->id][real_site] = state;
                    break;
                }
        } else {
            // intersection is empty
            subst++;
            for (state = 0; state < nstates; state++)
                if (x[state] & bit) {
                    // assign the first admissible state
                    sequences[node->id][real_site] = state;
                    break;
                }
            for (state = 0; state < nstates; state++)
                if (y[state] & bit) {
                    // assign the first admissible state
                    sequences[dad->id][real_site] = state;
                    break;
                }
        }
    }

    ASSERT(subst == branch_subst);
    sum_score += subst;
    fixOneNegativeBranch(correctBranchLengthF81(subst*persite, alpha), dad_branch, dad);
    
    
    // walking down the tree to assign node states
    for (int id = 1; id < nodes1.size(); id++) {
        // arrange such that states of dad are known
        if (done[nodes1[id]->id]) {
            dad = (PhyloNode*)nodes1[id];
            node = (PhyloNode*)nodes2[id];
        } else {
            ASSERT(done[nodes2[id]->id]);
            dad = (PhyloNode*)nodes2[id];
            node = (PhyloNode*)nodes1[id];
        }
        done[node->id] = true;
        // now determine states of node
        dad_branch = (PhyloNeighbor*)dad->findNeighbor(node);
        node_branch = (PhyloNeighbor*)node->findNeighbor(dad);
        subst = 0;
        for (site = 0, real_site = 0; site < nsites; site++) {
            size_t offset = nstates*site;
            UINT *x = dad_branch->partial_pars + offset;
            //UINT *y = node_branch->partial_pars + offset;
            int state;
            UINT bit = 1;
            for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites; s++, bit = bit << 1, real_site++) {
                StateType dad_state = sequences[dad->id][real_site];
                ASSERT(dad_state < nstates);
                //ASSERT(y[dad_state] & bit);
                if (x[dad_state] & bit) {
                    // same state as dad
                    sequences[node->id][real_site] = dad_state;
                } else {
                    // different state from dad
                    subst++;
                    for (state = 0; state < nstates; state++)
                        if (x[state] & bit) {
                            // assign the first admissible state
                            sequences[node->id][real_site] = state;
                            break;
                        }
                }
            }
        }
        fixOneNegativeBranch(correctBranchLengthF81(subst*persite, alpha), dad_branch, dad);
//        computeParsimonyBranchFast(dad_branch, dad, &branch_subst);
//        ASSERT(subst <= branch_subst);
        sum_score += subst;
    }
    
    ASSERT(pars_score == sum_score);
    
    return nodes1.size();
}


/****************************************************************************
 Sankoff parsimony function
 ****************************************************************************/


void PhyloTree::initCostMatrix(CostMatrixType cost_type) {
    if(cost_matrix){
        aligned_free(cost_matrix);
        cost_matrix = NULL;
    }
    ASSERT(aln);
    int cost_nstates = aln->num_states;
    // allocate memory for cost_matrix
    cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);
    
    switch (cost_type) {
        case CM_LINEAR:
            for(int i = 0; i < cost_nstates; i++){
                for(int j = 0; j < cost_nstates; j++)
                    cost_matrix[i * cost_nstates + j] = abs(i-j);
            }
            break;
        case CM_UNIFORM:
            for(int i = 0; i < cost_nstates; i++){
                for(int j = 0; j < cost_nstates; j++)
                    cost_matrix[i * cost_nstates + j] = ((i==j) ? 0 : 1);
            }
            break;
    }
    clearAllPartialLH();
}

void PhyloTree::loadCostMatrixFile(char * file_name){
    if(cost_matrix){
        aligned_free(cost_matrix);
        cost_matrix = NULL;
    }
    //    if(strcmp(file_name, "fitch") == 0)
    ////    if(file_name == NULL)
    //        cost_matrix = new SankoffCostMatrix(aln->num_states);
    //    else
    //        cost_matrix = new SankoffCostMatrix(file_name);
    
    int cost_nstates;
    if(strcmp(file_name, "fitch") == 0 || strcmp(file_name, "e") == 0) { // uniform cost
        cost_nstates = aln->num_states;
        cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);
        for(int i = 0; i < cost_nstates; i++)
            for(int j = 0; j < cost_nstates; j++){
                if(j == i) cost_matrix[i * cost_nstates + j] = 0;
                else cost_matrix[i * cost_nstates + j] = 1;
            }
    } else{ // Sankoff cost
        cout << "Loading cost matrix from " << file_name << "..." << endl;
        ifstream fin(file_name);
        if(!fin.is_open()){
            outError("Reading cost matrix file cannot perform. Please check your input file!");
        }
        fin >> cost_nstates;
        if (cost_nstates != aln->num_states)
            outError("Cost matrix file does not have the same size as number of alignment states");
        // allocate memory for cost_matrix
        cost_matrix = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);
        
        // read numbers from file
        for(int i = 0; i < cost_nstates; i++){
            for(int j = 0; j < cost_nstates; j++)
                fin >> cost_matrix[i * cost_nstates + j];
        }
        
        fin.close();
        
    }
    
    int i, j, k;
    bool changed = false;
    
    for (k = 0; k < cost_nstates; k++)
        for (i = 0; i < cost_nstates; i++)
            for (j = 0; j < cost_nstates; j++)
                if (cost_matrix[i*cost_nstates+j] > cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j]) {
                    changed = true;
                    cost_matrix[i*cost_nstates+j] = cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j];
                }
    
    if (changed) {
        cout << "WARING: Cost matrix does not satisfy triangular inenquality and is automatically fixed to:" << endl;
        cout << cost_nstates << endl;
        for (i = 0; i < cost_nstates; i++) {
            for (j = 0; j < cost_nstates; j++)
                cout << "  " << cost_matrix[i*cost_nstates+j];
            cout << endl;
        }
    } else {
        cout << "Cost matrix satisfies triangular inenquality" << endl;
    }
}

void PhyloTree::computeTipPartialParsimony() {
    if ((tip_partial_lh_computed & 2) != 0)
        return;
    tip_partial_lh_computed |= 2;
    
    int i, state, nstates = aln->num_states;
    ASSERT(tip_partial_pars);
    // ambiguous characters
    int ambi_aa[] = {
        4+8, // B = N or D
        32+64, // Z = Q or E
        512+1024 // U = I or L
    };
    
    memset(tip_partial_pars, 0, (aln->STATE_UNKNOWN+1)*nstates*sizeof(UINT));
    
    // initialize real states with cost_matrix
    memcpy(tip_partial_pars, cost_matrix, nstates*nstates*sizeof(UINT));

    UINT *this_tip_partial_pars;
    
    switch (aln->seq_type) {
        case SEQ_DNA:
            for (state = 4; state < 18; state++) {
                int cstate = state-nstates+1;
                this_tip_partial_pars = &tip_partial_pars[state*nstates];
                for (i = 0; i < nstates; i++) {
                    if ((cstate) & (1 << i))
                        this_tip_partial_pars[i] = 0;
                    else {
                        this_tip_partial_pars[i] = UINT_MAX;
                        for (int j = 0; j < nstates; j++)
                            if ((cstate) & (1 << j))
                                this_tip_partial_pars[i] = min(this_tip_partial_pars[i], cost_matrix[i*nstates+j]);
                    }
                }
            }
            break;
        case SEQ_PROTEIN:
            for (state = 0; state < sizeof(ambi_aa)/sizeof(int); state++) {
                this_tip_partial_pars = &tip_partial_pars[(state+20)*nstates];
                for (i = 0; i < nstates; i++) {
                    if (ambi_aa[state] & (1 << i))
                        this_tip_partial_pars[i] = 0;
                    else {
                        this_tip_partial_pars[i] = UINT_MAX;
                        for (int j = 0; j < nstates; j++)
                            if (ambi_aa[state] & (1 << j))
                                this_tip_partial_pars[i] = min(this_tip_partial_pars[i], cost_matrix[i*nstates+j]);
                    }
                }
            }
            break;
        case SEQ_POMO:
            ASSERT(0 && "POMO not handled with Sankoff parsimony");
            break;
        default:
            break;
    }
}
/**
 compute partial parsimony score of the subtree rooted at dad
 @param dad_branch the branch leading to the subtree
 @param dad its dad, used to direct the traversal
 */
void PhyloTree::computePartialParsimonySankoff(PhyloNeighbor *dad_branch, PhyloNode *dad){
    // don't recompute the parsimony
    if (dad_branch->partial_lh_computed & 2)
        return;
    
    Node *node = dad_branch->node;
    //assert(node->degree() <= 3);
    /*
     if(aln->num_states != cost_nstates){
     cout << "Your cost matrix is not compatible with the alignment"
     << " in terms of number of states. Please check!" << endl;
     exit(1);
     }
     */
    int nstates = aln->num_states;
    assert(dad_branch->partial_pars);
    
    int pars_block_size = getBitsBlockSize();
    

    // internal node
    UINT i, j, ptn, min_child_ptn_pars;
    
    UINT * partial_pars = dad_branch->partial_pars;
    memset(partial_pars, 0, sizeof(UINT)*pars_block_size);

    PhyloNeighbor *left = NULL, *right = NULL;
    
    FOR_NEIGHBOR_IT(node, dad, it)
        if ((*it)->node->name != ROOT_NAME) {
            if (!(*it)->node->isLeaf())
                computePartialParsimonySankoff((PhyloNeighbor*) (*it), (PhyloNode*) node);
            if (!left)
                left = ((PhyloNeighbor*)*it);
            else
                right = ((PhyloNeighbor*)*it);
        }
    
    if (!left->node->isLeaf() && right->node->isLeaf()) {
        // swap leaf and internal node
        PhyloNeighbor *tmp = left;
        left = right;
        right = tmp;
    }
    ASSERT(node->degree() >= 3);
    
    if (node->degree() > 3) {
        // multifurcating node
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++) {
            int ptn_start_index = ptn*nstates;
            UINT *partial_pars_ptr = &partial_pars[ptn_start_index];

            FOR_NEIGHBOR_IT(node, dad, it) if ((*it)->node->name != ROOT_NAME) {
                if ((*it)->node->isLeaf()) {
                    // leaf node
                    UINT *partial_pars_child_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][(*it)->node->id]*nstates];
                
                    for(i = 0; i < nstates; i++){
                        partial_pars_ptr[i] += partial_pars_child_ptr[i];
                    }
                } else {
                    // internal node
                    UINT *partial_pars_child_ptr = &((PhyloNeighbor*) (*it))->partial_pars[ptn_start_index];
                    UINT *cost_matrix_ptr = cost_matrix;
                    
                    for (i = 0; i < nstates; i++){
                        // min(j->i) from child_branch
                        min_child_ptn_pars = partial_pars_child_ptr[0] + cost_matrix_ptr[0];
                        for(j = 1; j < nstates; j++) {
                            UINT value = partial_pars_child_ptr[j] + cost_matrix_ptr[j];
                            min_child_ptn_pars = min(value, min_child_ptn_pars);
                        }
                        partial_pars_ptr[i] += min_child_ptn_pars;
                        cost_matrix_ptr += nstates;
                    }
                }
            }
        }
    } else if (left->node->isLeaf() && right->node->isLeaf()) {
        // tip-tip case
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            int ptn_start_index = ptn*nstates;
            
            UINT *left_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            UINT *right_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][right->node->id]*nstates];
            UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
            
            for (i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                partial_pars_ptr[i] = left_ptr[i] + right_ptr[i];
            }
        }
    } else if (left->node->isLeaf() && !right->node->isLeaf()) {
        // tip-inner case
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            int ptn_start_index = ptn*nstates;
            
            UINT *left_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            UINT *right_ptr = &right->partial_pars[ptn_start_index];
            UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
            UINT *cost_matrix_ptr = cost_matrix;
            UINT right_contrib;
            
            for(i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                for(j = 1; j < nstates; j++) {
                    right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
                }
                partial_pars_ptr[i] = left_ptr[i] + right_contrib;
                cost_matrix_ptr += nstates;
            }
        }
    } else {
        // inner-inner case
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            int ptn_start_index = ptn*nstates;
            
            UINT *left_ptr = &left->partial_pars[ptn_start_index];
            UINT *right_ptr = &right->partial_pars[ptn_start_index];
            UINT *partial_pars_ptr = &partial_pars[ptn_start_index];
            UINT *cost_matrix_ptr = cost_matrix;
            UINT left_contrib, right_contrib;
            
            for(i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                left_contrib = left_ptr[0] + cost_matrix_ptr[0];
                right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                for(j = 1; j < nstates; j++) {
                    left_contrib = min(left_ptr[j] + cost_matrix_ptr[j], left_contrib);
                    right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
                }
                partial_pars_ptr[i] = left_contrib+right_contrib;
                cost_matrix_ptr += nstates;
            }
        }
        
    }
    
    dad_branch->partial_lh_computed |= 2;
}

/**
 compute tree parsimony score based on a particular branch
 @param dad_branch the branch leading to the subtree
 @param dad its dad, used to direct the traversal
 @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
 @return parsimony score of the tree
 */
int PhyloTree::computeParsimonyBranchSankoff(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {

    if ((tip_partial_lh_computed & 2) == 0)
        computeTipPartialParsimony();

    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    
    if (!central_partial_pars)
        initializeAllPartialPars();
    
    // DTH: I don't really understand what this is for. ###########
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //        cout << "swapped\n";
    }
    
    //int nptn = aln->size();
    //    if(!_pattern_pars) _pattern_pars = aligned_alloc<BootValTypePars>(nptn+VCSIZE_USHORT);
    //    memset(_pattern_pars, 0, sizeof(BootValTypePars) * (nptn+VCSIZE_USHORT));
    
    if ((dad_branch->partial_lh_computed & 2) == 0 && !node->isLeaf())
        computePartialParsimonySankoff(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0 && !dad->isLeaf())
        computePartialParsimonySankoff(node_branch, node);
    
    // now combine likelihood at the branch
    UINT tree_pars = 0;
    int nstates = aln->num_states;
    UINT i, j, ptn;
    UINT branch_pars = 0;
    
    if (dad->isLeaf()) {
        // external node
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            int ptn_start_index = ptn * nstates;
            UINT *node_branch_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][dad->id]*nstates];
            UINT *dad_branch_ptr = &dad_branch->partial_pars[ptn_start_index];
            UINT min_ptn_pars = node_branch_ptr[0] + dad_branch_ptr[0];
            UINT br_ptn_pars = node_branch_ptr[0];
            for (i = 1; i < nstates; i++){
                // min(j->i) from node_branch
                UINT min_score = node_branch_ptr[i] + dad_branch_ptr[i];
                if (min_score < min_ptn_pars) {
                    min_ptn_pars = min_score;
                    br_ptn_pars = node_branch_ptr[i];
                }
            }
            //_pattern_pars[ptn] = min_ptn_pars;
            tree_pars += min_ptn_pars * aln->ordered_pattern[ptn].frequency;
            branch_pars += br_ptn_pars * aln->ordered_pattern[ptn].frequency;
        }
    }  else {
        // internal node
        for (ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            int ptn_start_index = ptn * nstates;
            UINT *node_branch_ptr = &node_branch->partial_pars[ptn_start_index];
            UINT *dad_branch_ptr = &dad_branch->partial_pars[ptn_start_index];
            UINT *cost_matrix_ptr = cost_matrix;
            UINT min_ptn_pars = UINT_MAX;
            UINT br_ptn_pars = UINT_MAX;
            for(i = 0; i < nstates; i++){
                // min(j->i) from node_branch
                UINT min_score = node_branch_ptr[0] + cost_matrix_ptr[0];
                UINT branch_score = cost_matrix_ptr[0];
                for(j = 1; j < nstates; j++) {
                    UINT value = node_branch_ptr[j] + cost_matrix_ptr[j];
                    if (value < min_score) {
                        min_score = value;
                        branch_score = cost_matrix_ptr[j];
                    }
                    
                }
                min_score = min_score + dad_branch_ptr[i];
                if (min_score < min_ptn_pars) {
                    min_ptn_pars = min_score;
                    br_ptn_pars = branch_score;
                }
                cost_matrix_ptr += nstates;
            }
            //_pattern_pars[ptn] = min_ptn_pars;
            tree_pars += min_ptn_pars * aln->ordered_pattern[ptn].frequency;
            branch_pars += br_ptn_pars * aln->ordered_pattern[ptn].frequency;
        }
    }
    if (branch_subst)
        *branch_subst = branch_pars;
    //    cout << endl;
    return tree_pars;
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
    if (verbose_mode >= VB_MAX)
        cout << "computeParsimony: " << computeParsimony() << endl;
    size_t index = (2*leafNum-3)*2;
    size_t pars_block_size = getBitsBlockSize();

    UINT *tmp_partial_pars;
    tmp_partial_pars = newBitsBlock();
    if (size == 3)
        best_pars_score = computeParsimony();

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
