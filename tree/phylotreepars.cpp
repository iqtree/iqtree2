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
    if (dad_branch->isParsimonyComputed()) {
        return;
    }
    PhyloNode* node    = dad_branch->getNode();
    int        nstates = aln->getMaxNumStates();
    int        site    = 0;

    dad_branch->setParsimonyComputed(true);

    vector<Alignment*> *partitions = NULL;
    if (aln->isSuperAlignment())
        partitions = &((SuperAlignment*)aln)->partitions;
    else {
        partitions = new vector<Alignment*>;
        partitions->push_back(aln);
    }

    if (node->name == ROOT_NAME) {
        ASSERT(dad);
        int pars_size = getBitsBlockSize();
        memset(dad_branch->partial_pars, 255, pars_size*sizeof(UINT));
        size_t nsites = (aln->num_parsimony_sites+UINT_BITS-1)/UINT_BITS;
        dad_branch->partial_pars[nstates*nsites] = 0;
    } else if (node->isLeaf() && dad) {
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
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, pit) {
            if (pit->node->name != ROOT_NAME && !pit->isParsimonyComputed()) {
                computePartialParsimonyFast(pit, node);
            }
            if (!left) left = pit; else right = pit;
        }
        computePartialParsimonyOutOfTreeFast(left->partial_pars, right->partial_pars, dad_branch->partial_pars);
    }
    if (!aln->isSuperAlignment()) {
        delete partitions;
    }
}

void PhyloTree::computePartialParsimonyOutOfTreeFast(const UINT* left_partial_pars,
                                                    const UINT* right_partial_pars,
                                                    UINT* dad_partial_pars) const {

    int  nstates = aln->getMaxNumStates();
    UINT score   = 0;
    int  nsites  = aln->num_parsimony_sites;
    
    nsites = (nsites+UINT_BITS-1)/UINT_BITS;

    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>200)
        #endif
        for (int site = 0; site<nsites; ++site) {
            UINT w;
            size_t offset = nstates*site;
            const UINT* x = left_partial_pars + offset;
            const UINT* y = right_partial_pars + offset;
            UINT*       z = dad_partial_pars + offset;
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
        #pragma omp parallel for reduction(+: score) if(nsites > 800/nstates)
        #endif
        for (int site = 0; site<nsites; ++site) {
            UINT   w = 0;
            size_t offset = nstates*site;
            const UINT*  x = left_partial_pars  + offset;
            const UINT*  y = right_partial_pars + offset;
            UINT*        z = dad_partial_pars   + offset;
            
            for (int i = 0; i < nstates; i++) {
                z[i] = x[i] & y[i];
                w |= z[i];
            }
            w = ~w;
            score += vml_popcnt(w);
            for (int i = 0; i < nstates; i++) {
                z[i] |= w & (x[i] | y[i]);
            }
        }
        break;
    }
    auto total = nstates*nsites;
    dad_partial_pars[total] = score + left_partial_pars[total] + right_partial_pars[total];
}


int PhyloTree::computeParsimonyBranchFast(PhyloNeighbor *dad_branch,
                                          PhyloNode *dad, int *branch_subst) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    ASSERT(node_branch);
    if (central_partial_pars==nullptr) {
        initializeAllPartialPars();
    }
    if (!dad_branch->isParsimonyComputed()) {
        computePartialParsimonyFast(dad_branch, dad);
    }
    if (!node_branch->isParsimonyComputed()) {
        computePartialParsimonyFast(node_branch, node);
    }
    return computeParsimonyOutOfTreeFast(dad_branch->partial_pars,
                                         node_branch->partial_pars,
                                         branch_subst );
}

int PhyloTree::computeParsimonyOutOfTreeFast(const UINT* dad_partial_pars,
                                             const UINT* node_partial_pars,
                                             int*        branch_subst) const {
    int nsites = (aln->num_parsimony_sites + UINT_BITS-1) / UINT_BITS;
    int nstates = aln->getMaxNumStates();

    int scoreid = ((aln->num_parsimony_sites+UINT_BITS-1)/UINT_BITS)*nstates;
    UINT sum_end_node = (dad_partial_pars[scoreid] + node_partial_pars[scoreid]);
    UINT score = sum_end_node;

    UINT lower_bound = best_pars_score;
    if (branch_subst) lower_bound = INT_MAX;
    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>200)
        #endif
        for (int site = 0; site < nsites; ++site) {
            size_t offset = 4*site;
            const UINT *x = dad_partial_pars + offset;
            const UINT *y = node_partial_pars + offset;
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
        #pragma omp parallel for reduction(+: score) if(nsites > 800/nstates)
        #endif
        for (int site = 0; site < nsites; ++site) {
            size_t offset = nstates * site;
            const UINT *x = dad_partial_pars + offset;
            const UINT *y = node_partial_pars + offset;
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
    if (branch_subst) {
            *branch_subst = score - sum_end_node;
    }
    //    score += sum_end_node;
        return score;
}


void PhyloTree::computeAllPartialPars(PhyloNode *node, PhyloNode *dad) {
	if (!node) node = getRoot();
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        
        //Note: This was already checking whether likelihood had
        //      been computed, rather than parsimony! (James B. 14-Sep-2020).
        //      This might be a bug.  It certainly looks like one!
        if (!nei->isLikelihoodComputed()) {
            computePartialParsimony(nei, node);
        }
        PhyloNeighbor *rev = nei->getNode()->findNeighbor(node);
        
        //Note: this, too, was a "likelihood computed?" check.
        if (!rev->isLikelihoodComputed()) {
            computePartialParsimony(rev, nei->getNode());
        }
        computeAllPartialPars(nei->getNode(), node);
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
    PhyloNodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    clearAllPartialLH();
    clearAllPartialParsimony(false);
    
    int sum_score = 0;
    double persite = 1.0/getAlnNSite();
    double alpha = (site_rate) ? site_rate->getGammaShape() : 1.0;
//    int pars_score;
    //int i, state;

    PhyloNode*     dad         = nodes1[0];
    PhyloNeighbor* dad_branch  = dad->findNeighbor(nodes2[0]);
    PhyloNode*     node        = nodes2[0];
    PhyloNeighbor* node_branch = node->findNeighbor(nodes1[0]);

    // determine state of the root
    int branch_subst = 0;
    int pars_score   = computeParsimonyBranchFast(dad_branch, dad, &branch_subst);
    int nsites       = (aln->num_parsimony_sites + UINT_BITS-1) / UINT_BITS;
    int nstates      = aln->getMaxNumStates();

    vector<vector<StateType> > sequences;
    sequences.resize(nodeNum, vector<StateType>(aln->num_parsimony_sites, aln->STATE_UNKNOWN));
    BoolVector done(nodeNum, false);
    done[node->id] = done[dad->id]  = true;

    int subst = 0; //count of substitutions
    
    for (int site = 0, real_site = 0; site < nsites; site++) {
        size_t offset = nstates*site;
        UINT *x = dad_branch->partial_pars  + offset;
        UINT *y = node_branch->partial_pars + offset;
        UINT w  = x[0] & y[0];
        for (int state = 1; state < nstates; state++) {
            w |= x[state] & y[state];
        }
        UINT bit = 1;
        StateType* dadSeq  = sequences[dad->id].data();
        StateType* nodeSeq = sequences[node->id].data();

        for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites; s++, bit = bit << 1, real_site++)
        if (w & bit) {
            // intersection is non-empty
            for (int state = 0; state < nstates; state++) {
                if ((x[state] & bit) && (y[state] & bit)) {
                    // assign the first state in the intersection
                    nodeSeq[real_site] = dadSeq[real_site] = state;
                    break;
                }
            }
        } else {
            // intersection is empty
            subst++;
            for (int state = 0; state < nstates; state++) {
                if (x[state] & bit) {
                    // assign the first admissible state
                    nodeSeq[real_site] = state;
                    break;
                }
            }
            for (int state = 0; state < nstates; state++) {
                if (y[state] & bit) {
                    // assign the first admissible state
                    dadSeq[real_site] = state;
                    break;
                }
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
            dad  = nodes1[id];
            node = nodes2[id];
        } else {
            ASSERT(done[nodes2[id]->id]);
            dad  = nodes2[id];
            node = nodes1[id];
        }
        done[node->id] = true;
        // now determine states of node
        dad_branch  = dad->findNeighbor(node);
        node_branch = node->findNeighbor(dad);
        subst = 0;
        for (int site = 0, real_site = 0; site < nsites; site++) {
            size_t offset = nstates*site;
            UINT*  x      = dad_branch->partial_pars + offset;
            //UINT* y     = node_branch->partial_pars + offset;
            UINT bit      = 1;
            const StateType* dadSeq  = sequences[dad->id].data();
            StateType*       nodeSeq = sequences[node->id].data();
            for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites; s++, bit = bit << 1, real_site++) {
                StateType dad_state = dadSeq[real_site];
                ASSERT(dad_state < nstates);
                //ASSERT(y[dad_state] & bit);
                if (x[dad_state] & bit) {
                    // same state as dad
                    nodeSeq[real_site] = dad_state;
                } else {
                    // different state from dad
                    subst++;
                    for (int state = 0; state < nstates; state++) {
                        if (x[state] & bit) {
                            // assign the first admissible state
                            nodeSeq[real_site] = state;
                            break;
                        }
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
        cost_matrix = nullptr; //(though... aligned_free set it to null already).
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
    clearAllPartialParsimony(false);
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
    if ((tip_partial_lh_computed & 2) != 0) {
        return;
    }
    tip_partial_lh_computed |= 2;
    
    const int    nstates = aln->num_states;
    const size_t nptn    = aln->ordered_pattern.size();
    const size_t maxptn  = get_safe_upper_limit_float(nptn);

    if (ptn_freq_pars == nullptr) {
        ptn_freq_pars = aligned_alloc<UINT>(get_safe_upper_limit_float(getAlnNPattern()));
    }
    for (int ptn = 0; ptn < nptn; ptn++) {
        ptn_freq_pars[ptn] = aln->ordered_pattern[ptn].frequency;
    }
    for (int ptn = nptn; ptn < maxptn; ptn++) {
        ptn_freq_pars[ptn] = 0;
    }

    ASSERT(tip_partial_pars);
    memset(tip_partial_pars, 0, (aln->STATE_UNKNOWN+1)*nstates*sizeof(UINT));
    
    // initialize real states with cost_matrix
    memcpy(tip_partial_pars, cost_matrix, nstates*nstates*sizeof(UINT));

    UINT *this_tip_partial_pars;
    
    switch (aln->seq_type) {
        case SEQ_DNA:
            for (int state = 4; state < 18; state++) {
                int cstate = state-nstates+1;
                this_tip_partial_pars = &tip_partial_pars[state*nstates];
                for (int i = 0; i < nstates; i++) {
                    if ((cstate) & (1 << i))
                        this_tip_partial_pars[i] = 0;
                    else {
                        this_tip_partial_pars[i] = UINT_MAX;
                        for (int j = 0; j < nstates; j++) {
                            if ((cstate) & (1 << j)) {
                                this_tip_partial_pars[i] = min(this_tip_partial_pars[i], cost_matrix[i*nstates+j]);
                            }
                        }
                    }
                }
            }
            break;
        case SEQ_PROTEIN:
            {
                // ambiguous characters
                int ambi_aa[] = {
                    4+8, // B = N or D
                    32+64, // Z = Q or E
                    512+1024 // U = I or L
                };
                for (int state = 0; state < sizeof(ambi_aa)/sizeof(int); state++) {
                    this_tip_partial_pars = &tip_partial_pars[(state+20)*nstates];
                    for (int i = 0; i < nstates; i++) {
                        if (ambi_aa[state] & (1 << i)) {
                            this_tip_partial_pars[i] = 0;
                        }
                        else {
                            this_tip_partial_pars[i] = UINT_MAX;
                            for (int j = 0; j < nstates; j++) {
                                if (ambi_aa[state] & (1 << j)) {
                                    this_tip_partial_pars[i] = min(this_tip_partial_pars[i],
                                                                   cost_matrix[i*nstates+j]);
                                }
                            }
                        }
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
    if (dad_branch->isParsimonyComputed()) {
        return;
    }
    PhyloNode *node = dad_branch->getNode();
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
    UINT min_child_ptn_pars;
    
    UINT * partial_pars = dad_branch->partial_pars;

    PhyloNeighbor* left  = nullptr;
    PhyloNeighbor* right = nullptr;
    
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        if (nei->node->name != ROOT_NAME) {
            if (!nei->node->isLeaf())
                computePartialParsimonySankoff(nei, node);
            if (!left) {
                left = nei;
            }
            else {
                right = nei;
            }
        }
    }
    dad_branch->setParsimonyComputed(true);

    if (left==nullptr && right==nullptr && 0<=node->id && node->id<aln->getNSeq() ) {
        //
        //James B. This calculates a partial parsimony vector oriented
        //         at a leaf (as these are needed during parsimony placement,
        //         when CandidateTaxon's constructor is calculating parsimony
        //         for new_interior->findNeighbor(new_leaf).
        //
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (UINT ptn = 0; ptn < aln->ordered_pattern.size(); ++ptn){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            int         ptn_start_index  = ptn*nstates;
            const UINT* leaf_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][node->id]*nstates];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            for (UINT i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                partial_pars_ptr[i] = leaf_ptr[i];
            }
        }
        return;
    }
    
    ASSERT(node->degree() >= 3);
    if (node->degree() > 3) {
        memset(partial_pars, 0, sizeof(UINT)*pars_block_size);
        // multifurcating node
        #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
        #pragma omp parallel for
        #endif
        for (UINT ptn = 0; ptn < aln->ordered_pattern.size(); ptn++) {
            int   ptn_start_index  = ptn*nstates;
            UINT* partial_pars_ptr = &partial_pars[ptn_start_index];

            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) if (nei->node->name != ROOT_NAME) {
                PhyloNode* child = nei->getNode();
                if (child->isLeaf()) {
                    // leaf node
                    const UINT *partial_pars_child_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][child->id]*nstates];
                
                    for (UINT i = 0; i < nstates; i++) {
                        partial_pars_ptr[i] += partial_pars_child_ptr[i];
                    }
                } else {
                    // internal node
                    const UINT* partial_pars_child_ptr = & nei->partial_pars[ptn_start_index];
                    const UINT* cost_matrix_ptr        = cost_matrix;
                    
                    for (UINT i = 0; i < nstates; i++){
                        // min(j->i) from child_branch
                        min_child_ptn_pars = partial_pars_child_ptr[0] + cost_matrix_ptr[0];
                        for(UINT j = 1; j < nstates; j++) {
                            UINT value = partial_pars_child_ptr[j] + cost_matrix_ptr[j];
                            min_child_ptn_pars = min(value, min_child_ptn_pars);
                        }
                        partial_pars_ptr[i] += min_child_ptn_pars;
                        cost_matrix_ptr     += nstates;
                    }
                }
            }
        }
        return;
    }
    
    if (!left->node->isLeaf() && right->node->isLeaf()) {
        std::swap(left, right);
    }
    
    if (left->node->isLeaf() && right->node->isLeaf()) {
        // tip-tip case
        #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
        #pragma omp parallel for
        #endif
        for (UINT ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            int         ptn_start_index  = ptn*nstates;
            const UINT* left_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            const UINT* right_ptr        = &tip_partial_pars[aln->ordered_pattern[ptn][right->node->id]*nstates];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            
            for (UINT i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                partial_pars_ptr[i] = left_ptr[i] + right_ptr[i];
            }
        }
        return;
    }
    
    if (left->node->isLeaf() && !right->node->isLeaf()) {
        // tip-inner case
        // Can't use computePartialParsimonyOutOfTreeSankoff, because
        // reading left state from...
        //   &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
        // ...rather than from...
        //   &some_partial_pars[ptn_start_index].
        
        #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
        #pragma omp parallel for
        #endif
        for (UINT ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
            int         ptn_start_index  = ptn*nstates;
            const UINT* left_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            const UINT* right_ptr        = &right->partial_pars[ptn_start_index];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            const UINT* cost_matrix_ptr  = cost_matrix;
            UINT        right_contrib;
            
            for (UINT i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                for (UINT j = 1; j < nstates; j++) {
                    right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
                }
                partial_pars_ptr[i] = left_ptr[i] + right_contrib;
                cost_matrix_ptr += nstates;
            }
        }
        return;
    }
    // inner-inner case
    computePartialParsimonyOutOfTreeSankoff(left->partial_pars, right->partial_pars, partial_pars );
}

void PhyloTree::computePartialParsimonyOutOfTreeSankoff(const UINT* left_partial_pars,
                                                        const UINT* right_partial_pars,
                                                        UINT* dad_partial_pars) const {
    int nstates  = aln->num_states;
    int ptnCount = aln->ordered_pattern.size();
    #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
    #pragma omp parallel for
    #endif
    for (UINT ptn = 0; ptn < ptnCount; ++ptn){
        // ignore const ptn because it does not affect pars score
        //if (aln->at(ptn).isConst()) continue;
        int ptn_start_index          = ptn * nstates;
        const UINT* left_ptr         = &left_partial_pars[ptn_start_index];
        const UINT* right_ptr        = &right_partial_pars[ptn_start_index];
        UINT*       partial_pars_ptr = &dad_partial_pars[ptn_start_index];
        UINT*       cost_matrix_ptr  = cost_matrix;
        
        for ( UINT i = 0; i < nstates; ++i ){
            // min(j->i) from child_branch
            UINT left_contrib  = left_ptr[0]  + cost_matrix_ptr[0];
            UINT right_contrib = right_ptr[0] + cost_matrix_ptr[0];
            for(UINT j = 1; j < nstates; j++) {
                left_contrib  = min(left_ptr[j] + cost_matrix_ptr[j], left_contrib);
                right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
            }
            partial_pars_ptr[i] = left_contrib + right_contrib;
            cost_matrix_ptr    += nstates;
        }
    }
}

/**
 compute tree parsimony score based on a particular branch
 @param dad_branch the branch leading to the subtree
 @param dad its dad, used to direct the traversal
 @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
 @return parsimony score of the tree
 */
int PhyloTree::computeParsimonyBranchSankoff(PhyloNeighbor *dad_branch,
                                             PhyloNode *dad, int *branch_subst) {
    if ((tip_partial_lh_computed & 2) == 0) {
        computeTipPartialParsimony();
    }
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    assert(node_branch);
    
    if (!central_partial_pars) {
        initializeAllPartialPars();
    }
    
    // if node is a leaf, swap node and dad
    // (so there are two cases to worry about later,
    //  rather than three).
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    
    if (!dad_branch->isParsimonyComputed() && !node->isLeaf()) {
        computePartialParsimonySankoff(dad_branch, dad);
    }
    if (!node_branch->isParsimonyComputed() && !dad->isLeaf()) {
        computePartialParsimonySankoff(node_branch, node);
    }
    
    // now combine likelihood at the branch
    if (dad->isLeaf()) {
        // one of the nodes is external node
        int nstates = aln->num_states;
        UINT tree_pars = 0;
        UINT branch_pars = 0;
        size_t ptnCount = aln->ordered_pattern.size();
        #ifdef _OPENMP //James B. Parallelized this loop
        #pragma omp parallel for reduction(+:tree_pars,branch_pars)
        #endif
        for (int ptn = 0; ptn < ptnCount ; ptn++){
            int         ptn_start_index = ptn * nstates;
            const UINT* node_branch_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][dad->id]*nstates];
            const UINT* dad_branch_ptr  = &dad_branch->partial_pars[ptn_start_index];
            UINT        min_ptn_pars    = node_branch_ptr[0] + dad_branch_ptr[0];
            UINT        br_ptn_pars     = node_branch_ptr[0];
            for (UINT i = 1; i < nstates; i++){
                // min(j->i) from node_branch
                UINT min_score = node_branch_ptr[i] + dad_branch_ptr[i];
                if (min_score < min_ptn_pars) {
                    min_ptn_pars = min_score;
                    br_ptn_pars  = node_branch_ptr[i];
                }
            }
            tree_pars   += min_ptn_pars * aln->ordered_pattern[ptn].frequency;
            branch_pars += br_ptn_pars  * aln->ordered_pattern[ptn].frequency;
        }
        if (branch_subst) {
            *branch_subst = branch_pars;
        }
        return tree_pars;
    }  else {
        // internal node
        return computeParsimonyOutOfTreeSankoff(dad_branch->partial_pars, node_branch->partial_pars, branch_subst);
    }
}

int PhyloTree::computeParsimonyOutOfTreeSankoff(const UINT* dad_partial_pars,
                                                const UINT* node_partial_pars,
                                                int* branch_subst) const {
    int  nstates     = aln->num_states;
    UINT tree_pars   = 0;
    UINT branch_pars = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:tree_pars,branch_pars)
    #endif
    for (int ptn = 0; ptn < aln->ordered_pattern.size(); ptn++){
        int         ptn_start_index = ptn * nstates;
        const UINT* node_branch_ptr = &node_partial_pars[ptn_start_index];
        const UINT* dad_branch_ptr  = &dad_partial_pars[ptn_start_index];
        UINT*       cost_matrix_ptr = cost_matrix;
        UINT        min_ptn_pars    = UINT_MAX;
        UINT        br_ptn_pars     = UINT_MAX;
        for(UINT i = 0; i < nstates; i++){
            // min(j->i) from node_branch
            UINT min_score    = node_branch_ptr[0] + cost_matrix_ptr[0];
            UINT branch_score = cost_matrix_ptr[0];
            for(UINT j = 1; j < nstates; j++) {
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
        tree_pars   += min_ptn_pars * aln->ordered_pattern[ptn].frequency;
        branch_pars += br_ptn_pars  * aln->ordered_pattern[ptn].frequency;
    }
    if (branch_subst) {
        *branch_subst = branch_pars;
    }
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

void PhyloTree::create3TaxonTree(IntVector &taxon_order, int *rand_stream) {
    freeNode();
    size_t nseq = aln->getNSeq();
    taxon_order.resize(nseq);
    for (size_t i = 0; i < nseq; ++i)
        taxon_order[i] = i;
    // randomize the addition order
    my_random_shuffle(taxon_order.begin(), taxon_order.end(), rand_stream);
    
    root = newNode(nseq);
    
    // create star tree
    for (leafNum = 0; leafNum < 3; ++leafNum) {
        if (leafNum < 3 && verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(taxon_order[leafNum]) << " to the tree" << endl;
        Node *new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
        root->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(root, -1.0);
    }
    root = root->neighbors[0]->node;
}

void PhyloTree::copyConstraintTree(MTree *tree, IntVector &taxon_order, int *rand_stream) {
    MTree::copyTree(tree);
    // assign proper taxon IDs
    NodeVector nodes;
    NodeVector::iterator it;
    getTaxa(nodes);
    leafNum = nodes.size();
    vector<int> pushed;
    pushed.resize(aln->getNSeq(), 0);
    
    // name map for fast lookup
    StrVector seq_names = aln->getSeqNames();
    StringIntMap name2id;
    for (auto sit = seq_names.begin(); sit != seq_names.end(); sit++)
        name2id[*sit] = sit - seq_names.begin();
    
    // reindex taxon ID from alignment
    for (it = nodes.begin(); it != nodes.end(); it++) {
        (*it)->id = name2id[(*it)->name];
        ASSERT((*it)->id >= 0);
        taxon_order.push_back((*it)->id);
        pushed[(*it)->id] = 1;
    }
    ASSERT(taxon_order.size() == constraintTree.leafNum);

    // reindex internal nodes properly
    nodes.clear();
    getInternalNodes(nodes);
    for (it = nodes.begin(); it != nodes.end(); it++)
        (*it)->id = aln->getNSeq() + (it - nodes.begin());

    // add the remaining taxa
    for (size_t i = 0; i < aln->getNSeq(); ++i)
        if (!pushed[i]) {
            taxon_order.push_back(i);
        }
    // randomize the addition order
    my_random_shuffle(taxon_order.begin()+constraintTree.leafNum, taxon_order.end(), rand_stream);
}

/**
 get all neighboring branches to a removed node
 */
void getNeiBranches(NeighborVec &removed_nei, NodeVector &attached_node, NodeVector &added_nodes, int i,
                    NodeVector &nodes1, NodeVector &nodes2)
{
    // get target branches surrounding attached_node
    FOR_NEIGHBOR_IT(attached_node[i], NULL, it) {
        if (attached_node[i]->id < (*it)->node->id) {
            nodes1.push_back(attached_node[i]);
            nodes2.push_back((*it)->node);
        } else {
            nodes2.push_back(attached_node[i]);
            nodes1.push_back((*it)->node);
        }
    }
    // get target branches surrounding previous added_nodes
    int j;
    for (j = i-1; j >= 0; j--) {
        if (attached_node[j] != attached_node[i]) {
            break;
        }
        Node *node = added_nodes[j];
        FOR_NEIGHBOR_IT(node, NULL, it) {
            if (node->id < (*it)->node->id) {
                bool present = false;
                for (int k = 0; k < nodes1.size(); k++)
                    if (node == nodes1[k] && (*it)->node == nodes2[k]) {
                        present = true;
                        break;
                    }
                if (present) continue;
                nodes1.push_back(node);
                nodes2.push_back((*it)->node);
            } else {
                bool present = false;
                for (int k = 0; k < nodes1.size(); k++)
                    if (node == nodes2[k] && (*it)->node == nodes1[k]) {
                        present = true;
                        break;
                    }
                if (present) continue;
                nodes2.push_back(node);
                nodes1.push_back((*it)->node);
            }
        }
        // check that exactly two branches are added
    }
    ASSERT(nodes1.size() == 3 + (i-j-1)*2);
}

void PhyloTree::insertNode2Branch(PhyloNode* added_node, PhyloNode* target_node, PhyloNode* target_dad) {
    target_node->updateNeighbor(target_dad, added_node, -1.0);
    target_dad->updateNeighbor(target_node, added_node, -1.0);
    added_node->updateNeighbor(DUMMY_NODE_1, target_node, -1.0);
    added_node->updateNeighbor(DUMMY_NODE_2, target_dad, -1.0);
    added_node->findNeighbor(target_node)->partial_pars =
    target_dad->findNeighbor(added_node)->partial_pars;
    added_node->findNeighbor(target_dad)->partial_pars =
    target_node->findNeighbor(added_node)->partial_pars;
    
    added_node->findNeighbor(target_node)->partial_lh_computed =
    target_dad->findNeighbor(added_node)->partial_lh_computed;
    added_node->findNeighbor(target_dad)->partial_lh_computed =
    target_node->findNeighbor(added_node)->partial_lh_computed;

    PhyloNode* ass_node = added_node->firstNeighbor()->getNode();
    ass_node->findNeighbor(added_node)->clearPartialLh();
    ass_node->clearReversePartialParsimony(added_node);
}

int PhyloTree::computeParsimonyTree(const char *out_prefix, Alignment *alignment, int *rand_stream) {
    aln = alignment;
    size_t nseq = aln->getNSeq();
    if (nseq < 3) {
        outError(ERR_FEW_TAXA);
    }

    IntVector taxon_order;
    taxon_order.reserve(aln->getNSeq());
    
    NeighborVec removed_nei; // removed Neighbor
    PhyloNodeVector attached_node; // node attached to removed Neighbor
    PhyloNodeVector added_nodes; // newly added nodes
    int newNodeID;
    size_t index;
    size_t pars_block_size = getBitsBlockSize();

    if (constraintTree.empty()) {
        create3TaxonTree(taxon_order, rand_stream);
        ASSERT(leafNum == 3);
        initializeAllPartialPars();
        index = (2*leafNum-3)*2;
        newNodeID = nseq + leafNum - 2;
    } else {
        // first copy the constraint tree
        copyConstraintTree(&constraintTree, taxon_order, rand_stream);
        newNodeID = nodeNum - leafNum + nseq;
        index     = (branchNum)*2;
        
        // initialize partial_pars to reuse later
        initializeAllPartialPars();
        
        // extract a bifurcating subtree and get removed nodes to insert later
        extractBifurcatingSubTree(removed_nei, attached_node, rand_stream);
        
        added_nodes.reserve(removed_nei.size());
    }
    if (verbose_mode >= VB_MAX) {
        cout << "computeParsimony: " << computeParsimony() << endl;
    }

    //UINT *tmp_partial_pars;
    //tmp_partial_pars = newBitsBlock();
    if (nseq == 3) {
        best_pars_score = computeParsimony();
    }

    best_pars_score = 0;
    if (leafNum == nseq) {
        outWarning("Constraint tree has all taxa and is bifurcating, which strictly enforces final tree!");
    }
    
    // stepwise adding the next taxon for the remaining taxa
    initProgress(leafNum*2 + nseq*(nseq+1) - leafNum*(leafNum+1),
                 "Constructing parsimony tree", "", "");
    for (int step = 0; leafNum < nseq; ++step) {
        PhyloNodeVector nodes1, nodes2;
        PhyloNode*      target_node = nullptr;
        PhyloNode*      target_dad  = nullptr;
        best_pars_score = UINT_MAX;
        
        // create a new node attached to new taxon or removed node
        PhyloNode *added_node = newNode(newNodeID++);
        PhyloNode *new_taxon;

        if (step < removed_nei.size()) {
            // add the removed_nei (from constraint tree) back to the tree
            getNeiBranches(removed_nei, attached_node, added_nodes, step, nodes1, nodes2);
            new_taxon = (PhyloNode*)removed_nei[step]->node;
            added_node->neighbors.push_back(removed_nei[step]);
            new_taxon->updateNeighbor(attached_node[step], added_node);
            added_nodes.push_back(added_node);
        } else {
            // add new taxon to the tree
            if (verbose_mode >= VB_MAX) {
                hideProgress();
                cout << "Adding " << aln->getSeqName(taxon_order[leafNum]) << " to the tree..." << endl;
                showProgress();
            }
            getBranches(nodes1, nodes2);

            // allocate a new taxon
            new_taxon = newNode(taxon_order[leafNum], aln->getSeqName(taxon_order[leafNum]).c_str());
            
            // link new_taxon and added_node
            added_node->addNeighbor(new_taxon, -1.0);
            new_taxon->addNeighbor(added_node, -1.0);
            
            // allocate memory
            new_taxon->findNeighbor(added_node)->partial_pars = central_partial_pars + ((index++) * pars_block_size);
            added_node->findNeighbor(new_taxon)->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        }
        // preserve two neighbors
        added_node->addNeighbor(DUMMY_NODE_1, -1.0);
        added_node->addNeighbor(DUMMY_NODE_2, -1.0);

        for (int branch_num = 0; branch_num < nodes1.size(); branch_num++) {
            int score = addTaxonMPFast(new_taxon, added_node, nodes1[branch_num], nodes2[branch_num]);
            if (score < best_pars_score) {
                best_pars_score = score;
                target_node = nodes1[branch_num];
                target_dad  = nodes2[branch_num];
            }
        }
        
        if (verbose_mode >= VB_MAX) {
            hideProgress();
            cout << ", score = " << best_pars_score << endl;
            showProgress();
        }
        // now insert the new node in the middle of the branch node-dad
        insertNode2Branch(added_node, target_node, target_dad);

        // assign partial_pars storage
        target_dad->findNeighbor(added_node)->clearComputedFlags();
        target_dad->findNeighbor(added_node)->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        target_node->findNeighbor(added_node)->clearComputedFlags();
        target_node->findNeighbor(added_node)->partial_pars = central_partial_pars + ((index++) * pars_block_size);

        target_dad->clearReversePartialParsimony(added_node);
        target_node->clearReversePartialParsimony(added_node);

        // increase number of taxa
        leafNum += getNumTaxa(new_taxon, added_node);
        trackProgress(nodes1.size());
    }
    doneProgress();
    
    ASSERT(index == 4*leafNum-6);

    nodeNum = 2 * leafNum - 2;
    initializeTree();
    // parsimony tree is always unrooted
    bool orig_rooted = rooted;
    rooted = false;
    setAlignment(alignment);
    fixNegativeBranch(true);
    // convert to rooted tree if originally so
    if (orig_rooted) {
        convertToRooted();
    }
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        printTree(file_name.c_str(), WT_NEWLINE + WT_BR_LEN);
    }
//    if (isSuperTree())
//        ((PhyloSuperTree*)this)->mapTrees();
    
    return best_pars_score;
}

int PhyloTree::addTaxonMPFast(PhyloNode *added_taxon, PhyloNode* added_node, PhyloNode* node, PhyloNode* dad) {

    // now insert the new node in the middle of the branch node-dad
    insertNode2Branch(added_node, node, dad);

    // compute the likelihood
    int score = computeParsimonyBranch(added_taxon->findNeighbor(added_node), added_taxon);

    // remove the added node
    node->updateNeighbor(added_node, dad);
    dad->updateNeighbor(added_node, node);
    added_node->updateNeighbor(node, DUMMY_NODE_1);
    added_node->updateNeighbor(dad, DUMMY_NODE_2);

    // set partial_pars to COMPUTED
    node->findNeighbor(dad)->setParsimonyComputed(true);
    dad->findNeighbor(node)->setParsimonyComputed(true);

    // now tranverse the tree downwards

//    FOR_NEIGHBOR_IT(node, dad, it){
//        addTaxonMPFast(added_node, target_node, target_dad, target_partial_pars, (*it)->node, node);
//    }
    return score;

}

void PhyloTree::extractBifurcatingSubTree(NeighborVec &removed_nei, NodeVector &attached_node, int *rand_stream) {
    NodeVector nodes;
    getMultifurcatingNodes(nodes);
    if (nodes.empty())
        return;
    
    int i;
    
    computeBranchDirection();
    
    // firstly make bifurcating tree
    for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        Node *node = (*it);
        int id[3];
        id[0] = -1;
        // find the neighbor toward root to preserve root
        for (i = 0; i < node->neighbors.size(); i++)
            if (((PhyloNeighbor*)node->neighbors[i])->direction == TOWARD_ROOT) {
                id[0] = i;
                break;
            }
        ASSERT(id[0] >= 0);
        // randomly choose 2 neighbors to reserve
        do {
            id[1] = random_int(node->degree(), rand_stream);
        } while (id[1] == id[0]);
        
        do {
            id[2] = random_int(node->degree(), rand_stream);
        } while (id[2] == id[0] || id[2] == id[1]);
        
        std::sort(id, id+3);
        
        // remove taxa
        int cur_size = removed_nei.size();
        for (i = 0; i < node->degree(); i++)
            if (i != id[0] && i != id[1] && i != id[2]) {
                removed_nei.push_back(node->neighbors[i]);
                attached_node.push_back(node);
            }
        // randomize removed_nei
        my_random_shuffle(removed_nei.begin() + cur_size, removed_nei.end(), rand_stream);
        
        // remove neigbors to make bifurcating tree
        node->neighbors[0] = node->neighbors[id[0]];
        node->neighbors[1] = node->neighbors[id[1]];
        node->neighbors[2] = node->neighbors[id[2]];
        node->neighbors.resize(3);
    }
    
    leafNum = getNumTaxa();
    
}
