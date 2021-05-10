/*
 * phylotreepars.cpp
 *
 * Fast implementation of parsimony kernel
 *
 *  Created on: May 18, 2015
 *      Author: minh
 */

#include "phylotree.h"
#include "phylotreethreadingcontext.h"
#include "phylosupertree.h"
#include <placement/taxontoplace.h>
#include <placement/blockallocator.h>
#include <placement/parallelparsimonycalculator.h>
#include <utils/timeutil.h> //for getRealTime

#if defined (__GNUC__) || defined(__clang__) && !defined(CLANG_UNDER_VS)
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
static inline uint32_t vml_popcntl (uint64_t ab) 
{
    union {
        uint64_t value64;
        struct {
            //endianness does not matter here, which is why these are 
            //first32 and second32, rather than low32 and high32.
            uint32_t first32;
            uint32_t second32;
        };
    } bitness_converter;
    bitness_converter.value64 = ab;
    return vml_popcnt(bitness_converter.first32) + vml_popcnt(bitness_converter.second32);
}
#endif

#if defined (_MSC_VER)
#if defined ( __SSE4_2__ ) || defined (__AVX__)
#include <nmmintrin.h>
#define __builtin_popcount  _mm_popcnt_u32
#define __builtin_popcountl _mm_popcnt_u64
#else
#define __builtin_popcount  vml_popcnt
#define __builtin_popcountl vml_popcntl
#endif
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
        memset(dad_branch->partial_pars, 255, pars_block_size*sizeof(UINT));
        size_t bits_per_state  = max(aln->size(), (size_t)aln->num_variant_sites);
        size_t uints_per_state = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
        size_t total           = aln->getMaxNumStates() * uints_per_state;
        dad_branch->partial_pars[total]=0; //Todo: don't we want to count
    } else if (node->isLeaf() && dad) {
        // external node
        int leafid = node->id;
        memset(dad_branch->partial_pars, 0, pars_block_size*sizeof(UINT));
        int max_sites = ((aln->num_parsimony_sites+UINT_BITS-1)/UINT_BITS)*UINT_BITS;
        int ambi_aa[] = {2, 3, 5, 6, 9, 10}; // {4+8, 32+64, 512+1024};
//        if (aln->ordered_pattern.empty())
//            aln->orderPatternByNumChars();
        ASSERT(!aln->ordered_pattern.empty());
        int start_pos = 0;
        for (vector<Alignment*>::iterator alnit = partitions->begin(); alnit != partitions->end(); alnit++) {
            int end_pos = start_pos + static_cast<int>((*alnit)->ordered_pattern.size());
            switch ((*alnit)->seq_type) {
            case SeqType::SEQ_DNA:
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
            case SeqType::SEQ_PROTEIN:
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
                    if (aln->seq_type == SeqType::SEQ_POMO && state >= static_cast<int>((*alnit)->num_states) 
                        && state < static_cast<int>((*alnit)->STATE_UNKNOWN)) {
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
        if (site < max_sites) {
            dad_branch->partial_pars[(site/UINT_BITS)*nstates] |= ~((1<<(site%UINT_BITS)) - 1);
        }
        {
            //Count sites where the taxon, indicated by leafid, is the
            //taxon that has a "singleton" state (is the only one in the
            //pattern with this state, when all of the others that have
            //known states have the one known state).
            //
            //(Such sites are parsimony-uninformative but tracking how
            // many of them there are, per taxon, potentially makes
            // parsimony branch lengths more accurate).
            //
            size_t bits_per_state  = max(aln->size(), (size_t)aln->num_variant_sites);
            size_t uints_per_state = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
            size_t total           = aln->getMaxNumStates() * uints_per_state;
            if ( 0 <= leafid && leafid < aln->singleton_parsimony_states.size() ) {
                dad_branch->partial_pars[total] = aln->singleton_parsimony_states[leafid];
            }
        }
    } else {
        // internal node
        ASSERT(node->degree() == 3 || (dad==nullptr && 1<node->degree())  );  // it works only for strictly bifurcating tree
        PhyloNeighbor *left  = nullptr;
        PhyloNeighbor *right = nullptr; // left & right are two neighbors leading to 2 subtrees
        
        //
        //Note: This was running out of stack, in deep trees. So it has
        //      been rewritten to use a vector of things to do, and a
        //      second vector (layers) indicating which of those are in
        //      the same "layer" of the tree (can be done in parallel!):
        //      The content of things_to_do is calculated breadth-first
        //      to make the parallelization easier.
        //
        intptr_t start_of_layer = 0;
        std::vector< std::pair< PhyloNeighbor*, PhyloNode* > > things_to_do;
        std::vector< intptr_t> layers;
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, node_nei) {
            if (node_nei->node->name != ROOT_NAME && !node_nei->isParsimonyComputed()) {
                things_to_do.emplace_back(node_nei, node);
            }
            if (!left) left = node_nei; else right = node_nei;
        }
        while (start_of_layer<static_cast<intptr_t>(things_to_do.size())) {
            layers.push_back(start_of_layer);
            intptr_t start_of_last_layer = start_of_layer;
            start_of_layer = things_to_do.size();
            for (intptr_t i = things_to_do.size()-1; start_of_last_layer<=i; --i ) {
                PhyloNode* node_below = things_to_do[i].first->getNode();
                PhyloNode* node_here  = things_to_do[i].second;
                FOR_EACH_PHYLO_NEIGHBOR(node_below, node_here, it, nei_below) {
                    if (nei_below->node->name != ROOT_NAME && !nei_below->isParsimonyComputed()) {
                        //LOG_LINE(VerboseMode::VB_MIN, "To do " << things_to_do.size()
                        //         << " is child of to do " << i);
                        things_to_do.emplace_back(nei_below, node_below);
                    }
                }
            }
        }
        if (!layers.empty()) {
            //Work up from the bottom layer, computing all the
            //partial parsimonies in each layer, in parallel
            layers.push_back(things_to_do.size());
            for (intptr_t layer_no = layers.size()-2; 0<=layer_no; --layer_no) {
                intptr_t start_index = layers[layer_no];
                intptr_t stop_index  = layers[layer_no+1];
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for (intptr_t i=start_index; i<stop_index; ++i) {
                    PhyloNeighbor* stack_nei  = things_to_do[i].first;
                    PhyloNode*     stack_node = things_to_do[i].second;
                    computePartialParsimonyFast(stack_nei, stack_node);
                    //LOG_LINE(VerboseMode::VB_MIN, "To do " << i
                    //         << " set score " << getSubTreeParsimonyFast(stack_nei));
                }
                things_to_do.resize(start_index);
            }
        }
        ASSERT(left!=nullptr && right!=nullptr);
        if (left!=nullptr && right!=nullptr) {
            computePartialParsimonyOutOfTreeFast(left->partial_pars,
                                                 right->partial_pars,
                                                 dad_branch->partial_pars);
        }
    }
    if (!aln->isSuperAlignment()) {
        delete partitions;
    }
}

double PhyloTree::computePartialParsimonyOutOfTreeFast(const UINT* left_partial_pars,
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
    size_t bits_per_state   = max(aln->size(), (size_t)aln->num_variant_sites);
    size_t uints_per_state  = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
    size_t total            = aln->getMaxNumStates() * uints_per_state;
    dad_partial_pars[total] = score + left_partial_pars[total] + right_partial_pars[total];
    return  dad_partial_pars[total];
}

int PhyloTree::getSubTreeParsimonyFast(PhyloNeighbor* dad_branch) const {
    //Must agree with how computePartialParsimonyOutOfTreeFast
    //records total subtree parsimony
    size_t bits_per_state  = max(aln->size(), (size_t)aln->num_variant_sites);
    size_t uints_per_state = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
    size_t total           = aln->getMaxNumStates() * uints_per_state;
    return dad_branch->partial_pars[total];
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

    size_t bits_per_state  = max(aln->size(), (size_t)aln->num_variant_sites);
    size_t uints_per_state = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
    size_t total           = aln->getMaxNumStates() * uints_per_state;
    
    UINT   sum_end_node    = (dad_partial_pars[total] + node_partial_pars[total]);
    UINT   score           = sum_end_node;

    UINT lower_bound = best_pars_score;
    if (branch_subst) {
        lower_bound = INT_MAX;
    }
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
        if (!nei->isParsimonyComputed()) {
            computePartialParsimony(nei, node);
        }
        PhyloNeighbor *rev = nei->getNode()->findNeighbor(node);
        
        //Note: this, too, was a "likelihood computed?" check.
        if (!rev->isParsimonyComputed()) {
            computePartialParsimony(rev, nei->getNode());
        }
        computeAllPartialPars(nei->getNode(), node);
    }
}

double PhyloTree::jukesCantorCorrection(double dist, double alpha) const {
    double z = (double) aln->num_states / (aln->num_states - 1);
    double x = 1.0 - (z * dist);
    if (x > 0) {
        if (alpha <= 0.0) {
            dist = -log(x) / z;
        } else {
            //if (verbose_mode >= VerboseMode::VB_MAX) cout << "alpha: " << alpha << endl;
            dist = alpha * (pow(x, -1.0/alpha) - 1) / z;
        }
    }
    // Branch lengths under PoMo are #events, which is ~N^2 * #substitutions
    if (aln->seq_type == SeqType::SEQ_POMO)
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
    
    int    sum_score = 0;
    double persite   = 1.0/getAlnNSite();
    double alpha     = (site_rate) ? site_rate->getGammaShape() : 1.0;

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

        for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites;
             s++, bit = bit << 1, real_site++)
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
    double branch_length = correctBranchLengthF81(subst*persite, alpha);
    if (branch_length <= 0.0) {
        branch_length = params->min_branch_length;
    }
    fixOneNegativeBranch(branch_length, dad_branch, dad);
    
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
        const StateType* dadSeq  = sequences[dad->id].data();
        StateType*       nodeSeq = sequences[node->id].data();
        subst = 0;
        for (int site = 0, real_site = 0; site < nsites; site++) {
            size_t offset = nstates*site;
            UINT*  x      = dad_branch->partial_pars + offset;
            UINT   bit    = 1;
            for (int s = 0; s < UINT_BITS && real_site < aln->num_parsimony_sites;
                 s++, bit = bit << 1, real_site++) {
                StateType dad_state = dadSeq[real_site];
                ASSERT(static_cast<int>(dad_state) < nstates);
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
        if (dad->isLeaf()) {
            if (dad->id < aln->singleton_parsimony_states.size()) {
                subst += aln->singleton_parsimony_states[dad->id];
            }
        }
        if (node->isLeaf()) {
            if (node->id < aln->singleton_parsimony_states.size()) {
                subst += aln->singleton_parsimony_states[node->id];
            }
        }
        double branch_length = correctBranchLengthF81(subst*persite, alpha);
        if (branch_length <= 0.0) {
            branch_length = params->min_branch_length;
        }
        fixOneNegativeBranch(branch_length, dad_branch, dad);
        sum_score += subst;
    }
    if (pars_score!=sum_score) {
        hideProgress();
        LOG_LINE(VerboseMode::VB_MIN, "pars_score " << pars_score << " but sum_score " << sum_score);
        showProgress();
    }
    ASSERT(pars_score == sum_score);
    return static_cast<int>(nodes1.size());
}

/****************************************************************************
 Sankoff parsimony function
 ****************************************************************************/

bool PhyloTree::isUsingSankoffParsimony() const {
    return cost_matrix != nullptr;
}

void PhyloTree::stopUsingSankoffParsimony() {
    if (isUsingSankoffParsimony()) {
        deleteAllPartialParsimony();
        aligned_free(cost_matrix);
        setParsimonyKernel(params->SSE);
    }
}

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
        cost_matrix  = aligned_alloc<unsigned int>(cost_nstates * cost_nstates);
        for (int i = 0; i < cost_nstates; i++) {
            for (int j = 0; j < cost_nstates; j++) {
                if (j == i) {
                    cost_matrix[i * cost_nstates + j] = 0;
                }
                else {
                    cost_matrix[i * cost_nstates + j] = 1;
                }
            }
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
        for (int i = 0; i < cost_nstates; i++) {
            for (int j = 0; j < cost_nstates; j++) {
                fin >> cost_matrix[i * cost_nstates + j];
            }
        }
        fin.close();
    }
    
    bool changed = false;
    
    for (int k = 0; k < cost_nstates; ++k) {
        for (int i = 0; i < cost_nstates; ++i) {
            for (int j = 0; j < cost_nstates; ++j) {
                if (cost_matrix[i*cost_nstates+j] > cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j]) {
                    changed = true;
                    cost_matrix[i*cost_nstates+j] = cost_matrix[i*cost_nstates+k] + cost_matrix[k*cost_nstates+j];
                }
            }
        }
    }
    if (changed) {
        cout << "WARING: Cost matrix does not satisfy triangular inenquality and is automatically fixed to:" << endl;
        cout << cost_nstates << endl;
        for (int i = 0; i < cost_nstates; ++i) {
            for (int j = 0; j < cost_nstates; ++j) {
                cout << "  " << cost_matrix[i*cost_nstates+j];
            }
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
    const intptr_t nptn = aln->ordered_pattern.size();
    const intptr_t maxptn = get_safe_upper_limit_float(nptn);

    if (ptn_freq_pars == nullptr) {
        ptn_freq_pars = aligned_alloc<UINT>(get_safe_upper_limit_float(getAlnNPattern()));
    }
    for (intptr_t ptn = 0; ptn < nptn; ptn++) {
        ptn_freq_pars[ptn] = aln->ordered_pattern[ptn].frequency;
    }
    for (intptr_t ptn = nptn; ptn < maxptn; ptn++) {
        ptn_freq_pars[ptn] = 0;
    }

    ASSERT(tip_partial_pars);
    memset(tip_partial_pars, 0, (aln->STATE_UNKNOWN+1)*nstates*sizeof(UINT));
    
    // initialize real states with cost_matrix
    memcpy(tip_partial_pars, cost_matrix, nstates*nstates*sizeof(UINT));

    UINT *this_tip_partial_pars;
    
    switch (aln->seq_type) {
        case SeqType::SEQ_DNA:
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
        case SeqType::SEQ_PROTEIN:
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
        case SeqType::SEQ_POMO:
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
void PhyloTree::computePartialParsimonySankoff(PhyloNeighbor *dad_branch,
                                               PhyloNode *dad){
    // don't recompute the parsimony
    if (dad_branch->isParsimonyComputed()) {
        return;
    }
    PhyloNode* node        = dad_branch->getNode();
    intptr_t   ptnCount    = aln->ordered_pattern.size();
    int        nstates     = aln->num_states;
    intptr_t   total_index = ptnCount*nstates;
    assert(dad_branch->partial_pars);
        
    // internal node
    UINT min_child_ptn_pars;
    UINT * partial_pars = dad_branch->partial_pars;
    
    computeTipPartialParsimony();
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
    UINT score = 0;
    if (left==nullptr && right==nullptr && 0<=node->id && node->id<aln->getNSeq() ) {
        //
        //James B. This calculates a partial parsimony vector oriented
        //         at a leaf (as these are needed during parsimony placement,
        //         when TaxonToPlace's constructor is calculating parsimony
        //         for new_interior->findNeighbor(new_leaf).
        //
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:score)
        #endif
        for (intptr_t ptn = 0; ptn < ptnCount; ++ptn){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            intptr_t    ptn_start_index  = ptn*nstates;
            const UINT* leaf_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][node->id]*nstates];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            for (int i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                partial_pars_ptr[i] = leaf_ptr[i];
            }
            UINT here = partial_pars_ptr[0];
            for ( int i = 0; i < nstates; ++i ){
                UINT there = partial_pars_ptr[i] ;
                here = ( here < there ) ? here : there;
            }
            score += here;
        }
        partial_pars[total_index] = score;
        return;
    }
    ASSERT(node->degree() >= 3);
    if (node->degree() > 3) {
        memset(partial_pars, 0, sizeof(UINT)*pars_block_size);
        // multifurcating node
        #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
        #pragma omp parallel for reduction(+:score)
        #endif
        for (intptr_t ptn = 0; ptn < ptnCount; ptn++) {
            intptr_t ptn_start_index  = ptn*nstates;
            UINT* partial_pars_ptr = &partial_pars[ptn_start_index];
            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
                if (nei->node->name != ROOT_NAME) {
                    PhyloNode* child = nei->getNode();
                    if (child->isLeaf()) {
                        // leaf node
                        const UINT *partial_pars_child_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][child->id]*nstates];
                        
                        for (int i = 0; i < nstates; i++) {
                            partial_pars_ptr[i] += partial_pars_child_ptr[i];
                        }
                    } else {
                        // internal node
                        const UINT* partial_pars_child_ptr = & nei->partial_pars[ptn_start_index];
                        const UINT* cost_matrix_ptr        = cost_matrix;
                        
                        for (int i = 0; i < nstates; i++){
                            // min(j->i) from child_branch
                            min_child_ptn_pars = partial_pars_child_ptr[0] + cost_matrix_ptr[0];
                            for(int j = 1; j < nstates; j++) {
                                UINT value         = partial_pars_child_ptr[j] + cost_matrix_ptr[j];
                                min_child_ptn_pars = min(value, min_child_ptn_pars);
                            }
                            partial_pars_ptr[i] += min_child_ptn_pars;
                            cost_matrix_ptr     += nstates;
                        }
                    }
                }
            }
            UINT here = partial_pars_ptr[0];
            for ( int i = 0; i < nstates; ++i ){
                UINT there = partial_pars_ptr[i] ;
                here = ( here < there ) ? here : there;
            }
            score += here;
        }
        partial_pars[total_index] = score;
        return;
    }
    ASSERT(left!=nullptr);
    ASSERT(right!=nullptr);
    if (left==nullptr || right==nullptr) {
        return;
    }
    if (!left->node->isLeaf() && right->node->isLeaf()) {
        std::swap(left, right);
    }
    if (left->node->isLeaf() && right->node->isLeaf()) {
        // tip-tip case
        #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
        #pragma omp parallel for reduction(+:score)
        #endif
        for (intptr_t ptn = 0; ptn < ptnCount; ++ptn){
            intptr_t    ptn_start_index  = ptn*nstates;
            const UINT* left_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            const UINT* right_ptr        = &tip_partial_pars[aln->ordered_pattern[ptn][right->node->id]*nstates];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            
            for (int i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                UINT sum            = left_ptr[i] + right_ptr[i];
                partial_pars_ptr[i] = sum;
            }
            UINT here = partial_pars_ptr[0];
            for (int i = 0; i < nstates; ++i ){
                UINT there = partial_pars_ptr[i] ;
                here = ( here < there ) ? here : there;
            }
            score += here;
        }
        partial_pars[total_index] = score;
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
        #pragma omp parallel for reduction(+:score)
        #endif
        for (intptr_t ptn = 0; ptn < ptnCount; ptn++){
            intptr_t    ptn_start_index  = ptn*nstates;
            const UINT* left_ptr         = &tip_partial_pars[aln->ordered_pattern[ptn][left->node->id]*nstates];
            const UINT* right_ptr        = &right->partial_pars[ptn_start_index];
            UINT*       partial_pars_ptr = &partial_pars[ptn_start_index];
            const UINT* cost_matrix_ptr  = cost_matrix;
            UINT        right_contrib;
            
            for (int i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                for (int j = 1; j < nstates; j++) {
                    right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
                }
                UINT sum = left_ptr[i] + right_contrib;
                partial_pars_ptr[i] = sum;
                cost_matrix_ptr += nstates;
            }
            UINT here = partial_pars_ptr[0];
            for (int i = 1; i < nstates; ++i) {
                UINT there = partial_pars_ptr[i];
                here = (here<there) ? here : there;
            }
            score += here;
        }
        partial_pars[total_index] = score;
        return;
    }
    // inner-inner case
    computePartialParsimonyOutOfTreeSankoff(left->partial_pars, right->partial_pars, partial_pars );
}

double PhyloTree::computePartialParsimonyOutOfTreeSankoff(const UINT* left_partial_pars,
                                                          const UINT* right_partial_pars,
                                                          UINT* dad_partial_pars) const {
    int      nstates  = aln->num_states;
    intptr_t ptnCount = aln->ordered_pattern.size();
    UINT     score    = 0;
    #ifdef _OPENMP //James B. This for-loop parallelized, 18-Sep-2020.
    #pragma omp parallel for reduction(+:score)
    #endif
    for (intptr_t ptn = 0; ptn < ptnCount; ++ptn){
        // ignore const ptn because it does not affect pars score
        //if (aln->at(ptn).isConst()) continue;
        intptr_t    ptn_start_index  = ptn * nstates;
        const UINT* left_ptr         = &left_partial_pars  [ ptn_start_index ];
        const UINT* right_ptr        = &right_partial_pars [ ptn_start_index ];
        UINT*       partial_pars_ptr = &dad_partial_pars   [ ptn_start_index ];
        UINT*       cost_matrix_ptr  = cost_matrix;
        
        for (int i = 0; i < nstates; ++i ){
            // min(j->i) from child_branch
            UINT left_contrib  = left_ptr[0]  + cost_matrix_ptr[0];
            UINT right_contrib = right_ptr[0] + cost_matrix_ptr[0];
            for(int j = 1; j < nstates; j++) {
                left_contrib  = min(left_ptr[j]  + cost_matrix_ptr[j], left_contrib);
                right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
            }
            partial_pars_ptr[i] = left_contrib + right_contrib;
            cost_matrix_ptr    += nstates;
        }
        UINT here = partial_pars_ptr[0];
        for (int i = 0; i < nstates; ++i ){
            UINT there = partial_pars_ptr[i] ;
            here = ( here < there ) ? here : there;
        }
        score += here;
    }
    size_t totalOffset            = nstates * ptnCount;
    dad_partial_pars[totalOffset] = score;
    return score;
}

int PhyloTree::getSubTreeParsimonySankoff(PhyloNeighbor* dad_branch) const {
    if (dad_branch->partial_pars==nullptr) {
        return 0;
    }
    size_t nstates     = aln->num_states;
    intptr_t ptnCount  = aln->ordered_pattern.size();
    size_t totalOffset = nstates * ptnCount;
    return dad_branch->partial_pars[totalOffset];
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
        int      nstates     = aln->num_states;
        intptr_t ptnCount    = aln->ordered_pattern.size();
        UINT     tree_pars   = 0;
        UINT     branch_pars = 0;
        #ifdef _OPENMP //James B. Parallelized this loop
        #pragma omp parallel for reduction(+:tree_pars,branch_pars)
        #endif
        for (int ptn = 0; ptn < ptnCount ; ptn++){
            int         ptn_start_index = ptn * nstates;
            const UINT* node_branch_ptr = &tip_partial_pars[aln->ordered_pattern[ptn][dad->id]*nstates];
            const UINT* dad_branch_ptr  = &dad_branch->partial_pars[ptn_start_index];
            UINT        min_ptn_pars    = node_branch_ptr[0] + dad_branch_ptr[0];
            UINT        br_ptn_pars     = node_branch_ptr[0];
            for (int i = 1; i < nstates; i++){
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
        return computeParsimonyOutOfTreeSankoff(dad_branch->partial_pars,
                                                node_branch->partial_pars,
                                                branch_subst);
    }
}

int PhyloTree::computeParsimonyOutOfTreeSankoff(const UINT* dad_partial_pars,
                                                const UINT* node_partial_pars,
                                                int* branch_subst) const {
    int    nstates     = aln->num_states;
    intptr_t ptnCount    = aln->ordered_pattern.size();
    UINT   tree_pars   = 0;
    UINT   branch_pars = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:tree_pars,branch_pars)
    #endif
    for (intptr_t ptn = 0; ptn < ptnCount; ptn++){
        intptr_t    ptn_start_index = ptn * nstates;
        const UINT* node_branch_ptr = &node_partial_pars[ptn_start_index];
        const UINT* dad_branch_ptr  = &dad_partial_pars[ptn_start_index];
        UINT*       cost_matrix_ptr = cost_matrix;
        UINT        min_ptn_pars    = UINT_MAX;
        UINT        br_ptn_pars     = UINT_MAX;
        for(int i = 0; i < nstates; i++){
            // min(j->i) from node_branch
            UINT min_score    = node_branch_ptr[0] + cost_matrix_ptr[0];
            UINT branch_score = cost_matrix_ptr[0];
            for(int j = 1; j < nstates; j++) {
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
        auto freq    = aln->ordered_pattern[ptn].frequency;
        tree_pars   += min_ptn_pars * freq;
        branch_pars += br_ptn_pars  * freq;
    }
    if (branch_subst) {
        *branch_subst = branch_pars;
    }
    return tree_pars;
}

/****************************************************************************
 Stepwise addition (greedy) by maximum parsimony
 ****************************************************************************/

void PhyloTree::randomizeTaxonOrder(int* rand_stream,
                                    IntVector& taxon_order) {
    int nseq = aln->getNSeq32();
    taxon_order.resize(nseq);
    for (int i = 0; i < nseq; ++i) {
        taxon_order[i] = static_cast<int>(i);
    }
    // randomize the addition order
    my_random_shuffle(taxon_order.begin(), taxon_order.end(),
                      rand_stream);
}

void PhyloTree::create3TaxonTree(IntVector &taxon_order,
                                 int *rand_stream) {
    freeNode();
    deleteAllPartialLhAndParsimony();
    randomizeTaxonOrder(rand_stream, taxon_order);
    root = newNode(aln->getNSeq32());

    // create star tree
    for (leafNum = 0; leafNum < 3; ++leafNum) {
        auto name = aln->getSeqName(taxon_order[leafNum]);
        if (leafNum < 3 && verbose_mode >= VerboseMode::VB_MAX) {
            cout << "Add " << name << " to the tree" << endl;
        }
        Node *new_taxon = newNode(taxon_order[leafNum], name.c_str());
        root->addNeighbor(new_taxon, -1.0);
        new_taxon->addNeighbor(root, -1.0);
    }
    root = root->neighbors[0]->node;
}

void PhyloTree::copyConstraintTree(MTree *tree, IntVector &taxon_order,
                                   int *rand_stream) {
    MTree::copyTree(tree);
    // assign proper taxon IDs
    NodeVector nodes;
    NodeVector::iterator it;
    getTaxa(nodes);
    leafNum = static_cast<int>(nodes.size());
    vector<int> pushed;
    pushed.resize(aln->getNSeq(), 0);
    
    // name map for fast lookup
    StrVector seq_names = aln->getSeqNames();
    StringIntMap name2id;
    for (auto sit = seq_names.begin(); sit != seq_names.end(); sit++) {
        name2id[*sit] = static_cast<int>(sit - seq_names.begin());
    }
    
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
        (*it)->id = aln->getNSeq32() + static_cast<int>(it - nodes.begin());

    // add the remaining taxa
    for (int i = 0; i < aln->getNSeq(); ++i) {
        if (!pushed[i]) {
            taxon_order.push_back(i);
        }
    }
    // randomize the addition order
    my_random_shuffle(taxon_order.begin()+constraintTree.leafNum, taxon_order.end(), rand_stream);
}

/**
 get all neighboring branches to a removed node
 */
void getNeiBranches(PhyloNeighborVec &removed_nei,
                    NodeVector &attached_node, NodeVector &added_nodes, int i,
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

void PhyloTree::insertNode2Branch(PhyloNode* added_node,
                                  PhyloNode* target_node,
                                  PhyloNode* target_dad) {
    target_node->updateNeighbor(target_dad,   added_node,  -1.0);
    target_dad->updateNeighbor (target_node,  added_node,  -1.0);
    added_node->updateNeighbor (DUMMY_NODE_1, target_node, -1.0);
    added_node->updateNeighbor (DUMMY_NODE_2, target_dad,  -1.0);
    
    PhyloNeighbor* down       = added_node->findNeighbor (target_node);
    PhyloNeighbor* old_down   = target_dad->findNeighbor (added_node);
    PhyloNeighbor* up         = added_node->findNeighbor (target_dad);
    PhyloNeighbor* old_up     = target_node->findNeighbor(added_node);
    
    down->partial_pars        = old_down->partial_pars;
    up->partial_pars          = old_up->partial_pars;
    
    down->partial_lh_computed = old_down->partial_lh_computed;
    up->partial_lh_computed   = old_up->partial_lh_computed;

    PhyloNode* leaf = added_node->firstNeighbor()->getNode();
    leaf->findNeighbor(added_node)->setParsimonyComputed(false);
    leaf->clearReversePartialParsimony(added_node);
}

int PhyloTree::computeParsimonyTree(Alignment* alignment,
                                    int* random_number_stream, 
                                    const char* out_prefix,
                                    const char*& doing_what) {
    int score;
    aln = alignment;
    if (params->use_batch_parsimony_addition) {
        doing_what = "Generalized Batched Parsimony";
        score = computeParsimonyTreeBatch(random_number_stream);
    }
    else if (params->use_compute_parsimony_tree_new) {
        doing_what = "Experimental Step-wise Parsimony";
        score = computeParsimonyTreeNew(random_number_stream);
    }
    else {
        doing_what = "Step-wise Parsimony";
        score = computeParsimonyTreeOld(random_number_stream);
    }
    nodeNum = 2 * leafNum - 2;
    initializeTree();
    
    // parsimony tree is always unrooted
    bool orig_rooted = rooted;
    rooted = false;
    setAlignment(alignment);
    if (params->compute_likelihood) {
        fixNegativeBranch(true);
    } else {
        double score = best_pars_score;
        setAllBranchLengthsFromParsimony(true, score);
    }
    if (orig_rooted) {
        // convert to rooted tree if originally so
        convertToRooted();
    }
    if (out_prefix!=nullptr && 0<strlen(out_prefix)) {
        string file_name = out_prefix;
        file_name += ".parstree";
        printTree(file_name.c_str(), WT_NEWLINE + WT_BR_LEN);
    }
    //    if (isSuperTree())
    //        ((PhyloSuperTree*)this)->mapTrees();
    return score;
}

int PhyloTree::computeParsimonyTreeBatch(int *rand_stream) {
    freeNode();
    deleteAllPartialLhAndParsimony();
    IntVector taxon_order;
    create3TaxonTree(taxon_order, rand_stream);
    int parsimony_score = computeParsimony();

    //Decide on sizes for sample trees.
    std::vector<intptr_t> sizes;
    intptr_t n = taxon_order.size();
    while (3 < n ) {
        sizes.push_back(n);
        n = (intptr_t)floor(pow(n,1.0/3.0)); //pow(n,1.0/3.0)
    }

    n = 3;
    intptr_t taxon_count = taxon_order.size();
    do {
        intptr_t p = sizes.back();
        sizes.pop_back();
        IntVector taxa_this_time;
        taxa_this_time.resize(p-n);
        for (intptr_t i=n; i<p; ++i) {
            taxa_this_time[i-n]=taxon_order[i];
        }        
        std::stringstream desc;
        desc << "Expanding tree to " << p << " taxa";
        addNewTaxaToTree(taxa_this_time, desc.str().c_str(), "C{MP}", true);
        if (p < taxon_count) {
            doParsimonySPR(3, false, 3, true );
        }
        n = p;
    }
    while ( n < taxon_count );
    deleteAllPartialParsimony();
    int score = computeParsimony("Determining parsimony of tree" 
                                 " generated by batch parsimony");
    return score;
}

int PhyloTree::computeParsimonyTreeNew(int *rand_stream) {
    if (!constraintTree.empty()) {
        NodeVector multis;
        getMultifurcatingNodes(multis);
        if (!multis.empty()) {
            //This function doesn't (yet) handle constraint trees
            //that aren't bifurcated.
            return computeParsimonyTreeOld(rand_stream);
        }
    }
    size_t nseq = aln->getNSeq();
    if (nseq < 3) {
        outError(ERR_FEW_TAXA);
    }
    freeNode();
    deleteAllPartialLhAndParsimony();
    IntVector taxon_order;
    taxon_order.resize(nseq);
    if (!constraintTree.empty()) {
        copyConstraintTree(&constraintTree, taxon_order, rand_stream);
    } else {
        create3TaxonTree(taxon_order, rand_stream);
    }
    for (size_t i=0; i<leafNum; ++i) {
        LOG_LINE(VerboseMode::VB_DEBUG, "Initial node " << aln->getSeqName(taxon_order[i]));
    }
    
    PhyloTreeThreadingContext context(*this, params->parsimony_uses_max_threads);
    setParsimonyKernel(params->SSE);
    ensureCentralPartialParsimonyIsAllocated(num_threads);
    int index_parsimony = 0;
    initializeAllPartialPars(index_parsimony);
    BlockAllocator block_allocator(*this, index_parsimony);

    intptr_t count_to_add = taxon_order.size()-leafNum;

    //On each iteration 4*leafNum-4 parsimony vectors
    //are recalculated, so the sum is 2*(nseq*nseq-1)
    double work_estimate = (double)nseq * ((double)nseq - 1.0) * 2.0
        - (double)leafNum * ((double)leafNum - 1.0) * 2.0
        + 4 * (double)leafNum - 4 /*initial tree's parsimony vector count*/
        + (double)count_to_add * (double)num_threads; /*threading overhead*/

    initProgress(work_estimate,
                 "Constructing parsimony tree", "", "");
        
    TypedTaxaToPlace<TaxonToPlace> candidates;
    candidates.resize(count_to_add);
    int first_new_interior_id = renumberInternalNodes();
    for (int i=0; i<count_to_add; ++i) {
        int         taxonId   = taxon_order[i+leafNum];
        std::string taxonName = aln->getSeqName(taxonId);
        candidates[i].initialize(&block_allocator, first_new_interior_id+i,
                                 taxonId, taxonName, true);
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t i=0; i<count_to_add; ++i) {
        candidates[i].computeParsimony(this);
        if (i%100==99) {
            trackProgress(100);
        }
    }
    trackProgress(static_cast<double>(count_to_add%100));

    std::vector<UINT*> buffer;
    block_allocator.allocateVectorOfParsimonyBlocks(num_threads, buffer);

    int parsimony_score;
    {
        ParallelParsimonyCalculator ppc(*this, true);
        parsimony_score = ppc.computeAllParsimony(getRoot()->firstNeighbor(), getRoot());
        LOG_LINE(VerboseMode::VB_DEBUG, "Parsimony score for first "
                 << leafNum << " taxa is " << parsimony_score);
    }
    IntVector scores;
    scores.resize(nseq * 2 - 3);
    PhyloBranchVector branches;
    getBranches(branches);
    for (intptr_t i=0; i<count_to_add; ++i) {
        intptr_t branch_count = branches.size();
        auto pars = candidates[i].getParsimonyBlock();
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(num_threads)
        #endif
        for (intptr_t j=0; j<branch_count; ++j) {
            PhyloBranch b(branches[j]);
            PhyloNeighbor* nei1 = b.getLeftNeighbor();
            PhyloNeighbor* nei2 = b.getRightNeighbor();

            int    t    = context.getThreadNumber();
            double link = computePartialParsimonyOutOfTree(nei1->partial_pars,
                                                           nei2->partial_pars,
                                                           buffer[t]);
            int  branch = 0;
            computeParsimonyOutOfTree(pars, buffer[t], &branch);
            scores[j] = static_cast<int>(link + branch); 
                //cost to link anything there
                //plus cost to link this there
            if ((j%100)==99) {
                trackProgress(100.0);
            }
        }
        trackProgress(static_cast<double>(branch_count % 100));
        intptr_t bestJ = 0;
        double   bestScore = scores[0];
        for (intptr_t j=1; j<branch_count; ++j) {
            if (scores[j]<bestScore) {
                bestScore = scores[j];
                bestJ = j;
            }
        }
        trackProgress(num_threads);
        
#if (0)
        std::stringstream s;
        s << candidates[i].new_leaf->name;
        for (size_t j=0; j<branch_count; ++j) {
            s << " " << scores[j];
        }
        LOG_LINE(VerboseMode::VB_DEBUG, s.str());
#endif
        
        PhyloBranch b(branches[bestJ]);
        PhyloNode*  newInterior = candidates[i].new_interior;
        PhyloNode*  leaf        = candidates[i].new_leaf;
        
        newInterior->addNeighbor(b.first, -1);
        PhyloNeighbor* newNeiLeft  = newInterior->findNeighbor(b.first);
        block_allocator.allocateMemoryFor(newNeiLeft);
        PhyloNeighbor* oldNeiLeft  = b.second->findNeighbor(b.first);
        std::swap(newNeiLeft->partial_pars,  oldNeiLeft->partial_pars);
        newNeiLeft->setParsimonyComputed(true);
        oldNeiLeft->node = newInterior;
        oldNeiLeft->setParsimonyComputed(false);
    
        newInterior->addNeighbor(b.second, -1);
        PhyloNeighbor* newNeiRight = newInterior->findNeighbor(b.second);
        block_allocator.allocateMemoryFor(newNeiRight);
        PhyloNeighbor* oldNeiRight = b.first->findNeighbor(b.second);
        std::swap(newNeiRight->partial_pars,  oldNeiRight->partial_pars);
        newNeiRight->setParsimonyComputed(true);
        oldNeiRight->node = newInterior;
        oldNeiRight->setParsimonyComputed(false);

        newInterior->is_floating_interior = false;

        PhyloNeighbor* leafNei = leaf->findNeighbor(newInterior);
        block_allocator.allocateMemoryFor(leafNei);
        leafNei->setParsimonyComputed(false);
        newInterior->clearReversePartialParsimony(leaf);
        leaf->clearReversePartialParsimony(newInterior);
        
        branches.emplace_back(b.second, newInterior);
        branches.emplace_back(newInterior, leaf);
        branches[bestJ] = PhyloBranch(b.first, newInterior);

        ParallelParsimonyCalculator ppc(*this, true);
        ppc.computeReverseParsimony(newInterior, leaf);
        LOG_LINE(VerboseMode::VB_DEBUG, "After adding " << candidates[i].taxonName
                 << " parsimony score was " << bestScore);
        ++leafNum;
        parsimony_score = (int)floor(bestScore);
    }
    doneProgress();
    parsimony_score = computeParsimony("Computing parsimony score after stepwise addition");
    return parsimony_score; 
}

int PhyloTree::computeParsimonyTreeOld(int *rand_stream) {
    int nseq = aln->getNSeq32();
    if (nseq < 3) {
        outError(ERR_FEW_TAXA);
    }
    IntVector taxon_order;
    taxon_order.reserve(aln->getNSeq());
    
    PhyloNeighborVec removed_nei; // removed Neighbor
    PhyloNodeVector  attached_node; // node attached to removed Neighbor
    PhyloNodeVector  added_nodes; // newly added nodes
    int    newNodeID;
    size_t index;

    PhyloTreeThreadingContext context(*this, params->parsimony_uses_max_threads);
    setParsimonyKernel(params->SSE);
    
    if (constraintTree.empty()) {
        create3TaxonTree(taxon_order, rand_stream);
        ASSERT(leafNum == 3);
        initializeAllPartialPars();
        index     = (2*leafNum-3)*2;
        newNodeID = nseq + leafNum - 2;
    } else {
        // first copy the constraint tree
        double copyStart = getRealTime();
        copyConstraintTree(&constraintTree, taxon_order, rand_stream);
        LOG_LINE(VerboseMode::VB_MED, "Time to copy constraint tree "
                 << (getRealTime()-copyStart) << " secs");
        
        newNodeID = nodeNum + static_cast<int>(nseq-leafNum);
        index     = (branchNum)*2;
        
        // initialize partial_pars to reuse later
        initializeAllPartialPars();
        
        // extract a bifurcating subtree and get removed nodes to insert later
        double extractStart = getRealTime();
        extractBifurcatingSubTree(removed_nei, attached_node, rand_stream);
        LOG_LINE(VerboseMode::VB_MED, "Time to extract bifurcating subtree "
                 << (getRealTime()-extractStart) << " secs");

        added_nodes.reserve(removed_nei.size());
    }
    if (verbose_mode >= VerboseMode::VB_MAX) {
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
    TimeKeeper usefulTime("inserting nodes");
    double fudge      = (double)leafNum * (double)leafNum * 0.01;
    double work_to_do = fudge + (double)leafNum * 2
        + (double)nseq * (double)(nseq + 1)
        - (double)leafNum * (double)(leafNum + 1);
    initProgress(work_to_do,
                 "Constructing parsimony tree", "", "");
    
    for (int step = 0; static_cast<int>(leafNum) < nseq; ++step) {
        PhyloNodeVector nodes1, nodes2;
        PhyloNode*      target_node = nullptr;
        PhyloNode*      target_dad  = nullptr;
        best_pars_score = UINT_MAX;
        
        // create a new node attached to new taxon or removed node
        PhyloNode *added_node = newNode(newNodeID++);
        PhyloNode *new_taxon;
        auto seq_name = aln->getSeqName(taxon_order[leafNum]);

        if (step < static_cast<intptr_t>(removed_nei.size())) {
            // add the removed_nei (from constraint tree) back to the tree
            getNeiBranches(removed_nei, attached_node, added_nodes, step, nodes1, nodes2);
            new_taxon = removed_nei[step]->getNode();
            added_node->neighbors.push_back(removed_nei[step]);
            new_taxon->updateNeighbor(attached_node[step], added_node);
            added_nodes.push_back(added_node);
        } else {
            // add new taxon to the tree
            getBranches(nodes1, nodes2);

            // allocate a new taxon
            new_taxon = newNode(taxon_order[leafNum], seq_name.c_str());
            
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

        usefulTime.start();
        for (int branch_num = 0; branch_num < nodes1.size(); branch_num++) {
            UINT score = addTaxonMPFast(new_taxon, added_node,
                                        nodes1[branch_num], nodes2[branch_num]);
            if (branch_num == 0 || score < best_pars_score) {
                best_pars_score = score;
                target_node = nodes1[branch_num];
                target_dad  = nodes2[branch_num];
            }
        }
        usefulTime.stop();
        
        LOG_LINE(VerboseMode::VB_MAX, "Added " << seq_name << " to the tree"
                 << ", score " << best_pars_score );

        // now insert the new node in the middle of the branch node-dad
        insertNode2Branch(added_node, target_node, target_dad);

        // assign partial_pars storage
        PhyloNeighbor* neiDown = target_dad->findNeighbor(added_node);
        neiDown->clearComputedFlags();
        neiDown->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        
        PhyloNeighbor* neiUp = target_node->findNeighbor(added_node);
        neiUp->clearComputedFlags();
        neiUp->partial_pars = central_partial_pars + ((index++) * pars_block_size);
        
        ASSERT(neiUp->partial_pars < tip_partial_pars);

        target_dad->clearReversePartialParsimony(added_node);
        target_node->clearReversePartialParsimony(added_node);

        // increase number of taxa
        leafNum += getNumTaxa(new_taxon, added_node);
        trackProgress(nodes1.size() + ((step==0) ? fudge : 0));
    }
    doneProgress();
    ASSERT(index == 4*leafNum-6);
    return best_pars_score;
}

UINT PhyloTree::addTaxonMPFast(PhyloNode *added_taxon, PhyloNode* added_node,
                               PhyloNode* node, PhyloNode* dad) {

    // now insert the new node in the middle of the branch node-dad
    insertNode2Branch(added_node, node, dad);

    // compute the parsimony score
    PhyloNeighbor* nei_down = added_taxon->findNeighbor(added_node);
    UINT score = computeParsimonyBranch(nei_down, added_taxon);

    //Would have preferred...
    //  ParallelParsimonyCalculator c(*this);
    //  int score = computeParsimonyBranch(added_taxon->findNeighbor(added_node), added_taxon);
    //...
    //but... on average, three parsimony vector calculations are needed (at different levels)
    //so there is nothing to gain by parallelizing the computation.

    // remove the added node
    node->updateNeighbor      (added_node, dad);
    dad->updateNeighbor       (added_node, node);
    added_node->updateNeighbor(node,       DUMMY_NODE_1);
    added_node->updateNeighbor(dad,        DUMMY_NODE_2);

    // set partial_pars to COMPUTED
    //node->findNeighbor(dad)->setParsimonyComputed(true);
    //dad->findNeighbor(node)->setParsimonyComputed(true);

    return score;
}

void PhyloTree::extractBifurcatingSubTree(PhyloNeighborVec& removed_nei,
                                          PhyloNodeVector&  attached_node,
                                          int *rand_stream) {
    PhyloNodeVector nodes;
    getMultifurcatingNodes(nodes);
    if (nodes.empty()) {
        return;
    }
    int i;
    
    computeBranchDirection();
    
    // firstly make bifurcating tree
    for (auto it = nodes.begin(); it != nodes.end(); it++)
    {
        PhyloNode *node = *it;
        int id[3];
        id[0] = -1;
        // find the neighbor toward root to preserve root
        for (i = 0; i < node->neighbors.size(); i++)
            if (node->getNeighborByIndex(i)->direction == TOWARD_ROOT) {
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
        int cur_size = static_cast<int>(removed_nei.size());
        for (i = 0; i < node->degree(); i++)
            if (i != id[0] && i != id[1] && i != id[2]) {
                removed_nei.push_back(node->getNeighborByIndex(i));
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

int PhyloTree::setAllBranchLengthsFromParsimony(bool    recalculate_parsimony,
                                                double& tree_parsimony) {
    if (recalculate_parsimony) {
        initializeAllPartialPars();
        tree_parsimony = computeParsimony("Recalculating tree parsimony"
                                          " to determine branch lengths",
                                          true, false);
    }
    PhyloNodeVector nodes1;
    PhyloNodeVector nodes2;
    getBranches(nodes1, nodes2);
    intptr_t branch_count = nodes1.size();
    int fixes = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:fixes)
    #endif
    for (intptr_t b = 0; b < branch_count; ++b) {
        PhyloNode* left  = nodes1[b];
        PhyloNode* right = nodes2[b];
        fixes += setOneBranchLengthFromParsimony(tree_parsimony,
                                                 left, right);
    }
    return fixes;
}

int PhyloTree::setOneBranchLengthFromParsimony(double tree_parsimony,
                                               PhyloNode* node1,
                                               PhyloNode* node2 ) const {
    PhyloNeighbor* leftNei  = node1->findNeighbor(node2);
    PhyloNeighbor* rightNei = node2->findNeighbor(node1);
    ASSERT(leftNei->isParsimonyComputed());
    ASSERT(rightNei->isParsimonyComputed());
    int one_if_corrected = (leftNei->length <= 0 ||
                            rightNei->length<=0) ? 1 : 0;
    double branch_cost = tree_parsimony
                       - getSubTreeParsimony(leftNei)
                       - getSubTreeParsimony(rightNei);
    if (branch_cost<1) {
        branch_cost = 1;
    }
    double branch_length    = branch_cost / getAlnNSite();
    double alpha            = (site_rate) ? site_rate->getGammaShape() : 1.0;
    double corrected_length = correctBranchLengthF81(branch_length, alpha);
    if (corrected_length < params->min_branch_length) {
        corrected_length = params->min_branch_length;
    }
    leftNei->length = rightNei->length = corrected_length;
    return one_if_corrected;
}
