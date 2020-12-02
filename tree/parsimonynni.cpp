//
//  parsimonynni.cpp
//  alignment
//
//  Created by James Barbetti on 30/11/20.
//

#include <stdio.h>
#include <placement/blockallocator.h>
#include <placement/targetbranch.h>
#include <placement/parallelparsimonycalculator.h>
#include <utils/timekeeper.h>

class PossibleExchange {
public:
    int         branch_number;
    PhyloNode*  left;
    PhyloBranch middle;
    PhyloNode*  right;
    double delta; //change in parsimony score.  The lower the better!
    PossibleExchange(): left(nullptr), middle(nullptr, nullptr)
                      , right(nullptr), delta(0) {
    }
    void consider(int num, PhyloNode* leftNode,
                  const PhyloBranch& middleBranch,
                  PhyloNode* rightNode, double score_delta) {
        if (score_delta < delta) {
            branch_number = num;
            left          = leftNode;
            middle        = middleBranch;
            right         = rightNode;
            delta         = score_delta;
        }
    }
    bool operator < (const PossibleExchange& rhs) const {
        return delta < rhs.delta;
    }
    bool operator <= (const PossibleExchange& rhs) const {
        return delta <= rhs.delta;
    }
    PhyloNode* getOtherLeftNode() const {
        FOR_EACH_ADJACENT_PHYLO_NODE(middle.first, left, it, node) {
            if (node!=middle.second) {
                return node;
            }
        }
        return nullptr;
    }
    PhyloNode* getOtherRightNode() const {
        FOR_EACH_ADJACENT_PHYLO_NODE(middle.second, right, it, node) {
            if (node!=middle.first) {
                return node;
            }
        }
        return nullptr;
    }
};

#define GET_OTHER_ADJACENT_PHYLO_NODES(node, dad, child_one, child_two) \
    if (1) { \
        child_one = child_two = nullptr; \
        FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) { \
            if (child_one==nullptr) child_one = child; else child_two = child; \
        } \
    } else 0

void PhyloTree::clearInwardParsimonyFromNeighbors(PhyloNode* node1, PhyloNode* node2) {
    //Todo: Do I need this?
    FOR_EACH_ADJACENT_PHYLO_NODE(node1, nullptr, it, node_X)
    {
        PhyloNeighbor* nei = node_X->findNeighbor(node1);
        nei->setParsimonyComputed(false);
        
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node2, nullptr, it, node_Y) {
        PhyloNeighbor* nei = node_Y->findNeighbor(node2);
        nei->setParsimonyComputed(false);
    }
}

void PhyloTree::doParsimonyNNI() {
    size_t max_iterations = params->parsimony_nni_iterations;
    PhyloBranchVector branches;
    getBranches(branches);
    size_t branch_count = branches.size();
    std::vector<TargetBranch> targets;
    
    if (0<params->num_threads) {
        num_threads = params->num_threads;
    } else {
        ensureNumberOfThreadsIsSet(params, true);
    }
    deleteAllPartialParsimony();
    setParsimonyKernel(params->SSE);
    ensureCentralPartialParsimonyIsAllocated( num_threads * 2 + 2 );
    int parsimony_index = 0;
    initializeAllPartialPars(parsimony_index);
    std::vector<UINT*> buffer1, buffer2;
    for (int i=0; i<=num_threads; ++i) {
        size_t offset = parsimony_index * pars_block_size;
        buffer1.push_back( central_partial_pars + offset );
        buffer2.push_back( buffer1[i] + pars_block_size );
        ASSERT( buffer2[i] < tip_partial_pars );
        parsimony_index += 2;
    }

    TimeKeeper rescoring ("rescoring parsimony");
    TimeKeeper evaluating("evaluating NNI moves");
    TimeKeeper sorting   ("sorting NNI moves");
    TimeKeeper applying  ("applying NNI moves");
    
    initProgress(branch_count*max_iterations,
                 "Looking for parsimony nearest neighbor exchanges", "", "");
    for (size_t iteration = 1; iteration <= max_iterations; ++iteration ) {
        rescoring.start();
        double parsimony_score = computeParsimony("Determining two-way parsimony", true);
        rescoring.stop();
        LOG_LINE(VB_MIN, "Parsimony score before parsimony NNI"
                 << " iteration " << iteration
                 << " was " << parsimony_score);

        evaluating.start();
        std::vector<double> branch_cost;
        branch_cost.resize(branch_count);
        
        if ( 1 < iteration ) {
            branches.clear();
            getBranches(branches);
        }
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t i=0; i<branch_count; ++i) {
            PhyloBranch    tb(branches[i]);
            PhyloNeighbor* nei1 = tb.first->findNeighbor(tb.second);
            PhyloNeighbor* nei2 = tb.second->findNeighbor(tb.first);
            int middleParsimony = 0;
            computeParsimonyOutOfTree(nei1->partial_pars,
                                      nei2->partial_pars,
                                      &middleParsimony);
            int rightParsimony = getSubTreeParsimony(nei1, tb.first);
            int leftParsimony  = getSubTreeParsimony(nei2, tb.second);
            branch_cost[i] = leftParsimony + middleParsimony + rightParsimony;
            if (i<10) {
                LOG_LINE(VB_DEBUG, "branch=" << i
                         << ", left=" << leftParsimony
                         << ", mid=" << middleParsimony
                         << ", right=" << rightParsimony);
            }
        }
        std::vector<PossibleExchange> best_exchange;
        best_exchange.resize(branch_count);
        #ifdef _OPENMP
        //#pragma omp parallel for
        #endif
        for (int i=0; i<branch_count; ++i) {
            PhyloBranch tb(branches[i]);
            if ( tb.first->degree()==3 && tb.second->degree()==3 ) {
                int t;
                #ifdef _OPENMP
                    t = omp_get_thread_num();
                #else
                    t = 0;
                #endif
                
                PhyloNode* left1;
                PhyloNode* left2;
                GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                               left1, left2);
                PhyloNode* right1;
                PhyloNode* right2;
                GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                               right1, right2);
                
                ParallelParsimonyCalculator ppc(*this);
                PossibleExchange& best = best_exchange[i];
                
                double cost1 = ppc.parsimonyLink4Cost(left1, right1,
                                                      tb.first, tb.second,
                                                      left2, right2,
                                                      buffer1[t], buffer2[t]);
                LOG_LINE(VB_DEBUG, "for " << i << " cost1 " << cost1 << ","
                         << " oldcost " << branch_cost[i] );

                best.consider(i, left1, tb, right2, cost1 - branch_cost[i]);
                
                double cost2 = ppc.parsimonyLink4Cost(left1, right2, tb.first,
                                                      tb.second, left2, right1,
                                                      buffer1[t], buffer2[t]);
                LOG_LINE(VB_DEBUG, "for " << i << " cost2 " << cost2 << ","
                         << " oldcost " << branch_cost[i] );
                best.consider(i, left1, tb, right1, cost2 - branch_cost[i]);
            }
        }
        evaluating.stop();

        sorting.start();
        std::sort(best_exchange.begin(), best_exchange.end());
        sorting.stop();
        
        applying.start();
        size_t count_considered = 0;
        size_t count_exchanged  = 0;
        size_t expected_delta   = 0;
        size_t actual_delta     = 0;
        for (size_t i=0; i<branch_count; ++i) {
            PossibleExchange& x = best_exchange[i];
            if (0 <= x.delta) {
                break;
            }
            if ( x.middle.first->hasNeighbor(x.left)
                && x.middle.first->hasNeighbor(x.middle.second)
                && x.middle.second->hasNeighbor(x.right) ) {
                ++count_considered;
                expected_delta -= x.delta;
                PhyloNode* left2  = x.getOtherLeftNode();
                PhyloNode* right2 = x.getOtherRightNode();
                ParallelParsimonyCalculator ppc(*this);
                double swap_cost;
                swap_cost = ppc.parsimonyLink4Cost(x.right, left2, x.middle.first,
                                                   x.middle.second, x.left, right2,
                                                   buffer1[0], buffer2[0]);
                double noswap_cost;
                noswap_cost = ppc.parsimonyLink4Cost(x.right, x.left, x.middle.first,
                                                     x.middle.second, left2, right2,
                                                     buffer1[1], buffer2[1]);
                LOG_LINE(VB_DEBUG, "branch " << x.branch_number << ","
                         << " swap cost "    << swap_cost << ","
                         << " no-swap cost " << noswap_cost );
                if (swap_cost < noswap_cost) {
                    LOG_LINE(VB_MAX, "branch " << x.branch_number
                             << " had predicted delta " << x.delta
                             << ", swap cost " << swap_cost
                             << ", and no-swap cost " << noswap_cost);
                    actual_delta += noswap_cost - swap_cost;
                    ++count_exchanged;
                    //Swap inward neighbours
                    x.left->updateNeighbor (x.middle.first,  x.middle.second);
                    x.right->updateNeighbor(x.middle.second, x.middle.first);

                    //Swap outward neighbours and views (still up-to-date!)
                    x.middle.first->updateNeighbor (x.left,  x.right);
                    x.middle.second->updateNeighbor(x.right, x.left);
                    std::swap(x.middle.first->findNeighbor (x.right)->partial_pars,
                              x.middle.second->findNeighbor(x.left )->partial_pars);
                    std::swap(x.middle.first->findNeighbor (x.right)->id,
                              x.middle.second->findNeighbor(x.left )->id);
                    //Swap branch ids too.

                    //Replace views on the middle branch with those we calculated
                    //when figuring out swap_cost.
                    PhyloNeighbor* nei_to_left = x.middle.second->findNeighbor(x.middle.first);
                    std::swap(nei_to_left->partial_pars, buffer1[0]);
                    nei_to_left->setParsimonyComputed(true);
                    
                    PhyloNeighbor* nei_to_right = x.middle.first->findNeighbor(x.middle.second);
                    std::swap(nei_to_right->partial_pars, buffer2[0]);
                    nei_to_right->setParsimonyComputed(true);
                                        
                    //Mark inward views as out of date
                    x.middle.first->clearReversePartialParsimony (x.middle.second);
                    x.middle.second->clearReversePartialParsimony(x.middle.first);
                }
            }
            trackProgress(1.0);
        }
        applying.stop();
        if (count_exchanged==0) {
            break;
        }
        LOG_LINE(VB_MIN, "Parsimony-NNI Iteration " << iteration << ","
                 << " considered " << count_considered << " NNIs,"
                 << " and did " << count_exchanged << " interchanges,"
                 << " with expected delta " << expected_delta << ","
                 << " and actual delta " << actual_delta);
    }
    doneProgress();
    rescoring.start();
    double parsimony_score = computeParsimony("Recalculating one-way parsimony");
    LOG_LINE(VB_MIN, "Parsimony score after parsimony NNI"
             << " was " << parsimony_score);
    rescoring.stop();
    
    if (VB_MED <= verbose_mode) {
        hideProgress();
        std::cout.precision(4);
        rescoring.report();
        evaluating.report();
        sorting.report();
        applying.report();
        showProgress();
    }
    
    deleteAllPartialParsimony();
}
