//
//  parsimonyspr.cpp
//  Parsimony SPR implementation
//
//  Created by James Barbetti on 8/12/20.
//

#include "phylotree.h"
#include "iqtree.h"
#include <placement/targetbranch.h>            //for TargetBranchRange
#include <placement/placementcostcalculator.h> //for ParsimonyCostCalculator
#include <utils/timekeeper.h>

PhyloNodeVector PhyloTree::getTaxaNodesInIDOrder() const {
    PhyloNodeVector taxa;
    getTaxa(taxa);
    PhyloNodeVector result;
    result.resize(leafNum, nullptr);
    for (auto it=taxa.begin(); it!=taxa.end(); ++it) {
        auto node = *it;
        ASSERT(0<=node->id && node->id<static_cast<int>(leafNum));
        ASSERT(result[node->id] == nullptr);
        result[node->id] = node;
    }
    for (auto itResult=result.begin(); itResult!=result.end(); ++itResult) {
        ASSERT((*itResult)!=nullptr);
    }
    return result;
}

namespace {
    struct ParsimonySPRMove {
    public:
        double   benefit;
        intptr_t source_branch_id;
        intptr_t target_branch_id;
        bool     isForward;        //branch.first moves, branch.second does not
        
        //Members used for checking whether an SPR is still valid
        //after other changes might have been made to the tree.
        PhyloNode* source_first;
        PhyloNode* source_second;
        PhyloNode* target_first;
        PhyloNode* target_second;
        
        int64_t positions_considered;
        
        ParsimonySPRMove(const ParsimonySPRMove& rhs) = default;
        
        ParsimonySPRMove& operator=(const ParsimonySPRMove& rhs) = default;
        
        bool operator < (const ParsimonySPRMove& rhs) const {
            return benefit < rhs.benefit;
        }
        
        bool operator <= (const ParsimonySPRMove& rhs) const {
            return benefit <= rhs.benefit;
        }

        struct SPRSearch {
            public:
            const PhyloTree&         tree;
            const TargetBranchRange& branches;
            const TargetBranch&      source;
            double                   discon;
            ParsimonySPRMove&        put_answer_here;
            
            SPRSearch(const PhyloTree& phylo_tree,
                      const TargetBranchRange& target_branches,
                      double disconnection_benefit,
                      ParsimonySPRMove& output)
                : tree(phylo_tree), branches(target_branches)
                , source(target_branches[output.source_branch_id])
                , discon(disconnection_benefit)
                , put_answer_here(output) {
            }
            void searchForForwardsSPR(PhyloNode* current, PhyloNode* prev,
                                      int radius) {
                --radius;
                FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
                    if (0<radius) {
                        searchForForwardsSPR(next, current, radius);
                    }
                    int target_branch_id = (*it)->id;
                    const TargetBranch& target = branches[target_branch_id];
                    double cost    = source.getForwardConnectionCost(tree, target);
                    double benefit = discon - cost;
                    if (put_answer_here.benefit<benefit) {
                        put_answer_here.benefit          = benefit;
                        put_answer_here.target_branch_id = target_branch_id;
                        put_answer_here.isForward        = true;
                    }
                    ++put_answer_here.positions_considered;
                }
            }
            void searchForBackwardsSPR(PhyloNode* current, PhyloNode* prev,
                                      int radius) {
                --radius;
                FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
                    if (0<radius) {
                        searchForBackwardsSPR(next, current, radius);
                    }
                    int target_branch_id = (*it)->id;
                    const TargetBranch& target = branches[target_branch_id];
                    double cost    = source.getBackwardConnectionCost(tree, target);
                    double benefit = discon - cost;
                    if (put_answer_here.benefit<benefit) {
                        put_answer_here.benefit          = benefit;
                        put_answer_here.target_branch_id = target_branch_id;
                        put_answer_here.isForward        = false;
                    }
                    ++put_answer_here.positions_considered;
                }
            }
        };
    public:
        ParsimonySPRMove() {
            initialize(0);
        }
        void initialize(intptr_t source_branch) {
            benefit          = -1.0;
            source_branch_id = source_branch;
            target_branch_id = -1;
            isForward        = false;
            source_first     = nullptr;
            source_second    = nullptr;
            target_first     = nullptr;
            target_second    = nullptr;
            positions_considered = 0;
        }
        void findForwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                            int radius, double disconnection_benefit) {
            SPRSearch s(tree, branches, disconnection_benefit, *this);
            const TargetBranch& tb = branches[source_branch_id];
            PhyloNode* left;
            PhyloNode* right;
            GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                           left, right);
            s.searchForForwardsSPR(left,  tb.first, radius);
            s.searchForForwardsSPR(right, tb.first, radius);
        }
        void findBackwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                             int radius, double disconnection_benefit) {
            SPRSearch s(tree, branches, disconnection_benefit, *this);
            const TargetBranch& tb = branches[source_branch_id];
            PhyloNode* left;
            PhyloNode* right;
            GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                           left, right);
            s.searchForBackwardsSPR(left,  tb.second, radius);
            s.searchForBackwardsSPR(right, tb.second, radius);
        }
        double getBenefit() const {
            return benefit;
        }
        bool isAConnectedThroughBToC(PhyloNode* a, PhyloNode* b,
                                     PhyloNode* c, PhyloBranchVector& path) const {
            //Note: This uses a breadth-first search, because the
            //SPR radius is likely to be low (and a depth-first
            //search would have an expected running time proportional
            //to the size of the tree) (though it would be prettier!).
            std::vector<PhyloBranch> this_layer;
            std::vector<PhyloBranch> previous_layers;
            this_layer.push_back(PhyloBranch(a,b));
            do
            {
                std::vector<PhyloBranch> next_layer;
                for (auto it=this_layer.begin(); it!=this_layer.end(); ++it) {
                    PhyloBranch branch(*it);
                    FOR_EACH_ADJACENT_PHYLO_NODE(branch.second, branch.first, itNode, x) {
                        if (x==c) {
                            //Reconstruct the pathway, all the way back to a,
                            //from what is in previous_layers.
                            path.emplace_back(branch.second, x);
                            path.emplace_back(branch);
                            PhyloNode* hunting = branch.first;
                            for (auto rev=previous_layers.rbegin();
                                 rev!=previous_layers.rend(); ++rev) {
                                PhyloBranch way_back(*rev);
                                if (way_back.second==hunting) {
                                    path.emplace_back(way_back);
                                    hunting = way_back.first;
                                }
                            }
                            path.reverseAll();
                            return true;
                        }
                        else {
                            next_layer.emplace_back(branch.second, x);
                        }
                    }
                    previous_layers.emplace_back(branch);
                }
                std::swap(this_layer, next_layer);
            } while ( ! this_layer.empty());
            return false;
        }
        bool isStillPossible(const TargetBranchRange& branches,
                             PhyloBranchVector& path) const {
            path.clear();
            const TargetBranch& source = branches[source_branch_id];
            if (source.first  != source_first)   return false;
            if (source.second != source_second) return false;
            const TargetBranch& target = branches[target_branch_id];
            if (target.first  != target_first)   return false;
            if (target.second != target_second) return false;
            if (isForward) {
                return isAConnectedThroughBToC(source.second, source.first, target.first, path);
            } else {
                return isAConnectedThroughBToC(source.first, source.second, target.first, path);
            }
        }
        bool isNoLongerPossible(const TargetBranchRange& branches,
                                PhyloBranchVector& path) const {
            return !isStillPossible(branches, path);
        }
        double recalculateBenefit(PhyloTree& tree, TargetBranchRange& branches,
                                  LikelihoodBlockPairs &blocks) {
            TargetBranch& source = branches[source_branch_id];
            TargetBranch& target = branches[target_branch_id];
            source.computeState(tree, source_branch_id, blocks);
            target.computeState(tree, target_branch_id, blocks);
            BenefitPair benefitPair = source.getPartialDisconnectionBenefit(tree, branches);
            double   benefit = isForward ? benefitPair.forwardBenefit : benefitPair.backwardBenefit;
            double   cost    = isForward
                             ? source.getForwardConnectionCost(tree, target)
                             : source.getBackwardConnectionCost(tree, target);
            return benefit - cost;
        }
        void apply(PhyloTree& tree, LikelihoodBlockPairs blocks, TargetBranchRange& branches) {
            TargetBranch& source     = branches[source_branch_id];
            TargetBranch& target     = branches[target_branch_id];
            PhyloNode*    moved_node = isForward ? source.first : source.second;
            PhyloNode*    other_node = isForward ? source.second : source.first;
            PhyloNode*    new_left   = target.first;
            PhyloNode*    new_right  = target.second;
            PhyloNode*    snip_left;
            PhyloNode*    snip_right;
            GET_OTHER_ADJACENT_PHYLO_NODES(moved_node, other_node, snip_left, snip_right);
            int           snip_left_id  = snip_left->findNeighbor(moved_node)->id;
            int           snip_right_id = snip_right->findNeighbor(moved_node)->id;
            
            //Update node linkage in the tree (note, the neighbors
            //of moved_node have to be done in two phases, steps 3+4
            //and 5+6, just in case there's overlap, with one of
            //new_left, new_right being the same as one of snip_left,
            //and snip_right).
            //
            snip_left ->updateNeighbor(moved_node, snip_right);
            snip_right->updateNeighbor(moved_node, snip_left);
            moved_node->updateNeighbor(snip_left,  DUMMY_NODE_1);
            moved_node->updateNeighbor(snip_right, DUMMY_NODE_2);
            moved_node->updateNeighbor(DUMMY_NODE_1, new_left);
            moved_node->updateNeighbor(DUMMY_NODE_2, new_right);
            new_left  ->updateNeighbor(new_right,  moved_node);
            new_right ->updateNeighbor(new_left,   moved_node);
            
            TargetBranch& left_branch  = branches[snip_left_id];
            left_branch.updateMapping(snip_left_id, new_left, moved_node);
                    
            TargetBranch& right_branch = branches[snip_right_id];
            right_branch.updateMapping(snip_right_id, new_right, moved_node);
            
            target.updateMapping(target_branch_id, snip_left, snip_right);
            
            //Todo: return if logging level isn't high
            
            if (VB_MAX<=verbose_mode) {
                double score;
                left_branch   .computeState (tree, snip_left_id,     blocks);
                right_branch  .computeState (tree, snip_right_id,    blocks);
                target        .computeState (tree, target_branch_id, blocks);
                score = source.computeState (tree, source_branch_id, blocks);

                TREE_LOG_LINE(tree, VB_MAX, "Updated parsimony score"
                              << " after applying SPR move was " << score);
            }
        }
    };
};

void PhyloTree::doParsimonySPR() {
    size_t   max_iterations  = params->parsimony_spr_iterations; //assumed >0
    intptr_t branch_count    = leafNum * 2 - 3;                  //assumed > 3
    int      index_parsimony = 0;

    TimeKeeper initializing ("initializing");
    TimeKeeper rescoring    ("rescoring parsimony");
    TimeKeeper evaluating   ("evaluating SPR moves");
    TimeKeeper sorting      ("sorting SPR moves");
    TimeKeeper applying     ("applying SPR moves");

    double work_estimate = (double)branch_count * ((double)max_iterations * 2.0 + 1.0);
    initProgress(work_estimate,
                 "Looking for parsimony SPR moves", "", "");

    initializing.start();
    bool zeroNumThreadsWhenDone = false;
    if (num_threads<1) {
        zeroNumThreadsWhenDone = true;
        num_threads = params->num_threads;
        if (num_threads==0) {
            #ifdef _OPENMP
                num_threads = omp_get_max_threads();
            #else
                num_threads = 1;
            #endif
        }
        ASSERT(0 < num_threads);
    }
    deleteAllPartialLhAndParsimony();
    initializeTree(); //to ensure all branches are properly numbered
    ensureCentralPartialParsimonyIsAllocated(branch_count);
    initializeAllPartialPars(index_parsimony);
    BlockAllocator          block_allocator(*this, index_parsimony);
    ParsimonyCostCalculator calculator(isUsingSankoffParsimony());
    TargetBranchRange       targets(*this, &block_allocator, &calculator, true);
    initializing.stop();
    
    size_t  spr_moves_applied    = 0;
    size_t  spr_moves_considered = 0;
    int64_t positions_considered = 0;
    for (size_t iteration=1; iteration<=max_iterations;++iteration) {
        rescoring.start();
        int parsimony_score = computeParsimony("Determining two-way parsimony", true, true );
        rescoring.stop();
        if (iteration==1) {
        LOG_LINE(VB_DEBUG, "Parsimony score before parsimony SPR"
                 << " iteration " << iteration
                 << " was " << parsimony_score);
        } else {
            LOG_LINE(VB_MIN, "Applied " << spr_moves_applied << " move"
                         << ((1==spr_moves_applied) ? "" : "s")
                         << " (out of " << spr_moves_considered << ")"
                         << " in iteration " << (iteration-1)
                         << " (parsimony now " << parsimony_score << ")");
        }

        evaluating.start();
        LikelihoodBlockPairs dummyBlocks;

        LOG_LINE(VB_DEBUG, "Computing branch initial states");
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=0; i<branch_count; ++i) {
            TargetBranch&     tb   = targets[i];
            tb.computeState(*this, i, dummyBlocks);
            LOG_LINE(VB_DEBUG, "Branch " << i
                     << " has branch cost " << tb.getBranchCost()
                     << " and connection_cost " << tb.getConnectionCost() );
        }
        LOG_LINE(VB_DEBUG, "finding best SPR move for each branch");
        std::vector<ParsimonySPRMove> moves;
        moves.resize(branch_count);
                
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(num_threads) reduction(+:positions_considered)
        #endif
        for (intptr_t i=0; i<branch_count; ++i) {
            TargetBranch&     tb   = targets[i];
            ParsimonySPRMove& move = moves[i];
            move.initialize(i);
            BenefitPair benefit = tb.getPartialDisconnectionBenefit(*this, targets);
            //LOG_LINE(VB_MIN, "for s=" << i << " bf=" << benefit.forwardBenefit
            //         << ", bb=" << benefit.backwardBenefit);
            
            if (tb.first->isInterior()) {
                move.findForwardSPR(*this, targets, params->sprDist,
                                    benefit.forwardBenefit);
            }
            if (tb.second->isInterior()) {
                move.findBackwardSPR(*this, targets, params->sprDist,
                                     benefit.backwardBenefit);
            }
            LOG_LINE(VB_DEBUG, "move s=" << i << ", d=" << move.target_branch_id
                     << ", f=" << move.isForward << ", b=" << move.benefit);
            move.source_first  = tb.first;
            move.source_second = tb.second;
            const TargetBranch& target = targets[move.target_branch_id];
            move.target_first  = (0<move.benefit) ? target.first  : nullptr;
            move.target_second = (0<move.benefit) ? target.second : nullptr;
            if (i%100 == 99) {
                trackProgress(100.0);
            }
            positions_considered += move.positions_considered;
        }
        trackProgress(static_cast<double>(branch_count % 100));
        evaluating.stop();
    
        LOG_LINE(VB_DEBUG, "sorting SPR moves");
        sorting.start();
        auto first         = moves.begin();
        auto firstNegative = std::partition(first, moves.end(),
                                            [](const ParsimonySPRMove& move){ return 0 < move.benefit; } );
        std::sort(first, firstNegative);
        sorting.stop();
        
        applying.start();
        spr_moves_considered=firstNegative-first;
        LOG_LINE(VB_MAX, "Considering " << spr_moves_considered << " potentially beneficial SPR moves");
        spr_moves_applied=0;

        size_t i=spr_moves_considered;
        while (0<i) {
            --i;
            ParsimonySPRMove& move = moves[i];
            LOG_LINE(VB_DEBUG, "considering SPR move " << i
                     << " with benefit " << move.benefit
                     << " linking branch " << move.source_branch_id
                     << ((move.isForward) ? " forward" : " backward" )
                     << " to branch " << move.target_branch_id);
            if ( move.getBenefit() <= 0 ) {
                break;
            }
            PhyloBranchVector path;
            if ( move.isNoLongerPossible(targets, path) ) {
                //The tree topology has changed and this SPR move
                //is no longer legal (source doesn't exist, target
                //branch doesn't exist, or target branch is now on
                //"the wrong side"of the source branch (so connecting
                //to it would create a disconnected subtree and introduce
                //a cycle containing the source branch).
                LOG_LINE(VB_DEBUG, "Best move for branch " << move.source_branch_id
                         << " is no longer possible.");
                continue;
            }
            double benefit = move.recalculateBenefit(*this, targets, dummyBlocks) ;
            if ( benefit <= 0) {
                LOG_LINE(VB_DEBUG, "Best move for branch " << move.source_branch_id
                         << " is no longer beneficial (net cost delta now " << benefit  << ").");
                continue;
            }
            LOG_LINE(VB_MAX, "Applying SPR move " << i
                     << " with original benefit " << move.benefit
                     << " and current benefit " << benefit
                     << " linking branch " << move.source_branch_id
                     << ((move.isForward) ? " forward" : " backward" )
                     << " to branch " << move.target_branch_id);
            move.apply(*this, dummyBlocks, targets);
            ++spr_moves_applied;
        }
        applying.stop();
        if (spr_moves_applied==0) {
            break;
        }
    }
    doneProgress();
    


    initializing.start();
    deleteAllPartialParsimony();
    initializing.stop();

    rescoring.start();
    double parsimony_score = computeParsimony("Recalculating one-way parsimony", false, true);
    
    LOG_LINE(VB_MIN, "Applied " << spr_moves_applied << " move"
                 << ((1==spr_moves_applied) ? "" : "s")
                 << " (out of " << spr_moves_considered << ")"
                 << " in last iteration "
                 << " (parsimony now " << parsimony_score << ")"
                << " (total SPR moves examined " << positions_considered << ")");
    rescoring.stop();
    if (zeroNumThreadsWhenDone) {
        num_threads = 0;
    }

    if (VB_MED <= verbose_mode) {
        hideProgress();
        std::cout.precision(4);
        initializing.report();
        rescoring.report();
        evaluating.report();
        sorting.report();
        applying.report();
        showProgress();
    }
}

void IQTree::doPLLParsimonySPR() {
    double init_start = getRealTime();
    StrVector oldNames;
    bool areNamesDummied = false;
    if (!isInitializedPLL()) {
        //Names containing '/' characters give PLL trouble.
        //Simplest way to avoid the risk... is to dummy all
        //the names for the (PLL) duration.
        areNamesDummied = true;
        oldNames        = aln->getSeqNames();
        PhyloNodeVector taxa ( getTaxaNodesInIDOrder() );
        ASSERT(taxa.size() == aln->getSeqNames().size());
        intptr_t taxa_count = taxa.size();
        for (intptr_t i=0; i<taxa_count; ++i) {
            std::stringstream dummy;
            dummy << "S" << i;
            std::string dummyName = dummy.str();
            taxa[i]->name = dummyName;
            aln->setSeqName(static_cast<int>(i), dummyName);
        }
        initializePLL(*params);
    }
    string constructedTreeString = getTreeString();
    //LOG_LINE(VB_MAX, constructedTreeString);
    pllReadNewick(constructedTreeString);
    double opt_start = getRealTime();
    int iterations = pllOptimizeWithParsimonySPR(pllInst, pllPartitions,
        params->parsimony_spr_iterations,  
        params->sprDist);
    double read_start = getRealTime();
    pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
        pllInst->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE,
        PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
    string pllTreeString = string(pllInst->tree_string);
    PhyloTree::readTreeString(pllTreeString);
    
    //LOG_LINE(VB_MAX, pllTreeString);x
    if (areNamesDummied) {
        PhyloNodeVector taxa ( getTaxaNodesInIDOrder() );
        intptr_t taxa_count     = taxa.size();
        intptr_t seq_name_count = aln->getSeqNames().size();
        ASSERT(taxa_count == seq_name_count);
        for (intptr_t i=0; i<taxa_count; ++i) {
            taxa[i]->name = oldNames[i];
            aln->setSeqName(static_cast<int>(i), oldNames[i]);
        }
        pllDestroyInstance(pllInst);
        pllInst = nullptr;
    }

    double pars_start = getRealTime();
    initializeAllPartialPars();
    auto optimized_parsimony = computeParsimony("Computing post optimization parsimony");
    LOG_LINE(VB_MIN, "After " << iterations << " rounds of Parsimony SPR,"
        << " parsimony score was " << optimized_parsimony);
    double pll_overhead = (opt_start - init_start) + (pars_start - read_start);
    double spr_time     = (read_start - opt_start);
    LOG_LINE(VB_MED, "PLL overhead was " << pll_overhead << " wall-clock seconds");
    LOG_LINE(VB_MED, "SPR optimization took " << spr_time << " wall-clock seconds");
    
    double fix_start = getRealTime();
    fixNegativeBranch();
    double fix_time  = getRealTime()-fix_start;
    LOG_LINE(VB_MED, "Fixing -ve branches tool " << fix_time << " wall-clock seconds");
}
