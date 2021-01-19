//
//  parsimonyspr.cpp
//  Parsimony SPR implementation
//  Created by James Barbetti on 8-Dec-2020.
//

#include "parsimonymove.h"
#include <placement/targetbranch.h>            //for TargetBranchRange
#include <placement/placementcostcalculator.h> //for ParsimonyCostCalculator
#include <utils/timekeeper.h>                  //for TimeKeeper

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
    struct ParsimonyLazySPRMove : public ParsimonyMove {
    public:
        typedef  ParsimonyLazySPRMove this_type;
        typedef  ParsimonyMove        super;

        intptr_t target_branch_id;
        bool     isForward;        //branch.first moves, branch.second does not
        
        //Members used for checking whether an SPR is still valid
        //after other changes might have been made to the tree.
    private:
        PhyloNode* source_first;
        PhyloNode* source_second;
        PhyloNode* target_first;
        PhyloNode* target_second;
        
    public:
        
        ParsimonyLazySPRMove(const this_type& rhs) = default;
        
        ParsimonyLazySPRMove& operator=(const this_type& rhs) = default;
        
        struct LazySPRSearch {
            public:
            const PhyloTree&         tree;
            const TargetBranchRange& branches;
            const TargetBranch&      source;
            double                   discon;
            ParsimonyLazySPRMove&    put_answer_here;
            
            LazySPRSearch(const PhyloTree& phylo_tree,
                      const TargetBranchRange& target_branches,
                      double disconnection_benefit,
                      ParsimonyLazySPRMove& output)
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
        ParsimonyLazySPRMove(): super() {
            initialize(0, true);
        }
        virtual void initialize(intptr_t source_branch, bool beLazy) {
            lazy             = beLazy;
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
        virtual void finalize(PhyloTree& tree,
                      const TargetBranchRange& branches) {
            if (target_branch_id < 0) {
                return;
            }
            if (0 < benefit) {
                TREE_LOG_LINE(tree, VB_DEBUG, "move s=" << source_branch_id
                    << ",d=" << target_branch_id
                    << ", f=" << isForward
                    << ", b=" << benefit);
            }
            auto source   = branches[source_branch_id];
            source_first  = source.first;
            source_second = source.second;
            auto target   = branches[target_branch_id];
            target_first  = (0<benefit) ? target.first  : nullptr;
            target_second = (0<benefit) ? target.second : nullptr;
        }
        void findForwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                            int radius, double disconnection_benefit) {
            LazySPRSearch s(tree, branches, disconnection_benefit, *this);
            const TargetBranch& tb = branches[source_branch_id];
            PhyloNode* left;
            PhyloNode* right;
            GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                           left, right);
            s.searchForForwardsSPR(left,  tb.first, radius);
            s.searchForForwardsSPR(right, tb.first, radius);
        }
        void findBackwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                             int radius, double disconnection_benefit) {
            LazySPRSearch s(tree, branches, disconnection_benefit, *this);
            const TargetBranch& tb = branches[source_branch_id];
            PhyloNode* left;
            PhyloNode* right;
            GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                           left, right);
            s.searchForBackwardsSPR(left,  tb.second, radius);
            s.searchForBackwardsSPR(right, tb.second, radius);
        }
        virtual bool isStillPossible(const TargetBranchRange& branches,
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
        virtual double recalculateBenefit(PhyloTree& tree, TargetBranchRange& branches,
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
        //Note: callers may assume that ParsimonyLazySPRMove::apply
        //      is its own inverse (that calling it a second time,
        //      immediately after it is called the first time, will
        //      reverse the changes that calling it the first time
        //      made).
        //
        virtual double apply(PhyloTree& tree, LikelihoodBlockPairs blocks,
                     TargetBranchRange& branches) {
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
            left_branch.updateMapping(snip_left_id, new_left, moved_node, true);
                    
            TargetBranch& right_branch = branches[snip_right_id];
            right_branch.updateMapping(snip_right_id, new_right, moved_node, true);
            
            target.updateMapping(target_branch_id, snip_left, snip_right, true);
            //Note: The branch that target_branch_id now refers to,
            //      is the branch that, were we reversing the SPR,
            //      would be the "new" target branch (it's the branch
            //      that came to be when we "snipped out" moved_node).
            //      It has to be this way, if apply() is to be its own
            //      inverse.
            
            double score;
            left_branch   .computeState (tree, snip_left_id,     blocks);
            right_branch  .computeState (tree, snip_right_id,    blocks);
            target        .computeState (tree, target_branch_id, blocks);
            score = source.computeState (tree, source_branch_id, blocks);
            
            left_branch.setParsimonyLength(tree);
            right_branch.setParsimonyLength(tree);
            target.setParsimonyLength(tree);
            source.setParsimonyLength(tree);

            TREE_LOG_LINE(tree, VB_MAX, "Updated parsimony score"
                          << " after applying SPR move was " << score);
            return score;
        } //ParsimonyLazySPRMove::apply
    }; //ParsimonyLazySPRMove

    class ParsimonySPRMove: public ParsimonyLazySPRMove {
        public:
            typedef  ParsimonySPRMove this_type;
            typedef  ParsimonyLazySPRMove super;
            ParsimonySPRMove(const this_type& rhs) = default;
            ParsimonySPRMove& operator=(const this_type& rhs) = default;
            ParsimonySPRMove() = default;
        
            struct ProperSPRSearch: public LazySPRSearch {
                public:
                typedef LazySPRSearch super;
                std::vector<UINT*>& path_parsimony;
                    //This is a vector of size at least radius + 2, for any
                    //radius passed to searchForForwardsSPR or to
                    //searchForBackwardsSPR.  The [tree->params->sprDist+1]
                    //entry in the vector can be a copy (when disconnecting
                    //x from left and right), it is either the parsimony
                    //vector of the view, from x, towards the subtree
                    //containing left, OR the parsimony vector of the view
                    //from x, towards the subtree containing right.
                    //The others are all "private" to the search and will
                    //be calculated (and repeatedly overwritten) during the
                    //search.                
                
                ProperSPRSearch(const PhyloTree& phylo_tree,
                          const TargetBranchRange& target_branches,
                          double disconnection_benefit,
                          std::vector<UINT*>& path_parsimony_to_use,
                          ParsimonyLazySPRMove& output)
                    : super(phylo_tree, target_branches, disconnection_benefit, output)
                    , path_parsimony(path_parsimony_to_use) /*proper*/ {
                }
                PhyloNode* other_adjacent_node(PhyloNode*a, PhyloNode*b, PhyloNode*c) {
                    //Return the other node, adjacent to a, that isn't b or c
                    FOR_EACH_ADJACENT_PHYLO_NODE(a, b, it, x) {
                        if (x!=c) return x;
                    }
                    ASSERT(0 && "could not find other adjacent node");
                    return nullptr;
                }
                void prepareToSearch(PhyloNode* left, PhyloNode* right,
                                     PhyloNode* snipped, int radius) {
                    //Sets up path_parsimony[radius+1] to the view from
                    //snipped, of the subtree containing right.
                    path_parsimony.resize(radius+2, nullptr);
                    path_parsimony[radius+1] = snipped->findNeighbor(right)->get_partial_pars();
                    discon = tree.getSubTreeParsimony(left->findNeighbor(snipped), left)
                           - tree.getSubTreeParsimony(snipped->findNeighbor(right), snipped);
                }
                void searchForForwardsSPR(PhyloNode* current, PhyloNode* prev,
                                          int radius, double parsimony) {
                    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
                        //path_parsimony[radius+1] has already been set "above"
                        //(in the callstack), or in prepareToSearch.
                        UINT*      on_path_vector  = path_parsimony[radius+1];
                        PhyloNode* off_path_node   = other_adjacent_node(current, next, prev);
                        UINT*      off_path_vector = current->findNeighbor(off_path_node)->get_partial_pars();
                        tree.computePartialParsimonyOutOfTree
                            ( on_path_vector, off_path_vector
                            , path_parsimony[radius] );
                        if (1<radius) {
                            searchForForwardsSPR(next, current, radius-1, parsimony);
                        }
                        int target_branch_id = (*it)->id;
                        
                        double pruned_tree_score = tree.computePartialParsimonyOutOfTree
                                                   ( current->findNeighbor(next)->get_partial_pars()
                                                   , path_parsimony[radius], path_parsimony[0] );
                        int new_branch_cost = 0;
                        tree.computeParsimonyOutOfTree
                            ( path_parsimony[0]
                            , source.first->findNeighbor(source.second)->get_partial_pars()
                            , &new_branch_cost );
                        double subtree_cost = tree.getSubTreeParsimony
                                              ( source.first->findNeighbor(source.second),
                                                source.first );
                        double benefit = parsimony       - pruned_tree_score
                                       - new_branch_cost - subtree_cost;
                        
                        if (put_answer_here.benefit<benefit) {
                            put_answer_here.benefit          = benefit;
                            put_answer_here.target_branch_id = target_branch_id;
                            put_answer_here.isForward        = true;
                        }
                        ++put_answer_here.positions_considered;
                    }
                }
                void searchForBackwardsSPR(PhyloNode* current, PhyloNode* prev,
                                          int radius, double parsimony) {
                    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
                        UINT*      on_path_vector  = path_parsimony[radius+1];
                        PhyloNode* off_path_node   = other_adjacent_node(current, next, prev);
                        UINT*      off_path_vector = current->findNeighbor(off_path_node)->get_partial_pars();
                        tree.computePartialParsimonyOutOfTree
                            ( on_path_vector, off_path_vector
                            , path_parsimony[radius] );
                        if (1<radius) {
                            searchForBackwardsSPR(next, current, radius-1, parsimony);
                        }
                        int target_branch_id = (*it)->id;
                        
                        double pruned_tree_score = tree.computePartialParsimonyOutOfTree
                                                   ( current->findNeighbor(next)->get_partial_pars()
                                                   , path_parsimony[radius], path_parsimony[0] );
                        int new_branch_cost = 0;
                        tree.computeParsimonyOutOfTree
                            ( path_parsimony[0]
                            , source.second->findNeighbor(source.first)->get_partial_pars()
                            , &new_branch_cost );
                        double subtree_cost = tree.getSubTreeParsimony
                                              ( source.second->findNeighbor(source.first),
                                                source.second );
                        double benefit = parsimony       - pruned_tree_score
                                       - new_branch_cost - subtree_cost;
                        
                        if (put_answer_here.benefit<benefit) {
                            put_answer_here.benefit          = benefit;
                            put_answer_here.target_branch_id = target_branch_id;
                            put_answer_here.isForward        = false;
                        }
                        ++put_answer_here.positions_considered;
                    }
                }
            };
            void findForwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                                int radius, double disconnection_benefit,
                                std::vector<UINT*> &path_parsimony,
                                double parsimony_score) {
                if (lazy) {
                    super::findForwardLazySPR(tree, branches, radius, disconnection_benefit);
                    return;
                }
                ProperSPRSearch s(tree, branches, disconnection_benefit,
                                  path_parsimony, *this);
                const TargetBranch& tb = branches[source_branch_id];
                PhyloNode* left;
                PhyloNode* right;
                GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                               left, right);
                s.prepareToSearch     (left,  right,    tb.first, radius);
                s.searchForForwardsSPR(left,  tb.first, radius,   parsimony_score);
                s.prepareToSearch     (right, left,     tb.first, radius);
                s.searchForForwardsSPR(right, tb.first, radius,   parsimony_score);
            }
            void findBackwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                                 int radius, double disconnection_benefit,
                                 std::vector<UINT*> &path_parsimony,
                                 double parsimony_score) {
                if (lazy) {
                    super::findBackwardLazySPR(tree, branches, radius, disconnection_benefit);
                    return;
                }
                ProperSPRSearch s(tree, branches, disconnection_benefit,
                                  path_parsimony, *this);
                const TargetBranch& tb = branches[source_branch_id];
                PhyloNode* left;
                PhyloNode* right;
                GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                               left, right);
                s.prepareToSearch      (left,  right,     tb.second, radius);
                s.searchForBackwardsSPR(left,  tb.second, radius,    parsimony_score);
                s.prepareToSearch      (right, left,      tb.second, radius);
                s.searchForBackwardsSPR(right, tb.second, radius,    parsimony_score);
            }
    }; //ParsimonySPRMove
}; //namespace

void PhyloTree::doParsimonySPR() {
    size_t   max_iterations  = params->parsimony_spr_iterations; //assumed >0
    bool     lazy_mode       = params->use_lazy_parsimony_spr;
    auto     radius          = params->sprDist;
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
    size_t path_overhead = num_threads * (radius+1);
    std::vector< std::vector<UINT*> > per_thread_path_parsimony;
    
    deleteAllPartialLhAndParsimony();
    initializeTree(); //to ensure all branches are properly numbered
    ensureCentralPartialParsimonyIsAllocated(branch_count + path_overhead);
    initializeAllPartialPars(index_parsimony);
    BlockAllocator          block_allocator(*this, index_parsimony);
    ParsimonyCostCalculator calculator(isUsingSankoffParsimony());
    TargetBranchRange       targets(*this, &block_allocator, &calculator, true);
    
    //Allocate per-thread parsimony vector work areas used to calculate
    //modified parsimony scores along the path between the
    //pruning and regrafting points.
    per_thread_path_parsimony.resize(num_threads);
    for (int thread=0; thread<num_threads; ++thread) {
        block_allocator.allocateVectorOfParsimonyBlocks
            ( radius+1, per_thread_path_parsimony[thread]);
    }
    
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
            TargetBranch&     source = targets[i];
            ParsimonySPRMove& move   = moves[i];
            BenefitPair benefit = source.getPartialDisconnectionBenefit(*this, targets);
            //LOG_LINE(VB_MIN, "for s=" << i << " bf=" << benefit.forwardBenefit
            //         << ", bb=" << benefit.backwardBenefit);
#ifdef _OPENMP
            int thread = omp_get_thread_num();
            ASSERT(0<=thread && thread<num_threads);
#else
            int thread = 0;
#endif
            move.initialize(i, lazy_mode);
            if (source.first->isInterior()) {
                move.findForwardSPR(*this, targets, radius,
                                    benefit.forwardBenefit,
                                    per_thread_path_parsimony[thread],
                                    parsimony_score);
            }
            if (source.second->isInterior()) {
                move.findBackwardSPR(*this, targets, radius,
                                     benefit.backwardBenefit,
                                     per_thread_path_parsimony[thread],
                                     parsimony_score);
            }
            move.finalize(*this, targets);
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
                                            [](const ParsimonyLazySPRMove& move) { return 0 < move.benefit; } );
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
            double benefit = move.getBenefit();
            if (lazy_mode) {
                benefit = move.recalculateBenefit(*this, targets, dummyBlocks) ;
                if ( benefit <= 0) {
                    LOG_LINE(VB_DEBUG, "Best move for branch " << move.source_branch_id
                             << " is no probably longer beneficial"
                             << " (net cost delta now " << benefit  << ").");
                    continue;
                }
            }
            else {
                //If we're not running in lazy mode, we'll just apply the move
                //(since figuring out what the benefit would really be takes
                //so long we might as well do it and see what is after we've
                //done it!). If it were re-evaluated properly, and then done,
                //that'd mean calculating new parsimony views twice (and we'd
                //never back out).
                //But if it is *not* re-evaluated properly, but just done, parsimony
                //views are calculated ONCE unless the move is not an improvement
                //(in which case it gets re-evaluated when backing out).
                //Admittedly more parsimony views get marked as out of date,
                //for SPR moves that turned out not to be beneficial, this
                //way, but after the first iteration, almost all SPR moves we
                //try here turn out to be beneficial.
                //So: not only is not re-evaluating the move simpler, it's more
                //efficient (in terms of CPU cost), on average.
            }
            LOG_LINE(VB_MAX, "Applying SPR move " << i
                     << " with original benefit " << move.getBenefit()
                     << " and current benefit " << benefit
                     << " linking branch " << move.source_branch_id
                     << ((move.isForward) ? " forward" : " backward" )
                     << " to branch " << move.target_branch_id);
            double revised_score = move.apply(*this, dummyBlocks, targets);
            if (parsimony_score <= revised_score) {
                const char* same_or_worse = (parsimony_score < revised_score)
                    ? " a worse " : " the same ";
                LOG_LINE(VB_MAX, "Reverting SPR move; as it resulted in"
                         << same_or_worse << "parsimony score"
                         << " (" << revised_score << ")");
                //ParsimonyMove::apply() is its own inverse.  Calling it again
                //with the same parameters, reverses what it did.
                revised_score = move.apply(*this, dummyBlocks, targets);
                ASSERT( revised_score == parsimony_score );
            } else {
                parsimony_score = revised_score;
                ++spr_moves_applied;
            }
        }
        applying.stop();
        if (spr_moves_applied==0) {
            break;
        }
    }
    initializing.start();
    deleteAllPartialParsimony();
    initializing.stop();

    rescoring.start();
    double parsimony_score = computeParsimony();
    rescoring.stop();

    LOG_LINE(VB_MIN, "Applied " << spr_moves_applied << " move"
                 << ((1==spr_moves_applied) ? "" : "s")
                 << " (out of " << spr_moves_considered << ")"
                 << " in last iteration "
                 << " (parsimony now " << parsimony_score << ")"
                 << " (total SPR moves examined " << positions_considered << ")");

    doneProgress();

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
