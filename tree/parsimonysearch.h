//
//  parsimonysearch.h
//  Created by James Barbetti on 19/1/21.
//

#include "phylotree.h"
#include "phylotreethreadingcontext.h"
#include "parsimonysearchparameters.h"
#include <utils/timekeeper.h>                  //for TimeKeeper
#include <placement/targetbranch.h>            //for TargetBranchRange
#include <placement/placementcostcalculator.h> //for ParsimonyCostCalculator

template <class Move>
int PhyloTree::doParsimonySearch(ParsimonySearchParameters& s) {
    //
    //The "big picture" steps, after initialization is done
    //  1. Evaluation (scoring possible moves) (in parallel).
    //     Moves for each "starting point") (e.g. the "connecting"
    //     branch of a subtree to be pruned and regrafted, in SPR)
    //     are considered *sequentially*, but every thread can run
    //     a sequential search for a different "starting point",
    //     independently of all the other threads. Each search
    //     finds the most beneficial move (if any) with a given
    //     starting point.
    //  2. Sorting (of beneficial moves) (sequential).
    //     Moves that aren't beneficial are filtered out.
    //     Moves that are beneficial are sorted by *ascending*
    //     benefit.
    //  3. Applying (of the most beneficial moves) (sequential).
    //     The beneficial moves are considered from right to left
    //     (order of decreasing benefit).  They are (perhaps
    //     re-evaluated) and applied (and tree parsimony
    //     recalculated). Those that increased the parsimony score
    //     (required more evolution) and backed out of
    //     (taking advantage of the way that immediately calling
    //     ParsimonyMove::apply again is *required* to back out of
    //     the move that was just applied).
    //
    
    best_pars_score          = UINT_MAX;
    intptr_t branch_count    = leafNum * 2 - 3; //assumed leafNum > 3
    int      index_parsimony = 0;
    
    std::string task_name    = "Looking for parsimony " + s.name + " moves";
    if (s.overall.getActivity().empty()) {
        s.overall.setActivity(task_name.c_str());
    }
    if (s.be_quiet) {
        task_name.clear();
    }

    double work_estimate = (double)branch_count * ((double)s.iterations * 2.5 + 1.0);
    //assumes that:
    // A. rescoring parsimony costs 1 (per source branch)
    // B. searching for possible moves costs 1 (per source branch)
    // C. applying moves costs 0.5 (per source branch) (probably a high estimate)
    //    (although only a small fraction of moves are applied, the
    //     partial rescoring necessary to calculate updated parsimony
    //     costs is surprisingly expensive).
    //
    
    initProgress(work_estimate, task_name.c_str(), "", "");

    s.overall.start();
    s.initializing.start();
    
    PhyloTreeThreadingContext context(*this, params->parsimony_uses_max_threads);
    deleteAllPartialLhAndParsimony();
    initializeTree(); //to ensure all branches are properly numbered
    
    intptr_t blocks_per_thread = Move::getParsimonyVectorSize(s.radius);
    intptr_t min_threads       = Move::getMinimumPathVectorCount();
    intptr_t thread_count      = num_threads;
    
    ParsimonyPathVector per_thread_path_parsimony(blocks_per_thread, min_threads, thread_count);
    
    intptr_t num_path_vectors  = per_thread_path_parsimony.getNumberOfPathsRequired();
    intptr_t path_overhead     = per_thread_path_parsimony.getTotalNumberOfBlocksRequired();
        
    auto count_of_branch_vectors_needed =
        ( s.calculate_connection_costs ? branch_count : 0 )
        + path_overhead;
    ensureCentralPartialParsimonyIsAllocated(count_of_branch_vectors_needed);
    initializeAllPartialPars(index_parsimony);
    BlockAllocator           block_allocator(*this, index_parsimony);
    PlacementCostCalculator  dumb_calculator;
    ParsimonyCostCalculator  fancy_calculator(isUsingSankoffParsimony());
    PlacementCostCalculator* calculator = &dumb_calculator;
    if (s.calculate_connection_costs) {
        calculator = &fancy_calculator;
    }
    TargetBranchRange       targets(*this, &block_allocator,
                                    calculator, true);
        
    //Allocate per-thread parsimony vector work areas used to calculate
    //modified parsimony scores along the path between the
    //pruning and regrafting points.
    per_thread_path_parsimony.resize(num_path_vectors);
    for (int vector=0; vector<num_path_vectors; ++vector) {
        block_allocator.allocateVectorOfParsimonyBlocks
            (per_thread_path_parsimony.getBlocksPerThread(),
            per_thread_path_parsimony[vector]);
    }
        
    s.initializing.stop();
    int rc = optimizeSubtreeParsimony<Move>
           (s, targets, per_thread_path_parsimony, context, true);
    doneProgress();
    return rc;
}

template <class Move>
int PhyloTree::optimizeSubtreeParsimony(ParsimonySearchParameters& s,
                                        TargetBranchRange& targets,
                                        ParsimonyPathVector& per_thread_path_parsimony,
                                        PhyloTreeThreadingContext& context,
                                        bool rescore_when_done) {

    TimeKeeper timeSpent("optimizing");
    intptr_t branch_count        = targets.size();
    
    size_t  moves_considered     = 0;
    size_t  moves_still_possible = 0;
    size_t  moves_applied        = 0;
    int64_t positions_considered = 0;
    double parsimony_score;
    timeSpent.start();
    for (intptr_t iteration=1; iteration<=s.iterations; ++iteration) {
        s.rescoring.start();
        parsimony_score = computeParsimony("Determining two-way parsimony", true, true,
                                           targets[0].first->findNeighbor(targets[0].second),
                                           targets[0].first);
        s.rescoring.stop();
        if (!s.be_quiet) {
            if (iteration==1) {
                LOG_LINE(VerboseMode::VB_MIN, "Before doing (up to) "
                         << s.iterations << " rounds of parsimony " << s.name << ","
                         << " parsimony score was " << parsimony_score);
            } else if (!s.be_quiet) {
                LOG_LINE(VerboseMode::VB_MIN, 
                         "Applied " << moves_applied << " move"
                         << ((1==moves_applied) ? "" : "s")
                         << " (out of " << moves_considered << ")"
                         << " (" << moves_still_possible << " still possible)"
                         << " in iteration " << (iteration-1)
                         << " (parsimony now " << parsimony_score << ")"
                         << " after " << timeSpent.getElapsedDescription() );
            }
        }
        s.evaluating.start();
        LikelihoodBlockPairs dummyBlocks(0);
        
        LOG_LINE(VerboseMode::VB_DEBUG, "Computing branch initial states");
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (intptr_t i=0; i<branch_count; ++i) {
            TargetBranch&     tb   = targets[i];
            tb.computeState(*this, parsimony_score, i, dummyBlocks);
            LOG_LINE(VerboseMode::VB_DEBUG, "Branch " << i
                     << " has branch cost " << tb.getBranchCost()
                     << " and connection_cost " << tb.getConnectionCost() );
        }
        LOG_LINE(VerboseMode::VB_DEBUG, "finding best " << s.name << " move for each branch");
        std::vector<Move> moves;
        moves.resize(branch_count);
        
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) reduction(+:positions_considered) schedule(dynamic)
#endif
        for (intptr_t i=0; i<branch_count; ++i) {
            Move&         move   = moves[i];
            int thread = context.getThreadNumber();
            move.initialize(i, s.lazy_mode);
            move.findMove(*this, targets, s.radius,
                          per_thread_path_parsimony[thread],
                          parsimony_score);
            
            move.finalize(*this, targets);
            if (i%100 == 99) {
                trackProgress(100.0);
            }
            positions_considered += move.positions_considered;
        }
        trackProgress(static_cast<double>(branch_count % 100));
        s.evaluating.stop();
        
        LOG_LINE(VerboseMode::VB_DEBUG, "sorting " << s.name << " moves");
        s.sorting.start();
        auto first         = moves.begin();
        auto firstNegative = std::partition(first, moves.end(),
                                            [](const Move& move) { return 0 < move.benefit; } );
        std::sort(first, firstNegative);
        s.sorting.stop();
        
        s.applying.start();
        moves_considered=firstNegative-first;
        LOG_LINE(VerboseMode::VB_MAX, "Considering " << moves_considered
                 << " potentially beneficial " << s.name << " moves");
        moves_applied=0;
        moves_still_possible=0;
        
        size_t i=moves_considered;
        double work_step = 0;
        if (0<moves_considered) {
            work_step = (double) branch_count / (double) moves_considered * 0.1;
        }
        while (0<i) {
            --i;
            Move& move = moves[i];
            LOG_LINE(VerboseMode::VB_DEBUG, "considering " << s.name << " move " << i
                     << " with benefit " << move.benefit
                     << move.getDescription() );
            
            PhyloBranchVector path;
            if ( move.isNoLongerPossible(targets, path) ) {
                //The tree topology has changed and this move
                //is no longer legal (source doesn't exist, target
                //branch doesn't exist, or target branch is now on
                //"the wrong side"of the source branch (so connecting
                //to it would create a disconnected subtree and introduce
                //a cycle containing the source branch).
                LOG_LINE(VerboseMode::VB_DEBUG, "Best move" 
                         " for branch " << move.source_branch_id
                         << " is no longer possible.");
                continue;
            }
            ++moves_still_possible;
            double benefit = move.getBenefit();
            if (s.lazy_mode) {
                benefit = move.recalculateBenefit(*this, parsimony_score,
                                                  targets, dummyBlocks,
                                                  per_thread_path_parsimony) ;
                if ( benefit <= 0) {
                    LOG_LINE(VerboseMode::VB_DEBUG, 
                             "Best move for branch " << move.source_branch_id
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
                //for moves that turned out not to be beneficial, this
                //way, but after the first iteration, almost all moves we
                //try here turn out to be beneficial.
                //So: not only is not re-evaluating the move simpler, it's more
                //efficient (in terms of CPU cost), on average.
            }
            LOG_LINE(VerboseMode::VB_MAX, 
                     "Applying " << s.name << " move " << i
                     << " with original benefit " << move.getBenefit()
                     << " and current benefit " << benefit
                     << move.getDescription());
            double revised_score = move.apply(*this, parsimony_score,
                                              targets, dummyBlocks,
                                              per_thread_path_parsimony);
            if (parsimony_score <= revised_score) {
                const char* same_or_worse = (parsimony_score < revised_score)
                    ? " a worse " : " the same ";
                LOG_LINE(VerboseMode::VB_MAX, "Reverting " << s.name 
                         << " move; as it resulted in"
                         << same_or_worse << "parsimony score"
                         << " (" << revised_score << ")");
                //ParsimonyMove::apply() is its own inverse.  Calling it again
                //with the same parameters, reverses what it did.
                parsimony_score = move.apply(*this, parsimony_score,
                                             targets, dummyBlocks,
                                             per_thread_path_parsimony);
                //Was:
                //    double reverted_score = move.apply ... etc ...
                //    ASSERT( reverted_score == parsimony_score );
                //But:
                //    can't use that for now.
            } else {
                parsimony_score = revised_score;
                ++moves_applied;
            }
            if (!s.be_quiet) {
                trackProgress(work_step);
            }
        }
        s.applying.stop();
        if (moves_applied==0) {
            break;
        }
    }
    
    if (rescore_when_done) {
        s.initializing.start();
        clearAllPartialParsimony(false);
        s.initializing.stop();
        
        s.rescoring.start();
        parsimony_score = computeParsimony("Determining one-way parsimony", false, true,
                                           targets[0].first->findNeighbor(targets[0].second),
                                           targets[0].first);
        s.rescoring.stop();
    }
    
    if (!s.be_quiet) {
        LOG_LINE(VerboseMode::VB_MIN, 
                 "Applied " << moves_applied << " move"
                 << ((1==moves_applied) ? "" : "s")
                 << " (out of " << moves_considered << ")"
                 << " (" << moves_still_possible << " still possible)"
                 << " in last iteration "
                 << " (parsimony now " << parsimony_score << ")"
                 << " (total " << s.name << " moves examined "
                 << positions_considered << ")");
    }
    s.overall.stop();
    hideProgress();
    s.report();
    showProgress();

    return (int)floor(parsimony_score);
}
