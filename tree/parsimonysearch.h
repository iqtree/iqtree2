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
void PhyloTree::doParsimonySearch(const ParsimonySearchParameters& s) {
    intptr_t branch_count    = leafNum * 2 - 3; //assumed > 3
    int      index_parsimony = 0;
    
    TimeKeeper initializing ("initializing");
    TimeKeeper rescoring    ("rescoring parsimony");
    TimeKeeper evaluating   ("evaluating " + s.name + " moves");
    TimeKeeper sorting      ("sorting " + s.name + " moves");
    TimeKeeper applying     ("applying " + s.name + " moves");

    double work_estimate = (double)branch_count * ((double)s.iterations * 2.0 + 1.0);
    initProgress(work_estimate,
                 "Looking for parsimony " + s.name + " moves", "", "");

    initializing.start();
    
    PhyloTreeThreadingContext context(*this, params->parsimony_uses_max_threads);
    std::vector< std::vector<UINT*> > per_thread_path_parsimony;
    intptr_t pv_per_thread = Move::getParsimonyVectorSize(s.radius);
    intptr_t path_overhead = num_threads * pv_per_thread;
        
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
        (pv_per_thread, per_thread_path_parsimony[thread]);
    }
        
    initializing.stop();
    
    size_t  moves_applied    = 0;
    size_t  moves_considered = 0;
    int64_t positions_considered = 0;
    for (size_t iteration=1; iteration<=s.iterations;++iteration) {
        rescoring.start();
        int parsimony_score = computeParsimony("Determining two-way parsimony", true, true );
        rescoring.stop();
        if (iteration==1) {
            LOG_LINE(VB_DEBUG, "Parsimony score before parsimony " << s.name
                     << " iteration " << iteration
                     << " was " << parsimony_score);
        } else {
            LOG_LINE(VB_MIN, "Applied " << moves_applied << " move"
                     << ((1==moves_applied) ? "" : "s")
                     << " (out of " << moves_considered << ")"
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
        LOG_LINE(VB_DEBUG, "finding best " << s.name << " move for each branch");
        std::vector<Move> moves;
        moves.resize(branch_count);
        
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) reduction(+:positions_considered)
#endif
        for (intptr_t i=0; i<branch_count; ++i) {
            TargetBranch& source = targets[i];
            Move&         move   = moves[i];
            BenefitPair   benefit = source.getPartialDisconnectionBenefit(*this, targets);
            //LOG_LINE(VB_MIN, "for s=" << i << " bf=" << benefit.forwardBenefit
            //         << ", bb=" << benefit.backwardBenefit);
            int thread = context.getThreadNumber();
            move.initialize(i, s.lazy_mode);
            move.findMove(*this, targets, s.radius,
                          benefit.forwardBenefit,
                          per_thread_path_parsimony[thread],
                          parsimony_score);
            
            move.finalize(*this, targets);
            if (i%100 == 99) {
                trackProgress(100.0);
            }
            positions_considered += move.positions_considered;
        }
        trackProgress(static_cast<double>(branch_count % 100));
        evaluating.stop();
        
        LOG_LINE(VB_DEBUG, "sorting " << s.name << " moves");
        sorting.start();
        auto first         = moves.begin();
        auto firstNegative = std::partition(first, moves.end(),
                                            [](const Move& move) { return 0 < move.benefit; } );
        std::sort(first, firstNegative);
        sorting.stop();
        
        applying.start();
        moves_considered=firstNegative-first;
        LOG_LINE(VB_MAX, "Considering " << moves_considered
                 << " potentially beneficial " << s.name << " moves");
        moves_applied=0;
        
        size_t i=moves_considered;
        while (0<i) {
            --i;
            Move& move = moves[i];
            LOG_LINE(VB_DEBUG, "considering " << s.name << " move " << i
                     << " with benefit " << move.benefit
                     << move.getDescription() );
            if ( move.getBenefit() <= 0 ) {
                break;
            }
            PhyloBranchVector path;
            if ( move.isNoLongerPossible(targets, path) ) {
                //The tree topology has changed and this move
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
            if (s.lazy_mode) {
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
                //for moves that turned out not to be beneficial, this
                //way, but after the first iteration, almost all moves we
                //try here turn out to be beneficial.
                //So: not only is not re-evaluating the move simpler, it's more
                //efficient (in terms of CPU cost), on average.
            }
            LOG_LINE(VB_MAX, "Applying " << s.name << " move " << i
                     << " with original benefit " << move.getBenefit()
                     << " and current benefit " << benefit
                     << move.getDescription());
            double revised_score = move.apply(*this, targets, dummyBlocks);
            if (parsimony_score <= revised_score) {
                const char* same_or_worse = (parsimony_score < revised_score)
                ? " a worse " : " the same ";
                LOG_LINE(VB_MAX, "Reverting " << s.name << " move; as it resulted in"
                         << same_or_worse << "parsimony score"
                         << " (" << revised_score << ")");
                //ParsimonyMove::apply() is its own inverse.  Calling it again
                //with the same parameters, reverses what it did.
                revised_score = move.apply(*this, targets, dummyBlocks);
                ASSERT( revised_score == parsimony_score );
            } else {
                parsimony_score = revised_score;
                ++moves_applied;
            }
        }
        applying.stop();
        if (moves_applied==0) {
            break;
        }
    }
    initializing.start();
    deleteAllPartialParsimony();
    initializing.stop();
    
    rescoring.start();
    double parsimony_score = computeParsimony();
    rescoring.stop();
    
    LOG_LINE(VB_MIN, "Applied " << moves_applied << " move"
             << ((1==moves_applied) ? "" : "s")
             << " (out of " << moves_considered << ")"
             << " in last iteration "
             << " (parsimony now " << parsimony_score << ")"
             << " (total " << s.name << " moves examined "
             << positions_considered << ")");
    
    doneProgress();
    
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
