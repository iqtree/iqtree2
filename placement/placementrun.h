//
// placementrun.h
// Does resource acquisition and release, for a parsimony
// (or likelihood) placement calculation, on a PhyloTree.
//
// Created by James Barbetti on 9/10/20.
//

#ifndef placementrun_h
#define placementrun_h

#include <tree/phylotree.h>          //for PhyloTree
#include "blockallocator.h"          //for BlockAllocator
#include "placement.h"               //for Placement::CostFunction
#include "searchheuristic.h"         //for SearchHeuristic
#include "placementoptimizer.h"      //for TaxonPlacementOptimizer (and others)
#include "placementcostcalculator.h" //for PlacementCostCalculator

#define NEW_TAXON_MAJOR     (0)
#define TARGET_BRANCH_MAJOR (1)

struct PlacementRun {
public:
    PhyloTree&                phylo_tree;
    IntVector                 taxa_ids_to_add;
    size_t                    taxa_per_batch;    //Must be 1 or more
    size_t                    inserts_per_batch; //Must be 1 or more
    BlockAllocator*           block_allocator;
    Placement::CostFunction   costFunction; //parsimony or likelihood
    SearchHeuristic*          heuristic; //global? or localized?
    TaxonPlacementOptimizer*  taxon_placement_optimizer;
    BatchPlacementOptimizer*  batch_placement_optimizer;
    GlobalPlacementOptimizer* global_placement_optimizer;
    PlacementCostCalculator*  calculator;
    
    size_t taxa_inserted_this_batch;
    size_t taxa_inserted_in_total;
    size_t taxa_inserted_nearby;

    PlacementRun(PhyloTree& tree, const IntVector& taxaIdsToAdd);
    void setUpAllocator(int extra_parsimony_blocks, bool trackLikelihood=false,
                        int extra_lh_blocks=0) ;
    void prepareForPlacementRun();
    void prepareForBatch();
    
    template <class T> double doBatchPlacement(TaxaToPlace<T>& candidates,
                                             size_t batchStart, size_t batchStop,
                                             TargetBranchRange& targets) {
        ASSERT( !candidates.empty() );
        T* candidateStart = candidates.data() + batchStart;
        T* candidateStop  = candidates.data() + batchStop;
        ASSERT( !targets.empty() );
        TargetBranch* pointStart = targets.data();
        TargetBranch* pointStop  = pointStart + targets.size();
        
        double refreshTime = 0;
#if (NEW_TAXON_MAJOR)
        for (TargetBranch* point = pointStart; point<pointStop; ++point) {
            point->computeState(phylo_tree);
        }
        for (T* c = candidateStart; c<candidateStop; ++c) {
            TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scoring ... " << c->taxonName);
            c->findPlacement(phylo_tree, targets,
                             heuristic, calculator);
            const PossiblePlacement& p = c->getBestPlacement();
            TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scored " << p.score << " for placement"
                     << " of " << c->taxonName << " with lengths "
                     << p.lenToNode1 << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon);
        }
#else //TARGET_BRANCH_MAJOR
        for (TargetBranch* point = pointStart; point<pointStop; ++point) {
            TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scoring target branch " << (point-pointStart) << " of " << (pointStop-pointStart));
            double computeStart = getRealTime();
            point->computeState(phylo_tree);
            refreshTime += getRealTime() - computeStart;
            point->costPlacementOfTaxa(phylo_tree, &targets, point-pointStart,
                                       candidateStart, candidateStop,
                                       heuristic, calculator,
                                       point==pointStart );
            double forgetStart = getRealTime();
            point->forgetState();
            refreshTime += getRealTime() - forgetStart;
        }
#endif
        return refreshTime;
    }
    
    template <class T>
    void selectPlacementsForInsertion(TaxaToPlace<T>& candidates,
                                      size_t  batchStart, size_t batchStop,
                                      size_t& insertStop) {
        inserts_per_batch = Placement::getInsertsPerBatch(taxa_ids_to_add.size(), batchStop-batchStart);
        insertStop        = batchStart + inserts_per_batch;
        std::sort( candidates.begin() + batchStart, candidates.begin() + batchStop);
        if (batchStop <= insertStop) {
            insertStop = batchStop; //Want them all
        }
    }
    
    void startBatchInsert();
    void insertTaxon(TaxonToPlace& c, TargetBranchRange& targets);
    template <class T>
    void doneBatch(TaxaToPlace<T>& candidates, size_t batchStart, size_t batchStop) {
        if ( 1 < batchStop - batchStart ) {
            TREE_LOG_LINE ( phylo_tree, VB_MED,  "Inserted " << (taxa_inserted_this_batch)
                      << " out of a batch of " << (batchStop - batchStart) << "." );
        }
        batch_placement_optimizer->cleanUpAfterBatch(candidates, batchStart, batchStop, &phylo_tree);
        if (calculator->usesLikelihood()) {
            phylo_tree.fixNegativeBranch();
        }
        if (taxa_inserted_this_batch == 0) {
            outError("No taxa inserted in batch");
        }
    }
    template <class T>
    void donePass(TaxaToPlace<T>& candidates, size_t batchStart, TargetBranchRange& targets) {
        targets.removeUsed();
        //Remove all the candidates that we were able to place
        std::vector<T> oldCandidates;
        std::swap(oldCandidates, candidates);
        //1. Any candidates not considered this time go to the
        //   first batch to consider in the next pass.
        size_t newTaxaCount = candidates.size();
        for (size_t r=batchStart; r<newTaxaCount; ++r) {
            candidates.emplace_back(oldCandidates[r]);
        }
        //2. Any candidates that were considered, but were not
        //   inserted, are to be considered in the next pass.
        for (size_t r=0; r<batchStart; ++r) {
            if (!oldCandidates[r].inserted) {
                //Keep this one to be considered next time
                candidates.emplace_back(oldCandidates[r]);
            }
        }
        inserts_per_batch = Placement::getInsertsPerBatch(taxa_ids_to_add.size(), taxa_per_batch);
        TREE_LOG_LINE ( phylo_tree, VB_MAX, "At the end of this pass, index_lhs was "
                  << block_allocator->getLikelihoodBlockCount() << ", index_pars was "
                  << block_allocator->getParsimonyBlockCount());
    }


    void logSubtreesNearAddedTaxa() const;
    ~PlacementRun();
};

#include <stdio.h>

#endif /* placementrun_h */
