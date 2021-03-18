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
#include "placement.h"               //for PlacementParameters
#include "searchheuristic.h"         //for SearchHeuristic
#include "placementoptimizer.h"      //for TaxonPlacementOptimizer (and others)
#include "placementcostcalculator.h" //for PlacementCostCalculator
#include <utils/timekeeper.h>        //for TimeKeeper

#define NEW_TAXON_MAJOR     (0)
#define TARGET_BRANCH_MAJOR (1)

class TaxonPlacement {
    //Used to group inserts by target branch (the cache
    //works better, if they're grouped in this way!).
    //This class declared here because only PlacementRun uses it.
public:
    size_t                   candidate_index; //which one
    size_t                   target_index;    //goes where
    const PossiblePlacement* placement;       //with what score
    bool operator <  ( const TaxonPlacement& rhs) const;
    bool operator <= ( const TaxonPlacement& rhs) const;
    TaxonPlacement();
    TaxonPlacement(const TaxonPlacement& rhs) = default;
};


struct PlacementRun {
public:
    PhyloTree&                 phylo_tree;
    const PlacementParameters& placement_params;
    IntVector                  taxa_ids_to_add;   //id numbers of taxa to add
    bool                       be_quiet;
    size_t                     taxa_per_batch;    //Must be 1 or more
    size_t                     inserts_per_batch; //Must be 1 or more
    BlockAllocator*            block_allocator;   //
    SearchHeuristic*           heuristic;         //global? or localized?
    PlacementCostCalculator*   calculator;        //
    bool                       use_likelihood;    //true if heuristic or calculator do

    TaxonPlacementOptimizer*   taxon_placement_optimizer;
    BatchPlacementOptimizer*   batch_placement_optimizer;
    GlobalPlacementOptimizer*  global_placement_optimizer;

    size_t                     taxa_inserted_this_batch;
    size_t                     taxa_inserted_in_total; //for a given pass
    size_t                     taxa_inserted_nearby;   //taxa that couldn't be inserted at their
                                                       //preferred branch, but were inserted nearby

    TimeKeeper overall;
    TimeKeeper initializing;
    TimeKeeper refreshTime;
    TimeKeeper searchTime;
    TimeKeeper rankingTime;
    TimeKeeper insertTime;
    TimeKeeper optoTime;

    PlacementRun(PhyloTree& tree, const PlacementParameters& parameters,
                 const IntVector& taxaIdsToAdd, bool be_silent);
    
    void setUpAllocator(int extra_parsimony_blocks, bool trackLikelihood=false,
                        int extra_lh_blocks=0) ;
    void prepareForPlacementRun();
    
    /** called before each batch */
    void prepareForBatch();
    
    /** cost potential placements for a batch
     @param candidates the taxa being added to the tree
     @param batchStart  the taxon index (index into candidates) for the first taxon in this batch
     @param batchStop one more than the taxon index for the last taxon in this batch
     @param refreshTime time spent refreshing parsimony (or likelihood)
                        information will be recorded against this TimeKeeper instance
     @param evaluationTime time spent evaluating placements will be recorded against this one */
    void doBatchPlacementCosting(TaxaToPlace& candidates,
                                   size_t batchStart, size_t batchStop,
                                   TargetBranchRange& targets );
    
    /** decide which, of the taxa in a batch are to be inserted, and move them
       to the front of the batch.
     @param candidates the taxa being added to the tree
     @param batchStart  the taxon index (index into candidates) for the first taxon in this batch
     @param batchStop one more than the taxon index for the last taxon in this batch
     @param[out] insertStop one more than the taon index for the last taxon to be inserted
     */
    void selectPlacementsForInsertion(TaxaToPlace& candidates,
                                     size_t  batchStart, size_t batchStop,
                                     size_t& insertStop);
    void startBatchInsert();
    
    void doBatchInsert(TaxaToPlace& candidates,
                       size_t  batchStart, size_t insertStop,
                       LikelihoodBlockPairs& spare_blocks,
                       double estimate_per_placement,
                       TargetBranchRange& targets,
                       ParsimonySearchParameters &s,
                       ParsimonyPathVector &pv);
    
    void doPartOfBatchInsert(TaxaToPlace& candidates,
                             std::vector<TaxonPlacement>& inserts,
                             size_t  insertStart, size_t insertStop,
                             LikelihoodBlockPairs& spare_blocks,
                             double estimate_per_placement,
                             TargetBranchRange& targets,
                             ParsimonySearchParameters &s,
                             ParsimonyPathVector &pv);
    
    /** insert a single taxon
     @param taxa the taxa being added to the tree
     @param taxon_index  the taxon index (index into taxa) for the taxon to be inserted
     @param targets the target branches. This is passed as a parameter, because when
     a taxon is inserted, that marks an existing target branch as used up, and adds three
         new target branches (for the new branches, on either side of the new interior node,
     and between the new interior node and the taxon's leaf node.
     @param blocks [todo: explain]
     */
    void insertTaxon(TaxaToPlace& taxa, size_t taxon_index,
                     TargetBranchRange& targets,
                     LikelihoodBlockPairs& blocks);
    
    /** called when a batch has been processed.
     @param taxa  the taxa being added to the tree
     @param start_taxon_index  the taxon index (index into candidates) for the first taxon in the batch
     @param stop_taxon_index  one more than the taxon index for the last taxon in the batch*/
    void doneBatch  (TaxaToPlace& taxa,
                     intptr_t start_taxon_index, intptr_t stop_taxon_index,
                     TargetBranchRange& targets);
    
    /** Remove taxa that have been inserted, from a TypedTaxaToPlace<T> container,
            and update them, remove target branches that no longer exist,
            from a TargetBranchRange&.
            @param candidates - the taxa (including some that have been inserted)
             @param batchStart - the index of the first taxon, in candidates, that
                        wasn't included in a batch, in the last pass (these go to
                        the front in the next pass, to ensure they will be included in
                        a batch during the next pass.  May be candidates.size().
             @param targets - the target branches (including some that no longer
                        exist, because taxa have been inserted into them).  */
    template <class T>
    void donePass(TypedTaxaToPlace<T>& candidates, size_t batchStart,
                  TargetBranchRange& targets) {
        optoTime.start();
        targets.removeUsed();
        //Remove all the candidates that we were able to place
        TypedTaxaToPlace<T> oldCandidates;
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
        inserts_per_batch = placement_params.getInsertsPerBatch(taxa_ids_to_add.size(), taxa_per_batch);
        TREE_LOG_LINE ( phylo_tree, VB_MAX, "At the end of this pass, index_lhs was "
                  << block_allocator->getLikelihoodBlockCount() << ", index_pars was "
                  << block_allocator->getParsimonyBlockCount());
        optoTime.stop();
    }

    /** Log tree structures *near* the newly added taxa*/
    void logSubtreesNearAddedTaxa() const;

    void donePlacement();
    void reportActivity() const;

    ~PlacementRun();
};

#endif /* placementrun_h */
