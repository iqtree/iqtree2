//
// placementrun.cpp
// Implementation of the PlacementRun class.
//
// Created by James Barbetti on 09-Oct-2020.
//

#include "placementrun.h"
#include "placementcostcalculator.h" //for

PlacementRun::PlacementRun(PhyloTree& tree, const IntVector& taxaIdsToAdd)
    : phylo_tree(tree), taxa_ids_to_add(taxaIdsToAdd)
    , taxa_per_batch(Placement::getTaxaPerBatch(taxaIdsToAdd.size()))
    , inserts_per_batch(Placement::getInsertsPerBatch(taxa_ids_to_add.size(), taxa_per_batch))
    , block_allocator(nullptr)
    , heuristic(SearchHeuristic::getSearchHeuristic())
    , taxon_placement_optimizer(TaxonPlacementOptimizer::getNewTaxonPlacementOptimizer())
    , batch_placement_optimizer(BatchPlacementOptimizer::getNewBatchPlacementOptimizer())
    , global_placement_optimizer(GlobalPlacementOptimizer::getNewGlobalPlacementOptimizer())
    , calculator(PlacementCostCalculator::getNewCostCalculator())
    , taxa_inserted_this_batch(0), taxa_inserted_in_total(0), taxa_inserted_nearby(0) {
    if (calculator->usesLikelihood()) {
        phylo_tree.prepareToComputeDistances(); //Set up state look-up vectors
    }
}

void PlacementRun::setUpAllocator(int extra_parsimony_blocks, bool trackLikelihood,
                    int extra_lh_blocks) {
    int      index_parsimony        = 0;
    int      index_lh               = 0;
    phylo_tree.deleteAllPartialLh();
    phylo_tree.ensurePartialLHIsAllocated(extra_parsimony_blocks, extra_lh_blocks);
    phylo_tree.initializeAllPartialLh(index_parsimony, index_lh);
    block_allocator = trackLikelihood
        ? new LikelihoodBlockAllocator(phylo_tree, index_parsimony, index_lh)
        : new BlockAllocator(phylo_tree, index_parsimony);
}
                     
void PlacementRun::prepareForPlacementRun() {
    TREE_LOG_LINE ( phylo_tree, VB_MED, "After overallocating lh blocks, index_lh was "
              << block_allocator->getLikelihoodBlockCount() );
    if (VB_MED <= verbose_mode) {
        double curScore = phylo_tree.computeLikelihood();
        TREE_LOG_LINE ( phylo_tree, VB_MED, "Likelihood score before insertions was " << curScore );
        #if (0)
            curScore = phylo_tree.optimizeAllBranches(2);
            LOG_LINE ( phylo_tree, VB_MED, "Optimized likelihood score before insertions was " << curScore);
        #endif
    }
    TREE_LOG_LINE ( phylo_tree, VB_MED, "Batch size is " << taxa_per_batch
              << " and the number of inserts per batch is " << inserts_per_batch);
    if (calculator->usesSankoffParsimony()) {
        phylo_tree.computeTipPartialParsimony();
    }
}

void PlacementRun::prepareForBatch() {
    if (calculator->usesLikelihood()) {
        phylo_tree.clearAllPartialLH(false);
        phylo_tree.clearAllScaleNum(false);
        double likelihoodScore = phylo_tree.computeLikelihood();
        TREE_LOG_LINE( phylo_tree, VB_MIN, "Log-likelihood is currently " << likelihoodScore);
    }
    phylo_tree.clearAllPartialParsimony(false);
}

 double PlacementRun::doBatchPlacementCosting(TaxaToPlace& candidates,
                                              size_t batchStart, size_t batchStop,
                                              TargetBranchRange& targets) {
    ASSERT( !candidates.isEmpty() );
    ASSERT( !targets.empty() );
    
    double refreshTime = 0;
    size_t targetCount = targets.size();
    heuristic->prepareToFilter(phylo_tree, targets, 0, targetCount ,
                               candidates, batchStart, batchStop);

#if (NEW_TAXON_MAJOR)
     
    double refreshStart = getRealTime();
    LikelihoodBlockPairs blocks(2);
    for (int t=0; t<targets; ++t) {
        TargetBranch* target = targets.getTargetBranch(t);
        target->computeState(phylo_tree, t, blocks);
    }
    refreshTime += getRealTime() - refreshStart;

    for (int i = batchStart; i<batchStop; ++i) {
        TaxonToPlace& c = candidates.getTaxonByIndex(i);
        TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scoring ... " << c,taxonName);
        c.findPlacement(phylo_tree, i, targets,
                         heuristic, calculator);
        const PossiblePlacement& p = c.getBestPlacement();
        TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scored " << p.score << " for placement"
                      << " of " << c.taxonName << " with lengths "
                 << p.lenToNode1 << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon);
    }
     
#else //TARGET_BRANCH_MAJOR
    LikelihoodBlockPairs blocks(2);
    for (size_t t = 0; t < targetCount; ++t ) {
        TargetBranch* target = targets.getTargetBranch(t);
        TREE_LOG_LINE(phylo_tree, VB_DEBUG,
                      "Scoring target branch " << t << " of " << targetCount);
        
        double computeStart = getRealTime();
        target->computeState(phylo_tree, t, blocks);
        refreshTime += getRealTime() - computeStart;
        
        target->costPlacementOfTaxa(phylo_tree,
                                   targets,    t,
                                   candidates, batchStart, batchStop,
                                   heuristic,  calculator,
                                   t == 0 );
        
        double forgetStart = getRealTime();
        target->forgetState();
        refreshTime += getRealTime() - forgetStart;
    }
     
#endif
    heuristic->doneFiltering();
    return refreshTime;
}

void PlacementRun::selectPlacementsForInsertion(TaxaToPlace& candidates,
                                  size_t  batchStart, size_t batchStop,
                                  size_t& insertStop) {
    inserts_per_batch = Placement::getInsertsPerBatch(taxa_ids_to_add.size(), batchStop-batchStart);
    insertStop        = batchStart + inserts_per_batch;
    candidates.sortBatch( batchStart, batchStop);
    if (batchStop <= insertStop) {
        insertStop = batchStop; //Want them all
    }
}

void PlacementRun::startBatchInsert() {
    taxa_inserted_this_batch = 0;
}

void PlacementRun::insertTaxon(TaxaToPlace& taxa, size_t taxon_index,
                               TargetBranchRange& targets,
                               LikelihoodBlockPairs& blocks) {
    TaxonToPlace& c = taxa.getTaxonByIndex(taxon_index);
    const char* verb = "inserted";
    const char* where ;
    if (c.canInsert()) {
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        c.insertIntoTree(phylo_tree, *block_allocator, blocks, targets, *calculator);
        where = "on its desired branch";
    } else {
        //Another candidate taxon has gotten there first
        ++taxa_inserted_nearby;
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        c.insertNearby(phylo_tree, *block_allocator, blocks, targets, *calculator);
        where = "near its desired branch";
    }
    
    if (( verbose_mode >= VB_MIN && !phylo_tree.params->suppress_list_of_sequences)
        || verbose_mode >= VB_MED ) {
        const PossiblePlacement& p = c.getBestPlacement();
        stringstream s;
        s << taxa_inserted_in_total << ". " << verb << " "
            << c.taxonName << " " << where
        << " (branch index " << p.getTargetIndex() << "). It had ";
        if (!calculator->usesLikelihood()) {
            s << "parsimony score " << (int)(p.score);
        } else {
            s << "likelihood score " << p.score;
        }
        s << " (and path lengths " << p.lenToNode1
            << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon << ")";
        phylo_tree.logLine(s.str());
    }    
    taxon_placement_optimizer->cleanUpAfterTaxonPlacement(taxa, taxon_index, targets, phylo_tree);
}

void PlacementRun::doneBatch(TaxaToPlace& candidates,
                             size_t batchStart, size_t batchStop,
                             TargetBranchRange& targets) {
    if ( 1 < batchStop - batchStart ) {
        TREE_LOG_LINE ( phylo_tree, VB_MED,  "Inserted " << (taxa_inserted_this_batch)
                  << " out of a batch of " << (batchStop - batchStart) << "." );
    }
    batch_placement_optimizer->cleanUpAfterBatch(candidates, batchStart, batchStop, targets, phylo_tree);
    if (calculator->usesLikelihood()) {
        phylo_tree.fixNegativeBranch();
    }
    if (taxa_inserted_this_batch == 0) {
        outError("No taxa inserted in batch");
    }
}

void PlacementRun::logSubtreesNearAddedTaxa() const {
    PhyloNodeVector idToNode;
    this->phylo_tree.getArrayOfTaxaNodesById(nullptr, nullptr, idToNode);
    for (size_t i=0; i<taxa_ids_to_add.size(); ++i) {
        PhyloNode*     leaf      = idToNode[taxa_ids_to_add[i]];
        if ( leaf != nullptr ) {
            PhyloNeighbor* upLink    = leaf->firstNeighbor();
            PhyloNode*     upNode    = upLink->getNode();
            PhyloNeighbor* leftLink  = nullptr;
            PhyloNeighbor* rightLink = nullptr;
            double leftLength  = -1;
            double rightLength = -1;
            FOR_EACH_PHYLO_NEIGHBOR(upNode, leaf, itNei, nei) {
                if (leftLink==nullptr) {
                    leftLink   = nei;
                    leftLength = nei->length;
                } else {
                    rightLink   = nei;
                    rightLength = nei->length;
                }
            }
            auto length = leaf->firstNeighbor()->length;
            std::cout << (i+1) << "." << "Node [" << leaf->id << "]=" << leaf->name
                << " now has branch length " << length
                << " (interior left branch " << leftLength
                << " , and right branch " << rightLength << ")"
                << std::endl;
        }
    }
}
                     
PlacementRun::~PlacementRun() {
    if (calculator->usesLikelihood()) {
        phylo_tree.doneComputingDistances();
    }
    delete global_placement_optimizer;
    delete batch_placement_optimizer;
    delete taxon_placement_optimizer;
    delete heuristic;
    delete calculator;
    delete block_allocator;
}

