//
// placementrun.cpp
// Implementation of the PlacementRun class.
//
// Created by James Barbetti on 09-Oct-2020.
//

#include "placementrun.h"
#include "placementcostcalculator.h"
#include <utils/timeutil.h>                 //for getRealTime()
#include <tree/phylotreethreadingcontext.h> //

TaxonPlacement::TaxonPlacement(): candidate_index(0),
    target_index(0), placement(nullptr) {
}

bool TaxonPlacement::operator < ( const TaxonPlacement& rhs) const {
    if (target_index < rhs.target_index) return true;
    if (rhs.target_index < target_index) return false;
    return (*placement < *rhs.placement);
}

bool TaxonPlacement::operator <= ( const TaxonPlacement& rhs) const {
    if (target_index < rhs.target_index) return true;
    if (rhs.target_index < target_index) return false;
    return (*placement <= *rhs.placement);
}

PlacementRun::PlacementRun(PhyloTree& tree, const PlacementParameters& parameters,
                           const IntVector& taxaIdsToAdd, bool be_silent)
    : phylo_tree(tree), placement_params(parameters)
    , taxa_ids_to_add(taxaIdsToAdd), be_quiet(be_silent)
    , taxa_per_batch(placement_params.getTaxaPerBatch(taxaIdsToAdd.size()))
    , inserts_per_batch(placement_params.getInsertsPerBatch(taxa_ids_to_add.size(), taxa_per_batch))
    , block_allocator(nullptr)
    , heuristic(SearchHeuristic::getSearchHeuristic(placement_params))
    , calculator(PlacementCostCalculator::getNewCostCalculator(placement_params))
    , use_likelihood(heuristic->usesLikelihood() || calculator->usesLikelihood())
    , taxon_placement_optimizer(TaxonPlacementOptimizer::getNewTaxonPlacementOptimizer())
    , batch_placement_optimizer(BatchPlacementOptimizer::getNewBatchPlacementOptimizer(be_silent))
    , global_placement_optimizer(GlobalPlacementOptimizer::getNewGlobalPlacementOptimizer(use_likelihood, be_silent))
    , taxa_inserted_this_batch(0), taxa_inserted_in_total(0), taxa_inserted_nearby(0)
    , overall      ("adding new taxa"), initializing ("initializing")
    , refreshTime  ("Refresh"), searchTime ("Search"), rankingTime  ("Ranking")
    , insertTime   ("Insert"),  optoTime   ("Post-Batch Optimization") {
    if (use_likelihood) {
        phylo_tree.prepareToComputeDistances(); //Set up state look-up vectors
    }
}

void PlacementRun::setUpAllocator(int extra_parsimony_blocks,
                                  bool trackLikelihood,
                                  int extra_lh_blocks) {
    int      index_parsimony = 0;
    int      index_lh        = 0;
    phylo_tree.deleteAllPartialLhAndParsimony(); 
    if (trackLikelihood) {
        phylo_tree.ensurePartialLHIsAllocated(extra_parsimony_blocks, extra_lh_blocks);
        phylo_tree.initializeAllPartialLh(index_parsimony, index_lh);
        block_allocator = new LikelihoodBlockAllocator(phylo_tree, index_parsimony, index_lh);
    } else {
        phylo_tree.ensureCentralPartialParsimonyIsAllocated(extra_parsimony_blocks);
        phylo_tree.initializeAllPartialPars(index_parsimony);
        block_allocator = new BlockAllocator(phylo_tree, index_parsimony);
    }
}
                     
void PlacementRun::prepareForPlacementRun() {
    if (use_likelihood) {
        TREE_LOG_LINE(phylo_tree, VB_MED, "After overallocating lh blocks, index_lh was "
            << block_allocator->getLikelihoodBlockCount());
    }
    if (VB_MED <= verbose_mode &&
        (phylo_tree.params->compute_likelihood || use_likelihood) &&
        phylo_tree.hasModel()) {
        phylo_tree.configureLikelihoodKernel(*phylo_tree.params, true);
        double curScore = phylo_tree.computeLikelihood();
        TREE_LOG_LINE(phylo_tree, VB_MED, "Likelihood score before insertions was " << curScore);
        if (use_likelihood) {
            curScore = phylo_tree.optimizeAllBranches(2);
            TREE_LOG_LINE(phylo_tree, VB_MED, "Optimized likelihood score"
                          " before insertions was " << curScore);
        }
    }
    if (!use_likelihood) {
        phylo_tree.deleteAllPartialLh();
    }
    if (!be_quiet) {
        TREE_LOG_LINE ( phylo_tree, VB_MED, "Batch size is " << taxa_per_batch
                        << " and the number of inserts per batch is " << inserts_per_batch);
    }
    if (calculator->usesSankoffParsimony()) {
        phylo_tree.computeTipPartialParsimony();
    }
}

void PlacementRun::prepareForBatch() {
    refreshTime.start();
    if (calculator->usesLikelihood()) {
        phylo_tree.clearAllPartialLH(false);
        phylo_tree.clearAllScaleNum(false);
        double likelihoodScore = phylo_tree.computeLikelihood();
        TREE_LOG_LINE( phylo_tree, VB_MIN, "Log-likelihood is currently " << likelihoodScore);
    }
    phylo_tree.clearAllPartialParsimony(false);
    refreshTime.stop();
}

void PlacementRun::doBatchPlacementCosting(TaxaToPlace& candidates,
                                             size_t batchStart, size_t batchStop,
                                             TargetBranchRange& targets) {
    ASSERT( !candidates.isEmpty() );
    ASSERT( !targets.empty() );
    
    size_t targetCount = targets.size();
    heuristic->prepareToFilter(phylo_tree, targets, 0, targetCount ,
                               candidates, batchStart, batchStop);
    double parsimony_score = -1;

#if (NEW_TAXON_MAJOR)

    refreshTime.start();
    LikelihoodBlockPairs blocks(2);
    for (int t=0; t<targets; ++t) {
        TargetBranch* target = targets.getTargetBranch(t);
        target->computeState(phylo_tree, parsimony_score, t, blocks);
    }
    refreshTime.stop();
    
    evaluationTime.start();
    for (int i = batchStart; i<batchStop; ++i) {
        TaxonToPlace& c = candidates.getTaxonByIndex(i);
        TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scoring ... " << c,taxonName);
        c.findPlacement(phylo_tree, i, targets,
                         heuristic, calculator);
        const PossiblePlacement& p = c.getBestPlacement();
        TREE_LOG_LINE(phylo_tree, VB_DEBUG, "Scored " << p.score << " for placement"
                      << " of " << c.taxonName << " with lengths "
                 << p.lenToNode1 << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon);
        if ((i-batchStart) % 1000 ) == 0 ) {
            phylo_tree.trackProgress(1000.0);
        }
    }
    evaluationTime.stop();
#else //TARGET_BRANCH_MAJOR
    LikelihoodBlockPairs blocks(2);
    size_t tStep = 1;
    while (tStep * 100 <= targetCount) {
        tStep *= 10;
    }
    double workStep = (double)( batchStop - batchStart ) / (double)targetCount;
    double workDone = 0;
    
    for (size_t t = 0; t < targetCount; ++t ) {
        TargetBranch* target = targets.getTargetBranch(t);
        if ((t & tStep) == 0) {
            TREE_LOG_LINE(phylo_tree, VB_DEBUG,
                "Scoring target branch " << t << " of " << targetCount);
        }
        refreshTime.start();
        target->computeState(phylo_tree, parsimony_score, t, blocks);
        refreshTime.stop();
        
        searchTime.start();
        target->costPlacementOfTaxa(phylo_tree,
                                   targets,    t,
                                   candidates, batchStart, batchStop,
                                   heuristic,  calculator,
                                   t == 0 );
        searchTime.stop();
        
        workDone += workStep;
        if (workDone > 1000.0) {
            workDone -= 1000.0;
            phylo_tree.trackProgress(1000.0);
        }
    }
#endif
    heuristic->doneFiltering();
}

void PlacementRun::selectPlacementsForInsertion(TaxaToPlace& candidates,
                                  size_t  batchStart, size_t batchStop,
                                  size_t& insertStop) {
    rankingTime.start();
    inserts_per_batch = placement_params.getInsertsPerBatch(taxa_ids_to_add.size(), batchStop-batchStart);
    insertStop        = batchStart + inserts_per_batch;
    candidates.sortBatch( batchStart, batchStop);
    if (batchStop <= insertStop) {
        insertStop = batchStop; //Want them all
    }
    rankingTime.stop();
}

void PlacementRun::startBatchInsert() {
    taxa_inserted_this_batch = 0;
}

void PlacementRun::doBatchInsert
        (TaxaToPlace& candidates, size_t  insertStart, size_t insertStop,
         LikelihoodBlockPairs& spare_blocks, double estimate_per_placement,
         TargetBranchRange& targets,
         ParsimonySearchParameters &s, ParsimonyPathVector &pv) {
    insertTime.start();
    std::vector<TaxonPlacement> inserts;
    inserts.resize(insertStop-insertStart);
    for ( size_t i = insertStart; i < insertStop; ++i) {
        TaxonPlacement& insert = inserts[i-insertStart];
        insert.candidate_index = i;
        insert.placement       = &candidates.getTaxonByIndex(i).getBestPlacement();
        insert.target_index    = insert.placement->getTargetIndex();
    }
    insertTime.stop(); //includes optimisation time

    rankingTime.start();
    std::sort(inserts.begin(), inserts.end() );
    rankingTime.stop();

    doPartOfBatchInsert(candidates, inserts, 0, inserts.size(),
                        spare_blocks, estimate_per_placement,
                        targets, s, pv);
}

void PlacementRun::doPartOfBatchInsert
        (TaxaToPlace& candidates, std::vector<TaxonPlacement>& inserts,
         size_t  insertStart, size_t insertStop,
         LikelihoodBlockPairs& spare_blocks, double estimate_per_placement,
         TargetBranchRange& targets,
         ParsimonySearchParameters &s, ParsimonyPathVector &pv) {
    size_t j = insertStart;
    while ( j < insertStop ) {
        size_t h = j;
        size_t t = inserts[h].target_index;
        for (++j; j<insertStop; ++j) {
            if (t != inserts[j].target_index) {
                break;
            }
        }
        
        if (targets[t].isOutOfDate()) {
#if (0)
            TimeKeeper icky2("rescoring target");
            icky2.start();
#endif
            double parsimony_score = -1.0;
            targets[t].computeState(phylo_tree, parsimony_score,
                                    t, spare_blocks);
            ASSERT(targets[t].getLeftNeighbor()->isParsimonyComputed());
            ASSERT(targets[t].getRightNeighbor()->isParsimonyComputed());
#if (0)
            icky2.stop();
            TREE_LOG_LINE(phylo_tree, VB_MIN,
                          "Rescoring target " << t
                          << " for insert " << h << " took "
                          << icky2.elapsed_wallclock_time << " wall, "
                          << icky2.elapsed_cpu_time << " cpu");
        } else {
            TREE_LOG_LINE(phylo_tree, VB_MIN, "Target " << t
                          << " for insert " << h << " already up to date ");
#endif
        }
        
        size_t insertCount = j - h;
        if (insertCount<256) {
            //Insert candidates h through j-1
            insertTime.start();
#if (0)
            TimeKeeper icky("inserting");
            icky.start();
#endif
            for (size_t i=h; i<j; ++i) {
                TaxonPlacement& insert    = inserts[i];
#if (0)
                TaxonToPlace&   candidate = candidates.getTaxonByIndex
                                            (insert.candidate_index);
                auto check_index = candidate.getBestPlacement().getTargetIndex();
#endif
                insertTaxon(candidates, insert.candidate_index,
                            targets, spare_blocks);
                if ((taxa_inserted_in_total % 1000) == 0) {
                    phylo_tree.trackProgress(1000.0*estimate_per_placement);
                }
#if (0)
                TREE_LOG_LINE(phylo_tree, VB_MIN, "Insert " << i
                              << ", of Candidate " << insert.candidate_index
                              << ", at or near target branch " << insert.target_index
                              << " (check " << check_index << ") done");
#endif
            }
#if (0)
            icky.stop();
            TREE_LOG_LINE(phylo_tree, VB_MIN, "Inserting I " << h
                          << " thru " << (j-1)
                          << " at T " << inserts[h].target_index << " took "
                          << icky.elapsed_wallclock_time << " wall, "
                          << icky.elapsed_cpu_time << " cpu");
#endif
            insertTime.stop();
        } else {
            //The recursive bit!
            size_t sampleCount = (size_t)floor(sqrt(insertCount));
            size_t sampleStop  = h + sampleCount;
            double inside_estimate = (double)insertCount
                                   / (double)(sampleCount + insertCount)
                                   * estimate_per_placement;
            
            doPartOfBatchInsert(candidates, inserts, h, sampleStop,
                                spare_blocks, inside_estimate,
                                targets, s, pv);
            for (size_t i=sampleStop; i<j; ++i) {
                TaxonPlacement& insert    = inserts[i];
                TaxonToPlace&   candidate = candidates.getTaxonByIndex
                                            (insert.candidate_index);
                candidate.findNewPlacement(phylo_tree, *block_allocator,
                                           spare_blocks, targets, *calculator);
                insert.placement          = &candidate.getBestPlacement();
                insert.target_index       = insert.placement->getTargetIndex();
#if (0)
                TREE_LOG_LINE(phylo_tree, VB_MIN, "New place for insert " << i
                              << " of candidate " << insert.candidate_index
                              << " was target branch " << insert.target_index
                              << " with pars score " << insert.placement->parsimony_score);
#endif
            }
            rankingTime.start();
            std::sort(inserts.begin() + sampleStop, inserts.begin() + j );
#if (0)
            for (size_t i=sampleStop; i<j; ++i) {
                TREE_LOG_LINE(phylo_tree, VB_MIN, "Insert " << i
                              << " candidate " << inserts[i].candidate_index
                              << " has best placement t.b. " << inserts[i].target_index);
            }
#endif
            rankingTime.stop();
            doPartOfBatchInsert(candidates, inserts, sampleStop, j,
                                spare_blocks, inside_estimate,
                                targets, s, pv);
        }
        if (1<j-h) {
            optoTime.start();
            //Run SPR optimization on all of the taxa that
            //targeted the same branch (between front and back)
            auto max_out_threads = phylo_tree.params->parsimony_uses_max_threads;
            PhyloTreeThreadingContext context(phylo_tree, max_out_threads);
            phylo_tree.optimizePlacementRegion(s, targets, t, pv,
                                               context, spare_blocks);
            optoTime.stop();
        }
    }
}

void PlacementRun::insertTaxon(TaxaToPlace& taxa, size_t taxon_index,
                               TargetBranchRange& targets,
                               LikelihoodBlockPairs& blocks) {
    TaxonToPlace& c = taxa.getTaxonByIndex(taxon_index);
    const char* verb = "inserted";
    const char* where ;
    if (c.canInsert()) {
        c.insertIntoTree(phylo_tree, *block_allocator,
                         blocks, targets, *calculator);
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        where = "on its desired branch";
    } else {
        //Another candidate taxon has gotten there first
        c.insertNearby(phylo_tree, *block_allocator,
                       blocks, targets, *calculator);
        ++taxa_inserted_nearby;
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        where = "near its desired branch";
    }
    
    if (( verbose_mode >= VB_MIN
          && !phylo_tree.params->suppress_list_of_sequences)
          || verbose_mode >= VB_MAX ) {
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
    taxon_placement_optimizer->optimizeAfterTaxonPlacement(taxa, taxon_index,
                                                           targets, phylo_tree);
}

void PlacementRun::doneBatch(TaxaToPlace& candidates,
                             intptr_t batchStart, intptr_t batchStop,
                             TargetBranchRange& targets) {
    optoTime.start();
    intptr_t batch_size = batchStop - batchStart;
    if ( static_cast<intptr_t>(taxa_inserted_this_batch) < batch_size && !be_quiet ) {
        TREE_LOG_LINE ( phylo_tree, VB_MED,  "Inserted " << (taxa_inserted_this_batch)
                  << " out of a batch of " << batch_size << "." );
    }
    batch_placement_optimizer->optimizeAfterBatch(candidates, batchStart, batchStop, targets, phylo_tree);
    if (calculator->usesLikelihood()) {
        phylo_tree.fixNegativeBranch();
    }
    if (taxa_inserted_this_batch == 0) {
        outError("No taxa inserted in batch");
    }
    optoTime.stop();
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
            double leftLength  = -1;
            double rightLength = -1;
            FOR_EACH_PHYLO_NEIGHBOR(upNode, leaf, itNei, nei) {
                if (leftLink==nullptr) {
                    leftLink   = nei;
                    leftLength = nei->length;
                } else {
                    rightLength = nei->length;
                }
            }
            auto length = leaf->firstNeighbor()->length;
            if (VB_MED <= verbose_mode) {
                std::cout << (i+1) << "." << "Node [" << leaf->id << "]=" << leaf->name
                    << " now has branch length " << length
                    << " (interior left branch " << leftLength
                    << " , and right branch " << rightLength << ")"
                    << std::endl;
            }
        }
    }
}

void PlacementRun::donePlacement() {
    optoTime.start();
    if (VB_MED <= verbose_mode && use_likelihood) {
        logSubtreesNearAddedTaxa();
    }
    if (!be_quiet) {
        TREE_LOG_LINE(phylo_tree, VB_MED, "Tidying up tree after inserting taxa.");
    }
    global_placement_optimizer->optimizeAfterPlacement(phylo_tree);
    if (VB_MED <= verbose_mode && use_likelihood) {
        logSubtreesNearAddedTaxa();
    }
    optoTime.stop();
}

void PlacementRun::reportActivity() const {
    phylo_tree.hideProgress();
    std::cout.precision(4);
    refreshTime.report();
    searchTime.report();
    insertTime.report();
    optoTime.report();
    phylo_tree.showProgress();
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

