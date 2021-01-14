//
// placementoptimizer.cpp
// Implementations of the TaxonPlacementOptimizer,
// BatchPlacementOptimizer, and GlobalPlacementOptimizer classes.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "placement.h"
#include "placementoptimizer.h"

TaxonPlacementOptimizer::TaxonPlacementOptimizer() = default;
TaxonPlacementOptimizer::~TaxonPlacementOptimizer() = default;
void TaxonPlacementOptimizer::optimizeAfterTaxonPlacement(TaxaToPlace& taxa,
                                                         size_t taxon_index,
                                                         TargetBranchRange& targets,
                                                         PhyloTree& tree) {}

BatchPlacementOptimizer::BatchPlacementOptimizer() = default;
void BatchPlacementOptimizer::optimizeAfterBatch(TaxaToPlace& taxa,
                                                intptr_t firstTaxon, intptr_t lastTaxon,
                                                TargetBranchRange& targets,
                                                PhyloTree& tree) {
    if (VB_MIN <= verbose_mode) {
        std::stringstream s;
        s << "Processed batch of "
          << (lastTaxon - firstTaxon) << " taxa";
        tree.logLine(s.str() );
    }
}

BatchPlacementOptimizer::~BatchPlacementOptimizer() = default;

GlobalPlacementOptimizer::GlobalPlacementOptimizer() = default;
GlobalPlacementOptimizer::~GlobalPlacementOptimizer() = default;

void GlobalPlacementOptimizer::optimizeAfterPlacement(PhyloTree& tree) {
    tree.initializeTree(); //Make sure branch numbers et cetera are set.
    tree.deleteAllPartialLh();

    if (tree.isUsingSankoffParsimony() && !tree.params->sankoff_cost_file) {
        tree.deleteAllPartialParsimony();
        tree.stopUsingSankoffParsimony();
        tree.setParsimonyKernel(tree.params->SSE);
    }
    tree.determineBlockSizes();

    TREE_LOG_LINE(tree, VB_MED, 
        "Number of leaves " << tree.leafNum
        << ", of nodes "    << tree.nodeNum
        << ", of branches " << tree.branchNum);

    //First, recompute parsimony
    int parsimony_score = tree.computeParsimony("Computing parsimony (after adding taxa)");
    TREE_LOG_LINE(tree, VB_MIN, "Parsimony score after adding taxa was " << parsimony_score);

    //Then, fix any negative branches
    tree.initializeAllPartialLh();
    double negativeStart = getRealTime();
    tree.fixNegativeBranch();
    double negativeElapsed = getRealTime() - negativeStart;
    TREE_LOG_LINE(tree, VB_MED, "Fixing negative branches took " << negativeElapsed << " wall-clock seconds");

    //And, at last, compute the tree likelihood
    double likelyStart   = getRealTime();
    double likelihood    = tree.computeLikelihood();
    double likelyElapsed = getRealTime() - likelyStart;
    TREE_LOG_LINE(tree, VB_MIN, "Likelihood score after adding taxa was " << likelihood
        << " (and took " << likelyElapsed << " wall-clock seconds to calculate)");
}

GlobalLikelihoodPlacementOptimizer::GlobalLikelihoodPlacementOptimizer()  = default;

GlobalLikelihoodPlacementOptimizer::~GlobalLikelihoodPlacementOptimizer() = default;

void GlobalLikelihoodPlacementOptimizer::optimizeAfterPlacement(PhyloTree& tree) {
    super::optimizeAfterPlacement(tree);
    //Optimize
    double optimizeStart   = getRealTime();
    double score           = tree.optimizeAllBranches();
    double optimizeElapsed = getRealTime() - optimizeStart;
    TREE_LOG_LINE(tree, VB_MIN, "After optimizing"
        << " (which took " << optimizeElapsed << " wall-clock seconds)"
        << ", likelihood score was " << score);
}

TaxonPlacementOptimizer* TaxonPlacementOptimizer::getNewTaxonPlacementOptimizer() {
    //auto localCleanup = Placement::getLocalOptimizationAlgorithm();
    return new TaxonPlacementOptimizer();
}

BatchPlacementOptimizer* BatchPlacementOptimizer::getNewBatchPlacementOptimizer() {
    //auto batchCleanup = Placement::getBatchOptimizationAlgorithm();
    return new BatchPlacementOptimizer();
}

GlobalPlacementOptimizer* GlobalPlacementOptimizer::getNewGlobalPlacementOptimizer(bool useLikelihood) {
    //auto globalCleanup = Placement::getGlobalOptimizationAlgorithm();
    return useLikelihood ? new GlobalLikelihoodPlacementOptimizer()
                         : new GlobalPlacementOptimizer();
}


