//
// placementoptimizer.cpp
// Implementations of the TaxonPlacementOptimizer,
// BatchPlacementOptimizer, and GlobalPlacementOptimizer classes.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "placement.h"
#include "placementoptimizer.h"
#include <utils/timeutil.h>

TaxonPlacementOptimizer::TaxonPlacementOptimizer() = default;
TaxonPlacementOptimizer::~TaxonPlacementOptimizer() = default;
void TaxonPlacementOptimizer::optimizeAfterTaxonPlacement(TaxaToPlace& taxa,
                                                         size_t taxon_index,
                                                         TargetBranchRange& targets,
                                                         PhyloTree& tree) {}

BatchPlacementOptimizer::BatchPlacementOptimizer(bool shut_up): be_quiet(shut_up) {}
void BatchPlacementOptimizer::optimizeAfterBatch(TaxaToPlace& taxa,
                                                intptr_t firstTaxon, intptr_t lastTaxon,
                                                TargetBranchRange& targets,
                                                PhyloTree& tree) {
    if (VB_MIN <= verbose_mode && !be_quiet) {
        std::stringstream s;
        s << "Processed batch of "
          << (lastTaxon - firstTaxon) << " taxa";
        tree.logLine(s.str() );
    }
}

BatchPlacementOptimizer::~BatchPlacementOptimizer() = default;

GlobalPlacementOptimizer::GlobalPlacementOptimizer(bool be_silent) : be_quiet(be_silent) {}
GlobalPlacementOptimizer::~GlobalPlacementOptimizer() = default;

void GlobalPlacementOptimizer::optimizeAfterPlacement(PhyloTree& tree) {
    tree.initializeTree(); //Make sure branch numbers et cetera are set.
    tree.deleteAllPartialLh();

    bool turnOffSankoff = tree.isUsingSankoffParsimony()
                          && !tree.params->sankoff_cost_file;
    tree.deleteAllPartialParsimony();
    if (turnOffSankoff) {
        tree.stopUsingSankoffParsimony();
        tree.setParsimonyKernel(tree.params->SSE);
    }
    tree.determineBlockSizes();

    if (tree.params->compute_ml_dist && tree.hasModel() ) {
        //First, recompute parsimony
        int parsimony_score = tree.computeParsimony("Computing parsimony (after adding taxa)",
                                                    true, false);
        if (!be_quiet) {
            TREE_LOG_LINE(tree, VB_MIN, "Parsimony score after adding taxa was " << parsimony_score);
        }
        //Then, fix any negative branches
        tree.initializeAllPartialLh();
        double negativeStart = getRealTime();
        auto   fixed = tree.fixNegativeBranch();
        double negativeElapsed = getRealTime() - negativeStart;
        TREE_LOG_LINE(tree, VB_MED, "Fixing " << fixed
                      << " negative branches took " << negativeElapsed
                      << " wall-clock seconds");

        //And, at last, compute the tree likelihood
        double likelyStart   = getRealTime();
        double likelihood    = tree.computeLikelihood();
        double likelyElapsed = getRealTime() - likelyStart;
        TREE_LOG_LINE(tree, VB_MIN, "Likelihood score after adding taxa was " << likelihood
            << " (and took " << likelyElapsed << " wall-clock seconds to calculate)");
    } else {
        double score;
        tree.setAllBranchLengthsFromParsimony(true, score);
        if (!be_quiet) {
            TREE_LOG_LINE(tree, VB_MIN, "Parsimony score after adding taxa was " << score);
        }
    }
}

GlobalLikelihoodPlacementOptimizer::GlobalLikelihoodPlacementOptimizer
    (bool be_silent): super(be_silent) {}

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

BatchPlacementOptimizer* BatchPlacementOptimizer::getNewBatchPlacementOptimizer(bool be_quiet) {
    //auto batchCleanup = Placement::getBatchOptimizationAlgorithm();
    return new BatchPlacementOptimizer(be_quiet);
}

GlobalPlacementOptimizer* GlobalPlacementOptimizer::getNewGlobalPlacementOptimizer
    (bool useLikelihood, bool be_silent) {
    return useLikelihood ? new GlobalLikelihoodPlacementOptimizer(be_silent)
                         : new GlobalPlacementOptimizer(be_silent);
}


