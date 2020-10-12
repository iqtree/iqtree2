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
void TaxonPlacementOptimizer::cleanUpAfterTaxonPlacement(const TaxonToPlace& taxon,
                                              PhyloTree* tree) {}

BatchPlacementOptimizer::BatchPlacementOptimizer() = default;
BatchPlacementOptimizer::~BatchPlacementOptimizer() = default;

GlobalPlacementOptimizer::GlobalPlacementOptimizer() = default;
GlobalPlacementOptimizer::~GlobalPlacementOptimizer() = default;
void GlobalPlacementOptimizer::cleanUpAfterPlacement(PhyloTree* tree) {
}

TaxonPlacementOptimizer* TaxonPlacementOptimizer::getTaxonPlacementOptimizer() {
    auto localCleanup = Placement::getLocalCleanupAlgorithm();
    switch (localCleanup) {
        default:
            break;
    }
    return new TaxonPlacementOptimizer();
}

BatchPlacementOptimizer* BatchPlacementOptimizer::getBatchPlacementOptimizer() {
    auto batchCleanup = Placement::getBatchCleanupAlgorithm();
    switch (batchCleanup) {
        default:
            break;
    }
    return new BatchPlacementOptimizer();
}

GlobalPlacementOptimizer* GlobalPlacementOptimizer::getGlobalPlacementOptimizer() {
    auto globalCleanup = Placement::getGlobalCleanupAlgorithm();
    switch (globalCleanup) {
        default:
            break;
    }
    return new GlobalPlacementOptimizer();
}


