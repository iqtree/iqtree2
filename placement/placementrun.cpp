//
// placementrun.cpp
// Implementation of the PlacementRun class.
//
// Created by James Barbetti on 09-Oct-2020.
//

#include "placementrun.h"
#include "placementcostcalculator.h" //for

PlacementRun::PlacementRun(PhyloTree& tree): phylo_tree(tree)
    , block_allocator(nullptr)
    , costFunction(Placement::getCostFunction())
    , heuristic(SearchHeuristic::getSearchHeuristic())
    , taxon_placement_optimizer(TaxonPlacementOptimizer::getTaxonPlacementOptimizer())
    , batch_placement_optimizer(BatchPlacementOptimizer::getBatchPlacementOptimizer())
    , global_placement_optimizer(GlobalPlacementOptimizer::getGlobalPlacementOptimizer())
    , calculator(PlacementCostCalculator::getCostCalculator(costFunction)) {
}
void PlacementRun::setUpAllocator(int extra_parsimony_blocks, bool trackLikelihood,
                    int extra_lh_blocks) {
    int      index_parsimony        = 0;
    int      index_lh               = 0;
    phylo_tree.ensurePartialLHIsAllocated(extra_parsimony_blocks, extra_lh_blocks);
    phylo_tree.initializeAllPartialLh(index_parsimony, index_lh);
    block_allocator = trackLikelihood
        ? new LikelihoodBlockAllocator(phylo_tree, index_parsimony, index_lh)
        : new BlockAllocator(phylo_tree, index_parsimony);

}

PlacementRun::~PlacementRun() {
    delete global_placement_optimizer;
    delete batch_placement_optimizer;
    delete taxon_placement_optimizer;
    delete heuristic;
    delete calculator;
    delete block_allocator;
}

