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

struct PlacementRun {
public:
    PhyloTree&                phylo_tree;
    BlockAllocator*           block_allocator;
    Placement::CostFunction   costFunction; //parsimony or likelihood
    SearchHeuristic*          heuristic; //global? or localized?
    TaxonPlacementOptimizer*  taxon_placement_optimizer;
    BatchPlacementOptimizer*  batch_placement_optimizer;
    GlobalPlacementOptimizer* global_placement_optimizer;
    PlacementCostCalculator*  calculator;

    PlacementRun(PhyloTree& tree);
    void setUpAllocator(int extra_parsimony_blocks, bool trackLikelihood=false,
                        int extra_lh_blocks=0) ;
    ~PlacementRun();
};

#include <stdio.h>

#endif /* placementrun_h */
