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
    , costFunction(Placement::getCostFunction())
    , heuristic(SearchHeuristic::getSearchHeuristic())
    , taxon_placement_optimizer(TaxonPlacementOptimizer::getTaxonPlacementOptimizer())
    , batch_placement_optimizer(BatchPlacementOptimizer::getBatchPlacementOptimizer())
    , global_placement_optimizer(GlobalPlacementOptimizer::getGlobalPlacementOptimizer())
    , calculator(PlacementCostCalculator::getCostCalculator(costFunction))
    , taxa_inserted_this_batch(0), taxa_inserted_in_total(0), taxa_inserted_nearby(0) {
    if (calculator->usesLikelihood()) {
        phylo_tree.prepareToComputeDistances(); //Set up state look-up vectors
    }
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
    if (costFunction == Placement::SANKOFF_PARSIMONY) {
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

void PlacementRun::startBatchInsert() {
    taxa_inserted_this_batch = 0;
}

void PlacementRun::insertTaxon(TaxonToPlace& c, TargetBranchRange& targets) {
    const char* verb = "inserted";
    const char* where ;
    if (c.canInsert()) {
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        c.insertIntoTree(phylo_tree, block_allocator, targets, *calculator);
        where = "on its desired branch";
    } else {
        //Another candidate taxon has gotten there first
        ++taxa_inserted_nearby;
        ++taxa_inserted_this_batch;
        ++taxa_inserted_in_total;
        c.insertNearby(phylo_tree, block_allocator, targets, *calculator);
        where = "near its desired branch";
    }
    
    if (( verbose_mode >= VB_MIN && !phylo_tree.params->suppress_list_of_sequences)
        || verbose_mode >= VB_MED ) {
        stringstream s;
        s << taxa_inserted_in_total << ". " << verb << " "
            << c.taxonName << " " << where << ". It had ";
        const PossiblePlacement& p = c.getBestPlacement();
        if (costFunction==Placement::MAXIMUM_PARSIMONY ||
            costFunction==Placement::SANKOFF_PARSIMONY) {
            s << "parsimony score " << (int)(p.score);
        } else {
            s << "likelihood score " << p.score;
        }
        s << " (and path lengths " << p.lenToNode1
            << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon << ")";
        phylo_tree.logLine(s.str());
    }    
    taxon_placement_optimizer->cleanUpAfterTaxonPlacement(c, &phylo_tree);
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

