//
// taxontoplace.cpp
// Implementation of the TaxonToPlace and LessFussyTaxon
// classes.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "taxontoplace.h"

#include "blockallocator.h"
#include "placementcostcalculator.h"
#include "searchheuristic.h"

TaxonToPlace::TaxonToPlace() : bestPlacement(), taxonId(-1), taxonName()
               , inserted(false), new_leaf(nullptr), new_interior(nullptr)
               , partial_pars(nullptr) {
}

TaxonToPlace::TaxonToPlace(const TaxonToPlace& rhs) = default; //copy everything!

TaxonToPlace::TaxonToPlace(BlockAllocator* ba, int interior_id, int leaf_id,
                           std::string name, bool delay) {
    initialize(ba, interior_id, leaf_id, name, delay);
}

void TaxonToPlace::initialize(BlockAllocator* ba, int interior_id, int leaf_id,
                              std::string name, bool delayCompute) {
    PhyloTree& phylo_tree = ba->getTree();
    inserted     = false;
    partial_lh   = nullptr;
    scale_num    = nullptr;
    taxonId      = leaf_id;
    taxonName    = name;
    new_leaf     = phylo_tree.newNode(taxonId, taxonName.c_str());
    new_interior = phylo_tree.newNode(interior_id);
    new_interior->is_floating_interior = true;
    new_interior->addNeighbor(new_leaf, -1 );
    new_leaf->addNeighbor(new_interior, -1 );
    PhyloNeighbor* nei = new_interior->firstNeighbor();
    ba->allocateMemoryFor(nei); //The allocator knows if partial_lh & scale_num wanted
    partial_pars = nei->partial_pars;
    partial_lh   = nei->partial_lh;
    scale_num    = nei->scale_num;
    if (!delayCompute) {
        phylo_tree.computePartialParsimony(nei, new_interior);
    }
}

void TaxonToPlace::computeParsimony(PhyloTree* tree) {
    //Assumes new_leaf is the first neighbor of new_interior
    PhyloNeighbor* nei = new_interior->firstNeighbor();
    tree->computePartialParsimony(nei, new_interior);
}

TaxonToPlace::~TaxonToPlace() = default;
    
const UINT* TaxonToPlace::getParsimonyBlock() const {
    return partial_pars;
}

const double* TaxonToPlace::getLikelihoodBlock() const {
    return partial_lh;
}

const UBYTE* TaxonToPlace::getScaleNumBlock() const {
    return scale_num;
}

int TaxonToPlace::getTaxonId() const {
    return taxonId;
}

size_t TaxonToPlace::considerPlacements(const PossiblePlacement* placements,
                                  size_t count) {
    size_t bestI = 0;
    for (size_t i=1; i<count; ++i) {
        if (placements[bestI].score < placements[i].score ) {
            bestI = i;
        }
    }
    bestPlacement = placements[bestI];
    return bestI;
}

bool TaxonToPlace::considerAdditionalPlacement(const PossiblePlacement& placement) {
    bool best = (!bestPlacement.canStillUse()
                 || placement.score < bestPlacement.score);
    #if (0)
        std::cout << "  xplace " << placement.getTargetIndex()
                  << " score "   << placement.score << std::endl;
    #endif
    if (best) {
        bestPlacement = placement;
    }
    return best;
}

const PossiblePlacement& TaxonToPlace::getBestPlacement() const {
    return bestPlacement;
}

bool TaxonToPlace::canInsert() const {
    return bestPlacement.canStillUse();
}

void TaxonToPlace::findPlacement(PhyloTree& phylo_tree,
                                 size_t taxonIndex,
                                 TargetBranchRange& range,
                                 SearchHeuristic* heuristic,
                                 const PlacementCostCalculator* calculator) {
    //
    //Note: for now, heuristic is ignored...
    //
    size_t target_branch_count = range.size();
    std::vector<PossiblePlacement> placements;
    placements.resize(target_branch_count);
    PossiblePlacement* placementArray = placements.data();
    
    if ( verbose_mode >= VB_DEBUG ) {
        std::stringstream s1;
        s1 << "Scoring " << this->taxonName;
        phylo_tree.logLine(s1.str());
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<target_branch_count; ++i) {
        TargetBranchRef target(&range, i);
        if (heuristic->isPlacementWorthTrying(*this, taxonIndex, target)) {
            PossiblePlacement* p = placementArray + i;
            p->setTargetBranch(&range, i);
            calculator->assessPlacementCost(phylo_tree, *this, *p);
        }
    }
    phylo_tree.trackProgress(static_cast<double>(target_branch_count));
    
    auto bestI = considerPlacements(placements.data(), placements.size());
    
    if ( verbose_mode >= VB_MED ) {
        std::stringstream s3;
        s3 << "Best (lowest) score for " << this->taxonName 
            << " was " << bestPlacement.score << " at place " << bestI;
        phylo_tree.logLine(s3.str());
    }
    inserted      = false;
}
bool TaxonToPlace::operator < (const TaxonToPlace& rhs) const {
    return bestPlacement.score < rhs.bestPlacement.score; //lowest score is best
}
bool TaxonToPlace::operator <= (const TaxonToPlace& rhs) const {
    return bestPlacement.score <= rhs.bestPlacement.score; //lowest score is best
}
void TaxonToPlace::insertIntoTree(PhyloTree& phylo_tree, BlockAllocator& b,
                    LikelihoodBlockPairs& blocks, TargetBranchRange& dest,
                    PlacementCostCalculator& calculator) {
    //
    //Assumes, canInsert() returned true, and the tree has not
    //been modified in the meantime.
    //
    TargetBranch* target = bestPlacement.getTarget();
    PhyloNode* node_1    = const_cast<PhyloNode*>(target->first);
    PhyloNode* node_2    = const_cast<PhyloNode*>(target->second);
    PhyloNeighbor* down  = new_interior->findNeighbor(new_leaf);
    down->length         = bestPlacement.lenToNewTaxon;
    PhyloNeighbor* up    = new_leaf->findNeighbor(new_interior);
    up->length           = bestPlacement.lenToNewTaxon;
    target->handOverComputedStateTo(up);

    new_interior->addNeighbor(node_2, bestPlacement.lenToNode2 );
    b.handOverComputedState  (node_1->findNeighbor(node_2), new_interior->findNeighbor(node_2) );
    node_1->updateNeighbor   (node_2, new_interior, bestPlacement.lenToNode1);
    node_1->findNeighbor     (new_interior)->clearComputedFlags();

    new_interior->addNeighbor(node_1, bestPlacement.lenToNode1 );
    b.handOverComputedState  (node_2->findNeighbor(node_1), new_interior->findNeighbor(node_1) );
    node_2->updateNeighbor   (node_1, new_interior, bestPlacement.lenToNode2);
    node_1->findNeighbor     (new_interior)->clearComputedFlags();

    inserted                    = true;
    double score                = -1.0;
    bool   lh_needed            = calculator.usesLikelihood();
    ReplacementBranchList* reps = new ReplacementBranchList;
    reps->emplace_back ( dest.addNewRef(b, blocks, new_interior, node_1  , score, lh_needed) );
    reps->emplace_back ( dest.addNewRef(b, blocks, new_interior, node_2  , score, lh_needed) );
    reps->emplace_back ( dest.addNewRef(b, blocks, new_interior, new_leaf, score, lh_needed) );
    bestPlacement.target_branch.getTarget()->takeOwnershipOfReplacementVector(reps);
    //Set this information in case it's logged higher up
    //in the call stack (because... it might well be).
    bestPlacement.lenToNode1    = new_interior->findNeighbor(node_1)->length;
    bestPlacement.lenToNode2    = new_interior->findNeighbor(node_2)->length;
    bestPlacement.lenToNewTaxon = new_interior->findNeighbor(new_leaf)->length;
    
    FOR_EACH_ADJACENT_PHYLO_NODE(new_interior, nullptr, it, node) {
        ASSERT(node->findNeighbor(new_interior)->isParsimonyComputed());
        ASSERT(new_interior->findNeighbor(node)->isParsimonyComputed());
    }
    
    ++phylo_tree.leafNum;
    phylo_tree.branchNum += 2;
    phylo_tree.nodeNum   += 2;
}
void TaxonToPlace::forgetGazumpedPlacements() {
    bestPlacement.forget();
}
void TaxonToPlace::findNewPlacement(PhyloTree& phylo_tree, BlockAllocator& b,
                                     LikelihoodBlockPairs& blocks,
                                     TargetBranchRange& dest,
                                     PlacementCostCalculator& calculator ) {
    const TargetBranch* blocked_target = bestPlacement.getTarget();
    forgetGazumpedPlacements();
    std::vector<PossiblePlacement> placements;
    assessNewTargetBranches(phylo_tree, calculator,
                            blocked_target, placements);
    for (auto it = placements.begin(); it != placements.end(); ++it) {
        considerAdditionalPlacement(*it);
    }
}
bool TaxonToPlace::insertNearby(PhyloTree& phylo_tree, BlockAllocator& b,
                                LikelihoodBlockPairs& blocks,
                                TargetBranchRange& dest,
                                PlacementCostCalculator& calculator ) {
    findNewPlacement(phylo_tree, b, blocks, dest, calculator);
    if (!canInsert()) {
        return false;
    }
    insertIntoTree(phylo_tree, b, blocks, dest, calculator);
    return true;
}
void TaxonToPlace::assessNewTargetBranches(PhyloTree& phylo_tree,
                                           PlacementCostCalculator& calculator,
                                           const TargetBranch* tb,
                                           std::vector<PossiblePlacement>& scores) {
    if (tb==nullptr  ) {
        return;
    }
    auto reps = tb->getReplacements();
    if (reps==nullptr) {
        return;
    }
    std::vector< ReplacementBranchList* > stack;
    stack.push_back(reps);
    LikelihoodBlockPairs blocks(2);
    double parsimony_score = -1.0;
    while (!stack.empty()) {
        reps = stack.back();
        stack.pop_back();
        for (auto it=reps->begin(); it!=reps->end(); ++it) {
            if (it->isUsedUp()) {
                stack.push_back(it->getTarget()->getReplacements());
            } else {
                PossiblePlacement p;
                p.setTargetBranch(*it);
                TargetBranch* tb = it->getTarget();
                ASSERT(tb->stillExists());
                tb->computeState(phylo_tree, parsimony_score,
                                              p.getTargetIndex(), blocks);
                scores.emplace_back(p);
            }
        }
    }
    intptr_t score_count = scores.size();
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t t = 0 ; t < score_count; ++ t) {
        calculator.assessPlacementCost(phylo_tree, *this, scores[t]);
    }
}

LessFussyTaxon::LessFussyTaxon() : super() { }
LessFussyTaxon::LessFussyTaxon(const LessFussyTaxon& rhs) = default; //copy everything!
LessFussyTaxon::LessFussyTaxon(BlockAllocator* ba, int interior_id, int leaf_id,
                               std::string name, bool delay)
                             : super(ba, interior_id, leaf_id, name, delay) {
}
size_t LessFussyTaxon::considerPlacements(const PossiblePlacement* placements,
                                          size_t count) {
    placementStore.clear();
    size_t            bestI = 0;
    for (int i=0; i < count; ++i) {
        considerAdditionalPlacement(placements[i]);
        if (i==0 || placements[i] < placements[bestI]) {
            bestI    = i;
        }
    }
    return bestI;
}

bool LessFussyTaxon::considerAdditionalPlacement(const PossiblePlacement& placement) {
    size_t placement_count = placementStore.size();
    if (max_placements_to_keep <= placement_count) {
        if (placementStore.back() < placement) {
            return false;
        }
        else {
            placementStore.pop_back();
        }
    }
    placementStore.emplace_back(PossiblePlacement());
    size_t sweep = placementStore.size() - 1;
    while ( 0 < sweep && placement < placementStore[sweep-1]) {
        placementStore[sweep] = placementStore[sweep - 1];
        --sweep;
    }
    placementStore[sweep] = placement;
    if (sweep == 0) {
        bestPlacement = placement;
        return true;
    }
    return false;
}
void LessFussyTaxon::forgetGazumpedPlacements() {
    size_t w = 0;
    for (size_t r = 0; r < placementStore.size(); ++r) {
        if (placementStore[r].canStillUse()) {
            if (w < r) {
                placementStore[w] = placementStore[r];
            }
            ++w;
        }
    }
    placementStore.resize(w);
    if (w == 0) {
        bestPlacement.forget();
    }
    else {
        bestPlacement = placementStore[0];
    }
}
