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

TaxonToPlace::TaxonToPlace(BlockAllocator* ba, int id, std::string name)
    : taxonId(id), taxonName(name), inserted(false)
    , partial_lh(nullptr), scale_num(nullptr) {
    PhyloTree& phylo_tree = ba->getTree();
    new_leaf     = phylo_tree.newNode(taxonId, taxonName.c_str());
    new_interior = phylo_tree.newNode();
    new_interior->addNeighbor(new_leaf, -1 );
    new_leaf->addNeighbor(new_interior, -1 );
    PhyloNeighbor* nei = new_interior->firstNeighbor();
    ba->allocateMemoryFor(nei); //The allocator knows if partial_lh & scale_num wanted
    phylo_tree.computePartialParsimony(nei, new_interior);
    partial_pars = nei->partial_pars;
    partial_lh = nei->partial_lh;
    scale_num  = nei->scale_num;
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
        if (heuristic->isPlacementWorthTrying(this, target)) {
            PossiblePlacement* p = placementArray + i;
            p->setTargetBranch(&range, i);
            calculator->assessPlacementCost(phylo_tree, this, p);
        }
    }
    phylo_tree.trackProgress(target_branch_count);
    
    auto bestI = considerPlacements(placements.data(), placements.size());
    
    if ( verbose_mode >= VB_MED ) {
        std::stringstream s3;
        s3 << "Best (lowest) score for " << this->taxonName << " was " << bestPlacement.score << " at place " << bestI;
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
void TaxonToPlace::insertIntoTree(PhyloTree& phylo_tree, BlockAllocator* b,
                    TargetBranchRange& dest,
                    PlacementCostCalculator& calculator) {
    //
    //Assumes, canInsert() returned true, and the tree has not
    //been modified in the meantime.
    //
    PhyloNode* node_1    = const_cast<PhyloNode*>(bestPlacement.getTarget()->first);
    PhyloNode* node_2    = const_cast<PhyloNode*>(bestPlacement.getTarget()->second);
    TargetBranch* target = bestPlacement.target_branch.getTarget();
    
    PhyloNeighbor* down  = new_interior->findNeighbor(new_leaf);
    down->length         = bestPlacement.lenToNewTaxon;

    PhyloNeighbor* up    = new_leaf->findNeighbor(new_interior);
    up->length = bestPlacement.lenToNewTaxon;
    up->clearComputedFlags();
    target->handOverComputedStateTo( up);

    new_interior->addNeighbor(node_2, bestPlacement.lenToNode2 );
    b->handOverComputedState (node_1->findNeighbor(node_2), new_interior->findNeighbor(node_2) );
    node_1->updateNeighbor   (node_2, new_interior, bestPlacement.lenToNode1);
    node_1->findNeighbor     (new_interior)->clearComputedFlags();

    new_interior->addNeighbor(node_1, bestPlacement.lenToNode1 );
    b->handOverComputedState (node_2->findNeighbor(node_1), new_interior->findNeighbor(node_1) );
    node_2->updateNeighbor   (node_1, new_interior, bestPlacement.lenToNode2);
    node_1->findNeighbor     (new_interior)->clearComputedFlags();
            
    //Todo: redo likelihood too
    inserted = true;
    bool  likelihood_needed = calculator.usesLikelihood();
    ReplacementBranchList* reps = new ReplacementBranchList;
    reps->emplace_back ( dest.addNewRef(b, new_interior, node_1  , likelihood_needed ) );
    reps->emplace_back ( dest.addNewRef(b, new_interior, node_2  , likelihood_needed ) );
    reps->emplace_back ( dest.addNewRef(b, new_interior, new_leaf, likelihood_needed ) );
    bestPlacement.target_branch.getTarget()->takeOwnershipOfReplacementVector(reps);
}
void TaxonToPlace::forgetGazumpedPlacements() {
    bestPlacement.forget();
}
bool TaxonToPlace::insertNearby(PhyloTree& phylo_tree, BlockAllocator* b,
                  TargetBranchRange& dest,
                  PlacementCostCalculator& calculator ) {
    TargetBranch* blocked_target = bestPlacement.getTarget();
    forgetGazumpedPlacements();
    std::vector<PossiblePlacement> placements;
    assessNewTargetBranches(phylo_tree, calculator,
                            blocked_target, placements);
    for (auto it = placements.begin(); it != placements.end(); ++it) {
        considerAdditionalPlacement(*it);
    }
    if (!canInsert()) {
        return false;
    }
    insertIntoTree(phylo_tree, b, dest, calculator);
    return true;
}
void TaxonToPlace::assessNewTargetBranches(PhyloTree& phylo_tree,
                                           PlacementCostCalculator& calculator,
                                           TargetBranch* tb,
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
    while (!stack.empty()) {
        reps = stack.back();
        stack.pop_back();
        for (auto it=reps->begin(); it!=reps->end(); ++it) {
            if (it->isUsedUp()) {
                stack.push_back(it->getTarget()->getReplacements());
            } else {
                PossiblePlacement p;
                p.setTargetBranch(*it);
                p.getTarget()->computeState(phylo_tree);
                calculator.assessPlacementCost(phylo_tree, this, &p);
                scores.emplace_back(p);
            }
        }
    }
}


LessFussyTaxon::LessFussyTaxon() : super() { }
LessFussyTaxon::LessFussyTaxon(const LessFussyTaxon& rhs) = default; //copy everything!
LessFussyTaxon::LessFussyTaxon(BlockAllocator* ba, int id, std::string name)
    : super(ba, id, name) {
}
size_t LessFussyTaxon::considerPlacements(const PossiblePlacement* placements, size_t count) {
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
