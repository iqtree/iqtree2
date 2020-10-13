//
// targetbranch.cpp
// Implementation of the TargetBranch class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "targetbranch.h"
#include "blockallocator.h"
#include "parallelparsimonycalculator.h"
#include "placementcostcalculator.h"
#include "placementtraversalinfo.h"
#include "possibleplacement.h"
#include "searchheuristic.h"
#include "taxontoplace.h"

TargetBranch::TargetBranch() : super(nullptr, nullptr), blocker(nullptr)
               , partial_pars(nullptr), partial_lh(nullptr)
               , scale_num(nullptr), used(false)
               , replacements(nullptr) {}

TargetBranch::~TargetBranch() {
    delete replacements;
}

TargetBranch::TargetBranch(const TargetBranch& rhs): super(rhs) {
    blocker      = rhs.blocker;
    partial_pars = rhs.partial_pars;
    partial_lh   = rhs.partial_lh;
    scale_num    = rhs.scale_num;
    used         = rhs.used;
    replacements = (rhs.replacements==nullptr) ? nullptr
        : new ReplacementBranchList(*rhs.replacements);
}

TargetBranch& TargetBranch::operator= (const TargetBranch& rhs) {
    if (this!=&rhs) {
        super::operator=(rhs);
        blocker      = rhs.blocker;
        partial_pars = rhs.partial_pars;
        partial_lh   = rhs.partial_lh;
        scale_num    = rhs.scale_num;
        used         = rhs.used;
        delete         replacements;
        replacements = (rhs.replacements==nullptr) ? nullptr
        : new ReplacementBranchList(*rhs.replacements);
    }
    return *this;
}

TargetBranch::TargetBranch(BlockAllocator* allocator,
               PhyloNode* node1, PhyloNode* node2,
               bool likelihood_wanted)
    : super(node1, node2), blocker(allocator)
    , partial_pars(nullptr), partial_lh(nullptr)
    , scale_num(nullptr), used(false), replacements(nullptr) {
    blocker->allocateParsimonyBlock(partial_pars);
    blocker->allocateLikelihoodBlocks(partial_lh, scale_num);
    used  = false;
}

void TargetBranch::computeState(PhyloTree& phylo_tree) const {
    PhyloNeighbor* neigh1 = first->findNeighbor(second);
    PhyloNeighbor* neigh2 = second->findNeighbor(first);
    ParallelParsimonyCalculator c(phylo_tree);
    c.computePartialParsimony(neigh1, first);
    c.computePartialParsimony(neigh2, second);
    c.calculate();
    phylo_tree.computePartialParsimonyOutOfTree
        ( neigh1->partial_pars, neigh2->partial_pars, partial_pars );
    if (partial_lh != nullptr) {
        blocker->makeTreeReady(first, second);
        LikelihoodBufferSet localBuffers(phylo_tree.tree_buffers);
        double score = -phylo_tree.computeLikelihoodBranch(neigh1, first, localBuffers);
        TREE_LOG_LINE(phylo_tree, VB_DEBUG,
                      "Preliminary likelihood score for place "
                      << pointer_to_hex(this)
                      << " was " << score);
        
        double old_length            = neigh1->length;
        double half_old_length       = 0.5 * old_length;
        PlacementTraversalInfo info(phylo_tree, nullptr, nullptr);
        
        PhyloNode fakeExterior;
        PhyloNode fakeInterior;
        fakeInterior.addNeighbor(second,        half_old_length);
        fakeInterior.addNeighbor(first,         half_old_length);
        fakeInterior.addNeighbor(&fakeExterior, 0.0);
        fakeExterior.addNeighbor(&fakeInterior, 0.0);
        
        info.computePartialLikelihood(neigh1, first);
        PhyloNeighbor* fakeNeigh1 = fakeInterior.findNeighbor(second);
        fakeNeigh1->partial_lh = neigh1->partial_lh;
        fakeNeigh1->scale_num  = neigh1->scale_num;
        fakeNeigh1->setLikelihoodComputed(true);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "fakeNeigh1 " << pointer_to_hex(fakeNeigh1)
        //              << " copied partial_lh " << pointer_to_hex(neigh1->partial_lh)
        //              << " from neigh1 " << pointer_to_hex(neigh1));

        info.computePartialLikelihood(neigh2, second);
        PhyloNeighbor* fakeNeigh2 = fakeInterior.findNeighbor(first);
        fakeNeigh2->partial_lh    = neigh2->partial_lh;
        fakeNeigh2->scale_num     = neigh2->scale_num;
        fakeNeigh2->setLikelihoodComputed(true);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "fakeNeigh2 " << pointer_to_hex(fakeNeigh2)
        //              << " copied partial_lh " << pointer_to_hex(neigh2->partial_lh)
        //              << " from neigh2 " << pointer_to_hex(neigh2));

        PhyloNeighbor* fakeNeigh3 = fakeExterior.findNeighbor(&fakeInterior);
        fakeNeigh3->partial_lh    = partial_lh;
        fakeNeigh3->scale_num     = scale_num;
        fakeNeigh3->setLikelihoodComputed(false);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "neigh3(Up) " << pointer_to_hex(fakeNeigh3)
        //              << " using placement's partial_lh " << pointer_to_hex(partial_lh));
        info.computePartialLikelihood(fakeNeigh3, &fakeExterior);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "neigh3(Up) " << pointer_to_hex(fakeNeigh3)
        //              << " now has partial_lh " << fakeNeigh3->partial_lh);
                                
        //PhyloNeighbor* fakeNeigh4 = fakeInterior.findNeighbor(&fakeExterior);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "neigh4(Down) " << pointer_to_hex(fakeNeigh4));

        neigh1->length = old_length;
        neigh1->setLikelihoodComputed(false);
        neigh2->length = old_length;
        neigh2->setLikelihoodComputed(false);
    }
}

void TargetBranch::forgetState() const {
    //Todo: Ditch partial_lh and scale_num, if this instance owns them.
}

bool TargetBranch::isUsedUp() const {
    return used;
}

void TargetBranch::handOverComputedStateTo(PhyloNeighbor* nei) {
    nei->partial_pars = partial_pars;
    nei->partial_lh   = partial_lh;
    nei->scale_num    = scale_num;
    partial_pars      = nullptr;
    partial_lh        = nullptr;
    scale_num         = nullptr;
    nei->setParsimonyComputed  ( true );
    nei->setLikelihoodComputed ( nei->partial_lh != nullptr );
    used  = true;
}

const UINT* TargetBranch::getParsimonyBlock() const {
    return partial_pars;
}

const double* TargetBranch::getLikelihoodBlock() const {
    return partial_lh;
}

const UBYTE* TargetBranch::getScaleNumBlock() const {
    return scale_num;
}

void TargetBranch::takeOwnershipOfReplacementVector(ReplacementBranchList* branches) {
    replacements = branches;
}

ReplacementBranchList* TargetBranch::getReplacements() {
    return replacements;
}

TargetBranchRange::TargetBranchRange(PhyloTree& phylo_tree, BlockAllocator* b,
               PlacementCostCalculator* calculator): super() {
    PhyloNodeVector v1, v2;
    phylo_tree.getBranches(v1, v2);
    reserve(v1.size());
    if ( verbose_mode >= VB_DEBUG ) {
        std::stringstream s1;
        s1 << "TargetBranchRange will have " << v1.size() << " entries";
        phylo_tree.logLine(s1.str());
    }
    for (int i=0; i<v1.size(); ++i) {
        emplace_back(b, v1[i], v2[i], calculator->usesLikelihood());
    }
}

void TargetBranchRange::removeUsed() {
    int w=0;
    for (int r=0; r<size(); ++r) {
        if (!at(r).isUsedUp()) {
            at(w) = at(r);
            ++w;
        }
    }
    resize(w);
}

TargetBranchRef::TargetBranchRef(): target_range(nullptr), target_index(0) {}

TargetBranchRef::TargetBranchRef(const TargetBranchRef& r)
    : target_range(r.target_range), target_index(r.target_index) {}

TargetBranchRef::TargetBranchRef(TargetBranchRange* range, size_t index)
    : target_range(range), target_index(index) {}

TargetBranchRef& TargetBranchRef::operator=(const TargetBranchRef& r) = default;

bool TargetBranchRef::isUsedUp() const {
    if (target_range == nullptr) {
        return true;
    }
    return target_range->at(target_index).isUsedUp();
}

PhyloNode* TargetBranchRef::getFirst() const {
    if (target_range == nullptr) {
        return nullptr;
    }
    return target_range->at(target_index).first;
}

PhyloNode* TargetBranchRef::getSecond() const {
    if (target_range == nullptr) {
        return nullptr;
    }
    return target_range->at(target_index).second;
}

TargetBranch* TargetBranchRef::getTarget() const {
    if (target_range == nullptr) {
        return nullptr;
    }
    return &target_range->at(target_index);
}

size_t TargetBranchRef::getTargetIndex() const {
    return target_index;
}

TargetBranchRef TargetBranchRange::addNewRef(BlockAllocator* allocator,
                                             PhyloNode* node1, PhyloNode* node2,
                                             bool likelihood_wanted) {
    emplace_back(allocator, node1, node2, likelihood_wanted);
    back().computeState(allocator->getTree());
    return TargetBranchRef(this, size()-1);
}

template<>
void TargetBranch::costPlacementOfTaxa
    (PhyloTree&         phylo_tree,
     TargetBranchRange* targets,
     size_t             targetNumber,
     TaxonToPlace*      candidateStart,
     TaxonToPlace*      candidateStop,
     SearchHeuristic*   heuristic,
     PlacementCostCalculator* calculator,
     bool               isFirstTargetBranch) const {
        
    TargetBranchRef here(targets, targetNumber);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (TaxonToPlace* candidate = candidateStart;
         candidate < candidateStop; ++candidate) {
        if (heuristic->isPlacementWorthTrying(candidate, here)) {
            PossiblePlacement p;
            p.setTargetBranch(here);
            calculator->assessPlacementCost(phylo_tree, candidate, &p);
            candidate->considerAdditionalPlacement(p);
        }
    }
    phylo_tree.trackProgress(candidateStop - candidateStart);
}

template<>
void TargetBranch::costPlacementOfTaxa
    (PhyloTree&         phylo_tree,
     TargetBranchRange* targets,
     size_t             targetNumber,
     LessFussyTaxon*    candidateStart,
     LessFussyTaxon*    candidateStop,
     SearchHeuristic*   heuristic,
     PlacementCostCalculator* calculator,
     bool               isFirstTargetBranch) const {
        
    double candidateCount = candidateStop - candidateStart;
    TargetBranchRef here(targets, targetNumber);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (LessFussyTaxon* candidate = candidateStart;
         candidate < candidateStop; ++candidate) {
        if (heuristic->isPlacementWorthTrying(candidate, here)) {
            PossiblePlacement p;
            p.setTargetBranch(targets, targetNumber);
            calculator->assessPlacementCost(phylo_tree, candidate, &p);
            candidate->considerAdditionalPlacement(p);
        }
    }
    phylo_tree.trackProgress(candidateCount);
}
