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

void TargetBranch::dumpNeighbor(VerboseMode level, const char* name,
                                PhyloTree& phylo_tree, PhyloNeighbor* nei) const
{
    PhyloNode* node = nei->getNode();
    std::stringstream s;
    s << name << " (scale factor " << nei->lh_scale_factor << ")";
    if (node->isLeaf()) {
        s << " points to exterior " << node->id
          << " (" << node->name << ")\n";
    } else {
        s << " points to interior\n";
    }
    if (nei->partial_lh!=nullptr) {
        s.precision(6);
        std::string sep  = name;
        sep += " ";
        std::string sep2 = "\n" + sep;
        for (int i=0; i<phylo_tree.lh_block_size; ++i) {
            s << (((i&7)==0) ? sep : " ");
            s << nei->partial_lh[i];
            sep = sep2;
        }
        s << "\n";
        s << "Scale";
        for (int i=0; i<phylo_tree.scale_block_size; ++i) {
            s << " " << nei->scale_num[i];
        }
        s << "\n";
    } else {
        s << name << " does not have a partial likelihood block";
    }
    TREE_LOG_LINE(phylo_tree, level, s.str() );
}

void TargetBranch::computeState(PhyloTree& phylo_tree,
                                LikelihoodBlockPairs &blocks) const {
    PhyloNeighbor* neigh1 = first->findNeighbor(second);
    PhyloNeighbor* neigh2 = second->findNeighbor(first);
    ParallelParsimonyCalculator c(phylo_tree);
    c.computePartialParsimony(neigh1, first);
    c.computePartialParsimony(neigh2, second);
    c.calculate();
    double pars = phylo_tree.computePartialParsimonyOutOfTree
        ( neigh1->partial_pars, neigh2->partial_pars, partial_pars );
    TREE_LOG_LINE(phylo_tree, VB_MAX, "TB Parsimony was " << pars
                  << ", neigh1->length was " << neigh1->length);
    if (partial_lh != nullptr) {
        blocker->makeTreeReady(first, second);
        LikelihoodBufferSet localBuffers(phylo_tree.tree_buffers);
        double score = -phylo_tree.computeLikelihoodBranch(neigh1, first, localBuffers);
        TREE_LOG_LINE(phylo_tree, VB_MIN,
                      "Preliminary likelihood score for place "
                      << pointer_to_hex(this)
                      << " was " << score);
        double old_length      = neigh1->length;
        double half_old_length = 0.5 * old_length;
        PlacementTraversalInfo info(phylo_tree, localBuffers, nullptr, nullptr);
        
        PhyloNode fakeInterior;
        fakeInterior.addNeighbor(second,        half_old_length);
        fakeInterior.addNeighbor(first,         half_old_length);

        PhyloNode fakeExterior;
        fakeInterior.addNeighbor(&fakeExterior, half_old_length);
        fakeExterior.addNeighbor(&fakeInterior, half_old_length);

        neigh1->length = half_old_length;
        neigh2->length = half_old_length;
        
        blocks.resize(2);

        PhyloNeighbor* fakeNeigh1 = fakeInterior.findNeighbor(second);
        fakeNeigh1->partial_lh    = neigh1->partial_lh;
        fakeNeigh1->scale_num     = neigh1->scale_num;
        fakeNeigh1->lh_scale_factor = neigh1->lh_scale_factor;
        if (neigh1->partial_lh != nullptr ) {
            blocks[0].allocate(phylo_tree);
            blocks[0].lendTo(neigh1);
            neigh1->setLikelihoodComputed(false);
            //info.computePartialLikelihood(neigh1, first);
        }

        PhyloNeighbor* fakeNeigh2 = fakeInterior.findNeighbor(first);
        fakeNeigh2->partial_lh    = neigh2->partial_lh;
        fakeNeigh2->scale_num     = neigh2->scale_num;
        fakeNeigh2->lh_scale_factor = neigh2->lh_scale_factor;
        if (neigh2->partial_lh != nullptr ) {
            blocks[1].allocate(phylo_tree);
            blocks[1].lendTo(neigh2);
            neigh2->setLikelihoodComputed(false);
            //info.computePartialLikelihood(neigh2, second);
        }
        
        score = -phylo_tree.computeLikelihoodBranch(neigh1, first, localBuffers);
        TREE_LOG_LINE(phylo_tree, VB_MIN,
                      "Second likelihood score for place "
                      << pointer_to_hex(this)
                      << " was " << score);
        
        std::swap(neigh1->partial_lh, fakeNeigh1->partial_lh );
        std::swap(neigh1->scale_num,  fakeNeigh1->scale_num  );
        std::swap(neigh2->partial_lh, fakeNeigh2->partial_lh );
        std::swap(neigh2->scale_num,  fakeNeigh2->scale_num  );
        fakeNeigh1->setLikelihoodComputed(true);
        fakeNeigh2->setLikelihoodComputed(true);

        //phylo_tree.tracing_lh = true;
        PhyloNeighbor* fakeNeigh3 = fakeExterior.findNeighbor(&fakeInterior);
        fakeNeigh3->partial_lh    = partial_lh;
        fakeNeigh3->scale_num     = scale_num;
        fakeNeigh3->setLikelihoodComputed(false);
        info.computePartialLikelihood(fakeNeigh3, &fakeExterior);

        neigh1->length = old_length;
        neigh2->length = old_length;

#if (0)
        dumpNeighbor(VB_MIN, "LO", phylo_tree, neigh1);
        dumpNeighbor(VB_MIN, "RO", phylo_tree, neigh2);
        dumpNeighbor(VB_MIN, "LF", phylo_tree, fakeNeigh1);
        dumpNeighbor(VB_MIN, "RF", phylo_tree, fakeNeigh2);
        dumpNeighbor(VB_MIN, "MF", phylo_tree, fakeNeigh3);
#endif
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
    int w = 0;
    for (int r = 0; r < size(); ++r ) {
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

TargetBranchRef TargetBranchRange::addNewRef(BlockAllocator& allocator,
                                             LikelihoodBlockPairs& blocks,
                                             PhyloNode* node1, PhyloNode* node2,
                                             bool likelihood_wanted) {
    emplace_back(&allocator, node1, node2, likelihood_wanted);
    back().computeState(allocator.getTree(), blocks);
    return TargetBranchRef(this, size()-1);
}

void TargetBranch::costPlacementOfTaxa
    (PhyloTree&         phylo_tree,
     TargetBranchRange& targets,
     size_t             targetNumber,
     TaxaToPlace&       candidates,
     size_t             candidateStartIndex,
     size_t             candidateStopIndex,
     SearchHeuristic*   heuristic,
     PlacementCostCalculator* calculator,
     bool               isFirstTargetBranch) const {

    TargetBranchRef here(&targets, targetNumber);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = candidateStartIndex; i< candidateStopIndex; ++i) {
        TaxonToPlace& candidate = candidates.getTaxonByIndex(i);
        if (heuristic->isPlacementWorthTrying(candidate, i, here)) {
            PossiblePlacement p;
            p.setTargetBranch(here);
            calculator->assessPlacementCost(phylo_tree, candidate, p);
            candidate.considerAdditionalPlacement(p);
        }
    }
    phylo_tree.trackProgress(candidateStopIndex - candidateStartIndex);
}
