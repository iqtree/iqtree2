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
               , partial_pars(nullptr)
               , connection_cost(0), branch_cost(0), parsimony_dirtiness(1)
               , partial_lh(nullptr)
               , scale_num(nullptr), branch_lh_scale_factor(0)
               , used(false)
               , replacements(nullptr) {}

TargetBranch::~TargetBranch() {
    delete replacements;
}

void TargetBranch::copyComputedState(const TargetBranch& rhs) {
    partial_pars           = rhs.partial_pars;
    connection_cost        = rhs.connection_cost;
    branch_cost            = rhs.branch_cost;
    parsimony_dirtiness    = rhs.parsimony_dirtiness;
    partial_lh             = rhs.partial_lh;
    scale_num              = rhs.scale_num;
    branch_lh_scale_factor = rhs.branch_lh_scale_factor;
}

TargetBranch::TargetBranch(const TargetBranch& rhs): super(rhs) {
    blocker                = rhs.blocker;
    used                   = rhs.used;
    replacements           = (rhs.replacements==nullptr) ? nullptr
                           : new ReplacementBranchList(*rhs.replacements);
    copyComputedState(rhs);
}

TargetBranch& TargetBranch::operator= (const TargetBranch& rhs) {
    if (this!=&rhs) {
        super::operator=(rhs);
        blocker                = rhs.blocker;
        used                   = rhs.used;
        delete replacements;
        replacements           = (rhs.replacements==nullptr) ? nullptr
                               : new ReplacementBranchList(*rhs.replacements);
        copyComputedState(rhs);
    }
    return *this;
}

TargetBranch::TargetBranch(BlockAllocator* allocator,
               PhyloNode* node1, PhyloNode* node2,
               bool parsimony_wanted, bool likelihood_wanted)
    : super(node1, node2), blocker(allocator), partial_pars(nullptr)
    , connection_cost(0), branch_cost(0), parsimony_dirtiness(1)
    , partial_lh(nullptr), scale_num(nullptr)
    , branch_lh_scale_factor(0), used(false)
    , replacements(nullptr) {
    if (parsimony_wanted) {
        blocker->allocateParsimonyBlock(partial_pars);
    }
    if (likelihood_wanted) {
        blocker->allocateLikelihoodBlocks(partial_lh, scale_num);
    }
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
        s << name << " does not have a partial likelihood block.\n";
    }
    TREE_LOG_LINE(phylo_tree, level, s.str() );
}

double TargetBranch::computeState(PhyloTree& phylo_tree,
                                  double& tree_parsimony_score,
                                  intptr_t target_branch_index,
                                  LikelihoodBlockPairs &blocks) {
    PhyloNeighbor* neigh1   = getLeftNeighbor();
    PhyloNeighbor* neigh2   = getRightNeighbor();
    ParallelParsimonyCalculator c(phylo_tree, false);
    parsimony_dirtiness    += c.schedulePartialParsimony(neigh1, first);
    parsimony_dirtiness    += c.schedulePartialParsimony(neigh2, second);
    c.calculate();
    int second_subtree_cost = phylo_tree.getSubTreeParsimony(neigh1);
    int first_subtree_cost  = phylo_tree.getSubTreeParsimony(neigh2);
    if ( partial_pars != nullptr ) {
        if (0<parsimony_dirtiness) {
            connection_cost = phylo_tree.computePartialParsimonyOutOfTree
                              (neigh1->partial_pars, neigh2->partial_pars,
                               partial_pars);
        }
    } else {
        if (tree_parsimony_score==-1) {
            int branch_cost_dummy = 0;
            tree_parsimony_score  = phylo_tree.computeParsimonyOutOfTree
                                    ( neigh1->partial_pars, neigh2->partial_pars,
                                      &branch_cost_dummy);
        }
        connection_cost = tree_parsimony_score;
    }
    parsimony_dirtiness = 0;
    branch_cost         = connection_cost - first_subtree_cost - second_subtree_cost;
    if (phylo_tree.leafNum < 50) {
        TREE_LOG_LINE(phylo_tree, VB_MAX, "TB Parsimony (net) was "
                      << connection_cost << " (" << branch_cost << ")"
                      << ", neigh1->length was " << neigh1->length);
    }
    if (partial_lh != nullptr) {
        blocker->makeTreeReady(first, second);
        LikelihoodBufferSet localBuffers(phylo_tree.tree_buffers);
        double score = -phylo_tree.computeLikelihoodBranch(neigh1, first, localBuffers);
        if (phylo_tree.leafNum < 50) {
            TREE_LOG_LINE(phylo_tree, VB_MAX, "TB Initial Likelihood was " << score);
        }

        double old_length      = neigh1->length;
        double half_old_length = 0.5 * old_length;
        PlacementTraversalInfo info(phylo_tree, localBuffers, nullptr, nullptr);
        
        PhyloNode fakeInterior;
        fakeInterior.addNeighbor(second,        half_old_length);
        fakeInterior.addNeighbor(first,         half_old_length);

        PhyloNode fakeExterior;
        fakeInterior.addNeighbor(&fakeExterior, half_old_length);
        fakeExterior.addNeighbor(&fakeInterior, half_old_length);

        blocks.resize(2);

        PhyloNeighbor* fakeNeigh1   = fakeInterior.findNeighbor(second);
        if (!second->isLeaf()) {
            blocks[0].allocate(phylo_tree);
            blocks[0].lendTo(fakeNeigh1);
            info.computePartialLikelihood(fakeNeigh1, first);
            //dumpNeighbor(VB_MIN, "LO", phylo_tree, neigh1);
            //dumpNeighbor(VB_MIN, "LF", phylo_tree, fakeNeigh1);
        }

        PhyloNeighbor* fakeNeigh2   = fakeInterior.findNeighbor(first);
        if (!first->isLeaf()) {
            blocks[1].allocate(phylo_tree);
            blocks[1].lendTo(fakeNeigh2);
            info.computePartialLikelihood(fakeNeigh2, second);
            //dumpNeighbor(VB_MIN, "RO", phylo_tree, neigh2);
            //dumpNeighbor(VB_MIN, "RF", phylo_tree, fakeNeigh2);
        }

        PhyloNeighbor* fakeNeigh3   = fakeExterior.findNeighbor(&fakeInterior);
        fakeNeigh3->partial_lh      = partial_lh;
        fakeNeigh3->scale_num       = scale_num;
        fakeNeigh3->lh_scale_factor = branch_lh_scale_factor;
        fakeNeigh3->setLikelihoodComputed(false);
        info.computePartialLikelihood(fakeNeigh3, &fakeExterior);
        branch_lh_scale_factor      = fakeNeigh3->lh_scale_factor;
        //dumpNeighbor(VB_MIN, "MF", phylo_tree, fakeNeigh3);
    }
    return connection_cost;
}

void TargetBranch::updateMapping(intptr_t branch_id,
                                 PhyloNode* updated_first,
                                 PhyloNode* updated_second,
                                 bool forceClearOfReversePartialParsimony) {
    first  = updated_first;
    second = updated_second;
    PhyloNeighbor* first_nei  = first->findNeighbor(second);
    PhyloNeighbor* second_nei = second->findNeighbor(first);
    first_nei ->setParsimonyComputed(false);
    second_nei->setParsimonyComputed(false);
    first_nei ->id = static_cast<int>(branch_id);
    second_nei->id = static_cast<int>(branch_id);
    if (forceClearOfReversePartialParsimony) {
        first->clearReversePartialParsimony(second);
        second->clearReversePartialParsimony(first);
    }
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

UINT* TargetBranch::getParsimonyBlock() const {
    return partial_pars;
}

double* TargetBranch::getLikelihoodBlock() const {
    return partial_lh;
}

UBYTE* TargetBranch::getScaleNumBlock() const {
    return scale_num;
}

double TargetBranch::getLhScaleFactor() const {
    return branch_lh_scale_factor;
}

void TargetBranch::setLhScaleFactor(double v) {
    branch_lh_scale_factor = v;
}

void TargetBranch::takeOwnershipOfReplacementVector(ReplacementBranchList* branches) {
    replacements = branches;
}

ReplacementBranchList* TargetBranch::getReplacements() const {
    return replacements;
}

TargetBranchRange::TargetBranchRange(PhyloTree& phylo_tree, BlockAllocator* b,
                                     PlacementCostCalculator* calculator,
                                     bool match_branch_numbers): super() {
    PhyloNodeVector v1, v2;
    if (match_branch_numbers) {
        phylo_tree.getBranchesInIDOrder(v1, v2);
    }
    else {
        phylo_tree.getBranches(v1, v2);
    }
    reserve(v1.size());
    if ( verbose_mode >= VB_DEBUG ) {
        std::stringstream s1;
        s1 << "TargetBranchRange will have " << v1.size() << " entries";
        phylo_tree.logLine(s1.str());
    }
    for (int i=0; i<v1.size(); ++i) {
        emplace_back(b, v1[i], v2[i],
                     calculator->usesParsimony(),
                     calculator->usesLikelihood());
    }
}

TargetBranchRange::TargetBranchRange(const TargetBranchRange& tbr,
                                     const std::vector<size_t>& indicesOfSubset) {
    for (size_t i : indicesOfSubset) {
        emplace_back(tbr.at(i));
    }
    for (TargetBranch& b : *this) {
        if (b.replacements != nullptr ) {
            //it's a pointer to a vector of TargetBranchRef
            //instances, which will return to tbr, which we
            //don't want.  Besides which, the indices are into tbr
            //rather than into (this).  Ditch all that!  Or, better yet,
            //throw.  Since it shouldn't ever happen.
            ASSERT(0 && "replaced branches shouldn't be copied between ranges");
            delete b.replacements;
            b.replacements = nullptr;
        }
        b.used = false;
    }
}


void TargetBranchRange::getNodes(NodeVector& vec) const {
    std::set<PhyloNode*> seen;
    for (int r = 0; r < size(); ++r ) {
        PhyloNode* first = at(r).first;
        if (seen.insert(first).second) {
            vec.push_back(first);
        }
        PhyloNode* second = at(r).second;
        if (seen.insert(second).second) {
            vec.push_back(second);
        }
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

void TargetBranchRange::reload(const PhyloTree& phylo_tree) {
    PhyloNodeVector first_node;
    PhyloNodeVector second_node;
    phylo_tree.getBranchesInIDOrder(first_node, second_node);
    size_t branch_count = first_node.size();
    resize(branch_count);
    for (int b=0; b<branch_count; ++b) {
        TargetBranch& tb(at(b));
        tb.first  = first_node[b];
        tb.second = second_node[b];
        tb.connection_cost = 0;
        tb.branch_cost = 0;
        tb.used = false;
        if (tb.replacements!=nullptr) {
            delete tb.replacements;
        }
    }
}

void TargetBranchRange::getFinalReplacementBranchIndexes(intptr_t top_index,
                                                     std::vector<size_t> &ids) const {
    ids.clear();
    const TargetBranch& top = at(top_index);
    if (top.getReplacements() == nullptr) {
        return;
    }
    ReplacementBranchList next_layer;
    for (TargetBranchRef& branch : *(top.getReplacements())) {
        next_layer.push_back(branch);
    }
    while (!next_layer.empty()) {
        ReplacementBranchList current_layer;
        std::swap(current_layer, next_layer);
        for (TargetBranchRef& branch : current_layer) {
            TargetBranch* b = branch.getTarget();
            ReplacementBranchList* reps = b->getReplacements();
            if (reps==nullptr) {
                ids.push_back(branch.getTargetIndex());
            } else {
                ASSERT(!reps->empty());
                for (TargetBranchRef& rep_branch : *reps) {
                    next_layer.push_back(rep_branch);
                }
            }
        }
    }
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

const TargetBranch* TargetBranchRef::getTarget() const {
    if (target_range == nullptr) {
        return nullptr;
    }
    return &target_range->at(target_index);
}

TargetBranch* TargetBranchRef::getTarget() {
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
                                             double& parsimony_score,
                                             bool likelihood_wanted) {
    size_t index /*of added target branch*/ = size();
    emplace_back(&allocator, node1, node2, true,
                 likelihood_wanted);
    back().computeState(allocator.getTree(), parsimony_score, index, blocks);
    back().setParsimonyLength(allocator.getTree());
    return TargetBranchRef(this, index);
}

void TargetBranch::costPlacementOfTaxa
    (PhyloTree&         phylo_tree,
     TargetBranchRange& targets,
     size_t             targetNumber,
     TaxaToPlace&       candidates,
     intptr_t           candidateStartIndex,
     intptr_t           candidateStopIndex,
     SearchHeuristic*   heuristic,
     PlacementCostCalculator* calculator,
     bool               isFirstTargetBranch) const {

    TargetBranchRef here(&targets, targetNumber);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t i = candidateStartIndex; i< candidateStopIndex; ++i) {
        TaxonToPlace& candidate = candidates.getTaxonByIndex(i);
        if (heuristic->isPlacementWorthTrying(candidate, i, here)) {
            PossiblePlacement p;
            p.setTargetBranch(here);
            calculator->assessPlacementCost(phylo_tree, candidate, p);
            candidate.considerAdditionalPlacement(p);
        }
    }
}

bool TargetBranch::isExternalBranch() const {
    return first->isLeaf() || second->isLeaf();
}

double TargetBranch::getBranchCost() const {
    return branch_cost;
}

double TargetBranch::getConnectionCost() const {
    return connection_cost;
}

double TargetBranch::getFullDisconnectionBenefit(const PhyloTree& phylo_tree,
                                                 const TargetBranchRange& branches) const {
    ASSERT(!isExternalBranch());
    ASSERT(first->degree()  == 3);
    ASSERT(second->degree() == 3);
    //If this branch is CD, and C's other neighbors are A and B,
    //and D's other neighbors are E and F, then the parsimony benefit
    //is cost(AC)+cost(BC)+cost(CD)+cost(DE)+cost(DF) - cost(AB) - cost(EF)
    //
    //    From:     To:
    //
    //    A   B    A---B
    //     \ /
    //      C
    //      |
    //      D
    //     / \
    //    E   F    E---F
    //
    //Since it's so expensive to calculate, the output of this function
    //should probably be cached.  This is the bisection benefit in a TBR
    //iteration.
    //
    
    double disconnection_benefit = branch_cost;
    PhyloNeighborVec orphans;
    FOR_EACH_PHYLO_NEIGHBOR(first, second, it, x_to_c) {
        PhyloNeighbor* c_to_x = first->findNeighbor(x_to_c->getNode());
        int costCX = 0;
        phylo_tree.computeParsimonyOutOfTree(c_to_x->partial_pars,
                                             x_to_c->partial_pars,
                                             &costCX);
        disconnection_benefit += costCX;
        orphans.push_back(c_to_x);
    }
    FOR_EACH_PHYLO_NEIGHBOR(second, first, it, x_to_d) {
        PhyloNeighbor* d_to_x = second->findNeighbor(x_to_d->getNode());
        int costDX = 0;
        phylo_tree.computeParsimonyOutOfTree(d_to_x->partial_pars,
                                             x_to_d->partial_pars,
                                             &costDX);
        disconnection_benefit += costDX;
        orphans.push_back(d_to_x);
    }
    for (size_t i=0; i<orphans.size(); i+=2) {
        int costXY = 0;
        phylo_tree.computeParsimonyOutOfTree(orphans[i]->partial_pars,
                                             orphans[i+1]->partial_pars,
                                             &costXY);
        disconnection_benefit -= costXY;
    }
    return disconnection_benefit;
}

double TargetBranch::getFullConnectionCost(const PhyloTree& phylo_tree,
                                       const TargetBranch& other_branch) const {
    //Letter-H connection.
    //If this branch is AB, and the other is EF, returns
    //cost(AC)+cost(BC)+cost(CD)+cost(DE)+cost(DF) - cost(AB) - cost(EF)
    //(the cost of adding the cross-bar of the H).
    //
    //    From:     To:
    //
    //    A---B     A   B
    //               \ /
    //                C
    //                |
    //                D
    //               / \
    //    E---F     E   F
    //
    //This is the reconnection cost in a TBR iteration.
    //
    
    int    costCD              = 0;
    phylo_tree.computeParsimonyOutOfTree(partial_pars,
                                         other_branch.partial_pars,
                                         &costCD);
    return costCD;
}

BenefitPair TargetBranch::getPartialDisconnectionBenefit(const PhyloTree& phylo_tree,
                                                         const TargetBranchRange& branches) const {
    //
    //Letter-T disconnection
    //If this branch is CD, and C is connected to A and B, determine the
    //benefit (in terms of reduced parsimony score) of carrying out
    //disconnection A     B of the subtree rooted with C:    A-----B
    //               \   /
    //                 C                                        C
    //                 |                                        |
    //                 D                                        D
    //             (whatever)                               (whatever)
    //
    //or , from:                       to:              
    //             (whatever)                               (whatever)
    //                 C                                        C
    //                 |                                        |
    //                 D                                        D
    //                / \
    //               E   F                                   E-----F
    //
    //This is used, in conjunction with getTeeConnectionCost, when evaluating
    //possible SPR (subtree prune and regraft) moves.  
    //Calling this function is expensive (6 out-of-tree branch length 
    //calculations), so the result should probably be cached.
    //
    BenefitPair result;
    if (first->isInterior()) {
        ASSERT(first->degree() == 3);
        result.hasForwardBenefit = true;
        result.forwardBenefit = branch_cost;
    }
    if (second->isInterior()) {
        ASSERT(second->degree() == 3);
        result.hasBackwardBenefit = true;
        result.backwardBenefit = branch_cost;
    }
    return result;
}

double TargetBranch::getForwardConnectionCost(const PhyloTree& phylo_tree,
                                              const TargetBranch& other_branch) const {
    //Letter-T connection
    //If this branch is CD, and the other is AB, returns the cost of
    //connecting A     B the subtree C-D by linking A and B to C.
    //            \   /    (if C were *not* connected to anything).
    //              C
    //              |
    //              D
    //It is assumed that A-B is, at present, not connected to C via D
    //(it *might* be connected to D via C).
    //It is also assumed that C is an interior node!
    //
    auto   view_from_D = other_branch.partial_pars;   //what the updated view from D would be
    auto   view_to_D   = first->findNeighbor(second)->partial_pars;
           //what the existing view from C would be
    int    updatedCDCost;
    phylo_tree.computeParsimonyOutOfTree(view_from_D, view_to_D,
                                         &updatedCDCost);
    
    return (double)updatedCDCost;
}

double TargetBranch::getBackwardConnectionCost(const PhyloTree& phylo_tree,
                                               const TargetBranch& other_branch) const {
    //Letter-T connection
    //If this branch is CD, and the other is EF, returns the cost of
    //connecting    C the subtree D-C by linking E and F to D (if 
    //              |         D were *not* connected to anything).
    //              D
    //             / \
    //            E   F
    //It is assumed that E-F is, at present, not connected to D via C
    //(it *might* be connected to C via D).
    //It is also assumed that D is an interior node!
    //
    auto view_from_C = other_branch.partial_pars;
    auto view_to_C   = second->findNeighbor(first)->partial_pars;
    int  updatedCDCost;
    phylo_tree.computeParsimonyOutOfTree(view_from_C, view_to_C, &updatedCDCost);
    
    return (double)updatedCDCost;
}

void TargetBranch::setParsimonyLength(PhyloTree& tree) {
    auto nei        = first->findNeighbor(second);
    auto backnei    = second->findNeighbor(first);
    double parsimony_length = branch_cost / tree.getAlnNSite();
    if (parsimony_length == 0.0) {
        parsimony_length = tree.params->min_branch_length;
    }
    nei->length     = parsimony_length;
    backnei->length = parsimony_length;
}

bool TargetBranch::isOutOfDate() {
    return getLeftNeighbor()->isParsimonyComputed() &&
           getRightNeighbor()->isParsimonyComputed();
}

TargetBranch& TargetBranch::clearReverseParsimony() {
    first->clearReversePartialParsimony(second);
    second->clearReversePartialParsimony(first);
    return *this;
}

