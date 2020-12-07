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
               , connection_cost(0), branch_cost(0), net_connection_cost(0)
               , partial_lh(nullptr)
               , scale_num(nullptr), branch_lh_scale_factor(0)
               , used(false)
               , replacements(nullptr) {}

TargetBranch::~TargetBranch() {
    delete replacements;
}

TargetBranch::TargetBranch(const TargetBranch& rhs): super(rhs) {
    blocker                = rhs.blocker;
    partial_pars           = rhs.partial_pars;
    connection_cost        = rhs.connection_cost;
    branch_cost            = rhs.branch_cost;
    net_connection_cost    = rhs.net_connection_cost;
    partial_lh             = rhs.partial_lh;
    scale_num              = rhs.scale_num;
    branch_lh_scale_factor = rhs.branch_lh_scale_factor;
    used                   = rhs.used;
    replacements           = (rhs.replacements==nullptr) ? nullptr
                           : new ReplacementBranchList(*rhs.replacements);
}

TargetBranch& TargetBranch::operator= (const TargetBranch& rhs) {
    if (this!=&rhs) {
        super::operator=(rhs);
        blocker                = rhs.blocker;
        partial_pars           = rhs.partial_pars;
        connection_cost        = rhs.connection_cost;
        branch_cost            = rhs.branch_cost;
        net_connection_cost    = rhs.net_connection_cost;
        partial_lh             = rhs.partial_lh;
        scale_num              = rhs.scale_num;
        branch_lh_scale_factor = rhs.branch_lh_scale_factor;
        used                   = rhs.used;
        delete replacements;
        replacements           = (rhs.replacements==nullptr) ? nullptr
                               : new ReplacementBranchList(*rhs.replacements);
    }
    return *this;
}

TargetBranch::TargetBranch(BlockAllocator* allocator,
               PhyloNode* node1, PhyloNode* node2,
               bool likelihood_wanted)
    : super(node1, node2), blocker(allocator), partial_pars(nullptr)
    , connection_cost(0), branch_cost(0), net_connection_cost(0)
    , partial_lh(nullptr), scale_num(nullptr)
    , branch_lh_scale_factor(0), used(false)
    , replacements(nullptr) {
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
        s << name << " does not have a partial likelihood block.\n";
    }
    TREE_LOG_LINE(phylo_tree, level, s.str() );
}

void TargetBranch::computeState(PhyloTree& phylo_tree,
                                size_t target_branch_index,
                                LikelihoodBlockPairs &blocks) {
    PhyloNeighbor* neigh1 = first->findNeighbor(second);
    PhyloNeighbor* neigh2 = second->findNeighbor(first);
    ParallelParsimonyCalculator c(phylo_tree, false);
    c.schedulePartialParsimony(neigh1, first);
    c.schedulePartialParsimony(neigh2, second);
    c.calculate();
    connection_cost = phylo_tree.computePartialParsimonyOutOfTree
        ( neigh1->partial_pars, neigh2->partial_pars, partial_pars );
    branch_cost =  phylo_tree.computeParsimonyOutOfTree
                            ( neigh1->partial_pars, neigh2->partial_pars);
    net_connection_cost = connection_cost - branch_cost;
    if (phylo_tree.leafNum < 50) {
        TREE_LOG_LINE(phylo_tree, VB_MAX, "TB Parsimony (net) was "
                      << connection_cost << " (" << net_connection_cost << ")"
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
    size_t index /*of added target branch*/ = size();
    emplace_back(&allocator, node1, node2, likelihood_wanted);
    back().computeState(allocator.getTree(), index, blocks);
    return TargetBranchRef(this, index);
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

bool TargetBranch::isExternalBranch() const {
    return first->isLeaf() || second->isLeaf();
}

double TargetBranch::getBranchCost() const {
    return branch_cost;
}

double TargetBranch::getFullDisconnectionBenefit(const PhyloTree& phylo_tree) const {
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
    //This isn't as expensive as it sounds, because the net cost of connecting
    //C to A and B (minus the cost of disconnecting A and B), 
    //is stored in this->net_connection_cost, and the net cost of
    //connecting D to E and F is stored in other_branch.net_connection_cost.
    //So only the cost of the CD branch needs to be calculated.
    //
    //This is the reconnection cost in a TBR iteration.
    //
    
    double costACPlusBCMinusAB = net_connection_cost; //Connecting C to AB
    double costDEPlusDFMinusEF = other_branch.net_connection_cost; //D to EF
    int    costCD              = 0;
    phylo_tree.computeParsimonyOutOfTree(partial_pars,
                                         other_branch.partial_pars,
                                         &costCD);
    
    return costACPlusBCMinusAB + costDEPlusDFMinusEF + costCD;
}

CostPair TargetBranch::getPartialDisconnectionBenefit(const PhyloTree& phylo_tree) const {
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
    CostPair result;
    if (!first->isLeaf()) {
        ASSERT(first->degree() == 3);
        result.hasForwardCost = true;
        result.forwardCost    = branch_cost;
        PhyloNeighborVec orphans;
        FOR_EACH_PHYLO_NEIGHBOR(first, second, it, x_to_c) {
            PhyloNeighbor* c_to_x = first->findNeighbor(x_to_c->getNode());
            int costCX = 0;
            phylo_tree.computeParsimonyOutOfTree(c_to_x->partial_pars,
                                                 x_to_c->partial_pars,
                                                 &costCX);
            result.forwardCost += costCX;
            orphans.push_back(c_to_x);
        }
        for (size_t i=0; i<orphans.size(); i+=2) {
            int costAB = 0;
            phylo_tree.computeParsimonyOutOfTree(orphans[i]->partial_pars,
                                                 orphans[i+1]->partial_pars,
                                                 &costAB);
            result.forwardCost -= costAB;
        }
    }
    if (!second->isLeaf()) {
        ASSERT(second->degree() == 3);
        result.hasBackwardCost = true;
        result.backwardCost    = branch_cost;
        PhyloNeighborVec orphans;
        FOR_EACH_PHYLO_NEIGHBOR(second, first, it, x_to_d) {
            PhyloNeighbor* d_to_x = second->findNeighbor(x_to_d->getNode());
            int costDX = 0;
            phylo_tree.computeParsimonyOutOfTree(d_to_x->partial_pars,
                                                 x_to_d->partial_pars,
                                                 &costDX);
            result.backwardCost += costDX;
            orphans.push_back(d_to_x);
        }
        for (size_t i=0; i<orphans.size(); i+=2) {
            int costEF = 0;
            phylo_tree.computeParsimonyOutOfTree(orphans[i]->partial_pars,
                                                 orphans[i+1]->partial_pars,
                                                 &costEF);
            result.backwardCost -= costEF;
        }
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
    phylo_tree.computeParsimonyOutOfTree(view_from_D,
                                         view_to_D,
                                         &updatedCDCost);
    double cost = other_branch.net_connection_cost
                + updatedCDCost - branch_cost;
    return cost;
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
    double cost = other_branch.net_connection_cost
                + updatedCDCost - branch_cost;
    return cost;
}
