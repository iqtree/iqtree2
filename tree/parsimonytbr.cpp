//
//  parsimonytbr.cpp
//  (Draft) Parsimony TBR implementation
//  Created by James Barbetti on 08-Dec-2020 (as a stub).
//  (First implementation using lazy TBR, 18-Jan-2021).
//  Note: ParsimonyLazyTBRMove's "benefit estimates"
//  are so bad that it isn't of any real use.  But a proper
//  TBR implementation could be based on it (or perhaps if
//  its estimates were better, it would be worth using).
//  (ParsimonyLazyTBRMove::apply() works, all it needs is
//  better benefit estimates).
//

#include "parsimonytbr.h"

#include "phylotree.h"
#include "parsimonysearch.h"
#include <placement/placementcostcalculator.h>     //for ParsimonyCostCalculator
#include <utils/timekeeper.h>                      //for TimeKeeper


ParsimonyLazyTBRMove::ParsimonyLazyTBRMove(): super() {
    initialize(0, true);
}
    
/*static*/ intptr_t ParsimonyLazyTBRMove::getParsimonyVectorSize(intptr_t radius) {
    return 0;
}
    
void ParsimonyLazyTBRMove::initialize(intptr_t id_of_source_branch, bool beLazy) {
    lazy                    = beLazy;
    benefit                 = -1.0;
    source_branch_id        = id_of_source_branch;
    first_target_branch_id  = -1;
    second_target_branch_id = -1;
    copy_of_source          = PhyloBranch(nullptr, nullptr);
    copy_of_first_target    = PhyloBranch(nullptr, nullptr);
    copy_of_second_target   = PhyloBranch(nullptr, nullptr);
    positions_considered    = 0;
    better_positions        = 0;
}
    
std::string ParsimonyLazyTBRMove::getDescription() const {
    std::stringstream s;
    s << " linking branch " << first_target_branch_id
      << " to branch " << second_target_branch_id
      << " at the expense of branch " << source_branch_id;
    return s.str();
}
    
void ParsimonyLazyTBRMove::finalize(PhyloTree& tree,
                                    const TargetBranchRange& branches) {
    if (first_target_branch_id < 0 || second_target_branch_id < 0) {
        return;
    }
    if (0<benefit) {
        TREE_LOG_LINE(tree, VerboseMode::VB_DEBUG, 
            "move s=" << source_branch_id
            << ", t1=" << first_target_branch_id
            << ", t1=" << second_target_branch_id
            << ", b=" << benefit);
        copy_of_first_target  = branches[first_target_branch_id];
        copy_of_source        = branches[source_branch_id];
        copy_of_second_target = branches[second_target_branch_id];
    }
}
    
bool ParsimonyLazyTBRMove::doBranchesTouch(const TargetBranchRange& branches,
                     intptr_t id_1, intptr_t id_2) const {
    const TargetBranch branch_1 = branches[id_1];
    const TargetBranch branch_2 = branches[id_2];
    FOR_EACH_ADJACENT_PHYLO_NODE(branch_1.first, nullptr, it, node) {
        if (node == branch_2.first || node == branch_2.second ) {
            return true;
        }
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(branch_1.second, nullptr, it, node) {
        if (node == branch_2.first || node == branch_2.second ) {
            return true;
        }
    }
    return false;
}
    
bool ParsimonyLazyTBRMove::isStillPossible(const TargetBranchRange& branches,
                                           PhyloBranchVector& path) const {
    path.clear();
    if (branches[source_branch_id] != copy_of_source) {
        return false;
    }
    if (branches[first_target_branch_id] != copy_of_first_target) {
        return false;
    }
    if (branches[second_target_branch_id] != copy_of_second_target) {
        return false;
    }
    if (doBranchesTouch(branches, source_branch_id, first_target_branch_id)) {
        return false;
    }
    if (doBranchesTouch(branches, source_branch_id, second_target_branch_id)) {
        return false;
    }
    if (doBranchesTouch(branches, first_target_branch_id, second_target_branch_id)) {
        return false;
    }
    PhyloNode* front       = copy_of_source.first;  //front of moving branch
                                                    //(looking forward)
    PhyloNode* back        = copy_of_source.second; //back of moving branch
    PhyloNode* x           = copy_of_first_target.first;
    PhyloBranchVector        pathLeft;  //path to first target branch
    PhyloBranchVector        pathRight; //path to second target branch
    bool x_front_connected = isAConnectedThroughBToC(back,front,x, pathLeft);
    PhyloNode* y           = copy_of_second_target.first;
    bool y_back_connected  = isAConnectedThroughBToC(front, back, y, pathRight);
    if (x_front_connected!=y_back_connected) {
        return false;
    }
    if (x_front_connected) {
        std::swap(pathLeft, path);
    } else {
        ASSERT(isAConnectedThroughBToC(front, back,  x, pathLeft));
        ASSERT(isAConnectedThroughBToC(back,  front, y, path));
    }
    path.reverseAll();
    path.pop_back();
    for (auto branch : pathLeft) {
        path.push_back(branch);
    }
    return true;
}
double ParsimonyLazyTBRMove::recalculateBenefit
               ( PhyloTree& tree, double parsimony_score,
                 TargetBranchRange& branches,
                 LikelihoodBlockPairs &blocks,
                 ParsimonyPathVector& parsimony_path_vectors) const {
    
    auto t1 = branches[first_target_branch_id];
    auto t2 = branches[second_target_branch_id];
    int reconnect_cost = 0;
    tree.computeParsimonyOutOfTree(t1.getParsimonyBlock(),
                                   t2.getParsimonyBlock(),
                                   &reconnect_cost);
    return branches[source_branch_id].getBranchCost() - reconnect_cost;
}
    
    
        
ParsimonyLazyTBRMove::LazyTBRSearch::LazyTBRSearch
    ( const PhyloTree& t, const TargetBranchRange& b, int r,
      std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m )
    : tree(t), branches(b), max_radius(r),
      path_parsimony(p), parsimony_score(s), move(m) {
    auto source = branches[m.source_branch_id];
    front       = source.first;
    back        = source.second;
}
void ParsimonyLazyTBRMove::LazyTBRSearch::searchPart1
    ( PhyloNode* from, PhyloNode* prev, int depth ) {
    if (1<depth) {
        first_id    = prev->findNeighbor(from)->id;
        first_depth = depth;
        FOR_EACH_ADJACENT_PHYLO_NODE(back, front , it, node) {
            searchPart2(node, back, depth+1);
        }
    }
    if (depth<max_radius-3) {
        FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
            searchPart1(node, from, depth+1);
        }
    }
}
void ParsimonyLazyTBRMove::LazyTBRSearch::searchPart2
    ( PhyloNode* from, PhyloNode* prev, int depth) {
    if (first_depth+2<depth) {
        intptr_t second_id  = from->findNeighbor(prev)->id;
        auto branch_one     = branches[first_id];
        auto branch_two     = branches[second_id];
        int  reconnect_cost = 0;
        tree.computeParsimonyOutOfTree(branch_one.getParsimonyBlock(),
                                       branch_two.getParsimonyBlock(),
                                       &reconnect_cost);
        double gain = move.disconnection_benefit - reconnect_cost;
        considerMove(first_id, second_id, gain, depth);
    }
    if (depth<max_radius-1) {
        FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
            searchPart2(node, from, depth+1);
        }
    }
}
void ParsimonyLazyTBRMove::LazyTBRSearch::considerMove
     ( intptr_t first_id, intptr_t second_id,
       double gain, int depth) {
    ++move.positions_considered;
    if (gain > move.benefit ||
        (gain==move.benefit && depth<move.depth)) {
        move.first_target_branch_id  = first_id;
        move.second_target_branch_id = second_id;
        move.benefit                 = gain;
        move.depth                   = depth;
    }
}

void ParsimonyLazyTBRMove::findMove(const PhyloTree& tree,
                      const TargetBranchRange& branches,
                      int radius,
                      std::vector<UINT*>& path_parsimony,
                      double parsimony_score) {
    auto source_branch    = branches[source_branch_id];
    PhyloNode* front      = source_branch.first;
    PhyloNode* back       = source_branch.second;
    
    if (front->isLeaf()) {
        return;
    }
    if (back->isLeaf()) {
        return;
    }
    disconnection_benefit = source_branch.getBranchCost();
    depth                 = radius;
    LazyTBRSearch s(tree, branches, radius,
                path_parsimony, parsimony_score, *this);
    FOR_EACH_ADJACENT_PHYLO_NODE(front, back, it, node) {
        s.searchPart1(node, front, 0);
    }
}

void ParsimonyLazyTBRMove::getOtherNeighbors(PhyloNode* of, PhyloNode* but_not,
                       PhyloNode** put_here, intptr_t* branch_ids) {
    ASSERT( of->degree() == 3);
    put_here[0] = nullptr;
    put_here[1] = nullptr;
    FOR_EACH_ADJACENT_PHYLO_NODE(of, but_not, it, node) {
        *put_here   = node;
        ++put_here;
        *branch_ids = (*it)->id;
        ++branch_ids;
    }
}

void ParsimonyLazyTBRMove::getBranchNodes(const TargetBranch& b, PhyloNode** put_here) {
    put_here[0] = b.first;
    put_here[1] = b.second;
}

void ParsimonyLazyTBRMove::disconnect(PhyloNode* first, PhyloNode* second, PhyloNode* third) {
    first->updateNeighbor  ( third, second );
    second->updateNeighbor ( third, first  );
}

void ParsimonyLazyTBRMove::reconnect(PhyloNode* first, PhyloNode* second,
                                     PhyloNode* third, PhyloNode* fourth) {
    first->updateNeighbor  ( second, third );
    second->updateNeighbor ( first,  third );
    PhyloNode* old_first  = nullptr;
    PhyloNode* old_second = nullptr;
    FOR_EACH_ADJACENT_PHYLO_NODE(third, fourth, it, node) {
        if (old_first==nullptr) {
            old_first = node;
        } else {
            old_second = node;
        }
    }
    ASSERT(old_first  != nullptr );
    ASSERT(old_second != nullptr );
    third->updateNeighbor(old_first,  first);
    third->updateNeighbor(old_second, second);
}

void ParsimonyLazyTBRMove::updateBranch
    ( TargetBranchRange& branches, intptr_t id,
      PhyloNode* left, PhyloNode* right) {
    TargetBranch& branch = branches[id];
    branch.updateMapping(id, left, right, false);
}

double ParsimonyLazyTBRMove::apply
    ( PhyloTree& tree, double parsimony_score,
      TargetBranchRange& branches, LikelihoodBlockPairs blocks,
      ParsimonyPathVector& parsimony_path_vectors) {
    //
    //Apply a TBR move (letters are indicative):
    //  A   B  G-H       A-B G   H
    //   \ /                  \ /
    //    C                    C
    //    |        -->         |
    //    D                    D
    //   / \                  / \
    //  E   F  I-J       E-F I   J
    //
    //Six branch IDs are reassigned as follows (numbers
    //are indexes into the branch_ids array, see below).
    // 0. what was GH, becomes AB
    // 1. what was IJ, becomes EF
    // 2. what was AC, becomes GC
    // 3. what was BC, becomes HC
    // 4. what was ED, becomes ID
    // 5. what was FD, becomes JD
    //
    TargetBranch& t1    = branches[first_target_branch_id].clearReverseParsimony();
    TargetBranch& t2    = branches[second_target_branch_id].clearReverseParsimony();
    TargetBranch& moved = branches[source_branch_id].clearReverseParsimony();
    
    PhyloNode* nodes[10];     //nodes that are invoved (A through J)
    intptr_t   branch_ids[7]; //branches that get messed with (element [6] is moved)
    branch_ids[0] = first_target_branch_id;
    branch_ids[1] = second_target_branch_id;
    getOtherNeighbors(moved.first, moved.second, &nodes[0], &branch_ids[2]); //A, B
    getBranchNodes   (moved, &nodes[2]);                                     //C, D
    getOtherNeighbors(moved.second, moved.first, &nodes[4], &branch_ids[4]); //E, F
    getBranchNodes   (t1,    &nodes[6]);                                     //G, H
    getBranchNodes   (t2,    &nodes[8]);                                     //I, J
    branch_ids[6] = source_branch_id;
    
    disconnect( nodes[0], nodes[1], nodes[2] );//A<->B rather than linking C
    disconnect( nodes[4], nodes[5], nodes[3] );//E<->F rather than link D
    reconnect ( nodes[6], nodes[7], nodes[2]/*C*/, nodes[3] );//Link C to G and H and vice versa
    reconnect ( nodes[8], nodes[9], nodes[3]/*D*/, nodes[2] );//Link D to I and J and vice versa
    
    updateBranch ( branches, branch_ids[0], nodes[0], nodes[1]); //t1 now AB
    updateBranch ( branches, branch_ids[1], nodes[4], nodes[5]); //t2 now EF
    updateBranch ( branches, branch_ids[2], nodes[6], nodes[2]); //AC becomes GC
    updateBranch ( branches, branch_ids[3], nodes[7], nodes[2]); //BC becomes HC
    updateBranch ( branches, branch_ids[4], nodes[8], nodes[3]); //ED becomes ID
    updateBranch ( branches, branch_ids[5], nodes[9], nodes[3]); //FD becomes JD

    parsimony_score = -1;
    for (int i=0; i<7; ++i) {
        auto id = branch_ids[i];
        TargetBranch& branch = branches[id];
        branch.computeState(tree, parsimony_score, id, blocks);
        branch.setParsimonyLength(tree);
    }
    TREE_LOG_LINE(tree, VerboseMode::VB_MAX, 
                  "Updated parsimony score"
                  << " after applying TBR move was " << parsimony_score);
    return parsimony_score;
}
/*static*/ intptr_t ProperParsimonyTBRMove::getParsimonyVectorSize(intptr_t radius) {
    return radius+1;
}
ProperParsimonyTBRMove::ProperTBRSearch::ProperTBRSearch(const PhyloTree& t, const TargetBranchRange& b, int r,
          std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m )
    : super(t, b, r, p, s, m), front_parsimony(nullptr)
    , back_parsimony(nullptr) {
    path.resize(max_radius, nullptr);
}
UINT* ProperParsimonyTBRMove::ProperTBRSearch::offPathParsimony(PhyloNode* a, PhyloNode* b, PhyloNode* c ) {
    //parsimony, viewed from a, to its non-b, non-c neighbor
    //ASSERT(a!=nullptr && b!=nullptr && c!=nullptr && b!=c && a!=b && a!=c);
    FOR_EACH_PHYLO_NEIGHBOR(a, b, it, nei) {
        if (nei->getNode() != c) {
            return nei->get_partial_pars();
        }
    }
    ASSERT(false && "could not find third adjacent node");
    return nullptr;
}

void ProperParsimonyTBRMove::ProperTBRSearch::searchPart1(PhyloNode* from, PhyloNode* prev, int depth) {
    path[depth] = prev;
    if (0==depth) {
        front_parsimony = offPathParsimony(front, from, back);
    }
    else {
        UINT* on_path     = (depth>1) ? path_parsimony[depth-2]: front_parsimony;
        UINT* off_path    = offPathParsimony ( prev, from, path[depth-1] );
        tree.computePartialParsimonyOutOfTree(off_path, on_path,
                                                  path_parsimony[depth-1]);
        if (1<depth) {
            PhyloNeighbor* nei = prev->findNeighbor(from);
            front_score        = tree.computePartialParsimonyOutOfTree
                                 ( path_parsimony[depth-1],
                                   nei->get_partial_pars(),
                                   path_parsimony[max_radius-1]);
            first_id           = nei->id;
            first_depth        = depth;
            FOR_EACH_ADJACENT_PHYLO_NODE(back, front , it, node) {
                searchPart2(node, back, depth+1);
            }
        }
    }
    if (depth<max_radius-3) {
        FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
            searchPart1(node, from, depth+1);
        }
    }
}
void ProperParsimonyTBRMove::ProperTBRSearch::searchPart2(PhyloNode* from, PhyloNode* prev,
                 int depth) {
    path[depth] = prev;
    if ( first_depth + 1 == depth ) {
        back_parsimony = offPathParsimony(back, from, front);
    }
    else {
        UINT* on_path     = (first_depth+2<depth) ? path_parsimony[depth-2]: back_parsimony;
        UINT* off_path    = offPathParsimony(prev, from, path[depth-1]);
        tree.computePartialParsimonyOutOfTree(off_path, on_path,
                                              path_parsimony[depth-1]);
        if ( first_depth + 2 < depth ) {
            PhyloNeighbor* nei = prev->findNeighbor(from);
            back_score         = tree.computePartialParsimonyOutOfTree
                                 ( path_parsimony[depth-1],
                                   nei->get_partial_pars(),
                                   path_parsimony[max_radius] );
            int reconnect_cost = 0;
            tree.computeParsimonyOutOfTree(path_parsimony[max_radius-1],
                                           path_parsimony[max_radius],
                                           &reconnect_cost );
            double gain        = parsimony_score - front_score
                               - back_score - reconnect_cost;
            
            considerMove(first_id, nei->id, gain, depth);
        }
    }
    if ( depth+1 < max_radius ) {
        FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
            searchPart2(node, from, depth+1);
        }
    }
}
    
void ProperParsimonyTBRMove::findMove(const PhyloTree& tree,
                      const TargetBranchRange& branches,
                      int radius,
                      std::vector<UINT*>& path_parsimony,
                      double parsimony_score) {
    if (lazy) {
        super::findMove(tree, branches, radius,
                        path_parsimony, parsimony_score);
        return;
    }
    auto source_branch    = branches[source_branch_id];
    PhyloNode* front      = source_branch.first;
    PhyloNode* back       = source_branch.second;
    if (front->isLeaf()) {
        return;
    }
    if (back->isLeaf()) {
        return;
    }
    disconnection_benefit = source_branch.getBranchCost();
    depth                 = radius;
    ProperTBRSearch s(tree, branches, radius,
                path_parsimony, parsimony_score, *this);
    FOR_EACH_ADJACENT_PHYLO_NODE(front, back, it, node) {
        s.searchPart1(node, front, 0);
    }
}

int PhyloTree::doParsimonyTBR(VerboseMode how_loud) {
    if (leafNum<6) {
        return computeParsimony();
    }
    ParsimonySearchParameters s("TBR");
        
    s.iterations                 = params->parsimony_tbr_iterations;
    s.lazy_mode                  = params->use_lazy_parsimony_tbr;
    s.radius                     = params->tbr_radius;
    s.calculate_connection_costs = s.lazy_mode;
    s.be_quiet                   = verbose_mode < how_loud;

    if (s.lazy_mode) {
        return doParsimonySearch<ParsimonyLazyTBRMove>(s);
    } else {
        return doParsimonySearch<ProperParsimonyTBRMove>(s);
    }
}
