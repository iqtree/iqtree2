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

#include "phylotree.h"
#include "parsimonymove.h"
#include "parsimonysearch.h"
#include <placement/placementcostcalculator.h>     //for ParsimonyCostCalculator
#include <utils/timekeeper.h>                      //for TimeKeeper

namespace {
struct ParsimonyLazyTBRMove : public ParsimonyMove {
public:
    typedef  ParsimonyLazyTBRMove this_type;
    typedef  ParsimonyMove        super;

    double   disconnection_benefit; //benefit of snipping the source branch
    intptr_t first_target_branch_id;
    intptr_t second_target_branch_id;
    intptr_t better_positions;
    
    PhyloBranch copy_of_source;
    PhyloBranch copy_of_first_target;
    PhyloBranch copy_of_second_target;
    
    int depth; //search depth is used for tie breaks
               //when two TBR moves have equal benefit. Short range
               //tbr moves mess up less of the tree, and are to
               //be preferred over long-range moves, for that reason.

    ParsimonyLazyTBRMove(const this_type& rhs)            = default;
    ParsimonyLazyTBRMove& operator=(const this_type& rhs) = default;
    ParsimonyLazyTBRMove(): super() {
        initialize(0, true);
    }
    
    static intptr_t getParsimonyVectorSize(intptr_t radius) {
        return 0;
    }
    
    virtual void initialize(intptr_t id_of_source_branch, bool beLazy) {
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
    
    virtual std::string getDescription() const {
        std::stringstream s;
        s << " linking branch " << first_target_branch_id
          << " to branch " << second_target_branch_id
          << " at the expense of branch " << source_branch_id;
        return s.str();
    }
    
    virtual void finalize(PhyloTree& tree,
                  const TargetBranchRange& branches) {
        if (first_target_branch_id < 0 || second_target_branch_id < 0) {
            return;
        }
        if (0<benefit) {
            TREE_LOG_LINE(tree, VB_DEBUG, "move s=" << source_branch_id
                << ", t1=" << first_target_branch_id
                << ", t1=" << second_target_branch_id
                << ", b=" << benefit);
            copy_of_first_target  = branches[first_target_branch_id];
            copy_of_source        = branches[source_branch_id];
            copy_of_second_target = branches[second_target_branch_id];
        }
    }
    
    bool doBranchesTouch(const TargetBranchRange& branches,
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
    
    virtual bool isStillPossible(const TargetBranchRange& branches,
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
    virtual double recalculateBenefit(PhyloTree& tree, TargetBranchRange& branches,
                              LikelihoodBlockPairs &blocks) const {
        
        auto t1 = branches[first_target_branch_id];
        auto t2 = branches[second_target_branch_id];
        int reconnect_cost = 0;
        tree.computeParsimonyOutOfTree(t1.getParsimonyBlock(),
                                       t2.getParsimonyBlock(),
                                       &reconnect_cost);
        return branches[source_branch_id].getBranchCost() - reconnect_cost;
    }
    
    class LazyTBRSearch {
    public:
        const PhyloTree&         tree;
        const TargetBranchRange& branches;
        int                      max_radius;
        std::vector<UINT*>&      path_parsimony;
        double                   parsimony_score; //of the tree as it is *before* any TBR move
        ParsimonyLazyTBRMove&    move;
        PhyloNode*               front;
        PhyloNode*               back;
        intptr_t                 first_id;
        int                      first_depth;
        
        LazyTBRSearch(const PhyloTree& t, const TargetBranchRange& b, int r,
                  std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m )
            : tree(t), branches(b), max_radius(r)
            , path_parsimony(p), parsimony_score(s), move(m) {
            auto source = branches[m.source_branch_id];
            front       = source.first;
            back        = source.second;
        }
         /**
          * @param from    where we are
          * @param prev    where we were
          * @param depth  how far we have gone (0 if we are on a branch adjacent
          *              to the source branch).
          * @note  this should not be declared virtual, as it has the same name
          *        as a member function of one of its subclasses and we want
          *        LazyTBRSearch::findMove to see *this* function, not the version
          *        in the subclass.
          */
        void searchPart1(PhyloNode* from, PhyloNode* prev, int depth) {
            if (depth<max_radius-3) {
                FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
                    searchPart1(node, from, depth+1);
                }
            }
            if (1<depth) {
                first_id    = prev->findNeighbor(from)->id;
                first_depth = depth;
                FOR_EACH_ADJACENT_PHYLO_NODE(back, front , it, node) {
                    searchPart2(node, back, depth+1);
                }
            }
        }
        /**
         * @param from    where we are
         * @param prev    where we were
         * @param depth  how far we have searched on both sides of the
         *              source branch (==first_depth+1, if we are looking at
         *              a branch adjacent to the "back" end of the source branch).
         * @note  this should not be declared virtual, as it has the same name
         *        as a member function of one of its subclasses and we want
         *        LazyTBRSearch::findMove to see *this* function, not the version
         *        in the subclass.
         */
        void searchPart2(PhyloNode* from, PhyloNode* prev,
                         int depth) {
            if (depth<max_radius-1) {
                FOR_EACH_ADJACENT_PHYLO_NODE(from, prev, it, node) {
                    searchPart2(node, from, depth+1);
                }
            }
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
        }
        void considerMove(intptr_t first_id, intptr_t second_id,
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
    };
    
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius, double disconnection_benefit,
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
    
    void getOtherNeighbors(PhyloNode* of, PhyloNode* but_not,
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
    
    void getBranchNodes(const TargetBranch& b, PhyloNode** put_here) {
        put_here[0] = b.first;
        put_here[1] = b.second;
    }
    
    void disconnect(PhyloNode* first, PhyloNode* second, PhyloNode* third) {
        first->updateNeighbor  ( third, second );
        second->updateNeighbor ( third, first  );
    }
    
    void reconnect(PhyloNode* first, PhyloNode* second,
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
    
    void updateBranch(TargetBranchRange& branches, intptr_t id,
                      PhyloNode* left, PhyloNode* right) {
        TargetBranch& branch = branches[id];
        branch.updateMapping(id, left, right, false);
    }
    
    virtual double apply(PhyloTree& tree,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks) {
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

        double score;
        for (int i=0; i<7; ++i) {
            auto id = branch_ids[i];
            TargetBranch& branch = branches[id];
            score = branch.computeState(tree, id, blocks);
            branch.setParsimonyLength(tree);
        }
        TREE_LOG_LINE(tree, VB_MAX, "Updated parsimony score"
                      << " after applying TBR move was " << score);
        return score;
    }
}; //ParsimonyLazyTBRMove


struct ProperParsimonyTBRMove : public ParsimonyLazyTBRMove {
public:
    typedef ParsimonyLazyTBRMove super;
    static intptr_t getParsimonyVectorSize(intptr_t radius) {
        return radius+1;
    }
    class ProperTBRSearch: public LazyTBRSearch {
    public:
        typedef LazyTBRSearch super;
        std::vector<PhyloNode*> path;
        UINT*  front_parsimony;
        double front_score;
        UINT*  back_parsimony;
        double back_score;
        ProperTBRSearch(const PhyloTree& t, const TargetBranchRange& b, int r,
                  std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m )
            : super(t, b, r, p, s, m), front_parsimony(nullptr)
            , back_parsimony(nullptr) {
            path.resize(max_radius, nullptr);
        }
        inline UINT* offPathParsimony(PhyloNode* a, PhyloNode* b, PhyloNode* c ) {
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
        
        /**
         * @param from    where we are
         * @param prev    where we were
         * @param depth  how far we have gone (0 if we are on a branch adjacent
         *              to the source branch).
         * @note  path_parsimony[max_radius-1] is the calculated vector
         *        for the view of the subtree from the first target branch.
         *        path[0..depth-1] are the nodes visited in the path to "prev".
         */
        void searchPart1(PhyloNode* from, PhyloNode* prev, int depth) {
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
        /**
         * @param from    where we are (on the back side of the source branch)
         * @param prev    where we were
         * @param depth  how far we have searched on both sides of the
         *              source branch (==first_depth+1, if we are looking at
         *              a branch adjacent to the "back" end of the source branch).
         * @note  path_parsimony[max_radius] is the calculated vector
         *        for the view of the subtree from the second target branch.
         *        path[0..first_depth]] are the nodes visited in the path from
         *        front to the first target branch.
         *        path[first_depth+1..depth-1] are the nodes visited in the path
         *        from back to where we are search now.
         */
        void searchPart2(PhyloNode* from, PhyloNode* prev,
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
                    double move_score  = tree.computeParsimonyOutOfTree
                                         ( path_parsimony[max_radius-1],
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
    };
    
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius, double disconnection_benefit,
                          std::vector<UINT*>& path_parsimony,
                          double parsimony_score) {
        if (lazy) {
            super::findMove(tree, branches, radius,
                            disconnection_benefit,
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
}; //ProperParsimonyTBRMove
}; //namespace

void PhyloTree::doParsimonyTBR() {
    if (leafNum<6) {
        return;
    }
    ParsimonySearchParameters s;
        
    s.name                      = "TBR";
    s.iterations                = params->parsimony_tbr_iterations;
    s.lazy_mode                 = params->use_lazy_parsimony_tbr;
    s.radius                    = params->tbr_radius; 

    if (s.lazy_mode) {
        doParsimonySearch<ParsimonyLazyTBRMove>(s);
    } else {
        doParsimonySearch<ProperParsimonyTBRMove>(s);
    }
}
