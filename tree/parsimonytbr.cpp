//
//  parsimonytbr.cpp
//  (Draft) Parsimony TBR implementation
//  Created by James Barbetti on 08-Dec-2020 (as a stub).
//  (First implementation using lazy TBR, 18-Jan-2021).
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

    intptr_t first_target_branch_id;
    intptr_t second_target_branch_id;
    intptr_t better_positions;
    
    PhyloBranch copy_of_source;
    PhyloBranch copy_of_first_target;
    PhyloBranch copy_of_second_target;
    
    ParsimonyLazyTBRMove(const this_type& rhs)            = default;
    ParsimonyLazyTBRMove& operator=(const this_type& rhs) = default;
    ParsimonyLazyTBRMove(): super() {
        initialize(0, true);
    }
    
    static intptr_t getParsimonyVectorSize(intptr_t radius) {
        return 0;
    }
    
    virtual void initialize(intptr_t first_branch_id, bool beLazy) {
        lazy                    = beLazy;
        benefit                 = -1.0;
        source_branch_id        = -1;
        first_target_branch_id  = first_branch_id;
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
        if (source_branch_id < 0 || second_target_branch_id < 0) {
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
        PhyloNode* a  = copy_of_first_target.first;
        PhyloNode* b  = copy_of_source.first;
        PhyloNode* c  = copy_of_second_target.first;
        bool connected = isAConnectedThroughBToC(a,b,c, path);
        if (!connected) {
            return false;
        }
        //There's more to it! The path must go *through* the source
        //branch (the branch that will move), not just touch it.
        PhyloNode* b2 = copy_of_source.second;
        for (PhyloBranch branch : path) {
            if (branch.first==b && branch.second==b2) return true;
            if (branch.first==b2 && branch.second==b) return true;
        }
        path.clear();
        return false;
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
    
    struct TBRSnip {
    public:
        intptr_t bisection_branch_id;
        double   bisection_benefit;
        TBRSnip(): bisection_branch_id(-1), bisection_benefit(0.0) {}
        TBRSnip& operator = (const TBRSnip& rhs) = default;
        bool operator < (const TBRSnip& rhs) const {
            return bisection_benefit < rhs.bisection_benefit;
        }
    };
    typedef std::vector<TBRSnip> SnipList;
    
    struct TBRState: public TBRSnip {
    public:
        intptr_t reconnection_branch_id;
        double   reconnection_cost;
        TBRState(): reconnection_branch_id(-1), reconnection_cost(0.0) {}
        bool operator < (const TBRState& rhs) const {
            return bisection_benefit - reconnection_cost
            <  rhs.bisection_benefit - rhs.reconnection_cost;
        }
    };
    
    TBRState searchForMove(const PhyloTree& tree,
                           const TargetBranchRange& branches,
                           PhyloNode* here, PhyloNode* prev,
                           SnipList& list, int radius, int max_radius) {
        ++positions_considered;
        bool messed_with_path = false;
        if (4<radius) {
            if (list[radius-3] < list[radius-4]) {
                messed_with_path = true;
                std::swap(list[radius-4], list[radius-3]);
            }
        }
        TBRState            best;
        if (1<radius) {
            intptr_t            branch_id    = here->findNeighbor(prev)->id;
            const TargetBranch& branch       = branches[branch_id];
            list[radius].bisection_branch_id = branch_id;
            list[radius].bisection_benefit   = branch.getBranchCost();
            if (4<radius) {
                int  connection_cost         = 0;
                auto start_branch            = branches[first_target_branch_id];
                tree.computeParsimonyOutOfTree(start_branch.getParsimonyBlock(),
                                               branch.getParsimonyBlock(),
                                               &connection_cost);
                TBRSnip& snip                = list[radius-3];
                best.bisection_branch_id     = snip.bisection_branch_id;
                best.bisection_benefit       = snip.bisection_benefit;
                best.reconnection_branch_id  = branch_id;
                best.reconnection_cost       = connection_cost;
                connection_cost = 0;
                better_positions += (best.bisection_benefit > connection_cost ) ? 1 : 0;
            }
        }
        if (radius<max_radius) {
            FOR_EACH_ADJACENT_PHYLO_NODE(here, prev, it, next) {
                auto move = searchForMove(tree, branches, next, here,
                                          list, radius+1, max_radius);
                if (best < move) {
                    best = move;
                }
            }
        }
        if (messed_with_path) {
            std::swap(list[radius-4], list[radius-3]);
        }
        return best;
    }
    
    void findBestMoveWithFirstTarget(const PhyloTree& tree,
                                     const TargetBranchRange& branches,
                                     int max_radius, intptr_t first_target_id,
                                     bool lazy_mode) {
        auto t1 = branches[first_target_id];
        SnipList list;
        list.resize(max_radius+1);
        TBRState best;
        FOR_EACH_ADJACENT_PHYLO_NODE(t1.first, t1.second, it, node) {
            TBRState state = searchForMove(tree, branches, t1.first,
                                           node, list, 0, max_radius);
            if (best<state) {
                best=state;
            }
        }
        FOR_EACH_ADJACENT_PHYLO_NODE(t1.second, t1.first, it, node) {
            TBRState state = searchForMove(tree, branches, t1.second,
                                           node, list, 0, max_radius);
            if (best<state) {
                best=state;
            }
        }
        auto best_score = best.bisection_benefit - best.reconnection_cost;
        if ( 0 < best_score ) {
            second_target_branch_id = best.reconnection_branch_id;
            source_branch_id        = best.bisection_branch_id;
            benefit                 = best_score;
        }
    }
    
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius, double disconnection_benefit,
                          std::vector<UINT*> &path_parsimony,
                          double parsimony_score) {
        findBestMoveWithFirstTarget(tree, branches, radius,
                                    first_target_branch_id, lazy);
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
            TargetBranch& branch = branches[branch_ids[i]];
            score = branch.computeState(tree, id, blocks);
            branch.setParsimonyLength(tree);
        }
        TREE_LOG_LINE(tree, VB_MAX, "Updated parsimony score"
                      << " after applying TBR move was " << score);
        return score;
    }
    
}; //ParsimonyLazyTBRMove
}; //namespace


void PhyloTree::doParsimonyTBR() {
    if (leafNum<6) {
        return;
    }
    ParsimonySearchParameters s;
        
    s.name                      = "TBR";
    s.iterations                = params->parsimony_tbr_iterations;
    s.lazy_mode                 = params->use_lazy_parsimony_tbr;
    s.radius                    = params->sprDist; //no tbrDist member (yet)

    doParsimonySearch<ParsimonyLazyTBRMove>(s);
}
