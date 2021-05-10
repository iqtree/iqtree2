//
//  parsimonymove.cpp
//  Class for:
//  (a) finding a tree rearrangement (lazy NNI, SPR, or TBR)
//      that has a fixed branch (e.g. for NNI, the branch 
//      "around" which the interchange happens, for SPR, 
//      the branch at the "top" of the subtree that is to 
//      be pruned and regrafted, for TBR, the branch that
//      will be removed), without modifying the tree,
//  (b) checking whether a tree rearrangement is still valid
//  (c) making (or reversing) a tree rearrangement.
//
//  Created by James Barbetti on 18-Jan-2021.
//

#include "parsimonymove.h"

ParsimonyPathVector::ParsimonyPathVector
    ( intptr_t blocks, intptr_t min_threads, intptr_t threads_to_use )
    : pv_per_thread(blocks), min_thread_count(min_threads), thread_count(threads_to_use) {
    ASSERT(0<=pv_per_thread);
    ASSERT(1<=min_threads);
    ASSERT(1<=threads_to_use);
}

intptr_t ParsimonyPathVector::getBlocksPerThread() const {
    return pv_per_thread;
}

intptr_t ParsimonyPathVector::getNumberOfPathsRequired() const {
    return (thread_count < min_thread_count) ? min_thread_count : thread_count;
}

intptr_t ParsimonyPathVector::getTotalNumberOfBlocksRequired() const {
    return getNumberOfPathsRequired() * pv_per_thread;
}

ParsimonyMove::ParsimonyMove(): lazy(false), benefit(-1.0)
                              , source_branch_id(-1)
                              , positions_considered(0) {
}

intptr_t ParsimonyMove::getParsimonyVectorSize(intptr_t radius) {
    ASSERT(0 && "getParsimonyVectorSize implementation not provided");
    return 0;
}

intptr_t ParsimonyMove::getMinimumPathVectorCount() {
    return 1;
}

bool ParsimonyMove::operator < (const ParsimonyMove& rhs) const {
    if (benefit < rhs.benefit) return true;
    if (rhs.benefit < benefit) return false;
    return source_branch_id < rhs.source_branch_id;
}

bool ParsimonyMove::operator <= (const ParsimonyMove& rhs) const {
    if (benefit < rhs.benefit) return true;
    if (rhs.benefit < benefit) return false;
    return source_branch_id <= rhs.source_branch_id;
}

double ParsimonyMove::getBenefit() const {
    return benefit;
}

bool ParsimonyMove::isNoLongerPossible(const TargetBranchRange& branches,
                                       PhyloBranchVector& path) const {
    return !isStillPossible(branches, path);
}

bool ParsimonyMove::isAConnectedThroughBToC
     ( PhyloNode* a, PhyloNode* b,
       PhyloNode* c, PhyloBranchVector& path) const {
    //Note: This uses a breadth-first search, because the
    //SPR radius is likely to be low (and a depth-first
    //search would have an expected running time proportional
    //to the size of the tree) (though it would be prettier!).
         
    //
    //Note 1: it would be faster, to search from *both* ends
    //        with a second search radiating outward from c
    //        if (c connects to a before b connects to c, return
    //        false, no such path). Though, that would
    //        be more complicated and you'd need a BoolVector
    //        to keep track of which nodes had been reached, to know
    //        when a successful path.
    //        (first half:  from a, through b to x)
    //        (second half: from c to x)
    //Note 2: there doesn't seem to be much point making it
    //        faster, because checking nodes are still connected
    //        is ony a small part of the running time for parsimony
    //        searches.
    //
         
    std::vector<PhyloBranch> this_layer;
    std::vector<PhyloBranch> previous_layers;
    this_layer.push_back(PhyloBranch(a,b));
    do
    {
        std::vector<PhyloBranch> next_layer;
        for (auto it=this_layer.begin(); it!=this_layer.end(); ++it) {
            PhyloBranch branch(*it);
            FOR_EACH_ADJACENT_PHYLO_NODE(branch.second, branch.first, itNode, x) {
                if (x==c) {
                    //Reconstruct the pathway, all the way back to a,
                    //from what is in previous_layers.
                    path.emplace_back(branch.second, x);
                    path.emplace_back(branch);
                    PhyloNode* hunting = branch.first;
                    for (auto rev=previous_layers.rbegin();
                         rev!=previous_layers.rend(); ++rev) {
                        PhyloBranch way_back(*rev);
                        if (way_back.second==hunting) {
                            path.emplace_back(way_back);
                            hunting = way_back.first;
                        }
                    }
                    path.reverseAll();
                    return true;
                }
                else {
                    next_layer.emplace_back(branch.second, x);
                }
            }
            previous_layers.emplace_back(branch);
        }
        std::swap(this_layer, next_layer);
    } while ( ! this_layer.empty());
    return false;
}
