//
//  parsimonymove.cpp
//  Created by James Barbetti on 18-Jan-2021.
//

#include "parsimonymove.h"


ParsimonyMove::ParsimonyMove(): lazy(false), benefit(-1.0)
                              , source_branch_id(-1)
                              , positions_considered(0) {
}

bool ParsimonyMove::operator < (const ParsimonyMove& rhs) const {
    return benefit < rhs.benefit;
}

bool ParsimonyMove::operator <= (const ParsimonyMove& rhs) const {
    return benefit <= rhs.benefit;
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
