//
//  parsimonytbr.cpp
//
//  Created by James Barbetti on 8/12/20.
//

#include "phylotree.h"

void PhyloTree::doParsimonyTBR() {
    size_t max_iterations = params->parsimony_tbr_iterations;
    if (max_iterations<1) {
        return;
    }
    PhyloBranchVector branches;
    getBranches(branches);
    size_t branch_count = branches.size();
    if (branch_count<4) {
        return;
    }
}

