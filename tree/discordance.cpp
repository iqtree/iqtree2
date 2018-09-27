//
//  discordance.cpp
//  tree
//
//  Created by Minh Bui on 24/9/18.
//

#include "phylotree.h"

void PhyloTree::computeSiteConcordanceFactor(map<int,BranchSupportInfo> &branch_supports) {
    Branches branches;
    getInnerBranches(branches);
    for (auto it = branches.begin(); it != branches.end(); it++) {
        double num_sites;
        double sup = computeSiteConcordanceFactor(it->second, num_sites);
        string sup_str = convertDoubleToString(round(sup*1000)/10);
        Node *node = it->second.second;

        auto brsup = branch_supports.find(node->id);
        ASSERT(brsup != branch_supports.end());
        brsup->second.siteCF = sup;
        brsup->second.siteN = num_sites;
        brsup->second.length = node->findNeighbor(it->second.first)->length;

        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += "[&sCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",!sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }
}

double PhyloTree::computeSiteConcordanceFactor(Branch &branch, double &num_sites) {
    vector<IntVector> taxa;
    taxa.resize(4);

    // extract the taxa from the two left subtrees
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        getTaxaID(taxa[id], (*it)->node, branch.first);
        id++;
        if (id > 2)
            outError(__func__, " only work with bifurcating tree");
    }

    // extract the taxa from the two right subtrees
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        getTaxaID(taxa[id], (*it)->node, branch.second);
        id++;
        if (id > 4)
            outError(__func__, " only work with bifurcating tree");
    }
    
    double sum_support = 0.0;
    int sum_sites = 0;
    for (int i = 0; i < Params::getInstance().site_concordance; i++) {
        int j;
        // get a random quartet
        IntVector quartet;
        quartet.resize(taxa.size());
        for (j = 0; j < taxa.size(); j++) {
            quartet[j] = taxa[j][random_int(taxa[j].size())];
        }

        int support[3] = {0, 0, 0};
        for (auto pat = aln->begin(); pat != aln->end(); pat++) {
            if (!pat->isInformative()) continue;
            bool informative = true;
            for (j = 0; j < quartet.size(); j++)
                if (pat->at(quartet[j]) >= aln->num_states) {
                    informative = false;
                    break;
                }
            if (!informative) continue;
            if (pat->at(quartet[0]) == pat->at(quartet[1]) && pat->at(quartet[2]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[2]))
                support[0] += pat->frequency;
            if (pat->at(quartet[0]) == pat->at(quartet[2]) && pat->at(quartet[1]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[1]))
                support[1] += pat->frequency;
            if (pat->at(quartet[0]) == pat->at(quartet[3]) && pat->at(quartet[1]) == pat->at(quartet[2]) && pat->at(quartet[0]) != pat->at(quartet[1]))
                support[2] += pat->frequency;
        }
        int sum = support[0] + support[1] + support[2];
        sum_sites += sum;
        if (sum > 0)
            sum_support += ((double)support[0]) / sum;
    }
    num_sites = (double)sum_sites / Params::getInstance().site_concordance;
    return sum_support / Params::getInstance().site_concordance;
}
