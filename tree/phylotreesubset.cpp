#include "phylotree.h"

PhyloNodeVector PhyloTree::getAllNodesInTree() const {
    auto root = getRoot();
    PhyloNodeVector answer;
    getAllNodesInSubtree(root, nullptr, answer);
    return answer;
}

void PhyloTree::setLeafSubsetNumbersFromAlignment() {
    for (PhyloNode* leaf: getTaxaNodesInIDOrder()) {
        int subset_number = aln->getSequenceSubset(leaf->id);
        leaf->setSubsetNumber(subset_number);
    }
}

#define SUBSET_UNKNOWN (-1)

void PhyloTree::computeSubsetNumbersForInternalNodes() {
    setLeafSubsetNumbersFromAlignment();
    PhyloNodeVector all_nodes (getAllNodesInTree());
    PhyloNodeVector layer; //starts out with leaf nodes
    for (PhyloNode* visit: all_nodes) {
        //Mark all interior node subsets as uncalculated
        if (!visit->isLeaf()) {
            visit->setSubsetNumber(SUBSET_UNKNOWN);
        } else {
            layer.push_back(visit);
        }
    }
    do {
        PhyloNodeVector next_layer; //interior nodes with subset number
                                    //set, in this do-loop iteration.
        for (PhyloNode* visited : layer ) {
            int subset = visited->getSubsetNumber();
            FOR_EACH_ADJACENT_PHYLO_NODE
                (visited, nullptr, it, interior) {
                if (interior->getSubsetNumber()==SUBSET_UNKNOWN) {
                    FOR_EACH_ADJACENT_PHYLO_NODE
                        (interior, visited, it2, next_door ) {
                        if (next_door->getSubsetNumber() == subset) {
                            //A second node, adjacent to interior, is
                            //in the same subset that visited was.
                            interior->setSubsetNumber(subset);
                            next_layer.push_back(interior);
                            break;
                        }
                    }
                }
            }
        }
        std::swap(layer, next_layer);
        //Next layer of interior nodes (that just had their
        //subset numbers calculated)... will be the layer of
        //most recently visited nodes, in the next iteration 
        //of the do-loop.
    }
    while (!layer.empty());

    for (PhyloNode* visit: all_nodes) {
        int subset_number = visit->getSubsetNumber();
        ASSERT(subset_number != SUBSET_UNKNOWN);
    }
}
