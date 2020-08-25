//
//  placement.cpp
//  Evolutionary placement: adding taxa to a phylogenetic tree
//  This file created by James Barbetti on 25/8/20, but:
//  1. addNewTaxToTree was formerly in phylotree.cpp;
//  2. addTaxonML likewise (and was the work of BUI Quang Minh)
//
#include "phylotree.h"

/****************************************************************************
 Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::addTaxonML(Node *added_node, Node* &target_node,
                             Node* &target_dad, Node *node, Node *dad) {
    Neighbor *dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    double halfLen = 0.5 * len;
    node->updateNeighbor(dad, added_node, halfLen);
    dad->updateNeighbor(node, added_node, halfLen);
    added_node->updateNeighbor((Node*) 1, node, halfLen);
    added_node->updateNeighbor((Node*) 2, dad, halfLen);
    // compute the likelihood

    double best_score = optimizeChildBranches((PhyloNode*) added_node);
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*) 1, len);
    added_node->updateNeighbor(dad, (Node*) 2, len);

    // now traverse the tree downwards
    FOR_NEIGHBOR_IT(node, dad, it){
        Node *target_node2 = nullptr;
        Node *target_dad2  = nullptr;
        double score = addTaxonML(added_node, target_node2,
                                  target_dad2, (*it)->node, node);
        if (score > best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad  = target_dad2;
        }
    }
    return best_score;
}

void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd) {
    //
    //Assumes: The tree is rooted.
    //Notes: This is a fixed-insertion-order, one-at-a-time, no-search-heuristic
    //       global likelihood implementation (suitable for adding only a few
    //       taxa at a time).
    //Todo:  Support options requested via params->incremental_method
    //       1. Batching of taxa (more than 1 at a time)
    //       2. Parsimony rather than likelihood
    //       3. The use of search heuristics
    //       4. Extra restructuring after each batch of insertions
    //
    initProgress(taxaIdsToAdd.size(),
                 "Adding new taxa to tree",
                 "added", "taxon");
    for (size_t i=0; i<taxaIdsToAdd.size(); ++i) {
        int         taxonId   = taxaIdsToAdd[i];
        std::string taxonName = aln->getSeqName(taxonId);
        
        Node* new_taxon  = newNode(taxonId, taxonName.c_str());
        Node* added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        //add dummy neighbors (need this, because of how addTaxonML works)
        added_node->addNeighbor((Node*) 1, 1.0);
        added_node->addNeighbor((Node*) 2, 1.0);
        
        //Todo: see note 3. Search heuristics
        //      (for restricting the subtree to search)
        Node* search_start    = root->neighbors[0]->node;
        Node* search_backstop = root;

        Node *target_node = nullptr;
        Node *target_dad  = nullptr;
        addTaxonML(added_node, target_node, target_dad, search_start, search_backstop);
        
        //
        // now insert the new node in the middle of the branch node-dad
        // (question: why the middle?  Why didn't addTaxonML say *where* it
        //  wanted to put the new node, and return three lengths:
        //  target_dad:added_node, added_node:target_node, added_node:new_taxon).
        // surely it calculated them!   When figuring out the score?)
        //
        double len     = target_dad->findNeighbor(target_node)->length;
        double halfLen = len * 0.5;
        target_node->updateNeighbor(target_dad, added_node, halfLen);
        target_dad->updateNeighbor(target_node, added_node, halfLen);
        
        //replace dummy neighbours (also needed, due to how addTaxonML works)
        added_node->updateNeighbor((Node*) 1, target_node, halfLen);
        added_node->updateNeighbor((Node*) 2, target_dad, halfLen);
        
        trackProgress(1.0);
    }
    doneProgress();
    clearAllPartialLH(true);
    clearAllScaleNum();
}
