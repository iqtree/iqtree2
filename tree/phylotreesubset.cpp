#include "phylotree.h"
#include <utils/gzstream.h>
#include <model/modeldivergent.h>

void PhyloTree::setUpSubtreesForDivergentModels(ModelSubst* top_model) {
    DivergentModels div_models;
    top_model->getDivergentModels(div_models);
    if (div_models.empty()) {
        //No divergent models?  Nothing to do!
        return;
    }
    std::string div_file_path = params->divergence_graph_file_path;
    if (div_file_path.empty()) {
        outError("Cannot use divergent models without a divergence graph");
    }

    //Load up divergence graph
    MTree div_tree;
    bool  is_div_tree_rooted = false;
    std::cout << "Loading divergence graph"
              << " from file: " << div_file_path 
              << "." << std::endl;
    div_tree.readTree(div_file_path.c_str(), is_div_tree_rooted);
    Node* div_root = div_tree.getRoot();

    //Set up a mapping from taxon name to taxon id.
    NameToIDMap name_to_id;
    intptr_t taxa_in_alignment = aln->getMapFromNameToID(name_to_id);
    intptr_t taxa_in_div_graph = div_tree.setNodeIdsFromMap(name_to_id);
    if (taxa_in_alignment!=taxa_in_div_graph) {
        std::stringstream complaint;
        complaint << "There were " << taxa_in_alignment << " taxa"
                  << " in the alignment, but only " 
                  << taxa_in_div_graph
                  << " of those taxa were mentioned,"
                  << " in the divergence graph.";
        outWarning(complaint.str());
    }

    //For each divergent model, identify the way that it
    //partitions the taxa in the tree
    std::vector<IntVector> subsets_for_all_models;
    ASSERT(div_models.size()==1);
    for (ModelDivergent* div_model: div_models) {
        std::vector<IntVector> subsets_for_model;        
        div_model->identifyTaxonSubsets(div_root,
                                        subsets_for_model);
        //Todo: splitSubsets(subsets_for_model, subsets_for_all_models);
        std::swap(subsets_for_all_models, subsets_for_model);
    }

    //Set up an array that maps taxon id to subset,
    //and an array that maps taxon id to taxon node
    //in the tree
    PhyloNodeVector taxa(getTaxaNodesInIDOrder());
    IntVector       taxon_to_subset(aln->getNSeq());
    int subset_number=0;
    for (IntVector& subset : subsets_for_all_models) {
        for (int taxon_id : subset) {
            taxon_to_subset[taxon_id] = subset_number;
            if (taxa.size()<=taxon_id) {
                continue;
            }
            if (taxa[taxon_id]==nullptr) {
                //If taxon is not currently in the tree, never
                //mind.  Perhaps it will be added later.
                continue;
            }
            taxa[taxon_id]->setSubsetNumber(subset_number);
        }
        ++subset_number;
    }

    computeSubsetNumbersForInternalNodes();

    for (ModelDivergent* div_model: div_models) {
        div_model->mapTaxonSubsetsToModels
            (div_root, subset_number, taxon_to_subset);
    }
}

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
