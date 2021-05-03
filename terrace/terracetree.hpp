/*
 *  terracetree.hpp
 *
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#ifndef terracetree_hpp
#define terracetree_hpp

#include <stdio.h>
#include "tree/mtree.h"
#include "terracenode.hpp"

class TerraceTree: public MTree {

public:
    /**
     constructor
     */
    TerraceTree();
    
    /**
     constructor
     */
    //TerraceTree(SuperAlignment *alignment);
    
    /**
     destructor
     */
    ~TerraceTree();
    
    /**
     map taxon name to pointer to node
     */
    unordered_map<string, Node*> leafNodes;
    void fillLeafNodes();
    void fillbrNodes();
    
    unordered_map<int, NodeVector> brNodes;
    
    /**
     read the tree from the input file in newick format
     @param infile the input file file.
     @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(const char *infile, bool &is_rooted);
    
    /**
     read the tree from the ifstream in newick format
     @param in the input stream.
     @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(istream &in, bool &is_rooted);
    
    /**
     allocate a new node. Override this if you have an inherited Node class.
     @param node_id node ID
     @param node_name node name
     @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = NULL);
    
    /**
     allocate a new node. Override this if you have an inherited Node class.
     @param node_id node ID
     @param node_name node name issued by an interger
     @return a new node
     */
    virtual Node* newNode(int node_id, int node_name);
    
    /**
     copy the tree given a list of taxon names that should remain on the tree (not yet a 0-1 vector)
     */
    void copyTree_byTaxonNames(MTree *tree, vector<string> taxon_names);
    
    /**
     *  Clean all info about link neighbours and taxa
     */
    
    void cleanAllLinkINFO(bool clean_induced_part_maps = false, TerraceNode *node = nullptr, TerraceNode *dad = nullptr);
    
    /**
     *  Insert a new taxon on given branch 
     */
    TerraceNode* insertNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, bool update_leafNode = false,bool update_brNodes=false);
    
    /**
     *  Remove one taxon
     */
    void remove_taxon(string taxon_name,bool update_leafNode = false,bool update_brNodes=false);
    
    /**
     *  print a tree, but taking into account that it can be empty, with one or two taxa or many taxa
     */
    void print_terrace_tree(bool draw = false, ostream &out = cout);
    
};

/**
 * Get a tree topology as a string
 */

string getTreeTopologyString(MTree* tree);

#endif /* terracetree_hpp */
