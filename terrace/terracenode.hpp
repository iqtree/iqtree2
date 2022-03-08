/*
 *  terracenode.hpp
 *
 *  Created on: Sep 9, 2020
 *      Author: Olga
 */

#ifndef terracenode_hpp
#define terracenode_hpp

#include <stdio.h>
#include "tree/node.h"

/**
 A neighbor in a tree, which belongs to some terrace
 
 @author Olga Chernomor <olga.chernomor@univie.ac.at>
 */
class TerraceNeighbor : public Neighbor {
    
    friend class TerraceNode;
    friend class TerraceTree;
    
public:
    /**
     construct class with a node and length
     @param anode the other end of the branch
     @param alength length of branch
     */
    TerraceNeighbor(Node *anode, double alength) : Neighbor(anode, alength) {
    }
    
    /**
     construct class with a node and length
     @param anode the other end of the branch
     @param alength length of branch
     @param aid branch ID
     */
    TerraceNeighbor(Node *anode, double alength, int aid) : Neighbor(anode, alength, aid) {
    }
    
    ~TerraceNeighbor();
    
    // TODECIDE: for simplicity you might want to have just vector of nodes, oder?
    
    /**
        vector of size m (m = #partitions) for representative tree and x for a back map from induced partition tree to the representative (for different partitions x is a different number)
     */
    NeighborVec link_neighbors;
    
    /**
     * vector to save backward maps from induced partition subtrees (low level) to induced partition trees (top level)
     */
    NeighborVec link_neighbors_lowtop_back;
    
    /**
        vector of taxa, which were collapsed onto branch. This vector will be used for partition trees.
     */
    vector<Node*> taxa_to_insert;
    
    
    /**
     * print information about link_neighbours and about taxa_to_insert
     */
    void printInfo(Node *dad);
    
    void delete_ptr_members();
    
};

/**
 Node of a tree on a terrace
 
 @author Olga Chernomor <olga.chernomor@univie.ac.at>
 */
class TerraceNode : public Node
{
    friend class TerraceTree;
    
public:
    /**
     constructor
     */
    TerraceNode();
    
    /**
     constructor
     @param aid id of this node
     */
    TerraceNode(int aid);
    
    /**
     constructor
     @param aid id of this node
     @param aname name of this nodexw
     */
    TerraceNode(int aid, int aname);
    
    /**
     constructor
     @param aid id of this node
     @param aname name of this node
     */
    TerraceNode(int aid, const char *aname);
    
    /**
     initialization
     */
    void init();
    
    /**
     add a neighbor
     @param node the neighbor node
     @param length branch length
     @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);
    
    /**
     used for the destructor
     */
    virtual void deleteNode();
    
    virtual ~TerraceNode();
    
    /*
     * For each node, save
     *  - branches below the node, that have empty images;
     *  - taxa below the node, that are absent on partition tree.
     */
    
    
    /* WARNING: there is no appropriate way to get the correct neighbors from the branch id or just the nodes. You need to pay attention to the direction!
     * Solution: get a pair of nodes, where first is always a dad and the other one is node
     */
   
    NeighborVec empty_br_node_nei;
    NeighborVec empty_br_dad_nei;
    
};


#endif /* terracenode_hpp */
