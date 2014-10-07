/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef NODE_H
#define NODE_H

#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

//#include <sys/time.h>
//#include <time.h>
#include <math.h>
#include "ncl/ncl.h"

#include "tools.h"

using namespace std;

/*--------------------------------------------------------------*/

class Node;

/**
    Neighbor list of a node in the tree
 */
class Neighbor {
public:

    /**
        the other end of the branch
     */
    Node *node;

    /**
        branch length
     */
    double length;

    /**
        branch ID
     */
    int id;

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
     */
    Neighbor(Node *anode, double alength) {
        node = anode;
        length = alength;
        id = -1;
    }

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param id branch ID
     */
    Neighbor(Node *anode, double alength, int aid) {
        node = anode;
        length = alength;
        id = aid;
    }

    /**
        construct class with another Neighbor
        @param nei another Neighbor
     */
    Neighbor(Neighbor *nei) {
        node = nei->node;
        length = nei->length;
        id = nei->id;
    }

    /**
        destructor
     */
    virtual ~Neighbor() {
    }
};

/**
    Neighbor vector
 */
typedef vector<Neighbor*> NeighborVec;

/**
    Node vector
 */
typedef vector<Node*> NodeVector;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
A Node in the tree
@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
 */
class Node {
public:
    /**
        node id.
     */
    int id;

    /**
        node name
     */
    string name;

    /**
        list of neighbors
     */
    NeighborVec neighbors;

    /**
        the height of subtree rooted at this node, used for greedy algorithm
     */
    double height;

    /**
        child of maximal height of subtree rooted at this node, used for greedy algorithm
     */
    Neighbor *highestNei;


    /**
     *      List of closest leaves to the current node.
     */
    NodeVector closestLeaves;

    /**
        constructor
     */
    Node() {
        id = -1;
        height = -1;
    };


    /**
        constructor
        @param aid id of this node
     */
    Node(int aid);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    Node(int aid, int aname);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    Node(int aid, const char *aname);

    /**
        destructor
     */
    virtual ~Node();

    /**
        used for the destructor
     */
    virtual void deleteNode();


    /**
        @return true of this is a leaf
     */
    bool isLeaf();

    /**
     *  @return TRUE if this node is a leaf in a cherry
     */
    bool isInCherry();

    /**
        @return TRUE if this node is a cherry, FALSE otherwise
     */
    bool isCherry();

    /**
        @return the number of adjacent nodes
     */
    int degree();

    /** calculate the height of the subtree rooted at this node,
        given the dad. Also return the lowestLeaf.
        @param dad the dad of this node
        @return the leaf at the lowest level. Also modify the height, highestNei of this class.
     */
    Node *calcHeight(Node *dad = NULL);


    /**
     * Calculate the distance between 2 nodes. Only for binary tree.
     * @param parner the other node
     * @return the distance
     */
    int calDist(Node *parner, Node *dad = NULL, int curLen = 0);

    /** calculate the longest path in the subtree (version 2: more efficient)
        @param node1 the returned node1 of the one end of the path
        @param node2 the returned node2 of the one end of the path
        @return the length of the longest path
     */
    double longestPath2(Node* &node1, Node* &node2);

    /**
        @param node the target node
        @return the iterator to the neighbor that has the node. If not found, return NULL
     */
    Neighbor *findNeighbor(Node *node);

    /**
     * @brief check whether the two nodes are neighbors
     * @param[in] node the other node
     */
    bool isNeighbor(Node *node);

    /**
        @param node the target node
        @return the iterator to the neighbor that has the node. If not found, return neighbors.end()
     */
    NeighborVec::iterator findNeighborIt(Node *node);

    /**
        update the neighbor node with the newnode
        @param node old neighbor node
        @param newnode new neighbor node
        @param newlen new length applied for the corresponding branch
     */
    void updateNeighbor(Node* node, Node *newnode, double newlen);

    /**
        update the neighbor node with the newnode
        @param node old neighbor node
        @param newnode new neighbor node
        @return length applied for the corresponding branch
     */
    double updateNeighbor(Node* node, Node *newnode);

    /**
        update the neighbor node with the newnode
        @param nei_it iterator to the neighbor
        @param newnei new neighbor
     */
    void updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei);

    /**
        update the neighbor node with the newnode
        @param nei_it iterator to the neighbor
        @param newnei new neighbor
        @param newlen new branch length
     */
    void updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei, double newlen);

    /**
        update the neighbor node with the newnode
        @param node old neighbor node
        @param newnei new neighbor
     */
    void updateNeighbor(Node *node, Neighbor *newnei);

    /**
        update the neighbor node with the newnode
        @param node old neighbor node
        @param newnei new neighbor
        @param newlen new branch length
     */
    void updateNeighbor(Node *node, Neighbor *newnei, double newlen);

    /**
        add a neighbor
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);
};

class Branch {
public:
    Node* node1;

    Node* node2;

    string key;

    Branch(Node* node1, Node* node2) {
        assert(node1->isNeighbor(node2));
        assert(node2->isNeighbor(node1));

        if (node1->id < node2->id) {
            this->node1 = node1;
            this->node2 = node2;
        } else {
            this->node1 = node2;
            this->node2 = node1;
        }

        key = convertIntToString(this->node1->id) + convertIntToString(this->node2->id);
    }

    inline string getKey() {
        return key;
    }
};

/*
    some macros to transverse neighbors of a node
 */
#define FOR_NEIGHBOR(mynode, mydad, it) \
	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_NEIGHBOR_IT(mynode, mydad, it) \
	for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_NEIGHBOR_DECLARE(mynode, mydad, it) \
	NeighborVec::iterator it; \
	for (it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))




/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
    nodecmp, for pruning algorithm
 */
struct nodecmp {

    /**
        nodecmp, for pruning algorithm
     */
    bool operator()(const Node* s1, const Node* s2) const {
        return (s1->neighbors[0]->length) < (s2->neighbors[0]->length);
    }
};

inline int nodenamecmp(const Node* a, const Node* b) {
    return (a->name < b->name);
}

/**
    set of leaves, sorted in ascending order by the length of the incident branch.
    For pruning algorithm
 */
typedef multiset<Node*, nodecmp> LeafSet;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
    map from leaf name to node class
 */
typedef map<const string, Node*> LeafMapName;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
    neighborcmp, for greedy algorithm
 */
struct neighborcmp {

    /**
        neighborcmp, for greedy algorithm
     */
    bool operator()(const Neighbor* s1, const Neighbor* s2) const {
        return ((s1->length + s1->node->height) > (s2->length + s2->node->height));
    }
};

/**
    set of branches, sorted in descending order by the height of the corresponding subtree.
    For greedy algorithm.
 */
typedef multiset<Neighbor*, neighborcmp> NeighborSet;


#endif
