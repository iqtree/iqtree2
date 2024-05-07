/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include <cmath>
#include "ncl/ncl.h"
#include "utils/tools.h"
#include "pda/split.h"
#include "model/modelsubst.h"
#include "alignment/sequence.h"

using namespace std;

/*--------------------------------------------------------------*/

class Node;

#define BA_BOOTSTRAP "B"
#define BA_CERTAINTY "C"


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
    *   The set of taxa underneath the neighbor
    */
    Split* split;

    /**
        Flexible attributes of the branch as key-value pairs (2018-10-08)
     */
    map<string,string> attributes;
    
    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
     */
    Neighbor(Node *anode, double alength) {
        node = anode;
        length = alength;
        id = -1;
        split = NULL;
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
        split = NULL;
    }

    /**
        construct class with another Neighbor
        @param nei another Neighbor
     */
    Neighbor(Neighbor *nei) {
        node = nei->node;
        length = nei->length;
        id = nei->id;
        split = NULL;
        attributes = nei->attributes;
    }

    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual Neighbor* newNeighbor() {
        return (new Neighbor(this));
    }
    

    /**
        destructor
     */
    virtual ~Neighbor() {
    }

    /**
        get branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual double getLength(int c) { return length; }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void getLength(DoubleVector &vec) {
        vec.resize(1);
        vec[0] = length;
//        vec.resize(1, length);
    }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec (OUT) destination branch length vector
        @param start_pos starting position in vec to copy to
    */
    virtual void getLength(DoubleVector &vec, int start_pos) { 
        ASSERT(start_pos < vec.size());
        vec[start_pos] = length;
    }

    /**
        set branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual void setLength(int c, double len) { length = len; }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void setLength(DoubleVector &vec) { 
        ASSERT(vec.size() == 1);
        length = vec[0];
    }

    /**
        set branch length by length of a Neighbor, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param nei source neigbor to copy branch lengths
        @return branch length for class c
    */
    virtual void setLength(Neighbor *nei) { length = nei->length; }

    /**
        set branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec source branch length vector
        @param start_pos starting position in vec to copy from
    */
    virtual void setLength(DoubleVector &vec, int start_pos, int num_elem) { 
        ASSERT(start_pos < vec.size());
        ASSERT(num_elem == 1);
        length = vec[start_pos];
    }

    /************** attribute processing ***************/
    /**
     @param key key name
     @param[out] value value for key
     @return true if key exists, false otherwise
     */
    template<class T>
    bool getAttr(string key, T& value) {
        map<string,string>::iterator it = attributes.find(key);
        if (it == attributes.end())
            return false;
        stringstream ss(it->second);
        ss >> value;
        return true;
    }

    /**
     put pair of (key,value) to checkpoint
     @param key key name
     @param value value
     */
    template<class T>
    void putAttr(string key, T value) {
        stringstream ss;
        ss.precision(10);
        ss << value;
        attributes[key] = ss.str();
    }

};

#define PUT_ATTR(branch, value) branch->putAttr(#value, value)
#define GET_ATTR(branch, value) branch->getAttr(#value, value)

/**
    Neighbor vector
 */
typedef vector<Neighbor*> NeighborVec;

/**
    Node vector
 */
typedef vector<Node*> NodeVector;

typedef pair<Node*, Node*> Branch;
typedef vector<Branch> BranchVector;
typedef map<int, Branch> Branches;

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
        sequence
     */
    Sequence* sequence = NULL;

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
        sequence = NULL;
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

    /**
        add a neighbor for heterotachy model
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, DoubleVector &length, int id = -1);

};
/*
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
*/

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
