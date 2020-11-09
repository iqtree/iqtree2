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
#include "node.h"

//#include <sys/time.h>
//#include <time.h>
#include <cmath>

//#define INFINITY 1000000000

/*********************************************
        class Node
 *********************************************/

Node::Node() {
    id = -1;
    height = -1;
    is_floating_interior = false;
};

Node::Node(int aid) {
    id = aid;
    height = -1;
    is_floating_interior = false;
}

Node::Node(int aid, int aname) {
    id = aid;
    char str[20];
    sprintf(str, "%d", aname);
    name = str;
    height = -1;
    is_floating_interior = false;
}

Node::Node(int aid, const char *aname) {
    id = aid;
    if (aname != nullptr ) {
        name = aname;
    }
    height = -1;
    is_floating_interior = false;
}

bool Node::isLeaf() {
    return neighbors.size() <= 1 && ! is_floating_interior;
}

bool Node::isInCherry() {
    if (this->isLeaf()) {
        if (neighbors[0]->node->isCherry()) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool Node::isCherry() {
    int num_leaves = 0;
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
        num_leaves += (*it)->node->isLeaf() ? 1 : 0;
    }
    return (num_leaves > 1);
}

int Node::degree() {
    return neighbors.size();
}

/** calculate the height of the subtree rooted at this node,
        given the dad. Also return the lowestLeaf.
        @param dad the dad of this node
        @return the leaf at the lowest level. Also modify the height, highestNei of this class.
 */
Node *Node::calcHeight(Node *dad) {
    if (isLeaf() && dad != NULL) {
        // if a leaf, but not the root
        height = 0;
        highestNei = NULL;
        return this;
    }
    // scan through all children
    height = -INFINITY;
    Node *lowestLeaf = NULL;
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if ((*it)->node != dad) {
            Node *leaf = (*it)->node->calcHeight(this);
            if ((*it)->node->height + (*it)->length > height) {
                height = (*it)->node->height + (*it)->length;
                highestNei = (*it);
                lowestLeaf = leaf;
            }
        }
    return lowestLeaf;
}

int Node::calDist(Node* partner, Node* dad, int curLen) {
    if ( this->isLeaf() && this != partner && dad != NULL )
        return 0;
    if ( this->isLeaf() && dad == NULL ) {
        return this->neighbors[0]->node->calDist(partner, this, 1);
    } else {
        Node* left = nullptr;
        Node* right = nullptr;
        for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
            if ((*it)->node != dad) {
                if (left == NULL)
                    left = (*it)->node;
                else
                    right = (*it)->node;
            }
        }
        curLen++;
        int sumLeft = 0;
        int sumRight = 0;
        if (left!=nullptr) { //JB 10-Jul-2020 Don't segfault if null
            if (left->isLeaf()) {
                if (left == partner)
                {
                    //cout << " I found you baby" << endl;
                    return curLen;
                }
            }
            else {
                sumLeft = left->calDist(partner, this, curLen);
            }
        }
        if (right != nullptr) { //JB 10-Jul-2020 Don't segfault if null
            if (right->isLeaf()) {
                if (right == partner) {
                    //cout << " I found you baby" << endl;
                    return curLen;
                }
            }
            else {
                sumRight = right->calDist(partner, this, curLen);
            }
        }
        return sumRight + sumLeft;
    }
}

/**
        efficient longest path algorithm
 */
double Node::longestPath2(Node* &node1, Node* &node2) {
    // step 1: find the farthest leaf from this node (as a leaf)
    ASSERT(isLeaf());
    node1 = calcHeight();
    // step 2: find the farthest leaf from node1
    node2 = node1->calcHeight();
    return node1->height;
}

Neighbor *Node::findNeighbor(Node *node) {
	int size = neighbors.size();
    for (int i = 0; i < size; i++) {
        if (neighbors[i]->node == node) {
            return neighbors[i];
        }
    }
    cout << "ERROR : Could not find neighbors of node " << node->id << endl;
    ASSERT(0);
    return NULL;
}

bool Node::isNeighbor(Node* node) {
    int size = neighbors.size();
    for (int i = 0; i < size; i++)
        if (neighbors[i]->node == node) return true;
    return false;
}

NeighborVec::iterator Node::findNeighborIt(Node *node) {
    NeighborVec::iterator it;
    for ( it = neighbors.begin(); it != neighbors.end(); it++) {
        if ((*it)->node == node) {
            return it;
        }
    }
    ASSERT(0);
    return it;
}

void Node::addNeighbor(Node *node, double length, int branch_id) {
    neighbors.push_back(new Neighbor(node, length, branch_id));
}

void Node::addNeighbor(Node *node, DoubleVector &length, int branch_id) {
//	assert(!length.empty());
    if (length.empty())
        addNeighbor(node, -1.0, branch_id);
    else
        addNeighbor(node, length[0], branch_id);
}

void Node::updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei) {
    ASSERT(nei_it != neighbors.end());
    *nei_it = newnei;
}

void Node::updateNeighbor(NeighborVec::iterator nei_it, Neighbor *newnei, double newlen) {
    ASSERT(nei_it != neighbors.end());
    *nei_it = newnei;
    newnei->length = newlen;
}

void Node::updateNeighbor(Node *node, Neighbor *newnei) {
    NeighborVec::iterator nei_it = findNeighborIt(node);
    ASSERT(nei_it != neighbors.end());
    *nei_it = newnei;
}

void Node::updateNeighbor(Node *node, Neighbor *newnei, double newlen) {
    NeighborVec::iterator nei_it = findNeighborIt(node);
    ASSERT(nei_it != neighbors.end());
    *nei_it = newnei;
    newnei->length = newlen;
}

void Node::updateNeighbor(Node* node, Node *newnode, double newlen) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if ((*it)->node == node) {
            (*it)->node = newnode;
            (*it)->length = newlen;
            break;
        }
}

double Node::updateNeighbor(Node* node, Node *newnode) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        if ((*it)->node == node) {
            (*it)->node = newnode;
            return (*it)->length;
        }
    return -1;
}

Node::~Node() {
    NeighborVec::reverse_iterator it;
    for (it = neighbors.rbegin(); it != neighbors.rend(); it++) {
        delete (*it);
    }
    neighbors.clear();
}

