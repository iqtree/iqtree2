/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
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
#ifndef SUPERNODE_H
#define SUPERNODE_H

#include "phylonode.h"

#define FOR_EACH_ADJACENT_SUPER_NODE(mynode, mydad, it, mychild) \
    for (SuperNode* mychild=nullptr, *child2x=(SuperNode*)(mynode); child2x!=nullptr; child2x=nullptr) \
        for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
            if ((mychild = (SuperNode*)(*it)->node ) && mychild != mydad )

#define FOR_EACH_SUPER_NEIGHBOR(mynode, mydad, it, nei) \
    for (size_t ncx = (mynode)->neighbors.size(); ncx!=0; ncx=0) \
        for (SuperNeighbor* nei=nullptr, *nei2x=(SuperNeighbor*)(mynode)->neighbors[0]; nei2x!=nullptr ; nei2x=nullptr) \
            for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
                if ((nei = (SuperNeighbor*)(*it)) && nei->getNode() != (mydad) )

#define FOR_SUPER_NEIGHBOR(mynode, mydad, it) \
	for (it = SuperNeighborVec::iterator((mynode)->neighbors.begin()); \
		 it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_SUPER_NEIGHBOR_IT(mynode, mydad, it) \
	for (SuperNeighborVec::iterator it((mynode)->neighbors.begin()); \
		 it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))

#define FOR_SUPER_NEIGHBOR_DECLARE(mynode, mydad, it) \
	SuperNeighborVec::iterator it((mynode)->neighbors.begin()); \
	for (; it != (mynode)->neighbors.end(); it++) \
		if ((*it)->node != (mydad))


class SuperNode;

/**
A neighbor in a phylogenetic SUPER tree

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class SuperNeighbor : public PhyloNeighbor {

	friend class SuperNode;
	friend class PhyloSuperTree;

public:
	/**
		construct class with a node and length
		@param anode the other end of the branch
		@param alength length of branch
	*/
    SuperNeighbor(Node *anode, double alength);

	/**
		construct class with a node and length
		@param anode the other end of the branch
		@param alength length of branch
		@param aid branch ID
	*/
    SuperNeighbor(Node *anode, double alength, int aid);

    /**
     construct class with another Neighbor
     @param nei another Neighbor
     */
    SuperNeighbor(SuperNeighbor *nei);
    
    virtual SuperNode* getNode();
    
    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual SuperNeighbor* newNeighbor();

	/**
		vector of size m (m = #partitions)
	*/
	PhyloNeighborVec link_neighbors;

};

/**
Node of a super tree

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class SuperNode : public PhyloNode
{
	friend class PhyloSuperTree;

public:
	/**
		constructor 
	*/
    SuperNode();

	/**
		constructor 
		@param aid id of this node
	*/
	SuperNode(int aid);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	SuperNode(int aid, int aname);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	SuperNode(int aid, const char *aname);

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

    virtual SuperNeighbor* findNeighbor(Node* node);

	virtual SuperNeighbor* firstNeighbor();

	virtual SuperNeighbor* getNeighborByIndex(size_t index);
    
    ~SuperNode();

};

typedef SubclassPointerVector<SuperNode,     PhyloNodeVector>  SuperNodeVector;
typedef SubclassPointerVector<SuperNeighbor, PhyloNeighborVec> SuperNeighborVec;

class SuperBranch : public pair<SuperNode*, SuperNode*> {
public:
    typedef pair<SuperNode*, SuperNode*> super;
    SuperBranch(Node* a, Node* b) : super( (SuperNode*)a, (SuperNode*)b ) {
    }
    inline SuperNeighbor* lookingLeft()  const { return second->findNeighbor(first); }
    inline SuperNeighbor* lookingRight() const { return first->findNeighbor(second); }
};
#endif
