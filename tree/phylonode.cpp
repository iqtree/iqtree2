//
// C++ Implementation: phylonode
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "phylonode.h"


void PhyloNeighbor::clearForwardPartialLh(Node *dad) {
	clearPartialLh();
	for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNeighbor*)*it)->clearForwardPartialLh(node);
}

void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
//	PhyloNeighbor *node_nei = (PhyloNeighbor*)findNeighbor(dad);
//	assert(node_nei);
//	node_nei->partial_lh_computed = 0;
	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
		if ((*it)->node != dad) {
            PhyloNeighbor *nei = (PhyloNeighbor*)(*it)->node->findNeighbor(this);
			nei->partial_lh_computed = 0;
            nei->size = 0;
			((PhyloNode*)(*it)->node)->clearReversePartialLh(this);
		}
}

void PhyloNode::clearAllPartialLh(bool make_null, PhyloNode* dad) {
	PhyloNeighbor* node_nei = (PhyloNeighbor*)findNeighbor(dad);
	node_nei->partial_lh_computed = 0;
	if (make_null) node_nei->partial_lh = NULL;


	if (Params::getInstance().lh_mem_save == LM_MEM_SAVE)
		node_nei->size = 0;

	node_nei = (PhyloNeighbor*)dad->findNeighbor(this);
	node_nei->partial_lh_computed = 0;
	if (make_null) {
		node_nei->partial_lh = NULL;
	}
	if (Params::getInstance().lh_mem_save == LM_MEM_SAVE) {
		node_nei->size = 0;
	}
	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
		if ((*it)->node != dad) {
			((PhyloNode*)(*it)->node)->clearAllPartialLh(make_null, this);
		}
	}
}


PhyloNode::PhyloNode()
 : Node()
{
	init();
}


PhyloNode::PhyloNode(int aid) : Node(aid)
{
	init();
}

PhyloNode::PhyloNode(int aid, int aname) : Node (aid, aname) {
	init();
}


PhyloNode::PhyloNode(int aid, const char *aname) : Node(aid, aname) {
	init();
}

void PhyloNode::init() {
	//partial_lh = NULL;
}


void PhyloNode::addNeighbor(Node *node, double length, int id) {
	neighbors.push_back(new PhyloNeighbor(node, length, id));
}


int PhyloNode::computeSize(Node *dad) {
    PhyloNeighbor *nei = (PhyloNeighbor*)dad->findNeighbor(this);
    if (nei->size > 0)
        return nei->size;

    if (isLeaf()) {
        nei->size = 1;
        return nei->size;
    }
    nei->size = 0;
    FOR_NEIGHBOR_IT(this, dad, it) {
        nei->size += ((PhyloNode*)(*it)->node)->computeSize(this);
    }
    return nei->size;
}

