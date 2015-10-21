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

void PhyloNeighbor::reorientPartialLh(Node *dad) {
    if (partial_lh)
        return;
    bool done = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighbor *backnei = (PhyloNeighbor*)(*it)->node->findNeighbor(node);
        if (backnei->partial_lh) {
            partial_lh = backnei->partial_lh;
            scale_num = backnei->scale_num;
            backnei->partial_lh = NULL;
            backnei->scale_num = NULL;
            backnei->partial_lh_computed &= ~1; // clear bit
            done = true;
            break;
        }
    }
    assert(done && "partial_lh is not re-oriented");
}


void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
//	PhyloNeighbor *node_nei = (PhyloNeighbor*)findNeighbor(dad);
//	assert(node_nei);
//	node_nei->partial_lh_computed = 0;
	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
		if ((*it)->node != dad) {
			((PhyloNeighbor*)(*it)->node->findNeighbor(this))->partial_lh_computed = 0;
			((PhyloNode*)(*it)->node)->clearReversePartialLh(this);
		}
}

void PhyloNode::clearAllPartialLh(bool make_null, PhyloNode *dad) {
	PhyloNeighbor *node_nei = (PhyloNeighbor*)findNeighbor(dad);
	node_nei->partial_lh_computed = 0;
	if (make_null) node_nei->partial_lh = NULL;

	node_nei = (PhyloNeighbor*)dad->findNeighbor(this);
	node_nei->partial_lh_computed = 0;
	if (make_null) node_nei->partial_lh = NULL;

	for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNode*)(*it)->node)->clearAllPartialLh(make_null, this);
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
