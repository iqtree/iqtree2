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

std::string pointer_to_hex(void *ptr) {
    uintptr_t p = reinterpret_cast<uintptr_t>(ptr);
    std::stringstream s;
    s << "0x"
        << std::setfill ('0') << std::setw(sizeof(p)*2)
        << std::hex << p;
    return s.str();
}

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
    bool zeroSize = (Params::getInstance().lh_mem_save == LM_MEM_SAVE);
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
        PhyloNeighbor* nei = (PhyloNeighbor*)(*it);
        nei->partial_lh_computed = 0;
        if (make_null) {
#if (0)
            std::cout << "zapping partial_lh block " << pointer_to_hex(nei->partial_lh)
                << " to front-neighbour " << pointer_to_hex(nei)
                << " of node " << pointer_to_hex(this) << std::endl;
#endif
            nei->partial_lh = nullptr;
        }
        if (zeroSize) {
            nei->size = 0;
        }
        if (nei->node != dad) {
            ((PhyloNode*)(*it)->node)->clearAllPartialLh(make_null, this);
        }
    }
    if (dad==nullptr) {
        return;
    }
    PhyloNeighbor* backnei = dad->findNeighbor(this);
    if (backnei==nullptr) {
        return;
    }
    backnei->partial_lh_computed = 0;
    if (make_null) {
#if (0)
        std::cout << "zapping partial_lh block " << pointer_to_hex(backnei->partial_lh)
            << " to back-neighbour " << pointer_to_hex(backnei)
            << " of dad " << pointer_to_hex(dad) << std::endl;
#endif
        backnei->partial_lh = nullptr;
    }
    if (zeroSize) {
        backnei->size = 0;
    }
}

void PhyloNode::clearAllScaleNum(PhyloNode* dad) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
        PhyloNeighbor* nei = (PhyloNeighbor*)(*it);
        nei->scale_num = nullptr;
        if (nei->node != dad) {
            ((PhyloNode*)(*it)->node)->clearAllScaleNum(this);
        } else {
            PhyloNeighbor* backnei = nei->getNode()->findNeighbor(this);
            if (backnei!=nullptr) {
                backnei->scale_num = nullptr;
            }
        }
    }
}

void PhyloNode::clearAllPartialParsimony(PhyloNode* dad) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
        PhyloNeighbor* nei = (PhyloNeighbor*)(*it);
        nei->partial_pars = nullptr;
        if ((*it)->node != dad) {
            ((PhyloNode*)(*it)->node)->clearAllPartialParsimony(this);
        } else {
            PhyloNeighbor* backnei = nei->getNode()->findNeighbor(this);
            if (backnei!=nullptr) {
                backnei->partial_pars = nullptr;
            }
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

int PhyloNode::computeSize(PhyloNode *dad) {
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

PhyloNeighbor* PhyloNode::findNeighbor(Node* node)
{
    return (PhyloNeighbor*) Node::findNeighbor(node);
}


