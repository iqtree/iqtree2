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
	setLikelihoodComputed(false);
	for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it ++)
		if ((*it)->node != dad)
			((PhyloNeighbor*)*it)->clearForwardPartialLh(node);
}

void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
    //	PhyloNeighbor *node_nei = findNeighbor(dad);
    //	assert(node_nei);
    //	node_nei->setLikelihoodComputed(false);
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it ++) {
        PhyloNode* node = (PhyloNode*)(*it)->node;
        if (node != dad) {
            PhyloNeighbor* nei  = node->findNeighbor(this);
            nei->setLikelihoodComputed(false);
            nei->size = 0;
            node->clearReversePartialLh(this);
        }
    }
}

void PhyloNode::clearReversePartialParsimony(PhyloNode* dad) {
    FOR_EACH_PHYLO_NEIGHBOR(this, dad, it, nei) {
        PhyloNode* node = nei->getNode();
        if (node != dad ) {
            PhyloNeighbor* reverseNei = node->findNeighbor(this);
            reverseNei->setParsimonyComputed(false);
            node->clearReversePartialParsimony(this);
        }
    }
}

void PhyloNode::clearAllPartialLh(bool set_to_null, PhyloNode* dad) {
    bool zeroSize = (Params::getInstance().lh_mem_save == LM_MEM_SAVE);
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
        PhyloNeighbor* nei   = (PhyloNeighbor*)(*it);
        PhyloNode*     child = nei->getNode();
        nei->setLikelihoodComputed(false);
        if (set_to_null) {
            nei->partial_lh = nullptr;
        }
        if (zeroSize) {
            nei->size = 0;
        }
        if (child != dad) {
            child->clearAllPartialLh(set_to_null, this);
        }
    }
    if (dad==nullptr) {
        return;
    }
    PhyloNeighbor* backnei = dad->findNeighbor(this);
    if (backnei==nullptr) {
        return;
    }
    backnei->setLikelihoodComputed(false);
    if (set_to_null) {
        backnei->partial_lh = nullptr;
    }
    if (zeroSize) {
        backnei->size = 0;
    }
}

void PhyloNode::clearAllScaleNum(bool set_to_null, PhyloNode* dad) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
        PhyloNeighbor* nei   = (PhyloNeighbor*)(*it);
        PhyloNode*     child = nei->getNode();
        if (set_to_null) {
            nei->scale_num = nullptr;
        }
        if (child != dad) {
            child->clearAllScaleNum(set_to_null, this);
        } else {
            PhyloNeighbor* backnei = child->findNeighbor(this);
            if (backnei!=nullptr && set_to_null) {
                backnei->scale_num = nullptr;
            }
        }
    }
}

void PhyloNode::clearAllPartialParsimony(bool set_to_null, PhyloNode* dad) {
    for (NeighborVec::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
        PhyloNeighbor* nei = (PhyloNeighbor*)(*it);
        nei->setParsimonyComputed(false);
        if (set_to_null) {
            nei->partial_pars = nullptr;
        }
        if (nei->getNode() != dad) {
            nei->getNode()->clearAllPartialParsimony(set_to_null, this);
        } else {
            PhyloNeighbor* backnei = nei->getNode()->findNeighbor(this);
            if (backnei!=nullptr) {
                backnei->setParsimonyComputed(false);
                if (set_to_null) {
                    backnei->partial_pars = nullptr;
                }
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
    PhyloNeighbor *nei = dad->findNeighbor(this);
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

PhyloNeighbor* PhyloNode::firstNeighbor() {
    if (neighbors.empty()) {
        return nullptr;
    }
    return (PhyloNeighbor*) neighbors[0];
}

