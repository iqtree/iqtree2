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

std::string pointer_to_hex(const void *ptr) {
    uintptr_t p = reinterpret_cast<uintptr_t>(ptr);
    std::stringstream s;
    s << "0x"
      << std::hex << p;
    return s.str();
}

void PhyloNeighbor::clearForwardPartialLh(PhyloNode* dad) {
    setLikelihoodComputed(false);
    FOR_EACH_PHYLO_NEIGHBOR(getNode(), dad, it, nei) {
        nei->clearForwardPartialLh(getNode());
    }
}

void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
    //	PhyloNeighbor *node_nei = findNeighbor(dad);
    //	assert(node_nei);
    //	node_nei->setLikelihoodComputed(false);
    FOR_EACH_ADJACENT_PHYLO_NODE(this, dad, it, node) {
        PhyloNeighbor* backNei = node->findNeighbor(this);
        backNei->setLikelihoodComputed(false);
        backNei->size = 0;
        node->clearReversePartialLh(this);
    }
}

void PhyloNode::clearReversePartialParsimony(PhyloNode* dad) {
    FOR_EACH_ADJACENT_PHYLO_NODE(this, dad, it, node) {
        PhyloNeighbor* reverseNei = node->findNeighbor(this);
        reverseNei->setParsimonyComputed(false);
        node->clearReversePartialParsimony(this);
    }
}

void PhyloNode::clearAllPartialLh(bool set_to_null, PhyloNode* dad) {
    bool zeroSize = (Params::getInstance().lh_mem_save == LM_MEM_SAVE);
    FOR_EACH_PHYLO_NEIGHBOR(this, nullptr, it, nei) {
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
    FOR_EACH_PHYLO_NEIGHBOR(this, nullptr, it, nei) {
        PhyloNode*     child = nei->getNode();
        if (set_to_null) {
            //if (nei->scale_num) std::cout << " Nulling " << pointer_to_hex(backnei) << std::endl;
            nei->scale_num = nullptr;
        }
        if (child != dad) {
            child->clearAllScaleNum(set_to_null, this);
        } else {
            PhyloNeighbor* backnei = child->findNeighbor(this);
            if (backnei!=nullptr && set_to_null) {
                //std::cout << " Nulling " << pointer_to_hex(backnei) << std::endl;
                backnei->scale_num = nullptr;
            }
        }
    }
}

void PhyloNode::clearAllPartialParsimony(bool set_to_null, PhyloNode* dad) {
    FOR_EACH_PHYLO_NEIGHBOR(this, nullptr, it, nei) {
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
    FOR_EACH_ADJACENT_PHYLO_NODE(this, dad, it, child) {
        nei->size += child->computeSize(this);
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

PhyloNeighbor* PhyloNode::getNeighborByIndex(size_t index) {
    ASSERT(index < neighbors.size());
    return (PhyloNeighbor*)neighbors[index];
}

