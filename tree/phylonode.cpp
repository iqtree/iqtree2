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

bool PhyloNeighbor::isLikelihoodComputed() const {
    return ( partial_lh_computed & LIKELIHOOD_IS_COMPUTED ) != 0;
}

void PhyloNeighbor::setLikelihoodComputed(bool set) {
    if (set) {
        partial_lh_computed |= LIKELIHOOD_IS_COMPUTED;
    } else {
        partial_lh_computed &= ~LIKELIHOOD_IS_COMPUTED;
    }
}

bool PhyloNeighbor::isParsimonyComputed() const {
   return ( partial_lh_computed & PARSIMONY_IS_COMPUTED ) != 0;
}
       
void PhyloNeighbor::setParsimonyComputed(bool set) {
    if (set) {
        partial_lh_computed |= PARSIMONY_IS_COMPUTED;
    } else {
        partial_lh_computed &= ~PARSIMONY_IS_COMPUTED;
    }
}

UINT* PhyloNeighbor::get_partial_pars() {
    return partial_pars;
}

void PhyloNeighbor::clearComputedFlags() {
    partial_lh_computed = 0;
}

void PhyloNeighbor::copyComputedState(const PhyloNeighbor* donor) {
    this->partial_lh          = donor->partial_lh;
    this->scale_num           = donor->scale_num;
    this->partial_lh_computed = donor->partial_lh_computed;
    this->partial_pars        = donor->partial_pars;
    this->length              = donor->length;
}

void PhyloNode::clearReversePartialLh(PhyloNode *dad) {
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

bool PhyloNode::hasNeighbor(PhyloNode* node) {
    FOR_EACH_ADJACENT_PHYLO_NODE(this, nullptr, it, child) {
        if (child==node) {
            return true;
        }
    }
    return false;
}

PhyloNeighbor* PhyloNode::firstNeighbor() const {
    if (neighbors.empty()) {
        return nullptr;
    }
    return (PhyloNeighbor*) neighbors[0];
}

PhyloNeighbor* PhyloNode::lastNeighbor() const {
    if (neighbors.empty()) {
        return nullptr;
    }
    return (PhyloNeighbor*) neighbors[neighbors.size()-1];
}

PhyloNeighbor* PhyloNode::getNeighborByIndex(size_t index) {
    ASSERT(index < neighbors.size());
    return (PhyloNeighbor*)neighbors[index];
}

PhyloBranch::PhyloBranch() : super(nullptr, nullptr) {}

PhyloBranch::PhyloBranch(const Branch& copyMe)
    : super((PhyloNode*)(copyMe.first), (PhyloNode*)(copyMe.second)) {}

PhyloBranch::PhyloBranch(PhyloNode* left, PhyloNode* right)
    : super(left, right) {}

int PhyloBranch::getBranchID() const {
    FOR_EACH_PHYLO_NEIGHBOR(first, nullptr, it, nei) {
        if (nei->getNode() == second) {
            return nei->id;
        }
    }
    ASSERT(false);
    return -1;
}

PhyloNeighbor* PhyloBranch::getLeftNeighbor() const {
    return first->findNeighbor(second);
}

PhyloNeighbor* PhyloBranch::getRightNeighbor() const {
    return second->findNeighbor(first);
}

bool PhyloBranch::stillExists() const {
    return first->hasNeighbor(second) &&
           second->hasNeighbor(first);
}



