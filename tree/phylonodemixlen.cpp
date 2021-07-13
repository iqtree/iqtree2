//
//  phylonodemixlen.cpp
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//

#include "phylonodemixlen.h"

PhyloNeighborMixlen::PhyloNeighborMixlen
	( Node *anode, double alength, int aid ) 
	: PhyloNeighbor(anode, alength, aid) {
	lengths.clear();
}

PhyloNeighborMixlen::PhyloNeighborMixlen
	( Node *anode, DoubleVector &alength, int aid) 
	: PhyloNeighbor(anode, -1.0, aid) {
	lengths = alength;
	if (!lengths.empty()) {
		length = 0.0;
		for (int i = 0; i < lengths.size(); i++) {
			length += lengths[i];
		}
		length /= lengths.size();
	}
}

PhyloNeighborMixlen::PhyloNeighborMixlen
	( const PhyloNeighborMixlen &nei ) : PhyloNeighbor(nei) {
	lengths = nei.lengths;
}
    
PhyloNeighborMixlen* PhyloNeighborMixlen::newNeighbor() const {
	return (new PhyloNeighborMixlen(*this));
}

double PhyloNeighborMixlen::getLength(int c) { 
	if (lengths.empty()) {
		return length; 
	}
	ASSERT(c < lengths.size());
	return lengths[c];
}

void PhyloNeighborMixlen::getLength
	( DoubleVector &vec ) {
	vec = lengths;
}

void PhyloNeighborMixlen::getLength(DoubleVector &vec, 
                                    int start_pos) { 
	ASSERT(start_pos+lengths.size() <= vec.size());
	for (int i = 0; i < lengths.size(); i++)
		vec[start_pos+i] = lengths[i];
}

void PhyloNeighborMixlen::setLength(int c, double len) {
	if (lengths.empty()) {
		length = len;
		return;
	}
	ASSERT(c < lengths.size());
	lengths[c] = len;
}

void PhyloNeighborMixlen::setLength(const DoubleVector &vec) { 
	lengths = vec;
}

void PhyloNeighborMixlen::setLength(Neighbor* nei) { 
	length  = nei->length; 
	auto cast_nei = dynamic_cast<PhyloNeighborMixlen*>(nei);
	if (cast_nei!=nullptr) {
		lengths = cast_nei->lengths;
	}
}
    
void PhyloNeighborMixlen::setLength
	( const DoubleVector &vec, int start_pos, int num_elem ) { 
	ASSERT(start_pos+num_elem <= vec.size());
	lengths.clear();
	lengths.insert(lengths.begin(), vec.begin()+start_pos, 
	               vec.begin()+start_pos+num_elem);
}

PhyloNodeMixlen* PhyloNeighborMixlen::getNode() const {
	return (PhyloNodeMixlen*)node;
}
    
PhyloNodeMixlen::PhyloNodeMixlen() 
	: PhyloNode() {}

PhyloNodeMixlen::PhyloNodeMixlen(int aid) 
	: PhyloNode(aid) {}

PhyloNodeMixlen::PhyloNodeMixlen(int aid, int aname) 
	: PhyloNode(aid, aname) {}

PhyloNodeMixlen::PhyloNodeMixlen(int aid, const char *aname) 
	: PhyloNode(aid, aname) {}

void PhyloNodeMixlen::addNeighbor
	(Node *node, double length, int id) {
	auto nei = new PhyloNeighborMixlen(node, length, id);
	neighbors.push_back(nei);
}

void PhyloNodeMixlen::addNeighbor
	(Node *node, DoubleVector &length, int id) {
	if (length.empty()) {
		addNeighbor(node, -1.0, id);
	}
	else if (length.size() == 1) {
		addNeighbor(node, length[0], id);
	}
	else {
		auto nei = new PhyloNeighborMixlen(node, length, id);
		neighbors.push_back(nei);
	}
}

PhyloNeighborMixlen* PhyloNodeMixlen::firstNeighbor() const {
	if (neighbors.empty()) {
		return nullptr;
	}
	return (PhyloNeighborMixlen*)neighbors[0];
}

 PhyloNeighborMixlen* PhyloNodeMixlen::findNeighbor(PhyloNodeMixlen* node) const {
	 for (Neighbor* nei : neighbors) {
		 if (nei->getNode()==node) {
			 return dynamic_cast<PhyloNeighborMixlen*>(nei);
		 }
	 }
     cout << "ERROR : Could not find node " << node->id
          << " as a neighbor of node " << id << endl;
	 ASSERT(0);
	 return nullptr;
 }
