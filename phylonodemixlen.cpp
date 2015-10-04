//
//  phylonodemixlen.cpp
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#include "phylonodemixlen.h"

void PhyloNodeMixlen::addNeighbor(Node *node, double length, int id) {
	neighbors.push_back(new PhyloNeighborMixlen(node, length, id));
}

void PhyloNodeMixlen::addNeighbor(Node *node, DoubleVector &length, int id) {
	if (length.empty())
		addNeighbor(node, -1.0, id);
	else if (length.size() == 1)
		addNeighbor(node, length[0], id);
	else
		neighbors.push_back(new PhyloNeighborMixlen(node, length, id));
}
