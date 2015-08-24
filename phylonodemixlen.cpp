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
