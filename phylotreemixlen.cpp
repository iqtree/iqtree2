//
//  phylotreemixlen.cpp
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#include "phylotreemixlen.h"

PhyloTreeMixlen::PhyloTreeMixlen() : IQTree() {
	mixlen = 1;
}

PhyloTreeMixlen::PhyloTreeMixlen(Alignment *aln, int mixlen) : IQTree(aln) {
	cout << "Initializing heterotachy model with " << mixlen << " mixture branch lengths" << endl;
	this->mixlen = mixlen;
}

void PhyloTreeMixlen::setMixlen(int mixlen) {

}


