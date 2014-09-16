/*
 * phylotreeeigen.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: minh
 */




#include "phylotree.h"
#include "gtrmodel.h"

/**
 * this version uses RAxML technique that stores the product of partial likelihoods and eigenvectors at node
 * for faster branch length optimization
 */
void PhyloTree::computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {

}
