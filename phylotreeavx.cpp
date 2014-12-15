/*
 * phylotreeavx.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */


#include "phylokernel.h"
#include "vectorclass/vectorclass.h"

#ifndef __AVX__
#error "You must compile this file with AVX enabled!"
#endif

double PhyloTree::computeLikelihoodFromBufferEigenAVX_DNA() {
	return computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>();
}

void PhyloTree::computePartialLikelihoodEigenAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	computePartialLikelihoodEigenSIMD<Vec4d, 4, 4>(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranchEigenAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	return computeLikelihoodBranchEigenSIMD<Vec4d, 4, 4>(dad_branch, dad);
}

void PhyloTree::computeLikelihoodDervEigenAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
	computeLikelihoodDervEigenSIMD<Vec4d, 4, 4>(dad_branch, dad, df, ddf);
}

double PhyloTree::computeLikelihoodFromBufferEigenAVX_PROT() {
	return computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>();
}

void PhyloTree::computePartialLikelihoodEigenAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	computePartialLikelihoodEigenSIMD<Vec4d, 4, 20>(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranchEigenAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	return computeLikelihoodBranchEigenSIMD<Vec4d, 4, 20>(dad_branch, dad);
}

void PhyloTree::computeLikelihoodDervEigenAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
	computeLikelihoodDervEigenSIMD<Vec4d, 4, 20>(dad_branch, dad, df, ddf);
}
