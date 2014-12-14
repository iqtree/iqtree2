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
	return computeLikelihoodFromBufferEigenSSE<Vec4d, 4, 4>();
}

void PhyloTree::computePartialLikelihoodEigenTipAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	computePartialLikelihoodEigenTipSSE<Vec4d, 4, 4>(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranchEigenTipAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	return computeLikelihoodBranchEigenTipSSE<Vec4d, 4, 4>(dad_branch, dad);
}

void PhyloTree::computeLikelihoodDervEigenTipAVX_DNA(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
	computeLikelihoodDervEigenTipSSE<Vec4d, 4, 4>(dad_branch, dad, df, ddf);
}

double PhyloTree::computeLikelihoodFromBufferEigenAVX_PROT() {
	return computeLikelihoodFromBufferEigenSSE<Vec4d, 4, 20>();
}

void PhyloTree::computePartialLikelihoodEigenTipAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	computePartialLikelihoodEigenTipSSE<Vec4d, 4, 20>(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranchEigenTipAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	return computeLikelihoodBranchEigenTipSSE<Vec4d, 4, 20>(dad_branch, dad);
}

void PhyloTree::computeLikelihoodDervEigenTipAVX_PROT(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
	computeLikelihoodDervEigenTipSSE<Vec4d, 4, 20>(dad_branch, dad, df, ddf);
}
