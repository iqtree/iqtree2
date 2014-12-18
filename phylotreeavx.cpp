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

void PhyloTree::setDotProductAVX() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec8f, 8>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec4d, 4>;
#endif

}

void PhyloTree::setLikelihoodKernelAVX() {
	switch(aln->num_states) {
	case 4:
		computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
		computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
		computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
		computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
		break;
	case 20:
		computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
		computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
		computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
		computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
		break;
	default:
		assert(0);
		break;
	}
}

