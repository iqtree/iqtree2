/*
 * phylotreeavx.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */


#include "phylokernel.h"
//#include "phylokernelmixture.h"
//#include "phylokernelmixrate.h"
#include "phylokernelsitemodel.h"
#include "vectorclass/vectorclass.h"

#ifndef __AVX__
#error "You must compile this file with AVX enabled!"
#endif

void PhyloTree::setParsimonyKernelAVX() {
	computeParsimonyBranchPointer = &PhyloTree::computeParsimonyBranchFastSIMD<Vec8ui>;
    computePartialParsimonyPointer = &PhyloTree::computePartialParsimonyFastSIMD<Vec8ui>;
}

void PhyloTree::setDotProductAVX() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec8f, 8>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec4d, 4>;
#endif

        dotProductDouble = &PhyloTree::dotProductSIMD<double, Vec4d, 4>;
}

void PhyloTree::setLikelihoodKernelAVX() {
    setParsimonyKernelAVX();
    if (model_factory && model_factory->model->isSiteSpecificModel()) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeSitemodelLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeSitemodelLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computeSitemodelPartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeSitemodelLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeSitemodelLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
            computeLikelihoodDervPointer = &PhyloTree::computeSitemodelLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
            computePartialLikelihoodPointer = &PhyloTree::computeSitemodelPartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeSitemodelLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeSitemodelLikelihoodBranchEigen;
            computeLikelihoodDervPointer = &PhyloTree::computeSitemodelLikelihoodDervEigen;
            computePartialLikelihoodPointer = &PhyloTree::computeSitemodelPartialLikelihoodEigen;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeSitemodelLikelihoodFromBufferEigen;
            break;        
        }
        return;
    }

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
	case 64:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 64>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 64>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 64>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 64>;
		break;
	default:
		assert(0);
		break;
	}
}

