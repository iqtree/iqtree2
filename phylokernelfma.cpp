/*
 * phylokernelfma.cpp
 *
 *  Created on: Sept 25, 2016
 *      Author: minh
 */


#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_exp.h"
#include "phylokernel.h"
//#include "phylokernelsafe.h"
//#include "phylokernelmixture.h"
//#include "phylokernelmixrate.h"
#include "phylokernelsitemodel.h"

#include "phylokernelnew.h"
#define KERNEL_FIX_STATES
#include "phylokernelnew.h"

#if !defined(__AVX2__) && !defined(__FMA__)
#error "You must compile this file with AVX2 or FMA enabled!"
#endif

void PhyloTree::setDotProductFMA() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec8f>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec4d>;
#endif
        dotProductDouble = &PhyloTree::dotProductSIMD<double, Vec4d>;
}

void PhyloTree::setLikelihoodKernelFMA() {
    vector_size = 4;
//    setParsimonyKernelAVX();
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

    if (params->lk_safe_scaling) {
	switch(aln->num_states) {
        case 2:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 2, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 2, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 2, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 2, true>;
            break;
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 4, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 4, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 4, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 4, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 20, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 20, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 20, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 20, true>;
            break;
        case 64:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 64, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 64, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 64, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 64, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec4d, SAFE_LH, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec4d, SAFE_LH, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec4d, SAFE_LH, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, SAFE_LH, true>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 2, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 2, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 2, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 2, true>;
		break;
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 4, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 4, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 4, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 4, true>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 20, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 20, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 20, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 20, true>;
		break;
	case 64:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 64, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 64, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 64, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 64, true>;
		break;
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec4d, NORM_LH, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec4d, NORM_LH, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec4d, NORM_LH, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, NORM_LH, true>;
		break;
	}
}

