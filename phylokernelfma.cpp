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
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 2>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 2>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 2>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 2>;
            break;
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 4>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 4>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 20>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 20>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 20>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 20>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 2>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 2>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 2>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 2>;
		break;
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 4>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 4>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 4>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 4>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 20>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 20>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 20>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 20>;
		break;
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH>;
		break;
	}
}

