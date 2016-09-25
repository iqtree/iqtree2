/*
 * phylokernelavx512.cpp
 *
 *  Created on: Sept 25, 2016
 *      Author: minh
 */


#define MAX_VECTOR_SIZE 512 // for VectorClass

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


#if !defined ( __AVX512F__ ) && !defined ( __AVX512__ )
#error "You must compile this file with AVX512 enabled!"
#endif

void PhyloTree::setDotProductAVX512() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec16f>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec8d>;
#endif
        dotProductDouble = &PhyloTree::dotProductSIMD<double, Vec8d>;
}

void PhyloTree::setLikelihoodKernelAVX512() {
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
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 2>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 2>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 2>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 2>;
            break;
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 4>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 4>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 20>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 20>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 20>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 20>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 2>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 2>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 2>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 2>;
		break;
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 4>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 4>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 4>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 4>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 20>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 20>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 20>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 20>;
		break;
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH>;
		break;
	}
}

