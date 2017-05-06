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
//#include "phylokernelsitemodel.h"

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
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec8d, NORM_LH, 4, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec8d, NORM_LH, 4, true, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec8d, NORM_LH, 4, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 4, true, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec8d, NORM_LH, 20, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec8d, NORM_LH, 20, true, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec8d, NORM_LH, 20, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 20, true, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec8d, NORM_LH, true, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec8d, NORM_LH, true, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec8d, NORM_LH, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec8d, NORM_LH, true, true>;
            break;
        }
        return;
    }

    if (params->lk_safe_scaling) {
	switch(aln->num_states) {
        case 2:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 2, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 2, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 2, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 2, true>;
            break;
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 4, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 4, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 4, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 4, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, SAFE_LH, 20, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, SAFE_LH, 20, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, SAFE_LH, 20, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, SAFE_LH, 20, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec8d, SAFE_LH, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec8d, SAFE_LH, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec8d, SAFE_LH, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec8d, SAFE_LH, true>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 2, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 2, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 2, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 2, true>;
		break;
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 4, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 4, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 4, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 4, true>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec8d, NORM_LH, 20, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec8d, NORM_LH, 20, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec8d, NORM_LH, 20, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec8d, NORM_LH, 20, true>;
		break;
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec8d, NORM_LH, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec8d, NORM_LH, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec8d, NORM_LH, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec8d, NORM_LH, true>;
		break;
	}
}

