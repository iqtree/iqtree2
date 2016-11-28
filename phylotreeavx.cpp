/*
 * phylotreeavx.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */


#include "vectorclass/vectormath_exp.h"
#include "vectorclass/vectorclass.h"
#include "phylokernel.h"
//#include "phylokernelsafe.h"
//#include "phylokernelmixture.h"
//#include "phylokernelmixrate.h"
//#include "phylokernelsitemodel.h"

#include "phylokernelnew.h"
#include "phylokernelnonrev.h"
#define KERNEL_FIX_STATES
#include "phylokernelnew.h"
#include "phylokernelnonrev.h"

#ifndef __AVX__
#error "You must compile this file with AVX enabled!"
#endif

void PhyloTree::setParsimonyKernelAVX() {
	computeParsimonyBranchPointer = &PhyloTree::computeParsimonyBranchFastSIMD<Vec8ui>;
    computePartialParsimonyPointer = &PhyloTree::computePartialParsimonyFastSIMD<Vec8ui>;
}

void PhyloTree::setDotProductAVX() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec8f>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec4d>;
#endif
        dotProductDouble = &PhyloTree::dotProductSIMD<double, Vec4d>;
}

void PhyloTree::setLikelihoodKernelAVX() {
    vector_size = 4;
    setParsimonyKernelAVX();

    if (model_factory && model_factory->model->isSiteSpecificModel() && (params->lk_safe_scaling || leafNum >= params->numseq_safe_scaling)) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, SAFE_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, SAFE_LH, 4, false, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec4d, SAFE_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, SAFE_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, SAFE_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec4d, SAFE_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 20, false, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec4d, SAFE_LH, false, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec4d, SAFE_LH, false, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec4d, SAFE_LH, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, SAFE_LH, false, true>;
            break;
        }
        return;
    }

    if (model_factory && model_factory->model->isSiteSpecificModel()) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, NORM_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, NORM_LH, 4, false, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec4d, NORM_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, NORM_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, NORM_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec4d, NORM_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 20, false, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec4d, NORM_LH, false, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec4d, NORM_LH, false, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec4d, NORM_LH, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, NORM_LH, false, true>;
            break;
        }
        return;
    }

    if ((model_factory && !model_factory->model->isReversible()) || params->kernel_nonrev) {
        // if nonreversible model
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeNonrevLikelihoodBranchSIMD<Vec4d, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeNonrevLikelihoodDervSIMD<Vec4d, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodSIMD<Vec4d, 4>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeNonrevLikelihoodBranchGenericSIMD<Vec4d>;
            computeLikelihoodDervPointer = &PhyloTree::computeNonrevLikelihoodDervGenericSIMD<Vec4d>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodGenericSIMD<Vec4d>;
            break;
        }
        computeLikelihoodFromBufferPointer = NULL;
        return;        
    }

    if (params->lk_safe_scaling || leafNum >= params->numseq_safe_scaling) {
	switch(aln->num_states) {
        /*
        case 2:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 2>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 2>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 2>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 2>;
            break;
        */
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
        /*
        case 64:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 64>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 64>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 64>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 64>;
            break;
        */
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec4d, SAFE_LH>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec4d, SAFE_LH>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec4d, SAFE_LH>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, SAFE_LH>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
    /*
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 2>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 2>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 2>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 2>;
		break;
    */
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
    /*
	case 64:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 64>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 64>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 64>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 64>;
		break;
    */
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec4d, NORM_LH>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec4d, NORM_LH>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec4d, NORM_LH>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, NORM_LH>;
		break;
	}
}

