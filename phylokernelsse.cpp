/*
 * phylokernelavx.cpp
 *
 *  Created on: Sept 25, 2016
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


#if !defined ( __SSE2__ ) && !defined ( __x86_64__ )
#error "You must compile this file with SSE2 enabled!"
#endif

void PhyloTree::setParsimonyKernelSSE() {
	computeParsimonyBranchPointer = &PhyloTree::computeParsimonyBranchFastSIMD<Vec4ui>;
    computePartialParsimonyPointer = &PhyloTree::computePartialParsimonyFastSIMD<Vec4ui>;
}

void PhyloTree::setDotProductSSE() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec4f>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec2d>;
#endif
        dotProductDouble = &PhyloTree::dotProductSIMD<double, Vec2d>;
}

void PhyloTree::setLikelihoodKernelSSE() {
    vector_size = 2;
    setParsimonyKernelSSE();

    if (model_factory && model_factory->model->isSiteSpecificModel() && (params->lk_safe_scaling || leafNum >= params->numseq_safe_scaling)) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 4, false, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec2d, SAFE_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, SAFE_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 20, false, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec2d, SAFE_LH, false, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec2d, SAFE_LH, false, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec2d, SAFE_LH, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d, SAFE_LH, false, true>;
            break;
        }
        return;
    }

    if (model_factory && model_factory->model->isSiteSpecificModel()) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 4, false, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec2d, NORM_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, NORM_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 20, false, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec2d, NORM_LH, false, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec2d, NORM_LH, false, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec2d, NORM_LH, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d, NORM_LH, false, true>;
            break;
        }
        return;
    }

    if ((model_factory && !model_factory->model->isReversible()) || params->kernel_nonrev) {
        // if nonreversible model
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeNonrevLikelihoodBranchSIMD<Vec2d, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeNonrevLikelihoodDervSIMD<Vec2d, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodSIMD<Vec2d, 4>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeNonrevLikelihoodBranchGenericSIMD<Vec2d>;
            computeLikelihoodDervPointer = &PhyloTree::computeNonrevLikelihoodDervGenericSIMD<Vec2d>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodGenericSIMD<Vec2d>;
            break;
        }
        computeLikelihoodFromBufferPointer = NULL;
        return;        
    }

    if (params->lk_safe_scaling || leafNum >= params->numseq_safe_scaling) {
	switch(aln->num_states) {
        case 2:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, SAFE_LH, 2>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, SAFE_LH, 2>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, SAFE_LH, 2>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 2>;
            break;
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, SAFE_LH, 4>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, SAFE_LH, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, SAFE_LH, 4>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 4>;
            break;
        case 20:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, SAFE_LH, 20>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, SAFE_LH, 20>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, SAFE_LH, 20>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 20>;
            break;
        case 64:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, SAFE_LH, 64>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, SAFE_LH, 64>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, SAFE_LH, 64>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, SAFE_LH, 64>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec2d, SAFE_LH>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec2d, SAFE_LH>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec2d, SAFE_LH>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d, SAFE_LH>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, NORM_LH, 2>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, NORM_LH, 2>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, NORM_LH, 2>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 2>;
		break;
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, NORM_LH, 4>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, NORM_LH, 4>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, NORM_LH, 4>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 4>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, NORM_LH, 20>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, NORM_LH, 20>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, NORM_LH, 20>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 20>;
		break;
	case 64:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec2d, NORM_LH, 64>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec2d, NORM_LH, 64>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec2d, NORM_LH, 64>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, NORM_LH, 64>;
		break;
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec2d, NORM_LH>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec2d, NORM_LH>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec2d, NORM_LH>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d, NORM_LH>;
		break;
	}
}

