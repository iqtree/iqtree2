/*
 * phylokernelfma.cpp
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

    if (model_factory && model_factory->model->isSiteSpecificModel() && (params->lk_safe_scaling || leafNum >= params->numseq_safe_scaling)) {
    	// safe site-specific model
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, SAFE_LH, 4, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, SAFE_LH, 4, true, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec4d, SAFE_LH, 4, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 4, true, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, SAFE_LH, 20, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, SAFE_LH, 20, true, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec4d, SAFE_LH, 20, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 20, true, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec4d, SAFE_LH, true, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec4d, SAFE_LH, true, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec4d, SAFE_LH, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, SAFE_LH, true, true>;
            break;
        }
        return;
    }

    if (model_factory && model_factory->model->isSiteSpecificModel()) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, NORM_LH, 4, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, NORM_LH, 4, true, true>;
            computePartialLikelihoodPointer    =  &PhyloTree::computePartialLikelihoodSIMD  <Vec4d, NORM_LH, 4, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 4, true, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec4d, NORM_LH, 20, true, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec4d, NORM_LH, 20, true, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec4d, NORM_LH, 20, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 20, true, true>;
            break;
        default:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD        <Vec4d, NORM_LH, true, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD            <Vec4d, NORM_LH, true, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD      <Vec4d, NORM_LH, true, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, NORM_LH, true, true>;
            break;
        }
        return;
    }

    if ((model_factory && !model_factory->model->isReversible()) || params->kernel_nonrev) {
        // if nonreversible model
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer = &PhyloTree::computeNonrevLikelihoodBranchSIMD<Vec4d, 4, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeNonrevLikelihoodDervSIMD<Vec4d, 4, true>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodSIMD<Vec4d, 4, true>;
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
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 2, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 2, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 2, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 2, true>;
            break;
        */
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
        /*
        case 64:
            computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, SAFE_LH, 64, true>;
            computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, SAFE_LH, 64, true>;
            computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, SAFE_LH, 64, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, SAFE_LH, 64, true>;
            break;
        */
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
    /*
	case 2:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 2, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 2, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 2, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 2, true>;
		break;
    */
	case 4:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 4, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 4, true>;
        computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec4d, NORM_LH, 4, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 4, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 4, true>;
		break;
	case 20:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 20, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 20, true>;
        computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec4d, NORM_LH, 20, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 20, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 20, true>;
		break;
    /*
	case 64:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSIMD<Vec4d, NORM_LH, 64, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSIMD<Vec4d, NORM_LH, 64, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSIMD<Vec4d, NORM_LH, 64, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec4d, NORM_LH, 64, true>;
		break;
    */
	default:
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchGenericSIMD<Vec4d, NORM_LH, true>;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervGenericSIMD<Vec4d, NORM_LH, true>;
        computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenGenericSIMD<Vec4d, NORM_LH, true>;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodGenericSIMD<Vec4d, NORM_LH, true>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec4d, NORM_LH, true>;
		break;
	}
}

