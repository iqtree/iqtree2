/*
 * phylokernelavx.cpp
 *
 *  Created on: Sept 25, 2016
 *      Author: minh
 */

#include <vectorclass/vectorclass.h>
#include <vectorclass/vectormath_exp.h>
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
    if (cost_matrix) {
        // Sankoff kernel
        computeParsimonyBranchPointer           = &PhyloTree::computeParsimonyBranchSankoffSIMD<Vec4ui>;
        computeParsimonyOutOfTreePointer        = &PhyloTree::computeParsimonyOutOfTreeSankoffSIMD<Vec4ui>;
        computePartialParsimonyPointer          = &PhyloTree::computePartialParsimonySankoffSIMD<Vec4ui>;
        computePartialParsimonyOutOfTreePointer = &PhyloTree::computePartialParsimonyOutOfTreeSankoffSIMD<Vec4ui>;

        return;
    }
    // Fitch kernel
	computeParsimonyBranchPointer           = &PhyloTree::computeParsimonyBranchFastSIMD<Vec4ui>;
    computeParsimonyOutOfTreePointer        = &PhyloTree::computeParsimonyOutOfTreeSIMD<Vec8ui>;
    computePartialParsimonyPointer          = &PhyloTree::computePartialParsimonyFastSIMD<Vec4ui>;
    computePartialParsimonyOutOfTreePointer = &PhyloTree::computePartialParsimonyOutOfTreeSIMD<Vec8ui>;
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
    bool site_model = model_factory && model_factory->model->isSiteSpecificModel();

    if (site_model && ((model_factory && !model_factory->model->isReversible()) || params->kernel_nonrev))
        outError("Site-specific model is not yet supported for nonreversible models");
    
    setParsimonyKernelSSE();
    computeLikelihoodDervMixlenPointer = NULL;

    if (site_model && safe_numeric) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 4, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD  <Vec2d, SAFE_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, SAFE_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 20, false, true>;
            break;
        default:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchGenericSIMD    <Vec2d, SAFE_LH, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervGenericSIMD      <Vec2d, SAFE_LH, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodGenericSIMD   <Vec2d, SAFE_LH, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d, false, true>;
            break;
        }
        return;
    }

    if (site_model) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 4, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 4, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD  <Vec2d, NORM_LH, 4, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 4, false, true>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 20, false, true>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 20, false, true>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, NORM_LH, 20, false, true>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 20, false, true>;
            break;
        default:
            ASSERT(0);
            break;
        }
        return;
    }

    if ((model_factory && !model_factory->model->isReversible()) || params->kernel_nonrev) {
        // if nonreversible model
        if (safe_numeric) {
        switch (aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer  = &PhyloTree::computeNonrevLikelihoodBranchSIMD <Vec2d, SAFE_LH, 4>;
            computeLikelihoodDervPointer    = &PhyloTree::computeNonrevLikelihoodDervSIMD   <Vec2d, SAFE_LH, 4>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodSIMD<Vec2d, SAFE_LH, 4>;
            break;
        default:
            computeLikelihoodBranchPointer  = &PhyloTree::computeNonrevLikelihoodBranchGenericSIMD <Vec2d, SAFE_LH>;
            computeLikelihoodDervPointer    = &PhyloTree::computeNonrevLikelihoodDervGenericSIMD   <Vec2d, SAFE_LH>;
            computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodGenericSIMD<Vec2d, SAFE_LH>;
            break;
        }
        } else {
            switch (aln->num_states) {
                case 4:
                    computeLikelihoodBranchPointer  = &PhyloTree::computeNonrevLikelihoodBranchSIMD <Vec2d, NORM_LH, 4>;
                    computeLikelihoodDervPointer    = &PhyloTree::computeNonrevLikelihoodDervSIMD   <Vec2d, NORM_LH, 4>;
                    computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodSIMD<Vec2d, NORM_LH, 4>;
                    break;
                default:
                    computeLikelihoodBranchPointer  = &PhyloTree::computeNonrevLikelihoodBranchGenericSIMD <Vec2d, NORM_LH>;
                    computeLikelihoodDervPointer    = &PhyloTree::computeNonrevLikelihoodDervGenericSIMD   <Vec2d, NORM_LH>;
                    computePartialLikelihoodPointer = &PhyloTree::computeNonrevPartialLikelihoodGenericSIMD<Vec2d, NORM_LH>;
                    break;
            }
        }
        computeLikelihoodFromBufferPointer = NULL;
        return;        
    }

    if (safe_numeric) {
	switch(aln->num_states) {
        case 4:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 4>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 4>;
            computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec2d, SAFE_LH, 4>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, SAFE_LH, 4>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 4>;
            break;
        case 20:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, SAFE_LH, 20>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, SAFE_LH, 20>;
            computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec2d, SAFE_LH, 20>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, SAFE_LH, 20>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 20>;
            break;
        default:
            computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchGenericSIMD    <Vec2d, SAFE_LH>;
            computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervGenericSIMD      <Vec2d, SAFE_LH>;
            computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenGenericSIMD<Vec2d, SAFE_LH>;
            computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodGenericSIMD   <Vec2d, SAFE_LH>;
            computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferGenericSIMD<Vec2d>;
            break;
        }
        return;
    }

	switch(aln->num_states) {
	case 4:
        computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 4>;
        computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 4>;
        computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec2d, NORM_LH, 4>;
        computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, NORM_LH, 4>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 4>;
		break;
	case 20:
        computeLikelihoodBranchPointer     = &PhyloTree::computeLikelihoodBranchSIMD    <Vec2d, NORM_LH, 20>;
        computeLikelihoodDervPointer       = &PhyloTree::computeLikelihoodDervSIMD      <Vec2d, NORM_LH, 20>;
        computeLikelihoodDervMixlenPointer = &PhyloTree::computeLikelihoodDervMixlenSIMD<Vec2d, NORM_LH, 20>;
        computePartialLikelihoodPointer    = &PhyloTree::computePartialLikelihoodSIMD   <Vec2d, NORM_LH, 20>;
        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferSIMD<Vec2d, 20>;
		break;
	default:
        ASSERT(0);
		break;
	}
}

