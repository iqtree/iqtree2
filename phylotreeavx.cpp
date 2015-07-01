/*
 * phylotreeavx.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */


#include "phylokernel.h"
#include "phylokernelmixture.h"
#include "phylokernelmixrate.h"
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

}

void PhyloTree::setLikelihoodKernelAVX() {
    setParsimonyKernelAVX();
	switch(aln->num_states) {
	case 4:
		if (model_factory && model_factory->model->isMixture()) {
			if (model_factory->fused_mix_rate) {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixrateLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
//		        cout << "Fast-AVX-semi-mixture" << endl;
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
//		        cout << "Fast-AVX-mixture" << endl;
			}
		} else {
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
//	        cout << "Fast-AVX" << endl;
		}
		break;
	case 20:
		if (model_factory && model_factory->model->isMixture()) {
			if (model_factory->fused_mix_rate) {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixrateLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
//		        cout << "Fast-AVX-semi-mixture" << endl;
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
//		        cout << "Fast-AVX-mixture" << endl;
			}
		} else {
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
//	        cout << "Fast-AVX" << endl;
		}
		break;
	case 64: // Codon or PoMo N=11
		if (model_factory && model_factory->model->isMixture()) {
			if (model_factory->fused_mix_rate) {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigenSIMD<Vec4d, 4, 64>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigenSIMD<Vec4d, 4, 64>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigenSIMD<Vec4d, 4, 64>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixrateLikelihoodFromBufferEigenSIMD<Vec4d, 4, 64>;
//		        cout << "Fast-AVX-semi-mixture" << endl;
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec4d, 4, 64>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec4d, 4, 64>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec4d, 4, 64>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec4d, 4, 64>;
//		        cout << "Fast-AVX-mixture" << endl;
			}
		} else {
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 64>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 64>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 64>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 64>;
//	        cout << "Fast-AVX" << endl;
		}
		break;
    case 28:// PoMo N=5
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 28>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 28>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 28>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 28>;
            break;
    case 40:// PoMo N=7
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 40>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 40>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 40>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 40>;
            break;
    case 52:// PoMo N=9
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 52>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 52>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 52>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 52>;
            break;
    case 76:// PoMo N=13
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 76>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 76>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 76>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 76>;
            break;
    case 88:// PoMo N=15
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 88>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 88>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 88>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 88>;
            break;
    case 100:// PoMo N=17
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 100>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 100>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 100>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 100>;
            break;
    case 112:// PoMo N=19
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 112>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 112>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 112>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 112>;
            break;
	default:
		assert(0);
		break;
	}
}

