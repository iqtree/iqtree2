/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "phylotree.h"
#include "phylokernel.h"
#include "phylokernelmixture.h"
#include "phylokernelmixrate.h"
#include "model/modelgtr.h"


/* BQM: to ignore all-gapp subtree at an alignment site */
//#define IGNORE_GAP_LH

//#define USING_SSE

void PhyloTree::setLikelihoodKernel(LikelihoodKernel lk) {
	if (instruction_set >= 7) {
		setDotProductAVX();
	} else {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec4f, 4>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec2d, 2>;
#endif
	}

	sse = lk;
    if (!aln || lk == LK_NORMAL) {
        computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchNaive;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervNaive;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodNaive;
        computeLikelihoodFromBufferPointer = NULL;
        return;
    }

	switch(aln->num_states) {
	case 4:
		switch(sse) {
		case LK_SSE:
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSSE<4>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSSE<4>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSSE<4>;
	        computeLikelihoodFromBufferPointer = NULL;
			break;
		case LK_EIGEN:
			if (model_factory && model_factory->model->isMixture()) {
				if (model_factory->fused_mix_rate) {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigen<4>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigen<4>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigen<4>;
					computeLikelihoodFromBufferPointer = NULL;
				} else {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigen<4>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigen<4>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigen<4>;
					computeLikelihoodFromBufferPointer = NULL;
				}
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigen<4>;
				computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigen<4>;
				computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigen<4>;
				computeLikelihoodFromBufferPointer = NULL;
			}
			break;
		case LK_EIGEN_SSE:
			if (instruction_set >= 7) {
				// CPU supports AVX
				setLikelihoodKernelAVX();
			} else {
				// CPU does not support AVX
				if (model_factory && model_factory->model->isMixture() && !model_factory->fused_mix_rate) {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec2d, 2, 4>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec2d, 2, 4>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec2d, 2, 4>;
					computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec2d, 2, 4>;
				} else {
					computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec2d, 2, 4>;
					computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec2d, 2, 4>;
					computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec2d, 2, 4>;
					computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec2d, 2, 4>;
				}
			}
			break;
		case LK_NORMAL:
			break;
		}
		break;
	case 20:
		switch(sse) {
		case LK_SSE:
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSSE<20>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSSE<20>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSSE<20>;
	        computeLikelihoodFromBufferPointer = NULL;
			break;
		case LK_EIGEN:
			if (model_factory && model_factory->model->isMixture()) {
				if (model_factory->fused_mix_rate) {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigen<20>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigen<20>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigen<20>;
					computeLikelihoodFromBufferPointer = NULL;
				} else {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigen<20>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigen<20>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigen<20>;
					computeLikelihoodFromBufferPointer = NULL;
				}
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigen<20>;
				computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigen<20>;
				computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigen<20>;
				computeLikelihoodFromBufferPointer = NULL;
			}
			break;
		case LK_EIGEN_SSE:
			if (instruction_set >= 7) {
				setLikelihoodKernelAVX();
			} else {
				if (model_factory && model_factory->model->isMixture() && !model_factory->fused_mix_rate) {
					computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec2d, 2, 20>;
					computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec2d, 2, 20>;
					computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec2d, 2, 20>;
					computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec2d, 2, 20>;
				} else {
					computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec2d, 2, 20>;
					computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec2d, 2, 20>;
					computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec2d, 2, 20>;
					computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec2d, 2, 20>;
				}
			}
			break;
		case LK_NORMAL:
			break;
		}
		break;

	case 2:
		switch(sse) {
		case LK_SSE:
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchSSE<2>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervSSE<2>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodSSE<2>;
	        computeLikelihoodFromBufferPointer = NULL;
			break;
		case LK_EIGEN:
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigen<2>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigen<2>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigen<2>;
	        computeLikelihoodFromBufferPointer = NULL;
			break;
		case LK_EIGEN_SSE:
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec2d, 2, 2>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec2d, 2, 2>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec2d, 2, 2>;
	        computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec2d, 2, 2>;
			break;
		case LK_NORMAL:
			break;
		}
		break;

	default:
		computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchNaive;
        computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervNaive;
        computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodNaive;
        computeLikelihoodFromBufferPointer = NULL;
		break;
	}
}

void PhyloTree::changeLikelihoodKernel(LikelihoodKernel lk) {
	if (sse == lk) return;
	if ((sse == LK_EIGEN || sse == LK_EIGEN_SSE) && (lk == LK_NORMAL || lk == LK_SSE)) {
		// need to increase the memory usage when changing from new kernel to old kernel
		setLikelihoodKernel(lk);
		deleteAllPartialLh();
		initializeAllPartialLh();
		clearAllPartialLH();
	} else {
		// otherwise simply assign variable sse
		setLikelihoodKernel(lk);
	}
}

/*******************************************************
 *
 * master function: wrapper for other optimized functions
 *
 ******************************************************/

void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	(this->*computePartialLikelihoodPointer)(dad_branch, dad);
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
	return (*this.*computeLikelihoodBranchPointer)(dad_branch, dad);

}

void PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
	(this->*computeLikelihoodDervPointer)(dad_branch, dad, df, ddf);
}


double PhyloTree::computeLikelihoodFromBuffer() {
	assert(current_it && current_it_back);

	if (computeLikelihoodFromBufferPointer)
		return (this->*computeLikelihoodFromBufferPointer)();
	else
		return (this->*computeLikelihoodBranchPointer)(current_it, (PhyloNode*)current_it_back->node);

}

void PhyloTree::computeTipPartialLikelihood() {
	if (tip_partial_lh_computed)
		return;
	tip_partial_lh_computed = true;
	int m, i, x, state, nstates = aln->num_states, nmixtures = model->getNMixtures();
	double *all_inv_evec = model->getInverseEigenvectors();
	assert(all_inv_evec);
	assert(tip_partial_lh);

	for (state = 0; state < nstates; state++) {
		double *this_tip_partial_lh = &tip_partial_lh[state*nstates*nmixtures];
		for (m = 0; m < nmixtures; m++) {
			double *inv_evec = &all_inv_evec[m*nstates*nstates];
			for (i = 0; i < nstates; i++)
				this_tip_partial_lh[m*nstates + i] = inv_evec[i*nstates+state];
		}
	}
	// special treatment for unknown char
	for (i = 0; i < nstates; i++) {
		double *this_tip_partial_lh = &tip_partial_lh[aln->STATE_UNKNOWN*nstates*nmixtures];
		for (m = 0; m < nmixtures; m++) {
			double *inv_evec = &all_inv_evec[m*nstates*nstates];
			double lh_unknown = 0.0;
			for (x = 0; x < nstates; x++)
				lh_unknown += inv_evec[i*nstates+x];
			this_tip_partial_lh[m*nstates + i] = lh_unknown;
		}
	}

	double lh_ambiguous;
	// ambiguous characters
	int ambi_aa[2] = {4+8, 32+64};
	switch (aln->seq_type) {
	case SEQ_DNA:
		for (state = 4; state < 18; state++) {
			int cstate = state-nstates+1;
			double *this_tip_partial_lh = &tip_partial_lh[state*nstates*nmixtures];
			for (m = 0; m < nmixtures; m++) {
				double *inv_evec = &all_inv_evec[m*nstates*nstates];
				for (i = 0; i < nstates; i++) {
					lh_ambiguous = 0.0;
					for (x = 0; x < nstates; x++)
						if ((cstate) & (1 << x))
							lh_ambiguous += inv_evec[i*nstates+x];
					this_tip_partial_lh[m*nstates+i] = lh_ambiguous;
				}
			}
		}
		break;
	case SEQ_PROTEIN:
		//map[(unsigned char)'B'] = 4+8+19; // N or D
		//map[(unsigned char)'Z'] = 32+64+19; // Q or E
		for (state = 0; state < 2; state++) {
			double *this_tip_partial_lh = &tip_partial_lh[(state+20)*nstates*nmixtures];
			for (m = 0; m < nmixtures; m++) {
				double *inv_evec = &all_inv_evec[m*nstates*nstates];
				for (i = 0; i < nstates; i++) {
					lh_ambiguous = 0.0;
					for (x = 0; x < 7; x++)
						if (ambi_aa[state] & (1 << x))
							lh_ambiguous += inv_evec[i*nstates+x];
					this_tip_partial_lh[m*nstates+i] = lh_ambiguous;
				}
			}
		}
		break;
	default:
		break;
	}


	//-------------------------------------------------------
	// initialize ptn_freq and ptn_invar
	//-------------------------------------------------------

	size_t nptn = aln->getNPattern();
	size_t maxptn = get_safe_upper_limit(nptn+model_factory->unobserved_ptns.size());
	int ptn;
	for (ptn = 0; ptn < nptn; ptn++)
		ptn_freq[ptn] = (*aln)[ptn].frequency;
	for (ptn = nptn; ptn < maxptn; ptn++)
		ptn_freq[ptn] = 0.0;

	// for +I model
	computePtnInvar();
}

void PhyloTree::computePtnInvar() {
	size_t nptn = aln->getNPattern(), ptn;
	size_t maxptn = get_safe_upper_limit(nptn+model_factory->unobserved_ptns.size());
	int nstates = aln->num_states;

    double *state_freq = aligned_alloc<double>(nstates);
    model->getStateFrequency(state_freq);
	memset(ptn_invar, 0, maxptn*sizeof(double));
	double p_invar = site_rate->getPInvar();
	if (p_invar != 0.0) {
		for (ptn = 0; ptn < nptn; ptn++)
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
					ptn_invar[ptn] = p_invar * state_freq[(int) (*aln)[ptn][0]];
			}
		// ascertmain bias correction
		for (ptn = 0; ptn < model_factory->unobserved_ptns.size(); ptn++)
			ptn_invar[nptn+ptn] = p_invar * state_freq[(int)model_factory->unobserved_ptns[ptn]];
	}
	aligned_free(state_freq);
}

/*******************************************************
 *
 * non-vectorized likelihood functions.
 * this version uses Alexis' technique that stores the
 * dot product of partial likelihoods and eigenvectors at node
 * for faster branch length optimization
 *
 ******************************************************/

template <const int nstates>
void PhyloTree::computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;

    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    PhyloNode *node = (PhyloNode*)(dad_branch->node);

	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;

		if (!tip_partial_lh_computed)
			computeTipPartialLikelihood();
		return;
	}

    size_t ptn, c;
    size_t orig_ntn = aln->size();
    size_t ncat = site_rate->getNRate();
    const size_t nstatesqr=nstates*nstates;
    size_t i, x;
    size_t block = nstates * ncat;

	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();
	assert(inv_evec && evec);
	double *eval = model->getEigenvalues();

    dad_branch->lh_scale_factor = 0.0;

	// internal node
	assert(node->degree() == 3); // it works only for strictly bifurcating tree
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
	}

	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
	}
	if ((left->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigen<nstates>(left, node);
	if ((right->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigen<nstates>(right, node);
	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;
	double partial_lh_tmp[nstates];
	double *eleft = new double[block*nstates], *eright = new double[block*nstates];

	// precompute information buffer
	for (c = 0; c < ncat; c++) {
		double *expleft = new double[nstates];
		double *expright = new double[nstates];
		double len_left = site_rate->getRate(c) * left->length;
		double len_right = site_rate->getRate(c) * right->length;
		for (i = 0; i < nstates; i++) {
			expleft[i] = exp(eval[i]*len_left);
			expright[i] = exp(eval[i]*len_right);
		}
		for (x = 0; x < nstates; x++)
			for (i = 0; i < nstates; i++) {
				eleft[c*nstatesqr+x*nstates+i] = evec[x*nstates+i] * expleft[i];
				eright[c*nstatesqr+x*nstates+i] = evec[x*nstates+i] * expright[i];
			}
		delete [] expright;
		delete [] expleft;
	}

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = new double[(aln->STATE_UNKNOWN+1)*block];
		double *partial_lh_right = new double[(aln->STATE_UNKNOWN+1)*block];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			for (x = 0; x < block; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_left[state*block+x] = vleft;
			}
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
			for (x = 0; x < block; x++) {
				double vright = 0.0;
				for (i = 0; i < nstates; i++) {
					vright += eright[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_right[state*block+x] = vright;
			}
		}

		for (x = 0; x < block; x++) {
			size_t addr = aln->STATE_UNKNOWN * block;
			partial_lh_left[addr+x] = 1.0;
			partial_lh_right[addr+x] = 1.0;
		}


		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i, partial_lh_tmp)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			int state_left = (ptn < orig_ntn) ? (aln->at(ptn))[left->node->id] : model_factory->unobserved_ptns[ptn-orig_ntn];
			int state_right = (ptn < orig_ntn) ? (aln->at(ptn))[right->node->id] : model_factory->unobserved_ptns[ptn-orig_ntn];
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				double *left = partial_lh_left + (state_left*block+c*nstates);
				double *right = partial_lh_right + (state_right*block+c*nstates);
				for (x = 0; x < nstates; x++) {
					partial_lh_tmp[x] = left[x] * right[x];
				}

				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					}
					partial_lh[c*nstates+i] = res;
				}
			}
		}
		delete [] partial_lh_right;
		delete [] partial_lh_left;
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

		// pre compute information for left tip
		double *partial_lh_left = new double[(aln->STATE_UNKNOWN+1)*block];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			for (x = 0; x < block; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_left[state*block+x] = vleft;
			}
		}
		for (x = 0; x < block; x++) {
			size_t addr = aln->STATE_UNKNOWN * block;
			partial_lh_left[addr+x] = 1.0;
		}


		double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i, partial_lh_tmp)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_right = right->partial_lh + ptn*block;
			int state_left = (ptn < orig_ntn) ? (aln->at(ptn))[left->node->id] : model_factory->unobserved_ptns[ptn-orig_ntn];
            double lh_max = 0.0;
            
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
					size_t addr = c*nstatesqr+x*nstates;
					vleft = partial_lh_left[state_left*block+c*nstates+x];
					for (i = 0; i < nstates; i++) {
						vright += eright[addr+i] * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft * (vright);
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					}
					partial_lh[c*nstates+i] = res;
                    lh_max = max(fabs(res), lh_max);
				}
			}
            if (lh_max < SCALING_THRESHOLD) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
					partial_lh[i] *= SCALING_THRESHOLD_INVER;
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
				dad_branch->scale_num[ptn] += 1;
            }


		}
		dad_branch->lh_scale_factor += sum_scale;
		delete [] partial_lh_left;

	} else {
		// both left and right are internal node

		double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i, partial_lh_tmp)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = left->partial_lh + ptn*block;
			double *partial_lh_right = right->partial_lh + ptn*block;
            double lh_max = 0.0;
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
					size_t addr = c*nstatesqr+x*nstates;
					for (i = 0; i < nstates; i++) {
						vleft += eleft[addr+i] * partial_lh_left[c*nstates+i];
						vright += eright[addr+i] * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft*vright;
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					}
					partial_lh[c*nstates+i] = res;
                    lh_max = max(lh_max, fabs(res));
				}
			}
            if (lh_max < SCALING_THRESHOLD) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
                    partial_lh[i] *= SCALING_THRESHOLD_INVER;
				}
				// unobserved const pattern will never have underflow
                sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
				dad_branch->scale_num[ptn] += 1;
            }

		}
		dad_branch->lh_scale_factor += sum_scale;

	}

	delete [] eright;
	delete [] eleft;
}

template <const int nstates>
void PhyloTree::computeLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen<nstates>(node_branch, node);
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    double *eval = model->getEigenvalues();
    assert(eval);

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
				double *theta = theta_all + ptn*block;
				double *lh_tip = tip_partial_lh + ((int)((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]))*nstates;
				for (i = 0; i < block; i++) {
					theta[i] = lh_tip[i%nstates] * partial_lh_dad[i];
				}

			}
			// ascertainment bias correction
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
		    double *partial_lh_dad = dad_branch->partial_lh;

	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    	for (i = 0; i < all_entries; i++) {
				theta_all[i] = partial_lh_node[i] * partial_lh_dad[i];
			}
	    }
		theta_computed = true;
	}

    double *val0 = new double[block];
    double *val1 = new double[block];
    double *val2 = new double[block];
	for (c = 0; c < ncat; c++) {
		double prop = site_rate->getProp(c);
		for (i = 0; i < nstates; i++) {
			double cof = eval[i]*site_rate->getRate(c);
			double val = exp(cof*dad_branch->length) * prop;
			double val1_ = cof*val;
			val0[c*nstates+i] = val;
			val1[c*nstates+i] = val1_;
			val2[c*nstates+i] = cof*val1_;
		}
	}


    double my_df = 0.0, my_ddf = 0.0, prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
		double lh_ptn = ptn_invar[ptn], df_ptn = 0.0, ddf_ptn = 0.0;
		double *theta = theta_all + ptn*block;
		for (i = 0; i < block; i++) {
			lh_ptn += val0[i] * theta[i];
			df_ptn += val1[i] * theta[i];
			ddf_ptn += val2[i] * theta[i];
		}

        assert(lh_ptn > 0.0);
        
        if (ptn < orig_nptn) {
			double df_frac = df_ptn / lh_ptn;
			double ddf_frac = ddf_ptn / lh_ptn;
			double freq = ptn_freq[ptn];
			double tmp1 = df_frac * freq;
			double tmp2 = ddf_frac * freq;
			my_df += tmp1;
			my_ddf += tmp2 - tmp1 * df_frac;
		} else {
			// ascertainment bias correction
			prob_const += lh_ptn;
			df_const += df_ptn;
			ddf_const += ddf_ptn;
		}
    }
	df = my_df;
	ddf = my_ddf;

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
    }


    delete [] val2;
    delete [] val1;
    delete [] val0;
}

template <const int nstates>
double PhyloTree::computeLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen<nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *val = new double[block];
	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		double prop = site_rate->getProp(c);
		for (i = 0; i < nstates; i++)
			val[c*nstates+i] = exp(eval[i]*len) * prop;
	}

	double prob_const = 0.0;
	memset(_pattern_lh_cat, 0, nptn*ncat*sizeof(double));

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
    	double *partial_lh_node = new double[(aln->STATE_UNKNOWN+1)*block];
    	IntVector states_dad = aln->seq_states[dad->id];
    	states_dad.push_back(aln->STATE_UNKNOWN);
    	// precompute information from one tip
    	for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
    		double *lh_node = partial_lh_node +(*it)*block;
    		double *lh_tip = tip_partial_lh + (*it)*nstates;
    		double *val_tmp = val;
			for (c = 0; c < ncat; c++) {
				for (i = 0; i < nstates; i++) {
					  lh_node[i] = val_tmp[i] * lh_tip[i];
				}
				lh_node += nstates;
				val_tmp += nstates;
			}
    	}

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c)
#endif
    	for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = ptn_invar[ptn];
			double *lh_cat = _pattern_lh_cat + ptn*ncat;
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
			int state_dad = (ptn < orig_nptn) ? (aln->at(ptn))[dad->id] : model_factory->unobserved_ptns[ptn-orig_nptn];
			double *lh_node = partial_lh_node + state_dad*block;
			for (c = 0; c < ncat; c++) {
				for (i = 0; i < nstates; i++) {
					*lh_cat += lh_node[i] * partial_lh_dad[i];
				}
				lh_node += nstates;
				partial_lh_dad += nstates;
				lh_ptn += *lh_cat;
				lh_cat++;
			}
			assert(lh_ptn > 0.0);
			if (ptn < orig_nptn) {
				lh_ptn = log(lh_ptn);
				_pattern_lh[ptn] = lh_ptn;
				tree_lh += lh_ptn * ptn_freq[ptn];
			} else {
				prob_const += lh_ptn;
			}
		}
		delete [] partial_lh_node;
    } else {
    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c)
#endif
    	for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = ptn_invar[ptn];
			double *lh_cat = _pattern_lh_cat + ptn*ncat;
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
			double *partial_lh_node = node_branch->partial_lh + ptn*block;
			double *val_tmp = val;
			for (c = 0; c < ncat; c++) {
				for (i = 0; i < nstates; i++) {
					*lh_cat +=  val_tmp[i] * partial_lh_node[i] * partial_lh_dad[i];
				}
				lh_ptn += *lh_cat;
				partial_lh_node += nstates;
				partial_lh_dad += nstates;
				val_tmp += nstates;
				lh_cat++;
			}

			assert(lh_ptn > 0.0);
            if (ptn < orig_nptn) {
				lh_ptn = log(lh_ptn);
				_pattern_lh[ptn] = lh_ptn;
				tree_lh += lh_ptn * ptn_freq[ptn];
			} else {
				prob_const += lh_ptn;
			}
		}
    }


    if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		assert(!isnan(tree_lh) && !isinf(tree_lh));
    }

	assert(!isnan(tree_lh) && !isinf(tree_lh));

    delete [] val;
    return tree_lh;
}


/************************************************************************************************
 *
 *   SSE vectorized functions of the Naive implementation
 *
 *************************************************************************************************/

template<const int NSTATES>
inline double PhyloTree::computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    PhyloNode *node = (PhyloNode*) dad_branch->node; // Node A
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad); // Node B
    assert(node_branch);
    if (!central_partial_lh)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ptn, cat, state1, state2;
    double *partial_lh_site;
    double *partial_lh_child;
    double *trans_state;
    double p_invar = site_rate->getPInvar();
    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();
    int block = numStates * numCat;

    double p_var_cat = (1.0 - p_invar) / (double) numCat;

    EIGEN_ALIGN16 double *trans_mat_orig = new double[numCat * tranSize + 1];
    double *trans_mat = trans_mat_orig;
    if (((intptr_t) trans_mat) % 16 != 0)
        trans_mat = trans_mat + 1;
    EIGEN_ALIGN16 double state_freq[NSTATES];
    model->getStateFrequency(state_freq);
    for (cat = 0; cat < numCat; cat++) {
        double *trans_cat = trans_mat + (cat * tranSize);
        model_factory->computeTransMatrix(dad_branch->length * site_rate->getRate(cat), trans_cat);
        for (state1 = 0; state1 < NSTATES; state1++) {
            double *trans_mat_state = trans_cat + (state1 * NSTATES);
            for (state2 = 0; state2 < NSTATES; state2++)
                trans_mat_state[state2] *= state_freq[state1];
        }
    }

    double prob_const = 0.0; // probability of unobserved const patterns

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, cat)
#endif
    for (ptn = 0; ptn < alnSize; ++ptn) {
        double lh_ptn = 0.0; // likelihood of the pattern
        for (cat = 0; cat < numCat; cat++) {
            partial_lh_site = node_branch->partial_lh + (ptn * block + cat * NSTATES);
            partial_lh_child = dad_branch->partial_lh + (ptn * block + cat * NSTATES);
            trans_state = trans_mat + cat * tranSize;
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_child(&partial_lh_child[0]);
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_site(&partial_lh_site[0]);
            Map<Matrix<double, NSTATES, NSTATES>, Aligned> eigen_trans_state(&trans_state[0]);
            lh_ptn += (eigen_partial_lh_child * eigen_trans_state).dot(eigen_partial_lh_site);
        }
        if (ptn < orig_alnSize) {
			lh_ptn *= p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}
			lh_ptn = log(lh_ptn);
			tree_lh += lh_ptn * (aln->at(ptn).frequency);
			_pattern_lh[ptn] = lh_ptn;
			// BQM: pattern_lh contains the LOG-likelihood, not likelihood
        } else {
			lh_ptn = lh_ptn*p_var_cat + p_invar*state_freq[(int)model_factory->unobserved_ptns[ptn-orig_alnSize]];
			prob_const += lh_ptn;

        }
    }
    if (orig_alnSize < alnSize) {
    	// ascertainment bias correction
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_alnSize; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
    }

    delete[] trans_mat_orig;
    return tree_lh;
}

template<int NSTATES>
void PhyloTree::computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    //double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    dad_branch->lh_scale_factor = 0.0;

    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();
    int block = numStates * numCat;
    size_t lh_size = alnSize * block;
    memset(dad_branch->scale_num, 0, alnSize * sizeof(UBYTE));

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_lh, 0, lh_size * sizeof(double));
        for (ptn = 0; ptn < alnSize; ++ptn) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);

            if (node->name == ROOT_NAME) {
                state = aln->STATE_UNKNOWN;
            } else if (ptn < orig_alnSize){
                state = (aln->at(ptn))[node->id];
            } else {
            	state = model_factory->unobserved_ptns[ptn-orig_alnSize];
            }

            if (state == aln->STATE_UNKNOWN) {
#ifndef KEEP_GAP_LH
                dad_branch->scale_num[ptn] = -1;
#endif
                for (int state2 = 0; state2 < block; state2++) {
                    partial_lh_site[state2] = 1.0;
                }
            } else if (state < NSTATES) {
                double *_par_lh_site = partial_lh_site + state;
                for (cat = 0; cat < numCat; cat++) {
                    *_par_lh_site = 1.0;
                    _par_lh_site += NSTATES;
                }
            } else if (aln->seq_type == SEQ_DNA) {
                // ambiguous character, for DNA, RNA
                state = state - (NSTATES - 1);
                for (int state2 = 0; state2 < NSTATES; state2++)
                    if (state & (1 << state2)) {
                        for (cat = 0; cat < numCat; cat++)
                            partial_lh_site[cat * NSTATES + state2] = 1.0;
                    }
            } else if (aln->seq_type == SEQ_PROTEIN) {
                // ambiguous character, for DNA, RNA
                state = state - (NSTATES);
                assert(state < 2);
                int state_map[2] = {4+8,32+64};
                for (int state2 = 0; state2 <= 6; state2++)
                    if (state_map[(int)state] & (1 << state2)) {
                        for (cat = 0; cat < numCat; cat++)
                            partial_lh_site[cat * NSTATES + state2] = 1.0;
                    }
            } else {
            	outError("Internal error ", __func__);
            }
        }
    } else {
        // internal node
        EIGEN_ALIGN16 double *trans_mat_orig = new double[numCat * tranSize + 2];
        double *trans_mat = trans_mat_orig;
        if (((intptr_t) trans_mat) % 16 != 0)
            trans_mat = trans_mat + 1;
        for (ptn = 0; ptn < lh_size; ++ptn)
            dad_branch->partial_lh[ptn] = 1.0;
#ifndef KEEP_GAP_LH
        for (ptn = 0; ptn < alnSize; ptn++)
            dad_branch->scale_num[ptn] = -1;
#endif
        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodSSE<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node);
            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;
            for (cat = 0; cat < numCat; cat++) {
                model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), &trans_mat[cat * tranSize]);
            }
            partial_lh_site = dad_branch->partial_lh;
            partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh;
            double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, cat, partial_lh_site, partial_lh_child)
#endif
            for (ptn = 0; ptn < alnSize; ++ptn)
#ifndef KEEP_GAP_LH
            if (((PhyloNeighbor*) (*it))->scale_num[ptn] < 0) {
#ifndef _OPENMP
                partial_lh_site += NSTATES * numCat;
                partial_lh_child += NSTATES * numCat;
#endif
            } else
#endif
            {
#ifndef KEEP_GAP_LH
                if (dad_branch->scale_num[ptn] < 0)
                dad_branch->scale_num[ptn] = 0;
#endif
#ifdef _OPENMP
                int lh_offset = ptn*block;
                partial_lh_site = dad_branch->partial_lh + lh_offset;
                partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh + lh_offset;
#endif
                dad_branch->scale_num[ptn] += ((PhyloNeighbor*) (*it))->scale_num[ptn];
                double *partial_lh_block = partial_lh_site;
                double *trans_state = trans_mat;
                bool do_scale = true;
                for (cat = 0; cat < numCat; cat++)
                {
                    MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                    MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                    MappedMat(NSTATES) ei_trans_state(trans_state);
                    ei_partial_lh_site.array() *= (ei_partial_lh_child * ei_trans_state).array();
                    partial_lh_site += NSTATES;
                    partial_lh_child += NSTATES;
                    trans_state += tranSize;
                }
                for (cat = 0; cat < block; cat++)
                if (partial_lh_block[cat] > SCALING_THRESHOLD) {
                    do_scale = false;
                    break;
                }
                if (do_scale) {
                    // unobserved const pattern will never have underflow
                    Map<VectorXd, Aligned> ei_lh_block(partial_lh_block, block);
                    ei_lh_block *= SCALING_THRESHOLD_INVER;
                    sum_scale += LOG_SCALING_THRESHOLD *  (*aln)[ptn].frequency;
                    dad_branch->scale_num[ptn] += 1;
                }
            }
            dad_branch->lh_scale_factor += sum_scale;
        }
        delete[] trans_mat_orig;
    }

    dad_branch->partial_lh_computed |= 1;
}

/****************************************************************************
 computing derivatives of likelihood function
 ****************************************************************************/
template<int NSTATES>
inline void PhyloTree::computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    //assert(node_branch);
    // swap node and dad if node is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(node_branch, node);
    df = ddf = 0.0;
    int cat = 0;
    double *partial_lh_site = node_branch->partial_lh;
    double *partial_lh_child = dad_branch->partial_lh;
    double lh_ptn; // likelihood of the pattern
    double lh_ptn_derv1;
    double lh_ptn_derv2;
    double derv1_frac;
    double derv2_frac;
    double *trans_state;
    double *derv1_state;
    double *derv2_state;
    double p_invar = site_rate->getPInvar();

    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();

    double p_var_cat = (1.0 - p_invar) / (double) numCat;
    double state_freq[NSTATES];
    model->getStateFrequency(state_freq);
    double *trans_mat_orig  = new double[numCat * tranSize + 1];
    double *trans_derv1_orig  = new double[numCat * tranSize + 1];
    double *trans_derv2_orig  = new double[numCat * tranSize + 1];
    // make alignment 16
    double *trans_mat = trans_mat_orig, *trans_derv1 = trans_derv1_orig, *trans_derv2 = trans_derv2_orig;
    if (((intptr_t) trans_mat) % 16 != 0)
        trans_mat = trans_mat + 1;
    if (((intptr_t) trans_derv1) % 16 != 0)
        trans_derv1 = trans_derv1 + 1;
    if (((intptr_t) trans_derv2) % 16 != 0)
        trans_derv2 = trans_derv2 + 1;

    int discrete_cat = site_rate->getNDiscreteRate();
    if (!site_rate->isSiteSpecificRate())
        for (cat = 0; cat < discrete_cat; cat++) {
            double *trans_cat = trans_mat + (cat * tranSize);
            double *derv1_cat = trans_derv1 + (cat * tranSize);
            double *derv2_cat = trans_derv2 + (cat * tranSize);
            double rate_val = site_rate->getRate(cat);
            model_factory->computeTransDervFreq(dad_branch->length, rate_val, state_freq, trans_cat, derv1_cat,
                    derv2_cat);
        }
    int dad_state = aln->STATE_UNKNOWN;
    double my_df = 0.0;
    double my_ddf = 0.0;
    double prob_const = 0.0, prob_const_derv1 = 0.0, prob_const_derv2 = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf,prob_const, prob_const_derv1, prob_const_derv2) \
	private(cat, partial_lh_child, partial_lh_site,\
	lh_ptn, lh_ptn_derv1, lh_ptn_derv2, derv1_frac, derv2_frac, dad_state, trans_state, derv1_state, derv2_state)
#endif
    for (int ptn = 0; ptn < alnSize; ++ptn) {
#ifdef _OPENMP
        int lh_offset = ptn*numCat*numStates;
        partial_lh_site = node_branch->partial_lh + lh_offset;
        partial_lh_child = dad_branch->partial_lh + lh_offset;
#endif
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        int padding = 0;
        dad_state = aln->STATE_UNKNOWN; // FOR TUNG: This is missing in your codes!
        if (dad->isLeaf()) {
        	if (ptn < orig_alnSize)
        		dad_state = (*aln)[ptn][dad->id];
        	else
        		dad_state = model_factory->unobserved_ptns[ptn-orig_alnSize];
        }
        padding = dad_state * NSTATES;
        if (dad_state < NSTATES) {
            //external node
            trans_state = trans_mat + padding;
            derv1_state = trans_derv1 + padding;
            derv2_state = trans_derv2 + padding;
            for (cat = 0; cat < numCat; cat++) {
                MappedVec(NSTATES)ei_partial_lh_child(partial_lh_child);
                MappedVec(NSTATES) ei_trans_state(trans_state);
                MappedVec(NSTATES) ei_derv1_state(derv1_state);
                MappedVec(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += ei_partial_lh_child.dot(ei_trans_state);
                lh_ptn_derv1 += ei_partial_lh_child.dot(ei_derv1_state);
                lh_ptn_derv2 += ei_partial_lh_child.dot(ei_derv2_state);
                partial_lh_child += NSTATES;
                partial_lh_site += NSTATES;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        } else {
            // internal node, or external node but ambiguous character
            trans_state = trans_mat;
            derv1_state = trans_derv1;
            derv2_state = trans_derv2;
            for (cat = 0; cat < numCat; cat++) {
                MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                MappedMat(NSTATES) ei_trans_state(trans_state);
                MappedMat(NSTATES) ei_derv1_state(derv1_state);
                MappedMat(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += (ei_partial_lh_child * ei_trans_state).dot(ei_partial_lh_site);
                lh_ptn_derv1 += (ei_partial_lh_child * ei_derv1_state).dot(ei_partial_lh_site);
                lh_ptn_derv2 += (ei_partial_lh_child * ei_derv2_state).dot(ei_partial_lh_site);
                partial_lh_site += NSTATES;
                partial_lh_child += NSTATES;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        }
        if (ptn < orig_alnSize) {
			lh_ptn = lh_ptn * p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}
			double pad = p_var_cat / lh_ptn;
			if (std::isinf(pad)) {
				lh_ptn_derv1 *= p_var_cat;
				lh_ptn_derv2 *= p_var_cat;
				derv1_frac = lh_ptn_derv1 / lh_ptn;
				derv2_frac = lh_ptn_derv2 / lh_ptn;
			} else {
				derv1_frac = lh_ptn_derv1 * pad;
				derv2_frac = lh_ptn_derv2 * pad;
			}
	        double freq = aln->at(ptn).frequency;
			double tmp1 = derv1_frac * freq;
			double tmp2 = derv2_frac * freq;
			my_df += tmp1;
			my_ddf += tmp2 - tmp1 * derv1_frac;
        } else {
        	lh_ptn = lh_ptn*p_var_cat + p_invar*state_freq[(int)model_factory->unobserved_ptns[ptn-orig_alnSize]];
        	prob_const += lh_ptn;
        	prob_const_derv1 += lh_ptn_derv1 * p_var_cat;
        	prob_const_derv2 += lh_ptn_derv2 * p_var_cat;
        }
    }
    if (orig_alnSize < alnSize) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	derv1_frac = prob_const_derv1 / prob_const;
    	derv2_frac = prob_const_derv2 / prob_const;
    	int nsites = aln->getNSite();
    	my_df += nsites * derv1_frac;
    	my_ddf += nsites *(derv2_frac + derv1_frac*derv1_frac);
    }

    delete[] trans_derv2_orig;
    delete[] trans_derv1_orig;
    delete[] trans_mat_orig;
    df = my_df;
    ddf = my_ddf;
}

