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
#include "gtrmodel.h"

/* BQM: to ignore all-gapp subtree at an alignment site */
//#define IGNORE_GAP_LH

#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_exp.h"
#ifdef __AVX
#define VectorClass Vec4d
#define VCSIZE 4
#pragma message "Using AVX instructions"
#else
#define VectorClass Vec2d
#define VCSIZE 2
//#pragma message "Using SS3 instructions"
#endif

//#define USING_SSE

void PhyloTree::computeTipPartialLikelihood() {
	if (tip_partial_lh_computed)
		return;
	tip_partial_lh_computed = true;
	int i, x, state, nstates = aln->num_states;

	double *evec = aligned_alloc_double(nstates*nstates);
    double *inv_evec = aligned_alloc_double(nstates*nstates);
	double **_evec = model->getEigenvectors(), **_inv_evec = model->getInverseEigenvectors();
	assert(_inv_evec && _evec);
	for (i = 0; i < nstates; i++) {
		memcpy(evec+i*nstates, _evec[i], nstates*sizeof(double));
		memcpy(inv_evec+i*nstates, _inv_evec[i], nstates*sizeof(double));
	}

	assert(tip_partial_lh);
	for (state = 0; state < nstates; state++)
		for (i = 0; i < nstates; i++)
			tip_partial_lh[state*nstates + i] = inv_evec[i*nstates+state];
	// special treatment for unknown char
	for (i = 0; i < nstates; i++) {
		double lh_unknown = 0.0;
		for (x = 0; x < nstates; x++)
			lh_unknown += inv_evec[i*nstates+x];
		tip_partial_lh[aln->STATE_UNKNOWN * nstates + i] = lh_unknown;
	}

	double lh_ambiguous;
	// ambiguous characters
	int ambi_aa[2] = {4+8, 32+64};
	switch (aln->seq_type) {
	case SEQ_DNA:
		for (state = 4; state < 18; state++) {
			int cstate = state-nstates+1;
			for (i = 0; i < nstates; i++) {
				lh_ambiguous = 0.0;
				for (x = 0; x < nstates; x++)
					if ((cstate) & (1 << x))
						lh_ambiguous += inv_evec[i*nstates+x];
				tip_partial_lh[state*nstates+i] = lh_ambiguous;
			}
		}
		break;
	case SEQ_PROTEIN:
		//map[(unsigned char)'B'] = 4+8+19; // N or D
		//map[(unsigned char)'Z'] = 32+64+19; // Q or E
		for (state = 0; state < 2; state++) {
			for (i = 0; i < nstates; i++) {
				lh_ambiguous = 0.0;
				for (x = 0; x < 7; x++)
					if (ambi_aa[state] & (1 << x))
						lh_ambiguous += inv_evec[i*nstates+x];
				tip_partial_lh[(state+20)*nstates+i] = lh_ambiguous;
			}
		}
		break;
	default:
		break;
	}
	aligned_free(inv_evec);
	aligned_free(evec);
}

/**
 * this version uses Alexis' technique that stores the dot product of partial likelihoods and eigenvectors at node
 * for faster branch length optimization
 */
template <const int nstates>
void PhyloTree::computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);
	if (node->isLeaf()) {
		if (!tip_partial_lh_computed)
			computeTipPartialLikelihood();
		return;
	}

    size_t ptn, c;
    size_t nptn = aln->size();
    size_t ncat = site_rate->getNRate();
    //size_t nstates = aln->num_states;
    //const size_t ncat = 4;
    //const size_t nstates = 4;
    const size_t nstatesqr=nstates*nstates;
    size_t i, x;
    //const size_t block = nstates * ncat;
    size_t block = nstates * ncat;

    double *partial_lh = dad_branch->partial_lh;

    double *evec = new double[nstates*nstates];
    double *inv_evec = new double[nstates*nstates];
	double **_evec = model->getEigenvectors(), **_inv_evec = model->getInverseEigenvectors();
	assert(_inv_evec && _evec);
	for (i = 0; i < nstates; i++) {
		memcpy(evec+i*nstates, _evec[i], nstates*sizeof(double));
		memcpy(inv_evec+i*nstates, _inv_evec[i], nstates*sizeof(double));
	}
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
		computePartialLikelihoodEigen<nstates>(left, node, pattern_scale);
	if ((right->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigen<nstates>(right, node, pattern_scale);
	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;
	double *partial_lh_tmp = new double[nstates];
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

	MappedMat(nstates) ei_inv_evec(inv_evec);
	MappedRowVec(nstates) ei_partial_lh_tmp(partial_lh_tmp);

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = new double[(aln->STATE_UNKNOWN+1)*block];
		double *partial_lh_right = new double[(aln->STATE_UNKNOWN+1)*block];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
#ifdef USING_SSE
			double *eleft_tmp = eleft;
			double *partial_lh_left_tmp = partial_lh_left + state*block;
			MappedRowVec(nstates) ei_tip_partial_lh(tip_partial_lh+state*nstates);
			for (c = 0; c < ncat; c++) {
				MappedMat(nstates) ei_eleft(eleft_tmp);
				MappedRowVec(nstates) ei_partial_lh_left(partial_lh_left_tmp);
				ei_partial_lh_left = ei_tip_partial_lh * ei_eleft;
				eleft_tmp += nstatesqr;
				partial_lh_left_tmp += nstates;

			}
#else
			for (x = 0; x < block; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_left[state*block+x] = vleft;
			}
#endif
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
#ifdef USING_SSE
			double *eright_tmp = eright;
			double *partial_lh_right_tmp = partial_lh_right + state*block;
			MappedRowVec(nstates) ei_tip_partial_lh(tip_partial_lh+state*nstates);
			for (c = 0; c < ncat; c++) {
				MappedMat(nstates) ei_eright(eright_tmp);
				MappedRowVec(nstates) ei_partial_lh_right(partial_lh_right_tmp);
				ei_partial_lh_right = ei_tip_partial_lh * ei_eright;
				eright_tmp += nstatesqr;
				partial_lh_right_tmp += nstates;

			}
#else
			for (x = 0; x < block; x++) {
				double vright = 0.0;
				for (i = 0; i < nstates; i++) {
					vright += eright[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_right[state*block+x] = vright;
			}
#endif
		}

		for (x = 0; x < block; x++) {
			size_t addr = aln->STATE_UNKNOWN * block;
			partial_lh_left[addr+x] = 1.0;
			partial_lh_right[addr+x] = 1.0;
		}


		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

		for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];
			int state_right = (aln->at(ptn))[right->node->id];
#ifdef USING_SSE

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				double *left = partial_lh_left + (state_left*block+c*nstates);
				double *right = partial_lh_right + (state_right*block+c*nstates);
				MappedRowVec(nstates) ei_partial_lh(partial_lh);
				MappedRowVec(nstates) ei_left(left);
				MappedRowVec(nstates) ei_right(right);

				ei_partial_lh_tmp.array() = ei_left.array() * ei_right.array();
				ei_partial_lh =  ei_partial_lh_tmp * ei_inv_evec;

				partial_lh = partial_lh + nstates;
			}
#else
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
			partial_lh += block;
#endif
		}
		delete [] partial_lh_right;
		delete [] partial_lh_left;
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

		// pre compute information for left tip
		double *partial_lh_left = new double[(aln->STATE_UNKNOWN+1)*block];
		double *partial_lh_right = right->partial_lh;

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
#ifdef USING_SSE
			double *eleft_tmp = eleft;
			double *partial_lh_left_tmp = partial_lh_left + state*block;
			MappedRowVec(nstates) ei_tip_partial_lh(tip_partial_lh+state*nstates);
			for (c = 0; c < ncat; c++) {
				MappedMat(nstates) ei_eleft(eleft_tmp);
				MappedRowVec(nstates) ei_partial_lh_left(partial_lh_left_tmp);
				ei_partial_lh_left = ei_tip_partial_lh * ei_eleft;
				eleft_tmp += nstatesqr;
				partial_lh_left_tmp += nstates;

			}
#else
			for (x = 0; x < block; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
				}
				partial_lh_left[state*block+x] = vleft;
			}
#endif
		}
		for (x = 0; x < block; x++) {
			size_t addr = aln->STATE_UNKNOWN * block;
			partial_lh_left[addr+x] = 1.0;
		}



		double sum_scale = 0.0;
		for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];
#ifdef USING_SSE
			double *partial_lh_block = partial_lh;
			double *eright_tmp = eright;
			double *partial_lh_left_tmp = partial_lh_left + state_left*block;
			for (c = 0; c < ncat; c++) {
				MappedMat(nstates) ei_eright(eright_tmp);
				MappedRowVec(nstates) ei_partial_lh_right(partial_lh_right);
				MappedRowVec(nstates) ei_partial_lh_left(partial_lh_left_tmp);
				MappedRowVec(nstates) ei_partial_lh(partial_lh);
				ei_partial_lh_tmp.array() = (ei_partial_lh_right * ei_eright).array() * ei_partial_lh_left.array();
				ei_partial_lh =  ei_partial_lh_tmp * ei_inv_evec;
				partial_lh_right += nstates;
				partial_lh += nstates;
				eright_tmp += nstatesqr;
				partial_lh_left_tmp += nstates;
			}

            // check if one should scale partial likelihoods
			bool do_scale = true;
            for (i = 0; i < block; i++)
				if (fabs(partial_lh_block[i]) > SCALING_THRESHOLD) {
					do_scale = false;
					break;
				}
            if (do_scale) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
					partial_lh_block[i] /= SCALING_THRESHOLD;
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

#else
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
				}
			}
            // check if one should scale partial likelihoods
			bool do_scale = true;
            for (i = 0; i < block; i++)
				if (fabs(partial_lh[i]) > SCALING_THRESHOLD) {
					do_scale = false;
					break;
				}
            if (do_scale) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
					partial_lh[i] /= SCALING_THRESHOLD;
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

			partial_lh += block;
			partial_lh_right += block;
#endif

		}
		dad_branch->lh_scale_factor += sum_scale;
		delete [] partial_lh_left;

	} else {
		// both left and right are internal node
		double *partial_lh_left = left->partial_lh;
		double *partial_lh_right = right->partial_lh;

		double sum_scale = 0.0;
		for (ptn = 0; ptn < nptn; ptn++) {
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];
#ifdef USING_SSE
			double *partial_lh_block = partial_lh;
			double *eleft_tmp = eleft;
			double *eright_tmp = eright;
			for (c = 0; c < ncat; c++) {
				MappedMat(nstates) ei_eleft(eleft_tmp);
				MappedMat(nstates) ei_eright(eright_tmp);
				MappedRowVec(nstates) ei_partial_lh_left(partial_lh_left);
				MappedRowVec(nstates) ei_partial_lh_right(partial_lh_right);
				MappedRowVec(nstates) ei_partial_lh(partial_lh);
				ei_partial_lh_tmp.array() = (ei_partial_lh_left * ei_eleft).array() * (ei_partial_lh_right * ei_eright).array();
				ei_partial_lh =  ei_partial_lh_tmp * ei_inv_evec;

				partial_lh_left += nstates;
				partial_lh_right += nstates;
				partial_lh += nstates;
				eleft_tmp += nstatesqr;
				eright_tmp += nstatesqr;
			}
            // check if one should scale partial likelihoods
			bool do_scale = true;
            for (i = 0; i < block; i++)
				if (fabs(partial_lh_block[i]) > SCALING_THRESHOLD) {
					do_scale = false;
					break;
				}
            if (do_scale) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
					partial_lh_block[i] /= SCALING_THRESHOLD;
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

#else
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
				}
			}

            // check if one should scale partial likelihoods
			bool do_scale = true;
            for (i = 0; i < block; i++)
				if (fabs(partial_lh[i]) > SCALING_THRESHOLD) {
					do_scale = false;
					break;
				}
            if (do_scale) {
				// now do the likelihood scaling
				for (i = 0; i < block; i++) {
					partial_lh[i] /= SCALING_THRESHOLD;
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

			partial_lh += block;
			partial_lh_left += block;
			partial_lh_right += block;
#endif
		}
		dad_branch->lh_scale_factor += sum_scale;

	}

	delete [] eright;
	delete [] eleft;
	delete [] partial_lh_tmp;
	delete [] inv_evec;
	delete [] evec;
}

template <const int nstates>
double PhyloTree::computeLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
    df = ddf = 0.0;
    size_t ncat = site_rate->getNRate();
//    const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
//    const size_t block = ncat * nstates;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization
	    double *partial_lh_dad = dad_branch->partial_lh;
    	double *theta = theta_all;

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
			for (ptn = 0; ptn < nptn; ptn++) {
				int state_dad = (aln->at(ptn))[dad->id];
#ifdef USING_SSE
				MappedRowVec(nstates) ei_tip_partial_lh(tip_partial_lh + state_dad*nstates);
				for (c = 0; c < ncat; c++) {
					MappedRowVec(nstates) ei_theta(theta);
					MappedRowVec(nstates) ei_partial_lh_dad(partial_lh_dad);
					ei_theta.array() = ei_tip_partial_lh.array() * ei_partial_lh_dad.array();
					theta += nstates;
					partial_lh_dad += nstates;
				}

#else
				for (i = 0; i < block; i++) {
					theta[i] = tip_partial_lh[state_dad*nstates+i%nstates] * partial_lh_dad[i];
				}

				//partial_lh_node += nstates;
				partial_lh_dad += block;
				theta += block;
#endif
			}
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
	    	size_t all_entries = nptn*block;
#ifdef USING_SSE
	    	MappedArrDyn ei_theta(theta, all_entries);
	    	MappedArrDyn ei_partial_lh_node(partial_lh_node, all_entries);
	    	MappedArrDyn ei_partial_lh_dad(partial_lh_dad, all_entries);
	    	ei_theta = ei_partial_lh_node * ei_partial_lh_dad;
#else
			for (i = 0; i < all_entries; i++) {
				theta[i] = partial_lh_node[i] * partial_lh_dad[i];
			}
#endif
	    }
		theta_computed = true;
	}

    double *val0 = new double[block];
    double *val1 = new double[block];
    double *val2 = new double[block];
	for (c = 0; c < ncat; c++) {
		for (i = 0; i < nstates; i++) {
			double cof = eval[i]*site_rate->getRate(c);
			double val = exp(cof*dad_branch->length);
			double val1_ = cof*val;
			val0[c*nstates+i] = val;
			val1[c*nstates+i] = val1_;
			val2[c*nstates+i] = cof*val1_;
		}
	}


    double *theta = theta_all;
	for (ptn = 0; ptn < nptn; ptn++) {
		double lh_ptn = 0.0, df_ptn = 0.0, ddf_ptn = 0.0;
#ifdef USING_SSE
		double *val0_tmp = val0;
		double *val1_tmp = val1;
		double *val2_tmp = val2;
		for (c = 0; c < ncat; c++) {
			MappedVec(nstates) ei_theta(theta);
			MappedVec(nstates) ei_val0(val0_tmp);
			MappedVec(nstates) ei_val1(val1_tmp);
			MappedVec(nstates) ei_val2(val2_tmp);

			lh_ptn += ei_val0.dot(ei_theta);
			df_ptn += ei_val1.dot(ei_theta);
			ddf_ptn += ei_val2.dot(ei_theta);
			theta += nstates;
			val0_tmp += nstates;
			val1_tmp += nstates;
			val2_tmp += nstates;
		}
#else
		for (i = 0; i < block; i++) {
			lh_ptn += val0[i] * theta[i];
			df_ptn += val1[i] * theta[i];
			ddf_ptn += val2[i] * theta[i];
		}
		theta += block;

#endif
		lh_ptn *= p_var_cat;
		df_ptn *= p_var_cat;
		ddf_ptn *= p_var_cat;
		double df_frac = df_ptn / lh_ptn;
		double ddf_frac = ddf_ptn / lh_ptn;
		double freq = (*aln)[ptn].frequency;
		double tmp1 = df_frac * freq;
		double tmp2 = ddf_frac * freq;
		df += tmp1;
		ddf += tmp2 - tmp1 * df_frac;
		assert(lh_ptn > 0.0);
		lh_ptn = log(lh_ptn);
		tree_lh += lh_ptn * freq;
		_pattern_lh[ptn] = lh_ptn;
    }

//	cout.unsetf(ios::fixed);
//	for (ptn = 0; ptn < nptn; ptn++)
//		cout << _pattern_lh[ptn] << " ";
//	cout << endl;
//	abort();


    delete [] val2;
    delete [] val1;
    delete [] val0;
    return tree_lh;
}

template <const int nstates>
double PhyloTree::computeLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
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
//    const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
//    const size_t block = ncat * nstates;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *partial_lh_dad = dad_branch->partial_lh;
    double *partial_lh_node = node_branch->partial_lh;
    double *val = new double[block];
	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		for (i = 0; i < nstates; i++)
			val[c*nstates+i] = exp(eval[i]*len);
	}
    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0;
			int state_dad = (aln->at(ptn))[dad->id];
#ifdef USING_SSE
			double *val_tmp = val;
			MappedVec(nstates) ei_tip_partial_lh(tip_partial_lh+state_dad*nstates);
			for (c = 0; c < ncat; c++){
				MappedVec(nstates) ei_val(val_tmp);
				MappedVec(nstates) ei_partial_lh_dad(partial_lh_dad);
				lh_ptn += (ei_val.array() * ei_tip_partial_lh.array() * ei_partial_lh_dad.array()).sum();
				partial_lh_dad += nstates;
				val_tmp += nstates;
			}
#else
			for (c = 0; c < ncat; c++)
				for (i = 0; i < nstates; i++) {
					lh_ptn +=  val[c*nstates+i] * tip_partial_lh[state_dad*nstates+i] * partial_lh_dad[c*nstates+i];
				}
			partial_lh_dad += block;
#endif
			lh_ptn *= p_var_cat;
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
		}
    } else {
    	// both dad and node are internal nodes
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0;
#ifdef USING_SSE
			double *val_tmp = val;
			for (c = 0; c < ncat; c++){
				MappedVec(nstates) ei_val(val_tmp);
				MappedVec(nstates) ei_partial_lh_dad(partial_lh_dad);
				MappedVec(nstates) ei_partial_lh_node(partial_lh_node);
				lh_ptn += (ei_val.array() * ei_partial_lh_node.array() * ei_partial_lh_dad.array()).sum();
				partial_lh_dad += nstates;
				partial_lh_node += nstates;
				val_tmp += nstates;
			}

#else
			for (i = 0; i < block; i++)
				lh_ptn +=  val[i] * partial_lh_node[i] * partial_lh_dad[i];
			partial_lh_node += block;
			partial_lh_dad += block;
#endif
			lh_ptn *= p_var_cat;
			assert(lh_ptn > 0.0);
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
		}
    }
    if (pattern_lh)
        memmove(pattern_lh, _pattern_lh, aln->size() * sizeof(double));
    delete [] val;
    return tree_lh;
}

/************************************************************************************************
 *
 *   SSE vectorized versions of above functions
 *
 *************************************************************************************************/

static inline Vec2d horizontal_add(Vec2d const & a, Vec2d const & b) {
#if  INSTRSET >= 3  // SSE3
    return _mm_hadd_pd(a,b);
#else
#error "Unsupported yet"
#endif
}

static inline double horizontal_max(Vec2d const &a) {
	double t[2]  __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)));
	a.store_a(t);
	if (t[0] > t[1]) return t[0];
	return t[1];
}

/***
 * AVX support codes
 */
#ifdef __AVX

static inline Vec2d horizontal_add(Vec4d const & a, Vec4d const & b) {
	__m256d t1 = _mm256_hadd_pd(a,b);
	__m128d low = _mm256_castpd256_pd128(t1);
	__m128d high = _mm256_extractf128_pd(t1,1);
	return _mm_add_pd(low, high);
}

/*
static inline Vec4d horizontal_add(Vec4d const & a, Vec4d const & b, Vec4d const & c, Vec4d const & d) {
	__m256d t1 = _mm256_hadd_pd(a,b);
	__m256d t2 = _mm256_hadd_pd(c,d);
	__m128d low = _mm256_castpd256_pd128(t1);
	__m128d high = _mm256_extractf128_pd(t1,1);
	return _mm_add_pd(low, high);
}*/

static inline double horizontal_max(Vec4d const &a) {
	__m128d low = _mm256_castpd256_pd128(a);
	__m128d high = _mm256_extractf128_pd(a,1);
	low = max(low, high);
	return horizontal_max(low);
}

#endif

/* required to compute the absolute values of double precision numbers with SSE3 */

const union __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)))
{
  uint64_t i[2];
  __m128d m;
} absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};

template <const int nstates>
void PhyloTree::computePartialLikelihoodEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);
	if (node->isLeaf()) {
		if (!tip_partial_lh_computed)
			computeTipPartialLikelihood();
		return;
	}

    size_t ptn, c;
    size_t nptn = aln->size();
    size_t ncat = site_rate->getNRate();
    //size_t nstates = aln->num_states;
    //const size_t ncat = 4;
    //const size_t nstates = 4;
    const size_t nstatesqr=nstates*nstates;
    size_t i, x;
    //const size_t block = nstates * ncat;
    size_t block = nstates * ncat;

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
		computePartialLikelihoodEigenTipSSE<nstates>(left, node, pattern_scale);
	if ((right->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigenTipSSE<nstates>(right, node, pattern_scale);

    double *partial_lh = dad_branch->partial_lh;

    double *evec = aligned_alloc_double(nstates*nstates);
    double *inv_evec = aligned_alloc_double(nstates*nstates);
	double **_evec = model->getEigenvectors(), **_inv_evec = model->getInverseEigenvectors();

	VectorClass vc_inv_evec[nstates*nstates/VCSIZE];
	assert(_inv_evec && _evec);
	for (i = 0; i < nstates; i++) {
		memcpy(evec+i*nstates, _evec[i], nstates*sizeof(double));
		memcpy(inv_evec+i*nstates, _inv_evec[i], nstates*sizeof(double));
		for (x = 0; x < nstates/VCSIZE; x++)
			vc_inv_evec[i*nstates/VCSIZE+x].load_a(&inv_evec[i*nstates+x*VCSIZE]);
	}
	double *eval = model->getEigenvalues();

#ifdef __AVX
	Vec2d v2d_inv_evec[nstates*nstates/2];
	for (i = 0; i < nstates; i++) {
		for (x = 0; x < nstates/2; x++)
			v2d_inv_evec[i*nstates/2+x].load_a(&inv_evec[i*nstates+x*2]);
	}
#else
#define v2d_inv_evec vc_inv_evec
#endif

    dad_branch->lh_scale_factor = 0.0;
	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;

//	double *partial_lh_tmp = aligned_alloc_double(nstates);
	double *eleft = aligned_alloc_double(block*nstates);
	double *eright = aligned_alloc_double(block*nstates);
	VectorClass vleft, vleft2, vright, vright2, res1, res2;


	// precompute information buffer
	for (c = 0; c < ncat; c++) {
		double *expleft = aligned_alloc_double(nstates);
		double *expright = aligned_alloc_double(nstates);
		double len_left = site_rate->getRate(c) * left->length;
		double len_right = site_rate->getRate(c) * right->length;
		for (i = 0; i < nstates; i++) {
			expleft[i] = exp(eval[i]*len_left);
			expright[i] = exp(eval[i]*len_right);
		}
		for (x = 0; x < nstates; x++)
			for (i = 0; i < nstates; i+=VCSIZE) {
				res1.load_a(&evec[x*nstates+i]);
				(res1 * VectorClass().load_a(&expleft[i])) . store_a(&eleft[c*nstatesqr+x*nstates+i]);
				(res1 * VectorClass().load_a(&expright[i])) . store_a(&eright[c*nstatesqr+x*nstates+i]);
			}
		aligned_free(expright);
		aligned_free(expleft);
	}

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = aligned_alloc_double((aln->STATE_UNKNOWN+1)*block);
		double *partial_lh_right = aligned_alloc_double((aln->STATE_UNKNOWN+1)*block);

		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			for (i = 0; i < nstates/VCSIZE; i++)
				vc_partial_lh_tmp[i].load_a(&tip_partial_lh[state*nstates+i*VCSIZE]);
			for (x = 0; x < block; x+=2) {
				vleft = VectorClass().load_a(&eleft[x*nstates]) * vc_partial_lh_tmp[0];
				vleft2 = VectorClass().load_a(&eleft[(x+1)*nstates]) * vc_partial_lh_tmp[0];
				for (i = VCSIZE; i < nstates; i+=VCSIZE) {
					vleft += VectorClass().load_a(&eleft[x*nstates+i]) * vc_partial_lh_tmp[i/VCSIZE];
					vleft2 += VectorClass().load_a(&eleft[(x+1)*nstates+i]) * vc_partial_lh_tmp[i/VCSIZE];
				}
				horizontal_add(vleft, vleft2).store_a(&partial_lh_left[state*block+x]);
			}
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
			for (i = 0; i < nstates/VCSIZE; i++)
				vc_partial_lh_tmp[i].load_a(&tip_partial_lh[state*nstates+i*VCSIZE]);
			for (x = 0; x < block; x+=2) {
				vright = VectorClass().load_a(&eright[x*nstates]) * vc_partial_lh_tmp[0];
				vright2 = VectorClass().load_a(&eright[(x+1)*nstates]) * vc_partial_lh_tmp[0];
				for (i = VCSIZE; i < nstates; i+=VCSIZE) {
					vright += VectorClass().load_a(&eright[x*nstates+i]) * vc_partial_lh_tmp[i/VCSIZE];
					vright2 += VectorClass().load_a(&eright[(x+1)*nstates+i]) * vc_partial_lh_tmp[i/VCSIZE];
				}
				horizontal_add(vright, vright2).store_a(&partial_lh_right[state*block+x]);
			}
		}

		size_t addr_unknown = aln->STATE_UNKNOWN * block;
		for (x = 0; x < block; x++) {
			partial_lh_left[addr_unknown+x] = 1.0;
			partial_lh_right[addr_unknown+x] = 1.0;
		}


		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

	    /*
	    if (nstates == 4 && ncat == 0) {
	    	// Codes taken and adapted from PLL
	    	double *x3, *uX1, *uX2;
	    	int state_left, state_right;
	    	size_t j, k;
	    	__m128d EVV[8];
	    	  for(k = 0; k < 8; k++)
	    	    EVV[k] = _mm_load_pd(&inv_evec[k * 2]);

	    	for (i = 0; i < nptn; i++)
	    	{
	    		x3 = &partial_lh[i * 16];
	    		state_left = (aln->at(i))[left->node->id];
	    		state_right = (aln->at(i))[right->node->id];

	    		uX1 = &partial_lh_left[16 * state_left];
	    		uX2 = &partial_lh_right[16 * state_right];

	    		for (j = 0; j < 4; j++)
	    		{
	    			__m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	    			__m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	    			__m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
	    			__m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );

	    			//
	    			// multiply left * right
	    			//
	    			__m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
	    			__m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );

	    			//
	    			// multiply with EV matrix (!?)
	    			//
	    			__m128d EV_t_l0_k0 = EVV[0];
	    			__m128d EV_t_l0_k2 = EVV[1];
	    			__m128d EV_t_l1_k0 = EVV[2];
	    			__m128d EV_t_l1_k2 = EVV[3];
	    			__m128d EV_t_l2_k0 = EVV[4];
	    			__m128d EV_t_l2_k2 = EVV[5];
	    			__m128d EV_t_l3_k0 = EVV[6];
	    			__m128d EV_t_l3_k2 = EVV[7];

	    			EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    			EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    			EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

	    			EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    			EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

	    			EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    			EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

	    			EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    			EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    			EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

	    			EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    			EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    			EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

	    			EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

	    			_mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
	    			_mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
	    		}
	    	}

	    	// END OF DNA with 4 categories
	    } else*/
	    for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];
			int state_right = (aln->at(ptn))[right->node->id];

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				double *left = partial_lh_left + (state_left*block+c*nstates);
				double *right = partial_lh_right + (state_right*block+c*nstates);
				for (x = 0; x < nstates; x+=VCSIZE) {
					vc_partial_lh_tmp[x/VCSIZE] = (VectorClass().load_a(&left[x]) * VectorClass().load_a(&right[x]));
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=2) {
					res1 = vc_partial_lh_tmp[0] * vc_inv_evec[i*nstates/VCSIZE];
					res2 = vc_partial_lh_tmp[0] * vc_inv_evec[(i+1)*nstates/VCSIZE];
					for (x = VCSIZE; x < nstates; x+=VCSIZE) {
						res1 += vc_partial_lh_tmp[x/VCSIZE] * vc_inv_evec[(i*nstates+x)/VCSIZE];
						res2 += vc_partial_lh_tmp[x/VCSIZE] * vc_inv_evec[((i+1)*nstates+x)/VCSIZE];
					}
					horizontal_add(res1,res2).store_a(&partial_lh[c*nstates+i]);
				}
			}
			partial_lh += block;
		}
		aligned_free(partial_lh_right);
		aligned_free(partial_lh_left);
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

		// pre compute information for left tip
		double *partial_lh_left = aligned_alloc_double((aln->STATE_UNKNOWN+1)*block);
		double *partial_lh_right = right->partial_lh;

		VectorClass vc_tip_lh[nstates/VCSIZE];
		VectorClass vc_lh_right[nstates/VCSIZE];
		Vec2d v2d_partial_lh_tmp[nstates/2];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			for (i = 0; i < nstates/VCSIZE; i++)
				vc_tip_lh[i].load_a(&tip_partial_lh[state*nstates+i*VCSIZE]);
			for (x = 0; x < block; x+=2) {
				vleft = VectorClass().load_a(&eleft[x*nstates]) * vc_tip_lh[0];
				vleft2 = VectorClass().load_a(&eleft[(x+1)*nstates]) * vc_tip_lh[0];
				for (i = VCSIZE; i < nstates; i+=VCSIZE) {
					vleft += VectorClass().load_a(&eleft[x*nstates+i]) * vc_tip_lh[i/VCSIZE];
					vleft2 += VectorClass().load_a(&eleft[(x+1)*nstates+i]) * vc_tip_lh[i/VCSIZE];
				}
				horizontal_add(vleft, vleft2).store_a(&partial_lh_left[state*block+x]);
			}
		}

		for (x = 0; x < block; x++) {
			size_t addr = aln->STATE_UNKNOWN * block;
			partial_lh_left[addr+x] = 1.0;
		}


		double sum_scale = 0.0;

		/*
		if (nstates == 4 && ncat == 0) {
	    	// Codes taken and adapted from PLL
	    	double *x2, *x3, *uX1;
	    	int state_left;
	    	size_t j, k;
	    	__m128d values[8], EVV[8];
	    	double maxima[2] __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)));
	    	  for(k = 0; k < 8; k++)
	    	    EVV[k] = _mm_load_pd(&inv_evec[k * 2]);

			for (i = 0; i < nptn; i++)
			{
				state_left = (aln->at(i))[left->node->id];
				__m128d maxv =_mm_setzero_pd();

				x2 = &partial_lh_right[i * 16];
				x3 = &partial_lh[i * 16];

				uX1 = &partial_lh_left[16 * state_left];

				for (j = 0; j < 4; j++)
				{

					//
					// multiply/add right side
					//
					double *x2_p = &x2[j*4];
					double *right_k0_p = &eright[j*16];
					double *right_k1_p = &eright[j*16 + 1*4];
					double *right_k2_p = &eright[j*16 + 2*4];
					double *right_k3_p = &eright[j*16 + 3*4];
					__m128d x2_0 = _mm_load_pd( &x2_p[0] );
					__m128d x2_2 = _mm_load_pd( &x2_p[2] );

					__m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
					__m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
					__m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
					__m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
					__m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
					__m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
					__m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
					__m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );



					right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
					right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

					right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
					right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

					right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
					right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
					right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);


					right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
					right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


					right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
					right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

					right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
					right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
					right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

					{
						//
						// load left side from tip vector
						//

						__m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
						__m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );


						//
						// multiply left * right
						//

						__m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
						__m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );


						//
						// multiply with EV matrix (!?)
						//

						__m128d EV_t_l0_k0 = EVV[0];
						__m128d EV_t_l0_k2 = EVV[1];
						__m128d EV_t_l1_k0 = EVV[2];
						__m128d EV_t_l1_k2 = EVV[3];
						__m128d EV_t_l2_k0 = EVV[4];
						__m128d EV_t_l2_k2 = EVV[5];
						__m128d EV_t_l3_k0 = EVV[6];
						__m128d EV_t_l3_k2 = EVV[7];


						EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
						EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
						EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

						EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
						EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

						EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
						EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

						EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
						EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
						EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

						EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
						EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
						EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

						EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

						values[j * 2]     = EV_t_l0_k0;
						values[j * 2 + 1] = EV_t_l2_k0;

						maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
						maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
					}
				}


				_mm_store_pd(maxima, maxv);

				if(PLL_MAX(maxima[0], maxima[1]) < PLL_MINLIKELIHOOD)
				{
					__m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

					_mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));
					_mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
					_mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
					_mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
					_mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));
					_mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
					_mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
					_mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));

					sum_scale += LOG_SCALING_THRESHOLD * (*aln)[i].frequency;
					dad_branch->scale_num[i] += 1;
					if (pattern_scale)
						pattern_scale[i] += LOG_SCALING_THRESHOLD;

				}
				else
				{
					_mm_store_pd(&x3[0], values[0]);
					_mm_store_pd(&x3[2], values[1]);
					_mm_store_pd(&x3[4], values[2]);
					_mm_store_pd(&x3[6], values[3]);
					_mm_store_pd(&x3[8], values[4]);
					_mm_store_pd(&x3[10], values[5]);
					_mm_store_pd(&x3[12], values[6]);
					_mm_store_pd(&x3[14], values[7]);
				}
			}

		} else*/

		for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (i = 0; i < nstates/VCSIZE; i++)
					vc_lh_right[i].load_a(&partial_lh_right[c*nstates+i*VCSIZE]);
				for (x = 0; x < nstates; x+=2) {
					size_t addr = c*nstatesqr+x*nstates;
					vright = VectorClass().load_a(&eright[addr]) * vc_lh_right[0];
					vright2 = VectorClass().load_a(&eright[addr+nstates]) * vc_lh_right[0];
					for (i = VCSIZE; i < nstates; i+=VCSIZE) {
						vright += VectorClass().load_a(&eright[addr+i]) * vc_lh_right[i/VCSIZE];
						vright2 += VectorClass().load_a(&eright[addr+i+nstates]) * vc_lh_right[i/VCSIZE];
					}
					v2d_partial_lh_tmp[x/2] = Vec2d().load(&partial_lh_left[state_left*block+c*nstates+x])
							* horizontal_add(vright, vright2);
				}
				// compute dot-product with inv_eigenvector
				// NOTE THAT it does not work with AVX yet
				for (i = 0; i < nstates; i+=2) {
					res1 = v2d_partial_lh_tmp[0] * v2d_inv_evec[i*nstates/2];
					res2 = v2d_partial_lh_tmp[0] * v2d_inv_evec[(i+1)*nstates/2];
					for (x = 2; x < nstates; x+=2) {
						res1 += v2d_partial_lh_tmp[x/2] * v2d_inv_evec[(i*nstates+x)/2];
						res2 += v2d_partial_lh_tmp[x/2] * v2d_inv_evec[((i+1)*nstates+x)/2];
					}
					horizontal_add(res1,res2).store_a(&partial_lh[c*nstates+i]);
				}
			}
            // check if one should scale partial likelihoods
			VectorClass vc_max = abs(VectorClass().load_a(partial_lh));
			for (i = VCSIZE; i < block; i+=VCSIZE) {
				vc_max = max(vc_max, abs(VectorClass().load_a(&partial_lh[i])));
			}
			double vmax = horizontal_max(vc_max);
            if (vmax < SCALING_THRESHOLD) {
            	// now do the likelihood scaling
            	VectorClass scale_thres(SCALING_THRESHOLD_INVER);
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

			partial_lh += block;
			partial_lh_right += block;
		}
		dad_branch->lh_scale_factor += sum_scale;
		aligned_free(partial_lh_left);

	} else {
		// both left and right are internal node
		double *partial_lh_left = left->partial_lh;
		double *partial_lh_right = right->partial_lh;

		double sum_scale = 0.0;
		Vec2d v2d_partial_lh_tmp[nstates/2];
		VectorClass vc_lh_left[nstates/VCSIZE], vc_lh_right[nstates/VCSIZE];

/*
		if (nstates == 4 && ncat == 0) {
	    	// Codes taken and adapted from PLL
			double  *x1, *x2, *x3;
			size_t j, k;
			__m128d values[8], EVV[8];
			double maxima[2] __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)));
			for(k = 0; k < 8; k++)
				EVV[k] = _mm_load_pd(&inv_evec[k * 2]);

			for (i = 0; i < nptn; i++)
			{
				__m128d maxv =_mm_setzero_pd();


				x1 = &partial_lh_left[i * 16];
				x2 = &partial_lh_right[i * 16];
				x3 = &partial_lh[i * 16];

				for (j = 0; j < 4; j++)
				{

					double *x1_p = &x1[j*4];
					double *left_k0_p = &eleft[j*16];
					double *left_k1_p = &eleft[j*16 + 1*4];
					double *left_k2_p = &eleft[j*16 + 2*4];
					double *left_k3_p = &eleft[j*16 + 3*4];

					__m128d x1_0 = _mm_load_pd( &x1_p[0] );
					__m128d x1_2 = _mm_load_pd( &x1_p[2] );

					__m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
					__m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
					__m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
					__m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
					__m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
					__m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
					__m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
					__m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );

					left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
					left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

					left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
					left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

					left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
					left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
					left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

					left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
					left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

					left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
					left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

					left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
					left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
					left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);


					//
					// multiply/add right side
					//
					double *x2_p = &x2[j*4];
					double *right_k0_p = &eright[j*16];
					double *right_k1_p = &eright[j*16 + 1*4];
					double *right_k2_p = &eright[j*16 + 2*4];
					double *right_k3_p = &eright[j*16 + 3*4];
					__m128d x2_0 = _mm_load_pd( &x2_p[0] );
					__m128d x2_2 = _mm_load_pd( &x2_p[2] );

					__m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
					__m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
					__m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
					__m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
					__m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
					__m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
					__m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
					__m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );

					right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
					right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

					right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
					right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

					right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
					right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
					right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

					right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
					right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


					right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
					right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

					right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
					right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
					right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

					//
					// multiply left * right
					//

					__m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
					__m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );


					//
					// multiply with EV matrix (!?)
					//

					__m128d EV_t_l0_k0 = EVV[0];
					__m128d EV_t_l0_k2 = EVV[1];
					__m128d EV_t_l1_k0 = EVV[2];
					__m128d EV_t_l1_k2 = EVV[3];
					__m128d EV_t_l2_k0 = EVV[4];
					__m128d EV_t_l2_k2 = EVV[5];
					__m128d EV_t_l3_k0 = EVV[6];
					__m128d EV_t_l3_k2 = EVV[7];


					EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
					EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
					EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

					EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
					EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

					EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
					EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

					EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
					EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
					EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );


					EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
					EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
					EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

					EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );


					values[j * 2] = EV_t_l0_k0;
					values[j * 2 + 1] = EV_t_l2_k0;

					maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
					maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
				}


				_mm_store_pd(maxima, maxv);

				if(PLL_MAX(maxima[0], maxima[1]) < PLL_MINLIKELIHOOD)
				{
					__m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

					_mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));
					_mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
					_mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
					_mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
					_mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));
					_mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
					_mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
					_mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));

					sum_scale += LOG_SCALING_THRESHOLD * (*aln)[i].frequency;
					dad_branch->scale_num[i] += 1;
					if (pattern_scale)
						pattern_scale[i] += LOG_SCALING_THRESHOLD;
				}
				else
				{
					_mm_store_pd(&x3[0], values[0]);
					_mm_store_pd(&x3[2], values[1]);
					_mm_store_pd(&x3[4], values[2]);
					_mm_store_pd(&x3[6], values[3]);
					_mm_store_pd(&x3[8], values[4]);
					_mm_store_pd(&x3[10], values[5]);
					_mm_store_pd(&x3[12], values[6]);
					_mm_store_pd(&x3[14], values[7]);
				}
			}
		} else*/
		for (ptn = 0; ptn < nptn; ptn++) {
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (i = 0; i < nstates/VCSIZE; i++) {
					vc_lh_left[i].load_a(&partial_lh_left[c*nstates+i*VCSIZE]);
					vc_lh_right[i].load_a(&partial_lh_right[c*nstates+i*VCSIZE]);
				}
				for (x = 0; x < nstates; x+=2) {
					size_t addr = c*nstatesqr+x*nstates;
					vleft = VectorClass().load_a(&eleft[addr]) * vc_lh_left[0];
					vleft2 = VectorClass().load_a(&eleft[addr+nstates]) * vc_lh_left[0];
					vright = VectorClass().load_a(&eright[addr]) * vc_lh_right[0];
					vright2 = VectorClass().load_a(&eright[addr+nstates]) * vc_lh_right[0];
					for (i = VCSIZE; i < nstates; i+=VCSIZE) {
						vleft += VectorClass().load_a(&eleft[addr+i]) * vc_lh_left[i/VCSIZE];
						vleft2 += VectorClass().load_a(&eleft[addr+i+nstates]) * vc_lh_left[i/VCSIZE];
						vright += VectorClass().load_a(&eright[addr+i]) * vc_lh_right[i/VCSIZE];
						vright2 += VectorClass().load_a(&eright[addr+i+nstates]) * vc_lh_right[i/VCSIZE];
					}
					v2d_partial_lh_tmp[x/2] = horizontal_add(vleft, vleft2) * horizontal_add(vright, vright2);
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=2) {
					res1 = v2d_partial_lh_tmp[0] * v2d_inv_evec[i*nstates/2];
					res2 = v2d_partial_lh_tmp[0] * v2d_inv_evec[(i+1)*nstates/2];
					for (x = 2; x < nstates; x+=2) {
						res1 += v2d_partial_lh_tmp[x/2] * v2d_inv_evec[(i*nstates+x)/2];
						res2 += v2d_partial_lh_tmp[x/2] * v2d_inv_evec[((i+1)*nstates+x)/2];
					}
					horizontal_add(res1,res2).store_a(&partial_lh[c*nstates+i]);
				}
			}

            // check if one should scale partial likelihoods
			VectorClass vc_max = abs(VectorClass().load_a(partial_lh));
			for (i = VCSIZE; i < block; i+=VCSIZE) {
				vc_max = max(vc_max, abs(VectorClass().load_a(&partial_lh[i])));
			}
			double vmax = horizontal_max(vc_max);
            if (vmax < SCALING_THRESHOLD) {
				// now do the likelihood scaling
            	VectorClass scale_thres(SCALING_THRESHOLD_INVER);
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
				dad_branch->scale_num[ptn] += 1;
				if (pattern_scale)
					pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }

			partial_lh += block;
			partial_lh_left += block;
			partial_lh_right += block;
		}
		dad_branch->lh_scale_factor += sum_scale;

	}

	aligned_free(eright);
	aligned_free(eleft);
//	aligned_free(partial_lh_tmp);
	aligned_free(inv_evec);
	aligned_free(evec);
}

template <const int nstates>
double PhyloTree::computeLikelihoodDervEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computePartialLikelihoodEigenTipSSE<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenTipSSE<nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    size_t ncat = site_rate->getNRate();
//    const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    //assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
//    const size_t block = ncat * nstates;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *val0 = aligned_alloc_double(block);
    double *val1 = aligned_alloc_double(block);
    double *val2 = aligned_alloc_double(block);
    double *state_freq = aligned_alloc_double(nstates);
    model->getStateFrequency(state_freq);

	for (c = 0; c < ncat; c++) {
		for (i = 0; i < nstates; i++) {
			double cof = eval[i]*site_rate->getRate(c);
			double val = exp(cof*dad_branch->length);
			double val1_ = cof*val;
			val0[c*nstates+i] = val;
			val1[c*nstates+i] = val1_;
			val2[c*nstates+i] = cof*val1_;
		}
	}
	VectorClass vc_val0[block/VCSIZE];
	VectorClass vc_val1[block/VCSIZE];
	VectorClass vc_val2[block/VCSIZE];
	for (i = 0; i < block/VCSIZE; i++) {
		vc_val0[i].load_a(&val0[i*VCSIZE]);
		vc_val1[i].load_a(&val1[i*VCSIZE]);
		vc_val2[i].load_a(&val2[i*VCSIZE]);
	}

	assert(theta_all);
	if (!theta_computed) {
		theta_computed = true;
		// precompute theta for fast branch length optimization

	    double *partial_lh_dad = dad_branch->partial_lh;
		double *theta = theta_all;

		if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
			for (ptn = 0; ptn < nptn; ptn++) {
				int state_dad = (aln->at(ptn))[dad->id];
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&tip_partial_lh[state_dad*nstates+i%nstates]) * VectorClass().load_a(&partial_lh_dad[i]))
							.store_a(&theta[i]);
				}
				partial_lh_dad += block;
				theta += block;
			}
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
	    	size_t all_entries = nptn*block;
			for (i = 0; i < all_entries; i+=VCSIZE) {
				(VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i]))
						.store_a(&theta[i]);
			}
	    }
		if (nptn % 2 != 0) {
			theta = &theta_all[nptn*block];
			VectorClass vc1(p_var_cat);
			for (i = 0; i < block/VCSIZE; i++) {
				(vc1 / vc_val0[i]).store_a(&theta[i*VCSIZE]);
			}
	//		memset(&theta[nptn*block], 0, block*sizeof(double));
		}
	}



    double *theta = theta_all;

	VectorClass vc_ptn, vc_df, vc_ddf, vc_theta, vc_ptn2, vc_df2, vc_ddf2, vc_theta2;
	Vec2d v2d_var_cat(p_var_cat, p_var_cat);
	Vec2d v2d_1(1.0, 1.0);
	Vec2d lh_final(0.0, 0.0), df_final(0.0, 0.0), ddf_final(0.0, 0.0);

	// allocate 1 more element
	double *ptn_freq = aligned_alloc_double(nptn+1);
	for (ptn = 0; ptn < nptn; ptn++)
		ptn_freq[ptn] = (*aln)[ptn].frequency;
	ptn_freq[nptn] = 0.0;

	assert(p_invar == 0.0);

	// perform 2 sites at the same time for SSE
	for (ptn = 0; ptn < nptn; ptn+=2) {
		// these stores values of 2 consecutive patterns
		Vec2d lh_ptn, df_ptn, ddf_ptn, inv_lh_ptn;

		vc_theta.load_a(theta);
		vc_theta2.load_a(theta+block);

		vc_ptn = vc_val0[0] * vc_theta;
		vc_df = vc_val1[0] * vc_theta;
		vc_ddf = vc_val2[0] * vc_theta;

		vc_ptn2 = vc_val0[0] * vc_theta2;
		vc_df2 = vc_val1[0] * vc_theta2;
		vc_ddf2 = vc_val2[0] * vc_theta2;

		for (i = 1; i < block/VCSIZE; i++) {
			vc_theta.load_a(&theta[i*VCSIZE]);
			vc_theta2.load_a(&theta[i*VCSIZE+block]);
			vc_ptn = mul_add(vc_theta, vc_val0[i], vc_ptn);
			vc_df = mul_add(vc_theta, vc_val1[i], vc_df);
			vc_ddf = mul_add(vc_theta, vc_val2[i], vc_ddf);
			vc_ptn2 = mul_add(vc_theta2, vc_val0[i], vc_ptn2);
			vc_df2 = mul_add(vc_theta2, vc_val1[i], vc_df2);
			vc_ddf2 = mul_add(vc_theta2, vc_val2[i], vc_ddf2);
		}
		theta += block*2;

		lh_ptn = horizontal_add(vc_ptn, vc_ptn2) * v2d_var_cat;

		// TODO: account for +I model
//		if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
//			lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
//		}
		inv_lh_ptn = v2d_1/lh_ptn;
		df_ptn = (horizontal_add(vc_df, vc_df2) * v2d_var_cat * inv_lh_ptn);
		ddf_ptn = (horizontal_add(vc_ddf, vc_ddf2) * v2d_var_cat * inv_lh_ptn);
		lh_ptn = log(lh_ptn);

		Vec2d freq;
		freq.load_a(&ptn_freq[ptn]);

		lh_final += lh_ptn * freq;
		df_final += df_ptn * freq;
		ddf_final += (ddf_ptn - df_ptn * df_ptn) * freq;

		lh_ptn.store_a(&_pattern_lh[ptn]);

    }
	tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor + horizontal_add(lh_final);
	df = horizontal_add(df_final);
	ddf = horizontal_add(ddf_final);

	aligned_free(ptn_freq);
	aligned_free(state_freq);
    aligned_free(val2);
    aligned_free(val1);
    aligned_free(val0);
    return tree_lh;
}


template <const int nstates>
double PhyloTree::computeLikelihoodBranchEigenTipSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
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
        computePartialLikelihoodEigenTipSSE<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenTipSSE<nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();
//    const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    //assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
//    const size_t block = ncat * nstates;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *partial_lh_dad = dad_branch->partial_lh;
    double *partial_lh_node = node_branch->partial_lh;
    double *val = aligned_alloc_double(block);
    double *state_freq = aligned_alloc_double(nstates);
    model->getStateFrequency(state_freq);

	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		for (i = 0; i < nstates; i++)
			val[c*nstates+i] = exp(eval[i]*len);
	}
    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
    	VectorClass vc_val, vc_tip_partial_lh, vc_partial_lh_dad;
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn;
			VectorClass vc_ptn = 0.0;
			int state_dad = (aln->at(ptn))[dad->id];
			for (c = 0; c < ncat; c++)
				for (i = 0; i < nstates; i+=VCSIZE) {
					vc_val.load_a(&val[c*nstates+i]);
					vc_tip_partial_lh.load_a(&tip_partial_lh[state_dad*nstates+i]);
					vc_partial_lh_dad.load_a(&partial_lh_dad[c*nstates+i]);
					vc_ptn += vc_val * vc_tip_partial_lh * vc_partial_lh_dad;
				}
			lh_ptn = horizontal_add(vc_ptn) * p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}

			assert(lh_ptn > 0.0);
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
			//partial_lh_node += nstates;
			partial_lh_dad += block;
		}
    } else {
    	// both dad and node are internal nodes
    	VectorClass vc_val, vc_partial_lh_node, vc_partial_lh_dad;

		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0;
			VectorClass vc_ptn = 0.0;
			for (i = 0; i < block; i+=VCSIZE) {
				vc_val.load_a(&val[i]);
				vc_partial_lh_node.load_a(&partial_lh_node[i]);
				vc_partial_lh_dad.load_a(&partial_lh_dad[i]);
				vc_ptn += vc_val * vc_partial_lh_node * vc_partial_lh_dad;
			}
			lh_ptn = horizontal_add(vc_ptn) * p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}

			assert(lh_ptn > 0.0);
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
			partial_lh_node += block;
			partial_lh_dad += block;
		}
    }
    if (pattern_lh)
        memmove(pattern_lh, _pattern_lh, aln->size() * sizeof(double));
    aligned_free(state_freq);
    aligned_free(val);
    return tree_lh;
}


/************************************************************************************************
 *
 *   SSE vectorized functions of the Naive implementation
 *
 *************************************************************************************************/

template<const int NSTATES>
inline double PhyloTree::computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
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
        //cout << "swapped\n";
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

    if (pattern_lh) {
        memmove(pattern_lh, _pattern_lh, orig_alnSize * sizeof(double));
    }
    delete[] trans_mat_orig;
    return tree_lh;
}

template<int NSTATES>
void PhyloTree::computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    //double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    //double *partial_lh_block;
    //bool do_scale = true;
    //double freq;
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
        //double *partial_lh_site;
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

//            } else {
//                // ambiguous character, for DNA, RNA
//                state = state - (NSTATES - 1);
//                for (int state2 = 0; state2 < NSTATES && state2 <= 6; state2++)
//                    if (state & (1 << state2)) {
//                        double *_par_lh_site = partial_lh_site + state2;
//                        for (cat = 0; cat < numCat; cat++) {
//                            *_par_lh_site = 1.0;
//                            _par_lh_site += NSTATES;
//                        }
//                    }
//            }
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
            computePartialLikelihoodSSE<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);
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
                    //ei_partial_lh_site.noalias() = (ei_partial_lh_child * ei_trans_state).cwiseProduct(ei_partial_lh_site);
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
                    if (pattern_scale)
                    pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
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
inline double PhyloTree::computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
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
    double *trans_mat_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
    double *trans_derv1_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
    double *trans_derv2_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
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
#pragma omp parallel for reduction(+: tree_lh, my_df, my_ddf,prob_const, prob_const_derv1, prob_const_derv2) \
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
			lh_ptn = log(lh_ptn);
			tree_lh += lh_ptn * freq;
			_pattern_lh[ptn] = lh_ptn;
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
    	prob_const = log(prob_const);
    	tree_lh -= nsites * prob_const;
    	for (int ptn = 0; ptn < orig_alnSize; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    }

    delete[] trans_derv2_orig;
    delete[] trans_derv1_orig;
    delete[] trans_mat_orig;
    df = my_df;
    ddf = my_ddf;
    return tree_lh;
}

template<int NSTATES>
void PhyloTree::computeThetaSSE(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    cout << "This has not been implemented yet" << endl;
    exit(1);
}

void PhyloTree::computeTheta(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeThetaSSE<2>(dad_branch, dad);
        case 4:
            return computeThetaSSE<4>(dad_branch, dad);
        case 20:
            return computeThetaSSE<20>(dad_branch, dad);
        default:
            computeThetaNaive(dad_branch, dad);
            break;
        }
    } else {
        computeThetaNaive(dad_branch, dad);
    }
}

void PhyloTree::computeThetaNaive(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    //cout << "Computing theta vector " << endl;
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
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
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);

    double *partial_lh_site = node_branch->partial_lh;
    double *partial_lh_child = dad_branch->partial_lh;
    double *theta_ptn = theta_all;
    int alnSize = aln->getNPattern();
    // BQM's question: why state_freq has size = alnSize?
    double* state_freq = new double[alnSize];
    model->getStateFrequency(state_freq);
    int num_cat = site_rate->getNRate();
    GTRModel* gtr_model = reinterpret_cast<GTRModel *>(model);
    //double* eigen_coff = gtr_model->getEigenCoeff();
    double** inv_eigen_vector = gtr_model->getInverseEigenvectors();
	double* partial_lh_child_ptn = partial_lh_child;
	double* partial_lh_site_ptn = partial_lh_site;

	int numStates = aln->num_states;
	for (int ptn = 0; ptn < alnSize; ++ptn) {
	    for (int i = 0; i < numStates; ++i) {
	    	partial_lh_child_ptn = partial_lh_child + ptn * numStates * num_cat;
	    	partial_lh_site_ptn = partial_lh_site + ptn * numStates * num_cat;
		    for (int c = 0; c < num_cat; ++c) {
				double term1 = 0;
				double term2 = 0;
				for (int x = 0; x < numStates; ++x) {
					// Compute Sigma_x pi_x u_xi L^h_a(x,c)
					term1 += inv_eigen_vector[i][x] * partial_lh_site_ptn[x];
					term2 += inv_eigen_vector[i][x] * partial_lh_child_ptn[x];
				}
	            partial_lh_child_ptn += numStates;
	            partial_lh_site_ptn += numStates;
				*theta_ptn = term1 * term2;
				//cout << "theta_ptn : " << *theta_ptn << endl;
				theta_ptn++;
			}
		}
	}
	delete[] state_freq;
}


void PhyloTree::initiateMyEigenCoeff() {
    assert(model);
    GTRModel* gtr_model = reinterpret_cast<GTRModel *>(model);
    double* eigen_coff = gtr_model->getEigenCoeff();
    int numStates = aln->num_states;
    if (!myEigenCoeff) {
        myEigenCoeff = new double[numStates * numStates * numStates];
    }
    int i = 0;
    for (int j = 0; j < numStates; ++j)
        for (int xa = 0; xa < numStates; ++xa)
            for (int xb = 0; xb < numStates; ++xb) {
                myEigenCoeff[i] = eigen_coff[xa * numStates * numStates + xb * numStates + j];
                ++i;
            }
}


void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
	switch(aln->num_states) {
	case 4:
		switch(sse) {
		case LK_SSE: computePartialLikelihoodSSE<4>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN: computePartialLikelihoodEigen<4>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN_TIP_SSE: computePartialLikelihoodEigenTipSSE<4>(dad_branch, dad, pattern_scale); break;
		case LK_NORMAL: computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
		}
		break;
	case 20:
		switch(sse) {
		case LK_SSE: computePartialLikelihoodSSE<20>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN: computePartialLikelihoodEigen<20>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN_TIP_SSE: computePartialLikelihoodEigenTipSSE<20>(dad_branch, dad, pattern_scale); break;
		case LK_NORMAL: computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: computePartialLikelihoodSSE<2>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN: computePartialLikelihoodEigen<2>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN_TIP_SSE: computePartialLikelihoodEigenTipSSE<2>(dad_branch, dad, pattern_scale); break;
		case LK_NORMAL: computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
		}
		break;

	default:
		computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
	}
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
	switch(aln->num_states) {
	case 4:
		switch(sse) {
		case LK_SSE: return computeLikelihoodBranchSSE<4>(dad_branch, dad, pattern_lh);
		case LK_EIGEN: return computeLikelihoodBranchEigen<4>(dad_branch, dad, pattern_lh);
		case LK_EIGEN_TIP_SSE: return computeLikelihoodBranchEigenTipSSE<4>(dad_branch, dad, pattern_lh);
		case LK_NORMAL: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
		}
		break;
	case 20:
		switch(sse) {
		case LK_SSE: return computeLikelihoodBranchSSE<20>(dad_branch, dad, pattern_lh);
		case LK_EIGEN: return computeLikelihoodBranchEigen<20>(dad_branch, dad, pattern_lh);
		case LK_EIGEN_TIP_SSE: return computeLikelihoodBranchEigenTipSSE<20>(dad_branch, dad, pattern_lh);
		case LK_NORMAL: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: return computeLikelihoodBranchSSE<2>(dad_branch, dad, pattern_lh);
		case LK_EIGEN: return computeLikelihoodBranchEigen<2>(dad_branch, dad, pattern_lh);
		case LK_EIGEN_TIP_SSE: return computeLikelihoodBranchEigenTipSSE<2>(dad_branch, dad, pattern_lh);
		case LK_NORMAL: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
		}
		break;

	default:
		return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
	}
	return 0.0;
}

/*
 * This function is called millions times. So it is not a good idea to
 * have a if and switch here.
 */
double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    switch (aln->num_states) {
    case 4:
    	switch(sse) {
    	case LK_SSE: return computeLikelihoodDervSSE<4>(dad_branch, dad, df, ddf);
    	case LK_EIGEN: return computeLikelihoodDervEigen<4>(dad_branch, dad, df, ddf);
    	case LK_EIGEN_TIP_SSE: return computeLikelihoodDervEigenTipSSE<4>(dad_branch, dad, df, ddf);
    	case LK_NORMAL: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
    	}
    	break;
	case 20:
		switch(sse) {
		case LK_SSE: return computeLikelihoodDervSSE<20>(dad_branch, dad, df, ddf);
		case LK_EIGEN: return computeLikelihoodDervEigen<20>(dad_branch, dad, df, ddf);
		case LK_EIGEN_TIP_SSE: return computeLikelihoodDervEigenTipSSE<20>(dad_branch, dad, df, ddf);
		case LK_NORMAL: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: return computeLikelihoodDervSSE<2>(dad_branch, dad, df, ddf);
		case LK_EIGEN: return computeLikelihoodDervEigen<2>(dad_branch, dad, df, ddf);
		case LK_EIGEN_TIP_SSE: return computeLikelihoodDervEigenTipSSE<2>(dad_branch, dad, df, ddf);
		case LK_NORMAL: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
		}
		break;
	default:
		return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);

    }
    return 0.0;
}

double PhyloTree::computeLikelihoodDervFast(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeLikelihoodDervFastSSE<2>(dad_branch, dad, df, ddf);
        case 4:
            return computeLikelihoodDervFastSSE<4>(dad_branch, dad, df, ddf);
        case 20:
            return computeLikelihoodDervFastSSE<20>(dad_branch, dad, df, ddf);
        default:
            return computeLikelihoodDervFastNaive(dad_branch, dad, df, ddf);
            //cout << "Bad number of states: " << aln->num_states << endl;
            //exit(1);
        }
    } else {
        return computeLikelihoodDervFastNaive(dad_branch, dad, df, ddf);
    }
}

double PhyloTree::computeLikelihoodDervFastNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if dad is a leaf
    // NEW: swap if root_state is given
    if (node->isLeaf() || (node->name == ROOT_NAME && root_state != aln->STATE_UNKNOWN)) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double lh_ptn = 0.0;
    double lh_ptn_derv1 = 0.0;
    double lh_ptn_derv2 = 0.0;
    double my_df = 0.0;
    double my_ddf = 0.0;
    double p_invar = site_rate->getPInvar();
    int ncat = site_rate->getNRate();
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    double derv1_frac;
    double derv2_frac;
    int nstates = aln->num_states;
    double* state_freq = new double[nstates];
    model->getStateFrequency(state_freq);
    int discrete_cat = site_rate->getNDiscreteRate();
    double *theta_ptn = theta_all;
    size_t nptn = aln->getNPattern();
    double t = dad_branch->length;

    double* rates = new double[discrete_cat];
    for (int i = 0; i < discrete_cat; i++) {
    	rates[i] = site_rate->getRate(i);
    }
    //double* rates = site_rate->getRates();
    double* lambda = dynamic_cast<GTRModel*>(model)->getEigenvalues();
    int block = discrete_cat * nstates;
    // array containing lambda_i * r_c
    double* lambda_r = new double[ block ];
    // array containing exp(lambda_i * r_c * t)
    double* exp_part = new double[ block ];
    // array containing square of lambda_i * r_c
    double* lambda_r_sqr = new double[ block ];

    // now initialize all the arrays (pre-computation before coming to the big loop)
    for (int i = 0; i < nstates; i++) {
    	for (int c = 0; c < discrete_cat; c++) {
    		lambda_r[ i * discrete_cat + c ] = lambda[i] * rates[c];
    		lambda_r_sqr[ i * discrete_cat + c ] = lambda_r[ i * discrete_cat + c ] * lambda_r[ i * discrete_cat + c ];
    		exp_part [ i * discrete_cat + c ] = exp( lambda_r[ i * discrete_cat + c ] * t);
    	}
    }

    int pointer_jump = nstates * discrete_cat;
    for (int ptn = 0; ptn < nptn; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        for (int i = 0; i < nstates; i++) {
        	for (int c = 0; c < discrete_cat; c++) {
        		double base = exp_part[ i * discrete_cat + c ] * theta_ptn[ i * discrete_cat + c ];
        		lh_ptn += base;
        		lh_ptn_derv1 +=  lambda_r[ i * discrete_cat + c ] * base;
        		lh_ptn_derv2 +=  lambda_r_sqr[ i * discrete_cat + c ] * base;
        	}
        }
        theta_ptn += pointer_jump;
        lh_ptn = lh_ptn * p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
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
        double freq = (*aln)[ptn].frequency;
        double tmp1 = derv1_frac * freq;
        double tmp2 = derv2_frac * freq;
        my_df += tmp1;
        my_ddf += tmp2 - tmp1 * derv1_frac;
        lh_ptn = log(lh_ptn);
        //cout << lh_ptn << endl;
        //cout << "lh_ptn = " << lh_ptn << endl;
        tree_lh += lh_ptn * aln->at(ptn).frequency;
        //cout << tree_lh << endl;
        _pattern_lh[ptn] = lh_ptn;
    }
    delete [] lambda_r;
    delete [] lambda_r_sqr;
    delete [] exp_part;
    delete [] rates;
    delete [] state_freq;
    df = my_df;
    ddf = my_ddf;
    return tree_lh;
}

template<int NSTATES>
double PhyloTree::computeLikelihoodDervFastSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if node is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double lh_ptn = 0.0;
    double lh_ptn_derv1 = 0.0;
    double lh_ptn_derv2 = 0.0;
    double p_invar = site_rate->getPInvar();
    int numCat = site_rate->getNRate();
    double p_var_cat = (1.0 - p_invar) / (double) numCat;
    double derv1_frac;
    double derv2_frac;
    double state_freq[NSTATES];
    int discrete_cat = site_rate->getNDiscreteRate();
    model->getStateFrequency(state_freq);
    double* rates = new double[discrete_cat];
    for (int i = 0; i < numCat; i++) {
    	rates[i] = site_rate->getRate(i);
    }
    Map<Matrix<double, Dynamic, 1>, Aligned> ei_rates(rates, discrete_cat);
    Map<Matrix<double, 1, NSTATES>, Aligned> ei_eigenvalues(dynamic_cast<GTRModel*>(model)->getEigenvalues());
    Matrix(NSTATES)ei_rates_times_lambdas = ei_rates * ei_eigenvalues;
    Matrix(NSTATES)expo_time =
    (ei_rates_times_lambdas * dad_branch->length).array().exp();
    Matrix(NSTATES)expo_time_derv1 = expo_time.cwiseProduct(
            ei_rates_times_lambdas);
    Matrix(NSTATES)expo_time_derv2 = expo_time_derv1.cwiseProduct(
            ei_rates_times_lambdas);
    double *theta_ptn = theta_all;
    int num_patterns = aln->getNPattern();
    for (int ptn = 0; ptn < num_patterns; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        /*
         for (int k = 0; k < discrete_cat; k++) {
         for (int i = 0; i < num_states; i++) {
         double x =theta_ptn[k*num_states+i] *
         exp(site_rate->getRates()[k]*dad_branch->length*dynamic_cast<GTRModel*> (model)->getEigenvalues()[i]);
         //cout << exp(site_rate->getRates()[k]*dad_branch->length*dynamic_cast<GTRModel*> (model)->getEigenvalues()[i]) << " ";
         lh_ptn += x;
         //cout << x << " ";
         //cout << theta_ptn[k*num_states+i] << " ";
         double y = x * site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i];
         //cout << y << " ";
         //cout << site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i] << " ";
         lh_ptn_derv1 += y;
         lh_ptn_derv2 += y*site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i];

         }
         }
         */

        Map<Matrix<double, Dynamic, Dynamic, RowMajor>, Aligned> ei_theta_ptn(theta_ptn, discrete_cat, NSTATES);
        //Map<Matrix<double, NSTATES, NSTATES, RowMajor>, Aligned> ei_theta_ptn(theta_ptn);
        //cout << "ei_theta_ptn" << endl;
        //cout << ei_theta_ptn << endl;
        //ArrayXXd ei_theta_times_expo = ei_theta_ptn * expo_time;

        lh_ptn = ei_theta_ptn.cwiseProduct(expo_time).sum();
        lh_ptn_derv1 = ei_theta_ptn.cwiseProduct(expo_time_derv1).sum();
        lh_ptn_derv2 = ei_theta_ptn.cwiseProduct(expo_time_derv2).sum();

        /*
         for (int cat = 0; cat < discrete_cat; ++cat) {
         MappedRowVec(NSTATES) ei_theta_ptn_cat(theta_ptn);
         lh_ptn += ei_theta_ptn_cat.dot(expo_time.row(cat));
         lh_ptn_derv1 += ei_theta_ptn_cat.dot(expo_time_derv1.row(cat));
         lh_ptn_derv2 += ei_theta_ptn_cat.dot(expo_time_derv2.row(cat));
         theta_ptn += NSTATES;
         }
         */

        /*
         ArrayXXd ei_theta_times_expo_derv1 = ei_rates_times_lambdas * ei_theta_times_expo;
         ArrayXXd ei_theta_times_expo_derv1 = ei_theta_ptn.cwiseProduct(expo_time_derv1);

         ArrayXXd ei_theta_times_expo_derv2 = ei_theta_ptn.cwiseProduct(expo_time_derv2);
         lh_ptn = ei_theta_times_expo.sum();
         lh_ptn_derv1 = ei_theta_times_expo_derv1.sum();
         lh_ptn_derv2 = (ei_theta_times_expo_derv1 * ei_rates_times_lambdas).sum();
         */
        theta_ptn += NSTATES * discrete_cat;

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
        double tmp1 = derv1_frac * aln->at(ptn).frequency;
        double tmp2 = derv2_frac * aln->at(ptn).frequency;
        df += tmp1;
        ddf += tmp2 - tmp1 * derv1_frac;
        lh_ptn = log(lh_ptn);
        tree_lh += lh_ptn * aln->at(ptn).frequency;
        _pattern_lh[ptn] = lh_ptn;
    }
    return tree_lh;

}
