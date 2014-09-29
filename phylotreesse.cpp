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

void PhyloTree::computeTipPartialLikelihood() {
	if (tip_partial_lh_computed)
		return;
	tip_partial_lh_computed = true;
	int i, x, state, nstates = aln->num_states;

	double *evec = new double[nstates*nstates];
    double *inv_evec = new double[nstates*nstates];
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
	delete [] inv_evec;
	delete [] evec;
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
    size_t nstatesqr=nstates*nstates, i, x;
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

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = new double[(aln->STATE_UNKNOWN+1)*block];
		double *partial_lh_right = new double[(aln->STATE_UNKNOWN+1)*block];

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			for (x = 0; x < block; x++) {
				double vleft1 = 0.0, vleft2 = 0.0;
				for (i = 0; i < nstates; i+=2) {
					vleft1 += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
					vleft2 += eleft[x*nstates+i+1] * tip_partial_lh[state*nstates+i+1];
				}
				partial_lh_left[state*block+x] = vleft1+vleft2;
			}
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
			for (x = 0; x < block; x++) {
				double vright1 = 0.0, vright2 = 0.0;
				for (i = 0; i < nstates; i+=2) {
					vright1 += eright[x*nstates+i] * tip_partial_lh[state*nstates+i];
					vright2 += eright[x*nstates+i+1] * tip_partial_lh[state*nstates+i+1];
				}
				partial_lh_right[state*block+x] = vright1+vright2;
			}
		}

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];
			int state_right = (aln->at(ptn))[right->node->id];
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				double *left = partial_lh_left + (state_left*block+c*nstates);
				double *right = partial_lh_right + (state_right*block+c*nstates);
				for (x = 0; x < nstates; x++) {
					partial_lh_tmp[x] = left[x] * right[x];
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res1 = 0.0, res2=0.0;
					for (x = 0; x < nstates; x+=2) {
						res1 += partial_lh_tmp[x]*inv_evec[i*nstates+x];
						res2 += partial_lh_tmp[x+1]*inv_evec[i*nstates+x+1];
					}
					partial_lh[c*nstates+i] = res1+res2;
				}
			}
			partial_lh += block;
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
			for (x = 0; x < block; x++) {
				double vleft1 = 0.0, vleft2 = 0.0;
				for (i = 0; i < nstates; i+=2) {
					vleft1 += eleft[x*nstates+i] * tip_partial_lh[state*nstates+i];
					vleft2 += eleft[x*nstates+i+1] * tip_partial_lh[state*nstates+i+1];
				}
				partial_lh_left[state*block+x] = vleft1+vleft2;
			}
		}

		double sum_scale = 0.0;
		for (ptn = 0; ptn < nptn; ptn++) {
			int state_left = (aln->at(ptn))[left->node->id];

			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright1 = 0.0, vright2 = 0.0;
					size_t addr = c*nstatesqr+x*nstates;
					vleft = partial_lh_left[state_left*block+c*nstates+x];
					for (i = 0; i < nstates; i+=2) {
						vright1 += eright[addr+i] * partial_lh_right[c*nstates+i];
						vright2 += eright[addr+i+1] * partial_lh_right[c*nstates+i+1];
					}
					partial_lh_tmp[x] = vleft * (vright1+vright2);
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res1 = 0.0, res2=0.0;
					for (x = 0; x < nstates; x+=2) {
						res1 += partial_lh_tmp[x]*inv_evec[i*nstates+x];
						res2 += partial_lh_tmp[x+1]*inv_evec[i*nstates+x+1];
					}
					partial_lh[c*nstates+i] = res1+res2;
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
			//partial_lh_left += nstates;
			partial_lh_right += block;
		}
		dad_branch->lh_scale_factor += sum_scale;
		delete [] partial_lh_left;

	} else {
		// both left and right are internal node
		double *partial_lh_left = left->partial_lh;
		double *partial_lh_right = right->partial_lh;

		double sum_scale = 0.0;
//		memset(partial_lh, 0, block*nptn*sizeof(double));
		for (ptn = 0; ptn < nptn; ptn++) {
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
//					prod = vleft * vright;
//					for (i = 0; i < nstates; i++)
//						partial_lh[c*nstates+i] += prod * inv_evec[i*nstates+x];
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					double res1 = 0.0, res2=0.0;
					for (x = 0; x < nstates; x+=2) {
						res1 += partial_lh_tmp[x]*inv_evec[i*nstates+x];
						res2 += partial_lh_tmp[x+1]*inv_evec[i*nstates+x+1];
					}
					partial_lh[c*nstates+i] = res1+res2;
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
    //const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
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
				for (i = 0; i < block; i++) {
					//theta[i] = partial_lh_node[i%nstates] * partial_lh_dad[i];
					theta[i] = tip_partial_lh[state_dad*nstates+i%nstates] * partial_lh_dad[i];
				}
				//partial_lh_node += nstates;
				partial_lh_dad += block;
				theta += block;
			}
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
	    	size_t all_entries = nptn*block;
			for (i = 0; i < all_entries; i++) {
				theta[i] = partial_lh_node[i] * partial_lh_dad[i];
			}
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
		for (i = 0; i < block; i++) {
			lh_ptn += val0[i] * theta[i];
			df_ptn += val1[i] * theta[i];
			ddf_ptn += val2[i] * theta[i];
		}
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
		//assert(lh_ptn > 0.0);
		lh_ptn = log(lh_ptn);
		tree_lh += lh_ptn * freq;
		_pattern_lh[ptn] = lh_ptn;
		theta += block;
    }
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
    //const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
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
			for (c = 0; c < ncat; c++)
				for (i = 0; i < nstates; i++) {
//					lh_ptn +=  val[c*nstates+i] * partial_lh_node[i] * partial_lh_dad[c*nstates+i];
					lh_ptn +=  val[c*nstates+i] * tip_partial_lh[state_dad*nstates+i] * partial_lh_dad[c*nstates+i];
				}
			lh_ptn *= p_var_cat;
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
			//partial_lh_node += nstates;
			partial_lh_dad += block;
		}
    } else {
    	// both dad and node are internal nodes
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0;
			for (i = 0; i < block; i++)
				lh_ptn +=  val[i] * partial_lh_node[i] * partial_lh_dad[i];
			lh_ptn *= p_var_cat;
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
    delete [] val;
    return tree_lh;
}

#include "vectorclass/vectorclass.h"

/************************************************************************************************
 *
 *   SSE vectorized versions of above functions
 *
 *************************************************************************************************/

#include "pll/mem_alloc.h"

#ifdef __AVX
#define VectorClass Vec4d
#pragma message "Using AVX instructions"
#else
#define VectorClass Vec2d
//#pragma message "Using SS3 instructions"
#endif

double *aligned_alloc_double(size_t size) {
	void *res;
	rax_posix_memalign(&res, MEM_ALIGNMENT, size*sizeof(double));
	return (double*)res;
}

void aligned_free(void *mem) {
	free(mem);
}

template<const int nstates>
void PhyloTree::computePartialLikelihoodEigenSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    size_t ptn, c;
    size_t nptn = aln->size();
    size_t ncat = site_rate->getNRate();
    //const size_t ncat = 4;

    //size_t nstates = aln->num_states ;
    //const size_t nstates = 4;
    size_t nstatesqr = nstates*nstates, i, x;
    size_t block = nstates * ncat;

    PhyloNode *node = (PhyloNode*)(dad_branch->node);
    double *partial_lh = dad_branch->partial_lh;

    double *evec = aligned_alloc_double(nstates*nstates);
    double *inv_evec = aligned_alloc_double(nstates*nstates);
	double **_evec = model->getEigenvectors(), **_inv_evec = model->getInverseEigenvectors();
	assert(_inv_evec && _evec);
	for (i = 0; i < nstates; i++) {
		memcpy(evec+i*nstates, _evec[i], nstates*sizeof(double));
		memcpy(inv_evec+i*nstates, _inv_evec[i], nstates*sizeof(double));
	}
	double *eval = model->getEigenvalues();

    dad_branch->lh_scale_factor = 0.0;

	if (node->isLeaf()) {
		// external node
		assert(node->id < aln->getNSeq());
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		for (ptn = 0; ptn < nptn; ptn++) {
			int state = (aln->at(ptn))[node->id];
			if (state < nstates) {
				// simple state
				for (i = 0; i < nstates; i++)
					partial_lh[i] = inv_evec[i*nstates+state];
			} else if (state == aln->STATE_UNKNOWN) {
				// gap or unknown state
				//dad_branch->scale_num[ptn] = -1;
				memset(partial_lh, 0, nstates*sizeof(double));
				for (i = 0; i < nstates; i++) {
					for (x = 0; x < nstates; x++) {
						partial_lh[i] += inv_evec[i*nstates+x];
					}
				}
			} else {
				// ambiguous state
				memset(partial_lh, 0, nstates*sizeof(double));
				state -= (nstates-1);
				for (i = 0; i < nstates; i++) {
					for (x = 0; x < nstates; x++)
						if (state & (1 << x))
							partial_lh[i] += inv_evec[i*nstates+x];
				}

			}
			partial_lh += nstates;
		} // for loop over ptn
		aligned_free(inv_evec);
		aligned_free(evec);
		return;
	}

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
		computePartialLikelihoodEigenSSE<nstates>(left, node, pattern_scale);
	if ((right->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigenSSE<nstates>(right, node, pattern_scale);
	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;
	double *partial_lh_left = left->partial_lh, *partial_lh_right = right->partial_lh;
	double *partial_lh_tmp = aligned_alloc_double(nstates);
	double *eleft = aligned_alloc_double(block*nstates);
	double *eright = aligned_alloc_double(block*nstates);
	double *expleft = aligned_alloc_double(nstates);
	double *expright = aligned_alloc_double(nstates);

	for (c = 0; c < ncat; c++) {
		VectorClass vc_a, vc_b, vc_c;
		VectorClass vc_left;
		VectorClass vc_right;

		double len_left = site_rate->getRate(c) * left->length;
		double len_right = site_rate->getRate(c) * right->length;
		for (i = 0; i < nstates; i++) {
			expleft[i] = exp(eval[i]*len_left);
			expright[i] = exp(eval[i]*len_right);
		}
		for (x = 0; x < nstates; x++) {
			size_t addr = c*nstatesqr+x*nstates;
			for (i = 0; i < nstates; i+=VectorClass::size()) {
				vc_a.load_a(evec+x*nstates+i);
				vc_b.load_a(expleft+i);
				vc_c.load_a(expright+i);
				vc_left = vc_a * vc_b;
				vc_right = vc_a * vc_c;
				vc_left.store_a(eleft+addr+i);
				vc_right.store_a(eright+addr+i);
			}
			/*
			for (i = 0; i < nstates; i++) {
				eleft[addr+i] = evec[x*nstates+i] * expleft[i];
				eright[addr+i] = evec[x*nstates+i] * expright[i];
			}*/
		}
	}
	aligned_free(expright);
	aligned_free(expleft);

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case
		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		VectorClass vc_left, vc_right, vc_res;
		for (ptn = 0; ptn < nptn; ptn++) {
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					size_t addr = c*nstatesqr+x*nstates;
					vc_left=VectorClass().load_a(eleft+addr) * VectorClass().load_a(partial_lh_left);
					vc_right = VectorClass().load_a(eright+addr) * VectorClass().load_a(partial_lh_right);
					for (i = VectorClass::size(); i < nstates; i+=VectorClass::size()) {
						vc_left += VectorClass().load_a(eleft+addr+i) * VectorClass().load_a(partial_lh_left+i);
						vc_right += VectorClass().load_a(eright+addr+i) * VectorClass().load_a(partial_lh_right+i);
					}
					partial_lh_tmp[x] = horizontal_add(vc_left) * horizontal_add(vc_right);
				}
				for (i = 0; i < nstates; i++) {
					vc_res = VectorClass().load_a(partial_lh_tmp) * VectorClass().load_a(inv_evec+i*nstates);
					for (x = VectorClass::size(); x < nstates; x+=VectorClass::size()) {
						vc_res += VectorClass().load_a(partial_lh_tmp+x) * VectorClass().load_a(inv_evec+i*nstates+x);
					}
					partial_lh[c*nstates+i] = horizontal_add(vc_res);
				}
			}
			partial_lh += block;
			partial_lh_left += nstates;
			partial_lh_right += nstates;
		}
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));
		double sum_scale = 0.0;
		VectorClass vc_left, vc_right, vc_res;
		for (ptn = 0; ptn < nptn; ptn++) {
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					size_t addr = c*nstatesqr+x*nstates;
					vc_left = VectorClass().load_a(eleft+addr) * VectorClass().load_a(partial_lh_left);
					vc_right = VectorClass().load_a(eright+addr) * VectorClass().load_a(partial_lh_right+c*nstates);
					for (i = VectorClass::size(); i < nstates; i+=VectorClass::size()) {
						vc_left += VectorClass().load_a(eleft+addr+i) * VectorClass().load_a(partial_lh_left+i);
						vc_right += VectorClass().load_a(eright+addr+i) * VectorClass().load_a(partial_lh_right+c*nstates+i);
					}
					partial_lh_tmp[x] = horizontal_add(vc_left) * horizontal_add(vc_right);

				}
				for (i = 0; i < nstates; i++) {
					vc_res = VectorClass().load_a(partial_lh_tmp) * VectorClass().load_a(inv_evec+i*nstates);
					for (x = VectorClass::size(); x < nstates; x+=VectorClass::size()) {
						vc_res += VectorClass().load_a(partial_lh_tmp+x) * VectorClass().load_a(inv_evec+i*nstates+x);
					}
					partial_lh[c*nstates+i] = horizontal_add(vc_res);
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
			partial_lh_left += nstates;
			partial_lh_right += block;
		}
		dad_branch->lh_scale_factor += sum_scale;

	} else {
		// both left and right are internal node
		double sum_scale = 0.0;
		VectorClass vc_left, vc_right, vc_res;
		for (ptn = 0; ptn < nptn; ptn++) {
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					size_t addr = c*nstatesqr+x*nstates;
					vc_left = VectorClass().load_a(eleft+addr) * VectorClass().load_a(partial_lh_left+c*nstates);
					vc_right = VectorClass().load_a(eright+addr) * VectorClass().load_a(partial_lh_right+c*nstates);
					for (i = VectorClass::size(); i < nstates; i+=VectorClass::size()) {
						vc_left += VectorClass().load_a(eleft+addr+i) * VectorClass().load_a(partial_lh_left+c*nstates+i);
						vc_right += VectorClass().load_a(eright+addr+i) * VectorClass().load_a(partial_lh_right+c*nstates+i);
					}
					partial_lh_tmp[x] = horizontal_add(vc_left) * horizontal_add(vc_right);
				}
				for (i = 0; i < nstates; i++) {
					vc_res = VectorClass().load_a(partial_lh_tmp) * VectorClass().load_a(inv_evec+i*nstates);
					for (x = VectorClass::size(); x < nstates; x+=VectorClass::size()) {
						vc_res += VectorClass().load_a(partial_lh_tmp+x) * VectorClass().load_a(inv_evec+i*nstates+x);
					}
					partial_lh[c*nstates+i] = horizontal_add(vc_res);
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
		}
		dad_branch->lh_scale_factor += sum_scale;

	}

	aligned_free(eright);
	aligned_free(eleft);
	aligned_free(partial_lh_tmp);
	aligned_free(inv_evec);
	aligned_free(evec);
}

template<const int nstates>
double PhyloTree::computeLikelihoodDervEigenSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computePartialLikelihoodEigenSSE<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSSE<nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    //size_t ncat = site_rate->getNRate();
    const size_t ncat = 4;

    //double p_invar = site_rate->getPInvar();
    //assert(p_invar == 0.0); // +I model not supported yet
    //double p_var_cat = (1.0 - p_invar) / (double) ncat;
    const double p_var_cat = 0.25;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
    const size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *partial_lh_dad = dad_branch->partial_lh;
    double *partial_lh_node = node_branch->partial_lh;
    double *val0 = aligned_alloc_double(block);
    double *val1 = aligned_alloc_double(block);
    double *val2 = aligned_alloc_double(block);
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
	VectorClass vc_res;
    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn, df_ptn, ddf_ptn;
			VectorClass vc_ptn = 0.0;
			VectorClass vc_df_ptn = 0.0;
			VectorClass vc_ddf_ptn = 0.0;
			for (i = 0; i < block; i+=VectorClass::size()) {
				vc_res = VectorClass().load_a(partial_lh_node+i%nstates) * VectorClass().load_a(partial_lh_dad+i);
				vc_ptn += vc_res * VectorClass().load_a(val0+i);
				vc_df_ptn += vc_res * VectorClass().load_a(val1+i);
				vc_ddf_ptn += vc_res * VectorClass().load_a(val2+i);
			}
			lh_ptn = horizontal_add(vc_ptn);
			df_ptn = horizontal_add(vc_df_ptn);
			ddf_ptn = horizontal_add(vc_ddf_ptn);
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
//			assert(lh_ptn > 0.0);
	        lh_ptn = log(lh_ptn);
	        tree_lh += lh_ptn * freq;
	        _pattern_lh[ptn] = lh_ptn;
			partial_lh_node += nstates;
			partial_lh_dad += block;
		}
    } else {
    	// both dad and node are internal nodes
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn, df_ptn, ddf_ptn;
			VectorClass vc_ptn = 0.0;
			VectorClass vc_df_ptn = 0.0;
			VectorClass vc_ddf_ptn = 0.0;
			for (i = 0; i < block; i+=VectorClass::size()) {
				vc_res = VectorClass().load_a(partial_lh_node+i) * VectorClass().load_a(partial_lh_dad+i);
				vc_ptn += vc_res * VectorClass().load_a(val0+i);
				vc_df_ptn += vc_res * VectorClass().load_a(val1+i);
				vc_ddf_ptn += vc_res * VectorClass().load_a(val2+i);
			}
			lh_ptn = horizontal_add(vc_ptn);
			df_ptn = horizontal_add(vc_df_ptn);
			ddf_ptn = horizontal_add(vc_ddf_ptn);

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
			partial_lh_node += block;
			partial_lh_dad += block;
		}
    }
    aligned_free(val2);
    aligned_free(val1);
    aligned_free(val0);
    return tree_lh;
}

template<const int nstates>
double PhyloTree::computeLikelihoodBranchEigenSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
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
        computePartialLikelihoodEigenSSE<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSSE<nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();
    //const size_t ncat = 4;

    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    //size_t nstates = aln->num_states;
    //const size_t nstates = 4;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *partial_lh_dad = dad_branch->partial_lh;
    double *partial_lh_node = node_branch->partial_lh;
    double *val = aligned_alloc_double(block);
	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		for (i = 0; i < nstates; i++)
			val[c*nstates+i] = exp(eval[i]*len);
	}
    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn;
			VectorClass vc_ptn = 0.0;
			for (i = 0; i < block; i+=VectorClass::size()) {
				vc_ptn += VectorClass().load_a(partial_lh_node+i%nstates) * VectorClass().load_a(partial_lh_dad+i) *
						VectorClass().load_a(val+i);
			}
			lh_ptn = horizontal_add(vc_ptn);
			lh_ptn *= p_var_cat;
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
			partial_lh_node += nstates;
			partial_lh_dad += block;
		}
    } else {
    	// both dad and node are internal nodes
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn;
			/*
			for (i = 0; i < block; i++)
				lh_ptn +=  val[i] * partial_lh_node[i] * partial_lh_dad[i];*/
			VectorClass vc_ptn = 0.0;
			for (i = 0; i < block; i+=VectorClass::size()) {
				vc_ptn += VectorClass().load_a(partial_lh_node+i) * VectorClass().load_a(partial_lh_dad+i) * VectorClass().load_a(val+i);
			}
			lh_ptn = horizontal_add(vc_ptn);
			lh_ptn *= p_var_cat;
			//assert(lh_ptn > 0.0);
			lh_ptn = log(lh_ptn);
			_pattern_lh[ptn] = lh_ptn;
			tree_lh += lh_ptn * aln->at(ptn).frequency;
			partial_lh_node += block;
			partial_lh_dad += block;
		}
    }
    if (pattern_lh)
        memmove(pattern_lh, _pattern_lh, aln->size() * sizeof(double));
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
		case LK_EIGEN_SSE: computePartialLikelihoodEigenSSE<4>(dad_branch, dad, pattern_scale); break;
		case LK_NORMAL: computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
		}
		break;
	case 20:
		switch(sse) {
		case LK_SSE: computePartialLikelihoodSSE<20>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN: computePartialLikelihoodEigen<20>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN_SSE: computePartialLikelihoodEigenSSE<20>(dad_branch, dad, pattern_scale); break;
		case LK_NORMAL: computePartialLikelihoodNaive(dad_branch, dad, pattern_scale); break;
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: computePartialLikelihoodSSE<2>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN: computePartialLikelihoodEigen<2>(dad_branch, dad, pattern_scale); break;
		case LK_EIGEN_SSE: computePartialLikelihoodEigenSSE<2>(dad_branch, dad, pattern_scale); break;
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
		case LK_EIGEN_SSE: return computeLikelihoodBranchEigenSSE<4>(dad_branch, dad, pattern_lh);
		case LK_NORMAL: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
		}
		break;
	case 20:
		switch(sse) {
		case LK_SSE: return computeLikelihoodBranchSSE<20>(dad_branch, dad, pattern_lh);
		case LK_EIGEN: return computeLikelihoodBranchEigen<20>(dad_branch, dad, pattern_lh);
		case LK_EIGEN_SSE: return computeLikelihoodBranchEigenSSE<20>(dad_branch, dad, pattern_lh);
		case LK_NORMAL: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: return computeLikelihoodBranchSSE<2>(dad_branch, dad, pattern_lh);
		case LK_EIGEN: return computeLikelihoodBranchEigen<2>(dad_branch, dad, pattern_lh);
		case LK_EIGEN_SSE: return computeLikelihoodBranchEigenSSE<2>(dad_branch, dad, pattern_lh);
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
    	case LK_EIGEN_SSE: return computeLikelihoodDervEigenSSE<4>(dad_branch, dad, df, ddf);
    	case LK_NORMAL: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
    	}
    	break;
	case 20:
		switch(sse) {
		case LK_SSE: return computeLikelihoodDervSSE<20>(dad_branch, dad, df, ddf);
		case LK_EIGEN: return computeLikelihoodDervEigen<20>(dad_branch, dad, df, ddf);
		case LK_EIGEN_SSE: return computeLikelihoodDervEigenSSE<20>(dad_branch, dad, df, ddf);
		case LK_NORMAL: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
		}
		break;
	case 2:
		switch(sse) {
		case LK_SSE: return computeLikelihoodDervSSE<2>(dad_branch, dad, df, ddf);
		case LK_EIGEN: return computeLikelihoodDervEigen<2>(dad_branch, dad, df, ddf);
		case LK_EIGEN_SSE: return computeLikelihoodDervEigenSSE<2>(dad_branch, dad, df, ddf);
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
