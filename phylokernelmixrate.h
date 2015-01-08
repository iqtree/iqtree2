/*
 * phylokernelmixrate.h
 *
 *  Created on: Jan 7, 2015
 *      Author: minh
 */

#ifndef PHYLOKERNELMIXRATE_H_
#define PHYLOKERNELMIXRATE_H_

#include "model/modelmixture.h"

template <const int nstates>
void PhyloTree::computeMixratePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {
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
    assert(ncat == model->getNMixtures());
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
		computeMixratePartialLikelihoodEigen<nstates>(left, node);
	if ((right->partial_lh_computed & 1) == 0)
		computeMixratePartialLikelihoodEigen<nstates>(right, node);
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
			expleft[i] = exp(eval[c*nstates+i]*len_left);
			expright[i] = exp(eval[c*nstates+i]*len_right);
		}
		for (x = 0; x < nstates; x++)
			for (i = 0; i < nstates; i++) {
				eleft[c*nstatesqr+x*nstates+i] = evec[c*nstatesqr+x*nstates+i] * expleft[i];
				eright[c*nstatesqr+x*nstates+i] = evec[c*nstatesqr+x*nstates+i] * expright[i];
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
			for (c = 0; c < ncat; c++)
			for (x = 0; x < nstates; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[c*nstatesqr+x*nstates+i] * tip_partial_lh[state*block+c*nstates+i];
				}
				partial_lh_left[state*block+c*nstates+x] = vleft;
			}
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
			for (c = 0; c < ncat; c++)
			for (x = 0; x < nstates; x++) {
				double vright = 0.0;
				for (i = 0; i < nstates; i++) {
					vright += eright[c*nstatesqr+x*nstates+i] * tip_partial_lh[state*block+c*nstates+i];
				}
				partial_lh_right[state*block+c*nstates+x] = vright;
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
						res += partial_lh_tmp[x]*inv_evec[c*nstatesqr+i*nstates+x];
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
			for (c = 0; c < ncat; c++)
			for (x = 0; x < nstates; x++) {
				double vleft = 0.0;
				for (i = 0; i < nstates; i++) {
					vleft += eleft[c*nstatesqr+x*nstates+i] * tip_partial_lh[state*block+c*nstates+i];
				}
				partial_lh_left[state*block+c*nstates+x] = vleft;
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
						res += partial_lh_tmp[x]*inv_evec[c*nstatesqr+i*nstates+x];
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
						res += partial_lh_tmp[x]*inv_evec[c*nstatesqr+i*nstates+x];
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
void PhyloTree::computeMixrateLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computeMixratePartialLikelihoodEigen<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeMixratePartialLikelihoodEigen<nstates>(node_branch, node);
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
				double *lh_tip = tip_partial_lh + ((int)((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]))*nstates*ncat;
				for (i = 0; i < block; i++) {
					theta[i] = lh_tip[i] * partial_lh_dad[i];
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
			double cof = eval[c*nstates+i]*site_rate->getRate(c);
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
double PhyloTree::computeMixrateLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {
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
        computeMixratePartialLikelihoodEigen<nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeMixratePartialLikelihoodEigen<nstates>(node_branch, node);
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
			val[c*nstates+i] = exp(eval[c*nstates+i]*len) * prop;
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
    		double *lh_tip = tip_partial_lh + (*it)*block;
    		for (i = 0; i < block; i++)
    			lh_node[i] = val[i]*lh_tip[i];
//    		double *val_tmp = val;
//			for (c = 0; c < ncat; c++) {
//				for (i = 0; i < nstates; i++) {
//					  lh_node[i] = val_tmp[i] * lh_tip[i];
//				}
//				lh_node += nstates;
//				val_tmp += nstates;
//			}
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


#endif /* PHYLOKERNELMIXRATE_H_ */
