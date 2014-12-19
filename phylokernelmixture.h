/*
 * phylokernelmixture.h
 *
 *  Created on: Dec 19, 2014
 *      Author: minh
 */

#ifndef PHYLOKERNELMIXTURE_H_
#define PHYLOKERNELMIXTURE_H_


/************************************************************************************************
 *
 *   Highly optimized vectorized versions of likelihood functions
 *
 *************************************************************************************************/


template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeMixturePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;

    size_t nptn = aln->size() + model_factory->unobserved_ptns.size();
    PhyloNode *node = (PhyloNode*)(dad_branch->node);

	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;
	    //memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

		if (!tip_partial_lh_computed)
			computeTipPartialLikelihood();
		return;
	}

    size_t ptn, c;
    size_t orig_ntn = aln->size();

    size_t ncat = site_rate->getNRate();
    assert(nstates == aln->num_states && nstates >= VCSIZE && VCSIZE == VectorClass().size());
    assert(model->isReversible()); // only works with reversible model!
    const size_t nstatesqr=nstates*nstates;
    size_t i, x, j;
    size_t block = nstates * ncat;

	// internal node
	assert(node->degree() == 3); // it works only for strictly bifurcating tree
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
	}

	if (!left->node->isLeaf() && right->node->isLeaf()) {
		// swap left and right
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
	}
	if ((left->partial_lh_computed & 1) == 0)
		computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(left, node);
	if ((right->partial_lh_computed & 1) == 0)
		computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(right, node);

	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();

	VectorClass vc_inv_evec[nstates*nstates/VCSIZE];
	assert(inv_evec && evec);
	for (i = 0; i < nstates; i++) {
		for (x = 0; x < nstates/VCSIZE; x++)
			// inv_evec is not aligned!
			vc_inv_evec[i*nstates/VCSIZE+x].load(&inv_evec[i*nstates+x*VCSIZE]);
	}
	double *eval = model->getEigenvalues();

	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;

	VectorClass *eleft = (VectorClass*)aligned_alloc<double>(block*nstates);
	VectorClass *eright = (VectorClass*)aligned_alloc<double>(block*nstates);

	// precompute information buffer
	for (c = 0; c < ncat; c++) {
		VectorClass vc_evec;
		VectorClass expleft[nstates/VCSIZE];
		VectorClass expright[nstates/VCSIZE];
		double len_left = site_rate->getRate(c) * left->length;
		double len_right = site_rate->getRate(c) * right->length;
		for (i = 0; i < nstates/VCSIZE; i++) {
			// eval is not aligned!
			expleft[i] = exp(VectorClass().load(&eval[i*VCSIZE]) * VectorClass(len_left));
			expright[i] = exp(VectorClass().load(&eval[i*VCSIZE]) * VectorClass(len_right));
		}
		for (x = 0; x < nstates; x++)
			for (i = 0; i < nstates/VCSIZE; i++) {
				// evec is not be aligned!
				vc_evec.load(&evec[x*nstates+i*VCSIZE]);
				eleft[c*nstatesqr/VCSIZE+x*nstates/VCSIZE+i] = (vc_evec * expleft[i]);
				eright[c*nstatesqr/VCSIZE+x*nstates/VCSIZE+i] = (vc_evec * expright[i]);
			}
	}

	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);
		double *partial_lh_right = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);

		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
			VectorClass vleft[VCSIZE];
			size_t addr = state*nstates;
			for (i = 0; i < nstates/VCSIZE; i++)
				vc_partial_lh_tmp[i].load_a(&tip_partial_lh[addr+i*VCSIZE]);
			for (x = 0; x < block; x+=VCSIZE) {
				addr = x*nstates/VCSIZE;
				for (j = 0; j < VCSIZE; j++)
					vleft[j] = eleft[addr+j*nstates/VCSIZE] * vc_partial_lh_tmp[0];
				for (i = 1; i < nstates/VCSIZE; i++) {
					for (j = 0; j < VCSIZE; j++)
						vleft[j] = mul_add(eleft[addr+j*nstates/VCSIZE+i], vc_partial_lh_tmp[i], vleft[j]);
				}
				horizontal_add(vleft).store_a(&partial_lh_left[state*block+x]);
			}
		}

		for (it = aln->seq_states[right->node->id].begin(); it != aln->seq_states[right->node->id].end(); it++) {
			int state = (*it);
			VectorClass vright[VCSIZE];
			VectorClass vc_partial_lh_tmp[nstates/VCSIZE];

			for (i = 0; i < nstates/VCSIZE; i++)
				vc_partial_lh_tmp[i].load_a(&tip_partial_lh[state*nstates+i*VCSIZE]);
			for (x = 0; x < block; x+=VCSIZE) {
				for (j = 0; j < VCSIZE; j++)
					vright[j] = eright[(x+j)*nstates/VCSIZE] * vc_partial_lh_tmp[0];
				for (i = 1; i < nstates/VCSIZE; i++) {
					for (j = 0; j < VCSIZE; j++)
						vright[j] = mul_add(eright[(x+j)*nstates/VCSIZE+i], vc_partial_lh_tmp[i], vright[j]);
				}
				horizontal_add(vright).store_a(&partial_lh_right[state*block+x]);
			}
		}

		size_t addr_unknown = aln->STATE_UNKNOWN * block;
		for (x = 0; x < block; x++) {
			partial_lh_left[addr_unknown+x] = 1.0;
			partial_lh_right[addr_unknown+x] = 1.0;
		}

		// assign pointers for left and right partial_lh
		double **lh_left_ptr = aligned_alloc<double*>(nptn);
		double **lh_right_ptr = aligned_alloc<double*>(nptn);
		for (ptn = 0; ptn < orig_ntn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
			lh_right_ptr[ptn] = &partial_lh_right[block * (aln->at(ptn))[right->node->id]];
		}
		for (ptn = orig_ntn; ptn < nptn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
			lh_right_ptr[ptn] = &partial_lh_right[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
		}

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass res[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i, j, vc_partial_lh_tmp, res)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;

	        double *lh_left = lh_left_ptr[ptn];
	        double *lh_right = lh_right_ptr[ptn];
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector

				for (x = 0; x < nstates/VCSIZE; x++) {
					vc_partial_lh_tmp[x] = (VectorClass().load_a(&lh_left[x*VCSIZE]) * VectorClass().load_a(&lh_right[x*VCSIZE]));
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++)
						for (j = 0; j < VCSIZE; j++) {
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec[(i+j)*nstates/VCSIZE+x], res[j]);
						}
					horizontal_add(res).store_a(&partial_lh[i]);
				}

				lh_left += nstates;
				lh_right += nstates;
				partial_lh += nstates;
			}
		}

	    aligned_free(lh_left_ptr);
	    aligned_free(lh_right_ptr);
		aligned_free(partial_lh_right);
		aligned_free(partial_lh_left);
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

		// pre compute information for left tip
		double *partial_lh_left = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);


		vector<int>::iterator it;
		for (it = aln->seq_states[left->node->id].begin(); it != aln->seq_states[left->node->id].end(); it++) {
			int state = (*it);
			VectorClass vc_tip_lh[nstates/VCSIZE];
			VectorClass vleft[VCSIZE];
			for (i = 0; i < nstates/VCSIZE; i++)
				vc_tip_lh[i].load_a(&tip_partial_lh[state*nstates+i*VCSIZE]);
			for (x = 0; x < block; x+=VCSIZE) {
				for (j = 0; j < VCSIZE; j++)
					vleft[j] = eleft[(x+j)*nstates/VCSIZE] * vc_tip_lh[0];
				for (i = 1; i < nstates/VCSIZE; i++) {
					for (j = 0; j < VCSIZE; j++)
						vleft[j] = mul_add(eleft[(x+j)*nstates/VCSIZE+i], vc_tip_lh[i], vleft[j]);
				}
				horizontal_add(vleft).store_a(&partial_lh_left[state*block+x]);
			}
		}

		size_t addr_unknown = aln->STATE_UNKNOWN * block;
		for (x = 0; x < block; x++) {
			partial_lh_left[addr_unknown+x] = 1.0;
		}

		// assign pointers for partial_lh_left
		double **lh_left_ptr = aligned_alloc<double*>(nptn);
		for (ptn = 0; ptn < orig_ntn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
		}
		for (ptn = orig_ntn; ptn < nptn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
		}

		double sum_scale = 0.0;
		VectorClass vc_lh_right[nstates/VCSIZE];
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass res[VCSIZE];
		VectorClass vc_max; // maximum of partial likelihood, for scaling check
		VectorClass vright[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private (ptn, c, x, i, j, vc_lh_right, vc_partial_lh_tmp, res, vc_max, vright)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;
	        double *partial_lh_right = right->partial_lh + ptn*block;

	        double *lh_left = lh_left_ptr[ptn];
			vc_max = 0.0;
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (i = 0; i < nstates/VCSIZE; i++)
					vc_lh_right[i].load_a(&partial_lh_right[i*VCSIZE]);

				for (x = 0; x < nstates/VCSIZE; x++) {
					size_t addr = c*nstatesqr/VCSIZE+x*nstates;
					for (j = 0; j < VCSIZE; j++) {
						vright[j] = eright[addr+nstates*j/VCSIZE] * vc_lh_right[0];
					}
					for (i = 1; i < nstates/VCSIZE; i++)
						for (j = 0; j < VCSIZE; j++) {
							vright[j] = mul_add(eright[addr+i+nstates*j/VCSIZE], vc_lh_right[i], vright[j]);
						}
					vc_partial_lh_tmp[x] = VectorClass().load_a(&lh_left[x*VCSIZE])
							* horizontal_add(vright);
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++) {
						for (j = 0; j < VCSIZE; j++) {
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec[(i+j)*nstates/VCSIZE+x], res[j]);
						}
					}
					VectorClass sum_res = horizontal_add(res);
					sum_res.store_a(&partial_lh[i]);
					vc_max = max(vc_max, abs(sum_res)); // take the maximum for scaling check
				}
				lh_left += nstates;
				partial_lh_right += nstates;
				partial_lh += nstates;
			}
            // check if one should scale partial likelihoods
			double lh_max = horizontal_max(vc_max);
            if (lh_max < SCALING_THRESHOLD) {
            	// now do the likelihood scaling
            	partial_lh -= block; // revert its pointer
            	VectorClass scale_thres(SCALING_THRESHOLD_INVER);
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
				dad_branch->scale_num[ptn] += 1;
				partial_lh += block; // increase the pointer again
            }

		}
		dad_branch->lh_scale_factor += sum_scale;

	    aligned_free(lh_left_ptr);
		aligned_free(partial_lh_left);

	} else {
		// both left and right are internal node

		double sum_scale = 0.0;
		VectorClass vc_max; // maximum of partial likelihood, for scaling check
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass vc_lh_left[nstates/VCSIZE], vc_lh_right[nstates/VCSIZE];
		VectorClass res[VCSIZE];
		VectorClass vleft[VCSIZE], vright[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for reduction (+: sum_scale) private(ptn, c, x, i, j, vc_max, vc_partial_lh_tmp, vc_lh_left, vc_lh_right, res, vleft, vright)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = left->partial_lh + ptn*block;
			double *partial_lh_right = right->partial_lh + ptn*block;

			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];
			vc_max = 0.0;
			for (c = 0; c < ncat; c++) {
				// compute real partial likelihood vector
				for (i = 0; i < nstates/VCSIZE; i++) {
					vc_lh_left[i].load_a(&partial_lh_left[i*VCSIZE]);
					vc_lh_right[i].load_a(&partial_lh_right[i*VCSIZE]);
				}

				for (x = 0; x < nstates/VCSIZE; x++) {
					size_t addr = c*nstatesqr/VCSIZE+x*nstates;
					for (j = 0; j < VCSIZE; j++) {
						size_t addr_com = addr+j*nstates/VCSIZE;
						vleft[j] = eleft[addr_com] * vc_lh_left[0];
						vright[j] = eright[addr_com] * vc_lh_right[0];
					}
					for (i = 1; i < nstates/VCSIZE; i++) {
						for (j = 0; j < VCSIZE; j++) {
							size_t addr_com = addr+i+j*nstates/VCSIZE;
							vleft[j] = mul_add(eleft[addr_com], vc_lh_left[i], vleft[j]);
							vright[j] = mul_add(eright[addr_com], vc_lh_right[i], vright[j]);
						}
					}
					vc_partial_lh_tmp[x] = horizontal_add(vleft) * horizontal_add(vright);
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++)
						for (j = 0; j < VCSIZE; j++)
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec[(i+j)*nstates/VCSIZE+x], res[j]);

					VectorClass sum_res = horizontal_add(res);
					sum_res.store_a(&partial_lh[i]);
					vc_max = max(vc_max, abs(sum_res)); // take the maximum for scaling check
				}
				partial_lh += nstates;
				partial_lh_left += nstates;
				partial_lh_right += nstates;
			}

            // check if one should scale partial likelihoods
			double lh_max = horizontal_max(vc_max);
            if (lh_max < SCALING_THRESHOLD) {
				// now do the likelihood scaling
            	partial_lh -= block; // revert its pointer
            	VectorClass scale_thres(SCALING_THRESHOLD_INVER);
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
				}
				// unobserved const pattern will never have underflow
				sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
				dad_branch->scale_num[ptn] += 1;
				partial_lh += block; // increase the pointer again
            }

		}
		dad_branch->lh_scale_factor += sum_scale;

	}

	aligned_free(eright);
	aligned_free(eleft);
}

template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeMixtureLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node);
    df = ddf = 0.0;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    double *eval = model->getEigenvalues();
    assert(eval);

	VectorClass *vc_val0 = (VectorClass*)aligned_alloc<double>(block);
	VectorClass *vc_val1 = (VectorClass*)aligned_alloc<double>(block);
	VectorClass *vc_val2 = (VectorClass*)aligned_alloc<double>(block);

	VectorClass vc_len = dad_branch->length;
	for (c = 0; c < ncat; c++) {
		VectorClass vc_rate = site_rate->getRate(c);
		VectorClass vc_prop = site_rate->getProp(c);
		for (i = 0; i < nstates/VCSIZE; i++) {
			VectorClass cof = VectorClass().load(&eval[i*VCSIZE]) * vc_rate;
			VectorClass val = exp(cof*vc_len) * vc_prop;
			VectorClass val1_ = cof*val;
			vc_val0[c*nstates/VCSIZE+i] = val;
			vc_val1[c*nstates/VCSIZE+i] = val1_;
			vc_val2[c*nstates/VCSIZE+i] = cof*val1_;
		}
	}

	assert(theta_all);
	if (!theta_computed) {
		theta_computed = true;
		// precompute theta for fast branch length optimization

		if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i)
#endif
			for (ptn = 0; ptn < orig_nptn; ptn++) {
			    double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
				double *theta = theta_all + ptn*block;
				double *lh_dad = &tip_partial_lh[(aln->at(ptn))[dad->id] * nstates];
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&lh_dad[i%nstates]) * VectorClass().load_a(&partial_lh_dad[i])).store_a(&theta[i]);
				}
			}
			// ascertainment bias correction
			for (ptn = orig_nptn; ptn < nptn; ptn++) {
			    double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
				double *theta = theta_all + ptn*block;
				double *lh_dad = &tip_partial_lh[model_factory->unobserved_ptns[ptn-orig_nptn] * nstates];
				for (i = 0; i < block; i+=VCSIZE) {
					(VectorClass().load_a(&lh_dad[i%nstates]) * VectorClass().load_a(&partial_lh_dad[i])).store_a(&theta[i]);
				}
			}
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
		    double *partial_lh_dad = dad_branch->partial_lh;
	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	    	for (i = 0; i < all_entries; i+=VCSIZE) {
				(VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i]))
						.store_a(&theta_all[i]);
			}
	    }
		if (nptn < maxptn) {
			// copy dummy values
			for (ptn = nptn; ptn < maxptn; ptn++)
				memcpy(&theta_all[ptn*block], theta_all, block*sizeof(double));
		}
	}



	VectorClass vc_ptn[VCSIZE], vc_df[VCSIZE], vc_ddf[VCSIZE], vc_theta[VCSIZE];
	VectorClass vc_unit = 1.0;
	VectorClass vc_freq;
	VectorClass df_final = 0.0, ddf_final = 0.0;
	// these stores values of 2 consecutive patterns
	VectorClass lh_ptn, df_ptn, ddf_ptn, inv_lh_ptn;

	// perform 2 sites at the same time for SSE/AVX efficiency

#ifdef _OPENMP
#pragma omp parallel private (ptn, i, j, vc_freq, vc_ptn, vc_df, vc_ddf, vc_theta, inv_lh_ptn, lh_ptn, df_ptn, ddf_ptn)
	{
	VectorClass df_final_th = 0.0;
	VectorClass ddf_final_th = 0.0;
#pragma omp for nowait
#endif
	for (ptn = 0; ptn < orig_nptn; ptn+=VCSIZE) {
		double *theta = theta_all + ptn*block;
		// initialization
		for (i = 0; i < VCSIZE; i++) {
			vc_theta[i].load_a(theta+i*block);
			vc_ptn[i] = vc_val0[0] * vc_theta[i];
			vc_df[i] = vc_val1[0] * vc_theta[i];
			vc_ddf[i] = vc_val2[0] * vc_theta[i];
		}

		for (i = 1; i < block/VCSIZE; i++) {
			for (j = 0; j < VCSIZE; j++) {
				vc_theta[j].load_a(&theta[i*VCSIZE+j*block]);
				vc_ptn[j] = mul_add(vc_theta[j], vc_val0[i], vc_ptn[j]);
				vc_df[j] = mul_add(vc_theta[j], vc_val1[i], vc_df[j]);
				vc_ddf[j] = mul_add(vc_theta[j], vc_val2[i], vc_ddf[j]);
			}
		}
		lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

		inv_lh_ptn = vc_unit / lh_ptn;

		vc_freq.load_a(&ptn_freq[ptn]);

		df_ptn = horizontal_add(vc_df) * inv_lh_ptn;
		ddf_ptn = horizontal_add(vc_ddf) * inv_lh_ptn;
		ddf_ptn = nmul_add(df_ptn, df_ptn, ddf_ptn);

#ifdef _OPENMP
		df_final_th = mul_add(df_ptn, vc_freq, df_final_th);
		ddf_final_th = mul_add(ddf_ptn, vc_freq, ddf_final_th);
#else
		df_final = mul_add(df_ptn, vc_freq, df_final);
		ddf_final = mul_add(ddf_ptn, vc_freq, ddf_final);
#endif

	}

#ifdef _OPENMP
#pragma omp critical
	{
		df_final += df_final_th;
		ddf_final += ddf_final_th;
	}
}
#endif
	df = horizontal_add(df_final);
	ddf = horizontal_add(ddf_final);

//	assert(isnormal(tree_lh));
	if (orig_nptn < nptn) {
		// ascertaiment bias correction
		VectorClass lh_final = 0.0;
		df_final = 0.0;
		ddf_final = 0.0;
		lh_ptn = 0.0;
		df_ptn = 0.0;
		ddf_ptn = 0.0;
		double prob_const, df_const, ddf_const;
		double *theta = &theta_all[orig_nptn*block];
		for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
			lh_final += lh_ptn;
			df_final += df_ptn;
			ddf_final += ddf_ptn;

			// initialization
			for (i = 0; i < VCSIZE; i++) {
				vc_theta[i].load_a(theta+i*block);
				vc_ptn[i] = vc_val0[0] * vc_theta[i];
				vc_df[i] = vc_val1[0] * vc_theta[i];
				vc_ddf[i] = vc_val2[0] * vc_theta[i];
			}

			for (i = 1; i < block/VCSIZE; i++) {
				for (j = 0; j < VCSIZE; j++) {
					vc_theta[j].load_a(&theta[i*VCSIZE+j*block]);
					vc_ptn[j] = mul_add(vc_theta[j], vc_val0[i], vc_ptn[j]);
					vc_df[j] = mul_add(vc_theta[j], vc_val1[i], vc_df[j]);
					vc_ddf[j] = mul_add(vc_theta[j], vc_val2[i], vc_ddf[j]);
				}
			}
			theta += block*VCSIZE;

			// ptn_invar[ptn] is not aligned
			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);

		}
		switch ((nptn-orig_nptn) % VCSIZE) {
		case 0:
			prob_const = horizontal_add(lh_final+lh_ptn);
			df_const = horizontal_add(df_final+df_ptn);
			ddf_const = horizontal_add(ddf_final+ddf_ptn);
			break;
		case 1:
			prob_const = horizontal_add(lh_final)+lh_ptn[0];
			df_const = horizontal_add(df_final)+df_ptn[0];
			ddf_const = horizontal_add(ddf_final)+ddf_ptn[0];
			break;
		case 2:
			prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1];
			df_const = horizontal_add(df_final)+df_ptn[0]+df_ptn[1];
			ddf_const = horizontal_add(ddf_final)+ddf_ptn[0]+ddf_ptn[1];
			break;
		case 3:
			prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2];
			df_const = horizontal_add(df_final)+df_ptn[0]+df_ptn[1]+df_ptn[2];
			ddf_const = horizontal_add(ddf_final)+ddf_ptn[0]+ddf_ptn[1]+ddf_ptn[2];
			break;
		default:
			assert(0);
			break;
		}
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
	}

    aligned_free(vc_val2);
    aligned_free(vc_val1);
    aligned_free(vc_val0);
}


template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeMixtureLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
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
        computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeMixturePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    double *eval = model->getEigenvalues();
    assert(eval);

    VectorClass *vc_val = (VectorClass*)aligned_alloc<double>(block);


	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		VectorClass vc_len(len);
		VectorClass vc_prop(site_rate->getProp(c));
		for (i = 0; i < nstates/VCSIZE; i++) {
			// eval is not aligned!
			vc_val[c*nstates/VCSIZE+i] = exp(VectorClass().load(&eval[i*VCSIZE]) * vc_len) * vc_prop;
		}
	}

	double prob_const = 0.0;

	if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
    	VectorClass vc_tip_partial_lh[nstates];
    	VectorClass vc_partial_lh_dad[VCSIZE], vc_ptn[VCSIZE];
    	VectorClass lh_final(0.0), vc_freq;
		VectorClass lh_ptn; // store likelihoods of VCSIZE consecutive patterns

    	double **lh_states_dad = aligned_alloc<double*>(maxptn);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		lh_states_dad[ptn] = &tip_partial_lh[(aln->at(ptn))[dad->id] * nstates];
    	for (ptn = orig_nptn; ptn < nptn; ptn++)
    		lh_states_dad[ptn] = &tip_partial_lh[model_factory->unobserved_ptns[ptn-orig_nptn] * nstates];
    	// initialize beyond #patterns for efficiency
    	for (ptn = nptn; ptn < maxptn; ptn++)
    		lh_states_dad[ptn] = &tip_partial_lh[aln->STATE_UNKNOWN * nstates];

		// copy dummy values because VectorClass will access beyond nptn
		for (ptn = nptn; ptn < maxptn; ptn++)
			memcpy(&dad_branch->partial_lh[ptn*block], dad_branch->partial_lh, block*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel private(ptn, i, j, vc_tip_partial_lh, vc_partial_lh_dad, vc_ptn, vc_freq, lh_ptn)
    {
    	VectorClass lh_final_th = 0.0;
#pragma omp for nowait
#endif
   		// main loop over all patterns with a step size of VCSIZE
		for (ptn = 0; ptn < orig_nptn; ptn+=VCSIZE) {
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;

			// initialize vc_tip_partial_lh
			for (j = 0; j < VCSIZE; j++) {
				double *lh_dad = lh_states_dad[ptn+j];
				for (i = 0; i < nstates/VCSIZE; i++) {
					vc_tip_partial_lh[j*(nstates/VCSIZE)+i].load_a(&lh_dad[i*VCSIZE]);
				}
				vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block]);
				vc_ptn[j] = vc_val[0] * vc_tip_partial_lh[j*(nstates/VCSIZE)] * vc_partial_lh_dad[j];
			}

			// compute vc_ptn
			for (i = 1; i < block/VCSIZE; i++)
				for (j = 0; j < VCSIZE; j++) {
					vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block+i*VCSIZE]);
					vc_ptn[j] = mul_add(vc_val[i] * vc_tip_partial_lh[j*(nstates/VCSIZE)+i%(nstates/VCSIZE)],
							vc_partial_lh_dad[j], vc_ptn[j]);
				}

			vc_freq.load_a(&ptn_freq[ptn]);
			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
			lh_ptn = log(lh_ptn);
			lh_ptn.store_a(&_pattern_lh[ptn]);

			// multiply with pattern frequency
#ifdef _OPENMP
			lh_final_th = mul_add(lh_ptn, vc_freq, lh_final_th);
#else
			lh_final = mul_add(lh_ptn, vc_freq, lh_final);
#endif
		}

#ifdef _OPENMP
#pragma omp critical
		{
			lh_final += lh_final_th;
    	}
    }
#endif
		tree_lh += horizontal_add(lh_final);
		assert(!isnan(tree_lh) && !isinf(tree_lh));

		// ascertainment bias correction
		if (orig_nptn < nptn) {
			lh_final = 0.0;
			lh_ptn = 0.0;
			for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
				double *partial_lh_dad = &dad_branch->partial_lh[ptn*block];
				lh_final += lh_ptn;

				// initialize vc_tip_partial_lh
				for (j = 0; j < VCSIZE; j++) {
					double *lh_dad = lh_states_dad[ptn+j];
					for (i = 0; i < nstates/VCSIZE; i++) {
						vc_tip_partial_lh[j*(nstates/VCSIZE)+i].load_a(&lh_dad[i*VCSIZE]);
					}
					vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block]);
					vc_ptn[j] = vc_val[0] * vc_tip_partial_lh[j*(nstates/VCSIZE)] * vc_partial_lh_dad[j];
				}

				// compute vc_ptn
				for (i = 1; i < block/VCSIZE; i++)
					for (j = 0; j < VCSIZE; j++) {
						vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block+i*VCSIZE]);
						vc_ptn[j] = mul_add(vc_val[i] * vc_tip_partial_lh[j*(nstates/VCSIZE)+i%(nstates/VCSIZE)],
								vc_partial_lh_dad[j], vc_ptn[j]);
					}
				// ptn_invar[ptn] is not aligned
				lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);
			}
			switch ((nptn-orig_nptn)%VCSIZE) {
			case 0: prob_const = horizontal_add(lh_final+lh_ptn); break;
			case 1: prob_const = horizontal_add(lh_final)+lh_ptn[0]; break;
			case 2: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]; break;
			case 3: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2]; break;
			default: assert(0); break;
			}
		}
		aligned_free(lh_states_dad);
    } else {
    	// both dad and node are internal nodes
    	VectorClass vc_partial_lh_node[VCSIZE];
    	VectorClass vc_partial_lh_dad[VCSIZE], vc_ptn[VCSIZE];
    	VectorClass lh_final(0.0), vc_freq;
		VectorClass lh_ptn;

		// copy dummy values because VectorClass will access beyond nptn
		for (ptn = nptn; ptn < maxptn; ptn++) {
			memcpy(&dad_branch->partial_lh[ptn*block], dad_branch->partial_lh, block*sizeof(double));
			memcpy(&node_branch->partial_lh[ptn*block], node_branch->partial_lh, block*sizeof(double));
		}

#ifdef _OPENMP
#pragma omp parallel private(ptn, i, j, vc_partial_lh_node, vc_partial_lh_dad, vc_ptn, vc_freq, lh_ptn)
		{
		VectorClass lh_final_th = 0.0;
#pragma omp for nowait
#endif
		for (ptn = 0; ptn < orig_nptn; ptn+=VCSIZE) {
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
			double *partial_lh_node = node_branch->partial_lh + ptn*block;

			for (j = 0; j < VCSIZE; j++)
				vc_ptn[j] = 0.0;

			for (i = 0; i < block; i+=VCSIZE) {
				for (j = 0; j < VCSIZE; j++) {
					vc_partial_lh_node[j].load_a(&partial_lh_node[i+j*block]);
					vc_partial_lh_dad[j].load_a(&partial_lh_dad[i+j*block]);
					vc_ptn[j] = mul_add(vc_val[i/VCSIZE] * vc_partial_lh_node[j], vc_partial_lh_dad[j], vc_ptn[j]);
				}
			}

			vc_freq.load_a(&ptn_freq[ptn]);

			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

			lh_ptn = log(lh_ptn);
			lh_ptn.store_a(&_pattern_lh[ptn]);
#ifdef _OPENMP
			lh_final_th = mul_add(lh_ptn, vc_freq, lh_final_th);
#else
			lh_final = mul_add(lh_ptn, vc_freq, lh_final);
#endif
		}
#ifdef _OPENMP
#pragma omp critical
		{
			lh_final += lh_final_th;
		}
	}
#endif

		tree_lh += horizontal_add(lh_final);
		assert(!isnan(tree_lh) && !isinf(tree_lh));

		if (orig_nptn < nptn) {
			// ascertainment bias correction
			lh_final = 0.0;
			lh_ptn = 0.0;
			double *partial_lh_node = &node_branch->partial_lh[orig_nptn*block];
			double *partial_lh_dad = &dad_branch->partial_lh[orig_nptn*block];

			for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
				lh_final += lh_ptn;

				for (j = 0; j < VCSIZE; j++)
					vc_ptn[j] = 0.0;

				for (i = 0; i < block; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						vc_partial_lh_node[j].load_a(&partial_lh_node[i+j*block]);
						vc_partial_lh_dad[j].load_a(&partial_lh_dad[i+j*block]);
						vc_ptn[j] = mul_add(vc_val[i/VCSIZE] * vc_partial_lh_node[j], vc_partial_lh_dad[j], vc_ptn[j]);
					}
				}

				// ptn_invar[ptn] is not aligned
				lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);
				partial_lh_node += block*VCSIZE;
				partial_lh_dad += block*VCSIZE;
			}
			switch ((nptn-orig_nptn)%VCSIZE) {
			case 0: prob_const = horizontal_add(lh_final+lh_ptn); break;
			case 1: prob_const = horizontal_add(lh_final)+lh_ptn[0]; break;
			case 2: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]; break;
			case 3: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2]; break;
			default: assert(0); break;
			}
		}
    }

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
    }

    aligned_free(vc_val);
    return tree_lh;
}

template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD() {


	assert(theta_all && theta_computed);

	double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;

    size_t ncat = site_rate->getNRate();
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    double *eval = model->getEigenvalues();
    assert(eval);

	VectorClass *vc_val0 = (VectorClass*)aligned_alloc<double>(block);

	VectorClass vc_len = current_it->length;
	for (c = 0; c < ncat; c++) {
		VectorClass vc_rate = site_rate->getRate(c);
		VectorClass vc_prop = site_rate->getProp(c);
		for (i = 0; i < nstates/VCSIZE; i++) {
			VectorClass cof = VectorClass().load(&eval[i*VCSIZE]) * vc_rate;
			VectorClass val = exp(cof*vc_len) * vc_prop;
			vc_val0[c*nstates/VCSIZE+i] = val;
		}
	}

	VectorClass vc_ptn[VCSIZE];
	VectorClass vc_freq;
	VectorClass lh_final = 0.0;
	// these stores values of 2 consecutive patterns
	VectorClass lh_ptn;

	// perform 2 sites at the same time for SSE/AVX efficiency

#ifdef _OPENMP
#pragma omp parallel private (ptn, i, j, vc_freq, vc_ptn, lh_ptn)
	{
	VectorClass lh_final_th = 0.0;
#pragma omp for nowait
#endif
	for (ptn = 0; ptn < orig_nptn; ptn+=VCSIZE) {
		double *theta = theta_all + ptn*block;
		// initialization
		for (i = 0; i < VCSIZE; i++) {
			vc_ptn[i] = vc_val0[0] * VectorClass().load_a(theta+i*block);
		}

		for (i = 1; i < block/VCSIZE; i++) {
			for (j = 0; j < VCSIZE; j++) {
				vc_ptn[j] = mul_add(VectorClass().load_a(&theta[i*VCSIZE+j*block]), vc_val0[i], vc_ptn[j]);
			}
		}
		lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
		lh_ptn = log(lh_ptn);
		lh_ptn.store_a(&_pattern_lh[ptn]);
		vc_freq.load_a(&ptn_freq[ptn]);

#ifdef _OPENMP
		lh_final_th = mul_add(lh_ptn, vc_freq, lh_final_th);
#else
		lh_final = mul_add(lh_ptn, vc_freq, lh_final);
#endif

	}

#ifdef _OPENMP
#pragma omp critical
	{
		lh_final += lh_final_th;
	}
}
#endif
	tree_lh += horizontal_add(lh_final);

	assert(isnormal(tree_lh));
	if (orig_nptn < nptn) {
		// ascertaiment bias correction
		lh_final = 0.0;
		lh_ptn = 0.0;
		double prob_const;// df_const, ddf_const;
		double *theta = &theta_all[orig_nptn*block];
		for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
			lh_final += lh_ptn;

			// initialization
			for (i = 0; i < VCSIZE; i++) {
				vc_ptn[i] = vc_val0[0] * VectorClass().load_a(theta+i*block);
			}

			for (i = 1; i < block/VCSIZE; i++) {
				for (j = 0; j < VCSIZE; j++) {
					vc_ptn[j] = mul_add(VectorClass().load_a(&theta[i*VCSIZE+j*block]), vc_val0[i], vc_ptn[j]);
				}
			}
			theta += block*VCSIZE;

			// ptn_invar[ptn] is not aligned
			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);

		}
		switch ((nptn-orig_nptn) % VCSIZE) {
		case 0:
			prob_const = horizontal_add(lh_final+lh_ptn);
			break;
		case 1:
			prob_const = horizontal_add(lh_final)+lh_ptn[0];
			break;
		case 2:
			prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1];
			break;
		case 3:
			prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2];
			break;
		default:
			assert(0);
			break;
		}
    	prob_const = log(1.0 - prob_const);
    	tree_lh -= aln->getNSite() * prob_const;
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
	}

    aligned_free(vc_val0);

    return tree_lh;
}




#endif /* PHYLOKERNELMIXTURE_H_ */
