/*
 * phylokernelsitemodel.cpp
 * likelihood kernel site-specific frequency model
 *
 *  Created on: Jan 9, 2016
 *      Author: minh
 */



#include "tree/phylotree.h"
#include "model/modelset.h"

void PhyloTree::computeSitemodelPartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {

    // don't recompute the likelihood
	ASSERT(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);


    size_t nstates = aln->num_states;
    size_t nptn = aln->size(), tip_block_size = get_safe_upper_limit(nptn)*nstates;
    size_t ptn, c;
    size_t ncat = site_rate->getNRate();
    size_t i, x;
    size_t block = nstates * ncat;
    ModelSet *models = (ModelSet*) model;
    ASSERT(models->size() == nptn);


	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;
		// scale number must be ZERO
//	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

		if (!tip_partial_lh_computed)
			computeTipPartialLikelihood();
   
		return;
	}
    
    dad_branch->lh_scale_factor = 0.0;

	// internal node
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
	FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)*it;
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
        if ((nei->partial_lh_computed & 1) == 0)
            computePartialLikelihood(nei, node);
        dad_branch->lh_scale_factor += nei->lh_scale_factor;
	}

    if (params->lh_mem_save == LM_PER_NODE && !dad_branch->partial_lh) {
        // re-orient partial_lh
        bool done = false;
        FOR_NEIGHBOR_IT(node, dad, it2) {
            PhyloNeighbor *backnei = ((PhyloNeighbor*)(*it2)->node->findNeighbor(node));
            if (backnei->partial_lh) {
                dad_branch->partial_lh = backnei->partial_lh;
                dad_branch->scale_num = backnei->scale_num;
                backnei->partial_lh = NULL;
                backnei->scale_num = NULL;
                backnei->partial_lh_computed &= ~1; // clear bit
                done = true;
                break;
            }
        }
        ASSERT(done && "partial_lh is not re-oriented");
    }


    double sum_scale = 0.0;
        
	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
	}
    
    if (node->degree() > 3) {

        /*--------------------- multifurcating node ------------------*/
    
        // now for-loop computing partial_lh over all site-patterns
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn++) {
            double partial_lh_all[block];
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            dad_branch->scale_num[ptn] = 0;

            double expchild[nstates];
            double *eval = models->at(ptn)->getEigenvalues();
            double *evec = models->at(ptn)->getEigenvectors();
            double *inv_evec = models->at(ptn)->getInverseEigenvectors();

            FOR_NEIGHBOR_IT(node, dad, it) {
                PhyloNeighbor *child = (PhyloNeighbor*)*it;
                if (child->node->isLeaf()) {
                    // external node
                    double *tip_partial_lh_child = tip_partial_lh + (child->node->id*tip_block_size)+ptn*nstates;
                    double *partial_lh = partial_lh_all;
                    for (c = 0; c < ncat; c++) {
                        double len_child = site_rate->getRate(c) * child->length;
                        for (i = 0; i < nstates; i++) {
                            expchild[i] = exp(eval[i]*len_child) * tip_partial_lh_child[i];
                        }
                        // compute real partial likelihood vector
                        for (x = 0; x < nstates; x++) {
                            double vchild = 0.0;
                            double *this_evec = evec + x*nstates;
                            for (i = 0; i < nstates; i++) {
                                vchild += this_evec[i] * expchild[i];
                            }
                            partial_lh[x] *= vchild;
                        }
                        partial_lh += nstates;
                    }
                } else {
                    // internal node
                    dad_branch->scale_num[ptn] += child->scale_num[ptn];
                    
                    double *partial_lh_child = child->partial_lh + ptn*block;
                    double *partial_lh = partial_lh_all;
                    for (c = 0; c < ncat; c++) {
                        double len_child = site_rate->getRate(c) * child->length;
                        for (i = 0; i < nstates; i++) {
                            expchild[i] = exp(eval[i]*len_child) * partial_lh_child[i];
                        }
                        // compute real partial likelihood vector
                        for (x = 0; x < nstates; x++) {
                            double vchild = 0.0;
                            double *this_evec = evec + x*nstates;
                            for (i = 0; i < nstates; i++) {
                                vchild += this_evec[i] * expchild[i];
                            }
                            partial_lh[x] *= vchild;
                        }
                        partial_lh_child += nstates;                        
                        partial_lh += nstates;
                    }
                } // if
                
                
            } // FOR_NEIGHBOR
            
        
            // compute dot-product with inv_eigenvector
            double lh_max = 0.0;
            double *partial_lh_tmp = partial_lh_all;
            double *partial_lh = dad_branch->partial_lh + ptn*block;
            for (c = 0; c < ncat; c++) {
                double *inv_evec_ptr = inv_evec;
                for (i = 0; i < nstates; i++) {
                    double res = 0.0;
                    for (x = 0; x < nstates; x++) {
                        res += partial_lh_tmp[x]*inv_evec_ptr[x];
                    }
                    inv_evec_ptr += nstates;
                    partial_lh[i] = res;
                    lh_max = max(lh_max, fabs(res));
                }
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }
            // check if one should scale partial likelihoods
            if (lh_max < SCALING_THRESHOLD) {
                partial_lh = dad_branch->partial_lh + ptn*block;
                if (lh_max == 0.0) {
                    // for very shitty data
                    for (c = 0; c < ncat; c++)
                        memcpy(&partial_lh[c*nstates], &tip_partial_lh[ptn*nstates], nstates*sizeof(double));
                    sum_scale += LOG_SCALING_THRESHOLD* 4 * ptn_freq[ptn];
                    //sum_scale += log(lh_max) * ptn_freq[ptn];
                    dad_branch->scale_num[ptn] += 4;
                    int nsite = aln->getNSite();
                    for (i = 0, x = 0; i < nsite && x < ptn_freq[ptn]; i++)
                        if (aln->getPatternID(i) == ptn) {
                            outWarning((string)"Numerical underflow for site " + convertIntToString(i+1));
                            x++;
                        }
                } else {
                    // now do the likelihood scaling
                    for (i = 0; i < block; i++) {
                        partial_lh[i] *= SCALING_THRESHOLD_INVER;
                        //partial_lh[i] /= lh_max;
                    }
                    // unobserved const pattern will never have underflow
                    sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
                    //sum_scale += log(lh_max) * ptn_freq[ptn];
                    dad_branch->scale_num[ptn] += 1;
                }
            }

        } // for ptn
        dad_branch->lh_scale_factor += sum_scale;               
                
        // end multifurcating treatment
    } else 
    if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *tip_partial_lh_left = tip_partial_lh + (left->node->id * tip_block_size);
        double *tip_partial_lh_right = tip_partial_lh + (right->node->id * tip_block_size);

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double partial_lh_tmp[nstates];
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = tip_partial_lh_left + ptn*nstates;
			double *partial_lh_right = tip_partial_lh_right + ptn*nstates;

            double expleft[nstates];
            double expright[nstates];
            double *eval = models->at(ptn)->getEigenvalues();
            double *evec = models->at(ptn)->getEigenvectors();
            double *inv_evec = models->at(ptn)->getInverseEigenvectors();

			for (c = 0; c < ncat; c++) {
                double len_left = site_rate->getRate(c) * left->length;
                double len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }

				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
                    double *this_evec = evec + x*nstates;
					for (i = 0; i < nstates; i++) {
						vleft += this_evec[i] * expleft[i];
						vright += this_evec[i] * expright[i];
					}
					partial_lh_tmp[x] = vleft*vright;
				}

                // do not increase partial_lh_left and right for tips
                
				// compute dot-product with inv_eigenvector
                double *inv_evec_ptr = inv_evec;
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec_ptr[x];
					}
                    inv_evec_ptr += nstates;
					partial_lh[c*nstates+i] = res;
				}
			}

		}

	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

        double *tip_partial_lh_left = tip_partial_lh + (left->node->id * tip_block_size);

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double partial_lh_tmp[nstates];
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = tip_partial_lh_left + ptn*nstates;
			double *partial_lh_right = right->partial_lh + ptn*block;
            double lh_max = 0.0;

            double expleft[nstates];
            double expright[nstates];
            double *eval = models->at(ptn)->getEigenvalues();
            double *evec = models->at(ptn)->getEigenvectors();
            double *inv_evec = models->at(ptn)->getInverseEigenvectors();

			for (c = 0; c < ncat; c++) {
                double len_left = site_rate->getRate(c) * left->length;
                double len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
                    double *this_evec = evec + x*nstates;
					for (i = 0; i < nstates; i++) {
						vleft += this_evec[i] * expleft[i];
						vright += this_evec[i] * expright[i];
					}
					partial_lh_tmp[x] = vleft*vright;
				}
                // do not increase partial_lh_left for left tip
                partial_lh_right += nstates;
                
				// compute dot-product with inv_eigenvector
                double *inv_evec_ptr = inv_evec;
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec_ptr[x];
					}
                    inv_evec_ptr += nstates;
					partial_lh[c*nstates+i] = res;
                    lh_max = max(lh_max, fabs(res));
				}
			}

            // check if one should scale partial likelihoods
            if (lh_max < SCALING_THRESHOLD) {
            	if (lh_max == 0.0) {
            		// for very shitty data
            		for (c = 0; c < ncat; c++)
            			memcpy(&partial_lh[c*nstates], &tip_partial_lh[ptn*nstates], nstates*sizeof(double));
					sum_scale += LOG_SCALING_THRESHOLD* 4 * ptn_freq[ptn];
					//sum_scale += log(lh_max) * ptn_freq[ptn];
					dad_branch->scale_num[ptn] += 4;
					int nsite = aln->getNSite();
					for (i = 0, x = 0; i < nsite && x < ptn_freq[ptn]; i++)
						if (aln->getPatternID(i) == ptn) {
							outWarning((string)"Numerical underflow for site " + convertIntToString(i+1));
							x++;
						}
            	} else {
					// now do the likelihood scaling
					for (i = 0; i < block; i++) {
						partial_lh[i] *= SCALING_THRESHOLD_INVER;
	                    //partial_lh[i] /= lh_max;
					}
					// unobserved const pattern will never have underflow
					sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
					//sum_scale += log(lh_max) * ptn_freq[ptn];
					dad_branch->scale_num[ptn] += 1;
            	}
            }

		}
		dad_branch->lh_scale_factor += sum_scale;
        
	} else {
        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/

#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			double partial_lh_tmp[nstates];
			double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = left->partial_lh + ptn*block;
			double *partial_lh_right = right->partial_lh + ptn*block;
            double lh_max = 0.0;

            double expleft[nstates];
            double expright[nstates];
            double *eval = models->at(ptn)->getEigenvalues();
            double *evec = models->at(ptn)->getEigenvectors();
            double *inv_evec = models->at(ptn)->getInverseEigenvectors();

			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];

			for (c = 0; c < ncat; c++) {
                double len_left = site_rate->getRate(c) * left->length;
                double len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
                    double *this_evec = evec + x*nstates;
					for (i = 0; i < nstates; i++) {
						vleft += this_evec[i] * expleft[i];
						vright += this_evec[i] * expright[i];
					}
					partial_lh_tmp[x] = vleft*vright;
				}
                partial_lh_left += nstates;
                partial_lh_right += nstates;
                
				// compute dot-product with inv_eigenvector
                double *inv_evec_ptr = inv_evec;
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++) {
						res += partial_lh_tmp[x]*inv_evec_ptr[x];
					}
                    inv_evec_ptr += nstates;
					partial_lh[c*nstates+i] = res;
                    lh_max = max(lh_max, fabs(res));
				}
			}

            // check if one should scale partial likelihoods
            if (lh_max < SCALING_THRESHOLD) {
            	if (lh_max == 0.0) {
            		// for very shitty data
            		for (c = 0; c < ncat; c++)
            			memcpy(&partial_lh[c*nstates], &tip_partial_lh[ptn*nstates], nstates*sizeof(double));
					sum_scale += LOG_SCALING_THRESHOLD* 4 * ptn_freq[ptn];
					//sum_scale += log(lh_max) * ptn_freq[ptn];
					dad_branch->scale_num[ptn] += 4;
					int nsite = aln->getNSite();
					for (i = 0, x = 0; i < nsite && x < ptn_freq[ptn]; i++)
						if (aln->getPatternID(i) == ptn) {
							outWarning((string)"Numerical underflow for site " + convertIntToString(i+1));
							x++;
						}
            	} else {
					// now do the likelihood scaling
					for (i = 0; i < block; i++) {
						partial_lh[i] *= SCALING_THRESHOLD_INVER;
	                    //partial_lh[i] /= lh_max;
					}
					// unobserved const pattern will never have underflow
					sum_scale += LOG_SCALING_THRESHOLD * ptn_freq[ptn];
					//sum_scale += log(lh_max) * ptn_freq[ptn];
					dad_branch->scale_num[ptn] += 1;
            	}
            }

		}
		dad_branch->lh_scale_factor += sum_scale;

	}

}

//template <const int nstates>
void PhyloTree::computeSitemodelLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);
        
    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();

	ASSERT(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
            
            double *tip_partial_lh_node = tip_partial_lh + (dad->id * get_safe_upper_limit(nptn)*nstates);
            
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
				double *theta = theta_all + ptn*block;
				double *lh_tip = tip_partial_lh_node + ptn*nstates;
                for (c = 0; c < ncat; c++) {
                    for (i = 0; i < nstates; i++) {
                        theta[i] = lh_tip[i] * partial_lh_dad[i];
                    }
                    partial_lh_dad += nstates;
                    theta += nstates;
                }

			}
	    } else 
        {
	    	// both dad and node are internal nodes

//	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				double *theta = theta_all + ptn*block;
			    double *partial_lh_node = node_branch->partial_lh + ptn*block;
			    double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
	    		for (i = 0; i < block; i++) {
	    			theta[i] = partial_lh_node[i] * partial_lh_dad[i];
	    		}
			}
	    }
		theta_computed = true;
	}

    ModelSet *models = (ModelSet*)model;
    double my_df = 0.0, my_ddf = 0.0;
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf) private(ptn, i, c) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
		double lh_ptn = ptn_invar[ptn], df_ptn = 0.0, ddf_ptn = 0.0;
		double *theta = theta_all + ptn*block;
        
        double *eval = models->at(ptn)->getEigenvalues();
        
        for (c = 0; c < ncat; c++) {
            double lh_cat = 0.0, df_cat = 0.0, ddf_cat = 0.0;
            for (i = 0; i < nstates; i++) {
                double cof = eval[i]*site_rate->getRate(c);
                double val = exp(cof*dad_branch->length) * theta[i];
                double val1 = cof*val;
                lh_cat += val;
                df_cat += val1;
                ddf_cat += cof*val1;
            }
            double prop = site_rate->getProp(c);
            lh_ptn += prop * lh_cat;
            df_ptn += prop * df_cat;
            ddf_ptn += prop * ddf_cat;
            theta += nstates;
        }

        lh_ptn = 1.0/fabs(lh_ptn);
        
        double df_frac = df_ptn * lh_ptn;
        double ddf_frac = ddf_ptn * lh_ptn;
        double freq = ptn_freq[ptn];
        double tmp1 = df_frac * freq;
        double tmp2 = ddf_frac * freq;
        my_df += tmp1;
        my_ddf += tmp2 - tmp1 * df_frac;
    }
	df = my_df;
	ddf = my_ddf;
    if (isnan(df) || isinf(df)) {
        df = 0.0;
        ddf = 0.0;
//        outWarning("Numerical instability (some site-likelihood = 0)");
    }

}

//template <const int nstates>
double PhyloTree::computeSitemodelLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad) {
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
//        computePartialLikelihoodEigen(dad_branch, dad);
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihoodEigen(node_branch, node);
        computePartialLikelihood(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();


	memset(_pattern_lh_cat, 0, sizeof(double)*nptn*ncat);
    ModelSet *models = (ModelSet*)model;

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
        double *tip_partial_lh_node = tip_partial_lh + (dad->id * get_safe_upper_limit(nptn)*nstates);
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = ptn_invar[ptn];
			double *lh_cat = _pattern_lh_cat + ptn*ncat;
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
			double *partial_lh_node = tip_partial_lh_node + ptn*nstates;
            double *eval = models->at(ptn)->getEigenvalues();

			for (c = 0; c < ncat; c++) {
                double len = site_rate->getRate(c)*dad_branch->length;
                double prop = site_rate->getProp(c);
				for (i = 0; i < nstates; i++) {
					*lh_cat +=  (exp(eval[i]*len) * partial_lh_node[i] * partial_lh_dad[i]);
				}
                *lh_cat *= prop;
				lh_ptn += *lh_cat;
                // don't increase partial_lh_node pointer
				partial_lh_dad += nstates;
				lh_cat++;
			}

            lh_ptn = log(fabs(lh_ptn));
            _pattern_lh[ptn] = lh_ptn;
            tree_lh += lh_ptn * ptn_freq[ptn];
		}
    } else 
    {
    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = ptn_invar[ptn];
			double *lh_cat = _pattern_lh_cat + ptn*ncat;
			double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
			double *partial_lh_node = node_branch->partial_lh + ptn*block;
            double *eval = models->at(ptn)->getEigenvalues();

			for (c = 0; c < ncat; c++) {
                double len = site_rate->getRate(c)*dad_branch->length;
                double prop = site_rate->getProp(c);
				for (i = 0; i < nstates; i++) {
					*lh_cat +=  (exp(eval[i]*len) * partial_lh_node[i] * partial_lh_dad[i]);
				}
                *lh_cat *= prop;
				lh_ptn += *lh_cat;
				partial_lh_node += nstates;
				partial_lh_dad += nstates;
				lh_cat++;
			}

            lh_ptn = log(fabs(lh_ptn));
            _pattern_lh[ptn] = lh_ptn;
            tree_lh += lh_ptn * ptn_freq[ptn];
		}
    }

    if (isnan(tree_lh) || isinf(tree_lh)) {
        cout << "WARNING: Numerical underflow caused by alignment sites";
        i = aln->getNSite();
        int j;
        for (j = 0, c = 0; j < i; j++) {
            ptn = aln->getPatternID(j);
            if (isnan(_pattern_lh[ptn]) || isinf(_pattern_lh[ptn])) {
                cout << " " << j+1;
                c++;
                if (c >= 10) {
                    cout << " ...";
                    break;
                }
            }
        }
        cout << endl;
        tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;
        for (ptn = 0; ptn < nptn; ptn++) {
            if (isnan(_pattern_lh[ptn]) || isinf(_pattern_lh[ptn])) {
                _pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += _pattern_lh[ptn] * ptn_freq[ptn];
        }
    }

	ASSERT(!isnan(tree_lh) && !isinf(tree_lh));

    return tree_lh;
}


double PhyloTree::computeSitemodelLikelihoodFromBufferEigen() {
	ASSERT(theta_all && theta_computed);

    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();

    ModelSet *models = (ModelSet*)model;
    
    double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
		double lh_ptn = ptn_invar[ptn];
		double *theta = theta_all + ptn*block;
        
        double *eval = models->at(ptn)->getEigenvalues();
        
        for (c = 0; c < ncat; c++) {
            double lh_cat = 0.0;
            double len = site_rate->getRate(c)*current_it->length;
            for (i = 0; i < nstates; i++) {
                lh_cat += exp(eval[i]*len) * theta[i];
            }
            lh_ptn += lh_cat * site_rate->getProp(c);
            theta += nstates;
        }

        lh_ptn = log(fabs(lh_ptn));
        _pattern_lh[ptn] = lh_ptn;
        tree_lh += lh_ptn * ptn_freq[ptn];
    }
    return tree_lh;
}

