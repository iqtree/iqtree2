/*
 * phylokernelsitemodel.h
 * optimize SIMD likelihood function for site-specific model
 *  Created on: Jan 12, 2016
 *      Author: minh
 */

#ifndef PHYLOKERNELSITEMODEL_H_
#define PHYLOKERNELSITEMODEL_H_

#include "phylotree.h"
#include "model/modelset.h"

inline double horizontal_add(double x) {
    return x;
}

template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeSitemodelPartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {

    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);

    size_t nptn = aln->size(), tip_block_size = get_safe_upper_limit(nptn)*nstates;
    size_t ptn, c;
    size_t ncat = site_rate->getNRate();
    size_t i, x, j;
    size_t block = nstates * ncat;
    ModelSet *models = (ModelSet*) model;
    assert(models->size() == nptn);


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
            computeSitemodelPartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(nei, node);
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
        assert(done && "partial_lh is not re-oriented");
    }


    double sum_scale = 0.0;
        
	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
	}
    
    assert(node->degree() == 3); // does not work with multifurcating tree yet
    
    if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *tip_partial_lh_left = tip_partial_lh + (left->node->id * tip_block_size);
        double *tip_partial_lh_right = tip_partial_lh + (right->node->id * tip_block_size);

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i, j) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			VectorClass partial_lh_tmp[nstates/VCSIZE];
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(tip_partial_lh_left + ptn*nstates);
			VectorClass *partial_lh_right = (VectorClass*)(tip_partial_lh_right + ptn*nstates);

            VectorClass expleft[nstates/VCSIZE];
            VectorClass expright[nstates/VCSIZE];
            VectorClass *eval = (VectorClass*)(models->at(ptn)->getEigenvalues());
            VectorClass *evec = (VectorClass*)(models->at(ptn)->getEigenvectors());
            VectorClass *inv_evec = (VectorClass*)(models->at(ptn)->getInverseEigenvectors());
            VectorClass vleft[VCSIZE];
            VectorClass vright[VCSIZE];
            VectorClass res[VCSIZE];
            VectorClass len_left, len_right;
            VectorClass *this_evec;

			for (c = 0; c < ncat; c++) {
                len_left = site_rate->getRate(c) * left->length;
                len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates/VCSIZE; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }
				// compute real partial likelihood vector
                this_evec = evec;
				for (x = 0; x < nstates/VCSIZE; x++) {
                    for (j = 0; j < VCSIZE; j++) {
                        vleft[j] = 0.0;
                        vright[j] = 0.0;
                        for (i = 0; i < nstates/VCSIZE; i++) {
                            vleft[j] += this_evec[i] * expleft[i];
                            vright[j] += this_evec[i] * expright[i];
                        }
                        this_evec += nstates/VCSIZE;
                    }
					partial_lh_tmp[x] = horizontal_add(vleft)*horizontal_add(vright);
				}
                
				// compute dot-product with inv_eigenvector
                this_evec = inv_evec;
				for (i = 0; i < nstates/VCSIZE; i++) {
                    for (j = 0; j < VCSIZE; j++) {
                        res[j] = 0.0;
                        for (x = 0; x < nstates/VCSIZE; x++) {
                            res[j] += partial_lh_tmp[x]*this_evec[x];
                        }
                        this_evec += nstates/VCSIZE;
                    }
                    partial_lh[i] = horizontal_add(res);
				}

//                partial_lh_left += nstates/VCSIZE;
//                partial_lh_right += nstates/VCSIZE;
                partial_lh += nstates/VCSIZE;
			}
		}

	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

        double *tip_partial_lh_left = tip_partial_lh + (left->node->id * tip_block_size);

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i, j) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			VectorClass partial_lh_tmp[nstates/VCSIZE];
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(tip_partial_lh_left + ptn*nstates);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;

            VectorClass expleft[nstates/VCSIZE];
            VectorClass expright[nstates/VCSIZE];
            VectorClass *eval = (VectorClass*)(models->at(ptn)->getEigenvalues());
            VectorClass *evec = (VectorClass*)(models->at(ptn)->getEigenvectors());
            VectorClass *inv_evec = (VectorClass*)(models->at(ptn)->getInverseEigenvectors());
            VectorClass vleft[VCSIZE];
            VectorClass vright[VCSIZE];
            VectorClass res[VCSIZE];
            VectorClass len_left, len_right;
            VectorClass *this_evec;
            
			for (c = 0; c < ncat; c++) {
                len_left = site_rate->getRate(c) * left->length;
                len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates/VCSIZE; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }
				// compute real partial likelihood vector
                this_evec = evec;
				for (x = 0; x < nstates/VCSIZE; x++) {
                    for (j = 0; j < VCSIZE; j++) {
                        vleft[j] = 0.0;
                        vright[j] = 0.0;
                        for (i = 0; i < nstates/VCSIZE; i++) {
                            vleft[j] += this_evec[i] * expleft[i];
                            vright[j] += this_evec[i] * expright[i];
                        }
                        this_evec += nstates/VCSIZE;
                    }
					partial_lh_tmp[x] = horizontal_add(vleft)*horizontal_add(vright);
				}
                
				// compute dot-product with inv_eigenvector
                this_evec = inv_evec;
				for (i = 0; i < nstates/VCSIZE; i++) {
                    for (j = 0; j < VCSIZE; j++) {
                        res[j] = 0.0;
                        for (x = 0; x < nstates/VCSIZE; x++) {
                            res[j] += partial_lh_tmp[x]*this_evec[x];
                        }
                        this_evec += nstates/VCSIZE;
                    }
                    lh_max = max(lh_max, abs(partial_lh[i] = horizontal_add(res)));
				}

//                partial_lh_left += nstates/VCSIZE;
                partial_lh_right += nstates/VCSIZE;
                partial_lh += nstates/VCSIZE;
			}

            // check if one should scale partial likelihoods
            double dmax = horizontal_max(lh_max);
            if (dmax < SCALING_THRESHOLD) {
                partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            	if (dmax == 0.0) {
            		// for very shitty data
            		for (c = 0; c < ncat; c++)
            			memcpy(&partial_lh[c*nstates/VCSIZE], &tip_partial_lh[aln->STATE_UNKNOWN*nstates], nstates*sizeof(double));
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
					for (i = 0; i < block/VCSIZE; i++) {
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
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i, j) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
			VectorClass partial_lh_tmp[nstates/VCSIZE];
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(left->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;

            VectorClass expleft[nstates/VCSIZE];
            VectorClass expright[nstates/VCSIZE];
            VectorClass *eval = (VectorClass*)(models->at(ptn)->getEigenvalues());
            VectorClass *evec = (VectorClass*)(models->at(ptn)->getEigenvectors());
            VectorClass *inv_evec = (VectorClass*)(models->at(ptn)->getInverseEigenvectors());
            VectorClass vleft[VCSIZE];
            VectorClass vright[VCSIZE];
            VectorClass res[VCSIZE];
            VectorClass len_left, len_right;
            VectorClass *this_evec;
            
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];

			for (c = 0; c < ncat; c++) {
                len_left = site_rate->getRate(c) * left->length;
                len_right = site_rate->getRate(c) * right->length;
                for (i = 0; i < nstates/VCSIZE; i++) {
                    expleft[i] = exp(eval[i]*len_left) * partial_lh_left[i];
                    expright[i] = exp(eval[i]*len_right) * partial_lh_right[i];
                }
				// compute real partial likelihood vector
                this_evec = evec;
				for (x = 0; x < nstates/VCSIZE; x++) {
                    for (j = 0; j < VCSIZE; j++) {
                        vleft[j] = 0.0;
                        vright[j] = 0.0;
//                        this_evec = evec + (x*VCSIZE+j)*nstates/VCSIZE;
                        for (i = 0; i < nstates/VCSIZE; i++) {
                            vleft[j] += this_evec[i] * expleft[i];
                            vright[j] += this_evec[i] * expright[i];
                        }
                        this_evec += nstates/VCSIZE;
                    }
					partial_lh_tmp[x] = horizontal_add(vleft)*horizontal_add(vright);
				}
                
				// compute dot-product with inv_eigenvector
                this_evec = inv_evec;
				for (i = 0; i < nstates/VCSIZE; i++) {
                    for (j = 0; j < VCSIZE; j++) {
                        res[j] = 0.0;
                        for (x = 0; x < nstates/VCSIZE; x++) {
                            res[j] += partial_lh_tmp[x]*this_evec[x];
                        }
                        this_evec += nstates/VCSIZE;
                    }
                    lh_max = max(lh_max, abs(partial_lh[i] = horizontal_add(res)));
				}

                partial_lh_left += nstates/VCSIZE;
                partial_lh_right += nstates/VCSIZE;
                partial_lh += nstates/VCSIZE;
			}

            // check if one should scale partial likelihoods
            double dmax = horizontal_max(lh_max);
            if (dmax < SCALING_THRESHOLD) {
                partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            	if (dmax == 0.0) {
            		// for very shitty data
            		for (c = 0; c < ncat; c++)
            			memcpy(&partial_lh[c*nstates/VCSIZE], &tip_partial_lh[aln->STATE_UNKNOWN*nstates], nstates*sizeof(double));
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
					for (i = 0; i < block/VCSIZE; i++) {
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

template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeSitemodelLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
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
        computeSitemodelPartialLikelihoodEigenSIMD<VectorClass,VCSIZE,nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeSitemodelPartialLikelihoodEigenSIMD<VectorClass,VCSIZE,nstates>(node_branch, node);
        
//    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t nptn = aln->size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
            
            double *tip_partial_lh_node = tip_partial_lh + (dad->id * get_safe_upper_limit(nptn)*nstates);
            
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
				VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
				VectorClass *lh_tip = (VectorClass*)(tip_partial_lh_node + ptn*nstates);
                for (c = 0; c < ncat; c++) {
                    for (i = 0; i < nstates/VCSIZE; i++) {
                        theta[i] = lh_tip[i] * partial_lh_dad[i];
                    }
                    partial_lh_dad += nstates/VCSIZE;
                    theta += nstates/VCSIZE;
                }

			}
	    } else 
        {
	    	// both dad and node are internal nodes

            size_t block_VCSIZE = block/VCSIZE;

//	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
				VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
			    VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
	    		for (i = 0; i < block_VCSIZE; i++) {
	    			theta[i] = partial_lh_node[i] * partial_lh_dad[i];
	    		}
			}
	    }
		if (nptn < maxptn) {
			// copy dummy values
			for (ptn = nptn; ptn < maxptn; ptn++)
				memcpy(&theta_all[ptn*block], theta_all, block*sizeof(double));
		}
		theta_computed = true;
	}

    ModelSet *models = (ModelSet*)model;
    VectorClass my_df = 0.0, my_ddf = 0.0;
    VectorClass dad_length = dad_branch->length;
    VectorClass unit = 1.0;
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf) private(ptn, i, c, j) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VCSIZE) {
        VectorClass lh_ptn[VCSIZE];
        VectorClass df_ptn[VCSIZE];
        VectorClass ddf_ptn[VCSIZE];
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);        
        VectorClass* eval;
        
        for (j = 0; j < VCSIZE; j++) {
            lh_ptn[j] = 0.0;
            df_ptn[j] = 0.0;
            ddf_ptn[j] = 0.0;
            if (ptn+j < nptn) {
                eval = (VectorClass*)models->at(ptn+j)->getEigenvalues();
            } else {
                eval = (VectorClass*)models->at(nptn-1)->getEigenvalues();
            }                
            for (c = 0; c < ncat; c++) {
                VectorClass cat_rate = site_rate->getRate(c);
                VectorClass lh_cat = 0.0, df_cat = 0.0, ddf_cat = 0.0;
                for (i = 0; i < nstates/VCSIZE; i++) {
                    VectorClass cof = eval[i]*cat_rate;
                    VectorClass val = exp(cof*dad_length) * theta[i];
                    VectorClass val1 = cof*val;
                    lh_cat += val;
                    df_cat += val1;
                    ddf_cat += cof*val1;
                }
                VectorClass prop = site_rate->getProp(c);
                lh_ptn[j] += prop * lh_cat;
                df_ptn[j] += prop * df_cat;
                ddf_ptn[j] += prop * ddf_cat;
                theta += nstates/VCSIZE;
            }
        }

        VectorClass inv_lh_ptn = horizontal_add(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
        inv_lh_ptn = unit / abs(inv_lh_ptn);
        VectorClass freq;
        freq.load_a(&ptn_freq[ptn]);
        
        VectorClass df_ptn_sum = horizontal_add(df_ptn) * inv_lh_ptn;
        VectorClass ddf_ptn_sum = horizontal_add(ddf_ptn) * inv_lh_ptn;
        ddf_ptn_sum = ddf_ptn_sum - df_ptn_sum*df_ptn_sum;
        
        my_df += df_ptn_sum * freq;
        my_ddf += ddf_ptn_sum * freq;
    }
	df = horizontal_add(my_df);
	ddf = horizontal_add(my_ddf);
    if (isnan(df) || isinf(df)) {
        df = 0.0;
        ddf = 0.0;
        outWarning("Numerical instability (some site-likelihood = 0)");
    }

}

template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeSitemodelLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
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
        computeSitemodelPartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computeSitemodelPartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node);
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t nptn = aln->size();
    size_t maxptn = get_safe_upper_limit(nptn);

    ModelSet *models = (ModelSet*)model;
    VectorClass tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    VectorClass dad_length = dad_branch->length;

    if (dad->isLeaf()) {
		// copy dummy values because VectorClass will access beyond nptn
		for (ptn = nptn; ptn < maxptn; ptn++)
			memcpy(&dad_branch->partial_lh[ptn*block], &dad_branch->partial_lh[(ptn-1)*block], block*sizeof(double));

    	// special treatment for TIP-INTERNAL NODE case
        double *tip_partial_lh_node = tip_partial_lh + (dad->id * get_safe_upper_limit(nptn)*nstates);
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c, j) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn+=VCSIZE) {
            VectorClass lh_ptn[VCSIZE];
            VectorClass* eval;
			VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_node = (VectorClass*)(tip_partial_lh_node + ptn*nstates);
            
            for (j = 0; j < VCSIZE; j++) {
                lh_ptn[j] = 0.0;
                if (ptn+j < nptn) {
                    eval = (VectorClass*)models->at(ptn+j)->getEigenvalues();
                } else {
                    eval = (VectorClass*)models->at(nptn-1)->getEigenvalues();
                }                
                for (c = 0; c < ncat; c++) {
                    VectorClass cat_rate = site_rate->getRate(c);
                    VectorClass lh_cat = 0.0;
                    for (i = 0; i < nstates/VCSIZE; i++) {
                        VectorClass cof = eval[i]*cat_rate;
                        VectorClass val = exp(cof*dad_length) * partial_lh_dad[i] * partial_lh_node[i];
                        lh_cat += val;
                    }
                    VectorClass prop = site_rate->getProp(c);
                    lh_ptn[j] += prop * lh_cat;
                    partial_lh_dad += nstates/VCSIZE;
//                    partial_lh_node += nstates/VCSIZE;
                }
                partial_lh_node += nstates/VCSIZE;
            }

            VectorClass freq;
            freq.load_a(&ptn_freq[ptn]);
            VectorClass lh_ptn_sum = horizontal_add(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
            lh_ptn_sum = log(abs(lh_ptn_sum));
            lh_ptn_sum.store_a(&_pattern_lh[ptn]);
            tree_lh += lh_ptn_sum * freq;
        }
    } else 
    {
    	// both dad and node are internal nodes
		// copy dummy values because VectorClass will access beyond nptn
		for (ptn = nptn; ptn < maxptn; ptn++) {
			memcpy(&dad_branch->partial_lh[ptn*block], &dad_branch->partial_lh[(ptn-1)*block], block*sizeof(double));
			memcpy(&node_branch->partial_lh[ptn*block], &node_branch->partial_lh[(ptn-1)*block], block*sizeof(double));
        }
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c, j) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn+=VCSIZE) {
            VectorClass lh_ptn[VCSIZE];
            VectorClass* eval;
			VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
            
            for (j = 0; j < VCSIZE; j++) {
                lh_ptn[j] = 0.0;
                if (ptn+j < nptn) {
                    eval = (VectorClass*)models->at(ptn+j)->getEigenvalues();
                } else {
                    eval = (VectorClass*)models->at(nptn-1)->getEigenvalues();
                }                
                for (c = 0; c < ncat; c++) {
                    VectorClass cat_rate = site_rate->getRate(c);
                    VectorClass lh_cat = 0.0;
                    for (i = 0; i < nstates/VCSIZE; i++) {
                        VectorClass cof = eval[i]*cat_rate;
                        VectorClass val = exp(cof*dad_length) * partial_lh_dad[i] * partial_lh_node[i];
                        lh_cat += val;
                    }
                    VectorClass prop = site_rate->getProp(c);
                    lh_ptn[j] += prop * lh_cat;
                    partial_lh_dad += nstates/VCSIZE;
                    partial_lh_node += nstates/VCSIZE;
                }
            }

            VectorClass freq;
            freq.load_a(&ptn_freq[ptn]);
            VectorClass lh_ptn_sum = horizontal_add(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
            lh_ptn_sum = log(abs(lh_ptn_sum));
            lh_ptn_sum.store_a(&_pattern_lh[ptn]);
            tree_lh += lh_ptn_sum * freq;
        }

    }

    double tree_lh_final = horizontal_add(tree_lh);

    if (isnan(tree_lh_final) || isinf(tree_lh_final)) {
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

	assert(!isnan(tree_lh_final) && !isinf(tree_lh_final));

    return tree_lh_final;
}


template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeSitemodelLikelihoodFromBufferEigenSIMD() {
	assert(theta_all && theta_computed);

//    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();

    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t nptn = aln->size();

    ModelSet *models = (ModelSet*)model;
    
    VectorClass dad_length = current_it->length;
    VectorClass tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) private(ptn, i, c, j) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VCSIZE) {
        VectorClass lh_ptn[VCSIZE];
        VectorClass* eval;
        VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
        
        for (j = 0; j < VCSIZE; j++) {
            lh_ptn[j] = 0.0;
            if (ptn+j < nptn) {
                eval = (VectorClass*)models->at(ptn+j)->getEigenvalues();
            } else {
                eval = (VectorClass*)models->at(nptn-1)->getEigenvalues();
            }                
            for (c = 0; c < ncat; c++) {
                VectorClass cat_rate = site_rate->getRate(c);
                VectorClass lh_cat = 0.0;
                for (i = 0; i < nstates/VCSIZE; i++) {
                    VectorClass cof = eval[i]*cat_rate;
                    VectorClass val = exp(cof*dad_length) * theta[i];
                    lh_cat += val;
                }
                VectorClass prop = site_rate->getProp(c);
                lh_ptn[j] += prop * lh_cat;
                theta += nstates/VCSIZE;
            }
        }

        VectorClass freq;
        freq.load_a(&ptn_freq[ptn]);
        VectorClass lh_ptn_sum = horizontal_add(lh_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
        lh_ptn_sum = log(abs(lh_ptn_sum));
        lh_ptn_sum.store_a(&_pattern_lh[ptn]);
        tree_lh += lh_ptn_sum * freq;
    }
    double tree_lh_final = horizontal_add(tree_lh);
    return tree_lh_final;
}


#endif /* PHYLOKERNELSITEMODEL_H_ */
