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

#endif /* PHYLOKERNELSITEMODEL_H_ */
