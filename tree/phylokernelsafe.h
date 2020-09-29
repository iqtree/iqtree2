/*
 * phylokernelsafe.h
 * Safe likelihood kernel that scales likelihood per category
 *
 *  Created on: Sept 23, 2016
 *      Author: minh
 */

#ifndef PHYLOKERNELSAFE_H_
#define PHYLOKERNELSAFE_H_

#include "phylotree.h"
//#include "vectorclass/vectorclass.h"
//#include "vectorclass/vectormath_exp.h"
#include "superalignment.h"


#ifdef __SSE__
inline Vec2d horizontal_add(Vec2d x[2]) {
#if  INSTRSET >= 3  // SSE3
    return _mm_hadd_pd(x[0],x[1]);
#elif INSTRSET >= 2
    Vec2d help0 = _mm_shuffle_pd(x[0], x[1], _MM_SHUFFLE2(0,0));
    Vec2d help1 = _mm_shuffle_pd(x[0], x[1], _MM_SHUFFLE2(1,1));
    return _mm_add_pd(help0, help1);
#else
#error "You must compile with SSE2 enabled!"
#endif
}

inline double horizontal_max(Vec2d const &a) {
    double x[2];
    a.store(x);
    return max(x[0],x[1]);
}

#endif

#ifdef __AVX__

inline Vec4d horizontal_add(Vec4d x[4]) {
	// {a[0]+a[1], b[0]+b[1], a[2]+a[3], b[2]+b[3]}
	__m256d sumab = _mm256_hadd_pd(x[0], x[1]);
	// {c[0]+c[1], d[0]+d[1], c[2]+c[3], d[2]+d[3]}
	__m256d sumcd = _mm256_hadd_pd(x[2], x[3]);

	// {a[0]+a[1], b[0]+b[1], c[2]+c[3], d[2]+d[3]}
	__m256d blend = _mm256_blend_pd(sumab, sumcd, 12/* 0b1100*/);
	// {a[2]+a[3], b[2]+b[3], c[0]+c[1], d[0]+d[1]}
	__m256d perm = _mm256_permute2f128_pd(sumab, sumcd, 0x21);

	return _mm256_add_pd(perm, blend);
}

inline double horizontal_max(Vec4d const &a) {
	__m128d high = _mm256_extractf128_pd(a,1);
	__m128d m = _mm_max_pd(_mm256_castpd256_pd128(a), high);
    double x[2];
    _mm_storeu_pd(x, m);
    return max(x[0],x[1]);
}

#endif // __AVX__

template <class Numeric, class VectorClass, const int VCSIZE>
Numeric PhyloTree::dotProductSIMD(Numeric *x, Numeric *y, int size) {
	VectorClass res = VectorClass().load_a(x) * VectorClass().load_a(y);
	for (int i = VCSIZE; i < size; i += VCSIZE)
		res = mul_add(VectorClass().load_a(&x[i]), VectorClass().load_a(&y[i]), res);
	return horizontal_add(res);
}

/************************************************************************************************
 *
 *   Highly optimized vectorized versions of likelihood functions
 *
 *************************************************************************************************/


template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {

    // don't recompute the likelihood
	ASSERT(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;

    num_partial_lh_computations++;

    size_t nptn = aln->size() + model_factory->unobserved_ptns.size();
    PhyloNode *node = dad_branch->getNode();

    if (!tip_partial_lh_computed)
        computeTipPartialLikelihood();

	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;
	    //memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		return;
	}

    size_t ptn, c;
    size_t orig_nptn = aln->size();

    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    ASSERT(nstates == aln->num_states && nstates >= VCSIZE && VCSIZE == VectorClass().size());
    ASSERT(model->isReversible()); // only works with reversible model!
    const size_t nstatesqr=nstates*nstates;
    size_t i, x, j;
    size_t block = nstates * ncat_mix;
    size_t tip_block = nstates * model->getNMixtures();
    size_t scale_size = nptn * ncat_mix;

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = m*nstatesqr;
    }

	// internal node
    dad_branch->lh_scale_factor = 0.0;
	PhyloNeighbor *left = nullptr, *right = nullptr; // left & right are two neighbors leading to 2 subtrees
    int num_leaves = 0;
	FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
		if (!left) left = nei; else right = nei;
        if ((nei->partial_lh_computed & 1) == 0)
            computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(nei, node);
        dad_branch->lh_scale_factor += nei->lh_scale_factor;
        if ((*it)->node->isLeaf()) num_leaves++;
	}

    if (params->lh_mem_save == LM_PER_NODE && !dad_branch->partial_lh) {
        // re-orient partial_lh
        bool done = false;
        FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it2, child) {
            PhyloNeighbor *backnei = child->findNeighbor(node);
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

	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();

	ASSERT(inv_evec && evec);
//	for (i = 0; i < tip_block; i++) {
//		for (x = 0; x < nstates/VCSIZE; x++)
//			// inv_evec is not aligned!
//			vc_inv_evec[i*nstates/VCSIZE+x].load_a(&inv_evec[i*nstates+x*VCSIZE]);
//	}
	double *eval = model->getEigenvalues();


    VectorClass *echildren = aligned_alloc<VectorClass>(block*nstates/VCSIZE*(node->degree()-1));
    double *partial_lh_leaves = NULL;
    if (num_leaves > 0)
        partial_lh_leaves = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block*num_leaves);
    VectorClass *echild = echildren;
    double *partial_lh_leaf = partial_lh_leaves;
    
    
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        VectorClass expchild[nstates/VCSIZE];
        VectorClass *echild_ptr = echild;
        // precompute information buffer
        for (c = 0; c < ncat_mix; c++) {
            VectorClass len_child = site_rate->getRate(c%ncat) * nei->length;
            double *eval_ptr = eval + mix_addr_nstates[c];
            double *evec_ptr = evec + mix_addr[c];
            for (i = 0; i < nstates/VCSIZE; i++) {
                // eval is not aligned!
                expchild[i] = exp(VectorClass().load_a(&eval_ptr[i*VCSIZE]) * len_child);
            }
            for (x = 0; x < nstates; x++) {
                for (i = 0; i < nstates/VCSIZE; i++) {
                    // evec is not be aligned!
                    echild_ptr[i] = (VectorClass().load_a(&evec_ptr[x*nstates+i*VCSIZE]) * expchild[i]);
                }
                echild_ptr += nstates/VCSIZE;
            }
        }

        // pre compute information for tip
        PhyloNode* child = nei->getNode();
        if (child->node->isLeaf()) {
            vector<int>::iterator it;
            for (it = aln->seq_states[child->id].begin()
                ; it != aln->seq_states[child->id].end(); it++) {
                int state = (*it);
                double *this_partial_lh_leaf = partial_lh_leaf + state*block;
                VectorClass *echild_ptr = echild;
                for (c = 0; c < ncat_mix; c++) {
                    VectorClass *this_tip_partial_lh = (VectorClass*)(tip_partial_lh + state*tip_block + mix_addr_nstates[c]);
                    for (x = 0; x < nstates; x++) {
                        VectorClass vchild = 0.0;
                        for (i = 0; i < nstates/VCSIZE; i++) {
                            vchild += echild_ptr[i] * this_tip_partial_lh[i];
                        }
                        this_partial_lh_leaf[x] = horizontal_add(vchild);
                        echild_ptr += nstates/VCSIZE;
                    }
                    this_partial_lh_leaf += nstates;
                }
            }
            size_t addr = aln->STATE_UNKNOWN * block;
            for (x = 0; x < block; x++) {
                partial_lh_leaf[addr+x] = 1.0;
            }
            partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
        }
        echild += block*nstates/VCSIZE;
    }
    
    VectorClass *eleft = echildren, *eright = echildren + block*nstates/VCSIZE;
    
	if (!left->node->isLeaf() && right->node->isLeaf()) {
        std::swap(left, right);
        std::swap(eleft, eright);
	}
    
    
    if (node->degree() > 3) {

        /*--------------------- multifurcating node ------------------*/
        // now for-loop computing partial_lh over all site-patterns
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn++) {
            double partial_lh_all[block];
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
            memset(scale_dad, 0, sizeof(UBYTE)*ncat_mix);

            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = (double*)echildren;

            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
                PhyloNode* child = nei->getNode();
                UBYTE *scale_child = nei->scale_num + ptn*ncat_mix;
                if (child->isLeaf()) {
                    // external node
                    int state_child = (ptn < orig_nptn) ? (aln->at(ptn))[child->id] : model_factory->unobserved_ptns[ptn-orig_nptn];
                    double *child_lh = partial_lh_leaf + state_child*block;
                    for (c = 0; c < block; c++) {
                        // compute real partial likelihood vector
                        partial_lh_all[c] *= child_lh[c];
                    }
                    partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                } else {
                    // internal node
                    double* partial_lh       = partial_lh_all;
                    double* partial_lh_child = nei->partial_lh + ptn*block;

                    double *echild_ptr = echild;
                    for (c = 0; c < ncat_mix; c++) {
                        scale_dad[c] += scale_child[c];
                        // compute real partial likelihood vector
                        for (x = 0; x < nstates; x++) {
                            double vchild = 0.0;
//                            double *echild_ptr = echild + (c*nstatesqr+x*nstates);
                            for (i = 0; i < nstates; i++) {
                                vchild += echild_ptr[i] * partial_lh_child[i];
                            }
                            echild_ptr += nstates;
                            partial_lh[x] *= vchild;
                        }
                        partial_lh += nstates;
                        partial_lh_child += nstates;
                    }
                } // if
                echild += block*nstates;
            } // FOR_EACH_PHYLO_NEIGHBOR
            
        
            // compute dot-product with inv_eigenvector
            double* partial_lh_tmp = partial_lh_all;
            double* partial_lh     = dad_branch->partial_lh + ptn*block;
            for (c = 0; c < ncat_mix; c++) {
                double lh_max = 0.0;
                double *inv_evec_ptr = inv_evec + mix_addr[c];
                for (i = 0; i < nstates; i++) {
                    double res = 0.0;
                    for (x = 0; x < nstates; x++) {
                        res += partial_lh_tmp[x]*inv_evec_ptr[x];
                    }
                    inv_evec_ptr += nstates;
                    partial_lh[i] = res;
                    lh_max = max(lh_max, fabs(res));
                }
                // check if one should scale partial likelihoods
                if (lh_max < SCALING_THRESHOLD && lh_max != 0.0) {
                    if (ptn_invar[ptn] == 0.0) {
                        // now do the likelihood scaling
                        for (i = 0; i < nstates; i++)
                            partial_lh[i] *= SCALING_THRESHOLD_INVER;
                        scale_dad[c] += 1;
                    }
                }
                partial_lh     += nstates;
                partial_lh_tmp += nstates;
            }

        } // for ptn

        // end multifurcating treatment
    } else if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case

		// pre compute information for both tips
		double *partial_lh_left = partial_lh_leaves;
		double *partial_lh_right = partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;

		// assign pointers for left and right partial_lh
        /*
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
        */

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, scale_size * sizeof(UBYTE));
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass res[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i, j, vc_partial_lh_tmp, res)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;

	        double *lh_left;
	        double *lh_right;
            if (ptn < orig_nptn) {
                lh_left = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
                lh_right = &partial_lh_right[block *  (aln->at(ptn))[right->node->id]];
            } else {
                lh_left = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_nptn]];
                lh_right = &partial_lh_right[block * model_factory->unobserved_ptns[ptn-orig_nptn]];
            }
			for (c = 0; c < ncat_mix; c++) {
                VectorClass *vc_inv_evec_ptr = (VectorClass*)(inv_evec + mix_addr[c]);
				// compute real partial likelihood vector

				for (x = 0; x < nstates/VCSIZE; x++) {
					vc_partial_lh_tmp[x] = (VectorClass().load_a(&lh_left[x*VCSIZE]) * VectorClass().load_a(&lh_right[x*VCSIZE]));
				}
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec_ptr[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++)
						for (j = 0; j < VCSIZE; j++) {
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec_ptr[(i+j)*nstates/VCSIZE+x], res[j]);
						}
					horizontal_add(res).store_a(&partial_lh[i]);
				}

				lh_left += nstates;
				lh_right += nstates;
				partial_lh += nstates;
			}
		}

	    //aligned_free(lh_right_ptr);
	    //aligned_free(lh_left_ptr);
	} else if (left->node->isLeaf() && !right->node->isLeaf()) {
		// special treatment to TIP-INTERNAL NODE case
		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, scale_size * sizeof(UBYTE));

		// pre compute information for left tip
		double *partial_lh_left = partial_lh_leaves;


		// assign pointers for partial_lh_left
        /*
		double **lh_left_ptr = aligned_alloc<double*>(nptn);
		for (ptn = 0; ptn < orig_ntn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
		}
		for (ptn = orig_ntn; ptn < nptn; ptn++) {
			lh_left_ptr[ptn] = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
		}
        */
		VectorClass vc_lh_right[nstates/VCSIZE];
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass res[VCSIZE];
		VectorClass vc_max; // maximum of partial likelihood, for scaling check
		VectorClass vright[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for private (ptn, c, x, i, j, vc_lh_right, vc_partial_lh_tmp, res, vc_max, vright)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;
	        double *partial_lh_right = right->partial_lh + ptn*block;

	        double *lh_left;
            if (ptn < orig_nptn) {
                lh_left = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
            } else {
                lh_left = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_nptn]];
            }

			for (c = 0; c < ncat_mix; c++) {
                vc_max = 0.0;
                VectorClass *vc_inv_evec_ptr = (VectorClass*)(inv_evec + mix_addr[c]);
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
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec_ptr[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++) {
						for (j = 0; j < VCSIZE; j++) {
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec_ptr[(i+j)*nstates/VCSIZE+x], res[j]);
						}
					}
					VectorClass sum_res = horizontal_add(res);
					sum_res.store_a(&partial_lh[i]);
					vc_max = max(vc_max, abs(sum_res)); // take the maximum for scaling check
				}
                // check if one should scale partial likelihoods
                double lh_max = horizontal_max(vc_max);
                if (lh_max < SCALING_THRESHOLD && ptn_invar[ptn] == 0.0 && lh_max != 0.0) {
                    // now do the likelihood scaling
                    VectorClass scale_thres(SCALING_THRESHOLD_INVER);
                    for (i = 0; i < block; i+=VCSIZE) {
                        (VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
                    }
                    dad_branch->scale_num[ptn*ncat_mix+c] += 1;
                }
				lh_left += nstates;
				partial_lh_right += nstates;
				partial_lh += nstates;
			}

		}

	} else {
		// both left and right are internal node

		VectorClass vc_max; // maximum of partial likelihood, for scaling check
		VectorClass vc_partial_lh_tmp[nstates/VCSIZE];
		VectorClass vc_lh_left[nstates/VCSIZE], vc_lh_right[nstates/VCSIZE];
		VectorClass res[VCSIZE];
		VectorClass vleft[VCSIZE], vright[VCSIZE];

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i, j, vc_max, vc_partial_lh_tmp, vc_lh_left, vc_lh_right, res, vleft, vright)
#endif
		for (ptn = 0; ptn < nptn; ptn++) {
	        double *partial_lh = dad_branch->partial_lh + ptn*block;
			double *partial_lh_left = left->partial_lh + ptn*block;
			double *partial_lh_right = right->partial_lh + ptn*block;
            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
            UBYTE *scale_left = left->scale_num + ptn*ncat_mix;
            UBYTE *scale_right = right->scale_num + ptn*ncat_mix; 

			for (c = 0; c < ncat_mix; c++) {
                scale_dad[c] = scale_left[c] + scale_right[c];
                vc_max = 0.0;
                VectorClass *vc_inv_evec_ptr = (VectorClass*)(inv_evec + mix_addr[c]);
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
						res[j] = vc_partial_lh_tmp[0] * vc_inv_evec_ptr[(i+j)*nstates/VCSIZE];
					}
					for (x = 1; x < nstates/VCSIZE; x++)
						for (j = 0; j < VCSIZE; j++)
							res[j] = mul_add(vc_partial_lh_tmp[x], vc_inv_evec_ptr[(i+j)*nstates/VCSIZE+x], res[j]);

					VectorClass sum_res = horizontal_add(res);
					sum_res.store_a(&partial_lh[i]);
					vc_max = max(vc_max, abs(sum_res)); // take the maximum for scaling check
				}
                // check if one should scale partial likelihoods
                double lh_max = horizontal_max(vc_max);
                if (lh_max < SCALING_THRESHOLD && ptn_invar[ptn] == 0.0 && lh_max != 0.0) {
                    // now do the likelihood scaling
                    VectorClass scale_thres(SCALING_THRESHOLD_INVER);
                    for (i = 0; i < block; i+=VCSIZE) {
                        (VectorClass().load_a(&partial_lh[i]) * scale_thres).store_a(&partial_lh[i]);
                    }
                    // unobserved const pattern will never have underflow
                    scale_dad[c] += 1;
                }
				partial_lh += nstates;
				partial_lh_left += nstates;
				partial_lh_right += nstates;
			}


		}

	}

	if (partial_lh_leaves)
        aligned_free(partial_lh_leaves);
	aligned_free(echildren);
}

template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeLikelihoodDervEigenSIMD
    ( PhyloNeighbor *dad_branch, PhyloNode *dad, 
      double &df, double &ddf) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node);
    df = ddf = 0.0;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    maxptn = max(maxptn, aln->size()+((model_factory->unobserved_ptns.size()+VCSIZE-1)/VCSIZE)*VCSIZE);

    size_t mix_addr_nstates[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    double *eval = model->getEigenvalues();
    ASSERT(eval);

	VectorClass *vc_val0 = (VectorClass*)aligned_alloc<double>(block);
	VectorClass *vc_val1 = (VectorClass*)aligned_alloc<double>(block);
	VectorClass *vc_val2 = (VectorClass*)aligned_alloc<double>(block);

	VectorClass vc_len = dad_branch->length;
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        size_t mycat = c%ncat;
        double *eval_ptr = eval + m*nstates;
		VectorClass vc_rate = site_rate->getRate(mycat);
		VectorClass vc_prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
		for (i = 0; i < nstates/VCSIZE; i++) {
			VectorClass cof = VectorClass().load_a(&eval_ptr[i*VCSIZE]) * vc_rate;
			VectorClass val = exp(cof*vc_len) * vc_prop;
			VectorClass val1_ = cof*val;
			vc_val0[c*nstates/VCSIZE+i] = val;
			vc_val1[c*nstates/VCSIZE+i] = val1_;
			vc_val2[c*nstates/VCSIZE+i] = cof*val1_;
		}
	}

	ASSERT(theta_all);
	if (!theta_computed) {
		theta_computed = true;
        double scale_all = 0.0;
		// precompute theta for fast branch length optimization

		if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) reduction(+: scale_all)
#endif
			for (ptn = 0; ptn < nptn; ptn++) {
			    double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
                UBYTE *scale_dad = dad_branch->scale_num+ptn*ncat_mix;
				double *theta = theta_all + ptn*block;
                double *this_tip_partial_lh = tip_partial_lh + tip_block*((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]);
                UBYTE min_scale = scale_dad[0];
                for (c = 1; c < ncat_mix; c++)
                    min_scale = min(min_scale, scale_dad[c]);

                scale_all += (double)min_scale;

                for (c = 0; c < ncat_mix; c++) {
                    double *lh_dad = this_tip_partial_lh + mix_addr_nstates[c];
                    if (scale_dad[c] == min_scale) {
                        for (i = 0; i < nstates; i+=VCSIZE) {
                            (VectorClass().load_a(&lh_dad[i]) * VectorClass().load_a(&partial_lh_dad[i])).store_a(&theta[i]);
                        }
                    } else if (scale_dad[c] == min_scale+1) {
                        for (i = 0; i < nstates; i+=VCSIZE) {
                            (VectorClass().load_a(&lh_dad[i]) * VectorClass().load_a(&partial_lh_dad[i]) * VectorClass(SCALING_THRESHOLD)).store_a(&theta[i]);
                        }
                    } else {
                        memset(theta, 0, sizeof(double)*nstates);
                    }
                    partial_lh_dad += nstates;
                    theta += nstates;
                }
			}
	    } else {
	    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for private(i, c) reduction(+: scale_all)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				double *theta = theta_all + ptn*block;
			    double *partial_lh_node = node_branch->partial_lh + ptn*block;
			    double *partial_lh_dad = dad_branch->partial_lh + ptn*block;

                size_t ptn_ncat = ptn*ncat_mix; 
                UBYTE *scale_dad = dad_branch->scale_num + ptn_ncat;
                UBYTE *scale_node = node_branch->scale_num + ptn_ncat;
                UBYTE sum_scale[ncat_mix];
                UBYTE min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
                for (c = 1; c < ncat_mix; c++) {
                    sum_scale[c] = scale_dad[c] + scale_node[c];
                    min_scale = min(min_scale, sum_scale[c]);
                }
                scale_all += (double)min_scale;

                for (c = 0; c < ncat_mix; c++) {
                    if (sum_scale[c] == min_scale) {
                        for (i = 0; i < nstates; i++) {
                            (VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i])).store_a(&theta[i]);
                        }
                    } else if (sum_scale[c] == min_scale+1) {
                        for (i = 0; i < nstates; i++) {
                            (VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i]) * VectorClass(SCALING_THRESHOLD)).store_a(&theta[i]);
                        }
                    } else {
                        memset(theta, 0, sizeof(double)*nstates);
                    }
                    theta += nstates;
                    partial_lh_dad += nstates;
                    partial_lh_node += nstates;
                }
			}
	    }
		if (nptn < maxptn) {
			// copy dummy values
			for (ptn = nptn; ptn < maxptn; ptn++)
				memcpy(&theta_all[ptn*block], theta_all, block*sizeof(double));
		}
        buffer_scale_all = scale_all*LOG_SCALING_THRESHOLD;
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

		inv_lh_ptn = vc_unit / abs(lh_ptn);

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
    
    ASSERT(!isnan(df) && !isinf(df) && "Numerical underflow for SIMD lh-derivative");

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
			df_ptn = horizontal_add(vc_df);
			ddf_ptn = horizontal_add(vc_ddf);

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
			ASSERT(0);
			break;
		}
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
	}
    ASSERT(!isnan(df));
    aligned_free(vc_val2);
    aligned_free(vc_val1);
    aligned_free(vc_val0);
}


template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh) {
        initializeAllPartialLh();
    }
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    size_t mix_addr_nstates[ncat_mix];

    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    maxptn = max(maxptn, aln->size()+((model_factory->unobserved_ptns.size()+VCSIZE-1)/VCSIZE)*VCSIZE);
    double *eval = model->getEigenvalues();
    ASSERT(eval);

    VectorClass *vc_val = (VectorClass*)aligned_alloc<double>(block);


	for (c = 0; c < ncat_mix; c++) {
        size_t mycat = c%ncat;
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        double *eval_ptr = eval + mix_addr_nstates[c];
		VectorClass vc_len(site_rate->getRate(mycat)*dad_branch->length);
		VectorClass vc_prop(site_rate->getProp(c) * model->getMixtureWeight(m));
		for (i = 0; i < nstates/VCSIZE; i++) {
			// eval is not aligned!
			vc_val[c*nstates/VCSIZE+i] = exp(VectorClass().load_a(&eval_ptr[i*VCSIZE]) * vc_len) * vc_prop;
		}
	}

	double prob_const = 0.0;

	if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case

    	// precompute information from one tip        
    	double *partial_lh_node = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);
    	IntVector states_dad = aln->seq_states[dad->id];
    	states_dad.push_back(aln->STATE_UNKNOWN);
    	for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
    		double *lh_node = partial_lh_node + (*it)*block;
    		double *lh_tip = tip_partial_lh + (*it)*tip_block;
    		VectorClass *vc_val_tmp = vc_val;
            for (c = 0; c < ncat_mix; c++) {
                double *this_lh_tip = lh_tip + mix_addr_nstates[c];
                for (i = 0; i < nstates; i+=VCSIZE) {
                    (vc_val_tmp[i/VCSIZE] * VectorClass().load_a(&this_lh_tip[i])).store_a(&lh_node[i]);
                }
                lh_node += nstates;
                vc_val_tmp += nstates/VCSIZE;
            }
    	}


    	//VectorClass vc_tip_partial_lh[nstates];
    	//VectorClass vc_partial_lh_dad[VCSIZE]
        VectorClass vc_ptn[VCSIZE];
    	VectorClass lh_final(0.0), vc_freq;
		VectorClass lh_ptn; // store likelihoods of VCSIZE consecutive patterns

		int *ptn_states_dad = aligned_alloc<int>(maxptn);
		for (ptn = 0; ptn < orig_nptn; ptn++)
			ptn_states_dad[ptn] = (aln->at(ptn))[dad->id];
		for (ptn = orig_nptn; ptn < nptn; ptn++)
			ptn_states_dad[ptn] = model_factory->unobserved_ptns[ptn-orig_nptn];
		// initialize beyond #patterns for efficiency
		for (ptn = nptn; ptn < maxptn; ptn++)
			ptn_states_dad[ptn] = aln->STATE_UNKNOWN;

		// copy dummy values because VectorClass will access beyond nptn
		for (ptn = nptn; ptn < maxptn; ptn++)
			memcpy(&dad_branch->partial_lh[ptn*block], dad_branch->partial_lh, block*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel private(ptn, i, j, vc_ptn, vc_freq, lh_ptn)
    {
    	VectorClass lh_final_th = 0.0;
#pragma omp for nowait
#endif
   		// main loop over all patterns with a step size of VCSIZE
		for (ptn = 0; ptn < orig_nptn; ptn+=VCSIZE) {
			//double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
            VectorClass vc_scale;
			for (j = 0; j < VCSIZE; j++) {
                vc_ptn[j] = 0.0;
				double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                UBYTE *scale_dad = dad_branch->scale_num + (ptn+j)*ncat_mix;
                // determine the min scaling
                UBYTE min_scale = *min_element(scale_dad, scale_dad+ncat_mix);
                vc_scale.insert(j, (double)min_scale);

				double *lh_node = &partial_lh_node[ptn_states_dad[ptn+j]*block];

                for (c = 0; c < ncat_mix; c++) {
                    VectorClass this_vc_ptn = 0.0;
                    if (scale_dad[c] <= min_scale+1) {
                        for (i = 0; i < nstates; i+=VCSIZE) {
                            this_vc_ptn = mul_add(VectorClass().load_a(&lh_node[i]), VectorClass().load_a(&partial_lh_dad[i]), this_vc_ptn);
                        }
                        if (scale_dad[c] == min_scale)
                            vc_ptn[j] += this_vc_ptn;
                        else
                            vc_ptn[j] += this_vc_ptn * VectorClass(SCALING_THRESHOLD);
                    }
                    lh_node += nstates;
                    partial_lh_dad += nstates;
                }
			}

			vc_freq.load_a(&ptn_freq[ptn]);
			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
			lh_ptn = log(abs(lh_ptn)) + vc_scale*LOG_SCALING_THRESHOLD;
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

        ASSERT(!isnan(tree_lh) & !isinf(tree_lh) && "Numerical underflow for SIMD lh-branch");

        // ascertainment bias correction
		if (orig_nptn < nptn) {
			lh_final = 0.0;
			lh_ptn = 0.0;
			for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
//				double *partial_lh_dad = &dad_branch->partial_lh[ptn*block];
                VectorClass vc_scale;
				lh_final += lh_ptn;
				for (j = 0; j < VCSIZE; j++) {
					vc_ptn[j] = 0.0;
					double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                    UBYTE *scale_dad = dad_branch->scale_num + (ptn+j)*ncat_mix;
                    // determine the min scaling
                    UBYTE min_scale = *min_element(scale_dad, scale_dad+ncat_mix);
                    vc_scale.insert(j, min_scale);

					int state_dad = ptn_states_dad[ptn+j];
					double *lh_node = &partial_lh_node[state_dad*block];

                    for (c = 0; c < ncat_mix; c++) {
                        VectorClass this_vc_ptn = 0.0;
                        if (scale_dad[c] <= min_scale+1) {
                            for (i = 0; i < nstates; i+=VCSIZE) {
                                this_vc_ptn = mul_add(VectorClass().load_a(&lh_node[i]), VectorClass().load_a(&partial_lh_dad[i]), this_vc_ptn);
                            }
                            if (scale_dad[c] == min_scale)
                                vc_ptn[j] += this_vc_ptn;
                            else
                                vc_ptn[j] += this_vc_ptn * VectorClass(SCALING_THRESHOLD);
                        }
                        lh_node += nstates;
                        partial_lh_dad += nstates;
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    if (min_scale >= 1)
                        vc_ptn[j] = vc_ptn[j] * VectorClass(SCALING_THRESHOLD);
				}
				// ptn_invar[ptn] is not aligned
				lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);
			}
			switch ((nptn-orig_nptn)%VCSIZE) {
			case 0: prob_const = horizontal_add(lh_final+lh_ptn); break;
			case 1: prob_const = horizontal_add(lh_final)+lh_ptn[0]; break;
			case 2: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]; break;
			case 3: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2]; break;
			default: ASSERT(0); break;
			}
		}
		aligned_free(ptn_states_dad);
		aligned_free(partial_lh_node);

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

            VectorClass vc_scale;
			for (j = 0; j < VCSIZE; j++) {
                vc_ptn[j] = 0.0;
				double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                double *partial_lh_node = node_branch->partial_lh + (ptn+j)*block;
                VectorClass *val_tmp = vc_val;
                UBYTE *scale_dad = dad_branch->scale_num + (ptn+j)*ncat_mix;
                UBYTE *scale_node = node_branch->scale_num + (ptn+j)*ncat_mix;
                // determine the min scaling
                UBYTE sum_scale[ncat_mix];
                UBYTE min_scale = sum_scale[0] = scale_dad[0]+scale_node[0];
                for (c = 1; c < ncat_mix; c++) {
                    sum_scale[c] = scale_dad[c] + scale_node[c];
                    min_scale = min(min_scale, sum_scale[c]);
                }
                vc_scale.insert(j, min_scale);

                for (c = 0; c < ncat_mix; c++) {
                    if (sum_scale[c] <= min_scale+1) {
                        VectorClass this_vc_ptn = 0.0;
                        for (i = 0; i < nstates; i+=VCSIZE) {
                            this_vc_ptn = mul_add(VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i]), val_tmp[i/VCSIZE], this_vc_ptn);
                        }
                        if (sum_scale[c] == min_scale)
                            vc_ptn[j] += this_vc_ptn;
                        else
                            vc_ptn[j] += this_vc_ptn * VectorClass(SCALING_THRESHOLD);
                    }
                    partial_lh_node += nstates;
                    partial_lh_dad += nstates;
                    val_tmp += nstates/VCSIZE;
                }
			}

			vc_freq.load_a(&ptn_freq[ptn]);

			lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);

			lh_ptn = log(abs(lh_ptn)) + vc_scale*LOG_SCALING_THRESHOLD;
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
		ASSERT(!isnan(tree_lh) && !isinf(tree_lh));

		if (orig_nptn < nptn) {
			// ascertainment bias correction
			lh_final = 0.0;
			lh_ptn = 0.0;

			for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
				lh_final += lh_ptn;

                VectorClass vc_scale;
                for (j = 0; j < VCSIZE; j++) {
                    vc_ptn[j] = 0.0;
                    double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                    double *partial_lh_node = node_branch->partial_lh + (ptn+j)*block;
                    VectorClass *val_tmp = vc_val;
                    UBYTE *scale_dad = dad_branch->scale_num + (ptn+j)*ncat_mix;
                    UBYTE *scale_node = node_branch->scale_num + (ptn+j)*ncat_mix;
                    // determine the min scaling
                    UBYTE sum_scale[ncat_mix];
                    UBYTE min_scale = sum_scale[0] = scale_dad[0]+scale_node[0];
                    for (c = 1; c < ncat_mix; c++) {
                        sum_scale[c] = scale_dad[c] + scale_node[c];
                        min_scale = min(min_scale, sum_scale[c]);
                    }
                    vc_scale.insert(j, min_scale);

                    for (c = 0; c < ncat_mix; c++) {
                        if (sum_scale[c] <= min_scale+1) {
                            VectorClass this_vc_ptn = 0.0;
                            for (i = 0; i < nstates; i+=VCSIZE) {
                                this_vc_ptn = mul_add(VectorClass().load_a(&partial_lh_node[i]) * VectorClass().load_a(&partial_lh_dad[i]), val_tmp[i/VCSIZE], this_vc_ptn);
                            }
                            if (sum_scale[c] == min_scale)
                                vc_ptn[j] += this_vc_ptn;
                            else
                                vc_ptn[j] += this_vc_ptn * VectorClass(SCALING_THRESHOLD);
                        }
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                        val_tmp += nstates/VCSIZE;
                    }
                    if (min_scale >= 1)
                        vc_ptn[j] *= VectorClass(SCALING_THRESHOLD);
                }

                /*
				for (j = 0; j < VCSIZE; j++)
					vc_ptn[j] = 0.0;

				for (i = 0; i < block; i+=VCSIZE) {
					for (j = 0; j < VCSIZE; j++) {
						vc_partial_lh_node[j].load_a(&partial_lh_node[i+j*block]);
						vc_partial_lh_dad[j].load_a(&partial_lh_dad[i+j*block]);
						vc_ptn[j] = mul_add(vc_val[i/VCSIZE] * vc_partial_lh_node[j], vc_partial_lh_dad[j], vc_ptn[j]);
					}
				}


                // bugfix 2016-01-21, prob_const can be rescaled
                for (j = 0; j < VCSIZE; j++)
                    if (dad_branch->scale_num[ptn+j] + node_branch->scale_num[ptn+j] >= 1)
                        vc_ptn[j] = vc_ptn[j] * SCALING_THRESHOLD;
                */
				// ptn_invar[ptn] is not aligned
				lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);
			}
			switch ((nptn-orig_nptn)%VCSIZE) {
			case 0: prob_const = horizontal_add(lh_final+lh_ptn); break;
			case 1: prob_const = horizontal_add(lh_final)+lh_ptn[0]; break;
			case 2: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]; break;
			case 3: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2]; break;
			default: ASSERT(0); break;
			}
		}
    }

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
        ASSERT(prob_const < 1.0 && prob_const >= 0.0);
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
    }

    aligned_free(vc_val);
    return tree_lh;
}

template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeLikelihoodFromBufferEigenSIMD() {


	ASSERT(theta_all && theta_computed);

	double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;

    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    size_t block = ncat_mix * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, j;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
//    size_t maxptn = ((nptn+VCSIZE-1)/VCSIZE)*VCSIZE;
    double *eval = model->getEigenvalues();
    ASSERT(eval);

	VectorClass *vc_val0 = (VectorClass*)aligned_alloc<double>(block);

	VectorClass vc_len = current_it->length;
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        double *eval_ptr = eval + (m)*nstates;
        size_t mycat = c%ncat;
		VectorClass vc_rate = site_rate->getRate(mycat);
		VectorClass vc_prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
		for (i = 0; i < nstates/VCSIZE; i++) {
			VectorClass cof = VectorClass().load_a(&eval_ptr[i*VCSIZE]) * vc_rate;
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
		lh_ptn = log(abs(lh_ptn));
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
	tree_lh += horizontal_add(lh_final) + buffer_scale_all;

    ASSERT(!isnan(tree_lh) && !isinf(tree_lh) && "Numerical underflow for SIMD lh-FromBuffer");

	if (orig_nptn < nptn) {
		// ascertaiment bias correction
		lh_final = 0.0;
		lh_ptn = 0.0;
		double prob_const;// df_const, ddf_const;
		double *theta = &theta_all[orig_nptn*block];

        UBYTE sum_scale_num[(nstates+VCSIZE)*ncat_mix];
        memset(sum_scale_num, 0, sizeof(UBYTE)*(nstates+VCSIZE));
        if (current_it->node->isLeaf())
            memcpy(sum_scale_num, current_it_back->scale_num+orig_nptn*ncat_mix, sizeof(UBYTE)*(nptn-orig_nptn)*ncat_mix);
        else if (current_it_back->node->isLeaf())
            memcpy(sum_scale_num, current_it->scale_num+orig_nptn*ncat_mix, sizeof(UBYTE)*(nptn-orig_nptn)*ncat_mix);
        else {
            UBYTE *cur_scale_num = current_it->scale_num + orig_nptn*ncat_mix;
            UBYTE *back_scale_num = current_it_back->scale_num + orig_nptn*ncat_mix;
            c = (nptn-orig_nptn)*ncat_mix;
            for (i = 0; i < c; i++)
                sum_scale_num[i] = cur_scale_num[i] + back_scale_num[i];
        }

        for (ptn = orig_nptn; ptn < nptn; ptn++) {
			//lh_final += lh_ptn;
			// initialization
            VectorClass this_vc_ptn = vc_val0[0] * VectorClass().load_a(theta);

            UBYTE *this_sum_scale = sum_scale_num + (ptn-orig_nptn)*ncat_mix;
            UBYTE min_scale = *min_element(this_sum_scale, this_sum_scale + ncat_mix);

			for (i = 1; i < block/VCSIZE; i++) {
                this_vc_ptn = mul_add(VectorClass().load_a(&theta[i*VCSIZE]), vc_val0[i], this_vc_ptn);
			}
			theta += block;

            // bugfix 2016-01-21, prob_const can be rescaled
            if (min_scale >= 1)
                this_vc_ptn *= VectorClass(SCALING_THRESHOLD);

			// no +I for +ASC!
			prob_const = horizontal_add(this_vc_ptn);

		}

        /*
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
        */
    	prob_const = log(1.0 - prob_const);
    	tree_lh -= aln->getNSite() * prob_const;
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
	}

    aligned_free(vc_val0);

    return tree_lh;
}

/****************************************************************************
        Highly optimized Parsimony function
 ****************************************************************************/

#ifdef _MSC_VER
	#define MEM_ALIGN_BEGIN __declspec(align(32))
	#define MEM_ALIGN_END
#else
	#define MEM_ALIGN_BEGIN
	#define MEM_ALIGN_END __attribute__((aligned(32)))
#endif

inline UINT fast_popcount(Vec4ui &x) {
    MEM_ALIGN_BEGIN UINT vec[4] MEM_ALIGN_END;
    x.store_a(vec);
    return popcount_lauradoux(vec, 4);
}

inline UINT fast_popcount(Vec8ui &x) {
#if defined (__GNUC__) || defined(__clang__)
    MEM_ALIGN_BEGIN uint64_t vec[4] MEM_ALIGN_END;
    MEM_ALIGN_BEGIN uint64_t res[4] MEM_ALIGN_END;
    Vec8ui y;
    x.store_a(vec);
    __asm("popcntq %1, %0" : "=r"(res[0]) : "r"(vec[0]) : );
    __asm("popcntq %1, %0" : "=r"(res[1]) : "r"(vec[1]) : );
    __asm("popcntq %1, %0" : "=r"(res[2]) : "r"(vec[2]) : );
    __asm("popcntq %1, %0" : "=r"(res[3]) : "r"(vec[3]) : );
    y.load_a(res);
    return horizontal_add(y);
#else
    MEM_ALIGN_BEGIN uint64_t vec[4] MEM_ALIGN_END;
    MEM_ALIGN_BEGIN int res[4] MEM_ALIGN_END;
    Vec4ui y;
    x.store_a(vec);
    res[0] = _mm_popcnt_u64(vec[0]);
    res[1] = _mm_popcnt_u64(vec[1]);
    res[2] = _mm_popcnt_u64(vec[2]);
    res[3] = _mm_popcnt_u64(vec[3]);
    y.load_a(res);
    return horizontal_add(y);
#endif

}


inline void horizontal_popcount(Vec4ui &x) {
    MEM_ALIGN_BEGIN UINT vec[4] MEM_ALIGN_END;
    x.store_a(vec);
    vec[0] = vml_popcnt(vec[0]);
    vec[1] = vml_popcnt(vec[1]);
    vec[2] = vml_popcnt(vec[2]);
    vec[3] = vml_popcnt(vec[3]);
    x.load_a(vec);
}

inline void horizontal_popcount(Vec8ui &x) {
    MEM_ALIGN_BEGIN UINT vec[8] MEM_ALIGN_END;
    x.store_a(vec);
    vec[0] = vml_popcnt(vec[0]);
    vec[1] = vml_popcnt(vec[1]);
    vec[2] = vml_popcnt(vec[2]);
    vec[3] = vml_popcnt(vec[3]);
    vec[4] = vml_popcnt(vec[4]);
    vec[5] = vml_popcnt(vec[5]);
    vec[6] = vml_popcnt(vec[6]);
    vec[7] = vml_popcnt(vec[7]);
    x.load_a(vec);
}

template<class VectorClass>
void PhyloTree::computePartialParsimonyFastSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (dad_branch->partial_lh_computed & 2) {
        return;
    }
    PhyloNode *node    = dad_branch->getNode();
    int nstates        = aln->getMaxNumStates();
    int site           = 0;
    const int VCSIZE   = VectorClass::size();
    const int NUM_BITS = VectorClass::size() * UINT_BITS;

    dad_branch->partial_lh_computed |= 2;

    if (node->isLeaf() && dad) {
        // external node
        vector<Alignment*> *partitions = NULL;
        if (aln->isSuperAlignment())
            partitions = &((SuperAlignment*)aln)->partitions;
        else {
            partitions = new vector<Alignment*>;
            partitions->push_back(aln);
        }
        if (aln->ordered_pattern.empty())
            aln->orderPatternByNumChars();
        int leafid = node->id;
        int pars_size = getBitsBlockSize();
        memset(dad_branch->partial_pars, 0, pars_size*sizeof(UINT));
    	int ambi_aa[] = {2, 3, 5, 6, 9, 10}; // {4+8, 32+64, 512+1024};
        UINT *x = dad_branch->partial_pars;
        int start_pos = 0;

        for (vector<Alignment*>::iterator alnit = partitions->begin(); alnit != partitions->end(); alnit++) {
            int end_pos = start_pos + (*alnit)->ordered_pattern.size();
            switch ((*alnit)->seq_type) {
            case SEQ_DNA:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (state < 4) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            x[state*VCSIZE + site/UINT_BITS] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            UINT bit1 = (1 << (site%UINT_BITS));
                            UINT *p = x+(site/UINT_BITS);
                            p[0] |= bit1;
                            p[VCSIZE] |= bit1;
                            p[2*VCSIZE] |= bit1;
                            p[3*VCSIZE] |= bit1;
                        }
                    } else {
                        state -= 3;
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            UINT *p = x + ((site/UINT_BITS));
                            
                            UINT bit1 = (1 << (site%UINT_BITS));
                            for (int i = 0; i < 4; i++)
                                if (state & (1<<i))
                                    p[i*VCSIZE] |= bit1;
                        }
                    }
                }
                break;
            case SEQ_PROTEIN:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (state < 20) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            x[state*VCSIZE + site/UINT_BITS] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            UINT bit1 = (1 << (site%UINT_BITS));
                            UINT *p = x+(site/UINT_BITS);
                            for (int i = 0; i < 20; i++)
                                p[i*VCSIZE] |= bit1;
                        }
                    } else {
                        ASSERT(state < 23);
                        state = (state-20)*2;
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            UINT *p = x + ((site/UINT_BITS));
                            UINT bit1 = (1 << (site%UINT_BITS));

                            p[ambi_aa[state]*VCSIZE] |= bit1;
                            p[ambi_aa[state+1]*VCSIZE] |= bit1;
                        }
                    }
                }
                break;
            default:
                for (int patid = start_pos; patid != end_pos; patid++) {
                    Alignment::iterator pat = aln->ordered_pattern.begin()+ patid;
                    int state = pat->at(leafid);
                    int freq = pat->frequency;
                    if (state < (*alnit)->num_states) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            x[state*VCSIZE + site/UINT_BITS] |= (1 << (site % UINT_BITS));
                        }
                    } else if (state == (*alnit)->STATE_UNKNOWN) {
                        for (int j = 0; j < freq; j++, site++) {
                            if (site == NUM_BITS) {
                                x += nstates*VCSIZE;
                                site = 0;
                            }
                            UINT bit1 = (1 << (site%UINT_BITS));
                            UINT *p = x+(site/UINT_BITS);
                            for (int i = 0; i < (*alnit)->num_states; i++)
                                p[i*VCSIZE] |= bit1;
                        }
                    } else {
                        ASSERT(0);
                    }
                }
                break;
            } // end of switch
            start_pos = end_pos;
        } // of end FOR LOOP

        ASSERT(start_pos == aln->ordered_pattern.size());
//        assert(site == aln->num_informative_sites % NUM_BITS);
        // add dummy states
        if (site > 0 && site < NUM_BITS) {
            x += site/UINT_BITS;
        	*x |= ~((1<<(site%UINT_BITS)) - 1);
            x++;
            int max_sites = ((site+UINT_BITS-1)/UINT_BITS);
            memset(x, 255, (VCSIZE - max_sites)*sizeof(UINT));
        }
        if (!aln->isSuperAlignment())
            delete partitions;
    } else {
        // internal node
        ASSERT(node->degree() == 3); // it works only for strictly bifurcating tree
        PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
            if (nei->node->name != ROOT_NAME && (nei->partial_lh_computed & 2) == 0) {
                computePartialParsimonyFastSIMD<VectorClass>(nei, node);
            }
            if (!left) left = pit; else right = pit;
        }
//        VectorClass score = 0;
        UINT score = 0;
        int nsites = (aln->num_informative_sites+NUM_BITS-1)/NUM_BITS;
        int entry_size = nstates * VCSIZE;
        
        switch (nstates) {
        case 4:
            #ifdef _OPENMP
            #pragma omp parallel for private (site) reduction(+: score) if(nsites>200)
            #endif
			for (site = 0; site<nsites; site++) {
                size_t offset = entry_size*site;
                VectorClass *x = (VectorClass*)(left->partial_pars + offset);
                VectorClass *y = (VectorClass*)(right->partial_pars + offset);
                VectorClass *z = (VectorClass*)(dad_branch->partial_pars + offset);
                z[0] = x[0] & y[0];
                z[1] = x[1] & y[1];
                z[2] = x[2] & y[2];
                z[3] = x[3] & y[3];
                VectorClass w = z[0] | z[1] | z[2] | z[3];
				w = ~w;
                z[0] |= w & (x[0] | y[0]);
                z[1] |= w & (x[1] | y[1]);
                z[2] |= w & (x[2] | y[2]);
                z[3] |= w & (x[3] | y[3]);
//				horizontal_popcount(w);
//                score += w;
                score += fast_popcount(w);
//                x += 4;
//                y += 4;
//                z += 4;
			}

			break;
        default:
            #ifdef _OPENMP
            #pragma omp parallel for private (site) reduction(+: score) if(nsites > 800/nstates)
            #endif
			for (site = 0; site<nsites; site++) {
                size_t offset = entry_size*site;
                VectorClass *x = (VectorClass*)(left->partial_pars + offset);
                VectorClass *y = (VectorClass*)(right->partial_pars + offset);
                VectorClass *z = (VectorClass*)(dad_branch->partial_pars + offset);
				int i;
				VectorClass w = 0;
				for (i = 0; i < nstates; i++) {
                    z[i] = x[i] & y[i];
                    w |= z[i];
				}
				w = ~w;
				for (i = 0; i < nstates; i++) {
                    z[i] |= w & (x[i] | y[i]);
				}
//				horizontal_popcount(w);
//                score += w;
                score += fast_popcount(w);
                x += nstates;
                y += nstates;
                z += nstates;
			}
			break;
        }
//        UINT sum_score = horizontal_add(score); 
//        UINT *zscore = (UINT*)z;
//        UINT *xscore = (UINT*)x;
//        UINT *yscore = (UINT*)y;
        dad_branch->partial_pars[nstates*VCSIZE*nsites] = score + left->partial_pars[nstates*VCSIZE*nsites] + right->partial_pars[nstates*VCSIZE*nsites];
    }
}

template<class VectorClass>
int PhyloTree::computeParsimonyBranchFastSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    ASSERT(node_branch);
    if (!central_partial_pars)
        initializeAllPartialPars();
    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFastSIMD<VectorClass>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimonyFastSIMD<VectorClass>(node_branch, node);
    int site;
    int nstates = aln->getMaxNumStates();

//    VectorClass score = 0;
//    VectorClass w;

    const int NUM_BITS = VectorClass::size() * UINT_BITS;
    int nsites = (aln->num_informative_sites + NUM_BITS - 1)/NUM_BITS;
    int entry_size = nstates * VectorClass::size();
    
    int scoreid = nsites*entry_size;
    UINT sum_end_node = (dad_branch->partial_pars[scoreid] + node_branch->partial_pars[scoreid]);
    UINT score = sum_end_node;
    UINT lower_bound = best_pars_score;
    if (branch_subst) lower_bound = INT_MAX;
    
    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for private (site) reduction(+: score) if(nsites>200)
        #endif
		for (site = 0; site < nsites; site++) {
            size_t offset = entry_size*site;
            VectorClass *x = (VectorClass*)(dad_branch->partial_pars + offset);
            VectorClass *y = (VectorClass*)(node_branch->partial_pars + offset);
            VectorClass w = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
			w = ~w;
//			horizontal_popcount(w);
//            score += w;
            score += fast_popcount(w);
            #ifndef _OPENMP
            if (score >= lower_bound) 
                break;
            #endif
		}
		break;
    default:
        #ifdef _OPENMP
        #pragma omp parallel for private (site) reduction(+: score) if(nsites > 800/nstates)
        #endif
		for (site = 0; site < nsites; site++) {
            size_t offset = entry_size*site;
            VectorClass *x = (VectorClass*)(dad_branch->partial_pars + offset);
            VectorClass *y = (VectorClass*)(node_branch->partial_pars + offset);
            VectorClass w = x[0] & y[0];
			for (int i = 1; i < nstates; i++) {
                w |= x[i] & y[i];
			}
			w = ~w;
//			horizontal_popcount(w);
//            score += w;
            score += fast_popcount(w);
            #ifndef _OPENMP
            if (score >= lower_bound) 
                break;
            #endif
		}
		break;
    }
//    UINT sum_score = horizontal_add(score);
//    if (branch_subst)
//        *branch_subst = sum_score;
    if (branch_subst)
        *branch_subst = score - sum_end_node;
//    UINT *xscore = (UINT*)x;
//    UINT *yscore = (UINT*)y;
//    sum_score += *xscore + *yscore;
//    score += *xscore + *yscore;
//    return sum_score;
    return score;
}


#endif /* PHYLOKERNELSAFE_H_ */
