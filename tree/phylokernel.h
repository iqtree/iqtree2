/*
 * phylokernel.h
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */

#ifndef PHYLOKERNEL_H_
#define PHYLOKERNEL_H_

#include "phylotree.h"
#include "alignment/superalignment.h"
#if defined(CLANG_UNDER_VS)
#   define _mm_popcnt_u64 __popcnt64
#endif
#include <vectorclass/vectorclass.h>

#ifdef __SSE2__
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

template <class Numeric, class VectorClass>
Numeric PhyloTree::dotProductSIMD(Numeric *x, Numeric *y, int size) {
    VectorClass res = VectorClass().load_a(x) * VectorClass().load_a(y);
    for (int i = VectorClass::size(); i < size; i += VectorClass::size())
        res = mul_add(VectorClass().load_a(&x[i]), VectorClass().load_a(&y[i]), res);
    return horizontal_add(res);
}

/************************************************************************************************
 *
 *   Highly optimized vectorized versions of likelihood functions
 *
 *************************************************************************************************/

/*
template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                  LikelihoodBufferSet& buffers) {

    // don't recompute the likelihood
    assert(dad);
    if (dad_branch->isLikelihoodComputed()) {
        return;
    }
    dad_branch->setLikelihoodComputed(true);

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
    assert(nstates == aln->num_states && nstates >= VCSIZE && VCSIZE == VectorClass().size());
    assert(model->isReversible()); // only works with reversible model!
    const size_t nstatesqr=nstates*nstates;
    size_t i, x, j;
    size_t block = nstates * ncat_mix;
    size_t tip_block = nstates * model->getNMixtures();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = m*nstatesqr;
    }

    // internal node
    dad_branch->lh_scale_factor = 0.0;
    PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
    int num_leaves = 0;
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        if (!left) left = nei; else right = nei;
        if (!nei->isLikelihoodComputed()) {
            computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(nei, node, buffers);
        }
        dad_branch->lh_scale_factor += nei->lh_scale_factor;
        if ((*it)->node->isLeaf()) num_leaves++;
    }

    if (params->lh_mem_save == LM_PER_NODE && !dad_branch->partial_lh) {
        // re-orient partial_lh
        bool done = false;
        FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it2, child) {
            PhyloNeighbor *backnei = child->findNeighbor(node);
            if (backnei->partial_lh) {
                std::swap(dad_branch->partial_lh, backnei->partial_lh);
                std::swap(dad_branch->scale_num,  backnei->scale_num);
                backnei->clearComputedFlags();
                done = true;
                break;
            }
        }
        assert(done && "partial_lh is not re-oriented");
    }

    double *evec = model->getEigenvectors();
    double *inv_evec = model->getInverseEigenvectors();

    assert(inv_evec && evec);
//    for (i = 0; i < tip_block; i++) {
//        for (x = 0; x < nstates/VCSIZE; x++)
//            // inv_evec is not aligned!
//            vc_inv_evec[i*nstates/VCSIZE+x].load_a(&inv_evec[i*nstates+x*VCSIZE]);
//    }
    double *eval = model->getEigenvalues();


    VectorClass *echildren = aligned_alloc<VectorClass>(block*nstates/VCSIZE*(node->degree()-1));
    double *partial_lh_leaves = NULL;
    if (num_leaves > 0)
        partial_lh_leaves = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block*num_leaves);
    VectorClass *echild = echildren;
    double *partial_lh_leaf = partial_lh_leaves;
    
    
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
        VectorClass expchild[nstates/VCSIZE];
        VectorClass *echild_ptr = echild;
        // precompute information buffer
        for (c = 0; c < ncat_mix; c++) {
            VectorClass len_child = site_rate->getRate(c%ncat) * child->length;
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
        if (child->node->isLeaf()) {
            vector<int>::iterator it;
            for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
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

        //--------------------- multifurcating node ------------------//
        double sum_scale = 0.0;
        // now for-loop computing partial_lh over all site-patterns
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, c, x, i) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn++) {
            double partial_lh_all[block];
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            dad_branch->scale_num[ptn] = 0;
                
            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = (double*)echildren;

            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, child) {
                if (child->node->isLeaf()) {
                    // external node
                    int state_child = (ptn < orig_nptn) ? (aln->at(ptn))[child->node->id] : model_factory->unobserved_ptns[ptn-orig_nptn];
                    double *child_lh = partial_lh_leaf + state_child*block;
                    for (c = 0; c < block; c++) {
                        // compute real partial likelihood vector
                        partial_lh_all[c] *= child_lh[c];
                    }
                    partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                } else {
                    // internal node
                    double *partial_lh = partial_lh_all;
                    double *partial_lh_child = child->partial_lh + ptn*block;
                    dad_branch->scale_num[ptn] += child->scale_num[ptn];

                    double *echild_ptr = echild;
                    for (c = 0; c < ncat_mix; c++) {
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
            } // FOR_NEIGHBOR
            
        
            // compute dot-product with inv_eigenvector
            double lh_max = 0.0;
            double *partial_lh_tmp = partial_lh_all;
            double *partial_lh = dad_branch->partial_lh + ptn*block;
            for (c = 0; c < ncat_mix; c++) {
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
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }
            // check if one should scale partial likelihoods
            if (lh_max < SCALING_THRESHOLD) {
                partial_lh = dad_branch->partial_lh + ptn*block;
                if (lh_max == 0.0) {
                    // for very shitty data
                    for (c = 0; c < ncat_mix; c++)
                        memcpy(&partial_lh[c*nstates], &tip_partial_lh[aln->STATE_UNKNOWN*nstates], nstates*sizeof(double));
                    sum_scale += LOG_SCALING_THRESHOLD* 4 * ptn_freq[ptn];
                    //sum_scale += log(lh_max) * ptn_freq[ptn];
                    dad_branch->scale_num[ptn] += 4;
                    size_t nsite = aln->getNSite();
                    for (i = 0, x = 0; i < nsite && x < ptn_freq[ptn]; i++)
                        if (aln->getPatternID(i) == ptn) {
                            outWarning((string)"Numerical underflow for site " + convertIntToString(i+1));
                            x++;
                        }
                } else if (ptn_invar[ptn] == 0.0) {
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
    } else if (left->node->isLeaf() && right->node->isLeaf()) {
        // special treatment for TIP-TIP (cherry) case

        // pre compute information for both tips
        double *partial_lh_left = partial_lh_leaves;
        double *partial_lh_right = partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;

        // assign pointers for left and right partial_lh
//        double **lh_left_ptr = aligned_alloc<double*>(nptn);
//        double **lh_right_ptr = aligned_alloc<double*>(nptn);
//        for (ptn = 0; ptn < orig_ntn; ptn++) {
//            lh_left_ptr[ptn] = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
//            lh_right_ptr[ptn] = &partial_lh_right[block * (aln->at(ptn))[right->node->id]];
//        }
//        for (ptn = orig_ntn; ptn < nptn; ptn++) {
//            lh_left_ptr[ptn] = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
//            lh_right_ptr[ptn] = &partial_lh_right[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
//        }

        // scale number must be ZERO
        memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
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
        memcpy(dad_branch->scale_num, right->scale_num, nptn * sizeof(UBYTE));

        // pre compute information for left tip
        double *partial_lh_left = partial_lh_leaves;


        // assign pointers for partial_lh_left
//        double **lh_left_ptr = aligned_alloc<double*>(nptn);
//        for (ptn = 0; ptn < orig_ntn; ptn++) {
//            lh_left_ptr[ptn] = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
//        }
//        for (ptn = orig_ntn; ptn < nptn; ptn++) {
//            lh_left_ptr[ptn] = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_ntn]];
//        }
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

            double *lh_left;
            if (ptn < orig_nptn) {
                lh_left = &partial_lh_left[block *  (aln->at(ptn))[left->node->id]];
            } else {
                lh_left = &partial_lh_left[block * model_factory->unobserved_ptns[ptn-orig_nptn]];
            }
            vc_max = 0.0;
            for (c = 0; c < ncat_mix; c++) {
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
                lh_left += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
            }
            // check if one should scale partial likelihoods
            double lh_max = horizontal_max(vc_max);
            if (lh_max < SCALING_THRESHOLD && ptn_invar[ptn] == 0.0) {
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

        //aligned_free(lh_left_ptr);

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
            for (c = 0; c < ncat_mix; c++) {
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
                partial_lh += nstates;
                partial_lh_left += nstates;
                partial_lh_right += nstates;
            }

            // check if one should scale partial likelihoods
            double lh_max = horizontal_max(vc_max);
            if (lh_max < SCALING_THRESHOLD && ptn_invar[ptn] == 0.0) {
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

    if (partial_lh_leaves)
        aligned_free(partial_lh_leaves);
    aligned_free(echildren);
}

template <class VectorClass, const int VCSIZE, const int nstates>
void PhyloTree::computeLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                               double &df, double &ddf,
                                               LikelihoodBufferSet& buffers) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    if (!dad_branch->isLikelihoodComputed())
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad), buffers;
    if (!node_branch->isLikelihoodComputed())
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node, buffers);
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
    assert(eval);

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

    assert(theta_all);
    if (!theta_computed) {
        theta_computed = true;
        // precompute theta for fast branch length optimization

        if (dad->isLeaf()) {
            // special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c)
#endif
            for (ptn = 0; ptn < nptn; ptn++) {
                double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
                double *theta = theta_all + ptn*block;
                double *this_tip_partial_lh = tip_partial_lh + tip_block*((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]);
                for (c = 0; c < ncat_mix; c++) {
                    double *lh_dad = this_tip_partial_lh + mix_addr_nstates[c];
                    for (i = 0; i < nstates; i+=VCSIZE) {
                        (VectorClass().load_a(&lh_dad[i]) * VectorClass().load_a(&partial_lh_dad[i])).store_a(&theta[i]);
                    }
                    partial_lh_dad += nstates;
                    theta += nstates;
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
    if (isnan(df) || isinf(df)) {
        df = 0.0;
        ddf = 0.0;
//        outWarning("Numerical instability (some site-likelihood = 0)");
    }


//    assert(isnormal(tree_lh));
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
            assert(0);
            break;
        }
        prob_const = 1.0 - prob_const;
        double df_frac = df_const / prob_const;
        double ddf_frac = ddf_const / prob_const;
        size_t nsites = aln->getNSite();
        df += nsites * df_frac;
        ddf += nsites *(ddf_frac + df_frac*df_frac);
    }
    assert(!isnan(df));
    aligned_free(vc_val2);
    aligned_free(vc_val1);
    aligned_free(vc_val0);
}


template <class VectorClass, const int VCSIZE, const int nstates>
double PhyloTree::computeLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                   LikelihoodBufferSet& buffers) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    if (!dad_branch->isLikelihoodComputed()) {
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(dad_branch, dad, buffers);
    }
    if (!node_branch->isLikelihoodComputed()) {
        computePartialLikelihoodEigenSIMD<VectorClass, VCSIZE, nstates>(node_branch, node, buffers);
    }
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
    assert(eval);

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
//        double **lh_states_dad = aligned_alloc<double*>(maxptn);
//        for (ptn = 0; ptn < orig_nptn; ptn++)
//            lh_states_dad[ptn] = &tip_partial_lh[(aln->at(ptn))[dad->id] * tip_block];
//        for (ptn = orig_nptn; ptn < nptn; ptn++)
//            lh_states_dad[ptn] = &tip_partial_lh[model_factory->unobserved_ptns[ptn-orig_nptn] * tip_block];
//        // initialize beyond #patterns for efficiency
//        for (ptn = nptn; ptn < maxptn; ptn++)
//            lh_states_dad[ptn] = &tip_partial_lh[aln->STATE_UNKNOWN * tip_block];
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

            for (j = 0; j < VCSIZE; j++) {
                vc_ptn[j] = 0.0;
                double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                int state_dad = ptn_states_dad[ptn+j];
                double *lh_node = &partial_lh_node[state_dad*block];
                for (i = 0; i < block; i+=VCSIZE) {
                    vc_ptn[j] = mul_add(VectorClass().load_a(&lh_node[i]),
                            VectorClass().load_a(&partial_lh_dad[i]), vc_ptn[j]);
                }
            }

            // initialize vc_tip_partial_lh
//            for (j = 0; j < VCSIZE; j++) {
//                double *lh_dad = lh_states_dad[ptn+j];
//                for (i = 0; i < nstates/VCSIZE; i++) {
//                    vc_tip_partial_lh[j*(nstates/VCSIZE)+i].load_a(&lh_dad[i*VCSIZE]);
//                }
//                vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block]);
//                vc_ptn[j] = vc_val[0] * vc_tip_partial_lh[j*(nstates/VCSIZE)] * vc_partial_lh_dad[j];
//            }
//
//            // compute vc_ptn
//            for (i = 1; i < block/VCSIZE; i++)
//                for (j = 0; j < VCSIZE; j++) {
//                    vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block+i*VCSIZE]);
//                    vc_ptn[j] = mul_add(vc_val[i] * vc_tip_partial_lh[j*(nstates/VCSIZE)+i%(nstates/VCSIZE)],
//                            vc_partial_lh_dad[j], vc_ptn[j]);
//                }

            vc_freq.load_a(&ptn_freq[ptn]);
            lh_ptn = horizontal_add(vc_ptn) + VectorClass().load_a(&ptn_invar[ptn]);
            lh_ptn = log(abs(lh_ptn));
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
        if (isnan(tree_lh) || isinf(tree_lh)) {
            cout << "WARNING: Numerical underflow caused by alignment sites there";
            i = aln->getNSite();
            for (j = 0; j < i; j++) {
                ptn = aln->getPatternID(j);
                if (isnan(_pattern_lh[ptn]) || isinf(_pattern_lh[ptn])) {
                    cout << " " << j+1;
                }
            }
            tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
            for (ptn = 0; ptn < orig_nptn; ptn++) {
                if (isnan(_pattern_lh[ptn]) || isinf(_pattern_lh[ptn])) {
                    _pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
                }
                tree_lh += _pattern_lh[ptn] * ptn_freq[ptn];
            }
            cout << endl;
            if (verbose_mode >= VB_MED) {
                printTree(cout);
                cout << endl;
            }
//            cout << "WARNING: Tree log-likelihood is set to " << tree_lh << endl;
        }

        if (orig_nptn < nptn) {
            lh_final = 0.0;
            lh_ptn = 0.0;
            for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
//                double *partial_lh_dad = &dad_branch->partial_lh[ptn*block];
                lh_final += lh_ptn;
                for (j = 0; j < VCSIZE; j++) {
                    vc_ptn[j] = 0.0;
                    double *partial_lh_dad = dad_branch->partial_lh + (ptn+j)*block;
                    int state_dad = ptn_states_dad[ptn+j];
                    double *lh_node = &partial_lh_node[state_dad*block];
                    for (i = 0; i < block; i+=VCSIZE) {
                        vc_ptn[j] = mul_add(VectorClass().load_a(&lh_node[i]),
                                VectorClass().load_a(&partial_lh_dad[i]), vc_ptn[j]);
                    }
                }

                // bugfix 2016-01-21, prob_const can be rescaled
                for (j = 0; j < VCSIZE; j++)
                    if (dad_branch->scale_num[ptn+j] >= 1)
                        vc_ptn[j] = vc_ptn[j] * SCALING_THRESHOLD;

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
        aligned_free(ptn_states_dad);
        aligned_free(partial_lh_node);
        


        // ascertainment bias correction
//        if (orig_nptn < nptn) {
//            lh_final = 0.0;
//            lh_ptn = 0.0;
//            for (ptn = orig_nptn; ptn < nptn; ptn+=VCSIZE) {
//                double *partial_lh_dad = &dad_branch->partial_lh[ptn*block];
//                lh_final += lh_ptn;
//
//                // initialize vc_tip_partial_lh
//                for (j = 0; j < VCSIZE; j++) {
//                    double *lh_dad = lh_states_dad[ptn+j];
//                    for (i = 0; i < nstates/VCSIZE; i++) {
//                        vc_tip_partial_lh[j*(nstates/VCSIZE)+i].load(&lh_dad[i*VCSIZE]); // lh_dad is not aligned!
//                    }
//                    vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block]);
//                    vc_ptn[j] = vc_val[0] * vc_tip_partial_lh[j*(nstates/VCSIZE)] * vc_partial_lh_dad[j];
//                }
//
//                // compute vc_ptn
//                for (i = 1; i < block/VCSIZE; i++)
//                    for (j = 0; j < VCSIZE; j++) {
//                        vc_partial_lh_dad[j].load_a(&partial_lh_dad[j*block+i*VCSIZE]);
//                        vc_ptn[j] = mul_add(vc_val[i] * vc_tip_partial_lh[j*(nstates/VCSIZE)+i%(nstates/VCSIZE)],
//                                vc_partial_lh_dad[j], vc_ptn[j]);
//                    }
//
//                // bugfix 2016-01-21, prob_const can be rescaled
//                for (j = 0; j < VCSIZE; j++)
//                    if (dad_branch->scale_num[ptn+j] >= 1)
//                        vc_ptn[j] = vc_ptn[j] * SCALING_THRESHOLD;
//
//                // ptn_invar[ptn] is not aligned
//                lh_ptn = horizontal_add(vc_ptn) + VectorClass().load(&ptn_invar[ptn]);
//            }
//            switch ((nptn-orig_nptn)%VCSIZE) {
//            case 0: prob_const = horizontal_add(lh_final+lh_ptn); break;
//            case 1: prob_const = horizontal_add(lh_final)+lh_ptn[0]; break;
//            case 2: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]; break;
//            case 3: prob_const = horizontal_add(lh_final)+lh_ptn[0]+lh_ptn[1]+lh_ptn[2]; break;
//            default: assert(0); break;
//            }
//        }
//        aligned_free(lh_states_dad);
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

            lh_ptn = log(abs(lh_ptn));
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

                // bugfix 2016-01-21, prob_const can be rescaled
                for (j = 0; j < VCSIZE; j++)
                    if (dad_branch->scale_num[ptn+j] + node_branch->scale_num[ptn+j] >= 1)
                        vc_ptn[j] = vc_ptn[j] * SCALING_THRESHOLD;

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
        assert(prob_const < 1.0 && prob_const >= 0.0);
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


    assert(theta_all && theta_computed);

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
    assert(eval);

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
    tree_lh += horizontal_add(lh_final);
    if (isnan(tree_lh) || isinf(tree_lh)) {
        cout << "WARNING: Numerical underflow caused by alignment sites here";
        i = aln->getNSite();
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
        for (ptn = 0; ptn < orig_nptn; ptn++) {
            if (isnan(_pattern_lh[ptn]) || isinf(_pattern_lh[ptn])) {
                _pattern_lh[ptn] = LOG_SCALING_THRESHOLD*4; // log(2^(-1024))
            }
            tree_lh += _pattern_lh[ptn] * ptn_freq[ptn];
        }
    }

    if (orig_nptn < nptn) {
        // ascertaiment bias correction
        lh_final = 0.0;
        lh_ptn = 0.0;
        double prob_const;// df_const, ddf_const;
        double *theta = &theta_all[orig_nptn*block];

        UBYTE sum_scale_num[nstates+VCSIZE];
        memset(sum_scale_num, 0, sizeof(UBYTE)*(nstates+VCSIZE));
        if (current_it->node->isLeaf())
            memcpy(sum_scale_num, current_it_back->scale_num+orig_nptn, sizeof(UBYTE)*(nptn-orig_nptn));
        else if (current_it_back->node->isLeaf())
            memcpy(sum_scale_num, current_it->scale_num+orig_nptn, sizeof(UBYTE)*(nptn-orig_nptn));
        else {
            for (ptn = orig_nptn; ptn < nptn; ptn++)
                sum_scale_num[ptn-orig_nptn] = current_it->scale_num[ptn] + current_it_back->scale_num[ptn];
        }

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

            // bugfix 2016-01-21, prob_const can be rescaled
            for (j = 0; j < VCSIZE; j++)
                if (sum_scale_num[ptn+j-orig_nptn] >= 1)
                    vc_ptn[j] = vc_ptn[j] * SCALING_THRESHOLD;

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
*/
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
    MEM_ALIGN_BEGIN uint64_t vec[4] MEM_ALIGN_END;
    MEM_ALIGN_BEGIN uint64_t res[4] MEM_ALIGN_END;
    x.store_a(vec);
    #if (defined (__GNUC__) || defined(__clang__)) && !defined(CLANG_UNDER_VS)
        __asm("popcntq %1, %0" : "=r"(res[0]) : "r"(vec[0]) : );
        __asm("popcntq %1, %0" : "=r"(res[1]) : "r"(vec[1]) : );
        __asm("popcntq %1, %0" : "=r"(res[2]) : "r"(vec[2]) : );
        __asm("popcntq %1, %0" : "=r"(res[3]) : "r"(vec[3]) : );
    #else
        res[0] = _mm_popcnt_u64(vec[0]);
        res[1] = _mm_popcnt_u64(vec[1]);
        res[2] = _mm_popcnt_u64(vec[2]);
        res[3] = _mm_popcnt_u64(vec[3]);
    #endif
    Vec8ui y;
    y.load_a(res);
    return horizontal_add(y);
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
    if (dad_branch->isParsimonyComputed()) {
        return;
    }
    PhyloNode* node     = dad_branch->getNode();
    int        nstates  = aln->getMaxNumStates();
    int        site     = 0;
    const int  VCSIZE   = VectorClass::size();
    const int  NUM_BITS = VectorClass::size() * UINT_BITS;

    dad_branch->setParsimonyComputed(true);

    if (node->name == ROOT_NAME) {
        ASSERT(dad);
        // special treatment for root node
//        if (aln->ordered_pattern.empty())
//            aln->orderPatternByNumChars();
//        ASSERT(!aln->ordered_pattern.empty());
        size_t pars_size = getBitsBlockSize();
        memset(dad_branch->partial_pars, 255, pars_size*sizeof(UINT));
        size_t nsites = (aln->num_parsimony_sites+NUM_BITS-1)/NUM_BITS;
        dad_branch->partial_pars[nstates*VCSIZE*nsites] = 0;
    } else if (node->isLeaf() && dad) {
        // external node
        vector<Alignment*> *partitions = NULL;
        if (aln->isSuperAlignment())
            partitions = &((SuperAlignment*)aln)->partitions;
        else {
            partitions = new vector<Alignment*>;
            partitions->push_back(aln);
        }
//        if (aln->ordered_pattern.empty())
//            aln->orderPatternByNumChars();
//        ASSERT(!aln->ordered_pattern.empty());
        int leafid = node->id;
        size_t pars_size = getBitsBlockSize();
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
                if (aln->seq_type == SEQ_POMO && state >= nstates 
                    && state < aln->STATE_UNKNOWN) {
                    state -= nstates;
                    ASSERT(state < aln->pomo_sampled_states.size());
                    int id1 = aln->pomo_sampled_states[state] & 3;
                    int id2 = (aln->pomo_sampled_states[state] >> 16) & 3;
                    int value1 = (aln->pomo_sampled_states[state] >> 2) & 16383;
                    int value2 = aln->pomo_sampled_states[state] >> 18;
                    double weight1 = ((double)value1)/(value1+value2);
//                    int N = aln->virtual_pop_size;
//                    int M = value1 + value2;

                    // 2016-09-30: resolving polymorphic states to fixed states

                    // value1 = value1*N/(value1+value2);
                    int real_state;

                    if (weight1 < 1.0/4)
                        real_state = id2;
                    else if (weight1 > 3.0/4)
                        real_state = id1;
                    else
                        real_state = (*alnit)->STATE_UNKNOWN;
                    /*
                    if (value1 == 0) 
                        real_state = id2;
                    else if (value1 >= N)
                        real_state = id1;
                    else {
                        int j;
                        if (id1 == 0) j = id2 - 1;
                        else j = id1 + id2;
                        real_state = 4 + j*(N-2) + j + value1 - 1;
                    }
                    */
                    state = real_state;
                    ASSERT(state < 4 || state == (*alnit)->STATE_UNKNOWN);
//                    assert(state < nstates);
                }
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
            } // FOR loop
            break; // of switch
            } // end of switch
            start_pos = end_pos;
        } // of end FOR LOOP

        ASSERT(start_pos == aln->ordered_pattern.size());
//        assert(site == aln->num_parsimony_sites % NUM_BITS);
        // add dummy states
        if (site > 0 && site < NUM_BITS) {
            x += site/UINT_BITS;
            if (site % UINT_BITS != 0) {
                *x |= ~((1<<(site%UINT_BITS)) - 1);
                x++;
            }
            int max_sites = ((site+UINT_BITS-1)/UINT_BITS);
            memset(x, 255, (VCSIZE - max_sites)*sizeof(UINT));
        }
        if (!aln->isSuperAlignment())
            delete partitions;
    } else {
        // internal node
        ASSERT(node->degree() == 3); // it works only for strictly bifurcating tree
        PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, pit) {
            if (!pit->isParsimonyComputed()) {
                computePartialParsimonyFastSIMD<VectorClass>(pit, node);
            }
            if (!left) left = pit; else right = pit;
        }
        
        computePartialParsimonyOutOfTreeSIMD<VectorClass>(left->partial_pars,
                                                          right->partial_pars,
                                                          dad_branch->partial_pars);
    }
}

template<class VectorClass>
void PhyloTree::computePartialParsimonyOutOfTreeSIMD(const UINT* left_partial_pars,
                                                     const UINT* right_partial_pars,
                                                     UINT* dad_partial_pars) const {

    const int NUM_BITS   = VectorClass::size() * UINT_BITS;
    int       nstates    = aln->getMaxNumStates();
    UINT      score      = 0;
    size_t    nsites     = (aln->num_parsimony_sites+NUM_BITS-1)/NUM_BITS;
    const int VCSIZE     = VectorClass::size();
    int       entry_size = nstates * VCSIZE;
    
    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>num_threads*10)
        #endif
        for (int site = 0; site<nsites; ++site) {
            size_t  offset = entry_size*site;
            VectorClass* x = (VectorClass*)(left_partial_pars  + offset);
            VectorClass* y = (VectorClass*)(right_partial_pars + offset);
            VectorClass* z = (VectorClass*)(dad_partial_pars   + offset);
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
            score += fast_popcount(w);
        }
        break;
            
    default:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>num_threads*10)
        #endif
        for (int site = 0; site<nsites; ++site) {
            size_t offset = entry_size*site;
            VectorClass* x = (VectorClass*)(left_partial_pars + offset);
            VectorClass* y = (VectorClass*)(right_partial_pars + offset);
            VectorClass* z = (VectorClass*)(dad_partial_pars + offset);
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
            score += fast_popcount(w);
        }
        break;
    }
    auto total = nstates*VCSIZE*nsites;
    dad_partial_pars[total] = score + left_partial_pars[total] + right_partial_pars[total];
}

template<class VectorClass>
int PhyloTree::computeParsimonyBranchFastSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    ASSERT(node_branch);
    if (central_partial_pars==nullptr) {
        initializeAllPartialPars();
    }
    if (!dad_branch->isParsimonyComputed()) {
        computePartialParsimonyFastSIMD<VectorClass>(dad_branch, dad);
    }
    if (!node_branch->isParsimonyComputed()) {
        computePartialParsimonyFastSIMD<VectorClass>(node_branch, node);
    }
    return computeParsimonyOutOfTreeSIMD<VectorClass>(dad_branch->partial_pars,
                                                      node_branch->partial_pars,
                                                      branch_subst);
}

template<class VectorClass>
int PhyloTree::computeParsimonyOutOfTreeSIMD(const UINT* dad_partial_pars,
                                             const UINT* node_partial_pars,
                                             int* branch_subst) const {
    int nstates = aln->getMaxNumStates();

//    VectorClass score = 0;
//    VectorClass w;

    const int NUM_BITS = VectorClass::size() * UINT_BITS;
    int nsites         = (aln->num_parsimony_sites + NUM_BITS - 1)/NUM_BITS;
    int entry_size     = nstates * VectorClass::size();
    
    int  scoreid        = nsites*entry_size;
    UINT sum_end_node  = (dad_partial_pars[scoreid] + node_partial_pars[scoreid]);
    UINT score         = sum_end_node;
    UINT lower_bound   = best_pars_score;
    if (branch_subst) {
        lower_bound = INT_MAX;
    }
    switch (nstates) {
    case 4:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>num_threads*10)
        #endif
        for (int offset = 0; offset < scoreid; offset+=entry_size) {
            VectorClass* x = (VectorClass*)(dad_partial_pars + offset);
            VectorClass* y = (VectorClass*)(node_partial_pars + offset);
            VectorClass  w = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
            w = ~w;
            score += fast_popcount(w);
            #ifndef _OPENMP
            if (score >= lower_bound) 
                break;
            #endif
        }
        break;
    default:
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+: score) if(nsites>num_threads*10)
        #endif
        for (int offset = 0; offset < scoreid; offset+=entry_size) {
            VectorClass *x = (VectorClass*)(dad_partial_pars + offset);
            VectorClass *y = (VectorClass*)(node_partial_pars + offset);
            VectorClass w  = x[0] & y[0];
            for (int i = 1; i < nstates; i++) {
                w |= x[i] & y[i];
            }
            w = ~w;
            score += fast_popcount(w);
            #ifndef _OPENMP
            if (score >= lower_bound) 
                break;
            #endif
        }
        break;
    }
    if (branch_subst) {
        *branch_subst = score - sum_end_node;
    }
    return score;
}

/****************************************************************************
 Sankoff parsimony function
 ****************************************************************************/

template<class VectorClass>
void PhyloTree::computePartialParsimonySankoffSIMD(PhyloNeighbor *dad_branch,
                                                   PhyloNode *dad){
    // don't recompute the parsimony
    if (dad_branch->isParsimonyComputed()) {
        return;
    }
    PhyloNode *node = dad_branch->getNode();
    //assert(node->degree() <= 3);
    /*
     if(aln->num_states != cost_nstates){
     cout << "Your cost matrix is not compatible with the alignment"
     << " in terms of number of states. Please check!" << endl;
     exit(1);
     }
     */
    int nstates = aln->num_states;
    assert(dad_branch->partial_pars);
    
    size_t pars_block_size = getBitsBlockSize();
    
    // internal node
    
    UINT * partial_pars = dad_branch->partial_pars;
    
    PhyloNeighbor *left = nullptr, *right = nullptr;
    
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        PhyloNode* child = nei->getNode();
        if (child->name != ROOT_NAME) {
            if (!child->isLeaf()) {
                computePartialParsimonySankoffSIMD<VectorClass>(nei, node);
            }
            if (!left) {
                left = nei;
            }
            else {
                right = nei;
            }
        }
    }
    dad_branch->setParsimonyComputed(true);
    computeTipPartialParsimony();
    memset(partial_pars, 0, sizeof(UINT)*pars_block_size);
    if (left==nullptr && right==nullptr && 0<=node->id && node->id<aln->getNSeq()) {
        //
        //James B. This calculates a partial parsimony vector oriented
        //         at a leaf (as these are needed during parsimony placement,
        //         when TaxonToPlace's constructor is calculating parsimony
        //         for new_interior->findNeighbor(new_leaf).
        //
        //Note: Instead of "write out-of-order into a tip buffer" and then copy to
        //      partial_pars, this code writes out-of-order into partial_pars_ptr.
        //      (so it can be parallelized).
        //
        int ptnCount = aln->ordered_pattern.size();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t ptn = 0; ptn < ptnCount; ptn+=VectorClass::size()) {
            int ptn_start_index = ptn*nstates;
            VectorClass *partial_pars_ptr = (VectorClass*)&partial_pars[ptn_start_index];
            for (int i = 0; i < VectorClass::size(); i++) {
                UINT*       tip_buffer_ptr         = ((UINT*)partial_pars_ptr) + i;
                size_t      offset                 = aln->ordered_pattern[ptn+i][node->id]*nstates;
                const UINT* partial_pars_child_ptr = &tip_partial_pars[offset];
                for (int j = 0; j < nstates; j++, tip_buffer_ptr += VectorClass::size()) {
                    *tip_buffer_ptr += partial_pars_child_ptr[j];
                }
            }
        }
        return;
    }
    
    if (!left->node->isLeaf() && right->node->isLeaf()) {
        // swap leaf and internal node so that left will be the leaf
        std::swap(left, right);
    }
    
    ASSERT(node->degree() >= 3);
    
    if (node->degree() > 3) {
        // multifurcating node
        int ptnCount = aln->ordered_pattern.size();
        #ifdef _OPENMP //Can now be parallelized, because tip buffer no longer used.
        #pragma omp parallel for
        #endif
        for (size_t ptn = 0; ptn < ptnCount; ptn+=VectorClass::size()) {
            int ptn_start_index = ptn*nstates;
            VectorClass *partial_pars_ptr = (VectorClass*)&partial_pars[ptn_start_index];
            FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) if ((*it)->node->name != ROOT_NAME) {
                PhyloNode* child = nei->getNode();
                if (child->isLeaf()) {
                    // leaf node
                    // Note: James B. Rewrote this 18-Sep-2020, so that it
                    //       doesn't have to use a tip buffer.
                    for (int i = 0; i < VectorClass::size(); i++) {
                        UINT* partial_pars_child_ptr = &tip_partial_pars[aln->ordered_pattern[ptn+i][child->id]*nstates];
                        UINT* tip_buffer_ptr         = ((UINT*)partial_pars_ptr) + i;
                        for (int j = 0; j < nstates; ++j, tip_buffer_ptr += VectorClass::size()) {
                            *tip_buffer_ptr += partial_pars_child_ptr[j];
                        }
                    }
                } else {
                    // internal node
                    VectorClass* partial_pars_child_ptr = (VectorClass*) & nei->partial_pars[ptn_start_index];
                    UINT*        cost_matrix_ptr        = cost_matrix;
                    for (int i = 0; i < nstates; ++i){
                        // min(j->i) from child_branch
                        VectorClass min_child_ptn_pars = partial_pars_child_ptr[0] + cost_matrix_ptr[0];
                        for (int j = 1; j < nstates; j++) {
                            min_child_ptn_pars = min(partial_pars_child_ptr[j] + cost_matrix_ptr[j], min_child_ptn_pars);
                        }
                        partial_pars_ptr[i] += min_child_ptn_pars;
                        cost_matrix_ptr += nstates;
                    }
                }
            }
        }
        return;
    }
    
    if (left->node->isLeaf() && right->node->isLeaf()) {
        // tip-tip case
        // Note: James B. Rewrote this 18-Sep-2020, so that it
        //       doesn't use a tip buffer (and can be parallelized).
        size_t ptnCount = aln->ordered_pattern.size();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t ptn = 0; ptn < ptnCount; ptn+=VectorClass::size()){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            size_t       ptn_start_index  = ptn*nstates;
            VectorClass* partial_pars_ptr = (VectorClass*)&partial_pars[ptn_start_index];

            // load data for tip
            for (int i = 0; i < VectorClass::size(); i++) {
                UINT* left_ptr             = &tip_partial_pars[aln->ordered_pattern[ptn+i][left->node->id]*nstates];
                UINT* right_ptr            = &tip_partial_pars[aln->ordered_pattern[ptn+i][right->node->id]*nstates];
                UINT* tip_buffer_ptr       = ((UINT*)partial_pars_ptr) + i;
                for (int j = 0; j < nstates; j++) {
                    *tip_buffer_ptr       += (left_ptr[j] + right_ptr[j]);
                    tip_buffer_ptr        += VectorClass::size();
                }
            }
        }
        return;
    }
    
    if (left->node->isLeaf() && !right->node->isLeaf()) {
        // tip-inner case
        // Note: this still needs to use a tip buffer, but that's created inside
        //       the loop, so the loop can be parallelized.
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t ptn = 0; ptn < aln->ordered_pattern.size(); ptn+=VectorClass::size()){
            // ignore const ptn because it does not affect pars score
            //if (aln->at(ptn).isConst()) continue;
            size_t      ptn_start_index = ptn*nstates;
            VectorClass tip_buffer[nstates];
            
            for (int i = 0; i < VectorClass::size(); i++) {
                const UINT *left_ptr = &tip_partial_pars[aln->ordered_pattern[ptn+i][left->node->id]*nstates];
                UINT *tip_buffer_ptr = (UINT*)tip_buffer + i;
                for (int j = 0; j < nstates; j++) {
                    *tip_buffer_ptr = left_ptr[j];
                    tip_buffer_ptr += VectorClass::size();
                }
            }
            
            const VectorClass* right_ptr        = (VectorClass*)&right->partial_pars[ptn_start_index];
            VectorClass*       partial_pars_ptr = (VectorClass*)&partial_pars[ptn_start_index];
            const UINT*        cost_matrix_ptr  = cost_matrix;
            VectorClass right_contrib;
            
            for(int i = 0; i < nstates; i++){
                // min(j->i) from child_branch
                right_contrib = right_ptr[0] + cost_matrix_ptr[0];
                for(int j = 1; j < nstates; j++) {
                    right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
                }
                partial_pars_ptr[i] = tip_buffer[i] + right_contrib;
                cost_matrix_ptr += nstates;
            }
        }
        return;
    }
    
    // inner-inner case
    computePartialParsimonyOutOfTreeSankoffSIMD<VectorClass>
        ( left->partial_pars, right->partial_pars, partial_pars );
}

template<class VectorClass>
void PhyloTree::computePartialParsimonyOutOfTreeSankoffSIMD
        (const UINT* left_partial_pars, const UINT* right_partial_pars,
         UINT*       dad_partial_pars) const
{
    size_t ptnCount = aln->ordered_pattern.size();
    size_t ptnStep  = VectorClass::size();
    int nstates = aln->num_states;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t ptn = 0; ptn < ptnCount; ptn+=ptnStep){
        // ignore const ptn because it does not affect pars score
        //if (aln->at(ptn).isConst()) continue;
        int          ptn_start_index  = ptn*nstates;
        VectorClass* left_ptr         = (VectorClass*)&left_partial_pars[ptn_start_index];
        VectorClass* right_ptr        = (VectorClass*)&right_partial_pars[ptn_start_index];
        VectorClass* partial_pars_ptr = (VectorClass*)&dad_partial_pars[ptn_start_index];
        UINT *cost_matrix_ptr         = cost_matrix;
        VectorClass left_contrib, right_contrib;
        
        for (int i = 0; i < nstates; i++){
            // min(j->i) from child_branch
            left_contrib  =  left_ptr[0] + cost_matrix_ptr[0];
            right_contrib = right_ptr[0] + cost_matrix_ptr[0];
            for (int j = 1; j < nstates; j++) {
                left_contrib  = min( left_ptr[j] + cost_matrix_ptr[j],  left_contrib);
                right_contrib = min(right_ptr[j] + cost_matrix_ptr[j], right_contrib);
            }
            partial_pars_ptr[i] = left_contrib + right_contrib;
            cost_matrix_ptr    += nstates;
        }
    }
}

template<class VectorClass>
int PhyloTree::computeParsimonyBranchSankoffSIMD(PhyloNeighbor *dad_branch,
                                                 PhyloNode *dad, int *branch_subst) {
    if ((tip_partial_lh_computed & 2) == 0) {
        computeTipPartialParsimony();
    }
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    assert(node_branch);
    
    if (!central_partial_pars) {
        initializeAllPartialPars();
    }
    
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    
    //int nptn = aln->size();
    //    if(!_pattern_pars) _pattern_pars = aligned_alloc<BootValTypePars>(nptn+VCSIZE_USHORT);
    //    memset(_pattern_pars, 0, sizeof(BootValTypePars) * (nptn+VCSIZE_USHORT));
    
    if (!dad_branch->isParsimonyComputed() && !node->isLeaf()) {
        computePartialParsimonySankoffSIMD<VectorClass>(dad_branch, dad);
    }
    
    if (!node_branch->isParsimonyComputed() && !dad->isLeaf()) {
        computePartialParsimonySankoffSIMD<VectorClass>(node_branch, node);
    }
    
    // now combine likelihood at the branch
    int nstates = aln->num_states;
    
    if (dad->isLeaf()) {
        VectorClass* tip_buffer = aligned_alloc<VectorClass>(nstates);
        VectorClass  tree_pars   = 0;
        VectorClass  branch_pars = 0;
        // external node
        for (size_t ptn = 0; ptn < aln->ordered_pattern.size(); ptn+=VectorClass::size()){
            int ptn_start_index = ptn * nstates;
            for (int  i = 0; i < VectorClass::size(); i++) {
                UINT *node_branch_ptr = &tip_partial_pars[aln->ordered_pattern[ptn+i][dad->id]*nstates];
                UINT *tip_buffer_ptr = (UINT*)tip_buffer + i;
                for (int j = 0; j < nstates; j++, tip_buffer_ptr += VectorClass::size()) {
                    *tip_buffer_ptr = node_branch_ptr[j];
                }
            }
            VectorClass* dad_branch_ptr = (VectorClass*)&dad_branch->partial_pars[ptn_start_index];
            VectorClass  min_ptn_pars   = tip_buffer[0] + dad_branch_ptr[0];
            VectorClass  br_ptn_pars    = tip_buffer[0];
            for (int i = 1; i < nstates; i++){
                // min(j->i) from node_branch
                VectorClass min_score = tip_buffer[i] + dad_branch_ptr[i];
                br_ptn_pars  = select(min_score < min_ptn_pars, tip_buffer[i], br_ptn_pars);
                min_ptn_pars = min(min_ptn_pars, min_score);
            }
            //_pattern_pars[ptn] = min_ptn_pars;
            tree_pars   += min_ptn_pars * VectorClass().load_a(&ptn_freq_pars[ptn]);
            branch_pars += br_ptn_pars  * VectorClass().load_a(&ptn_freq_pars[ptn]);
        }
        aligned_free(tip_buffer);
        if (branch_subst != nullptr) {
            *branch_subst = horizontal_add(branch_pars);
        }
        return horizontal_add(tree_pars);
    }  else {
        // internal node
        return computeParsimonyOutOfTreeSankoffSIMD<VectorClass>
               ( dad_branch->partial_pars, node_branch->partial_pars, branch_subst);
    }
}

template<class VectorClass>
int PhyloTree::computeParsimonyOutOfTreeSankoffSIMD(const UINT* dad_partial_pars,
                                                    const UINT* node_partial_pars,
                                                    int*        branch_subst) const {
    size_t ptnCount    = aln->ordered_pattern.size();
    size_t ptnStep     = VectorClass::size();
    int    nstates     = aln->num_states;
    UINT   tree_pars   = 0; //James B. tree_pars and branch_pars now UINT
    UINT   branch_pars = 0; //so that the for-loop can be parallelized.

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:tree_pars,branch_pars)
    #endif
    for (size_t ptn = 0; ptn < ptnCount; ptn += ptnStep){
        int          ptn_start_index = ptn * nstates;
        VectorClass* node_branch_ptr = (VectorClass*)&node_partial_pars[ptn_start_index];
        VectorClass* dad_branch_ptr  = (VectorClass*)&dad_partial_pars [ptn_start_index];
        UINT*        cost_matrix_ptr = cost_matrix;
        VectorClass  min_ptn_pars    = UINT_MAX;
        VectorClass  br_ptn_pars     = UINT_MAX;
        for(int i = 0; i < nstates; ++i){
            // min(j->i) from node_branch
            VectorClass min_score    = node_branch_ptr[0] + cost_matrix_ptr[0];
            VectorClass branch_score = cost_matrix_ptr[0];
            for(int j = 1; j < nstates; ++j) {
                VectorClass value = node_branch_ptr[j] + cost_matrix_ptr[j];
                branch_score      = select(value < min_score, cost_matrix_ptr[j], branch_score);
                min_score         = min(value, min_score);
            }
            min_score       += dad_branch_ptr[i];
            br_ptn_pars      = select(min_score < min_ptn_pars, branch_score, br_ptn_pars);
            min_ptn_pars     = min(min_score, min_ptn_pars);
            cost_matrix_ptr += nstates;
        }
        //_pattern_pars[ptn] = min_ptn_pars;
        tree_pars   += horizontal_add(min_ptn_pars * VectorClass().load_a(&ptn_freq_pars[ptn]));
        branch_pars += horizontal_add(br_ptn_pars  * VectorClass().load_a(&ptn_freq_pars[ptn]));
    }
    if (branch_subst != nullptr) {
        *branch_subst = branch_pars;
    }
    return tree_pars;
}

#endif /* PHYLOKERNEL_H_ */
