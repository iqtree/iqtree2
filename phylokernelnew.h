/*
 * phylokernelnew.h
 * Newly revised kernel based on vectorizing over alignment patterns
 *
 *  Created on: Sept 23, 2016
 *      Author: minh
 */


#if !defined(PHYLOKERNELNEW_H_) || !defined(PHYLOKERNELNEW_STATE_H_)

#ifdef KERNEL_FIX_STATES
#   define PHYLOKERNELNEW_STATE_H_
#else
#   define PHYLOKERNELNEW_H_
#endif


#include "phylotree.h"
#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_exp.h"
#include "vectorf64.h"



/*******************************************************
 *
 * NEW! highly-vectorized likelihood functions.
 *
 ******************************************************/


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates>
void PhyloTree::computePartialLikelihoodNewSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#else
template <class VectorClass, const bool SAFE_NUMERIC>
void PhyloTree::computePartialLikelihoodSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#endif
{

    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    PhyloNode *node = (PhyloNode*)(dad_branch->node);


#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    if (!tip_partial_lh_computed)
        computeTipPartialLikelihood();

	if (node->isLeaf()) {
	    dad_branch->lh_scale_factor = 0.0;
		return;
	}
    
    size_t ptn, c;
    size_t orig_ntn = aln->size();
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }
    size_t i, x;
    size_t block = nstates * ncat_mix;
    size_t tip_block = nstates * model->getNMixtures();
    size_t scale_size = nptn * ncat_mix;
    
	double *evec = model->getEigenvectors();
	double *inv_evec = model->getInverseEigenvectors();
	assert(inv_evec && evec);
	double *eval = model->getEigenvalues();

    dad_branch->lh_scale_factor = 0.0;

	// internal node
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
    int num_leaves = 0;
	FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)*it;
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
        if ((nei->partial_lh_computed & 1) == 0)
            computePartialLikelihood(nei, node);
        dad_branch->lh_scale_factor += nei->lh_scale_factor;
        if (nei->node->isLeaf())
            num_leaves ++;
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

    // precompute buffer to save times
    double *echildren = aligned_alloc<double>(block*nstates*(node->degree()-1));
    double *partial_lh_leaves = NULL;
    if (num_leaves > 0)
        partial_lh_leaves = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block*num_leaves);
    double *echild = echildren;
    double *partial_lh_leaf = partial_lh_leaves;

    FOR_NEIGHBOR_IT(node, dad, it) {
        double expchild[nstates];
        PhyloNeighbor *child = (PhyloNeighbor*)*it;
        // precompute information buffer
        double *echild_ptr = echild;
        for (c = 0; c < ncat_mix; c++) {
            double len_child = site_rate->getRate(c%ncat) * child->length;
            double *eval_ptr = eval + mix_addr_nstates[c];
            double *evec_ptr = evec + mix_addr[c];
            for (i = 0; i < nstates; i++) {
                expchild[i] = exp(eval_ptr[i]*len_child);
            }
            for (x = 0; x < nstates; x++) {
                for (i = 0; i < nstates; i++) {
                    echild_ptr[i] = evec_ptr[x*nstates+i] * expchild[i];
                }
                echild_ptr += nstates;
            }
        }

        // pre compute information for tip
        if (child->node->isLeaf()) {
            vector<int>::iterator it;
            for (it = aln->seq_states[child->node->id].begin(); it != aln->seq_states[child->node->id].end(); it++) {
                int state = (*it);
                double *this_partial_lh_leaf = partial_lh_leaf + state*block;
                double *echild_ptr = echild;
                for (c = 0; c < ncat_mix; c++) {
                    double *this_tip_partial_lh = tip_partial_lh + state*tip_block + mix_addr_nstates[c];
                    for (x = 0; x < nstates; x++) {
                        double vchild = 0.0;
                        for (i = 0; i < nstates; i++) {
                            vchild += echild_ptr[i] * this_tip_partial_lh[i];
                        }
                        this_partial_lh_leaf[x] = vchild;
                        echild_ptr += nstates;
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
        echild += block*nstates;
    }
    
    
    double *eleft = echildren, *eright = echildren + block*nstates;
    
	if (!left->node->isLeaf() && right->node->isLeaf()) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
        double *etmp = eleft;
        eleft = eright;
        eright = etmp;
	}
    
    if (node->degree() > 3) {
        /*--------------------- multifurcating node ------------------*/

        double *partial_lh_all_dbl = aligned_alloc<double>(block*VectorClass::size()*2);
        double *vec_tip = partial_lh_all_dbl + block*VectorClass::size();

        // now for-loop computing partial_lh over all site-patterns
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
        for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
            VectorClass *partial_lh_all = (VectorClass*)partial_lh_all_dbl;
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
//            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
//            memset(scale_dad, 0, sizeof(UBYTE)*ncat_mix);

            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = echildren;

            FOR_NEIGHBOR_IT(node, dad, it) {
                PhyloNeighbor *child = (PhyloNeighbor*)*it;
//                UBYTE *scale_child = child->scale_num + ptn*ncat_mix;
                if (child->node->isLeaf()) {
                    // external node
                    // load data for tip
                    for (i = 0; i < VectorClass::size(); i++) {
                        double *child_lh;
                        if (ptn+i < orig_ntn)
                            child_lh = partial_lh_leaf + block*(aln->at(ptn+i))[child->node->id];
                        else if (ptn+i < nptn)
                            child_lh = partial_lh_leaf + block*model_factory->unobserved_ptns[ptn+i-orig_ntn];
                        else
                            child_lh = partial_lh_leaf + block*aln->STATE_UNKNOWN;
                        double *this_vec_tip = vec_tip+i;
                        for (c = 0; c < block; c++)
                            this_vec_tip[c*VectorClass::size()] = child_lh[c];
                    }
                    VectorClass *vtip = (VectorClass*)vec_tip;
                    for (c = 0; c < block; c++) {
                        // compute real partial likelihood vector
                        partial_lh_all[c] *= vtip[c];
                    }
                    partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                } else {
                    // internal node
                    VectorClass *partial_lh = partial_lh_all;
                    VectorClass *partial_lh_child = (VectorClass*)(child->partial_lh + ptn*block);

                    double *echild_ptr = echild;
                    for (c = 0; c < ncat_mix; c++) {
//                        scale_dad[c] += scale_child[c];
                        // compute real partial likelihood vector
                        for (x = 0; x < nstates; x++) {
                            VectorClass vchild = echild_ptr[0] * partial_lh_child[0];
//                            double *echild_ptr = echild + (c*nstatesqr+x*nstates);
                            for (i = 1; i < nstates; i++) {
                                vchild = mul_add(echild_ptr[i], partial_lh_child[i], vchild);
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
            VectorClass *partial_lh_tmp = partial_lh_all;
            VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            for (c = 0; c < ncat_mix; c++) {
                VectorClass lh_max = 0.0;
                double *inv_evec_ptr = inv_evec + mix_addr[c];
                for (i = 0; i < nstates; i++) {
					VectorClass res = partial_lh_tmp[0]*inv_evec_ptr[0];
					for (x = 1; x < nstates; x++) {
						res = mul_add(partial_lh_tmp[x], inv_evec_ptr[x], res);
					}
                    inv_evec_ptr += nstates;
                    partial_lh[i] = res;
                    res = abs(res);
                    lh_max = max(lh_max, res);
                }
                // check if one should scale partial likelihoods
//                if (lh_max < SCALING_THRESHOLD && lh_max != 0.0) {
//                    //assert(lh_max != 0.0 && "Numerical underflow for multifurcation node");
//                    if (ptn_invar[ptn] == 0.0) {
//                        // now do the likelihood scaling
//                        for (i = 0; i < nstates; i++)
//                            partial_lh[i] *= SCALING_THRESHOLD_INVER;
//                        scale_dad[c] += 1;
//                    }
//                }
                partial_lh += nstates;
                partial_lh_tmp += nstates;
            }

        } // for ptn
//        dad_branch->lh_scale_factor += sum_scale;               
        aligned_free(partial_lh_all_dbl);

        // end multifurcating treatment
    } else if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *partial_lh_left = partial_lh_leaves;
        double *partial_lh_right = partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;
        // TODO: this is not thread-safe
        double *vec_left = aligned_alloc<double>(block*VectorClass::size()*2 + nstates*VectorClass::size());
        double *vec_right = vec_left + block*VectorClass::size();
        double *partial_lh_dbl = vec_right + block*VectorClass::size();

		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, scale_size * sizeof(UBYTE));
#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass* partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            VectorClass *vright = (VectorClass*)vec_right;

            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip_left, *tip_right;
                if (ptn+x < orig_ntn) {
                    tip_left  = partial_lh_left  + block * (aln->at(ptn+x))[left->node->id];
                    tip_right = partial_lh_right + block * (aln->at(ptn+x))[right->node->id];
                } else if (ptn+x < nptn) {
                    tip_left  = partial_lh_left  + block * model_factory->unobserved_ptns[ptn+x-orig_ntn];
                    tip_right = partial_lh_right + block * model_factory->unobserved_ptns[ptn+x-orig_ntn];
                } else {
                    tip_left  = partial_lh_left  + block * aln->STATE_UNKNOWN;
                    tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                }
                double *this_vec_left = vec_left+x;
                double *this_vec_right = vec_right+x;
                for (i = 0; i < block; i++) {
                    this_vec_left[i*VectorClass::size()] = tip_left[i];
                    this_vec_right[i*VectorClass::size()] = tip_right[i];
                }
            }


			for (c = 0; c < ncat_mix; c++) {
                double *inv_evec_ptr = inv_evec + mix_addr[c];
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					partial_lh_tmp[x] = vleft[x] * vright[x];
				}

				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					VectorClass res = partial_lh_tmp[0]*inv_evec_ptr[0];
					for (x = 1; x < nstates; x++) {
						res = mul_add(partial_lh_tmp[x], inv_evec_ptr[x], res);
					}
                    inv_evec_ptr += nstates;
					partial_lh[i] = res;
				}

                // increase pointer
                vleft += nstates;
                vright += nstates;
                partial_lh += nstates;
			}
		}

        aligned_free(vec_left);

	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num, right->scale_num, scale_size * sizeof(UBYTE));


        double *partial_lh_left = partial_lh_leaves;

        // TODO: this is not thread-safe
        double *vec_left = aligned_alloc<double>(block*VectorClass::size() + nstates*VectorClass::size());
        double *partial_lh_dbl = vec_left + block*VectorClass::size();


#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass *partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip;
                if (ptn+x < orig_ntn) {
                    tip = partial_lh_left + block*(aln->at(ptn+x))[left->node->id];
                } else if (ptn+x < nptn) {
                    tip = partial_lh_left + block*model_factory->unobserved_ptns[ptn+x-orig_ntn];
                } else {
                    tip = partial_lh_left + block*aln->STATE_UNKNOWN;
                }
                double *this_vec_left = vec_left+x;
                for (i = 0; i < block; i++)
                    this_vec_left[i*VectorClass::size()] = tip[i];
            }

            double *eright_ptr = eright;
			for (c = 0; c < ncat_mix; c++) {
                VectorClass lh_max = 0.0;
                double *inv_evec_ptr = inv_evec + mix_addr[c];
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					VectorClass vright = VectorClass(eright_ptr[0]) * partial_lh_right[0];
					for (i = 1; i < nstates; i++) {
						vright = mul_add(VectorClass(eright_ptr[i]), partial_lh_right[i], vright);
					}
                    eright_ptr += nstates;
					partial_lh_tmp[x] = vleft[x] * (vright);
				}
                
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					VectorClass res = partial_lh_tmp[0]*inv_evec_ptr[0];
					for (x = 1; x < nstates; x++) {
						res = mul_add(partial_lh_tmp[x], inv_evec_ptr[x], res);
					}
                    inv_evec_ptr += nstates;
					partial_lh[i] = res;
                    res = abs(res);
                    lh_max = max(lh_max, res);
				}
                // check if one should scale partial likelihoods
//                auto small_lh = ((lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(ptn_invar+ptn) == 0.0));
//                for (x = 0; x < VectorClass::size(); x++)
//                if (small_lh[x]) {
//                    // this only happens occasionally, no need for vectorization
//                    // BQM 2016-05-03: only scale for non-constant sites
//                    // now do the likelihood scaling
//                    double *this_partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
//                    for (i = 0; i < nstates; i++)
//                        this_partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
//                    dad_branch->scale_num[(ptn+x)*ncat_mix+c] += 1;
//                }
                vleft += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
			}

		}
        aligned_free(vec_left);
//		dad_branch->lh_scale_factor += sum_scale;
//		delete [] partial_lh_left;

	} else {

        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/
        double *partial_lh_dbl = aligned_alloc<double>(nstates*VectorClass::size());

#ifdef _OPENMP
#pragma omp parallel for private(ptn, c, x, i) schedule(static)
#endif
		for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass *partial_lh_tmp = (VectorClass*)partial_lh_dbl;
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(left->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
//            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
//            UBYTE *scale_left = left->scale_num + ptn*ncat_mix;
//            UBYTE *scale_right = right->scale_num + ptn*ncat_mix; 

            double *eleft_ptr = eleft;
            double *eright_ptr = eright;

			for (c = 0; c < ncat_mix; c++) {
//                for (x = 0; x < VectorClass::size(); x++)
//                    scale_dad[x*ncat_mix+c] = scale_left[x*ncat_mix+c] + scale_right[x*ncat_mix+c];
                VectorClass lh_max = 0.0;
                double *inv_evec_ptr = inv_evec + mix_addr[c];
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					VectorClass vleft =  eleft_ptr[0]  * partial_lh_left[0];
                    VectorClass vright = eright_ptr[0] * partial_lh_right[0];
//					size_t addr = c*nstatesqr+x*nstates;
					for (i = 1; i < nstates; i++) {
						vleft  = mul_add(eleft_ptr[i],  partial_lh_left[i],  vleft);
						vright = mul_add(eright_ptr[i], partial_lh_right[i], vright);
					}
                    eleft_ptr += nstates;
                    eright_ptr += nstates;
					partial_lh_tmp[x] = vleft*vright;
//                    assert(partial_lh_tmp[x] != 0.0);
				}
                
				// compute dot-product with inv_eigenvector
				for (i = 0; i < nstates; i++) {
					VectorClass res = partial_lh_tmp[0]*inv_evec_ptr[0];
					for (x = 1; x < nstates; x++) {
						res = mul_add(partial_lh_tmp[x], inv_evec_ptr[x], res);
					}
                    inv_evec_ptr += nstates;
					partial_lh[i] = res;
                    res = abs(res);
                    lh_max = max(lh_max, res);
				}
                // check if one should scale partial likelihoods
//                auto small_lh = ((lh_max < SCALING_THRESHOLD) & (lh_max != 0.0) & (VectorClass().load_a(ptn_invar+ptn) == 0.0));
//                for (x = 0; x < VectorClass::size(); x++)
//                if (small_lh[x]) {
//                    // this only happens occasionally, no need for vectorization
//                    // BQM 2016-05-03: only scale for non-constant sites
//                    // now do the likelihood scaling
//                    double *this_partial_lh = dad_branch->partial_lh + (ptn*block + c*nstates*VectorClass::size() + x);
//                    for (i = 0; i < nstates; i++)
//                        this_partial_lh[i*VectorClass::size()] *= SCALING_THRESHOLD_INVER;
//                    scale_dad[x*ncat_mix+c] += 1;
//                }
                partial_lh_left += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
			}

		}
//		dad_branch->lh_scale_factor += sum_scale;

        aligned_free(partial_lh_dbl);

	}

    if (partial_lh_leaves)
        aligned_free(partial_lh_leaves);
    aligned_free(echildren);
}


#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates>
void PhyloTree::computeLikelihoodDervNewSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
#else
template <class VectorClass, const bool SAFE_NUMERIC>
void PhyloTree::computeLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf)
#endif
{
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

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    assert(eval);

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case

            double *vec_tip = aligned_alloc<double>(tip_block*VectorClass::size());

#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
				VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
//                UBYTE *scale_dad = dad_branch->scale_num+ptn*ncat_mix;
				VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
                //load tip vector
                for (i = 0; i < VectorClass::size(); i++) {
                    double *this_tip_partial_lh = tip_partial_lh + tip_block*((ptn+i < orig_nptn) ? (aln->at(ptn+i))[dad->id] :  aln->STATE_UNKNOWN);
//                double *this_tip_partial_lh = tip_partial_lh + tip_block*((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]);
                    double *this_vec_tip = vec_tip+i;
                    for (c = 0; c < tip_block; c++)
                        this_vec_tip[c*VectorClass::size()] = this_tip_partial_lh[c];

                }
//                UBYTE min_scale = scale_dad[0];
//                for (c = 1; c < ncat_mix; c++)
//                    min_scale = min(min_scale, scale_dad[c]);
                for (c = 0; c < ncat_mix; c++) {
                    VectorClass *lh_tip = (VectorClass*)(vec_tip + mix_addr_nstates[c]*VectorClass::size());
//                    if (scale_dad[c] == min_scale) {
                        for (i = 0; i < nstates; i++) {
                            theta[i] = lh_tip[i] * partial_lh_dad[i];
                        }
//                    } else if (scale_dad[c] == min_scale+1) {
//                        for (i = 0; i < nstates; i++) {
//                            theta[i] = lh_tip[i] * partial_lh_dad[i] * SCALING_THRESHOLD;
//                        }
//                    } else {
//                        memset(theta, 0, sizeof(double)*nstates);
//                    }
                    partial_lh_dad += nstates;
                    theta += nstates;
                }

			}
            aligned_free(vec_tip);
			// ascertainment bias correction
	    } else {
	    	// both dad and node are internal nodes

//	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static)
#endif
	    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
				VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
			    VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
			    VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);

//                size_t ptn_ncat = ptn*ncat_mix; 
//                UBYTE *scale_dad = dad_branch->scale_num + ptn_ncat;
//                UBYTE *scale_node = node_branch->scale_num + ptn_ncat;
//                UBYTE sum_scale[ncat_mix];
//                UBYTE min_scale = sum_scale[0] = scale_dad[0] + scale_node[0];
//                for (c = 1; c < ncat_mix; c++) {
//                    sum_scale[c] = scale_dad[c] + scale_node[c];
//                    min_scale = min(min_scale, sum_scale[c]);
//                }
                for (c = 0; c < ncat_mix; c++) {
//                    if (sum_scale[c] == min_scale) {
                        for (i = 0; i < nstates; i++) {
                            theta[i] = partial_lh_node[i] * partial_lh_dad[i];
                        }
//                    } else if (sum_scale[c] == min_scale+1) {
//                        for (i = 0; i < nstates; i++) {
//                            theta[i] = partial_lh_node[i] * partial_lh_dad[i] * SCALING_THRESHOLD;
//                        }
//                    } else {
//                        memset(theta, 0, sizeof(double)*nstates);
//                    }
                    theta += nstates;
                    partial_lh_dad += nstates;
                    partial_lh_node += nstates;
                }
			}
	    }
		theta_computed = true;
	}

    double *val0 = aligned_alloc<double>(block);
    double *val1 = aligned_alloc<double>(block);
    double *val2 = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        double *eval_ptr = eval + mix_addr_nstates[c];
        size_t mycat = c%ncat;
        double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        size_t addr = c*nstates;
		for (i = 0; i < nstates; i++) {
			double cof = eval_ptr[i]*site_rate->getRate(mycat);
			double val = exp(cof*dad_branch->length) * prop;
			double val1_ = cof*val;
			val0[addr+i] = val;
			val1[addr+i] = val1_;
			val2[addr+i] = cof*val1_;
		}
	}


    VectorClass my_df = 0.0, my_ddf = 0.0, vc_prob_const = 0.0, vc_df_const = 0.0, vc_ddf_const = 0.0;
    double prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
		VectorClass lh_ptn;
        lh_ptn.load_a(&ptn_invar[ptn]);
        VectorClass df_ptn = 0.0;
        VectorClass ddf_ptn = 0.0;
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
		for (i = 0; i < block; i++) {
			lh_ptn  = mul_add(val0[i], theta[i], lh_ptn);
			df_ptn  = mul_add(val1[i], theta[i], df_ptn);
			ddf_ptn = mul_add(val2[i], theta[i], ddf_ptn);
		}

//        assert(lh_ptn > 0.0);
        lh_ptn = abs(lh_ptn);
        
        if (ptn < orig_nptn) {
            lh_ptn = 1.0 / lh_ptn;
			VectorClass df_frac = df_ptn * lh_ptn;
			VectorClass ddf_frac = ddf_ptn * lh_ptn;
			VectorClass freq;
            freq.load_a(&ptn_freq[ptn]);
			VectorClass tmp1 = df_frac * freq;
			VectorClass tmp2 = ddf_frac * freq;
			my_df += tmp1;
//			my_ddf += tmp2 - tmp1 * df_frac;
            my_ddf += nmul_add(tmp1, df_frac, tmp2);
		} else {
			// ascertainment bias correction
			vc_prob_const += lh_ptn;
			vc_df_const += df_ptn;
			vc_ddf_const += ddf_ptn;
		}
    }
	df = horizontal_add(my_df);
	ddf = horizontal_add(my_ddf);
    prob_const = horizontal_add(vc_prob_const);
    df_const = horizontal_add(vc_df_const);
    ddf_const = horizontal_add(vc_ddf_const);
    
    assert(!isnan(df) && !isinf(df) && "Numerical underflow for lh-derivative");

    if (isnan(df) || isinf(df)) {
        df = 0.0;
        ddf = 0.0;
//        outWarning("Numerical instability (some site-likelihood = 0)");
    }

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
    }


    aligned_free(val2);
    aligned_free(val1);
    aligned_free(val0);
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates>
double PhyloTree::computeLikelihoodBranchNewSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#else
template <class VectorClass, const bool SAFE_NUMERIC>
double PhyloTree::computeLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad)
#endif
{
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
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    VectorClass vc_tree_lh = 0.0;
#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    double *eval = model->getEigenvalues();
    assert(eval);

    double *val = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t mycat = c%ncat;
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
        double *eval_ptr = eval + mix_addr_nstates[c];
		double len = site_rate->getRate(mycat)*dad_branch->length;
		double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        double *this_val = val + c*nstates;
		for (i = 0; i < nstates; i++)
			this_val[i] = exp(eval_ptr[i]*len) * prop;
	}

	double prob_const = 0.0;
    VectorClass vc_prob_const = 0.0;
	memset(_pattern_lh_cat, 0, sizeof(double)*nptn*ncat_mix);

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
    	double *partial_lh_node = aligned_alloc<double>((aln->STATE_UNKNOWN+1)*block);
    	IntVector states_dad = aln->seq_states[dad->id];
    	states_dad.push_back(aln->STATE_UNKNOWN);
    	// precompute information from one tip
    	for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
    		double *lh_node = partial_lh_node +(*it)*block;
    		double *val_tmp = val;
            double *this_tip_partial_lh = tip_partial_lh + (*it)*tip_block;
			for (c = 0; c < ncat_mix; c++) {
                double *lh_tip = this_tip_partial_lh + mix_addr_nstates[c];
				for (i = 0; i < nstates; i++) {
					  lh_node[i] = val_tmp[i] * lh_tip[i];
				}
				lh_node += nstates;
				val_tmp += nstates;
			}
    	}

        double *vec_tip = aligned_alloc<double>(block*VectorClass::size());

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass lh_ptn;
            lh_ptn.load_a(&ptn_invar[ptn]);
            VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
//            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
            VectorClass *lh_node = (VectorClass*)vec_tip;

            //load tip vector
            for (i = 0; i < VectorClass::size(); i++) {
                double *lh_tip = partial_lh_node + block*(((ptn+i) < orig_nptn) ? (aln->at(ptn+i))[dad->id] : aln->STATE_UNKNOWN);
//                double *lh_tip = partial_lh_node + block*(((ptn+i) < orig_nptn) ? (aln->at(ptn+i))[dad->id] : model_factory->unobserved_ptns[ptn+i-orig_nptn]);
                double *this_vec_tip = vec_tip+i;
                for (c = 0; c < block; c++)
                    this_vec_tip[c*VectorClass::size()] = lh_tip[c];

            }
            // determine the min scaling
//            UBYTE min_scale = scale_dad[0];
//            for (c = 1; c < ncat_mix; c++) 
//                min_scale = min(min_scale, scale_dad[c]);

            for (c = 0; c < ncat_mix; c++) {
//                if (scale_dad[c] <= min_scale+1) {
                    // only compute for least scale category
                    *lh_cat = (lh_node[0] * partial_lh_dad[0]);
                    for (i = 1; i < nstates; i++) {
                        *lh_cat = mul_add(lh_node[i], partial_lh_dad[i], *lh_cat);
                    }
//                    if (scale_dad[c] != min_scale)
//                        *lh_cat *= SCALING_THRESHOLD;
                    lh_ptn += *lh_cat;
//                }
                lh_node += nstates;
                partial_lh_dad += nstates;
                lh_cat++;
            }
//			assert(lh_ptn > -1e-10);
			if (ptn < orig_nptn) {
				//lh_ptn = log(abs(lh_ptn)) + LOG_SCALING_THRESHOLD*min_scale;
                lh_ptn = log(abs(lh_ptn));
				lh_ptn.store_a(&_pattern_lh[ptn]);
				vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
			} else {
                // bugfix 2016-01-21, prob_const can be rescaled
//                if (min_scale >= 1)
//                    lh_ptn *= SCALING_THRESHOLD;
//				_pattern_lh[ptn] = lh_ptn;
				vc_prob_const += lh_ptn;
			}
		}
        aligned_free(vec_tip);
		aligned_free(partial_lh_node);
    } else {
    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, i, c) schedule(static)
#endif
    	for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
			VectorClass lh_ptn;
            lh_ptn.load_a(&ptn_invar[ptn]);
            VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
            VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
            double *val_tmp = val;
//            UBYTE *scale_dad = dad_branch->scale_num + ptn*ncat_mix;
//            UBYTE *scale_node = node_branch->scale_num + ptn*ncat_mix;
//            UBYTE sum_scale[ncat_mix*VectorClass::size()];
//            UBYTE min_scale[VectorClass::size()];
//            VectorClass vc_min_scale;
//
//            for (i = 0; i < VectorClass::size(); i++) {
//                min_scale[i] = sum_scale[i] = scale_dad[i] + scale_node[i];
//                for (c = 1; c < ncat_mix; c++) {
//                    sum_scale[c*VectorClass::size()+i] = scale_dad[c*VectorClass::size()+i] + scale_node[c*VectorClass::size()+i];
//                    if (min_scale[i] > sum_scale[c*VectorClass::size()+i])
//                        min_scale[i] = sum_scale[c*VectorClass::size()+i];
//                }
//            }
            for (c = 0; c < ncat_mix; c++) {
                VectorClass value = val_tmp[0] * partial_lh_node[0] * partial_lh_dad[0];
                for (i = 1; i < nstates; i++) {
                    value =  mul_add(val_tmp[i] * partial_lh_node[i], partial_lh_dad[i], value);
                }
                *lh_cat = value;
                lh_ptn += value;
                /*
                // TODO: scaling not working yet
                for (i = 0; i < VectorClass::size(); i++)
                    if (sum_scale[c*VectorClass::size()+i] == min_scale[i]) {
                        lh_cat[i] = value[i];
                        lh_ptn += value;
                    } else if (sum_scale[c*VectorClass::size()+i] == min_scale[i]+1) {
                        // only compute for least scale category
                        value *= VectorClass(SCALING_THRESHOLD);
                        *lh_cat = value;
                        lh_ptn += value;
                    }
                */
                partial_lh_node += nstates;
                partial_lh_dad += nstates;
                val_tmp += nstates;
                lh_cat ++;
            }

//			assert(lh_ptn > 0.0);
            if (ptn < orig_nptn) {
//				lh_ptn = log(abs(lh_ptn)) + LOG_SCALING_THRESHOLD*vc_min_scale;
                lh_ptn = log(abs(lh_ptn));
				lh_ptn.store_a(&_pattern_lh[ptn]);
				vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
			} else {
                // bugfix 2016-01-21, prob_const can be rescaled
                // TODO
//                if (min_scale >= 1)
//                    lh_ptn *= SCALING_THRESHOLD;
//				_pattern_lh[ptn] = lh_ptn;
				vc_prob_const += lh_ptn;
			}
		}
    }

    tree_lh += horizontal_add(vc_tree_lh);
    assert(!isnan(tree_lh) && !isinf(tree_lh) && "Numerical underflow for lh-branch");

    if (orig_nptn < nptn) {
    	// ascertainment bias correction
        prob_const = horizontal_add(vc_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        assert(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		assert(!isnan(tree_lh) && !isinf(tree_lh));
    }

    aligned_free(val);
    return tree_lh;
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const bool SAFE_NUMERIC, const int nstates>
double PhyloTree::computeLikelihoodFromBufferNewSIMD()
#else
template <class VectorClass, const bool SAFE_NUMERIC>
double PhyloTree::computeLikelihoodFromBufferSIMD()
#endif
{

	assert(theta_all && theta_computed);

	double tree_lh = current_it->lh_scale_factor + current_it_back->lh_scale_factor;
    VectorClass vc_tree_lh = 0.0;

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();

    size_t block = ncat_mix * nstates;
//    size_t tip_block = nstates * model->getNMixtures();
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();

    size_t mix_addr_nstates[ncat_mix], mix_addr[ncat_mix];
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;
    for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        mix_addr_nstates[c] = m*nstates;
        mix_addr[c] = mix_addr_nstates[c]*nstates;
    }

    double *eval = model->getEigenvalues();
    assert(eval);

    double *val0 = aligned_alloc<double>(block);
	for (c = 0; c < ncat_mix; c++) {
        size_t m = c/denom;
        double *eval_ptr = eval + mix_addr_nstates[c];
        size_t mycat = c%ncat;
        double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        size_t addr = c*nstates;
		for (i = 0; i < nstates; i++) {
			double cof = eval_ptr[i]*site_rate->getRate(mycat);
			double val = exp(cof*current_it->length) * prop;
			val0[addr+i] = val;
		}
	}


    VectorClass vc_prob_const = 0.0;
    double prob_const = 0.0;
//    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i) schedule(static)
#endif
    for (ptn = 0; ptn < nptn; ptn+=VectorClass::size()) {
		VectorClass lh_ptn;
        lh_ptn.load_a(&ptn_invar[ptn]);
		VectorClass *theta = (VectorClass*)(theta_all + ptn*block);
		for (i = 0; i < block; i++) {
			lh_ptn  = mul_add(val0[i], theta[i], lh_ptn);
		}

//        assert(lh_ptn > 0.0);
        if (ptn < orig_nptn) {
            //lh_ptn = log(abs(lh_ptn)) + LOG_SCALING_THRESHOLD*min_scale;
            lh_ptn = log(abs(lh_ptn));
            lh_ptn.store_a(&_pattern_lh[ptn]);
            vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
        } else {
            // bugfix 2016-01-21, prob_const can be rescaled
//                if (min_scale >= 1)
//                    lh_ptn *= SCALING_THRESHOLD;
//				_pattern_lh[ptn] = lh_ptn;
            vc_prob_const += lh_ptn;
        }
    }

    tree_lh += horizontal_add(vc_tree_lh);
    assert(!isnan(tree_lh) && !isinf(tree_lh) && "Numerical underflow for lh-FromBuffer");

    if (orig_nptn < nptn) {
    	// ascertainment bias correction
        prob_const = horizontal_add(vc_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        assert(prob_const < 1.0 && prob_const >= 0.0);

        // BQM 2015-10-11: fix this those functions using _pattern_lh_cat
//        double inv_const = 1.0 / (1.0-prob_const);
//        size_t nptn_cat = orig_nptn*ncat;
//    	for (ptn = 0; ptn < nptn_cat; ptn++)
//            _pattern_lh_cat[ptn] *= inv_const;
        
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
		assert(!isnan(tree_lh) && !isinf(tree_lh));
    }

    aligned_free(val0);
    return tree_lh;
}


#endif //PHYLOKERNELNEW_H_
