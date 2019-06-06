/*
 * phylokernelnonrev.h
 * Kernel based on vectorizing over alignment patterns for non-reversible models
 *
 *  Created on: Nov 4, 2016
 *      Author: minh
 */


#if !defined(PHYLOKERNELNONREV_H_) || !defined(PHYLOKERNELNONREV_STATE_H_)

#ifdef KERNEL_FIX_STATES
    #define PHYLOKERNELNONREV_STATE_H_
#else
    #define PHYLOKERNELNONREV_H_
#endif

#include "phylotree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#include <thread>

using namespace std;



#ifdef KERNEL_FIX_STATES
template <class VectorClass, const int nstates, const bool FMA>
void PhyloTree::computeNonrevPartialLikelihoodSIMD(TraversalInfo &info, size_t ptn_lower, size_t ptn_upper, int thread_id) {
#else
template <class VectorClass, const bool FMA>
void PhyloTree::computeNonrevPartialLikelihoodGenericSIMD(TraversalInfo &info, size_t ptn_lower, size_t ptn_upper, int thread_id) {
#endif

    PhyloNeighbor *dad_branch = info.dad_branch;
    PhyloNode *dad = info.dad;

	ASSERT(dad);
    PhyloNode *node = (PhyloNode*)(dad_branch->node);

//    assert(dad_branch->direction != UNDEFINED_DIRECTION);

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif

	if (node->isLeaf()) {
		return;
	}
    
    ASSERT(node->degree() >= 3);
    
    size_t ptn, c;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t i, x;
    size_t block = nstates * ncat_mix;

	// internal node
	PhyloNeighbor *left = NULL, *right = NULL; // left & right are two neighbors leading to 2 subtrees
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!left) left = (PhyloNeighbor*)(*it); else right = (PhyloNeighbor*)(*it);
	}

    // precomputed buffer to save times
    double *buffer_partial_lh_ptr = buffer_partial_lh + (getBufferPartialLhSize() - 2*block*VectorClass::size()*num_threads);
    double *echildren = info.echildren;
    double *partial_lh_leaves = info.partial_lh_leaves;

    double *eleft = echildren, *eright = echildren + block*nstates;
    
	if ((!left->node->isLeaf() && right->node->isLeaf())) {
		PhyloNeighbor *tmp = left;
		left = right;
		right = tmp;
        double *etmp = eleft;
        eleft = eright;
        eright = etmp;
    }

    if (node->degree() > 3) {

        /*--------------------- multifurcating node ------------------*/
        double *vec_tip = buffer_partial_lh_ptr + (block*2)*VectorClass::size()*thread_id;
        VectorClass *vtip = (VectorClass*)vec_tip;

        // now for-loop computing partial_lh over all site-patterns
        for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
            VectorClass *partial_lh_all = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            for (i = 0; i < block; i++)
                partial_lh_all[i] = 1.0;
            memset(&dad_branch->scale_num[ptn], 0, sizeof(UBYTE)*VectorClass::size());
                
            double *partial_lh_leaf = partial_lh_leaves;
            double *echild = echildren;

            FOR_NEIGHBOR_IT(node, dad, it) {
                PhyloNeighbor *child = (PhyloNeighbor*)*it;
                if (child->node->isLeaf()) {
                    // external node
                    // load data for tip
                    for (x = 0; x < VectorClass::size(); x++) {
                        double *tip_child;
                        if (isRootLeaf(child->node))
                            tip_child = partial_lh_leaf;
                        else if (ptn+x < orig_nptn)
                            tip_child = partial_lh_leaf + block * (aln->at(ptn+x))[child->node->id];
                        else if (ptn+x < max_orig_nptn)
                            tip_child = partial_lh_leaf + block * aln->STATE_UNKNOWN;
                        else if (ptn+x < nptn)
                            tip_child = partial_lh_leaf + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                        else
                            tip_child = partial_lh_leaf + block * aln->STATE_UNKNOWN;
                        double *this_vec_tip = vec_tip+x;
                        for (i = 0; i < block; i++) {
                            *this_vec_tip = tip_child[i];
                            this_vec_tip += VectorClass::size();
                        }
                    }
                    for (c = 0; c < block; c++) {
                        // compute real partial likelihood vector
                        partial_lh_all[c] *= vtip[c];
                    }
                    partial_lh_leaf += (aln->STATE_UNKNOWN+1)*block;
                } else {
                    // internal node
                    VectorClass *partial_lh = partial_lh_all;
                    VectorClass *partial_lh_child = (VectorClass*)(child->partial_lh + ptn*block);
                    for (i = 0; i < VectorClass::size(); i++)
                        dad_branch->scale_num[ptn+i] += child->scale_num[ptn+i];

                    double *echild_ptr = echild;
                    for (c = 0; c < ncat_mix; c++) {
                        // compute real partial likelihood vector
                        for (x = 0; x < nstates; x++) {
                            VectorClass vchild;
//                            for (i = 0; i < nstates; i++) {
//                                vchild += echild_ptr[i] * partial_lh_child[i];
//                            }
#ifdef KERNEL_FIX_STATES
                            dotProductVec<VectorClass, double, nstates, FMA>(echild_ptr, partial_lh_child, vchild);
#else
                            dotProductVec<VectorClass, double, FMA>(echild_ptr, partial_lh_child, vchild, nstates);
#endif
                            echild_ptr += nstates;
                            partial_lh[x] *= vchild;
                        }
                        partial_lh += nstates;
                        partial_lh_child += nstates;
                    }
                } // if
                echild += block*nstates;
            } // FOR_NEIGHBOR
            
        
            VectorClass lh_max = partial_lh_all[0];
            for (i = 1; i < block; i++)
                lh_max = max(lh_max, partial_lh_all[i]);

            // check if one should scale partial likelihoods
            auto underflown = (lh_max < SCALING_THRESHOLD);
            if (horizontal_or(underflown)) {
                // now do the likelihood scaling
                for (x = 0; x < VectorClass::size(); x++)
                if (underflown[x]) {
                    double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                    // now do the likelihood scaling
                    for (i = 0; i < block; i++) {
                        partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                    }
                    dad_branch->scale_num[ptn+x] += 1;
                }
            }

        } // for ptn

        // end multifurcating treatment
    } else if (left->node->isLeaf() && right->node->isLeaf()) {

        /*--------------------- TIP-TIP (cherry) case ------------------*/

        double *partial_lh_left = partial_lh_leaves;
        double *partial_lh_right = partial_lh_leaves + (aln->STATE_UNKNOWN+1)*block;
        double *vec_left = buffer_partial_lh_ptr + (block*2)*VectorClass::size()*thread_id;
        double *vec_right =  &vec_left[block*VectorClass::size()];

        if (isRootLeaf(right->node)) {
            // swap so that left node is the root
            PhyloNeighbor *tmp = left;
            left = right;
            right = tmp;
            double *etmp = eleft;
            eleft = eright;
            eright = etmp;
            etmp = partial_lh_left;
            partial_lh_left = partial_lh_right;
            partial_lh_right = etmp;
        }
    
		// scale number must be ZERO
	    memset(dad_branch->scale_num + ptn_lower, 0, (ptn_upper-ptn_lower) * sizeof(UBYTE));

        if (isRootLeaf(left->node)) {
            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                double *vright = dad_branch->partial_lh + ptn*block;
                VectorClass *partial_lh = (VectorClass*)vright;
                // load data for tip
                for (x = 0; x < VectorClass::size(); x++) {
                    double *tip_right;
                    if (ptn+x < orig_nptn)
                        tip_right = partial_lh_right + block * (aln->at(ptn+x))[right->node->id];
                    else if (ptn+x < max_orig_nptn)
                        tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                    else if (ptn+x < nptn)
                        tip_right = partial_lh_right + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                    else
                        tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                    double *this_vec_right = vright+x;
                    for (i = 0; i < block; i++) {
                        *this_vec_right = tip_right[i];
                        this_vec_right += VectorClass::size();
                    }
                }
                for (i = 0; i < block; i++)
                    partial_lh[i] *= partial_lh_left[i];
            }
        } else
		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            VectorClass *vright = (VectorClass*)vec_right;
            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip_left, *tip_right;
                if (ptn+x < orig_nptn) {
                    tip_left  = partial_lh_left  + block * (aln->at(ptn+x))[left->node->id];
                    tip_right = partial_lh_right + block * (aln->at(ptn+x))[right->node->id];
                } else if (ptn+x < max_orig_nptn) {
                    tip_left  = partial_lh_left  + block * aln->STATE_UNKNOWN;
                    tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                } else if (ptn+x < nptn) {
                    tip_left  = partial_lh_left  + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                    tip_right = partial_lh_right + block * model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                } else {
                    tip_left  = partial_lh_left  + block * aln->STATE_UNKNOWN;
                    tip_right = partial_lh_right + block * aln->STATE_UNKNOWN;
                }
                double *this_vec_left = vec_left+x;
                double *this_vec_right = vec_right+x;
                for (i = 0; i < block; i++) {
                    *this_vec_left = tip_left[i];
                    *this_vec_right = tip_right[i];
                    this_vec_left += VectorClass::size();
                    this_vec_right += VectorClass::size();
                }
            }


            for (i = 0; i < block; i++)
                partial_lh[i] = vleft[i] * vright[i];
		}
	} else if (isRootLeaf(left->node) && !right->node->isLeaf()) {
        // left is root node
        /*--------------------- ROOT-INTERNAL NODE case ------------------*/

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num + ptn_lower, right->scale_num + ptn_lower, (ptn_upper-ptn_lower) * sizeof(UBYTE));

        double *partial_lh_left = partial_lh_leaves;

		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            double *eright_ptr = eright;
            double *lh_left = partial_lh_left;

			for (c = 0; c < ncat_mix; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					VectorClass vright;
#ifdef KERNEL_FIX_STATES
                    dotProductVec<VectorClass, double, nstates, FMA>(eright_ptr, partial_lh_right, vright);
#else
                    dotProductVec<VectorClass, double, FMA>(eright_ptr, partial_lh_right, vright, nstates);
#endif
                    eright_ptr += nstates;
                    partial_lh[x] = lh_left[x]*vright;
				}
                partial_lh_right += nstates;
                lh_left += nstates;
                partial_lh += nstates;
			}
		}

	} else if (left->node->isLeaf() && !right->node->isLeaf()) {

        /*--------------------- TIP-INTERNAL NODE case ------------------*/

		// only take scale_num from the right subtree
		memcpy(dad_branch->scale_num + ptn_lower, right->scale_num + ptn_lower, (ptn_upper-ptn_lower) * sizeof(UBYTE));


        double *partial_lh_left = partial_lh_leaves;
        double *vec_left = buffer_partial_lh_ptr + (block*2)*VectorClass::size()*thread_id;

		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass *vleft = (VectorClass*)vec_left;
            // load data for tip
            for (x = 0; x < VectorClass::size(); x++) {
                double *tip;
                if (ptn+x < orig_nptn)
                    tip = partial_lh_left + block*(aln->at(ptn+x))[left->node->id];
                else if (ptn+x < max_orig_nptn)
                    tip = partial_lh_left + block*aln->STATE_UNKNOWN;
                else if (ptn+x < nptn)
                    tip = partial_lh_left + block*model_factory->unobserved_ptns[ptn+x-max_orig_nptn];
                else
                    tip = partial_lh_left + block*aln->STATE_UNKNOWN;
                double *this_vec_left = vec_left+x;
                for (i = 0; i < block; i++) {
                    *this_vec_left = tip[i];
                    this_vec_left += VectorClass::size();
                }
            }

            VectorClass lh_max = 0.0;
            
            double *eright_ptr = eright;

			for (c = 0; c < ncat_mix; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
					VectorClass vright;
#ifdef KERNEL_FIX_STATES
                    dotProductVec<VectorClass, double, nstates, FMA>(eright_ptr, partial_lh_right, vright);
#else
                    dotProductVec<VectorClass, double, FMA>(eright_ptr, partial_lh_right, vright, nstates);
#endif
                    eright_ptr += nstates;
                    lh_max = max(lh_max, (partial_lh[x] = vleft[x]*vright));
				}
                vleft += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
			}
            // check if one should scale partial likelihoods
            auto underflown = (lh_max < SCALING_THRESHOLD);
            if (horizontal_or(underflown)) {
                // now do the likelihood scaling
                for (x = 0; x < VectorClass::size(); x++)
                if (underflown[x]) {
                    double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                    // now do the likelihood scaling
                    for (i = 0; i < block; i++) {
                        partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                    }
                    dad_branch->scale_num[ptn+x] += 1;
                }
            }
		}

	} else {

        /*--------------------- INTERNAL-INTERNAL NODE case ------------------*/

		for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
			VectorClass *partial_lh = (VectorClass*)(dad_branch->partial_lh + ptn*block);
			VectorClass *partial_lh_left = (VectorClass*)(left->partial_lh + ptn*block);
			VectorClass *partial_lh_right = (VectorClass*)(right->partial_lh + ptn*block);
            VectorClass lh_max = 0.0;
            for (i = 0; i < VectorClass::size(); i++)
                dad_branch->scale_num[ptn+i] = left->scale_num[ptn+i] + right->scale_num[ptn+i];

            double *eleft_ptr = eleft;
            double *eright_ptr = eright;

			for (c = 0; c < ncat_mix; c++) {
				// compute real partial likelihood vector
				for (x = 0; x < nstates; x++) {
#ifdef KERNEL_FIX_STATES
                    dotProductDualVec<VectorClass, double, nstates, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh[x]);
#else
                    dotProductDualVec<VectorClass, double, FMA>(eleft_ptr, partial_lh_left, eright_ptr, partial_lh_right, partial_lh[x], nstates);
#endif
                    eleft_ptr += nstates;
                    eright_ptr += nstates;
					lh_max=max(lh_max, partial_lh[x]);
				}
                partial_lh_left += nstates;
                partial_lh_right += nstates;
                partial_lh += nstates;
            }

            // check if one should scale partial likelihoods
            auto underflown = (lh_max < SCALING_THRESHOLD);
            if (horizontal_or(underflown)) {
                // now do the likelihood scaling
                for (x = 0; x < VectorClass::size(); x++)
                if (underflown[x]) {
                    double *partial_lh = dad_branch->partial_lh + (ptn*block + x);
                    // now do the likelihood scaling
                    for (i = 0; i < block; i++) {
                        partial_lh[i*VectorClass::size()] = ldexp(partial_lh[i*VectorClass::size()], SCALING_THRESHOLD_EXP);
                    }
                    dad_branch->scale_num[ptn+x] += 1;
                }
            }

		}

	}
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const int nstates, const bool FMA>
void PhyloTree::computeNonrevLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double *df, double *ddf) {
#else
template <class VectorClass, const bool FMA>
void PhyloTree::computeNonrevLikelihoodDervGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double *df, double *ddf) {
#endif

//    assert(rooted);

    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf() || (dad_branch->direction == AWAYFROM_ROOT && !isRootLeaf(dad))) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, false);
#endif

#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t nstatesqr = nstates*nstates;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    size_t block = ncat_mix * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    bool isASC = model_factory->unobserved_ptns.size() > 0;

//    double *trans_mat = new double[block*nstates*3];
    double *trans_mat = buffer_partial_lh;
    double *trans_derv1 = buffer_partial_lh + block*nstates;
    double *trans_derv2 = trans_derv1 + block*nstates;
    double *buffer_partial_lh_ptr = buffer_partial_lh + get_safe_upper_limit(3*block*nstates);

	for (c = 0; c < ncat_mix; c++) {
        size_t mycat = c%ncat;
        size_t m = c/denom;
        double cat_rate = site_rate->getRate(mycat);
		double len = cat_rate * dad_branch->length;
		double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        double *this_trans_mat = &trans_mat[c*nstatesqr];
        double *this_trans_derv1 = &trans_derv1[c*nstatesqr];
        double *this_trans_derv2 = &trans_derv2[c*nstatesqr];
        model->computeTransDerv(len, this_trans_mat, this_trans_derv1, this_trans_derv2, m);
        double prop_rate = prop * cat_rate;
        double prop_rate_2 = prop_rate * cat_rate;
		for (i = 0; i < nstatesqr; i++) {
			this_trans_mat[i] *= prop;
            this_trans_derv1[i] *= prop_rate;
            this_trans_derv2[i] *= prop_rate_2;
        }
        if (!rooted) {
            // for unrooted tree, multiply with state_freq
            double state_freq[nstates];
            model->getStateFrequency(state_freq, m);
            for (i = 0; i < nstates; i++) {
                for (size_t x = 0; x < nstates; x++) {
                    this_trans_mat[x] *= state_freq[i];
                    this_trans_derv1[x] *= state_freq[i];
                    this_trans_derv2[x] *= state_freq[i];
                }
                this_trans_mat += nstates;
                this_trans_derv1 += nstates;
                this_trans_derv2 += nstates;
            }
        }
	}

    VectorClass all_df(0.0), all_ddf(0.0);
    VectorClass all_prob_const(0.0), all_df_const(0.0), all_ddf_const(0.0);

    vector<size_t> limits;
    computeBounds<VectorClass>(num_threads, nptn, limits);
//    double *buffer_partial_lh_ptr = buffer_partial_lh;

    if (dad->isLeaf()) {
         // make sure that we do not estimate the virtual branch length from the root
        ASSERT(!isRootLeaf(dad));
    	// special treatment for TIP-INTERNAL NODE case
//    	double *partial_lh_node = new double[(aln->STATE_UNKNOWN+1)*block*3];
        double *partial_lh_node = buffer_partial_lh_ptr;
        double *partial_lh_derv1 = partial_lh_node + (aln->STATE_UNKNOWN+1)*block;
        double *partial_lh_derv2 = partial_lh_derv1 + (aln->STATE_UNKNOWN+1)*block;
        buffer_partial_lh_ptr += get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block*3);
        IntVector states_dad = aln->seq_states[dad->id];
        states_dad.push_back(aln->STATE_UNKNOWN);
        // precompute information from one tip
        for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
            double *lh_node = partial_lh_node +(*it)*block;
            double *lh_derv1 = partial_lh_derv1 +(*it)*block;
            double *lh_derv2 = partial_lh_derv2 +(*it)*block;
            double *lh_tip = tip_partial_lh + (*it)*nstates;
            double *trans_mat_tmp = trans_mat;
            double *trans_derv1_tmp = trans_derv1;
            double *trans_derv2_tmp = trans_derv2;
            for (c = 0; c < ncat_mix; c++) {
                for (i = 0; i < nstates; i++) {
                    lh_node[i] = 0.0;
                    lh_derv1[i] = 0.0;
                    lh_derv2[i] = 0.0;
                    for (size_t x = 0; x < nstates; x++) {
                        lh_node[i] += trans_mat_tmp[x] * lh_tip[x];
                        lh_derv1[i] += trans_derv1_tmp[x] * lh_tip[x];
                        lh_derv2[i] += trans_derv2_tmp[x] * lh_tip[x];
                    }
                    trans_mat_tmp += nstates;
                    trans_derv1_tmp += nstates;
                    trans_derv2_tmp += nstates;
                }
                lh_node += nstates;
                lh_derv1 += nstates;
                lh_derv2 += nstates;
            }
        }

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static,1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0), vc_df_const(0.0), vc_ddf_const(0.0);
            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];
            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            double *vec_tip = buffer_partial_lh_ptr + block*3*VectorClass::size()*thread_id;

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn, df_ptn, ddf_ptn;
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);

                //load tip vector
                for (i = 0; i < VectorClass::size(); i++) {
                    size_t state_dad;
                    if (ptn+i < orig_nptn)
                        state_dad = block * (aln->at(ptn+i))[dad->id];
                    else if (ptn+i < max_orig_nptn)
                        state_dad = block * aln->STATE_UNKNOWN;
                    else if (ptn+i < nptn)
                        state_dad = block * model_factory->unobserved_ptns[ptn+i-max_orig_nptn];
                    else
                        state_dad = block * aln->STATE_UNKNOWN;

                    double *lh_tip = partial_lh_node + state_dad;
                    double *lh_derv1 = partial_lh_derv1 + state_dad;
                    double *lh_derv2 = partial_lh_derv2 + state_dad;

                    double *this_vec_tip = vec_tip+i;
                    double *this_derv1 = this_vec_tip + block*VectorClass::size();
                    double *this_derv2 = this_derv1 + block*VectorClass::size();
                    for (c = 0; c < block; c++) {
                        *this_vec_tip = lh_tip[c];
                        *this_derv1 = lh_derv1[c];
                        *this_derv2 = lh_derv2[c];
                        this_vec_tip += VectorClass::size();
                        this_derv1 += VectorClass::size();
                        this_derv2 += VectorClass::size();
                    }

                }

                VectorClass *lh_node = (VectorClass*)vec_tip;
                VectorClass *lh_derv1 = (VectorClass*)vec_tip + block;
                VectorClass *lh_derv2 = (VectorClass*)lh_derv1 + block;
#ifdef KERNEL_FIX_STATES
                dotProductTriple<VectorClass, VectorClass, nstates, FMA, false>(lh_node, lh_derv1, lh_derv2, partial_lh_dad, lh_ptn, df_ptn, ddf_ptn, block);
#else
                dotProductTriple<VectorClass, VectorClass, FMA, false>(lh_node, lh_derv1, lh_derv2, partial_lh_dad, lh_ptn, df_ptn, ddf_ptn, block, nstates);
#endif
                lh_ptn = (lh_ptn + VectorClass().load_a(&ptn_invar[ptn]));

                if (ptn < orig_nptn) {
                    lh_ptn = 1.0 / lh_ptn;
                    VectorClass df_frac = df_ptn * lh_ptn;
                    VectorClass ddf_frac = ddf_ptn * lh_ptn;
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass tmp1 = df_frac * freq;
                    VectorClass tmp2 = ddf_frac * freq;
                    my_df += tmp1;
                    my_ddf += nmul_add(tmp1, df_frac, tmp2);
                } else {
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                        df_ptn.cutoff(nptn-ptn);
                        ddf_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    double *lh_ptn_ptr = (double*)&lh_ptn;
                    for (i = 0; i < VectorClass::size(); i++)
                        if (dad_branch->scale_num[ptn+i] >= 1)
                            lh_ptn_ptr[i] *= SCALING_THRESHOLD;
                    vc_prob_const += lh_ptn;
                    vc_df_const += df_ptn;
                    vc_ddf_const += ddf_ptn;
                }
            } // FOR ptn
        #ifdef _OPENMP
        #pragma omp critical
        #endif
            {
                all_df += my_df;
                all_ddf += my_ddf;
                if (isASC) {
                    all_prob_const += vc_prob_const;
                    all_df_const += vc_df_const;
                    all_ddf_const += vc_ddf_const;
                }
            }
        } // FOR thread_id

//		delete [] partial_lh_node;
    } else {

    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static,1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            VectorClass my_df(0.0), my_ddf(0.0), vc_prob_const(0.0), vc_df_const(0.0), vc_ddf_const(0.0);
            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];
            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn(0.0), df_ptn(0.0), ddf_ptn(0.0);
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
                double *trans_mat_tmp = trans_mat;
                double *trans_derv1_tmp = trans_derv1;
                double *trans_derv2_tmp = trans_derv2;
                for (c = 0; c < ncat_mix; c++) {
                    for (i = 0; i < nstates; i++) {
                        VectorClass lh_state;
                        VectorClass lh_derv1;
                        VectorClass lh_derv2;
#ifdef KERNEL_FIX_STATES
                        dotProductTriple<VectorClass, double, nstates, FMA, false>(trans_mat_tmp, trans_derv1_tmp, trans_derv2_tmp, partial_lh_node, lh_state, lh_derv1, lh_derv2, nstates);
#else
                        dotProductTriple<VectorClass, double, FMA, false>(trans_mat_tmp, trans_derv1_tmp, trans_derv2_tmp, partial_lh_node, lh_state, lh_derv1, lh_derv2, nstates, nstates);
#endif
                        lh_ptn = mul_add(partial_lh_dad[i], lh_state, lh_ptn);
                        df_ptn = mul_add(partial_lh_dad[i], lh_derv1, df_ptn);
                        ddf_ptn = mul_add(partial_lh_dad[i], lh_derv2, ddf_ptn);
                        trans_mat_tmp += nstates;
                        trans_derv1_tmp += nstates;
                        trans_derv2_tmp += nstates;
                    }
                    partial_lh_node += nstates;
                    partial_lh_dad += nstates;
                }
                lh_ptn = (lh_ptn + VectorClass().load_a(&ptn_invar[ptn]));

                if (ptn < orig_nptn) {
                    lh_ptn = 1.0 / lh_ptn;
                    VectorClass df_frac = df_ptn * lh_ptn;
                    VectorClass ddf_frac = ddf_ptn * lh_ptn;
                    VectorClass freq;
                    freq.load_a(&ptn_freq[ptn]);
                    VectorClass tmp1 = df_frac * freq;
                    VectorClass tmp2 = ddf_frac * freq;
                    my_df += tmp1;
                    my_ddf += nmul_add(tmp1, df_frac, tmp2);
                } else {
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                        df_ptn.cutoff(nptn-ptn);
                        ddf_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    double *lh_ptn_ptr = (double*)&lh_ptn;
                    for (i = 0; i < VectorClass::size(); i++)
                        if (dad_branch->scale_num[ptn+i] >= 1)
                            lh_ptn_ptr[i] *= SCALING_THRESHOLD;
                    vc_prob_const += lh_ptn;
                    vc_df_const += df_ptn;
                    vc_ddf_const += ddf_ptn;
                }
            } // FOR ptn
        #ifdef _OPENMP
        #pragma omp critical
        #endif
            {
                all_df += my_df;
                all_ddf += my_ddf;
                if (isASC) {
                    all_prob_const += vc_prob_const;
                    all_df_const += vc_df_const;
                    all_ddf_const += vc_ddf_const;
                }
            }
        } // FOR thread
    }

	*df = horizontal_add(all_df);
	*ddf = horizontal_add(all_ddf);
    ASSERT(std::isfinite(*df) && "Numerical underflow for non-rev lh-derivative");

	if (isASC) {
        double prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;
        prob_const = horizontal_add(all_prob_const);
        df_const = horizontal_add(all_df_const);
        ddf_const = horizontal_add(all_ddf_const);
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	*df += nsites * df_frac;
    	*ddf += nsites *(ddf_frac + df_frac*df_frac);
    }
}

#ifdef KERNEL_FIX_STATES
template <class VectorClass, const int nstates, const bool FMA>
double PhyloTree::computeNonrevLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
#else
template <class VectorClass, const bool FMA>
double PhyloTree::computeNonrevLikelihoodBranchGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad) {
#endif

//    assert(rooted);

    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf() || (dad_branch->direction == AWAYFROM_ROOT && !isRootLeaf(dad))) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }

#ifdef KERNEL_FIX_STATES
    computeTraversalInfo<VectorClass, nstates>(node, dad, false);
#else
    computeTraversalInfo<VectorClass>(node, dad, false);
#endif

    double tree_lh = 0.0;
#ifndef KERNEL_FIX_STATES
    size_t nstates = aln->num_states;
#endif
    size_t nstatesqr = nstates*nstates;
    size_t ncat = site_rate->getNRate();
    size_t ncat_mix = (model_factory->fused_mix_rate) ? ncat : ncat*model->getNMixtures();
    size_t denom = (model_factory->fused_mix_rate) ? 1 : ncat;

    size_t block = ncat_mix * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t orig_nptn = aln->size();
    size_t max_orig_nptn = ((orig_nptn+VectorClass::size()-1)/VectorClass::size())*VectorClass::size();
    size_t nptn = max_orig_nptn+model_factory->unobserved_ptns.size();
    bool isASC = model_factory->unobserved_ptns.size() > 0;

    vector<size_t> limits;
    computeBounds<VectorClass>(num_threads, nptn, limits);

//    double *trans_mat = new double[block*nstates];
    double *trans_mat = buffer_partial_lh;
    double *buffer_partial_lh_ptr = buffer_partial_lh + block*nstates;
	for (c = 0; c < ncat_mix; c++) {
        size_t mycat = c%ncat;
        size_t m = c/denom;
		double len = site_rate->getRate(mycat) * dad_branch->length;
		double prop = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        double *this_trans_mat = &trans_mat[c*nstatesqr];
        model->computeTransMatrix(len, this_trans_mat, m);
		for (i = 0; i < nstatesqr; i++)
			this_trans_mat[i] *= prop;
        if (!rooted) {
            // if unrooted tree, multiply with frequency
            double state_freq[nstates];
            model->getStateFrequency(state_freq, m);
            for (i = 0; i < nstates; i++) {
                for (size_t x = 0; x < nstates; x++)
                    this_trans_mat[x] *= state_freq[i];
                this_trans_mat += nstates;
            }
        }
	}

    VectorClass all_tree_lh(0.0);
    VectorClass all_prob_const(0.0);

    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
//    	double *partial_lh_node = new double[(aln->STATE_UNKNOWN+1)*block];
        double *partial_lh_node = buffer_partial_lh_ptr;
        buffer_partial_lh_ptr += get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block);

        if (isRootLeaf(dad)) {
            for (c = 0; c < ncat_mix; c++) {
                double *lh_node = partial_lh_node + c*nstates;
                size_t m = c/denom;
                model->getStateFrequency(lh_node, m);
                double prop = site_rate->getProp(c%ncat) * model->getMixtureWeight(m);
                for (i = 0; i < nstates; i++)
                    lh_node[i] *= prop;
            }
        } else {
            IntVector states_dad = aln->seq_states[dad->id];
            states_dad.push_back(aln->STATE_UNKNOWN);
            // precompute information from one tip
            for (IntVector::iterator it = states_dad.begin(); it != states_dad.end(); it++) {
                double *lh_node = partial_lh_node +(*it)*block;
                double *lh_tip = tip_partial_lh + (*it)*nstates;
                double *trans_mat_tmp = trans_mat;
                for (c = 0; c < ncat_mix; c++) {
                    for (i = 0; i < nstates; i++) {
                        lh_node[i] = 0.0;
                        for (size_t x = 0; x < nstates; x++)
                            lh_node[i] += trans_mat_tmp[x] * lh_tip[x];
                        trans_mat_tmp += nstates;
                    }
                    lh_node += nstates;
                }
            }
        }

    	// now do the real computation
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static,1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            VectorClass vc_tree_lh(0.0), vc_prob_const(0.0);
            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];
            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            // reset memory for _pattern_lh_cat
//            memset(_pattern_lh_cat+ptn_lower*ncat_mix, 0, (ptn_upper-ptn_lower)*ncat_mix*sizeof(double));

            double *vec_tip = buffer_partial_lh_ptr + block*VectorClass::size()*thread_id;

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn;
                lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                VectorClass *lh_node = (VectorClass*)vec_tip;

                //load tip vector
                for (i = 0; i < VectorClass::size(); i++) {
                    double *lh_tip;
                    if (isRootLeaf(dad))
                        lh_tip = partial_lh_node;
                    else if (ptn+i < orig_nptn)
                        lh_tip = partial_lh_node + block*(aln->at(ptn+i))[dad->id];
                    else if (ptn+i < max_orig_nptn)
                        lh_tip = partial_lh_node + block*aln->STATE_UNKNOWN;
                    else if (ptn+i < nptn)
                        lh_tip = partial_lh_node + block*model_factory->unobserved_ptns[ptn+i-max_orig_nptn];
                    else
                        lh_tip = partial_lh_node + block*aln->STATE_UNKNOWN;

                    double *this_vec_tip = vec_tip+i;
                    for (c = 0; c < block; c++) {
                        *this_vec_tip = lh_tip[c];
                        this_vec_tip += VectorClass::size();
                    }

                }

                if (_pattern_lh_cat_state) {
                    // naively compute pattern_lh per category per state
                    VectorClass *lh_state = (VectorClass*)(_pattern_lh_cat_state + ptn*block);
                    for (c = 0; c < ncat_mix; c++) {
                        for (i=0; i < nstates; i++) {
                            lh_cat[c] += (lh_state[i] = lh_node[i]*partial_lh_dad[i]);
                        }
                        lh_node += nstates;
                        partial_lh_dad += nstates;
                        lh_state += nstates;
                        lh_ptn += lh_cat[c];
                    }
                } else {
                    for (c = 0; c < ncat_mix; c++) {
    #ifdef KERNEL_FIX_STATES
                        dotProductVec<VectorClass, VectorClass, nstates, FMA>(lh_node, partial_lh_dad, lh_cat[c]);
    #else
                        dotProductVec<VectorClass, VectorClass, FMA>(lh_node, partial_lh_dad, lh_cat[c], nstates);
    #endif
                        lh_node += nstates;
                        partial_lh_dad += nstates;
                        lh_ptn += lh_cat[c];
                    }
                }
                VectorClass vc_min_scale;
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                for (i = 0; i < VectorClass::size(); i++) {
                    vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i];
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;

//                lh_ptn = abs(lh_ptn);
//                assert(horizontal_and(lh_ptn > 0));
                if (ptn < orig_nptn) {
                    lh_ptn = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&_pattern_lh[ptn]);
                    vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    if (horizontal_or(vc_min_scale != 0.0)) {
                        // some entries are rescaled
                        double *lh_ptn_dbl = (double*)&lh_ptn;
                        for (i = 0; i < VectorClass::size(); i++)
                            if (vc_min_scale_ptr[i] != 0.0)
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                    }
                    vc_prob_const += lh_ptn;
                }
            } // FOR ptn
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                all_tree_lh += vc_tree_lh;
                if (isASC)
                    all_prob_const += vc_prob_const;
            }
        } // FOR thread_id
    } else {

    	// both dad and node are internal nodes
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, c) schedule(static,1) num_threads(num_threads)
#endif
        for (int thread_id = 0; thread_id < num_threads; thread_id++) {
            VectorClass vc_tree_lh(0.0), vc_prob_const(0.0);
            size_t ptn_lower = limits[thread_id];
            size_t ptn_upper = limits[thread_id+1];
            // first compute partial_lh
            for (vector<TraversalInfo>::iterator it = traversal_info.begin(); it != traversal_info.end(); it++)
                computePartialLikelihood(*it, ptn_lower, ptn_upper, thread_id);

            // reset memory for _pattern_lh_cat
            memset(_pattern_lh_cat+ptn_lower*ncat_mix, 0, (ptn_upper-ptn_lower)*ncat_mix*sizeof(double));

            for (ptn = ptn_lower; ptn < ptn_upper; ptn+=VectorClass::size()) {
                VectorClass lh_ptn;
                lh_ptn.load_a(&ptn_invar[ptn]);
                VectorClass *lh_cat = (VectorClass*)(_pattern_lh_cat + ptn*ncat_mix);
                VectorClass *partial_lh_dad = (VectorClass*)(dad_branch->partial_lh + ptn*block);
                VectorClass *partial_lh_node = (VectorClass*)(node_branch->partial_lh + ptn*block);
                double *trans_mat_tmp = trans_mat;
                if (_pattern_lh_cat_state) {
                    VectorClass *lh_state = (VectorClass*)(_pattern_lh_cat_state + ptn*block);
                    for (c = 0; c < ncat_mix; c++) {
                        for (i = 0; i < nstates; i++) {
    #ifdef KERNEL_FIX_STATES
                            dotProductVec<VectorClass, double, nstates, FMA>(trans_mat_tmp, partial_lh_node, lh_state[i]);
    #else
                            dotProductVec<VectorClass, double, FMA>(trans_mat_tmp, partial_lh_node, lh_state[i], nstates);
    #endif
                            lh_cat[c] += (lh_state[i] *= partial_lh_dad[i]);
                            trans_mat_tmp += nstates;
                        }
                        lh_ptn += lh_cat[c];
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                        lh_state += nstates;
                    }
                } else {
                    for (c = 0; c < ncat_mix; c++) {
                        for (i = 0; i < nstates; i++) {
                            VectorClass lh_state;
    #ifdef KERNEL_FIX_STATES
                            dotProductVec<VectorClass, double, nstates, FMA>(trans_mat_tmp, partial_lh_node, lh_state);
    #else
                            dotProductVec<VectorClass, double, FMA>(trans_mat_tmp, partial_lh_node, lh_state, nstates);
    #endif
                            lh_cat[c] = mul_add(partial_lh_dad[i], lh_state, lh_cat[c]);
                            trans_mat_tmp += nstates;
                        }
                        lh_ptn += lh_cat[c];
                        partial_lh_node += nstates;
                        partial_lh_dad += nstates;
                    }
                }
                VectorClass vc_min_scale;
                double* vc_min_scale_ptr = (double*)&vc_min_scale;
                for (i = 0; i < VectorClass::size(); i++) {
                    vc_min_scale_ptr[i] = dad_branch->scale_num[ptn+i] + node_branch->scale_num[ptn+i];
                }
                vc_min_scale *= LOG_SCALING_THRESHOLD;
//                lh_ptn = abs(lh_ptn);
                ASSERT(horizontal_and(lh_ptn > 0));
                if (ptn < orig_nptn) {
                    lh_ptn = log(lh_ptn) + vc_min_scale;
                    lh_ptn.store_a(&_pattern_lh[ptn]);
                    vc_tree_lh = mul_add(lh_ptn, VectorClass().load_a(&ptn_freq[ptn]), vc_tree_lh);
                } else {
                    // ascertainment bias correction
                    if (ptn+VectorClass::size() > nptn) {
                        // cutoff the last entries if going beyond
                        lh_ptn.cutoff(nptn-ptn);
                    }
                    // bugfix 2016-01-21, prob_const can be rescaled
                    if (horizontal_or(vc_min_scale != 0.0)) {
                        // some entries are rescaled
                        double *lh_ptn_dbl = (double*)&lh_ptn;
                        for (i = 0; i < VectorClass::size(); i++)
                            if (vc_min_scale_ptr[i] != 0.0)
                                lh_ptn_dbl[i] *= SCALING_THRESHOLD;
                    }
                    vc_prob_const += lh_ptn;
                }
            } // FOR ptn
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                all_tree_lh += vc_tree_lh;
                if (isASC)
                    all_prob_const += vc_prob_const;
            }
        } // FOR thread_id
    }

    tree_lh = horizontal_add(all_tree_lh);

    if (!std::isfinite(tree_lh)) {
        model->writeInfo(cout);
        site_rate->writeInfo(cout);
        ASSERT(0 && "Numerical underflow for non-rev lh-branch");
    }

    if (isASC) {
    	// ascertainment bias correction
        double prob_const = horizontal_add(all_prob_const);
        if (prob_const >= 1.0 || prob_const < 0.0) {
            printTree(cout, WT_TAXON_ID + WT_BR_LEN + WT_NEWLINE);
            model->writeInfo(cout);
        }
        ASSERT(prob_const < 1.0 && prob_const >= 0.0);

    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_nptn; ptn+=VectorClass::size())
            (VectorClass().load_a(&_pattern_lh[ptn])-prob_const).store_a(&_pattern_lh[ptn]);
    	tree_lh -= aln->getNSite()*prob_const;
		ASSERT(std::isfinite(tree_lh));
    }

    return tree_lh;
}


#endif
