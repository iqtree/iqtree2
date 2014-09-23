/*
 * phylotreeeigen.cpp
 *
 *  Created on: Sep 15, 2014
 *      Author: minh
 */




#include "phylotree.h"
#include "gtrmodel.h"

/**
 * this version uses RAxML technique that stores the product of partial likelihoods and eigenvectors at node
 * for faster branch length optimization
 */
void PhyloTree::computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
	assert(dad);
    if (dad_branch->partial_lh_computed & 1)
        return;
    dad_branch->partial_lh_computed |= 1;
    size_t ptn, c;
    size_t nptn = aln->size(), ncat = site_rate->getNRate();
    size_t nstates = aln->num_states, nstatesqr = nstates*nstates, i, x;
    size_t block = nstates * ncat;

    PhyloNode *node = (PhyloNode*)(dad_branch->node);
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
			} else if (state == STATE_UNKNOWN) {
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
		delete [] inv_evec;
		delete [] evec;
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
		computePartialLikelihoodEigen(left, node, pattern_scale);
	if ((right->partial_lh_computed & 1) == 0)
		computePartialLikelihoodEigen(right, node, pattern_scale);
	dad_branch->lh_scale_factor = left->lh_scale_factor + right->lh_scale_factor;
	double *partial_lh_left = left->partial_lh, *partial_lh_right = right->partial_lh;
	double *partial_lh_tmp = new double[nstates];
	double *eleft = new double[block*nstates], *eright = new double[block*nstates];
	for (c = 0; c < ncat; c++) {
		double len_left = site_rate->getRate(c) * left->length;
		double len_right = site_rate->getRate(c) * right->length;
		for (x = 0; x < nstates; x++)
			for (i = 0; i < nstates; i++) {
				eleft[c*nstatesqr+x*nstates+i] = evec[x*nstates+i] * exp(eval[i]*len_left);
				eright[c*nstatesqr+x*nstates+i] = evec[x*nstates+i] * exp(eval[i]*len_right);
			}
	}

	//memset(partial_lh, 0, nptn*block*sizeof(double));
	if (left->node->isLeaf() && right->node->isLeaf()) {
		// special treatment for TIP-TIP (cherry) case
		// scale number must be ZERO
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		for (ptn = 0; ptn < nptn; ptn++) {
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
					for (i = 0; i < nstates; i++) {
						vleft += eleft[c*nstatesqr+x*nstates+i] * partial_lh_left[i];
						vright += eright[c*nstatesqr+x*nstates+i] * partial_lh_right[i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					partial_lh[c*nstates+i] = res;
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
		for (ptn = 0; ptn < nptn; ptn++) {
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
					for (i = 0; i < nstates; i++) {
						vleft += eleft[c*nstatesqr+x*nstates+i] * partial_lh_left[i];
						vright += eright[c*nstatesqr+x*nstates+i] * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					partial_lh[c*nstates+i] = res;
				}
			}
			partial_lh += block;
			partial_lh_left += nstates;
			partial_lh_right += block;
		}
	} else {
		// both left and right are internal node
		for (ptn = 0; ptn < nptn; ptn++) {
			dad_branch->scale_num[ptn] = left->scale_num[ptn] + right->scale_num[ptn];
			for (c = 0; c < ncat; c++) {
				for (x = 0; x < nstates; x++) {
					double vleft = 0.0, vright = 0.0;
					for (i = 0; i < nstates; i++) {
						vleft += eleft[c*nstatesqr+x*nstates+i] * partial_lh_left[c*nstates+i];
						vright += eright[c*nstatesqr+x*nstates+i] * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i*nstates+x];
					partial_lh[c*nstates+i] = res;
				}
			}
			partial_lh += block;
			partial_lh_left += block;
			partial_lh_right += block;
		}
	}

	delete [] eright;
	delete [] eleft;
	delete [] partial_lh_tmp;
	delete [] inv_evec;
	delete [] evec;
}


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
        computePartialLikelihoodEigen(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    size_t ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    size_t nstates = aln->num_states;
    size_t block = ncat * nstates;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i;
    size_t nptn = aln->size();
    double *eval = model->getEigenvalues();
    assert(eval);

    double *partial_lh_dad = dad_branch->partial_lh;
    double *partial_lh_node = node_branch->partial_lh;
    double *val = new double[block], *cof = new double[block];
	for (c = 0; c < ncat; c++) {
		double len = site_rate->getRate(c)*dad_branch->length;
		for (i = 0; i < nstates; i++) {
			val[c*nstates+i] = exp(eval[i]*len);
			cof[c*nstates+i] = eval[i]*site_rate->getRate(c);
		}
	}
    if (dad->isLeaf()) {
    	// special treatment for TIP-INTERNAL NODE case
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0, df_ptn = 0.0, ddf_ptn = 0.0;
			for (c = 0; c < ncat; c++)
				for (i = 0; i < nstates; i++) {
					double res = val[c*nstates+i] * partial_lh_node[i] * partial_lh_dad[c*nstates+i];
					lh_ptn += res;
					res *= cof[c*nstates+i];
					df_ptn += res;
					res *= cof[c*nstates+i];
					ddf_ptn += res;
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
	        lh_ptn = log(lh_ptn);
	        tree_lh += lh_ptn * freq;
	        _pattern_lh[ptn] = lh_ptn;
			partial_lh_node += nstates;
			partial_lh_dad += block;
		}
    } else {
    	// both dad and node are internal nodes
		for (ptn = 0; ptn < nptn; ptn++) {
			double lh_ptn = 0.0, df_ptn = 0.0, ddf_ptn = 0.0;
			for (i = 0; i < block; i++) {
				double res = val[i] * partial_lh_node[i] * partial_lh_dad[i];
				lh_ptn += res;
				res *= cof[i];
				df_ptn += res;
				res *= cof[i];
				ddf_ptn += res;
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
	        lh_ptn = log(lh_ptn);
	        tree_lh += lh_ptn * freq;
	        _pattern_lh[ptn] = lh_ptn;
			partial_lh_node += block;
			partial_lh_dad += block;
		}
    }
    delete [] cof;
    delete [] val;
    return tree_lh;
}

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
        computePartialLikelihoodEigen(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodEigen(node_branch, node);
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    size_t ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    assert(p_invar == 0.0); // +I model not supported yet
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    size_t nstates = aln->num_states;
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
			for (c = 0; c < ncat; c++)
				for (i = 0; i < nstates; i++)
					lh_ptn +=  val[c*nstates+i] * partial_lh_node[i] * partial_lh_dad[c*nstates+i];
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
			double lh_ptn = 0.0;
			for (i = 0; i < block; i++)
				lh_ptn +=  val[i] * partial_lh_node[i] * partial_lh_dad[i];
			lh_ptn *= p_var_cat;
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
