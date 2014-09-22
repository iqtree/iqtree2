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
    size_t nstates = aln->num_states, i, x;
    size_t block = nstates * ncat;

    PhyloNode *node = (PhyloNode*)(dad_branch->node);
    double *partial_lh = dad_branch->partial_lh;
	assert(node->id < aln->getNSeq());

	double **evec = model->getEigenvectors(), **inv_evec = model->getInverseEigenvectors();
	double *eval = model->getEigenvalues();
	assert(inv_evec && evec);

    dad_branch->lh_scale_factor = 0.0;

	if (node->isLeaf()) {
		// external node
	    memset(dad_branch->scale_num, 0, nptn * sizeof(UBYTE));
		for (ptn = 0; ptn < nptn; ptn++) {
			int state = (aln->at(ptn))[node->id];
			if (state < nstates) {
				// simple state
				for (i = 0; i < nstates; i++)
					partial_lh[i] = inv_evec[i][state];
			} else if (state == STATE_UNKNOWN) {
				// gap or unknown state
				//dad_branch->scale_num[ptn] = -1;
				memset(partial_lh, 0, nstates*sizeof(double));
				for (i = 0; i < nstates; i++) {
					for (x = 0; x < nstates; x++) {
						partial_lh[i] += inv_evec[i][x];
					}
				}
			} else {
				// ambiguous state
				memset(partial_lh, 0, nstates*sizeof(double));
				state -= (nstates-1);
				for (i = 0; i < nstates; i++) {
					for (x = 0; x < nstates; x++)
						if (state & (1 << x))
							partial_lh[i] += inv_evec[i][x];
				}

			}
			partial_lh += nstates;
		} // for loop over ptn
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
						vleft += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*left->length) * partial_lh_left[i];
						vright += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*right->length) * partial_lh_right[i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i][x];
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
						vleft += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*left->length) * partial_lh_left[i];
						vright += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*right->length) * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i][x];
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
						vleft += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*left->length) * partial_lh_left[c*nstates+i];
						vright += evec[x][i] * exp(eval[i]*site_rate->getRate(c)*right->length) * partial_lh_right[c*nstates+i];
					}
					partial_lh_tmp[x] = vleft * vright;
				}
				for (i = 0; i < nstates; i++) {
					double res = 0.0;
					for (x = 0; x < nstates; x++)
						res += partial_lh_tmp[x]*inv_evec[i][x];
					partial_lh[c*nstates+i] = res;
				}
			}
			partial_lh += block;
			partial_lh_left += block;
			partial_lh_right += block;
		}
	}

	delete [] partial_lh_tmp;
}
