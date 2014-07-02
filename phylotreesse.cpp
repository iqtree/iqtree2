/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "phylotree.h"
#include "gtrmodel.h"

/* BQM: to ignore all-gapp subtree at an alignment site */
//#define IGNORE_GAP_LH

template<const int NSTATES>
inline double PhyloTree::computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    PhyloNode *node = (PhyloNode*) dad_branch->node; // Node A
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad); // Node B
    assert(node_branch);
    if (!central_partial_lh)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ptn, cat, state1, state2;
    double *partial_lh_site;
    double *partial_lh_child;
    double *trans_state;
    double p_invar = site_rate->getPInvar();
    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();
    int block = numStates * numCat;

    double p_var_cat = (1.0 - p_invar) / (double) numCat;

    EIGEN_ALIGN16 double *trans_mat_orig = new double[numCat * tranSize + 1];
    double *trans_mat = trans_mat_orig;
    if (((intptr_t) trans_mat) % 16 != 0)
        trans_mat = trans_mat + 1;
    EIGEN_ALIGN16 double state_freq[NSTATES];
    model->getStateFrequency(state_freq);
    for (cat = 0; cat < numCat; cat++) {
        double *trans_cat = trans_mat + (cat * tranSize);
        model_factory->computeTransMatrix(dad_branch->length * site_rate->getRate(cat), trans_cat);
        for (state1 = 0; state1 < NSTATES; state1++) {
            double *trans_mat_state = trans_cat + (state1 * NSTATES);
            for (state2 = 0; state2 < NSTATES; state2++)
                trans_mat_state[state2] *= state_freq[state1];
        }
    }

    double prob_const = 0.0; // probability of unobserved const patterns

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, prob_const) private(ptn, cat)
#endif
    for (ptn = 0; ptn < alnSize; ++ptn) {
        double lh_ptn = 0.0; // likelihood of the pattern
        for (cat = 0; cat < numCat; cat++) {
            partial_lh_site = node_branch->partial_lh + (ptn * block + cat * NSTATES);
            partial_lh_child = dad_branch->partial_lh + (ptn * block + cat * NSTATES);
            trans_state = trans_mat + cat * tranSize;
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_child(&partial_lh_child[0]);
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_site(&partial_lh_site[0]);
            Map<Matrix<double, NSTATES, NSTATES>, Aligned> eigen_trans_state(&trans_state[0]);
            lh_ptn += (eigen_partial_lh_child * eigen_trans_state).dot(eigen_partial_lh_site);
        }
        if (ptn < orig_alnSize) {
			lh_ptn *= p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}
			lh_ptn = log(lh_ptn);
			tree_lh += lh_ptn * (aln->at(ptn).frequency);
			_pattern_lh[ptn] = lh_ptn;
			// BQM: pattern_lh contains the LOG-likelihood, not likelihood
        } else {
			lh_ptn = lh_ptn*p_var_cat + p_invar*state_freq[(int)model_factory->unobserved_ptns[ptn-orig_alnSize]];
			prob_const += lh_ptn;

        }
    }
    if (orig_alnSize < alnSize) {
    	// ascertainment bias correction
    	prob_const = log(1.0 - prob_const);
    	for (ptn = 0; ptn < orig_alnSize; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    	tree_lh -= aln->getNSite()*prob_const;
    }

    if (pattern_lh) {
        memmove(pattern_lh, _pattern_lh, orig_alnSize * sizeof(double));
    }
    delete[] trans_mat_orig;
    return tree_lh;
}

template<int NSTATES>
void PhyloTree::computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    //double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    //double *partial_lh_block;
    //bool do_scale = true;
    //double freq;
    dad_branch->lh_scale_factor = 0.0;

    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();
    int block = numStates * numCat;
    size_t lh_size = alnSize * block;
    memset(dad_branch->scale_num, 0, alnSize * sizeof(UBYTE));

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_lh, 0, lh_size * sizeof(double));
        //double *partial_lh_site;
        for (ptn = 0; ptn < alnSize; ++ptn) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);

            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else if (ptn < orig_alnSize){
                state = (aln->at(ptn))[node->id];
            } else {
            	state = model_factory->unobserved_ptns[ptn-orig_alnSize];
            }

            if (state == STATE_UNKNOWN) {
#ifndef KEEP_GAP_LH
                dad_branch->scale_num[ptn] = -1;
#endif
                for (int state2 = 0; state2 < block; state2++) {
                    partial_lh_site[state2] = 1.0;
                }
            } else if (state < NSTATES) {
                cat = 0;
                double *_par_lh_site = partial_lh_site + state;
                while (true) {
                    *_par_lh_site = 1.0;
                    ++cat;
                    if (cat == numCat)
                        break;
                    _par_lh_site += NSTATES;
                }
            } else {
                // ambiguous character, for DNA, RNA
                state = state - (NSTATES - 1);
                for (int state2 = 0; state2 < NSTATES && state2 <= 6; state2++)
                    if (state & (1 << state2)) {
                        cat = 0;
                        double *_par_lh_site = partial_lh_site + state2;
                        while (true) {
                            *_par_lh_site = 1.0;
                            ++cat;
                            if (cat == numCat)
                                break;
                            _par_lh_site += NSTATES;
                        }
                    }
            }
        }
    } else {
        // internal node
        EIGEN_ALIGN16 double *trans_mat_orig = new double[numCat * tranSize + 2];
        double *trans_mat = trans_mat_orig;
        if (((intptr_t) trans_mat) % 16 != 0)
            trans_mat = trans_mat + 1;
        for (ptn = 0; ptn < lh_size; ++ptn)
            dad_branch->partial_lh[ptn] = 1.0;
#ifndef KEEP_GAP_LH
        for (ptn = 0; ptn < alnSize; ptn++)
            dad_branch->scale_num[ptn] = -1;
#endif
        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodSSE<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);
            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;
            for (cat = 0; cat < numCat; cat++) {
                model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), &trans_mat[cat * tranSize]);
            }
            partial_lh_site = dad_branch->partial_lh;
            partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh;
            double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, cat, partial_lh_site, partial_lh_child)
#endif
            for (ptn = 0; ptn < alnSize; ++ptn)
#ifndef KEEP_GAP_LH
            if (((PhyloNeighbor*) (*it))->scale_num[ptn] < 0) {
#ifndef _OPENMP
                partial_lh_site += NSTATES * numCat;
                partial_lh_child += NSTATES * numCat;
#endif
            } else
#endif
            {
#ifndef KEEP_GAP_LH
                if (dad_branch->scale_num[ptn] < 0)
                dad_branch->scale_num[ptn] = 0;
#endif
#ifdef _OPENMP
                int lh_offset = ptn*block;
                partial_lh_site = dad_branch->partial_lh + lh_offset;
                partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh + lh_offset;
#endif
                dad_branch->scale_num[ptn] += ((PhyloNeighbor*) (*it))->scale_num[ptn];
                double *partial_lh_block = partial_lh_site;
                double *trans_state = trans_mat;
                cat = 0;
                bool do_scale = true;
                while (true) {
                    ++cat;
                    MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                    MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                    MappedMat(NSTATES) ei_trans_state(trans_state);
                    //ei_partial_lh_site.noalias() = (ei_partial_lh_child * ei_trans_state).cwiseProduct(ei_partial_lh_site);
                    ei_partial_lh_site.array() *= (ei_partial_lh_child * ei_trans_state).array();
                    partial_lh_site += NSTATES;
                    partial_lh_child += NSTATES;
                    if (cat == numCat)
                    break;
                    else
                    trans_state += tranSize;
                }
                for (cat = 0; cat < block; cat++)
                if (partial_lh_block[cat] > SCALING_THRESHOLD) {
                    do_scale = false;
                    break;
                }
                if (do_scale) {
                    // unobserved const pattern will never have underflow
                    Map<VectorXd, Aligned> ei_lh_block(partial_lh_block, block);
                    ei_lh_block *= SCALING_THRESHOLD_INVER;
                    sum_scale += LOG_SCALING_THRESHOLD *  (*aln)[ptn].frequency;
                    dad_branch->scale_num[ptn] += 1;
                    if (pattern_scale)
                    pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
                }
            }
            dad_branch->lh_scale_factor += sum_scale;
        }
        delete[] trans_mat_orig;
    }

    dad_branch->partial_lh_computed |= 1;
}

/****************************************************************************
 computing derivatives of likelihood function
 ****************************************************************************/
template<int NSTATES>
inline double PhyloTree::computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    //assert(node_branch);
    // swap node and dad if node is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES>(node_branch, node);
    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    int cat = 0;
    double *partial_lh_site = node_branch->partial_lh;
    double *partial_lh_child = dad_branch->partial_lh;
    double lh_ptn; // likelihood of the pattern
    double lh_ptn_derv1;
    double lh_ptn_derv2;
    double derv1_frac;
    double derv2_frac;
    double *trans_state;
    double *derv1_state;
    double *derv2_state;
    double p_invar = site_rate->getPInvar();

    int numCat = site_rate->getNRate();
    int numStates = model->num_states;
    int tranSize = numStates * numStates;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();

    double p_var_cat = (1.0 - p_invar) / (double) numCat;
    double state_freq[NSTATES];
    model->getStateFrequency(state_freq);
    double *trans_mat_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
    double *trans_derv1_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
    double *trans_derv2_orig EIGEN_ALIGN16 = new double[numCat * tranSize + 1];
    // make alignment 16
    double *trans_mat = trans_mat_orig, *trans_derv1 = trans_derv1_orig, *trans_derv2 = trans_derv2_orig;
    if (((intptr_t) trans_mat) % 16 != 0)
        trans_mat = trans_mat + 1;
    if (((intptr_t) trans_derv1) % 16 != 0)
        trans_derv1 = trans_derv1 + 1;
    if (((intptr_t) trans_derv2) % 16 != 0)
        trans_derv2 = trans_derv2 + 1;

    int discrete_cat = site_rate->getNDiscreteRate();
    if (!site_rate->isSiteSpecificRate())
        for (cat = 0; cat < discrete_cat; cat++) {
            double *trans_cat = trans_mat + (cat * tranSize);
            double *derv1_cat = trans_derv1 + (cat * tranSize);
            double *derv2_cat = trans_derv2 + (cat * tranSize);
            double rate_val = site_rate->getRate(cat);
            model_factory->computeTransDervFreq(dad_branch->length, rate_val, state_freq, trans_cat, derv1_cat,
                    derv2_cat);
        }
    int dad_state = STATE_UNKNOWN;
    double my_df = 0.0;
    double my_ddf = 0.0;
    double prob_const = 0.0, prob_const_derv1 = 0.0, prob_const_derv2 = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh, my_df, my_ddf,prob_const, prob_const_derv1, prob_const_derv2) \
	private(cat, partial_lh_child, partial_lh_site,\
	lh_ptn, lh_ptn_derv1, lh_ptn_derv2, derv1_frac, derv2_frac, dad_state, trans_state, derv1_state, derv2_state)
#endif
    for (int ptn = 0; ptn < alnSize; ++ptn) {
#ifdef _OPENMP
        int lh_offset = ptn*numCat*numStates;
        partial_lh_site = node_branch->partial_lh + lh_offset;
        partial_lh_child = dad_branch->partial_lh + lh_offset;
#endif
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        int padding = 0;
        dad_state = STATE_UNKNOWN; // FOR TUNG: This is missing in your codes!
        if (dad->isLeaf()) {
        	if (ptn < orig_alnSize)
        		dad_state = (*aln)[ptn][dad->id];
        	else
        		dad_state = model_factory->unobserved_ptns[ptn-orig_alnSize];
        }
        padding = dad_state * NSTATES;
        if (dad_state < NSTATES) {
            //external node
            trans_state = trans_mat + padding;
            derv1_state = trans_derv1 + padding;
            derv2_state = trans_derv2 + padding;
            for (cat = 0; cat < numCat; cat++) {
                MappedVec(NSTATES)ei_partial_lh_child(partial_lh_child);
                MappedVec(NSTATES) ei_trans_state(trans_state);
                MappedVec(NSTATES) ei_derv1_state(derv1_state);
                MappedVec(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += ei_partial_lh_child.dot(ei_trans_state);
                lh_ptn_derv1 += ei_partial_lh_child.dot(ei_derv1_state);
                lh_ptn_derv2 += ei_partial_lh_child.dot(ei_derv2_state);
                partial_lh_child += NSTATES;
                partial_lh_site += NSTATES;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        } else {
            // internal node, or external node but ambiguous character
            trans_state = trans_mat;
            derv1_state = trans_derv1;
            derv2_state = trans_derv2;
            for (cat = 0; cat < numCat; cat++) {
                MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                MappedMat(NSTATES) ei_trans_state(trans_state);
                MappedMat(NSTATES) ei_derv1_state(derv1_state);
                MappedMat(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += (ei_partial_lh_child * ei_trans_state).dot(ei_partial_lh_site);
                lh_ptn_derv1 += (ei_partial_lh_child * ei_derv1_state).dot(ei_partial_lh_site);
                lh_ptn_derv2 += (ei_partial_lh_child * ei_derv2_state).dot(ei_partial_lh_site);
                partial_lh_site += NSTATES;
                partial_lh_child += NSTATES;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        }
        if (ptn < orig_alnSize) {
			lh_ptn = lh_ptn * p_var_cat;
			if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
				lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
			}
			double pad = p_var_cat / lh_ptn;
			if (std::isinf(pad)) {
				lh_ptn_derv1 *= p_var_cat;
				lh_ptn_derv2 *= p_var_cat;
				derv1_frac = lh_ptn_derv1 / lh_ptn;
				derv2_frac = lh_ptn_derv2 / lh_ptn;
			} else {
				derv1_frac = lh_ptn_derv1 * pad;
				derv2_frac = lh_ptn_derv2 * pad;
			}
	        double freq = aln->at(ptn).frequency;
			double tmp1 = derv1_frac * freq;
			double tmp2 = derv2_frac * freq;
			my_df += tmp1;
			my_ddf += tmp2 - tmp1 * derv1_frac;
			lh_ptn = log(lh_ptn);
			tree_lh += lh_ptn * freq;
			_pattern_lh[ptn] = lh_ptn;
        } else {
        	lh_ptn = lh_ptn*p_var_cat + p_invar*state_freq[(int)model_factory->unobserved_ptns[ptn-orig_alnSize]];
        	prob_const += lh_ptn;
        	prob_const_derv1 += lh_ptn_derv1 * p_var_cat;
        	prob_const_derv2 += lh_ptn_derv2 * p_var_cat;
        }
    }
    if (orig_alnSize < alnSize) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	derv1_frac = prob_const_derv1 / prob_const;
    	derv2_frac = prob_const_derv2 / prob_const;
    	int nsites = aln->getNSite();
    	my_df += nsites * derv1_frac;
    	my_ddf += nsites *(derv2_frac + derv1_frac*derv1_frac);
    	prob_const = log(prob_const);
    	tree_lh -= nsites * prob_const;
    	for (int ptn = 0; ptn < orig_alnSize; ptn++)
    		_pattern_lh[ptn] -= prob_const;
    }

    delete[] trans_derv2_orig;
    delete[] trans_derv1_orig;
    delete[] trans_mat_orig;
    df = my_df;
    ddf = my_ddf;
    return tree_lh;
}

template<int NSTATES>
void PhyloTree::computeThetaSSE(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    cout << "This has not been implemented yet" << endl;
    exit(1);
}

void PhyloTree::computeTheta(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeThetaSSE<2>(dad_branch, dad);
        case 4:
            return computeThetaSSE<4>(dad_branch, dad);
        case 20:
            return computeThetaSSE<20>(dad_branch, dad);
        default:
            computeThetaNaive(dad_branch, dad);
            break;
        }
    } else {
        computeThetaNaive(dad_branch, dad);
    }
}

void PhyloTree::computeThetaNaive(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    //cout << "Computing theta vector " << endl;
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if node is a leaf
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

    double *partial_lh_site = node_branch->partial_lh;
    double *partial_lh_child = dad_branch->partial_lh;
    double *theta_ptn = theta_all;
    int alnSize = aln->getNPattern();
    // BQM's question: why state_freq has size = alnSize?
    double* state_freq = new double[alnSize];
    model->getStateFrequency(state_freq);
    int num_cat = site_rate->getNRate();
    GTRModel* gtr_model = reinterpret_cast<GTRModel *>(model);
    //double* eigen_coff = gtr_model->getEigenCoeff();
    double** inv_eigen_vector = gtr_model->getInverseEigenvectors();
	double* partial_lh_child_ptn = partial_lh_child;
	double* partial_lh_site_ptn = partial_lh_site;

	int numStates = aln->num_states;
	for (int ptn = 0; ptn < alnSize; ++ptn) {
	    for (int i = 0; i < numStates; ++i) {
	    	partial_lh_child_ptn = partial_lh_child + ptn * numStates * num_cat;
	    	partial_lh_site_ptn = partial_lh_site + ptn * numStates * num_cat;
		    for (int c = 0; c < num_cat; ++c) {
				double term1 = 0;
				double term2 = 0;
				for (int x = 0; x < numStates; ++x) {
					// Compute Sigma_x pi_x u_xi L^h_a(x,c)
					term1 += inv_eigen_vector[i][x] * partial_lh_site_ptn[x];
					term2 += inv_eigen_vector[i][x] * partial_lh_child_ptn[x];
				}
	            partial_lh_child_ptn += numStates;
	            partial_lh_site_ptn += numStates;
				*theta_ptn = term1 * term2;
				//cout << "theta_ptn : " << *theta_ptn << endl;
				theta_ptn++;
			}
		}
	}
	delete[] state_freq;
}

void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    if (sse) {
    	if (params->fastSSE && getModel()->isReversible() && getModel()->getEigenvectors()) {
			switch (aln->num_states) {
			case 2:
				return computePartialLikelihoodFast<2>(dad_branch, dad, pattern_scale);
			case 4:
				return computePartialLikelihoodFast<4>(dad_branch, dad, pattern_scale);
			case 20:
				return computePartialLikelihoodFast<20>(dad_branch, dad, pattern_scale);
			/*
			case 61:
				return computePartialLikelihoodSSE<61>(dad_branch, dad, pattern_scale);
			*/
			default:
				return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
			}
    	} else {
			switch (aln->num_states) {
			case 2:
				return computePartialLikelihoodSSE<2>(dad_branch, dad, pattern_scale);
			case 4:
				return computePartialLikelihoodSSE<4>(dad_branch, dad, pattern_scale);
			case 20:
				return computePartialLikelihoodSSE<20>(dad_branch, dad, pattern_scale);
			/*
			case 61:
				return computePartialLikelihoodSSE<61>(dad_branch, dad, pattern_scale);
			*/
			default:
				return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
			}
    	}
    } else {
        return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
    }
}

void PhyloTree::initiateMyEigenCoeff() {
    assert(model);
    GTRModel* gtr_model = reinterpret_cast<GTRModel *>(model);
    double* eigen_coff = gtr_model->getEigenCoeff();
    int numStates = aln->num_states;
    if (!myEigenCoeff) {
        myEigenCoeff = new double[numStates * numStates * numStates];
    }
    int i = 0;
    for (int j = 0; j < numStates; ++j)
        for (int xa = 0; xa < numStates; ++xa)
            for (int xb = 0; xb < numStates; ++xb) {
                myEigenCoeff[i] = eigen_coff[xa * numStates * numStates + xb * numStates + j];
                ++i;
            }
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeLikelihoodBranchSSE<2>(dad_branch, dad, pattern_lh);
        case 4:
            return computeLikelihoodBranchSSE<4>(dad_branch, dad, pattern_lh);
        case 20:
            return computeLikelihoodBranchSSE<20>(dad_branch, dad, pattern_lh);
        /*
        case 61:
            return computeLikelihoodBranchSSE<61>(dad_branch, dad, pattern_lh);
        */
        default:
            return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
        }
    } else {
        return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
    }

}

/*
 * This function is called millions times. So it is not a good idea to
 * have a if and switch here.
 */
double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeLikelihoodDervSSE<2>(dad_branch, dad, df, ddf);
        case 4:
            return computeLikelihoodDervSSE<4>(dad_branch, dad, df, ddf);
        case 20:
            return computeLikelihoodDervSSE<20>(dad_branch, dad, df, ddf);
        /*
        case 61:
            return computeLikelihoodDervSSE<61>(dad_branch, dad, df, ddf);
        */
        default:
            return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
            //cout << "Wrong number of states: " << aln->num_states << endl;
            //return -1.0;
        }
    } else {
        return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
    }
}

double PhyloTree::computeLikelihoodDervFast(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    if (sse) {
        switch (aln->num_states) {
        case 2:
            return computeLikelihoodDervFastSSE<2>(dad_branch, dad, df, ddf);
        case 4:
            return computeLikelihoodDervFastSSE<4>(dad_branch, dad, df, ddf);
        case 20:
            return computeLikelihoodDervFastSSE<20>(dad_branch, dad, df, ddf);
        default:
            return computeLikelihoodDervFastNaive(dad_branch, dad, df, ddf);
            //cout << "Bad number of states: " << aln->num_states << endl;
            //exit(1);
        }
    } else {
        return computeLikelihoodDervFastNaive(dad_branch, dad, df, ddf);
    }
}

double PhyloTree::computeLikelihoodDervFastNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if dad is a leaf
    // NEW: swap if root_state is given
    if (node->isLeaf() || (node->name == ROOT_NAME && root_state != STATE_UNKNOWN)) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double lh_ptn = 0.0;
    double lh_ptn_derv1 = 0.0;
    double lh_ptn_derv2 = 0.0;
    double my_df = 0.0;
    double my_ddf = 0.0;
    double p_invar = site_rate->getPInvar();
    int ncat = site_rate->getNRate();
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    double derv1_frac;
    double derv2_frac;
    int nstates = aln->num_states;
    double* state_freq = new double[nstates];
    model->getStateFrequency(state_freq);
    int discrete_cat = site_rate->getNDiscreteRate();
    double *theta_ptn = theta_all;
    size_t nptn = aln->getNPattern();
    double t = dad_branch->length;

    double* rates = new double[discrete_cat];
    for (int i = 0; i < discrete_cat; i++) {
    	rates[i] = site_rate->getRate(i);
    }
    //double* rates = site_rate->getRates();
    double* lambda = dynamic_cast<GTRModel*>(model)->getEigenvalues();
    int block = discrete_cat * nstates;
    // array containing lambda_i * r_c
    double* lambda_r = new double[ block ];
    // array containing exp(lambda_i * r_c * t)
    double* exp_part = new double[ block ];
    // array containing square of lambda_i * r_c
    double* lambda_r_sqr = new double[ block ];

    // now initialize all the arrays (pre-computation before coming to the big loop)
    for (int i = 0; i < nstates; i++) {
    	for (int c = 0; c < discrete_cat; c++) {
    		lambda_r[ i * discrete_cat + c ] = lambda[i] * rates[c];
    		lambda_r_sqr[ i * discrete_cat + c ] = lambda_r[ i * discrete_cat + c ] * lambda_r[ i * discrete_cat + c ];
    		exp_part [ i * discrete_cat + c ] = exp( lambda_r[ i * discrete_cat + c ] * t);
    	}
    }

    int pointer_jump = nstates * discrete_cat;
    for (int ptn = 0; ptn < nptn; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        for (int i = 0; i < nstates; i++) {
        	for (int c = 0; c < discrete_cat; c++) {
        		double base = exp_part[ i * discrete_cat + c ] * theta_ptn[ i * discrete_cat + c ];
        		lh_ptn += base;
        		lh_ptn_derv1 +=  lambda_r[ i * discrete_cat + c ] * base;
        		lh_ptn_derv2 +=  lambda_r_sqr[ i * discrete_cat + c ] * base;
        	}
        }
        theta_ptn += pointer_jump;
        lh_ptn = lh_ptn * p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
            lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
        }
        double pad = p_var_cat / lh_ptn;
        if (std::isinf(pad)) {
            lh_ptn_derv1 *= p_var_cat;
            lh_ptn_derv2 *= p_var_cat;
            derv1_frac = lh_ptn_derv1 / lh_ptn;
            derv2_frac = lh_ptn_derv2 / lh_ptn;
        } else {
            derv1_frac = lh_ptn_derv1 * pad;
            derv2_frac = lh_ptn_derv2 * pad;
        }
        double freq = (*aln)[ptn].frequency;
        double tmp1 = derv1_frac * freq;
        double tmp2 = derv2_frac * freq;
        my_df += tmp1;
        my_ddf += tmp2 - tmp1 * derv1_frac;
        lh_ptn = log(lh_ptn);
        //cout << lh_ptn << endl;
        //cout << "lh_ptn = " << lh_ptn << endl;
        tree_lh += lh_ptn * aln->at(ptn).frequency;
        //cout << tree_lh << endl;
        _pattern_lh[ptn] = lh_ptn;
    }
    delete [] lambda_r;
    delete [] lambda_r_sqr;
    delete [] exp_part;
    delete [] rates;
    delete [] state_freq;
    df = my_df;
    ddf = my_ddf;
    return tree_lh;
}

template<int NSTATES>
double PhyloTree::computeLikelihoodDervFastSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if node is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    double lh_ptn = 0.0;
    double lh_ptn_derv1 = 0.0;
    double lh_ptn_derv2 = 0.0;
    double p_invar = site_rate->getPInvar();
    int numCat = site_rate->getNRate();
    double p_var_cat = (1.0 - p_invar) / (double) numCat;
    double derv1_frac;
    double derv2_frac;
    double state_freq[NSTATES];
    int discrete_cat = site_rate->getNDiscreteRate();
    model->getStateFrequency(state_freq);
    double* rates = new double[discrete_cat];
    for (int i = 0; i < numCat; i++) {
    	rates[i] = site_rate->getRate(i);
    }
    Map<Matrix<double, Dynamic, 1>, Aligned> ei_rates(rates, discrete_cat);
    Map<Matrix<double, 1, NSTATES>, Aligned> ei_eigenvalues(dynamic_cast<GTRModel*>(model)->getEigenvalues());
    Matrix(NSTATES)ei_rates_times_lambdas = ei_rates * ei_eigenvalues;
    Matrix(NSTATES)expo_time =
    (ei_rates_times_lambdas * dad_branch->length).array().exp();
    Matrix(NSTATES)expo_time_derv1 = expo_time.cwiseProduct(
            ei_rates_times_lambdas);
    Matrix(NSTATES)expo_time_derv2 = expo_time_derv1.cwiseProduct(
            ei_rates_times_lambdas);
    double *theta_ptn = theta_all;
    int num_patterns = aln->getNPattern();
    for (int ptn = 0; ptn < num_patterns; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        /*
         for (int k = 0; k < discrete_cat; k++) {
         for (int i = 0; i < num_states; i++) {
         double x =theta_ptn[k*num_states+i] *
         exp(site_rate->getRates()[k]*dad_branch->length*dynamic_cast<GTRModel*> (model)->getEigenvalues()[i]);
         //cout << exp(site_rate->getRates()[k]*dad_branch->length*dynamic_cast<GTRModel*> (model)->getEigenvalues()[i]) << " ";
         lh_ptn += x;
         //cout << x << " ";
         //cout << theta_ptn[k*num_states+i] << " ";
         double y = x * site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i];
         //cout << y << " ";
         //cout << site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i] << " ";
         lh_ptn_derv1 += y;
         lh_ptn_derv2 += y*site_rate->getRates()[k] * dynamic_cast<GTRModel*> (model)->getEigenvalues()[i];

         }
         }
         */

        Map<Matrix<double, Dynamic, Dynamic, RowMajor>, Aligned> ei_theta_ptn(theta_ptn, discrete_cat, NSTATES);
        //Map<Matrix<double, NSTATES, NSTATES, RowMajor>, Aligned> ei_theta_ptn(theta_ptn);
        //cout << "ei_theta_ptn" << endl;
        //cout << ei_theta_ptn << endl;
        //ArrayXXd ei_theta_times_expo = ei_theta_ptn * expo_time;

        lh_ptn = ei_theta_ptn.cwiseProduct(expo_time).sum();
        lh_ptn_derv1 = ei_theta_ptn.cwiseProduct(expo_time_derv1).sum();
        lh_ptn_derv2 = ei_theta_ptn.cwiseProduct(expo_time_derv2).sum();

        /*
         for (int cat = 0; cat < discrete_cat; ++cat) {
         MappedRowVec(NSTATES) ei_theta_ptn_cat(theta_ptn);
         lh_ptn += ei_theta_ptn_cat.dot(expo_time.row(cat));
         lh_ptn_derv1 += ei_theta_ptn_cat.dot(expo_time_derv1.row(cat));
         lh_ptn_derv2 += ei_theta_ptn_cat.dot(expo_time_derv2.row(cat));
         theta_ptn += NSTATES;
         }
         */

        /*
         ArrayXXd ei_theta_times_expo_derv1 = ei_rates_times_lambdas * ei_theta_times_expo;
         ArrayXXd ei_theta_times_expo_derv1 = ei_theta_ptn.cwiseProduct(expo_time_derv1);

         ArrayXXd ei_theta_times_expo_derv2 = ei_theta_ptn.cwiseProduct(expo_time_derv2);
         lh_ptn = ei_theta_times_expo.sum();
         lh_ptn_derv1 = ei_theta_times_expo_derv1.sum();
         lh_ptn_derv2 = (ei_theta_times_expo_derv1 * ei_rates_times_lambdas).sum();
         */
        theta_ptn += NSTATES * discrete_cat;

        lh_ptn = lh_ptn * p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
            lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
        }

        double pad = p_var_cat / lh_ptn;
        if (std::isinf(pad)) {
            lh_ptn_derv1 *= p_var_cat;
            lh_ptn_derv2 *= p_var_cat;
            derv1_frac = lh_ptn_derv1 / lh_ptn;
            derv2_frac = lh_ptn_derv2 / lh_ptn;
        } else {
            derv1_frac = lh_ptn_derv1 * pad;
            derv2_frac = lh_ptn_derv2 * pad;
        }
        double tmp1 = derv1_frac * aln->at(ptn).frequency;
        double tmp2 = derv2_frac * aln->at(ptn).frequency;
        df += tmp1;
        ddf += tmp2 - tmp1 * derv1_frac;
        lh_ptn = log(lh_ptn);
        tree_lh += lh_ptn * aln->at(ptn).frequency;
        _pattern_lh[ptn] = lh_ptn;
    }
    return tree_lh;

}

template<int NSTATES>
void PhyloTree::computePartialLikelihoodFast(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    //double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    //double *partial_lh_block;
    //bool do_scale = true;
    //double freq;
    dad_branch->lh_scale_factor = 0.0;
    double **eigen_vec = model->getEigenvectors();
    double **inv_eigen_vec = model->getInverseEigenvectors();
    double *eigen_val = model->getEigenvalues();

    int numCat = site_rate->getNRate();
    //int numStates = model->num_states;
    const int tranSize = NSTATES * NSTATES;
    int alnSize = aln->size() + model_factory->unobserved_ptns.size();
    int orig_alnSize = aln->size();
    int block = NSTATES * numCat;
    size_t lh_size = alnSize * block;
    memset(dad_branch->scale_num, 0, alnSize * sizeof(UBYTE));

    double *sum_inv_eigen_vec = new double[NSTATES];
    memset(sum_inv_eigen_vec, 0, NSTATES*sizeof(double));
    {
    for (int state2 = 0; state2 < NSTATES; state2++) {
    	for (int state3 = 0; state3 < NSTATES; state3++) {
    		sum_inv_eigen_vec[state2] += inv_eigen_vec[state2][state3];
    	}
    }
    }

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_lh, 0, lh_size * sizeof(double));
        //double *partial_lh_site;
        for (ptn = 0; ptn < alnSize; ++ptn) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);

            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else if (ptn < orig_alnSize){
                state = (aln->at(ptn))[node->id];
            } else {
            	state = model_factory->unobserved_ptns[ptn-orig_alnSize];
            }

            if (state == STATE_UNKNOWN) {
#ifndef KEEP_GAP_LH
                dad_branch->scale_num[ptn] = -1;
#endif
                for (cat = 0; cat < numCat; cat++) {
                	memcpy(partial_lh_site + cat*NSTATES, sum_inv_eigen_vec, NSTATES*sizeof(double));
                }
            } else if (state < NSTATES) {
            	for (int state2 = 0; state2 < NSTATES; state2++)
            		partial_lh_site[state2] = inv_eigen_vec[state2][(int)state];
                for (cat = 1; cat < numCat; cat++) {
                	memcpy(partial_lh_site + cat*NSTATES, partial_lh_site, NSTATES*sizeof(double));
                }
            } else {
                // ambiguous character, for DNA, RNA
                state = state - (NSTATES - 1);
                for (int state2 = 0; state2 < NSTATES; state2++) {
                	partial_lh_site[state2] = 0.0;
					for (int state3 = 0; state3 < NSTATES && state3 <= 6; state3++)
						if (state & (1 << state3)) {
							partial_lh_site[state2] += inv_eigen_vec[state2][state3];
						}
                }
                for (cat = 1; cat < numCat; cat++) {
                	memcpy(partial_lh_site + cat*NSTATES, partial_lh_site, NSTATES*sizeof(double));
                }
            }
        }
    } else {
        // internal node
        EIGEN_ALIGN16 double *trans_mat_orig = new double[numCat * tranSize + 2];
        double *trans_mat = trans_mat_orig;
        if (((intptr_t) trans_mat) % 16 != 0)
            trans_mat = trans_mat + 1;
        for (ptn = 0; ptn < lh_size; ++ptn)
            dad_branch->partial_lh[ptn] = 1.0;
#ifndef KEEP_GAP_LH
        for (ptn = 0; ptn < alnSize; ptn++)
            dad_branch->scale_num[ptn] = -1;
#endif
        FOR_NEIGHBOR_IT(node, dad, it)if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodFast<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);
            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;
            for (cat = 0; cat < numCat; cat++) {
                //model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), &trans_mat[cat * tranSize]);
            	for (int state = 0; state < NSTATES; state++)
            		for (int state2 = 0; state2 < NSTATES; state2++)
            			trans_mat[cat*tranSize + state*NSTATES + state2] =
            				eigen_vec[state][state2] * exp(eigen_val[state2] * site_rate->getRate(cat) * (*it)->length);
            }
            partial_lh_site = dad_branch->partial_lh;
            partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh;
            double sum_scale = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sum_scale) private(ptn, cat, partial_lh_site, partial_lh_child)
#endif
            for (ptn = 0; ptn < alnSize; ++ptn)
#ifndef KEEP_GAP_LH
            if (((PhyloNeighbor*) (*it))->scale_num[ptn] < 0) {
	#ifndef _OPENMP
                partial_lh_site += NSTATES * numCat;
                partial_lh_child += NSTATES * numCat;
	#endif // _OPENMP
            } else
#endif // KEEP_GAP_LH
            {
#ifndef KEEP_GAP_LH
                if (dad_branch->scale_num[ptn] < 0)
                dad_branch->scale_num[ptn] = 0;
#endif
#ifdef _OPENMP
                int lh_offset = ptn*block;
                partial_lh_site = dad_branch->partial_lh + lh_offset;
                partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh + lh_offset;
#endif
                dad_branch->scale_num[ptn] += ((PhyloNeighbor*) (*it))->scale_num[ptn];
                double *partial_lh_block = partial_lh_site;
                double *trans_state = trans_mat;
                bool do_scale = true;
                for (cat=0; cat < numCat; cat++) {
                    MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                    MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                    MappedMat(NSTATES) ei_trans_state(trans_state);
                    //ei_partial_lh_site.noalias() = (ei_partial_lh_child * ei_trans_state).cwiseProduct(ei_partial_lh_site);
                    ei_partial_lh_site.array() *= (ei_partial_lh_child * ei_trans_state).array();
                    partial_lh_site += NSTATES;
                    partial_lh_child += NSTATES;
                    trans_state += tranSize;
                }
                for (cat = 0; cat < block; cat++)
                if (partial_lh_block[cat] > SCALING_THRESHOLD) {
                    do_scale = false;
                    break;
                }
                if (do_scale) {
                    // unobserved const pattern will never have underflow
                    Map<VectorXd, Aligned> ei_lh_block(partial_lh_block, block);
                    ei_lh_block *= SCALING_THRESHOLD_INVER;
                    sum_scale += LOG_SCALING_THRESHOLD *  (*aln)[ptn].frequency;
                    dad_branch->scale_num[ptn] += 1;
                    if (pattern_scale)
                    pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
                }
            }
            dad_branch->lh_scale_factor += sum_scale;
        }
        delete[] trans_mat_orig;
    }

    dad_branch->partial_lh_computed |= 1;
    delete [] sum_inv_eigen_vec;
}
