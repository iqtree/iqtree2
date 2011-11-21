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
//#include "acml_mv/acml_mv.h"

template<int NSTATES>
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
        computePartialLikelihoodSSE<NSTATES> (dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES> (node_branch, node);
    // now combine likelihood at the branch

    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ptn, cat, state1, state2;
    double *partial_lh_site;
    double *partial_lh_child;
    double *trans_state;
    double p_invar = site_rate->getPInvar();
    double p_var_cat = (1.0 - p_invar) / (double) numCat;    

    EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
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

        lh_ptn *= p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
            lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
        }                
        lh_ptn += p_invar_ptns[ptn];
        //lh_ptns[ptn] = lh_ptn;
        tree_lh += log(lh_ptn) * ptn_freqs[ptn];
    }
    //vrda_log(alnSize, lh_ptns.data(), lh_ptns_log.data());
    //tree_lh += (lh_ptns_log * ptn_freqs).sum();
    if (pattern_lh) {
		for (ptn = 0; ptn < alnSize; ++ptn)
			pattern_lh[ptn] = lh_ptns_log[ptn];
	}
    return tree_lh;
}

template<int NSTATES>
inline void PhyloTree::computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    double *partial_lh_block;
    bool do_scale = true;
    double freq;
    dad_branch->lh_scale_factor = 0.0;

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_lh, 0, lh_size * sizeof (double));
        //double *partial_lh_site;
        for (ptn = 0; ptn < alnSize; ++ptn) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);

            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else {                
                state = (aln->at(ptn))[node->id];
            }

            if (state == STATE_UNKNOWN) {
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
                for (int state2 = 0; state2 < NSTATES; state2++)
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
        EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
        for (ptn = 0; ptn < lh_size; ++ptn)
            dad_branch->partial_lh[ptn] = 1.0;

        FOR_NEIGHBOR_IT(node, dad, it)
        if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodSSE<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);
            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;
            for (cat = 0; cat < numCat; cat++) {
                model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), &trans_mat[cat * tranSize]);
            }
            partial_lh_site = dad_branch->partial_lh;
            partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh;            
            for (ptn = 0; ptn < alnSize; ++ptn) {
                partial_lh_block = partial_lh_site;
                freq = ptn_freqs[ptn];
                trans_state = trans_mat;
                cat = 0;
                do_scale = true;
                while (true) {
                    ++cat;
                    MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                    MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                    MappedMat(NSTATES) ei_trans_state(trans_state);
                    ei_partial_lh_site.noalias() = (ei_partial_lh_child * ei_trans_state).cwiseProduct(ei_partial_lh_site);
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
                    Map<ArrayXd, Aligned> ei_lh_block(partial_lh_block, block);
                    ei_lh_block *= SCALING_THRESHOLD_INVER;                    
                    dad_branch->lh_scale_factor += LOG_SCALING_THRESHOLD * freq;
                    if (pattern_scale)
                        pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
                }                
            }
        }
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
        computePartialLikelihoodSSE<NSTATES > (dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES > (node_branch, node);
    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    int ptn = 0;
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
    p_invar = site_rate->getPInvar();
    p_var_cat = (1.0 - p_invar) / (double) numCat;
    double state_freq[NSTATES];
    model->getStateFrequency(state_freq);
    EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
    EIGEN_ALIGN16 double trans_derv1[numCat * tranSize];
    EIGEN_ALIGN16 double trans_derv2[numCat * tranSize];
    int discrete_cat = site_rate->getNDiscreteRate();
    if (!site_rate->isSiteSpecificRate())
        for (cat = 0; cat < discrete_cat; cat++) {
            //trans_mat[cat] = model->newTransMatrix();
            double *trans_cat = trans_mat + (cat * tranSize);
            double *derv1_cat = trans_derv1 + (cat * tranSize);
            double *derv2_cat = trans_derv2 + (cat * tranSize);
            double rate_val = site_rate->getRate(cat);
            //double rate_sqr = rate_val * rate_val;
            model_factory->computeTransDervFreq(dad_branch->length, rate_val, state_freq, trans_cat, derv1_cat, derv2_cat);
        }
    bool not_ptn_cat = (site_rate->getPtnCat(0) < 0);
    int dad_state = STATE_UNKNOWN;
    for ( ; ptn < alnSize; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        double freq = ptn_freqs[ptn];
        int padding = 0;
        if (dad->isLeaf()) {
            dad_state = (*aln)[ptn][dad->id];
            padding = dad_state * NSTATES;
        }
        if (dad_state < NSTATES) {
            //external node
            trans_state = trans_mat + padding;
            derv1_state = trans_derv1 + padding;
            derv2_state = trans_derv2 + padding;
            cat = 0;
            while (true) {
                ++cat;
                MappedVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                MappedVec(NSTATES) ei_trans_state(trans_state);
                MappedVec(NSTATES) ei_derv1_state(derv1_state);
                MappedVec(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += ei_partial_lh_child.dot(ei_trans_state);
                lh_ptn_derv1 += ei_partial_lh_child.dot(ei_derv1_state);
                lh_ptn_derv2 += ei_partial_lh_child.dot(ei_derv2_state);
                partial_lh_child += NSTATES;
                partial_lh_site += NSTATES;
                if ( cat == numCat )
                    break;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;                
            }
        } else {
            // internal node, or external node but ambiguous character
            cat = 0;
            trans_state = trans_mat;
            derv1_state = trans_derv1;
            derv2_state = trans_derv2;            
            while (true) {
                ++cat;
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
                if (cat == numCat)
                    break;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        }
        lh_ptn = lh_ptn * p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < NSTATES) {
            lh_ptn += p_invar * state_freq[(int) (*aln)[ptn][0]];
        }        
        double tmp = p_var_cat / lh_ptn;
        derv1_frac = lh_ptn_derv1 * tmp;
        derv2_frac = lh_ptn_derv2 * tmp;
        df += derv1_frac * freq;
        ddf += (derv2_frac - derv1_frac * derv1_frac) * freq;
        //lh_ptns[ptn] = lh_ptn;
        tree_lh += log(lh_ptn) * freq;
    }
    //vrda_log(alnSize, lh_ptns.data(), lh_ptns_log.data());q
    //tree_lh += (lh_ptns_log * ptn_freqs).sum();
    if (!isnormal(tree_lh)){
    	outError("Log-likelihood is NaN. Please contact the author");
    }
    return tree_lh;
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    if (sse) {
        switch (aln->num_states) {            
            case 2: return computeLikelihoodBranchSSE <2> (dad_branch, dad, pattern_lh);
            case 4: return computeLikelihoodBranchSSE <4> (dad_branch, dad, pattern_lh);
            case 20: return computeLikelihoodBranchSSE <20> (dad_branch, dad, pattern_lh);
            default: return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
        }
    } else {
        return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
    }

}


void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    if (sse) {
        switch (aln->num_states) {
            case 2: return computePartialLikelihoodSSE < 2 > (dad_branch, dad, pattern_scale);
            case 4: return computePartialLikelihoodSSE < 4 > (dad_branch, dad, pattern_scale);
            //case 4: return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
            case 20:return computePartialLikelihoodSSE < 20 > (dad_branch, dad, pattern_scale);
            default: return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
        }
    } else {
        return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
    }    
}

double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    if (sse) {
        switch (aln->num_states) {
            case 2: return computeLikelihoodDervSSE < 2 > (dad_branch, dad, df, ddf);
            case 4: return computeLikelihoodDervSSE < 4 > (dad_branch, dad, df, ddf);
            case 20:return computeLikelihoodDervSSE < 20 > (dad_branch, dad, df, ddf);
            default: return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
        }
    } else {
        return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
    }
}
