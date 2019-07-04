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
#include "alignmentpairwise.h"
#include "tree/phylosupertree.h"

AlignmentPairwise::AlignmentPairwise()
        : Alignment(), Optimization()
{
    pair_freq = NULL;
}

AlignmentPairwise::AlignmentPairwise(PhyloTree *atree, int seq1, int seq2) : Alignment(), Optimization() {
    tree = atree;
    seq_id1 = seq1;
    seq_id2 = seq2;
    num_states = tree->aln->num_states;
    STATE_UNKNOWN = tree->aln->STATE_UNKNOWN;
    pair_freq = NULL;

    if (tree->getRate()->isSiteSpecificRate() || tree->getModel()->isSiteSpecificModel()) return;

    // categorized rates
    if (tree->getRate()->getPtnCat(0) >= 0) {
        int size_sqr = num_states * num_states;
        int total_size = size_sqr * tree->getRate()->getNDiscreteRate();
        pair_freq = new double[total_size];
        memset(pair_freq, 0, sizeof(double)*total_size);
        int i = 0;
        for (Alignment::iterator it = tree->aln->begin(); it != tree->aln->end(); it++, i++) {
            int state1 = tree->aln->convertPomoState((*it)[seq_id1]);
            int state2 = tree->aln->convertPomoState((*it)[seq_id2]);
            addPattern(state1, state2, it->frequency, tree->getRate()->getPtnCat(i));
            /*
            	if (state1 < num_states && state2 < num_states)
            		pair_freq[tree->getRate()->getPtnCat(i)*size_sqr + state1*num_states + state2] += it->frequency;*/
        }
        return;
    }

    pair_freq = new double[num_states * num_states];
    memset(pair_freq, 0, sizeof(double) * num_states * num_states);
    for (Alignment::iterator it = tree->aln->begin(); it != tree->aln->end(); it++) {
        int state1 = tree->aln->convertPomoState((*it)[seq_id1]);
        int state2 = tree->aln->convertPomoState((*it)[seq_id2]);
        addPattern(state1, state2, it->frequency);
        /*		if (state1 < num_states && state2 < num_states)
        			pair_freq[state1 * num_states + state2] += it->frequency;*/
    }
}

bool AlignmentPairwise::addPattern(int state1, int state2, int freq, int cat) {
    int i;
    if (state1 == STATE_UNKNOWN || state2 == STATE_UNKNOWN) return true;

    double *pair_pos = pair_freq + (cat*num_states*num_states);
    // unambiguous case
    if (state1 < num_states && state2 < num_states) {
        pair_pos[state1*num_states + state2] += freq;
        return false;
    }

    return true;

    if (state1 < num_states) {
        // ambiguous character, for DNA, RNA
        state2 = state2 - (num_states - 1);
        for (i = 0; i < num_states; i++)
            if (state2 & (1 << i))
                pair_pos[state1*num_states + i] += freq;
        return false;
    }

    if (state2 < num_states) {
        // ambiguous character, for DNA, RNA
        state1 = state1 - (num_states - 1);
        for (i = 0; i < num_states; i++)
            if (state1 & (1 << i))
                pair_pos[i*num_states + state2] += freq;
        return false;
    }

    return true;
}

double AlignmentPairwise::computeFunction(double value) {

    RateHeterogeneity *site_rate = tree->getRate();
    int ncat = site_rate->getNDiscreteRate();
    ModelSubst *model = tree->getModel();
    int trans_size = tree->getModel()->getTransMatrixSize();
    int cat, i;
    int nptn = tree->aln->getNPattern();
    double lh = 0.0;

    // site-specific rates
    if (site_rate->isSiteSpecificRate()) {
        for (i = 0; i < nptn; i++) {
            int state1 = tree->aln->at(i)[seq_id1];
            int state2 = tree->aln->at(i)[seq_id2];
            if (state1 >= num_states || state2 >= num_states) continue;
            double trans = tree->getModelFactory()->computeTrans(value * site_rate->getPtnRate(i), state1, state2);
            lh -= log(trans) * tree->aln->at(i).frequency;

        }
        return lh;
    }

    if (tree->getModel()->isSiteSpecificModel()) {
        for (i = 0; i < nptn; i++) {
            int state1 = tree->aln->at(i)[seq_id1];
            int state2 = tree->aln->at(i)[seq_id2];
            if (state1 >= num_states || state2 >= num_states) continue;
            double trans = tree->getModel()->computeTrans(value, model->getPtnModelID(i), state1, state2);
            lh -= log(trans) * tree->aln->at(i).frequency;

        }
		return lh;
	}
    
    double *trans_mat = new double[trans_size];

    // categorized rates
    if (site_rate->getPtnCat(0) >= 0) {
        for (cat = 0; cat < ncat; cat++) {
            tree->getModelFactory()->computeTransMatrix(value*site_rate->getRate(cat), trans_mat);
            double *pair_pos = pair_freq + cat*trans_size;
            for (i = 0; i < trans_size; i++) if (pair_pos[i] > Params::getInstance().min_branch_length) {
                    if (trans_mat[i] <= 0) throw "Negative transition probability";
                    lh -= pair_pos[i] * log(trans_mat[i]);
                }
        }
        delete [] trans_mat;
        return lh;
    }

    double *sum_trans_mat = new double[trans_size];

    if (tree->getModelFactory()->site_rate->getGammaShape() == 0.0)
        tree->getModelFactory()->computeTransMatrix(value, sum_trans_mat);
    else {
        tree->getModelFactory()->computeTransMatrix(value * site_rate->getRate(0), sum_trans_mat);
        for (cat = 1; cat < ncat; cat++) {
            tree->getModelFactory()->computeTransMatrix(value * site_rate->getRate(cat), trans_mat);
            for (i = 0; i < trans_size; i++)
                sum_trans_mat[i] += trans_mat[i];
        }
    }
    for (i = 0; i < trans_size; i++) {
        lh -= pair_freq[i] * log(sum_trans_mat[i]);
    }
    delete [] sum_trans_mat;
    delete [] trans_mat;
    // negative log-likelihood (for minimization)
    return lh;
}

void AlignmentPairwise::computeFuncDerv(double value, double &df, double &ddf) {
    RateHeterogeneity *site_rate = tree->getRate();
    int ncat = site_rate->getNDiscreteRate();
    ModelSubst *model = tree->getModel();
    int trans_size = tree->getModel()->getTransMatrixSize();
    int cat, i;
    int nptn = tree->aln->getNPattern();
//    double lh = 0.0;
    df = 0.0;
    ddf = 0.0;

    if (site_rate->isSiteSpecificRate()) {
        for (i = 0; i < nptn; i++) {
            int state1 = tree->aln->at(i)[seq_id1];
            int state2 = tree->aln->at(i)[seq_id2];
            if (state1 >= num_states || state2 >= num_states) continue;
            double rate_val = site_rate->getPtnRate(i);
            double rate_sqr = rate_val * rate_val;
            double derv1, derv2;
            double trans = tree->getModelFactory()->computeTrans(value * rate_val, state1, state2, derv1, derv2);
//            lh -= log(trans) * tree->aln->at(i).frequency;
            double d1 = derv1 / trans;
            df -= rate_val * d1 * tree->aln->at(i).frequency;
            ddf -= rate_sqr * (derv2/trans - d1*d1) * tree->aln->at(i).frequency;

        }
//        return lh;
        return;
    }

    
    if (tree->getModel()->isSiteSpecificModel()) {
        for (i = 0; i < nptn; i++) {
            int state1 = tree->aln->at(i)[seq_id1];
            int state2 = tree->aln->at(i)[seq_id2];
            if (state1 >= num_states || state2 >= num_states) continue;
            double rate_val = site_rate->getPtnRate(i);
            double rate_sqr = rate_val * rate_val;
            double derv1, derv2;
            double trans = tree->getModel()->computeTrans(value * rate_val,model->getPtnModelID(i), state1, state2, derv1, derv2);
//            lh -= log(trans) * tree->aln->at(i).frequency;
            double d1 = derv1 / trans;
            df -= rate_val * d1 * tree->aln->at(i).frequency;
            ddf -= rate_sqr * (derv2/trans - d1*d1) * tree->aln->at(i).frequency;

        }
//        return lh;
        return;
    }

    double *trans_mat = new double[trans_size];
	double *trans_derv1 = new double[trans_size];
	double *trans_derv2 = new double[trans_size];

    // categorized rates
    if (site_rate->getPtnCat(0) >= 0) {
        for (cat = 0; cat < ncat; cat++) {
            double rate_val = site_rate->getRate(cat);
            double derv1 = 0.0, derv2 = 0.0;
            tree->getModelFactory()->computeTransDerv(value*rate_val, trans_mat, trans_derv1, trans_derv2);
            double *pair_pos = pair_freq + cat*trans_size;
            for (i = 0; i < trans_size; i++) if (pair_pos[i] > 0) {
                    if (trans_mat[i] <= 0) throw "Negative transition probability";
                    double d1 = trans_derv1[i] / trans_mat[i];
                    derv1 += pair_pos[i] * d1;
                    derv2 += pair_pos[i] * (trans_derv2[i]/trans_mat[i] - d1 * d1);
//                    lh -= pair_pos[i] * log(trans_mat[i]);
                }
            df -= derv1 * rate_val;
            ddf -= derv2 * rate_val * rate_val;
        }
        delete [] trans_derv2;
		delete [] trans_derv1;
		delete [] trans_mat;
//        return lh;
        return;
    }


    double *sum_trans = new double[trans_size];
	double *sum_derv1 = new double[trans_size];
	double *sum_derv2 = new double[trans_size];
    memset(sum_trans, 0, sizeof(double) * trans_size);
    memset(sum_derv1, 0, sizeof(double) * trans_size);
    memset(sum_derv2, 0, sizeof(double) * trans_size);

    for (cat = 0; cat < ncat; cat++) {
        double rate_val = site_rate->getRate(cat);
        double prop_val = site_rate->getProp(cat);
        if (tree->getModelFactory()->site_rate->getGammaShape() == 0.0)
            rate_val = 1.0;

        double coeff1 = rate_val * prop_val;
        double coeff2 = rate_val * coeff1;
        tree->getModelFactory()->computeTransDerv(value * rate_val, trans_mat, trans_derv1, trans_derv2);
        for (i = 0; i < trans_size; i++) {
            sum_trans[i] += trans_mat[i] * prop_val;
            sum_derv1[i] += trans_derv1[i] * coeff1;
            sum_derv2[i] += trans_derv2[i] * coeff2;
        }
    }
    
    // 2019-07-03: incorporate p_invar
    double p_invar = site_rate->getPInvar();
    if (p_invar > 0.0)
        for (i = 0; i < num_states; i++)
            sum_trans[i*num_states+i] += p_invar;
    
    for (i = 0; i < trans_size; i++)
        if (pair_freq[i] > Params::getInstance().min_branch_length && sum_trans[i] > 0.0) {
//            lh -= pair_freq[i] * log(sum_trans[i]);
            double d1 = sum_derv1[i] / sum_trans[i];
            df -= pair_freq[i] * d1;
            ddf -= pair_freq[i] * (sum_derv2[i]/sum_trans[i] - d1 * d1);
        }
    delete [] sum_derv2;
	delete [] sum_derv1;
	delete [] sum_trans;
	delete [] trans_derv2;
	delete [] trans_derv1;
	delete [] trans_mat;
    // negative log-likelihood (for minimization)
//    return lh;
    return;
}

double AlignmentPairwise::optimizeDist(double initial_dist, double &d2l) {
    // initial guess of the distance using Juke-Cantor correction
    double dist = initial_dist;

    d2l = -1.0;

    // if no model or rate is specified, return the JC distance and set variance to const
    if (!tree->getModelFactory() || !tree->getRate()) return dist;

    double negative_lh, ferror;
    double max_genetic_dist = MAX_GENETIC_DIST;
    if (tree->aln->seq_type == SEQ_POMO) {
        int N = tree->aln->virtual_pop_size;
        max_genetic_dist *= N*N;
    }
    if (tree->optimize_by_newton) // Newton-Raphson method
        dist = minimizeNewton(Params::getInstance().min_branch_length, dist, max_genetic_dist, Params::getInstance().min_branch_length, d2l);
    else // Brent method
        dist = minimizeOneDimen(Params::getInstance().min_branch_length, dist, max_genetic_dist, Params::getInstance().min_branch_length, &negative_lh, &ferror);

    return dist;
}

double AlignmentPairwise::optimizeDist(double initial_dist) {
	double d2l;
	return optimizeDist(initial_dist, d2l);
}


AlignmentPairwise::~AlignmentPairwise()
{
    if (pair_freq) delete [] pair_freq;
}


