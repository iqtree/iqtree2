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

AlignmentPairwise::AlignmentPairwise()
 : Alignment(), Optimization()
{
	pair_freq = NULL;
}

AlignmentPairwise::AlignmentPairwise(PhyloTree *atree, int seq_id1, int seq_id2) : Alignment(), Optimization() {
	tree = atree;
	num_states = tree->aln->num_states;
	pair_freq = new int[num_states * num_states];
	memset(pair_freq, 0, sizeof(int) * num_states * num_states);
	for (Alignment::iterator it = tree->aln->begin(); it != tree->aln->end(); it++) {
		int state1 = (*it)[seq_id1];
		int state2 = (*it)[seq_id2];
		if (state1 < num_states && state2 < num_states)
			pair_freq[state1 * num_states + state2] += it->frequency;
	}
}

double AlignmentPairwise::computeFunction(double value) {

	RateHeterogeneity *site_rate = tree->getRate();
	int ncat = site_rate->getNRate();
	int trans_size = num_states * num_states;
	int cat, i;

	double trans_mat[trans_size];
	double sum_trans_mat[trans_size];

	tree->getModelFactory()->computeTransMatrix(value * site_rate->getRate(0), sum_trans_mat);
	for (cat = 1; cat < ncat; cat++) {
		tree->getModelFactory()->computeTransMatrix(value * site_rate->getRate(cat), trans_mat);
		for (i = 0; i < trans_size; i++)
			sum_trans_mat[i] += trans_mat[i];
	}
	double lh = 0.0;
	for (i = 0; i < trans_size; i++) {
		lh += pair_freq[i] * log(sum_trans_mat[i]);
	}
	// negative log-likelihood (for minimization)
	return -lh;
}

double AlignmentPairwise::computeFuncDerv(double value, double &df, double &ddf) {
	RateHeterogeneity *site_rate = tree->getRate();
	int ncat = site_rate->getNRate();
	int trans_size = num_states * num_states;
	int cat, i;

	double trans_mat[trans_size], trans_derv1[trans_size], trans_derv2[trans_size];
	double sum_trans[trans_size],sum_derv1[trans_size], sum_derv2[trans_size];

	tree->getModelFactory()->computeTransDerv(value * site_rate->getRate(0), sum_trans, sum_derv1, sum_derv2);
	for (cat = 1; cat < ncat; cat++) {
		tree->getModelFactory()->computeTransDerv(value * site_rate->getRate(cat), trans_mat, trans_derv1, trans_derv2);
		for (i = 0; i < trans_size; i++) {
			sum_trans[i] += trans_mat[i];
			sum_derv1[i] += trans_derv1[i];
			sum_derv2[i] += trans_derv2[i];
		}
	}
	double lh = 0.0;
	df = 0.0; ddf = 0.0;
	for (i = 0; i < trans_size; i++) {
		lh += pair_freq[i] * log(sum_trans[i]);
		double d1 = sum_derv1[i] / sum_trans[i];
		df += pair_freq[i] * d1;
		ddf += pair_freq[i] * (sum_derv2[i]/sum_trans[i] - d1 * d1);
	}
	// negative log-likelihood (for minimization)
	df = -df; ddf = -ddf;
	return -lh;
}

double AlignmentPairwise::optimizeDist(double initial_dist) {
	// initial guess of the distance using Juke-Cantor correction
	double dist = initial_dist; 

	// if no model or rate is specified, return the JC distance
	if (!tree->getModelFactory() || !tree->getRate()) return dist;

	double negative_lh, ferror;
	//if (tree->optimize_by_newton) // Newton-Raphson method
		//dist = minimizeNewton(1e-6, dist, MAX_GENETIC_DIST, 1e-6, negative_lh);
	//else // Brent method
		dist = minimizeOneDimen(1e-6, dist, MAX_GENETIC_DIST, 1e-6, &negative_lh, &ferror); 
	return dist;
}

AlignmentPairwise::~AlignmentPairwise()
{
	if (pair_freq) delete [] pair_freq;
}


