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

AlignmentPairwise::AlignmentPairwise(Alignment *aln, int seq_id1, int seq_id2) : Alignment(), Optimization() {
	num_states = aln->num_states;
	pair_freq = new int[num_states * num_states];
	memset(pair_freq, 0, sizeof(int) * num_states * num_states);
	for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
		int state1 = (*it)[seq_id1];
		int state2 = (*it)[seq_id2];
		if (state1 < num_states && state2 < num_states)
			pair_freq[state1 * num_states + state2] += it->frequency;
	}
	model_factory = NULL;
	site_rate = NULL;
}

double AlignmentPairwise::computeFunction(double value) {

	int ncat = site_rate->getNRate();
	int trans_size = num_states * num_states;
	int cat, i;

	double trans_mat[trans_size];
	double sum_trans_mat[trans_size];

	for (cat = 0; cat < ncat; cat++) {
		model_factory->computeTransMatrix(value * site_rate->getRate(cat), trans_mat);
		if (cat == 0) 
			memcpy(sum_trans_mat, trans_mat, sizeof(double)*trans_size);
		else
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

double AlignmentPairwise::optimizeDist(double initial_dist) {
	// initial guess of the distance using Juke-Cantor correction
	double dist = initial_dist; 

	// if no model or rate is specified, return the JC distance
	if (!model_factory || !site_rate) return dist;

	double negative_lh, ferror;
	dist = minimizeOneDimen(1e-6, dist, MAX_GENETIC_DIST, 1e-6, &negative_lh, &ferror); 
	return dist;
}

AlignmentPairwise::~AlignmentPairwise()
{
	if (pair_freq) delete [] pair_freq;
}


