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
#include "ratemeyerhaeseler.h"
#include "phylotree.h"

const double MIN_SITE_RATE = 1e-6;
const double MAX_SITE_RATE = 100.0;
const double TOL_SITE_RATE = 1e-6;


RateMeyerHaeseler::RateMeyerHaeseler()
 : RateHeterogeneity()
{
	name = "MvH";
	full_name = "Meyer & von Haeseler (2003)";
	dist_mat = NULL;
}


RateMeyerHaeseler::~RateMeyerHaeseler()
{
	if (dist_mat) delete [] dist_mat;
}


double RateMeyerHaeseler::optimizeParameters() {
	assert(phylo_tree);
	if (!dist_mat) {
		dist_mat = new double[phylo_tree->leafNum * phylo_tree->leafNum];
		// compute the distance based on the path lengths between taxa of the tree
		phylo_tree->calcDist(dist_mat);
	}
	if (empty()) resize(phylo_tree->aln->getNPattern(), 1.0);
	IntVector ok_ptn;
	ok_ptn.resize(size(), 0);
	double sum = 0.0;
	int i;
	int ok_sites = 0;
	int saturated_sites = 0;
	int invar_sites = 0;
	int ambiguous_sites = 0;
	int nseq = phylo_tree->leafNum;
	int nstates = phylo_tree->aln->num_states;
	for (i = 0; i < size(); i++) {
		int freq = phylo_tree->aln->at(i).frequency;
		if (phylo_tree->aln->at(i).computeAmbiguousChar(nstates) < nseq-2) {
			optimizeSiteRate(i);
			if (at(i) < MAX_SITE_RATE) {
				if(at(i) > MIN_SITE_RATE) {
					sum += at(i) * freq;
					ok_sites += freq;
					ok_ptn[i] = 1;
				} else  invar_sites += freq;
			} else  saturated_sites += freq; 
		} else { ambiguous_sites += freq; at(i) = 0; }
	} 

	// now scale such that the mean of rates is 1
	for (i = 0; i < size(); i++) {
		if (ok_ptn[i]) at(i) = at(i) * ok_sites / sum;
	}

	if (ambiguous_sites) {
		stringstream str;
		str << ambiguous_sites << " sites were ignored due to too many gaps or ambiguous characters";
		outWarning(str.str());
	}
	if (saturated_sites) {
		stringstream str;
		str << saturated_sites << " sites show saturated (too high) rates";
		outWarning(str.str());
	}
	cout << invar_sites << " sites have zero rate" << endl;
}

double RateMeyerHaeseler::optimizeSiteRate(int site) {
	optimizing_site = site;
	double negative_lh;
	double current_rate = at(site);
	double ferror, optx;
	optx = minimizeOneDimen(MIN_SITE_RATE, current_rate, MAX_SITE_RATE, TOL_SITE_RATE, &negative_lh, &ferror);
	if (optx > MAX_SITE_RATE*0.99) optx = MAX_SITE_RATE;
	if (optx < MIN_SITE_RATE*1.01) optx = 0;
	at(site) = optx;
	return optx;
}

double RateMeyerHaeseler::computeFunction(double value) {
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	int i, j;
	double lh = 0.0;
	SubstModel *model = phylo_tree->getModel();
	Pattern *pat = & phylo_tree->aln->at(optimizing_site);
	
	for (i = 0; i < nseq-1; i++)
		for (j = i+1; j < nseq; j++) {
			int state1 = pat->at(i);
			int state2 = pat->at(j);
			if (state1 >= nstate || state2 >= nstate) continue;
			lh += log(model->computeTrans(value * dist_mat[i*nseq + j], state1, state2));
		}
	return -lh;
}

void RateMeyerHaeseler::writeSiteRates(const char *file_name) {
	int nsite = phylo_tree->aln->getNSite();
	int i;
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name);
		for (i = 0; i < nsite; i++) 
			out << i+1 << "\t" << at(phylo_tree->aln->getPatternID(i)) << endl;
		out.close();
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}
}
