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
#include "ratemeyerhaeseler.h"



RateMeyerHaeseler::RateMeyerHaeseler()
 : RateHeterogeneity()
{
	name = "MH";
	full_name = "Meyer & von Haeseler (2003)";
	dist_mat = NULL;
}


RateMeyerHaeseler::~RateMeyerHaeseler()
{
	if (dist_mat) delete [] dist_mat;
}


double RateMeyerHaeseler::getRate(int category) {
	if (category < size())
		return at(category);

	return 1.0;
}

void RateMeyerHaeseler::getRates(DoubleVector &rates) {
	rates.clear();
	if (empty()) {
		rates.resize(phylo_tree->aln->size(), 1.0);
	} else {
		rates.insert(rates.begin(), begin(), end());
	} 
}

void RateMeyerHaeseler::setRates(DoubleVector &rates) {
	clear();
	insert(begin(), rates.begin(), rates.end());
}

void RateMeyerHaeseler::initializeRates() {

	int i, j, rate_id = 0;
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;

	resize(phylo_tree->aln->getNPattern(), 1.0);

	for (Alignment::iterator pat = phylo_tree->aln->begin(); pat != phylo_tree->aln->end(); pat++, rate_id++) {
		double diff = 0, total = 0;
		for (i = 0; i < nseq-1; i++)
			for (j = i+1; j < nseq; j++) {
				int state1 = pat->at(i);
				int state2 = pat->at(j);
				if (state1 >= nstate || state2 >= nstate) continue;
				total += dist_mat[state1 * nstate + state2];
				if (state1 != state2) diff += dist_mat[state1 * nstate + state2];
		}
		double obs_diff = double(diff) / total;
		double tolog = 1.0 - obs_diff*nstate/(nstate-1);
		if (tolog > 0.0) {
			at(rate_id) = -log(tolog) * (nstate-1) / nstate;
		} else at(rate_id) = 10.0;
		
	}
}

double RateMeyerHaeseler::optimizeParameters() {
	assert(phylo_tree);
	if (!dist_mat) {
		dist_mat = new double[phylo_tree->leafNum * phylo_tree->leafNum];
	}
	// compute the distance based on the path lengths between taxa of the tree
	phylo_tree->calcDist(dist_mat);
	if (empty()) 
		initializeRates();
	IntVector ok_ptn;
	ok_ptn.resize(size(), 0);
	double sum = 0.0;
	int i;
	int ok_sites = 0;
	int saturated_sites = 0, saturated_ptn = 0;
	int invar_sites = 0;
	int ambiguous_sites = 0;
	int nseq = phylo_tree->leafNum;
	int nstates = phylo_tree->aln->num_states;
	for (i = 0; i < size(); i++) {
		int freq = phylo_tree->aln->at(i).frequency;
		if (phylo_tree->aln->at(i).computeAmbiguousChar(nstates) < nseq-2) {
			optimizeSiteRate(i);
			if (at(i) == 0.0) invar_sites += freq; 
			if (at(i) == MAX_SITE_RATE) {
				saturated_sites += freq; 
				saturated_ptn ++;
			}
		} else { ambiguous_sites += freq; }
		if (at(i) < MAX_SITE_RATE) {
			sum += at(i) * freq;
			ok_ptn[i] = 1;
			ok_sites += freq;
		}
	} 

	// now scale such that the mean of rates is 1
	for (i = 0; i < size(); i++) {
		if (ok_ptn[i]) at(i) = at(i) * ok_sites / sum;
	}

	if (ambiguous_sites) {
		stringstream str;
		str << ambiguous_sites << " sites contain too many gaps or ambiguous characters";
		outWarning(str.str());
	}
	if (saturated_sites) {
		stringstream str;
		str << saturated_sites << " sites (" << saturated_ptn << " patterns) show too high rates (>=" << MAX_SITE_RATE << ")";
		outWarning(str.str());
	}
	cout << invar_sites << " sites have zero rate" << endl;
	phylo_tree->clearAllPartialLh();
	return phylo_tree->computeLikelihood();
}

double RateMeyerHaeseler::optimizeSiteRate(int site) {
	optimizing_site = site;

	double max_rate = MAX_SITE_RATE;

		double minf = INFINITY, minx;
	double negative_lh;
	double current_rate = at(site);
	double ferror, optx, optx2;
    if (phylo_tree->optimize_by_newton) // Newton-Raphson method 
	{
    	optx = minimizeNewton(MIN_SITE_RATE, current_rate, max_rate, TOL_SITE_RATE, negative_lh);
    }
    else 
		optx = minimizeOneDimen(MIN_SITE_RATE, current_rate, max_rate, TOL_SITE_RATE, &negative_lh, &ferror);
	//negative_lh = brent(MIN_SITE_RATE, current_rate, max_rate, 1e-3, &optx);
	if (optx > max_rate*0.99) optx = MAX_SITE_RATE;
	if (optx < MIN_SITE_RATE*1.01) optx = 0;
	at(site) = optx;
//#ifndef NDEBUG		
	if (optx == MAX_SITE_RATE) {
		ofstream out;
	
		if (verbose_mode >= VB_MED)  {
			cout << "Checking pattern " << site << " (" << current_rate << ")" << endl;
			out.open("x", ios::app);
			out << site;
		}
		for (double val=0.1; val <= 30; val += 0.1) {
			double f = computeFunction(val);
			
			if (verbose_mode >= VB_MED) out << " " << f;
			if (f < minf) { minf = f; minx = val; }
			if (verbose_mode < VB_MED && minf < negative_lh) break;
		}
		if (verbose_mode >= VB_MED) { 
			out << endl;
			out.close();
		}
		//cout << "minx: " << minx << " " << minf << endl;
		if (negative_lh > minf+1e-3) {
			optx = minimizeOneDimen(MIN_SITE_RATE, minx, max_rate, 1e-3, &negative_lh, &ferror);
			at(site) = optx;
			if (verbose_mode >= VB_MED)
				cout << "FIX rate: " << minx << " , " << optx << endl;
		}
	}
//#endif

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

double RateMeyerHaeseler::computeFuncDerv(double value, double &df, double &ddf) {
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	int i, j;
	double lh = 0.0;
	double trans, derv1, derv2;
	SubstModel *model = phylo_tree->getModel();
	Pattern *pat = & phylo_tree->aln->at(optimizing_site);
	df = ddf = 0.0;
	for (i = 0; i < nseq-1; i++)
		for (j = i+1; j < nseq; j++) {
			int state1 = pat->at(i);
			int state2 = pat->at(j);
			if (state1 >= nstate || state2 >= nstate) continue;
			trans = model->computeTrans(value * dist_mat[i*nseq + j], state1, state2, derv1, derv2);
			lh += log(trans);
			double t = derv1 / trans;
			df -= t;
			ddf -= (derv2/trans - t*t);
		}
	return -lh;
}


void RateMeyerHaeseler::runIterativeProc(Params &params, IQPTree &tree) {
	int i;
	if (tree.leafNum < 25) 
		cout << "WARNING: Method not recommended for less than 25 sequences" << endl;
	ofstream out("x");
	out.close();
	setTree(&tree);
	RateHeterogeneity *backup_rate = tree.getRate();
	if (backup_rate->getGammaShape() > 0 ) {
		backup_rate->computePatternRates(*this);
		double sum = 0.0;
		for (i = 0; i < size(); i++)
			sum += at(i) * phylo_tree->aln->at(i).frequency;
		sum /=  phylo_tree->aln->getNSite();
		if (fabs(sum - 1.0) > 0.0001) {
			if (verbose_mode >= VB_MED)
				cout << "Normalizing Gamma rates (" << sum << ")" << endl;
			for (i = 0; i < size(); i++)
				at(i) /= sum;
		}
	}
	tree.getModelFactory()->site_rate = this;
	tree.setRate(this);

	
	//if  (empty()) initializeRates();

	DoubleVector prev_rates;
	getRates(prev_rates);
	//setRates(prev_rates);
	string rate_file = params.out_prefix;
	rate_file += ".mhrate";
	double prev_lh = tree.getBestScore();
	string dist_file = params.out_prefix;
	dist_file += ".tdist";

	for (i = 2; i < 100; i++) {
		optimizeParameters();
		writeSiteRates(*this, rate_file.c_str());
		phylo_tree->aln->printDist(dist_file.c_str(), dist_mat);
		tree.curScore = tree.optimizeAllBranches();
		if (params.min_iterations) tree.curScore = tree.optimizeNNI();
		cout << "Current Log-likelihood: " << tree.curScore << endl;
		if (tree.curScore <= prev_lh + TOL_LIKELIHOOD) {
			setRates(prev_rates);
			tree.curScore = tree.optimizeAllBranches();
			tree.setBestScore(tree.curScore);
			break;
		}
		getRates(prev_rates);
		prev_lh = tree.curScore;
		tree.setBestScore(tree.curScore);
	}
	cout << "Optimization took " << i-1 << " rounds to finish" << endl;
	//tree.getModelFactory()->site_rate = backup_rate;
	//tree.setRate(backup_rate);
}