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



RateMeyerHaeseler::RateMeyerHaeseler(double mean_rate_val)
 : RateHeterogeneity()
{
	name = "+MH";
	full_name = "Meyer & von Haeseler (2003)";
	dist_mat = NULL;
	mean_rate = mean_rate_val;
}


RateMeyerHaeseler::~RateMeyerHaeseler()
{
	if (dist_mat) delete [] dist_mat;
}

int RateMeyerHaeseler::getNDim() {
	if (phylo_tree) 
		return phylo_tree->aln->getNPattern()-1;
	if (empty()) return 0;
	return size()-1; 
}

/*
double RateMeyerHaeseler::getRate(int category) {
	if (category < size())
		return at(category);

	return 1.0;
}*/

double RateMeyerHaeseler::getPtnRate(int ptn) {
	if (ptn < size())
		return at(ptn);

	return 1.0;
}

void RateMeyerHaeseler::getRates(DoubleVector &rates) {
	rates.clear();
	if (empty()) {
		rates.resize(phylo_tree->aln->size(), mean_rate);
	} else {
		rates.insert(rates.begin(), begin(), end());
	} 
}

void RateMeyerHaeseler::setRates(DoubleVector &rates) {
	clear();
	insert(begin(), rates.begin(), rates.end());
}

void RateMeyerHaeseler::initializeRates() {

	int i, j, rate_id = 0, state1, state2;
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;

	if (nseq < 25) 
		outWarning("Meyer & von Haeseler model is not recommended for < 25 sequences\n");

	resize(phylo_tree->aln->getNPattern(), mean_rate);

	for (Alignment::iterator pat = phylo_tree->aln->begin(); pat != phylo_tree->aln->end(); pat++, rate_id++) {
		double diff = 0, total = 0;
		for (i = 0; i < nseq-1; i++) if ((state1 = pat->at(i)) < nstate)
			for (j = i+1; j < nseq; j++) if ((state2 = pat->at(j)) < nstate) {
				//total += dist_mat[state1 * nstate + state2];
				//if (state1 != state2) diff += dist_mat[state1 * nstate + state2];
				total++;
				if (state1 != state2) diff++;
		}
		double obs_diff = double(diff) / total;
		double tolog = 1.0 - obs_diff*nstate/(nstate-1);
		if (tolog > 0.0) {
			at(rate_id) = -log(tolog) * (nstate-1) / nstate;
		} else at(rate_id) = obs_diff;
		
	}
}


double RateMeyerHaeseler::optimizeRate(int pattern) {
	optimizing_pattern = pattern;

	double max_rate = MAX_SITE_RATE;

	double minf = INFINITY, minx = 0;
	double negative_lh;
	double current_rate = at(pattern);
	double ferror, optx;
	/* constant site alway have ZERO rates */
	if (phylo_tree->aln->at(pattern).is_const) {
		return (at(pattern) = MIN_SITE_RATE);
	}
    if (phylo_tree->optimize_by_newton) // Newton-Raphson method 
	{
    	optx = minimizeNewtonSafeMode(MIN_SITE_RATE, current_rate, max_rate, TOL_SITE_RATE, negative_lh);
		if (optx > MAX_SITE_RATE*0.99 || (optx < MIN_SITE_RATE*2 && !phylo_tree->aln->at(pattern).is_const)) 
		{
			double optx2, negative_lh2;
			optx2 = minimizeOneDimen(MIN_SITE_RATE, current_rate, max_rate, TOL_SITE_RATE, &negative_lh2, &ferror);
			if (negative_lh2 < negative_lh - 1e-4) {
				cout << "+++NEWTON IS WRONG for pattern " << pattern << ": " << optx2 << " " << 
				negative_lh2 << " (Newton: " << optx << " " << negative_lh <<")" << endl;
			}
			if (negative_lh < negative_lh2 - 1e-4 && verbose_mode >= VB_MED) {
				cout << "Brent is wrong for pattern " << pattern << ": " << optx2 << " " << 
				negative_lh2 << " (Newton: " << optx << " " << negative_lh <<")" << endl;
			}
		}
    }
    else {
		optx = minimizeOneDimen(MIN_SITE_RATE, current_rate, max_rate, TOL_SITE_RATE, &negative_lh, &ferror);
		double fnew;
		if ((optx < max_rate) && (fnew = computeFunction(max_rate)) <= negative_lh+TOL_SITE_RATE) {
			optx = max_rate;
			negative_lh = fnew;
		}
		if ((optx > MIN_SITE_RATE) && (fnew = computeFunction(MIN_SITE_RATE)) <= negative_lh+TOL_SITE_RATE) {
			optx = MIN_SITE_RATE;
			negative_lh = fnew;
		}
	}
	//negative_lh = brent(MIN_SITE_RATE, current_rate, max_rate, 1e-3, &optx);
	if (optx > max_rate*0.99) optx = MAX_SITE_RATE;
	if (optx < MIN_SITE_RATE*2) optx = MIN_SITE_RATE;
	at(pattern) = optx;
//#ifndef NDEBUG		
	if (optx == MAX_SITE_RATE || (optx == MIN_SITE_RATE && !phylo_tree->aln->at(pattern).is_const)) {
		ofstream out;
	
		if (verbose_mode >= VB_MED)  {
			cout << "Checking pattern " << pattern << " (" << current_rate << ", " << optx << ")" << endl;
			out.open("x", ios::app);
			out << pattern;
		}
		for (double val=0.1; val <= 100; val += 0.1) {
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
			at(pattern) = optx;
			if (verbose_mode >= VB_MED)
				cout << "FIX rate: " << minx << " , " << optx << endl;
		}
	}
//#endif

	return optx;
}


void RateMeyerHaeseler::optimizeRates() {
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
			optimizeRate(i);
			if (at(i) == MIN_SITE_RATE) invar_sites += freq; 
			if (at(i) == MAX_SITE_RATE) {
				saturated_sites += freq; 
				saturated_ptn ++;
			}
		} else { at(i) = 1.0; ambiguous_sites += freq; }
		if (at(i) < MAX_SITE_RATE) 
		{
			if (at(i) > MIN_SITE_RATE) sum += at(i) * freq;
			ok_ptn[i] = 1;
			ok_sites += freq;
		}
	} 

	// now scale such that the mean of rates is 1
	double scale_f = ok_sites * mean_rate / sum;
	if (mean_rate > 0.0)
	for (i = 0; i < size(); i++) {
		if (ok_ptn[i] && at(i) > MIN_SITE_RATE) at(i) = at(i) * scale_f;
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

}

double RateMeyerHaeseler::optimizeParameters() {
	assert(phylo_tree);
	double tree_lh = phylo_tree->computeLikelihood();
	//double tree_lh = phylo_tree->optimizeAllBranches();
	DoubleVector prev_rates;
	getRates(prev_rates);

	optimizeRates();

// DEBUG
/*
	writeSiteRates(*this, "xrate");
	double pattern_lh[size()];
	double new_tree_lh1 = phylo_tree->computeLikelihood(pattern_lh);
	ofstream out("xlogl");
	double check = 0.0;
	for (i = 0; i < size(); i++) { 
		out << pattern_lh[i] << endl;
		check += pattern_lh[i] * phylo_tree->aln->at(i).frequency;
	}
	out.close();
	cout << "Site rates and pattern log-likelihood written " << new_tree_lh1 << " " << check << endl;
	exit(1);*/
// END DEBUG
/*
	double new_tree_lh = phylo_tree->computeLikelihood();
	if (new_tree_lh < tree_lh - TOL_LIKELIHOOD) {
		cout << "Worse likelihood (" << new_tree_lh << "), rolling back site rates..." << endl;
		setRates(prev_rates);
		phylo_tree->clearAllPartialLh();
		new_tree_lh = phylo_tree->computeLikelihood();
		cout << "Back up likelihood: " << new_tree_lh << endl;
	}
*/
	
	phylo_tree->clearAllPartialLh();
	stringstream best_tree_string;
	phylo_tree->printTree(best_tree_string, WT_BR_LEN + WT_TAXON_ID);
	double new_tree_lh = phylo_tree->optimizeAllBranches(1);
	if (new_tree_lh < tree_lh + 1e-5) {
		cout << "Worse likelihood (" << new_tree_lh << "), roll back site rates..." << endl;
		setRates(prev_rates);
		phylo_tree->rollBack(best_tree_string);
		//phylo_tree->clearAllPartialLh();
		new_tree_lh = phylo_tree->computeLikelihood();
		cout << "Backup log-likelihood: " << new_tree_lh << endl;
		return classifyRates(new_tree_lh);
	}
	
	return new_tree_lh;
}


double RateMeyerHaeseler::computeFunction(double value) {
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	int i, j, state1, state2;
	double lh = 0.0;
	SubstModel *model = phylo_tree->getModel();
	Pattern *pat = & phylo_tree->aln->at(optimizing_pattern);
	
	for (i = 0; i < nseq-1; i++) if ((state1 = pat->at(i)) < nstate) 
		for (j = i+1; j < nseq; j++) if ((state2 = pat->at(j)) < nstate) 
			lh -= log(model->computeTrans(value * dist_mat[i*nseq + j], state1, state2));
	return lh;
}

double RateMeyerHaeseler::computeFuncDerv(double value, double &df, double &ddf) {
	int nseq = phylo_tree->leafNum;
	int nstate = phylo_tree->getModel()->num_states;
	int i, j, state1, state2;
	double lh = 0.0;
	double trans, derv1, derv2;
	SubstModel *model = phylo_tree->getModel();
	Pattern *pat = & phylo_tree->aln->at(optimizing_pattern);
	df = ddf = 0.0;
	for (i = 0; i < nseq-1; i++) if ((state1 = pat->at(i)) < nstate) 
		for (j = i+1; j < nseq; j++) if ((state2 = pat->at(j)) < nstate) {
			double dist = dist_mat[i*nseq + j];
			trans = model->computeTrans(value * dist, state1, state2, derv1, derv2);
			lh -= log(trans);
			double t1 = derv1 / trans;
			double t2 = derv2 / trans;
			df -= t1 * dist;
			ddf -= dist * dist * (t2 - t1*t1);
		}
	return lh;
}


void RateMeyerHaeseler::runIterativeProc(Params &params, IQPTree &tree) {
	int i;
	if (verbose_mode >= VB_MED) {
		ofstream out("x");
		out.close();
	}
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

	//setRates(prev_rates);
	string rate_file = params.out_prefix;
	rate_file += ".mhrate";
	double prev_lh = tree.getBestScore();
	string dist_file = params.out_prefix;
	dist_file += ".tdist";
	tree.getModelFactory()->stopStoringTransMatrix();

	for (i = 2; i < 100; i++) {
		DoubleVector prev_rates;
		getRates(prev_rates);
		writeSiteRates(prev_rates, rate_file.c_str());
		tree.curScore = optimizeParameters();
		//phylo_tree->aln->printDist(dist_file.c_str(), dist_mat);
		tree.curScore = tree.optimizeAllBranches(i);
		cout << "Current Log-likelihood: " << tree.curScore << endl;
		if (tree.curScore <= prev_lh + 1e-4) {
			tree.setBestScore(tree.curScore);
			break;
		}
		prev_lh = tree.curScore;
		tree.setBestScore(tree.curScore);
	}
	cout << "Optimization took " << i-1 << " rounds to finish" << endl;
	tree.getModelFactory()->startStoringTransMatrix();
	//tree.getModelFactory()->site_rate = backup_rate;
	//tree.setRate(backup_rate);
}