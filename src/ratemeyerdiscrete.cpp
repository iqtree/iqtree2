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
#include "ratemeyerdiscrete.h"
#include "kmeans/KMeans.h"

RateMeyerDiscrete::RateMeyerDiscrete(int ncat)
 : RateMeyerHaeseler()
{
	ncategory = ncat;
	rates = NULL;
	if (ncat > 0) {
		rates = new double[ncategory];
		memset(rates, 0, sizeof(double) * ncategory);
	}
	name = "+M";
	name += convertIntToString(ncategory);
	full_name += " with " + convertIntToString(ncategory) + " categories";
}


RateMeyerDiscrete::~RateMeyerDiscrete()
{
	if (rates) delete [] rates;
}

void RateMeyerDiscrete::optimizeRates() {
	RateMeyerHaeseler::optimizeRates();
}

void RateMeyerDiscrete::classifyRatesKMeans() {

	assert(ncategory > 0);
	int nsites = phylo_tree->aln->getNSite();

	// clustering the rates with k-means
	//AddKMeansLogging(&cout, false);
	double points[nsites];
	int assignments[nsites];
	int attempts = nsites, i;
	for (i = 0; i < nsites; i++) {
		if (at(phylo_tree->aln->getPatternID(i)) == MAX_SITE_RATE) 
			points[i] = log(1e6);
		else 
			points[i] = log(at(phylo_tree->aln->getPatternID(i)));
	}
	memset(rates, 0, sizeof(double) * ncategory);

	double cost = RunKMeansPlusPlus(nsites, ncategory, 1, points, attempts, rates, assignments);
	cout << "Rates are classified by k-means (" << attempts << " runs) with cost " << cost << endl;
	// assign the categorized rates
	double sum = 0.0, ok = 0.0;
	for (i = 0; i < ncategory; i++) rates[i] = exp(rates[i]);
	for (i = 0; i < nsites; i++) {
		double site_r = rates[assignments[i]];
		if (site_r < 2*MIN_SITE_RATE) site_r = MIN_SITE_RATE;
		if (site_r > MAX_SITE_RATE*0.99) site_r = MAX_SITE_RATE;
		at(phylo_tree->aln->getPatternID(i)) = site_r;
		if (site_r < MAX_SITE_RATE) { sum += site_r; ok += mean_rate; }
	}

	if (fabs(sum - ok) > 1e-3) {
		//cout << "Normalizing " << sum << " / " << ok << endl;
		double scale_f = ok / sum;
		for (i = 0; i < size(); i++) {
			if (at(i) > MIN_SITE_RATE && at(i) < MAX_SITE_RATE) at(i) = at(i) * scale_f;
	}

	}
	std::sort(rates, rates + ncategory);
	if (rates[0] < MIN_SITE_RATE) rates[0] = MIN_SITE_RATE;
	if (verbose_mode >= VB_MED) {
		//cout << "K-means cost: " << cost << endl;
		for (i = 0; i < ncategory; i++) cout << rates[i] << " ";
		cout << endl;
	}

}


double RateMeyerDiscrete::classifyRates(double tree_lh) {
	double new_tree_lh;
	if (ncategory > 0) {
		classifyRatesKMeans();
		phylo_tree->clearAllPartialLh();
		new_tree_lh = phylo_tree->computeLikelihood();
		return new_tree_lh;
	}

	// identifying proper number of categories
	int nptn = phylo_tree->aln->getNPattern();
	rates = new double[nptn];
	DoubleVector continuous_rates;
	getRates(continuous_rates);

	for (ncategory = 4; ; ncategory++) {
		classifyRatesKMeans();
		phylo_tree->clearAllPartialLh();
		new_tree_lh = phylo_tree->optimizeAllBranches();
		cout << "For " << ncategory << " categories, LogL = " << new_tree_lh << endl;
		if (new_tree_lh > tree_lh - 3.0) break;
		setRates(continuous_rates);
	}

	cout << endl << "Number of categories is set to " << ncategory << endl;
	return new_tree_lh;
}

