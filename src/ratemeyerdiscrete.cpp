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
	rates = new double[ncategory];
	memset(rates, 0, sizeof(double) * ncategory);
	name = "+M";
	name += convertIntToString(ncategory);
	full_name += " with " + convertIntToString(ncategory) + " categories";
}


RateMeyerDiscrete::~RateMeyerDiscrete()
{
	if (rates) delete [] rates;
}

void RateMeyerDiscrete::optimizeRates() {
	assert(ncategory > 0);

	RateMeyerHaeseler::optimizeRates();
	int nsites = phylo_tree->aln->getNSite();

	// clustering the rates with k-means
	//AddKMeansLogging(&cout, false);
	double points[nsites];
	int assignments[nsites];
	int attempts = 10, i;
	for (i = 0; i < nsites; i++) {
		points[i] = at(phylo_tree->aln->getPatternID(i));
		if (points[i] == MIN_SITE_RATE) points[i] = -MAX_SITE_RATE;
	}
	memset(rates, 0, sizeof(double) * ncategory);

	double cost = RunKMeansPlusPlus(nsites, ncategory, 1, points, attempts, rates, assignments);
	// assign the categorized rates
	double sum = 0.0, ok = 0.0;
	for (i = 0; i < nsites; i++) {
		double site_r = rates[assignments[i]];
		if (site_r < MIN_SITE_RATE) site_r = MIN_SITE_RATE;
		at(phylo_tree->aln->getPatternID(i)) = site_r;
		if (site_r < MAX_SITE_RATE) { sum += site_r; ok += mean_rate; }
	}

	if (fabs(sum - ok) > 1e-3) 
		cout << "Wrong normalization " << sum << " " << ok << endl;

	std::sort(rates, rates + ncategory);
	if (rates[0] < MIN_SITE_RATE) rates[0] = MIN_SITE_RATE;
	if (verbose_mode >= VB_MED) {
		cout << "K-means cost: " << cost << endl;
		for (i = 0; i < ncategory; i++) cout << rates[i] << " ";
		cout << endl;
	}

}
