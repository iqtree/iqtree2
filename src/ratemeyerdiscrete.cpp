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
#include "modeltest_wrapper.h"

/**
	Solve k-means problem for one-dimensional data with dynamic programming
	@param n number of data points
	@param ncat number of clusters
	@param data data point of size n: x[0..n-1]
	@param center (OUT) output k centers of k clusters: center[0...k-1] will be filled
	@param cluster (OUT) cluster assignments for each data point: cluster[0...n-1] will be filled
	@return the minimum sum of squares over all k clusters
*/
double kMeansOneDim(int n, int ncat, double *data, double *center, int *cluster) {
	int i, j, m, k = ncat;
	if (ncat == 0) k = n;
	/**
		dynamic programming cost matrix, c[i][j] = cost of i clusters for {x1...xj}
	*/
	double *c[k]; 
	/**
		id is used to trace back the minimal solution
	*/
	double *id[k]; 
	/** 
		c1[i][j] = cost of 1 cluster for {xi...xj}
	*/
	double *c1[n];
	/** 
		m1[i][j] = mean of {xi...xj}
	*/
	double *m1[n];
	
	double x[n]; // sorted data points

	double h[n]; // Hartigan index

	// allocate memory 
	for (i = 0; i < k; i++) c[i] = new double[n];
	for (i = 0; i < k; i++) id[i] = new double[n];
	for (i = 0; i < n; i++) c1[i] = new double[n];
	for (i = 0; i < n; i++) m1[i] = new double[n];

	// first sort data into x
	memmove(x, data, sizeof(double)*n);
	std::sort(x, x+n);
	// first compute c1 matrix
	for (i = 0; i < n; i++) {
		double sum = 0.0;
		for (j = i; j < n; j++) {
			sum += x[j];
			double mean = sum / (j-i+1);
			m1[i][j] = mean;
			double ss = 0; 
			for (m = i; m <= j; m++) 
				ss += (x[m]-mean)*(x[m]-mean); // sum of squared difference
				//ss += fabs(x[m]-mean); // sum of absolute difference
			c1[i][j] = ss;
		}
	}

	/* now compute dynamic programming matrix */
	// initialize the 1st row
	for (j = 0; j < n; j++) {
		c[0][j] = c1[0][j];
		id[0][j] = -1;
	}
	for (i = 1; i < k; i++) {
		// no i clusters exist for less than i data points
		for (j = 0; j < i; j++) { c[i][j] = INFINITY; id[i][j] = -1; }
		for (j = i; j < n; j++) {
			c[i][j] = INFINITY;
			for (m = i-1; m < j; m++)
				if (c[i][j] > c[i-1][m] + c1[m+1][j]) {
					c[i][j] = c[i-1][m] + c1[m+1][j];
					id[i][j] = m;
				}
		}
		// compute Hartigan index
		h[i-1] = (n-i-1)*(c[i-1][n-1]-c[i][n-1]) / c[i][n-1];
		//cout << i << " clusters " << h[i-1] << endl;
	}

	double min_cost = c[k-1][n-1];
	int bound[k+1];
	// now trace back
	bound[k] = n-1;
	for (i = k-1; i >= 0; i--) {
		bound[i] = id[i][bound[i+1]];
	}

	for (i = 0; i < k; i++) {
		center[i] = m1[bound[i]+1][bound[i+1]];
		for (j = 0; j < n; j++)
			if (data[j] <= x[bound[i+1]] && data[j] >= x[bound[i]+1])
				cluster[j] = i;
	}

	// free memory
	for (i = n-1; i >= 0; i--) delete [] m1[i];
	for (i = n-1; i >= 0; i--) delete [] c1[i];
	for (i = k-1; i >= 0; i--) delete [] id[i];
	for (i = k-1; i >= 0; i--) delete [] c[i];
	return min_cost;
}

/************************************************
	RateMeyerDiscrete
************************************************/
RateMeyerDiscrete::RateMeyerDiscrete(int ncat)
 : RateMeyerHaeseler()
{
	ncategory = ncat;
	rates = NULL;
	ptn_cat = NULL;
	is_categorized = false;
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

bool RateMeyerDiscrete::isSiteSpecificRate() { 
	return !is_categorized; 
}

int RateMeyerDiscrete::getNDiscreteRate() { 
	if (ncategory < 1) return 1;
	return ncategory; 
}

double RateMeyerDiscrete::getRate(int category) {
	if (isSiteSpecificRate()) return RateMeyerHaeseler::getRate(category);
	assert(category < ncategory); 
	return rates[category]; 
}

double RateMeyerDiscrete::getPtnCat(int ptn) {
	if (!ptn_cat) return -1;
	return ptn_cat[ptn];
}

void RateMeyerDiscrete::optimizeRates() {
	RateMeyerHaeseler::optimizeRates();
}

void RateMeyerDiscrete::classifyRatesKMeans() {

	assert(ncategory > 0);
	int nptn = size();

	// clustering the rates with k-means
	//AddKMeansLogging(&cout, false);
	double points[nptn];
	int i;
	if (!ptn_cat) ptn_cat = new int[nptn];
	for (i = 0; i < nptn; i++) {
		if (at(i) == MAX_SITE_RATE) 
			points[i] = log(1e6);
		else 
			points[i] = log(at(i));
	}
	memset(rates, 0, sizeof(double) * ncategory);

	//double cost = RunKMeansPlusPlus(nsites, ncategory, 1, points, attempts, rates, assignments);
	double cost = kMeansOneDim(nptn, ncategory, points, rates, ptn_cat);
	//cout << "Rates are classified by k-means (" << attempts << " runs) with cost " << cost << endl;
	//cout << "Minimum cost by dynamic programming is " << cost1 << endl;
	//if (cost > cost1+1e-6)
		//cout << "k-means++ is stuck in local minimum, global is " << cost1 << endl;
	// assign the categorized rates
	double sum = 0.0, ok = 0.0;
	for (i = 0; i < ncategory; i++) rates[i] = exp(rates[i]);
	if (rates[0] < MIN_SITE_RATE) rates[0] = MIN_SITE_RATE;
	if (rates[ncategory-1] > MAX_SITE_RATE) rates[ncategory-1] = MAX_SITE_RATE;
	if (verbose_mode >= VB_MED) {
		cout << "K-means cost: " << cost << endl;
		for (i = 0; i < ncategory; i++) cout << rates[i] << " ";
		cout << endl;
	}

	for (i = 0; i < nptn; i++) {
		at(i) = rates[ptn_cat[i]];
		if (at(i) < MAX_SITE_RATE) { 
			sum += at(i) * phylo_tree->aln->at(i).frequency; 
			ok += mean_rate * phylo_tree->aln->at(i).frequency; 
		}
	}

	if (fabs(sum - ok) > 1e-3) {
		//cout << "Normalizing " << sum << " / " << ok << endl;
		double scale_f = ok / sum;
		for (i = 0; i < size(); i++) {
			if (at(i) > MIN_SITE_RATE && at(i) < MAX_SITE_RATE) at(i) = at(i) * scale_f;
	}

	}

}


double RateMeyerDiscrete::classifyRates(double tree_lh) {
	double new_tree_lh;
	is_categorized = true;
	if (ncategory > 0) {
		cout << endl << "Classifying rates into " << ncategory << " categories..." << endl;
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
		cout << "For " << ncategory << " categories, LogL = " << new_tree_lh;
		double lh_diff = 2*(tree_lh - new_tree_lh);
		double df = (nptn - ncategory);
		double pval = computePValueChiSquare(lh_diff, df);
		cout << ", p-value = " << pval;
		cout << endl;
		//if (new_tree_lh > tree_lh - 3.0) break;
		if (pval > 0.05) break;
		setRates(continuous_rates);
	}

	cout << endl << "Number of categories is set to " << ncategory << endl;
	return new_tree_lh;
}



