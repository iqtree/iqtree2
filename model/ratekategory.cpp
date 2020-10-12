/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  BUI Quang Minh <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "tree/phylotree.h"
#include "ratekategory.h"

RateKategory::RateKategory(int ncat, PhyloTree *tree)
{
	ncategory = ncat;
	phylo_tree = tree;
	rates = new double[ncategory];
	name = "+K";
	name += convertIntToString(ncategory);
	full_name = "KAT";
	full_name += " with " + convertIntToString(ncategory) + " categories";
	if (ncategory == 1) { rates[0] = 1.0; return; } 
	int i;
	for (i = 0; i < ncategory; i++) do { rates[i] = random_double(); } while (rates[i]<0.1 || rates[i] > 0.9);
	//for (i = 0; i < ncategory; i++) rates[i] = 1.0 + i;
	double sum = 0.0;
	for (i = 0; i < ncategory; i++) sum += rates[i];
	for (i = 0; i < ncategory; i++) rates[i] = rates[i]*ncategory/sum;
}

RateKategory::~RateKategory()
{
	delete [] rates;
	rates = NULL;
}

double RateKategory::targetFunk(double x[])
{
	getVariables(x);
	if (rates[ncategory-1] < 1e-4) return 1.0e+12;
	ASSERT(phylo_tree);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateKategory::optimizeParameters(double gradient_epsilon)
{
	int ndim = getNDim();
	
	// return if nothing to be optimized
	if (ndim == 0) return 0.0;

	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters..." << endl;

	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	int i;
	double score;
	
	// by BFGS algorithm
	setVariables(variables);
	for (i = 1; i <= ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = 1e-4;
		upper_bound[i] = ncategory;
		bound_check[i] = false;
	}

	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, 1e-6));

	getVariables(variables);
	//sort(rates, rates+ncategory);
	phylo_tree->clearAllPartialLH();
    score = phylo_tree->computeLikelihood();
	
	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}

int RateKategory::computePatternRates(DoubleVector& pattern_rates, IntVector& pattern_cat)
{
	cout << "Computing site rates by empirical Bayes..." << endl;

	phylo_tree->computePatternLhCat(WSL_RATECAT);

	int npattern = phylo_tree->aln->getNPattern();
	pattern_rates.resize(npattern);
	pattern_cat.resize(npattern);

    double *lh_cat = phylo_tree->tree_buffers._pattern_lh_cat;
	for (int i = 0; i < npattern; i++) {
		double sum_rate = 0.0, sum_lh = 0.0;
		int best = 0;
		for (int c = 0; c < ncategory; c++) {
			sum_rate += rates[c] * lh_cat[c];
			sum_lh += lh_cat[c];
			if (lh_cat[c] > lh_cat[best]) best = c;
		}
		pattern_rates[i] = sum_rate / sum_lh;
		pattern_cat[i] = best;
		lh_cat += ncategory;
	}

    return ncategory;

//	int npattern = phylo_tree->aln->getNPattern();
//	double *ptn_rates = new double[npattern];
//	phylo_tree->computeLikelihoodBranchNaive(phylo_tree->getRoot()->firstNeighbor(),
//		phylo_tree->getRoot(), nullptr, ptn_rates);
//
//	pattern_rates.clear();
//	pattern_rates.insert(pattern_rates.begin(), ptn_rates, ptn_rates + npattern);
//	pattern_cat.resize(npattern, 0);
//	for (int i = 0; i < npattern; i++)
//		for (int j = 1; j < ncategory; j++)
//			if (fabs(rates[j] - ptn_rates[i]) < fabs(rates[pattern_cat[i]] - ptn_rates[i]))
//				pattern_cat[i] = j;
//	delete [] ptn_rates;
}

bool RateKategory::getVariables(double* variables)
{
	if (ncategory == 1) return false;
    bool changed = (rates[0] != 1.0);
	rates[0] = 1.0;
	changed |= memcmpcpy(rates, variables+1, (ncategory-1) * sizeof(double));
	double sum = 0.0;
	int i;
	for (i = 0; i < ncategory-1; i++) 
		sum += rates[i];
	/*
	for (i = 0; i < ncategory; i++) 
		rates[i] = rates[i]*ncategory/sum;*/
    changed |= (rates[ncategory-1] != ncategory - sum);
	rates[ncategory-1] = ncategory - sum;
    return changed;
}

void RateKategory::setVariables(double* variables)
{
	if (ncategory == 1) return;
	memcpy(variables+1, rates, (ncategory-1) * sizeof(double));
}


void RateKategory::writeInfo(ostream& out)
{
	out << "Rates: ";
	for (int i = 0; i < ncategory; i++)
		out << " " << rates[i];
	out << endl;
	out << "BIC: " << -2 * phylo_tree->computeLikelihood() + getNDim() * log(phylo_tree->getAlnNSite()) << endl;
}

void RateKategory::writeParameters(ostream& out)
{
}
