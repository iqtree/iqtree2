/*
 * ratefree.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: minh
 */

#include "phylotree.h"
#include "ratefree.h"

RateFree::RateFree(int ncat, PhyloTree *tree) : RateHeterogeneity() {
	ncategory = ncat;
	phylo_tree = tree;

	rates = new double[ncategory];
	prop  = new double[ncategory];
	for (int i = 0; i < ncategory; i++) {
		prop[i] = 1.0 / ncategory;
		rates[i] = 1.0;
	}
	name = "+FR";
	name += convertIntToString(ncategory);
	full_name = "FreeRate";
	full_name += " with " + convertIntToString(ncategory) + " categories";

}

RateFree::~RateFree() {
	delete [] rates;
}

string RateFree::getNameParams() {
	stringstream str;
	str << "FR" << ncategory << "{";
	for (int i = 0; i < ncategory; i++) {
		if (i > 0) str << ",";
		str << prop[i]<< "," << rates[i];
	}
	str << "}";
	return str.str();
}


double RateFree::targetFunk(double x[]) {
	return 0.0;
}

/**
	optimize parameters. Default is to optimize gamma shape
	@return the best likelihood
*/
double RateFree::optimizeParameters(double epsilon) {
	return phylo_tree->computeLikelihood();
}


void RateFree::setVariables(double *variables) {
	if (getNDim() == 0) return;
	int i;
	for (i = 0; i < ncategory-1; i++)
		variables[i+1] = rates[i];
}

void RateFree::getVariables(double *variables) {
	if (getNDim() == 0) return;
}
