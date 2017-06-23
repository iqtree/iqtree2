/*
 * rateheterotachyinvar.cpp
 *
 *  Created on: Nov 11, 2016
 *      Author: minh
 */

#include "rateheterotachyinvar.h"

RateHeterotachyInvar::RateHeterotachyInvar(int ncat, string params, double p_invar_sites, PhyloTree *tree)
: RateInvar(p_invar_sites, tree), RateHeterotachy(ncat, params, tree)
{
	cur_optimize = 0;
	name = "+I" + name;
	full_name = "Invar+" + full_name;
    setNCategory(ncat);
}

void RateHeterotachyInvar::startCheckpoint() {
    checkpoint->startStruct("RateHeterotachyInvar" + convertIntToString(ncategory));
}

void RateHeterotachyInvar::saveCheckpoint() {
    RateInvar::saveCheckpoint();
    RateHeterotachy::saveCheckpoint();
}

void RateHeterotachyInvar::restoreCheckpoint() {
    RateInvar::restoreCheckpoint();
    RateHeterotachy::restoreCheckpoint();
}

double RateHeterotachyInvar::getRate(int category) {
    return 1.0;
    return 1.0 / (1.0 - p_invar);
}

void RateHeterotachyInvar::setNCategory(int ncat) {
	RateHeterotachy::setNCategory(ncat);
	name = "+I" + name;
	full_name = "Invar+" + full_name;
}

double RateHeterotachyInvar::computeFunction(double value) {
	p_invar = value;
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateHeterotachyInvar::targetFunk(double x[]) {
	return RateHeterotachy::targetFunk(x);
}

void RateHeterotachyInvar::writeInfo(ostream &out) {
	RateInvar::writeInfo(out);
	RateHeterotachy::writeInfo(out);

}

/**
	write parameters, used with modeltest
	@param out output stream
*/
void RateHeterotachyInvar::writeParameters(ostream &out) {
	RateInvar::writeParameters(out);
	RateHeterotachy::writeParameters(out);
}

void RateHeterotachyInvar::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	RateHeterotachy::setBounds(lower_bound, upper_bound, bound_check);
	if (RateInvar::getNDim() == 0) return;
	int ndim = getNDim()-1;
	RateInvar::setBounds(lower_bound+ndim, upper_bound+ndim, bound_check+ndim);
}

/**
	optimize parameters
	@return the best likelihood
*/
double RateHeterotachyInvar::optimizeParameters(double gradient_epsilon) {
	double tree_lh;
	tree_lh = RateHeterotachy::optimizeParameters(gradient_epsilon);
	return tree_lh;
}

void RateHeterotachyInvar::setVariables(double *variables) {
	RateHeterotachy::setVariables(variables);
	if (RateInvar::getNDim() == 0) return;
	variables[getNDim()] = p_invar;
}

/**
	this function is served for the multi-dimension optimization. It should assign the model parameters
	from a vector of variables that is index from 1 (NOTE: not from 0)
	@param variables vector of variables, indexed from 1
*/
bool RateHeterotachyInvar::getVariables(double *variables) {
	bool changed = RateHeterotachy::getVariables(variables);
	if (RateInvar::getNDim() == 0) return changed;
    changed |= (p_invar != variables[getNDim()]);
	p_invar = variables[getNDim()];
    return changed;
}

