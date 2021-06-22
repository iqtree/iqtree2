/*
 * ratefreeinvar.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: minh
 */

#include "ratefreeinvar.h"

RateFreeInvar::RateFreeInvar(int ncat,
                             PhyloTree* tree, 
                             PhyloTree* report_to_tree) 
	: RateFree(ncat, tree, report_to_tree), invar(0.1, tree) {
}

RateFreeInvar::RateFreeInvar(int ncat, double start_alpha, string params,
                             bool sorted_rates, double p_invar_sites,
                             string opt_alg, PhyloTree *tree)
    : RateFree(ncat, start_alpha, params, sorted_rates, opt_alg, tree)
	, invar(p_invar_sites, tree)
{
	cur_optimize = 0;
	name         = "+I" + name;
	full_name    = "Invar+" + full_name;
    setNCategory(ncat);
}

void RateFreeInvar::startCheckpoint() {
    checkpoint->startStruct("RateFreeInvar" + convertIntToString(ncategory));
}

void RateFreeInvar::saveCheckpoint() {
    invar.saveCheckpoint();
    super::saveCheckpoint();
}

void RateFreeInvar::restoreCheckpoint() {
    invar.restoreCheckpoint();
    super::restoreCheckpoint();
}

int RateFreeInvar::getNDim() const { 
	return invar.getNDim() + super::getNDim(); 
}

double RateFreeInvar::getProp(int category) const {
	return prop[category]; 
}

double RateFreeInvar::getRate(int category) const {
	return super::getRate(category); 
}

std::string RateFreeInvar::getNameParams() const {
	return invar.getNameParams() + super::getNameParams();
}


void RateFreeInvar::setNCategory(int ncat) {
	super::setNCategory(ncat);
	name = "+I" + name;
	full_name = "Invar+" + full_name;
}

double RateFreeInvar::computeFunction(double value) {
	invar.setPInvar(value);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateFreeInvar::targetFunk(double x[]) {
	return super::targetFunk(x);
}

void RateFreeInvar::writeInfo(ostream &out) {
	invar.writeInfo(out);
	super::writeInfo(out);

}

/**
	write parameters, used with modeltest
	@param out output stream
*/
void RateFreeInvar::writeParameters(ostream &out) {
	invar.writeParameters(out);
	super::writeParameters(out);
}

void RateFreeInvar::setBounds(double *lower_bound, double *upper_bound,
                              bool *bound_check) {
	super::setBounds(lower_bound, upper_bound, bound_check);
	if (invar.getNDim() == 0) {
		return;
	}
	int ndim = getNDim()-1;
	invar.setBounds(lower_bound+ndim, upper_bound+ndim, bound_check+ndim);
}

/**
	optimize parameters
	@return the best likelihood
*/
double RateFreeInvar::optimizeParameters(double gradient_epsilon,
                                         PhyloTree* report_to_tree) {
	double tree_lh;
	tree_lh = super::optimizeParameters(gradient_epsilon,
                                           report_to_tree);
	return tree_lh;
}

void RateFreeInvar::setVariables(double *variables) {
	super::setVariables(variables);
	if (invar.getNDim() == 0) {
		return;
	}
	variables[getNDim()] = invar.getPInvar();
}

/**
	this function is served for the multi-dimension optimization. It should assign the model parameters
	from a vector of variables that is index from 1 (NOTE: not from 0)
	@param variables vector of variables, indexed from 1
*/
bool RateFreeInvar::getVariables(double *variables) {
	bool changed = super::getVariables(variables);
	if (invar.getNDim() == 0) {
		return changed;
	}
    changed |= (invar.getPInvar() != variables[getNDim()]);
	invar.setPInvar(variables[getNDim()]);
    return changed;
}

