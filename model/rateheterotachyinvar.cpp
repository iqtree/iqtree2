/*
 * rateheterotachyinvar.cpp
 *
 *  Created on: Nov 11, 2016
 *      Author: minh
 */

#include "rateheterotachyinvar.h"

RateHeterotachyInvar::RateHeterotachyInvar(int ncat, PhyloTree* tree, 
                                           PhyloTree* report_to_tree)
: super(ncat, tree, report_to_tree), invar(0, tree) {
	cur_optimize = 0;
	name         = "+I"     + name;
	full_name    = "Invar+" + full_name;
    setNCategory(ncat);
}

RateHeterotachyInvar::RateHeterotachyInvar(int ncat, const string& params, 
                                           double p_invar_sites, PhyloTree* tree)
: super(ncat, params, tree), invar(p_invar_sites, tree)
{
	cur_optimize = 0;
	name         = "+I"     + name;
	full_name    = "Invar+" + full_name;
    setNCategory(ncat);
}

void RateHeterotachyInvar::startCheckpoint() {
    checkpoint->startStruct("RateHeterotachyInvar" + convertIntToString(ncategory));
}

void RateHeterotachyInvar::saveCheckpoint() {
	invar.setCheckpoint(checkpoint);
    invar.saveCheckpoint();
    super::saveCheckpoint();
}

void RateHeterotachyInvar::restoreCheckpoint() {
	invar.setCheckpoint(checkpoint);
    invar.restoreCheckpoint();
    super::restoreCheckpoint();
}

double RateHeterotachyInvar::getRate(int category) const {
    return 1.0;
    return 1.0 / (1.0 - invar.getPInvar());
}

void RateHeterotachyInvar::setNCategory(int ncat) {
	super::setNCategory(ncat);
	name = "+I" + name;
	full_name = "Invar+" + full_name;
}

double RateHeterotachyInvar::computeFunction(double value) {
	invar.setPInvar(value);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateHeterotachyInvar::targetFunk(double x[]) {
	return super::targetFunk(x);
}

void RateHeterotachyInvar::writeInfo(ostream &out) {
	invar.writeInfo(out);
	super::writeInfo(out);

}

/**
	write parameters, used with modeltest
	@param out output stream
*/
void RateHeterotachyInvar::writeParameters(ostream &out) {
	invar.writeParameters(out);
	super::writeParameters(out);
}

void RateHeterotachyInvar::setBounds(double *lower_bound, double *upper_bound,
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
double RateHeterotachyInvar::optimizeParameters(double gradient_epsilon,
                                                PhyloTree* report_to_tree) {
	double tree_lh;
	tree_lh = super::optimizeParameters(gradient_epsilon,
                                        report_to_tree);
	return tree_lh;
}

void RateHeterotachyInvar::setVariables(double *variables) {
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
bool RateHeterotachyInvar::getVariables(const double *variables) {
	bool changed = super::getVariables(variables);
	if (invar.getNDim() == 0) { 
		return changed;
	}
    changed |= (invar.getPInvar() != variables[getNDim()]);
	invar.setPInvar(variables[getNDim()]);
    return changed;
}

