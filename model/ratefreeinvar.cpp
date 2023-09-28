/*
 * ratefreeinvar.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: minh
 */

#include "ratefreeinvar.h"

RateFreeInvar::RateFreeInvar(int ncat, double start_alpha, string params, bool sorted_rates, double p_invar_sites, string opt_alg, PhyloTree *tree)
: RateInvar(p_invar_sites, tree), RateFree(ncat, start_alpha, params, sorted_rates, opt_alg, tree)
{
	cur_optimize = 0;
	name = "+I" + name;
	full_name = "Invar+" + full_name;
    setNCategory(ncat);
}

void RateFreeInvar::startCheckpoint() {
    checkpoint->startStruct("RateFreeInvar" + convertIntToString(ncategory));
}

void RateFreeInvar::saveCheckpoint() {
    RateInvar::saveCheckpoint();
    RateFree::saveCheckpoint();
}

void RateFreeInvar::restoreCheckpoint() {
    RateInvar::restoreCheckpoint();
    RateFree::restoreCheckpoint();
}

/**
    return the number of dimensions
*/
int RateFreeInvar::getNDim() {
    if (optimizing_params == 1) {
        // only rates
        return RateFree::getNDim();
    } else {
        return RateInvar::getNDim() + RateFree::getNDim();
    }
}

void RateFreeInvar::setNCategory(int ncat) {
    /* NHANLT: bug fixed
     we should NOT re-initialize probs equally
     instead, we should renormalize props and rates according to invar_prop
     */
	// RateFree::setNCategory(ncat);
    double sum = 0;
    // normalize prop
    for (int i = 0; i < ncategory; i++) {
        prop[i] *= (1.0-getPInvar());
        sum += rates[i] * prop[i];
    }
    // normalize rates
    for (int i = 0; i < ncategory; i++) {
        rates[i] /= sum;
    }
    
    // the following two lines should be commented, as the same lines appear inside the constructor
	// name = "+I" + name;
	// full_name = "Invar+" + full_name;
}

double RateFreeInvar::computeFunction(double value) {
	p_invar = value;
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateFreeInvar::targetFunk(double x[]) {
    getVariables(x);
    phylo_tree->clearAllPartialLH();
    return -phylo_tree->computeLikelihood();
}

void RateFreeInvar::writeInfo(ostream &out) {
	RateInvar::writeInfo(out);
	RateFree::writeInfo(out);

}

/**
	write parameters, used with modeltest
	@param out output stream
*/
void RateFreeInvar::writeParameters(ostream &out) {
	RateInvar::writeParameters(out);
	RateFree::writeParameters(out);
}

void RateFreeInvar::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	RateFree::setBounds(lower_bound, upper_bound, bound_check);
	if (RateInvar::getNDim() == 0) return;
    if (optimizing_params != 1) {
        int ndim = getNDim()-1;
        RateInvar::setBounds(lower_bound+ndim, upper_bound+ndim, bound_check+ndim);
    }
}

/**
	optimize parameters
	@return the best likelihood
*/
double RateFreeInvar::optimizeParameters(double gradient_epsilon) {
	double tree_lh;
    //optimize_alg = "1-BFGS";
	tree_lh = RateFree::optimizeParameters(gradient_epsilon);
	return tree_lh;
}

void RateFreeInvar::setVariables(double *variables) {
    int n = getNDim();
    if (n == 0) return;
    int i;

    if (optimizing_params == 2) {
        // proportions
        for (i = 0; i < ncategory-1; i++)
            variables[i+1] = prop[i] / prop[ncategory-1];
        if (RateInvar::getNDim() > 0)
            variables[n] = getPInvar() / prop[ncategory-1];
    } else if (optimizing_params == 1) {
        // rates
        for (i = 0; i < ncategory-1; i++)
            variables[i+1] = rates[i] / rates[ncategory-1];
    } else {
        // both rates and weights
        for (i = 0; i < ncategory-1; i++) {
            variables[i+1] = prop[i] / prop[ncategory-1];
        }
        for (i = 0; i < ncategory-1; i++) {
            variables[i+ncategory] = rates[i] / rates[ncategory-1];
        }
        if (RateInvar::getNDim() > 0)
            variables[n] = getPInvar() / prop[ncategory-1];
    }
}

/**
	this function is served for the multi-dimension optimization. It should assign the model parameters
	from a vector of variables that is index from 1 (NOTE: not from 0)
	@param variables vector of variables, indexed from 1
*/
bool RateFreeInvar::getVariables(double *variables) {
    int n = getNDim();
    if (n == 0) return false;
    int i;
    bool changed = false;

    double sum = 1.0;
    if (optimizing_params == 2) {
        // proportions
        for (i = 0; i < ncategory-1; i++) {
            sum += variables[i+1];
        }
        if (RateInvar::getNDim() > 0)
            sum += variables[n];
        else {
            // the proportion of invariable sites are fixed
            sum /= (1.0 - getPInvar());
        }
        for (i = 0; i < ncategory-1; i++) {
            changed |= (prop[i] != variables[i+1] / sum);
            prop[i] = variables[i+1] / sum;
        }
        changed |= (prop[ncategory-1] != 1.0 / sum);
        prop[ncategory-1] = 1.0 / sum;
        if (RateInvar::getNDim() > 0) {
            changed |= (p_invar != variables[n] / sum);
            setPInvar(variables[n] / sum);
        }
        RateFree::rescaleRates();
    } else if (optimizing_params == 1) {
        // rates
        sum = prop[ncategory-1];
        for (i = 0; i < ncategory-1; i++) {
            sum += prop[i] * variables[i+1];
        }
        for (i = 0; i < ncategory-1; i++) {
            changed |= (rates[i] != variables[i+1] / sum);
            rates[i] = variables[i+1] / sum;
        }
        changed |= (rates[ncategory-1] != 1.0 / sum);
        rates[ncategory-1] = 1.0 / sum;
    } else {
        // both weights and rates
        for (i = 0; i < ncategory-1; i++) {
            sum += variables[i+1];
        }
        if (RateInvar::getNDim() > 0)
            sum += variables[n];
        else {
            // the proportion of invariable sites are fixed
            sum /= (1.0 - getPInvar());
        }
        for (i = 0; i < ncategory-1; i++) {
            changed |= (prop[i] != variables[i+1] / sum);
            prop[i] = variables[i+1] / sum;
        }
        changed |= (prop[ncategory-1] != 1.0 / sum);
        prop[ncategory-1] = 1.0 / sum;
        if (RateInvar::getNDim() > 0) {
            changed |= (p_invar != variables[n] / sum);
            setPInvar(variables[n] / sum);
        }

        // then rates
        sum = prop[ncategory-1];
        for (i = 0; i < ncategory-1; i++) {
            sum += prop[i] * variables[i+ncategory];
        }
        for (i = 0; i < ncategory-1; i++) {
            changed |= (rates[i] != variables[i+ncategory] / sum);
            rates[i] = variables[i+ncategory] / sum;
        }
        changed |= (rates[ncategory-1] != 1.0 / sum);
        rates[ncategory-1] = 1.0 / sum;
    }
    return changed;
}
