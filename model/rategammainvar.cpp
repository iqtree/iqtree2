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
#include "rategammainvar.h"

RateGammaInvar::RateGammaInvar(int ncat, double shape, bool median,
		double p_invar_sites, bool simultaneous, PhyloTree *tree) :
		RateInvar(p_invar_sites, tree), RateGamma(ncat, shape, median, tree) {
	name = "+I" + name;
	full_name = "Invar+" + full_name;
	joint_optimize = simultaneous;
    cur_optimize = 0;
	computeRates();
}

void RateGammaInvar::saveCheckpoint() {
    checkpoint->startStruct("RateGammaInvar");
//    CKP_SAVE(joint_optimize);
    checkpoint->endStruct();
    RateInvar::saveCheckpoint();
    RateGamma::saveCheckpoint();
}

void RateGammaInvar::restoreCheckpoint() {
    // should restore p_invar first before gamma, because RateGamma will call computeRates()
    RateInvar::restoreCheckpoint();
    RateGamma::restoreCheckpoint();
    checkpoint->startStruct("RateGammaInvar");
//    CKP_RESTORE(joint_optimize);
    checkpoint->endStruct();
}

void RateGammaInvar::setNCategory(int ncat) {
	RateGamma::setNCategory(ncat);
	name = "+I" + name;
	full_name = "Invar+" + full_name;
	computeRates();
}

string RateGammaInvar::getNameParams() {
	return RateInvar::getNameParams() + RateGamma::getNameParams();
}

double RateGammaInvar::computeFunction(double value) {
	if (cur_optimize == 0)
		gamma_shape = value;
	else
		p_invar = value;
	// need to compute rates again if p_inv or Gamma shape changes!
	computeRates();
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

void RateGammaInvar::writeInfo(ostream &out) {
	RateInvar::writeInfo(out);
	RateGamma::writeInfo(out);
}

void RateGammaInvar::writeParameters(ostream &out) {
	RateInvar::writeParameters(out);
	RateGamma::writeParameters(out);
}

void RateGammaInvar::setVariables(double *variables) {
	RateGamma::setVariables(variables);
	int gid = RateGamma::getNDim();
	RateInvar::setVariables(variables+gid);
}

bool RateGammaInvar::getVariables(double *variables) {
	int gid = RateGamma::getNDim();
	bool changed = RateGamma::getVariables(variables);
	changed |= RateInvar::getVariables(variables+gid);
    return changed;
}

double RateGammaInvar::targetFunk(double x[]) {
	assert(phylo_tree);
	getVariables(x);
	// need to compute rates again if p_inv or Gamma shape changes!
	RateGamma::computeRates();
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

void RateGammaInvar::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int gid = RateGamma::getNDim();
	RateGamma::setBounds(lower_bound, upper_bound, bound_check);
	RateInvar::setBounds(lower_bound+gid, upper_bound+gid, bound_check+gid);
}

double RateGammaInvar::optimizeParameters(double gradient_epsilon) {

	int ndim = getNDim();

	// return if nothing to be optimized
	if (ndim == 0)
		return phylo_tree->computeLikelihood();


	if (!joint_optimize) {
		double lh = phylo_tree->computeLikelihood();
		cur_optimize = 0;
		double gamma_lh;
		if (Params::getInstance().testAlpha) {
			gamma_lh = RateGamma::optimizeParameters(gradient_epsilon, 0.05, 10);
		} else {
            gamma_lh = RateGamma::optimizeParameters(gradient_epsilon);
        }
		assert(gamma_lh >= lh-0.1);
		cur_optimize = 1;
		double invar_lh = -DBL_MAX;
        invar_lh = RateInvar::optimizeParameters(gradient_epsilon);
		assert(invar_lh >= gamma_lh-0.1);
		//lh = tree_lh;

		//assert(gamma_lh >= invar_lh - 0.1);
//		phylo_tree->clearAllPartialLH();
//		return gamma_lh;
        cur_optimize = 0;
        return invar_lh;
	}

/*
	if (!joint_optimize) {
//		double lh = phylo_tree->computeLikelihood();
		cur_optimize = 1;
		double invar_lh = -DBL_MAX;
		invar_lh = RateInvar::optimizeParameters(gradient_epsilon);
//		assert(tree_lh >= lh-0.1);
//		lh = tree_lh;
		cur_optimize = 0;
		double gamma_lh;
		if (Params::getInstance().testAlpha) {
			gamma_lh = RateGamma::optimizeParameters(gradient_epsilon, 0.05, 10);
		} else {
			gamma_lh = RateGamma::optimizeParameters(gradient_epsilon);
		}
		//assert(gamma_lh >= invar_lh - 0.1);
		phylo_tree->clearAllPartialLH();
		return gamma_lh;
	}
*/

	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters by BFGS..." << endl;

	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	double score;

	// by BFGS algorithm
	setVariables(variables);
	setBounds(lower_bound, upper_bound, bound_check);

	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_GAMMA_SHAPE));

	getVariables(variables);

	phylo_tree->clearAllPartialLH();
    score = phylo_tree->computeLikelihood();

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}


int RateGammaInvar::computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat) {
	//cout << "Computing Gamma site rates by empirical Bayes..." << endl;

	phylo_tree->computePatternLhCat(WSL_RATECAT);

	int npattern = phylo_tree->aln->getNPattern();
	pattern_rates.resize(npattern);
	pattern_cat.resize(npattern);

	double *lh_cat = phylo_tree->_pattern_lh_cat;
	for (int i = 0; i < npattern; i++) {
		double sum_rate = 0.0, sum_lh = phylo_tree->ptn_invar[i];
		int best = 0;
        double best_lh = phylo_tree->ptn_invar[i];
		for (int c = 0; c < ncategory; c++) {
			sum_rate += rates[c] * lh_cat[c];
			sum_lh += lh_cat[c];
			if (lh_cat[c] > best_lh || (lh_cat[c] == best_lh && random_double()<0.5)) { // break tie at random
                best = c+1;
                best_lh = lh_cat[c];
            }
		}
		pattern_rates[i] = sum_rate / sum_lh;
		pattern_cat[i] = best;
		lh_cat += ncategory;
	}
    return ncategory+1;
}

