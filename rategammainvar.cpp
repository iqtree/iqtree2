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

RateGammaInvar::RateGammaInvar(int ncat, double shape, bool median, double p_invar_sites, bool simultaneous, PhyloTree *tree)
: RateInvar(p_invar_sites, tree), RateGamma(ncat, shape, median, tree)
{
	name = "+I" + name;
	full_name = "Invar+" + full_name;
	optimize_gamma_invar_by_bfgs = simultaneous;
}


double RateGammaInvar::computeFunction(double value) {
	if (cur_optimize == 0)
		gamma_shape = value;
	else 
		p_invar = value;
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
	int gid = RateGamma::getNDim();
	if (gid == 1) {
		variables[1] = gamma_shape;
	}
	if (RateInvar::getNDim() == 1) {
		variables[gid+1] = p_invar;
	}
}

void RateGammaInvar::getVariables(double *variables) {
	int gid = RateGamma::getNDim();
	if (gid == 1) {
		gamma_shape = variables[1];
	}
	if (RateInvar::getNDim() == 1) {
		p_invar = variables[gid+1];
	}
}

double RateGammaInvar::targetFunk(double x[]) {
	assert(phylo_tree);
	/*
	int gid = RateGamma::getNDim();
	if (gid == 1 && gamma_shape != x[1]) {
		gamma_shape = x[1];
		RateGamma::computeRates();
		phylo_tree->clearAllPartialLH();
	}
	if (RateInvar::getNDim() == 1) {
		p_invar = x[gid+1];
	}*/
	getVariables(x);
	RateGamma::computeRates();
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}


double RateGammaInvar::optimizeParameters() {


	if (!optimize_gamma_invar_by_bfgs) {
		double tree_lh;
		cur_optimize = 1;
		tree_lh = RateInvar::optimizeParameters();
		cur_optimize = 0;
		tree_lh = RateGamma::optimizeParameters();
		phylo_tree->clearAllPartialLH();
		return tree_lh;
	}

	int ndim = getNDim();

	// return if nothing to be optimized
	if (ndim == 0)
		return phylo_tree->computeLikelihood();

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
	int gid = RateGamma::getNDim();
	if (gid == 1) {
		lower_bound[1] = MIN_GAMMA_SHAPE;
		upper_bound[1] = MAX_GAMMA_SHAPE;
		bound_check[1] = false;
	}
	if (RateInvar::getNDim() == 1) {
		lower_bound[gid+1] = 1e-6;
		upper_bound[gid+1] = phylo_tree->aln->frac_const_sites;
		bound_check[gid+1] = false;
	}
		//packData(variables, lower_bound, upper_bound, bound_check);
	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, TOL_GAMMA_SHAPE);

	getVariables(variables);

	phylo_tree->clearAllPartialLH();

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



