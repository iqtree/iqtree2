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
	joint_optimize = simultaneous;
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

void RateGammaInvar::getVariables(double *variables) {
	int gid = RateGamma::getNDim();
	RateGamma::getVariables(variables);
	RateInvar::getVariables(variables+gid);
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

double RateGammaInvar::optimizeParameters(double epsilon) {


	if (!joint_optimize) {
		double tree_lh;
		cur_optimize = 0;
		tree_lh = RateGamma::optimizeParameters(epsilon);
		cur_optimize = 1;
		tree_lh = RateInvar::optimizeParameters(epsilon);
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
	setBounds(lower_bound, upper_bound, bound_check);

	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(epsilon, TOL_GAMMA_SHAPE));

	getVariables(variables);

	phylo_tree->clearAllPartialLH();

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



