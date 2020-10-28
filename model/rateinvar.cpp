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
#include "rateinvar.h"

RateInvar::RateInvar(double p_invar_sites, PhyloTree *tree)
 : RateHeterogeneity()
{
	if (tree) {
        if (tree->aln->frac_const_sites == 0.0)
            p_invar = 0.0;
        else
            p_invar = max(tree->aln->frac_const_sites/2.0, MIN_PINVAR);
//		p_invar = MIN_PINVAR;
	} else
		p_invar = MIN_PINVAR;
	fix_p_invar = false;
//    optimize_p_invar = true;
	phylo_tree = tree;
	name = "+I";
	full_name = "Invar";
	if (p_invar_sites >= 0) {
		p_invar = p_invar_sites;
		// true unless -optfromgiven cmd line option
		fix_p_invar = !(Params::getInstance().optimize_from_given_params);
	}
}

void RateInvar::startCheckpoint() {
    checkpoint->startStruct("RateInvar");
}

void RateInvar::saveCheckpoint() {
    startCheckpoint();
    CKP_SAVE(p_invar);
//    CKP_SAVE(fix_p_invar);
//    CKP_SAVE(optimize_p_invar);
    endCheckpoint();
    RateHeterogeneity::saveCheckpoint();
}

void RateInvar::restoreCheckpoint() {
    RateHeterogeneity::restoreCheckpoint();
    startCheckpoint();
    CKP_RESTORE(p_invar);
//    CKP_RESTORE(fix_p_invar);
//    CKP_RESTORE(optimize_p_invar);
    endCheckpoint();
}

string RateInvar::getNameParams() {
	ostringstream str;
	str << "+I{" << p_invar << '}';
	return str.str();
}

double RateInvar::computeFunction(double p_invar_value) {
	p_invar = p_invar_value;
	// fix bug: computeTip... will update ptn_invar vector
//	phylo_tree->computePtnInvar();
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

double RateInvar::targetFunk(double x[]) {
	getVariables(x);
	// fix bug: computeTip... will update ptn_invar vector
	phylo_tree->computePtnInvar();
	return -phylo_tree->computeLikelihood();
}

void RateInvar::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	if (getNDim() == 0) return;
	lower_bound[1] = MIN_PINVAR;
	upper_bound[1] = phylo_tree->aln->frac_const_sites;
	bound_check[1] = false;
}

double RateInvar::optimizeParameters(double gradient_epsilon) {
    if (phylo_tree->aln->frac_const_sites == 0.0)
        return -computeFunction(0.0);
//	if (fix_p_invar || !optimize_p_invar)
	if (fix_p_invar) {
		return -computeFunction(p_invar);
	}
	TREE_LOG_LINE(*phylo_tree, VB_MAX, "Optimizing proportion of invariable sites...");

	double negative_lh;
	double ferror;
	p_invar = minimizeOneDimen(MIN_PINVAR, p_invar, min(phylo_tree->aln->frac_const_sites, 1.0-MIN_PINVAR), max(gradient_epsilon, TOL_PINVAR), &negative_lh, &ferror);
	//p_invar = minimizeOneDimen(MIN_PINVAR, p_invar, 1.0 - MIN_PINVAR, TOL_PINVAR, &negative_lh, &ferror);
//    phylo_tree->clearAllPartialLH();
//	phylo_tree->computePtnInvar();
//	return -negative_lh;
    return -computeFunction(p_invar);
}

void RateInvar::writeInfo(ostream &out) {
	out << "Proportion of invariable sites: " << p_invar << endl;
}

void RateInvar::writeParameters(ostream &out) {
	out << "\t" << p_invar;
}

void RateInvar::setVariables(double *variables) {
	if (RateInvar::getNDim() == 0) return;
	variables[1] = p_invar;
}

bool RateInvar::getVariables(double *variables) {
	if (RateInvar::getNDim() == 0) return false;
    bool changed = (p_invar != variables[1]);
	p_invar = variables[1];
    return changed;
}
