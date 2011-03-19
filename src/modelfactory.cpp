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
#include "modelfactory.h"
#include "rategamma.h"
#include "rategammainvar.h"
#include "gtrmodel.h"
#include "modelnonrev.h"
#include "modeldna.h"
#include "modelprotein.h"
#include "ratemeyerhaeseler.h"
#include "ratemeyerdiscrete.h"
#include "ngs.h"

ModelFactory::ModelFactory() { 
	model = NULL; 
	site_rate = NULL;
	store_trans_matrix = false;
	is_storing = false;
}

ModelFactory::ModelFactory(Params &params, PhyloTree *tree) { 
	store_trans_matrix = params.store_trans_matrix;
	is_storing = false;

	string model_str = params.model_name;
	if (model_str == "") {
		if (tree->aln->num_states == 4) model_str = "HKY";
		else if (tree->aln->num_states == 20) model_str = "WAG";
		else model_str = "JC";
	}
	string::size_type pos = model_str.find('+');

	/* create site-rate heterogeneity */
	if (pos != string::npos) {
		string rate_str = model_str.substr(pos);
		if (rate_str == "+I") {
			site_rate = new RateInvar(params.p_invar_sites, tree);
		} else if (rate_str.substr(0,4) == "+I+G" || rate_str.substr(0,4) == "+G+I") {
			if (rate_str.length() > 4) {
				params.num_rate_cats = convert_int(rate_str.substr(4).c_str());
				if (params.num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateGammaInvar(params.num_rate_cats, params.gamma_shape, params.gamma_median, params.p_invar_sites, tree);
		} else if (rate_str.substr(0,2) == "+G") {
			if (rate_str.length() > 2) {
				params.num_rate_cats = convert_int(rate_str.substr(2).c_str());
				if (params.num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateGamma(params.num_rate_cats, params.gamma_shape, params.gamma_median, tree);
		} else if (rate_str.substr(0,2) == "+M") {
			tree->sse = false;
			params.rate_mh_type = true;
			if (rate_str.length() > 2) {
				params.num_rate_cats = convert_int(rate_str.substr(2).c_str());
				if (params.num_rate_cats < 0) outError("Wrong number of rate categories");
			} else params.num_rate_cats = -1;
			if (params.num_rate_cats >= 0)
				site_rate = new RateMeyerDiscrete(params.num_rate_cats, params.mcat_type, 
					params.rate_file, tree, params.rate_mh_type);
			else
				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
			site_rate->setTree(tree);
		} else if (rate_str.substr(0,4) == "+CAT") {
			tree->sse = false;
			params.rate_mh_type = false;
			if (rate_str.length() > 4) {
				params.num_rate_cats = convert_int(rate_str.substr(4).c_str());
				if (params.num_rate_cats < 0) outError("Wrong number of rate categories");
			} else params.num_rate_cats = -1;
			if (params.num_rate_cats >= 0)
				site_rate = new RateMeyerDiscrete(params.num_rate_cats, params.mcat_type, 
					params.rate_file, tree, params.rate_mh_type);
			else
				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
			site_rate->setTree(tree);
		} else if (rate_str.substr(0,2) == "+F") {
			tree->sse = false;
			if (rate_str.length() > 2) {
				params.num_rate_cats = convert_int(rate_str.substr(2).c_str());
				if (params.num_rate_cats < 0) outError("Wrong number of rate categories");
			} else params.num_rate_cats = -1;
			site_rate = new NGSRate(tree);
			site_rate->setTree(tree);
		} else
			outError("Invalid rate heterogeneity type");
		model_str = model_str.substr(0, pos);
	} else {
		site_rate = new RateHeterogeneity();
		site_rate->setTree(tree);
	} 	

	/* create substitution model */

	if (model_str == "JC" || model_str == "Poisson"/*&& (params.freq_type == FREQ_UNKNOWN || params.freq_type == FREQ_EQUAL)*/) {
		 model = new SubstModel(tree->aln->num_states);
	} else if (model_str == "GTR") {
		model = new GTRModel(tree);
		((GTRModel*)model)->init(params.freq_type);
	} else if (model_str == "UNREST") {
		params.freq_type = FREQ_EQUAL;
		params.optimize_by_newton = false;
		tree->optimize_by_newton = false;
		model = new ModelNonRev(tree);
		((ModelNonRev*)model)->init(params.freq_type);
	} else if (tree->aln->num_states == 4) {
		model = new ModelDNA(model_str.c_str(), params.freq_type, tree);
	} else if (tree->aln->num_states == 20) {
		model = new ModelProtein(model_str.c_str(), params.freq_type, tree);
	} else {
		outError("Unsupported model type");
	}

	tree->discardSaturatedSite(params.discard_saturated_site);

}


double ModelFactory::optimizeParameters(bool fixed_len, bool write_info) {
	assert(model);
	assert(site_rate);

	time_t begin_time, cur_time;
	time(&begin_time);

	double cur_lh;
	PhyloTree *tree = site_rate->getTree();
	assert(tree);

	stopStoringTransMatrix();
	if (fixed_len) 
		cur_lh = tree->computeLikelihood();
	else {
		cur_lh = tree->optimizeAllBranches(1);
	}
	if (verbose_mode >= VB_MED || write_info) 
		cout << "Initial log-likelihood: " << cur_lh << endl;
	int i;
	bool optimize_rate = true;
	for (i = 2; i < 500; i++) {
		double model_lh = model->optimizeParameters();
		//if (model_lh != 0.0) cur_lh = model_lh;
/*
		if (model_lh != 0.0 && !fixed_len)
			model_lh = optimizeAllBranches(3); */
		double rate_lh = 0.0;
		if (optimize_rate) {
			rate_lh = site_rate->optimizeParameters();
			if (rate_lh < model_lh+1e-6) optimize_rate = false;
		}
/*		if (rate_lh != 0.0 && !fixed_len)
			rate_lh = optimizeAllBranches(2);*/
		if (model_lh == 0.0 && rate_lh == 0.0) {
			if (!fixed_len) cur_lh = tree->optimizeAllBranches();
			break;
		}
		double new_lh = (rate_lh != 0.0) ? rate_lh : model_lh;
		
		if (verbose_mode >= VB_MED) {
			model->writeInfo(cout);
			site_rate->writeInfo(cout);
		}
		if (new_lh > cur_lh + 1e-3) {
			if (!fixed_len)
				cur_lh = tree->optimizeAllBranches(i<5 ? i : 5);  // loop only 5 times in total
			else
				cur_lh = new_lh;
			if (verbose_mode >= VB_MED || write_info)
				cout << "Current log-likelihood: " << cur_lh << endl;
		} else {
			site_rate->classifyRates(new_lh);
			if (!fixed_len) cur_lh = tree->optimizeAllBranches();
			break;
		}
	}
	if (verbose_mode == VB_MIN && write_info) {
		model->writeInfo(cout);
		site_rate->writeInfo(cout);
	}
	time(&cur_time);
	double elapsed_secs = difftime(cur_time,begin_time);
	if (write_info)
		cout << "Parameters optimization took " << i-1 << " rounds (" << elapsed_secs << " sec) to finish" << endl << endl;
	startStoringTransMatrix();
	return cur_lh;
}

void ModelFactory::startStoringTransMatrix() {
	if (!store_trans_matrix) return;
	is_storing = true;
}

void ModelFactory::stopStoringTransMatrix() {
	if (!store_trans_matrix) return;
	is_storing = false;
	if (!empty()) {
		for (iterator it = begin(); it != end(); it++)
			delete it->second;
		clear();
	}
}


double ModelFactory::computeTrans(double time, int state1, int state2) {
	return model->computeTrans(time, state1, state2);
}

double ModelFactory::computeTrans(double time, int state1, int state2, double &derv1, double &derv2) {
	return model->computeTrans(time, state1, state2, derv1, derv2);
}

void ModelFactory::computeTransMatrix(double time, double *trans_matrix) {
	if (!store_trans_matrix || !is_storing) {
		model->computeTransMatrix(time, trans_matrix);
		return;
	}
	int mat_size = model->num_states * model->num_states;
	iterator ass_it = find(round(time * 1e6));
	if (ass_it == end()) {
		// allocate memory for 3 matricies
		double *trans_entry = new double[mat_size * 3];
		trans_entry[mat_size] = trans_entry[mat_size+1] = 0.0;
		model->computeTransMatrix(time, trans_entry);
		ass_it = insert(value_type(round(time * 1e6), trans_entry)).first;
	} else {
		//if (verbose_mode >= VB_MAX) 
			//cout << "ModelFactory bingo" << endl;
	} 
	
	memcpy(trans_matrix, ass_it->second, mat_size * sizeof(double));
}

void ModelFactory::computeTransMatrixFreq(double time, double *state_freq, double *trans_matrix) {
	int nstates = model->num_states;
	computeTransMatrix(time, trans_matrix);
	for (int state1 = 0; state1 < nstates; state1++) {
		double *trans_mat_state = trans_matrix + (state1 * nstates);
		for (int state2 = 0; state2 < nstates; state2++)
			trans_mat_state[state2] *= state_freq[state1];
	}
}

void ModelFactory::computeTransDerv(double time, double *trans_matrix, 
	double *trans_derv1, double *trans_derv2) {
	if (!store_trans_matrix || !is_storing) {
		model->computeTransDerv(time, trans_matrix, trans_derv1, trans_derv2);
		return;
	}
	int mat_size = model->num_states * model->num_states;
	iterator ass_it = find(round(time * 1e6));
	if (ass_it == end()) {
		// allocate memory for 3 matricies
		double *trans_entry = new double[mat_size * 3];
		trans_entry[mat_size] = trans_entry[mat_size+1] = 0.0;
		model->computeTransDerv(time, trans_entry, trans_entry+mat_size, trans_entry+(mat_size*2));
		ass_it = insert(value_type(round(time * 1e6), trans_entry)).first;
	} else if (ass_it->second[mat_size] == 0.0 && ass_it->second[mat_size+1] == 0.0) {
		double *trans_entry = ass_it->second;
		model->computeTransDerv(time, trans_entry, trans_entry+mat_size, trans_entry+(mat_size*2));
	}
	memcpy(trans_matrix, ass_it->second, mat_size * sizeof(double));
	memcpy(trans_derv1, ass_it->second + mat_size, mat_size * sizeof(double));
	memcpy(trans_derv2, ass_it->second + (mat_size*2), mat_size * sizeof(double));
}

void ModelFactory::computeTransDervFreq(double time, double rate_val, double *state_freq, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2) 
{
	int nstates = model->num_states;	
	double rate_sqr = rate_val*rate_val;
	computeTransDerv(time * rate_val, trans_matrix, trans_derv1, trans_derv2);
	for (int state1 = 0; state1 < nstates; state1++) {
		double *trans_mat_state = trans_matrix + (state1 * nstates);
		double *trans_derv1_state = trans_derv1 + (state1 * nstates);
		double *trans_derv2_state = trans_derv2 + (state1 * nstates);

		for (int state2 = 0; state2 < nstates; state2++) {
			trans_mat_state[state2] *= state_freq[state1];
			trans_derv1_state[state2] *= (state_freq[state1] * rate_val);
			trans_derv2_state[state2] *= (state_freq[state1] * rate_sqr);
		}
	}
}

ModelFactory::~ModelFactory()
{
	for (iterator it = begin(); it != end(); it++)
		delete it->second;
	clear();
}
