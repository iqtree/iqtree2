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
#include "modeldna.h"
#include "modelprotein.h"

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
			site_rate = new RateGammaInvar(params.num_rate_cats, params.gamma_shape, params.p_invar_sites, tree);
		} else if (rate_str.substr(0,2) == "+G") {
			if (rate_str.length() > 2) {
				params.num_rate_cats = convert_int(rate_str.substr(2).c_str());
				if (params.num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateGamma(params.num_rate_cats, params.gamma_shape, tree);
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
	} else if (tree->aln->num_states == 4) {
		model = new ModelDNA(model_str.c_str(), params.freq_type, tree);
	} else if (tree->aln->num_states == 20) {
		model = new ModelProtein(model_str.c_str(), params.freq_type, tree);
	} else {
		outError("Unsupported model type");
	}
}


double ModelFactory::optimizeParameters(bool fixed_len, bool write_info) {
	assert(model);
	assert(site_rate);

	double cur_lh;
	PhyloTree *tree = site_rate->getTree();
	assert(tree);

	stopStoringTransMatrix();
	if (fixed_len) 
		cur_lh = tree->computeLikelihood();
	else {
		cur_lh = tree->optimizeAllBranches(1);
	}
	int i;
	for (i = 2; i < 100; i++) {
		double model_lh = model->optimizeParameters();
/*
		if (model_lh != 0.0 && !fixed_len)
			model_lh = optimizeAllBranches(3); */
		double rate_lh = site_rate->optimizeParameters();
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
		if (new_lh > cur_lh + 1e-4) {
			if (!fixed_len)
				cur_lh = tree->optimizeAllBranches(i);  // loop only 5 times in total
			else
				cur_lh = new_lh;
			if (verbose_mode > VB_MIN)
				cout << "Current Log-likelihood: " << cur_lh << endl;
		} else {
			if (!fixed_len) cur_lh = tree->optimizeAllBranches();
			break;
		}
	}
	if (verbose_mode == VB_MIN && write_info) {
		model->writeInfo(cout);
		site_rate->writeInfo(cout);
	}
	if (write_info)
		cout << "Optimization took " << i-1 << " rounds to finish" << endl;
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

ModelFactory::~ModelFactory()
{
	for (iterator it = begin(); it != end(); it++)
		delete it->second;
	clear();
}
