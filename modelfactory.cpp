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
#include "modelbin.h"
#include "modelcodon.h"
#include "modelset.h"
#include "ratemeyerhaeseler.h"
#include "ratemeyerdiscrete.h"
#include "ratekategory.h"
#include "ngs.h"
#include <string>
#include "timeutil.h"

const char OPEN_BRACKET = '{';
const char CLOSE_BRACKET = '}';

ModelFactory::ModelFactory() { 
	model = NULL; 
	site_rate = NULL;
	store_trans_matrix = false;
	is_storing = false;
	joint_optimize = false;
}

ModelSubst* ModelFactory::createModel(string model_str, StateFreqType freq_type, string freq_params,
		PhyloTree* tree, bool count_rates)
{
	ModelSubst *model = NULL;
	//cout << "Numstates: " << tree->aln->num_states << endl;
	string model_params;
	size_t pos = model_str.find(OPEN_BRACKET);
	if (pos != string::npos) {
		if (model_str.find(CLOSE_BRACKET) != model_str.length()-1)
			outError("Close bracket not found at the end of ", model_str);
		model_params = model_str.substr(pos+1, model_str.length()-pos-2);
		model_str = model_str.substr(0, pos);
	}

	if ((model_str == "JC" && tree->aln->seq_type == SEQ_DNA) ||
		(model_str == "POISSON" && tree->aln->seq_type == SEQ_PROTEIN) ||
		(model_str == "JC2" && tree->aln->seq_type == SEQ_BINARY) ||
		(model_str == "JCC" && tree->aln->seq_type == SEQ_CODON) ||
		(model_str == "MK" && tree->aln->seq_type == SEQ_MORPH))
	{
		model = new ModelSubst(tree->aln->num_states);
	} else 
	if ((model_str == "GTR" && tree->aln->seq_type == SEQ_DNA) ||
		(model_str == "GTR2" && tree->aln->seq_type == SEQ_BINARY) ||
		(model_str == "GTR20" && tree->aln->seq_type == SEQ_PROTEIN)) {
		model = new GTRModel(tree, count_rates);
		if (freq_params != "")
			((GTRModel*)model)->readStateFreq(freq_params);
		if (model_params != "")
			((GTRModel*)model)->readRates(model_params);
		((GTRModel*)model)->init(freq_type);
	} else if (model_str == "UNREST") {
		freq_type = FREQ_EQUAL;
		//params.optimize_by_newton = false;
		tree->optimize_by_newton = false;
		model = new ModelNonRev(tree, count_rates);
		((ModelNonRev*)model)->init(freq_type);
	} else if (tree->aln->seq_type == SEQ_BINARY) {
		model = new ModelBIN(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_DNA) {

		model = new ModelDNA(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_PROTEIN) {

		model = new ModelProtein(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_CODON) {
		model = new ModelCodon(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else {
		outError("Unsupported model type");
	}

	return model;
}


ModelFactory::ModelFactory(Params &params, PhyloTree *tree) { 
	store_trans_matrix = params.store_trans_matrix;
	is_storing = false;
	joint_optimize = params.optimize_model_rate_joint;

	string model_str = params.model_name;
	if (model_str == "") {
		if (tree->aln->seq_type == SEQ_DNA) model_str = "HKY";
		else if (tree->aln->seq_type == SEQ_PROTEIN) model_str = "WAG";
		else if (tree->aln->seq_type == SEQ_BINARY) model_str = "JC2";
		else if (tree->aln->seq_type == SEQ_CODON) model_str = "JCC";
		else if (tree->aln->seq_type == SEQ_MORPH) model_str = "MK";
		else model_str = "JC";
	}
	string::size_type posfreq;
	StateFreqType freq_type = params.freq_type;

	if (freq_type == FREQ_UNKNOWN) {
		switch (tree->aln->seq_type) {
		case SEQ_BINARY: freq_type = FREQ_ESTIMATE; break; // default for binary: optimized frequencies
		case SEQ_PROTEIN: freq_type = FREQ_USER_DEFINED; break; // default for protein: frequencies of the empirical AA matrix
		case SEQ_MORPH: freq_type = FREQ_EQUAL; break;
		default: freq_type = FREQ_EMPIRICAL; break; // default for DNA and others: counted frequencies from alignment
		}
	}

	size_t close_bracket;
	string freq_params;
	if ((posfreq = model_str.find("+F")) != string::npos) {
		if (model_str.length() > posfreq+2 && model_str[posfreq+2] == OPEN_BRACKET) {
			close_bracket = model_str.find(CLOSE_BRACKET, posfreq);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", model_str);
			if (close_bracket != model_str.length()-1)
				outError("Wrong close bracket position ", model_str);
			freq_type = FREQ_USER_DEFINED;
			freq_params = model_str.substr(posfreq+3, close_bracket-posfreq-3);
		} else if (model_str.substr(posfreq) == "+FC" || model_str.substr(posfreq) == "+Fc" || model_str.substr(posfreq) == "+F")
			freq_type = FREQ_EMPIRICAL;
		else if (model_str.substr(posfreq) == "+FU" || model_str.substr(posfreq) == "+Fu")
			freq_type = FREQ_USER_DEFINED;
		else if (model_str.substr(posfreq) == "+FQ" || model_str.substr(posfreq) == "+Fq")
			freq_type = FREQ_EQUAL;
		else if (model_str.substr(posfreq) == "+FO" || model_str.substr(posfreq) == "+Fo")
			freq_type = FREQ_ESTIMATE;
		else if (model_str.substr(posfreq) == "+F1x4")
			freq_type = FREQ_CODON_1x4;
		else if (model_str.substr(posfreq) == "+F3x4")
			freq_type = FREQ_CODON_3x4;
		else if (model_str.substr(posfreq) == "+F3x4C" || model_str.substr(posfreq) == "+F3x4c")
			freq_type = FREQ_CODON_3x4C;
		else outError("Unknown state frequency type ",model_str.substr(posfreq));
		model_str = model_str.substr(0, posfreq);
	}
	string::size_type posI = model_str.find("+I");
	string::size_type posG = model_str.find("+G");
	string::size_type posX;
	/* create site-rate heterogeneity */
	int num_rate_cats = params.num_rate_cats;
	double gamma_shape = params.gamma_shape;
	double p_invar_sites = params.p_invar_sites;
	if (posI != string::npos) {
		if (model_str.length() > posI+2 && model_str[posI+2] == OPEN_BRACKET) {
			close_bracket = model_str.find(CLOSE_BRACKET, posI);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", model_str);
			p_invar_sites = convert_double(model_str.substr(posI+3, close_bracket-posI-3).c_str());
			if (p_invar_sites <= 0 || p_invar_sites >= 1)
				outError("p_invar must be in (0,1)");
		} else if (model_str.length() > posI+2 && model_str[posI+2] != '+')
			outError("Wrong model name ", model_str);
	}
	if (posG != string::npos) {
		int end_pos = 0;
		if (model_str.length() > posG+2 && isdigit(model_str[posG+2])) {
			num_rate_cats = convert_int(model_str.substr(posG+2).c_str(), end_pos);
				if (num_rate_cats < 1) outError("Wrong number of rate categories");
			}
		if (model_str.length() > posG+2+end_pos && model_str[posG+2+end_pos] == OPEN_BRACKET) {
			close_bracket = model_str.find(CLOSE_BRACKET, posG);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", model_str);
			gamma_shape = convert_double(model_str.substr(posG+3+end_pos, close_bracket-posG-3-end_pos).c_str());
			if (gamma_shape < MIN_GAMMA_SHAPE || gamma_shape > MAX_GAMMA_SHAPE) {
				stringstream str;
				str << "Gamma shape parameter " << gamma_shape << "out of range ["
						<< MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]" << endl;
				outError(str.str());
			}
		} else if (model_str.length() > posG+2+end_pos && model_str[posG+2+end_pos] != '+')
			outError("Wrong model name ", model_str);
	}
	if (model_str.find('+') != string::npos) {
		//string rate_str = model_str.substr(pos);
		if (posI != string::npos && posG != string::npos) {
			site_rate = new RateGammaInvar(num_rate_cats, gamma_shape, params.gamma_median,
					p_invar_sites, params.optimize_model_rate_joint, tree);
		} else if (posI != string::npos) {
			site_rate = new RateInvar(p_invar_sites, tree);
		} else if (posG != string::npos) {
			site_rate = new RateGamma(num_rate_cats, gamma_shape, params.gamma_median, tree);
		} else if ((posX = model_str.find("+M")) != string::npos) {
			tree->sse = false;
			params.rate_mh_type = true;
			if (model_str.length() > posX+2 && isdigit(model_str[posX+2])) {
				num_rate_cats = convert_int(model_str.substr(posX+2).c_str());
				if (num_rate_cats < 0) outError("Wrong number of rate categories");
			} else num_rate_cats = -1;
			if (num_rate_cats >= 0)
				site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type, 
					params.rate_file, tree, params.rate_mh_type);
			else
				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
			site_rate->setTree(tree);
		} else if ((posX = model_str.find("+D")) != string::npos) {
			tree->sse = false;
			params.rate_mh_type = false;
			if (model_str.length() > posX+2 && isdigit(model_str[posX+2])) {
				num_rate_cats = convert_int(model_str.substr(posX+2).c_str());
				if (num_rate_cats < 0) outError("Wrong number of rate categories");
			} else num_rate_cats = -1;
			if (num_rate_cats >= 0)
				site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type, 
					params.rate_file, tree, params.rate_mh_type);
			else
				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
			site_rate->setTree(tree);
		} else if ((posX = model_str.find("+NGS")) != string::npos) {
			tree->sse = false;
			if (model_str.length() > posX+4 && isdigit(model_str[posX+4])) {
				num_rate_cats = convert_int(model_str.substr(posX+4).c_str());
				if (num_rate_cats < 0) outError("Wrong number of rate categories");
			} else num_rate_cats = -1;
			site_rate = new NGSRateCat(tree, num_rate_cats);
			site_rate->setTree(tree);
		} else if ((posX = model_str.find("+NGS")) != string::npos) {
			tree->sse = false;
			if (model_str.length() > posX+4 && isdigit(model_str[posX+4])) {
				num_rate_cats = convert_int(model_str.substr(posX+4).c_str());
				if (num_rate_cats < 0) outError("Wrong number of rate categories");
			} else num_rate_cats = -1;
			site_rate = new NGSRate(tree);
			site_rate->setTree(tree);
		} else if ((posX = model_str.find("+K")) != string::npos) {
			if (model_str.length() > posX+2 && isdigit(model_str[posX+2])) {
				num_rate_cats = convert_int(model_str.substr(posX+2).c_str());
				if (num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateKategory(num_rate_cats, tree);
		} else
			outError("Invalid rate heterogeneity type");
		model_str = model_str.substr(0, model_str.find('+'));
	} else {
		site_rate = new RateHeterogeneity();
		site_rate->setTree(tree);
	} 	

	/* create substitution model */

	if (!params.site_freq_file) {
		model = createModel(model_str, freq_type, freq_params, tree);
	} else { 
		// site-specific model
		if (model_str == "JC" || model_str == "POSSION") 
			outError("JC is not suitable for site-specific model");
		model = new ModelSet(model_str.c_str(), tree);
		ModelSet *models = (ModelSet*)model; // assign pointer for convenience
		models->init(params.freq_type);
		IntVector site_model;
		vector<double*> freq_vec;
		readSiteFreq(tree->aln, params.site_freq_file, site_model, freq_vec);
		tree->aln->regroupSitePattern(freq_vec.size(), site_model);
		//tree->aln->ungroupSitePattern();
		tree->setAlignment(tree->aln);
		int i;
		models->pattern_model_map.resize(tree->aln->getNPattern(), -1);
		for (i = 0; i < tree->aln->getNSite(); i++) {
			models->pattern_model_map[tree->aln->getPatternID(i)] = site_model[i];
			//cout << "site " << i << " ptn " << tree->aln->getPatternID(i) << " -> model " << site_model[i] << endl;
		}
		double *state_freq = new double[model->num_states];
		double *rates = new double[model->getNumRateEntries()];
		for (i = 0; i < freq_vec.size(); i++) {
			GTRModel *modeli;
			if (i == 0) {
				modeli = (GTRModel*)createModel(model_str, params.freq_type, "", tree, true);
				modeli->getStateFrequency(state_freq);
				modeli->getRateMatrix(rates);
			} else {
				modeli = (GTRModel*)createModel(model_str, FREQ_EQUAL, "", tree, false);
				modeli->setStateFrequency(state_freq);
				modeli->setRateMatrix(rates);
			}
			if (freq_vec[i])
				modeli->setStateFrequency (freq_vec[i]);

			modeli->init(FREQ_USER_DEFINED);
			models->push_back(modeli);
		}
		delete [] rates;
		delete [] state_freq;
		cout << "Alignment is divided into " << models->size() << " partitions with " << tree->aln->getNPattern() << " patterns" << endl;
		for (vector<double*>::reverse_iterator it = freq_vec.rbegin(); it != freq_vec.rend(); it++)
			if (*it) delete [] (*it);
	} 
	tree->discardSaturatedSite(params.discard_saturated_site);

}

int ModelFactory::getNParameters() {
	int df = model->getNDim() + site_rate->getNDim() + site_rate->phylo_tree->branchNum;
	if (model->freq_type == FREQ_EMPIRICAL) df += model->num_states-1;
	return df;
}
void ModelFactory::readSiteFreq(Alignment *aln, char* site_freq_file, IntVector &site_model, vector<double*> &freq_vec)
{
	cout << "Reading site-specific state frequency file " << site_freq_file << " ..." << endl;
	site_model.resize(aln->getNSite(), -1);
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(site_freq_file);
		double freq;
		string site_spec;
		int specified_sites = 0;
		in.exceptions(ios::badbit);
		for (int model_id = 0; !in.eof(); model_id++) {
			// remove the failbit
			in >> site_spec;
			if (in.eof()) break;
			IntVector site_id;
			extractSiteID(aln, site_spec.c_str(), site_id);
			specified_sites += site_id.size();
			if (site_id.size() == 0) throw "No site ID specified";
			for (IntVector::iterator it = site_id.begin(); it != site_id.end(); it++) {
				if (site_model[*it] != -1) throw "Duplicated site ID";
				site_model[*it] = model_id;
			}
			double *site_freq_entry = new double[aln->num_states];
			double sum = 0;
			for (int i = 0; i < aln->num_states; i++) {
				in >> freq;
				if (freq <= 0.0 || freq >= 1.0) throw "Invalid frequency entry";
				site_freq_entry[i] = freq;
				sum += freq;
			}
			if (fabs(sum-1.0) > 1e-4) throw "Frequencies do not sum up to 1";
			aln->convfreq(site_freq_entry); // regularize frequencies (eg if some freq = 0)
			freq_vec.push_back(site_freq_entry);
		}
		if (specified_sites < site_model.size()) {
			// there are some unspecified sites
			cout << site_model.size() - specified_sites << " unspecified sites will get default frequencies" << endl;
			for (int i = 0; i < site_model.size(); i++)
				if (site_model[i] == -1) 
					site_model[i] = freq_vec.size();
			freq_vec.push_back(NULL);
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}
}

double ModelFactory::optimizeParametersOnly(double epsilon) {
	if (!joint_optimize) {
		double model_lh = model->optimizeParameters(epsilon);
		double rate_lh = site_rate->optimizeParameters(epsilon);
		if (rate_lh == 0.0) return model_lh;
		return rate_lh;
	}

	int ndim = getNDim();

	// return if nothing to be optimized
	if (ndim == 0) return 0.0;

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	int i;
	double score;

	// setup the bounds for model
	setVariables(variables);
	int model_ndim = model->getNDim();
	for (i = 1; i <= model_ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = MIN_RATE;
		upper_bound[i] = MAX_RATE;
		bound_check[i] = false;
	}

	if (model->freq_type == FREQ_ESTIMATE) {
		for (i = model_ndim-model->num_states+2; i <= model_ndim; i++)
			upper_bound[i] = 1.0;
	}

	// setup the bounds for site_rate
	site_rate->setBounds(lower_bound+model_ndim, upper_bound+model_ndim, bound_check+model_ndim);

	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(epsilon, TOL_RATE));

	getVariables(variables);
	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(true);
	model->decomposeRateMatrix();
	site_rate->phylo_tree->clearAllPartialLH();

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}



double ModelFactory::optimizeParameters(bool fixed_len, bool write_info, double logl_epsilon) {
	assert(model);
	assert(site_rate);

	//time_t begin_time, cur_time;
	//time(&begin_time);
	double begin_time = getCPUTime();
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
		cout << "1. Initial log-likelihood: " << cur_lh << endl;
	int i;
	//bool optimize_rate = true;
	double param_epsilon = logl_epsilon; // epsilon for parameters starts at epsilon for logl
	for (i = 2; i < 100; i++, param_epsilon/=4.0) {
		/*
		double model_lh = model->optimizeParameters(param_epsilon);
		double rate_lh = 0.0;
		if (optimize_rate) {
			rate_lh = site_rate->optimizeParameters(param_epsilon);
			if (rate_lh < model_lh+1e-6 && model_lh != 0.0) optimize_rate = false;
		}
		if (model_lh == 0.0 && rate_lh == 0.0) {
			if (!fixed_len) cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
			break;
		}
		double new_lh = (rate_lh != 0.0) ? rate_lh : model_lh;
		*/
		double new_lh = optimizeParametersOnly(param_epsilon);
		if (new_lh == 0.0) {
			if (!fixed_len) cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
			break;
		}
		if (verbose_mode >= VB_MED) {
			model->writeInfo(cout);
			site_rate->writeInfo(cout);
		}
		if (!fixed_len)
			new_lh = tree->optimizeAllBranches(min(i,3), logl_epsilon);  // loop only 3 times in total (previously in v0.9.6 5 times)
		if (new_lh > cur_lh + logl_epsilon) {
			if (param_epsilon > (new_lh - cur_lh) * logl_epsilon)
				param_epsilon = (new_lh - cur_lh) * logl_epsilon;
			cur_lh = new_lh;
			if (verbose_mode >= VB_MED || write_info)
				cout << i << ". Current log-likelihood: " << cur_lh << endl;
		} else {
			site_rate->classifyRates(new_lh);
			if (!fixed_len) cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
				break;
		}
	}
	if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << cur_lh << endl;
	if (verbose_mode <= VB_MIN && write_info) {
		model->writeInfo(cout);
		site_rate->writeInfo(cout);
	}
	//time(&cur_time);
	//double elapsed_secs = difftime(cur_time,begin_time);
	double elapsed_secs = getCPUTime() - begin_time;
	if (write_info)
		cout << "Parameters optimization took " << i-1 << " rounds (" << elapsed_secs << " sec)" << endl << endl;
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
	if (!store_trans_matrix || !is_storing || model->isSiteSpecificModel()) {
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
	if (model->isSiteSpecificModel()) {
		model->computeTransMatrixFreq(time, trans_matrix);
		return;
	}
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
	if (!store_trans_matrix || !is_storing || model->isSiteSpecificModel()) {
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
	if (model->isSiteSpecificModel()) {
		model->computeTransDervFreq(time, rate_val, trans_matrix, trans_derv1, trans_derv2);
		return;
	}
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

/************* FOLLOWING SERVE FOR JOINT OPTIMIZATION OF MODEL AND RATE PARAMETERS *******/
int ModelFactory::getNDim()
{
	return model->getNDim() + site_rate->getNDim();
}

double ModelFactory::targetFunk(double x[]) {
	model->getVariables(x);
	// need to compute rates again if p_inv or Gamma shape changes!
	if (model->state_freq[model->num_states-1] < MIN_RATE) return 1.0e+12;
	model->decomposeRateMatrix();
	site_rate->phylo_tree->clearAllPartialLH();
	return site_rate->targetFunk(x + model->getNDim());
}

void ModelFactory::setVariables(double *variables) {
	model->setVariables(variables);
	site_rate->setVariables(variables + model->getNDim());
}

void ModelFactory::getVariables(double *variables) {
	model->getVariables(variables);
	site_rate->getVariables(variables + model->getNDim());
}

