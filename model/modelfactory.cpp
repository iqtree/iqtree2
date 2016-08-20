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
#include "modelgtr.h"
#include "modelnonrev.h"
#include "modeldna.h"
#include "modelprotein.h"
#include "modelbin.h"
#include "modelcodon.h"
#include "modelmorphology.h"
#include "modelset.h"
#include "modelmixture.h"
#include "ratemeyerhaeseler.h"
#include "ratemeyerdiscrete.h"
#include "ratekategory.h"
#include "ratefree.h"
#include "ratefreeinvar.h"
#include "ngs.h"
#include <string>
#include "timeutil.h"
#include "myreader.h"
#include <sstream>

ModelsBlock *readModelsDefinition(Params &params) {

	ModelsBlock *models_block = new ModelsBlock;

	try
	{
		// loading internal model definitions
		stringstream in(builtin_mixmodels_definition);
        assert(in && "stringstream is OK");
		NxsReader nexus;
		nexus.Add(models_block);
	    MyToken token(in);
	    nexus.Execute(token);
//	    int num_model = 0, num_freq = 0;
//	    for (ModelsBlock::iterator it = models_block->begin(); it != models_block->end(); it++)
//	    	if ((*it).flag & NM_FREQ) num_freq++; else num_model++;
//	    cout << num_model << " models and " << num_freq << " frequency vectors loaded" << endl;
	} catch (...) {
        assert(0 && "predefined mixture models initialized");
    }

	if (params.model_def_file) {
		cout << "Reading model definition file " << params.model_def_file << " ... ";
		MyReader nexus(params.model_def_file);
		nexus.Add(models_block);
	    MyToken token(nexus.inf);
	    nexus.Execute(token);
	    int num_model = 0, num_freq = 0;
	    for (ModelsBlock::iterator it = models_block->begin(); it != models_block->end(); it++)
	    	if ((*it).flag & NM_FREQ) num_freq++; else num_model++;
	    cout << num_model << " models and " << num_freq << " frequency vectors loaded" << endl;
	}
	return models_block;
}

ModelFactory::ModelFactory() : CheckpointFactory() { 
	model = NULL; 
	site_rate = NULL;
	store_trans_matrix = false;
	is_storing = false;
	joint_optimize = false;
	fused_mix_rate = false;
	unobserved_ptns = "";
}

size_t findCloseBracket(string &str, size_t start_pos) {
	int counter = 0;
	for (size_t pos = start_pos+1; pos < str.length(); pos++) {
		if (str[pos] == '{') counter++;
		if (str[pos] == '}') {
			if (counter == 0) return pos; else counter--;
		}
	}
	return string::npos;
}

ModelFactory::ModelFactory(Params &params, PhyloTree *tree, ModelsBlock *models_block) : CheckpointFactory() {
	store_trans_matrix = params.store_trans_matrix;
	is_storing = false;
	joint_optimize = params.optimize_model_rate_joint;
	fused_mix_rate = false;
    string model_str = params.model_name;
	string rate_str;

	try {


	if (model_str == "") {
		if (tree->aln->seq_type == SEQ_DNA) model_str = "HKY";
		else if (tree->aln->seq_type == SEQ_PROTEIN) model_str = "LG";
		else if (tree->aln->seq_type == SEQ_BINARY) model_str = "GTR2";
		else if (tree->aln->seq_type == SEQ_CODON) model_str = "GY";
		else if (tree->aln->seq_type == SEQ_MORPH) model_str = "MK";
		else model_str = "JC";
		outWarning("Default model "+model_str + " may be under-fitting. Use option '-m TEST' to determine the best-fit model.");
	}

	/********* preprocessing model string ****************/
	NxsModel *nxsmodel  = NULL;

    string new_model_str = "";
    size_t mix_pos;
    for (mix_pos = 0; mix_pos < model_str.length(); mix_pos++) {
        size_t next_mix_pos = model_str.find_first_of("+*", mix_pos);
        string sub_model_str = model_str.substr(mix_pos, next_mix_pos-mix_pos);
        nxsmodel = models_block->findMixModel(sub_model_str);
        if (nxsmodel) sub_model_str = nxsmodel->description;
        new_model_str += sub_model_str;
        if (next_mix_pos != string::npos)
            new_model_str += model_str[next_mix_pos];
        else 
            break;
        mix_pos = next_mix_pos;
    }
    if (new_model_str != model_str)
        cout << "Model " << model_str << " is alias for " << new_model_str << endl;
    model_str = new_model_str;
    
//	nxsmodel = models_block->findModel(model_str);
//	if (nxsmodel && nxsmodel->description.find_first_of("+*") != string::npos) {
//		cout << "Model " << model_str << " is alias for " << nxsmodel->description << endl;
//		model_str = nxsmodel->description;
//	}

	// decompose model string into model_str and rate_str string
	size_t spec_pos = model_str.find_first_of("{+*");
	if (spec_pos != string::npos) {
		if (model_str[spec_pos] == '{') {
			// scan for the corresponding '}'
			size_t pos = findCloseBracket(model_str, spec_pos);
			if (pos == string::npos)
				outError("Model name has wrong bracket notation '{...}'");
			rate_str = model_str.substr(pos+1);
			model_str = model_str.substr(0, pos+1);
		} else {
			rate_str = model_str.substr(spec_pos);
			model_str = model_str.substr(0, spec_pos);
		}
	}

//	nxsmodel = models_block->findModel(model_str);
//	if (nxsmodel && nxsmodel->description.find("MIX") != string::npos) {
//		cout << "Model " << model_str << " is alias for " << nxsmodel->description << endl;
//		model_str = nxsmodel->description;
//	}

	/******************** initialize state frequency ****************************/

	StateFreqType freq_type = params.freq_type;

	if (freq_type == FREQ_UNKNOWN) {
		switch (tree->aln->seq_type) {
		case SEQ_BINARY: freq_type = FREQ_ESTIMATE; break; // default for binary: optimized frequencies
		case SEQ_PROTEIN: freq_type = FREQ_USER_DEFINED; break; // default for protein: frequencies of the empirical AA matrix
		case SEQ_MORPH: freq_type = FREQ_EQUAL; break;
		case SEQ_CODON: freq_type = FREQ_UNKNOWN; break;
		default: freq_type = FREQ_EMPIRICAL; break; // default for DNA and others: counted frequencies from alignment
		}
	}

    // first handle mixture frequency
    string::size_type posfreq = rate_str.find("+FMIX");
	string freq_params;
    size_t close_bracket;

    if (posfreq != string::npos) {
		string freq_str;
		size_t last_pos = rate_str.find_first_of("+*", posfreq+1);
        
		if (last_pos == string::npos) {
			freq_str = rate_str.substr(posfreq);
			rate_str = rate_str.substr(0, posfreq);
		} else {
			freq_str = rate_str.substr(posfreq, last_pos-posfreq);
			rate_str = rate_str.substr(0, posfreq) + rate_str.substr(last_pos);
		}
        
        if (freq_str[5] != OPEN_BRACKET)
            outError("Mixture-frequency must start with +FMIX{");
        close_bracket = freq_str.find(CLOSE_BRACKET);
        if (close_bracket == string::npos)
            outError("Close bracket not found in ", freq_str);
        if (close_bracket != freq_str.length()-1)
            outError("Wrong close bracket position ", freq_str);
        freq_type = FREQ_MIXTURE;
        freq_params = freq_str.substr(6, close_bracket-6);
    }

    // then normal frequency
    if (rate_str.find("+FO") != string::npos)
        posfreq = rate_str.find("+FO");
    else if (rate_str.find("+Fo") != string::npos)
        posfreq = rate_str.find("+Fo");
    else
        posfreq = rate_str.find("+F");
        
    bool optimize_mixmodel_weight = params.optimize_mixmodel_weight;

	if (posfreq != string::npos) {
		string freq_str;
		size_t last_pos = rate_str.find_first_of("+*", posfreq+1);
		if (last_pos == string::npos) {
			freq_str = rate_str.substr(posfreq);
			rate_str = rate_str.substr(0, posfreq);
		} else {
			freq_str = rate_str.substr(posfreq, last_pos-posfreq);
			rate_str = rate_str.substr(0, posfreq) + rate_str.substr(last_pos);
		}

        if (freq_str.length() > 2 && freq_str[2] == OPEN_BRACKET) {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with user-defined frequency is not allowed");
			close_bracket = freq_str.find(CLOSE_BRACKET);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", freq_str);
			if (close_bracket != freq_str.length()-1)
				outError("Wrong close bracket position ", freq_str);
			freq_type = FREQ_USER_DEFINED;
			freq_params = freq_str.substr(3, close_bracket-3);
		} else if (freq_str == "+FC" || freq_str == "+Fc" || freq_str == "+F") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "empirical," + freq_params;
                optimize_mixmodel_weight = true;
            } else
                freq_type = FREQ_EMPIRICAL;
		} else if (freq_str == "+FU" || freq_str == "+Fu") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with user-defined frequency is not allowed");
            else
                freq_type = FREQ_USER_DEFINED;
		} else if (freq_str == "+FQ" || freq_str == "+Fq") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with equal frequency is not allowed");
            else
                freq_type = FREQ_EQUAL;
		} else if (freq_str == "+FO" || freq_str == "+Fo") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "optimize," + freq_params;
                optimize_mixmodel_weight = true;                
            } else
                freq_type = FREQ_ESTIMATE;
		} else if (freq_str == "+F1x4" || freq_str == "+F1X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + freq_str + " is not allowed");
            else
                freq_type = FREQ_CODON_1x4;
		} else if (freq_str == "+F3x4" || freq_str == "+F3X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + freq_str + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4;
		} else if (freq_str == "+F3x4C" || freq_str == "+F3x4c" || freq_str == "+F3X4C" || freq_str == "+F3X4c") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + freq_str + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4C;
		} else outError("Unknown state frequency type ",freq_str);
//		model_str = model_str.substr(0, posfreq);
	}

	/******************** initialize model ****************************/

	if (tree->aln->site_state_freq.empty()) {
		if (model_str.substr(0, 3) == "MIX" || freq_type == FREQ_MIXTURE) {
			string model_list;
			if (model_str.substr(0, 3) == "MIX") {
				if (model_str[3] != OPEN_BRACKET)
					outError("Mixture model name must start with 'MIX{'");
				if (model_str.rfind(CLOSE_BRACKET) != model_str.length()-1)
					outError("Close bracket not found at the end of ", model_str);
				model_list = model_str.substr(4, model_str.length()-5);
				model_str = model_str.substr(0, 3);
			}
			model = new ModelMixture(params.model_name, model_str, model_list, models_block, freq_type, freq_params, tree, optimize_mixmodel_weight);
		} else {
//			string model_desc;
//			NxsModel *nxsmodel = models_block->findModel(model_str);
//			if (nxsmodel) model_desc = nxsmodel->description;
			model = createModel(model_str, models_block, freq_type, freq_params, tree);
		}

//		fused_mix_rate &= model->isMixture() && site_rate->getNRate() > 1;
	} else {
		// site-specific model
		if (model_str == "JC" || model_str == "POISSON")
			outError("JC is not suitable for site-specific model");
		model = new ModelSet(model_str.c_str(), tree);
		ModelSet *models = (ModelSet*)model; // assign pointer for convenience
		models->init((params.freq_type != FREQ_UNKNOWN) ? params.freq_type : FREQ_EMPIRICAL);
		int i;
		models->pattern_model_map.resize(tree->aln->getNPattern(), -1);
		for (i = 0; i < tree->aln->getNSite(); i++) {
			models->pattern_model_map[tree->aln->getPatternID(i)] = tree->aln->site_model[i];
			//cout << "site " << i << " ptn " << tree->aln->getPatternID(i) << " -> model " << site_model[i] << endl;
		}
		double *state_freq = new double[model->num_states];
		double *rates = new double[model->getNumRateEntries()];
		for (i = 0; i < tree->aln->site_state_freq.size(); i++) {
			ModelGTR *modeli;
			if (i == 0) {
				modeli = (ModelGTR*)createModel(model_str, models_block, (params.freq_type != FREQ_UNKNOWN) ? params.freq_type : FREQ_EMPIRICAL, "", tree, true);
				modeli->getStateFrequency(state_freq);
				modeli->getRateMatrix(rates);
			} else {
				modeli = (ModelGTR*)createModel(model_str, models_block, FREQ_EQUAL, "", tree, false);
				modeli->setStateFrequency(state_freq);
				modeli->setRateMatrix(rates);
			}
			if (tree->aln->site_state_freq[i])
				modeli->setStateFrequency (tree->aln->site_state_freq[i]);

			modeli->init(FREQ_USER_DEFINED);
			models->push_back(modeli);
		}
		delete [] rates;
		delete [] state_freq;

        // delete information of the old alignment
//        tree->aln->ordered_pattern.clear();
//        tree->deleteAllPartialLh();
	}
    
//	if (model->isMixture())
//		cout << "Mixture model with " << model->getNMixtures() << " components!" << endl;

	/******************** initialize ascertainment bias correction model ****************************/

	string::size_type posasc;

	if ((posasc = rate_str.find("+ASC")) != string::npos) {
		// ascertainment bias correction
		unobserved_ptns = tree->aln->getUnobservedConstPatterns();
		// rebuild the seq_states to contain states of unobserved constant patterns
		tree->aln->buildSeqStates(true);
//		if (unobserved_ptns.size() <= 0)
//			outError("Invalid use of +ASC because all constant patterns are observed in the alignment");
		if (unobserved_ptns.size() < tree->aln->getNumNonstopCodons())
			outError("Invalid use of +ASC because constant patterns are observed in the alignment");
		cout << "Ascertainment bias correction: " << unobserved_ptns.size() << " unobservable constant patterns"<< endl;
		rate_str = rate_str.substr(0, posasc) + rate_str.substr(posasc+4);
	}


	/******************** initialize site rate heterogeneity ****************************/

	string::size_type posI = rate_str.find("+I");
	string::size_type posG = rate_str.find("+G");
	string::size_type posG2 = rate_str.find("*G");
    if (posG != string::npos && posG2 != string::npos) {
        cout << "NOTE: both +G and *G were specified, continue with " 
            << ((posG < posG2)? rate_str.substr(posG,2) : rate_str.substr(posG2,2)) << endl;
    }
    if (posG2 != string::npos && posG2 < posG) {
        posG = posG2;
        fused_mix_rate = true;
    }
//	if (posG == string::npos) {
//		posG = rate_str.find("*G");
//		if (posG != string::npos)
//			fused_mix_rate = true;
//	}
	string::size_type posR = rate_str.find("+R"); // FreeRate model
	string::size_type posR2 = rate_str.find("*R"); // FreeRate model
    if (posR != string::npos && posR2 != string::npos) {
        cout << "NOTE: both +R and *R were specified, continue with " 
            << ((posR < posR2)? rate_str.substr(posR,2) : rate_str.substr(posR2,2)) << endl;
    }
    if (posR2 != string::npos && posR2 < posR) {
        posR = posR2;
        fused_mix_rate = true;
    }
    
//	if (posR == string::npos) {
//		posR = rate_str.find("*R");
//		if (posR != string::npos)
//			fused_mix_rate = true;
//	}
	if (posG != string::npos && posR != string::npos) {
        if (posG == posG2 && posR != posR2) {
            outWarning("Both Gamma and FreeRate models were specified, continue with Gamma model because *G has higher priority than +R");
            posR = string::npos;
        } else {
            outWarning("Both Gamma and FreeRate models were specified, continue with FreeRate model");
            posG = string::npos;
        }
    }
	string::size_type posX;
	/* create site-rate heterogeneity */
	int num_rate_cats = params.num_rate_cats;
	if (fused_mix_rate) num_rate_cats = model->getNMixtures();
	double gamma_shape = params.gamma_shape;
	double p_invar_sites = params.p_invar_sites;
	string freerate_params = "";
	if (posI != string::npos) {
		// invariable site model
		if (rate_str.length() > posI+2 && rate_str[posI+2] == OPEN_BRACKET) {
			close_bracket = rate_str.find(CLOSE_BRACKET, posI);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", rate_str);
			p_invar_sites = convert_double(rate_str.substr(posI+3, close_bracket-posI-3).c_str());
			if (p_invar_sites < 0 || p_invar_sites >= 1)
				outError("p_invar must be in [0,1)");
		} else if (rate_str.length() > posI+2 && rate_str[posI+2] != '+')
			outError("Wrong model name ", rate_str);
	}
	if (posG != string::npos) {
		// Gamma rate model
		int end_pos = 0;
		if (rate_str.length() > posG+2 && isdigit(rate_str[posG+2])) {
			num_rate_cats = convert_int(rate_str.substr(posG+2).c_str(), end_pos);
			if (num_rate_cats < 1) outError("Wrong number of rate categories");
		}
		if (rate_str.length() > posG+2+end_pos && rate_str[posG+2+end_pos] == OPEN_BRACKET) {
			close_bracket = rate_str.find(CLOSE_BRACKET, posG);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", rate_str);
			gamma_shape = convert_double(rate_str.substr(posG+3+end_pos, close_bracket-posG-3-end_pos).c_str());
//			if (gamma_shape < MIN_GAMMA_SHAPE || gamma_shape > MAX_GAMMA_SHAPE) {
//				stringstream str;
//				str << "Gamma shape parameter " << gamma_shape << "out of range ["
//						<< MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]" << endl;
//				outError(str.str());
//			}
		} else if (rate_str.length() > posG+2+end_pos && rate_str[posG+2+end_pos] != '+')
			outError("Wrong model name ", rate_str);
	}
	if (posR != string::npos) {
		// FreeRate model
		int end_pos = 0;
		if (rate_str.length() > posR+2 && isdigit(rate_str[posR+2])) {
			num_rate_cats = convert_int(rate_str.substr(posR+2).c_str(), end_pos);
				if (num_rate_cats < 1) outError("Wrong number of rate categories");
			}
		if (rate_str.length() > posR+2+end_pos && rate_str[posR+2+end_pos] == OPEN_BRACKET) {
			close_bracket = rate_str.find(CLOSE_BRACKET, posR);
			if (close_bracket == string::npos)
				outError("Close bracket not found in ", rate_str);
			freerate_params = rate_str.substr(posR+3+end_pos, close_bracket-posR-3-end_pos).c_str();
		} else if (rate_str.length() > posR+2+end_pos && rate_str[posR+2+end_pos] != '+')
			outError("Wrong model name ", rate_str);
	}
	if (rate_str.find('+') != string::npos || rate_str.find('*') != string::npos) {
		//string rate_str = model_str.substr(pos);
		if (posI != string::npos && posG != string::npos) {
			site_rate = new RateGammaInvar(num_rate_cats, gamma_shape, params.gamma_median,
					p_invar_sites, params.optimize_alg_gammai, tree, false);
		} else if (posI != string::npos && posR != string::npos) {
			site_rate = new RateFreeInvar(num_rate_cats, gamma_shape, freerate_params, !fused_mix_rate, p_invar_sites, params.optimize_alg, tree);
		} else if (posI != string::npos) {
			site_rate = new RateInvar(p_invar_sites, tree);
		} else if (posG != string::npos) {
			site_rate = new RateGamma(num_rate_cats, gamma_shape, params.gamma_median, tree);
		} else if (posR != string::npos) {
			site_rate = new RateFree(num_rate_cats, gamma_shape, freerate_params, !fused_mix_rate, params.optimize_alg, tree);
//		} else if ((posX = rate_str.find("+M")) != string::npos) {
//			tree->setLikelihoodKernel(LK_NORMAL);
//			params.rate_mh_type = true;
//			if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
//				num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
//				if (num_rate_cats < 0) outError("Wrong number of rate categories");
//			} else num_rate_cats = -1;
//			if (num_rate_cats >= 0)
//				site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type, 
//					params.rate_file, tree, params.rate_mh_type);
//			else
//				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
//			site_rate->setTree(tree);
//		} else if ((posX = rate_str.find("+D")) != string::npos) {
//			tree->setLikelihoodKernel(LK_NORMAL);
//			params.rate_mh_type = false;
//			if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
//				num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
//				if (num_rate_cats < 0) outError("Wrong number of rate categories");
//			} else num_rate_cats = -1;
//			if (num_rate_cats >= 0)
//				site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type, 
//					params.rate_file, tree, params.rate_mh_type);
//			else
//				site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
//			site_rate->setTree(tree);
//		} else if ((posX = rate_str.find("+NGS")) != string::npos) {
//			tree->setLikelihoodKernel(LK_NORMAL);
//			if (rate_str.length() > posX+4 && isdigit(rate_str[posX+4])) {
//				num_rate_cats = convert_int(rate_str.substr(posX+4).c_str());
//				if (num_rate_cats < 0) outError("Wrong number of rate categories");
//			} else num_rate_cats = -1;
//			site_rate = new NGSRateCat(tree, num_rate_cats);
//			site_rate->setTree(tree);
//		} else if ((posX = rate_str.find("+NGS")) != string::npos) {
//			tree->setLikelihoodKernel(LK_NORMAL);
//			if (rate_str.length() > posX+4 && isdigit(rate_str[posX+4])) {
//				num_rate_cats = convert_int(rate_str.substr(posX+4).c_str());
//				if (num_rate_cats < 0) outError("Wrong number of rate categories");
//			} else num_rate_cats = -1;
//			site_rate = new NGSRate(tree);
//			site_rate->setTree(tree);
		} else if ((posX = rate_str.find("+K")) != string::npos) {
			if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
				num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
				if (num_rate_cats < 1) outError("Wrong number of rate categories");
			}
			site_rate = new RateKategory(num_rate_cats, tree);
		} else
			outError("Invalid rate heterogeneity type");
//		if (model_str.find('+') != string::npos)
//			model_str = model_str.substr(0, model_str.find('+'));
//		else
//			model_str = model_str.substr(0, model_str.find('*'));
	} else {
		site_rate = new RateHeterogeneity();
		site_rate->setTree(tree);
	} 	

	if (fused_mix_rate) {
		if (!model->isMixture())
			outError("Model is not a mixture model");
		if (model->getNMixtures() != site_rate->getNRate())
			outError("Mixture model and site rate model do not have the same number of categories");
		ModelMixture *mmodel = (ModelMixture*)model;
		// reset mixture model
		mmodel->fix_prop = true;
		for (ModelMixture::iterator it = mmodel->begin(); it != mmodel->end(); it++) {
			(*it)->total_num_subst = 1.0;
			mmodel->prop[it-mmodel->begin()] = 1.0;
		}
		mmodel->decomposeRateMatrix();
	}

	tree->discardSaturatedSite(params.discard_saturated_site);

	} catch (const char* str) {
		outError(str);
	}

}

void ModelFactory::setCheckpoint(Checkpoint *checkpoint) {
	CheckpointFactory::setCheckpoint(checkpoint);
	model->setCheckpoint(checkpoint);
	site_rate->setCheckpoint(checkpoint);
}

void ModelFactory::saveCheckpoint() {
    model->saveCheckpoint();
    site_rate->saveCheckpoint();
    checkpoint->startStruct("ModelFactory");
//    CKP_SAVE(fused_mix_rate);
//    CKP_SAVE(unobserved_ptns);
//    CKP_SAVE(joint_optimize);
    checkpoint->endStruct();
    CheckpointFactory::saveCheckpoint();
}

void ModelFactory::restoreCheckpoint() {
    model->restoreCheckpoint();
    site_rate->restoreCheckpoint();
    checkpoint->startStruct("ModelFactory");
//    CKP_RESTORE(fused_mix_rate);
//    CKP_RESTORE(unobserved_ptns);
//    CKP_RESTORE(joint_optimize);
    checkpoint->endStruct();
}

int ModelFactory::getNParameters() {
	int df = model->getNDim() + model->getNDimFreq() + site_rate->getNDim() + site_rate->phylo_tree->branchNum;
	return df;
}

double ModelFactory::initGTRGammaIParameters(RateHeterogeneity *rate, ModelSubst *model, double initAlpha,
                                           double initPInvar, double *initRates, double *initStateFreqs)  {

    RateHeterogeneity* rateGammaInvar = rate;
    ModelGTR* modelGTR = (ModelGTR*)(model);
    modelGTR->setRateMatrix(initRates);
    modelGTR->setStateFrequency(initStateFreqs);
    rateGammaInvar->setGammaShape(initAlpha);
    rateGammaInvar->setPInvar(initPInvar);
    modelGTR->decomposeRateMatrix();
    site_rate->phylo_tree->clearAllPartialLH();
    return site_rate->phylo_tree->computeLikelihood();
}

double ModelFactory::optimizeParametersOnly(double gradient_epsilon) {
	double logl;
	/* Optimize substitution and heterogeneity rates independently */
	if (!joint_optimize) {
		double model_lh = model->optimizeParameters(gradient_epsilon);
		double rate_lh = site_rate->optimizeParameters(gradient_epsilon);
		if (rate_lh == 0.0)
			logl = model_lh;
		else
			logl = rate_lh;
	} else {
		/* Optimize substitution and heterogeneity rates jointly using BFGS */
		logl = optimizeAllParameters(gradient_epsilon);
	}
	return logl;
}

double ModelFactory::optimizeAllParameters(double gradient_epsilon) {
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
        for (i = model_ndim- model->num_states+2; i <= model_ndim; i++)
            upper_bound[i] = 1.0;
    }

    // setup the bounds for site_rate
    site_rate->setBounds(lower_bound+model_ndim, upper_bound+model_ndim, bound_check+model_ndim);

    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));

    getVariables(variables);
    //if (freq_type == FREQ_ESTIMATE) scaleStateFreq(true);
    model->decomposeRateMatrix();
    site_rate->phylo_tree->clearAllPartialLH();

    score = site_rate->phylo_tree->computeLikelihood();

    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    return score;
}

double ModelFactory::optimizeParametersGammaInvar(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    if (!site_rate->isGammai())
        return optimizeParameters(fixed_len, write_info, logl_epsilon, gradient_epsilon);
        
	double begin_time = getRealTime();
        
    PhyloTree *tree = site_rate->getTree();
	double frac_const = tree->aln->frac_const_sites;
    tree->setCurScore(tree->computeLikelihood());

	/* Back up branch lengths and substitutional rates */
	DoubleVector initBranLens;
	DoubleVector bestLens;
	tree->saveBranchLengths(initBranLens);
    bestLens = initBranLens;
	int numRateEntries = tree->getModel()->getNumRateEntries();
	double *rates = new double[numRateEntries];
	double *bestRates = new double[numRateEntries];
	tree->getModel()->getRateMatrix(rates);
	int numStates = tree->aln->num_states;
	double *state_freqs = new double[numStates];
	tree->getModel()->getStateFrequency(state_freqs);

	/* Best estimates found */
	double *bestStateFreqs =  new double[numStates];
	double bestLogl = -DBL_MAX;
	double bestAlpha = 0.0;
	double bestPInvar = 0.0;

	double testInterval = (frac_const - MIN_PINVAR * 2) / 9;
	double initPInv = MIN_PINVAR;
	double initAlpha = site_rate->getGammaShape();

    if (Params::getInstance().opt_gammai_fast) {
        initPInv = frac_const/2;
        bool stop = false;
        while(!stop) {
            if (write_info) {
                cout << endl;
                cout << "Testing with init. pinv = " << initPInv << " / init. alpha = "  << initAlpha << endl;
            }

            vector<double> estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon, tree, site_rate, rates, state_freqs,
                                                                   initPInv, initAlpha, initBranLens);


            if (write_info) {
                cout << "Est. p_inv: " << estResults[0] << " / Est. gamma shape: " << estResults[1]
                << " / Logl: " << estResults[2] << endl;
            }

            if (estResults[2] > bestLogl) {
                bestLogl = estResults[2];
                bestAlpha = estResults[1];
                bestPInvar = estResults[0];
                bestLens.clear();
                tree->saveBranchLengths(bestLens);
                tree->getModel()->getRateMatrix(bestRates);
                tree->getModel()->getStateFrequency(bestStateFreqs);
                if (estResults[0] < initPInv) {
                    initPInv = estResults[0] - testInterval;
                    if (initPInv < 0.0)
                        initPInv = 0.0;
                } else {
                    initPInv = estResults[0] + testInterval;
                    if (initPInv > frac_const)
                        initPInv = frac_const;
                }
                //cout << "New initPInv = " << initPInv << endl;
            }  else {
                stop = true;
            }
        }
    } else {
        // Now perform testing different initial p_inv values
        while (initPInv <= frac_const) {
            if (write_info) {
                cout << endl;
                cout << "Testing with init. pinv = " << initPInv << " / init. alpha = "  << initAlpha << endl;
            }
            vector<double> estResults; // vector of p_inv, alpha and logl
            if (Params::getInstance().opt_gammai_keep_bran)
                estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon, tree, site_rate, rates, state_freqs,
                                                                          initPInv, initAlpha, bestLens);
            else
                estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon, tree, site_rate, rates, state_freqs,
                                                                      initPInv, initAlpha, initBranLens);
            if (write_info) {
                cout << "Est. p_inv: " << estResults[0] << " / Est. gamma shape: " << estResults[1]
                << " / Logl: " << estResults[2] << endl;
            }

            initPInv = initPInv + testInterval;

            if (estResults[2] > bestLogl) {
                bestLogl = estResults[2];
                bestAlpha = estResults[1];
                bestPInvar = estResults[0];
                bestLens.clear();
                tree->saveBranchLengths(bestLens);
                tree->getModel()->getRateMatrix(bestRates);
                tree->getModel()->getStateFrequency(bestStateFreqs);
            }
        }
    }

    site_rate->setGammaShape(bestAlpha);
    site_rate->setPInvar(bestPInvar);
	((ModelGTR*) tree->getModel())->setRateMatrix(bestRates);
	((ModelGTR*) tree->getModel())->setStateFrequency(bestStateFreqs);
	tree->restoreBranchLengths(bestLens);
	tree->getModel()->decomposeRateMatrix();

	tree->clearAllPartialLH();
	tree->setCurScore(tree->computeLikelihood());
    assert(fabs(tree->getCurScore() - bestLogl) < 1.0);
    if (write_info) {    
        cout << endl;
        cout << "Best p_inv: " << bestPInvar << " / best gamma shape: " << bestAlpha << " / ";
        cout << "Logl: " << tree->getCurScore() << endl;
    }

	delete [] rates;
	delete [] state_freqs;
	delete [] bestRates;
	delete [] bestStateFreqs;
    
	double elapsed_secs = getRealTime() - begin_time;
	if (write_info)
		cout << "Parameters optimization took " << elapsed_secs << " sec" << endl;
    
    // updating global variable is not safe!
//	Params::getInstance().testAlpha = false;
    
    // 2016-03-14: this was missing!
    return tree->getCurScore();
}

vector<double> ModelFactory::optimizeGammaInvWithInitValue(int fixed_len, double logl_epsilon, double gradient_epsilon,
                                                 PhyloTree *tree, RateHeterogeneity *site_rates, double *rates,
                                                 double *state_freqs, double initPInv, double initAlpha,
                                                 DoubleVector &lenvec) {
    tree->restoreBranchLengths(lenvec);
    ((ModelGTR*) tree->getModel())->setRateMatrix(rates);
    ((ModelGTR*) tree->getModel())->setStateFrequency(state_freqs);
    tree->getModel()->decomposeRateMatrix();
    site_rates->setPInvar(initPInv);
    site_rates->setGammaShape(initAlpha);
    tree->clearAllPartialLH();
    optimizeParameters(fixed_len, false, logl_epsilon, gradient_epsilon);

    vector<double> estResults;
    double estPInv = tree->getRate()->getPInvar();
    double estAlpha = tree->getRate()->getGammaShape();
    double logl = tree->getCurScore();
    estResults.push_back(estPInv);
    estResults.push_back(estAlpha);
    estResults.push_back(logl);
    return estResults;
}


double ModelFactory::optimizeParameters(int fixed_len, bool write_info,
                                        double logl_epsilon, double gradient_epsilon) {
	assert(model);
	assert(site_rate);

    double defaultEpsilon = logl_epsilon;

	double begin_time = getRealTime();
	double cur_lh;
	PhyloTree *tree = site_rate->getTree();
	assert(tree);

	stopStoringTransMatrix();
    // modified by Thomas Wong on Sept 11, 15
    // no optimization of branch length in the first round
    cur_lh = tree->computeLikelihood();
    tree->setCurScore(cur_lh);
	if (verbose_mode >= VB_MED || write_info) 
		cout << "1. Initial log-likelihood: " << cur_lh << endl;

	// For UpperBounds -----------
	//cout<<"MLCheck = "<<tree->mlCheck <<endl;
	if(tree->mlCheck == 0){
		tree->mlInitial = cur_lh;
	}
	// ---------------------------


	int i;
	//bool optimize_rate = true;
//	double gradient_epsilon = min(logl_epsilon, 0.01); // epsilon for parameters starts at epsilon for logl
	for (i = 2; i < tree->params->num_param_iterations; i++) {
        double new_lh;

        // changed to opimise edge length first, and then Q,W,R inside the loop by Thomas on Sept 11, 15
		if (fixed_len == BRLEN_OPTIMIZE)
			new_lh = tree->optimizeAllBranches(min(i,3), logl_epsilon);  // loop only 3 times in total (previously in v0.9.6 5 times)
        else if (fixed_len == BRLEN_SCALE) {
            double scaling = 1.0;
            new_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
        }
        new_lh = optimizeParametersOnly(gradient_epsilon);

		if (new_lh == 0.0) {
            if (fixed_len == BRLEN_OPTIMIZE)
                cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
            else if (fixed_len == BRLEN_SCALE) {
                double scaling = 1.0;
                cur_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
            }
			break;
		}
		if (verbose_mode >= VB_MED) {
			model->writeInfo(cout);
			site_rate->writeInfo(cout);
		}
		if (new_lh > cur_lh + logl_epsilon) {
			cur_lh = new_lh;
			if (verbose_mode >= VB_MED || write_info)
				cout << i << ". Current log-likelihood: " << cur_lh << endl;
		} else {
			site_rate->classifyRates(new_lh);
            if (fixed_len == BRLEN_OPTIMIZE)
                cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
            else if (fixed_len == BRLEN_SCALE) {
                double scaling = 1.0;
                cur_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
            }
            break;
		}
	}

	// normalize rates s.t. branch lengths are #subst per site
//    if (Params::getInstance().optimize_alg_gammai != "EM") 
    {
        double mean_rate = site_rate->rescaleRates();
        if (mean_rate != 1.0) {
            if (fixed_len == BRLEN_FIX)
                outError("Fixing branch lengths not supported under specified site rate model");
            tree->scaleLength(mean_rate);
            tree->clearAllPartialLH();
        }
    }
    

    
	if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << cur_lh << endl;

	// For UpperBounds -----------
	if(tree->mlCheck == 0)
		tree->mlFirstOpt = cur_lh;
	// ---------------------------

	if (verbose_mode <= VB_MIN && write_info) {
		model->writeInfo(cout);
		site_rate->writeInfo(cout);
	}
	double elapsed_secs = getRealTime() - begin_time;
	if (write_info)
		cout << "Parameters optimization took " << i-1 << " rounds (" << elapsed_secs << " sec)" << endl;
	startStoringTransMatrix();

	// For UpperBounds -----------
	tree->mlCheck = 1;
	// ---------------------------

	tree->setCurScore(cur_lh);
	return cur_lh;
}

/**
 * @return TRUE if parameters are at the boundary that may cause numerical unstability
 */
bool ModelFactory::isUnstableParameters() {
	if (model->isUnstableParameters()) return true;
	return false;
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

bool ModelFactory::getVariables(double *variables) {
	bool changed = model->getVariables(variables);
	changed |= site_rate->getVariables(variables + model->getNDim());
    return changed;
}


