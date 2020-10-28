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
#include "modelmarkov.h"
#include "modelliemarkov.h"
#include "modeldna.h"
#include "modelprotein.h"
#include "modelbin.h"
#include "modelcodon.h"
#include "modelmorphology.h"
#include "modelpomo.h"
#include "modelset.h"
#include "modelmixture.h"
#include "ratemeyerhaeseler.h"
#include "ratemeyerdiscrete.h"
#include "ratekategory.h"
#include "ratefree.h"
#include "ratefreeinvar.h"
#include "rateheterotachy.h"
#include "rateheterotachyinvar.h"
//#include "ngs.h"
#include <string>
#include "utils/timeutil.h"
#include "nclextra/myreader.h"
#include <sstream>

string::size_type findSubStr(string &name, string sub1, string sub2) {
    string::size_type pos1, pos2;
    for (pos1 = 0; pos1 != string::npos; pos1++) {
        pos1 = name.find(sub1, pos1);
        if (pos1 == string::npos)
            break;
        if (pos1+2 >= name.length() || !isalpha(name[pos1+2])) {
            break;
        }
    }
    
    for (pos2 = 0; pos2 != string::npos; pos2++) {
        pos2 = name.find(sub2, pos2);
        if (pos2 == string::npos)
            break;
        if (pos2+2 >= name.length() ||!isalpha(name[pos2+2]))
            break;
    }
    
    if (pos1 != string::npos && pos2 != string::npos) {
        return min(pos1, pos2);
    } else if (pos1 != string::npos)
        return pos1;
    else
        return pos2;
}

string::size_type posRateHeterotachy(string model_name) {
    return findSubStr(model_name, "+H", "*H");
}

string::size_type posRateFree(string &model_name) {
    return findSubStr(model_name, "+R", "*R");
}

string::size_type posPOMO(string &model_name) {
    return findSubStr(model_name, "+P", "*P");
}

ModelsBlock *readModelsDefinition(Params &params) {

    ModelsBlock *models_block = new ModelsBlock;

    try
    {
        // loading internal model definitions
        stringstream in(builtin_mixmodels_definition);
        ASSERT(in && "stringstream is OK");
        NxsReader nexus;
        nexus.Add(models_block);
        MyToken token(in);
        nexus.Execute(token);
    } catch (...) {
        ASSERT(0 && "predefined mixture models not initialized");
    }

    try
    {
        // loading internal protei model definitions
        stringstream in(builtin_prot_models);
        ASSERT(in && "stringstream is OK");
        NxsReader nexus;
        nexus.Add(models_block);
        MyToken token(in);
        nexus.Execute(token);
    } catch (...) {
        ASSERT(0 && "predefined protein models not initialized");
    }

    if (params.model_def_file) {
        cout << "Reading model definition file " << params.model_def_file << " ... ";
        MyReader nexus(params.model_def_file);
        nexus.Add(models_block);
        MyToken token(nexus.inf);
        nexus.Execute(token);
        int num_model = 0, num_freq = 0;
        for (ModelsBlock::iterator it = models_block->begin(); it != models_block->end(); it++)
            if (it->second.flag & NM_FREQ) num_freq++; else num_model++;
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
    ASC_type = ASC_NONE;
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

ModelFactory::ModelFactory(Params &params, string &model_name, PhyloTree *tree, ModelsBlock *models_block) : CheckpointFactory() {
    store_trans_matrix = params.store_trans_matrix;
    is_storing = false;
    joint_optimize = params.optimize_model_rate_joint;
    fused_mix_rate = false;
    ASC_type = ASC_NONE;
    string model_str = model_name;
    string rate_str;

    try {


    if (model_str == "") {
        if (tree->aln->seq_type == SEQ_DNA) model_str = "HKY";
        else if (tree->aln->seq_type == SEQ_PROTEIN) model_str = "LG";
        else if (tree->aln->seq_type == SEQ_BINARY) model_str = "GTR2";
        else if (tree->aln->seq_type == SEQ_CODON) model_str = "GY";
        else if (tree->aln->seq_type == SEQ_MORPH) model_str = "MK";
        else if (tree->aln->seq_type == SEQ_POMO) model_str = "HKY+P";
        else model_str = "JC";
        if (tree->aln->seq_type != SEQ_POMO && !params.model_joint)
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

    //    nxsmodel = models_block->findModel(model_str);
    //    if (nxsmodel && nxsmodel->description.find_first_of("+*") != string::npos) {
    //        cout << "Model " << model_str << " is alias for " << nxsmodel->description << endl;
    //        model_str = nxsmodel->description;
    //    }

    // Detect PoMo and throw error if sequence type is PoMo but +P is
    // not given.  This makes the model string cleaner and
    // compareable.
    string::size_type p_pos = posPOMO(model_str);
    bool pomo = (p_pos != string::npos);

    if ((p_pos == string::npos) &&
        (tree->aln->seq_type == SEQ_POMO))
        outError("Provided alignment is exclusively used by PoMo but model string does not contain, e.g., \"+P\".");

    // Decompose model string into model_str and rate_str string.
    size_t spec_pos = model_str.find_first_of("{+*");
    if (spec_pos != string::npos) {
        if (model_str[spec_pos] == '{') {
            // Scan for the corresponding '}'.
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
    
    // decompose +F from rate_str
    string freq_str = "";
    while ((spec_pos = rate_str.find("+F")) != string::npos) {
        size_t end_pos = rate_str.find_first_of("+*", spec_pos+1);
        if (end_pos == string::npos) {
            freq_str += rate_str.substr(spec_pos);
            rate_str = rate_str.substr(0, spec_pos);
        } else {
            freq_str += rate_str.substr(spec_pos, end_pos - spec_pos);
            rate_str = rate_str.substr(0, spec_pos) + rate_str.substr(end_pos);
        }
    }

    // set to model_joint if set
    if (Params::getInstance().model_joint) {
        model_str = Params::getInstance().model_joint;
        freq_str = "";
        while ((spec_pos = model_str.find("+F")) != string::npos) {
            size_t end_pos = model_str.find_first_of("+*", spec_pos+1);
            if (end_pos == string::npos) {
                freq_str += model_str.substr(spec_pos);
                model_str = model_str.substr(0, spec_pos);
            } else {
                freq_str += model_str.substr(spec_pos, end_pos - spec_pos);
                model_str = model_str.substr(0, spec_pos) + model_str.substr(end_pos);
            }
        }
    }
        
    // move error model +E from rate_str to model_str
    //string seqerr_str = "";
    while ((spec_pos = rate_str.find("+E")) != string::npos) {
        size_t end_pos = rate_str.find_first_of("+*", spec_pos+1);
        if (end_pos == string::npos) {
            model_str += rate_str.substr(spec_pos);
            rate_str = rate_str.substr(0, spec_pos);
        } else {
            model_str += rate_str.substr(spec_pos, end_pos - spec_pos);
            rate_str = rate_str.substr(0, spec_pos) + rate_str.substr(end_pos);
        }
    }

    // PoMo; +NXX and +W or +S because those flags are handled when
    // reading in the data.  Set PoMo parameters (heterozygosity).
    size_t n_pos_start = rate_str.find("+N");
    size_t n_pos_end   = rate_str.find_first_of("+", n_pos_start+1);
    if (n_pos_start != string::npos) {
        if (!pomo)
            outError("Virtual population size can only be set with PoMo.");
        if (n_pos_end != string::npos)
            rate_str = rate_str.substr(0, n_pos_start)
                + rate_str.substr(n_pos_end);
        else
            rate_str = rate_str.substr(0, n_pos_start);
    }

    size_t wb_pos = rate_str.find("+WB");
    if (wb_pos != string::npos) {
      if (!pomo)
        outError("Weighted binomial sampling can only be used with PoMo.");
      rate_str = rate_str.substr(0, wb_pos)
        + rate_str.substr(wb_pos+3);
    }
    size_t wh_pos = rate_str.find("+WH");
    if (wh_pos != string::npos) {
        if (!pomo)
            outError("Weighted hypergeometric sampling can only be used with PoMo.");
        rate_str = rate_str.substr(0, wh_pos)
            + rate_str.substr(wh_pos+3);
    }
    size_t s_pos = rate_str.find("+S");
    if ( s_pos != string::npos) {
        if (!pomo)
            outError("Binomial sampling can only be used with PoMo.");
        rate_str = rate_str.substr(0, s_pos)
            + rate_str.substr(s_pos+2);
    }

    // In case of PoMo, check that only supported flags are given.
    if (pomo) {
        if (rate_str.find("+ASC") != string::npos)
            // TODO DS: This is an important feature, because then,
            // PoMo can be applied to SNP data only.
            outError("PoMo does not yet support ascertainment bias correction (+ASC).");
        if (posRateFree(rate_str) != string::npos)
            outError("PoMo does not yet support free rate models (+R).");
        if (rate_str.find("+FMIX") != string::npos)
            outError("PoMo does not yet support frequency mixture models (+FMIX).");
        if (posRateHeterotachy(rate_str) != string::npos)
            outError("PoMo does not yet support heterotachy models (+H).");
    }

    // PoMo. The +P{}, +GXX and +I flags are interpreted during model creation.
    // This is necessary for compatibility with mixture models. If there is no
    // mixture model, move +P{}, +GXX and +I flags to model string. For mixture
    // models, the heterozygosity can be set separately for each model and the
    // +P{}, +GXX and +I flags should already be inside the model definition.
    if (model_str.substr(0, 3) != "MIX" && pomo) {
      // +P{} flag.
      p_pos = posPOMO(rate_str);
      if (p_pos != string::npos) {
        if (rate_str[p_pos+2] == '{') {
          string::size_type close_bracket = rate_str.find("}");
          if (close_bracket == string::npos)
            outError("No closing bracket in PoMo parameters.");
          else {
            string pomo_heterozygosity = rate_str.substr(p_pos+3,close_bracket-p_pos-3);
            rate_str = rate_str.substr(0, p_pos) + rate_str.substr(close_bracket+1);
            model_str += "+P{" + pomo_heterozygosity + "}";
          }
        }
        else {
          rate_str = rate_str.substr(0, p_pos) + rate_str.substr(p_pos + 2);
          model_str += "+P";
        }
      }

      // +G flag.
      size_t pomo_rate_start_pos;
      if ((pomo_rate_start_pos = rate_str.find("+G")) != string::npos) {
        string pomo_rate_str = "";
        size_t pomo_rate_end_pos = rate_str.find_first_of("+*", pomo_rate_start_pos+1);
        if (pomo_rate_end_pos == string::npos) {
          pomo_rate_str = rate_str.substr(pomo_rate_start_pos, rate_str.length() - pomo_rate_start_pos);
          rate_str = rate_str.substr(0, pomo_rate_start_pos);
          model_str += pomo_rate_str;
        } else {
          pomo_rate_str = rate_str.substr(pomo_rate_start_pos, pomo_rate_end_pos - pomo_rate_start_pos);
          rate_str = rate_str.substr(0, pomo_rate_start_pos) + rate_str.substr(pomo_rate_end_pos);
          model_str += pomo_rate_str;
        }
      }

      // // +I flag.
      // size_t pomo_irate_start_pos;
      // if ((pomo_irate_start_pos = rate_str.find("+I")) != string::npos) {
      //   string pomo_irate_str = "";
      //   size_t pomo_irate_end_pos = rate_str.find_first_of("+*", pomo_irate_start_pos+1);
      //   if (pomo_irate_end_pos == string::npos) {
      //     pomo_irate_str = rate_str.substr(pomo_irate_start_pos, rate_str.length() - pomo_irate_start_pos);
      //     rate_str = rate_str.substr(0, pomo_irate_start_pos);
      //     model_str += pomo_irate_str;
      //   } else {
      //     pomo_irate_str = rate_str.substr(pomo_irate_start_pos, pomo_irate_end_pos - pomo_irate_start_pos);
      //     rate_str = rate_str.substr(0, pomo_irate_start_pos) + rate_str.substr(pomo_irate_end_pos);
      //     model_str += pomo_irate_str;
      //   }
    }

    //    nxsmodel = models_block->findModel(model_str);
    //    if (nxsmodel && nxsmodel->description.find("MIX") != string::npos) {
    //        cout << "Model " << model_str << " is alias for " << nxsmodel->description << endl;
    //        model_str = nxsmodel->description;
    //    }

    /******************** initialize state frequency ****************************/

    StateFreqType freq_type = params.freq_type;

    if (freq_type == FREQ_UNKNOWN) {
        switch (tree->aln->seq_type) {
        case SEQ_BINARY: freq_type = FREQ_ESTIMATE; break; // default for binary: optimized frequencies
        case SEQ_PROTEIN: break; // let ModelProtein decide by itself
        case SEQ_MORPH: freq_type = FREQ_EQUAL; break;
        case SEQ_CODON: freq_type = FREQ_UNKNOWN; break;
            break;
        default: freq_type = FREQ_EMPIRICAL; break; // default for DNA, PoMo and others: counted frequencies from alignment
        }
    }

    // first handle mixture frequency
    string::size_type posfreq = freq_str.find("+FMIX");
    string freq_params;
    size_t close_bracket;

    if (posfreq != string::npos) {
        string fmix_str;
        size_t last_pos = freq_str.find_first_of("+*", posfreq+1);

        if (last_pos == string::npos) {
            fmix_str = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        } else {
            fmix_str = freq_str.substr(posfreq, last_pos-posfreq);
            freq_str = freq_str.substr(0, posfreq) + freq_str.substr(last_pos);
        }

        if (fmix_str[5] != OPEN_BRACKET)
            outError("Mixture-frequency must start with +FMIX{");
        close_bracket = fmix_str.find(CLOSE_BRACKET);
        if (close_bracket == string::npos)
            outError("Close bracket not found in ", fmix_str);
        if (close_bracket != fmix_str.length()-1)
            outError("Wrong close bracket position ", fmix_str);
        freq_type = FREQ_MIXTURE;
        freq_params = fmix_str.substr(6, close_bracket-6);
    }

    // then normal frequency
    if (freq_str.find("+FO") != string::npos)
        posfreq = freq_str.find("+FO");
    else if (freq_str.find("+Fo") != string::npos)
        posfreq = freq_str.find("+Fo");
    else
        posfreq = freq_str.find("+F");

    bool optimize_mixmodel_weight = params.optimize_mixmodel_weight;

    if (posfreq != string::npos) {
        string fstr;
        size_t last_pos = freq_str.find_first_of("+*", posfreq+1);
        if (last_pos == string::npos) {
            fstr = freq_str.substr(posfreq);
            freq_str = freq_str.substr(0, posfreq);
        } else {
            fstr = freq_str.substr(posfreq, last_pos-posfreq);
            freq_str = freq_str.substr(0, posfreq) + freq_str.substr(last_pos);
        }

        if (fstr.length() > 2 && fstr[2] == OPEN_BRACKET) {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with user-defined frequency is not allowed");
            close_bracket = fstr.find(CLOSE_BRACKET);
            if (close_bracket == string::npos)
                outError("Close bracket not found in ", fstr);
            if (close_bracket != fstr.length()-1)
                outError("Wrong close bracket position ", fstr);
            freq_type = FREQ_USER_DEFINED;
            freq_params = fstr.substr(3, close_bracket-3);
        } else if (fstr == "+FC" || fstr == "+Fc" || fstr == "+F") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "empirical," + freq_params;
                optimize_mixmodel_weight = true;
            } else
                freq_type = FREQ_EMPIRICAL;
        } else if (fstr == "+FU" || fstr == "+Fu") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with user-defined frequency is not allowed");
            else
                freq_type = FREQ_USER_DEFINED;
        } else if (fstr == "+FQ" || fstr == "+Fq") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with equal frequency is not allowed");
            else
            freq_type = FREQ_EQUAL;
        } else if (fstr == "+FO" || fstr == "+Fo") {
            if (freq_type == FREQ_MIXTURE) {
                freq_params = "optimize," + freq_params;
                optimize_mixmodel_weight = true;
            } else
                freq_type = FREQ_ESTIMATE;
    } else if (fstr == "+F1x4" || fstr == "+F1X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_1x4;
        } else if (fstr == "+F3x4" || fstr == "+F3X4") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4;
        } else if (fstr == "+F3x4C" || fstr == "+F3x4c" || fstr == "+F3X4C" || fstr == "+F3X4c") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_CODON_3x4C;
        } else if (fstr == "+FRY") {
        // MDW to Minh: I don't know how these should interact with FREQ_MIXTURE,
        // so as nearly everything else treats it as an error, I do too.
        // BQM answer: that's fine
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_RY;
        } else if (fstr == "+FWS") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_WS;
        } else if (fstr == "+FMK") {
            if (freq_type == FREQ_MIXTURE)
                outError("Mixture frequency with " + fstr + " is not allowed");
            else
                freq_type = FREQ_DNA_MK;
        } else {
            // might be "+F####" where # are digits
            try {
                freq_type = parseStateFreqDigits(fstr.substr(2)); // throws an error if not in +F#### format
            } catch (...) {
                outError("Unknown state frequency type ",fstr);
            }
        }
//          model_str = model_str.substr(0, posfreq);
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
            model = new ModelMixture(model_name, model_str, model_list, models_block, freq_type, freq_params, tree, optimize_mixmodel_weight);
        } else {
            //            string model_desc;
            //            NxsModel *nxsmodel = models_block->findModel(model_str);
            //            if (nxsmodel) model_desc = nxsmodel->description;
            model = createModel(model_str, models_block, freq_type, freq_params, tree);
        }
//        fused_mix_rate &= model->isMixture() && site_rate->getNRate() > 1;
    } else {
        // site-specific model
        if (model_str == "JC" || model_str == "POISSON")
            outError("JC is not suitable for site-specific model");
        model = new ModelSet(model_str.c_str(), tree);
        ModelSet *models = (ModelSet*)model; // assign pointer for convenience
        models->init((params.freq_type != FREQ_UNKNOWN) ? params.freq_type : FREQ_EMPIRICAL);
        models->pattern_model_map.resize(tree->aln->getNPattern(), -1);
        for (size_t i = 0; i < tree->aln->getNSite(); ++i) {
            models->pattern_model_map[tree->aln->getPatternID(i)] = tree->aln->site_model[i];
            //cout << "site " << i << " ptn " << tree->aln->getPatternID(i) << " -> model " << site_model[i] << endl;
        }
        double *state_freq = new double[model->num_states];
        double *rates = new double[model->getNumRateEntries()];
        for (size_t i = 0; i < tree->aln->site_state_freq.size(); ++i) {
            ModelMarkov *modeli;
            if (i == 0) {
                modeli = (ModelMarkov*)createModel(model_str, models_block, (params.freq_type != FREQ_UNKNOWN) ? params.freq_type : FREQ_EMPIRICAL, "", tree);
                modeli->getStateFrequency(state_freq);
                modeli->getRateMatrix(rates);
            } else {
                modeli = (ModelMarkov*)createModel(model_str, models_block, FREQ_EQUAL, "", tree);
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

        models->joinEigenMemory();
        models->decomposeRateMatrix();

        // delete information of the old alignment
//        tree->aln->ordered_pattern.clear();
//        tree->deleteAllPartialLh();
    }

//    if (model->isMixture())
//        cout << "Mixture model with " << model->getNMixtures() << " components!" << endl;

    /******************** initialize ascertainment bias correction model ****************************/

    string::size_type posasc;

    if ((posasc = rate_str.find("+ASC_INF")) != string::npos) {
        // ascertainment bias correction
        ASC_type = ASC_INFORMATIVE;
        tree->aln->getUnobservedConstPatterns(ASC_type, unobserved_ptns);
        
        // rebuild the seq_states to contain states of unobserved constant patterns
        //tree->aln->buildSeqStates(model->seq_states, true);
        if (tree->aln->num_informative_sites != tree->getAlnNSite()) {
            if (!params.partition_file) {
                string infsites_file = ((string)params.out_prefix + ".infsites.phy");
                tree->aln->printAlignment(params.aln_output_format, infsites_file.c_str(), false, NULL, EXCLUDE_UNINF);
                cerr << "For your convenience alignment with parsimony-informative sites printed to " << infsites_file << endl;
            }
            outError("Invalid use of +ASC_INF because of " + convertIntToString(tree->getAlnNSite() - tree->aln->num_informative_sites) +
                     " parsimony-uninformative sites in the alignment");
        }
        if (verbose_mode >= VB_MED)
            cout << "Ascertainment bias correction: " << unobserved_ptns.size() << " unobservable uninformative patterns"<< endl;
        rate_str = rate_str.substr(0, posasc) + rate_str.substr(posasc+8);
    } else if ((posasc = rate_str.find("+ASC_MIS")) != string::npos) {
        // initialize Holder's ascertainment bias correction model
        ASC_type = ASC_VARIANT_MISSING;
        tree->aln->getUnobservedConstPatterns(ASC_type, unobserved_ptns);
        // rebuild the seq_states to contain states of unobserved constant patterns
        //tree->aln->buildSeqStates(model->seq_states, true);
        if (tree->aln->frac_invariant_sites > 0) {
            if (!params.partition_file) {
                string varsites_file = ((string)params.out_prefix + ".varsites.phy");
                tree->aln->printAlignment(params.aln_output_format, varsites_file.c_str(), false, NULL, EXCLUDE_INVAR);
                cerr << "For your convenience alignment with variable sites printed to " << varsites_file << endl;
            }
            outError("Invalid use of +ASC_MIS because of " + convertIntToString(tree->aln->frac_invariant_sites*tree->aln->getNSite()) +
                     " invariant sites in the alignment");
        }
        if (verbose_mode >= VB_MED)
            cout << "Holder's ascertainment bias correction: " << unobserved_ptns.size() << " unobservable constant patterns" << endl;
        rate_str = rate_str.substr(0, posasc) + rate_str.substr(posasc+8);
    } else if ((posasc = rate_str.find("+ASC")) != string::npos) {
        // ascertainment bias correction
        ASC_type = ASC_VARIANT;
        tree->aln->getUnobservedConstPatterns(ASC_type, unobserved_ptns);
        
        // delete rarely observed state
        for (int i = unobserved_ptns.size()-1; i >= 0; i--)
            if (model->state_freq[(int)unobserved_ptns[i][0]] < 1e-8)
                unobserved_ptns.erase(unobserved_ptns.begin() + i);
                
        // rebuild the seq_states to contain states of unobserved constant patterns
        //tree->aln->buildSeqStates(model->seq_states, true);
//        if (unobserved_ptns.size() <= 0)
//            outError("Invalid use of +ASC because all constant patterns are observed in the alignment");
        if (tree->aln->frac_invariant_sites > 0) {
//            cerr << tree->aln->frac_invariant_sites*tree->aln->getNSite() << " invariant sites observed in the alignment" << endl;
//            for (Alignment::iterator pit = tree->aln->begin(); pit != tree->aln->end(); pit++)
//                if (pit->isInvariant()) {
//                    string pat_str = "";
//                    for (Pattern::iterator it = pit->begin(); it != pit->end(); it++)
//                        pat_str += tree->aln->convertStateBackStr(*it);
//                    cerr << pat_str << " is invariant site pattern" << endl;
//                }
            if (!params.partition_file) {
                string varsites_file = ((string)params.out_prefix + ".varsites.phy");
                tree->aln->printAlignment(params.aln_output_format, varsites_file.c_str(), false, NULL, EXCLUDE_INVAR);
                cerr << "For your convenience alignment with variable sites printed to " << varsites_file << endl;
            }
            outError("Invalid use of +ASC because of " + convertIntToString(tree->aln->frac_invariant_sites*tree->aln->getNSite()) +
                " invariant sites in the alignment");
        }
        if (verbose_mode >= VB_MED)
            cout << "Ascertainment bias correction: " << unobserved_ptns.size() << " unobservable constant patterns"<< endl;
		rate_str = rate_str.substr(0, posasc) + rate_str.substr(posasc+4);
    } else {
        //tree->aln->buildSeqStates(model->seq_states, false);
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

    string::size_type posR = rate_str.find("+R"); // FreeRate model
    string::size_type posR2 = rate_str.find("*R"); // FreeRate model

    if (posG != string::npos && (posR != string::npos || posR2 != string::npos)) {
        outWarning("Both Gamma and FreeRate models were specified, continue with FreeRate model");
        posG = string::npos;
        fused_mix_rate = false;
    }

    if (posR != string::npos && posR2 != string::npos) {
        cout << "NOTE: both +R and *R were specified, continue with "
            << ((posR < posR2)? rate_str.substr(posR,2) : rate_str.substr(posR2,2)) << endl;
    }

    if (posR2 != string::npos && posR2 < posR) {
        posR = posR2;
        fused_mix_rate = true;
    }

    string::size_type posH = rate_str.find("+H"); // heterotachy model
    string::size_type posH2 = rate_str.find("*H"); // heterotachy model

    if (posG != string::npos && (posH != string::npos || posH2 != string::npos)) {
        outWarning("Both Gamma and heterotachy models were specified, continue with heterotachy model");
        posG = string::npos;
        fused_mix_rate = false;
    }

    if (posR != string::npos && (posH != string::npos || posH2 != string::npos)) {
        outWarning("Both FreeRate and heterotachy models were specified, continue with heterotachy model");
        posR = string::npos;
        fused_mix_rate = false;
    }

    if (posH != string::npos && posH2 != string::npos) {
        cout << "NOTE: both +H and *H were specified, continue with "
            << ((posH < posH2)? rate_str.substr(posH,2) : rate_str.substr(posH2,2)) << endl;
    }
    if (posH2 != string::npos && posH2 < posH) {
        posH = posH2;
        fused_mix_rate = true;
    }

    string::size_type posX;
    /* create site-rate heterogeneity */
    int num_rate_cats = params.num_rate_cats;
    if (fused_mix_rate && model->isMixture()) num_rate_cats = model->getNMixtures();
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
        } else if (rate_str.length() > posI+2 && rate_str[posI+2] != '+' && rate_str[posI+2] != '*')
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
//            if (gamma_shape < MIN_GAMMA_SHAPE || gamma_shape > MAX_GAMMA_SHAPE) {
//                stringstream str;
//                str << "Gamma shape parameter " << gamma_shape << "out of range ["
//                        << MIN_GAMMA_SHAPE << ',' << MAX_GAMMA_SHAPE << "]" << endl;
//                outError(str.str());
//            }
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

    string heterotachy_params = "";
    if (posH != string::npos) {
        // Heterotachy model
        int end_pos = 0;
        if (rate_str.length() > posH+2 && isdigit(rate_str[posH+2])) {
            num_rate_cats = convert_int(rate_str.substr(posH+2).c_str(), end_pos);
                if (num_rate_cats < 1) outError("Wrong number of rate categories");
        } else {
            if (!model->isMixture() || !fused_mix_rate)
                outError("Please specify number of heterotachy classes (e.g., +H2)");
        }
        if (rate_str.length() > posH+2+end_pos && rate_str[posH+2+end_pos] == OPEN_BRACKET) {
            close_bracket = rate_str.find(CLOSE_BRACKET, posH);
            if (close_bracket == string::npos)
                outError("Close bracket not found in ", rate_str);
            heterotachy_params = rate_str.substr(posH+3+end_pos, close_bracket-posH-3-end_pos).c_str();
        } else if (rate_str.length() > posH+2+end_pos && rate_str[posH+2+end_pos] != '+')
            outError("Wrong model name ", rate_str);
    }


    if (rate_str.find('+') != string::npos || rate_str.find('*') != string::npos) {
        //string rate_str = model_str.substr(pos);
        if (posI != string::npos && posH != string::npos) {
            site_rate = new RateHeterotachyInvar(num_rate_cats, heterotachy_params, p_invar_sites, tree);
        } else if (posH != string::npos) {
            site_rate = new RateHeterotachy(num_rate_cats, heterotachy_params, tree);
        } else if (posI != string::npos && posG != string::npos) {
            site_rate = new RateGammaInvar(num_rate_cats, gamma_shape, params.gamma_median,
                    p_invar_sites, params.optimize_alg_gammai, tree, false);
        } else if (posI != string::npos && posR != string::npos) {
            site_rate = new RateFreeInvar(num_rate_cats, gamma_shape, freerate_params, !fused_mix_rate, p_invar_sites, params.optimize_alg_freerate, tree);
        } else if (posI != string::npos) {
            site_rate = new RateInvar(p_invar_sites, tree);
        } else if (posG != string::npos) {
            site_rate = new RateGamma(num_rate_cats, gamma_shape, params.gamma_median, tree);
        } else if (posR != string::npos) {
            site_rate = new RateFree(num_rate_cats, gamma_shape, freerate_params, !fused_mix_rate, params.optimize_alg_freerate, tree);
//        } else if ((posX = rate_str.find("+M")) != string::npos) {
//            tree->setLikelihoodKernel(LK_NORMAL);
//            params.rate_mh_type = true;
//            if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
//                num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
//                if (num_rate_cats < 0) outError("Wrong number of rate categories");
//            } else num_rate_cats = -1;
//            if (num_rate_cats >= 0)
//                site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type,
//                    params.rate_file, tree, params.rate_mh_type);
//            else
//                site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
//            site_rate->setTree(tree);
//        } else if ((posX = rate_str.find("+D")) != string::npos) {
//            tree->setLikelihoodKernel(LK_NORMAL);
//            params.rate_mh_type = false;
//            if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
//                num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
//                if (num_rate_cats < 0) outError("Wrong number of rate categories");
//            } else num_rate_cats = -1;
//            if (num_rate_cats >= 0)
//                site_rate = new RateMeyerDiscrete(num_rate_cats, params.mcat_type,
//                    params.rate_file, tree, params.rate_mh_type);
//            else
//                site_rate = new RateMeyerHaeseler(params.rate_file, tree, params.rate_mh_type);
//            site_rate->setTree(tree);
//        } else if ((posX = rate_str.find("+NGS")) != string::npos) {
//            tree->setLikelihoodKernel(LK_NORMAL);
//            if (rate_str.length() > posX+4 && isdigit(rate_str[posX+4])) {
//                num_rate_cats = convert_int(rate_str.substr(posX+4).c_str());
//                if (num_rate_cats < 0) outError("Wrong number of rate categories");
//            } else num_rate_cats = -1;
//            site_rate = new NGSRateCat(tree, num_rate_cats);
//            site_rate->setTree(tree);
//        } else if ((posX = rate_str.find("+NGS")) != string::npos) {
//            tree->setLikelihoodKernel(LK_NORMAL);
//            if (rate_str.length() > posX+4 && isdigit(rate_str[posX+4])) {
//                num_rate_cats = convert_int(rate_str.substr(posX+4).c_str());
//                if (num_rate_cats < 0) outError("Wrong number of rate categories");
//            } else num_rate_cats = -1;
//            site_rate = new NGSRate(tree);
//            site_rate->setTree(tree);
        } else if ((posX = rate_str.find("+K")) != string::npos) {
            if (rate_str.length() > posX+2 && isdigit(rate_str[posX+2])) {
                num_rate_cats = convert_int(rate_str.substr(posX+2).c_str());
                if (num_rate_cats < 1) outError("Wrong number of rate categories");
            }
            site_rate = new RateKategory(num_rate_cats, tree);
        } else
            outError("Invalid rate heterogeneity type");
//        if (model_str.find('+') != string::npos)
//            model_str = model_str.substr(0, model_str.find('+'));
//        else
//            model_str = model_str.substr(0, model_str.find('*'));
    } else {
        site_rate = new RateHeterogeneity();
        site_rate->setTree(tree);
    }

    if (fused_mix_rate) {
        if (!model->isMixture()) {
            if (verbose_mode >= VB_MED)
                cout << endl << "NOTE: Using mixture model with unlinked " << model_str << " parameters" << endl;
            string model_list = model_str;
            delete model;
            for (int i = 1; i < site_rate->getNRate(); i++)
                model_list += "," + model_str;
            model = new ModelMixture(model_name, model_str, model_list, models_block, freq_type, freq_params, tree, optimize_mixmodel_weight);
        }
        if (model->getNMixtures() != site_rate->getNRate()) {
            outError("Mixture model and site rate model do not have the same number of categories");
        }
        //if (!tree->isMixlen()) {
        // reset mixture model
        model->setFixMixtureWeight(true);
        int mix, nmix = model->getNMixtures();
        for (mix = 0; mix < nmix; mix++) {
            ((ModelMarkov*)model->getMixtureClass(mix))->total_num_subst = 1.0;
            model->setMixtureWeight(mix, 1.0);
        }
        model->decomposeRateMatrix();
//        } else {
//            site_rate->setFixParams(1);
//            int c, ncat = site_rate->getNRate();
//            for (c = 0; c < ncat; c++)
//                site_rate->setProp(c, 1.0);
//        }
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

void ModelFactory::startCheckpoint() {
    checkpoint->startStruct("ModelFactory");
}

void ModelFactory::saveCheckpoint() {
    model->saveCheckpoint();
    site_rate->saveCheckpoint();
    startCheckpoint();
//    CKP_SAVE(fused_mix_rate);
//    CKP_SAVE(unobserved_ptns);
//    CKP_SAVE(joint_optimize);
    endCheckpoint();
    CheckpointFactory::saveCheckpoint();
}

void ModelFactory::restoreCheckpoint() {
    model->restoreCheckpoint();
    site_rate->restoreCheckpoint();
    startCheckpoint();
//    CKP_RESTORE(fused_mix_rate);
//    CKP_RESTORE(unobserved_ptns);
//    CKP_RESTORE(joint_optimize);
    endCheckpoint();
}

int ModelFactory::getNParameters(int brlen_type) {
    int df = model->getNDim() + model->getNDimFreq() + site_rate->getNDim() +
        site_rate->getTree()->getNBranchParameters(brlen_type);

    return df;
}
/*
double ModelFactory::initGTRGammaIParameters(RateHeterogeneity *rate, ModelSubst *model, double initAlpha,
                                           double initPInvar, double *initRates, double *initStateFreqs)  {

    RateHeterogeneity* rateGammaInvar = rate;
    ModelMarkov* modelGTR = (ModelMarkov*)(model);
    modelGTR->setRateMatrix(initRates);
    modelGTR->setStateFrequency(initStateFreqs);
    rateGammaInvar->setGammaShape(initAlpha);
    rateGammaInvar->setPInvar(initPInvar);
    modelGTR->decomposeRateMatrix();
    site_rate->phylo_tree->clearAllPartialLH();
    return site_rate->phylo_tree->computeLikelihood();
}
*/

double ModelFactory::optimizeParametersOnly(int num_steps, double gradient_epsilon, double cur_logl) {
    double logl;
    /* Optimize substitution and heterogeneity rates independently */
    if (!joint_optimize) {
        // more steps for fused mix rate model
        int steps;
        if (false && fused_mix_rate && model->getNDim() > 0 && site_rate->getNDim() > 0) {
            model->setOptimizeSteps(1);
            site_rate->setOptimizeSteps(1);
            steps = max(model->getNDim()+site_rate->getNDim(), num_steps) * 3;
        } else {
            steps = 1;
        }
        double prev_logl = cur_logl;
        for (int step = 0; step < steps; step++) {
            double model_lh = 0.0;
            // only optimized if model is not linked
            model_lh = model->optimizeParameters(gradient_epsilon);

            double rate_lh = site_rate->optimizeParameters(gradient_epsilon);

            if (rate_lh == 0.0)
                logl = model_lh;
            else
                logl = rate_lh;
            if (logl <= prev_logl + gradient_epsilon)
                break;
            prev_logl = logl;
        }
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
    if (!site_rate->isGammai() || site_rate->isFixPInvar() || site_rate->isFixGammaShape() || site_rate->getTree()->aln->frac_const_sites == 0.0 || model->isMixture())
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
//    int numRateEntries = tree->getModel()->getNumRateEntries();
    Checkpoint *model_ckp = new Checkpoint;
    Checkpoint *best_ckp = new Checkpoint;
    Checkpoint *saved_ckp = model->getCheckpoint();
    *model_ckp = *saved_ckp;
//    double *rates = new double[numRateEntries];
//    double *bestRates = new double[numRateEntries];
//    tree->getModel()->getRateMatrix(rates);
//    int numStates = tree->aln->num_states;
//    double *state_freqs = new double[numStates];
//    tree->getModel()->getStateFrequency(state_freqs);

    /* Best estimates found */
//    double *bestStateFreqs =  new double[numStates];
    double bestLogl = -DBL_MAX;
    double bestAlpha = 0.0;
    double bestPInvar = 0.0;

    size_t numberOfStartValues = 10; //Was hardcoded
    
    double testInterval = (frac_const - MIN_PINVAR * 2.0) / ((double)numberOfStartValues - 1.0);
    double initPInv = MIN_PINVAR;
    double initAlpha = site_rate->getGammaShape();

    if (Params::getInstance().opt_gammai_fast) {
        initPInv = frac_const/2;
        bool stop = false;
        while(!stop) {
            if (write_info) {
                tree->hideProgress();
                cout << endl;
                cout << "Testing with init. pinv = " << initPInv
                    << " / init. alpha = "  << initAlpha << endl;
                tree->showProgress();
            }
            vector<double> estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon,
                                                                   initPInv, initAlpha, initBranLens, model_ckp);

            if (write_info) {
                tree->hideProgress();
                cout << "Est. p_inv: " << estResults[0] << " / Est. gamma shape: " << estResults[1]
                << " / Logl: " << estResults[2] << endl;
                tree->showProgress();
            }
            if (estResults[2] > bestLogl) {
                bestLogl = estResults[2];
                bestAlpha = estResults[1];
                bestPInvar = estResults[0];
                bestLens.clear();
                tree->saveBranchLengths(bestLens);
                model->setCheckpoint(best_ckp);
                model->saveCheckpoint();
                model->setCheckpoint(saved_ckp);
//                *best_ckp = *saved_ckp;

//                tree->getModel()->getRateMatrix(bestRates);
//                tree->getModel()->getStateFrequency(bestStateFreqs);
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
        std::stringstream whatIAmDoing;
        whatIAmDoing << "Thoroughly optimizing +I+G parameters from "
            << numberOfStartValues << " start values";
        if (write_info) {
            if (!progress_display::getProgressDisplay()) {
                cout << whatIAmDoing.str() << "..." << endl;
            }
        }
        tree->initProgress(numberOfStartValues,
                           whatIAmDoing.str(),
                           "tried", "start value");
        while (initPInv <= frac_const) {
            vector<double> estResults; // vector of p_inv, alpha and logl
            if (Params::getInstance().opt_gammai_keep_bran) {
                estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon,
                    initPInv, initAlpha, bestLens, model_ckp);
            }
            else {
                estResults = optimizeGammaInvWithInitValue(fixed_len, logl_epsilon, gradient_epsilon,
                    initPInv, initAlpha, initBranLens, model_ckp);
            }
            if (write_info) {
                tree->hideProgress();
                cout << "Init pinv, alpha: " << initPInv << ", "  << initAlpha
                     << " / Estimate: " << estResults[0] << ", " << estResults[1]
                     << " / LogL: " << estResults[2] << endl;
                tree->showProgress();
            }
            initPInv = initPInv + testInterval;
            if (estResults[2] > bestLogl) {
                bestLogl = estResults[2];
                bestAlpha = estResults[1];
                bestPInvar = estResults[0];
                bestLens.clear();
                tree->saveBranchLengths(bestLens);
                model->setCheckpoint(best_ckp);
                model->saveCheckpoint();
                model->setCheckpoint(saved_ckp);
//                *best_ckp = *saved_ckp;

//                tree->getModel()->getRateMatrix(bestRates);
//                tree->getModel()->getStateFrequency(bestStateFreqs);
            }
            tree->trackProgress(1.0);
        }
        tree->doneProgress();
    }
    site_rate->setGammaShape(bestAlpha);
    site_rate->setPInvar(bestPInvar);

    // -- Mon Apr 17 21:12:14 BST 2017
    // DONE Minh, merged correctly
    model->setCheckpoint(best_ckp);
    model->restoreCheckpoint();
    model->setCheckpoint(saved_ckp);
    // ((ModelGTR*) tree->getModel())->setRateMatrix(bestRates);
    // ((ModelGTR*) tree->getModel())->setStateFrequency(bestStateFreqs);
    // --

    tree->restoreBranchLengths(bestLens);
    // tree->getModel()->decomposeRateMatrix();

    tree->clearAllPartialLH();
    tree->setCurScore(tree->computeLikelihood());
    if (write_info) {
        tree->hideProgress();
        cout << "Optimal pinv,alpha: " << bestPInvar << ", " << bestAlpha << " / ";
        cout << "LogL: " << tree->getCurScore() << endl << endl;
        tree->showProgress();
    }
    if (!tree->params->ignore_any_errors) {
        ASSERT(fabs(tree->getCurScore() - bestLogl) < 1.0);
    }

//    delete [] rates;
//    delete [] state_freqs;
//    delete [] bestRates;
//    delete [] bestStateFreqs;

    delete model_ckp;
    delete best_ckp;

    double elapsed_secs = getRealTime() - begin_time;
    if (write_info) {
        tree->hideProgress();
        cout << "Parameters optimization took " << elapsed_secs << " sec" << endl;
        tree->showProgress();
    }
    // updating global variable is not safe!
//    Params::getInstance().testAlpha = false;

    // 2016-03-14: this was missing!
    return tree->getCurScore();
}

vector<double> ModelFactory::optimizeGammaInvWithInitValue(int fixed_len, double logl_epsilon, double gradient_epsilon,
                                                 double initPInv, double initAlpha,
                                                 DoubleVector &lenvec, Checkpoint *model_ckp) {
    PhyloTree *tree = site_rate->getTree();
    tree->restoreBranchLengths(lenvec);

    // -- Mon Apr 17 21:12:24 BST 2017
    // DONE Minh: merged correctly
    Checkpoint *saved_ckp = model->getCheckpoint();
    model->setCheckpoint(model_ckp);
    model->restoreCheckpoint();
    model->setCheckpoint(saved_ckp);
    site_rate->setPInvar(initPInv);
    site_rate->setGammaShape(initAlpha);
    // --

    tree->clearAllPartialLH();
    optimizeParameters(fixed_len, false, logl_epsilon, gradient_epsilon);

    vector<double> estResults;
    double estPInv = site_rate->getPInvar();
    double estAlpha = site_rate->getGammaShape();
    double logl = tree->getCurScore();
    estResults.push_back(estPInv);
    estResults.push_back(estAlpha);
    estResults.push_back(logl);
    return estResults;
}


double ModelFactory::optimizeParameters(int fixed_len, bool write_info,
                                        double logl_epsilon, double gradient_epsilon) {
    ASSERT(model);
    ASSERT(site_rate);

    double begin_time = getRealTime();
    PhyloTree *tree = site_rate->getTree();
    ASSERT(tree);

    double estimatedIterations = tree->params->num_param_iterations; //for now
    tree->initProgress(estimatedIterations, "Optimizing Model Parameters", "finished", "iteration", true );
    
    stopStoringTransMatrix();
    // modified by Thomas Wong on Sept 11, 15
    // no optimization of branch length in the first round
    double optimizeStartTime = getRealTime();
    double cur_lh = tree->computeLikelihood();
    tree->setCurScore(cur_lh);
    tree->trackProgress(1.0);
    
    if (verbose_mode >= VB_MED || write_info) {
        tree->hideProgress();
        int p = cout.precision(); //We'll restore it later
        if (verbose_mode >= VB_DEBUG) {
            cout.precision(17);
        }
        if (verbose_mode >= VB_MED) {
            cout << "1. Initial log-likelihood: " << cur_lh << " (took " <<
            (getRealTime() - optimizeStartTime) << " wall-clock sec)" << endl;
        } else {
            cout << "1. Initial log-likelihood: " << cur_lh << endl;
        }
        cout.precision(p);
        if (verbose_mode >= VB_MAX) {
            tree->printTree(cout);
            cout << endl;
        }
        tree->showProgress();
    }

    // For UpperBounds -----------
    //cout<<"MLCheck = "<<tree->mlCheck <<endl;
    if(tree->mlCheck == 0){
        tree->mlInitial = cur_lh;
    }
    // ---------------------------

    int i;
    for (i = 2; i < tree->params->num_param_iterations; i++) {
        double new_lh;

        // changed to optimise edge length first, and then Q,W,R inside the loop by Thomas on Sept 11, 15
        if (fixed_len == BRLEN_OPTIMIZE) {
            TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch lengths");
            new_lh = tree->optimizeAllBranches(min(i, 3), logl_epsilon);  // loop only 3 times in total (previously in v0.9.6 5 times)
        }
        else if (fixed_len == BRLEN_SCALE) {
            double scaling = 1.0;
            TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch scaling");
            new_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
        } else {
            new_lh = cur_lh;
        }
        TREE_LOG_LINE(*tree, VB_MAX, "Optimizing parameters");
        new_lh = optimizeParametersOnly(i, gradient_epsilon, new_lh);
        if (new_lh == 0.0) {
            if (fixed_len == BRLEN_OPTIMIZE) {
                TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch lengths (2nd time)");
                cur_lh = tree->optimizeAllBranches(tree->params->num_param_iterations, logl_epsilon);
            } else if (fixed_len == BRLEN_SCALE) {
                double scaling = 1.0;
                TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch scaling (2nd time)");
                cur_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
            }
            break;
        }
        if (verbose_mode >= VB_MED) {
            tree->hideProgress();
            model->writeInfo(cout);
            site_rate->writeInfo(cout);
            if (fixed_len == BRLEN_SCALE) {
                cout << "Scaled tree length: " << tree->treeLength() << endl;
            }
            tree->showProgress();
        }
        if (new_lh > cur_lh + logl_epsilon) {
            cur_lh = new_lh;
            if (write_info) {
                tree->hideProgress();
                if (verbose_mode >= VB_MED) {
                    cout << i << ". Current log-likelihood: " << cur_lh
                    << " (after " << (getRealTime() - optimizeStartTime) << " wall-clock sec)"
                    << endl;
                } else {
                    cout << i << ". Current log-likelihood: " << cur_lh << endl;
                }
                tree->showProgress();
            }
        } else {
            site_rate->classifyRates(new_lh);
            if (fixed_len == BRLEN_OPTIMIZE) {
                TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch lengths (3rd time)");
                cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
            } else if (fixed_len == BRLEN_SCALE) {
                double scaling = 1.0;
                TREE_LOG_LINE(*tree, VB_MAX, "Optimizing branch scaling (3rd time)");
                cur_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
            }
            break;
        }
        tree->trackProgress(1.0);
    }
    tree->doneProgress();
    
    //normalize rates s.t. branch lengths are #subst per site
    //if (Params::getInstance().optimize_alg_gammai != "EM")
    {
        double mean_rate = site_rate->rescaleRates();
        if (fabs(mean_rate-1.0) > 1e-6) {
            if (fixed_len == BRLEN_FIX)
            {
                outError("Fixing branch lengths not supported under specified site rate model");
            }
            tree->scaleLength(mean_rate);
            tree->clearAllPartialLH();
        }
    }
    if (Params::getInstance().root_find && tree->rooted && Params::getInstance().root_move_dist > 0) {
        cur_lh = tree->optimizeRootPosition(Params::getInstance().root_move_dist, write_info, logl_epsilon);
        if (verbose_mode >= VB_MED || write_info) {
            tree->hideProgress();
            cout << "Rooting log-likelihood: " << cur_lh << endl;
            tree->showProgress();
        }
    }
    if (verbose_mode >= VB_MED || write_info) {
        tree->hideProgress();
        cout << "Optimal log-likelihood: " << cur_lh << endl;
        tree->showProgress();
    }
    // For UpperBounds -----------
    if(tree->mlCheck == 0) {
        tree->mlFirstOpt = cur_lh;
    }
    // ---------------------------

    if (verbose_mode <= VB_MIN && write_info) {
        model->writeInfo(cout);
        site_rate->writeInfo(cout);
        if (fixed_len == BRLEN_SCALE) {
            tree->hideProgress();
            cout << "Scaled tree length: " << tree->treeLength() << endl;
            tree->showProgress();
        }
    }
    double elapsed_secs = getRealTime() - begin_time;
    if (write_info) {
        tree->hideProgress();
        cout << "Parameters optimization took " << i-1 << " rounds (" << elapsed_secs << " sec)" << endl;
        tree->showProgress();
    }
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

void ModelFactory::computeTransMatrix(double time, double *trans_matrix, int mixture) {
    if (!store_trans_matrix || !is_storing || model->isSiteSpecificModel()) {
        model->computeTransMatrix(time, trans_matrix, mixture);
        return;
    }
    int mat_size = model->num_states * model->num_states;
    iterator ass_it = find(round(time * 1e6));
    if (ass_it == end()) {
        // allocate memory for 3 matricies
        double *trans_entry = new double[mat_size * 3];
        trans_entry[mat_size] = trans_entry[mat_size+1] = 0.0;
        model->computeTransMatrix(time, trans_entry, mixture);
        ass_it = insert(value_type(round(time * 1e6), trans_entry)).first;
    } else {
        //if (verbose_mode >= VB_MAX)
            //cout << "ModelFactory bingo" << endl;
    }

    memcpy(trans_matrix, ass_it->second, mat_size * sizeof(double));
}

void ModelFactory::computeTransDerv(double time, double *trans_matrix,
    double *trans_derv1, double *trans_derv2, int mixture) {
    if (!store_trans_matrix || !is_storing || model->isSiteSpecificModel()) {
        model->computeTransDerv(time, trans_matrix, trans_derv1, trans_derv2, mixture);
        return;
    }
    int mat_size = model->num_states * model->num_states;
    iterator ass_it = find(round(time * 1e6));
    if (ass_it == end()) {
        // allocate memory for 3 matricies
        double *trans_entry = new double[mat_size * 3];
        trans_entry[mat_size] = trans_entry[mat_size+1] = 0.0;
        model->computeTransDerv(time, trans_entry, trans_entry+mat_size, trans_entry+(mat_size*2), mixture);
        ass_it = insert(value_type(round(time * 1e6), trans_entry)).first;
    } else if (ass_it->second[mat_size] == 0.0 && ass_it->second[mat_size+1] == 0.0) {
        double *trans_entry = ass_it->second;
        model->computeTransDerv(time, trans_entry, trans_entry+mat_size, trans_entry+(mat_size*2), mixture);
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
