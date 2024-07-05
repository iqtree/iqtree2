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
#include "modelbin.h"

ModelBIN::ModelBIN(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree)
: ModelMarkov(tree)
{
	init(model_name, model_params, freq, freq_params);
}

void ModelBIN::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
	ASSERT(num_states == 2); // make sure that you create model for Binary data
	StateFreqType def_freq = FREQ_UNKNOWN;
	name = model_name;
	full_name = model_name;
	if (name == "JC2") {
		def_freq = FREQ_EQUAL;
	} else if (name == "GTR2") {
		def_freq = FREQ_ESTIMATE;
	} else {
		readParameters(model_name);
	}
    if (freq_params != "") {
        readStateFreq(freq_params);
    }
    if (model_params != "") {
      readRates(model_params);
    }
	if (freq == FREQ_UNKNOWN || def_freq == FREQ_EQUAL) freq = def_freq;
	ModelMarkov::init(freq);
}

void ModelBIN::startCheckpoint() {
    checkpoint->startStruct("ModelBIN");
}

string ModelBIN::getNameParams(bool show_fixed_params) {
    //if (num_params == 0) return name;
    ostringstream retname;
    retname << name;
//    if (!fixed_parameters) {
//        retname << '{';
//        int nrates = getNumRateEntries();
//        for (int i = 0; i < nrates; i++) {
//            if (i>0) retname << ',';
//            retname << rates[i];
//        }
//        retname << '}';
//    }
//    getNameParamsFreq(retname);
    retname << freqTypeString(freq_type, phylo_tree->aln->seq_type, true);
    if (freq_type == FREQ_EMPIRICAL || freq_type == FREQ_ESTIMATE ||
        (freq_type == FREQ_USER_DEFINED)) {
        retname << "{" << state_freq[0];
        for (int i = 1; i < num_states; i++)
            retname << "," << state_freq[i];
        retname << "}";
    }
    return retname.str();
}

void ModelBIN::printMrBayesModelText(RateHeterogeneity* rate, ofstream& out, string partition, string charset, bool isSuperTree, bool inclParams) {
    // MrBayes does not support Invariable Modifier for Binary data
    if (rate->isFreeRate() || rate->getPInvar() > 0.0) {
        warnLogStream("MrBayes does not support Invariable Sites with Binary Data! +I has been ignored!", out);
    }

    // Lset Parameters
    out << "  lset applyto=(" << partition << ") rates=";

    // Free Rate should be substituted by +G (+I not supported)
    bool hasGamma = rate->getGammaShape() != 0.0 || rate->isFreeRate();
    if (hasGamma) {
        // Rate Categories + Gamma
        out << "gamma ngammacat=" << rate->getNRate();
    } else
        out << "equal";

    out << ";" << endl;

    if (!inclParams) {
        if (freq_type == FREQ_EQUAL) out << "  prset applyto=(" << partition << ") statefreqpr=fixed(equal);" << endl;
        return;
    }

    // Prset Parameters
    out << "  prset applyto=(" << partition << ")";

    // Freerate (+R)
    // Get replacement Gamma Shape
    if (rate->isFreeRate()) {
        printMrBayesFreeRateReplacement(rate, charset, out, false);
    }

    // Gamma Distribution (+G/+R)
    // Dirichlet is not available here, use fixed
    if (rate->getGammaShape() > 0.0)
        out << " shapepr=fixed(" << minValueCheckMrBayes(rate->getGammaShape()) << ")";

    // State Frequencies
    if (freq_type == FREQ_EQUAL)
        out << " statefreqpr=fixed(equal)";
    else {
        out << " statefreqpr=dirichlet(";
        for (int i = 0; i < num_states; ++i) {
            if (i != 0) out << ", ";
            out << minValueCheckMrBayes(state_freq[i]);
        }
        out << ")";
    }

    out << ";" << endl;
}
