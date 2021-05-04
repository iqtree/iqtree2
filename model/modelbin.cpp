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

ModelBIN::ModelBIN(PhyloTree *tree, PhyloTree* report_to_tree) 
    : ModelMarkov(tree, report_to_tree) {
}

ModelBIN::ModelBIN(const char *model_name, string model_params,
                   StateFreqType freq, string freq_params,
                   PhyloTree *tree, PhyloTree* report_to_tree)
: ModelMarkov(tree, report_to_tree)
{
	init(model_name, model_params, freq, freq_params, report_to_tree);
}

void ModelBIN::init(const char *model_name, string model_params,
                    StateFreqType freq, string freq_params,
                    PhyloTree* report_to_tree)
{
	ASSERT(num_states == 2); // make sure that you create model for Binary data
	StateFreqType def_freq = StateFreqType::FREQ_UNKNOWN;
	name = model_name;
	full_name = model_name;
	if (name == "JC2") {
		def_freq = StateFreqType::FREQ_EQUAL;
	} else if (name == "GTR2") {
		def_freq = StateFreqType::FREQ_ESTIMATE;
	} else {
		readParameters(model_name, true, report_to_tree);
	}
    if (freq_params != "") {
        readStateFreq(freq_params, report_to_tree);
    }
    if (model_params != "") {
        readRates(model_params);
    }
    if (freq == StateFreqType::FREQ_UNKNOWN || 
        def_freq == StateFreqType::FREQ_EQUAL) {
        freq = def_freq;
    }
    ModelMarkov::init(freq, report_to_tree);
}

void ModelBIN::startCheckpoint() {
    checkpoint->startStruct("ModelBIN");
}

string ModelBIN::getNameParams() {
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
    if (freq_type == StateFreqType::FREQ_EMPIRICAL || 
        freq_type == StateFreqType::FREQ_ESTIMATE ||
        freq_type == StateFreqType::FREQ_USER_DEFINED) {
        retname << "{" << state_freq[0];
        for (int i = 1; i < num_states; i++) {
            retname << "," << state_freq[i];
        }
        retname << "}";
    }
    return retname.str();
}
