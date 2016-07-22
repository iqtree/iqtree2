//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"

ModelPoMoMixture::ModelPoMoMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool is_reversible,
                     string pomo_params, bool count_rates)
	: ModelMixture(orig_model_name, model_name, model_list, models_block,
		freq, freq_params, tree, optimize_weights, count_rates)
{
}
