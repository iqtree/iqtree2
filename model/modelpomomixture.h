//
//  modelpomomixture.h
//  iqtree
//  Mixture PoMo models to include e.g. Gamma-rate heterogeneity
//
//  Created by Minh Bui on 7/22/16.
//
//

#ifndef modelpomomixture_h
#define modelpomomixture_h

#include <stdio.h>
#include "modelpomo.h"
#include "modelmixture.h"

/**
    Mixture PoMo models
*/
class ModelPoMoMixture : public ModelMixture {

public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelPoMoMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool is_reversible,
                     string pomo_params, bool count_rates = true);

    void initMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool count_rates = true);

};

#endif /* modelpomomixture_h */
