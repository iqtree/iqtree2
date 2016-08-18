//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"

ModelPoMoMixture::ModelPoMoMixture(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params, string pomo_rate_str)
	: ModelPoMo(tree)
{

    init(model_name, model_params, freq_type, freq_params, is_reversible, pomo_params);

    // initialize pomo_mixture
    pomo_mixture = new ModelMixture(tree, false);
}


ModelPoMoMixture::~ModelPoMoMixture() {
    delete pomo_mixture;
}

