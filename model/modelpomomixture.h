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
class ModelPoMoMixture : public ModelPoMo {

public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelPoMoMixture(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params,
                     string pomo_rate_str);

    ~ModelPoMoMixture();

protected:

    /** associated mixture model */
    ModelMixture *pomo_mixture;

};

#endif /* modelpomomixture_h */
