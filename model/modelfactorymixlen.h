/*
 * modelfactorymixlen.h
 *
 *  Created on: Sep 2, 2015
 *      Author: minh
 */


#include "modelfactory.h"

class ModelFactoryMixlen : public ModelFactory {

public:

	/**
		constructor
		create substitution model with possible rate heterogeneity. Create proper class objects
		for two variables: model and site_rate. It takes the following field of params into account:
			model_name, num_rate_cats, freq_type, store_trans_matrix
		@param params program parameters
		@param tree associated phylogenetic tree
	*/
	ModelFactoryMixlen(Params &params, PhyloTree *tree, ModelsBlock *models_block);

};
