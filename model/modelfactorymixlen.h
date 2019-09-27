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
	ModelFactoryMixlen(Params &params, string &model_name, PhyloTree *tree, ModelsBlock *models_block);

	/**
		optimize model parameters and tree branch lengths
		@param fixed_len TRUE to fix branch lengths, default is false
		@return the best likelihood 
	*/
	virtual double optimizeParameters(int fixed_len = BRLEN_OPTIMIZE, bool write_info = true,
                                      double logl_epsilon = 0.1, double gradient_epsilon = 0.0001);

    /**
        sort classes in ascending order of tree lengths
        @return tree string with sorted branch lengths
    */
    string sortClassesByTreeLength();


    /**
     * @param brlen_type either BRLEN_OPTIMIZE, BRLEN_FIX or BRLEN_SCALE
     * @return #parameters of the model + # branches
     */
    virtual int getNParameters(int brlen_type);

};
