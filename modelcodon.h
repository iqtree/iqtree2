/*
 * modelcodon.h
 *
 *  Created on: May 24, 2013
 *      Author: minh
 */

#ifndef MODELCODON_H_
#define MODELCODON_H_

#include "gtrmodel.h"

/**
 * Codon substitution models
 */
class ModelCodon: public GTRModel {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelCodon(const char *model_name, string model_params, StateFreqType freq, string freq_params,
    		PhyloTree *tree, bool count_rates = true);
	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

    /**
     * destructor
     */
	virtual ~ModelCodon();
};

#endif /* MODELCODON_H_ */
