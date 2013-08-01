/*
 * modelcodonsemiempirical.h
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#ifndef MODELCODONSEMIEMPIRICAL_H_
#define MODELCODONSEMIEMPIRICAL_H_

#include "modelcodonempirical.h"
#include "modelcodonparametric.h"

class ModelCodonSemiEmpirical: public ModelCodonEmpirical, public ModelCodonParametric {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodonSemiEmpirical(const char *model_name, string model_params, StateFreqType freq, string freq_params,
    		PhyloTree *tree, bool count_rates = true);


	/**
	 * destructor
	 */
	virtual ~ModelCodonSemiEmpirical();

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

};

#endif /* MODELCODONSEMIEMPIRICAL_H_ */
