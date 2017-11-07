/*
 * modelunrest.h
 *
 *  Created on: 24/05/2016
 *      Author: Michael Woodhams
 */

#ifndef MODELUNREST_H_
#define MODELUNREST_H_

#include "modelmarkov.h"

class ModelUnrest: public ModelMarkov {
public:

    /** constructor */
	ModelUnrest(PhyloTree *tree, string model_params);

    /**
     * true if model_name is the name of some known non-reversible model
     */
	static bool validModelName(string model_name);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
    
protected:

	/**
	    Model parameters - cached so we know when they change, and thus when
	    recalculations are needed.

	 */
	double *model_parameters;

	/**
	 * Called from getVariables to update the rate matrix for the new
	 * model parameters.
	 */
	virtual void setRates();
};

#endif /* MODELUNREST_H_ */
