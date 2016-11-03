/*
 * modelunrest.h
 *
 *  Created on: 24/05/2016
 *      Author: Michael Woodhams
 */

#ifndef MODELUNREST_H_
#define MODELUNREST_H_

#include "modelnonrev.h"

class ModelUnrest: public ModelNonRev {
public:
	ModelUnrest(PhyloTree *tree, string model_params, bool count_rates);
	static bool validModelName(string model_name);
	void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
protected:
	virtual void setRates();
};

#endif /* MODELUNREST_H_ */
