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
	ModelUnrest(PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params);

    /**
     * true if model_name is the name of some known non-reversible model
     */
	static bool validModelName(string model_name);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
    
    /**
     set the state frequency vector.
     @param state_freq state frequency vector. Assume state_freq has size of num_states
     */
    virtual void setStateFrequency(double *state_freq);

    /**
     start structure for checkpointing
     */
    virtual void startCheckpoint();
    
    /**
     save object into the checkpoint
     */
    virtual void saveCheckpoint();
    
    /**
     restore object from the checkpoint
     */
    virtual void restoreCheckpoint();

protected:
    virtual void setRates();
};

#endif /* MODELUNREST_H_ */
