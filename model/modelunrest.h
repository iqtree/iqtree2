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
    typedef ModelMarkov super;
    /** constructor */
	ModelUnrest(PhyloTree *tree, const string& model_params,
                PhyloTree* report_to_tree);

    /**
     * true if model_name is the name of some known non-reversible model
     */
	static bool validModelName(const string& model_name);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double* lower_bound, double* upper_bound, 
                           bool*   bound_check) override;
    
    /**
     set the state frequency vector.
     @param state_freq state frequency vector. Assume state_freq has size of num_states
     */
    virtual void setStateFrequency(double *state_freq) override;

    /**
     start structure for checkpointing
     */
    virtual void startCheckpoint() override;
    
    /**
     save object into the checkpoint
     */
    virtual void saveCheckpoint() override;
    
    /**
     restore object from the checkpoint
     */
    virtual void restoreCheckpoint() override;

protected:
};

#endif /* MODELUNREST_H_ */
