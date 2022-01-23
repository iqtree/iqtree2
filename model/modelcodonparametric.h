/*
 * modelcodonparametric.h
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#ifndef MODELCODONPARAMETRIC_H_
#define MODELCODONPARAMETRIC_H_

#include "modelcodon.h"

/**
 * parametric codon model (e.g., Goldman-Yang, Muse-Gaut)
 */
class ModelCodonParametric: public ModelCodon {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodonParametric(const char* model_name, const std::string& model_params, 
	                     StateFreqType freq,     const std::string& freq_params,
    		             PhyloTree* tree, bool count_rates = true);

    /**
     * destructor
     */
	virtual ~ModelCodonParametric();

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *  model_name, const string& model_params, 
	                  StateFreqType freq,       const string& freq_params,
					  PhyloTree* report_to_tree) override;

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) override;

protected:

	/** initialize Muse-Gaut 1994 model */
	void initMG94();

	/** initialize Goldman-Yang 1994 model (simplified version with 2 parameters omega and kappa */
	void initGY94();

	void setRateGroup(IntVector& upper_triangle_entry_to_rate);

	void setRateGroupConstraint(const char* constraint);

};

#endif /* MODELCODONPARAMETRIC_H_ */
