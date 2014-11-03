/*
 * modelcodonempirical.h
 *
 *  Created on: May 29, 2013
 *      Author: minh
 */

#ifndef MODELCODONEMPIRICAL_H_
#define MODELCODONEMPIRICAL_H_

#include "modelcodon.h"

/**
 * empirical codon model (e.g., Kosiol et al. 2007)
 */
class ModelCodonEmpirical: virtual public ModelCodon {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodonEmpirical(const char *model_name, string model_params, StateFreqType freq, string freq_params,
    		PhyloTree *tree, bool count_rates = true);

	/**
	 * destructor
	 */
	virtual ~ModelCodonEmpirical();

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);


	/**
	 * read codon model from a stream, modying rates and state_freq accordingly
	 * @param in input stream containing lower triangular matrix of rates, frequencies and list of codons
	 */
	void readCodonModel(istream &in);

};

#endif /* MODELCODONEMPIRICAL_H_ */
