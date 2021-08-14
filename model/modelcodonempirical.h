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
class ModelCodonEmpirical: public ModelCodon {
public:
	typedef ModelCodon super;
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodonEmpirical(const char *model_name, const string& model_params, 
	                    StateFreqType freq, const string& freq_params,
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
	virtual void init(const char *  model_name, const string& model_params, 
	                  StateFreqType freq,       const string& freq_params,
					  PhyloTree* report_to_tree) override;

	/**
	 * read codon model from a stream, modying rates and state_freq accordingly
	 * @param in input stream containing lower triangular matrix of rates, frequencies and list of codons
	 */
	void readCodonModel(istream &in);

		//Supporting functions
		double** readRateMatrix        (istream& in) const;
		void     forgetRateMatrix      (double** q)  const;

		void     readCodonsAndStateMap (istream &in, 
                                        StrVector& codons,
									    IntVector& state_map);
		void     calculateRatesAndFrequencies(const StrVector& codons, 
		                                      const IntVector& state_map,
											  double** q,        
											  double* f);

		double*  readFrequencyVector   (istream& in) const;
		void     forgetFrequencyVector (double*  f)  const;
};

#endif /* MODELCODONEMPIRICAL_H_ */
