/*
 * modelmixture.h
 *
 *  Created on: Nov 29, 2014
 *      Author: minh
 */

#ifndef MODELMIXTURE_H_
#define MODELMIXTURE_H_

#include "phylotree.h"
#include "modelsubst.h"
#include "modelgtr.h"


const char OPEN_BRACKET = '{';
const char CLOSE_BRACKET = '}';

/**
 * create a substitution model
 * @param model_str model nme
 * @param freq_type state frequency type
 * @param freq_params frequency parameters
 * @param tree associated phylo tree
 * @param count_rates TRUE to assign rates counted from alignment, FALSE to not initialize rates
 * @return substitution model created
 */
ModelSubst *createModel(string model_str, StateFreqType freq_type, string freq_params,
		PhyloTree *tree, bool count_rates = true);


/**
 * mixture model
 */
class ModelMixture: public ModelGTR, vector<ModelSubst*> {
public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelMixture(string model_name, string model_list, StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates = true);

    virtual ~ModelMixture();


	/**
	 * @return TRUE if this is a mixture model, FALSE otherwise
	 */
	virtual bool isMixture() { return true; }


	/**
	 * @return the number of mixture model components
	 */
	virtual int getNMixtures() {return size(); }

	/**
	 * proportion of sites for each sub-models
	 */
	double *prop;

};

#endif /* MODELMIXTURE_H_ */
