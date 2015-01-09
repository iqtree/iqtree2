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
class ModelMixture: public ModelGTR, vector<ModelGTR*> {
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
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);


	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out);

	/**
		rates of mixture components
	*/
//	double *mix_rates;

	/**
	 * weight of each sub-model (must sum to 1)
	 */
	double *prop;

	/**
	 * TRUE to fix model weights
	 */
	bool fix_prop;

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables);

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
	*/
	virtual void getVariables(double *variables);

};

#endif /* MODELMIXTURE_H_ */
