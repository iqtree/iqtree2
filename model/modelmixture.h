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
#include "modelsblock.h"


const char OPEN_BRACKET = '{';
const char CLOSE_BRACKET = '}';

extern const string builtin_mixmodels_definition;

/**
 * create a substitution model
 * @param model_str model nme
 * @param freq_type state frequency type
 * @param freq_params frequency parameters
 * @param tree associated phylo tree
 * @param count_rates TRUE to assign rates counted from alignment, FALSE to not initialize rates
 * @return substitution model created
 */
ModelSubst *createModel(string model_str, ModelsBlock *models_block, StateFreqType freq_type, string freq_params,
		PhyloTree *tree, bool count_rates = true, string pomo_rate_str = "");


/**
 * mixture model
 */
class ModelMixture: virtual public ModelGTR, public vector<ModelGTR*> {
public:
    
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool count_rates = true);

    void initMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool count_rates = true);

    void initMem();

    /**
		constructor
		@param tree associated tree for the model
	*/
    ModelMixture(PhyloTree *tree, bool count_rates = true);


    virtual ~ModelMixture();

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();


	/**
	 * @return TRUE if this is a mixture model, FALSE otherwise
	 */
	virtual bool isMixture() { return true; }


	/**
	 * @return the number of mixture model components
	 */
	virtual int getNMixtures() {return size(); }

 	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual double getMixtureWeight(int cat) { return prop[cat]; }

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual double getMixtureWeight(int cat) { return prop[cat]; }

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual void setMixtureWeight(int cat, double weight) { prop[cat] = weight; }

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual void setFixMixtureWeight(bool fix_prop) { this->fix_prop = fix_prop; }

	/**
	 * @param cat mixture class ID
	 * @return corresponding mixture model component
	 */
    virtual ModelSubst* getMixtureClass(int cat) { return at(cat); }

	/**
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
		@return the number of dimensions corresponding to state frequencies
	*/
	virtual int getNDimFreq();
	
	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

    /** 
        optimize mixture weights using EM algorithm 
        @return log-likelihood of optimized weights
    */
    double optimizeWeights();

    /** 
        optimize rate parameters using EM algorithm
        @param gradient_epsilon
        @return log-likelihood of optimized parameters
    */
    double optimizeWithEM(double gradient_epsilon);


	/**
		optimize model parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon);

	/**
	 * @return TRUE if parameters are at the boundary that may cause numerical unstability
	 */
	virtual bool isUnstableParameters();

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
	 * @return model name
	 */
	virtual string getName();

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams();

    /**
     * compute the memory size for the model, can be large for site-specific models
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired() {
    	uint64_t mem = ModelGTR::getMemoryRequired();
    	for (iterator it = begin(); it != end(); it++)
    		mem += (*it)->getMemoryRequired();
    	return mem;
    }

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

	bool optimizing_submodels;

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
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(double *variables);

};

#endif /* MODELMIXTURE_H_ */
