/*
 * modelmixture.h
 *
 *  Created on: Nov 29, 2014
 *      Author: minh
 */

#ifndef MODELMIXTURE_H_
#define MODELMIXTURE_H_

#include "tree/phylotree.h"
#include "modelsubst.h"
#include "modelmarkov.h"
#include "nclextra/modelsblock.h"


extern const char* builtin_mixmodels_definition;

/**
 * create a substitution model
 * @param model_str model nme
 * @param freq_type state frequency type
 * @param freq_params frequency parameters
 * @param seqerr sequencing error model
 * @param tree associated phylo tree
 * @param count_rates TRUE to assign rates counted from alignment, FALSE to not initialize rates
 * @return substitution model created
 */
ModelSubst *createModel(string model_str, ModelsBlock *models_block,
                        StateFreqType freq_type, string freq_params,
                        PhyloTree *tree);


/**
 * mixture model
 */
class ModelMixture: virtual public ModelMarkov, public vector<ModelMarkov*> {
public:
    
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights);

    void initMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
    		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights);

    void initMem();

    /**
		constructor
		@param tree associated tree for the model
	*/
    ModelMixture(PhyloTree *tree);


    virtual ~ModelMixture();

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

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
	virtual void setMixtureWeight(int cat, double weight) { prop[cat] = weight; }

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual void setFixMixtureWeight(bool fix_weight) { this->fix_prop = fix_weight; }

	/**
	 * @param cat mixture class ID
	 * @return corresponding mixture model component
	 */
    virtual ModelSubst* getMixtureClass(int cat) { return at(cat); }

	/**
	 * @param cat mixture class ID
	 * @param m mixture model class to set
	 */
    virtual void setMixtureClass(int cat, ModelSubst* m) { at(cat) = (ModelMarkov*)m; }

	/**
		compute the state frequency vector
        @param mixture (optional) class for mixture model. 
            -1 to get weighted sum of class state frequency vector
		@param state_freq (OUT) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq, int mixture = 0);

	/**
		compute the transition probability matrix. One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
        @param mixture (optional) class for mixture model
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, int mixture = 0);


	/**
		compute the transition probability matrix.and the derivative 1 and 2
		@param time time between two events
        @param mixture (optional) class for mixture model
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	virtual void computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2, int mixture = 0);

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
        set number of optimization steps
        @param opt_steps number of optimization steps
    */
    virtual void setOptimizeSteps(int steps) { this->optimize_steps = steps; }

    /** @return true if model is fused with site_rate */
    bool isFused();

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
    	uint64_t mem = ModelMarkov::getMemoryRequired();
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

    /** number of optimization steps, default: ncategory*2 */
    int optimize_steps;    

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
