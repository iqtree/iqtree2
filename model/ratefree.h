/*
 * ratefree.h
 *
 *  Created on: Nov 3, 2014
 *      Author: minh
 */

#ifndef RATEFREE_H_
#define RATEFREE_H_

#include "rategamma.h"

class ModelFactory;

class RateFree: public RateGamma {
public:
	typedef RateGamma super;
	/**
		constructor (used by, for example, YAMLRateModelWrapper)
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
		@param report_to_tree send any log messages to this tree.
	*/
	RateFree(int ncat, PhyloTree *tree, PhyloTree* report_to_tree);

	/**
		constructor
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
        @param opt_alg optimization algorithm (1-BFGS, 2-BFGS, EM)
	*/
    RateFree(int ncat, double start_alpha, string params, bool sorted_rates, string opt_alg, PhyloTree *tree);

	virtual ~RateFree();

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

	/**
		@return true if this is a Gamma model (default: false)
	*/	
    virtual int isGammaRate() const override { return 0; }

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams() const override;

	virtual const string& getOptimizationAlgorithm() const;

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) const override { return prop[category]; }

	virtual void   setFixProportions(bool fixed);

	virtual void   setFixRates(bool fixed);

	virtual void   setGammaShape(double shape) override;

	virtual void   setOptimizationAlgorithm(const std::string& algorithm);

	virtual bool   isOptimizingProportions() const;

	virtual bool   isOptimizingRates() const;

	virtual bool   isOptimizingShapes() const override;

	virtual bool   areProportionsFixed() const;

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) override;

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double* lower_bound, double* upper_bound, 
	                       bool*   bound_check) override;

	/**
		optimize parameters. Default is to optimize gamma shape
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon,
                                      PhyloTree* report_to_tree) override;

    /** 
        optimize rate parameters using EM algorithm 
        @return log-likelihood of optimized parameters
    */
    double optimizeWithEM(PhyloTree* report_to_tree);

	void   doEStep(intptr_t nptn, double* new_prop, size_t nmix);

	int    doMStep(double* new_prop, size_t nmix);

	bool regularizeProportions(double* new_prop, size_t nmix, 
                               size_t  maxpropid);

	void optimizeRatesOneByOne(PhyloTree* tree, intptr_t nptn,
	                           bool& converged, ModelFactory* model_fac);

	/**
		return the number of dimensions
	*/
	virtual int getNDim() const override;

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) override;

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out) override;

    /**
        set number of rate categories
        @param ncat #categories
    */
	virtual void setNCategory(int ncat) override;

    /**
        initialize from checkpoint rates and prop from rate model with #category-1
    */
    virtual void initFromCatMinusOne() override;

	/**
	 * used to normal branch lengths if mean rate is not equal to 1 (e.g. FreeRate model)
	 * @return mean rate, default = 1
	 */
	virtual double meanRates() const override;

	/**
	 * rescale rates s.t. mean rate is equal to 1, useful for FreeRate model
	 * @return rescaling factor
	 */
	virtual double rescaleRates() override;

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables) override;

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(const double *variables) override;

	/**
	    sort updated/re-normalized rates
	 */
	virtual void sortUpdatedRates();

	/**
	 * proportion of sites for each rate categories
	 */
	double *prop;

	/** 1 to fix weights, 2 to fix both weights and rates */
	int fix_params;
    
    /** true to sort rate in increasing order, false otherwise */
    bool sorted_rates;

    /** 0: no, 1: rates, 2: weights */
    int optimizing_params;

    string optimize_alg;

	double proportion_tolerance;
	double rate_tolerance;
	
	virtual void setProportionTolerance(double tol) override;
	virtual void setRateTolerance(double tol) override;
};

#endif /* RATEFREE_H_ */
