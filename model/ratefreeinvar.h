/*
 * ratefreeinvar.h
 *
 *  Created on: Nov 7, 2014
 *      Author: minh
 */

#ifndef RATEFREEINVAR_H_
#define RATEFREEINVAR_H_

#include "rateinvar.h"
#include "ratefree.h"

class RateFreeInvar: public RateFree {
public:
	typedef RateFree super;

 	/**
		constructor (called from YAMLRateFreeInvar).
		@param ncat           number of rate categories
		@param tree           associated phylogenetic tree
		@param report_to_tree tree or what-have-you to log to
		@note  pretends that the proportion of invariant sites is 10%
		       (it gts overridden later).
	*/
	RateFreeInvar(int ncat, PhyloTree* tree, PhyloTree* report_to_tree);

 	/**
		constructor
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
	*/
    RateFreeInvar(int ncat, double start_alpha, const string& params, 
	              bool sorted_rates, double p_invar_sites, 
				  string opt_alg, PhyloTree* tree);

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
		return the number of dimensions
	*/
	virtual int getNDim() const override;

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) const override;

	/**
		get the rate of a specified category. Default returns 1.0 since it is homogeneous model
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) const override;

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams() const override;

	/**
		override function from Optimization class, 
		used by the minimizeOneDimen() to optimize
		p_invar or gamma shape parameter.
		@param value value of p_invar (if cur_optimize == 1)
		             or gamma shape (if cur_optimize == 0).
	*/
	virtual double computeFunction(double value) override;

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double* lower_bound, double* upper_bound, 
	                       bool*   bound_check) override;

	/**
		optimize parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon,
                                      PhyloTree* report_to_tree) override;


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) override;

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

	virtual void setNCategory(int ncat) override;

#ifdef _MSC_VER
	//MSVC generates warning messages about these member functions
	//being inherited "via dominance". Explictly declaring them
	//instead shuts those warnings up.
	virtual int    getNRate()         const override { return super::getNRate(); }
	virtual void   setRate(int category, double value) override { 
		super::setRate(category, value); 
	}
	virtual int    getNDiscreteRate() const override { return super::getNDiscreteRate(); }
	virtual double getGammaShape()    const override { return super::getGammaShape(); }
	virtual void   setGammaShape(double gs) override { super::setGammaShape(gs); }
	virtual bool   isFixGammaShape()  const override { return super::isFixGammaShape(); }
	virtual void   setFixGammaShape(bool fixGammaShape) override { 
		super::setFixGammaShape(fixGammaShape); 
	}
	virtual int    isGammaRate()      const override { return super::isGammaRate(); }
	virtual int    computePatternRates(double* pattern_lh_cat, 
	                                   DoubleVector& pattern_rates, 
									   IntVector& pattern_cat) override {
		return super::computePatternRates(pattern_lh_cat, pattern_rates, pattern_cat);
	}

	virtual double getPInvar()        const override { return invar.getPInvar(); }
	virtual void   setPInvar(double pInvar) override { invar.setPInvar(pInvar); }
	virtual bool   isFixPInvar()      const override { return invar.isFixPInvar(); }
	void    setFixPInvar(bool fixPInvar)    override { invar.setFixPInvar(fixPInvar); }


	virtual double meanRates()        const override { return super::meanRates(); }
	virtual double rescaleRates()           override { return super::rescaleRates(); }
	virtual void initFromCatMinusOne()      override { super::initFromCatMinusOne(); }
#endif

protected:
	RateInvar invar;

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

private:

	/**
		current parameter to optimize. 0 if gamma shape or 1 if p_invar.
	*/
	int cur_optimize;

};

#endif /* RATEFREEINVAR_H_ */
