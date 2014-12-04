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

class RateFreeInvar: public RateInvar, public RateFree {
public:

 	/**
		constructor
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
	*/
    RateFreeInvar(int ncat, double p_invar_sites, PhyloTree *tree);


	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return RateInvar::getNDim() + RateFree::getNDim(); }

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) { return (1.0-p_invar)*prop[category]; }

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams() {
		return RateInvar::getNameParams() + RateFree::getNameParams();
	}

	/**
		override function from Optimization class, used by the minimizeOneDimen() to optimize
		p_invar or gamma shape parameter.
		@param value value of p_invar (if cur_optimize == 1) or gamma shape (if cur_optimize == 0).
	*/
	virtual double computeFunction(double value);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		optimize parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double epsilon);


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

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

private:

	/**
		current parameter to optimize. 0 if gamma shape or 1 if p_invar.
	*/
	int cur_optimize;

};

#endif /* RATEFREEINVAR_H_ */
