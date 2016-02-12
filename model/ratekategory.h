/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  BUI Quang Minh <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef RATEKATEGORY_H
#define RATEKATEGORY_H

#include "rateheterogeneity.h"

class PhyloTree;

/**
among-site-rate model that sites are categorized into K categories of equal proportion
where the K rates are optimized by ML instead of the Gamma distribution

*/
class RateKategory : virtual public RateHeterogeneity
{
public:
	/**
		constructor
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
	*/
    RateKategory(int ncat, PhyloTree *tree);

	/**
		destructor
	*/
    virtual ~RateKategory();

	
	/**
		@return the number of rate categories
	*/
	virtual int getNRate() { return ncategory; }

	/**
		get the number of rate categories for site-specific category model
		@return the number of rate categories
	*/
	virtual int getNDiscreteRate() { return ncategory; }

	/**
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) { return rates[category]; }

	/**
		Compute site-specific rates. Override this for Gamma model
		@param pattern_rates (OUT) pattern rates. Resizing if necesary
        @return total number of categories
	*/
	virtual int computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat);

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);


	/**
		optimize model parameters
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double epsilon);

	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return (ncategory-1); }

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
		number of rate categories
	*/
	int ncategory;

	/**
		rates, containing ncategory elements
	*/
	double *rates;

	
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

#endif // RATEKATEGORY_H
