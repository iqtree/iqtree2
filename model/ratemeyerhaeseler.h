/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef RATEMEYERHAESELER_H
#define RATEMEYERHAESELER_H

#include "rateheterogeneity.h"
#include "utils/tools.h"
#include "tree/iqtree.h"


/**
Implementation for site-specific rates of Meyer & von Haeseler (2003)
Inherited from Optimization and the double vector for storing site-specific rates

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateMeyerHaeseler : public RateHeterogeneity, public DoubleVector
{
public:
	/**
		constructor
	*/
    RateMeyerHaeseler(char *file_name, PhyloTree *tree, bool rate_type);

    RateMeyerHaeseler();

	/**
		destructor
	*/
    ~RateMeyerHaeseler();

	void readRateFile(char *rate_file);

	/**
		@return true 
	*/
	virtual bool isSiteSpecificRate() { return true; }

	/**
		get the number of rate categories. 
		@return the number of rate categories
	*/
	//virtual int getNRate() { return size(); }


	/**
		return the number of dimensions
	*/
	virtual int getNDim();

	/**
		get the rate of a specified category
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	//virtual double getRate(int category);

	/**
		get the rate of a specified site-pattern. Default returns 1.0 since it is homogeneous model
		@param ptn pattern ID 
		@return the rate of the specified site-pattern
	*/
	virtual double getPtnRate(int ptn);

	/**
		Compute site-specific rates. Override this for Gamma model
		@param pattern_rates (OUT) pattern rates. Resizing if necesary
        @return total number of categories
	*/
	virtual int computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat);

	void getRates(DoubleVector &rates);


	void setRates(DoubleVector &rates);

	void initializeRates();

	/**
		optimize parameters, the rates in this case
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double epsilon);

	/**
		optimize rate of site
		@param pattern target pattern
		@return the optimized rate value, also update the corresponding element of the vector
	*/
	double optimizeRate(int pattern);

	/**
		optimize rates of all site-patterns
	*/
	virtual void optimizeRates();


	/**
		This function is inherited from Optimization class for optimizting site rates 
		@param value x-value of the function
		@return f(value) of function f you want to minimize
	*/
	virtual double computeFunction(double value);

	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		@param value x-value of the function
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
		@return f(value) of function f you want to minimize
	*/
	virtual void computeFuncDerv(double value, double &df, double &ddf);


	void runIterativeProc(Params &params, IQTree &tree);

	/**
		distance matrix inferred from the path lengths of the tree (not from the sequences)
	*/
	double *dist_mat;


protected:

	char *rate_file;

	/**
		current pattern under optimization. Note that this is not thread-safe
	*/
	int optimizing_pattern;

	/**
		FALSE to use MH Model, FALSE for using tree-likelihood
	*/
	bool rate_mh;

	double cur_scale;

	PhyloTree *ptn_tree;

	void prepareRateML(IntVector &ptn_id);
	void completeRateML();
};

#endif
