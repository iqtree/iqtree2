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
#ifndef RATEMEYERDISCRETE_H
#define RATEMEYERDISCRETE_H

#include "ratemeyerhaeseler.h"

/**
The discrete version of Meyer & von Haeseler rate class

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateMeyerDiscrete : public RateMeyerHaeseler
{
public:
 	/**
		constructor
		@param ncat number of rate categories
		@param cat_type category type, bitwise, incl. CAT_LOG, CAT_MEAN, CAT_PATTERN
		@param file_name rate file name, NULL if not inputed
		@param tree associated phylo tree
   */
   RateMeyerDiscrete(int ncat, int cat_type, char *file_name, PhyloTree *tree, bool rate_type);

   RateMeyerDiscrete();

	/**
		destructor
	*/
    virtual ~RateMeyerDiscrete();

	/**
		get the number of rate categories for site-specific category model
		@return the number of rate categories
	*/
	virtual int getNDiscreteRate();

	/**
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category);

	/**
		get the rate of a specified site-pattern. Default returns 1.0 since it is homogeneous model
		@param ptn pattern ID 
		@return the rate of the specified site-pattern
	*/
	virtual double getPtnRate(int ptn);

	/**
		get rate category of a specified site-pattern. 
		@param ptn pattern ID 
		@return the rate category of the specified site-pattern
	*/
	virtual int getPtnCat(int ptn);

	/**
		Compute site-specific rates. Override this for Gamma model
		@param pattern_rates (OUT) pattern rates. Resizing if necesary
        @return total number of categories        
	*/
	virtual int computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat);

	virtual bool isSiteSpecificRate();

	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return ncategory; }

	/**
		optimize rates of all site-patterns
		compute categorized rates from the "continuous" rate of the original Meyer & von Haeseler model.
		The current implementation uses the k-means algorithm with k-means++ package.
	*/
	virtual double optimizeParameters(double epsilon);

	/**
		classify rates into categories.
		@param tree_lh the current tree log-likelihood
	*/
	virtual double classifyRates(double tree_lh);

	/**
		classify rates into categories using k-means++ method.
		@return tree likelihood
	*/
	double classifyRatesKMeans();


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

	double optimizeCatRate(int cat);

	void normalizeRates();

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

protected:

	/**
		number of rate categories
	*/
	int ncategory;

	/**
		category index for every pattern
	*/
	int *ptn_cat;

	/**
		rates, containing ncategory elements
	*/
	double *rates;

	/**
		false at beginning, true after continuous rates were optimized
	*/
	bool is_categorized;

	int mcat_type;
	
	/**
		current category under optimization. Note that this is not thread-safe
	*/
	int optimizing_cat;

};

#endif
