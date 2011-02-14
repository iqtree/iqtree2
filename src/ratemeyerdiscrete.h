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
   */
   RateMeyerDiscrete(int ncat);

	/**
		destructor
	*/
    virtual ~RateMeyerDiscrete();

	/**
		@return the number of rate categories
	*/
	//virtual int getNRate() { return ncategory; }

	/**
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	//virtual double getRate(int category) { assert(category < ncategory); return rates[category]; }

	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return ncategory; }

	/**
		optimize rates of all site-patterns
		compute categorized rates from the "continuous" rate of the original Meyer & von Haeseler model.
		The current implementation uses the k-means algorithm with k-means++ package.
	*/
	virtual void optimizeRates();

	/**
		classify rates into categories.
		@param tree_lh the current tree log-likelihood
	*/
	virtual double classifyRates(double tree_lh);

	/**
		classify rates into categories using k-means++ method.
	*/
	void classifyRatesKMeans();

protected:

	/**
		number of rate categories
	*/
	int ncategory;

	/**
		rates, containing ncategory elements
	*/
	double *rates;

};

#endif
