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
#include "tools.h"

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
    RateMeyerHaeseler();

	/**
		destructor
	*/
    ~RateMeyerHaeseler();


	/**
		get the number of rate categories. 
		@return the number of rate categories
	*/
	virtual int getNRate() { return size(); }

	/**
		get the rate of a specified category
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) { assert(category < size()); return at(category); }


	/**
		optimize parameters, the rates in this case
		@return the best likelihood 
	*/
	virtual double optimizeParameters();

	/**
		optimize rate of site
		@param site target site
		@return the optimized rate value, also update the corresponding element of the vector
	*/
	double optimizeSiteRate(int site);

	/**
		This function is inherited from Optimization class for optimizting site rates 
		@param value x-value of the function
		@return f(value) of function f you want to minimize
	*/
	virtual double computeFunction(double value);


	/**
		write site-rates to a file in the following format:
		1  rate_1
		2  rate_2
		....
	*/
	void writeSiteRates(const char *file_name);

	/**
		distance matrix inferred from the path lengths of the tree (not from the sequences)
	*/
	double *dist_mat;

	/**
		current site under optimization. Note that this is not thread-safe
	*/
	int optimizing_site;
};

#endif
