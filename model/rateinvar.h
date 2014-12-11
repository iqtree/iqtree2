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
#ifndef RATEINVAR_H
#define RATEINVAR_H

#include "phylotree.h"
#include "rateheterogeneity.h"

const double MIN_PINVAR = 1e-6;
const double TOL_PINVAR = 1e-6;

/**
class for rate heterogeneity with a fraction of invariable sites

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateInvar : virtual public RateHeterogeneity
{
	friend class RateGammaInvar;

public:
	/**
		constructor
		@param p_invar_sites proportion of invariable sites
		@param tree associated phylogenetic tree
	*/
	RateInvar(double p_invar_sites, PhyloTree *tree);

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams();

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) { return 1.0 - p_invar; }

	/**
		get the proportion of invariable sites
		@return the proportion of invariable sites
	*/
	virtual double getPInvar() { return p_invar; }

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
		override function from Optimization class, used by the minimizeOneDimen() to optimize
		p_invar parameter
	*/
	virtual double computeFunction(double p_invar_value);

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);


	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return !fix_p_invar; }
	

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
		proportion of invariable sites
	*/
	double p_invar;
	
	/**
		TRUE to fix the proportion of invariable sites
	*/
	bool fix_p_invar;

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

};

#endif
