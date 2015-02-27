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
#ifndef RATEGAMMAINVAR_H
#define RATEGAMMAINVAR_H

#include "rateinvar.h"
#include "rategamma.h"

/**
class for I+G rate heterogeneity

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateGammaInvar : public RateInvar, public RateGamma
{
public:
 	/**
		constructor
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
	*/
    RateGammaInvar(int ncat, double shape, bool median, double p_invar_sites, bool simultaneous, PhyloTree *tree);

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) { return (1.0-p_invar)/ncategory; }

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams();

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
		return the number of dimensions
	*/
	virtual int getNDim() { return RateInvar::getNDim() + RateGamma::getNDim(); }

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

	/** TRUE to jointly optimize gamma shape and p_invar using BFGS, default: FALSE */
	bool joint_optimize;

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

#endif
