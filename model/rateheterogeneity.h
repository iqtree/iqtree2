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
#ifndef RATEHETEROGENEITY_H
#define RATEHETEROGENEITY_H


#include "optimization.h"
#include <string>
#include "tools.h"

using namespace std;

class PhyloTree;

const double MIN_SITE_RATE = 1e-6;
const double MAX_SITE_RATE = 100.0;
const double TOL_SITE_RATE = 1e-6;


/**
class for among-site rate heterogeneity, the default is homogeneous (equal) rate model

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/

class RateHeterogeneity : public Optimization
{
	friend class ModelFactory;

public:
	/**
		constructor
	*/
    RateHeterogeneity();

	/**
		destructor
	*/
    virtual ~RateHeterogeneity();

	/**
		set phylogenetic tree
		@param tree associated phyogenetic tree
	*/
	void setTree(PhyloTree *tree);

	/**
		set phylogenetic tree
		@param tree associated phyogenetic tree
	*/
	PhyloTree *getTree() { return phylo_tree; }

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams() { return name; }

	/**
		@return false by default. True if rates are site-specific (Meyer and von Haeseler (2003) model)
	*/
	virtual bool isSiteSpecificRate() { return false; }

	/**
		get the number of rate categories. The default returns 1 category since it is homogeneous model
		@return the number of rate categories
	*/
	virtual int getNRate() { return 1; }

	/**
		get the number of rate categories for site-specific category model
		The default returns 1 category since it is homogeneous model
		@return the number of rate categories
	*/
	virtual int getNDiscreteRate() { return 1; }

	/**
		get the rate of a specified category. Default returns 1.0 since it is homogeneous model
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) { return 1.0; }

	/**
		get the proportion of a specified category. Default returns 1.0 since it is homogeneous model
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) { return 1.0; }

	/**
		get the rate of a specified site-pattern. Default returns 1.0 since it is homogeneous model
		@param ptn pattern ID 
		@return the rate of the specified site-pattern
	*/
	virtual double getPtnRate(int ptn) { return 1.0; }

	/**
		get rate category of a specified site-pattern. Default returns -1 since it is homogeneous model
		@param ptn pattern ID 
		@return the rate category of the specified site-pattern
	*/
	virtual int getPtnCat(int ptn) { return -1; }

	/**
		get the proportion of invariable sites. Default returns 0.0 since it is homogeneous model
		@return the proportion of invariable sites
	*/
	virtual double getPInvar() { return 0.0; }

	/**
		get the Gamma shape. Default returns 0.0 since it is homogeneous model
		@return Gamma shape
	*/	
	virtual double getGammaShape() { return 0.0; }

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {}

	/**
		optimize parameters. Default does nothing
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double epsilon) { return 0.0; }

	/**
		classify rates into categories, this is meant for the discrete MH model. 
		The default just return tree_lh
		@param tree_lh the current tree log-likelihood
	*/
	virtual double classifyRates(double tree_lh) { return tree_lh; }

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) {}

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out) {}

	/**
		Compute site-specific rates. Override this for Gamma model
		@param pattern_rates (OUT) pattern rates. Resizing if necesary
	*/
	virtual void computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat) {}

	/**
		write site-rates to a file in the following format:
		1  rate_1
		2  rate_2
		....
		This function will call computePatternRates()
		@param file_name target file to write rates
	*/
	void writeSiteRates(const char *file_name);

	void writeSiteRates(ostream &out);

	void writeSiteRates(ostream &out, DoubleVector &pattern_rates, IntVector &pattern_cat);

	/**
		name of the rate heterogeneity type
	*/
	string name;


	/**
		full name of the rate heterogeneity type
	*/
	string full_name;

	/**
		phylogenetic tree associated
	*/
	PhyloTree *phylo_tree;

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables) {}

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
	*/
	virtual void getVariables(double *variables) {}

	
};
#endif
