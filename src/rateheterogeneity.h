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
using namespace std;

class PhyloTree;

const double MIN_SITE_RATE = 1e-6;
const double MAX_SITE_RATE = 200.0;
const double TOL_SITE_RATE = 1e-6;


/**
class for among-site rate heterogeneity, the default is homogeneous (equal) rate model

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/

class RateHeterogeneity : public Optimization
{
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
		@return false by default. True if rates are site-specific (Meyer and von Haeseler (2003) model)
	*/
	virtual bool isSiteSpecificRate() { return false; }

	/**
		get the number of rate categories. The default returns 1 category since it is homogeneous model
		@return the number of rate categories
	*/
	virtual int getNRate() { return 1; }

	/**
		get the rate of a specified category. Default returns 1.0 since it is homogeneous model
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) { return 1.0; }

	/**
		get the proportion of invariable sites. Default returns 0.0 since it is homogeneous model
		@return the proportion of invariable sites
	*/
	virtual double getPInvar() { return 0.0; }

	virtual double getGammaShape() { return 0.0; }

	/**
		optimize parameters. Default does nothing
		@return the best likelihood 
	*/
	virtual double optimizeParameters() { return 0.0; }

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
		name of the rate heterogeneity type
	*/
	string name;


	/**
		full name of the rate heterogeneity type
	*/
	string full_name;

protected:

	/**
		phylogenetic tree associated
	*/
	PhyloTree *phylo_tree;
	

};
#endif
