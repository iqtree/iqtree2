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

/**
class for rate heterogeneity with a fraction of invariable sites

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateInvar : virtual public RateHeterogeneity
{
public:
	/**
		constructor
		@param p_invar_sites proportion of invariable sites
		@param tree associated phylogenetic tree
	*/
	RateInvar(double p_invar_sites, PhyloTree *tree);


	/**
		get the proportion of invariable sites
		@return the proportion of invariable sites
	*/
	virtual double getPInvar() { return p_invar; }

	/**
		optimize parameters
		@return the best likelihood 
	*/
	virtual double optimizeParameters();

	/**
		override function from Optimization class, used by the minimizeOneDimen() to optimize
		p_invar parameter
	*/
	virtual double computeFunction(double p_invar_value);


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

};

#endif
