/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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
#ifndef MPDABLOCK_H
#define MPDABLOCK_H

#include "ncl/ncl.h"
#include "utils/tools.h"

class SplitGraph;

/**
PdaBlock to read from nexus file

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class MPdaBlock : public NxsBlock
{
public:
	friend class SplitGraph;
	friend class PDNetwork;
	friend class CircularNetwork;

	/**
		constructor, assigning an associated splits graph
		@param asgraph a splits graph
	*/
    MPdaBlock(SplitGraph *asgraph);

	/**
		destructor
	*/
    virtual ~MPdaBlock();

	/**
		print info to an output stream
		@param out output stream, cout for output to screen
	*/
	virtual void Report(ostream &out);

	/**
		reset the block
	*/
	virtual void Reset();

	/**
		called when some commands are skipped
		@param commandName command name
	*/
	virtual void SkippingCommand(NxsString commandName);


	/**
		read the file containing total budget and taxa costs informations
		@param params program parameters
	*/
	void readBudgetFile(Params &params);

	/**
		read the file containing total budget and area costs informations
		@param params program parameters
	*/
	void readBudgetAreaFile(Params &params);


	/**
		@return total budget
	*/
	double getBudget() {
		return budget;
	}

	/**
		@return min budget
	*/
	double getMinBudget() {
		return min_budget;
	}

	/**
		@return size of PD set
	*/
	int getSubSize() {
		return sub_size;
	}

	/**
		@return cost of a taxon
	*/
	double getCost(int tax_id) {
		ASSERT(tax_id < (int) costs.size());
		return costs[tax_id];
	}


protected:

	/**
		the associated splits graph
	*/
	SplitGraph *sgraph;

	/**
		total budget
	*/
	double budget;

	/**
		min budget, to compute PD sets with preservation
		costs from min_budget to budget
	*/
	double min_budget;

	/**
		size of PD set
	*/
	int sub_size;

	/**
		true if cost constrained PD problem
	*/
	bool cost_constrained;

	/**
		cost of each taxon
	*/
	vector<double> costs;

	/**
		main method to read block from file
		@param token a token reader
	*/
	virtual void Read(NxsToken &token);

};

#endif
