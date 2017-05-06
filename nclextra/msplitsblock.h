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
#ifndef MSPLITSBLOCK_H
#define MSPLITSBLOCK_H

#include "ncl/ncl.h"
//#include "splitgraph.h"

class SplitGraph;

/**
SplitsBlock to read from nexus file

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class MSplitsBlock : public NxsBlock
{
public:

	friend class SplitGraph;
	friend class MTree;

	/**
		constructor, assigning an associated splits graph
		@param asgraph a splits graph
	*/
    MSplitsBlock(SplitGraph *asgraph);

	/**
		destructor
	*/
	virtual ~MSplitsBlock();


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
		parse a line containing split
		@param token a token reader
	*/
	void AddSplit(NxsToken &token);

	/**
		@return cycle
	*/
	inline vector<int> &getCycle() {
		return cycle;
	}

protected:

	/**
		number of taxa
	*/
	int ntaxa;

	/**
		number of splits
	*/
	int nsplits;

	/**
		the associated splits graph
	*/
	SplitGraph *sgraph;

	/**
		taxa index around circle, if it is a circular split graph
	*/
	vector<int> cycle;

	/**
		main method to read block from file
		@param token a token reader
	*/
	virtual void Read(NxsToken &token);

};

#endif
