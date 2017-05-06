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
#ifndef PDTREESET_H
#define PDTREESET_H

#include "tree/mtreeset.h"
#include "pdtree.h"

/**
Vector of PDTree

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class PDTreeSet : public MTreeSet
{
public:
    PDTreeSet();

	/**
		constructor, read trees from user file
		@param params program parameters
	*/
	PDTreeSet(Params &params);

	/**
		@return a new tree
	*/
	virtual MTree *newTree() { return  (new PDTree()); }

	/**
		@return true if trees are rooted
	*/
	bool isRootedTrees();

	/**
		@return number of taxa
	*/
	int getNTaxa();


/********************************************************
	INITIALZATION
********************************************************/
	/**
		initialize the tree from program parameters
		@param params program parameters
	*/
	void init(Params &params);

	/**
		read the parameter from the file
		@param params program parameters
	*/
	void readParams(Params &params);

	/**
		Identify the root node if specified, include it into the initial set
		@param root_name name of the root node
	*/
	void readRootNode(const char *root_name);

	/**
		read the initial set of taxa to be included into PD-tree
		@param params program parameters
	*/
	void readInitialSet(Params &params);

protected:
	
	/**
		name of initial taxa, to be included into PD set
	*/
	StrVector init_taxa;

};

#endif
