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
#ifndef MATREE_H
#define MATREE_H

//#include "tools.h"
#include "mtree.h"
#include "mtreeset.h"

/**
Minh Anh: extended tree to serve some statistics
*/
class MaTree : public MTree
{
public:
/********************************************************
	CONSTRUCTORs, INITIALIZATION AND DESTRUCTORs
********************************************************/

	/**
		constructor, read tree from user file
		@param userTreeFile the name of the user tree
		@param is_rooted (IN/OUT) true if tree is rooted
	*/
	MaTree(const char *userTreeFile, bool &is_rooted) : MTree(userTreeFile, is_rooted) {};

	/**
		constructor, get from another tree
		@param tree another MTree
	*/
	MaTree(MTree &tree) : MTree(tree) {};

	/**
		constructor
	*/
    MaTree() : MTree() {};

/***************************************************************
	OUTPUT INFORMATION ABOUT THE BRANCHES
***************************************************************/
	/**
		Output information about branches on the tree into an output stream
		The information contains:
		- Number of external branches (number of leaves), minimum/maximum/sum of the external branches
		- Number of internal branches, minimum/maximum/sum of the internal branches
		- Number of branches, minimum/maximum/sum of all the branches
		@param out the output stream
	*/
	void printBrInfo(ostream& out);
	
	/**
		Compare this tree with each tree in a given set of trees
		@param trees (IN) the trees to compare
		@param brLenMatrix (OUT) a matrix of double, each row is a vector of double. The size of this vector is the number of branches in the tree, i.e. from 0 to 2n-3 if the tree is unrooted or from 0 to 2n-2 if the tree is rooted.
		If branch i is contained in the other tree (atree), element i is the length of this branch 
		on the other tree. If not, element i is set to -1.
		@param RFs (OUT) the Robinson-Foulds distance between this tree and each of the given trees.
		@param BSDs (OUT) the branch score distance between this tree and each of the given trees.
	*/
	void comparedTo(MTreeSet &trees, DoubleMatrix &brLenMatrix, IntVector &RFs, DoubleVector &BSDs);

	/**
		convert the tree into SplitIntMap, the integer number is the nodeID of the corresponding branch
		@param sim (OUT) resulting splitIntMap
		@param taxonID (IN) the ID of an external node (taxon) to be presented in all splits
	*/
//	void convertSplitIntMap(SplitIntMap &sim, const int taxonID);

	/**
		convert the tree into SplitIntMap, iterative procedure
		@param sim (OUT) resulting splitIntMap
		@param resp (internal) set of taxa below node
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param taxonID (IN) the ID of an external node (taxon) to be presented in all splits
	*/
	void convertSplitIntMap(SplitIntMap &sg, Split *resp, const int taxonID, Node *node = NULL, Node *dad = NULL);
};

#endif
