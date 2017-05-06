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
#ifndef GREEDY_H
#define GREEDY_H

#include "pdtree.h"

/**
Implementation of greedy algorithm with complexity O(n*logk)
@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class Greedy : public PDTree
{
public:
	/**
		construct from program parameters
		@param params program parameters
	*/
    Greedy(Params &params) : 
		PDTree(params) {}

	/**
		construct from a tree
		@param tree a tree class
	*/
    Greedy(PDTree &tree) : 
		PDTree(tree) {}

	/**
		constructor
	*/
	Greedy() : PDTree() {};

	/**
		run the algorithm
		@param params program parameters
		@param taxa_set (OUT) vector of PD sets
	*/
	void run(Params &params, vector<PDTaxaSet> &taxa_set);

	/**
		update the ordered list based on the recent longest path
		@param node the starting node
		@param subtree (OUT) resulted subtree
		@param cur_set the current set
	*/
	void updateOnLongestPath(Node *node, NodeVector &subtree, PDTaxaSet &cur_set);

	/**
		build the initial subtree based on the initial set of taxa
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param subtree (OUT) resulted subtree
		@param nodestack (TEMP) stack of node, used only by function
	*/
	void buildOnInitialSet(NodeVector &subtree, NodeVector &nodestack, Node *node = NULL, Node *dad = NULL);

	/**
		initialize the ordered list based on the initial subtree structure
		@param subtree vector containing nodes in the subtree
		@return the subtree length
	*/
	double updateOnInitialSet(NodeVector &subtree);

	/**
		@return innodes.begin().
	*/
	NeighborSet::iterator highestNeighbor();

	/**
		add an edge into the NeighborSet
	*/
	void addNeighbor(Neighbor* neigh);

	//NodeSet innodes;

	/**
		neighbor set
	*/
	NeighborSet neighset;

	/**
		list of nodes in the subtree
	*/
	NodeVector subtree;

private:

	/**
		size of list of nodes, used internally during greedy search
	*/
	int list_size;
};

#endif
