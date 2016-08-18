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
#ifndef PRUNING_H
#define PRUNING_H

#include "pdtree.h"

/**
Implementation of Pruning algorithm with complexity O(n*log(n-k))

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class Pruning : public PDTree
{
public:
	/**
		construct from program parameters
		@param params program parameters
	*/
    Pruning(Params &params) : 
		PDTree(params) {}

	/**
		constructor, get from another tree
		@param tree another MTree
	*/
    Pruning(PDTree &tree) : 
		PDTree(tree) {}

	/**
		constructor
	*/
	Pruning() : PDTree() {};

	/**
		run the algorithm
		@param params program parameters
		@param taxa_set (OUT) vector of PD sets
	*/
	void run(Params &params, vector<PDTaxaSet> &taxa_set);

	/**
		delete an external node 
		@param pos the position of the node in the LeafSet
	*/
	void deleteExNode(LeafSet::iterator pos);

	/**
		build the list of all leaves into field leaves.
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
	*/
	void buildLeaves(Node *node = NULL, Node *dad = NULL);

	/**
		print all leaves into screen
	*/
	void printLeaves();

	/**
		find the iterator to the leaf node in leaves field
		@param node a leaf node.
		@return iterator to the leaf node in leaves field
	*/
	LeafSet::iterator findNode(Node *node);

	/**
		@return leaves.begin().
	*/
	LeafSet::iterator nearestLeaf();

	/**
		mark the node in the initial set to be not PRUNABLE
	*/
	void doInitialSet();

	/**
		insert a leaf into the LeafSet
		@param leaf leaf node to be inserted
	*/
	void addLeaf(Node* leaf);

	/**
		leaf set of the tree, used for pruning algorithm
	*/
	LeafSet leaves;

	/**
		maximum size of the ordered list of leaves
	*/
	int list_size;

};

#endif
