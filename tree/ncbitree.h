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
#ifndef NCBITREE_H
#define NCBITREE_H

#include "mtree.h"

const int MAX_TAXONOMY_ID = 2000000;

/**
Class for processing NCBI Taxonomy tree

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class NCBITree : public MTree
{
public:
    NCBITree();
	
    ~NCBITree();


/********************************************************
	READ TREE FROM FILE
********************************************************/

	/**
		read the tree in nodes.dmp file from NCBI taxonomy
		@param infile the input file file.
		@param root_id taxon ID of the root
		@param taxon_level e.g. "species", "genus"; NULL to take all taxa (incl. subspecies)
	*/
	Node* readNCBITree(const char *infile, int root_id, const char* taxon_level, const char *ignore_level); 
	Node* readNCBITree(istream &in, int root_id, const char* taxon_level, const char *ignore_level);

	/**
		read names.dmp file. You must call readNCBITree() before calling this function
		@param infile input file name (typically names.dmp from NCBI)
		@param name_type type of the name
	*/
	void readNCBINames(const char* infile, const char *name_type = "scientific name");
	void readNCBINames(ifstream &in, const char *name_type = "scientific name");


protected:

	/**
		taxonomy rank of the nodes
	*/
	StrVector node_levels;

	/**
		vector of all taxonomical nodes
	*/
	NodeVector nodes;

	/**
		prune subtree below the taxon_level
		@return number of nodes pruned
	*/
	int pruneTaxa(StrVector &node_levels, const char* taxon_level, Node *node, Node *dad);

	/**
		prune all nodes that have degree of 2
		@return number of nodes pruned
	*/
	int pruneBridgeNodes(Node *node, Node *dad);

	void countNodeNum(Node *node, Node *dad);

	/**
		release the nemory.
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
	*/
	int freeNode(Node *node = NULL, Node *dad = NULL);

};

#endif
