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
#ifndef PDTREE_H
#define PDTREE_H

#include "tree/mtree.h"
#include "pda/split.h"


/**
Specialized Tree for Phylogenetic Diversity Algorithms
@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class PDTree: public MTree{
public:
	/**
		construct from program parameters
		@param params program parameters
	*/
    PDTree(Params &params);

	/**
		constructor, get from another tree
		@param tree another MTree
	*/
    PDTree(PDTree &tree);

	/**
		constructor
	*/
	PDTree() : MTree() {};


/********************************************************
	INITIALZATION
********************************************************/
	/**
		initialize the tree from program parameters
		@param params program parameters
	*/
	void init(Params &params);

	/**
		initialize tree, get from another tree
		@param tree another MTree
	*/
	void init(PDTree &tree);

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

	/**
		incoporate the parameters to the tree
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param scale (OUT) branch scaling factor
		@param tax_weight (OUT) taxa weights
	*/
	void incoporateParams(double &scale, DoubleVector &tax_weight, Node *node = NULL, Node* dad = NULL);

	/**
		build a set of leaf name, return to lsn
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param lsn (OUT) leaf set name
	*/
	void buildLeafMapName(LeafMapName &lsn, Node *node = NULL, Node* dad = NULL);

	int findNearestTaxon(Node* &taxon, Node *node, Node *dad = NULL);

/********************************************************
	Computing PD of area (user-defined set of taxa)
********************************************************/

	/**
		compute the PD score of a given taxa set in filename
		@param params program parameters
		@param taxa_set (OUT) corresponding set of taxa
		@param pd_more (OUT) more computed PD measures will be stored here
	*/
	void computePD(Params &params, vector<PDTaxaSet> &taxa_set, PDRelatedMeasures &pd_more);

	/**
		compute the PD score of a given taxa set with name in taxa_name
		@param taxa_name vector of name of all taxa
		@param taxa_set (OUT) corresponding set of taxa
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
	*/
	void makeTaxaSet(set<string> &taxa_name, PDTaxaSet &taxa_set, Node *node = NULL, Node *dad = NULL);

	/**
		compute the PD score of a given taxa set with name in taxa_name
		@param id_set (IN/OUT) corresponding set of taxa
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param curlen current length so far
		@return TRUE if the below subtree contains taxon in id_set
	*/
	bool calcPD(Split &id_set, double curlen = 0.0, Node *node = NULL, Node *dad = NULL);

	/**
		compute the EXCLUSIVE PD score of a given taxa set with name in taxa_name
		@param id_set (IN/OUT) corresponding set of taxa IDs
	*/
	void calcExclusivePD(Split &id_set);

	/**
		compute the area's PD ENDEMISM of set of area
		@param area_set set of area
		@param pd_endem (OUT) corresponding PD endemism
	*/
	void calcPDEndemism(vector<PDTaxaSet> &area_set, DoubleVector &pd_endem);

	/**
		compute the area's PD complementarity given a specific area
		@param area_set set of area
		@param area_name given area names as string separated by commas
		@param pd_comp (OUT) corresponding PD endemism
	*/
	void calcPDComplementarity(vector<PDTaxaSet> &area_set, char *area_name, DoubleVector &pd_comp);


/********************************************************
	VARIABLES
********************************************************/

	/**
	 	initial set of taxa which must be included into the final PD set
	*/
	NodeVector initialset;

protected:


};


#endif
