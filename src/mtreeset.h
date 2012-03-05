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
#ifndef MTREESET_H
#define MTREESET_H

#include "mtree.h"
#include "splitgraph.h"
#include "alignment.h"


/**
Set of trees

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class MTreeSet : public vector<MTree*> {
public:

    MTreeSet();

	/**
		constructor, read trees from user file
		@param userTreeFile the name of the user trees
		@param is_rooted (IN/OUT) true if tree is rooted
		@param burnin the number of beginning trees to be discarded
	*/
	MTreeSet(const char *userTreeFile, bool &is_rooted, int burnin, const char *tree_weight_file = NULL);

	/**
		initialize the tree from a NEWICK tree file
		@param userTreeFile the name of the user tree
		@param is_rooted (IN/OUT) true if tree is rooted
		@param burnin the number of beginning trees to be discarded
	*/
	void init(const char *userTreeFile, bool &is_rooted, int burnin, const char *tree_weight_file = NULL);

	void init(StringIntMap &treels, bool &is_rooted);


	/**
		read the tree from the input file in newick format
		@param userTreeFile the name of the user trees
		@param is_rooted (IN/OUT) true if tree is rooted
		@param burnin the number of beginning trees to be discarded
	*/
	void readTrees(const char *userTreeFile, bool &is_rooted, int burnin);

	/**
		check the consistency of trees: taxa names between trees are matched, same rooted or unrooted
	*/
	void checkConsistency();

	/**
		@return true if trees are rooted
	*/
	bool isRooted();

	/**
		print the tree to the output file in newick format
		@param outfile the output file.
		@param brtype type of branch to print
	*/
	void printTrees(const char *outfile, int brtype = WT_BR_LEN);

	/**
		print the tree to the output file in newick format
		@param out the output stream.
		@param brtype type of branch to print
	*/
	void printTrees(ostream & out, int brtype = WT_BR_LEN);

	/**
		convert all trees into the split system
		@param taxname certain taxa name
		@param sg (OUT) resulting split graph
		@param hash_ss (OUT) hash split set
		@param lensum TRUE if summing split length, FALSE to increment only
		@param weight_threshold minimum weight cutoff
		@param sort_taxa TRUE to sort taxa alphabetically
	*/
	void convertSplits(vector<string> &taxname, SplitGraph &sg, SplitIntMap &hash_ss, 
		int weighting_type, double weight_threshold, bool sort_taxa = true);

	/**
		convert all trees into the split system
		@param sg (OUT) resulting split graph
		@param hash_ss (OUT) hash split set
		@param lensum TRUE to assign split weight as sum of corresponding branch lengths. 
			Otherwise just count the number of branches.
		@param weight_threshold minimum weight cutoff
	*/
	void convertSplits(SplitGraph &sg, SplitIntMap &hash_ss, 
		int weighting_type, double weight_threshold);

	/**
		convert all trees into the split system
		@param sg (OUT) resulting split graph
		@param split_threshold only keep those splits which appear more than this threshold 
		@param lensum TRUE to assign split weight as sum of corresponding branch lengths. 
			Otherwise just count the number of branches.
		@param weight_threshold minimum weight cutoff
	*/
	void convertSplits(SplitGraph &sg, double split_threshold, 
		int weighting_type, double weight_threshold);

	/**
		compute the Robinson-Foulds distance between trees
		@param rfdist (OUT) RF distance
		@param mode 
	*/
	void computeRFDist(int *rfdist, int mode = RF_ALL_PAIR);

	/**
		compute the Robinson-Foulds distance between trees
		@param rfdist (OUT) RF distance
	*/
	void computeRFDist(int *rfdist, MTreeSet *treeset2, 
		const char* info_file = NULL, const char *tree_file = NULL);

	int categorizeDistinctTrees(IntVector &category);

	int sumTreeWeights();

	/**
		destructor
	*/
    virtual ~MTreeSet();

	/**
		new tree allocator
		@return a new tree
	*/
	virtual MTree *newTree() { return new MTree(); }

	IntVector tree_weights;

};

#endif
