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
#ifndef PHYLOSUPERTREE_H
#define PHYLOSUPERTREE_H

#include "iqptree.h"
#include "supernode.h"


struct PartitionInfo {
	string name;
	string model_name;
	string aln_file;
	string sequence_type;
	string position_spec;
};

/**
Phylogenetic tree for partition model (multi-gene alignment)

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class PhyloSuperTree : public IQPTree, public vector<PhyloTree* >
{
public:
	/** 
		constructor
	*/
    PhyloSuperTree();

	/** 
		constructor
	*/
    PhyloSuperTree(Params &params);


    ~PhyloSuperTree();

	virtual bool isSuperTree() { return true; }

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name
            @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = NULL);

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name issued by an interger
            @return a new node
     */
    virtual Node* newNode(int node_id, int node_name);


    /**
            compute the distance between 2 sequences.
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @param initial_dist initial distance
            @return distance between seq1 and seq2
     */
    virtual double computeDist(int seq1, int seq2, double initial_dist);


	/**
		TODO: create sub-trees of the current super-tree
	*/
	void mapTrees();

	void linkTree(int part, NodeVector &part_taxa, SuperNode *node = NULL, SuperNode *dad = NULL);

    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int iterations = 100);

	/**
		partition information
	*/
	vector<PartitionInfo> part_info; 

};

#endif
