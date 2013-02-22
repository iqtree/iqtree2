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

#include "iqtree.h"
#include "supernode.h"
#include "superalignment.h"


struct PartitionInfo {
	string name; // partition name
	string model_name; // model name
	string aln_file; // alignment file associated
	string sequence_type; // sequence type (DNA/AA/BIN)
	string position_spec; // position specification, e.g., "1-100\1 1-100\2"

	double cur_score; // current log-likelihood 

	DoubleVector null_score; // log-likelihood of each branch collapsed to zero
	DoubleVector opt_score; // optimized log-likelihood for every branch
	DoubleVector nni1_score; // log-likelihood for 1st NNI for every branch
	DoubleVector nni2_score; // log-likelihood for 2nd NNI for every branch

	DoubleVector cur_brlen; // current branch lengths
	DoubleVector opt_brlen; // optimized branch lengths for every branch
	DoubleVector nni1_brlen; // branch length for 1st NNI for every branch
	DoubleVector nni2_brlen; // branch length for 2nd NNI for every branch
};

/**
Phylogenetic tree for partition model (multi-gene alignment)

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class PhyloSuperTree : public IQTree, public vector<PhyloTree* >
{
public:
	/** 
		constructor
	*/
    PhyloSuperTree();

	/** 
		constructor
	*/
    PhyloSuperTree(SuperAlignment *alignment, PhyloSuperTree *super_tree);

	/** 
		constructor
	*/
    PhyloSuperTree(Params &params);


    ~PhyloSuperTree();

		/**
	 * setup all necessary parameters  (declared as virtual needed for phylosupertree)
	 */
	virtual void setParams(Params& params);
	
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
	 *		@return number of alignment patterns
	*/
	virtual int getAlnNPattern();

	/**
	 *		@return number of alignment sites
	*/
	virtual int getAlnNSite();

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

	void linkTrees();

	void linkBranch(int part, SuperNeighbor *nei, SuperNeighbor *dad_nei);

    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
            compute pattern likelihoods only if the accumulated scaling factor is non-zero.
            Otherwise, copy the pattern_lh attribute
            @param pattern_lh (OUT) pattern log-likelihoods,
                            assuming pattern_lh has the size of the number of patterns
     */
	virtual void computePatternLikelihood(double *pattern_lh, double *cur_logl = NULL);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD);

    /**
            optimize one branch length by ML
            @param node1 1st end node of the branch
            @param node2 2nd end node of the branch
            @param clearLH true to clear the partial likelihood, otherwise false
            @return likelihood score
     */
    virtual double optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true);

    /**
            search the best swap for a branch
            @return NNIMove The best Move/Swap
            @param cur_score the current score of the tree before the swaps
            @param node1 1 of the 2 nodes on the branch
            @param node2 1 of the 2 nodes on the branch
     */
    virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, bool approx_nni, double lh_contribution = -1.0);

    /**
            Do an NNI
     */
    virtual double doNNI(NNIMove move);

    /**
     * 	 Restore the branch lengths from the saved values
     */
    virtual void restoreAllBranLen(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            reinsert the whole list of leaves back into the tree
            @param del_leaves the list of deleted leaves, returned by deleteLeaves() function
     */
    virtual void reinsertLeaves(PhyloNodeVector &del_leaves);

	/**
		compute the weighted average of branch lengths over partitions
	*/
	void computeBranchLengths();

	void printMapInfo();
	
	/**
	 * initialize partition information for super tree
	*/
	virtual void initPartitionInfo();

	/**
		partition information
	*/
	vector<PartitionInfo> part_info; 

	

};

#endif
