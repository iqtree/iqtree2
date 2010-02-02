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
#ifndef IQPTREE_H
#define IQPTREE_H

#include "phylotree.h"
#include "phylonode.h"
#include <set>


/**
	nodeheightcmp, for building k-representative leaf set
*/
struct nodeheightcmp
{
/**
	nodeheightcmp, for building k-representative leaf set
*/

  bool operator()(const Node* s1, const Node* s2) const
  {
    return (s1->height) < (s2->height);
  }
};

/**
	NNISwap, define a NNI Swap or Move
*/
struct NNIMove
{
	PhyloNode *node1;
	Neighbor *node1Nei;	
	PhyloNode *node2;
	Neighbor *node2Nei;
	double score;
	
	bool operator<(const NNIMove& rhs) const
	{ 
		return score > rhs.score; 
	}
	
};

/**
	Representative Leaf Set, stored as a multiset template of STL, 
	sorted in ascending order of leaf's height
*/
typedef multiset<Node*, nodeheightcmp> RepresentLeafSet;


/**
Important Quartet Puzzling

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class IQPTree : public PhyloTree
{
public:
	/**
		constructor
	*/
    IQPTree();

	/**
		destructor
	*/
    virtual ~IQPTree();

	/**
		set k-representative parameter
		@param k_rep k-representative
	*/
    void setRepresentNum(int k_rep);

	/**
		set the probability of deleteing sequences for IQP algorithm
		@param p_del probability of deleting sequences
	*/
    void setProbDelete(double p_del);

	/**
		set the number of iterations for the IQPNNI algorithm
		@param iterations the number of iterations
	*/
	void setIQPIterations(int iterations);

	/**
		find the k-representative leaves under the node
		@param node the node at which the subtree is rooted
		@param dad the dad node of the considered subtree, to direct the search
		@param leaves (OUT) the k-representative leaf set
	*/
	void findRepresentLeaves(RepresentLeafSet &leaves, PhyloNode *node = NULL, PhyloNode *dad = NULL);

	/**
		perform one IQPNNI iteration
		@return current likelihood
	*/
	double doIQP();

	/**
		perform all IQPNNI iterations
		@return best likelihood found
		@param tree_file_name name of the tree file to write the best tree found
	*/
	double doIQPNNI(string tree_file_name);

/****************************************************************************
	Fast Nearest Neighbor Interchange by maximum likelihood
****************************************************************************/

	/**
		This implement the fastNNI algorithm proposed in PHYML paper
		TUNG: this is a virtual function, so it will be called automatically by optimizeNNIBranches()
		@return best likelihood found
	*/
	virtual double optimizeNNI();
	
	/*
	 * 	Do Simple NNI (Slow NNI)
	 */
	double optimizeNNISimple();

	/**
		search all positive NNI move on the current tree and save them on the possilbleNNIMoves list
	*/
	void generateAllPositiveNNIMoves(double cur_score,PhyloNode *node = NULL, PhyloNode *dad = NULL);
    
	
	/**
		search the best swap for a branch 
		@return NNIMove The best Move/Swap
		@param cur_score the current score of the tree before the swaps
		@param node1 1 of the 2 nodes on the branch
		@param node2 1 of the 2 nodes on the branch
	*/	
	NNIMove getBestNNIMoveForBranch(double cur_score, PhyloNode *node1, PhyloNode *node2); 
	
	/**
		distance matrix, used for IQP algorithm
	*/
	double *dist_matrix;
	
	/**
		add a NNI move to the list of possible NNI moves;
	*/
	void addPossibleNNIMove(NNIMove myMove);
	
	/**
		Described in PhyML paper : apply change to branches that not correspond to a swap  l = l + lamda(la - l)
	*/
	void applyBranchLengthChanges();
	
	/**
		Do an NNI
	*/
	void swapNNIBranch(PhyloNode *node1, PhyloNode *node2);

protected:
  
  	/**
		The list of possible NNI moves for the current tree;
	*/
	vector<NNIMove> possibleNNIMoves;
	
	/**
		k-representative parameter
	*/
	int k_represent;

	/**
		probability to delete a leaf
	*/
	double p_delete;

	/**
		number of IQPNNI iterations
	*/
	int iqpnni_iterations;

	/**
		bonus values of all branches, used for IQP algorithm
	*/
	double *bonus_values;

	/**
		delete a leaf from the tree, assume tree is birfucating
		@param leaf the leaf node to remove
	*/
	void deleteLeaf(Node *leaf);

	/**
		delete a set of leaves from tree (with the probability p_delete), assume tree is birfucating
		@param del_leaves (OUT) the list of deleted leaves
		@param adjacent_nodes (OUT) the corresponding list of nodes adjacent to the deleted leaves
	*/
	void deleteLeaves(PhyloNodeVector &del_leaves, PhyloNodeVector &adjacent_nodes);

	/**
		reinsert one leaf back into the tree
		@param leaf the leaf to reinsert
		@param adjacent_node the node adjacent to the leaf, returned by deleteLeaves() funcrion
		@param node one end node of the reinsertion branch in the existing tree
		@param dad the other node of the reinsertion branch in the existing tree
	*/
	void reinsertLeaf(Node *leaf, Node *adjacent_node, Node *node, Node *dad);
		
	/**
		reinsert the whole list of leaves back into the tree
		@param del_leaves the list of deleted leaves, returned by deleteLeaves() function
		@param adjacent_nodes the corresponding list of nodes adjacent to the deleted leaves, returned by deleteLeaves() function
	*/
	void reinsertLeaves(PhyloNodeVector &del_leaves, PhyloNodeVector &adjacent_nodes);

	/**
		assess a quartet with four taxa. Current implementation uses the four-point condition 
		based on distance matrix for quick evaluation.
		@param leaf0 one of the leaf in the existing sub-tree
		@param leaf1 one of the leaf in the existing sub-tree
		@param leaf2 one of the leaf in the existing sub-tree
		@param del_leaf a leaf that was deleted (not in the existing sub-tree)
	*/
	int assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf);

	/**
		assess the important quartets around a virtual root of the tree. 
		This function will assign bonus points to branches by updating the variable 'bonus_values'
		@param cur_root the current virtual root
		@param del_leaf a leaf that was deleted (not in the existing sub-tree)
	*/
	void assessQuartets(PhyloNode *cur_root, PhyloNode *del_leaf);

	/**
		initialize the bonus points to ZERO
	*/
	void initializeBonus();

	/**
		raise the bonus points for all branches in the subtree rooted at a node
		@param node the root of the sub-tree
		@param dad dad of 'node', used to direct the recursion
	*/
	void raiseBonus(Node *node, Node *dad);

	/**
		find the branch with the best bonus points
		@param node the root of the sub-tree
		@param dad dad of 'node', used to direct the recursion
		@param best_node (OUT) one end of the branch with highest bonus point
		@param best_dad (OUT) the other end of the branch with highest bonus point
	*/
	double findBestBonus(Node *node, Node *dad, Node *&best_node, Node *&best_dad);

};

#endif
