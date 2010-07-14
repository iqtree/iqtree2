//
// C++ Interface: phylotree
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include "mtree.h"
#include "alignment.h"
#include "substmodel.h"
#include "phylonode.h"
#include "optimization.h"
#include "rateheterogeneity.h"

const int MAX_SPR_MOVES = 20;

/**
	an SPR move.
*/
struct SPRMove {
	PhyloNode *prune_dad;
	PhyloNode *prune_node;
	PhyloNode *regraft_dad;
	PhyloNode *regraft_node;
	double score;
};

struct SPR_compare
{
  bool operator()(SPRMove s1, SPRMove s2) const
  {
    return s1.score > s2.score;
  }
};

class SPRMoves : public set<SPRMove, SPR_compare> {
public:
	void add(PhyloNode *prune_node, PhyloNode *prune_dad, 
		PhyloNode *regraft_node, PhyloNode *regraft_dad, double score);
};


/*
left_node-----------dad-----------right_node
                     |
                     |
                     |
                    node
*/
struct PruningInfo {
	NeighborVec::iterator dad_it_left, dad_it_right, left_it, right_it;
	Neighbor *dad_nei_left, *dad_nei_right, *left_nei, *right_nei;
	Node *node, *dad, *left_node, *right_node;
	double left_len, right_len;
	double *dad_lh_left, *dad_lh_right;
	
};


/**
Phylogenetic Tree class

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class PhyloTree : public MTree, public Optimization
{
public:

	/**
		constructor
	*/
	PhyloTree();

	/**
		destructor
	*/
	virtual ~PhyloTree();


	/**
		copy the phylogenetic tree structure into this tree, override to take sequence names
		in the alignment into account
		@param tree the tree to copy
	*/
	virtual void copyTree(MTree *tree);

	
	/**
		copy the phylogenetic tree structure into this tree, designed specifically for PhyloTree. 
		So there is some distinction with copyTree.
		@param tree the tree to copy
	*/
	void copyPhyloTree(PhyloTree *tree);


	/**
		set the alignment, important to compute parsimony or likelihood score
		@param alignment associated alignment
	*/
    void setAlignment(Alignment *alignment);

	/**
		set the substitution model, important to compute the likelihood
		@param amodel associated substitution model
	*/
	void setModel(SubstModel *amodel);

	/**
		set rate heterogeneity, important to compute the likelihood
		@param rate associated rate heterogeneity class
	*/
	void setRate(RateHeterogeneity *rate);

	/**
		create substitution model with possible rate heterogeneity. Create proper class objects
		for two variables: model and site_rate. It takes the following field of params into account:
			model_name, num_rate_cats, freq_type
		@param params program parameters
	*/
	void createModel(Params &params);

	/**
		get the name of the model
	*/
	string getModelName();

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
		this function return the parsimony or likelihood score of the tree. Default is 
		to compute the parsimony score. Override this function if you define a new
		score function.
		@return the tree score
	*/
	//virtual double computeScore() { return -computeLikelihood(); }
	//virtual double computeScore() { return (double)computeParsimonyScore(); }

/****************************************************************************
	Parsimony function
****************************************************************************/

	/**
		compute the parsimony score of the tree, given the alignment
		@return the parsimony score
	*/
	int computeParsimonyScore();


	/**
		compute the parsimony score of the tree, given the alignment
		@return the parsimony score
		@param node the current node
		@param dad dad of the node, used to direct the search
		@param ptn pattern ID
		@param states set of admissible states at the current node (in binary code)
	*/
	int computeParsimonyScore(int ptn, int &states, PhyloNode *node = NULL, PhyloNode *dad = NULL);


/****************************************************************************
	likelihood function
****************************************************************************/

	/**
		initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh 
	*/
	void initializeAllPartialLh();

	/**
		initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh 
		@param node the current node
		@param dad dad of the node, used to direct the search
		@param index the index 
	*/
	void initializeAllPartialLh(int &index, PhyloNode *node = NULL, PhyloNode *dad = NULL);


	/**
		clear all partial likelihood for a clean computation again
	*/
	void clearAllPartialLh();

	/**
		allocate memory for a partial likelihood vector
	*/
	double *newPartialLh();

	/**
		compute the partial likelihood at a subtree
		@param dad_branch the branch leading to the subtree
		@param dad its dad, used to direct the tranversal
	*/
	void computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL);

	/**
		compute tree likelihood on a branch. used to optimize branch length
		@param dad_branch the branch leading to the subtree
		@param dad its dad, used to direct the tranversal	
		@return tree likelihood
	*/
	double computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);


	/**
		compute the tree likelihood 
		@return tree likelihood
	*/
	double computeLikelihood();

	/**
		optimize model parameters and tree branch lengths
		@param fixed_len TRUE to fix branch lengths, default is false
		@return the best likelihood 
	*/
	virtual double optimizeModel(bool fixed_len = false);

/****************************************************************************
	computing derivatives of likelihood function
****************************************************************************/
	/**
		compute tree likelihood and derivatives on a branch. used to optimize branch length
		@param dad_branch the branch leading to the subtree
		@param dad its dad, used to direct the tranversal	
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
		@return tree likelihood
	*/
	double computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);


/****************************************************************************
	Stepwise addition (greedy) by maximum parsimony
****************************************************************************/

	/**
		grow the tree by step-wise addition
		@param alignment input alignment
	*/
	void growTreeMP(Alignment *alignment);

	/**
		used internally by growTreeMP() to find the best target branch to add into the tree
		@param added_node node to add
		@param target_node (OUT) one end of the best branch found
		@param target_dad (OUT) the other end of the best branch found
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the parsimony score of the tree
	*/
	int  addTaxonMP(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad);


/****************************************************************************
	Nearest Neighbor Interchange with parsimony
****************************************************************************/
	/**
		search by a nearest neigbor interchange with parsimony 
	*/
	void searchNNI();

	/**
		search by a nearest neigbor interchange with parsimony
		@param node the current node
		@param dad dad of the node, used to direct the search
		@param cur_score current score
		@return best score
	*/
	double searchNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

	/**
		try to swap the tree with nearest neigbor interchange at the branch connecting node1-node2.
		If a swap shows better score, return the swapped tree and the score.
		@param cur_score current score
		@param node1 1st end node of the branch
		@param node2 2nd end node of the branch
		@return best score
	*/
	double swapNNI(double cur_score, PhyloNode *node1, PhyloNode *node2);

/****************************************************************************
	Branch length optimization by maximum likelihood
****************************************************************************/

	/**
		optimize one branch length by ML
		@param node1 1st end node of the branch
		@param node2 2nd end node of the branch
                @param clearLH true to clear the partial likelihood, otherwise false
		@return likelihood score
	*/
	double optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH=true);

	/**
		optimize all branch lengths of the children of node
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the likelihood of the tree
	*/
	double optimizeChildBranches(PhyloNode *node, PhyloNode *dad = NULL);

	/**
		optimize all branch lengths at the subtree rooted at node step-by-step.
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the likelihood of the tree
	*/
	double optimizeAllBranches(PhyloNode *node, PhyloNode *dad = NULL);

	/**
		optimize all branch lengths of the tree
		@param iterations number of iterations to loop through all branches
		@return the likelihood of the tree
	*/
	double optimizeAllBranches(int iterations = 100);

	/**
		inherited from Optimization class, to return to likelihood of the tree
		when the current branch length is set to value
		@param value current branch length
		@return negative of likelihood (for minimization)
	*/
	virtual double computeFunction(double value);

	/**
		Inherited from Optimization class.
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		used by Newton raphson method to minimize the function.
		@param value current branch length
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
		@return negative of likelihood (for minimization)
	*/
	virtual double computeFuncDerv(double value, double &df, double &ddf);


/****************************************************************************
	Nearest Neighbor Interchange by maximum likelihood
****************************************************************************/

	/**
		search by a nearest neigbor interchange, then optimize branch lengths. Do it 
		until tree does not improve
		@return the likelihood of the tree	
	*/
	double optimizeNNIBranches();

	/**
		search by a nearest neigbor interchange
		@return the likelihood of the tree
	*/
	virtual double optimizeNNI();

	/**
		search by a nearest neigbor interchange 
		@param cur_score current likelihood score
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the likelihood of the tree
	*/
	double optimizeNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

	/**
		This is for ML. try to swap the tree with nearest neigbor interchange at the branch connecting node1-node2.
		If a swap shows better score, return the swapped tree and the score.
		@param cur_score current likelihood score
		@param node1 1st end node of the branch
		@param node2 2nd end node of the branch
		@return the likelihood of the tree
	*/
	double swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2);


/****************************************************************************
	Stepwise addition (greedy) by maximum likelihood
****************************************************************************/

	/**
		grow the tree by step-wise addition
		@param alignment input alignment
	*/
	void growTreeML(Alignment *alignment);

	/**
		used internally by growTreeML() to find the best target branch to add into the tree
		@param added_node node to add
		@param target_node (OUT) one end of the best branch found
		@param target_dad (OUT) the other end of the best branch found
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the likelihood of the tree
	*/
	double addTaxonML(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad);

/****************************************************************************
	compute BioNJ tree, a more accurate extension of Neighbor-Joining
****************************************************************************/

	/**
		compute BioNJ tree
		@param params program parameters
		@param alignment input alignment
		@param dist_mat (OUT) distance matrix
	*/
	void computeBioNJ(Params &params, Alignment *alignment, double* &dist_mat);

	/**
		Neighbor-joining tree might contain negative branch length. This
		function will fix this.
		@param fixed_length fixed branch length to set to negative branch lengths
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return The number of branches that have no/negative length
	*/
	int fixNegativeBranch(double fixed_length, Node *node = NULL, Node *dad = NULL);


/****************************************************************************
	Subtree Pruning and Regrafting by maximum likelihood
	NOTE: NOT DONE YET
****************************************************************************/
	
	/**
		search by Subtree pruning and regrafting
		@return the likelihood of the tree		
	*/
	double optimizeSPR();

	/**
		search by Subtree pruning and regrafting, then optimize branch lengths. Iterative until 
		no tree improvement found.
		@return the likelihood of the tree		
	*/
	double optimizeSPRBranches();

	/**
		search by Subtree pruning and regrafting at a current subtree
		@param cur_score current likelihood score
		@param node the current node
		@param dad dad of the node, used to direct the search
		@return the likelihood of the tree		
	*/
	double optimizeSPR(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

	/**
		move the subtree (dad1-node1) to the branch (dad2-node2)
	*/
	double swapSPR(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1, 
		PhyloNode *orig_node1, PhyloNode *orig_node2,
		PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor*> &spr_path);

	double assessSPRMove(double cur_score, const SPRMove &spr);

	void pruneSubtree(PhyloNode *node, PhyloNode *dad, PruningInfo &info);

	void regraftSubtree(PruningInfo &info, 
		PhyloNode *in_node, PhyloNode *in_dad);

	/**
		associated alignment	
	*/
	Alignment *aln;

	/**
		TRUE if you want to optimize branch lengths by Newton-Raphson method
	*/
	bool optimize_by_newton;

protected:

	/**
		assign the leaf names with the alignment sequence names, using the leaf ID for assignment.
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
	*/
	void assignLeafNames(Node *node = NULL, Node *dad = NULL);


	/**
		associated substitution model
	*/
	SubstModel *model;

	/**
		among-site rates 
	*/
	RateHeterogeneity *site_rate;

	/**
		current branch iterator, used by computeFunction() to optimize branch lengths
	*/
	PhyloNeighbor *current_it;
	/**
		current branch iterator of the other end, used by computeFunction() to optimize branch lengths
	*/
	PhyloNeighbor *current_it_back;

	/**
		spr moves
	*/
	SPRMoves spr_moves;

	/**
		SPR radius
	*/
	int spr_radius;


	/**
		the main memory storing all partial likelihoods for all neighbors of the tree.
		The variable partial_lh in PhyloNeighbor will be assigned to a region inside this variable.
	*/
	double *central_partial_lh;
};

#endif
