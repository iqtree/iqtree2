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
#include <map>
#include <stack>
#include <vector>
#include "stoprule.h"


#define SPEED_UP_AFTER_ITER_NUM 100

typedef std::map< string, double > BranLenMap;
typedef std::multiset< double, std::less< double > > multiSetDB;
typedef std::multiset< int, std::less< int > > MultiSetInt;

/**
        nodeheightcmp, for building k-representative leaf set
 */


class RepLeaf {
public:
    Node *leaf;
    int height;

    RepLeaf(Node *aleaf, int aheight = 0) {
        leaf = aleaf;
        height = aheight;
    }
};

struct nodeheightcmp {

    bool operator()(const RepLeaf* s1, const RepLeaf * s2) const {
        return (s1->height) < (s2->height);
    }
};

/**
        Representative Leaf Set, stored as a multiset template of STL,
        sorted in ascending order of leaf's height
 */
typedef multiset<RepLeaf*, nodeheightcmp> RepresentLeafSet;

/**
Important Quartet Puzzling

        @author BUI Quang Minh <minh.bui@univie.ac.at>
 */
class IQPTree : public PhyloTree {
public:
    /**
            default constructor
     */
    IQPTree();

    IQPTree(Alignment *aln);
	void setParams(Params &params);

    /**
            destructor
     */
    virtual ~IQPTree();

    /**
            print tree to .treefile
            @param params program parameters, field root is taken
     */
    void printResultTree(Params &params);

    /**
            print tree to out
            @param params program parameters, field root is taken
            @param out (OUT) output stream
     */
    void printResultTree(Params &params, ostream &out);

	/**
		print intermediate tree
	*/
	void printIntermediateTree(int brtype, Params &params);

    void setRootNode(char *my_root);

    /**
            set k-representative parameter
            @param k_rep k-representative
     */
   // void setRepresentNum(int k_rep);

    /**
            set the probability of deleteing sequences for IQP algorithm
            @param p_del probability of deleting sequences
     */
    //void setProbDelete(double p_del);

    double getProbDelete();

    void resetKDelete();
    void increaseKDelete();

    /**
            set the number of iterations for the IQPNNI algorithm
            @param stop_condition stop condition (SC_FIXED_ITERATION, SC_STOP_PREDICT)
            @param min_iterations the min number of iterations
            @param max_iterations the maximum number of iterations
     */
    void setIQPIterations(STOP_CONDITION stop_condition, double stop_confidence, int min_iterations, int max_iterations);

    /**
            @param assess_quartet the quartet assessment, either IQP_DISTANCE or IQP_PARSIMONY
     */
    //void setIQPAssessQuartet(IQP_ASSESS_QUARTET assess_quartet);

    /**
            find the k-representative leaves under the node
            @param node the node at which the subtree is rooted
            @param dad the dad node of the considered subtree, to direct the search
            @param leaves (OUT) the k-representative leaf set
     */
    RepresentLeafSet* findRepresentLeaves(vector<RepresentLeafSet*> &leaves, int nei_id,
            PhyloNode *dad);

    /**
            clear representative leave sets iteratively, called once a leaf is re-inserted into the tree
            @param node the node at which the subtree is rooted
            @param dad the dad node of the considered subtree, to direct the search
            @param leaves (OUT) the k-representative leaf set
     */
    void clearRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, Node *node, Node *dad);

    /**
            perform one IQPNNI iteration
            @return current likelihood
     */
    double doIQP();

    /**
     * Perturb the tree for the next round of local search by swaping position of 2 random leaves
     * @param nbDist The minimum distance between the 2 nodes that are swapped
     * @param nbTimes Number of times that the swap operations are carried out
     * @return The new loglikelihood of the tree
     */
    double perturb(int times);

    /**
     * Carry out Iterated Local Search
     * @param numIter number of iteration
     * @param perturbLevel the level of perturbation
     * @return tree's score
     */
    double doILS(Params &params, int perturbLevel);

    /**
     * TODO
     * @param node1
     * @param node2
     * @return
     */
    double swapTaxa(PhyloNode *node1, PhyloNode *node2);

    /**
            perform all IQPNNI iterations
            @return best likelihood found
            @param tree_file_name name of the tree file to write the best tree found
     */
    double doIQPNNI(Params &params);

    /****************************************************************************
            Fast Nearest Neighbor Interchange by maximum likelihood
     ****************************************************************************/

    /**
            This implement the fastNNI algorithm proposed in PHYML paper
            TUNG: this is a virtual function, so it will be called automatically by optimizeNNIBranches()
            @return best likelihood found
            @param skipped (OUT) 1 if current iteration is skipped, otherwise 0
            @param nni_count (OUT) the number of single NNI moves proceeded so far
     */
    virtual double optimizeNNI(bool beginHeu=false, int *skipped = NULL, int *nni_count = NULL);

    /**
            search all positive NNI move on the current tree and save them on the possilbleNNIMoves list            
     */
    virtual void genNNIMoves(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            search all positive NNI move on the current tree and save them on the possilbleNNIMoves list            
     */
    void genNNIMovesSort();


    /**
            search the best swap for a branch
            @return NNIMove The best Move/Swap
            @param cur_score the current score of the tree before the swaps
            @param node1 1 of the 2 nodes on the branch
            @param node2 1 of the 2 nodes on the branch
     */
    virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, double lh_contribution = -1.0);

    /**
            distance matrix, used for IQP algorithm
     */
    double *dist_matrix;

    /**
            add a NNI move to the list of possible NNI moves;
     */
    void addPossibleNNIMove(NNIMove myMove);

    /**
     * Described in PhyML paper: apply changes to all branches that do not correspond to a swap
     * with the following formula  l = l + lamda(la - l)
     * @param node
     * @param dad
     */
    void applyAllBranchLengthChanges(PhyloNode *node, PhyloNode *dad = NULL);

    /**
     * 	Save all the current branch lengths
     */
    void saveBranLens(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * 	 Restore the branch lengths from the saved values
     */
    virtual void restoreAllBranLen(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * Get the branch length of the branch node1-node2
     * @param node1
     * @param node2
     * @return the branch length
     */
    double getBranLen(PhyloNode *node1, PhyloNode *node2);


    /**
            Described in PhyML paper: apply change to branch that does not correspond to a swap with the following formula l = l + lamda(la - l)
            @param node1 the first node of the branch
            @param node2 the second node of the branch
     */
    double applyBranchLengthChange(PhyloNode *node1, PhyloNode *node2, bool nonNNIBranch);

    /**
            TODO
     */
    void applyChildBranchChanges(PhyloNode *node, PhyloNode *dad);


    /**
     * TODO
     * @return
     */
    int estimateNNICount(void);

    /**
     * TODO
     * @return
     */
    double estimateNNIDelta(void);

    /**
     *
     * @return
     */
    double getCurScore(void);

    /**
     *
     * @return
     */
    double getBestScore(void) {
        return bestScore;
    }

    /**
     *
     */
    void setBestScore(double score) {
        bestScore = score;
    }

    /**
            current parsimony score of the tree
     */
    int cur_pars_score;

    bool enable_parsimony;
    /**
            stopping rule
     */
    StopRule stop_rule;

    /**
     * Current score of the tree;
     */
    double curScore;

    /**
     *      Parsimony scores, used for linear regression
     */
    double* pars_scores;

	/**
		Log-likelihood variance
	*/
	double logl_variance;

    /**
     *      The coressponding log-likelihood score from computed indendently from the parsimony
     *      scores
     */
    double* lh_scores;

    Linear* linRegModel;

    inline void disableHeuristic() {
        enableHeuris = false;
    }

    inline void setSpeed_conf(double speed_conf) {
        this->speed_conf = speed_conf;
    }

    inline double getSpeed_conf() const {
        return speed_conf;
    }

    inline void setStartLambda(double startLambda) {
        this->startLambda = startLambda;
    }

    inline double getStartLambda() const {
        return startLambda;
    }


	inline double getNNICutoff() { return nni_cutoff; }

protected:


    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;

    /**
       The lamda number for NNI process (described in PhyML Paper)
     */
    double startLambda;

    /**
     * Array that stores the frequency that each taxa has been choosen to be swapped
     */
    map<int, int> freqList;

    /**
     * Taxa set
     */
    NodeVector taxaSet;

    int nbNNI;

    /**
     * 95% confidence value for number of NNIs found in one iteration
     */
    unsigned int nni_count_est;

    /**
     * 95% confidence value for likelihood improvement made by one NNI
     */
    double nni_delta_est;

    /**
     * Enable/Disable speed-up heuristics
     */
    bool enableHeuris;

    /**
     * Confidence level for speed up heuristics
     */
    double speed_conf;

    /**
     *  Vector contains number of NNIs used at each iterations
     */
    vector<int> vecNumNNI;    

    /**
     *  Vector contains approximated improvement pro NNI at each iterations
     */
    vector<double> vecImpProNNI;    

    /**
     * The current best score found
     */
    double bestScore;

    /**
            The list of possible NNI moves for the current tree;
     */
    vector<NNIMove> posNNIs;


    /**
            List contains non-conflicting NNI moves for the current tree;
     */
    vector<NNIMove> indeNNIs;

    /**
            Data structure (of type Map) which stores all the optimal branch lengths for all branches in the tree
     */
    BranLenMap mapOptBranLens;

    /**
     * 	Data structure (of type Map) used to store the original branch lengths of the tree
     */
    BranLenMap savedBranLens;


    int k_delete, k_delete_min, k_delete_max, k_delete_stay;

    /**
            number of representative leaves for IQP step
     */
    int k_represent;

public:
	int save_all_trees;
	/**
		this keeps the list of intermediate trees. 
		it will be activated if params.avoid_duplicated_trees is TRUE.
	*/
	StringIntMap treels;

	vector<double* > treels_ptnlh;
	
	DoubleVector treels_logl;

protected:
	/**
	*/
	vector<NNIInfo> nni_info;

	int nni_round;

	bool estimate_nni_cutoff;

	double nni_cutoff;

	bool nni_sort;

	bool testNNI;

	ofstream outNNI;
protected:
	bool print_tree_lh;
	bool write_intermediate_trees;

	ofstream out_treels, out_treelh, out_sitelh;

	void estimateNNICutoff(Params &params);

	void saveCurrentTree(PhyloNode* node1, PhyloNode *node2); // save current tree
	


    /**
            number of IQPNNI iterations
     */
    //int iqpnni_iterations;

    /**
            bonus values of all branches, used for IQP algorithm
     */
    //double *bonus_values;

    /**
            delete a set of leaves from tree (with the probability p_delete), assume tree is birfucating
            @param del_leaves (OUT) the list of deleted leaves
     */
    void deleteLeaves(PhyloNodeVector &del_leaves);

    /**
            reinsert the whole list of leaves back into the tree
            @param del_leaves the list of deleted leaves, returned by deleteLeaves() function
     */
    virtual void reinsertLeaves(PhyloNodeVector &del_leaves);

    /**
            assess a quartet with four taxa. Current implementation uses the four-point condition
            based on distance matrix for quick evaluation.
            @param leaf0 one of the leaf in the existing sub-tree
            @param leaf1 one of the leaf in the existing sub-tree
            @param leaf2 one of the leaf in the existing sub-tree
            @param del_leaf a leaf that was deleted (not in the existing sub-tree)
            @return 0, 1, or 2 depending on del_leaf should be in subtree containing leaf0, leaf1, or leaf2, respectively
     */
    int assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf);

    /**
            assess a quartet with four taxa using parsimony
            @param leaf0 one of the leaf in the existing sub-tree
            @param leaf1 one of the leaf in the existing sub-tree
            @param leaf2 one of the leaf in the existing sub-tree
            @param del_leaf a leaf that was deleted (not in the existing sub-tree)
            @return 0, 1, or 2 depending on del_leaf should be in subtree containing leaf0, leaf1, or leaf2, respectively
     */
    int assessQuartetParsimony(Node *leaf0, Node *leaf1, Node *leaf2,
            Node *del_leaf);

    /**
            assess the important quartets around a virtual root of the tree.
            This function will assign bonus points to branches by updating the variable 'bonus_values'
            @param cur_root the current virtual root
            @param del_leaf a leaf that was deleted (not in the existing sub-tree)
     */
    void assessQuartets(vector<RepresentLeafSet*> &leaves_vec, PhyloNode *cur_root, PhyloNode *del_leaf);

    /**
            initialize the bonus points to ZERO
            @param node the root of the sub-tree
            @param dad dad of 'node', used to direct the recursion
     */
    void initializeBonus(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            raise the bonus points for all branches in the subtree rooted at a node
            @param node the root of the sub-tree
            @param dad dad of 'node', used to direct the recursion
     */
    void raiseBonus(Neighbor *nei, Node *dad, double bonus);

    /**
            Bonuses are stored in a partial fashion. This function will propagate the bonus at every branch
            into the subtree at this branch.
            @param node the root of the sub-tree
            @param dad dad of 'node', used to direct the recursion
            @return the partial bonus of the branch (node -> dad)
     */
    double computePartialBonus(Node *node, Node* dad);

    /**
            determine the list of branches with the same best bonus point
            @param best_bonus the best bonus determined by findBestBonus()
            @param best_nodes (OUT) vector of one ends of the branches with highest bonus point
            @param best_dads (OUT) vector of the other ends of the branches with highest bonus point
            @param node the root of the sub-tree
            @param dad dad of 'node', used to direct the recursion
     */
    void findBestBonus(double &best_score, NodeVector &best_nodes, NodeVector &best_dads, Node *node = NULL, Node *dad = NULL);

};

void estimateNNICutoff(Params &params);


#endif
