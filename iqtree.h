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

#include <set>
#include <map>
#include <stack>
#include <vector>
#include "phylotree.h"
#include "phylonode.h"
#include "stoprule.h"
#include "mtreeset.h"

#include "phylolib/axml.h"
#include "nnisearch.h"
#include "phylolib.h"

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

struct IntBranchInfo {
    PhyloNode *node1;
    PhyloNode *node2;
    double lh_contribution; // log-likelihood contribution of this branch: L(T)-L(T|e=0)
};

inline int int_branch_cmp(const IntBranchInfo a, const IntBranchInfo b) {
    return (a.lh_contribution < b.lh_contribution);
}


/**
        Representative Leaf Set, stored as a multiset template of STL,
        sorted in ascending order of leaf's height
 */
typedef multiset<RepLeaf*, nodeheightcmp> RepresentLeafSet;

/**
Important Quartet Puzzling

        @author BUI Quang Minh <minh.bui@univie.ac.at>
 */
class IQTree : public PhyloTree {
public:
    /**
            default constructor
     */
    IQTree();

    IQTree(Alignment *aln);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
            destructor
     */
    virtual ~IQTree();

    void init();

    /**
     *   Create a parsimony tree using phylolib
     */
    double computeParsimonyTreePhylolib();

    /**
     * setup all necessary parameters  (declared as virtual needed for phylosupertree)
     */
    virtual void setParams(Params& params);

    /**
            print tree to .treefile
            @param params program parameters, field root is taken
     */
    void printResultTree(string suffix = "");
    /**
            print tree to out
            @param params program parameters, field root is taken
            @param out (OUT) output stream
     */
    void printResultTree(ostream &out);

    /**
     * print phylolib tree to a file.
     * @param suffix suffix string for the tree file
     */
    void printPhylolibTree(const char* suffix);


    /**
     *  print model parameters of Phylolib(rates, base frequencies, alpha) to stdout and
     *  to file
     */
    void printPhylolibModelParams(const char* suffix);

    /**
        print intermediate tree
     */
    void printIntermediateTree(int brtype);

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
            @param paramters given through command line and others
            @return current likelihood
     */
    double doIQP();

    /**
     * 	   perform a variable neighborhood search using
     * 	   NNI and SPR as the 2 alternative neighborhood
     *
     */
    double doVNS();

    bool containPosNNI(vector<NNIMove> posNNIs);

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
     */
    double doIQPNNI();

    /**
     * 		Perform random restart heuristic
     */
    void doRandomRestart();

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
    double optimizeNNI(bool beginHeu = false, int *skipped = NULL, int *nni_count = NULL);

    /**
     * 		Do fastNNI using RAxML kernel
     * 		@param beginHeu whether the heuristic is started
     * 		@param skipped (OUT) 1 if current iteration is skipped, otherwise 0
     *      @param nni_count (OUT) the number of single NNI moves proceeded so far
     */
    double optimizeNNIRax(bool beginHeu = false, int *skipped = NULL, int *nni_count = NULL);


    /**
            search all positive NNI move on the current tree and save them on the possilbleNNIMoves list
     */
    void genNNIMoves(bool approx_nni, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            search all positive NNI move on the current tree and save them
            on the possilbleNNIMoves list
     */
    void genNNIMovesSort(bool approx_nni);

    /**
            apply nni2apply NNIs from the non-conflicting NNI list
            @param nni2apply number of NNIs to apply from the list
     */
    void applyNNIs(int nni2apply);

    /**
            generate non conflicting NNI moves.
            moves are saved in vec_nonconf_nni
     */
    void genNonconfNNIs();

    /**
       search for the best NNI move corresponding to this branch
       @return NNIMove the best NNI, this NNI could be worse than the current tree
       according to the evaluation scheme in use
       @param node1 1 of the 2 nodes on the branch
       @param node2 1 of the 2 nodes on the branch
     * @param approx_nni evaluate NNI based on "Bayes"
     * @param useLS evaluate NNI based on Least Square
     */
    virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, bool approx_nni = false, bool useLS = false, double lh_contribution = -1.0);

    /**
            add a NNI move to the list of possible NNI moves;
     */
    void addPositiveNNIMove(NNIMove myMove);

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
            Described in PhyML paper: apply change to branch that does not
            correspond to a swap with the following formula l = l + lamda(la - l)
            @param node1 the first node of the branch
            @param node2 the second node of the branch
     */
    void changeBranLen(PhyloNode *node1, PhyloNode *node2, double branLen);

    /**
            Change all branch length according to the computed values during
     * 		NNI evaluation. There might be branches that are not be affected
     * 		since tree topology is changed after doing NNI
     */
    void changeAllBranches(PhyloNode *node = NULL, PhyloNode *dad = NULL);


    /**
     * Estimate the 95% quantile of the distribution of N (see paper for more d
                                                           details)
     * @return the estimated value
     */
    inline double estN95(void);

    /**
     * Estimate the median of the distribution of N (see paper for more d
                                                           details)
     * @return the estimated value
     */
    inline double estNMedian(void);

    /**
     * Estimate the median of the distribution of N (see paper for more d
                                                          details)
     * @return the estimated value
     */
    inline double estDeltaMedian(void);

    /**
     * Estimate the 95% quantile of the distribution of DELTA (see paper for
                                                               more detail)
     * @return the estimated value
     */
    inline double estDelta95(void);

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
     *      Parsimony scores, used for linear regression
     */
    double* pars_scores;

    /**
        Log-likelihood variastring IQPTree::bran2string(PhyloNode* node1, PhyloNode* node2)nce
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

    inline double getNNICutoff() {
        return nni_cutoff;
    }

    /**
     *  Tree data structure for RAxML kernel
     */
    tree* phyloTree;

protected:

    /**
     *  Current IQPNNI iteration number
     */
    int curIQPIter;
    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;

    /**
       The lambda number for NNI search (described in PhyML Paper)
     */
    double startLambda;

    /**
     * current lambda value in use
     */
    double curLambda;

    /**
     * Array that stores the frequency that each taxa has been choosen to be swapped
     */
    map<int, int> freqList;

    /**
     * Taxa set
     */
    NodeVector taxaSet;

    //int nbNNI;

    /**
     * confidence value for number of NNIs found in one iteration
     */
    int nni_count_est;

    /**
     * confidence value for likelihood improvement made by one NNI
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
            The list of positive NNI moves for the current tree;
     */
    vector<NNIMove> posNNIs;

    /**
     *  data structure to store delta LH (in NNICUT heuristic)
     */
    NNICUT nnicut;


    /**
            List contains non-conflicting NNI moves for the current tree;
     */
    vector<NNIMove> vec_nonconf_nni;

    /**
     *      Data structure to store how many times a leaf has been removed.
     *      LeafFreq is a struct that contains leaf_id and leaf_frequency
     */
    vector<LeafFreq> leaf_freqs;

    /**
            Data structure (of type Map) which stores all the optimal
            branch lengths for all branches in the tree
     */
    BranLenMap mapOptBranLens;

    /**
     * 	Data structure (of type Map) used to store the original branch
        lengths of the tree
     */
    BranLenMap savedBranLens;


    int k_delete, k_delete_min, k_delete_max, k_delete_stay;

    /**
            number of representative leaves for IQP step
     */
    int k_represent;

    /**
     *  Initialize the node frequency list (node_freqs)
     */
    void initLeafFrequency(PhyloNode* node = NULL, PhyloNode* dad = NULL);

    void clearLeafFrequency();

public:
    /****** following variables are for ultra-fast bootstrap *******/
    /** 2 to save all trees, 1 to save intermediate trees */
    int save_all_trees;

    /** TRUE to save also branch lengths into treels_newick */
    bool save_all_br_lens;

    /**
        this keeps the list of intermediate trees.
        it will be activated if params.avoid_duplicated_trees is TRUE.
     */
    StringIntMap treels;

    /** pattern log-likelihood vector for each treels */
    vector<double* > treels_ptnlh;

    /** tree log-likelihood for each treels */
    DoubleVector treels_logl;

    /** NEWICK string for each treels */
    StrVector treels_newick;

    /** maximum number of distinct candidate trees (tau parameter) */
    int max_candidate_trees;

    /** log-likelihood threshold (l_min) */
    double logl_cutoff;

    /** vector of bootstrap alignments generated */
    vector<IntVector> boot_samples;

    /** newick string of corresponding bootstrap trees */
    IntVector boot_trees;

	/** number of multiple optimal trees per replicate */
	IntVector boot_counts;

    /** corresponding RELL log-likelihood */
    DoubleVector boot_logl;

    /** Set of splits occuring in bootstrap trees */
    vector<SplitGraph*> boot_splits;

    /** Corresponding map for set of splits occuring in bootstrap trees */
    //SplitIntMap boot_splits_map;

    /** summarize all bootstrap trees */
    void summarizeBootstrap(Params &params, MTreeSet &trees);

    void summarizeBootstrap(Params &params);

    /** summarize bootstrap trees into split set */
    void summarizeBootstrap(SplitGraph &sg);

    /** @return TRUE if stopping criterion is met */
    bool checkBootstrapStopping();
	int getDelete() const;
	void setDelete(int _delete);

    int nni_round;
protected:
    /**** NNI cutoff heuristic *****/
    /**
     */
    vector<NNIInfo> nni_info;


    bool estimate_nni_cutoff;

    double nni_cutoff;

    bool nni_sort;

    bool testNNI;

    ofstream outNNI;
protected:


    bool print_tree_lh;

    int write_intermediate_trees;

    ofstream out_treels, out_treelh, out_sitelh, out_treebetter;

    void estimateNNICutoff(Params* params);

    virtual void saveCurrentTree(double logl); // save current tree

    void saveNNITrees(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    int duplication_counter;

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

    void deleteNonTabuLeaves(PhyloNodeVector &del_leaves);

    /**
     * 		delete a set of leaves from tree
     * 		non-cherry leaves are selected first
     * 		@param del_leaves (OUT) the list of deleted leaves
     */
    void deleteNonCherryLeaves(PhyloNodeVector &del_leaves);

    /**
            reinsert the whole list of leaves back into the tree
            @param del_leaves the list of deleted leaves, returned by deleteLeaves() function
     */
    virtual void reinsertLeaves(PhyloNodeVector &del_leaves);

    void reinsertLeavesByParsimony(PhyloNodeVector &del_leaves);

    void doParsimonyReinsertion();


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

    void estDeltaMin();

};

void estimateNNICutoff(Params &params);


#endif
