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
#include "node.h"

#include "pll/pll.h"
#include "nnisearch.h"
#include "candidateset.h"

#define BOOT_VAL_FLOAT
#define BootValType float
//#define BootValType double

typedef std::map< string, double > mapString2Double;
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
    //void printPhylolibModelParams(const char* suffix);

    /**
        print intermediate tree
     */
    void printIntermediateTree(int brtype);

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
            remove a portion of leaves and reinsert them using the IQP algorithm
     */
    void doIQP();

    /**
     * 		Perform a series of random NNI moves
     * 		@param numNNI number of random NNIs
     */
    void doRandomNNIs(int numNNI);

    /**
     *   input model parameters from IQ-TREE to PLL
     */
    void inputModelIQTree2PLL();

    /**
     *  input model parameters from PLL to IQ-TREE
     */
    void inputModelPLL2IQTree();

    /**
     *  get the rate parameters from PLL
     *  @return double array containing the 6 rates
     */
    double* getModelRatesFromPLL();

    /**
     *  get the alpha parameter from PLL for the GAMMA distribution of rate heterogenity
     *  @return alpha parameter
     */
    double getAlphaFromPLL();

    /**
     *  print model parameters from PLL
     */
    void pllPrintModelParams();

    /**
     * input the tree string from IQTree kernel to PLL kernel
     * @return
     */
    double inputTree2PLL(string treestring, bool computeLH = true);

    //bool containPosNNI(vector<NNIMove> posNNIs);

    /**
     * Perturb the tree for the next round of local search by swaping position of 2 random leaves
     * @param nbDist The minimum distance between the 2 nodes that are swapped
     * @param nbTimes Number of times that the swap operations are carried out
     * @return The new loglikelihood of the tree
     */
    double perturb(int times);

    /**
     * TODO
     * @param node1
     * @param node2
     * @return
     */
    double swapTaxa(PhyloNode *node1, PhyloNode *node2);

    /**
            perform tree search
            @return best likelihood found
     */
    double doTreeSearch();

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
    double optimizeNNI(int &nni_count, int &nni_steps);

    /**
     * 		Do fastNNI using PLL
     *
     *      @param nniCount (OUT) number of NNIs applied
     * 		@param nniSteps (OUT) number of NNI steps done
     */
    double pllOptimizeNNI(int &nniCount, int &nniSteps, SearchInfo &searchinfo);

    /**
            @brief evaluate all NNIs and store them in possilbleNNIMoves list
            @param  node    evaluate all NNIs of the subtree rooted at node
            @param  dad     a neighbor of \p node which does not belong to the subtree
                            being considered (used for traverse direction)

     */
    void evalNNIs(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * @brief Evaluate all NNIs defined in \a brans
     * @param[in] brans contains branches whose NNIs need to be evaluated
     */
    void evalNNIs(map<string, Branch> brans);

    /**
            search all positive NNI move on the current tree and save them
            on the possilbleNNIMoves list
     */
    void evalNNIsSort(bool approx_nni);

    /**
            apply nni2apply NNIs from the non-conflicting NNI list
            @param nni2apply number of NNIs to apply from the list
            @param changeBran whether or not the computed branch lengths should be applied
     */
    virtual void doNNIs(int nni2apply, bool changeBran = true);

    /**
     *  Restore the old 5 branch lengths stored in the NNI move.
     *  This is called after an NNI is reverted.
     *  @param nnimove the NNI move currently in consideration
     */
    //void restoreNNIBranches(NNIMove nnimove);

    /**
            generate non conflicting NNI moves.
            moves are saved in vec_nonconf_nni
     */
    void genNonconfNNIs();

    /**
            add a NNI move to the list of possible NNI moves;
     */
    void addPositiveNNIMove(NNIMove myMove);

    /**
     * 	Save all the current branch lengths
     */
    void saveBranches(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * 	 Restore the branch lengths from the saved values
     */
    virtual void restoreAllBrans(PhyloNode *node = NULL, PhyloNode *dad = NULL);

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
    double getAvgNumNNI(void);

    /**
     * Estimate the median of the distribution of N (see paper for more d
                                                          details)
     * @return the estimated value
     */
    double estDeltaMedian(void);

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
     *  set the current tree as the best tree
     */
    void setBestTree(string tree, double logl);

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


    inline double getNNICutoff() {
        return nni_cutoff;
    }

    /*
     *  Contains a sorted list of all NNIs (2n-6) evaluated for the current best tree
     *  The last element (nni_for_pertub.end()) is the best NNI
     */
    vector<pllNNIMove> nniListOfBestTree;


    /**
     *  Instance of the phylogenetic likelihood library. This is basically the tree data strucutre in RAxML
     */
    pllInstance *pllInst;

    /**
     *	PLL data structure for alignment
     */
    pllAlignmentData *pllAlignment;

    /**
     *  PLL data structure for storing phylognetic analysis options
     */
    pllInstanceAttr pllAttr;

    /**
     *  PLL partition list
     */
    partitionList * pllPartitions;

    /**
     *  information and parameters for the tree search procedure
     */
    SearchInfo searchinfo;

    /**
     *  Vector contains number of NNIs used at each iterations
     */
    vector<int> vecNumNNI;

    int getCurIteration() { return curIt; }

    /**
     * Do memory allocation and initialize parameter for UFBoot to run with PLL
     */
    void pllInitUFBootData();

    /**
     * Do memory deallocation for UFBoot data (PLL mode)
     */
    void pllDestroyUFBootData();

    /**
     * DTH:
     * Substitute bases in seq according to PLL's rules
     * This function should be updated if PLL's rules change.
     * @param seq: data of some sequence to be substituted
     * @param dataType: PLL_DNA_DATA or PLL_AA_DATA
     */
   void pllBaseSubstitute (char *str, int dataType);

   /*
    * An array to map site index in pllAlignment into IQTree pattern index
    * Born due to the order difference of these two
    * Will be deallocated in pllDestroyUFBootData()
    */
   int * pll2iqtree_pattern_index;

   /*
    * Build pll2iqtree_pattern_index
    * Must be called AFTER initializing PLL model
    */
   void pllBuildIQTreePatternIndex();

   /**
    * FOR TESTING:
    * Write to log file the freq of pllAlignment sites, and
    * freq of bootstrap site stored in pllUFBootDataPtr->boot_samples
    */
   void pllLogBootSamples(int** pll_boot_samples, int nsamples, int npatterns);

   /**
    * Convert certain arrays in pllUFBootDataPtr
    * into IQTree data structures
    * to be usable in IQTree::summarizeBootstrap()
    */
   void pllConvertUFBootData2IQTree();

protected:

    /**
     *  Current IQPNNI iteration number
     */
    int curIt;
    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;


    /**
     * Taxa set
     */
    NodeVector taxaSet;

    /**
     * confidence value for number of NNIs found in one iteration
     */
    int nni_count_est;

    /**
     * confidence value for likelihood improvement made by one NNI
     */
    double nni_delta_est;


    /**
     *  Vector contains approximated improvement pro NNI at each iterations
     */
    vector<double> vecImpProNNI;

    /**
        List of positive NNI for the current tree;
     */
    vector<NNIMove> plusNNIs;

    /**
        List of non-conflicting NNIs for the current tree;
     */
    vector<NNIMove> nonConfNNIs;

    /**
     *  NNIs that have been applied in the previous step
     */
    vector<NNIMove> appliedNNIs;

    /**
        Optimal branch lengths
     */
    mapString2Double optBrans;

    /**
     *  Set of all internal branches whose induced NNIs need to be evaluate
     *  in the next NNI step
     */
    map<string, Branch> brans2Eval;

    /**
     *  Update \a brans2Eval after \a all NNIs in nnis have been performed
     */
    void updateBrans2Eval(vector<NNIMove> nnis);

    /**
     *  Use fastNNI heuristic
     */
    bool fastNNI;

    /**
            Original branch lengths
     */
    mapString2Double orgBrans;

    int k_delete, k_delete_min, k_delete_max, k_delete_stay;

    /**
            number of representative leaves for IQP step
     */
    int k_represent;

public:

    /**
     *  @brief: optimize model parameters on the current tree
     *  either IQ-TREE or PLL
     *  @param imd_tree the input tree or NULL
     *  @param printInfo to print model parameters to the screen or not
     */
    string optimizeModelParameters(bool printInfo=false);

    /**
     *  variable storing the current best tree topology
     */
    topol* pllBestTree;


    /**
     * The current best score found
     */
    double bestScore;

    /**
     *  the current best tree
     */
    string bestTreeString;

    CandidateSet candidateTrees;


    /****** following variables are for ultra-fast bootstrap *******/

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
    vector<BootValType* > boot_samples;

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
