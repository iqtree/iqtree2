/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include "utils/stoprule.h"
#include "mtreeset.h"
#include "node.h"
#include "candidateset.h"
#include "utils/pllnni.h"

typedef std::map< string, double > mapString2Double;
typedef std::multiset< double, std::less< double > > multiSetDB;
typedef std::multiset< int, std::less< int > > MultiSetInt;

class RepLeaf {
public:
    Node *leaf;
    int height;

    RepLeaf(Node *aleaf, int aheight = 0) {
        leaf = aleaf;
        height = aheight;
    }
};

/**
        nodeheightcmp, for building k-representative leaf set
 */
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
    Main class for tree search
 */
class IQTree : public PhyloTree {
public:
    /**
            default constructor
     */
    IQTree();

    IQTree(Alignment *aln);

//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
            destructor
     */
    virtual ~IQTree();

    void init();

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    /** 
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    /**
        save UFBoot_trees.
        For MPI workers only save from sample_start to sample_end
        @param checkpoint Checkpoint object
    */
    void saveUFBoot(Checkpoint *checkpoint);

    /**

        restore UFBoot_trees from sample_start to sample_end (MPI)
        @param checkpoint Checkpoint object
    */
    void restoreUFBoot(Checkpoint *checkpoint);

    /**
     * setup all necessary parameters  (declared as virtual needed for phylosupertree)
     */
    virtual void initSettings(Params& params);

    void createPLLPartition(Params &params, ostream &pllPartitionFileHandle);

    void initializePLL(Params &params);

    bool isInitializedPLL();
    
    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    /**
            print tree to .treefile
            @param params program parameters, field root is taken
     */
    virtual void printResultTree(string suffix = "");
    /**
            print tree to out
            @param params program parameters, field root is taken
            @param out (OUT) output stream
     */
    void printResultTree(ostream &out);

    void printBestCandidateTree();

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
//    void setIQPIterations(STOP_CONDITION stop_condition, double stop_confidence, int min_iterations, int max_iterations);

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
     *  @brief get non-tabu branches from a set of branches
     *
     *  @param
     *      allBranches[IN] the inital branches
     *  @param
     *      initTabuSplits[IN] the tabu splits
     *  @param
     *        nonTabuBranches[OUT] non-tabu branches from \a allBranches
     *    @param[OUT]
     *        tabuBranches branches that are tabu
     */
    void getNonTabuBranches(Branches& allBranches, SplitGraph& tabuSplits, Branches& nonTabuBranches, Branches* tabuBranches = NULL);

    /**
     * @brief remove all branches corresponding to nnis
     * @param nodes1 node vector containing one end of the branches
     * @param nodes2 node vector containing the other end of the branches
     * @param nnis
     * @return
     */
    int removeNNIBranches(NodeVector& nodes1, NodeVector& nodes2, unordered_map<string, NNIMove> nnis);

    /**
     *         Perform a series of random NNI moves
     *         @return the perturbed newick string
     */
    string doRandomNNIs(bool storeTabu = false);

    /**
     *  Do a random NNI on splits that are shared among all the candidate trees.
     *  @return the perturbed newick string
     */
    string perturbStableSplits(double supportValue);


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
    virtual double doTreeSearch();

    /**
     *  Wrapper function that uses either PLL or IQ-TREE to optimize the branch length
     *  @param maxTraversal
     *      maximum number of tree traversal for branch length optimization
     *  @return NEWICK tree string
     */
    string optimizeBranches(int maxTraversal = 100);

    /**
     *  Wrapper function to compute tree log-likelihood.
     *  This function with call either PLL or IQ-TREE to compute tree log-likelihood
     *  @return current score of tree
     */
    double computeLogL();

    /**
     *    Print scores of tree used for generating offsprings
     */
    void printBestScores();

    /****************************************************************************
            Fast Nearest Neighbor Interchange by maximum likelihood
     ****************************************************************************/


    /**
     *  Optimize current tree using NNI
     *
     *  @return
     *      <number of NNI steps, number of NNIs> done
     */
    virtual pair<int, int> optimizeNNI(bool speedNNI = true);

    /**
     *  Return the current best score found
     */
    double getBestScore();

    /**
     * @brief Generate a list of internal branches on which NNI moves will be evaluated
     * @param
     *      nonNNIBranches [OUT] Branches on which NNI evaluation will be skipped
     * @param
     *      tabuSplits [IN] A list of splits that are considered tabu
     * @param
     *      candidateSplitHash [IN] Lists that appear on the best 20 candidate trees
     * @param
     *      dad [IN] for navigation
     * @param
     *      node[IN] for navigation
     * @return A list of branches for evaluating NNIs
     */
    void getNNIBranches(SplitIntMap &tabuSplits, SplitIntMap &candidateSplitHash, Branches &nonNNIBranches, Branches &outBranches, Node *dad = NULL, Node *node = NULL);

    /**
     *  Return internal branches that appear in \a candidateSplitHash
     *  and has support value >= \a supportValue.
     *  @param
     *      candidateSplitHash [IN]   A set of splits with the number of occurences.
     *  @param
     *      supportValue [IN]  Only consider split whose support value is higher than this number
     *  @param
     *      dad [IN] for navigation
     *  @param
     *      node[IN] for navigation
     *  @return
     *      A list of branches fufilling the aforementioned conditions.
     */
    void getStableBranches(SplitIntMap &candSplits, double supportValue, Branches &outBranches, Node *dad = NULL, Node *node = NULL);


    /**
     *
     *  Determine whether to evaluate NNI moves on the branch corresponding to the current split
     *
     *  @param curSplit [IN] the split that correspond to the current branch
     *  @param tabuSplits [IN] tabu splits
     *  @param candSplits [IN] splits contained in all candidate trees
     *  @param nonNNIBranches [OUT] branches that are not inserted to nniBranches are store here
     *  @param nniBranches [OUT] if the split is neither stable nor tabu it is inserted in this list
     */
    bool shouldEvaluate(Split* curSplit, SplitIntMap &tabuSplits, SplitIntMap &candSplits);


    /**
     *  @brief Only select NNI branches that are 2 branches away from the previously
     *  appied NNIs
     *  @param
     *      appliedNNIs List of previously applied NNIs
     *  @return
     *      List of branches to be evaluated
     */
    void filterNNIBranches(vector<NNIMove> &appliedNNIs, Branches &outBranches);

    
    /**
     * @brief get branches that correspond to the splits in \a nniSplits
     */
    void getSplitBranches(Branches &branches, SplitIntMap &splits, Node *dad = NULL, Node *node = NULL);

    /**
     *         Do fastNNI using PLL
     *
     *      @param nniCount (OUT) number of NNIs applied
     *         @param nniSteps (OUT) number of NNI steps done
     */
    double pllOptimizeNNI(int &nniCount, int &nniSteps, SearchInfo &searchinfo);

    /**
     *         @brief Perform NNI search on the current tree topology
     *         @return <number_of_NNIs, number_of_NNI_steps>
     *         This function will automatically use the selected kernel (either PLL or IQ-TREE)
     */
    virtual pair<int, int> doNNISearch(bool write_info = false);

    /**
            @brief evaluate all NNIs
            @param  node    evaluate all NNIs of the subtree rooted at node
            @param  dad     a neighbor of \p node which does not belong to the subtree
                            being considered (used for traverse direction)

     */
    //void evalNNIs(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * @brief Evaluate all NNIs on branch defined by \a branches
     *
     * @param nniBranches [IN] branches the branches on which NNIs will be evaluated
     * @return list positive NNIs
     */
    void evaluateNNIs(Branches &nniBranches, vector<NNIMove> &outNNIMoves);

    double optimizeNNIBranches(Branches &nniBranches);

    /**
            search all positive NNI move on the current tree and save them
            on the possilbleNNIMoves list
     */
    void evalNNIsSort(bool approx_nni);

    /**
            apply  NNIs from the non-conflicting NNI list
            @param compatibleNNIs vector of all compatible NNIs
            @param changeBran whether or not the computed branch lengths should be applied
     */
    virtual void doNNIs(vector<NNIMove> &compatibleNNIs, bool changeBran = true);

    /**
     *  Restore the old 5 branch lengths stored in the NNI move.
     *  This is called after an NNI is reverted.
     *  @param nnimove the NNI move currently in consideration
     */
    //void restoreNNIBranches(NNIMove nnimove);


    /**
     *  @brief get a list of compatible NNIs from a list of NNIs
     *  @param nniMoves [IN] list of NNIs
     *  @return list of compatible NNIs
     */
    void getCompatibleNNIs(vector<NNIMove> &nniMoves, vector<NNIMove> &compatibleNNIs);

    /**
            add a NNI move to the list of possible NNI moves;
     */
    void addPositiveNNIMove(NNIMove &myMove);

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
            current parsimony score of the tree
     */
    int cur_pars_score;

//    bool enable_parsimony;
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
     *  information and parameters for the tree search procedure
     */
    SearchInfo searchinfo;

    /**
     *  Vector contains number of NNIs used at each iterations
     */
    vector<int> vecNumNNI;


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
    *  Splits corresponding to random NNIs
    */
    SplitIntMap initTabuSplits;

    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;


    /**
     * Taxa set
     */
    NodeVector taxaSet;

    /**
     *  Vector contains approximated improvement pro NNI at each iterations
     */
    vector<double> vecImpProNNI;

    /**
        Optimal branch lengths
     */
//    mapString2Double optBrans;

    /**
     *  @brief get branches, on which NNIs are evaluated for the next NNI step.
     *  @param[out] nodes1 one ends of the branches
     *  @param[out] nodes2 the other ends of the branches
     *  @param[in] nnis NNIs that have been previously applied
     */
    void generateNNIBranches(NodeVector& nodes1, NodeVector& nodes2, unordered_map<string, NNIMove>& nnis);

    int k_delete, k_delete_min, k_delete_max, k_delete_stay;

    /**
            number of representative leaves for IQP step
     */
    int k_represent;

public:

    /**
     *  Candidate tree set (the current best N (default N = 5)
     *  NNI-optimal trees
     */
    CandidateSet candidateTrees;

    /**
     *  Set of all intermediate trees (initial trees, tree generated by NNI steps,
     *  NNI-optimal trees)
     */
    CandidateSet intermediateTrees;


    /**
     *  Update the candidate set with a new NNI-optimal tree. The maximum size of the candidate set
     *  is fixed to the initial setting. Thus, if the size exceed the maximum number of trees, the worse
     *  tree will be removed.
     *
     *  @param treeString
     *      the new tree
     *  @param score
     *      the score of the new tree
     *  @param updateStopRule
     *      Whether or not to update the stop rule
     *  @return relative position of the new tree to the current best.
     *      -1 if duplicated
     *      -2 if the candidate set is not updated
     */
    int addTreeToCandidateSet(string treeString, double score, bool updateStopRule, int sourceProcID);

    /**
        MPI: synchronize candidate trees between all processes
        @param nTrees number of trees to broadcast
        @param updateStopRule true to update stopping rule, false otherwise
    */
    void syncCandidateTrees(int nTrees, bool updateStopRule);

    /**
        MPI: synchronize tree of current iteration with master
        will update candidateset_changed
        @param curTree current tree

    */
    void syncCurrentTree();

    /**
        MPI: Master sends stop message to all workers
    */
    void sendStopMessage();

    /**
     *  Generate the initial parsimony/random trees, called by initCandidateTreeSet
     *  @param nParTrees number of parsimony/random trees to generate
     */
    void createInitTrees(int nParTrees);

    /**
     *  Generate the initial candidate tree set
     *  @param nParTrees number of parsimony trees to generate
     *  @param nNNITrees number of NNI locally optimal trees to generate
     */
    void initCandidateTreeSet(int nParTrees, int nNNITrees);

    /**
     * Generate the initial tree (usually used for model parameter estimation)
     */
    void computeInitialTree(LikelihoodKernel kernel, istream* in = NULL);

    /**
     *  @brief: optimize model parameters on the current tree
     *  either IQ-TREE or PLL
     *  @param printInfo to print model parameters to the screen or not
     *  @param epsilon likelihood epsilon for optimization
     *
     */
    virtual string optimizeModelParameters(bool printInfo = false, double epsilon = -1);

    /**
     *  @brief: either optimize model parameters on the current tree
     *  or restore them from a checkpoint (this function exists because the
     *  same things need to be done in two different places, in runTreeReconstruction)
     *  @param initEpsilon likelihood epsilon for optimization
     */

    virtual string ensureModelParametersAreSet(double initEpsilon);
    
    /**
     *  variable storing the current best tree topology
     */
    topol* pllBestTree;

    /****** following variables are for ultra-fast bootstrap *******/

    /** TRUE to save also branch lengths into treels_newick */
//    bool save_all_br_lens;

    /**
        this keeps the list of intermediate trees.
        it will be activated if params.avoid_duplicated_trees is TRUE.
     */
//    StringIntMap treels;

    /** pattern log-likelihood vector for each treels */
//    vector<double* > treels_ptnlh;

    /** OBSOLETE: tree log-likelihood for each treels */
//    DoubleVector treels_logl;

    /** NEWICK string for each treels */
//    StrVector treels_newick;

    /** OBSOLETE: maximum number of distinct candidate trees (tau parameter) */
//    int max_candidate_trees;

    /** log-likelihood threshold (l_min) */
    double logl_cutoff;

    /** vector of bootstrap alignments generated */
    vector<BootValType* > boot_samples;

    /** starting sample for UFBoot, used for MPI */
    int sample_start;

    /** end sample for UFBoot, used for MPI */
    int sample_end;

    /** newick string of corresponding bootstrap trees */
    StrVector boot_trees;

    /** bootstrap tree strings with branch lengths, for -wbtl option */
//    StrVector boot_trees_brlen;

    /** number of multiple optimal trees per replicate */
    IntVector boot_counts;

    /** corresponding RELL log-likelihood */
    DoubleVector boot_logl;

    /** corresponding log-likelihood on original alignment */
    DoubleVector boot_orig_logl;

    /** Set of splits occurring in bootstrap trees */
    vector<SplitGraph*> boot_splits;

    /** log-likelihood of bootstrap consensus tree */
    double boot_consense_logl;

    /** Robinson-Foulds distance between contree and ML tree */
    int contree_rfdist;

    /** Corresponding map for set of splits occurring in bootstrap trees */
    //SplitIntMap boot_splits_map;

    /** summarize all bootstrap trees */
    void summarizeBootstrap(Params &params, MTreeSet &trees);

    /** summarize bootstrap trees */
    virtual void summarizeBootstrap(Params &params);

    /** summarize bootstrap trees into split set */
    void summarizeBootstrap(SplitGraph &sg);

    /**
     write .ufboot trees file
     */
    virtual void writeUFBootTrees(Params &params);

    /** @return bootstrap correlation coefficient for assessing convergence */
    double computeBootstrapCorrelation();

    /**
        compute rootstrap supports for rooted tree, 2021-01-19 for Suha's paper
        @param[in] trees set of rooted trees
        @param[in] use_taxid true if rooted_trees has taxa IDs, false for taxa names
     */
    void computeRootstrap(MTreeSet &trees, bool use_taxid);

    /**
        compute rootstrap supports for unrooted tree using outgroup, 2021-01-19 for Suha's paper
     */
    void computeRootstrapUnrooted(MTreeSet &trees, const char* outgroup, bool use_taxid);

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

    //bool print_tree_lh;

    //int write_intermediate_trees;

    ofstream out_treels, out_treelh, out_sitelh, out_treebetter;
    string treels_name, out_lh_file, site_lh_file;

    void estimateNNICutoff(Params* params);

    virtual void saveCurrentTree(double logl); // save current tree


    void saveNNITrees(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    int duplication_counter;

    // MPI: vector of size = num processes, true if master should send candidate set to worker
    BoolVector candidateset_changed;

    // true if best candidate tree is changed
    bool bestcandidate_changed;

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
     *         delete a set of leaves from tree
     *         non-cherry leaves are selected first
     *         @param del_leaves (OUT) the list of deleted leaves
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

    void convertNNI2Splits(SplitIntMap &nniSplits, int numNNIs, vector<NNIMove> &compatibleNNIs);

    string generateParsimonyTree(int randomSeed);

    double doTreePerturbation();

    void estimateLoglCutoffBS();

    //void estimateNNICutoff(Params &params);

public:
    /**
     *  Return best tree string from the candidate set
     *
     *  @param numTrees
     *      Number of best trees to return
     *  @return
     *      A string vector of trees
     */
    virtual vector<string> getBestTrees(int numTrees = 0);

    /**
     *  Print the iteration number and the tree score
     */
    void printIterationInfo(int sourceProcID);

    /**
     *  Return branches that are 2 branches away from the branches, on which NNIs were applied
     *  in the previous NNI steps.
     *  @param
     *      previousNNIBranches[IN] a set of branches on which NNIs were performed in the previous NNI step.
     *  @return
     *      a set of branches, on which NNIs should be evaluated for the current NNI steps
     */
    Branches getReducedListOfNNIBranches(Branches &previousNNIBranches);


    // Diep added for UFBoot2-Corr
    void refineBootTrees();
    bool on_refine_btree;
    Alignment* saved_aln_on_refine_btree;
    vector<IntVector> boot_samples_int;
};
#endif
