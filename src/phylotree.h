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
#define NDEBUG
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD 80
#define EIGEN_UNROLLING_LIMIT 1000
//#define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (512*256)
//#define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (8*512*512)
#include "Eigen/Core"
#include "mtree.h"
#include "alignment.h"
#include "substmodel.h"
#include "modelfactory.h"
#include "phylonode.h"
#include "optimization.h"
#include "rateheterogeneity.h"


const double MIN_BRANCH_LEN = 0.000001;
const double MAX_BRANCH_LEN = 9.0;
const double TOL_BRANCH_LEN = 0.00001;
const double TOL_LIKELIHOOD = 0.0001;
const static double SCALING_THRESHOLD = sqrt(DBL_MIN);
const static double SCALING_THRESHOLD_INVER = 1/SCALING_THRESHOLD;
const static double LOG_SCALING_THRESHOLD = log(SCALING_THRESHOLD);

using namespace Eigen;

/**
 *  Row Major Array For Eigen
 */
typedef Array<double, Dynamic, Dynamic, RowMajor> RowMajorArrayXXd;

#define MappedMat(NSTATES) Map<Matrix<double, NSTATES, NSTATES>, Aligned>
#define MappedArr2D(NSTATES) Map<Array<double, NSTATES, NSTATES>, Aligned>
#define MappedRowVec(NSTATES) Map<Matrix<double, 1, NSTATES>, Aligned>
#define MappedVec(NSTATES) Map<Matrix<double, NSTATES, 1>, Aligned>
#define Matrix(NSTATES) Matrix<double, NSTATES, NSTATES>
#define RowVector(NSTATES) Matrix<double, 1, NSTATES>
#define MappedRowArr2DDyn Map<Array<double, Dynamic, Dynamic, RowMajor>, Aligned>

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

struct SPR_compare {

    bool operator()(SPRMove s1, SPRMove s2) const {
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
class PhyloTree : public MTree, public Optimization {
public:
    /**
       default constructor ( everything is initialized to NULL)
     */
    PhyloTree();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Constructor with given alignment
     * @param alignment
     */
    PhyloTree(Alignment *alignment);

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
		copy the sub-tree structure into this tree
		@param tree the tree to copy
		@param taxa_set 0-1 string of length leafNum (1 to keep the leaf)
	*/
	virtual void copyTree(MTree *tree, string &taxa_set);


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
            set the model factory
            @param model_fac model factory
     */
    void setModelFactory(ModelFactory *model_fac);

    /**
            set rate heterogeneity, important to compute the likelihood
            @param rate associated rate heterogeneity class
     */
    void setRate(RateHeterogeneity *rate);

    /**
            get rate heterogeneity
            @return associated rate heterogeneity class
     */
    RateHeterogeneity *getRate();

	void discardSaturatedSite(bool val);

    /**
            get the name of the model
     */
    string getModelName();

    SubstModel *getModel() {
        return model;
    }

    ModelFactory *getModelFactory() {
        return model_factory;
    }

	virtual bool isSuperTree() { return false; }

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
            compute the tree parsimony score
            @return tree likelihood
     */
    int computeParsimony();

    /**
            TODO: Compute partial parsimony score of the subtree rooted at dad
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
     */
    void computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad);


    /**
            compute tree parsimony score on a branch
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @return tree likelihood
     */
    int computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);

	void printParsimonyStates(PhyloNeighbor *dad_branch = NULL, PhyloNode *dad = NULL);


    /**
            SLOW VERSION: compute the parsimony score of the tree, given the alignment
            @return the parsimony score
     */
    int computeParsimonyScore();


    /**
            SLOW VERSION: compute the parsimony score of the tree, given the alignment
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
    virtual void initializeAllPartialLh();

    /**
            initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
            @param node the current node
            @param dad dad of the node, used to direct the search
            @param index the index
     */
    virtual void initializeAllPartialLh(int &index, PhyloNode *node = NULL, PhyloNode *dad = NULL);


    /**
            clear all partial likelihood for a clean computation again
     */
    void clearAllPartialLH();

    /**
            allocate memory for a partial likelihood vector
     */
    double *newPartialLh();

    /**
            compute the partial likelihood at a subtree
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
     */
    void computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL, double *pattern_scale = NULL);

    void computePartialLikelihoodNaive(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL, 
    	double *pattern_scale = NULL);

    template<int NSTATES>
    inline void computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad = NULL, double *pattern_scale = NULL);

    /**
            compute tree likelihood on a branch. used to optimize branch length
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    double computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

    template<int NSTATES>
    inline double computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh = NULL);

    double computeLikelihoodBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, 
    	double *pattern_lh = NULL, double *pattern_rate = NULL);


    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
            Roll back the tree saved with only Taxon IDs and branch lengths.
            For this function to work, one must printTree before with WT_TAXON_ID + WT_BR_LEN
            @param best_tree_string input stream to read from
     */
    void rollBack(istream &best_tree_string);

    bool checkEqualScalingFactor(double &sum_scaling, PhyloNode *node = NULL, PhyloNode *dad = NULL);

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

    double computeLikelihoodDervNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template<int NSTATES>
    inline double computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template<int NSTATES>
    double computeLikelihoodDervSSE_Test(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template<int NSTATES>
    double computeLikelihoodDervSSE_Test2(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template<int NSTATES>
    inline void computeLikelihoodDervSSE_PTN(int ptn, double *partial_lh_site, double *partial_lh_child, double *trans_state, double *derv1_state, double *derv2_state);
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
    int addTaxonMP(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad);


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
    double optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true);

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
    virtual double optimizeAllBranches(int iterations = 100);

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
    double optimizeNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL, ostream *out = NULL, int brtype = 0);

    /**
            This is for ML. try to swap the tree with nearest neigbor interchange at the branch connecting node1-node2.
            If a swap shows better score, return the swapped tree and the score.
            @param cur_score current likelihood score
            @param node1 1st end node of the branch
            @param node2 2nd end node of the branch
            @return the likelihood of the tree
     */
    double swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, ostream *out = NULL, int brtype = 0);


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
            Distance function
     ****************************************************************************/

    /**
            compute the distance between 2 sequences.
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @param initial_dist initial distance
            @return distance between seq1 and seq2
     */
    virtual double computeDist(int seq1, int seq2, double initial_dist);


    /**
            compute distance matrix, assume dist_mat is allocated by memory of size num_seqs * num_seqs.
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @return the longest distance
     */
    double computeDist(double *dist_mat);
    
    /**
            compute observed distance matrix, assume dist_mat is allocated by memory of size num_seqs * num_seqs.
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @return the longest distance
     */
	double computeObsDist(double *dist_mat);

    /**
            compute distance matrix, allocating memory if necessary
            @param params program parameters
            @param alignment input alignment
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @param dist_file (OUT) name of the distance file
            @return the longest distance
     */
    double computeDist(Params &params, Alignment *alignment, double* &dist_mat, string &dist_file);

    /**
            compute observed distance matrix, allocating memory if necessary
            @param params program parameters
            @param alignment input alignment
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @param dist_file (OUT) name of the distance file
            @return the longest distance
	*/
    double computeObsDist(Params &params, Alignment *alignment, double* &dist_mat, string &dist_file);

	/**
		correct the distances to follow metric property of triangle inequalities.
		Using the Floyd alogrithm.
		@param dist_mat (IN/OUT) the shortest path between all pairs of taxa
        @return the longest distance
	*/
	double correctDist(double *dist_mat);

    /****************************************************************************
            compute BioNJ tree, a more accurate extension of Neighbor-Joining
     ****************************************************************************/

    /**
            compute BioNJ tree
            @param params program parameters
            @param alignment input alignment
            @param dist_file distance matrix file
     */
    void computeBioNJ(Params &params, Alignment *alignment, string &dist_file);
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

    /****************************************************************************
            Approximate Likelihood Ratio Test with SH-like interpretation
     ****************************************************************************/

    void computeNNIPatternLh(double cur_lh,
            double &lh2, double *pattern_lh2,
            double &lh3, double *pattern_lh3,
            PhyloNode *node1, PhyloNode *node2);

    /**
            Resampling estimated log-likelihood (RELL)
     */
    void resampleLh(double **pat_lh, double *lh_new);

    /**
            Test one branch of the tree with aLRT SH-like interpretation
     */
    double testOneBranch(
            double best_score, double *pattern_lh,
            int times, PhyloNode *node1, PhyloNode *node2);

    /**
            Test all branches of the tree with aLRT SH-like interpretation
     */
    int testAllBranches(int threshold,
            double best_score, double *pattern_lh,
            int times, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /****************************************************************************
            Collapse stable (highly supported) clades by one representative
     ****************************************************************************/

    /**
            delete a leaf from the tree, assume tree is birfucating
            @param leaf the leaf node to remove
     */
    void deleteLeaf(Node *leaf);

    /**
            reinsert one leaf back into the tree
            @param leaf the leaf to reinsert
            @param adjacent_node the node adjacent to the leaf, returned by deleteLeaves() function
            @param node one end node of the reinsertion branch in the existing tree
            @param dad the other node of the reinsertion branch in the existing tree
     */
    void reinsertLeaf(Node *leaf, Node *node, Node *dad);

    bool isSupportedNode(PhyloNode* node, int min_support);

    /**
            Collapse stable (highly supported) clades by one representative
            @return the number of taxa prunned
     */
    int collapseStableClade(int min_support, NodeVector &pruned_taxa, StrVector &linked_name, double* &dist_mat);

    int restoreStableClade(Alignment *original_aln, NodeVector &pruned_taxa, StrVector &linked_name);

    /****************************************************************************
            Public variables
     ****************************************************************************/

    /**
            associated alignment
     */
    Alignment *aln;

    /**
     *      size of the alignment (used to avoid calling aln->size())
     */
    int alnSize;

    /**
     *      number of states ( used to avoid calling aln->num_states() )
     */
    int numStates;

    /**
            TRUE if you want to optimize branch lengths by Newton-Raphson method
     */
    bool optimize_by_newton;

    /**
     *      TRUE if the loglikelihood is computed using SSE
     */
    bool sse;

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
            Model factory includes SubstModel and RateHeterogeneity
            stores transition matrices computed before for efficiency purpose, eps. AA or CODON model.
     */
    ModelFactory *model_factory;

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

    /**
            the main memory storing all partial parsimony states for all neighbors of the tree.
            The variable partial_pars in PhyloNeighbor will be assigned to a region inside this variable.
     */
    UINT *central_partial_pars;

	/**
		TRUE to discard saturated for Meyer & von Haeseler (2003) model
	*/
	bool discard_saturated_site;
    
    /**
     * Number of category
     */
    int numCat;

    /**
     * Size of likelihood block for 1 pattern
     */
    int block;

    /**
     * Size of the partial likelihood (multiple of double)
     */
    int lh_size;

    int block_size;

    /**
     * Vector to store pattern frequency
     */
    ArrayXd ptn_freqs;

    /**
     * Size of the transition matrix
     */
    int tranSize;

    /**
     * Temporary partial likelihood array: used when swapping branch and recalculate the
     * likelihood --> avoid calling malloc
     */
    double *tmp_partial_lh1;
    double *tmp_partial_lh2;

    /****************************************************************************
            Vector of bit blocks, used for parsimony function
     ****************************************************************************/

    /**
            @return size of the bits block vector for one node
     */
    int getBitsBlockSize();

    /**
            allocate new memory for a bit block vector
            @return the allocated memory
     */
    UINT *newBitsBlock();

    /**
            @return size of the bits entry (for storing num_states bits)
     */
    int getBitsEntrySize();

	/**
		@param bits_entry
		@return TRUE if bits_entry contains all 0s, FALSE otherwise
	*/
	bool isEmptyBitsEntry(UINT *bits_entry);

	/**
		@param bits_entry1
		@param bits_entry1
		@param bits_union (OUT) union of bits_entry1 and bits_entry2
	*/
	void unionBitsEntry(UINT *bits_entry1, UINT *bits_entry2, UINT* &bits_union);

	/**
		set a single bit to 1
		@param bits_entry
		@param id index of the bit in the entry to set to 1
	*/
	void setBitsEntry(UINT* &bits_entry, int id);

	/**
		get a single bit content
		@param bits_entry
		@param id index of the bit in the entry
		@return TRUE if bit ID is 1, FALSE otherwise
	*/
	bool getBitsEntry(UINT* &bits_entry, int id);

    /**
            get bit blocks, each block span num_state bits
            @param bit_vec bit block vector
            @param index block index
            @param bits_entry (OUT) content of the block at index
     */
    void getBitsBlock(UINT *bit_vec, int index, UINT* &bits_entry);

    /**
            set bit blocks, each block span num_state bits
            @param bit_vec (OUT) bit block vector
            @param index block index
            @param bits_entry the content of the block at index
     */
    void setBitsBlock(UINT* &bit_vec, int index, UINT *bits_entry);
};

#endif
