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

#ifndef PHYLOTREE_H
#define PHYLOTREE_H
//#define NDEBUG
// comented out this for Mac

// PLEASE DONT TOUCH THESE VARIABLES ANYMORE!
#define EIGEN_NO_AUTOMATIC_RESIZING
//#define EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD 32 // PLEASE DONT TOUCH THESE VARIABLES ANYMORE!
//#define EIGEN_UNROLLING_LIMIT 1000 // PLEASE DONT TOUCH THESE VARIABLES ANYMORE!

//#define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (512*256)
//#define EIGEN_TUNE_FOR_CPU_CACHE_SIZE (8*512*512)
//#include <Eigen/Core>

#define EIGEN_PARTIAL_LIKELIHOOD (0) //Set to 1 to turn it all back on

#include "mtree.h"
#include "alignment/alignment.h"
#include "alignment/alignmentsummary.h"
#include "model/modelsubst.h"
#include "model/modelfactory.h"
#include "phylonode.h"
#include "utils/optimization.h"
#include "model/rateheterogeneity.h"
#include "candidateset.h"
#include "pll/pll.h"
#include "utils/checkpoint.h"
#include "constrainttree.h"
#include "memslot.h"
#include "utils/progress.h"
#include "alignedalloc.h"
#include "likelihoodbufferset.h"

class AlignmentPairwise;

#define BOOT_VAL_FLOAT
#define BootValType float
//#define BootValType double

enum CostMatrixType {CM_UNIFORM, CM_LINEAR};

//extern int instruction_set;

#define SAFE_LH   true  // safe likelihood scaling to avoid numerical underflow for ultra large trees
#define NORM_LH  false // normal likelihood scaling

const double TOL_BRANCH_LEN = 0.000001; // NEVER TOUCH THIS CONSTANT AGAIN PLEASE!
const double TOL_LIKELIHOOD = 0.001; // NEVER TOUCH THIS CONSTANT AGAIN PLEASE!
const double TOL_LIKELIHOOD_PARAMOPT = 0.001; // BQM: newly introduced for ModelFactory::optimizeParameters
//const static double SCALING_THRESHOLD = sqrt(DBL_MIN);
//const static double SCALING_THRESHOLD = 1e-100;
//const static double SCALING_THRESHOLD_INVER = 1 / SCALING_THRESHOLD;
//const static double LOG_SCALING_THRESHOLD = log(SCALING_THRESHOLD);
//#include "pll/pll.h"
// 2^256
//#define SCALING_THRESHOLD_INVER 115792089237316195423570985008687907853269984665640564039457584007913129639936.0
#define SCALING_THRESHOLD_EXP 256
//#define SCALING_THRESHOLD (1.0/SCALING_THRESHOLD_INVER)
// 2^{-256}
#define SCALING_THRESHOLD 8.636168555094444625386e-78
//#define SCALING_THRESHOLD ldexp(1.0, -256)
//#define LOG_SCALING_THRESHOLD log(SCALING_THRESHOLD)
#define LOG_SCALING_THRESHOLD -177.4456782233459932741

const int SPR_DEPTH = 2;

//using namespace Eigen;


#define LOG_LINE(lev,text) \
    if (verbose_mode >= (lev)) { \
        std::stringstream s; \
        s << text; \
        logLine(s.str()); \
    } else 0

#define TREE_LOG_LINE(t, lev, text) \
    if (verbose_mode >= (lev)) { \
        std::stringstream s; \
        s << text; \
        (t).logLine(s.str()); \
    } else 0

/**
 *  Row Major Array For Eigen
 */
//typedef Array<double, Dynamic, Dynamic, RowMajor> RowMajorArrayXXd;


typedef std::map< string, double > StringDoubleMap;
typedef std::map< int, PhyloNode* > IntPhyloNodeMap;

#define DUMMY_NODE_1 ((PhyloNode*)1)
#define DUMMY_NODE_2 ((PhyloNode*)2)

/*
#define MappedMat(NSTATES) Map<Matrix<double, NSTATES, NSTATES> >
#define MappedArr2D(NSTATES) Map<Array<double, NSTATES, NSTATES> >
#define MappedRowVec(NSTATES) Map<Matrix<double, 1, NSTATES> >
#define MappedVec(NSTATES) Map<Matrix<double, NSTATES, 1> >
#define Matrix(NSTATES) Matrix<double, NSTATES, NSTATES>
#define RowVector(NSTATES) Matrix<double, 1, NSTATES>
#define MappedRowArr2DDyn Map<Array<double, Dynamic, Dynamic, RowMajor> >
#define MappedArrDyn Map<Array<double, Dynamic, 1> >
#define MappedVecDyn(NSTATES) Map<Matrix<double, Dynamic, NSTATES> >
*/

const int MAX_SPR_MOVES = 20;

struct NNIMove {

    // Two nodes representing the central branch
    PhyloNode *node1, *node2;

    // Roots of the two subtree that are swapped
    NeighborVec::iterator node1Nei_it, node2Nei_it;

    // log-likelihood of the tree after applying the NNI
    double newloglh;

    int swap_id;

    // new branch lengths of 5 branches corresponding to the NNI
    DoubleVector newLen[5];

    // pattern likelihoods
    double *ptnlh;

    NNIMove(): node1(nullptr), node2(nullptr)
        , newloglh(-DBL_MAX), ptnlh(nullptr) {
    }
    bool operator<(const NNIMove & rhs) const {
        return newloglh > rhs.newloglh;
    }
};

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
                     |inline
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
 * This Structure is used in PhyloSuperTreePlen.
 */
struct SwapNNIParam {
    double nni1_score;
    double nni1_brlen;
    double nni2_score;
    double nni2_brlen;
    Neighbor* node1_nei;
    Neighbor* node2_nei;
    double *nni1_ptnlh;
    double *nni2_ptnlh;
};


struct LeafFreq {
    int leaf_id;

    int freq;

    bool operator<(const LeafFreq & rhs) const {
        return ( freq < rhs.freq);
    }
};


// **********************************************
// BEGIN definitions for likelihood mapping (HAS)
// **********************************************

/* maximum exp difference, such that 1.0+exp(-TP_MAX_EXP_DIFF) == 1.0 */
const double TP_MAX_EXP_DIFF = 40.0;

/* Index definition for counter array needed in likelihood mapping analysis (HAS) */
#define LM_REG1 0
#define LM_REG2 1
#define LM_REG3 2
#define LM_REG4 3
#define LM_REG5 4
#define LM_REG6 5
#define LM_REG7 6
#define LM_AR1  7
#define LM_AR2  8
#define LM_AR3  9
#define LM_MAX  10

struct QuartetGroups{
    int numGroups;		// number of clusters:
				// 0:	not initialized, default -> 1
				// 1:	no clusters - any (a,b)|(c,d)
				// 2:	2 clusters  - (a,a')|(b,b')
				// 3:	3 clusters  - (a,a')|(b,c)	[rare]
				// 4:	4 clusters  - (a,b)|(c,d)
    int numSeqs;		// number of seqs in alignment (should be #A+#B+#C+#D+#X)
    int numQuartSeqs;		// number of seqs in analysis  (should be #A+#B+#C+#D)
    int numGrpSeqs[5];		// number of seqs in cluster A, B, C, D, and X (exclude)
    int64_t uniqueQuarts;	// number of existing unique quartets for this grouping
    string Name[5];		// seqIDs of cluster A
    vector<int> GroupA;		// seqIDs of cluster A
    vector<int> GroupB;		// seqIDs of cluster B
    vector<int> GroupC;		// seqIDs of cluster C
    vector<int> GroupD;		// seqIDs of cluster D
    vector<int> GroupX;		// seqIDs of cluster X
};

struct QuartetInfo {
    int seqID[4];
    double logl[3];    // log-lh for {0,1}|{2,3}  {0,2}|{1,3}  {0,3}|{1,4}
    double qweight[3]; // weight for {0,1}|{2,3}  {0,2}|{1,3}  {0,3}|{1,4}
    int corner;        // for the 3 corners of the simplex triangle (0:top, 1:right, 2:left)
    int area;          // for the 7 areas of the simplex triangle
			// corners (0:top, 1:right, 2:left), rectangles (3:right, 4:left, 5:bottom), 6:center
};

struct SeqQuartetInfo {
    int64_t countarr[LM_MAX]; // the 7 areas of the simplex triangle [0-6; corners (0:top, 1:right, 2:left), rectangles (3:right, 4:left, 5:bottom), 6:center] and the 3 corners [7-9; 7:top, 8:right, 9:left]
};

// ********************************************
// END definitions for likelihood mapping (HAS)
// ********************************************


// ********************************************
// BEGIN traversal information
// ********************************************

class TraversalInfo {
public:
    PhyloNeighbor* dad_branch;
    PhyloNode*     dad;
    double*        echildren;
    double*        partial_lh_leaves;

    TraversalInfo(PhyloNeighbor *dad_branch, PhyloNode *dad) {
        this->dad                  = dad;
        this->dad_branch           = dad_branch;
        echildren                  = nullptr;
        partial_lh_leaves          = nullptr;
    }
    TraversalInfo(const TraversalInfo& rhs)
        : dad_branch(rhs.dad_branch), dad(rhs.dad)
        , echildren(rhs.echildren)
        , partial_lh_leaves(rhs.partial_lh_leaves) {
     }
    TraversalInfo& operator=(const TraversalInfo& rhs) {
        TraversalInfo dummy(rhs);
        swapWith(dummy);
        return *this;
    }
    void swapWith(TraversalInfo& rhs) {
        std::swap(dad_branch,        rhs.dad_branch);
        std::swap(dad,               rhs.dad);
        std::swap(echildren,         rhs.echildren);
        std::swap(partial_lh_leaves, rhs.partial_lh_leaves);
    }
};

// ********************************************
// END traversal information
// ********************************************

/**
Phylogenetic Tree class

        @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloTree : public MTree, public Optimization, public CheckpointFactory {

	friend class PhyloSuperTree;
	friend class PhyloSuperTreePlen;
	friend class RateGamma;
	friend class RateGammaInvar;
	friend class RateKategory;
    friend class ModelMixture;
    friend class RateFree;
    friend class RateHeterotachy;
    friend class PhyloTreeMixlen;
    friend class ModelFactoryMixlen;
    friend class MemSlotVector;
    friend class ModelFactory;
    friend class CandidateSet;
    friend class TaxonToPlace;
    friend class TargetBranch;
    friend class BlockAllocator;
    friend class PlacementTraversalInfo;
    friend class LikelihoodBlockAllocator;

public:
    typedef MTree super;
    
    /**
       default constructor ( everything is initialized to NULL)
     */
    PhyloTree();

//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Constructor with given alignment
     * @param alignment
     */
    PhyloTree(Alignment *aln);

    /**
     *  Create a phylotree from the tree string and assign alignment.
     *  Taxa IDs are numbered according to their orders in the alignment.
     */
    PhyloTree(string& treeString, Alignment *aln, bool isRooted);

    void init();

    /**
            destructor (and supporting deallocation functions)
     */
    virtual ~PhyloTree();

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /** 
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    /**
            read the tree from the input file in newick format
            @param infile the input file file.
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(const char *infile, bool &is_rooted);

    /**
            read the tree from the ifstream in newick format
            @param in the input stream.
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    virtual void readTree(istream &in, bool &is_rooted);

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
     copy the constraint tree structure into this tree and reindex node IDs accordingly
     @param tree the tree to copy
     @param[out] order of taxa with first part being in constraint tree
     */
    void copyConstraintTree(MTree *tree, IntVector &taxon_order, int *rand_stream);

    /**
            copy the phylogenetic tree structure into this tree, designed specifically for PhyloTree.
            So there is some distinction with copyTree.
            @param tree the tree to copy
            @param borrowSummary true if the alignment summary of the other tree is to be "borrowed"
     */
    void copyPhyloTree(PhyloTree *tree, bool borrowSummary);

    /**
            copy the phylogenetic tree structure into this tree, designed specifically for PhyloTree.
            So there is some distinction with copyTree.
            @param tree the tree to copy
            @param mix mixture ID of branch lengths
            @param borrowSummary true if the alignment summary of the other tree is to be "borrowed"
     */
    virtual void copyPhyloTreeMixlen(PhyloTree *tree, int mix, bool borrowSummary);


    /**
            Set the alignment, important to compute parsimony or likelihood score
            Assing taxa ids according to their position in the alignment
            @param alignment associated alignment
     */
    virtual void setAlignment(Alignment* alignment);

    void configureLikelihoodKernel(const Params& params);
    
    void configureModel(Params& params);
    
    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block);
    
    /*
     Modify the tree by marking a subset of the taxa as
    removed (if requested to do so).*/
    void removeSampleTaxaIfRequested();
    
    /**
     Modify the tree to match an alignment, resolving differences
     (by deleting nodes for taxa not found in the alignment,
     by adding nodes for taxa found in the alignment, but not in the tree)
     Nodes are assigned taxa ids according to their position in the alignment
     @param alignment associated alignment
     @return true if the tree had to be updated (nodes added or removed))
     */
    virtual bool updateToMatchAlignment(Alignment *alignment);

    bool shouldPlacementUseSankoffParsimony() const;
    bool shouldPlacementUseLikelihood() const;
    /**
     Modify a tree, by adding nodes for taxa found in the tree's alignment,
     but not in the tree, and assigning them taxa ids according to their
     position in the alignment
     @param taxaIdsToAdd
     @param index_parsimony indicates how many partial parsimony blocks
      have been allocated to nodes already in the tree
     @param index_lh indicates how many partial likelhood blocks
     have been allocated to nodes already in the tree
     */
    void addNewTaxaToTree(const IntVector& taxaIdsToAdd);
    
    double taxaAdditionWorkEstimate(size_t newTaxaCount, size_t taxaPerBatch, size_t insertsPerBatch);
    
    void finishUpAfterTaxaAddition();
    
    virtual PhyloNode* findFarthestLeaf(PhyloNode *node = nullptr,
                                       PhyloNode *dad = nullptr);
    
    /** set the root by name
        @param my_root root node name
        @param multi_taxa TRUE if my_root is a comma-separated list of nodes
     */
    virtual void setRootNode(const char *my_root, bool multi_taxa = false);


    /**
            set the substitution model, important to compute the likelihood
            @param amodel associated substitution model
     */
    void setModel(ModelSubst *amodel);

    /**
            set the model factory
            @param model_fac model factory
     */
    virtual void setModelFactory(ModelFactory *model_fac);

    /**
            indicates whether there is a model factory
            @return true if there is an associated model factory
     */
    bool hasModelFactory() const;
    
    /**
            set rate heterogeneity, important to compute the likelihood
            @param rate associated rate heterogeneity class
     */
    void setRate(RateHeterogeneity *rate);

    /**
            indicates whether there is a rate heterogeneity
            @return true if there is an associated rate heterogeneity
     */
    bool hasRateHeterogeneity() const;
    
    /**
            get rate heterogeneity
            @return associated rate heterogeneity class
     */
    RateHeterogeneity *getRate();

    void discardSaturatedSite(bool val);

    /** get substitution matrix name */
    string getSubstName();
    
    /** get rate heterogeneity name */
    string getRateName();
    
    /**
            get the name of the model
     */
    virtual string getModelName();

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}+I{pinvar}+G{alpha}
	 */
	virtual string getModelNameParams();

    ModelSubst *getModel() {
        return model;
    }

    ModelFactory *getModelFactory() {
        return model_factory;
    }

    virtual bool isSuperTree() {
        return false;
    }

    virtual bool isSuperTreeUnlinked() {
        return false;
    }

    /**
        @return true if this is a tree with mixture branch lengths, default: false
    */
    virtual bool isMixlen() { return false; }

    /**
        @return number of mixture branch lengths, default: 1
    */
    virtual int getMixlen() { return 1; }
    
    virtual PhyloNode* getRoot();

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name
            @return a new node
     */
    virtual PhyloNode* newNode(int node_id = -1, const char* node_name = NULL);

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name issued by an interger
            @return a new node
     */
    virtual PhyloNode* newNode(int node_id, int node_name);

    /**
            indicates whether a given node is a dummy "place-holder" node,
            used only for testing the placement cost of adding a new taxon,
            to an existing tree, by modifying the tree (for example, in 
            addTaxonML and growTreeML, in reinsertLeavesByParsimony, and in
            insertNode2Branch, addTaxonMPFast, and computeParsimonyTree).  

            @param node the node
            @return true if the node is a dummy "place-holder", false otherwise.
     */
    virtual bool isDummyNode(PhyloNode* node) const;

    /**
     *		@return number of alignment patterns
     */
    virtual size_t getAlnNPattern() const {
        return aln->getNPattern();
    }

    /**
     *		@return number of alignment sites
     */
    virtual size_t getAlnNSite() {
        return aln->getNSite();
    }

    virtual void computeBranchLengths() {
    }
    
    /**
     * save branch lengths into a vector
     */
    virtual void saveBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = NULL, PhyloNode *dad = NULL);
    /**
     * restore branch lengths from a vector previously called with saveBranchLengths
     */
    virtual void restoreBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /****************************************************************************
            Dot product
     ****************************************************************************/
    template <class Numeric, class VectorClass>
    Numeric dotProductSIMD(Numeric *x, Numeric *y, int size);

    typedef BootValType (PhyloTree::*DotProductType)(BootValType *x, BootValType *y, int size);
    DotProductType dotProduct;

    typedef double (PhyloTree::*DotProductDoubleType)(double *x, double *y, int size);
    DotProductDoubleType dotProductDouble;

    double dotProductDoubleCall(double *x, double *y, int size);

#if defined(BINARY32) || defined(__NOAVX__)
    void setDotProductAVX() {}
    void setDotProductFMA() {}
#else
    void setDotProductAVX();
    void setDotProductFMA();
    void setDotProductAVX512();
#endif

    void setDotProductSSE();

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
     * 		Return the approximated branch length estimation using corrected parsimony branch length
     * 		This is usually used as the starting point before using Newton-Raphson
     */
//    double computeCorrectedParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
            initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
            @return the number of partial parsimony blocks that were used by existing nodes of the tree
     */
    virtual int initializeAllPartialPars();

    /**
            initialize partial_pars vector of all PhyloNeighbors, allocating central_partial_pars
            @param node the current node
            @param dad dad of the node, used to direct the search
            @param index the index
     */
    virtual void initializeAllPartialPars(int &index, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            compute the tree parsimony score
            @return parsimony score of the tree
     */
    int computeParsimony(const char* taskDescription = "");

    typedef void (PhyloTree::*ComputePartialParsimonyType)(PhyloNeighbor *, PhyloNode *);
    ComputePartialParsimonyType computePartialParsimonyPointer;
    
    typedef void (PhyloTree::*ComputePartialParsimonyOutOfTreeType)(const UINT* left_partial_pars,
                                                     const UINT* right_partial_pars,
                                                     UINT* dad_partial_pars) const;
    ComputePartialParsimonyOutOfTreeType computePartialParsimonyOutOfTreePointer;


    /**
            Compute partial parsimony score of the subtree rooted at dad
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
     */
    virtual void computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad);
//    void computePartialParsimonyNaive(PhyloNeighbor *dad_branch, PhyloNode *dad);
    void computePartialParsimonyFast(PhyloNeighbor *dad_branch, PhyloNode *dad);
    
    void computePartialParsimonyOutOfTreeFast(const UINT* left_partial_pars,
                                              const UINT* right_partial_pars,
                                              UINT* dad_partial_pars) const;

    template<class VectorClass>
    void computePartialParsimonyFastSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);
    
    
    void computePartialParsimonyOutOfTree(const UINT* left_partial_pars,
                                          const UINT* right_partial_pars,
                                          UINT* dad_partial_pars) const;
    template<class VectorClass>
    void computePartialParsimonyOutOfTreeSIMD(const UINT* left_partial_pars,
                                              const UINT* right_partial_pars,
                                              UINT* dad_partial_pars) const;

    template<class VectorClass>
    void computePartialParsimonySankoffSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);

    template<class VectorClass>
    void computePartialParsimonyOutOfTreeSankoffSIMD(const UINT* left_partial_pars,
                                                     const UINT* right_partial_pars,
                                                     UINT*       dad_partial_pars) const;

    void computeReversePartialParsimony(PhyloNode *node, PhyloNode *dad);

    typedef int (PhyloTree::*ComputeParsimonyBranchType)(PhyloNeighbor *, PhyloNode *, int *);
    ComputeParsimonyBranchType computeParsimonyBranchPointer;

    typedef int (PhyloTree::*ComputeParsimonyOutOfTreeType)(const UINT* dad_partial_pars,
                                             const UINT* node_partial_pars,
                                             int* branch_subst) const;
    ComputeParsimonyOutOfTreeType computeParsimonyOutOfTreePointer;
    /**
            compute tree parsimony score on a branch
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
            @return parsimony score of the tree
     */
    virtual int computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);
    
    virtual int computeParsimonyOutOfTree(const UINT* dad_partial_pars,
                                          const UINT* node_partial_pars,
                                          int* branch_subst = nullptr) const;

    
//    int computeParsimonyBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);
    int computeParsimonyBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);
    int computeParsimonyOutOfTreeFast(const UINT* dad_partial_pars,
                                      const UINT* node_partial_pars,
                                      int*        branch_subst) const;
    
    template<class VectorClass>
        int computeParsimonyBranchFastSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);
    
    template<class VectorClass>
    int computeParsimonyOutOfTreeSIMD(const UINT* dad_partial_pars,
                                      const UINT* node_partial_pars,
                                      int*        branch_subst) const;

    template<class VectorClass>
    int computeParsimonyBranchSankoffSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = NULL);

    template<class VectorClass>
    int computeParsimonyOutOfTreeSankoffSIMD(const UINT* dad_partial_pars,
                                             const UINT* node_partial_pars,
                                             int*        branch_subst) const;

    //    void printParsimonyStates(PhyloNeighbor *dad_branch = NULL, PhyloNode *dad = NULL);

    virtual void setParsimonyKernel(LikelihoodKernel lk);
#if defined(BINARY32) || defined(__NOAVX__)
    virtual void setParsimonyKernelAVX() {}
#else
    virtual void setParsimonyKernelAVX();
#endif

    virtual void setParsimonyKernelSSE();

    /****************************************************************************
     Sankoff Parsimony function
     ****************************************************************************/

    /**
     * initialise cost_matrix as linear
     * initialize for 'nstates' and 'columns'
     */
    void initCostMatrix(CostMatrixType cost_type);
        
    /**
     * read the cost matrix file
     * initialize for 'nstates' and 'columns'
     */
    void loadCostMatrixFile(char* file_name = NULL);

    /*
     * For a leaf character corresponding to an ambiguous state
     * set elements corresponding to possible states to 0, others to UINT_MAX
     */
    void initLeafSiteParsForAmbiguousState(char state, UINT * site_partial_pars);

    /**
     compute partial parsimony score of the subtree rooted at dad
     @param dad_branch the branch leading to the subtree
     @param dad its dad, used to direct the traversal
     */
    void computePartialParsimonySankoff(PhyloNeighbor *dad_branch, PhyloNode *dad);
    
    void computePartialParsimonyOutOfTreeSankoff(const UINT* left_partial_pars,
                                                 const UINT* right_partial_pars,
                                                 UINT* dad_partial_pars) const;
    
    /**
     compute tree parsimony score based on a particular branch
     @param dad_branch the branch leading to the subtree
     @param dad its dad, used to direct the traversal
     @param branch_subst (OUT) if not NULL, the number of substitutions on this branch
     @return parsimony score of the tree
     */
    int computeParsimonyBranchSankoff(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                      int *branch_subst = NULL);
    
    int computeParsimonyOutOfTreeSankoff(const UINT* dad_partial_pars,
                                         const UINT* node_partial_pars,
                                         int* branch_subst = nullptr) const;
    
    /****************************************************************************
            likelihood function
     ****************************************************************************/

    size_t getBufferPartialLhSize();

    /**
            initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
     */
    virtual void initializeAllPartialLh();

    /**
            de-allocate central_partial_lh
     */
    virtual void deleteAllPartialLh();

    virtual void allocateCentralBlocks(size_t extra_parsimony_block_count,
                                       size_t extra_lh_block_count);
    
    virtual void getBlockSizes(size_t& nptn, uint64_t& pars_block_size,
                               uint64_t& lh_block_size, uint64_t& scale_block_size );
    //
    //Note: This probably needs to be passed a partition number too, if it is to work
    //      properly for PhyloSuperTreePlen.
    //      See... PhyloSuperTreePlen::initializeAllPartialLh
    //
    
    void ensurePartialLHIsAllocated(size_t count_of_extra_parsimony_blocks,
                                    size_t count_of_extra_lh_blocks);
        
    /**
            initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
            @param node the current node
            @param dad dad of the node, used to direct the search
            @param index the index
     */
    virtual void initializeAllPartialLh(int &index, int &indexlh,
                                        PhyloNode *node = NULL, PhyloNode *dad = NULL);


    /**
            clear all partial likelihood for a clean computation again
            @param set_to_null true to make all partial_lh become NULL
     */
    virtual void clearAllPartialLH(bool set_to_null = false);
    
    /**
            clear all partial parsimony data for a clean computation again
     @param set_to_null true to make all partial_pars become NULL
     */
    virtual void clearAllPartialParsimony(bool set_to_null);

    /**
            clear all scale number data for a clean computation again
            @param set_to_null true to make all scale_num become NULL
     */
    virtual void clearAllScaleNum(bool set_to_null);

    /**
     * compute all partial likelihoods if not computed before
     */
    void computeAllPartialLh(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * compute all partial parsimony vector if not computed before
     */
    void computeAllPartialPars(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
            allocate memory for a partial likelihood vector
     */
    double *newPartialLh();

    /** get the number of bytes occupied by partial_lh */
    size_t getPartialLhBytes();
    size_t getPartialLhSize();

    /**
            allocate memory for a scale num vector
     */
    UBYTE *newScaleNum();

    /** get the number of bytes occupied by scale_num */
    size_t getScaleNumBytes();
    size_t getScaleNumSize();

    /**
     * this stores partial_lh for each state at the leaves of the tree because they are the same between leaves
     * e.g. (1,0,0,0) for A,  (0,0,0,1) for T
     */
    double *tip_partial_lh;
    int tip_partial_lh_computed;
    UINT *tip_partial_pars;

    bool ptn_freq_computed;

    /** site log-likelihood buffer for robust phylogeny idea */
    double *_site_lh;
    
    /** vector size used by SIMD kernel */
    size_t vector_size;

    /** true if using safe numeric for likelihood kernel */
    bool safe_numeric;

    /** number of threads used for likelihood kernel */
    int num_threads;

    /** number of packets used for likelihood kernel (typically more) */
    int num_packets;

    /****************************************************************************
            helper functions for computing tree traversal
     ****************************************************************************/


    /**
        compute traversal_info of a subtree
    */
    bool computeTraversalInfo(PhyloNeighbor *dad_branch, PhyloNode *dad, double* &buffer);


    /**
        compute traversal_info of both subtrees
    */
    template<class VectorClass, const int nstates>
    void computeTraversalInfo(PhyloNode *node, PhyloNode *dad, LikelihoodBufferSet& buffers, bool compute_partial_lh);
    template<class VectorClass>
    void computeTraversalInfo(PhyloNode *node, PhyloNode *dad, LikelihoodBufferSet& buffers, bool compute_partial_lh);

    /**
        precompute info for models
    */
    template<class VectorClass, const int nstates>
    void computePartialInfo(TraversalInfo &info, VectorClass* buffer);
    template<class VectorClass>
    void computePartialInfo(TraversalInfo &info, VectorClass* buffer);

    /** 
        sort neighbor in descending order of subtree size (number of leaves within subree)
        @param node the starting node, NULL to start from the root
        @param dad dad of the node, used to direct the search
    */
    void sortNeighborBySubtreeSize(PhyloNode *node, PhyloNode *dad);

    /****************************************************************************
            computing partial (conditional) likelihood of subtrees
     ****************************************************************************/

    /** transform _pattern_lh_cat from "interleaved" to "sequential", due to vector_size > 1 */
    void transformPatternLhCat();

    // Compute the partial likelihoods LH (OUT) at the leaves for an observed PoMo
    // STATE (IN). Use binomial sampling unless hyper is true, then use
    // hypergeometric sampling.
    void computeTipPartialLikelihoodPoMo(int state, double *lh, bool hypergeometric=false);
    void computeTipPartialLikelihood();
    void computeTipPartialParsimony();
    void computePtnInvar();
    void computePtnFreq();
    
    void computePatternPacketBounds(int vector_size, int threads, int packets,
                                    size_t elements, vector<size_t> &limits);

    /**
            compute the partial likelihood at a subtree
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
     */
    virtual void computePartialLikelihood(TraversalInfo &info,
                                          size_t ptn_lower, size_t ptn_upper, int thread_id,
                                          LikelihoodBufferSet& buffers);
    typedef void (PhyloTree::*ComputePartialLikelihoodType)(TraversalInfo &info,
                                                            size_t ptn_lower, size_t ptn_upper, int thread_id,
                                                            LikelihoodBufferSet& buffers);
    ComputePartialLikelihoodType computePartialLikelihoodPointer;


    #if (EIGEN_PARTIAL_LIKELIHOOD)
    template <const int nstates>
    void computePartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                       LikelihoodBufferSet& buffers);

    void computeSitemodelPartialLikelihoodEigen(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                LikelihoodBufferSet& buffers);

    template <class VectorClass, const int VCSIZE, const int nstates>
    void computePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                           LikelihoodBufferSet& buffers);
    #endif
    
    void computeNonrevPartialLikelihood(TraversalInfo &info,
                                        size_t ptn_lower, size_t ptn_upper, int thread_id,
                                        LikelihoodBufferSet& buffers);
    template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA = false>
    void computeNonrevPartialLikelihoodGenericSIMD(TraversalInfo &info,
                                                   size_t ptn_lower, size_t ptn_upper, int thread_id,
                                                   LikelihoodBufferSet& buffers);
    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA = false>
    void computeNonrevPartialLikelihoodSIMD(TraversalInfo &info,
                                            size_t ptn_lower, size_t ptn_upper, int thread_id,
                                            LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA = false, const bool SITE_MODEL = false>
    void computePartialLikelihoodSIMD(TraversalInfo &info,
                                      size_t ptn_lower, size_t ptn_upper, int thread_id,
                                      LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA = false, const bool SITE_MODEL = false>
    void computePartialLikelihoodGenericSIMD(TraversalInfo &info,
                                             size_t ptn_lower, size_t ptn_upper, int thread_id,
                                             LikelihoodBufferSet& buffers);

    #if (EIGEN_PARTIAL_LIKELIHOOD)
    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeMixratePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                  LikelihoodBufferSet& buffers);

    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeMixturePartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                  LikelihoodBufferSet& buffers);

    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeSitemodelPartialLikelihoodEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad
                                                    LikelihoodBufferSet& buffers);
    #endif

    /****************************************************************************
            computing likelihood on a branch
     ****************************************************************************/

    /**
            compute tree likelihood on a branch. used to optimize branch length
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @return tree likelihood
     */
    virtual double computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                           LikelihoodBufferSet& buffers);

    typedef double (PhyloTree::*ComputeLikelihoodBranchType)(PhyloNeighbor*, PhyloNode*,
                                                             LikelihoodBufferSet&);
    ComputeLikelihoodBranchType computeLikelihoodBranchPointer;

    /**
     * MINH: this implements the fast alternative strategy for reversible model (March 2013)
     * where partial likelihoods at nodes store real partial likelihoods times eigenvectors
     */
//    template<int NSTATES>
//    inline double computeLikelihoodBranchFast(PhyloNeighbor *dad_branch, PhyloNode *dad);

    //template <const int nstates>
//    double computeLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad);

//    double computeSitemodelLikelihoodBranchEigen(PhyloNeighbor *dad_branch, PhyloNode *dad);

//    template <class VectorClass, const int VCSIZE, const int nstates>
//    double computeLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);

    double computeNonrevLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                         LikelihoodBufferSet& buffers);
    template<class VectorClass, const bool SAFE_NUMERIC, const bool FMA = false>
    double computeNonrevLikelihoodBranchGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                    LikelihoodBufferSet& buffers);
    template<class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA = false>
    double computeNonrevLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                             LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates, const bool FMA = false, const bool SITE_MODEL = false>
    double computeLikelihoodBranchSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                       LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const bool FMA = false, const bool SITE_MODEL = false>
    double computeLikelihoodBranchGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                              LikelihoodBufferSet& buffers);

    /*
    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeMixrateLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);

    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeMixtureLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);

    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeSitemodelLikelihoodBranchEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad);
    */

    /****************************************************************************
            computing likelihood on a branch using buffer
     ****************************************************************************/

    /**
            quickly compute tree likelihood on branch current_it <-> current_it_back given buffer (theta_all).
           	Used after optimizing branch length
            @param pattern_lh (OUT) if not NULL, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihoodFromBuffer();
    typedef double (PhyloTree::*ComputeLikelihoodFromBufferType)(LikelihoodBufferSet&);
    ComputeLikelihoodFromBufferType computeLikelihoodFromBufferPointer;

//    template <class VectorClass, const int VCSIZE, const int nstates>
//    double computeLikelihoodFromBufferEigenSIMD();

    template <class VectorClass, const int nstates, const bool FMA = false, const bool SITE_MODEL = false>
    double computeLikelihoodFromBufferSIMD(LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool FMA = false, const bool SITE_MODEL = false>
    double computeLikelihoodFromBufferGenericSIMD(LikelihoodBufferSet& buffers);

    /*
    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeMixrateLikelihoodFromBufferEigenSIMD();

    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeMixtureLikelihoodFromBufferEigenSIMD();

    template <class VectorClass, const int VCSIZE, const int nstates>
    double computeSitemodelLikelihoodFromBufferEigenSIMD();

    double computeSitemodelLikelihoodFromBufferEigen();
    */

    /**
            compute tree likelihood when a branch length collapses to zero
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @return tree likelihood
     */
    virtual double computeLikelihoodZeroBranch(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
        compute likelihood of rooted tree with virtual root (FOR TINA)
        @param dad_branch the branch leading to the subtree
        @param dad its dad, used to direct the tranversal
        @return tree likelihood
     */
//    virtual double computeLikelihoodRooted(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not nullptr, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
     * @return number of elements per site lhl entry, used in conjunction with computePatternLhCat
     */
    virtual int getNumLhCat(SiteLoglType wsl);

    /**
     * compute _pattern_lh_cat for site-likelihood per category
     * @return tree log-likelihood
     */
    virtual double computePatternLhCat(SiteLoglType wsl);

    /**
        compute state frequency for each pattern (for Huaichun)
        @param[out] ptn_state_freq state frequency vector per pattern, 
            should be pre-allocated with size of num_patterns * num_states
    */
    void computePatternStateFreq(double *ptn_state_freq);

    /****************************************************************************
            ancestral sequence reconstruction
     ****************************************************************************/

    /**
        initialize computing ancestral sequence probability for an internal node by marginal reconstruction
    */
    virtual void initMarginalAncestralState(ostream &out, bool &orig_kernel_nonrev, double* &ptn_ancestral_prob, int* &ptn_ancestral_seq);

    /**
        compute ancestral sequence probability for an internal node by marginal reconstruction
        (Yang, Kumar and Nei 1995)
        @param dad_branch branch leading to an internal node where to obtain ancestral sequence
        @param dad dad of the target internal node
        @param[out] ptn_ancestral_prob pattern ancestral probability vector of dad_branch->node
    */
    virtual void computeMarginalAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad,
        double *ptn_ancestral_prob, int *ptn_ancestral_seq);

    virtual void writeMarginalAncestralState(ostream &out, PhyloNode *node, double *ptn_ancestral_prob, int *ptn_ancestral_seq);

    /**
        end computing ancestral sequence probability for an internal node by marginal reconstruction
    */
    virtual void endMarginalAncestralState(bool orig_kernel_nonrev, double* &ptn_ancestral_prob, int* &ptn_ancestral_seq);

    /**
     	 compute the joint ancestral states at a pattern (Pupko et al. 2000)
     */
    void computeJointAncestralSequences(int *ancestral_seqs);

    /**
     * compute max ancestral likelihood according to
     *  step 1-3 of the dynamic programming algorithm of Pupko et al. 2000, MBE 17:890-896
     *  @param dad_branch branch leading to an internal node where to obtain ancestral sequence
     *  @param dad dad of the target internal node
     *  @param[out] C array storing all information about max ancestral states
     */
    void computeAncestralLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad, int *C);

    /**
     * compute max ancestral states according to
     *  step 4-5 of the dynamic programming algorithm of Pupko et al. 2000, MBE 17:890-896
     *  @param dad_branch branch leading to an internal node where to obtain ancestral sequence
     *  @param dad dad of the target internal node
     *  @param C array storing all information about max ancestral states
     *  @param[out] ancestral_seqs array of size nptn*nnode for ancestral sequences at all internal nodes
     */
    void computeAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad, int *C, int *ancestral_seqs);

    /**
            compute pattern likelihoods only if the accumulated scaling factor is non-zero.
            Otherwise, copy the pattern_lh attribute
            @param pattern_lh (OUT) pattern log-likelihoods,
                            assuming pattern_lh has the size of the number of patterns
            @param cur_logl current log-likelihood (for sanity check)
            @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
     */
    virtual void computePatternLikelihood(double *pattern_lh, double *cur_logl = NULL,
    		double *pattern_lh_cat = NULL, SiteLoglType wsl = WSL_RATECAT);

    /**
            compute pattern posterior probabilities per rate/mixture category
            @param pattern_prob_cat (OUT) all pattern-probabilities per category
            @param wsl either WSL_RATECAT, WSL_MIXTURE or WSL_MIXTURE_RATECAT
     */
    virtual void computePatternProbabilityCategory(double *pattern_prob_cat, SiteLoglType wsl);

    vector<uint64_t> ptn_cat_mask;

    /**
        compute categories for each pattern, update ptn_cat_mask
        @return max number of categories necessary
    */
    virtual int computePatternCategories(IntVector *pattern_ncat = NULL);

    /**
            Compute the variance in tree log-likelihood
            (Kishino & Hasegawa 1989, JME 29:170-179)
            @param pattern_lh pattern log-likelihoods, will be computed if NULL
            @param tree_lh tree log-likelihood, will be computed if ZERO
     */
    double computeLogLVariance(double *pattern_lh = NULL, double tree_lh = 0.0);

    /**
            Compute the variance in log-likelihood difference
            between the current tree and another tree.
            (Kishino & Hasegawa 1989, JME 29:170-179)
            @param pattern_lh_other pattern log-likelihoods of the other tree
            @param pattern_lh pattern log-likelihoods of current tree, will be computed if NULL
     */
    double computeLogLDiffVariance(double *pattern_lh_other, double *pattern_lh = NULL);

    /**
     *  \brief Estimate the observed branch length between \a dad_branch and \a dad analytically.
     *	The ancestral states of the 2 nodes are first computed (Yang, 2006).
     *	Branch length are then computed using analytical formula.
     *	@param[in] dad_branch must be an internal node
     *	@param[in] dad must be an internal node
     *	@return estimated branch length or -1.0 if one of the 2 nodes is leaf
     */
    double computeBayesianBranchLength(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
     * \brief Approximate the branch legnth between \a dad_branch and \a dad using Least Square instead
     * of Newton Raphson
     * @param[in] dad_branch
     * @param[in] dad
     * @return approximated branch length
     */
    double computeLeastSquareBranLen(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
     * Update all subtree distances that are affect by doing an NNI on branch (node1-node2)
     * @param nni NNI move that is carried out
     */
    void updateSubtreeDists(NNIMove &nni);

    /**
     * Compute all pairwise distance of subtree rooted at \a source and other subtrees
     */
    void computeSubtreeDists();

    void getUnmarkedNodes(PhyloNodeVector& unmarkedNodes, PhyloNode* node = NULL, PhyloNode* dad = NULL);

    void computeAllSubtreeDistForOneNode(PhyloNode* source, PhyloNode* nei1, PhyloNode* nei2, PhyloNode* node, PhyloNode* dad);

    double correctBranchLengthF81(double observedBran, double alpha);

    double computeCorrectedBayesianBranchLength(PhyloNeighbor *dad_branch, PhyloNode *dad);

    /**
            Compute the variance in log-likelihood difference
            between the current tree and another tree.
            (Kishino & Hasegawa 1989, JME 29:170-179)
            @param other_tree the other tree to compare
            @param pattern_lh pattern log-likelihoods of current tree, will be computed if NULL
     */
    double computeLogLDiffVariance(PhyloTree *other_tree, double *pattern_lh = NULL);

    /**
            Roll back the tree saved with only Taxon IDs and branch lengths.
            For this function to work, one must printTree before with WT_TAXON_ID + WT_BR_LEN
            @param best_tree_string input stream to read from
     */
    void rollBack(istream &best_tree_string);

    /**
            refactored 2015-12-22: Taxon IDs instead of Taxon names to save space!
            Read the tree saved with Taxon IDs and branch lengths.
            @param tree_string tree string to read from
            @param updatePLL if true, tree is read into PLL
     */
    virtual void readTreeString(const string &tree_string);

    /**
            Read the tree saved with Taxon names and branch lengths.
            @param tree_string tree string to read from
            @param updatePLL if true, tree is read into PLL
     */
    virtual void readTreeStringSeqName(const string &tree_string);

    /**
            Read the tree saved with Taxon Names and branch lengths.
            @param tree_string tree string to read from
     */
    void readTreeFile(const string &file_name);
    
    /*
            refactored 2015-12-22: Taxon IDs instead of Taxon names to save space!
     * Return the tree string contining taxon IDs and branch lengths
     * @return
     * @param format (WT_TAXON_ID, WT_BR_LEN, ...)
     * @return the tree string with the specified format
     */
    virtual string getTreeString();

    /**
     * Assign branch lengths for branch that has no or negative length
     * With single model branch lengths are assigned using parsimony. With partition model
     * branch lengths are assigned randomly
     * @param force_change if true then force fixing also positive branch lengths
     * @return number of branches fixed
     */
    virtual int wrapperFixNegativeBranch(bool force_change);

    /**
     * Read the newick string into PLL kernel
     * @param newickTree
     */
    void pllReadNewick(string newickTree);

    /**
     *  Return the sorted topology without branch length, used to compare tree topology
     *  @param
     *      printBranchLength true/false
     */
    string getTopologyString(bool printBranchLength);


    bool checkEqualScalingFactor(double &sum_scaling, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /****************************************************************************
            computing derivatives of likelihood function
     ****************************************************************************/

    //template <const int nstates>
//    void computeLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad,
    //                                double &df, double &ddf,
    //                                LikelihoodBufferSet& buffers);

//    void computeSitemodelLikelihoodDervEigen(PhyloNeighbor *dad_branch, PhyloNode *dad,
    //                                         double &df, double &ddf,
    //                                         LikelihoodBufferSet& buffers);

//    template <class VectorClass, const int VCSIZE, const int nstates>
//    void computeLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch,
    //                                    PhyloNode *dad, double &df, double &ddf,
    //                                    LikelihoodBufferSet& buffers);

    void computeNonrevLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                     double *df, double *ddf, LikelihoodBufferSet&);
    template<class VectorClass, const bool SAFE_NUMERIC, const bool FMA = false>
    void computeNonrevLikelihoodDervGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                double *df, double *ddf,
                                                LikelihoodBufferSet& buffers);
    template<class VectorClass, const bool SAFE_NUMERIC,
             const int nstates, const bool FMA = false>
    void computeNonrevLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                         double *df, double *ddf,
                                         LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates,
              const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodBufferSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                     size_t ptn_lower, size_t ptn_upper, int thread_id,
                                     LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC,
            const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodBufferGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                            size_t ptn_lower, size_t ptn_upper, int thread_id,
                                            LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates,
            const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodDervSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                   double *df, double *ddf,
                                   LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC,
              const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodDervGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                          double *df, double *ddf,
                                          LikelihoodBufferSet& buffers);

    /** For Mixlen stuffs */
    virtual int getCurMixture() { return 0; }

    template <class VectorClass, const bool SAFE_NUMERIC, const int nstates,
              const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodDervMixlenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                         double &df, double &ddf,
                                         LikelihoodBufferSet& buffers);

    template <class VectorClass, const bool SAFE_NUMERIC,
              const bool FMA = false, const bool SITE_MODEL = false>
    void computeLikelihoodDervMixlenGenericSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad,
                                                double &df, double &ddf,
                                                LikelihoodBufferSet& buffers);


    /*
    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeMixrateLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeMixtureLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);

    template <class VectorClass, const int VCSIZE, const int nstates>
    void computeSitemodelLikelihoodDervEigenSIMD(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf);
    */

    /**
            compute tree likelihood and derivatives on a branch. used to optimize branch length
            @param dad_branch the branch leading to the subtree
            @param dad its dad, used to direct the tranversal
            @param df (OUT) first derivative
            @param ddf (OUT) second derivative
            @return tree likelihood
     */
    void computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad,
                               double *df, double *ddf,
                               LikelihoodBufferSet& buffers);

    typedef void (PhyloTree::*ComputeLikelihoodDervType)(PhyloNeighbor *, PhyloNode *, double *, double *,
                                                         LikelihoodBufferSet&);
    ComputeLikelihoodDervType computeLikelihoodDervPointer;

    typedef void (PhyloTree::*ComputeLikelihoodDervMixlenType)(PhyloNeighbor *, PhyloNode *,
                                                               double &, double &,
                                                               LikelihoodBufferSet&);
    ComputeLikelihoodDervMixlenType computeLikelihoodDervMixlenPointer;

    /****************************************************************************
            Stepwise addition (greedy) by maximum parsimony
     ****************************************************************************/

    /** constraint tree used to guide tree search */
    ConstraintTree constraintTree;

    /**
     insert a node to a target branch
     @param added_node node to add
     @param target_node left node of the target branch
     @param target_dad right node of the target branch
     */
    void insertNode2Branch(PhyloNode* added_node, PhyloNode* target_node, PhyloNode* target_dad);
        
    /**
            FAST VERSION: used internally by computeParsimonyTree() to find the best target branch to add into the tree
            @param added_node node to add
            @param target_node (OUT) one end of the best branch found
            @param target_dad (OUT) the other end of the best branch found
            @param target_partial_pars (OUT) copy of the partial_pars corresponding to best branch
            @param node the current node
            @param dad dad of the node, used to direct the search
            @return the parsimony score of the tree
     */
    int addTaxonMPFast(PhyloNode *added_taxon, PhyloNode *added_node, PhyloNode *node, PhyloNode *dad);

    /**
        create a 3-taxon tree and return random taxon order
        @param[out] taxon_order random taxon order
     */
    void create3TaxonTree(IntVector &taxon_order, int *rand_stream);

    /**
     Extract a bifurcating subtree and return randomly removed Neighbor
     If the tree is bifurcating, nothing change
     @param[out] removed_nei vector of removed Neighbor
     @param[out] attached_node vector of node attached to removed Neighbor
     */
    void extractBifurcatingSubTree(PhyloNeighborVec &removed_nei, PhyloNodeVector &attached_node, int *rand_stream);
    

    /**
     * FAST VERSION: compute parsimony tree by step-wise addition
     * @param out_prefix prefix for .parstree file
     * @param alignment input alignment
     * @param rand_stream random stream
     * @return parsimony score
     */
    virtual int computeParsimonyTree(const char *out_prefix, Alignment *alignment, int *rand_stream);
        
    /****************************************************************************
            Branch length optimization by maximum likelihood
     ****************************************************************************/


    /**
        @param brlen_type either BRLEN_OPTIMIZE, BRLEN_FIX or BRLEN_SCALE
        @return the number of free branch parameters 
    */
    int getNBranchParameters(int brlen_type);

    /**
     * IMPORTANT: semantic change: this function does not return score anymore, for efficiency purpose
            optimize one branch length by ML
            @param node1 1st end node of the branch
            @param node2 2nd end node of the branch
            @param clearLH true to clear the partial likelihood, otherwise false
            @param maxNRStep maximum number of Newton-Raphson steps
     */
    virtual void optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true, int maxNRStep = 100);

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
    virtual void optimizeAllBranches(PhyloNode *node, PhyloNode *dad = NULL, int maxNRStep = 100);

    /**
     * optimize all branch lengths at the subtree rooted at node step-by-step.
     * Using Least Squares instead of Newton Raphson.
     * @param node the current node
     * @param dad dad of the node, used to direct the search
     */
    void optimizeAllBranchesLS(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    void computeBestTraversal(NodeVector &nodes, NodeVector &nodes2);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);

    void moveRoot(Node *node1, Node *node2);

    /**
        Optimize root position for rooted tree
        @param root_dist maximum distance to move root
        @param write_info true to write information to cout
        @param logl_epsilon epsilon of log-likelihood to consider as better
    */
    virtual double optimizeRootPosition(int root_dist, bool write_info, double logl_epsilon);

    /**
     Test all root positions for rooted tree
     @param write_info true to write information to cout
     @param logl_epsilon epsilon of log-likelihood to consider as better
     */
    virtual double testRootPosition(bool write_info, double logl_epsilon);

    /**
            inherited from Optimization class, to return to likelihood of the tree
            when the current branceh length is set to value
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
    virtual void computeFuncDerv(double value, double &df, double &ddf);

    /**
        optimize the scaling factor for tree length, given all branch lengths fixed
        @param scaling (IN/OUT) start value of scaling factor, and as output the optimal value
        @param gradient_epsilon gradient epsilon
        @return optimal tree log-likelihood
    */
    double optimizeTreeLengthScaling(double min_scaling, double &scaling, double max_scaling, double gradient_epsilon);

    /**
        print tree length scaling to a file (requested by Rob Lanfear)
        @param filename output file name written in YAML format 
    */
    void printTreeLengthScaling(const char *filename);

    /**
     optimize pattern-specific rates by maxmimum likelihood given the tree with fixed branch lengths
     This function will call optimizeTreeLengthScaling
     */
    void optimizePatternRates(DoubleVector &pattern_rates);

     /****************************************************************************
            Branch length optimization by Least Squares
     ****************************************************************************/

    /**
     * Estimate the current branch using least squares
     * @param node1 first node of the branch
     * @param node2 second node of the branch
     * @return
     */
    double optimizeOneBranchLS(PhyloNode *node1, PhyloNode *node2);

    /****************************************************************************
            Auxilary functions and varialbes for speeding up branch length optimization (RAxML Trick)
     ****************************************************************************/

    LikelihoodBufferSet tree_buffers;

    /**
     * frequencies of alignment patterns, used as buffer for likelihood computation
     */
    double *ptn_freq;
    
    /**
     * frequencies of aln->ordered_pattern, used as buffer for parsimony computation
     */
    UINT *ptn_freq_pars;

    /**
     * used as buffer for faster likelihood computation
     * for const pattern: it stores product of p_invar and state frequency
     * for other pattern: zero
     */
    double *ptn_invar;

    vector<TraversalInfo> traversal_info;


    /****************************************************************************
            Nearest Neighbor Interchange by maximum likelihood
     ****************************************************************************/

    /**
            Deprecated
            search by a nearest neigbor interchange, then optimize branch lengths. Do it
            until tree does not improve
            @return the likelihood of the tree
     */
//    double optimizeNNIBranches();

    /**
            search by a nearest neigbor interchange
            @return the likelihood of the tree
     */
    //double optimizeNNI();

    /**
            search by a nearest neigbor interchange
            @param cur_score current likelihood score
            @param node the current node
            @param dad dad of the node, used to direct the search
            @return the likelihood of the tree
     */
//    double optimizeNNI(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL
//            /*,ostream *out = NULL, int brtype = 0, ostream *out_lh = NULL, ostream *site_lh = NULL,
//    StringIntMap *treels = NULL, vector<double*> *treels_ptnlh = NULL, DoubleVector *treels_logl = NULL,
//    int *max_trees = NULL, double *logl_cutoff = NULL*/
//            );


    /**
       search for the best NNI move corresponding to this branch
       @return NNIMove the best NNI, this NNI could be worse than the current tree
       according to the evaluation scheme in use
       @param node1 1 of the 2 nodes on the branch
       @param node2 1 of the 2 nodes on the branch
       @param nniMoves (IN/OUT) detailed information of the 2 NNIs, set .ptnlh to compute pattern likelihoods
     */
    virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves = nullptr);

    /**
            Do an NNI
            @param move reference to an NNI move object containing information about the move
            @param clearLH decides whether or not the partial likelihood should be cleared
     */
    virtual void doNNI(NNIMove &move, bool clearLH = true);

    /**
     * [DEPRECATED]
     * Randomly choose perform an NNI, out of the two defined by branch node1-node2.
     * This function also clear the corresponding partial likelihood vectors
     *
     * @param branch on which a random NNI is done
     */
//    void doOneRandomNNI(Branch branch);

    /**
    *   Get a random NNI from an internal branch, checking for consistency with constraintTree
    *   @param branch the internal branch
    *   @return an NNIMove, node1 and node2 are set to NULL if not consistent with constraintTree
    */
    NNIMove getRandomNNI(Branch& branch);


    /**
     *   Apply 5 new branch lengths stored in the NNI move
     *   @param nnimove the NNI move currently in consideration
     */
    virtual void changeNNIBrans(NNIMove &nnimove);

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
            @param added_taxon leaf node to add
            @param added_node interior node to add
            @param node the current node
            @param dad dad of the node, used to direct the search
            @param isAddedAtMidpoint true if the node is to be added at the middle of
                    the branch (if false, can be added anywhere in the branch)
            @param target_node (OUT) one end of the best branch found
            @param target_dad (OUT) the other end of the best branch found
            @param len_to_new_taxon length of link between added_taxon and added_node
            @param len_to_node length of link between added_node and target_node
            @param len_to_dad length of link between added_dode and target_dad
            @return the likelihood of the tree
     */
    double addTaxonML(PhyloNode* added_taxon,  PhyloNode *added_node,
                      PhyloNode *node,         PhyloNode *dad,
                      bool isAddedAtMidpoint,
                      PhyloNode* &target_node, PhyloNode* &target_dad,
                      double& lenToNewTaxon, double& lenToNode, double& lenToDad);
    
    /**
            used internally by addTaxonML() to get parsimony branch length estimates
            which are used as initial estimates for calculating ML branch lengths,
            for candidate taxa at possible placement sites.
            @param  fromNode
            @param  toNode
            @return the updated branch length
     */

    double recomputeParsimonyBranchLength(PhyloNode* fromNode, PhyloNode* toNode);

    
    /****************************************************************************
            Distance function
     ****************************************************************************/

    virtual void prepareToComputeDistances();
    
    /**
            compute the distance between 2 sequences.
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @param initial_dist initial distance
            @param (OUT) variance of distance between seq1 and seq2
            @return distance between seq1 and seq2
     */
    
    virtual bool hasMatrixOfConvertedSequences() const;

    virtual size_t getConvertedSequenceLength() const;
    
    virtual const char* getConvertedSequenceByNumber(int seq1) const;
    
    virtual const int* getConvertedSequenceFrequencies() const;
    
    virtual const int* getConvertedSequenceNonConstFrequencies() const;
    
    virtual int  getSumOfFrequenciesForSitesWithConstantState(int state) const;
    
    virtual void doneComputingDistances();
    
    virtual double computeDist(int seq1, int seq2, double initial_dist, double &var);

    virtual double computeDist(int seq1, int seq2, double initial_dist);

    void decideDistanceFilePath(Params& params);
    
    void printDistanceFile();
    
    /**
            compute distance and variance matrix, assume dist_mat and var_mat are allocated by memory of size num_seqs * num_seqs.
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @param var_mat (OUT) variance matrix for distance matrix
            @return the longest distance
     */
    double computeDist(double *dist_mat, double *var_mat);

    double computeDist_Experimental(double *dist_mat, double *var_mat);
    
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
            @return the longest distance
     */
    double computeDist(Params &params, Alignment *alignment, double* &dist_mat, double* &var_mat);

    /**
            compute observed distance matrix, allocating memory if necessary
            @param params program parameters
            @param alignment input alignment
            @param dist_mat (OUT) distance matrix between all pairs of sequences in the alignment
            @return the longest distance
     */
    double computeObsDist(Params &params, Alignment *alignment, double* &dist_mat);

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
     */
    void computeBioNJ(Params &params);

    /**
        called by fixNegativeBranch to fix one branch
        @param branch_length new branch length
        @param dad_branch dad branch
        @param dad dad node
    */
    virtual void fixOneNegativeBranch(double branch_length, Neighbor *dad_branch, Node *dad);

    /**
            Neighbor-joining/parsimony tree might contain negative branch length. This
            function will fix this.
            @param fixed_length fixed branch length to set to negative branch lengths
            @param node the current node
            @param dad dad of the node, used to direct the search
            @return The number of branches that have no/negative length
     */
    virtual int fixNegativeBranch(bool force = false, PhyloNode *node = nullptr, PhyloNode *dad = nullptr);

    /**
     Jukes-Cantor correction with alpha shape of Gamma distribution
     @param dist observed distance
     @param alpha shape of Gamma distribution
     @return corrected JC distance
     */
    double JukesCantorCorrection(double dist, double alpha);
                                            
    /**
     set all branch lengths using parsimony
     */
    int setParsimonyBranchLengths();

    /**
        set negative branch to a new len
    */
    int setNegativeBranch(bool force, double newlen, Node *node = NULL, Node *dad = NULL);

    // OBSOLETE: assignRandomBranchLengths no longer needed, use fixNegativeBranch instead!
//    int assignRandomBranchLengths(bool force = false, Node *node = NULL, Node *dad = NULL);

    /* compute Bayesian branch lengths based on ancestral sequence reconstruction */
    void computeAllBayesianBranchLengths(PhyloNode *node = nullptr, PhyloNode *dad = nullptr);

    /**
        generate random tree
    */
    void generateRandomTree(TreeGenType tree_type);


    /**
        test the best number of threads
    */
    virtual int testNumThreads();
    
    /**
        ensure that the number of threads is set (perhaps calling testNumThreads to do so),
        perhaps calling warnNumThreads if it might be set too high.
     @param params (parameters supplied on the command-line)
     @param suppressAnyThreadCountWarnings (self explanatory)
    */
    int ensureNumberOfThreadsIsSet(Params *params,
                                   bool suppressAnyThreadCountWarnings);

    /**
        print warning about too many threads for short alignments
    */
    void warnNumThreads() const;

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
     *  original implementation by Minh
     */
    double optimizeSPR_old(double cur_score, PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     *  original implementation by Minh
     */
    double swapSPR_old(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1,
            PhyloNode *orig_node1, PhyloNode *orig_node2,
            PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor*> &spr_path);

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
    void resampleLh(double **pat_lh, double *lh_new, int *rstream);

    /**
            Test one branch of the tree with aLRT SH-like interpretation
     */
    double testOneBranch(double best_score, double *pattern_lh, 
            int reps, int lbp_reps,
            PhyloNode *node1, PhyloNode *node2, 
            double &lbp_support, double &aLRT_support, double &aBayes_support);

    /**
            Test all branches of the tree with aLRT SH-like interpretation
     */
    virtual int testAllBranches(int threshold, double best_score, double *pattern_lh,
            int reps, int lbp_reps, bool aLRT_test, bool aBayes_test,
            PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /****************************************************************************
            Quartet functions
     ****************************************************************************/

    QuartetGroups LMGroups;
    /**
     * for doLikelihoodMapping reportLikelihoodMapping: likelihood mapping information by region
     */
    vector<QuartetInfo> lmap_quartet_info;
    int areacount[8];
    int cornercount[4];
    // int areacount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // int cornercount[4] = {0, 0, 0, 0};

    /**
     * for doLikelihoodMapping, reportLikelihoodMapping: likelihood mapping information by sequence
     */
    vector<SeqQuartetInfo> lmap_seq_quartet_info;

    /** generate a bunch of quartets and compute likelihood for 3 quartet trees for each replicate
        @param lmap_num_quartets number of quartets
        @param lmap_quartet_info (OUT) vector of quartet information
    */
    void computeQuartetLikelihoods(vector<QuartetInfo> &lmap_quartet_info, QuartetGroups &LMGroups);

    /** main function that performs likelihood mapping analysis (Strimmer & von Haeseler 1997) */
    void doLikelihoodMapping();

    /** output results of likelihood mapping analysis */
    void reportLikelihoodMapping(ofstream &out);

    /** read clusters for likelihood mapping analysis */
    void readLikelihoodMappingGroups(char *filename, QuartetGroups &LMGroups);

    /**
     compute site concordance factor and assign node names
     */
    void computeSiteConcordance(map<string,string> &meanings);

    /**
     compute site concordance factor
     @param branch target branch
     @param nquartets number of quartets
     @param[out] info concordance information
     @param rstream random stream
     */
    virtual void computeSiteConcordance(Branch &branch, int nquartets, int *rstream);

    /**
     Compute gene concordance factor
     for each branch, assign how many times this branch appears in the input set of trees.
     Work fine also when the trees do not have the same taxon set.
     */
    void computeGeneConcordance(MTreeSet &trees, map<string,string> &meanings);

    /**
     Compute quartet concordance factor and internode certainty, similar to Zhou et al. biorxiv
     Work fine also when the trees do not have the same taxon set.
     */
    void computeQuartetConcordance(MTreeSet &trees);

    /**
     compute quartet concordance factor and assign node names
     @param branch target branch
     @param trees input tree set
     @return quartet concordance factor
     */
    double computeQuartetConcordance(Branch &branch, MTreeSet &trees);

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

    /**
            randomize the neighbor orders of all nodes
     */
    void randomizeNeighbors(Node *node = NULL, Node *dad = NULL);

    virtual void changeLikelihoodKernel(LikelihoodKernel lk);

    virtual void setLikelihoodKernel(LikelihoodKernel lk);

    virtual void setNumThreads(int num_threads);

#if defined(BINARY32) || defined(__NOAVX__)
    void setLikelihoodKernelAVX() {}
    void setLikelihoodKernelFMA() {}
#else
    void setLikelihoodKernelAVX();
    void setLikelihoodKernelFMA();
    void setLikelihoodKernelAVX512();
#endif
    virtual void setLikelihoodKernelSSE();
    
    /****************************************************************************
            Public variables
     ****************************************************************************/

    /**
            associated alignment
     */
    Alignment *aln;

    /**
     * Distance matrix
     */
    double *dist_matrix;

    /**
     * Variance matrix
     */
    double *var_matrix;

    /** distance matrix file */
    string dist_file;
    
    
    /** becomes true, if and when the distance matrix has been read from a file */
    bool is_dist_file_read;
    
    /**
            TRUE if you want to optimize branch lengths by Newton-Raphson method
     */
    bool optimize_by_newton;

    /**
     *      TRUE if the loglikelihood is computed using SSE
     */
    LikelihoodKernel sse;

    /**
     * for UpperBounds: Initial tree log-likelihood
     */
    double mlInitial;

    /**
     * for UpperBounds: Log-likelihood after optimization of model parameters in the beginning of tree search
     */
    double mlFirstOpt;

    /**
    * for Upper Bounds: how many NNIs have UB < L curScore, that is NNIs for which we don't need to compute likelihood
    */
	int skippedNNIub;

	/**
	* for Upper Bounds: how many NNIs were considered in total
	*/
	int totalNNIub;

    /**
     * for Upper Bounds: min, mean and max UB encountered during the tree search, such that UB < L curScore
     */

    //double minUB, meanUB, maxUB;

    /*
     * for UpperBounds: mlCheck = 1, if previous two values were already saved.
     * Needed, because parameter optimization is done twice before and after tree search
     */

    int mlCheck;

    /*
     * for Upper Bounds: min base frequency
     */

	double minStateFreq;

    /** sequence names that were removed */
	StrVector removed_seqs;

	/** sequence that are identical to one of the removed sequences */
	StrVector twin_seqs;

	size_t num_partial_lh_computations;

	/** remove identical sequences from the tree */
    virtual void removeIdenticalSeqs(Params &params);

    /** reinsert identical sequences into the tree and reset original alignment */
    virtual void reinsertIdenticalSeqs(Alignment *orig_aln);


    /**
            assign the leaf names with the alignment sequence names, using the leaf ID for assignment.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void assignLeafNames(Node *node = NULL, Node *dad = NULL);

    /**
     * initialize partition information for super tree
     */
    virtual void initPartitionInfo() {
    }

    /**
     * print transition matrix for all branches
     *
     */
    void printTransMatrices(Node *node = NULL, Node *dad = NULL);

    /**
     * compute the memory size required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired(size_t ncategory = 1, bool full_mem = false);

    /**
     * compute the memory size for top partitions required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequiredThreaded(size_t ncategory = 1, bool full_mem = false);
    
    void getMemoryRequired(uint64_t &partial_lh_entries, uint64_t &scale_num_entries, uint64_t &partial_pars_entries);

    /****** following variables are for ultra-fast bootstrap *******/
    /** 2 to save all trees, 1 to save intermediate trees */
    int save_all_trees;

    set<int> computeNodeBranchDists(Node *node = NULL, Node *dad = NULL);

    /*
     * Manuel's approach for analytic approximation of branch length given initial guess
        b0: initial guess for the maximum
        @return approximted branch length
    */
    double approxOneBranch(PhyloNode *node, PhyloNode *dad, double b0);

    void approxAllBranches(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    double getCurScore() {
		return curScore;
	}

	void setCurScore(double curScore) {
		this->curScore = curScore;
	}

	/**
	 * This will invalidate curScore variable, used whenever reading a tree!
	 */
	void resetCurScore(double score = 0.0) {
        if (score != 0.0)
            curScore = score;
        else
		    curScore = -DBL_MAX;
        if (model) {
            initializeAllPartialLh();
        }
	}
    
    /**
     * LikelihoodCostCalculator needs to see this (and I don't want to declare it a friend - James).
     */
    inline double getGammaShape() const {
        return site_rate ? site_rate->getGammaShape() : 1.0;
    }

    void computeSeqIdentityAlongTree(Split &resp, Node *node = NULL, Node *dad = NULL);
    void computeSeqIdentityAlongTree();

    double *getPatternLhCatPointer() { return tree_buffers._pattern_lh_cat; }
    
    /**
     * for rooted tree update direction for all branches
     */
    void computeBranchDirection(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
     * clear branch direction for all branches
     */
    void clearBranchDirection(PhyloNode *node = NULL, PhyloNode *dad = NULL);

    /**
        convert from unrooted to rooted tree
    */
    void convertToRooted();

    /**
        convert from rooted to unrooted tree
    */
    void convertToUnrooted();


	/**
		write site-rates to a file in the following format:
		1  rate_1
		2  rate_2
		....
		This function will call computePatternRates()
		@param out output stream to write rates
        @param bayes TRUE to use empirical Bayesian, false for ML method
	*/
	virtual void writeSiteRates(ostream &out, bool bayes, int partid = -1);

    /**
        write site log likelihood to a output stream
        @param out output stream
        @param wsl write site-loglikelihood type
        @param partid partition ID as first column of the line. -1 to omit it
    */
    virtual void writeSiteLh(ostream &out, SiteLoglType wsl, int partid = -1);

    /**
        write branches into a csv file
        Feature requested by Rob Lanfear
        @param out output stream
     */
    virtual void writeBranches(ostream &out);
    
    /**
       @param out the file path to which a (newly generated) distance file has been written
        (or blank, if it hasn't)
     */
    const string& getDistanceFileWritten() const;

    /** disable progress reporting for this tree
     */
    void showNoProgress();

    /** start reporting progress on a task
     @param size how big the task is
     @param name the name of the task (e.g. "evaluating candidate trees")
     @param verb the (past-tense) verb used to describe progress (e.g. "evaluated")
     @param noun the noun used (e.g. "candidate tree")
     @param isUpperBound true if (size) was only an upper bound on how big the task might be
     */
    virtual void initProgress(double size, std::string name, const char* verb, const char* noun, bool isUpperBound = false);

    /** track progress made on a task*/
    virtual void trackProgress(double amount);
    
    /** hide the progress made on a task (e.g. before writing to cout)*/
    virtual void hideProgress();
    
    /** hide the progress made on a task (e.g. after writing to cout)*/
    virtual void showProgress();
    
    /** report that a task is complete*/
    virtual void doneProgress();
    
protected:
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
     *  is the subtree distance matrix need to be computed or updated
     */
    bool subTreeDistComputed;

    /**
     * Map data structure to store distance Candidate trees between subtree.
     * The key is a string which is constructed by concatenating IDs of
     * the 2 nodes, e.g. 15-16
     */
    StringDoubleMap subTreeDists;

    StringDoubleMap subTreeWeights;

    /** distance (# of branches) between 2 nodes */
    int *nodeBranchDists;

    /**
     * A list containing all the marked list. This is used in the dynamic programming
     * algorithm for compute inter subtree distances
     */
    IntPhyloNodeMap markedNodeList;

    /** converted root state, for Tina's zoombie domain */
    char root_state;

    /**
            internal pattern likelihoods per category per state
            will be computed if not NULL and using non-reversible kernel 
    */
    double *_pattern_lh_cat_state;

    /**
            associated substitution model
     */
    ModelSubst *model;

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
            and by computePatternLikelihood() to compute all pattern likelihoods
     */
    PhyloNeighbor *current_it;
    /**
            current branch iterator of the other end, used by computeFunction() to optimize branch lengths
            and by computePatternLikelihood() to compute all pattern likelihoods
     */
    PhyloNeighbor *current_it_back;

    bool is_opt_scaling;

    /** current scaling factor for optimizeTreeLengthScaling() */
    double current_scaling;

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
    double *nni_partial_lh; // used for NNI functions

    /**
            the main memory storing all scaling event numbers for all neighbors of the tree.
            The variable scale_num in PhyloNeighbor will be assigned to a region inside this variable.
     */
    UBYTE *central_scale_num;
    /**
            The total size (in bytes) of the memory block pointed to by central_scale_num
     */
    size_t central_scale_num_size_in_bytes;
    UBYTE *nni_scale_num; // used for NNI functions

    /**
            the main memory storing all partial parsimony states for all neighbors of the tree.
            The variable partial_pars in PhyloNeighbor will be assigned to a region inside this variable.
     */
    UINT *central_partial_pars;

    virtual void reorientPartialLh(PhyloNeighbor* dad_branch, PhyloNode *dad);
    
    //----------- memory saving technique ------//

    /** maximum number of partial_lh_slots */
    int64_t max_lh_slots;

    /** mapping from */
    MemSlotVector mem_slots;

    /**
            TRUE to discard saturated for Meyer & von Haeseler (2003) model
     */
    bool discard_saturated_site;

    /**
     * Temporary partial likelihood array: used when swapping branch and recalculate the
     * likelihood --> avoid calling malloc everytime
     */
//    double *tmp_partial_lh1;
//    double *tmp_partial_lh2;

    /**
     *  Temporary array containing anscentral states.
     *  Used to avoid calling malloc
     */

//    double *tmp_anscentral_state_prob1;
//    double *tmp_anscentral_state_prob2;
    /** pattern-specific rates */
    //double *tmp_ptn_rates;

    /**
     * Temporary scale num array: used when swapping branch and recalculate the
     * likelihood --> avoid calling malloc
     */
//    UBYTE *tmp_scale_num1;
//    UBYTE *tmp_scale_num2;

    /****************************************************************************
            Vector of bit blocks, used for parsimony function
     ****************************************************************************/

    /**
            @return size of the bits block vector for one node
     */
    size_t getBitsBlockSize();

    /**
            allocate new memory for a bit block vector
            @return the allocated memory
     */
    UINT *newBitsBlock();

    virtual void saveCurrentTree(double logl) {
    } // save current tree


    /**
     * Current score of the tree;
     */
    double curScore;

    /** current best parsimony score */
    UINT best_pars_score;

    /** cost_matrix for non-uniform parsimony */
    unsigned int * cost_matrix; // Sep 2016: store cost matrix in 1D array

    /** stateful AlignmentPairwise instances used for distance processing*/
    std::vector<AlignmentPairwise*> distanceProcessors;

    /** Summary information (especially sequence-major matrix of pattern states */
    AlignmentSummary* summary;

    /** Indicates if summary was "borrowed" from another PhyloTree instance
        (if summary isn't borrowed, and isn't null, it needs to be deleted in the
        destructor*/
    bool isSummaryBorrowed;
    
    string distanceFileWritten;
        /** Is set if/when a distance file has been written*/

    
    /** stack of tasks in progress (top of stack is innermost task) */
    progress_display* progress;
    int  progressStackDepth;
    bool isShowingProgressDisabled;
    
    /** becomes true if/when user is warned about the threadcount in use for this tree */
    bool warnedAboutThreadCount;

    /** becomes true if/when user is warned about numerical underflow, during
        processing for this tree */
    bool warnedAboutNumericalUnderflow;

};
        
#endif
