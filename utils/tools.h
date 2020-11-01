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

#ifndef TOOLS_H
#define TOOLS_H

#include <iqtree_config.h>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <stdint.h>
#include <string.h>
#include <sstream>

//#include <sys/time.h>
//#include <time.h>
#include <sys/stat.h>
//#include <math.h>
#include "ncl/ncl.h"
#include "nclextra/msetsblock.h"

#define SPRNG
#include "sprng/sprng.h"

// redefine assertion
inline void _my_assert(const char* expression, const char *func, const char* file, int line)
{
    char *sfile = (char*)strrchr(file, '/');
    if (!sfile) sfile = (char*)file; else sfile++;
    cerr << sfile << ":" << line << ": " << func << ": Assertion `" << expression << "' failed." << endl;
    abort();
}
 
#ifdef NDEBUG
#define ASSERT(EXPRESSION) ((void)0)
#else
    #if defined(__GNUC__) || defined(__clang__)
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : _my_assert(#EXPRESSION, __PRETTY_FUNCTION__, __FILE__, __LINE__))
    #else
        #define ASSERT(EXPRESSION) ((EXPRESSION) ? (void)0 : _my_assert(#EXPRESSION, __func__, __FILE__, __LINE__))
    #endif
#endif


#define USE_HASH_MAP

#if defined(__GNUC__) && !defined(GCC_VERSION)
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
//#else
//#define GCC_VERSION 0
#endif

// for MSVC
#ifndef __func__
#define __func__ __FUNCTION__
#endif

#if defined(USE_HASH_MAP)
//    #include <unordered_map>
//    #include <unordered_set>

	#if defined(_MSC_VER)
		#include <unordered_map>
		#include <unordered_set>
    #elif defined(__clang__)
        // libc++ detected:     _LIBCPP_VERSION
        // libstdc++ detected:  __GLIBCXX__
        #if __has_include(<unordered_map>) // defines _LIBCPP_VERSION
            #include <unordered_map>
            #include <unordered_set>
        #else
            #include <tr1/unordered_map>
            #include <tr1/unordered_set>
            using namespace std::tr1;    
        #endif
	#elif !defined(__GNUC__)
		#include <hash_map>
		#include <hash_set>
		using namespace stdext;
	#elif GCC_VERSION < 40300
		#include <ext/hash_map>
		#include <ext/hash_set>
		using namespace __gnu_cxx;
		#define unordered_map hash_map
		#define unordered_set hash_set
	#else
		#include <tr1/unordered_map>
		#include <tr1/unordered_set>
		using namespace std::tr1;
	#endif

#else
	#include <map>
	#include <set>
#endif

using namespace std;


#if	defined(USE_HASH_MAP) && GCC_VERSION < 40300 && !defined(_MSC_VER) && !defined(__clang__)
/*
        Define the hash function of Split
 */
#if !defined(__GNUC__)
namespace stdext {
#else
namespace __gnu_cxx {
#endif

    template<>
    struct hash<string> {

        size_t operator()(string str) const {
            hash<const char*> hash_str;
            return hash_str(str.c_str());
        }
    };
} // namespace
#endif // USE_HASH_MAP


class Linear {
public:

    Linear(int n, double *x, double *y) {

        // calculate the averages of arrays x and y
        double xa = 0, ya = 0;
        for (int i = 0; i < n; i++) {
            xa += x[i];
            ya += y[i];
        }
        xa /= n;
        ya /= n;

        // calculate auxiliary sums
        double xx = 0, yy = 0, xy = 0;
        for (int i = 0; i < n; i++) {
            double tmpx = x[i] - xa, tmpy = y[i] - ya;
            xx += tmpx * tmpx;
            yy += tmpy * tmpy;
            xy += tmpx * tmpy;
        }

        // calculate regression line parameters

        // make sure slope is not infinite
        ASSERT(fabs(xx) != 0);

        m_b = xy / xx;
        m_a = ya - m_b * xa;
        m_coeff = (fabs(yy) == 0) ? 1 : xy / sqrt(xx * yy);

    }

    double getValue(double x) {
        return m_a + m_b * x;
    }

    //! Returns the slope of the regression line

    double getSlope() {
        return m_b;
    }

    //! Returns the intercept on the Y axis of the regression line

    double getIntercept() {
        return m_a;
    }

    //! Returns the linear regression coefficient

    double getCoefficient() {
        return m_coeff;
    }

private:

    double m_a, m_b, m_coeff;
};

/**
        vector of double number
 */
typedef vector<double> DoubleVector;

/**
        vector of int
 */
typedef vector<int> IntList;


/**
        vector of int
 */
typedef vector<int> IntVector;

/**
        vector of bool
 */


template <class T, class S> class CastingVector: public S {
    //
    //A subclass of vector (or other container) class S,
    //where the entries in S are treated as though
    //they are of type T (via implicit casting).
    //
public:
    typedef S super;
    typedef typename S::size_type size_type;
    
    CastingVector(): super() {}
    explicit CastingVector(size_type initialSize): super(initialSize) {}
    CastingVector(size_type initialSize, const T initialValue)
        : super(initialSize, initialValue) {}
    CastingVector(const CastingVector& rhs): super(rhs) {}
    
    class const_iterator: public super::const_iterator {
        public:
            const_iterator( typename super::const_iterator i) : super::const_iterator(i) {};
            inline T operator*() { return super::const_iterator::operator*(); }
    };

    class iterator : public super::iterator {
        public:
            iterator( typename super::iterator i) : super::iterator(i) {};
            inline T operator*() { return super::iterator::operator*(); }
    };
    
    const_iterator begin() const { return const_iterator ( super::begin() ); }
    iterator       begin()       { return iterator( super::begin() ); }
    const_iterator end()   const { return const_iterator ( super::end() ); }
    iterator       end()         { return iterator( super::end() ); }
    inline T operator[] (typename super::size_type i) const {
        return super::operator[] (i);
    }
    
    class reference {
        private:
            S& to_vector;
            size_type at_index;
        public:
            reference(S& vector, size_type index):
                to_vector(vector), at_index(index) {}
            operator T() { return to_vector[at_index]; }
            reference& operator= (const T new_value) {
                to_vector[at_index] = new_value;
                return *this;
            }
            reference& operator= (const reference& elsewhere) {
                to_vector[at_index] = elsewhere.to_vector[elsewhere.at_index];
                return *this;
            }
    };
    inline reference operator[] (typename super::size_type i) {
        return reference(*(dynamic_cast<S*>(this)), i);
    }
};

typedef CastingVector<bool, std::vector<char>> BoolVector;

/**
        vector of char
 */
typedef vector<char> CharVector;

/**
        vector of string
 */
typedef vector<string> StrVector;


/**
        matrix of double number
 */
#define mmatrix(T) vector< vector<T> >

/**
        matrix of double
 */

typedef mmatrix(double) DoubleMatrix;

typedef unsigned int UINT;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        run mode of program
 */
enum RunMode {
    DETECTED, GREEDY, PRUNING, BOTH_ALG, EXHAUSTIVE, DYNAMIC_PROGRAMMING,
    CALC_DIST, PD_USER_SET, PRINT_TAXA, PRINT_AREA, SCALE_BRANCH_LEN,
    SCALE_NODE_NAME, PD_DISTRIBUTION, LINEAR_PROGRAMMING, STATS //, GBO, MPRO
}; //STATS and GBO added by MA (STATS for some statistics on tree, GBO = guided 'bootstrap'

/**
        type of generating trees or splits graphs
 */
enum TreeGenType {
    NONE, YULE_HARDING, UNIFORM, CATERPILLAR, BALANCED,
    CIRCULAR_SPLIT_GRAPH, TAXA_SET, STAR_TREE
};

/**
        when writing tree:
                WT_BR_LEN - output branch length
                WT_BR_CLADE - put branch length into internal node name
                WT_TAXON_ID - output taxon ID
                WT_INT_NODE - for draw tree, draw the internal node ID
                WT_BR_SCALE - for draw tree, draw the branch proportional to its length
                WT_SORT_TAXA - sort the taxa s.t. subtrees with least taxon ID come first
                WT_APPEND    - append the output file
                WT_NEWLINE   - print a newline after
                WT_BR_LEN_FIXED_WIDTH - print branch length in fixed number format
 */
const int WT_BR_LEN = 1;
const int WT_BR_CLADE = 2;
const int WT_TAXON_ID = 4;
const int WT_INT_NODE = 8;
const int WT_BR_SCALE = 16;
const int WT_SORT_TAXA = 32;
const int WT_APPEND = 64;
const int WT_NEWLINE = 128;
const int WT_BR_LEN_FIXED_WIDTH = 256;
const int WT_BR_ID = 512;
const int WT_BR_LEN_ROUNDING = 1024;
const int WT_BR_LEN_SHORT = 2048; // store only 6 digits after the comma for branch lengths
const int WT_BR_ATTR = 4096; // print branch attributes
#ifdef TRUE
#undef TRUE
#endif
const int TRUE = 1;
#ifdef FALSE
#undef FALSE
#endif
const int FALSE = 0;

/**
 *  Specify different ways of doing an NNI.
 *  TOPO_ONLY: only change the tree topology
 *  TOPO_UPDATE_LV: the same as above but the partial likelihoods are update in addition
 *  NNI1: optimize the central branch after changing the tree topology
 *  NNI5: optimized the 5 affected branches after changing the tree topology
 */
enum NNI_Type {
    TOPO_ONLY,
    TOPO_UPDATE_LV,
    NNI1,
    NNI5
};

/**
        when computing Robinson-Foulds distances
 */
const int RF_ADJACENT_PAIR = 1;
const int RF_ALL_PAIR = 2;
const int RF_TWO_TREE_SETS = 3;
const int RF_TWO_TREE_SETS_EXTENDED = 4; // work for trees with non-equal taxon sets

/**
        split weight summarization
 */
const int SW_COUNT = 1; // just counting the number of splits
const int SW_SUM = 2; // take the sum of all split weights
const int SW_AVG_ALL = 3; // take the split weight average over all trees
const int SW_AVG_PRESENT = 4; // take the split weight average over all trees that the split is present

/**
        search mode
 */
//enum SearchMode {EXHAUSTIVE, EXHAUSTIVE_CIRCULAR};

/**
        input type, tree or splits graph
 */
enum InputType {
    IN_NEWICK, IN_NEXUS, IN_FASTA, IN_PHYLIP, IN_COUNTS, IN_CLUSTAL, IN_MSF, IN_OTHER
};

  // TODO DS: SAMPLING_SAMPLED is DEPRECATED and it is not possible to run PoMo with SAMPLING_SAMPLED.
enum SamplingType {
  SAMPLING_WEIGHTED_BINOM, SAMPLING_WEIGHTED_HYPER, SAMPLING_SAMPLED
};

/**
        verbose mode, determine how verbose should the screen be printed.
 */
enum VerboseMode {
    VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG
};

/**
        verbose level on the screen
 */
extern VerboseMode verbose_mode;

/**
        consensus reconstruction type
 */
enum ConsensusType {
    CT_NONE, CT_CONSENSUS_TREE, CT_CONSENSUS_NETWORK,
    CT_ASSIGN_SUPPORT, CT_ASSIGN_SUPPORT_EXTENDED, COMPARE
};

enum TestType {
    TEST_NONE, TEST_COMPATIBLE, TEST_CIRCULAR, TEST_WEAKLY_COMPATIBLE, TEST_K_COMPATIBLE
};

/**
        State frequency type
 */
enum StateFreqType {
    FREQ_UNKNOWN, FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, FREQ_ESTIMATE,
    FREQ_CODON_1x4, FREQ_CODON_3x4, FREQ_CODON_3x4C, // special frequency for codon model
    FREQ_MIXTURE, // mixture-frequency model
    // FREQ_DNA_RY has pi_A+pi_G = 0.5 = pi_C+pi_T. Similarly WS pairs (AT)(CG),
    // MK pairs (AC)(GT) in same way.
    FREQ_DNA_RY, FREQ_DNA_WS, FREQ_DNA_MK,
    // in following, digits indicate which frequencies must equal each other
    // (in ACGT order), e.g. 2131 means pi_C=pi_T (pi_A, pi_G unconstrained)
    FREQ_DNA_1112, FREQ_DNA_1121, FREQ_DNA_1211, FREQ_DNA_2111,
    FREQ_DNA_1122, FREQ_DNA_1212, FREQ_DNA_1221, 
    FREQ_DNA_1123, FREQ_DNA_1213, FREQ_DNA_1231, 
    FREQ_DNA_2113, FREQ_DNA_2131, FREQ_DNA_2311, 
};

/*
    outfile file format
 FORMAT_NORMAL: usual file format used so far
 FORMAT_CSV: csv file format
 FORMAT_TSV: tab separated file format
 */
enum FileFormat {
    FORMAT_NORMAL, FORMAT_CSV, FORMAT_TSV
};

enum ModelTestCriterion {
    MTC_AIC, MTC_AICC, MTC_BIC, MTC_ALL
};

/**
 PartitionFinder merging algorithm
 */
enum PartitionMerge {
    MERGE_NONE, MERGE_GREEDY, MERGE_RCLUSTER, MERGE_RCLUSTERF, MERGE_KMEANS
};

/**
        Stopping condition type
 */
enum STOP_CONDITION {
    SC_FIXED_ITERATION, SC_WEIBULL, SC_UNSUCCESS_ITERATION, SC_BOOTSTRAP_CORRELATION, SC_REAL_TIME
};

enum IQP_ASSESS_QUARTET {
    IQP_DISTANCE, IQP_PARSIMONY, IQP_BOOTSTRAP
};

enum LEAST_SQUARE_VAR {
    OLS, WLS_FIRST_TAYLOR, WLS_FITCH_MARGOLIASH, WLS_SECOND_TAYLOR, WLS_PAUPLIN
};

enum START_TREE_TYPE {
	STT_BIONJ, STT_PARSIMONY, STT_PLL_PARSIMONY, STT_RANDOM_TREE, STT_USER_TREE,
    STT_PARSIMONY_JOINING,
};

const int MCAT_LOG = 1; // categorize by log(rate) for Meyer & von Haeseler model
const int MCAT_MEAN = 2; // take the mean of rates for each category for Meyer & von Haeseler model
const int MCAT_PATTERN = 4; // categorize site-patterns instead of sites for Meyer & von Haeseler model

/* TODO DS: For PoMo, this setting does not make sense.  At the
   moment, when using PoMo, MAX_GENETIC_DIST is amended, wherever it
   is used. */
const double MAX_GENETIC_DIST = 9.0;

struct NNIInfo {
    double lh_score[4]; // tree log-likelihood of zero-branch, current tree, NNI tree 1, NNI tree 2
    double br_len[4]; // length of current branch, optimized branch, NNI branch 1, NNI branch 2
    int nni_round;
    int iqpnni_iteration;
};

/*
    0           = 80386 instruction set
    1  or above = SSE (XMM) supported by CPU (not testing for O.S. support)
    2  or above = SSE2
    3  or above = SSE3
    4  or above = Supplementary SSE3 (SSSE3)
    5  or above = SSE4.1
    6  or above = SSE4.2
    7  or above = AVX supported by CPU and operating system
    8  or above = AVX2
    9  or above = AVX512F
*/
enum LikelihoodKernel {
	LK_386, LK_SSE, LK_SSE2, LK_SSE3, LK_SSSE3, LK_SSE41, LK_SSE42, LK_AVX, LK_AVX_FMA, LK_AVX512
};

enum LhMemSave {
	LM_PER_NODE, LM_MEM_SAVE
};

enum SiteLoglType {
    WSL_NONE, WSL_SITE, WSL_RATECAT, WSL_MIXTURE, WSL_MIXTURE_RATECAT
};

enum SiteFreqType {
    WSF_NONE, WSF_POSTERIOR_MEAN, WSF_POSTERIOR_MAX
};

enum MatrixExpTechnique { 
    MET_SCALING_SQUARING, 
    MET_EIGEN3LIB_DECOMPOSITION,
    MET_EIGEN_DECOMPOSITION, 
    MET_LIE_MARKOV_DECOMPOSITION
};

/** ascertainment bias correction type */
enum ASCType {
    ASC_NONE, // no ASC
    ASC_VARIANT, // Lewis's correction for variant sites
    ASC_VARIANT_MISSING, // Holder's correction for variant sites with missing data
    ASC_INFORMATIVE, // correction for parsimony-informative sites
    ASC_INFORMATIVE_MISSING // Holder's correction for informative sites with missing data
};

enum AncestralSeqType {
    AST_NONE, AST_MARGINAL, AST_JOINT
};

enum SymTest {
    SYMTEST_NONE, SYMTEST_BINOM, SYMTEST_MAXDIV
};

const int BRLEN_OPTIMIZE = 0; // optimize branch lengths
const int BRLEN_FIX      = 1; // fix branch lengths
const int BRLEN_SCALE    = 2; // scale branch lengths
const int TOPO_UNLINKED  = 3; // unlinked/separate tree topologies between partitions

const int OUT_LOG       = 1; // .log file written or not
const int OUT_TREEFILE  = 2; // .treefile file written or not
const int OUT_IQTREE    = 4; // .iqtree file written or not
const int OUT_UNIQUESEQ = 8; // .uniqueseq file written or not


const double MIN_GAMMA_RATE = 1e-6;
// change from 0.01 to 0.02 as 0.01 causes numerical problems
const double MIN_GAMMA_SHAPE = 0.02;
const double MAX_GAMMA_SHAPE = 1000.0;
const double TOL_GAMMA_SHAPE = 0.001;


/** maximum number of newton-raphson steps for NNI branch evaluation */
extern int NNI_MAX_NR_STEP;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        program parameters, everything is specified here
        Use singleton pattern to avoid using global variable or
        having to pass the params variable around
 */
class Params {
public:
    static Params& getInstance();
private:
    Params () {}; // Disable constructor
    // Temoprarily commented out because void PhyloSuperTree::readPartition(Params &params)
    // make a copy of params?
    //Params (Params const&) {}; // Disable copy constructor
    //void operator=(Params const&) {}; // Disable assignment
public:

    /**
    *  Fast and accurate optimiation for alpha and p_invar
    */
    bool fai;

    /**
     *  Option to check memory consumption only
     */
    bool memCheck;

    /**
     *  The support threshold for stable splits (Default = 0.9)
     */
    double stableSplitThreshold;

    /**
     *  Option for adaptive perturbation.
     *  Branches that are shared among all candidate trees will be perturbed
     */
    bool adaptPertubation;

	/**
	 *  Option to do mutlipe start for estimating alpha and p_invar
	 */
	bool testAlpha;

    /**
     *  Restart the optimization of alpha and pinvar from different starting
     *  pinv values (supercedes the option testAlpha
     */
    bool opt_gammai;

    /**
     *  A faster version of opt_gammai using a heuristic similar to binary search.
     *  Thus, one does not need to perform 10 independent trials as in opt_gammai.
     */
    bool opt_gammai_fast;

    bool opt_gammai_keep_bran;

    /**
     *  Automatic adjust the log-likelihood espilon using some heuristic
     */
    bool testAlphaEpsAdaptive;

    /**
     *  Use random starting points for alpha
     */
    bool randomAlpha;

    /**
     *  Logl epsilon to test for initial alpha and pinvar values.
     *  This does not need to be small (default value = 100)
     */
    double testAlphaEps;

    /**
     *  Perform exhaustive search for parameter alpha and p_invar
     */
    bool exh_ai;

	/**
	 *  Text file contain all pairs of alpha and p_invar to
	 *  evaluate.
	 *  TODO Remove this option and implement the exhaustive search
	 *  directly into IQ-TREE
	 */
	char* alpha_invar_file;

	/**
	 *  Enable tabu search for NNI
	 */
	bool tabu;

    /**
	 *  Use (5+5)-ES strategy
	 */
	bool five_plus_five;

	/**
	 * Turn on feature to identify stable splits and fix them during tree search
	 */
	bool fixStableSplits;

	/**
	 *  Number of best trees used to compute stable splits
	 */
	int numSupportTrees;

	/**
	 *  Maximum number of trees stored in the candidate tree set
	 */
	int maxCandidates;

	/**
	 *  Number of starting parsimony trees
	 */
	int numInitTrees;

	/**
	 *  SPR distance (radius) for parsimony tree
	 */
	int sprDist;

    /** cost matrix file for Sankoff parsimony */
    char *sankoff_cost_file;
    
	/**
	 *  Number of NNI locally optimal trees generated from the set of parsimony trees
	 *  Default = 20 (out of 100 parsimony trees)
	 */
	int numNNITrees;

	/**
	 *  Number of best trees in the candidate set used to generate perturbed trees
	 *  In term of evolutionary algorithm, this is the population size
	 */
	int popSize;


	/**
	 *  heuristics for speeding up NNI evaluation
	 */
	bool speednni;


	/**
	 *  portion of NNI used for perturbing the tree
	 */
	double initPS;

	/**
	 *  logl epsilon for model parameter optimization
	 */
	double modelEps;
    
    /**
     logl epsilon for ModelFinder
     */
    double modelfinder_eps;

	/**
	 *  New search heuristics (DEFAULT: ON)
	 */
	bool snni;

	/**
	 *  Specify how the branch lengths are optimzed after each NNI operation
	 *  (No optimization, 1 branch optimization, 5 branch optimization)
	 */
    NNI_Type nni_type;

    /**
     *  Different type of Least Square variances
     */
	LEAST_SQUARE_VAR ls_var_type;

	/**
	 *  Threshold (likelihood difference between NNI and current tree)
	 *  to start optimizing 5 branches
	 */
	double nniThresHold;

	/**
	 *  Optimize 5 branches on NNI tree
	 */
	bool nni5;

	/**
	 *  Number of steps for the loop evaluating 5 branches around NNI 
	 */
	int nni5_num_eval;

	/**
	 *  Number of traversal for all branch lengths optimization of the initial tree 
	 */
	int brlen_num_traversal;

    /**
     *  Number of branch length optimization rounds performed after
     *  each NNI step (DEFAULT: 1)
     */
    int numSmoothTree;

    /**
     *   compute least square branches for a given tree
     */
    bool leastSquareBranch;

    /** TRUE to apply Manuel's analytic approximation formulae for branch length */
    bool manuel_analytic_approx;

    /** TRUE to compute parsimony branch length of final tree */
    bool pars_branch_length;

    /** TRUE to compute bayesian branch length for the final tree */
    bool bayes_branch_length;

    /**
     *  use Least Square to evaluate NNI
     */
    bool leastSquareNNI;

    /**
     *  epsilon value used to compare log-likelihood between trees
     */
    double loglh_epsilon;

    /*
     *  reinsert leaves back to tree using parsimony
     */
    bool reinsert_par;

    /*
     *  Option to evaluate 10 different starting tree and take the best
     */
    bool bestStart;

    /**
     *  Maximum running time of the tree search in minutes
     */
    double maxtime;

    /**
     *  Turn on parsimony branch length estimation
     */
    bool parbran;

    /**
     *  option to turn on phylogenetic library
     */
    bool pll;

    /**
     *  OBSOLETE! Stopping rule for the tree search
     */
//    bool autostop;

    /**
     *  Number of maximum unsuccessful iterations after the search is stopped.
     *  Used for the automatic stopping rule
     */
    int unsuccess_iteration;

    char *binary_aln_file;

    /**
     *  the speed up heuristic will be used after
     *  speedup_iter iteration
     */
    int speedup_iter;

    /**
     *  starting CPU time of the program
     */
    double startCPUTime;

    /** starting real time of the program */
    double start_real_time;

    /**
     *  Number iteration = num_taxa * iteration_multiple
     */
    int iteration_multiple;
    /**
             input file name
     */
    char *user_file;

    /* type of starting tree */
    START_TREE_TYPE start_tree;
    std::string start_tree_subtype_name;

    /** TRUE to infer fast ML tree for ModelFinder */
    bool modelfinder_ml_tree;
    
    /** TRUE to perform final model optimization */
    bool final_model_opt;
    
    /** name of constraint tree file in NEWICK format */
    char *constraint_tree_file;

    /**
            prefix of the output file, default is the same as input file
     */
    char *out_prefix;

    /**
            alignment file name
     */
    char *aln_file;

    /** true if sequential phylip format is used, default: false (interleaved format) */
    bool phylip_sequential_format;

    /**
     SYMTEST_NONE to not perform test of symmetry of Jermiin et al. (default)
     SYMTEST_MAXDIV to perform symmetry test on the pair with maximum divergence
     SYMTEST_BINOM to perform binomial test of all pair p-values
    */
    SymTest symtest;
    
    /** TRUE to do symtest then exist */
    bool symtest_only;
    
    /**
     1 to remove bad loci by SymTest
     2 to remove good loci by SymTest
     */
    int symtest_remove;
    
    /** true to keep zero which may result in many taxon pairs not testable (default: false) */
    bool symtest_keep_zero;
    
    /**
        which test used when removing loci
        0 test of symmetry (default)
        1 test of marginal symmetry
        2 test of internal symmetry
     */
    int symtest_type;
    
    /** pvalue cutoff (default: 0.05) */
    double symtest_pcutoff;

    /** TRUE to print all pairwise statistics */
    double symtest_stat;

    /** Times to shuffle characters within columns of the alignment */
    int symtest_shuffle;

    /**
            file containing multiple trees to evaluate at the end
     */
    string treeset_file;

    /** number of bootstrap replicates for tree topology test */
    int topotest_replicates;

    /** TRUE to optimize model parameters for topology test,
     FALSE (default) to only optimize branch lengths */
    bool topotest_optimize_model;

    /** true to perform weighted SH and KH test */
    bool do_weighted_test;

    /** true to do the approximately unbiased (AU) test */
    bool do_au_test;

    /**
            file specifying partition model
     */
    char *partition_file;

    /**
     *      IMPORTANT (2012-12-21): refactor this variable as below
     * 		defines the relation between edge lengths in supertree and subtrees
     * 		BRLEN_OPTIMIZE (0) for separate edge length (default)
     * 		BRLEN_FIX (1) for joint edge length
     * 		BRLEN_SCALE (2) for proportional edge length
     */
    int partition_type;

    /** PartitionFinder algorithm, default MERGE_NONE */
    PartitionMerge partition_merge;
    
    /** percentage for rcluster algorithm like PartitionFinder */
    double partfinder_rcluster; 

    /** absolute limit on #partition pairs for rcluster algorithm */
    size_t partfinder_rcluster_max;

    /** set of models used for model merging phase */
    string merge_models;

    /** set of rate models used for model merging phase */
    string merge_rates;

    /** use logarithm of rates for clustering algorithm */
    bool partfinder_log_rate;
    
    /** remove all-gap sequences in partition model to account for terrace default: TRUE */
    bool remove_empty_seq;

    /** use terrace aware data structure for partition models, default: TRUE */
    bool terrace_aware;

    /** check if the tree lies on a terrace */
    bool terrace_analysis;

    /**
            B, D, or P for Binary, DNA, or Protein sequences
     */
    char *sequence_type;

    /**
            alignment output file name
     */
    char *aln_output;

    /**
            file containing site likelihood as input for 'guided bootstrap' (added by MA)
     */
    char *siteLL_file;

    /**
            alignment where the gappy patterns will be superimposed into the input alignment
     */
    char *gap_masked_aln;

    /**
            alignment to be concatenated into the input alignment
     */
    char *concatenate_aln;

    /**
            file containing list of sites posititon to keep, format:
            pos1 pos2
            ....
     */
    char *aln_site_list;

    /**
            name of the reference sequence where aln_site_list is based on,
            NULL to take alignment positions.
     */
    char *ref_seq_name;

    /**
            alignment output format
     */
    InputType aln_output_format;
    
    /**
        output file format
     */
    FileFormat output_format;

    /**
     tree in extended newick format with node label like [&label=""]
     */
    bool newick_extended_format;
    
    /**
            TRUE to discard all gappy positions
     */
    bool aln_nogaps;

    /**
     * TRUE to discard all constant sites
     */
    bool aln_no_const_sites;

    /** TRUE to print .alninfo file */
    bool print_aln_info;

    /**
            OBSOLETE compute parsimony score on trees
     */
//    bool parsimony;

    /**
            compute random step-wise addition parsimony tree instead of BIONJ
     */
//    bool parsimony_tree;

    /**
             output file name
     */
    char *out_file;

    /**
             size of the maximal PD-tree
     */
    int sub_size;

    /**
             min size of the maximal PD-tree
             used to calculate all PD-k trees from min_size to sub_size
     */
    int min_size;

    /**
            step_size when running from min_size to sub_size
     */
    int step_size;

    /**
            conservation proprotion, another way of input set size
     */
    double pd_proportion;

    /**
            min conservation proprotion
     */
    double min_proportion;

    /**
            step conservation proprotion
     */
    double step_proportion;

    /**
            sample size for computing PD distribution
     */
    int sample_size;


    /**
            TRUE if want to find all optimal PD-k set
            with the same maximal PD score
     */
    bool find_all;

    /**
             type of random tree to be generated
     */
    TreeGenType tree_gen;

    /**
            when generating random split graph, specify the number of
            splits here!
     */
    int num_splits;

    /**
             running mode: which algorithms to be applied
     */
    RunMode run_mode;

    /**
             real running mode if run_mode == DETECTED
     */
    RunMode detected_mode;

    /**
             parameter file
     */
    char *param_file;

    /**
            file containing taxa names to be included into the PD-tree
     */
    char *initial_file;

    /**
            file containing area names to be included into the PD set
     */
    char *initial_area_file;

    /**
            file containing a list of specific taxa sets which user wants
            to compute PD score on these sets only
     */
    char *pdtaxa_file;

    /**
            sets relation file, in form of a distance matrix file
     */
    char *areas_boundary_file;

    /**
            boundary length modifier
     */
    double boundary_modifier;

    /**
            output file to store the distance matrix
     */
    char *dist_file;

    /**
            true if processing is incremental (previous trees, or distance
            matrices, that were based on a subset of the sequences,
            are to be modified, rather than reconstructed from scratch)
     */
    bool incremental;
    
    /**
     if processing is incremental, a string (upper case) describing how
     new sequences are to be handled (not actually supported as yet)
     */
    std::string incremental_method;
    
    /**
            TRUE to compute the observed distances instead of Juke-Cantor distances, default: FALSE
     */
    bool compute_obs_dist;

    /**
            true to treat unknown characters as different rather than discounting them, default: false.
     */
    bool count_unknown_as_different;
    
    /**
            TRUE to compute the Juke-Cantor distances, default: FALSE
     */
    bool compute_jc_dist;
    
    
    /**
            TRUE to use experimental implementation to calculate observed (or Jukes-Cantor) distances
     */
    bool experimental;

    /**
            TRUE to compute the maximum-likelihood distances
     */
    bool compute_ml_dist;

    /**
            TRUE to compute the maximum-likelihood tree
     */
    bool compute_ml_tree;

    /**
            TRUE to compute *only* the maximum-likelihood tree
      (without using parsimony or observed distances to compute an NJ tree first)
     */
    bool compute_ml_tree_only;
    /**
            file containing budget information
     */
    char *budget_file;

    /**
            used when generating pair of taxa set with overlapping
     */
    int overlap;

    // private use
    /**
             number of times to repeat the algorithms
     */
    int repeated_time;

    /**
             print no tree to output
     */
    int nr_output;

    /**
            input type, tree or splits graph
     */
    InputType intype;

    /**
            total budget, for cost constrained PD problem
     */
    int budget;

    /**
            minimum budget, for cost constrained PD problem
     */
    int min_budget;

    /**
            step_budget when running from min_budget to budget
     */
    int step_budget;

    /**
            name of the root taxon
     */
    const char *root;

    /**
            true if tree is forced to be rooted
     */
    bool is_rooted;

    /**
        maximum distance to move root
     */
    int root_move_dist;

    /**
     TRUE to find best root when optimizing model
     */
    bool root_find;

    /**
     TRUE to test all rooting positions at the end of the run
     */
    bool root_test;

    /**
            min branch length, used to create random tree/network
     */
    double min_len;

    /**
            mean branch length, used to create random tree/network
     */
    double mean_len;

    /**
            max branch length, used to create random tree/network
     */
    double max_len;

    /**
            number of internal branches to set zero length
     */
    int num_zero_len;

    /**
            random number seed
     */
    int ran_seed;

    /**
            run time of the algorithm
     */
    double run_time;

    /**
            limit on the number of optimal PD sets
     */
    int pd_limit;

    /**
            TRUE if one wants to calculate the PD gain matrix in terms of delta_k^j = pd(PD_k \/ {j}) - pd_k
     */
    bool calc_pdgain;

    /**
            TRUE if tree file contains more than 1 tree
     */
    bool multi_tree;

    /**
            2nd user tree used in assignBootstrapSupport
     */
    char *second_tree;

    /** 
        tag each branch with the tree ID where it occurs; "ALL" to tag all branches
    */
    char *support_tag;

    /**
        number of quartets for site concordance factor
     */
    int site_concordance;

    /**
     TRUE to print concordant sites per partition
     */
    bool site_concordance_partition;
    
    /** TRUE to print sCF for all sampled quartet */
    bool print_cf_quartets;

    /** TRUE to print trees associated with discordance factor 1 (NNI-1 tree) */
    bool print_df1_trees;
    
    /** 1 to compute internode certainty */
    int internode_certainty;
    
    /**
            2nd alignment used in computing multinomialProb (Added by MA)
     */
    char *second_align;
    /**
            type of consensus building
     */
    ConsensusType consensus_type;

    /**
            file containing weights for every tree in the input tree file
     */
    char *tree_weight_file;

    /**
            set the TRUE if want to find the minimal PD set, instead of the default maximal PD set
     */
    bool find_pd_min;

    /**
            set TRUE to find area's endemic PD instead of regular PD
     */
    bool endemic_pd;

    /**
            set TRUE to find exclusive PD instead of regular PD
     */
    bool exclusive_pd;

    /**
            to find PD complementarity given this area
     */
    char *complement_area;

    /**
            used for likelihood mapping: for each branch, print the four cluster
     */
    int branch_cluster;

    /**
            file containing taxa order
     */
    char *taxa_order_file;

    /**
            to scale branch length or clade support with a factor
     */
    double scaling_factor;

    /**
            TRUE if always use binary linear programming
     */
    bool binary_programming;

    /**
            test the input split system in one of the TestType
     */
    TestType test_input;

    /**
            burnin value: number of beginning trees to be discarded
     */
    int tree_burnin;

    /**
            maximum number of trees to consider (for e.g. consensus tree construction)
     */
    int tree_max_count;

    /**
        threshold of split frequency, splits appear less than threshold will be discarded
     */
    double split_threshold;

    /**
        thresholds of split frequency with back-slash separator
     */
    char* split_threshold_str;

    /**
            threshold of split weight, splits with weight less than or equal to threshold will be discarded
     */
    double split_weight_threshold;

    /** TRUE to collapse zero branches, default FALSE */
    bool collapse_zero_branch;

    /**
            Way to summarize split weight in the consensus tree or network: SW_SUM, SW_AVG_ALL, or SW_AVG_PRESENT
     */
    double split_weight_summary;

    /**
            TRUE if use quadratic programming (for GUROBI)
     */
    bool quad_programming;

    /**
            true if one wants to optimize tree by subtree pruning and regrafting
     */
    bool tree_spr;

    /**
            true if printing out of optimal sets in NEXUS format
     */
    bool nexus_output;

    /**
            k-representative parameter, used for IQP algorithm
     */
    int k_representative;

    /**
            probability of deleting a leaf, used for IQP algorithm
     */
    double p_delete;

    /**
            min number of iqpnni iterations
     */
    int min_iterations;

    /**
            max number of iqpnni iterations
     */
    int max_iterations;

    /**
            stop condition, SC_FIXED_ITERATION or SC_STOP_PREDICT
     */
    STOP_CONDITION stop_condition;

    /**
            confidence value for stop rule
     */
    double stop_confidence;

    /** number iterations for parameter optimization, default: 100 */
    int num_param_iterations;

    /** number of independent runs (-nrun option) */
    int num_runs;
    
    /**
            name of the substitution model (e.g., HKY, GTR, TN+I+G, JC+G, etc.)
     */
    string model_name;

    /** model name to initialize GTR20 or NONREV protein model */
    char* model_name_init;

    /** number of steps for linked model optimisation, default: 1 */
    int model_opt_steps;
    
    /** set of models for testing */
    string model_set;

    /** set of models to be added into default set */
    char *model_extra_set;

    /** subset of models for testing, e.g. viral, mitochondrial */
    char *model_subset;

    /** set of state frequencies model for testing */
    char *state_freq_set;

    /** set of rate heterogeneity model for testing */
    string ratehet_set;

    /** all models with score worse than the best score + this threshold will be ignored */
    double score_diff_thres;
    
    /** model defition file */
    char *model_def_file;

    /** TRUE to perform ModelOMatic method of Whelan et al. 2015 */
    bool modelomatic;
    
    /** true to redo model testing even if .model file exists */
    bool model_test_again;

    /** 0: use the same tree for model testing 
        1: estimate tree for each model, but initialize the tree for next model 
           by the tree reconstructed from the previous model
        2: estimate tree for each model independently
        */
    short int model_test_and_tree;

    /** true to fist test equal rate model, then test rate heterogeneity (default: false) */
    bool model_test_separate_rate;

    /** TRUE to optimize mixture model weights */
    bool optimize_mixmodel_weight;

    /** number of mixture branch lengths, default 1 */
    int num_mixlen;
    /** TRUE to always optimize rate matrix even if user parameters are specified in e.g. GTR{1,2,3,4,5} */
    bool optimize_rate_matrix;

    /**
            TRUE to store transition matrix into a hash table for computation efficiency
     */
    bool store_trans_matrix;

    /**
            state frequency type
     */
    StateFreqType freq_type;

    /** FALSE to set zero state frequency to 1e-4.
        Default: FALSE (version <= 1.5.5), TRUE (ver >= 1.5.6) */
    bool keep_zero_freq;
    
    /** minimal state frequency for optimisation, default=0.0001 */
    double min_state_freq;


    /**
            the number of rate categories
     */
    int num_rate_cats;

    /**
            maximum number of rate categories
     */
    int min_rate_cats;

    /**
            maximum number of rate categories
     */
    int max_rate_cats;

    /**
            shape parameter (alpha) of the Gamma distribution for site rates
     */
    double gamma_shape;

    /**
            minimum shape parameter (alpha) of the Gamma distribution for site rates
     */
    double min_gamma_shape;

    /**
            TRUE to use median rate for discrete categories, FALSE to use mean rate instead
     */
    bool gamma_median;

    /**
            proportion of invariable sites
     */
    double p_invar_sites;

    /** TRUE to optimize all model and rate parameters jointly by BFGS, default: FALSE */
    bool optimize_model_rate_joint;

    /**
            TRUE if you want to optimize branch lengths by Newton-Raphson method
     */
    bool optimize_by_newton;

    /** optimization algorithm for free rate model: 1-BFGS, 2-BFGS, EM */
    string optimize_alg_freerate;

    /** optimization algorithm for mixture (heterotachy) branch length models */
    string optimize_alg_mixlen;

    /**
     *  Optimization algorithm for +I+G
     */
    string optimize_alg_gammai;

    /**
     * If given model parameters on command line (e.g. -m RY3.4{0.2,-0.4})
     * treat these as fixed model parameters (if false), or treat them as 
     * starting point for optimization search (if true)?
     */

    bool optimize_from_given_params;

    /**
            BRLEN_OPTIMIZE optimize branch lengths during model optimization
            BRLEN_FIX      fix branch lengths during model optimization
            BRLEN_SCALE    scale all branch lengths by the same factor during model optimization
     */
    int fixed_branch_length;

    /** minimum branch length for optimization, default 0.000001 */
    double min_branch_length;

    /** maximum branch length for optimization, default 100 */
    double max_branch_length;


    /**
            criterion to assess important quartet
     */
    IQP_ASSESS_QUARTET iqp_assess_quartet;

    /**
     *      Using IQP algorithm to do tree perturbation
     */
    bool iqp;

    /**
            the LP file is in gurobi format or not
     */
    bool gurobi_format;

    /**
            number of threads for gurobi call
     */
    bool gurobi_threads;

    /**
            TRUE if doing bootstrap on the input trees (good, bad, ugly)
     */
    int num_bootstrap_samples;

    /** bootstrap specification of the form "l1:b1,l2:b2,...,lk:bk"
        to randomly draw b1 sites from the first l1 sites, etc. Note that l1+l2+...+lk
        must equal m, where m is the alignment length. Otherwise, an error will occur.
        The default bootstrap_spec == NULL, a standard procedure is applied, i.e., randomly draw m sites.
    */
    char *bootstrap_spec;

    /** 1 or 2 to perform transfer boostrap expectation (TBE) */
    int transfer_bootstrap;
    
    /** subsampling some number of partitions / sites for analysis */
    int subsampling;

    /** random seed number for subsampling */
    int subsampling_seed;

    /**
            1 if output all intermediate trees (initial trees, NNI-optimal trees and trees after each NNI step)
            2 if output all intermediate trees + 1-NNI-away trees
     */
    int write_intermediate_trees;

    /**
     *  Write all distinct intermediate trees and there likelihoods
     *  Note: intermediate trees are trees that have been visited by the search. These include trees created by
     *  NNI-steps within each NNI iteration.
     */
    bool writeDistImdTrees;

    /**
     *  Write trees obtained at the end of each NNI search
     */
    bool write_candidate_trees;


    /**
        TRUE to avoid duplicated trees while writing intermediate trees
     */
//    bool avoid_duplicated_trees;

    /**
            Robinson-Foulds distance computation mode: RF_ADJACENT PAIR, RF_ALL_PAIR
     */
    int rf_dist_mode;

    /**
     true to compute distance between the same k-th tree in two sets
     */
    bool rf_same_pair;
    
    /**
     true to normalize tree distances, false otherwise
     */
    bool normalize_tree_dist;
    
    /**
            compute the site-specific rates by Meyer & von Haeseler method
     */
    bool mvh_site_rate;

    /**
            FALSE to use MH Model, FALSE for using tree-likelihood
     */
    bool rate_mh_type;

    /**
            TRUE to discard saturated for Meyer & von Haeseler (2003) model
     */
    bool discard_saturated_site;

    /**
            rates will be normalized to this mean value
     */
    double mean_rate;

    /**
            Percentage threshold to accept a branch of the approximate likelihood ratio test
            (aLRT) with SH-like interpretation. See Guindon et al. (2010) Syst. Biol. for details.
            Default: 90%.
     */
    int aLRT_threshold;

    /**
            number of replicates, default: 1000
     */
    int aLRT_replicates;

    /** true to perform aLRT branch test of Anisimova & Gascuel (2006) */
    bool aLRT_test;

    /** true to perform aBayes branch test of Anisimova et al (2011) */
    bool aBayes_test;

    /**
            number of replicates for local bootstrap probabilities method of Adachi & Hasegawa (1996) in MOLPHY
     */
    int localbp_replicates;

    /**
            SSE Option
     */
    LikelihoodKernel SSE;

    /** TRUE for safe numerical scaling (per category; used for large trees), default: FALSE */
    bool lk_safe_scaling;

    /** minimum number of sequences to always use safe scaling, default: 2000 */
    int numseq_safe_scaling;

    /** if true, stubbornly keep processing */
    bool ignore_any_errors;
    
    /** TRUE to force using non-reversible likelihood kernel */
    bool kernel_nonrev;

    /**
     	 	WSL_NONE: do not print anything
            WSL_SITE: print site log-likelihood
            WSL_RATECAT: print site log-likelihood per rate category
            WSL_MIXTURE: print site log-likelihood per mixture class
            WSL_MIXTURE_RATECAT: print site log-likelihood per mixture class per rate category
            WSL_STATE: print site log-likelihood per state
     */
    SiteLoglType print_site_lh;

    /** TRUE to print partition log-likelihood, default: FALSE */
    bool print_partition_lh;

    /**
        control printing posterior probability of each site belonging to a rate/mixture categories
        same meaning as print_site_lh, but results are printed to .siteprob file
        WSL_RATECAT: print site probability per rate category
        WSL_MIXTURE: print site probability per mixture class
        WSL_MIXTURE_RATECAT: print site probability per mixture class per rate category
    */
    SiteLoglType print_site_prob;

    /**
        AST_NONE: do not print ancestral sequences (default)
        AST_MARGINAL: print ancestral sequences by marginal reconstruction
        AST_JOINT: print ancestral sequences by joint reconstruction
    */
    AncestralSeqType print_ancestral_sequence;

    /** minimum probability to assign an ancestral state */
    double min_ancestral_prob;

    /**
        0: print nothing
        1: print site state frequency vectors
    */
    SiteFreqType print_site_state_freq;

    /**
     0 (default): do not print .rate file
     1: print site-specific rates by empirical Bayes
     2: site-specific rates by maximum-likelihood
     */
    int print_site_rate;

    /* 1: print site posterior probability for many trees during tree search */
    int print_trees_site_posterior;

    /**
            TRUE to print tree log-likelihood
     */
    bool print_tree_lh;

    bool print_branch_lengths;

    /****** adaptive NNI search heuristic ******/

    /**
     *  Output log-likelihood
     */
    bool nni_lh;

    /**
     *  The number of iqp iteration before the heuristics is applied
     */
    int speedUpFromIter;

    /**
     *  Lambda in PhyML algorithm
     */
    double lambda;

    /**
     * Confidence level for the speed up heuristics
     */
    double speed_conf;

    bool new_heuristic;

    /***** WH-test (Weiss & von Haeseler 2003) *****/

    /**
            Results of Weiss & Haeseler test of model homogeneity
     */
    double whtest_simulations;
    double whtest_delta;
    double whtest_delta_quantile;
    double whtest_p_value;


    /**
            bit-wise type including MCAT_LOG, MCAT_MEAN
     */
    int mcat_type;

    /**
            initial rate file in format:
            Site Rate
            1  f_1
            2  f_2
            ...
     */
    char *rate_file;

    /***** NGS stuffs   ********/

    /**
            next-generation sequencing input file for Fritz project
     */
    char *ngs_file;

    /**
            next-generation sequencing input file containing mapped reads to the reference genome
     */
    char *ngs_mapped_reads;

    bool ngs_ignore_gaps;

    bool do_pars_multistate;

    /**
            File containing p-values of the genes, for GSS project with Roland
     */
    char *gene_pvalue_file;

    /**
            scaling factor for the p-values
     */
    double gene_scale_factor;

    /**
            transforming pvalues to logarithms
     */
    bool gene_pvalue_loga;

    /***** variables for reading NCBI taxonomy tree *******/

    /**
            NCBI taxonomy ID, for processing nodes.dmp file
     */
    int ncbi_taxid;

    /**
            NCBI taxon rank, restricting the tree to that rank
     */
    const char *ncbi_taxon_level;

    /**
            rank to ingore, e.g., "no rank", branch length to such node will be set to zero
     */
    const char *ncbi_ignore_level;

    /**
            typically names.dmp from NCBI
     */
    const char *ncbi_names_file;

    /**********************************************/
    /******* variables for ECOpd analysis *********/

	/**
		eco_dag_file - contains the food web in matrix form (n species, nxn matrix), 0 for no connection, 1 for predation of j predator on i prey
	*/
	char *eco_dag_file;

    /**
		eco_detail_file - contains IDs of species present in the final set and/or species absent in the TREE or SPLIT system, but present in the food web
	*/
	const char *eco_detail_file;

	/*
	 * the type of the phylo input - tree or network
	 */
	const char *eco_type;

	/*
		k% - percent of species to be conserved
	 */
	int k_percent;

    /*
		diet - percent of species diet to be preserved for species survival
	*/
	int diet_min;
	int diet_max;
	int diet_step;

    /*
		eco_run - run number, used when random branch length is assigned to the edges of an input tree
	*/
	int eco_run;

    /*
		eco_weighted - indicates whether to treat the food web as weighted or not weighted
	*/
	bool eco_weighted;

    /**********************************************/
    /****** variables for upper bound tests *******/
	bool upper_bound;
	bool upper_bound_NNI;
	/*
	 * fraction of current likelihood by which UB will be increased.
	 * if UBincreased < L, ignore corresponding NNI Add a comment to this line
	 */
	double upper_bound_frac;


    /**********************************************/
    /**** variables for ultra-fast bootstrap ******/

    /**
            number of replicates for guided bootstrap
     */
    int gbo_replicates;

	/* interval (l-epsilon,l+epsilon) indicates tie for bootstrap tree
	 * in this case, one tree is picked up at random
	 */
	double ufboot_epsilon;

    /**
            TRUE to check with different max_candidate_trees
     */
    int check_gbo_sample_size;

    /**
            TRUE to use RELL method of Simodaira Hasegawa, FALSE otherwise
     */
    bool use_rell_method;

    /**
            TRUE to use ELW method of Strimmer & Rambaut for new bootstrap, FALSE otherwise
     */
    bool use_elw_method;

    /**
            TRUE to weight each bootstrap sample by its probability, FALSE otherwise
     */
    bool use_weighted_bootstrap;

    /**
            TRUE to use the single ML tree per bootstrap, FALSE to include several sup-optima
     */
    bool use_max_tree_per_bootstrap;

    /** maximum number of candidate trees to consider for new bootstrap */
    int max_candidate_trees;

    /** TRUE if user_file contains topologically distinct trees */
    bool distinct_trees;

    /** NEW: TRUE to update bootstrap trees during the search (do not store treels_ptnlh).
            FALSE to call runGuidedBootstrap() at the end */
    bool online_bootstrap;

    /** minimal correlation coefficient for bootstrap stopping rule */
    double min_correlation;

    /** number of iterations between bootstrap stopping rule check */
    int step_iterations;

    /** TRUE to store all candidate trees in memory */
//    bool store_candidate_trees;

	/** true to print all UFBoot trees to a file */
	int print_ufboot_trees;

    /**********************************************/
    /**** variables for jackknife ******************/

    /** proportion of sites to be dropped in jackknife */
    double jackknife_prop;
    
    /**********************************************/
    /* variables for robust phylogeny (Lanfear & Holland project */
    
    /** proportion of sites to keep in robust phylogeny idea */
    double robust_phy_keep;

    /** use median log-likelihood instead of sum log-likelihood */
    double robust_median;
    
    /****** variables for NNI cutoff heuristics ******/

    /**
            TRUE to empirically estimate nni_cutoff
     */
    bool estimate_nni_cutoff;

    /**
            logl difference with zero-branch tree, to cutoff before evaluating NNI
     */
    double nni_cutoff;

    /**
            sort the NNI before evaluating
     */
    bool nni_sort;

    /**
            Obsolete: TRUE to optimize 5 branches around NNI
     */
    //bool nni_opt_5branches;

    /** print some output info for NNI */
    bool testNNI;

    /** TRUE to do approximate NNIs with approximate branch lengths before a normal NNI */
    bool approximate_nni;


    /** TRUE to compress big file using zlib */
    bool do_compression;

    /**
            number of bootstrap samples for AvH curiosity
     */
    int avh_test;

    /**
            number of bootstrap samples for Arndt's bootstrap plot
     */
    int bootlh_test;

    /**
            partition definition for Arndt's bootstrap plot
     */
    char* bootlh_partitions;

    /** precision when printing out for floating-point number */
    int numeric_precision;

    /** file containing state-frequencies per site for site-specific state frequency model
     * each line has n+1 entries (n=number of states):
     * site_ID state1_freq state2_freq ... staten_freq
     * where site_ID is from 1 to m (m=number of sites)
     */
    char *site_freq_file;

    /**
        user tree file used to estimate site-specific state frequency model 
    */
    char *tree_freq_file;

    /** number of threads for OpenMP version     */
    int num_threads;
    
    /** maximum number of threads, default: #CPU scores  */
    int num_threads_max;
    
    /** true to parallel ModelFinder by models instead of sites */
    bool openmp_by_model;

    /** either MTC_AIC, MTC_AICc, MTC_BIC */
    ModelTestCriterion model_test_criterion;

    /** either MTC_AIC, MTC_AICc, MTC_BIC, or MTC_ALL to stop +R increasing categories */
//    ModelTestCriterion model_test_stop_rule;

    /** sample size for AICc and BIC */
    int model_test_sample_size;

    /** root state, for Tina's zoombie domain */
    char *root_state;

	/**
	 * TRUE to print bootstrap alignments, default: false
	 */
	bool print_bootaln;

    /** TRUE to print bootstrapped site frequency for e.g. PMSF */
    bool print_boot_site_freq;

	/** true to print sub alignments of super alignment, default: false */
	bool print_subaln;

	/** print partition information */
	bool print_partition_info;

	/** TRUE to print concatenated alignment, default: false */
	bool print_conaln;

	/** TRUE to link alpha among Gamma model over partitions */
	bool link_alpha;

    /** TRUE to link substitution models over partitions */
    bool link_model;

    /** name of the joint model across partitions */
    char* model_joint;
    
	/** true to count all distinct trees visited during tree search */
	bool count_trees;

    /// True if PoMo is run; otherwise false.
    bool pomo;

    /// True if sampled input method is used (-st CR..); otherwise false.
    bool pomo_random_sampling;

	/// Virtual population size for PoMo model.
	int pomo_pop_size;

	/* -1 (auto-detect): will be set to 0 if there is enough memory, 1 otherwise
	 * 0: store all partial likelihood vectors
	 * 1: only store 1 partial likelihood vector per node */
	LhMemSave lh_mem_save;
    
    /** true to save buffer, default: false */
    bool buffer_mem_save;

    /** maximum size of memory allowed to use */
    double max_mem_size;

	/* TRUE to print .splits file in star-dot format */
	bool print_splits_file;
    
    /* TRUE to print .splits.nex file in NEXUS format */
    bool print_splits_nex_file;

    
    /** TRUE (default) to ignore identical sequences and add them back at the end */
    bool ignore_identical_seqs;

    /** TRUE to write initial tree to a file (default: false) */
    bool write_init_tree;

    /** TRUE to write branch lengths of partition trees for each branch of supertree */
    bool write_branches;
    
    /** frequencies of const patterns to be inserted into alignment */
    char *freq_const_patterns;
    /** BQM 2015-02-25: true to NOT rescale Gamma+Invar rates by (1-p_invar) */
    bool no_rescale_gamma_invar;

    /** true to compute sequence identity along tree */
    bool compute_seq_identity_along_tree;
    
    /** true to compute sequence composition */
    bool compute_seq_composition;
    
    /** true to ignore checkpoint file */
    bool ignore_checkpoint;
    /** number of quartets for likelihood mapping */
    int64_t lmap_num_quartets;

    /**
            file containing the cluster information for clustered likelihood mapping
     */
    char *lmap_cluster_file;

    /** time (in seconds) between checkpoint dump */
    int checkpoint_dump_interval;
    /** TRUE to print quartet log-likelihoods to .quartetlh file */
    bool print_lmap_quartet_lh;

    /** true if ignoring the "finished" flag in checkpoint file */
    bool force_unfinished;

    /** control output files to be written
     * OUT_LOG
     * OUT_TREEFILE
     * OUT_IQTREE
     */
    int suppress_output_flags;

    /** matrix exponentiation technique for nonreversible models, either 
        MET_SCALING_SQUARING 
        MET_EIGEN_DECOMPOSITION 
        MET_LIE_MARKOV_DECOMPOSITION
    */
    MatrixExpTechnique matrix_exp_technique;

    /**
     * Diep:
     * Data members for UFBoot2-Corr
     */
	bool ufboot2corr; // to turn on the correction mode for UFBoot under model violations, enable by "-bb <nrep> -correct
	bool u2c_nni5; // to use NNI5 during Refinement Step of UFBoot2-Corr

    /** method for phylogenetic dating, currently only LSD is supported */
    string dating_method;

    /** extra commands passed to the dating method */
    string dating_options;

    /** date file that has several lines, each line with a taxon name and date in YYYY-MM-DD */
    string date_file;
    
    /** tip date, a real number or YYYY-MM-DD */
    string date_tip;
    
    /** root date, a real number or YYYY-MM-DD */
    string date_root;
    
    /** false to remove outgroup from the dated tree, default: true */
    bool date_with_outgroup;
    
    /** true to print internal date files for debugging purpose */
    bool date_debug;
    
    /** number of replicates to compute date confidence interval */
    int date_replicates;
    
    /** standard deviation of lognormal relaxed clock model for confidence interval estimate */
    double clock_stddev;

    /** z-score for detecting outlier nodes */
    double date_outlier;

    /** supress the list of sequences */
    double suppress_list_of_sequences;

    /** supress warnings about low or zero distances */
    double suppress_zero_distance_warnings;

    /** supress notes about duplicate sequences */
    double suppress_duplicate_sequence_warnings;
    
    /** format to use when writing (.mldist) distance matrix files */
    string dist_format;

    /** compression level to use when writing (.mldist) distance matrix files */
    int    dist_compression_level;
    
};

/**
        related measures for PD
 */
struct PDRelatedMeasures {
    /**
            names of areas
     */
    vector<string> setName;

    /**
            pd scores of areas
     */
    DoubleVector PDScore;

    /**
            exclusive PD scores
     */
    DoubleVector exclusivePD;

    /**
            endemic pd scores of an area given all other areas
     */
    DoubleVector PDEndemism;

    /**
            pd-complementarity scores of an area given some provided area
     */
    DoubleVector PDComplementarity;

};



/*--------------------------------------------------------------*/
    
inline size_t get_safe_upper_limit(size_t cur_limit) {
    if (Params::getInstance().SSE >= LK_AVX512)
        // AVX-512
        return ((cur_limit+7)/8)*8;
    else
        if (Params::getInstance().SSE >= LK_AVX)
            // AVX
            return ((cur_limit+3)/4)*4;
        else
            // SSE
            return ((cur_limit+1)/2)*2;
}

inline size_t get_safe_upper_limit_float(size_t cur_limit) {
    if (Params::getInstance().SSE >= LK_AVX512)
        // AVX-512
        return ((cur_limit+15)/16)*16;
    else
        if (Params::getInstance().SSE >= LK_AVX)
            // AVX
            return ((cur_limit+7)/8)*8;
        else
            // SSE
            return ((cur_limit+3)/4)*4;
}
/*--------------------------------------------------------------*/

/**
        @return TRUE of ch is a control character (ascii <= 32)
 */
inline bool controlchar(char ch) {
    return ch <= 32;
}

inline bool is_newick_token(char ch) {
    return ch == ':' || ch == ';' || ch == ',' || ch == ')' || ch == '(' || ch == '[' || ch == ']';
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        print error message then exit program
 */
//void outError(char *error);

/**
        print error message then exit program
 */
void outError(const char *error, bool quit = true);

/**
        print error message then exit program
 */
void outError(string error, bool quit = true);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
        print double error messages then exit program
 */
void outError(const char *error, const char *msg, bool quit = true);

/**
        print double error messages then exit program
 */
void outError(const char *error, string msg, bool quit = true);

/**
        Output a warning message to screen
        @param error warning message
 */
void outWarning(const char *warn);
void outWarning(string warn);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


/**
        generate a random branch length under an exponential distribution
        with mean params.mean_len. Also make sure that the resulting
        length is in the range (params.min_len, params.max_len)
        @return the random branch length
 */
double randomLen(Params &params);

/**
        Compute the logarithm of the factorial of an integer number
        @param num: the number
        @return logarithm of (num! = 1*2*...*num)
 */
double logFac(const int num);

/**
 * Function to randomly select an element in a C++ container
 *
 * @param begin
 * @param end
 * @return
 */
template <typename I>
I random_element(I begin, I end);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*
        Error messages
 */
const char ERR_NO_TAXON[] = "Find no taxon with name ";
const char ERR_NO_AREA[] = "Find no area with name ";
const char ERR_NO_ROOT[] = "Root taxon not found: ";
const char ERR_ROOT_NET[] = "-root option is not available for network";
const char ERR_CONFLICT_ROOT[] = "Tree is already rooted, -o <taxon> is not allowed.";
const char ERR_DUPLICATED_TAXA[] = "Duplicated taxa name in the tree.";
const char ERR_FEW_TAXA[] = "Number of taxa must be greater than 2.";
const char ERR_NO_SPLITS[] = "No splits found!";
const char ERR_FEW_SPLITS[] = "Number of splits must be at least equal to the number of taxa";
const char ERR_NEG_BRANCH[] = "Negative branch length not allowed.";
const char ERR_NO_MEMORY[] = "Not enough memory!";

const char ERR_READ_INPUT[] = "File not found or incorrect input, pls check it again.";
const char ERR_UNEXPECTED_EOF[] = "Unexpected end of file.";
const char ERR_READ_ANY[] = "Unidentified error while reading file, pls check it carefully again.";
const char ERR_WRITE_OUTPUT[] = "Cannot write to file ";

const char ERR_NO_K[] = "You must specify the number of taxa in the PD set.";
const char ERR_TOO_SMALL_K[] = "Size of PD-set must be at least the size of initial set.";
const char ERR_NO_BUDGET[] = "Total budget is not specified or less than zero.";
const char ERR_TOO_SMALL_BUDGET[] = "Not enough budget to conserve the initial set of taxa.";
const char ERR_INTERNAL[] = "Internal error, pls contact authors!";

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
 * convert int to string
 * @param int
 * @return string
 */
string convertIntToString(int number);
string convertInt64ToString(int64_t number);

string convertDoubleToString(double number);

/**
 case-insensitive comparison between two strings
 @return true if two strings are equal.
 */
bool iEquals(const string a, const string b);
    
/**
 *
 * @param SRC
 * @param DEST
 * @return bool
 */
bool copyFile(const char SRC[], const char DEST[]);

/**
 * Check if the file exists
 * @param strFilename
 * @return
 */
bool fileExists(string strFilename);

/**
    check that path is a directory
 */
int isDirectory(const char *path);

/**
    get all file names in a directory (append them to a StrVector)
    @param path directory name
    @param[out] filenames vector of file names
    return 0 if FAIL, otherwise count of files found
 */
size_t getFilesInDir(const char *path, StrVector &filenames);

/**
        convert string to int, with error checking
        @param str original string
        @return the number
 */
int convert_int(const char *str);

/**
       convert string to int, with error checking (but not throwing on error)
       @param str original string
       @param defaultValue value to return if the string isn't numeric
       @return the number
*/
int convert_int_nothrow(const char* str, int defaultValue) throw();

/**
    convert string to int64, with error checking
    @param str original string
    @return the number
 */
int64_t convert_int64(const char *str);

/**
        convert string to int, with error checking
        @param str original string
        @param end_pos end position
        @return the number
 */
int convert_int(const char *str, int &end_pos);

/**
        convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
 */
void convert_int_vec(const char *str, IntVector &vec);

/**
        convert string to int64_t, with error checking
        @param str original string
        @return the number
 */
int64_t convert_int64(const char *str);

/**
        convert string to int64_t, with error checking
        @param str original string
        @param end_pos end position
        @return the number
 */
int64_t convert_int64(const char *str, int &end_pos);

/**
        convert string to double, with error checking
        @param str original string
        @return the double
 */
double convert_double(const char *str);

/**
        convert string to double, with error checking
        @param str original string
        @param end_pos end position
        @return the double
 */
double convert_double(const char *str, int &end_pos);

/**
       convert string to double, with error checking (but not throwing on error)
       @param str original string
       @param defaultValue value to return if the string isn't numeric
       @return the number
*/
double convert_double_nothrow(const char* str, double defaultValue) throw();


/**
        convert comma-separated string to integer vector, with error checking
        @param str original string with integers separated by comma
        @param vec (OUT) integer vector
        @param separator char separating elements
 */
void convert_double_vec(const char *str, DoubleVector &vec, char separator = ',');

/**
 * Convert seconds to hour, minute, second
 * @param sec
 * @return string represent hour, minute, second
 */
string convert_time(const double sec);


/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, int &lower, int &upper, int &step_size);

/**
        convert a string to to range lower:upper:step_size with error checking
        @param str original string
        @param lower (OUT) lower bound of the range
        @param upper (OUT) upper bound of the range
        @param step_size (OUT) step size of the range
 */
void convert_range(const char *str, double &lower, double &upper, double &step_size);

void convert_string_vec(const char *str, StrVector &str_vec, char separator = ',');

/**
    change unusual character in names into underscore (_)
    @param[in/out] name string name
    @return true if renamed, false otherwise
 */
bool renameString(string &name);

/**
        read the file containing branch/split scaling factor and taxa weights
        @param params program parameters
        @param ntaxa total number of taxa
        @param scale (OUT) scaling factor
        @param tax_name (OUT) vector of taxa names
        @param tax_weight (OUT) vector of corresponding taxa weights
 */
void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight);

/**
        read the initial taxa set from the file
        @param params program parameters
        @param ntaxa number of taxa
        @param tax_name (OUT) vector of taxa names
 */
void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name);

/**
        read the initial area set from the file
        @param params program parameters
        @param nareas number of areas
        @param area_name (OUT) vector of area names
 */
void readInitAreaFile(Params &params, int nareas, StrVector &area_name);


/**
        read a list of taxa set from a file, not in nexus format but as follows:
        n1
        tax-name-1
        ...
        tax-name-n1

        n2
        tax-name-1
        ...
        tax-name-n2
        ....

        @param filename file name
        @param sets (OUT) the returned sets of taxa
 */
void readTaxaSets(char *filename, MSetsBlock *sets);

/**
        read areas shared boundary file, in form of a standard distance matrix
        @param file_name file name
        @param areas the read sets block
        @param areas_shared_boundary (OUT) shared boundary length between areas.
                Diagonal elements represent the boundary length of single areas
 */
void readAreasBoundary(char *file_name, MSetsBlock *areas, double *areas_shared_boundary);

/**
        parse program argument into params
        @param argc number of arguments
        @param argv list of arguments
        @param params (OUT) program parameters
 */
void parseArg(int argc, char *argv[], Params &params);

/**
        detect the format of input file
        @param input_file file name
        @return
                IN_NEWICK if file in newick format,
                IN_NEXUS if in nexus format,
                IN_FASTA if in fasta format,
                IN_PHYLIP if in phylip format,
		IN_COUNTSFILE if in counts format (PoMo),
                IN_OTHER if file format unknown.
 */
InputType detectInputFile(const char *input_file);

/**
        if file exists, ask user to overwrite it or not
        @param filename file name
        @return TRUE if agree to overwrite an existing file, or simply file does not exist
 */
bool overwriteFile(char *filename);

/**
        print usage information
        @param argv program arguments list
 */
void usage(char* argv[]);

/**
 *   Print a string into a file
 */
void printString2File(string myString, string filename);

/**
 * print usage for iq-tree
 * @param program arguments list
 * @param full_command TRUE to print all available commands, FALSE to print normal usage dialog
 */
void usage_iqtree(char* argv[], bool full_command);

/**
        parse area name string, where names are separated by commas
        @param area_names a string of name
        @param areas (OUT) a set of name string
 */
void parseAreaName(char *area_names, set<string> &areas);

/**
 * generate 2 different random integer numbers smaller than a specific integer threshold
 * @param size integer threshold
 * @param &first first random integer number
 * @param @second second random integer number
 */
void get2RandNumb(const int size, int &first, int &second);

/*
inline double getCPUTime(clock_t startTime) {
        return double(clock() - startTime) / CLOCKS_PER_SEC;
}*/


/**
 *  Fills the range [first, last) with sequentially increasing values,
 *  starting with value and repetitively evaluating ++value.
 *  Introduced in C++11 --> this is a reimplementation
 */
template<class ForwardIterator, class T>
void iota( ForwardIterator first, ForwardIterator last, T value );

/**
        compute p-value for a chi-square value
        @param chi_square chi-square value
        @param df degree of freedom
        @return p-value
 */
double computePValueChiSquare(double x, int df);

/*--------------------------------------------------------------*/
/* random number generator */
/*--------------------------------------------------------------*/

extern int *randstream;

/**
 * initialize the random number generator
 * @param seed seed for generator
 * @param write_info true to write information, false otherwise (default)
 */
int init_random(int seed, bool write_info = false, int** rstream = NULL);

/**
 * finalize random number generator (e.g. free memory
 */
int finish_random(int *rstream = NULL);

/**
 * returns a random integer in the range [0; n - 1]
 * @param n upper-bound of random number
 */
int random_int(int n, int *rstream = NULL);

/**
 *  return a random integer in the range [a,b]
 */
//int randint(int a, int b);

/**
 * returns a random integer in the range [0; RAND_MAX - 1]
 * = random_int(RAND_MAX)
 */
//int random_int(int *rstream = NULL);

/**
 * returns a random floating-point nuber in the range [0; 1)
 */
double random_double(int *rstream = NULL);

template <class T>
void my_random_shuffle (T first, T last, int *rstream = NULL)
{
	int n = last - first;
	for (int i=n-1; i>0; --i) {
		swap (first[i],first[random_int(i+1, rstream)]);
	}
}

/**
 random resampling according to bootstrap or jackknife
 @param n sample size
 @param[in/out] sample array of size n with frequency of resampling
 @param rstream random number generator stream
*/
void random_resampling(int n, IntVector &sample, int *rstream = NULL);

#define RESAMPLE_NAME ((Params::getInstance().jackknife_prop == 0.0) ? "bootstrap" : "jackknife")
#define RESAMPLE_NAME_I ((Params::getInstance().jackknife_prop == 0.0) ? "Bootstrap" : "Jackknife")
#define RESAMPLE_NAME_UPPER ((Params::getInstance().jackknife_prop == 0.0) ? "BOOTSTRAP" : "JACKKNIFE")

/**
 * generic function for sorting by index
 */
template <class T>
void quicksort_index(T* arr, int* index, int left, int right) {
    int i = left, j = right, tmp2;
    T tmp, pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (pivot < arr[j])
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            tmp2 = index[i];
            index[i] = index[j];
            index[j] = tmp2;
            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quicksort_index(arr, index, left, j);
    if (i < right)
        quicksort_index(arr, index, i, right);
}

/**
 * generic function for sorting by index preseving entries in [first,last)
 * @param first first element
 * @param last last element
 * @param index (OUT) ordered index of elements from smallest to largest
 */
template <class T>
void sort_index(T* first, T* last, int *index) {
    T* x;
    int i;
    T* arr = new T[last - first];
    for (x = first, i = 0; x != last; x++, i++) {
        index[i] = i;
        arr[i] = *x;
    }
    ASSERT(last - first == i);
    quicksort_index(arr, index, 0, (last - first) - 1);
    delete [] arr;
}

/**
 * print the header of summary file
 */
void summarizeHeader(ostream &out, Params &params, bool budget_constraint, InputType analysis_type);

/**
 * print footer of summary file
 */
void summarizeFooter(ostream &out, Params &params);


/**
    remove white space at the beginning and end of the string
    @param str (IN/OUT) string to be trimmed
*/
void trimString(string &str);

/**
    get number of processor cores
*/
int countPhysicalCPUCores();

void print_stacktrace(ostream &out, unsigned int max_frames = 63);

/**
    quicksort template
*/
template<class T1, class T2>
void quicksort(T1* arr, int left, int right, T2* arr2 = NULL) {
      ASSERT(left <= right);
      int i = left, j = right;
      T1 pivot = arr[(left + right) / 2];

      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  T1 tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  if (arr2) {
                      T2 tmp2 = arr2[i];
                      arr2[i] = arr2[j];
                      arr2[j] = tmp2;
                  }
                  i++;
                  j--;
            }
      };

      /* recursion */
      if (left < j)
            quicksort(arr, left, j, arr2);
      if (i < right)
            quicksort(arr, i, right, arr2);
}

/* An optimized version of C̩dric Lauradoux's 64-bit merging3 algorithm
   implemented by Kim Walisch, see:
   http://code.google.com/p/primesieve/source/browse/trunk/src/soe/bithacks.h
   Modified ever so slightly to maintain the same API. Note that
   it assumes the buffer is a multiple of 64 bits in length.
*/
inline uint32_t popcount_lauradoux(unsigned *buf, int n) {
  const uint64_t* data = (uint64_t*) buf;
  uint32_t size = n/(sizeof(uint64_t)/sizeof(int));
  const uint64_t m1  = (0x5555555555555555ULL);
  const uint64_t m2  = (0x3333333333333333ULL);
  const uint64_t m4  = (0x0F0F0F0F0F0F0F0FULL);
  const uint64_t m8  = (0x00FF00FF00FF00FFULL);
  const uint64_t m16 = (0x0000FFFF0000FFFFULL);
  const uint64_t h01 = (0x0101010101010101ULL);

  uint32_t bitCount = 0;
  uint32_t i, j;
  uint64_t count1, count2, half1, half2, acc;
  uint64_t x;
  uint32_t limit30 = size - size % 30;

  // 64-bit tree merging (merging3)
  for (i = 0; i < limit30; i += 30, data += 30) {
    acc = 0;
    for (j = 0; j < 30; j += 3) {
      count1  =  data[j];
      count2  =  data[j+1];
      half1   =  data[j+2];
      half2   =  data[j+2];
      half1  &=  m1;
      half2   = (half2  >> 1) & m1;
      count1 -= (count1 >> 1) & m1;
      count2 -= (count2 >> 1) & m1;
      count1 +=  half1;
      count2 +=  half2;
      count1  = (count1 & m2) + ((count1 >> 2) & m2);
      count1 += (count2 & m2) + ((count2 >> 2) & m2);
      acc    += (count1 & m4) + ((count1 >> 4) & m4);
    }
    acc = (acc & m8) + ((acc >>  8)  & m8);
    acc = (acc       +  (acc >> 16)) & m16;
    acc =  acc       +  (acc >> 32);
    bitCount += (uint32_t)acc;
  }

  // count the bits of the remaining bytes (MAX 29*8) using
  // "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
  // the code uses wikipedia's 64-bit popcount_3() implementation:
  // http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
  for (i = 0; i < size - limit30; i++) {
    x = data[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    bitCount += (uint32_t)((x * h01) >> 56);
  }
  return bitCount;
}

/**
 * combination of memcmp and memcpy.
 * @param destination destination memory to copy to
 * @param source code memory to copy from
 * @param num number of bytes to copy
 * @return TRUE of memory are different, FALSE if identical
 */
bool memcmpcpy(void * destination, const void * source, size_t num);

/**
 *  Generating a unique integer from a pair of 2 integer
 *  This method is called cantor pairing function (see wikepedia).
 *  @param int1 the first integer
 *  @param int2 the second integer
 *  @return the encoding of the 2 integer
 */
int pairInteger(int int1, int int2);

/*
 * Given a model name, look in it for "+F..." and 
 * determine the StateFreqType. Returns FREQ_UNKNOWN if
 * unable to find a good +F... specifier
 */
StateFreqType parseStateFreqFromPlusF(string model_name);

/*
 * Given a string of 4 digits, return a StateFreqType according to
 * equality constraints expressed by those digits.
 * E.g. "1233" constrains pi_G=pi_T (ACGT order, 3rd and 4th equal)
 * which results in FREQ_DNA_2311. "5288" would give the same result.
 */
StateFreqType parseStateFreqDigits(string digits);

/*
 * All params in range [0,1] 
 * returns true if base frequencies have changed as a result of this call
 */
bool freqsFromParams(double *freq_vec, double *params, StateFreqType freq_type);

/*
 * For given freq_type, derives frequency parameters from freq_vec
 * All parameters are in range [0,1] (assuming freq_vec is valid)
 */
void paramsFromFreqs(double *params, double *freq_vec, StateFreqType freq_type);

/* 
 * Given a DNA freq_type and a base frequency vector, alter the
 * base freq vector to conform with the constraints of freq_type
 */
void forceFreqsConform(double *base_freq, StateFreqType freq_type);

/*
 * For given freq_type, how many parameters are needed to
 * determine frequenc vector?
 * BQM 2017-04-28: works for DNA and other data types
 */
 int nFreqParams(StateFreqType freq_type);

/*
 * For freq_type, and given every base must have frequency >= min_freq, set upper
 * and lower bounds for parameters.
 */
 void setBoundsForFreqType(double *lower_bound, 
                           double *upper_bound, 
                           bool *bound_check, 
                           double min_freq, 
                           StateFreqType freq_type);

template <typename T>
string NumberToString ( T Number )
{
    ostringstream ss;
    ss << Number;
    return ss.str();
}

template <typename T>
T StringToNumber ( const string &Text )
{
    istringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

// Calculate logarithm of binomial coefficient N choose i.
double binomial_coefficient_log(unsigned int N, unsigned int i);

// Calculate probability of having k out of N successes when the probability of
// a success is p under the binomial distribution.
double binomial_dist(unsigned int k, unsigned int N, double p);

// Calculate probability of having k out of n successes when there are K out of
// N successes in the pool under the hypergeometric distribution.
double hypergeometric_dist(unsigned int k, unsigned int n, unsigned int K, unsigned int N);

// Calculate the Frobenius norm of an N x N matrix M (flattened, rows
// concatenated) and linearly scaled by SCALE.
double frob_norm (double m[], int n, double scale=1.0);


#endif
