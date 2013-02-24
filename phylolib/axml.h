/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

/** @file axml.h
  * @brief contains various important functions
  * @todo this file will at some point be gone
  */
#ifndef AXML_H_
#define AXML_H_

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>
#include <pmmintrin.h>

#define BYTE_ALIGNMENT 32

#else

#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>

#define BYTE_ALIGNMENT 16

#else
#define BYTE_ALIGNMENT 1
#endif
#endif


#include "genericParallelization.h"

#define MAX_TIP_EV     0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define MAX_LOCAL_SMOOTHING_ITERATIONS     32          /* maximum iterations of smoothings per insert in the */
#define iterations     10          /* maximum iterations of iterations per insert */
#define newzpercycle   1           /* iterations of makenewz per tree traversal */
#define nmlngth        256         /* number of characters in species name */
#define deltaz         0.00001     /* test of net branch length change in update */
#define defaultz       0.9         /* value of z assigned as starting point */
#define unlikely       -1.0E300    /* low likelihood for initialization */


#define SUMMARIZE_LENGTH -3
#define SUMMARIZE_LH     -2
#define NO_BRANCHES      -1

#define MASK_LENGTH 32
#define GET_BITVECTOR_LENGTH(x) ((x % MASK_LENGTH) ? (x / MASK_LENGTH + 1) : (x / MASK_LENGTH))

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */

#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0  
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood




/* 18446744073709551616.0 */

/*4294967296.0*/

/* 18446744073709551616.0 */

/*  2**64 (exactly)  */
/* 4294967296 2**32 */

#define badRear         -1

#define NUM_BRANCHES     16

#define TRUE             1
#define FALSE            0



#define LIKELIHOOD_EPSILON 0.0000001

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

/* ALPHA_MIN is critical -> numerical instability, eg for 4 discrete rate cats                    */
/* and alpha = 0.01 the lowest rate r_0 is                                                        */
/* 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487             */
/* which leads to numerical problems Table for alpha settings below:                              */
/*                                                                                                */
/* 0.010000 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487    */
/* 0.010000 yielded nasty numerical bugs in at least one case !                                   */
/* 0.020000 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875    */
/* 0.030000 0.00000000000000000000476844846859006690412039180149775802624789852441798419292220    */
/* 0.040000 0.00000000000000049522423236954066431210260930029681736928018820007024736185030633    */
/* 0.050000 0.00000000000050625351310359203371872643495343928538368616365517027588794007897377    */
/* 0.060000 0.00000000005134625283884191118711474021861409372524676086868566926568746566772461    */
/* 0.070000 0.00000000139080650074206434685544624965062437960128249869740102440118789672851562    */
/* 0.080000 0.00000001650681201563587066858709818343436959153791576682124286890029907226562500    */
/* 0.090000 0.00000011301977332931251259273962858978301859735893231118097901344299316406250000    */
/* 0.100000 0.00000052651925834844387815526344648331402709118265192955732345581054687500000000    */


#define ALPHA_MIN    0.02
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0000001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.001

/* 
   previous values between 0.001 and 0.000001

   TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA 
   previous value of FREQ_MIN was: 0.000001, but this seemed to cause problems with some 
   of the 7-state secondary structure models with some rather exotic small toy test datasets,
   on the other hand 0.001 caused problems with some of the 16-state secondary structure models

   For some reason the frequency settings seem to be repeatedly causing numerical problems
   
*/

#define ITMAX 100



#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))

#define FABS(x) fabs(x)

#ifdef _USE_FPGA_LOG
extern double log_approx (double input);
#define LOG(x)  log_approx(x)
#else
#define LOG(x)  log(x)
#endif


#ifdef _USE_FPGA_EXP
extern double exp_approx (double x);
#define EXP(x)  exp_approx(x)
#else
#define EXP(x)  exp(x)
#endif


#define LOGF(x) logf(x)


#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#define programName        "RAxML-Light"
#define programVersion     "1.0.5"
#define programDate        "June 2011"


#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  CALC_BIPARTITIONS          3
#define  SPLIT_MULTI_GENE           4
#define  CHECK_ALIGNMENT            5
#define  PER_SITE_LL                6
#define  PARSIMONY_ADDITION         7
#define  CLASSIFY_ML                9
#define  DISTANCE_MODE              11
#define  GENERATE_BS                12
#define  COMPUTE_ELW                13
#define  BOOTSTOP_ONLY              14
#define  COMPUTE_LHS                17
#define  COMPUTE_BIPARTITION_CORRELATION 18
#define  THOROUGH_PARSIMONY         19
#define  COMPUTE_RF_DISTANCE        20
#define  MORPH_CALIBRATOR           21
#define  CONSENSUS_ONLY             22
#define  MESH_TREE_SEARCH           23
#define  FAST_SEARCH                24
#define  MORPH_CALIBRATOR_PARSIMONY 25
#define  SH_LIKE_SUPPORTS           28

#define  GPU_BENCHMARK              29

#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10


#define DAYHOFF    0
#define DCMUT      1
#define JTT        2
#define MTREV      3
#define WAG        4
#define RTREV      5
#define CPREV      6
#define VT         7
#define BLOSUM62   8
#define MTMAM      9
#define LG         10
#define MTART      11
#define MTZOA      12
#define PMB        13
#define HIVB       14
#define HIVW       15
#define JTTDCMUT   16
#define FLU        17 
#define AUTO       18
#define GTR        19  /* GTR always needs to be the last one */

#define NUM_PROT_MODELS 20

/* bipartition stuff */

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4



/* bootstopping stuff */

#define BOOTSTOP_PERMUTATIONS 100
#define START_BSTOP_TEST      10

#define FC_THRESHOLD          99
#define FC_SPACING            50
#define FC_LOWER              0.99
#define FC_INIT               20

#define FREQUENCY_STOP 0
#define MR_STOP        1
#define MRE_STOP       2
#define MRE_IGN_STOP   3

#define MR_CONSENSUS 0
#define MRE_CONSENSUS 1
#define STRICT_CONSENSUS 2



/* bootstopping stuff end */


#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define MIN_MODEL        -1
#define BINARY_DATA      0
#define DNA_DATA         1
#define AA_DATA          2
#define SECONDARY_DATA   3
#define SECONDARY_DATA_6 4
#define SECONDARY_DATA_7 5
#define GENERIC_32       6
#define GENERIC_64       7
#define MAX_MODEL        8

#define SEC_6_A 0
#define SEC_6_B 1
#define SEC_6_C 2
#define SEC_6_D 3
#define SEC_6_E 4

#define SEC_7_A 5
#define SEC_7_B 6
#define SEC_7_C 7
#define SEC_7_D 8
#define SEC_7_E 9
#define SEC_7_F 10

#define SEC_16   11
#define SEC_16_A 12
#define SEC_16_B 13
#define SEC_16_C 14
#define SEC_16_D 15
#define SEC_16_E 16
#define SEC_16_F 17
#define SEC_16_I 18
#define SEC_16_J 19
#define SEC_16_K 20

#define ORDERED_MULTI_STATE 0
#define MK_MULTI_STATE      1
#define GTR_MULTI_STATE     2





#define CAT         0
#define GAMMA       1
#define GAMMA_I     2

/* recomp */
#define SLOT_UNUSED            -2  /* value to mark an available vector */
#define NODE_UNPINNED          -3  /* marks an inner node as not available in RAM */
#define INNER_NODE_INIT_STLEN  -1  /* initialization */

#define MIN_RECOM_FRACTION     0.1 /* at least this % of inner nodes will be allocated in RAM */
#define MAX_RECOM_FRACTION     1.0 /* always 1, just there for boundary checks */
#define MEM_APROX_OVERHEAD     1.3 /* TODOFER can we measure this empirically? */


typedef  int pl_boolean;

typedef struct
{
  int numVectors;      /* #inner vectors in RAM*/
  int *iVector;        /* size: numVectors, stores node id || SLOT_UNUSED  */
  int *iNode;          /* size: inner nodes, stores slot id || NODE_UNPINNED */
  int *stlen;          /* #tips behind the current orientation of the indexed inner node (subtree size/cost) */ 
  int *unpinnable;     /* size:numVectors , TRUE if we dont need the vector */
  int maxVectorsUsed;  
  pl_boolean allSlotsBusy; /*on if all slots contain an ancesctral node (the usual case after first full traversal) */ 
#ifdef _DEBUG_RECOMPUTATION
  double pinTime;
  double recomStraTime;
#endif
} recompVectors;
/* E recomp */




typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct ent *next;
};

typedef struct ent entry;

typedef unsigned int hashNumberType;



/*typedef uint_fast32_t parsimonyNumber;*/

#define PCF 32

/*
  typedef uint64_t parsimonyNumber;

  #define PCF 16


typedef unsigned char parsimonyNumber;

#define PCF 2
*/

typedef struct
{
  hashNumberType tableSize;
  entry **table;
  hashNumberType entryCount;
}
  pl_hashtable;


struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};


typedef struct stringEnt stringEntry;
 
typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;





typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;

/** @brief To be commented
  *
  * The rest is her
  * @todo remove this line when finished commenting
  */
typedef struct
{
  int tipCase;                  /**< @brief What is this? */
  int pNumber;                  /**< @brief Or this */
  int qNumber;
  int rNumber;
  double qz[NUM_BRANCHES];
  double rz[NUM_BRANCHES];
  /* recom */
  int slot_p;
  int slot_q;
  int slot_r;
  /* E recom */
} traversalInfo;

typedef struct
{
  traversalInfo *ti;
  int count;
  int functionType;
  pl_boolean traversalHasChanged;
  pl_boolean *executeModel;
  double  *parameterValues;
} traversalData;


struct noderec;



typedef struct
{
 

  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;








typedef struct
{
  pl_boolean valid;
  int partitions;
  int *partitionList;
}
  linkageData;

typedef struct
{
  int entries;
  linkageData* ld;
}
  linkageList;



  /* the data structure below is fundamental for representing trees 
     in the library!

     Inner nodes are represented by three instances of the nodeptr data structure that is linked 
     via a cyclic list using the next pointer.

     So for building an inner node of the tree we need to allocate three nodeptr 
     data structures and link them together, e.g.:

     assuming that we have allocated space for an inner node 
     for nodeptr pointers p1, p2, p3, 

     we would then link them like this:

     p1->next = p2;
     p2->next = p3;
     p3->next = p1;

     also note that the node number that identifies the inner node 
     needs to be set to the same value.

     for n taxa, tip nodes are enumarated/indexed from 1....n,
     and inner node inbdices start at n+1. Assuming that we have 10 taxa 
     and this is our first inner node, we'd initialize the number as follows:

     p1->number = 11;
     p2->number = 11;
     p3->number = 11;

     Note that the node number is important for indexing tip sequence data as well as inner likelihood vectors 
     and that it is this number (the index) that actually gets stored in the traversal descriptor.

     Tip nodes are non-cyclic nodes that simply consist of one instance/allocation of nodeptr.

     if we have allocated a tip data structure nodeptr t1, 
     we would initialize it as follows:

     t1->number = 1;

     t1->next = NULL;

     now let's assume that we want to build a four taxon tree with tips t1, t2, t3, t4 
     and inner nodes (p1,p2,p3) and (q1,q2,q3).

     we first build the tips:

     t1->number = 1;
     t1->next = NULL;
     
     t2->number = 2;
     t2->next = NULL;

     t3->number = 3;
     t3->next = NULL;

     t4->number = 4;
     t4->next = NULL;
     
     now the first inner node

     p1->next = p2;
     p2->next = p3;
     p3->next = p1;    

     p1->number = 5;
     p2->number = 5;
     p3->number = 5;

     and the second inner node.

     q1->next = q2;
     q2->next = q3;
     q3->next = q1;    

     q1->number = 6;
     q2->number = 6;
     q3->number = 6;
     
     now we need to link the nodes together such that they form a tree, let's assume we want ((t1,t2), (t3, t4));

     we will have to link the nodes via the so-called back pointer,
     i.e.:

     let's connect node p with t1 and t2

     t1->back = p1;
     t2->back = p2;

     and vice versa:

     p1->back = t1;
     p2->back = t2;

     let's connect node p with node q:

     p3->back = q3;

     and vice versa:

     q3->back = p3;

     and now let's connect node q with tips t3 and t4:

     q1->back = t3;
     q2->back = t4;

     and vice versa:

     t3->back = q1;
     t4->back = q2;

     What remains to be done is to set up the branch lengths.
     Using the data structure below, we always have to store the 
     branch length twice for each "topological branch" unfortunately.

     Assuming that we are only estimating a single branch across all partitions 
     we'd just set the first index of the branch length array z[NUM_BRANCHES].

     e.g., 

     t3->z[0] = q1->z[0] = 0.9;

     the above operation for connecting nodes is implemented in functions hookup() which will set 
     the back pointers of two nodes that are to be connected as well as the branch lengths.

     The branchInfo data field is a pointer to a data-structure that stores meta-data and requires 
     the tree not to change while it is being used.
     
     Also, this pointer needs to be set by doing a full tree traversal on the tree.

     Note that q1->bInf == t3->bInf in the above example.

     The hash number is used for mapping bipartitions to a hash table as described in the following paper:

     A. Aberer, N. Pattengale, A. Stamatakis: "Parallelized phylogenetic post-analysis on multi-core architectures". Journal of Computational Science 1, 107-114, 2010.
     
     The support data field stores the support value for the branch associated with each nodeptr structure.
     Note that support always refers to branches. 

     Thus for consistency, q3->support must be equal to p3->support;

     Finally, the three char fields x, xPars and xBips are very very important!

     They are used to denote the presence/absence or if you want, direction of the 
     parsimony, bipartition, or likelihood vector at a node with respect to the virtual root.

     Essentially, they are just used as single presence/absence bits and ONLY for inner nodes!

     When setting up new inner nodes, one of the three pointers in the cyclic list must 
     have x = 1 and the other two x = 0;

     in the above example we could set:

     p1->x = 0;
     p2->x = 0;
     p3->x = 1;

     q1->x = 0;
     q2->x = 0;
     q3->x = 1;

     This would mean that the virtual root is located at the inner branch of the four taxon tree ((t1,t2),(t3,t4));

     When we re-root the tree at some other branch we need to update the location of the x pointer that is set to 1.

     This means if we root the tree at the branch leading to t1 we would set 

     p1->x = 1;
     p2->x = 0;
     p3->x = 0;

     the values for q remaon unchanged since q3 is still pointing toward the root.

     When we re-locate the root to branch p1 <-> t1 the fact that we have to "rotate" the x value that is set to 1
     to another node of the cyclic list representing the abstract topological node p, also tells us that we 
     need to re-compute the conditional likelihood array for p. 

     Note that, only one likelihood or parsimony array is stored per inner node and the location of x essentially tells us which subtree 
     it summarizes, if p1->x == 1, it summarizes subtree (t2, (t3, t4)), if p3->x = 1 the likelihood vector associated with 
     node p summarizes subtree (t1, t2).

  */
    


typedef  struct noderec
{
 
  branchInfo      *bInf;
  double           z[NUM_BRANCHES];
#ifdef _BAYESIAN 
  double           z_tmp[NUM_BRANCHES];
#endif 
  struct noderec  *next;
  struct noderec  *back;
  hashNumberType   hash;
  int              support;
  int              number;
  char             x;
  char             xPars;
  char             xBips;
}
  node, *nodeptr;

typedef struct
  {
    double lh;
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;





typedef unsigned int parsimonyNumber;


typedef struct {
  /* ALIGNMENT DATA */
  /* This depends only on the type of data in this partition of the alignment */
  int     dataType;         /* e.g. DNA_DATA, AA_DATA, etc  within range (MIN_MODEL, MAX_MODEL)  */
  int     states;           /* Number of states in inner vectors */
  int     maxTipStates;     /* Number of undetermined states (Possible states at the tips) */
  /* These are the boundaries of the partition itself, sites within these boudaries must share the type of data */
  char   *partitionName;
  int     lower;            /* starting position of the partition within [1, tr->originalCrunchedLength] */
  int     upper;            /* ending position  */
  int     width;            /* upper - lower, possibly we dont need this, number of site patterns*/
  int    *wgt;              /* Number of occurencies of each site pattern */
  double *empiricalFrequencies;    /* empirical Frequency of each state according to this alignment partition */


  /* MODEL OF RATE HETEROGENETY, We use either GAMMA or PSR */
  /* Rate heterogenety: Per Site Categories (PSR) model aka CAT, see updatePerSiteRates() */
  /* Rate of site i is given by perSiteRates[rateCategory[i]] */
  double *perSiteRates;     /* Values of rates*/
  int    *rateCategory;     /* Category index for each site */
  int     numberOfCategories;/* size of the set of possible categories */
  /* Rate heterogenety: GAMMA model of rate heterogenety */
  double alpha;             /* parameter to be optimized */
  double *gammaRates;       /* 4 gamma categories (rates), computed given an alpha*/


  /* TRANSITION MODEL: We always assume General Time Reversibility */
  /* Transistion probability matrix: P(t) = exp(Qt)*/
  /* Branch length t is the expected number of substitutions per site */
  /* Pij(t) is the probability of going from state i to state j in a branch of length t */
  /* Relative substitution rates (Entries in the Q matrix) */
  /* In GTR we can write Q = S * D, where S is a symmetrical matrix and D a diagonal with the state frequencies */
  double *substRates;       /* Entries in S, e.g. 6 free parameters in DNA */   
  double *frequencies;      /* State frequencies, entries in D, are initialized as empiricalFrequencies */
  /* Matrix decomposition: Explanation of the mathematical background? */
  double *EIGN;             /* eigenvalues */
  double *EV;               /* eigenvectors */
  double *EI;
  double *left;  
  double *right;
  double *tipVector; 
  /* Protein specific ?? */
  int     protModels;
  int     autoProtModels;
  int     protFreqs;
  /* specific for secondary structures ?? */
  pl_boolean nonGTR;
  int    *symmetryVector;
  int    *frequencyGrouping;


  /* LIKELIHOOD VECTORS */
  /* partial LH Inner vectors / ancestral vectors, we have 2*tips - 3 inner nodes */
  double          **xVector;          /* Probability entries for inner nodes */
  unsigned char   **yVector;          /* Tip entries (sequence) for tip nodes */
  unsigned int     *globalScaler;     /* Counters for scaling operations done at node i */
  /* These are for the saveMemory option (tracking gaps to skip computations and memory) */
  size_t           *xSpaceVector;     
  int               gapVectorLength;
  unsigned int     *gapVector;
  double           *gapColumn; 

  /* Parsimony vectors at each node */
  size_t parsimonyLength;
  parsimonyNumber *parsVect; 

  /* These buffers of size width are used to pre-compute ?? */
  double *sumBuffer; 
  double *ancestralBuffer;
} pInfo;



typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;


#define REARR_SETTING 1
#define FAST_SPRS     2
#define SLOW_SPRS     3


typedef struct {
 
  int state;

  /*unsigned int vLength;*/
  double accumulatedTime;  
  int rearrangementsMax;
  int rearrangementsMin;
  int thoroughIterations;
  int fastIterations;
  int mintrav;
  int maxtrav;
  int bestTrav;
  double startLH; 
  double lh;
  double previousLh;
  double difference;
  double epsilon;  
  pl_boolean impr;
  pl_boolean cutoff;  
       
  double tr_startLH;
  double tr_endLH;
  double tr_likelihood;
  double tr_bestOfNode;  
  double tr_lhCutoff;
  double tr_lhAVG;
  double tr_lhDEC;
  int    tr_NumberOfCategories;
  int    tr_itCount;  
  int    tr_doCutoff;
  int    tr_thoroughInsertion;
  int    tr_optimizeRateCategoryInvocations;

 
  /* prevent users from doing stupid things */

 
  int searchConvergenceCriterion;
  int rateHetModel;
  int maxCategories;
  int NumberOfModels;
  int numBranches;
  int originalCrunchedLength;    
  int mxtips;
  char seq_file[1024];
} checkPointState;



/* recomp */
#ifdef _DEBUG_RECOMPUTATION
typedef struct {
  unsigned long int numTraversals;
  unsigned long int tt;
  unsigned long int ti;
  unsigned long int ii;
  unsigned int *travlenFreq;
} traversalCounter;
#endif
/* E recomp */


typedef  struct  {

  int *ti;

  /* recomp */
  recompVectors *rvec;            /* this data structure tracks which vectors store which nodes */
  float maxMegabytesMemory;         /* User says how many MB in main memory should be used */
  float vectorRecomFraction;      /* vectorRecomFraction ~= 0.8 * maxMegabytesMemory  */
  pl_boolean useRecom;               /* ON if we apply recomputation of ancestral vectors*/
#ifdef _DEBUG_RECOMPUTATION 
  traversalCounter *travCounter;
  double stlenTime;
#endif
  /* E recomp */

 
  pl_boolean saveMemory;
  int              startingTree;
  long             randomNumberSeed;

  double          *lhs;
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
  int             *rateCategory;
  int             *aliaswgt;    /* weight by pattern */ 
  pl_boolean    manyPartitions;

  pl_boolean grouped;
  pl_boolean constrained;
  int threadID;
  volatile int numberOfThreads;

#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))
    
  int *partitionAssignment;     
 
  unsigned char *y_ptr; 
  
  double lower_spacing;
  double upper_spacing; 

  double *ancestralVector;
#endif
  


 
  
 

  stringHashtable  *nameHash;

  pInfo            *partitionData;
  

  char             *secondaryStructureInput;

  pl_boolean          *executeModel;

  double           *perPartitionLH;

  traversalData    td[1];

  int              maxCategories;
  int              categories;

  double           coreLZ[NUM_BRANCHES];
  int              numBranches;                 /* Number of length values per branch. 
                                                   Currently can be only 1 or number of
                                                   partitions */
  
  
 
  branchInfo	   *bInf;

  int              multiStateModel;


  pl_boolean curvatOK[NUM_BRANCHES];
  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */

  
  double           *fracchanges;

  /* model stuff end */

  unsigned char             **yVector;
  int              secondaryStructureModel;
  int              originalCrunchedLength;
 
 
  int              *secondaryStructurePairs;


  double            *partitionContributions;
  double            fracchange;
  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;

  double           startLH;
  double           endLH;
  double           likelihood;
  
 
  node           **nodep;
  nodeptr          nodeBaseAddress;
  node            *start;
  int              mxtips;  

  int              *constraintVector;
  int              numberOfSecondaryColumns;
  pl_boolean          searchConvergenceCriterion;
  int              ntips;
  int              nextnode;  
  int              NumberOfModels;    

  pl_boolean          bigCutoff;
  pl_boolean          partitionSmoothed[NUM_BRANCHES];
  pl_boolean          partitionConverged[NUM_BRANCHES];
  pl_boolean          rooted;
  pl_boolean          doCutoff;
 
  double         gapyness;

  char **nameList;
  char *tree_string;
  char *tree0;
  char *tree1;
  int treeStringLength;
 
  unsigned int bestParsimony;
  unsigned int *parsimonyScore;
  
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;

  double zqr[NUM_BRANCHES];
  double currentZQR[NUM_BRANCHES];

  double currentLZR[NUM_BRANCHES];
  double currentLZQ[NUM_BRANCHES];
  double currentLZS[NUM_BRANCHES];
  double currentLZI[NUM_BRANCHES];
  double lzs[NUM_BRANCHES];
  double lzq[NUM_BRANCHES];
  double lzr[NUM_BRANCHES];
  double lzi[NUM_BRANCHES];


  unsigned int **bitVectors;

  unsigned int vLength;

  pl_hashtable *h;
 
  int optimizeRateCategoryInvocations;

  checkPointState ckp;
  pl_boolean thoroughInsertion;
  pl_boolean useMedian;

} tree;


/***************************************************************/

typedef struct {
  int partitionNumber;
  int partitionLength;
} partitionType;

typedef struct
{
  double z[NUM_BRANCHES];
  nodeptr p, q;
  int cp, cq;
}
  connectRELL, *connptrRELL;

typedef  struct
{
  connectRELL     *connect; 
  int             start;
  double          likelihood;
}
  topolRELL;


typedef  struct
{
  int max;
  topolRELL **t;
}
  topolRELL_LIST;





/**************************************************************/



typedef struct conntyp {
    double           z[NUM_BRANCHES];           /* branch length */
    node            *p, *q;       /* parent and child sectors */
    void            *valptr;      /* pointer to value of subtree */
    int              descend;     /* pointer to first connect of child */
    int              sibling;     /* next connect from same parent */
    } pl_connect, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    pl_connect         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */

    } topol;

/* small helper data structure for printing out/downstream use of marginal ancestral probability vectors */
/* it is allocated as an array that has the same length as the input alignment and can be used to 
   index the ancestral states for each position/site/pattern */

typedef struct {
  double *probs; /* marginal ancestral states */
  char c; /* most likely stated, i.e. max(probs[i]) above */
  int states; /* number of states for this position */
} ancestralState;


typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    pl_boolean          improved;
    } bestlist;

#define randomTree    0
#define givenTree     1 
#define parsimonyTree 2

typedef  struct {
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  pl_boolean          initialSet;
  int              mode; 
  pl_boolean        perGeneBranchLengths;
  pl_boolean        permuteTreeoptimize; 
  pl_boolean        compressPatterns;
  double         likelihoodEpsilon;
  pl_boolean        useCheckpoint;
 
#ifdef _BAYESIAN 
  pl_boolean       bayesian;
  int           num_generations;
#endif

} analdef;




typedef struct 
{
  int leftLength;         /* s^2 */
  int rightLength;/* s^2 */
  int eignLength;/* s */
  int evLength;
  int eiLength;
  int substRatesLength;   /* (s^2 - s)/2 free model parameters for matrix Q i.e. substitution rates */
  int frequenciesLength;  /* s frequency of each state */ 
  int tipVectorLength;    /* ASK */
  int symmetryVectorLength;
  int frequencyGroupingLength;

  pl_boolean nonGTR;

  int undetermined;

  const char *inverseMeaning;

  int states;   /* s */

  pl_boolean smoothFrequencies;

  const unsigned  int *bitVector;

} partitionLengths;

/****************************** FUNCTIONS ****************************************************/

#ifdef _BAYESIAN 
extern void mcmc(tree *tr, analdef *adef);
#endif

#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))
pl_boolean isThisMyPartition(tree *localTree, int tid, int model);
void printParallelTimePerRegion(); 
#endif

extern void computePlacementBias(tree *tr, analdef *adef);

extern int lookupWord(char *s, stringHashtable *h);

extern void getDataTypeString(tree *tr, int model, char typeOfData[1024]);

extern unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength);
extern int countTips(nodeptr p, int numsp);
extern entry *initEntry(void);
extern void computeRogueTaxa(tree *tr, char* treeSetFileName, analdef *adef);
extern unsigned int precomputed16_bitcount(unsigned int n, char *bits_in_16bits);


/* Handling branch lengths*/
extern double get_branch_length(tree *tr, nodeptr p, int partition_id);
extern void set_branch_length(tree *tr, nodeptr p, int partition_id, double bl);


extern size_t discreteRateCategories(int rateHetModel);

extern const partitionLengths * getPartitionLengths(pInfo *p);
extern pl_boolean getSmoothFreqs(int dataType);
extern const unsigned int *getBitVector(int dataType);
extern int getUndetermined(int dataType);
extern int getStates(int dataType);
extern char getInverseMeaning(int dataType, unsigned char state);
extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );

extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
extern pl_boolean whitechar ( int ch );
extern void printResult ( tree *tr, analdef *adef, pl_boolean finalPrint );
extern void printBootstrapResult ( tree *tr, analdef *adef, pl_boolean finalPrint );
extern void printBipartitionResult ( tree *tr, analdef *adef, pl_boolean finalPrint );
extern void printLog ( tree *tr);
extern void printStartingTree ( tree *tr, analdef *adef, pl_boolean finalPrint );
extern void writeInfoFile ( analdef *adef, tree *tr, double t );
/* extern int main ( int argc, char *argv[] ); */
extern void calcBipartitions ( tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName );
extern void initReversibleGTR (tree *tr, int model);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void makeGammaCats (double alpha, double *gammaRates, int K, pl_boolean useMedian);
extern void initModel ( tree *tr, double **empiricalFrequencies);
extern void doAllInOne ( tree *tr, analdef *adef );

extern void classifyML(tree *tr, analdef *adef);

extern void resetBranches ( tree *tr );
extern void modOpt ( tree *tr, double likelihoodEpsilon);



extern void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( tree *tr, analdef *adef );
extern void thoroughOptimization ( tree *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( tree *tr, int mintrav, int maxtrav);

extern int checker ( tree *tr, nodeptr p );
extern pl_boolean tipHomogeneityChecker ( tree *tr, nodeptr p, int grouping );
extern void makeRandomTree ( tree *tr);
extern void nodeRectifier ( tree *tr );
extern void makeParsimonyTreeFast(tree *tr);
extern void allocateParsimonyDataStructures(tree *tr);
extern void freeParsimonyDataStructures(tree *tr);
extern void parsimonySPR(nodeptr p, tree *tr);

extern FILE *myfopen(const char *path, const char *mode);


extern pl_boolean initrav ( tree *tr, nodeptr p );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern void update ( tree *tr, nodeptr p );
extern void smooth ( tree *tr, nodeptr p );
extern void smoothTree ( tree *tr, int maxtimes );
extern void localSmooth ( tree *tr, nodeptr p, int maxtimes );
extern pl_boolean localSmoothMulti(tree *tr, nodeptr p, int maxtimes, int model);

extern void smoothRegion ( tree *tr, nodeptr p, int region );
extern void regionalSmooth ( tree *tr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( tree *tr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( tree *tr, nodeptr p );
extern pl_boolean insertBIG ( tree *tr, nodeptr p, nodeptr q, int numBranches);
extern pl_boolean insertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern pl_boolean testInsertBIG ( tree *tr, nodeptr p, nodeptr q );
extern void addTraverseBIG ( tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( tree *tr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern pl_boolean testInsertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( tree *tr );
extern int determineRearrangementSetting ( tree *tr, analdef *adef, bestlist *bestT, bestlist *bt );
extern void computeBIGRAPID ( tree *tr, analdef *adef, pl_boolean estimateModel);
extern pl_boolean treeEvaluate ( tree *tr, int maxSmoothIterations );
extern pl_boolean treeEvaluatePartition ( tree *tr, double smoothFactor, int model );

extern void meshTreeSearch(tree *tr, analdef *adef, int thorough);

extern void initTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, tree *tr, int index );

extern int  saveBestTree (bestlist *bt, tree *tr);
extern int  recallBestTree (bestlist *bt, int rank, tree *tr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern pl_boolean freeBestTree ( bestlist *bt );


extern char *Tree2String ( char *treestr, tree *tr, nodeptr p, pl_boolean printBranchLengths, pl_boolean printNames, pl_boolean printLikelihood, 
			   pl_boolean rellTree, pl_boolean finalPrint, int perGene, pl_boolean branchLabelSupport, pl_boolean printSHSupport);
extern void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission);
void printTopology(tree *tr, pl_boolean printInner);



extern int treeReadLen (FILE *fp, tree *tr, pl_boolean readBranches, pl_boolean readNodeLabels, pl_boolean topologyOnly);
extern void treeReadTopologyString(char *treeString, tree *tr);
extern pl_boolean treeReadLenMULT ( FILE *fp, tree *tr, analdef *adef );

extern void getStartingTree ( tree *tr);
extern double treeLength(tree *tr, int model);

extern void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef);
extern pl_boolean bootStop(tree *tr, pl_hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength);
extern void computeConsensusOnly(tree *tr, char* treeSetFileName, analdef *adef);
extern double evaluatePartialGeneric (tree *, int i, double ki, int _model);
extern void evaluateGeneric (tree *tr, nodeptr p, pl_boolean fullTraversal);
extern void newviewGeneric (tree *tr, nodeptr p, pl_boolean masked);

extern void newviewGenericAncestral(tree *tr, nodeptr p);
extern void newviewAncestralIterative(tree *tr);
extern void printAncestralState(nodeptr p, pl_boolean printStates, pl_boolean printProbs, tree *tr);

extern void newviewGenericMulti (tree *tr, nodeptr p, int model);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, pl_boolean mask);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern double evaluateGenericVector (tree *tr, nodeptr p);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern pl_boolean isTip(int number, int maxTips);

/* recom functions */
extern void computeTraversal(tree *tr, nodeptr p, pl_boolean partialTraversal);
extern void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, pl_boolean partialTraversal, recompVectors *rvec, pl_boolean useRecom);
extern void allocRecompVectorsInfo(tree *tr);
extern void allocTraversalCounter(tree *tr);
extern pl_boolean getxVector(recompVectors *rvec, int nodenum, int *slot, int mxtips);
extern pl_boolean needsRecomp(pl_boolean recompute, recompVectors *rvec, nodeptr p, int mxtips);
extern void unpinNode(recompVectors *v, int nodenum, int mxtips);
extern void protectNode(recompVectors *rvec, int nodenum, int mxtips);

extern void computeTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec, int *count);
extern void computeFullTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec);
extern void printTraversalInfo(tree *tr);
extern void countTraversal(tree *tr);

extern void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, pl_boolean saveMem, int maxCat, const int states);

extern void newviewIterative(tree *tr, int startIndex);

extern void evaluateIterative(tree *tr);

extern void *malloc_aligned( size_t size);

extern void storeExecuteMaskInTraversalDescriptor(tree *tr);
extern void storeValuesInTraversalDescriptor(tree *tr, double *value);


extern void makenewzIterative(tree *);
extern void execCore(tree *, volatile double *dlnLdlz, volatile double *d2lnLdlz2);

extern void determineFullTraversal(nodeptr p, tree *tr);
/*extern void optRateCat(tree *, int i, double lower_spacing, double upper_spacing, double *lhs);*/





extern double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model);
extern void evaluateGenericVectorIterative(tree *, int startIndex, int endIndex);
extern void categorizeIterative(tree *, int startIndex, int endIndex);

extern void fixModelIndices(tree *tr, int endsite, pl_boolean fixRates);
extern void calculateModelOffsets(tree *tr);
extern void gammaToCat(tree *tr);
extern void catToGamma(tree *tr, analdef *adef);


extern nodeptr findAnyTip(nodeptr p, int numsp);

extern void parseProteinModel(analdef *adef);



extern void computeNextReplicate(tree *tr, long *seed, int *originalRateCategories, int *originalInvariant, pl_boolean isRapid, pl_boolean fixRates);
/*extern void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);*/

extern void putWAG(double *ext_initialRates);

extern void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant);
extern void parseSecondaryStructure(tree *tr, analdef *adef, int sites);
extern void printPartitions(tree *tr);
extern void compareBips(tree *tr, char *bootStrapFileName, analdef *adef);
extern void computeRF(tree *tr, char *bootStrapFileName, analdef *adef);


extern  unsigned int **initBitVector(int mxtips, unsigned int *vectorLength);
extern pl_hashtable *copyHashTable(pl_hashtable *src, unsigned int vectorLength);
extern pl_hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(pl_hashtable *h, int state);
extern double convergenceCriterion(pl_hashtable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
extern void freeHashTable(pl_hashtable *h);
extern stringHashtable *initStringHashTable(hashNumberType n);
extern void addword(char *s, stringHashtable *h, int nodeNumber);


extern void printBothOpen(const char* format, ... );
extern void initRateMatrix(tree *tr);

extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, pl_hashtable *h, int treeNumber, int function, branchInfo *bInf,
				    int *countBranches, int treeVectorLength, pl_boolean traverseOnly, pl_boolean computeWRF, int processID);


extern  unsigned int bitcount_32_bit(unsigned int i); 
/* extern inline unsigned int bitcount_64_bit(unsigned long i); */

extern FILE *getNumberOfTrees(tree *tr, char *fileName, analdef *adef);

extern void writeBinaryModel(tree *tr);
extern void readBinaryModel(tree *tr);
extern void treeEvaluateRandom (tree *tr, double smoothFactor);
extern void treeEvaluateProgressive(tree *tr);

extern void testGapped(tree *tr);

extern pl_boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen);
extern pl_boolean compatible(entry* e1, entry* e2, unsigned int bvlen);

extern void perSiteLogLikelihoods(tree *tr, double *logLikelihoods);

extern int *permutationSH(tree *tr, int nBootstrap, long _randomSeed);
extern void perSiteLogLikelihoodsPthreads(tree *tr, double *lhs, int n, int tid);
extern void updatePerSiteRates(tree *tr, pl_boolean scaleRates);

extern void restart(tree *tr);

extern double getBranchLength(tree *tr, int perGene, nodeptr p);

#ifdef _IPTOL
extern void writeCheckpoint();
#endif

#ifdef _WAYNE_MPI

extern pl_boolean computeBootStopMPI(tree *tr, char *bootStrapFileName, analdef *adef, double *pearsonAverage);

#endif


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS) )



/* work tags for parallel regions */

#define THREAD_NEWVIEW                0        
#define THREAD_EVALUATE               1
#define THREAD_MAKENEWZ               2 
#define THREAD_MAKENEWZ_FIRST         3
#define THREAD_RATE_CATS              4
#define THREAD_COPY_RATE_CATS         5
#define THREAD_COPY_INIT_MODEL        6
#define THREAD_INIT_PARTITION         7
#define THREAD_OPT_ALPHA              8
#define THREAD_OPT_RATE               9
#define THREAD_COPY_ALPHA             10
#define THREAD_COPY_RATES             11
#define THREAD_PER_SITE_LIKELIHOODS   12
#define THREAD_NEWVIEW_ANCESTRAL      13
#define THREAD_GATHER_ANCESTRAL       14
#define THREAD_EXIT_GRACEFULLY        15

void threadMakeVector(tree *tr, int tid);
void threadComputeAverage(tree *tr, int tid);
void threadComputePearson(tree *tr, int tid);

extern void masterBarrier(int jobType, tree *tr);

#endif

#if (defined(_FINE_GRAIN_MPI) || (_USE_PTHREADS))

pl_boolean workerTrap(tree *tr); 
void initMPI(int argc, char *argv[]); 
void initializePartitions(tree *tr, tree *localTree, int tid, int n); 
void multiprocessorScheduling(tree *tr, int tid); 
void computeFraction(tree *localTree, int tid, int n); 
void computeFractionMany(tree *localTree, int tid); 
void initializePartitionsMaster(tree *tr, tree *localTree, int tid, int n); 
void startPthreads(tree *tr); 

typedef struct
{
  tree *tr;
  int threadNumber;
}
  threadData;
extern void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid);
void allocNodex(tree *tr, int tid, int n);
#endif


#ifdef __AVX
void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
		       double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
		       unsigned char *tipX1, unsigned char *tipX2,
		       int n,  double *left, double *right, int *wgt, int *scalerIncrement);


void newviewGenericCATPROT_AVX(int tipCase, double *extEV,
			       int *cptr,
			       double *x1, double *x2, double *x3, double *tipVector,
			       unsigned char *tipX1, unsigned char *tipX2,
			       int n, double *left, double *right, int *wgt, int *scalerIncrement);


void newviewGTRGAMMA_AVX(int tipCase,
			 double *x1_start, double *x2_start, double *x3_start,
			 double *EV, double *tipVector,
			 unsigned char *tipX1, unsigned char *tipX2,
			 const int n, double *left, double *right, int *wgt, int *scalerIncrement
			 );

#endif

void reorder( double *x, int n, int span );
void reorder_back( double *x, int n, int span );

static int virtual_width( int n ) {
    const int global_vw = 2;
    return (n+1) / global_vw * global_vw;
}

pl_boolean modelExists(char *model, tree *tr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* AXML_H_ */


