/** 
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Derived from 
 * RAxML-HPC, a program for sequential and parallel estimation of phylogenetic
 * trees by Alexandros Stamatakis
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Tomas Flouri
 * Tomas.Flouri@h-its.org
 *
 * When publishing work that uses PLL please cite PLL
 *
 * ABSTRACT
 * 
 * PLL is a highly optimized, parallelized software library to ease the
 * development of new software tools dealing with phylogenetic inference. Among
 * the functions included in PLL are 
 *
 * DOCUMENTATION
 *
 * Extensive documentation for using PLL is available online at
 * 
 *                 http://www.libpll.org
 *
 *
 * USAGE
 *
 * To use PLL, 
 *
 * @file pll.h
 * @brief Data structures for tree and model 
 *
 * @author Tomas Flouri
 * @author Fernando Izquierdo-Carrasco
 * @author Andre Aberer
 * @author Alexandros Stamatakis
 */

#ifndef __pll__
#define __pll__

#include <stdint.h>
#include <stdio.h>
#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __MIC_NATIVE
#define PLL_BYTE_ALIGNMENT 64
#define PLL_VECTOR_WIDTH 8
#elif defined (__AVX)

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <immintrin.h>
#include <pmmintrin.h>
#endif

#define PLL_BYTE_ALIGNMENT 32
#define PLL_VECTOR_WIDTH 4

#elif defined (__SSE3)

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

#define PLL_BYTE_ALIGNMENT 16
#define PLL_VECTOR_WIDTH 2

#else
#define PLL_BYTE_ALIGNMENT 8
#define PLL_VECTOR_WIDTH 1
#endif

#ifdef _MSC_VER
	#define PLL_ALIGN_BEGIN __declspec(align(PLL_BYTE_ALIGNMENT))
	#define PLL_ALIGN_END
#else
	#define PLL_ALIGN_BEGIN
	#define PLL_ALIGN_END __attribute__((aligned(PLL_BYTE_ALIGNMENT)))
#endif


#include "stack.h"
#include "newick.h"
#include "queue.h"

#define PLL_MAX_TIP_EV                          0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define PLL_MAX_LOCAL_SMOOTHING_ITERATIONS      32          /** @brief maximum iterations of smoothings per insert in the */
#define PLL_ITERATIONS                          10          /* maximum iterations of iterations per insert */
#define PLL_NEWZPERCYCLE                        10           /* iterations of makenewz per tree traversal */
#define PLL_NMLNGTH                             256         /* number of characters in species name */
#define PLL_DELTAZ                              0.00001     /* test of net branch length change in update */
#define PLL_DEFAULTZ                            0.9         /* value of z assigned as starting point */
#define PLL_UNLIKELY                            -1.0E300    /* low likelihood for initialization */
#define PLL_SUMMARIZE_LENGTH                    -3
#define PLL_SUMMARIZE_LH                        -2
#define PLL_NO_BRANCHES                         -1
#define PLL_MASK_LENGTH                         32
#define PLL_ZMIN                                1.0E-15  /* max branch prop. to -log(PLL_ZMIN) (= 34) */
#define PLL_ZMAX                                (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */
#define PLL_TWOTOTHE256                         115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define PLL_MINLIKELIHOOD                       (1.0/PLL_TWOTOTHE256)
#define PLL_MINUSMINLIKELIHOOD                  -PLL_MINLIKELIHOOD


#define PLL_FORMAT_PHYLIP                       1 
#define PLL_FORMAT_FASTA                        2
#define PLL_FORMAT_NEWICK                       3

#define PLL_NNI_P_NEXT                          1       /**< Use p->next for the NNI move */
#define PLL_NNI_P_NEXTNEXT                      2       /**< Use p->next->next for the NNI move */

#define PLL_BADREAR                             -1

#define PLL_NUM_BRANCHES                        16

#define PLL_TRUE                                1
#define PLL_FALSE                               0

#define PLL_REARRANGE_SPR                       0
#define PLL_REARRANGE_TBR                       1
#define PLL_REARRANGE_NNI                       2

#define PLL_AA_SCALE                            10.0
#define PLL_AA_SCALE_PLUS_EPSILON               10.001

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

#define PLL_ALPHA_MIN                           0.02
#define PLL_ALPHA_MAX                           1000.0

#define PLL_RATE_MIN                            0.0000001
#define PLL_RATE_MAX                            1000000.0

#define PLL_LG4X_RATE_MIN                       0.0000001
#define PLL_LG4X_RATE_MAX                       1000.0

#define PLL_FREQ_MIN                            0.001

#define PLL_NUM_AA_STATES                       20
#define PLL_NUM_DNA_STATES                      4

/* 
   previous values between 0.001 and 0.000001

   TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA 
   previous value of FREQ_MIN was: 0.000001, but this seemed to cause problems with some 
   of the 7-state secondary structure models with some rather exotic small toy test datasets,
   on the other hand 0.001 caused problems with some of the 16-state secondary structure models

   For some reason the frequency settings seem to be repeatedly causing numerical problems
*/

#define PLL_ITMAX                               100    /* max number of iterations in brent's algorithm */

#define PLL_SHFT(a,b,c,d)                       (a)=(b);(b)=(c);(c)=(d);
#define PLL_SIGN(a,b)                           ((b) > 0.0 ? fabs(a) : -fabs(a))
#define PLL_ABS(x)                              (((x)<0)   ?  (-(x)) : (x))
#define PLL_MIN(x,y)                            (((x)<(y)) ?    (x)  : (y))
#define PLL_MAX(x,y)                            (((x)>(y)) ?    (x)  : (y))
#define PLL_SWAP(x,y)                           do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)
#define PLL_SWAP_PTR(x,y) do{ char* _t = x; x = y; y = _t; } while(0)
#define PLL_SWAP_INT(x,y) do{ int _t = x; x = y; y = _t; } while(0)

#define PLL_POINT_GAMMA(prob,alpha,beta)        PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#define PLL_LIB_NAME                            "PLL"
#define PLL_LIB_VERSION                         "1.0.1"
#define PLL_LIB_DATE                            "November 3 2014"

/* aminoacid substitution models */
#define PLL_DAYHOFF                             0
#define PLL_DCMUT                               1
#define PLL_JTT                                 2
#define PLL_MTREV                               3
#define PLL_WAG                                 4
#define PLL_RTREV                               5
#define PLL_CPREV                               6
#define PLL_VT                                  7
#define PLL_BLOSUM62                            8
#define PLL_MTMAM                               9
#define PLL_LG                                  10
#define PLL_MTART                               11
#define PLL_MTZOA                               12
#define PLL_PMB                                 13
#define PLL_HIVB                                14
#define PLL_HIVW                                15
#define PLL_JTTDCMUT                            16
#define PLL_FLU                                 17 
#define PLL_AUTO                                18
#define PLL_LG4M                                 19
#define PLL_LG4X                                20
#define PLL_GTR                                 21  /* GTR always needs to be the last one */
#define PLL_NUM_PROT_MODELS                     22

/* information criteria for auto protein model selection */
#define PLL_AUTO_ML   0
#define PLL_AUTO_BIC  1
#define PLL_AUTO_AIC  2
#define PLL_AUTO_AICC 3

/* bipartition stuff */
#define PLL_BIPARTITIONS_RF                     4

/* scenarios for likelihood computation */
#define PLL_TIP_TIP                             0
#define PLL_TIP_INNER                           1
#define PLL_INNER_INNER                         2


/* available data types in PLL */
#define PLL_MIN_MODEL                          -1
#define PLL_BINARY_DATA                         0
#define PLL_DNA_DATA                            1
#define PLL_AA_DATA                             2
#define PLL_SECONDARY_DATA                      3
#define PLL_SECONDARY_DATA_6                    4
#define PLL_SECONDARY_DATA_7                    5
#define PLL_GENERIC_32                          6
#define PLL_GENERIC_64                          7
#define PLL_MAX_MODEL                           8

#define PLL_SEC_6_A                             0
#define PLL_SEC_6_B                             1
#define PLL_SEC_6_C                             2
#define PLL_SEC_6_D                             3
#define PLL_SEC_6_E                             4

#define PLL_SEC_7_A                             5
#define PLL_SEC_7_B                             6
#define PLL_SEC_7_C                             7
#define PLL_SEC_7_D                             8
#define PLL_SEC_7_E                             9
#define PLL_SEC_7_F                             10

#define PLL_SEC_16                              11
#define PLL_SEC_16_A                            12
#define PLL_SEC_16_B                            13
#define PLL_SEC_16_C                            14
#define PLL_SEC_16_D                            15
#define PLL_SEC_16_E                            16
#define PLL_SEC_16_F                            17
#define PLL_SEC_16_I                            18
#define PLL_SEC_16_J                            19
#define PLL_SEC_16_K                            20

#define PLL_ORDERED_MULTI_STATE                 0
#define PLL_MK_MULTI_STATE                      1
#define PLL_GTR_MULTI_STATE                     2


/* available models of rate heterogeneity in PLL */
#define PLL_CAT                                 0
#define PLL_GAMMA                               1

/* recomp */
#define PLL_SLOT_UNUSED                        -2  /* value to mark an available vector */
#define PLL_NODE_UNPINNED                      -3  /* marks an inner node as not available in RAM */
#define PLL_INNER_NODE_INIT_STLEN              -1  /* initialization */

#define PLL_MIN_RECOM_FRACTION     0.1 /* at least this % of inner nodes will be allocated in RAM */
#define PLL_MAX_RECOM_FRACTION     1.0 /* always 1, just there for boundary checks */


typedef  int pllBoolean;

/* @brief PLL instance attribute structure */
typedef struct
{
  int rateHetModel;
  int fastScaling;
  int saveMemory;
  int useRecom;
  long randomNumberSeed;
  int numberOfThreads;
} pllInstanceAttr;

/** @brief Stores the recomputation-state of likelihood vectors  */
typedef struct
{
  int numVectors;      /**< Number of inner vectors allocated in RAM*/
  int *iVector;        /**< size: numVectors, stores node id || PLL_SLOT_UNUSED  */
  int *iNode;          /**< size: inner nodes, stores slot id || PLL_NODE_UNPINNED */
  int *stlen;          /**< Number of tips behind the current orientation of the indexed inner node (subtree size/cost) */ 
  int *unpinnable;     /**< size:numVectors , TRUE if we dont need the vector */
  int maxVectorsUsed;  
  pllBoolean allSlotsBusy; /**< on if all slots contain an ancesctral node (the usual case after first full traversal) */ 
} recompVectors;
/* E recomp */

/** @brief ???
 * @todo add explanation, is this ever used?  */
 
typedef unsigned int hashNumberType;



/*typedef uint_fast32_t parsimonyNumber;*/

#define PLL_PCF 32

/** @brief ???Hash tables 
 * @todo add explanation of all hash tables  */
typedef struct pllBipartitionEntry
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct pllBipartitionEntry *next;
} pllBipartitionEntry;

//typedef struct
//{
//  hashNumberType tableSize;
//  entry **table;
//  hashNumberType entryCount;
//}
//  hashtable;
//struct stringEnt
//{
//  int nodeNumber;
//  char *word;
//  struct stringEnt *next;
//};
//
//typedef struct stringEnt stringEntry;
//typedef struct
//{
//  hashNumberType tableSize;
//  stringEntry **table;
//}
//  stringHashtable;

typedef struct pllHashItem
{
  void * data;
  char * str;
  struct pllHashItem * next;
} pllHashItem;

typedef struct pllHashTable
{
  unsigned int size;
  struct pllHashItem ** Items;
  unsigned int entries;
} pllHashTable;




/** @brief Per-site Rate category entry: likelihood per-site and CAT rate applied ???
  *
  */
typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}rateCategorize;

/** @brief Traversal descriptor entry.
  * 
  * Contains the information required to execute an operation in a step of the tree traversal.
  * q   r
  *  \ /
  *   p
  *
  * The entry defines 2 input/parent nodes (q and r) and one output/child node (p)
  * qz represents the branch length(s) of the branch connecting q and p
  * rz represents the branch length(s) of the branch connecting r and p
  * PLL_TIP_TIP     Both p and r are tips
  * PLL_INNER_INNER Both p and r are inner nodes
  * @note PLL_TIP_INNER   q is a tip and r is an inner node (by convention, flip q and r if required)
  */
typedef struct
{
  int tipCase;                  /**< Type of entry, must be PLL_TIP_TIP PLL_TIP_INNER or PLL_INNER_INNER */
  int pNumber;                  /**< should exist in some nodeptr p->number */
  int qNumber;                  /**< should exist in some nodeptr q->number */
  int rNumber;                  /**< should exist in some nodeptr r->number */
  double qz[PLL_NUM_BRANCHES];
  double rz[PLL_NUM_BRANCHES];
  /* recom */
  int slot_p;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node p, otherwise unused */
  int slot_q;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node q, otherwise unused */
  int slot_r;                   /**< In recomputation mode, the RAM slot index for likelihood vector of node r, otherwise unused */
  /* E recom */
} traversalInfo;

/** @brief Traversal descriptor.
  * 
  * Describes the state of a traversal descriptor
  */
typedef struct
{
  traversalInfo *ti;              /**< list of traversal steps */
  int count;                      /**< number of traversal steps */
  int functionType;
  pllBoolean traversalHasChanged;   
  pllBoolean *executeModel;           
  double  *parameterValues;
} traversalData;

/** @brief Node record structure
  * 
  * Each inner node is a trifurcation in the tree represented as a circular list containing 3 node records. One node record uniquely identifies a subtree, and the orientation of the likelihood vector within a node
  *
  * p1 -------> p2 ----> to the next node
  * ^           |
  * |-----p3<---|          
  * 
  */
struct noderec;

/** @brief Branch length information.
  * 
  * @todo add relevant info on where this is used ???
  */
typedef struct
{
  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;





/** @brief Linkage of partitions.
  * 
  * @todo add relevant info on where this is used ???
  */
typedef struct
{
  pllBoolean valid;
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



  /** 
   *
   * the data structure below is fundamental for representing trees 
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
     we'd just set the first index of the branch length array z[PLL_NUM_BRANCHES].

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

     @todo I think we should rename the back pointer. It's not back, it can be forward depending on the orientation. We should renmae it to outer. Back is too confusing, I would assume it's the opposite of next, i.e. previous.

     @struct noderec

     @brief Tree node record

     A node in a tree is a structure which contains a cyclic list of pointers to 3 nodes which we call a \e roundabout. The first node is the structure itself, and the other two nodes are accessed via \a noderec->next and \a noderec->next->next. To access the outer node with which each of the 3 nodes forms an edge one has to use the \a back pointer

     @var noderec::next
     @brief Next node in the roundabout

     @var noderec::back
     @brief Outer node

     @var noderec::number
     @brief Node identifier

     In general, tips (i.e. leaves) are numbered from 1 to \e n where \e n is the number of taxa. Identifiers for internal nodes start from \e n + 1. Note
     that for a given inner node, the identifier must be the same for all 3 nodes that compose it.

     @var info::z
     @brief The branch lengths per partition for the main node in the roundabout

     @todo Append an image
  */
typedef  struct noderec
{
 
  branchInfo      *bInf;
  double           z[PLL_NUM_BRANCHES];
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

typedef unsigned int parsimonyNumber;

/* @brief Alignment, transition model, model of rate heterogenety and likelihood vectors for one partition.
  * 
  * @todo De-couple into smaller data structures
  *
  * ALIGNMENT DATA 
  * This depends only on the type of data in this partition of the alignment 
  *
  * MODEL OF RATE HETEROGENETY, We use either GAMMA or PSR 
  * Rate heterogenety: Per Site Categories (PSR) model aka CAT, 
  * Rate of site i is given by perSiteRates[rateCategory[i]]
  *
  * TRANSITION MODEL: We always assume General Time Reversibility 
  * Transistion probability matrix: P(t) = exp(Qt)
  * Branch length t is the expected number of substitutions per site 
  * Pij(t) is the probability of going from state i to state j in a branch of length t 
  * Relative substitution rates (Entries in the Q matrix) 
  * In GTR we can write Q = S * D, where S is a symmetrical matrix and D a diagonal with the state frequencies 

    @var protModels
    @brief Protein models

    @detail Detailed protein models descriptiopn

    @var autoProtModels
    @brief Auto prot models
    @detail Detailed auto prot models
  */
 


/** @struct pInfo
    
    @brief Partition information structure

    This data structure encapsulates all properties and auxiliary variables that together
    consist a partition.

    @var pInfo::dataType
    @brief Type of data this partition contains

    Can be DNA (\b PLL_DNA_DATA) or AminoAcid (\b PLL_AA_DATA) data

    @var pInfo::states
    @brief Number of states

    Number of states this type of data can consist of

    @var pInfo::maxTipStates
    @brief Number of undetermined states (possible states at the tips)

    This is the total number of possible states that can appear in the alignment. This includes degenerate (undetermined) bases

    @var pInfo::partitionName
    @brief Name of partition

    A null-terminated string describing the name of partition

    @var pInfo::lower
    @brief Position of the first site in the alignment that is part of this partition [1, tr->originalCrunchedLength]

    @var pInfo::upper
    @brief Position of the last site that is part of this partition plus one (i.e. position of the first site that is not part of this partition) 

    @var pInfo::width
    @brief Number of sites in the partition (i.e. \a upper - \a lower)

    @var pInfo::wgt
    @brief Weight of site

    Number of times this particular site appeared in the partition before the duplicates were removed and replaced by this weight

    @var pInfo::empiricalFrequencies
    @brief Empirical frequency of each state in the current partition

    @var pInfo::perSiteRates
    @brief Per Site Categories (PSR) or aka CAT values for each rate

    @var pInfo::rateCategory
    @brief CAT category index for each site

    @var pInfo::numberOfCategories
    @brief CAT size of the set of possible categories

    @var pInfo::alpha
    @brief Gamma parameter to be optimized
    
    @var pInfo::gammaRates
    @brief Values of the 4 gamma categories (rates) computed given an alpha

    @var pInfo::substRates
    @brief Entries of substitution matrix, e.g. 6 free parameters in DNA

    In GTR we can write \f$ Q = S * D \f$, where \f$ S \f$ is a symmetrical matrix and \f$ D \f$ a diagonal with the state frequencies,
    which is represented by the array \a frequencies. The symmetrical matrix is the array \a substRates

    @var pInfo::frequencies
    @brief State frequencies, entries in D are initialized as empiricalFrequencies
    
    In GTR we can write \f$ Q = S * D \f$, where \f$ S \f$ is a symmetrical matrix and \f$ D \f$ a diagonal with the state frequencies,
    which is represented by the array \a frequencies. The symmetrical matrix is the array \a substRates

    @var pInfo::freqExponents

    @var pInfo::EIGN
    @brief Eigenvalues of Q matrix

    @var pInfo::EV
    @brief Eigenvectors of Q matrix

    @var pInfo::EI
    @brief Inverse eigenvectors of Q matrix

    @var pInfo::left
    @brief P matrix for the left term of the conditional likelihood equation

    @var pInfo::right
    @brief P matrix for the right term of the conditional likelihood equation

    @var pInfo::tipVector
    @brief Precomputed (based on current P matrix) conditional likelihood vectors for every possible base 

    @var pInfo::EIGN_LG4
    @brief Eigenvalues of Q matrix for the LG4 model

    @var pInfo::EV_LG4
    @brief Eigenvectors of Q matrix for the LG4 model

    @var pInfo::EI_LG4
    @brief Inverse eigenvectors of Q matrix for the LG4 model
    
    @var pInfo::frequencies_LG4
    @brief State frequencies for the LG4 model

    @var pInfo::tipVector_LG4
    @brief Precomputed (based on current P matrix) conditional likelihood vectors for every possible base for the LG4 model

    @var pInfo::substRates_LG4
    @brief Entries of substitution matrix for the LG4 model

    @var pInfo::protModels
    @brief Protein model for current partition

    In case \a pInfo::dataType is set to \a PLL_AA_DATA then \a protModels indicates the index in the global array \a protModels
    of the protein model that the current partition uses.

    @var pInfo::autoProtModels
    @brief Best fitted protein model for the \b PLL_AUTO partitions

    If \a protModels is set to \b PLL_AUTO then \a autoProtModels holds the currently detected best fitting protein model for the partition

    @var pInfo::protUseEmpiricalFreqs

    @var pInfo::nonGTR

    @var pInfo::optimizeBaseFrequencies

    @var pInfo::optimizeAlphaParameter

    @var pInfo::optimizeSubstitutionRates

    @var pInfo::symmetryVector

    @var pInfo::frequencyGrouping


    @todo
      Document freqExponents

*/



typedef struct {
  int     dataType;
  int     states;
  int     maxTipStates;
  char   *partitionName;
  int     lower;
  int     upper;
  int     width;
  int    *wgt;
  double *empiricalFrequencies; 


  /* MODEL OF RATE HETEROGENETY, We use either GAMMA or PSR */
  /* Rate heterogenety: Per Site Categories (PSR) model aka CAT, see updatePerSiteRates() */
  /* Rate of site i is given by perSiteRates[rateCategory[i]] */
  double *perSiteRates;
  int    *rateCategory;
  int     numberOfCategories;
  /* Rate heterogenety: GAMMA model of rate heterogenety */
  double alpha;
  double *gammaRates;


  /* TRANSITION MODEL: We always assume General Time Reversibility */
  /* Transistion probability matrix: P(t) = exp(Qt)*/
  /* Branch length t is the expected number of substitutions per site */
  /* Pij(t) is the probability of going from state i to state j in a branch of length t */
  /* Relative substitution rates (Entries in the Q matrix) */
  /* In GTR we can write Q = S * D, where S is a symmetrical matrix and D a diagonal with the state frequencies */
  double *substRates;       /**< TRANSITION MODEL Entries in S, e.g. 6 free parameters in DNA */   
  double *frequencies;      /**< State frequencies, entries in D, are initialized as empiricalFrequencies */
  double *freqExponents;
  /* Matrix decomposition: @todo map this syntax to Explanation of the mathematical background */
  double *EIGN;
  double *EV;
  double *EI;
  double *left;
  double *right;
  double *tipVector;


  /* asc bias */
  pllBoolean       ascBias;
  int           ascOffset;
  int         * ascExpVector;
  double      * ascSumBuffer;
  double      * ascVector;
  double        ascScaler[64];
  
  /* LG4 */

  double *EIGN_LG4[4];
  double *EV_LG4[4];
  double *EI_LG4[4];

  double *frequencies_LG4[4];
  double *tipVector_LG4[4];
  double *substRates_LG4[4];
  
  /* LG4X */

  double lg4x_weights[4];
  double lg4x_weightExponents[4];
  double lg4x_weightsBuffer[4];
  double lg4x_weightExponentsBuffer[4];
  double lg4x_weightLikelihood;
  
  /* Protein specific */
  int     protModels;			/**< Empirical model matrix */
  int     autoProtModels;		/**< Model selected with "auto" protein model */
  int     protUseEmpiricalFreqs;	/**< Whether to use empirical frequencies for protein model */

  pllBoolean nonGTR;
  pllBoolean optimizeBaseFrequencies;	/**< Whether to optimize base frequencies */
  pllBoolean optimizeAlphaParameter;	/**< Whether to optimize alpha parameters and gamma rates */
  pllBoolean optimizeSubstitutionRates;	/**< Whether to optimize substitution rates */
  int    *symmetryVector;		/**< Specify linkage between substitution rate parameters */
  int    *frequencyGrouping;

  /* LIKELIHOOD VECTORS */

  /* partial LH Inner vectors  ancestral vectors, we have 2*tips - 3 inner nodes */
  double          **xVector;          /**< Conditional likelihood vectors for inner nodes */
  unsigned char   **yVector;          /**< Tip entries (sequence) for tip nodes */
  unsigned int     *globalScaler;     /**< Counters for scaling operations done at node i */

  /* data structures for conducting per-site likelihood scaling.
     this allows to compute the per-site log likelihood scores 
     needed for RELL-based bootstrapping and all sorts of statistical 
     tests for comparing trees ! */
  int              **expVector;     /**< @brief An entry per inner node. Each element is an array of size the number of sites in the current partition and represents how many times the respective site has been scaled in the subtree rooted at the current node */
  size_t           *expSpaceVector; /**< @brief Each entry represents an inner node and states the size of the corresponding element in \a expVector, which is the number of sites for the current partition */

  /* These are for the saveMemory option (tracking gaps to skip computations and memory) */
  size_t           *xSpaceVector;       /* Size of conditional likelihood vectors per inner node */
  int               gapVectorLength;    /** Length of \a gapVector bitvector in unsigned integers assuming that \a unsigned \a int is 32bits. It is set to partition size / 32 */
  unsigned int     *gapVector;          /** A bit vector of size \a gapVectorLength * 32 bits. A bit is set to 1 if the corresponding */
  double           *gapColumn; 

  /* Parsimony vectors at each node */
  size_t parsimonyLength;
  parsimonyNumber *parsVect; 
  parsimonyNumber *perSiteParsScores;

  /* This buffer of size width is used to store intermediate values for the branch length optimization under 
     newton-raphson. The data in here can be re-used for all iterations irrespective of the branch length.
   */
  double *sumBuffer; 

  /* Buffer to store the per-site log likelihoods */
  double *perSiteLikelihoods;

  /* This buffer of size width is used to store the ancestral state at a node of the tree. */
  double *ancestralBuffer;

  /* From tree */
  pllBoolean executeModel;
  double fracchange;
  double rawFracchange;
  double partitionContribution;
  double partitionWeight;
  double partitionLH;

// #if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
  int partitionAssignment;
// #endif

} pInfo;

typedef struct
 {
   pInfo **partitionData;
   int numberOfPartitions;
   pllBoolean perGeneBranchLengths;
   pllBoolean dirty;
   linkageList *alphaList;
   linkageList *rateList;
   linkageList *freqList;
 }  partitionList;



#define PLL_REARR_SETTING 1
#define PLL_FAST_SPRS     2
#define PLL_SLOW_SPRS     3


/** @brief Checkpointing states. 
 * 
 * @todo Raxml specific 
  */
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
  pllBoolean impr;
  pllBoolean cutoff;  
       
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


/** @brief Tree topology.
 * 
 * @todo Apart from the topology this structure contains several fields that act like global variables in raxml
  */
typedef  struct  {

  int *ti;

  /* recomp */
  recompVectors *rvec;            /**< this data structure tracks which vectors store which nodes */
  float maxMegabytesMemory;       /**< User says how many MB in main memory should be used */
  float vectorRecomFraction;      /**< vectorRecomFraction ~= 0.8 * maxMegabytesMemory  */
  pllBoolean useRecom;               /**< ON if we apply recomputation of ancestral vectors*/
#ifdef _DEBUG_RECOMPUTATION 
  traversalCounter *travCounter;
  double stlenTime;
#endif
  /* E recomp */
  
  pllBoolean fastScaling;
  pllBoolean saveMemory;
  int              startingTree;
  long             randomNumberSeed;

  double          *lhs;         /**< Array to store per-site log likelihoods of \a originalCrunchedLength (compressed) sites */
  double          *patrat;      /**< rates per pattern */
  double          *patratStored; 
  int             *rateCategory;
  int             *aliaswgt;    /**< weight by pattern */ 
  pllBoolean    manyPartitions;

  pllBoolean grouped;              /**< No idea what this is, but is always set to PLL_FALSE */
  pllBoolean constrained;          /**< No idea what this is, but is always set to PLL_FALSE */
  int threadID;
  volatile int numberOfThreads;

//#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
 
  unsigned char *y_ptr; 
  
  double lower_spacing;
  double upper_spacing; 

  double *ancestralVector;

//#endif
  
  pllHashTable     *nameHash;
  char           ** tipNames;

  char             *secondaryStructureInput;

  traversalData    td[1];

  int              maxCategories;
  int              categories;

  double           coreLZ[PLL_NUM_BRANCHES];
  
 
  branchInfo       *bInf;

  int              multiStateModel;


  pllBoolean curvatOK[PLL_NUM_BRANCHES];

  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */

  /* model stuff end */
  unsigned char    **yVector;        /**< list of raw sequences (parsed from the alignment)*/

  int              secondaryStructureModel;
  int              originalCrunchedLength; /**< Length of alignment after removing duplicate sites in each partition */

  int              *secondaryStructurePairs;

  double            fracchange;      /**< Average substitution rate */
  double            rawFracchange;
  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;

  double           startLH;
  double           endLH;
  double           likelihood;           /**< last likelihood value evaluated for the current topology */
 
  node           **nodep;                /**< pointer to the list of nodes, which describe the current topology */
  nodeptr          nodeBaseAddress;
  node            *start;                /**< starting node by default for full traversals (must be a tip contained in the tree we are operating on) */
  int              mxtips;  /**< Number of tips in the topology */

  int              *constraintVector;   /**< @todo What is this? */
  int              numberOfSecondaryColumns;
  pllBoolean          searchConvergenceCriterion;
  int              ntips;
  int              nextnode;  

  pllBoolean          bigCutoff;
  pllBoolean          partitionSmoothed[PLL_NUM_BRANCHES];
  pllBoolean          partitionConverged[PLL_NUM_BRANCHES];
  pllBoolean          rooted;
  pllBoolean          doCutoff;
 
  double         gapyness;

  char **nameList;     /**< list of tips names (read from the phylip file) */
  char *tree_string;   /**< the newick representaion of the topology */
  char *tree0;
  char *tree1;
  int treeStringLength;
 
  unsigned int bestParsimony;
  unsigned int *parsimonyScore;
  
  double bestOfNode;
  nodeptr removeNode;   /**< the node that has been removed. Together with \a insertNode represents an SPR move */
  nodeptr insertNode;   /**< the node where insertion should take place . Together with \a removeNode represents an SPR move*/

  double zqr[PLL_NUM_BRANCHES];
  double currentZQR[PLL_NUM_BRANCHES];

  double currentLZR[PLL_NUM_BRANCHES];
  double currentLZQ[PLL_NUM_BRANCHES];
  double currentLZS[PLL_NUM_BRANCHES];
  double currentLZI[PLL_NUM_BRANCHES];
  double lzs[PLL_NUM_BRANCHES];
  double lzq[PLL_NUM_BRANCHES];
  double lzr[PLL_NUM_BRANCHES];
  double lzi[PLL_NUM_BRANCHES];


  unsigned int **bitVectors;

  unsigned int vLength;

  pllHashTable *h;                 /**< hashtable for ML convergence criterion */
  //hashtable *h;
 
  int optimizeRateCategoryInvocations;

  checkPointState ckp;
  pllBoolean thoroughInsertion; /**< true if the neighbor branches should be optimized when a subtree is inserted (slower)*/
  pllBoolean useMedian;

  int autoProteinSelectionType;

  pllStack * rearrangeHistory;


  /* analdef defines */
  /* TODO: Do some initialization */
  int              bestTrav;            /**< best rearrangement radius */
  int              max_rearrange;       /**< max. rearrangemenent radius */
  int              stepwidth;           /**< step in rearrangement radius */
  int              initial;             /**< user defined rearrangement radius which also sets bestTrav if initialSet is set */
  pllBoolean          initialSet;          /**< set bestTrav according to initial */
  int              mode;                /**< candidate for removal */
  pllBoolean        perGeneBranchLengths;
  pllBoolean        permuteTreeoptimize;   /**< randomly select subtrees for SPR moves */
  pllBoolean        compressPatterns;
  double         likelihoodEpsilon;
  pllBoolean        useCheckpoint;

} pllInstance;

/** @brief Stores data related to a NNI move  */
typedef struct {
        pllInstance * tr;
        nodeptr p;
        int nniType;
        double z[PLL_NUM_BRANCHES]; // optimize branch lengths
        double z0[PLL_NUM_BRANCHES]; // unoptimized branch lengths
        double likelihood;
        double deltaLH;
} nniMove;

/***************************************************************/

typedef struct {
  int partitionNumber;
  int partitionLength;
} partitionType;

typedef struct
{
  double z[PLL_NUM_BRANCHES];
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

/** @brief Connection within a topology.
*   */
typedef struct conntyp {
    double           z[PLL_NUM_BRANCHES];           /**< branch length */
    node            *p, *q;       /**< parent and child sectors */
    void            *valptr;      /**< pointer to value of subtree */
    int              descend;     /**< pointer to first connect of child */
    int              sibling;     /**< next connect from same parent */
    } pllConnect, *connptr;

/** @brief Single Topology
*   */
typedef  struct {
    double           likelihood;
    int              initialTreeNumber;
    pllConnect         *links;       /**< pointer to first connect (start) */
    node            *start;
    int              nextlink;    /**< index of next available connect */
                                  /**< tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;    /**< next available inner node for tree parsing */
    int              scrNum;      /**< position in sorted list of scores */
    int              tplNum;      /**< position in sorted list of trees */
    } topol;

/** @brief small helper data structure for printing out/downstream use of marginal ancestral probability vectors.
*
* it is allocated as an array that has the same length as the input alignment and can be used to 
*   index the ancestral states for each position/site/pattern 
*   */
typedef struct {
  double *probs; /**< marginal ancestral states */
  char c; /**< most likely stated, i.e. max(probs[i]) above */
  int states; /**< number of states for this position */
} ancestralState;

/** @brief List of topologies
*
*   */
typedef struct {
    double           best;        /**< highest score saved */
    double           worst;       /**< lowest score saved */
    topol           *start;       /**< starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /**< maximum topologies to save */
    int              nvalid;      /**< number of topologies saved */
    int              ninit;       /**< number of topologies initialized */
    int              numtrees;    /**< number of alternatives tested */
    pllBoolean          improved;
    } bestlist;

/** @brief  This is used to look up some hard-coded data for each data type 
*   */
typedef struct 
{
  int leftLength;         /**< s^2 */
  int rightLength;/**< s^2 */
  int eignLength;/**<  s */
  int evLength;
  int eiLength;
  int substRatesLength;   /**< (s^2 - s)/2 free model parameters for matrix Q i.e. substitution rates */
  int frequenciesLength;  /**< s frequency of each state */ 
  int tipVectorLength;    /* ??? */
  int symmetryVectorLength;
  int frequencyGroupingLength;

  pllBoolean nonGTR;
  pllBoolean optimizeBaseFrequencies;

  int undetermined;

  const char *inverseMeaning;

  int states;   /* s */

  pllBoolean smoothFrequencies;

  const unsigned  int *bitVector;

} partitionLengths;

typedef struct
{
  int rearrangeType;
  double  likelihood;

  union {
    struct {
      double * zp;
      double * zpn;
      double * zpnn;
      double * zqr;
      nodeptr pn;
      nodeptr pnn;
      nodeptr r;
      nodeptr p;
      nodeptr q;
    } SPR;
    struct {
      nodeptr origin;
      int swapType;
      double z[PLL_NUM_BRANCHES];
    } NNI;
  };
} pllRollbackInfo;


/** @struct pllRearrangeAttr
 
    @brief Structure holding attributes for searching possible tree rearrangements
    
    Holds the attributes for performing tree rearrangements.

    @var pllRearrangeAttr
      The origin node where the search should start

    @var pllRearrangeAttr:mintrav
      The minimum radius around the origin node \a p for which nodes should be tested

    @var pllRearrangeAttr:maxtrav
      The maximum radius around the origin node \a p for which nodes should be tested

    @var pllRearrangeAttr:max
      Maximum number of results to be returned
*/
typedef struct
 {
   nodeptr p;
   int mintrav;
   int maxtrav;
 } pllRearrangeAttr;

/** @typedef pllRearrangeInfo
    
    @brief Tree rearrangement information structure

    Holds information for conducting tree arrangements. This structure
    is the result of a tree arrangement search under given search
    attributes.

    @var pllRearrangeInfo::rearrangeType
      Type of rearrangement. Can be \b PLL_REARRANGE_SPR, \b PLL_REARRANGE_NNI or
      \b PLL_REARRANGE_TBR
    
    @var pllRearrangeInfo::likelihood
      Holds the computed likelihood for the addressed rearrangement

    @var pllRearrangeInfo::SPR::removeNode
      Node where to perform subtree pruning

    @var pllRearrangeInfo::SPR::insertNode
      Node where to place the pruned subtree

    @var pllRearrangeInfo::zqr
      Holds the computed branch lengths after the SPR
*/
typedef struct
 {
   int rearrangeType;
   double  likelihood;
   union {
     struct {
       nodeptr removeNode;
       nodeptr insertNode;
       double  zqr[PLL_NUM_BRANCHES];
     } SPR;
     struct {
       nodeptr originNode;
       int     swapType;
     } NNI;
   };
 } pllRearrangeInfo;


typedef struct
 {
   int max_entries;
   int entries;
   pllRearrangeInfo * rearr;
 } pllRearrangeList;

/** @brief Generic structure for storing a multiple sequence alignment */
typedef struct
 {
   int              sequenceCount;      /**< @brief Number of sequences */
   int              sequenceLength;     /**< @brief Length of sequences */
   int              originalSeqLength;  /**< @brief Original length of sequences (not modified after removing duplicates) */
   char          ** sequenceLabels;     /**< @brief An array of where the \a i-th element is the name of the \a i-th sequence */
   unsigned char ** sequenceData;       /**< @brief The actual sequence data */
   int            * siteWeights;        /**< @brief An array where the \a i-th element indicates how many times site \a i appeared (prior to duplicates removal) in the alignment */
 } pllAlignmentData;


/******************** START OF API FUNCTION DESCRIPTIONS ********************/


#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI)) || defined(_IQTREE_MPI)
pllBoolean isThisMyPartition(partitionList* pr, int tid, int model);
void printParallelTimePerRegion(void);
#endif

#ifdef _FINE_GRAIN_MPI
extern void pllFinalizeMPI (void);
#endif



/**
 * @brief Create the main instance of PLL
 *   
 * Create an instance of the phylogenetic likelihood library
 *
 * @param rateHetModel   Rate heterogeneity model
 * @param fastScaling    TODO: explain fastScaling here
 * @param saveMemory     TODO: explain saveMemory here
 * @param useRecom       If set to \b PLL_TRUE, enables ancestral state recomputation
 * 
 * @todo                 Document fastScaling, rate heterogeneity and saveMemory and useRecom
 *
 * @note                 Do not set \a saveMemory to when using \a useRecom as memory saving 
 *                       techniques are not yet implemented for ancestral state recomputation. 
 * 
 * @return               On success returns an instance to PLL, otherwise \b NULL
 */
extern pllInstance * pllCreateInstance (pllInstanceAttr * pInst);

/** 
 *  @ingroup instanceLinkingGroup
 *  @brief Load alignment to the PLL instance
 *   
 *   Loads (copies) the parsed alignment \a alignmentData to the PLL instance
 *   as a deep copy.
 * 
 *    @param tr              The library instance
 *    @param alignmentData   The multiple sequence alignment
 *    @param partitions      List of partitions
 *
 *    @return Returns 1 in case of success, 0 otherwise.
 */
extern int pllLoadAlignment (pllInstance * tr, 
                             pllAlignmentData * alignmentData, 
                             partitionList * pList);

/**
 * @brief Compute the empirical base frequencies for all partitions
 *
 * Compute the empirical base frequencies for all partitions in the list \a pl.
 *
 * @param pl                Partition list
 * @param alignmentData     Multiple sequence alignment
 *
 * @return   A list of \a pl->numberOfPartitions arrays each of size
             \a pl->partitionData[i]->states, where \a i is the \a i-th partition
*/
extern double ** pllBaseFrequenciesAlignment (pllAlignmentData * alignmentData, partitionList * pl);
extern double ** pllBaseFrequenciesInstance (pllInstance * tr, partitionList * pl);

/* pthreads and MPI */
extern void pllStartPthreads (pllInstance *tr, partitionList *pr);
extern void pllStopPthreads (pllInstance * tr);
extern void pllLockMPI (pllInstance * tr);
extern void pllInitMPI(int * argc, char **argv[]);


/* handling branch lengths*/
extern double pllGetBranchLength (pllInstance *tr, nodeptr p, int partition_id);
extern void pllSetBranchLength (pllInstance *tr, nodeptr p, int partition_id, double bl);
extern int pllNniSearch(pllInstance * tr, partitionList *pr, int estimateModel);
extern void pllOptimizeBranchLengths ( pllInstance *tr, partitionList *pr, int maxSmoothIterations );


extern void pllEvaluateLikelihood (pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean fullTraversal, pllBoolean getPerSiteLikelihoods);
extern void pllUpdatePartials (pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean masked);
extern void pllUpdatePartialsAncestral(pllInstance *tr, partitionList *pr, nodeptr p);
extern void pllNewviewIterative(pllInstance *tr, partitionList *pr, int startIndex);
extern void pllEvaluateIterative(pllInstance *tr, partitionList *pr, pllBoolean getPerSiteLikelihoods);

/* newick parser declarations */
extern pllNewickTree * pllNewickParseString (const char * newick);
extern pllNewickTree * pllNewickParseFile (const char * filename);
extern int pllValidateNewick (pllNewickTree *);
extern void pllNewickParseDestroy (pllNewickTree **);
extern int pllNewickUnroot (pllNewickTree * t);
extern char * pllTreeToNewick ( char *treestr, pllInstance *tr, partitionList *pr, nodeptr p,
      pllBoolean printBranchLengths, pllBoolean printNames, pllBoolean printLikelihood,
      pllBoolean rellTree, pllBoolean finalPrint, int perGene,
      pllBoolean branchLabelSupport, pllBoolean printSHSupport);

/* partition parser declarations */
extern void  pllQueuePartitionsDestroy (pllQueue ** partitions);
extern pllQueue * pllPartitionParse (const char * filename);
extern pllQueue * pllPartitionParseString (const char * p);
extern void pllPartitionDump (pllQueue * partitions);
void pllBaseSubstitute (pllInstance * tr, partitionList * partitions);
//void pllBaseSubstitute (pllAlignmentData * tr, partitionList * partitions);
partitionList * pllPartitionsCommit (pllQueue * parts, pllAlignmentData * alignmentData);
int pllPartitionsValidate (pllQueue * parts, pllAlignmentData * alignmentData);
extern void pllAlignmentRemoveDups (pllAlignmentData * alignmentData, partitionList * pl);
void pllPartitionsDestroy (pllInstance *, partitionList **);

/* alignment data declarations */
extern void pllAlignmentDataDestroy (pllAlignmentData *);
extern int pllAlignmentDataDumpFile (pllAlignmentData *, int, const char *);
extern void pllAlignmentDataDumpConsole (pllAlignmentData * alignmentData);
extern pllAlignmentData * pllInitAlignmentData (int, int);
extern pllAlignmentData * pllParseAlignmentFile (int fileType, const char *);
extern pllAlignmentData *pllParsePHYLIPString (const char *rawdata, long filesize);


/* model management */
int pllInitModel (pllInstance *, partitionList *);
void pllInitReversibleGTR(pllInstance * tr, partitionList * pr, int model);
void pllMakeGammaCats(double alpha, double *gammaRates, int K, pllBoolean useMedian);
int pllLinkAlphaParameters(char *string, partitionList *pr);
int pllLinkFrequencies(char *string, partitionList *pr);
int pllLinkRates(char *string, partitionList *pr);
int pllSetSubstitutionRateMatrixSymmetries(char *string, partitionList * pr, int model);
void pllSetFixedAlpha(double alpha, int model, partitionList * pr, pllInstance *tr);
void pllSetFixedBaseFrequencies(double *f, int length, int model, partitionList * pr, pllInstance *tr);
int  pllSetOptimizeBaseFrequencies(int model, partitionList * pr, pllInstance *tr);
void pllSetSubstitutionMatrix(double *q, int length, int model, partitionList * pr,  pllInstance *tr);
void pllSetFixedSubstitutionMatrix(double *q, int length, int model, partitionList * pr,  pllInstance *tr);
int pllGetInstRateMatrix (partitionList * pr, int model, double * outBuffer);
int pllOptimizeModelParameters(pllInstance *tr, partitionList *pr, double likelihoodEpsilon);
double pllGetAlpha (partitionList * pr, int pid);
void pllGetGammaRates (partitionList * pr, int pid, double * outBuffer);
extern void pllGetBaseFrequencies(partitionList * pr, int model, double * outBuffer);
extern void pllGetSubstitutionMatrix (partitionList * pr, int model, double * outBuffer);
void pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models);
extern void pllOptRatesGeneric(pllInstance *tr, partitionList *pr, double modelEpsilon, linkageList *ll);
extern void pllOptBaseFreqs(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll);
extern void pllOptAlphasGeneric(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll);
extern void pllOptLG4X(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll, int numberOfModels);

/* tree topology */
void pllTreeInitTopologyNewick (pllInstance *, pllNewickTree *, int);
void pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList);
void pllTreeInitTopologyForAlignment (pllInstance * tr, pllAlignmentData * alignmentData);
extern void pllMakeRandomTree ( pllInstance *tr);
void pllMakeParsimonyTree(pllInstance *tr);
extern void pllMakeParsimonyTreeFast(pllInstance *tr, partitionList *pr, int sprDist);
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions, int sprDist);
nodeptr pllGetRandomSubtree(pllInstance *);
extern void pllFreeParsimonyDataStructures(pllInstance *tr, partitionList *pr);
void pllDestroyInstance (pllInstance *);
extern void pllGetAncestralState(pllInstance *tr, partitionList *pr, nodeptr p, double * outProbs, char * outSequence);
unsigned int pllEvaluateParsimony(pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean full, pllBoolean perSiteScores);
void pllInitParsimonyStructures(pllInstance *tr, partitionList *pr, pllBoolean perSiteScores);

/* rearrange functions (NNI and SPR) */
pllRearrangeList * pllCreateRearrangeList (int max);
void pllDestroyRearrangeList (pllRearrangeList ** bestList);
void pllRearrangeSearch (pllInstance * tr, partitionList * pr, int rearrangeType, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList);
void pllRearrangeCommit (pllInstance * tr, partitionList * pr, pllRearrangeInfo * rearr, int saveRollbackInfo);
int pllRearrangeRollback (pllInstance * tr, partitionList * pr);
void pllClearRearrangeHistory (pllInstance * tr);
int pllRaxmlSearchAlgorithm (pllInstance * tr, partitionList * pr, pllBoolean estimateModel);
int pllGetTransitionMatrix (pllInstance * tr, partitionList * pr, nodeptr p, int model, int rate, double * outBuffer);
void pllGetTransitionMatrix2 (pllInstance * tr, partitionList * pr, int model, nodeptr p, double * outBuffer);
int pllGetCLV (pllInstance * tr, partitionList * pr, nodeptr p, int partition, double * outProbs);
extern int pllTopologyPerformNNI(pllInstance * tr, nodeptr p, int swap);

/* hash functions */
unsigned int pllHashString (const char * s, unsigned int size);
int pllHashAdd  (pllHashTable * hTable, unsigned int hash, const char * s, void * item);
pllHashTable * pllHashInit (unsigned int n);
int pllHashSearch (struct pllHashTable * hTable, char * s, void ** item);
void pllHashDestroy (struct pllHashTable ** hTable, void (*cbDealloc)(void *));

/* node specific functions */
nodeptr pllGetOrientedNodePointer (pllInstance * pInst, nodeptr p);

/* other functions */
extern char * pllReadFile (const char *, long *);
extern int * pllssort1main (char ** x, int n);
extern node ** pllGetInnerBranchEndPoints (pllInstance * tr);

/* ---------------- */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
