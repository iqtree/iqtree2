/*
 * pllInternal.h
 *
 *  Created on: Feb 17, 2014
 *      Author: diego
 */

#ifndef PLLINTERNAL_H_
#define PLLINTERNAL_H_

#include "pll.h"
#include "genericParallelization.h"
#include "errcodes.h"
#include "lexer.h"
#include "parsePartition.h"
#include "mem_alloc.h"

//extern int lookupWord(char *s, stringHashtable *h);

extern void getDataTypeString(pllInstance *tr, pInfo *partitionInfo, char typeOfData[1024]);
extern int countTips(nodeptr p, int numsp);
extern unsigned int precomputed16_bitcount(unsigned int n, char *bits_in_16bits);

extern size_t discreteRateCategories(int rateHetModel);

extern const partitionLengths * getPartitionLengths(pInfo *p);
extern pllBoolean getSmoothFreqs(int dataType);
extern const unsigned int *getBitVector(int dataType);
extern int getUndetermined(int dataType);
extern int getStates(int dataType);
extern char getInverseMeaning(int dataType, unsigned char state);
extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );

extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupFull ( nodeptr p, nodeptr q, double *z);
extern void hookupDefault ( nodeptr p, nodeptr q);
extern pllBoolean whitechar ( int ch );
extern void printLog ( pllInstance *tr);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void initModel ( pllInstance *tr, double **empiricalFrequencies, partitionList * partitions);

extern void resetBranches ( pllInstance *tr );
extern void modOpt ( pllInstance *tr, partitionList *pr, double likelihoodEpsilon);

extern void initializePartitionData(pllInstance *localTree, partitionList * localPartitions);
extern void initMemorySavingAndRecom(pllInstance *tr, partitionList *pr);

extern void nodeRectifier ( pllInstance *tr );
extern void allocateParsimonyDataStructures(pllInstance *tr, partitionList *pr);

extern FILE *myfopen(const char *path, const char *mode);

extern pllBoolean initrav ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void initravPartition ( pllInstance *tr, nodeptr p, int model );
extern void update ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void smooth ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void smoothTree ( pllInstance *tr, partitionList *pr, int maxtimes );
extern void localSmooth ( pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes );
extern pllBoolean localSmoothMulti(pllInstance *tr, nodeptr p, int maxtimes, int model);

extern void smoothRegion ( pllInstance *tr, partitionList *pr, nodeptr p, int region );
extern void regionalSmooth ( pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( pllInstance *tr, partitionList *pr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p );
extern pllBoolean insertBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q);
extern pllBoolean insertRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern pllBoolean testInsertBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern int NNI(pllInstance * tr, nodeptr p, int swap);
extern void addTraverseBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern pllBoolean testInsertRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( pllInstance *tr, partitionList *pr );

extern void initTL ( topolRELL_LIST *rl, pllInstance *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, pllInstance *tr, int n, int numBranches );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, pllInstance *tr, int index );

extern topol  *setupTopol (int maxtips);
extern void saveTree (pllInstance *tr, topol *tpl, int numBranches);
extern pllBoolean restoreTree (topol *tpl, pllInstance *tr, partitionList *pr);




extern int  saveBestTree (bestlist *bt, pllInstance *tr, int numBranches);
extern int  recallBestTree (bestlist *bt, int rank, pllInstance *tr, partitionList *pr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern pllBoolean freeBestTree ( bestlist *bt );


/* extern int treeReadLen (FILE *fp, pllInstance *tr, pllBoolean readBranches, pllBoolean readNodeLabels, pllBoolean topologyOnly);
extern void getStartingTree (pllInstance *tr); 
extern void treeReadTopologyString(char *treeString, pllInstance *tr);
extern double treeLength (pllInstance *tr, int model);*/
extern double evaluatePartialGeneric (pllInstance *, partitionList *pr, int i, double ki, int _model);
extern void newviewAncestralIterative(pllInstance *tr, partitionList *pr);
extern void printAncestralState(nodeptr p, pllBoolean printStates, pllBoolean printProbs, pllInstance *tr, partitionList *pr);
extern void makenewzGeneric(pllInstance *tr, partitionList * pr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, pllBoolean mask);
extern void makenewzGenericDistance(pllInstance *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (pllInstance *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (pllInstance *tr, nodeptr p, int model);
extern double evaluateGenericVector (pllInstance *tr, nodeptr p);
extern void categorizeGeneric (pllInstance *tr, nodeptr p);
extern double makenewzPartitionGeneric(pllInstance *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern pllBoolean isTip(int number, int maxTips);

/* recom functions */
extern void computeTraversal(pllInstance *tr, nodeptr p, pllBoolean partialTraversal, int numBranches);
extern void allocRecompVectorsInfo(pllInstance *tr);
extern void allocTraversalCounter(pllInstance *tr);
extern pllBoolean getxVector(recompVectors *rvec, int nodenum, int *slot, int mxtips);
extern pllBoolean needsRecomp(pllBoolean recompute, recompVectors *rvec, nodeptr p, int mxtips);
extern void unpinNode(recompVectors *v, int nodenum, int mxtips);
extern void protectNode(recompVectors *rvec, int nodenum, int mxtips);

/* Handling branch lengths*/
extern void computeTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec, int *count);
extern void computeFullTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec);
extern void printTraversalInfo(pllInstance *tr);
extern void countTraversal(pllInstance *tr);
extern void storeExecuteMaskInTraversalDescriptor(pllInstance *tr, partitionList *pr);
extern void storeValuesInTraversalDescriptor(pllInstance *tr, partitionList *pr, double *value);
extern void makenewzIterative(pllInstance *, partitionList *pr);
extern void execCore(pllInstance *, partitionList *pr, volatile double *dlnLdlz, volatile double *d2lnLdlz2);
extern void makePermutation(int *perm, int n, pllInstance *tr);
extern nodeptr findAnyTip(nodeptr p, int numsp);
extern void putWAG(double *ext_initialRates);
extern  unsigned int **initBitVector(int mxtips, unsigned int *vectorLength);
//extern hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(pllHashTable * h, int state);
extern double convergenceCriterion(pllHashTable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
//extern void freeHashTable(hashtable *h);
//extern stringHashtable *initStringHashTable(hashNumberType n);
//extern void addword(char *s, stringHashtable *h, int nodeNumber);
extern void initRateMatrix(pllInstance *tr, partitionList *pr);
extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, pllHashTable *h, int treeNumber, int function, branchInfo *bInf,
                                    int *countBranches, int treeVectorLength, pllBoolean traverseOnly, pllBoolean computeWRF, int processID);
extern  unsigned int bitcount_32_bit(unsigned int i);
extern __inline unsigned int bitcount_64_bit(uint64_t i);
extern void perSiteLogLikelihoods(pllInstance *tr, partitionList *pr, double *logLikelihoods);
extern void updatePerSiteRates(pllInstance *tr, partitionList *pr, pllBoolean scaleRates);
extern void restart(pllInstance *tr, partitionList *pr);

//extern const unsigned int mask32[32];

/** @brief Check whether the position \a pos in bitvector \a x is a gap

    @param x
      A bitvector represented by unsigned integers

    @param pos
      Position to check in \a x if it is set (i.e. it is a gap)

    @return
      Returns the value of the bit vector (\b 1 if set, \b 0 if not)
*/
//#ifndef __clang__
//inline
//#endif
pllBoolean isGap(unsigned int *x, int pos);

/** @brief Check whether the position \a pos in bitvector \a x is \b NOT a gap

    @param x
      A bitvector represented by unsigned integers

    @param pos
      Position to check in \a x if it is \b NOT set (i.e. it is \b NOT a gap)

    @return
      Returns the value of the bit vector (\b 1 if set, \b 0 if not)
*/
//#ifndef __clang__
//inline
//#endif
pllBoolean noGap(unsigned int *x, int pos);

//#ifndef __clang__
//__inline
//#endif
//pllBoolean isGap(unsigned int *x, int pos);

//#ifndef __clang__
//__inline
//#endif
//pllBoolean noGap(unsigned int *x, int pos);

/* from utils.h */
linkageList* initLinkageList(int *linkList, partitionList *pr);

#if (defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI) || defined(_USE_PTHREADS) )
/* work tags for parallel regions */

#define PLL_THREAD_NEWVIEW                  0
#define PLL_THREAD_EVALUATE                 1
#define PLL_THREAD_MAKENEWZ                 2
#define PLL_THREAD_MAKENEWZ_FIRST           3
#define PLL_THREAD_RATE_CATS                4
#define PLL_THREAD_COPY_RATE_CATS           5
#define PLL_THREAD_COPY_INIT_MODEL          6
#define PLL_THREAD_INIT_PARTITION           7
#define PLL_THREAD_OPT_ALPHA                8
#define PLL_THREAD_OPT_RATE                 9
#define PLL_THREAD_OPT_LG4X_RATE            10
#define PLL_THREAD_COPY_ALPHA               11
#define PLL_THREAD_COPY_RATES               12
#define PLL_THREAD_COPY_LG4X_RATES          13
#define PLL_THREAD_PER_SITE_LIKELIHOODS     14
#define PLL_THREAD_NEWVIEW_ANCESTRAL        15
#define PLL_THREAD_GATHER_ANCESTRAL         16
#define PLL_THREAD_EXIT_GRACEFULLY          17
#define PLL_THREAD_EVALUATE_PER_SITE_LIKES  18


typedef struct
{
  pllInstance *tr;

  partitionList *pr;
  int threadNumber;
}
  threadData;
extern void optRateCatPthreads(pllInstance *tr, partitionList *pr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid);
extern void pllMasterBarrier(pllInstance *, partitionList *, int);
#endif


#ifdef __AVX

extern void newviewGTRGAMMAPROT_AVX_LG4(int tipCase,
                                        double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
                                        int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n,
                                        double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);


extern void newviewGTRCAT_AVX_GAPPED_SAVE(int tipCase,  double *EV,  int *cptr,
                                   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
                                   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

extern void newviewGTRCATPROT_AVX_GAPPED_SAVE(int tipCase, double *extEV,
                                       int *cptr,
                                       double *x1, double *x2, double *x3, double *tipVector,
                                       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                       int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
                                       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

extern void  newviewGTRGAMMA_AVX_GAPPED_SAVE(int tipCase,
                                      double *x1_start, double *x2_start, double *x3_start,
                                      double *extEV, double *tipVector,
                                      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                      const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
                                      unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                      double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
                                      );

extern void newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(int tipCase,
                                         double *x1_start, double *x2_start, double *x3_start, double *extEV, double *tipVector,
                                         int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n,
                                         double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
                                         unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                         double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn);

extern void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
    double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);


extern void newviewGenericCATPROT_AVX(int tipCase, double *extEV,
    int *cptr,
    double *x1, double *x2, double *x3, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);


extern void newviewGTRGAMMA_AVX(int tipCase,
    double *x1_start, double *x2_start, double *x3_start,
    double *EV, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);

extern void newviewGTRGAMMAPROT_AVX(int tipCase,
                             double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                             int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n,
                             double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);

extern void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
                           int *cptr,
                           double *x1, double *x2, double *x3, double *tipVector,
                           int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                           int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);

#endif

extern int virtual_width( int n );
extern void computeAllAncestralVectors(nodeptr p, pllInstance *tr, partitionList *pr);

#endif /* PLLINTERNAL_H_ */
