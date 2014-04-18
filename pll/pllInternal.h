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
#include "hash.h"
#include "lexer.h"
#include "parsePartition.h"
#include "mem_alloc.h"

extern int lookupWord(char *s, stringHashtable *h);

extern void getDataTypeString(pllInstance *tr, pInfo *partitionInfo, char typeOfData[1024]);
extern int countTips(nodeptr p, int numsp);
extern entry *initEntry(void);
extern unsigned int precomputed16_bitcount(unsigned int n, char *bits_in_16bits);

extern size_t discreteRateCategories(int rateHetModel);

extern const partitionLengths * getPartitionLengths(pInfo *p);
extern boolean getSmoothFreqs(int dataType);
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
extern boolean whitechar ( int ch );
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

extern boolean initrav ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void initravPartition ( pllInstance *tr, nodeptr p, int model );
extern void update ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void smooth ( pllInstance *tr, partitionList *pr, nodeptr p );
extern void smoothTree ( pllInstance *tr, partitionList *pr, int maxtimes );
extern void localSmooth ( pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes );
extern boolean localSmoothMulti(pllInstance *tr, nodeptr p, int maxtimes, int model);

extern void smoothRegion ( pllInstance *tr, partitionList *pr, nodeptr p, int region );
extern void regionalSmooth ( pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( pllInstance *tr, partitionList *pr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p );
extern boolean insertBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q);
extern boolean insertRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern boolean testInsertBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern int NNI(pllInstance * tr, nodeptr p, int swap);
extern void addTraverseBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern boolean testInsertRestoreBIG ( pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( pllInstance *tr, partitionList *pr );

extern void initTL ( topolRELL_LIST *rl, pllInstance *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, pllInstance *tr, int n, int numBranches );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, pllInstance *tr, int index );

extern int  saveBestTree (bestlist *bt, pllInstance *tr, int numBranches);
extern int  recallBestTree (bestlist *bt, int rank, pllInstance *tr, partitionList *pr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern boolean freeBestTree ( bestlist *bt );

extern int treeReadLen (FILE *fp, pllInstance *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly);
extern void treeReadTopologyString(char *treeString, pllInstance *tr);
extern void getStartingTree (pllInstance *tr);
extern double treeLength (pllInstance *tr, int model);
extern double evaluatePartialGeneric (pllInstance *, partitionList *pr, int i, double ki, int _model);
extern void newviewAncestralIterative(pllInstance *tr, partitionList *pr);
extern void printAncestralState(nodeptr p, boolean printStates, boolean printProbs, pllInstance *tr, partitionList *pr);
extern void makenewzGeneric(pllInstance *tr, partitionList * pr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask);
extern void makenewzGenericDistance(pllInstance *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (pllInstance *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (pllInstance *tr, nodeptr p, int model);
extern double evaluateGenericVector (pllInstance *tr, nodeptr p);
extern void categorizeGeneric (pllInstance *tr, nodeptr p);
extern double makenewzPartitionGeneric(pllInstance *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern boolean isTip(int number, int maxTips);

/* recom functions */
extern void computeTraversal(pllInstance *tr, nodeptr p, boolean partialTraversal, int numBranches);
extern void allocRecompVectorsInfo(pllInstance *tr);
extern void allocTraversalCounter(pllInstance *tr);
extern boolean getxVector(recompVectors *rvec, int nodenum, int *slot, int mxtips);
extern boolean needsRecomp(boolean recompute, recompVectors *rvec, nodeptr p, int mxtips);
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
extern hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(hashtable *h, int state);
extern double convergenceCriterion(hashtable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
extern void freeHashTable(hashtable *h);
extern stringHashtable *initStringHashTable(hashNumberType n);
extern void addword(char *s, stringHashtable *h, int nodeNumber);
extern void printBothOpen(const char* format, ... );
extern void initRateMatrix(pllInstance *tr, partitionList *pr);
extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf,
                                    int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF, int processID);
extern  unsigned int bitcount_32_bit(unsigned int i);
extern inline unsigned int bitcount_64_bit(unsigned long i);
extern void perSiteLogLikelihoods(pllInstance *tr, partitionList *pr, double *logLikelihoods);
extern void updatePerSiteRates(pllInstance *tr, partitionList *pr, boolean scaleRates);
extern void restart(pllInstance *tr, partitionList *pr);
inline boolean isGap(unsigned int *x, int pos);
inline boolean noGap(unsigned int *x, int pos);

/* from utils.h */
linkageList* initLinkageList(int *linkList, partitionList *pr);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS) )
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
#define PLL_THREAD_COPY_ALPHA               10
#define PLL_THREAD_COPY_RATES               11
#define PLL_THREAD_PER_SITE_LIKELIHOODS     12
#define PLL_THREAD_NEWVIEW_ANCESTRAL        13
#define PLL_THREAD_GATHER_ANCESTRAL         14
#define PLL_THREAD_EXIT_GRACEFULLY          15
#define PLL_THREAD_EVALUATE_PER_SITE_LIKES  16


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
                                        double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);


extern void newviewGTRCAT_AVX_GAPPED_SAVE(int tipCase,  double *EV,  int *cptr,
                                   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
                                   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

extern void newviewGTRCATPROT_AVX_GAPPED_SAVE(int tipCase, double *extEV,
                                       int *cptr,
                                       double *x1, double *x2, double *x3, double *tipVector,
                                       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
                                       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

extern void  newviewGTRGAMMA_AVX_GAPPED_SAVE(int tipCase,
                                      double *x1_start, double *x2_start, double *x3_start,
                                      double *extEV, double *tipVector,
                                      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                      const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
                                      unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                      double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
                                      );

extern void newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(int tipCase,
                                         double *x1_start, double *x2_start, double *x3_start, double *extEV, double *tipVector,
                                         int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n,
                                         double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
                                         unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                         double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn);

extern void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
    double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);


extern void newviewGenericCATPROT_AVX(int tipCase, double *extEV,
    int *cptr,
    double *x1, double *x2, double *x3, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);


extern void newviewGTRGAMMA_AVX(int tipCase,
    double *x1_start, double *x2_start, double *x3_start,
    double *EV, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

extern void newviewGTRGAMMAPROT_AVX(int tipCase,
                             double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                             int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n,
                             double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

extern void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
                           int *cptr,
                           double *x1, double *x2, double *x3, double *tipVector,
                           int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                           int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling);

#endif

extern int virtual_width( int n );
extern void computeAllAncestralVectors(nodeptr p, pllInstance *tr, partitionList *pr);

#endif /* PLLINTERNAL_H_ */
