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
 * @file utils.c
 *  
 * @brief Miscellaneous general utility and helper functions
 */

#include "systypes.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <errno.h>
#include "cycle.h"
#include "mem_alloc.h" //for rax_malloc_string_copy



#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#if (defined(__AVX) || defined(__SSE3))
#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#endif
#endif
/*
   special bug fix, enforces denormalized numbers to be flushed to zero,
   without this program is a tiny bit faster though.
#include <emmintrin.h> 
#define MM_DAZ_MASK    0x0040
#define MM_DAZ_ON    0x0040
#define MM_DAZ_OFF    0x0000
*/
#endif

#include "pll.h"
#include "pllInternal.h"

void rax_malloc_string_copy(const char* source, char** dest)
{
    size_t bufLen = (strlen(source) + 1);
    *dest = (char*)rax_malloc(bufLen);
    #ifdef CLANG_UNDER_VS
        strcpy_s(*dest, bufLen, source);
    #else
        strcpy(*dest, source);
    #endif
}

#define GLOBAL_VARIABLES_DEFINITION

#include "globalVariables.h"


/* mappings of BIN/DNA/AA alphabet to numbers */

static const char PLL_MAP_BIN[256] =
 {
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  3, -1, -1,
    1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  3,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

static const char PLL_MAP_NT[256] =
 {
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15,
   -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, 15,
   -1, -1,  5,  6,  8,  8,  7,  9, 15, 10, -1, -1, -1, -1, -1, -1,
   -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, 15,
   -1, -1,  5,  6,  8,  8,  7,  9, 15, 10, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
 };

static const char PLL_MAP_AA[256] =
 {
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 22, -1, -1, 22, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 22,
   -1,  0, 20,  4,  3,  6, 13,  7,  8,  9, -1, 11, 10, 12,  2, -1,
   14,  5,  1, 15, 16, -1, 19, 17, 22, 18, 21, -1, -1, -1, -1, -1,
   -1,  0, 20,  4,  3,  6, 13,  7,  8,  9, -1, 11, 10, 12,  2, -1,
   14,  5,  1, 15, 16, -1, 19, 17, 22, 18, 21, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
 };





static void pllTreeInitDefaults (pllInstance * tr, int tips);
static void getInnerBranchEndPointsRecursive (nodeptr p, int tips, int * i, node **nodes);
#if (!defined(_FINE_GRAIN_MPI) && !defined(_USE_PTHREADS))
static void initializePartitionsSequential(pllInstance *tr, partitionList *pr);
#endif

/** @defgroup instanceLinkingGroup Linking topology, partition scheme and alignment to the PLL instance
    
    This set of functions handles the linking of topology, partition scheme and multiple sequence alignment
    with the PLL instance
*/
/***************** UTILITY FUNCTIONS **************************/

char *my_strndup(const char *s, size_t n) {
    char *ret = (char *) rax_malloc(n+1);
    strncpy(ret, s, n);
    ret[n] = 0;
    return ret;
}

#ifndef HAVE_STRTOK_R
char *strtok_r (char * s, const char * delim, char **save_ptr)
{  
  char *token;
   
  /* Scan leading delimiters */
  if (s == NULL) {
      s = *save_ptr;
  }   
  s += strspn (s, delim);
  if (*s == '\0')
   {
     *save_ptr = s;
     return NULL;
   }
   
  /* Find the end of the token. */
  token = s;
  s = strpbrk (token, delim);
  if (!s)
    *save_ptr = strchr (token, '\0');
  else
   {
     /* Terminate the token and make *SAVE_PTR point past it */
     *s = '\0';
     *save_ptr = s + 1;
   }
   
  return token;
}
#endif


void storeExecuteMaskInTraversalDescriptor(pllInstance *tr, partitionList *pr)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].executeModel[model] = pr->partitionData[model]->executeModel;

}

void storeValuesInTraversalDescriptor(pllInstance *tr, partitionList *pr, double *value)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].parameterValues[model] = value[model];
}

const unsigned int *getBitVector(int dataType)
{
  assert(PLL_MIN_MODEL < dataType && dataType < PLL_MAX_MODEL);

  return pLengths[dataType].bitVector;
}

/*
int getStates(int dataType)
{
  assert(PLL_MIN_MODEL < dataType && dataType < PLL_MAX_MODEL);

  return pLengths[dataType].states;
}
*/

int getUndetermined(int dataType)
{
  assert(PLL_MIN_MODEL < dataType && dataType < PLL_MAX_MODEL);

  return pLengths[dataType].undetermined;
}

const partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
              states    = p->states,
              tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(PLL_MIN_MODEL < dataType && dataType < PLL_MAX_MODEL);

  /*pLength.leftLength = pLength.rightLength = states * states;
    pLength.eignLength = states;
    pLength.evLength   = states * states;
    pLength.eiLength   = states * states;
    pLength.substRatesLength = (states * states - states) / 2;
    pLength.frequenciesLength = states;
    pLength.tipVectorLength   = tipLength * states;
    pLength.symmetryVectorLength = (states * states - states) / 2;
    pLength.frequencyGroupingLength = states;
    pLength.nonGTR = PLL_FALSE;*/
  return (&pLengths[dataType]); 
}

size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
  {
    case PLL_CAT:
      result = 1;
      break;
    case PLL_GAMMA:
      result = 4;
      break;
    default:
      result = 0;
      assert(0);
  }

  return result;
}



double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}


/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


/** @brief Check whether a node is a tip.
    
    Checks whether the node with number \a number is a tip.
    
    @param number
     Node number to be checked
   
    @param maxTips
     Number of tips in the tree
   
    @return
      \b PLL_TRUE if tip, \b PLL_FALSE otherwise
  */
pllBoolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return PLL_TRUE;
  else
    return PLL_FALSE;
}

/** @brief Set the orientation of a node

    Sets the orientation of node \a p. That means, it will reset the orientation
    \a p->next->x and \a p->next->next->x to 0 and of \a p->x to 1, meaning that
    the conditional likelihood vector for that node is oriented on \a p, i.e.
    the conditional likelihood vector represents the subtree rooted at \a p and
    not any other of the two nodes.

    @param p
      Node which we want to orient
*/
void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
  {
    p->x = s->x;
    s->x = 0;
  }

  assert(p->x);
}


/** @brief Connect two nodes and assign branch lengths 
  * 
  * Connect the two nodes \a p and \a q in each partition \e i with a branch of
  * length \a z[i]
  *
  * @param p
  *   Node \a p
  * 
  * @param q
  *   Node \a q
  *
  * @param numBranches
  *   Number of partitions
  */
void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

/* connects node p with q and assigns the branch lengths z for the whole vector*/
void hookupFull (nodeptr p, nodeptr q, double *z)
{
  //int i;

  p->back = q;
  q->back = p;

  memcpy(p->z, z, PLL_NUM_BRANCHES*sizeof(double) );
  memcpy(q->z, z, PLL_NUM_BRANCHES*sizeof(double) );
  //for(i = 0; i < numBranches; i++)
  //  p->z[i] = q->z[i] = z[i];

}

/* connect node p with q and assign the default branch lengths */
void hookupDefault (nodeptr p, nodeptr q)
{
  p->back = q;
  q->back = p;

// TODO: fix: this make parsimony tree computation very slow with increasing PLL_NUM_BRANCHES
//  for(i = 0; i < PLL_NUM_BRANCHES; i++)
//    p->z[i] = q->z[i] = PLL_DEFAULTZ;
    p->z[0] = q->z[0] = PLL_DEFAULTZ;
}


/***********************reading and initializing input ******************/



pllBoolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}
/*
static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
      y = 362436069,
      z = 21288629,
      w = 14921776,
      c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}
*/

/** @brief Get a random subtree

    Returns the root node of a randomly picked subtree of the tree in PLL
    instance \a tr. The picked subtree is guaranteed to have height over
    1, that is, the direct descendents of the returned (root) node are not tips.

    @param tr
      PLL instance

    @return
      The root node of the randomly picked subtree
*/
nodeptr pllGetRandomSubtree(pllInstance *tr)
{
  nodeptr p;
  do
  {
    int exitDirection = rand() % 3; 
    p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];
    switch(exitDirection)
    {
      case 0:
        break;
      case 1:
        p = p->next;
        break;
      case 2:
        p = p->next->next;
        break;
      default:
        assert(0);
    }
  }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));
  return p;
}
/* small example program that executes ancestral state computations 
   on the entire subtree rooted at p.

   Note that this is a post-order traversal.
*/

  
void computeAllAncestralVectors(nodeptr p, pllInstance *tr, partitionList *pr)
{
  /* if this is not a tip, for which evidently it does not make sense 
     to compute the ancestral sequence because we have the real one ....
  */

  if(!isTip(p->number, tr->mxtips))
    {
      /* descend recursively to compute the ancestral states in the left and right subtrees */

      computeAllAncestralVectors(p->next->back, tr, pr);
      computeAllAncestralVectors(p->next->next->back, tr, pr);
      
      /* then compute the ancestral state at node p */

      pllUpdatePartialsAncestral(tr, pr, p);

      /* and print it to terminal, the two booleans that are set to PLL_TRUE here 
         tell the function to print the marginal probabilities as well as 
         a discrete inner sequence, that is, ACGT etc., always selecting and printing 
         the state that has the highest probability */

      printAncestralState(p, PLL_TRUE, PLL_TRUE, tr, pr);
    }
}



void initializePartitionData(pllInstance *localTree, partitionList * localPartitions)
{
  /* in ancestralVectorWidth we store the total length in bytes (!) of 
     one conditional likelihood array !
     we need to know this length such that in the pthreads version the master thread can actually 
     gather the scattered ancestral probabilities from the threads such that they can be printed to screen!
  */

  size_t 
    maxCategories = (size_t)localTree->maxCategories;

  size_t 
    ancestralVectorWidth = 0,
    model; 

  int 
    tid  = localTree->threadID,
    innerNodes = localTree->mxtips - 2;

  if(tid > 0)
      localTree->rateCategory    = (int *)    rax_calloc((size_t)localTree->originalCrunchedLength, sizeof(int));           

  for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
    {
      size_t 
        width = localPartitions->partitionData[model]->width;

      const partitionLengths 
        *pl = getPartitionLengths(localPartitions->partitionData[model]);

      /* 
         globalScaler needs to be 2 * localTree->mxtips such that scalers of inner AND tip nodes can be added without a case switch
         to this end, it must also be initialized with zeros -> calloc
      */

      localPartitions->partitionData[model]->globalScaler    = (unsigned int *)rax_calloc(2 *(size_t)localTree->mxtips, sizeof(unsigned int));

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->left),  PLL_BYTE_ALIGNMENT, (size_t)pl->leftLength * (maxCategories + 1) * sizeof(double));
      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->right), PLL_BYTE_ALIGNMENT, (size_t)pl->rightLength * (maxCategories + 1) * sizeof(double));
      localPartitions->partitionData[model]->EIGN              = (double*)rax_malloc((size_t)pl->eignLength * sizeof(double));
      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->EV),    PLL_BYTE_ALIGNMENT, (size_t)pl->evLength * sizeof(double));
      localPartitions->partitionData[model]->EI                = (double*)rax_malloc((size_t)pl->eiLength * sizeof(double));
      localPartitions->partitionData[model]->substRates        = (double *)rax_malloc((size_t)pl->substRatesLength * sizeof(double));
      localPartitions->partitionData[model]->frequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->freqExponents     = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->empiricalFrequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->tipVector), PLL_BYTE_ALIGNMENT, (size_t)pl->tipVectorLength * sizeof(double));
      //localPartitions->partitionData[model]->partitionName      = NULL;   // very imporatant since it is deallocated in pllPartitionDestroy
      
       if(localPartitions->partitionData[model]->dataType == PLL_AA_DATA
               && (localPartitions->partitionData[model]->protModels == PLL_LG4M || localPartitions->partitionData[model]->protModels == PLL_LG4X))
        {
          int 
            k;
          
          for(k = 0; k < 4; k++)
            {       
              localPartitions->partitionData[model]->EIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));
              rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->EV_LG4[k]), PLL_BYTE_ALIGNMENT, pl->evLength * sizeof(double));
              localPartitions->partitionData[model]->EI_LG4[k]                = (double*)rax_malloc(pl->eiLength * sizeof(double));
              localPartitions->partitionData[model]->substRates_LG4[k]        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));
              localPartitions->partitionData[model]->frequencies_LG4[k]       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
              rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->tipVector_LG4[k]), PLL_BYTE_ALIGNMENT, pl->tipVectorLength * sizeof(double));
            }
        }

      localPartitions->partitionData[model]->symmetryVector    = (int *)rax_malloc((size_t)pl->symmetryVectorLength  * sizeof(int));
      localPartitions->partitionData[model]->frequencyGrouping = (int *)rax_malloc((size_t)pl->frequencyGroupingLength  * sizeof(int));

      localPartitions->partitionData[model]->perSiteRates      = (double *)rax_malloc(sizeof(double) * maxCategories);

      localPartitions->partitionData[model]->nonGTR = PLL_FALSE;

      localPartitions->partitionData[model]->gammaRates = (double*)rax_malloc(sizeof(double) * 4);
      localPartitions->partitionData[model]->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char*) * ((size_t)localTree->mxtips + 1));


      localPartitions->partitionData[model]->xVector = (double **)rax_calloc(sizeof(double*), (size_t)localTree->mxtips);

      if (localPartitions->partitionData[model]->ascBias)
       {
         localPartitions->partitionData[model]->ascOffset    = 4 * localPartitions->partitionData[model]->states * localPartitions->partitionData[model]->states;
         localPartitions->partitionData[model]->ascVector    = (double *)rax_malloc(innerNodes * 
                                                                                    localPartitions->partitionData[model]->ascOffset * 
                                                                                    sizeof(double));
         localPartitions->partitionData[model]->ascExpVector = (int *)rax_calloc(innerNodes *
                                                                                 localPartitions->partitionData[model]->states,
                                                                                 sizeof(int));
         localPartitions->partitionData[model]->ascSumBuffer = (double *)rax_malloc(localPartitions->partitionData[model]->ascOffset * sizeof(double)); 
       }


      /* 
         Initializing the xVector array like this is absolutely required !!!!
         I don't know which programming genius removed this, but it must absolutely stay in here!!!!
      */
      
      {
        int k;
        
        for(k = 0; k < localTree->mxtips; k++)
              localPartitions->partitionData[model]->xVector[k] = (double*)NULL;       
      }


      localPartitions->partitionData[model]->xSpaceVector = (size_t *)rax_calloc((size_t)localTree->mxtips, sizeof(size_t));

      const size_t span = (size_t)(localPartitions->partitionData[model]->states) *
              discreteRateCategories(localTree->rateHetModel);

#ifdef __MIC_NATIVE

      // Alexey: sum buffer buffer padding for Xeon PHI
      const int aligned_width = width % PLL_VECTOR_WIDTH == 0 ? width : width + (PLL_VECTOR_WIDTH - (width % PLL_VECTOR_WIDTH));

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->sumBuffer), PLL_BYTE_ALIGNMENT, aligned_width *
                                                                                      span *
                                                                                      sizeof(double));

      // Alexey: fill padding entries with 1. (will be corrected with site weights, s. below)
      {
          int k;
          for (k = width*span; k < aligned_width*span; ++k)
              localPartitions->partitionData[model]->sumBuffer[k] = 1.;
      }

#else

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->sumBuffer), PLL_BYTE_ALIGNMENT, width *
                                              span *
                                              sizeof(double));
#endif

      /* Initialize buffers to store per-site log likelihoods */

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->perSiteLikelihoods), PLL_BYTE_ALIGNMENT, width * sizeof(double));

      /* initialize data structures for per-site likelihood scaling */

      if(localTree->fastScaling)
        {
           localPartitions->partitionData[model]->expVector      = (int **)NULL;
           localPartitions->partitionData[model]->expSpaceVector = (size_t *)NULL;
        }
      else
        {        
          localPartitions->partitionData[model]->expVector      = (int **)rax_malloc(sizeof(int*) * innerNodes);
           
          /* 
             Initializing the expVector array like this is absolutely required !!!!
             Not doing this can (and did) cause segmentation faults !!!!
          */
          
          {
            int k;

            for(k = 0; k < innerNodes; k++)
              localPartitions->partitionData[model]->expVector[k] = (int*)NULL; 
          }

          localPartitions->partitionData[model]->expSpaceVector = (size_t *)rax_calloc(innerNodes, sizeof(size_t));
        }

      /* data structure to store the marginal ancestral probabilities in the sequential version or for each thread */

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->ancestralBuffer), PLL_BYTE_ALIGNMENT, width *
                                                                                 (size_t)(localPartitions->partitionData[model]->states) *
                                                                                 sizeof(double));

      /* count and accumulate how many bytes we will need for storing a full ancestral vector. for this we addf over the per-partition space requirements in bytes */
      /* ancestralVectorWidth += ((size_t)(pr->partitionData[model]->upper - pr->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double)); */
      ancestralVectorWidth += ((size_t)(localPartitions->partitionData[model]->upper - localPartitions->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double));
      /* :TODO: do we have to use the original tree for that   */

#ifdef __MIC_NATIVE

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->wgt), PLL_BYTE_ALIGNMENT, aligned_width * sizeof(int));

      // Alexey: fill padding entries with 0.
      {
          int k;
          for (k = width; k < aligned_width; ++k)
              localPartitions->partitionData[model]->wgt[k] = 0;
      }
#else
      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->wgt), PLL_BYTE_ALIGNMENT, width * sizeof(int));
#endif

      /* rateCategory must be assigned using rax_calloc() at start up there is only one rate category 0 for all sites */

      localPartitions->partitionData[model]->rateCategory = (int *)rax_calloc(width, sizeof(int));

      if(width > 0 && localTree->saveMemory)
        {
          localPartitions->partitionData[model]->gapVectorLength = ((int)width / 32) + 1;
          assert(4 == sizeof(unsigned int));
          localPartitions->partitionData[model]->gapVector = (unsigned int*)rax_calloc((size_t)localPartitions->partitionData[model]->gapVectorLength * 2 * (size_t)localTree->mxtips, sizeof(unsigned int));
          rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->gapColumn),PLL_BYTE_ALIGNMENT, ((size_t)localTree->mxtips) *
                                                                               ((size_t)(localPartitions->partitionData[model]->states)) *
                                                                               discreteRateCategories(localTree->rateHetModel) * sizeof(double));
        }
      else
        {
          localPartitions->partitionData[model]->gapVectorLength = 0;
          localPartitions->partitionData[model]->gapVector = (unsigned int*)NULL;
          localPartitions->partitionData[model]->gapColumn = (double*)NULL;
        }              
    }
}

int virtual_width( int n ) {
    const int global_vw = 2;
    return (n+1) / global_vw * global_vw;
}


void initMemorySavingAndRecom(pllInstance *tr, partitionList *pr)
{
  pllInstance  
    *localTree = tr; 
  partitionList
    *localPartitions = pr;
  size_t model; 

  /* initialize gap bit vectors at tips when memory saving option is enabled */

  if(localTree->saveMemory)
    {
      for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
        {
          int        
            undetermined = getUndetermined(localPartitions->partitionData[model]->dataType);

          size_t
            i,
            j,
            width =  localPartitions->partitionData[model]->width;

          if(width > 0)
            {                                        
              for(j = 1; j <= (size_t)(localTree->mxtips); j++)
                for(i = 0; i < width; i++)
                  if(localPartitions->partitionData[model]->yVector[j][i] == undetermined)
                    localPartitions->partitionData[model]->gapVector[localPartitions->partitionData[model]->gapVectorLength * j + i / 32] |= mask32[i % 32];
            }     
        }
    }
  /* recom */
  if(localTree->useRecom)
    allocRecompVectorsInfo(localTree);
  else
    localTree->rvec = (recompVectors*)NULL;
  /* E recom */
}

/** @brief Get the length of a specific branch

    Get the length of the branch specified by node \a p and \a p->back
    of partition \a partition_id.
    The branch length is decoded from the PLL representation.

    @param tr
      PLL instance

    @param p
      Specifies one end-point of the branch. The other one is \a p->back

    @param partition_id
      Specifies the partition

    @return
      The branch length
*/
double pllGetBranchLength (pllInstance *tr, nodeptr p, int partition_id)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < PLL_NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z = p->z[partition_id];
  if(z < PLL_ZMIN) z = PLL_ZMIN;
  if(z > PLL_ZMAX) z = PLL_ZMAX;
  return (-log(z) * tr->fracchange);
}

/** @brief Set the length of a specific branch

    Set the length of the branch specified by node \a p and \a p->back
    of partition \a partition_id.
    The function encodes the branch length to the PLL representation.

    @param tr
      PLL instance

    @param p
      Specifies one end-point of the branch. The other one is \a p->back

    @param partition_id
      Specifies the partition

    @param bl
      Branch length
*/
void pllSetBranchLength (pllInstance *tr, nodeptr p, int partition_id, double bl)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < PLL_NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z;
  z = exp((-1 * bl)/tr->fracchange);
  if(z < PLL_ZMIN) z = PLL_ZMIN;
  if(z > PLL_ZMAX) z = PLL_ZMAX;
  p->z[partition_id] = z;
}

#if (!defined(_FINE_GRAIN_MPI) && !defined(_USE_PTHREADS))
static void initializePartitionsSequential(pllInstance* tr, partitionList* pr)
{
    for (size_t model = 0; model < (size_t)pr->numberOfPartitions; model++) {
        assert(pr->partitionData[model]->width == pr->partitionData[model]->upper - pr->partitionData[model]->lower);
    }
    initializePartitionData(tr, pr);

    /* figure in tip sequence data per-site pattern weights */
    for (size_t model = 0; model < (size_t)pr->numberOfPartitions; model++)
    {
        size_t lower = pr->partitionData[model]->lower;
        size_t width = pr->partitionData[model]->upper - lower;
        for (size_t j = 1; j <= (size_t)tr->mxtips; j++)
        {
            pr->partitionData[model]->yVector[j] = &(tr->yVector[j][pr->partitionData[model]->lower]);
        }
        memcpy((void*)(&(pr->partitionData[model]->wgt[0])), (void*)(&(tr->aliaswgt[lower])), sizeof(int) * width);
    }

    initMemorySavingAndRecom(tr, pr);
}
#endif


/* interface to outside  */
//void initializePartitions(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
//{
//#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
//  initializePartitionsMaster(tr,localTree,pr,localPr,tid,n);
//#else
//  initializePartitionsSequential(tr, pr);
//#endif
//}

static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    rax_free(ll->ld[i].partitionList);         

  rax_free(ll->ld);
  rax_free(ll);   
}

/** @brief free all data structures associated to a partition
    
    frees all data structures allocated for this partition

    @param partitions
      the pointer to the partition list

    @param tips  
       number of tips in the tree      
*/
void 
pllPartitionsDestroy (pllInstance * tr, partitionList ** partitions)
{
  int i, j, tips;
  partitionList * pl = *partitions;

#ifdef _USE_PTHREADS
  int tid = tr->threadID;
  if (MASTER_P) {
     pllMasterBarrier (tr, pl, PLL_THREAD_EXIT_GRACEFULLY);
     pllStopPthreads (tr);
    }
#endif

  tips = tr->mxtips;

#ifdef _USE_PTHREADS
  if (MASTER_P) {
#endif
#ifdef _FINE_GRAIN_MPI
if (MASTER_P) {
     pllMasterBarrier (tr, pl, PLL_THREAD_EXIT_GRACEFULLY);
#endif
  freeLinkageList(pl->alphaList);
  freeLinkageList(pl->freqList); 
  freeLinkageList(pl->rateList);
#ifdef _FINE_GRAIN_MPI
}
#endif

#ifdef _USE_PTHREADS
  }
#endif
  for (i = 0; i < pl->numberOfPartitions; ++ i)
   {
     rax_free (pl->partitionData[i]->gammaRates);
     rax_free (pl->partitionData[i]->perSiteRates);
     rax_free (pl->partitionData[i]->globalScaler);
     rax_free (pl->partitionData[i]->left);
     rax_free (pl->partitionData[i]->right);
     rax_free (pl->partitionData[i]->EIGN);
     rax_free (pl->partitionData[i]->EV);
     rax_free (pl->partitionData[i]->EI);
     rax_free (pl->partitionData[i]->substRates);
     rax_free (pl->partitionData[i]->frequencies);
     rax_free (pl->partitionData[i]->freqExponents);
     rax_free (pl->partitionData[i]->empiricalFrequencies);
     rax_free (pl->partitionData[i]->tipVector);
     rax_free (pl->partitionData[i]->symmetryVector);
     rax_free (pl->partitionData[i]->frequencyGrouping);
     for (j = 0; j < tips; ++ j)
       rax_free (pl->partitionData[i]->xVector[j]);
     rax_free (pl->partitionData[i]->xVector);
     rax_free (pl->partitionData[i]->yVector);
     rax_free (pl->partitionData[i]->xSpaceVector);
     rax_free (pl->partitionData[i]->sumBuffer);
     rax_free (pl->partitionData[i]->ancestralBuffer);
     rax_free (pl->partitionData[i]->wgt);
     rax_free (pl->partitionData[i]->rateCategory);
     rax_free (pl->partitionData[i]->gapVector);
     rax_free (pl->partitionData[i]->gapColumn);
     rax_free (pl->partitionData[i]->perSiteLikelihoods);
     rax_free (pl->partitionData[i]->partitionName);
     rax_free (pl->partitionData[i]->expSpaceVector);
     /*TODO: Deallocate all entries of expVector */
     if (pl->partitionData[i]->expVector)
      {
        for (j = 0; j < tips - 2; ++ j)
          rax_free (pl->partitionData[i]->expVector[j]);
      }
     rax_free (pl->partitionData[i]->expVector);
     rax_free (pl->partitionData[i]);
   }
  rax_free (pl->partitionData);
  rax_free (pl);

  *partitions = NULL;

#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
     rax_free (tr->y_ptr);
#endif
}

/** @ingroup instanceLinkingGroup
    @brief Correspondance check between partitions and alignment

    This function checks whether the partitions to be created and the given
    alignment correspond, that is, whether each site of the alignment is
    assigned to exactly one partition.

    @param parts
      A list of partitions suggested by the caller

    @param alignmentData
      The multiple sequence alignment
    
    @return
      Returns \a 1 in case of success, otherwise \a 0
*/
int
pllPartitionsValidate (pllQueue * parts, pllAlignmentData * alignmentData)
{
  int nparts;
  char * used;
  struct pllQueueItem * elm;
  struct pllQueueItem * regionItem;
  pllPartitionRegion * region;
  pllPartitionInfo * pi;
  int i;

  /* check if the list contains at least one partition */
  nparts = pllQueueSize (parts);
  if (!nparts)          
    return (0);   

  /* pllBoolean array for marking that a site was assigned a partition */
  used = (char *) rax_calloc (alignmentData->sequenceLength, sizeof (char));

  /* traverse all partitions and their respective regions and mark sites */
  for (elm = parts->head; elm; elm = elm->next)
   {
     pi = (pllPartitionInfo *) elm->item;
     
     for (regionItem = pi->regionList->head; regionItem; regionItem = regionItem->next)
      {
        region = (pllPartitionRegion *) regionItem->item;
        
        if (region->start < 1 || region->end > alignmentData->sequenceLength) 
         {
           rax_free (used);
           return (0);
         }

        for (i = region->start - 1; i < region->end; i += region->stride)
         {
           if (used[i])
            {
              rax_free (used);
              return (0);
            }
           used[i] = 1; 
         }
      }
   }

  /* check whether all sites were assigned a partition */
  for (i = 0; i < alignmentData->sequenceLength; ++ i)
    if (used[i] != 1)
     {
       rax_free (used);
       return (0);
     }

  rax_free (used);
  return (1);
}

/** @brief Swap two sites in a buffer
    
    Swaps sites \a s1 and \a s2 in buffer \a buf which consists of \a nTaxa + 1
    taxa (i.e. rows), and the first row contains no information, i.e. it is not
    accessed.

    @param buffer
      Memory buffer

    @param s1
      First site

    @param s2
      Second site

    @param nTaxa
      Number of taxa, i.e. size of site
*/
static __inline void
swapSite (unsigned char ** buf, int s1, int s2, int nTaxa)
{
  int i;
  int x;

  for (i = 1; i <= nTaxa; ++ i)
  {
    x = buf[i][s1];
    buf[i][s1] = buf[i][s2];
    buf[i][s2] = x;
  }
}

/** @brief Constructs the list of partitions according to the proposed partition scheme
    
    A static function that construcs the \a partitionList structure according to
    the partition scheme \b AFTER the sites have been repositioned in contiguous
    regions according to the partition scheme.

    @param bounds  An array of the new starting and ending posititons of sites
    in the alignment for each partition.  This array is of size 2 * \a nparts.
    The elements are always couples (lower,upper). The upper bounds is a site
    that is not included in the partition

    @param nparts The number of partitions to be created

    @todo Fix the bug in PLL 
*/
static partitionList * createPartitions (pllQueue * parts, int * bounds)
{
  partitionList * pl;
  pllPartitionInfo * pi;
  struct pllQueueItem * elm;
  int i, j;

  pl = (partitionList *) rax_malloc (sizeof (partitionList));
  
  // TODO: fix this
  pl->perGeneBranchLengths =      0;

  // TODO: change PLL_NUM_BRANCHES to number of partitions I guess
  pl->partitionData = (pInfo **) rax_calloc (PLL_NUM_BRANCHES, sizeof (pInfo *));
  
  for (i = 0, elm = parts->head; elm; elm = elm->next, ++i)
  {
      pi = (pllPartitionInfo*)elm->item;

      /* check whether the data type is valid, and in case it's not, deallocate
         and return NULL */
      if (pi->dataType <= PLL_MIN_MODEL || pi->dataType >= PLL_MAX_MODEL)
      {
          for (j = 0; j < i; ++j)
          {
              rax_free(pl->partitionData[j]->partitionName);
              rax_free(pl->partitionData[j]);
          }
          rax_free(pl->partitionData);
          rax_free(pl);
          return (NULL);
      }

      pl->partitionData[i] = (pInfo*)rax_malloc(sizeof(pInfo));

      pl->partitionData[i]->lower = bounds[i << 1];
      pl->partitionData[i]->upper = bounds[(i << 1) + 1];
      pl->partitionData[i]->width = bounds[(i << 1) + 1] - bounds[i << 1];
      pl->partitionData[i]->partitionWeight = 1.0 * (double)pl->partitionData[i]->width;

      //the two flags below are required to allow users to set 
      //alpha parameters and substitution rates in the Q matrix 
      //to fixed values. These parameters will then not be optimized 
      //in the model parameter optimization functions
      //by default we assume that all parameters are being optimized, i.e., 
      //this has to be explicitly set by the user 

      pl->partitionData[i]->optimizeAlphaParameter = PLL_TRUE;
      pl->partitionData[i]->optimizeSubstitutionRates = PLL_TRUE;
      pl->partitionData[i]->dataType = pi->dataType;
      pl->partitionData[i]->protModels = -1;
      pl->partitionData[i]->protUseEmpiricalFreqs = -1;
      pl->partitionData[i]->maxTipStates = pLengths[pi->dataType].undetermined + 1;
      pl->partitionData[i]->optimizeBaseFrequencies = pi->optimizeBaseFrequencies;
      pl->partitionData[i]->ascBias = pi->ascBias;
      pl->partitionData[i]->parsVect = NULL;



      if (pi->dataType == PLL_AA_DATA)
      {
          if (pl->partitionData[i]->protModels != PLL_GTR)
              pl->partitionData[i]->optimizeSubstitutionRates = PLL_FALSE;
          pl->partitionData[i]->protUseEmpiricalFreqs = pi->protUseEmpiricalFreqs;
          pl->partitionData[i]->protModels = pi->protModels;
      }

      pl->partitionData[i]->states = pLengths[pl->partitionData[i]->dataType].states;
      pl->partitionData[i]->numberOfCategories = 1;
      pl->partitionData[i]->autoProtModels = 0;
      pl->partitionData[i]->nonGTR = PLL_FALSE;
      pl->partitionData[i]->partitionContribution = -1.0;
      pl->partitionData[i]->partitionLH = 0.0;
      pl->partitionData[i]->fracchange = 1.0;
      pl->partitionData[i]->executeModel = PLL_TRUE;

      rax_malloc_string_copy(pi->partitionName, &(pl->partitionData[i]->partitionName));
   }

  return (pl);
}


/** @ingroup instanceLinkingGroup
    @brief Constructs the proposed partition scheme 

    This function constructs the proposed partition scheme. It assumes
    that the partition scheme is correct.

    @note This function \b does \b not validate the partition scheme.
    The user must manually call the ::pllPartitionsValidate function
    for validation
    
    @param parts
      A list of partitions suggested by the caller

    @param alignmentData
      The multiple sequence alignment

    @return
      Returns a pointer to \a partitionList structure of partitions in case of success, \b NULL otherwise
*/
partitionList * pllPartitionsCommit (pllQueue * parts, pllAlignmentData * alignmentData)
{
  int * oi;
  int i, j, dst;
  struct pllQueueItem * elm;
  struct pllQueueItem * regionItem;
  pllPartitionRegion * region;
  pllPartitionInfo * pi;
  partitionList * pl;
  int * newBounds;
  int k, nparts;
  int tmpvar;
 

  dst = k = 0;
  oi  = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i) oi[i] = i;

  nparts = pllQueueSize (parts);
  newBounds = (int *) rax_malloc (2 * nparts * sizeof (int));

  /* reposition the sites in the alignment */
  for (elm = parts->head; elm; elm = elm->next, ++ k)
   {
     pi = (pllPartitionInfo *) elm->item;
     
     newBounds[k << 1] = dst;   /* set the lower column for this partition */
     for (regionItem = pi->regionList->head; regionItem; regionItem = regionItem->next)
      {
        region = (pllPartitionRegion *) regionItem->item;

        for (i = region->start - 1; i < region->end && i < alignmentData->sequenceLength; i += region->stride)
         {
           if (oi[i] == i)
            {
              swapSite (alignmentData->sequenceData, dst, i, alignmentData->sequenceCount);
              tmpvar = oi[i];
              oi[i] = oi[dst];
              oi[dst++] = tmpvar;
            }
           else
            {
              j = i;
              while (oi[j] != i) j = oi[j];

              swapSite (alignmentData->sequenceData, dst, j, alignmentData->sequenceCount);
              tmpvar = oi[j];
              oi[j] = oi[dst];
              oi[dst++] = tmpvar;
            }
         }
      }
     newBounds[(k << 1) + 1] = dst;    /* set the uppwer limit for this partition */
   }
  if ((pl = createPartitions (parts, newBounds)))
   { 
     pl->numberOfPartitions = nparts;
     pl->dirty = PLL_FALSE;
   }
  
  rax_free (newBounds);
  rax_free (oi);

  return (pl);
}

/** @brief Copy a site to another buffer

    Copies site \a from from buffer \a src to \a to in buffer \a dst. Both buffers
    must consist of \a nTaxa + 1 taxa and the first row contains no information, i.e.
    it is not accessed.

    @param dst
      Destination buffer

    @param src
      Source buffer

    @param to
      At which position in \a dst to copy the site to

    @param from
      Which site from \a src to copy

    @param nTaxa
      Number of taxa, i.e. size of site
*/
static __inline void
copySite (unsigned char ** dst, unsigned char ** src, int to, int from, int nTaxa)
{
  int i;

  for (i = 1; i <= nTaxa; ++ i)
   {
     dst[i][to] = src[i][from];
   }
}

/** @brief Remove duplicate sites from alignment and update weights vector

    Removes duplicate sites from the alignment given the partitions list
    and updates the weight vector of the alignment and the boundaries
    (upper, lower, width) for each partition.

    @param alignmentData
      The multiple sequence alignment
    
    @param pl
      List of partitions

*/
void 
pllAlignmentRemoveDups (pllAlignmentData * alignmentData, partitionList * pl)
{
  int i, j, k, p;
  char *** sites;
  void ** memptr;
  int ** oi;
  int dups = 0;
  int lower;

  /* allocate space for the transposed alignments (sites) for every partition */
  sites  = (char ***) rax_malloc (pl->numberOfPartitions * sizeof (char **));
  memptr = (void **)  rax_malloc (pl->numberOfPartitions * sizeof (void *));
  oi     = (int **)   rax_malloc (pl->numberOfPartitions * sizeof (int *));

  /* transpose the sites by partition */
  for (p = 0; p < pl->numberOfPartitions; ++ p)
   {
     sites[p]  = (char **) rax_malloc (sizeof (char *) * pl->partitionData[p]->width);
     memptr[p] = rax_malloc (sizeof (char) * (alignmentData->sequenceCount + 1) * pl->partitionData[p]->width);

     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        sites[p][i] = (char *) ((char*)memptr[p] + sizeof (char) * i * (alignmentData->sequenceCount + 1));
      }

     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        for (j = 0; j < alignmentData->sequenceCount; ++ j)
         {
           sites[p][i][j] = alignmentData->sequenceData[j + 1][pl->partitionData[p]->lower + i]; 
         }
        sites[p][i][j] = 0;
      }

     oi[p] = pllssort1main (sites[p], pl->partitionData[p]->width);

     for (i = 0; i < pl->partitionData[p]->width; ++ i) oi[p][i] = 1;

     for (i = 1; i < pl->partitionData[p]->width; ++ i)
      {
        if (! strcmp (sites[p][i], sites[p][i - 1]))
         {
           ++ dups;
           oi[p][i] = 0;
         }
      }
   }

  /* allocate memory for the alignment without duplicates*/
  rax_free (alignmentData->sequenceData[1]);
  rax_free (alignmentData->siteWeights);

  alignmentData->sequenceLength = alignmentData->sequenceLength - dups;
  alignmentData->sequenceData[0] = (unsigned char *) rax_malloc ((alignmentData->sequenceLength + 1) * sizeof (unsigned char) * alignmentData->sequenceCount);
  for (i = 0; i < alignmentData->sequenceCount; ++ i)
   {
     alignmentData->sequenceData[i + 1] = (unsigned char *) (alignmentData->sequenceData[0] + sizeof (unsigned char) * i * (alignmentData->sequenceLength + 1));
     alignmentData->sequenceData[i + 1][alignmentData->sequenceLength] = 0;
   }

  alignmentData->siteWeights    = (int *) rax_malloc ((alignmentData->sequenceLength) * sizeof (int));
  alignmentData->siteWeights[0] = 1;

  /* transpose sites back to alignment */
  for (p = 0, k = 0; p < pl->numberOfPartitions; ++ p)
   {
     lower = k;
     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        if (!oi[p][i])
         {
           ++ alignmentData->siteWeights[k - 1];
         }
        else
         {
           alignmentData->siteWeights[k] = 1;
           for (j = 0; j < alignmentData->sequenceCount; ++ j)
            {
              alignmentData->sequenceData[j + 1][k] = sites[p][i][j];
            }
           ++ k;
         }
      }
     pl->partitionData[p]->lower = lower;
     pl->partitionData[p]->upper = k;
     pl->partitionData[p]->width = k - lower;
   }

  /* deallocate storage for transposed alignment (sites) */
  for (p = 0; p < pl->numberOfPartitions; ++ p)
   {
     rax_free (oi[p]);
     rax_free (memptr[p]);
     rax_free (sites[p]);
   }
  rax_free (oi);
  rax_free (sites);
  rax_free (memptr);
}


/** @brief Compute the empirical frequencies of a partition
  
    Compute the empirical frequencies of partition \a partition and store them in
    \a pfreqs.

    @param partition
      The partition for which to compute empirical frequencies

    @param alignmentData
      The multiple sequence alignment

    @param smoothFrequencies
      Not needed?

    @param bitMask
      The bitmask

    @param pfreqs
      Array of size \a partition->states where the empirical frequencies for this partition are stored
*/
static int genericBaseFrequenciesAlignment (pInfo * partition, 
                                              pllAlignmentData * alignmentData, 
                                              pllBoolean smoothFrequencies,
                                              const unsigned int * bitMask, 
                                              double * pfreqs)
{
  double 
    wj, 
    acc,
    sumf[64],   
    temp[64];
 
  int     
    i, 
    j, 
    k, 
    l,
    numFreqs,
    lower,
    upper;

  unsigned char  *yptr;  
  const char * map;
  
  switch (partition->dataType)
   {
     case PLL_BINARY_DATA:
       map = PLL_MAP_BIN;
     case PLL_DNA_DATA:
       map = PLL_MAP_NT;
       break;
     case PLL_AA_DATA:
       map = PLL_MAP_AA;
       break;
     default:
       assert(0);
   }

  numFreqs = partition->states;
  lower    = partition->lower;
  upper    = partition->upper;

  for (l = 0; l < numFreqs; l++) {
      pfreqs[l] = 1.0 / ((double)numFreqs);
  }          
  for (k = 1; k <= 8; k++) 
    {                                                   
      for(l = 0; l < numFreqs; l++)
        sumf[l] = 0.0;
              
      for (i = 1; i <= alignmentData->sequenceCount; i++) 
        {                
          yptr = alignmentData->sequenceData[i];
          
          for(j = lower; j < upper; j++) 
            {
              if (map[yptr[j]] < 0) return (0);
              unsigned int code = bitMask[(unsigned char)map[yptr[j]]];
              assert(code >= 1);
              
              for(l = 0; l < numFreqs; l++)
                {
                  if((code >> l) & 1)
                    temp[l] = pfreqs[l];
                  else
                    temp[l] = 0.0;
                }                             
              
              for(l = 0, acc = 0.0; l < numFreqs; l++)
                {
                  if(temp[l] != 0.0)
                    acc += temp[l];
                }
              
              wj = alignmentData->siteWeights[j] / acc;
              
              for(l = 0; l < numFreqs; l++)
                {
                  if(temp[l] != 0.0)                
                    sumf[l] += wj * temp[l];                                                                                               
                }
            }
        }                     
      
      for(l = 0, acc = 0.0; l < numFreqs; l++)
        {
          if(sumf[l] != 0.0)
            acc += sumf[l];
        }
              
      for(l = 0; l < numFreqs; l++)
        pfreqs[l] = sumf[l] / acc;           
    }

   /* TODO: What is that? */
/*
  if(smoothFrequencies)         
   {;
    smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));     
   }
  else    
    {
      pllBoolean
        zeroFreq = PLL_FALSE;

      char 
        typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);  

      for(l = 0; l < numFreqs; l++)
        {
          if(pfreqs[l] == 0.0)
            {
              printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
              printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
              zeroFreq = PLL_TRUE;
            }
        }

      if(zeroFreq)
        exit(-1);

      for(l = 0; l < numFreqs; l++)
        {
          assert(pfreqs[l] > 0.0);
          tr->partitionData[model].frequencies[l] = pfreqs[l];
        }     
    }  
*/
  return (1);
  
}

static void  genericBaseFrequenciesInstance (pInfo * partition, 
                                             pllInstance * tr, 
                                             pllBoolean smoothFrequencies,
                                             const unsigned int * bitMask, 
                                             double * pfreqs)
{
  double 
    wj, 
    acc,
    sumf[64],   
    temp[64];
 
  int     
    i, 
    j, 
    k, 
    l,
    numFreqs,
    lower,
    upper;

  unsigned char  *yptr;  

  numFreqs = partition->states;
  lower    = partition->lower;
  upper    = partition->upper;

  for(l = 0; l < numFreqs; l++)     
    pfreqs[l] = 1.0 / ((double)numFreqs);
          
  for (k = 1; k <= 8; k++) 
    {                                                   
      for(l = 0; l < numFreqs; l++)
        sumf[l] = 0.0;
              
      for (i = 1; i <= tr->mxtips; i++) 
        {                
          yptr = tr->yVector[i];
          
          for(j = lower; j < upper; j++) 
            {
              unsigned int code = bitMask[yptr[j]];
              assert(code >= 1);
              
              for(l = 0; l < numFreqs; l++)
                {
                  if((code >> l) & 1)
                    temp[l] = pfreqs[l];
                  else
                    temp[l] = 0.0;
                }                             
              
              for(l = 0, acc = 0.0; l < numFreqs; l++)
                {
                  if(temp[l] != 0.0)
                    acc += temp[l];
                }
              
              wj = tr->aliaswgt[j] / acc;
              
              for(l = 0; l < numFreqs; l++)
                {
                  if(temp[l] != 0.0)                
                    sumf[l] += wj * temp[l];                                                                                               
                }
            }
        }                     
      
      for(l = 0, acc = 0.0; l < numFreqs; l++)
        {
          if(sumf[l] != 0.0)
            acc += sumf[l];
        }
              
      for(l = 0; l < numFreqs; l++)
        pfreqs[l] = sumf[l] / acc;           
    }

   /* TODO: What is that? */
/*
  if(smoothFrequencies)         
   {;
    smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));     
   }
  else    
    {
      pllBoolean
        zeroFreq = PLL_FALSE;

      char 
        typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);  

      for(l = 0; l < numFreqs; l++)
        {
          if(pfreqs[l] == 0.0)
            {
              printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
              printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
              zeroFreq = PLL_TRUE;
            }
        }

      if(zeroFreq)
        exit(-1);

      for(l = 0; l < numFreqs; l++)
        {
          assert(pfreqs[l] > 0.0);
          tr->partitionData[model].frequencies[l] = pfreqs[l];
        }     
    }  
*/

  
}

/**  Compute the empirical base frequencies of an alignment

     Computes the empirical base frequencies per partition of an alignment \a alignmentData
     given the partition structure \a pl.

     @param alignmentData The alignment structure for which to compute the empirical base frequencies
     @param pl            List of partitions
     @return Returns a list of frequencies for each partition
*/
double ** pllBaseFrequenciesAlignment (pllAlignmentData * alignmentData, partitionList * pl)
{
  int
    i,
    model;

  double 
    **freqs = (double **) rax_malloc (pl->numberOfPartitions * sizeof (double *));

  for (model = 0; model < pl->numberOfPartitions; ++ model)
    {
      freqs[model] = (double *) rax_malloc (pl->partitionData[model]->states * sizeof (double));
      
      switch  (pl->partitionData[model]->dataType)
        {
        case PLL_BINARY_DATA:
        case PLL_AA_DATA:
        case PLL_DNA_DATA:
          if (!genericBaseFrequenciesAlignment (pl->partitionData[model], 
                                                alignmentData, 
                                                pLengths[pl->partitionData[model]->dataType].smoothFrequencies,
                                                pLengths[pl->partitionData[model]->dataType].bitVector,
                                                freqs[model]
                                               ))
            return (NULL);
          break;
        default:
          {
            errno = PLL_UNKNOWN_MOLECULAR_DATA_TYPE;
            for (i = 0; i <= model; ++ i) rax_free (freqs[i]);
            rax_free (freqs);
            return (double **)NULL;
          }
        }
    }
  
  return (freqs);
}

/**  Compute the empirical base frequencies of the alignment incorporated in the instance

     Computes the empirical base frequencies per partition of the alignment
     incorporated in the instance \a tr given the partition structure \a pl.

     @param tr The instance for which to compute the empirical base frequencies
     @param pl List of partitions
     @return Returns a list of frequencies for each partition
*/
double ** pllBaseFrequenciesInstance (pllInstance * tr, partitionList * pl)
{
  int
    i,
    model;

  double 
    **freqs = (double **) rax_malloc (pl->numberOfPartitions * sizeof (double *));

  for (model = 0; model < pl->numberOfPartitions; ++ model)
    {
      freqs[model] = (double *) rax_malloc (pl->partitionData[model]->states * sizeof (double));
      
      switch  (pl->partitionData[model]->dataType)
        {
        case PLL_AA_DATA:
        case PLL_DNA_DATA:
        case PLL_BINARY_DATA:
          genericBaseFrequenciesInstance (pl->partitionData[model], 
                                          tr, 
                                          pLengths[pl->partitionData[model]->dataType].smoothFrequencies,
                                          pLengths[pl->partitionData[model]->dataType].bitVector,
                                          freqs[model]
                                          );
          break;
        default:
          {
            errno = PLL_UNKNOWN_MOLECULAR_DATA_TYPE;
            for (i = 0; i <= model; ++ i) rax_free (freqs[i]);
            rax_free (freqs);
            return (double **)NULL;
          }
        }
    }
  
  return (freqs);
}

void
pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models)
{
  int i;

  for (i = 0; i < models; ++ i)
   {
     rax_free ((*empiricalFrequencies)[i]);
   }
  rax_free (*empiricalFrequencies);

  *empiricalFrequencies = NULL;
}

int pllLoadAlignment (pllInstance * tr, pllAlignmentData * alignmentData, partitionList * partitions)
{
  int i;
  nodeptr node;
  pllHashItem * hItem;

  if (tr->mxtips != alignmentData->sequenceCount) return (0);

  tr->aliaswgt = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  memcpy (tr->aliaswgt, alignmentData->siteWeights, alignmentData->sequenceLength * sizeof (int));

  tr->originalCrunchedLength = alignmentData->sequenceLength;
  tr->rateCategory           = (int *)   rax_calloc (tr->originalCrunchedLength, sizeof (int));
  tr->patrat                 = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored           = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->lhs                    = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));

  /* allocate memory for the alignment */
  tr->yVector    = (unsigned char **) rax_malloc ((alignmentData->sequenceCount + 1) * sizeof (unsigned char *));                                                                                                                                                                      

  tr->yVector[0] = (unsigned char *)  rax_malloc (sizeof (unsigned char) * (alignmentData->sequenceLength + 1) * alignmentData->sequenceCount);
  for (i = 1; i <= alignmentData->sequenceCount; ++ i) 
   {                     
     tr->yVector[i] = (unsigned char *) (tr->yVector[0] + (i - 1) * (alignmentData->sequenceLength + 1) * sizeof (unsigned char));
     tr->yVector[i][alignmentData->sequenceLength] = 0;
   }                     
                         
  /* place sequences to tips */                              
  for (i = 1; i <= alignmentData->sequenceCount; ++ i)                      
   {                     
     if (!pllHashSearch (tr->nameHash, alignmentData->sequenceLabels[i],(void **)&node)) 
      {
        //rax_free (tr->originalCrunchedLength);
        rax_free (tr->rateCategory);
        rax_free (tr->patrat);
        rax_free (tr->patratStored);
        rax_free (tr->lhs);
        rax_free (tr->yVector[0]);
        rax_free (tr->yVector);
        return (0);
      }
     memcpy (tr->yVector[node->number], alignmentData->sequenceData[i], alignmentData->sequenceLength);
   }

  /* Do the base substitution (from A,C,G....  ->   0,1,2,3....)*/
  pllBaseSubstitute (tr, partitions);

  /* Populate tipNames */
  tr->tipNames = (char **) rax_calloc(tr->mxtips + 1, sizeof (char *));
  for (i = 0; (unsigned int)i < tr->nameHash->size; ++ i)
   {
     hItem = tr->nameHash->Items[i];

     for (; hItem; hItem = hItem->next)
      {
        tr->tipNames[((nodeptr)hItem->data)->number] = hItem->str; 
      }
   }

  return (1);
}

pllInstance * pllCreateInstance (pllInstanceAttr * attr)
{
  pllInstance * tr;

  if (attr->rateHetModel != PLL_GAMMA && attr->rateHetModel != PLL_CAT) return NULL;

  tr = (pllInstance *) rax_calloc (1, sizeof (pllInstance));

  tr->threadID          = 0;
  tr->rateHetModel      = attr->rateHetModel;
  tr->fastScaling       = attr->fastScaling;
  tr->saveMemory        = attr->saveMemory;
  tr->useRecom          = attr->useRecom;
  tr->likelihoodEpsilon = 0.01;
  
  tr->randomNumberSeed = attr->randomNumberSeed;
  tr->parsimonyScore   = NULL;

  /* remove it from the library */
  tr->useMedian         = PLL_FALSE;

  tr->maxCategories     = (attr->rateHetModel == PLL_GAMMA) ? 4 : 25;

  tr->numberOfThreads   = attr->numberOfThreads;
  tr->rearrangeHistory  = NULL;

  /* Lock the slave processors at this point */
#ifdef _FINE_GRAIN_MPI
  pllLockMPI (tr);
#endif

  return (tr);
}

/** @brief Initialize PLL tree structure with default values
    
    Initialize PLL tree structure with default values and allocate 
    memory for its elements.

    @todo
      STILL NOT FINISHED
*/
static void pllTreeInitDefaults (pllInstance * tr, int tips)
{
  nodeptr p0, p, q;
  int i, j;
  int inner;

  

  /* TODO: make a proper static setupTree function */

  inner = tips - 1;

  tr->mxtips = tips;

  tr->bigCutoff = PLL_FALSE;
  tr->treeStringLength = tr->mxtips * (PLL_NMLNGTH + 128) + 256 + tr->mxtips * 2;
  tr->tree_string = (char *) rax_calloc ( tr->treeStringLength, sizeof(char));
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->constraintVector = (int *)rax_malloc((2 * tr->mxtips) * sizeof(int));
  
  p0 = (nodeptr) rax_malloc ((tips + 3 * inner) * sizeof (node));
  assert (p0);

  tr->nodeBaseAddress  = p0;

  tr->nameList         = (char **)   rax_malloc ((tips + 1) * sizeof (char *));
  tr->nodep            = (nodeptr *) rax_malloc ((2 * tips) * sizeof (nodeptr));

  tr->autoProteinSelectionType = PLL_AUTO_ML;

  assert (tr->nameList && tr->nodep);

  tr->nodep[0] = NULL;          


  /* TODO: The line below was commented... why? */
  tr->fracchange = -1;
  tr->rawFracchange = -1;

  for (i = 1; i <= tips; ++ i)
   {
     p = p0++;

     //p->hash      = KISS32();     
     p->x         = 0;
     p->xBips     = 0;
     p->number    = i;
     p->next      = p;
     p->back      = NULL;
     p->bInf      = NULL;
     tr->nodep[i]  = p;
   }

  for (i = tips + 1; i <= tips + inner; ++i)
   {
     q = NULL;
     for (j = 1; j <= 3; ++ j)
     {
       p = p0++;
       if (j == 1)
        {
          p->xBips = 1;
          p->x = 1; //p->x     = 1;
        }
       else
        {
          p->xBips = 0;
          p->x     = 0;
        }
       p->number = i;
       p->next   = q;
       p->bInf   = NULL;
       p->back   = NULL;
       p->hash   = 0;
       q         = p;
     }
    p->next->next->next = p;
    tr->nodep[i]         = p;
   }

  tr->likelihood  = PLL_UNLIKELY;
  tr->start       = NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;

  for (i = 0; i < PLL_NUM_BRANCHES; ++ i) tr->partitionSmoothed[i] = PLL_FALSE;

  tr->bitVectors = NULL;
  tr->vLength    = 0;
  //tr->h          = NULL;

  /* TODO: Fix hash type */
  tr->nameHash   = pllHashInit (10 * tr->mxtips);

  /* TODO: do these options really fit here or should they be put elsewhere? */
  tr->td[0].count            = 0;
  tr->td[0].ti               = (traversalInfo *) rax_malloc (sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].parameterValues  = (double *) rax_malloc(sizeof(double) * (size_t)PLL_NUM_BRANCHES);
  tr->td[0].executeModel     = (pllBoolean *) rax_malloc (sizeof(pllBoolean) * (size_t)PLL_NUM_BRANCHES);
  tr->td[0].executeModel[0]  = PLL_TRUE;                                                                                                                                                                                                                                    
  for (i = 0; i < PLL_NUM_BRANCHES; ++ i) tr->td[0].executeModel[i] = PLL_TRUE;
}


/* @brief Check a parsed tree for inclusion in the current tree
   
   Check whether the set of leaves (taxa) of the parsed tree \a nTree is a
   subset of the leaves of the currently loaded tree.

   @param pInst
     PLL instance

   @param nTree
     Parsed newick tree structure

   @return
     Returns \b PLL_TRUE in case it is a subset, otherwise \b PLL_FALSE
*/
static int
checkTreeInclusion (pllInstance * pInst, pllNewickTree * nTree)
{
  pllStack * sList;
  pllNewickNodeInfo * sItem;
  void * dummy;

  if (!pInst->nameHash) return (PLL_FALSE);

  for (sList = nTree->tree; sList; sList = sList->next)
   {
     sItem = (pllNewickNodeInfo *) sList->item;
     if (!sItem->rank)   /* leaf */
      {
        if (!pllHashSearch (pInst->nameHash, sItem->name, &dummy)) return (PLL_FALSE);
      }
   }

  return (PLL_TRUE);
}

static void
updateBranchLength (nodeptr p, double old_fracchange, double new_fracchange)
{
  double z;
  int j;

  for (j = 0; j < PLL_NUM_BRANCHES; ++ j)
   {
     z = exp ((log (p->z[j]) * old_fracchange) / new_fracchange);
     if (z < PLL_ZMIN) z = PLL_ZMIN;
     if (z > PLL_ZMAX) z = PLL_ZMAX;
     p->z[j] = p->back->z[j] = z;
   }
}

static void
updateAllBranchLengthsRecursive (nodeptr p, int tips, double old_fracchange, double new_fracchange)
{
  updateBranchLength (p, old_fracchange, new_fracchange);

  if (!isTip (p->number, tips))
   {
     updateAllBranchLengthsRecursive (p->next->back,       tips, old_fracchange, new_fracchange);
     updateAllBranchLengthsRecursive (p->next->next->back, tips, old_fracchange, new_fracchange);
   }
}

static void
updateAllBranchLengths (pllInstance * tr, double old_fracchange, double new_fracchange)
{
  nodeptr p;

  p = tr->start;
  assert (isTip(p->number, tr->mxtips));

  updateAllBranchLengthsRecursive (p->back, tr->mxtips, old_fracchange, new_fracchange);

}


/** @brief Relink the taxa
    
    Relink the taxa by performing a preorder traversal of the unrooted binary tree.
    We assume that the tree is rooted such that the root is the only node of
    out-degree 3 and in-degree 0, while all the other inner nodes have in-degree
    1 and out-degree 2. Finally, the leaves have in-degree 1 and out-degree 0.

    @param pInst
      PLL instance

    @param nTree
      Parsed newick tree structure

    @param taxaExist
      Is the set of taxa of \a nTree a subset of the taxa of the current tree

    @return
*/
static int
linkTaxa (pllInstance * pInst, pllNewickTree * nTree, int taxaExist)
{
  nodeptr 
    parent,
    child;
  pllStack 
    * nodeStack = NULL,
    * current;
  int
    i,
    j,
    inner = nTree->tips + 1,
    leaf  = 1;
  double z;
  pllNewickNodeInfo * nodeInfo;

  if (!taxaExist) pllTreeInitDefaults (pInst, nTree->tips);

  /* Place the ternary root node 3 times on the stack such that later on
     three nodes use it as their parent */
  current = nTree->tree;
  for (parent = pInst->nodep[inner], i  = 0; i < 3; ++ i, parent = parent->next)
    pllStackPush (&nodeStack, parent);
  ++ inner;

  /* now traverse the rest of the nodes */
  for (current = current->next; current; current = current->next)
   {
     parent   = (nodeptr) pllStackPop (&nodeStack);
     nodeInfo = (pllNewickNodeInfo *) current->item;

     /* if inner node place it twice on the stack (out-degree 2) */
     if (nodeInfo->rank)
      {
        child = pInst->nodep[inner ++];
        pllStackPush (&nodeStack, child->next);
        pllStackPush (&nodeStack, child->next->next);
      }
     else /* check if taxon already exists, i.e. we loaded another tree topology */
      {
        if (taxaExist)
         {
           assert (pllHashSearch (pInst->nameHash, nodeInfo->name, (void **) &child));
         }
        else
         {
           child = pInst->nodep[leaf];
           pInst->nameList[leaf] = strdup (nodeInfo->name);
           pllHashAdd (pInst->nameHash, pllHashString(pInst->nameList[leaf], pInst->nameHash->size), pInst->nameList[leaf], (void *) (pInst->nodep[leaf]));
           ++ leaf;
         }
      }
     assert (parent);
     /* link parent and child */
     parent->back = child;
     child->back  = parent;

     if (!taxaExist) pInst->fracchange = 1;

     /* set the branch length */
     z = exp ((-1 * atof (nodeInfo->branch)) / pInst->fracchange);
     if (z < PLL_ZMIN) z = PLL_ZMIN;
     if (z > PLL_ZMAX) z = PLL_ZMAX;
     for (j = 0; j < PLL_NUM_BRANCHES; ++ j)
       parent->z[j] = child->z[j] = z;
   }
  pllStackClear (&nodeStack);

  return PLL_TRUE;
}

/** @brief Get the instantaneous rate matrix
    
    Obtain the instantaneous rate matrix (Q) for partitionm \a model
    of the partition list \a pr, and store it in an array \a outBuffer.
    
    @param tr        PLL instance
    @param pr        List of partitions
    @param model     Index of partition to use
    @param outBuffer Where to store the instantaneous rate matrix 

    @todo Currently, the Q matrix can be only obtained for DNA GTR data.

    @return Returns \b PLL_TRUE in case of success, otherwise \b PLL_FALSE
*/
int pllGetInstRateMatrix (partitionList * pr, int model, double * outBuffer)
{
  if (pr->partitionData[model]->dataType != PLL_DNA_DATA) return (PLL_FALSE);

  int  i;
  double mean = 0;
  double * substRates = pr->partitionData[model]->substRates;
  double * freqs = pr->partitionData[model]->frequencies;
  
  /* normalize substitution rates */
  for (i = 0; i < 6; ++ i)  substRates[i] /= substRates[5];

  outBuffer[0 * 4 + 1] = (substRates[0] * freqs[1]);
  outBuffer[0 * 4 + 2] = (substRates[1] * freqs[2]);
  outBuffer[0 * 4 + 3] = (substRates[2] * freqs[3]);

  outBuffer[1 * 4 + 0] = (substRates[0] * freqs[0]);
  outBuffer[1 * 4 + 2] = (substRates[3] * freqs[2]);
  outBuffer[1 * 4 + 3] = (substRates[4] * freqs[3]);

  outBuffer[2 * 4 + 0] = (substRates[1] * freqs[0]);
  outBuffer[2 * 4 + 1] = (substRates[3] * freqs[1]);
  outBuffer[2 * 4 + 3] = (substRates[5] * freqs[3]);

  outBuffer[3 * 4 + 0] = (substRates[2] * freqs[0]);
  outBuffer[3 * 4 + 1] = (substRates[4] * freqs[1]);
  outBuffer[3 * 4 + 2] = (substRates[5] * freqs[2]);

  outBuffer[0 * 4 + 0] = -(substRates[0] * freqs[1] + substRates[1] * freqs[2] + substRates[2] * freqs[3]);
  outBuffer[1 * 4 + 1] = -(substRates[0] * freqs[0] + substRates[3] * freqs[2] + substRates[4] * freqs[3]);
  outBuffer[2 * 4 + 2] = -(substRates[1] * freqs[0] + substRates[3] * freqs[1] + substRates[5] * freqs[3]);
  outBuffer[3 * 4 + 3] = -(substRates[2] * freqs[0] + substRates[4] * freqs[1] + substRates[5] * freqs[2]);

  for (i = 0; i <  4; ++ i) mean         += freqs[i] * (-outBuffer[i * 4 + i]);
  for (i = 0; i < 16; ++ i) outBuffer[i] /= mean;

  return (PLL_TRUE);
}

/** @ingroup instanceLinkingGroup
    @brief Initializes the PLL tree topology according to a parsed newick tree

    Set the tree topology based on a parsed and validated newick tree

    @param tree
      The PLL instance

    @param nt
      The \a pllNewickTree wrapper structure that contains the parsed newick tree

    @param useDefaultz
      If set to \b PLL_TRUE then the branch lengths will be reset to the default
      value.
*/
void
pllTreeInitTopologyNewick (pllInstance * tr, pllNewickTree * newick, int useDefaultz)
{
  linkTaxa (tr, newick, tr->nameHash && checkTreeInclusion (tr, newick));

  tr->start = tr->nodep[1];

  if (useDefaultz == PLL_TRUE)
    resetBranches (tr);
}

/** @brief Get the node oriented pointer from a round-about node

    Returns the pointer of the round-about node $p$ that has the orientation, i.e.
    has the \a x flag set to 1. In case a tip is passed, then the returned pointer
    is the same as the input.

    @param pInst  PLL instance
    @param p      One of the three pointers of a round-about node

    @return  Returns the the pointer that has the orientation
*/
nodeptr pllGetOrientedNodePointer (pllInstance * pInst, nodeptr p)
{
  if (p->number <= pInst->mxtips || p->x) return p;

  if (p->next->x) return p->next;

  return p->next->next;
}


//void
//pllTreeInitTopologyNewick (pllInstance * tr, pllNewickTree * nt, int useDefaultz)
//{
//  pllStack * nodeStack = NULL;
//  pllStack * head;
//  pllNewickNodeInfo * item;
//  int i, j, k;
//  
///*
//  for (i = 0; i < partitions->numberOfPartitions; ++ i)
//   {
//     partitions->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));
//     partitions->partitionData[i]->partitionContribution = -1.0;
//     partitions->partitionData[i]->partitionLH           =  0.0;
//     partitions->partitionData[i]->fracchange            =  1.0;
//   }
//*/
// 
//
// if (tr->nameHash)
//  {
//    if (checkTreeInclusion (tr, nt))
//     {
//       printf ("It is a subset\n");
//     }
//    else
//     {
//       printf ("It is not a subset\n");
//     }
//  }
//  
//  pllTreeInitDefaults (tr, nt->tips);
//
//  i = nt->tips + 1;
//  j = 1;
//  nodeptr v;
//  
//  
//  for (head = nt->tree; head; head = head->next)
//  {
//    item = (pllNewickNodeInfo *) head->item;
//    if (!nodeStack)
//     {
//       pllStackPush (&nodeStack, tr->nodep[i]);
//       pllStackPush (&nodeStack, tr->nodep[i]->next);
//       pllStackPush (&nodeStack, tr->nodep[i]->next->next);
//       ++i;
//     }
//    else
//     {
//       v = (nodeptr) pllStackPop (&nodeStack);
//       if (item->rank)  /* internal node */
//        {
//          v->back           = tr->nodep[i];
//          tr->nodep[i]->back = v; //t->nodep[v->number]
//          pllStackPush (&nodeStack, tr->nodep[i]->next);
//          pllStackPush (&nodeStack, tr->nodep[i]->next->next);
//          double z = exp((-1 * atof(item->branch))/tr->fracchange);
//          if(z < PLL_ZMIN) z = PLL_ZMIN;
//          if(z > PLL_ZMAX) z = PLL_ZMAX;
//          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
//             v->z[k] = tr->nodep[i]->z[k] = z;
//
//          ++ i;
//        }
//       else             /* leaf */
//        {
//          v->back           = tr->nodep[j];
//          tr->nodep[j]->back = v; //t->nodep[v->number];
//
//          double z = exp((-1 * atof(item->branch))/tr->fracchange);
//          if(z < PLL_ZMIN) z = PLL_ZMIN;
//          if(z > PLL_ZMAX) z = PLL_ZMAX;
//          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
//            v->z[k] = tr->nodep[j]->z[k] = z;
//            
//          //t->nameList[j] = strdup (item->name);
//          tr->nameList[j] = (char *) rax_malloc ((strlen (item->name) + 1) * sizeof (char));
//          strcpy (tr->nameList[j], item->name);
//          
//          pllHashAdd (tr->nameHash, tr->nameList[j], (void *) (tr->nodep[j]));
//          ++ j;
//        }
//     }
//  }
//  
//  tr->start = tr->nodep[1];
//  
//  pllStackClear (&nodeStack);
//
//  if (useDefaultz == PLL_TRUE) 
//    resetBranches (tr);
//}

/** @brief Initialize PLL tree with a random topology

    Initializes the PLL tree with a randomly created topology

    @todo
      Perhaps pass a seed?

    @param tr
      The PLL instance

    @param tips
      Number of tips

    @param nameList
      A set of \a tips names representing the taxa labels
*/
void 
pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList)
{
  int i;
  pllTreeInitDefaults (tr, tips);

  for (i = 1; i <= tips; ++i) {
      rax_malloc_string_copy(nameList[i], &(tr->nameList[i]));
      pllHashAdd(tr->nameHash, pllHashString(tr->nameList[i], tr->nameHash->size), tr->nameList[i], (void*)(tr->nodep[i]));
  }

  pllMakeRandomTree (tr);
}


/** @brief Initialize a tree that corresponds to a given (already parsed) alignment 

    Initializes the PLL tree such that it corresponds to the given alignment

    @todo
      nothing 

    @param tr
      The PLL instance

    @param alignmentData
      Parsed alignment
*/
void 
pllTreeInitTopologyForAlignment(pllInstance* tr, pllAlignmentData* alignmentData)
{
    int tips = alignmentData->sequenceCount;
    char** nameList = alignmentData->sequenceLabels;

    pllTreeInitDefaults(tr, tips);

    for (int i = 1; i <= tips; ++i) {
        rax_malloc_string_copy(nameList[i], &(tr->nameList[i]));
        pllHashAdd(tr->nameHash, pllHashString(tr->nameList[i], tr->nameHash->size), tr->nameList[i], (void*)(tr->nodep[i]));
    }
}


/** @brief Compute a randomized stepwise addition oder parsimony tree

    Implements the RAxML randomized stepwise addition order algorithm 

    @todo
      check functions that are invoked for potential memory leaks!

    @param tr
      The PLL instance

    @param partitions
      The partitions

    @param sprDist
      SPR distance for the SPR search in parsimony
*/
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions, int sprDist)
{
  allocateParsimonyDataStructures(tr, partitions);
  pllMakeParsimonyTreeFast(tr, partitions, sprDist);
  pllFreeParsimonyDataStructures(tr, partitions);
}

/** @brief Encode the alignment data to the PLL numerical representation
    
    Transforms the alignment to the PLL internal representation by substituting each base 
    with a specific digit.

    @param alignmentData  Multiple sequence alignment
    @param partitions     List of partitions
*/
void pllBaseSubstitute (pllInstance * tr, partitionList * partitions)
{
  const char * d;
  int i, j, k;

  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     switch (partitions->partitionData[i]->dataType)
      {
        case PLL_DNA_DATA:
          d = PLL_MAP_NT;
          break;
        case PLL_BINARY_DATA:
          d = PLL_MAP_BIN;
          break;
        case PLL_AA_DATA:
          d = PLL_MAP_AA;
          break;
        default:
          assert(0);
      }
     
     for (j = 1; j <= tr->mxtips; ++ j)
      {
        for (k = partitions->partitionData[i]->lower; k < partitions->partitionData[i]->upper; ++ k)
         {
           tr->yVector[j][k] = d[tr->yVector[j][k]];
         }
      }
   }
}

/** Clears the rearrangements history from PLL instance
    
    Clears the rearrangements rollback information (history) from the PLL instance \a tr.

    @param tr
      PLL instance
*/
void pllClearRearrangeHistory (pllInstance * tr)
{
  pllRollbackInfo * ri;

  while ((ri = (pllRollbackInfo *)pllStackPop (&(tr->rearrangeHistory))))
   {
     rax_free (ri);
   }
}

/** @brief Deallocate the PLL instance

    Deallocates the library instance and all its elements.

    @param tr
      The PLL instance
*/
void
pllDestroyInstance (pllInstance * tr)
{
  int i;

  for (i = 1; i <= tr->mxtips; ++ i)
    rax_free (tr->nameList[i]);
  
  pllHashDestroy (&(tr->nameHash), NULL);
  if (tr->yVector)
   {
     if (tr->yVector[0]) rax_free (tr->yVector[0]);
     rax_free (tr->yVector);
   }
  rax_free (tr->aliaswgt);
  rax_free (tr->rateCategory);
  rax_free (tr->patrat);
  rax_free (tr->patratStored);
  rax_free (tr->lhs);
  rax_free (tr->td[0].parameterValues);
  rax_free (tr->td[0].executeModel);
  rax_free (tr->td[0].ti);
  rax_free (tr->nameList);
  rax_free (tr->nodep);
  rax_free (tr->nodeBaseAddress);
  rax_free (tr->tree_string);
  rax_free (tr->tree0);
  rax_free (tr->tree1);
  rax_free (tr->tipNames);
  rax_free (tr->constraintVector);
  pllClearRearrangeHistory (tr);

  rax_free (tr);

#ifdef _FINE_GRAIN_MPI
  pllFinalizeMPI ();
#endif

}

/* initializwe a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
*/



static int init_Q_MatrixSymmetries(char *linkageString, partitionList * pr, int model)
{
  int 
    states = pr->partitionData[model]->states,
    numberOfRates = ((states * states - states) / 2), 
    *list = (int *)rax_malloc(sizeof(int) * numberOfRates),
    j,
    max = -1;

  char
    *str1,
    *saveptr,
    *ch,
    *token;

  rax_malloc_string_copy(linkageString, &ch);

  for(j = 0, str1 = ch; ;j++, str1 = (char *)NULL) 
    {
      token = strtok_r(str1, ",", &saveptr);
      if(token == (char *)NULL)
        break;
      if(!(j < numberOfRates))
        {
          errno = PLL_SUBSTITUTION_RATE_OUT_OF_BOUNDS;
          return PLL_FALSE;
        }
      list[j] = atoi(token);     
    }
  
  rax_free(ch);

  for(j = 0; j < numberOfRates; j++)
    {
      if(!(list[j] <= j))
        {
          errno = PLL_INVALID_Q_MATRIX_SYMMETRY;
          return PLL_FALSE;
        }
      
      if(!(list[j] <= max + 1))
        {
          errno = PLL_Q_MATRIX_SYMMETRY_OUT_OF_BOUNDS;
          return PLL_FALSE;
        }
      
      if(list[j] > max)
        max = list[j];
    }  
  
  for(j = 0; j < numberOfRates; j++)  
    pr->partitionData[model]->symmetryVector[j] = list[j];    

  //less than the maximum possible number of rate parameters

  if(max < numberOfRates - 1)    
    pr->partitionData[model]->nonGTR = PLL_TRUE;

  pr->partitionData[model]->optimizeSubstitutionRates = PLL_TRUE;

  rax_free(list);

  return PLL_TRUE;
}

/** @brief Check parameter linkage across partitions for consistency
 *
 * Checks that linked alpha, substitution rate and frequency model parameters 
 * across several partitions are consistent. E.g., when two partitions are linked 
 * via the alpha parameter, the alpha parameter should either be set to the same 
 * fixed value or it should be estimated!
 *
 * @param pr
 *   List of partitions
 *
 * @todo
 *   Call this in more functions, right now it's only invoked in the wrapper 
 *   for modOpt() 
 */
static int checkLinkageConsistency(partitionList *pr)
{
  if(pr->dirty)
    {
      int 
        i;
      
      linkageList 
        *ll;

      /* first deal with rates */

      ll = pr->rateList;
        
      for(i = 0; i < ll->entries; i++)
        {
          int
            partitions = ll->ld[i].partitions,
            reference = ll->ld[i].partitionList[0];
          
          if(pr->partitionData[reference]->dataType == PLL_AA_DATA)
            {
              if(pr->partitionData[reference]->protModels == PLL_GTR || pr->partitionData[reference]->nonGTR)                             
                {
                  if(!(pr->partitionData[reference]->optimizeSubstitutionRates == PLL_TRUE))
                    {
                      errno = PLL_INCONSISTENT_SUBST_RATE_OPTIMIZATION_SETTING;
                      return PLL_FALSE;
                    }
                }
              else              
                {
                  if(!(pr->partitionData[reference]->optimizeSubstitutionRates == PLL_FALSE))
                    {
                      errno = PLL_INCONSISTENT_SUBST_RATE_OPTIMIZATION_SETTING;
                      return PLL_FALSE;
                    }
                }                 
            }

          if(partitions > 1)
            {
              int
                j,
                k;
              
              for(k = 1; k < partitions; k++)
                {
                  int 
                    index = ll->ld[i].partitionList[k];
                  
                  int
                    states = pr->partitionData[index]->states,
                    rates = ((states * states - states) / 2);
                  
                  if(!(pr->partitionData[reference]->nonGTR == pr->partitionData[index]->nonGTR))
                    {
                      errno = PLL_INCONSISTENT_SUBST_RATE_OPTIMIZATION_SETTING;
                      return PLL_FALSE;
                    }
                  if(!(pr->partitionData[reference]->optimizeSubstitutionRates == pr->partitionData[index]->optimizeSubstitutionRates))
                    {
                      errno = PLL_INCONSISTENT_SUBST_RATE_OPTIMIZATION_SETTING;
                      return PLL_FALSE;
                    }
                
                  
                  if(pr->partitionData[reference]->nonGTR)
                    {              
                      
                      for(j = 0; j < rates; j++)                        
                        {
                          if(!(pr->partitionData[reference]->symmetryVector[j] == pr->partitionData[index]->symmetryVector[j]))
                            {
                              errno = PLL_INCONSISTENT_Q_MATRIX_SYMMETRIES_ACROSS_LINKED_PARTITIONS;
                              return PLL_FALSE;
                            }
                        }
                    }
                  
                 
                  for(j = 0; j < rates; j++)
                    {
                      if(!(pr->partitionData[reference]->substRates[j] == pr->partitionData[index]->substRates[j]))
                        {
                          errno = PLL_INCONSISTENT_Q_MATRIX_ENTRIES_ACROSS_LINKED_PARTITIONS;
                          return PLL_FALSE;
                        }
                    }
                }           
            }
        }
      
      /* then deal with alpha parameters */

      ll = pr->alphaList;

      for(i = 0; i < ll->entries; i++)
        {
          int
            partitions = ll->ld[i].partitions;
          
          if(partitions > 1)
            {
              int
                k, 
                reference = ll->ld[i].partitionList[0];
              
              for(k = 1; k < partitions; k++)
                {
                  int 
                    index = ll->ld[i].partitionList[k];                          

                  if(!(pr->partitionData[reference]->optimizeAlphaParameter == pr->partitionData[index]->optimizeAlphaParameter))
                    {
                      errno = PLL_INCONSISTENT_ALPHA_STATES_ACROSS_LINKED_PARTITIONS;
                      return PLL_FALSE;
                    }
                  if(!(pr->partitionData[reference]->alpha == pr->partitionData[index]->alpha))
                    {
                      errno = PLL_INCONSISTENT_ALPHA_VALUES_ACROSS_LINKED_PARTITIONS;
                      return PLL_FALSE;
                    }
                }           
            }
        }

      /* and then deal with base frequencies */

      ll = pr->freqList;

      for(i = 0; i < ll->entries; i++)
        {
          int     
            partitions = ll->ld[i].partitions;
          
          if(partitions > 1)
            {
              int               
                k, 
                reference = ll->ld[i].partitionList[0];
              
              for(k = 1; k < partitions; k++)
                {
                  int
                    j,
                    index = ll->ld[i].partitionList[k],
                    states = pr->partitionData[index]->states;                           

                  if(!(pr->partitionData[reference]->optimizeBaseFrequencies == pr->partitionData[index]->optimizeBaseFrequencies))
                    {
                      errno = PLL_INCONSISTENT_FREQUENCY_STATES_ACROSS_LINKED_PARTITIONS;
                      return PLL_FALSE;
                    }

                  for(j = 0; j < states; j++)
                    {
                      if(!(pr->partitionData[reference]->frequencies[j] == pr->partitionData[index]->frequencies[j]))
                        {
                          errno = PLL_INCONSISTENT_FREQUENCY_VALUES_ACROSS_LINKED_PARTITIONS;
                          return PLL_FALSE;
                        }
                    }
                }           
            }
        }
      
      pr->dirty = PLL_FALSE;
    }

  return PLL_TRUE;
}
/** @brief Set symmetries among parameters in the Q matrix
    
    Allows to link some or all rate parameters in the Q-matrix 
    for obtaining simpler models than GTR

    @param string
      string describing the symmetry pattern among the rates in the Q matrix

    @param pr
      List of partitions
      
    @param model
      Index of the partition for which we want to set the Q matrix symmetries

    @todo
      nothing
*/
int pllSetSubstitutionRateMatrixSymmetries(char *string, partitionList * pr, int model)
{
  int 
    result = init_Q_MatrixSymmetries(string, pr, model);

  pr->dirty = PLL_TRUE;

  return result;
}

/** @defgroup modelParamsGroup Model parameters setup and retrieval
    
    This set of functions is responsible for setting, retrieving, and optimizing
    model parameters. It also contains functions for linking model parameters
    across partitions.
*/

/** @ingroup modelParamsGroups
    @brief Set the alpha parameter of the Gamma model to a fixed value for a partition
    
    Sets the alpha parameter of the gamma model of rate heterogeneity to a fixed value
    and disables the optimization of this parameter 

    @param alpha
      alpha value

    @param model
      Index of the partition for which we want to set the alpha value

    @param pr
      List of partitions
      
    @param tr
      Library instance for which we want to fix alpha 

    @todo
      test if this works with the parallel versions
*/
void pllSetFixedAlpha(double alpha, int model, partitionList * pr, pllInstance *tr)
{
  //make sure that we are swetting alpha for a partition within the current range 
  //of partitions
  double old_fracchange = tr->fracchange;

  assert(model >= 0 && model < pr->numberOfPartitions);

  assert(alpha >= PLL_ALPHA_MIN && alpha <= PLL_ALPHA_MAX);

  //set the alpha paremeter 
  
  pr->partitionData[model]->alpha = alpha;

  //do the discretization of the gamma curve

  pllMakeGammaCats(pr->partitionData[model]->alpha, pr->partitionData[model]->gammaRates, 4, tr->useMedian);

  //broadcast the changed parameters to all threads/MPI processes 

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier(tr, pr, PLL_THREAD_COPY_ALPHA);
#endif

  pr->partitionData[model]->optimizeAlphaParameter = PLL_FALSE;

  pr->dirty = PLL_FALSE;
  updateAllBranchLengths (tr, old_fracchange, tr->fracchange);
}

/** @ingroup modelParamsGroups
    @brief Get the rate categories of the Gamma model of a partition

    Gets the gamma rate categories of the Gamma model of rate heterogeneity
    of partition \a pid from partition list \a pr.

    @param pr   List of partitions
    @param pid  Index of partition to use
    @param outBuffer  Output buffer where to store the rates
*/
void pllGetGammaRates (partitionList * pr, int pid, double * outBuffer)
{
  /* TODO: Change the hardcoded 4 and also add a check that this partition
     really uses gamma. Currently, instance is also not required */
  memcpy (outBuffer, pr->partitionData[pid]->gammaRates, 4 * sizeof (double));
}

/** @ingroup modelParamsGroups
    @brief Get the alpha parameter of the Gamma model of a partition

    Returns the alpha parameter of the gamma model of rate heterogeneity
    of partition \a pid from partition list \a pr.

    @param pr   List of partitions
    @param pid  Index of partition to use

    @return
      Alpha parameter
*/
double pllGetAlpha (partitionList * pr, int pid)
{
  /* TODO: check if the partition uses gamma */
  return (pr->partitionData[pid]->alpha);
}


/** @ingroup modelParamsGroups
    @brief Get the base frequencies of a partition

    Gets the base frequencies of partition \a model from partition list
    \a partitionList and stores them in \a outBuffer. Note that \outBuffer
    must be of size s, where s is the number of states.

    @param  tr       PLL instance
    @param pr        List of partitions
    @param model     Index of the partition for which we want to get the base frequencies
    @param outBuffer Buffer where to store the base frequencies
*/
void pllGetBaseFrequencies(partitionList * pr, int model, double * outBuffer)
{
  memcpy (outBuffer, pr->partitionData[model]->frequencies, pr->partitionData[model]->states * sizeof (double));
}


/** @ingroup modelParamsGroups
    @brief Set all base frequencies to a fixed value for a partition
    
    Sets all base freuqencies of a partition to fixed values and disables 
    ML optimization of these parameters 

    @param f
      array containing the base frequencies

    @param  length
      length of array f, this needs to be as long as the number of 
      states in the model, otherwise an assertion will fail!

    @param model
      Index of the partition for which we want to set the frequencies 

    @param pr
      List of partitions
      
    @param tr
      Library instance for which we want to fix the base frequencies

    @todo
      test if this works with the parallel versions
*/
void pllSetFixedBaseFrequencies(double *f, int length, int model, partitionList * pr, pllInstance *tr)
{
  int 
    i;

  double 
    acc = 0.0,
    old_fracchange;

  old_fracchange = tr->fracchange;

  //make sure that we are setting the base frequencies for a partition within the current range 
  //of partitions
  assert(model >= 0 && model < pr->numberOfPartitions);

  //make sure that the length of the input array f containing the frequencies 
  //is as long as the number of states in the model 
  assert(length == pr->partitionData[model]->states);


  //make sure that the base frequencies sum approximately to 1.0
  
  for(i = 0; i < length; i++)
    acc += f[i];

  if(fabs(acc - 1.0) > 0.000001)
    assert(0);

  //copy the base frequencies 
  memcpy(pr->partitionData[model]->frequencies, f, sizeof(double) * length);

  //re-calculate the Q matrix 
  pllInitReversibleGTR(tr, pr, model);


  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  
  pr->partitionData[model]->optimizeBaseFrequencies = PLL_FALSE;

  pr->dirty = PLL_TRUE;
  updateAllBranchLengths (tr, old_fracchange, tr->fracchange);
}

/** @ingroup modelParamsGroups
    @brief Set that the base freuqencies are optimized under ML
    
    The base freuqencies for partition model will be optimized under ML    

    @param model
      Index of the partition for which we want to optimize base frequencies 

    @param pr
      List of partitions
      
    @param tr
      Library instance for which we want to fix the base frequencies

    @todo
      test if this works with the parallel versions
*/
int pllSetOptimizeBaseFrequencies(int model, partitionList * pr, pllInstance *tr)
{
  int
    states,
    i;

  double 
    initialFrequency,
    acc = 0.0;

  //make sure that we are setting the base frequencies for a partition within the current range 
  //of partitions
  if(!(model >= 0 && model < pr->numberOfPartitions))
    {
      errno = PLL_PARTITION_OUT_OF_BOUNDS;
      return PLL_FALSE;
    }

  //set the number of states/ferquencies in this partition 
  states = pr->partitionData[model]->states;

  //set all frequencies to 1/states
  
  initialFrequency = 1.0 / (double)states;

  for(i = 0; i < states; i++)
    pr->partitionData[model]->frequencies[i] = initialFrequency;

  //make sure that the base frequencies sum approximately to 1.0
  
  for(i = 0; i < states; i++)
    acc += pr->partitionData[model]->frequencies[i];

  if(fabs(acc - 1.0) > 0.000001)
    {
      errno = PLL_BASE_FREQUENCIES_DO_NOT_SUM_TO_1;
      return PLL_FALSE;
    }

  //re-calculate the Q matrix 
  pllInitReversibleGTR(tr, pr, model);

  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  
  pr->partitionData[model]->optimizeBaseFrequencies = PLL_TRUE;

  pr->dirty = PLL_TRUE;

  return PLL_TRUE;
}




/** @ingroup modelParamsGroups
    @brief Get the substitution rates for a specific partition

    Gets the substitution rates of partition \a model from partition list
    \a partitionList and stores them in \a outBuffer. Note that \outBuffer
    must be of size (2 * s - s) / 2, where s is the number of states, i.e.
    the number of upper diagonal entries of the Q matrix.

    @param tr        PLL instance
    @param pr        List of partitions
    @param model     Index of partition for which we want to get the substitution rates
    @param outBuffer Buffer where to store the substitution rates.
*/
void pllGetSubstitutionMatrix (partitionList * pr, int model, double * outBuffer)
{
  int 
    rates,
    states;
  
  states = pr->partitionData[model]->states;
  rates = (states * states - states) / 2;

  memcpy (outBuffer, pr->partitionData[model]->substRates, rates * sizeof (double));
}

/** @ingroup modelParamsGroups
     @brief Set all substitution rates for a specific partition and disable ML optimization for them
    
    Sets all substitution rates of a partition to fixed values and disables 
    ML optimization of these parameters. It will automatically re-scale the relative rates  
    such that the last rate is 1.0 

    @param f
      array containing the substitution rates

    @param length
      length of array f, this needs to be as long as: (s * s - s) / 2,
      i.e., the number of upper diagonal entries of the Q matrix

    @param model
      Index of the partition for which we want to set/fix the substitution rates

    @param pr
      List of partitions
      
    @param tr
      Library instance for which we want to fix the substitution rates 

    @todo
      test if this works with the parallel versions
*/
void pllSetFixedSubstitutionMatrix(double *q, int length, int model, partitionList * pr,  pllInstance *tr)
{
  pllSetSubstitutionMatrix(q, length, model, pr, tr);
  pr->partitionData[model]->optimizeSubstitutionRates = PLL_FALSE;
}

/** @ingroup modelParamsGroups
     @brief Set all substitution rates for a specific partition
    
    Sets all substitution rates of a partition to the given values.
    It will automatically re-scale the relative rates such that the last rate is 1.0 

    @param f
      array containing the substitution rates

    @param length
      length of array f, this needs to be as long as: (s * s - s) / 2,
      i.e., the number of upper diagonal entries of the Q matrix

    @param model
      Index of the partition for which we want to set/fix the substitution rates

    @param pr
      List of partitions
      
    @param tr
      Library instance for which we want to fix the substitution rates 

    @todo
      test if this works with the parallel versions
*/
void pllSetSubstitutionMatrix(double *q, int length, int model, partitionList * pr,  pllInstance *tr)
{
  int 
    i,
    numberOfRates; 

  double
    scaler,
    old_fracchange;

  old_fracchange = tr->fracchange;

  //make sure that we are setting the Q matrix for a partition within the current range 
  //of partitions
  assert(model >= 0 && model < pr->numberOfPartitions);

  numberOfRates = (pr->partitionData[model]->states * pr->partitionData[model]->states - pr->partitionData[model]->states) / 2;

  //  make sure that the length of the array containing the subsitution rates 
  //  corresponds to the number of states in the model

  assert(length == numberOfRates);

  //automatically scale the last rate to 1.0 if this is not already the case

  if(q[length - 1] != 1.0)    
    scaler = 1.0 / q[length - 1]; 
  else
    scaler = 1.0;

  //set the rates for the partition and make sure that they are within the allowed bounds 

  for(i = 0; i < length; i++)
    {
      double
        r = q[i] * scaler;
      
      assert(r >= PLL_RATE_MIN && r <= PLL_RATE_MAX);
      
      pr->partitionData[model]->substRates[i] = r;
    }

  //re-calculate the Q matrix 
  pllInitReversibleGTR(tr, pr, model);

  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  

  pr->dirty = PLL_TRUE;
  updateAllBranchLengths (tr, old_fracchange, tr->fracchange);
}




/* initialize a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
*/

/** @ingroup modelParamsGroups
*/
linkageList* initLinkageList(int *linkList, partitionList *pr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  
  linkageList 
    *ll = (linkageList*)rax_malloc(sizeof(linkageList));
    
  /* figure out how many distinct parameters we need to estimate 
     in total, if all parameters are linked the result will be 1 if all 
     are unlinked the result will be pr->numberOfPartitions */
  
  for(i = 0; i < pr->numberOfPartitions; i++)
    {
      if(!(linkList[i] >= 0 && linkList[i] < pr->numberOfPartitions))
        {
          errno = PLL_LINKAGE_LIST_OUT_OF_BOUNDS;
          return (linkageList*)NULL;
        }

      if(!(linkList[i] <= i && linkList[i] <= numberOfModels + 1))
        {
          errno = PLL_LINKAGE_LIST_OUT_OF_BOUNDS;
          return (linkageList*)NULL;
        }

      if(linkList[i] > numberOfModels)
        numberOfModels = linkList[i];

    }

  numberOfModels++;
  
  /* allocate the linkage list data structure that containes information which parameters of which partition are 
     linked with each other.

     Note that we need a separate invocation of initLinkageList() and a separate linkage list 
     for each parameter type */

  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)rax_malloc(sizeof(linkageData) * numberOfModels);

  /* noe loop over the number of free parameters and assign the corresponding partitions to each parameter */

  for(i = 0; i < numberOfModels; i++)
    {
      /* 
         the valid flag is used for distinguishing between DNA and protein data partitions.
         This can be used to enable/disable parameter optimization for the paremeter 
         associated to the corresponding partitions. This deature is used in optRatesGeneric 
         to first optimize all DNA GTR rate matrices and then all PROT GTR rate matrices */

      ll->ld[i].valid = PLL_TRUE;
      partitions = 0;

      /* now figure out how many partitions share this joint parameter */

      for(k = 0; k < pr->numberOfPartitions; k++)
        if(linkList[k] == i)
          partitions++;     

      /* assign a list to store the partitions that share the parameter */

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)rax_malloc(sizeof(int) * partitions);
      
      /* now store the respective partition indices in this list */
      
      for(k = 0, pos = 0; k < pr->numberOfPartitions; k++)
        if(linkList[k] == i)
          ll->ld[i].partitionList[pos++] = k;
    }

  /* return the linkage list for the parameter */

  return ll;
}



static linkageList* initLinkageListString(char* linkageString, partitionList* pr)
{
    int* list = (int*)rax_malloc(sizeof(int) * pr->numberOfPartitions);
    char* saveptr, * ch;

    rax_malloc_string_copy(linkageString, &ch);
    char* str1 = ch;
    for (int j = 0; ; j++, str1 = (char*)NULL) {
        char* token = strtok_r(str1, ",", &saveptr);
        if (token == (char*)NULL) {
            break;
        }
        assert(j < pr->numberOfPartitions);
        list[j] = atoi(token);
    }
    rax_free(ch);
    linkageList* l = initLinkageList(list, pr);
    rax_free(list);

    return l;
}

/** @ingroup modelParamsGroups
    @brief Link alpha parameters across partitions
    
    Links alpha paremeters across partitions (GAMMA model of rate heterogeneity)

    @param string
      string describing the linkage pattern    

    @param pr
      List of partitions

    @todo
      test behavior/impact/mem-leaks of this when PSR model is used 
      it shouldn't do any harm, but it would be better to check!
*/
int pllLinkAlphaParameters(char *string, partitionList *pr)
{
  //assumes that it has already been assigned once
  freeLinkageList(pr->alphaList);
  
  pr->alphaList = initLinkageListString(string, pr); 

  pr->dirty = PLL_TRUE;
  
  if(!pr->alphaList)
    return PLL_FALSE;
  else
    return PLL_TRUE;
}

/** @ingroup modelParamsGroups
    @brief Link base frequency parameters across partitions
    
    Links base frequency paremeters across partitions

    @param string
      string describing the linkage pattern    

    @param pr
      List of partitions

    @todo
      semantics of this function not clear yet: right now this only has an effect 
      when we do a ML estimate of base frequencies 
      when we use empirical or model-defined (protein data) base frequencies, one could 
      maybe average over the per-partition frequencies, but the averages would need to be weighted 
      accodring on the number of patterns per partition 
*/
int pllLinkFrequencies(char *string, partitionList *pr)
{
  //assumes that it has already been assigned once
  freeLinkageList(pr->freqList);

  pr->freqList = initLinkageListString(string, pr);

  pr->dirty = PLL_TRUE;

  if(!pr->freqList)
    return PLL_FALSE;
  else
    return PLL_TRUE;
}

/** @ingroup modelParamsGroups
    @brief Link Substitution matrices across partitions
    
    Links substitution matrices (Q matrices) across partitions

    @param string
      string describing the linkage pattern    

    @param pr
      List of partitions

    @todo
      re-think/re-design how this is done for protein
      models
*/
int pllLinkRates(char *string, partitionList *pr)
{
  //assumes that it has already been assigned once
  freeLinkageList(pr->rateList);
  
  pr->rateList = initLinkageListString(string, pr);
  
  pr->dirty = PLL_TRUE;  

  if(!pr->dirty)
    return PLL_FALSE;
  else
    return PLL_TRUE;
}




/** @ingroup modelParamsGroups
    @brief Initialize partitions according to model parameters
    
    Initializes partitions according to model parameters.

    @param tr              The PLL instance
    @param partitions      List of partitions
    @param alignmentData   The parsed alignment
    @return                Returns \b PLL_TRUE in case of success, otherwise \b PLL_FALSE
*/
int pllInitModel (pllInstance * tr, partitionList * partitions) 
{
  double ** ef;
  int
    i,
    *unlinked = (int *)rax_malloc(sizeof(int) * partitions->numberOfPartitions);
  double old_fracchange = tr->fracchange;

  ef = pllBaseFrequenciesInstance (tr, partitions);

  if(!ef)
    return PLL_FALSE;

  
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#if (defined(__AVX) || defined(__SSE3))
#if !defined(__ARM_NEON) // flush zero always on in NEON
  _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif
#endif
#endif 

#ifdef _USE_PTHREADS
  tr->threadID = 0;
#ifndef _PORTABLE_PTHREADS
  /* not very portable thread to core pinning if PORTABLE_PTHREADS is not defined
     by defualt the cod ebelow is deactivated */
  pinToCore(0);
#endif
#endif

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* 
     this main function is the master thread, so if we want to run RAxML with n threads,
     we use pllStartPthreads to start the n-1 worker threads */
  
#ifdef _USE_PTHREADS
  pllStartPthreads (tr, partitions);
#endif

  /* via pllMasterBarrier() we invoke parallel regions in which all Pthreads work on computing something, mostly likelihood 
     computations. Have a look at execFunction() in axml.c where we siwtch of the different types of parallel regions.

     Although not necessary, below we copy the info stored on tr->partitionData to corresponding copies in each thread.
     While this is shared memory and we don't really need to copy stuff, it was implemented like this to allow for an easier 
     transition to a distributed memory implementation (MPI).
     */
#ifdef _FINE_GRAIN_MPI
  //MPI_Bcast (&(partitions->numberOfPartitions), 1, MPI_INT, MPI_ROOT, MPI_COMM_WORLD);
  MPI_Bcast (&(partitions->numberOfPartitions), 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  /* mpi version now also uses the generic barrier */
  pllMasterBarrier (tr, partitions, PLL_THREAD_INIT_PARTITION);
#else  /* SEQUENTIAL */
  /* 
     allocate the required data structures for storing likelihood vectors etc 
     */

  //initializePartitions(tr, tr, partitions, partitions, 0, 0);
  initializePartitionsSequential (tr, partitions);
#endif
  
  //initializePartitions (tr, tr, partitions, partitions, 0, 0);
  
  initModel (tr, ef, partitions);

  pllEmpiricalFrequenciesDestroy (&ef, partitions->numberOfPartitions);

  for(i = 0; i < partitions->numberOfPartitions; i++)
    unlinked[i] = i;

  //by default everything is unlinked initially 
  partitions->alphaList = initLinkageList(unlinked, partitions);
  partitions->freqList  = initLinkageList(unlinked, partitions);
  partitions->rateList  = initLinkageList(unlinked, partitions);

  rax_free(unlinked);

  updateAllBranchLengths (tr, old_fracchange ? old_fracchange : 1,  tr->fracchange);
  pllEvaluateLikelihood (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

  return PLL_TRUE;
}
 
/** @ingroup modelParamsGroups
    @brief Optimize all free model parameters of the likelihood model
    
    Initializes partitions according to model parameters.

    @param tr
      The PLL instance

    @param pr
      List of partitions

    @param likelihoodEpsilon
      Specifies up to which epsilon in likelihood values the iterative routine will 
      be optimizing the parameters  
*/
int pllOptimizeModelParameters(pllInstance *tr, partitionList *pr, double likelihoodEpsilon)
{
  //force the consistency check

  pr->dirty = PLL_TRUE;

  if(!checkLinkageConsistency(pr))
    return PLL_FALSE;

  modOpt(tr, pr, likelihoodEpsilon);

  return PLL_TRUE;
}

/** @brief Read the contents of a file
    
    Reads the ile \a filename and return its content. In addition
    the size of the file is stored in the input variable \a filesize.
    The content of the variable \a filesize can be anything and will
    be overwritten.

    @param filename
      Name of the input file

    @param filesize
      Input parameter where the size of the file (in bytes) will be stored

    @return
      Contents of the file
*/
char * 
pllReadFile (const char * filename, long * filesize)
{
  FILE * fp;
  char * rawdata;

  // FIX BUG: opening with "r" does not work on Windows
//  fp = fopen (filename, "r");
  printf("[PLL] Reading file %s...\n", filename);
  fp = fopen (filename, "rb");
  printf("[PLL] Success!\n");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1)
   {
     fclose (fp);
     return (NULL);
   }

  *filesize = ftell (fp);

  if (*filesize == -1) 
   {
     fclose (fp);
     return (NULL);
   }
  rewind (fp);

  /* allocate buffer and read file contents */
  rawdata = (char *) rax_malloc (((*filesize) + 1) * sizeof (char));
  if (rawdata) 
   {
     if (fread (rawdata, sizeof (char), *filesize, fp) != (size_t) *filesize) 
      {
        rax_free (rawdata);
        rawdata = NULL;
      }
     else
      {
        rawdata[*filesize] = 0;
      }
   }

  fclose (fp);

  return (rawdata);
}

static void getInnerBranchEndPointsRecursive (nodeptr p, int tips, int * i, node **nodes)
{
  if (!isTip (p->next->back->number, tips))
   {
     nodes[(*i)++] = p->next;
     getInnerBranchEndPointsRecursive(p->next->back, tips, i, nodes);
   }
  if (!isTip (p->next->next->back->number, tips))
   {
     nodes[(*i)++] = p->next->next;
     getInnerBranchEndPointsRecursive(p->next->next->back, tips, i, nodes);
   }
}

node ** pllGetInnerBranchEndPoints (pllInstance * tr)
{
  node ** nodes;
  nodeptr p;
  int i = 0;

  nodes = (node **) rax_calloc(tr->mxtips - 3, sizeof(node *));

  p = tr->start;
  assert (isTip(p->number, tr->mxtips));

  getInnerBranchEndPointsRecursive(p->back, tr->mxtips, &i, nodes);

  return nodes;
}

#if defined WIN32 || defined _WIN32 || defined __WIN32__
void* rax_calloc(size_t count, size_t size) {
	void* res = rax_malloc(size * count);
	memset(res,0,size * count);
	return res;
}
#endif

