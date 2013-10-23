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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

/** @file utils.c
 *  
 *  @brief Miscellaneous general utility and helper functions
 */
#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

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


#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#if (defined(__AVX) || defined(__SSE3))
#include <xmmintrin.h>
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

#define GLOBAL_VARIABLES_DEFINITION

#include "globalVariables.h"

static void pllTreeInitDefaults (pllInstance * tr, int tips);
static void initializePartitionsSequential(pllInstance *tr, partitionList *pr);


/** @defgroup instanceLinkingGroup Linking topology, partition scheme and alignment to the PLL instance
    
    This set of functions handles the linking of topology, partition scheme and multiple sequence alignment
    with the PLL instance
*/
/***************** UTILITY FUNCTIONS **************************/

static char *
my_strtok_r (char * s, const char * delim, char **save_ptr)
{  
  char *token;
   
  /* Scan leading delimiters */
  if (s == NULL)
    s = *save_ptr;
   
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

#if (defined(_SVID_SOURCE) || defined(_BSD_SOURCE) || defined(_POSIX_C_SOURCE) || defined(_XOPEN_SOURCE) || defined(_POSIX_SOURCE))
#define STRTOK_R strtok_r
#else
#define STRTOK_R my_strtok_r
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

void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}

/* Marked for deletion 
boolean getSmoothFreqs(int dataType)
{
  assert(PLL_MIN_MODEL < dataType && dataType < PLL_MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}
*/

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

/* Marked for deletion 
static int filexists(char *filename)
{
  FILE 
    *fp = fopen(filename,"rb");

  int res; 

  if(fp)
  {
    res = 1;
    fclose(fp);
  }
  else
    res = 0;

  return res;
}
*/


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
  {
    if(fp)
      return fp;
    else
    {	  
      printf("\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n", path);
      exit(-1);
      return (FILE *)NULL;
    }
  }
  else
  {
    if(fp)
      return fp;
    else
    {	 
      printf("\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
          path, mode);
      exit(-1);
      return (FILE *)NULL;
    }
  }


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
boolean isTip(int number, int maxTips)
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
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < PLL_NUM_BRANCHES; i++)
    p->z[i] = q->z[i] = PLL_DEFAULTZ;

}


/***********************reading and initializing input ******************/



boolean whitechar (int ch)
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
/*********************************** *********************************************************/








/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/



/* Delete it at some point */
void printLog(pllInstance *tr)
{
  FILE *logFile;
  double t;


  t = gettime() - masterTime;

  logFile = myfopen(logFileName, "ab");

  fprintf(logFile, "%f %f\n", t, tr->likelihood);

  fclose(logFile);


}


/************************************************************************************/

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

      pllNewviewGenericAncestral(tr, pr, p);

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
      
       if(localPartitions->partitionData[model]->dataType == PLL_AA_DATA && localPartitions->partitionData[model]->protModels == PLL_LG4)      
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


      /* 
	 Initializing the xVector array like this is absolutely required !!!!
	 I don't know which programming genious removed this, but it must absolutely stay in here!!!!
      */
      
      {
	int k;
	
	for(k = 0; k < localTree->mxtips; k++)
	      localPartitions->partitionData[model]->xVector[k] = (double*)NULL;       
      }


      localPartitions->partitionData[model]->xSpaceVector = (size_t *)rax_calloc((size_t)localTree->mxtips, sizeof(size_t));

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->sumBuffer), PLL_BYTE_ALIGNMENT, width *
										      (size_t)(localPartitions->partitionData[model]->states) *
										      discreteRateCategories(localTree->rateHetModel) *
										      sizeof(double));

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

      rax_posix_memalign ((void **)&(localPartitions->partitionData[model]->wgt), PLL_BYTE_ALIGNMENT,width * sizeof(int));

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

static void initializePartitionsSequential(pllInstance *tr, partitionList *pr)
{ 
  size_t
    model;

  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    assert(pr->partitionData[model]->width == pr->partitionData[model]->upper - pr->partitionData[model]->lower);

  initializePartitionData(tr, pr);

  /* figure in tip sequence data per-site pattern weights */ 
  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
  {
    size_t
      j;
    size_t lower = pr->partitionData[model]->lower;
    size_t width = pr->partitionData[model]->upper - lower;

    for(j = 1; j <= (size_t)tr->mxtips; j++)
    {
      pr->partitionData[model]->yVector[j] = &(tr->yVector[j][pr->partitionData[model]->lower]);
    }

    memcpy((void*)(&(pr->partitionData[model]->wgt[0])),         (void*)(&(tr->aliaswgt[lower])),      sizeof(int) * width);
  }  

  initMemorySavingAndRecom(tr, pr);
}


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

  /* boolean array for marking that a site was assigned a partition */
  used = (char *) rax_calloc (alignmentData->sequenceLength, sizeof (char));

  /* traverse all partitions and their respective regions and mark sites */
  for (elm = parts->head; elm; elm = elm->next)
   {
     pi = (pllPartitionInfo *) elm->item;
     
     for (regionItem = pi->regionList->head; regionItem; regionItem = regionItem->next)
      {
        region = (pllPartitionRegion *) regionItem->item;

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
static inline void
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

    @param bounds
      An array of the new starting and ending posititons of sites in the alignment for each partition.
      This array is of size 2 * \a nparts. The elements are always couples (lower,upper). The upper
      bounds is a site that is not included in the partition

    @todo
      Fix the bug in PLL 

    @param nparts
      The number of partitions to be created
      
*/
static partitionList *
createPartitions (pllQueue * parts, int * bounds)
{
  partitionList * pl;
  pllPartitionInfo * pi;
  struct pllQueueItem * elm;
  int i;

  pl = (partitionList *) rax_malloc (sizeof (partitionList));
  
  // TODO: fix this
  pl->perGeneBranchLengths =      0;

  // TODO: change PLL_NUM_BRANCHES to number of partitions I guess
  pl->partitionData = (pInfo **) rax_calloc (PLL_NUM_BRANCHES, sizeof (pInfo *));
  
  for (i = 0, elm = parts->head; elm; elm = elm->next, ++ i)
   {
     pi = (pllPartitionInfo *) elm->item;
     pl->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));

     pl->partitionData[i]->lower = bounds[i << 1];
     pl->partitionData[i]->upper = bounds[(i << 1) + 1];
     pl->partitionData[i]->width = bounds[(i << 1) + 1] - bounds[i << 1];

     //the two flags below are required to allow users to set 
     //alpha parameters and substitution rates in the Q matrix 
     //to fixed values. These parameters will then not be optimized 
     //in the model parameter optimization functions
     //by default we assume that all parameters are being optimized, i.e., 
     //this has to be explicitly set by the user 
     
     pl->partitionData[i]->optimizeAlphaParameter    = PLL_TRUE;
     pl->partitionData[i]->optimizeSubstitutionRates = PLL_TRUE;

     if (pi->dataType == PLL_DNA_DATA)
      {
        pl->partitionData[i]->protModels                = -1;
        pl->partitionData[i]->protFreqs                 = -1;
        pl->partitionData[i]->dataType                  = PLL_DNA_DATA;
        pl->partitionData[i]->maxTipStates              = 16;
        pl->partitionData[i]->optimizeBaseFrequencies   = pi->optimizeBaseFrequencies;
      }
     else /* PLL_AA_DATA */
      {
        pl->partitionData[i]->dataType                  = PLL_AA_DATA; 
        pl->partitionData[i]->protModels                = pi->protModels;
	if(pl->partitionData[i]->protModels != PLL_GTR)
	  pl->partitionData[i]->optimizeSubstitutionRates = PLL_FALSE;
        pl->partitionData[i]->maxTipStates              = 23;
        pl->partitionData[i]->protFreqs                 = pi->protFreqs;
        pl->partitionData[i]->protModels                = pi->protModels;
        pl->partitionData[i]->optimizeBaseFrequencies   = pi->optimizeBaseFrequencies;
      }
     
     pl->partitionData[i]->states                = pLengths[pl->partitionData[i]->dataType].states;
     pl->partitionData[i]->numberOfCategories    =        1;
     pl->partitionData[i]->autoProtModels        =        0;
     pl->partitionData[i]->nonGTR                =        PLL_FALSE;
     pl->partitionData[i]->partitionContribution =     -1.0;
     pl->partitionData[i]->partitionLH           =      0.0;
     pl->partitionData[i]->fracchange            =      1.0;
     pl->partitionData[i]->executeModel          =     PLL_TRUE;


     pl->partitionData[i]->partitionName         = (char *) rax_malloc ((strlen (pi->partitionName) + 1) * sizeof (char));
     strcpy (pl->partitionData[i]->partitionName, pi->partitionName);
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

        for (i = region->start - 1; i < region->end; i += region->stride)
         {
           if (oi[i] == i)
            {
              swapSite (alignmentData->sequenceData, dst, i, alignmentData->sequenceCount);
              oi[dst++] = i;
            }
           else if (oi[i] < i)
            {
              j = oi[i];
              while (j < i) j = oi[j];

              swapSite (alignmentData->sequenceData, dst, j, alignmentData->sequenceCount);
              oi[dst++] = j;
            }
         }
      }
     newBounds[(k << 1) + 1] = dst;    /* set the uppwer limit for this partition */
   }
  pl = createPartitions (parts, newBounds);
  pl->numberOfPartitions = nparts;
  pl->dirty = PLL_FALSE;
  
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
static inline void
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
     sites[p]  = (char **) rax_malloc (pl->partitionData[p]->width * sizeof (char *));
     memptr[p] = rax_malloc ((alignmentData->sequenceCount + 1) * pl->partitionData[p]->width * sizeof (char));

     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        sites[p][i] = (char *) (memptr[p] + i * (alignmentData->sequenceCount + 1) * sizeof (char));
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
     alignmentData->sequenceData[i + 1] = (unsigned char *) (alignmentData->sequenceData[0] + i * (alignmentData->sequenceLength + 1) * sizeof (unsigned char));
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
static void
genericBaseFrequencies (pInfo * partition, pllAlignmentData * alignmentData, boolean smoothFrequencies, const unsigned int * bitMask, double * pfreqs)
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
	      
      for (i = 1; i <= alignmentData->sequenceCount; i++) 
	{		 
          yptr = alignmentData->sequenceData[i];
	  
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
      boolean 
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

/** @brief Compute the empirical base frequencies for all partitions

    Compute the empirical base frequencies for all partitions in the list \a pl.

    @param pl
      Partition list
    
    @param alignmentData
      Multiple sequence alignment

    @return
      A list of \a pl->numberOfPartitions arrays each of size \a pl->partitionData[i]->states,
      where \a i is the \a i-th partition
*/
double **
pllBaseFrequenciesGTR (partitionList * pl, pllAlignmentData * alignmentData)
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
          genericBaseFrequencies (pl->partitionData[model], 
                                  alignmentData, 
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

/** @ingroup instanceLinkingGroup
    @brief Load alignment to the PLL instance
    
    Loads (copies) the parsed alignment to the PLL instance. Depending
    on how the \a bDeep flag is set, the alignment in the PLL instance
    is a deep or shallow copy of \a alignmentData

    @param tr
      The library instance

    @param alignmentData 
      The multiple sequence alignment

    @param partitions
      List of partitions

    @param bDeep
      Controls how the alignment is used within the PLL instance.
      If it is set to \b PLL_DEEP_COPY, then new memory will be allocated
      and the alignment will be copied there (deep copy). On the other
      hand, if \b PLL_SHALLOW_COPY is specified, only the pointers will be
      copied and therefore, the alignment will be shared between the 
      alignment structure and the library instance (shallow copy).

    @return
      Returns 1 in case of success, 0 otherwise.
*/
int
pllLoadAlignment (pllInstance * tr, pllAlignmentData * alignmentData, partitionList * partitions, int bDeep)
{
  int i;
  nodeptr node;

  if (tr->mxtips != alignmentData->sequenceCount) return (0);

  /* Do the base substitution (from A,C,G....  ->   0,1,2,3....)*/
  pllBaseSubstitute (alignmentData, partitions);

  tr->aliaswgt = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  memcpy (tr->aliaswgt, alignmentData->siteWeights, alignmentData->sequenceLength * sizeof (int));

  tr->originalCrunchedLength = alignmentData->sequenceLength;
  tr->rateCategory           = (int *)   rax_calloc (tr->originalCrunchedLength, sizeof (int));
  tr->patrat                 = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored           = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->lhs                    = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));

  /* allocate memory for the alignment */
  tr->yVector    = (unsigned char **) rax_malloc ((alignmentData->sequenceCount + 1) * sizeof (unsigned char *));                                                                                                                                                                      
  tr->bDeep = bDeep;

  if (bDeep == PLL_DEEP_COPY)
   {
     tr->yVector[0] = (unsigned char *)  rax_malloc (sizeof (unsigned char) * (alignmentData->sequenceLength + 1) * alignmentData->sequenceCount);
     for (i = 1; i <= alignmentData->sequenceCount; ++ i) 
      {                     
        tr->yVector[i] = (unsigned char *) (tr->yVector[0] + (i - 1) * (alignmentData->sequenceLength + 1) * sizeof (unsigned char));
        tr->yVector[i][alignmentData->sequenceLength] = 0;
      }                     
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
        if (bDeep == PLL_DEEP_COPY) rax_free (tr->yVector[0]);
        rax_free (tr->yVector);
        return (0);
      }
     if (bDeep == PLL_DEEP_COPY)
       memcpy (tr->yVector[node->number], alignmentData->sequenceData[i], alignmentData->sequenceLength );
     else
       tr->yVector[node->number] = alignmentData->sequenceData[i];
   }

  return (1);
}

/** @brief Create the main instance of PLL
    
    Create an instance of the phylogenetic likelihood library

    @param rateHetModel
      Rate heterogeneity model

    @param fastScaling
      explain fastScaling here

    @param saveMemory
      explain saveMemory here

    @param useRecom
      If set to \b PLL_TRUE, enables ancestral state recomputation
    
    @todo
      Document fastScaling, rate heterogeneity and saveMemory and useRecom

    @note
      Do not set \a saveMemory to when using \a useRecom as memory saving techniques 
      are not yet implemented for ancestral state recomputation. 

    @return
      On success returns an instance to PLL, otherwise \b NULL
*/
pllInstance *
pllCreateInstance (pllInstanceAttr * attr)
{
  pllInstance * tr;

  if (attr->rateHetModel != PLL_GAMMA && attr->rateHetModel != PLL_CAT) return NULL;

  tr = (pllInstance *) rax_calloc (1, sizeof (pllInstance));

  tr->threadID     = 0;
  tr->rateHetModel = attr->rateHetModel;
  tr->fastScaling  = attr->fastScaling;
  tr->saveMemory   = attr->saveMemory;
  tr->useRecom     = attr->useRecom;
  
  tr->randomNumberSeed = attr->randomNumberSeed;

  /* remove it from the library */
  tr->useMedian    = PLL_FALSE;

  tr->maxCategories = (attr->rateHetModel == PLL_GAMMA) ? 4 : 25;

  tr->numberOfThreads  = attr->numberOfThreads;
  tr->rearrangeHistory = NULL;

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

  
  p0 = (nodeptr) rax_malloc ((tips + 3 * inner) * sizeof (node));
  assert (p0);

  tr->nodeBaseAddress  = p0;

  tr->nameList         = (char **)   rax_malloc ((tips + 1) * sizeof (char *));
  tr->nodep            = (nodeptr *) rax_malloc ((2 * tips) * sizeof (nodeptr));
  assert (tr->nameList && tr->nodep);

  tr->nodep[0] = NULL;          


  /* TODO: FIX THIS! */
  //tr->fracchange = -1;

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
  tr->h          = NULL;

  /* TODO: Fix hash type */
  tr->nameHash   = pllHashInit (10 * tr->mxtips);

  /* TODO: do these options really fit here or should they be put elsewhere? */
  tr->td[0].count            = 0;
  tr->td[0].ti               = (traversalInfo *) rax_malloc (sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].parameterValues  = (double *) rax_malloc(sizeof(double) * (size_t)PLL_NUM_BRANCHES);
  tr->td[0].executeModel     = (boolean *) rax_malloc (sizeof(boolean) * (size_t)PLL_NUM_BRANCHES);
  tr->td[0].executeModel[0]  = PLL_TRUE;                                                                                                                                                                                                                                    
  for (i = 0; i < PLL_NUM_BRANCHES; ++ i) tr->td[0].executeModel[i] = PLL_TRUE;
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
pllTreeInitTopologyNewick (pllInstance * tr, pllNewickTree * nt, int useDefaultz)
{
  pllStack * nodeStack = NULL;
  pllStack * head;
  struct item_t * item;
  int i, j, k;
  
/*
  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     partitions->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));
     partitions->partitionData[i]->partitionContribution = -1.0;
     partitions->partitionData[i]->partitionLH           =  0.0;
     partitions->partitionData[i]->fracchange            =  1.0;
   }
*/
  
  pllTreeInitDefaults (tr, nt->tips);

  i = nt->tips + 1;
  j = 1;
  nodeptr v;
  
  
  for (head = nt->tree; head; head = head->next)
  {
    item = (struct item_t *) head->item;
    if (!nodeStack)
     {
       pllStackPush (&nodeStack, tr->nodep[i]);
       pllStackPush (&nodeStack, tr->nodep[i]->next);
       pllStackPush (&nodeStack, tr->nodep[i]->next->next);
       ++i;
     }
    else
     {
       v = (nodeptr) pllStackPop (&nodeStack);
       if (item->rank)  /* internal node */
        {
          v->back           = tr->nodep[i];
          tr->nodep[i]->back = v; //t->nodep[v->number]
          pllStackPush (&nodeStack, tr->nodep[i]->next);
          pllStackPush (&nodeStack, tr->nodep[i]->next->next);
          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
             v->z[k] = tr->nodep[i]->z[k] = z;

          ++ i;
        }
       else             /* leaf */
        {
          v->back           = tr->nodep[j];
          tr->nodep[j]->back = v; //t->nodep[v->number];

          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
            v->z[k] = tr->nodep[j]->z[k] = z;
            
          //t->nameList[j] = strdup (item->name);
          tr->nameList[j] = (char *) rax_malloc ((strlen (item->name) + 1) * sizeof (char));
          strcpy (tr->nameList[j], item->name);
          
          pllHashAdd (tr->nameHash, tr->nameList[j], (void *) (tr->nodep[j]));
          ++ j;
        }
     }
  }
  
  tr->start = tr->nodep[1];
  
  pllStackClear (&nodeStack);

  if (useDefaultz == PLL_TRUE) 
    resetBranches (tr);
}

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

  for (i = 1; i <= tips; ++ i)
   {
     tr->nameList[i] = (char *) rax_malloc ((strlen (nameList[i]) + 1) * sizeof (char));
     strcpy (tr->nameList[i], nameList[i]);
     pllHashAdd (tr->nameHash, tr->nameList[i], (void *) (tr->nodep[i]));
   }
  

  makeRandomTree (tr);
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
pllTreeInitTopologyForAlignment (pllInstance * tr, pllAlignmentData * alignmentData)
{
  int
    tips = alignmentData->sequenceCount,
    i;

  char 
    **nameList = alignmentData->sequenceLabels;
  
  pllTreeInitDefaults (tr, tips);

  for (i = 1; i <= tips; ++ i)
   {
     tr->nameList[i] = (char *) rax_malloc ((strlen (nameList[i]) + 1) * sizeof (char));
     strcpy (tr->nameList[i], nameList[i]);
     pllHashAdd (tr->nameHash, tr->nameList[i], (void *) (tr->nodep[i]));
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
*/
void pllComputeRandomizedStepwiseAdditionParsimonyTree(pllInstance * tr, partitionList * partitions)
{
  allocateParsimonyDataStructures(tr, partitions);
  makeParsimonyTreeFast(tr, partitions);
  freeParsimonyDataStructures(tr, partitions);
}

/** @brief Encode the alignment data to the PLL numerical representation
    
    Transforms the alignment to the PLL internal representation by substituting each base 
    with a specific digit.

    @param alignmentData
      Multiple sequence alignment

    @param partitions
      List of partitions
*/
void
pllBaseSubstitute (pllAlignmentData * alignmentData, partitionList * partitions)
{
  char meaningDNA[256];
  char  meaningAA[256];
  char * d;
  int i, j, k;

  for (i = 0; i < 256; ++ i)
   {
     meaningDNA[i] = -1;
     meaningAA[i]  = -1;
   }

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['Y'] = 10;
  meaningDNA['a'] =  1;
  meaningDNA['b'] = 14;
  meaningDNA['c'] =  2;
  meaningDNA['d'] = 13;
  meaningDNA['g'] =  4;
  meaningDNA['h'] = 11;
  meaningDNA['k'] = 12;
  meaningDNA['m'] =  3;
  meaningDNA['r'] =  5;
  meaningDNA['s'] =  6;
  meaningDNA['t'] =  8;
  meaningDNA['u'] =  8;
  meaningDNA['v'] =  7;
  meaningDNA['w'] =  9;
  meaningDNA['y'] = 10;

  meaningDNA['N'] =
  meaningDNA['n'] =
  meaningDNA['O'] =
  meaningDNA['o'] =
  meaningDNA['X'] =
  meaningDNA['x'] =
  meaningDNA['-'] =
  meaningDNA['?'] = 15;
 
  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['a'] =  0;  /* alanine */
  meaningAA['r'] =  1;  /* arginine */
  meaningAA['n'] =  2;  /*  asparagine*/
  meaningAA['d'] =  3;  /* aspartic */
  meaningAA['c'] =  4;  /* cysteine */
  meaningAA['q'] =  5;  /* glutamine */
  meaningAA['e'] =  6;  /* glutamic */
  meaningAA['g'] =  7;  /* glycine */
  meaningAA['h'] =  8;  /* histidine */
  meaningAA['i'] =  9;  /* isoleucine */
  meaningAA['l'] =  10; /* leucine */
  meaningAA['k'] =  11; /* lysine */
  meaningAA['m'] =  12; /* methionine */
  meaningAA['f'] =  13; /* phenylalanine */
  meaningAA['p'] =  14; /* proline */
  meaningAA['s'] =  15; /* serine */
  meaningAA['t'] =  16; /* threonine */
  meaningAA['w'] =  17; /* tryptophan */
  meaningAA['y'] =  18; /* tyrosine */
  meaningAA['v'] =  19; /* valine */
  meaningAA['b'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
  meaningAA['x'] = 
  meaningAA['?'] = 
  meaningAA['*'] = 
  meaningAA['-'] = 22;

  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     d = (partitions->partitionData[i]->dataType == PLL_DNA_DATA) ? meaningDNA : meaningAA;
     
     for (j = 1; j <= alignmentData->sequenceCount; ++ j)
      {
        for (k = partitions->partitionData[i]->lower; k < partitions->partitionData[i]->upper; ++ k)
         {
           alignmentData->sequenceData[j][k] = d[alignmentData->sequenceData[j][k]];
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
  
  pllHashDestroy (&(tr->nameHash), PLL_FALSE);
  if (tr->yVector)
   {
     if (tr->bDeep == PLL_DEEP_COPY)
      {
        if (tr->yVector[0]) rax_free (tr->yVector[0]);
      }
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

  ch = (char *) rax_malloc (strlen (linkageString) + 1);
  strcpy (ch, linkageString);


  for(j = 0, str1 = ch; ;j++, str1 = (char *)NULL) 
    {
      token = STRTOK_R(str1, ",", &saveptr);
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

/** @brief Set the alpha parameter of the Gamma model to a fixed value for a partition
    
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

  assert(model >= 0 && model < pr->numberOfPartitions);

  assert(alpha >= PLL_ALPHA_MIN && alpha <= PLL_ALPHA_MAX);

  //set the alpha paremeter 
  
  pr->partitionData[model]->alpha = alpha;

  //do the discretization of the gamma curve

  makeGammaCats(pr->partitionData[model]->alpha, pr->partitionData[model]->gammaRates, 4, tr->useMedian);

  //broadcast the changed parameters to all threads/MPI processes 

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier(tr, pr, PLL_THREAD_COPY_ALPHA);
#endif

  pr->partitionData[model]->optimizeAlphaParameter = PLL_FALSE;

  pr->dirty = PLL_FALSE;
}

/** @brief Set all base freuqncies to a fixed value for a partition
    
    Sets all base freuqencies of a partition to fixed values and disables 
    ML optimization of these parameters 

    @param f
      array containing the base frequencies

    @param 
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
    acc = 0.0;

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
  initReversibleGTR(tr, pr, model);


  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  
  pr->partitionData[model]->optimizeBaseFrequencies = PLL_FALSE;

  pr->dirty = PLL_TRUE;
}

/** @brief Set that the base freuqencies are optimized under ML
    
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
  initReversibleGTR(tr, pr, model);

  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  
  pr->partitionData[model]->optimizeBaseFrequencies = PLL_TRUE;

  pr->dirty = PLL_TRUE;

  return PLL_TRUE;
}






/** @brief Set all substitution rates to a fixed value for a specific partition
    
    Sets all substitution rates of a partition to fixed values and disables 
    ML optimization of these parameters. It will automatically re-scale the relative rates  
    such that the last rate is 1.0 

    @param f
      array containing the substitution rates

    @param 
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
  int 
    i,
    numberOfRates; 

  double
    scaler;

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
  initReversibleGTR(tr, pr, model);

  //broadcast the new Q matrix to all threads/processes 
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif
  
  pr->partitionData[model]->optimizeSubstitutionRates = PLL_FALSE;

  pr->dirty = PLL_TRUE;
}




/* initializwe a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
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



static linkageList* initLinkageListString(char *linkageString, partitionList * pr)
{
  int 
    *list = (int*)rax_malloc(sizeof(int) * pr->numberOfPartitions),
    j;

  linkageList 
    *l;

  char
    *str1,
    *saveptr,
//    *ch = strdup(linkageString),
    *ch,
    *token;
  
  ch = (char *) rax_malloc (strlen (linkageString) + 1);
  strcpy (ch, linkageString);

  for(j = 0, str1 = ch; ;j++, str1 = (char *)NULL) 
    {
      token = STRTOK_R(str1, ",", &saveptr);
      if(token == (char *)NULL)
	break;
      assert(j < pr->numberOfPartitions);
      list[j] = atoi(token);
      //printf("%d: %s\n", j, token);
    }
  
  rax_free(ch);

  l = initLinkageList(list, pr);
  
  rax_free(list);

  return l;
}

/** @brief Link alpha parameters across partitions
    
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

/** @brief Link base frequency parameters across partitions
    
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

/** @brief Link Substitution matrices across partitions
    
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




/** @brief Initialize partitions according to model parameters
    
    Initializes partitions according to model parameters.

    @param tr
      The PLL instance

    @param partitions
      List of partitions

    @param alignmentData
      The parsed alignment
    
    @return
      Returns \b PLL_TRUE in case of success, otherwise \b PLL_FALSE
*/
int pllInitModel (pllInstance * tr, partitionList * partitions, pllAlignmentData * alignmentData) 
{
  double ** ef;
  int
    i,
    *unlinked = (int *)rax_malloc(sizeof(int) * partitions->numberOfPartitions);

  ef = pllBaseFrequenciesGTR (partitions, alignmentData);

  if(!ef)
    return PLL_FALSE;

  
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#if (defined(__AVX) || defined(__SSE3))
  _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif
#endif 

  masterTime = gettime();         
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

  pllEvaluateGeneric (tr, partitions, tr->start, PLL_TRUE, PLL_FALSE);

  return PLL_TRUE;
}

/** @brief Optimize all free model parameters of the likelihood model
    
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

  fp = fopen (filename, "r");
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

