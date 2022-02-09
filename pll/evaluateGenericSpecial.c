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
 * @file evaluateGenericSpecial.c
 *   
 * @brief Functions for computing the log likelihood at a given branch of the tree (i.e. a virtual root that is placed at this branch)
 */
#include "mem_alloc.h"
#include "systypes.h"
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif

/* the set of functions in here computes the log likelihood at a given branch (the virtual root of a tree) */

/* includes for using SSE3 intrinsics */

#ifdef __SSE3
#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif
#endif


/** @defgroup evaluateLikelihoodGroup Likelihood evaluation
    
    This set of functions deals with the evaluation of likelihood for the current topology
*/







/* below are the function headers for unreadeble highly optimized versions of the above functions 
   for DNA and protein data that also use SSE3 intrinsics and implement some memory saving tricks.
   The actual functions can be found at the end of this source file. 
   All other likelihood function implementation files:

   newviewGenericSpacial.c
   makenewzSpecial.c
   evaluatePartialGenericSpecial.c

   are also structured like this 

   To decide which set of function implementations to use you will have to undefine or define _OPTIMIZED_FUNCTIONS 
   in the Makefile 
   */
#if (defined(__SSE3) || defined(__AVX))

static double evaluateGTRGAMMAPROT_LG4(int *ex1, int *ex2, int *wptr,
                                       double *x1, double *x2,  
                                       double *tipVector[4], 
                                       unsigned char *tipX1, int n, double *diagptable, const pllBoolean fastScaling,
                                       double * lg4_weights);

/* GAMMA for proteins with memory saving */

static double evaluateGTRGAMMAPROT_GAPPED_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                                double *x1, double *x2,  
                                                double *tipVector, 
                                                unsigned char *tipX1, int n, double *diagptable, 
                                                double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);


/* GAMMA for proteins */

static double evaluateGTRGAMMAPROT (const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                    double *x1, double *x2,  
                                    double *tipVector, 
                                    unsigned char *tipX1, int n, double *diagptable);

/* CAT for proteins */

static double evaluateGTRCATPROT (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                  double *x1, double *x2, double *tipVector,
                                  unsigned char *tipX1, int n, double *diagptable_start);


/* CAT for proteins with memory saving */

static double evaluateGTRCATPROT_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                       double *x1, double *x2, double *tipVector,
                                       unsigned char *tipX1, int n, double *diagptable_start, 
                                       double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

/* analogous DNA fuctions */

static double evaluateGTRCAT_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                   double *x1_start, double *x2_start, double *tipVector,                     
                                   unsigned char *tipX1, int n, double *diagptable_start,
                                   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA_GAPPED_SAVE(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                           double *x1_start, double *x2_start, 
                                           double *tipVector, 
                                           unsigned char *tipX1, const int n, double *diagptable,
                                           double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                               double *x1_start, double *x2_start, 
                               double *tipVector, 
                               unsigned char *tipX1, const int n, double *diagptable);


static double evaluateGTRCAT (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                              double *x1_start, double *x2_start, double *tipVector,                  
                              unsigned char *tipX1, int n, double *diagptable_start);


#endif

#if (defined(__AVX) || defined(__SSE3))
static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
                                      double *x1_start, double *x2_start, 
                                      double *tipVector, 
                                      unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling);

static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
                                     double *x1_start, double *x2_start, double *tipVector,                   
                                     unsigned char *tipX1, int n, double *diagptable_start, const pllBoolean fastScaling);
#endif


/* 
   global variables of pthreads version, reductionBuffer is the global array 
   that is used for implementing deterministic reduction operations, that is,
   the total log likelihood over the partial log lieklihoods for the sites that each thread has computed 

   NumberOfThreads is just the number of threads.

   Note the volatile modifier here, that guarantees that the compiler will not do weird optimizations 
   rearraengements of the code accessing those variables, because it does not know that several concurrent threads 
   will access those variables simulatenously 

   UPDATE: reductionBuffer is now merged with globalResult
   */


/* a pre-computed 32-bit integer mask */

extern const unsigned int mask32[32];

/* the function below computes the P matrix from the decomposition of the Q matrix and the respective rate categories for a single partition */

/** @brief Compute the diagonal of P matrix for a specific edge

    This function computes the diagonal of P matrix for a branch of length \a z
    from the decomposition of the Q matrix specified in \a EIGN and the respective
    rate categories \a rptr for a single partition. The diagonal is then stored in
    \a diagptable. 

    @param z                  Length of edge
    @param states             Number of states
    @param numberOfCategories Number of categories in the rate heterogeneity rate arrays
    @param rptr               Rate heterogeneity rate arrays
    @param EIGN               Eigenvalues
    @param diagptable         Where to store the resulting P matrix
*/
static void calcDiagptable(const double z, const int states, const int numberOfCategories, const double *rptr, const double *EIGN, double *diagptable)
{
  int 
    i, 
    l;

  double 
    lz,
    *lza = (double *)rax_malloc(sizeof(double) * states);

  /* transform the root branch length to the log and check if it is not too small */

  if (z < PLL_ZMIN) 
    lz = log(PLL_ZMIN);
  else
    lz = log(z);

  /* do some pre-computations to avoid redundant computations further below */

  for(i = 1; i < states; i++)      
    lza[i] = EIGN[i] * lz; 

  /* loop over the number of per-site or discrete gamma rate categories */

  for(i = 0; i < numberOfCategories; i++)
  {                    
    /* 
       diagptable is a pre-allocated array of doubles that stores the P-Matrix 
       the first entry is always 1.0 
       */
    diagptable[i * states] = 1.0;

    /* compute the P matrix for all remaining states of the model */

    for(l = 1; l < states; l++)
      diagptable[i * states + l] = exp(rptr[i] * lza[l]);
  }

  rax_free(lza);
}

/** @brief Compute the diagonal of P matrix for a specific edge for the LG4 model

    This function computes the diagonal of P matrix for a branch of length \a z
    from the decomposition of the 4 LG4 Q matrices specified in \a EIGN and the respective
    rate categories \a rptr for a single partition. The diagonal is then stored in
    \a diagptable. 

    @param z
      Length of edge

    @param states
      Number of states

    @param numberOfCategories
      Number of categories in the rate heterogeneity rate arrays

    @param rptr
      Rate heterogeneity rate arrays

    @param EIGN
      Eigenvalues of the 4 Q matrices

    @param diagptable
      Where to store the resulting P matrix

    @param numStates
      Number of states
*/
static void calcDiagptableFlex_LG4(double z, int numberOfCategories, double *rptr, double *EIGN[4], double *diagptable, const int numStates)
{
  int 
    i, 
    l;
  
  double 
    lz;
  
  assert(numStates <= 64);
  
  if (z < PLL_ZMIN) 
    lz = log(PLL_ZMIN);
  else
    lz = log(z);

  for(i = 0; i <  numberOfCategories; i++)
    {                  
      diagptable[i * numStates + 0] = 1.0;

      for(l = 1; l < numStates; l++)
        diagptable[i * numStates + l] = exp(rptr[i] * EIGN[i][l] * lz);                   
    }        
}

static void ascertainmentBiasSequence(unsigned char tip[32], int numStates)
{ 
  assert(numStates <= 32 && numStates > 1);

  switch(numStates)
    {
    case 2:     
      tip[0] = 1;
      tip[1] = 2;
      break;
    case 4:
      tip[0] = 1;
      tip[1] = 2;
      tip[2] = 4;
      tip[3] = 8;
      break;
    default:
      {
	int 
	  i;
	for(i = 0; i < numStates; i++)
	  {
	    tip[i] = i;
	    //printf("%c ", inverseMeaningPROT[i]);
	  }
	//printf("\n");
      }
      break;
    }
}

static double evaluateCatAsc(int *ex1, int *ex2,
			     double *x1, double *x2,  
			     double *tipVector, 
			     unsigned char *tipX1, int n, double *diagptable, const int numStates)
{
  double
    exponent,
    sum = 0.0, 
    unobserved,
    term,
    *left, 
    *right;
  
  int     
    i,    
    l;   
         
  unsigned char 
    tip[32];

  ascertainmentBiasSequence(tip, numStates);
   
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[numStates * tip[i]]);	  	  
	  right = &(x2[i * numStates]);

	  term = 0.0;
	         	      
	  for(l = 0; l < numStates; l++)
	    term += left[l] * right[l] * diagptable[l];	      	 	 	  	 

	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page:
	       If result underflows, and is not representable,
	       a range error occurs and 0.0 is returned.
	 */

	  exponent = pow(PLL_MINLIKELIHOOD, (double)ex2[i]);

	  unobserved = fabs(term) * exponent;

#ifdef _DEBUG_ASC
	  if(ex2[i] > 0)
	    {
	      printf("s %d\n", ex2[i]);
	      assert(0);
	    }
#endif	  
	    
	  sum += unobserved;
	}              
    }              
  else
    {           
      for (i = 0; i < n; i++) 
	{	  	 
	  term = 0.0;
	  	 
	  left  = &(x1[i * numStates]);
	  right = &(x2[i * numStates]);	    
	      
	  for(l = 0; l < numStates; l++)
	    term += left[l] * right[l] * diagptable[l];		  
	  
	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page:
	       If result underflows, and is not representable,
	       a range error occurs and 0.0 is returned.
	  */

	  exponent = pow(PLL_MINLIKELIHOOD, (double)(ex1[i] + ex2[i]));

	  unobserved = fabs(term) * exponent;
	  
#ifdef _DEBUG_ASC
	  if(ex2[i] > 0 || ex1[i] > 0)
	    {
	      printf("s %d %d\n", ex1[i], ex2[i]);
	      assert(0);
	    }
#endif

	  sum += unobserved;
	}             
    }        

  return  sum;
}


static double evaluateGammaAsc(int* ex1, int* ex2,
    double* x1, double* x2,
    double* tipVector,
    unsigned char* tipX1, int n, double* diagptable, const int numStates)
{
    double
        exponent,
        sum = 0.0,
        unobserved,
        term;

    int
        i,
        j,
        l;

    const int
        gammaStates = numStates * 4;

    unsigned char
        tip[32];

    ascertainmentBiasSequence(tip, numStates);

    if (tipX1)
    {
        for (i = 0; i < n; i++) {
            double* left = tipVector + numStates * tip[i];
            for (j = 0, term = 0.0; j < 4; j++) {
                double* right = x2 + gammaStates * i + numStates * j;
                for (l = 0; l < numStates; l++) {
                    term += left[l] * right[l] * diagptable[j * numStates + l];
                }
            }

            /* assumes that pow behaves as expected/specified for underflows
               from the man page:
                 If result underflows, and is not representable,
                 a range error occurs and 0.0 is returned.
            */

            exponent = pow(PLL_MINLIKELIHOOD, (double)ex2[i]);

            unobserved = fabs(term) * exponent;

#ifdef _DEBUG_ASC
            if (ex2[i] > 0)
            {
                printf("s %d\n", ex2[i]);
                assert(0);
            }
#endif	  

            sum += unobserved;
        }
    }
    else
    {
        for (i = 0; i < n; i++) {
            for (j = 0, term = 0.0; j < 4; j++) {
                double* left  = x1 + gammaStates * i + numStates * j;
                double* right = x2 + gammaStates * i + numStates * j;
                for (l = 0; l < numStates; l++) {
                    term += left[l] * right[l] * diagptable[j * numStates + l];
                }
            }

            /* assumes that pow behaves as expected/specified for underflows
               from the man page:
                 If result underflows, and is not representable,
                 a range error occurs and 0.0 is returned.
            */

            exponent = pow(PLL_MINLIKELIHOOD, (double)(ex1[i] + ex2[i]));

            unobserved = fabs(term) * exponent;

#ifdef _DEBUG_ASC
            if (ex2[i] > 0 || ex1[i] > 0)
            {
                printf("s %d %d\n", ex1[i], ex2[i]);
                assert(0);
            }
#endif

            sum += unobserved;
        }
    }
    return  sum;
}


/** @ingroup evaluateLikelihoodGroup
    @brief A generic (and slow) implementation of log likelihood evaluation of a tree using the GAMMA model of rate heterogeneity
    
    Computes the log likelihood of the topology for a specific partition, assuming
    that the GAMMA model of rate heterogeneity is used. The likelihood is computed at
    a virtual root placed at an edge whose two end-points (nodes) have the conditional
    likelihood vectors \a x1 and \a x2. 
    Furthermore, if \a getPerSiteLikelihoods is set to \b PLL_TRUE, then the log
    likelihood for each site is also computed and stored at the corresponding position
    in the array \a perSiteLikelihoods.

    @param fastScaling
      If set to \b PLL_FALSE, then the likelihood of each site is also multiplied by \a log(PLL_MINLIKELIHOOD) times the number
      of times it has been scaled down

    @param ex1
      An array that holds how many times a site has been scaled and points at the entries for node \a p. This
      parameter is used if \a fastScaling is set to \b PLL_FALSE.

    @param ex2
      An array that holds how many times a site has been scaled and points at the entries for node \a q. This
      parameter is used if \a fastScaling is set to \b PLL_TRUE.

    @param wptr
      Array holding the weight for each site in the compressed partition alignment

    @param x1_start
      Conditional likelihood vectors for one of the two end-points of the specific edge for which we are evaluating the likelihood

    @param x2_start
      Conditional likelihood vectors for the other end-point of the specific edge for which we are evaluating the likelihood

    @param tipVector
      Precomputed table where the number of rows is equal to the number of possible basepair characters for the current data 
      type, i.e.16 for DNA and 23 for AA, and each rows contains \a states elements each of which contains transition
      probabilities computed from the eigenvectors of the decomposed Q matrix.

    @param tipX1
      If one of the two end-points (nodes) of the specific edge (for which we are evaluating the likelihood) is a tip, then
      this holds a pointer to the sequence data (basepairs) already converted in the internal integer representation, and \a x2
      holds the conditional likelihood vectors for the internal node.

    @param n
      Number of sites for which we are doing the evaluation. For the single-thread version this is the 
      number of sites in the current partition, for multi-threads this is the number of sites assigned
      to the running thread from the current partition.

    @param diagptable
      Start of the array that contains the P-Matrix diagonal of the specific edge for which we are
      evaluating the likehood, and for each category of the GAMMA model

    @param states
      Number of states (4 for DNA, 20 for AA)

    @param perSiteLikelihoods
      Array to store per-site log likelihoods if \a getPerSiteLikelihoods is set to \b PLL_TRUE

    @param getPerSiteLikelihoods
      If set to \b PLL_TRUE then per-site log likelihoods are also computed and stored in \a perSiteLikelihoods

    @return
      The evaluated log likelihood of the tree topology
*/
static double evaluateGAMMA_FLEX(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                 double *x1_start, double *x2_start, 
                                 double *tipVector, 
                                 unsigned char *tipX1, const int n, double *diagptable, const int states, double *perSiteLikelihoods, pllBoolean getPerSiteLikelihoods)
{
  double   
    sum = 0.0, 
    term,
    *x1,
    *x2;

  int     
    i, 
    j,
    k;

  /* span is the offset within the likelihood array at an inner node that gets us from the values 
     of site i to the values of site i + 1 */

  const int 
    span = states * 4;


  /* we distingusih between two cases here: one node of the two nodes defining the branch at which we put the virtual root is 
     a tip. Both nodes can not be tips because we do not allow for two-taxon trees ;-) 
     Nota that, if a node is a tip, this will always be tipX1. This is done for code simplicity and the flipping of the nodes
     is done before when we compute the traversal descriptor.     
     */

  /* the left node is a tip */
  if(tipX1)
  {             
    /* loop over the sites of this partition */
    for (i = 0; i < n; i++)
    {
      /* access pre-computed tip vector values via a lookup table */
      x1 = &(tipVector[states * tipX1[i]]);      
      /* access the other(inner) node at the other end of the branch */
      x2 = &(x2_start[span * i]);        

      /* loop over GAMMA rate categories, hard-coded as 4 in RAxML */
      for(j = 0, term = 0.0; j < 4; j++)
        /* loop over states and multiply them with the P matrix */
        for(k = 0; k < states; k++)
          term += x1[k] * x2[j * states + k] * diagptable[j * states + k];                                                        

      /* take the log of the likelihood and multiply the per-gamma rate likelihood by 1/4.
         Under the GAMMA model the 4 discrete GAMMA rates all have the same probability 
         of 0.25 */

      if(!fastScaling)
        term = log(0.25 * fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));

      /* if required get the per-site log likelihoods.
         note that these are the plain per site log-likes, not 
         multiplied with the pattern weight value */
      
      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;
    }     
  }
  else
  {        
    for (i = 0; i < n; i++) 
    {
      /* same as before, only that now we access two inner likelihood vectors x1 and x2 */

      x1 = &(x1_start[span * i]);
      x2 = &(x2_start[span * i]);                 

      for(j = 0, term = 0.0; j < 4; j++)
        for(k = 0; k < states; k++)
          term += x1[j * states + k] * x2[j * states + k] * diagptable[j * states + k];

      if(!fastScaling)
        term = log(0.25 * fabs(term)) + ((ex1[i] + ex2[i])*log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));
      
      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;
    }                           
  }

  return sum;
} 

#if (defined(__SSE3) || defined(__AVX))
/** @ingroup evaluateLikelihoodGroup
    @brief Memory saving version of the generic (and slow) implementation of log likelihood evaluation of a tree using the GAMMA model of rate heterogeneity

    Computes the log likelihood of the topology for a specific partition, assuming
    that the GAMMA model of rate heterogeneity is used and memory saving technique
    is enabled. The likelihood is computed at a virtual root placed at an edge whose
    two end-points (nodes) have the conditional likelihood vectors \a x1 and \a x2. 
    Furthermore, if \a getPerSiteLikelihoods is set to \b PLL_TRUE, then the log
    likelihood for each site is also computed and stored at the corresponding position
    in the array \a perSiteLikelihoods.

    @param fastScaling
      If set to \b PLL_FALSE, then the likelihood of each site is also multiplied by \a log(PLL_MINLIKELIHOOD) times the number
      of times it has been scaled down

    @param ex1
      An array that holds how many times a site has been scaled and points at the entries for node \a p. This
      parameter is used if \a fastScaling is set to \b PLL_FALSE.

    @param ex2
      An array that holds how many times a site has been scaled and points at the entries for node \a q. This
      parameter is used if \a fastScaling is set to \b PLL_TRUE.

    @param wptr
      Array holding the weight for each site in the compressed partition alignment

    @param x1_start
      Conditional likelihood vectors for one of the two end-points of the specific edge for which we are evaluating the likelihood

    @param x2_start
      Conditional likelihood vectors for the other end-point of the specific edge for which we are evaluating the likelihood

    @param tipVector
      Precomputed table where the number of rows is equal to the number of possible basepair characters for the current data 
      type, i.e.16 for DNA and 23 for AA, and each rows contains \a states elements each of which contains transition
      probabilities computed from the eigenvectors of the decomposed Q matrix.

    @param tipX1
      If one of the two end-points (nodes) of the specific edge (for which we are evaluating the likelihood) is a tip, then
      this holds a pointer to the sequence data (basepairs) already converted in the internal integer representation, and \a x2
      holds the conditional likelihood vectors for the internal node.

    @param n
      Number of sites for which we are doing the evaluation. For the single-thread version this is the 
      number of sites in the current partition, for multi-threads this is the number of sites assigned
      to the running thread from the current partition.

    @param diagptable
      Start of the array that contains the P-Matrix diagonal of the specific edge for which we are
      evaluating the likehood, and for each category of the GAMMA model

    @param states
      Number of states (4 for DNA, 20 for AA)

    @param perSiteLikelihoods
      Array to store per-site log likelihoods if \a getPerSiteLikelihoods is set to \b PLL_TRUE

    @param getPerSiteLikelihoods
      If set to \b PLL_TRUE then per-site log likelihoods are also computed and stored in \a perSiteLikelihoods

    @param x1_gapColumn

    @param x2_gapColumn

    @param x1_gap
      Gap bitvector for the left child node

    @param x2_gap
      Gap bitvector for the right child node

    @return
      The evaluated log likelihood of the tree topology

    @todo
      Document x1_gapColumn, x2_gapColumn, x1_gap, x2_gap and add a brief description of how this technique works
*/
static double evaluateGAMMA_FLEX_SAVE(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                      double *x1_start, double *x2_start, 
                                      double *tipVector, 
                                      unsigned char *tipX1, const int n, double *diagptable, const int states, double *perSiteLikelihoods, pllBoolean getPerSiteLikelihoods,
                                      double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
    term,
    *x1,
    *x2,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;
    
  int     
    i, 
    j,
    k;

  /* span is the offset within the likelihood array at an inner node that gets us from the values 
     of site i to the values of site i + 1 */

  const int 
    span = states * 4;


  /* we distingusih between two cases here: one node of the two nodes defining the branch at which we put the virtual root is 
     a tip. Both nodes can not be tips because we do not allow for two-taxon trees ;-) 
     Nota that, if a node is a tip, this will always be tipX1. This is done for code simplicity and the flipping of the nodes
     is done before when we compute the traversal descriptor.     
     */

  /* the left node is a tip */
  if(tipX1)
  {             
    /* loop over the sites of this partition */
    for (i = 0; i < n; i++)
    {
      /* access pre-computed tip vector values via a lookup table */
      x1 = &(tipVector[states * tipX1[i]]);      
      /* access the other(inner) node at the other end of the branch */

      if(x2_gap[i / 32] & mask32[i % 32])
        x2 = x2_gapColumn;
      else
        {
          x2 = x2_ptr;
          x2_ptr += span;
        }

      /* loop over GAMMA rate categories, hard-coded as 4 in RAxML */
      for(j = 0, term = 0.0; j < 4; j++)
        /* loop over states and multiply them with the P matrix */
        for(k = 0; k < states; k++)
          term += x1[k] * x2[j * states + k] * diagptable[j * states + k];                                                        

      /* take the log of the likelihood and multiply the per-gamma rate likelihood by 1/4.
         Under the GAMMA model the 4 discrete GAMMA rates all have the same probability 
         of 0.25 */

      if(!fastScaling)
        term = log(0.25 * fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));

      /* if required get the per-site log likelihoods.
         note that these are the plain per site log-likes, not 
         multiplied with the pattern weight value */
      
      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;
    }     
  }
  else
  {        
    for (i = 0; i < n; i++) 
    {
      /* same as before, only that now we access two inner likelihood vectors x1 and x2 */
      
      if(x1_gap[i / 32] & mask32[i % 32])
        x1 = x1_gapColumn;
      else
        {
          x1 = x1_ptr;
          x1_ptr += span;
        }    

      if(x2_gap[i / 32] & mask32[i % 32])
        x2 = x2_gapColumn;
      else
        {
          x2 = x2_ptr;
          x2_ptr += span;
        }                 

      for(j = 0, term = 0.0; j < 4; j++)
        for(k = 0; k < states; k++)
          term += x1[j * states + k] * x2[j * states + k] * diagptable[j * states + k];

      if(!fastScaling)
        term = log(0.25 * fabs(term)) + ((ex1[i] + ex2[i])*log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));
      
      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;
    }                           
  }

  return sum;
} 
#endif

/** @ingroup evaluateLikelihoodGroup
    @brief A generic (and slow) implementation of log likelihood evaluation of a tree using the CAT model of rate heterogeneity
    
    Computes the log likelihood of the topology for a specific partition, assuming
    that the CAT model of rate heterogeneity is used. The likelihood is computed at
    a virtual root placed at an edge whose two end-points (nodes) have the conditional
    likelihood vectors \a x1 and \a x2. 
    Furthermore, if \a getPerSiteLikelihoods is set to \b PLL_TRUE, then the log
    likelihood for each site is also computed and stored at the corresponding position
    in the array \a perSiteLikelihoods.

    @param fastScaling
      If set to \b PLL_FALSE, then the likelihood of each site is also multiplied by \a log(PLL_MINLIKELIHOOD) times the number
      of times it has been scaled down

    @param ex1
      An array that holds how many times a site has been scaled and points at the entries for node \a p. This
      parameter is used if \a fastScaling is set to \b PLL_FALSE.

    @param ex2
      An array that holds how many times a site has been scaled and points at the entries for node \a q. This
      parameter is used if \a fastScaling is set to \b PLL_TRUE.

    @param cptr
      Array holding the rate for each site in the compressed partition alignment

    @param wptr
      Array holding the weight for each site in the compressed partition alignment

    @param x1
      Conditional likelihood vectors for one of the two end-points of the specific edge for which we are evaluating the likelihood

    @param x2
      Conditional likelihood vectors for the other end-point of the specific edge for which we are evaluating the likelihood

    @param tipVector
      Precomputed table where the number of rows is equal to the number of possible basepair characters for the current data type, 
      i.e.16 for DNA and 23 for AA, and each rows contains \a states elements each of which contains transition probabilities 
      computed from the eigenvectors of the decomposed Q matrix.

    @param tipX1
      If one of the two end-points (nodes) of the specific edge (for which we are evaluating the likelihood) is a tip, then
      this holds a pointer to the sequence data (basepairs) already converted in the internal integer representation, and \a x2
      holds the conditional likelihood vectors for the internal node.

    @param n
      Number of sites for which we are doing the evaluation. For the single-thread version this is the number of sites in the
      current partition, for multi-threads this is the number of sites assigned to the running thread from the current partition.

    @param diagptable_start
      Start of the array that contains the P-Matrix diagonal of the specific edge for which we are evaluating the likehood,
      and for each category of the CAT model

    @param states
      Number of states (4 for DNA, 20 for AA)

    @param perSiteLikelihoods
      Array to store per-site log likelihoods if \a getPerSiteLikelihoods is set to \b PLL_TRUE

    @param getPerSiteLikelihoods
      If set to \b PLL_TRUE then per-site log likelihoods are also computed and stored in \a perSiteLikelihoods

    @return
      The evaluated log likelihood of the tree topology
*/
static double evaluateCAT_FLEX (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                double *x1, double *x2, double *tipVector,
                                unsigned char *tipX1, int n, double *diagptable_start, const int states, double *perSiteLikelihoods, pllBoolean getPerSiteLikelihoods)
{
  double   
    sum = 0.0, 
    term,
    *diagptable,  
    *left, 
    *right;

  int     
    i, 
    l;                           

  /* chosing between tip vectors and non tip vectors is identical in all flavors of this function ,regardless 
     of whether we are using CAT, GAMMA, DNA or protein data etc */

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {
      /* same as in the GAMMA implementation */
      left = &(tipVector[states * tipX1[i]]);
      right = &(x2[states * i]);

      /* important difference here, we do not have, as for GAMMA 
         4 P matrices assigned to each site, but just one. However those 
         P-Matrices can be different for the sites.
         Hence we index into the precalculated P-matrices for individual sites 
         via the category pointer cptr[i]
         */
      diagptable = &diagptable_start[states * cptr[i]];                  

      /* similar to gamma, with the only difference that we do not integrate (sum)
         over the discrete gamma rates, but simply compute the likelihood of the 
         site and the given P-matrix */

      for(l = 0, term = 0.0; l < states; l++)
        term += left[l] * right[l] * diagptable[l];                        

      /* take the log */
       if(!fastScaling)
         term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
       else
         term = log(fabs(term));

       /* if required get the per-site log likelihoods.
          note that these are the plain per site log-likes, not 
          multiplied with the pattern weight value */

       if(getPerSiteLikelihoods)
         perSiteLikelihoods[i] = term;

      /* 
         multiply the log with the pattern weight of this site. 
         The site pattern for which we just computed the likelihood may 
         represent several alignment columns sites that have been compressed 
         into one site pattern if they are exactly identical AND evolve under the same model,
         i.e., form part of the same partition.
         */                  

      sum += wptr[i] * term;
    }      
  }    
  else
  {    
    for (i = 0; i < n; i++) 
    {   
      /* as before we now access the likelihood arrayes of two inner nodes */
      left  = &x1[states * i];
      right = &x2[states * i];

      diagptable = &diagptable_start[states * cptr[i]];         

      for(l = 0, term = 0.0; l < states; l++)
        term += left[l] * right[l] * diagptable[l];
      
      if(!fastScaling)
        term = log(fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));  

      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 

#if (defined(__SSE3) || defined(__AVX))
/** @ingroup evaluateLikelihoodGroup
    @brief A generic (and slow) implementation of log likelihood evaluation of a tree using the CAT model of rate heterogeneity with memory saving
    
    This is the same as ::evaluateCAT_FLEX but with the memory saving technique enabled.
    Please check ::evaluateCAT_FLEX for more information and a description of the common
    input parameters
    
    @param x1_gapColumn

    @param x2_gapColumn

    @param x1_gap
      Gap bitvector for the left child node

    @param x2_gap
      Gap bitvector for the right child node
    
    @todo
      Comment on x1_gapColumn and x2_gapColumn
*/
static double evaluateCAT_FLEX_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                     double *x1, double *x2, double *tipVector,
                                     unsigned char *tipX1, int n, double *diagptable_start, const int states, double *perSiteLikelihoods, pllBoolean getPerSiteLikelihoods,
                                     double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
    term,
    *diagptable,  
    *left, 
    *right,
    *left_ptr = x1,
    *right_ptr = x2;

  int     
    i, 
    l;                           

  /* chosing between tip vectors and non tip vectors is identical in all flavors of this function ,regardless 
     of whether we are using CAT, GAMMA, DNA or protein data etc */

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {
      /* same as in the GAMMA implementation */
      left = &(tipVector[states * tipX1[i]]);
   
      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
        {
          right = right_ptr;
          right_ptr += states;
        }         
      /* important difference here, we do not have, as for GAMMA 
         4 P matrices assigned to each site, but just one. However those 
         P-Matrices can be different for the sites.
         Hence we index into the precalculated P-matrices for individual sites 
         via the category pointer cptr[i]
         */
      diagptable = &diagptable_start[states * cptr[i]];                  

      /* similar to gamma, with the only difference that we do not integrate (sum)
         over the discrete gamma rates, but simply compute the likelihood of the 
         site and the given P-matrix */

      for(l = 0, term = 0.0; l < states; l++)
        term += left[l] * right[l] * diagptable[l];                        

      /* take the log */
       if(!fastScaling)
         term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
       else
         term = log(fabs(term));

       /* if required get the per-site log likelihoods.
          note that these are the plain per site log-likes, not 
          multiplied with the pattern weight value */

       if(getPerSiteLikelihoods)
         perSiteLikelihoods[i] = term;

      /* 
         multiply the log with the pattern weight of this site. 
         The site pattern for which we just computed the likelihood may 
         represent several alignment columns sites that have been compressed 
         into one site pattern if they are exactly identical AND evolve under the same model,
         i.e., form part of the same partition.
         */                  

      sum += wptr[i] * term;
    }      
  }    
  else
  {    
    for (i = 0; i < n; i++) 
    {   
      /* as before we now access the likelihood arrayes of two inner nodes */     

      if(isGap(x1_gap, i))
        left = x1_gapColumn;
      else
        {
          left = left_ptr;
          left_ptr += states;
        }       

      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
        {
          right = right_ptr;
          right_ptr += states;
        }       

      diagptable = &diagptable_start[states * cptr[i]];         

      for(l = 0, term = 0.0; l < states; l++)
        term += left[l] * right[l] * diagptable[l];
      
      if(!fastScaling)
        term = log(fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));  

      if(getPerSiteLikelihoods)
        perSiteLikelihoods[i] = term;

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 
#endif


/* This is the core function for computing the log likelihood at a branch */
/** @ingroup evaluateLikelihoodGroup
    @brief Evaluate the log likelihood of a specific branch of the topology
    
    Evaluates the likelihood of the tree topology assuming a virtual root is
    placed at the edge whose end-points are node with number \a pNumber and \a
    qNumber in the first slot of the traversal descriptor. The function first
    computes the conditional likelihoods for all necessary nodes (the ones in
    the traversal descriptor list) by calling the function \a pllNewviewIterative
    and then evaluates the likelihood at the root. In addition, if \a
    getPerSiteLikelihoods is set to \b PLL_TRUE, the per-site likelihoods are
    stored in \a tr->lhs.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param getPerSiteLikelihoods
      If set to \b PLL_TRUE, compute the log likelihood for each site. 

    @note
      This is an internal function and should not be called by the user. It assumes
      that a valid traversal descriptor has already been computed. It also assumes
      that the edge we are referring to is an edge that leads to a tip, i.e. either
      p or q of the first entry of traversal descriptor are tips.
*/
void pllEvaluateIterative(pllInstance *tr, partitionList *pr, pllBoolean getPerSiteLikelihoods)
{
  /* the branch lengths and node indices of the virtual root branch are always the first one that 
     are stored in the very important traversal array data structure that describes a partial or full tree traversal */

  /* get the branch length at the root */
  double 
    *pz = tr->td[0].ti[0].qz;   

  /* get the node number of the node to the left and right of the branch that defines the virtual rooting */

  int    
    pNumber = tr->td[0].ti[0].pNumber, 
    qNumber = tr->td[0].ti[0].qNumber, 
    p_slot,
    q_slot,
    model;
  
  pllBoolean
    fastScaling = tr->fastScaling;

  /* the slots are the entries in xVector where the LH vector is available */
  if(tr->useRecom)
    {
      p_slot = tr->td[0].ti[0].slot_p;
      q_slot = tr->td[0].ti[0].slot_q;
    }
  else
    {
      p_slot = pNumber - tr->mxtips - 1;
      q_slot = qNumber - tr->mxtips - 1;
    }
  
  /* before we can compute the likelihood at the virtual root, we need to do a partial or full tree traversal to compute 
     the conditional likelihoods of the vectors as specified in the traversal descriptor. Maintaining this tarversal descriptor consistent 
     will unfortunately be the responsibility of users. This is tricky, if as planned for here, we use a rooted view (described somewhere in Felsenstein's book)
     for the conditional vectors with respect to the tree
     */

  /* iterate over all valid entries in the traversal descriptor */

  pllNewviewIterative(tr, pr, 1);

  /* after the above call we are sure that we have properly and consistently computed the 
     conditionals to the right and left of the virtual root and we can now invoke the 
     the log likelihood computation */

  /* we need to loop over all partitions. Note that we may have a mix of DNA, protein binary data etc partitions */

  for(model = 0; model < pr->numberOfPartitions; model++)
    {    
      /* whats' the number of sites of this partition (at the current thread) */
      int           
        width = pr->partitionData[model]->width;
      
      /* 
         Important part of the tarversal descriptor: 
         figure out if we need to recalculate the likelihood of this 
         partition: 
         
         The reasons why this is important in terms of performance are given in this paper 
         here which you should actually read:
         
         A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009, accepted for publication, Vienna, Austria, September 2009
         
         The width > 0 check is for checking if under the cyclic data distribution of per-partition sites to threads this thread does indeed have a site 
         of the current partition.
         
      */

      if(tr->td[0].executeModel[model] && width > 0)
        {       
          int 
#if (defined(__SSE3) || defined(__AVX))
            rateHet = (int)discreteRateCategories(tr->rateHetModel),
#endif
            categories,
            ascWidth = pr->partitionData[model]->states,
            
            /* get the number of states in the partition, e.g.: 4 = DNA, 20 = Protein */
            
            states = pr->partitionData[model]->states,
            *ex1 = NULL,
            *ex2 = NULL,
            *ex1_asc = NULL,
            *ex2_asc = NULL;
          
          double 
            *rateCategories = (double*)NULL,
            z, 
            partitionLikelihood = 0.0,
            *x1_start           = NULL,
            *x2_start           = NULL,
            *diagptable         = NULL,
            *x1_start_asc       = NULL,
            *x2_start_asc       = NULL;

#if (defined(__SSE3) || defined(__AVX))
          double
            *x1_gapColumn = (double*)NULL,
            *x2_gapColumn = (double*)NULL;
#endif
          
#if (defined(__SSE3) || defined(__AVX))
          unsigned int
            *x1_gap = (unsigned int*)NULL,
            *x2_gap = (unsigned int*)NULL;       
#endif
          
          unsigned char 
            *tip = (unsigned char*)NULL;          
          
          /* 
             figure out if we are using the CAT or GAMMA model of rate heterogeneity 
             and set pointers to the rate heterogeneity rate arrays and also set the 
             number of distinct rate categories appropriately.
             
             Under GAMMA this is constant and hard-coded as 4, weheras under CAT 
             the number of site-wise rate categories can vary in the course of computations 
             up to a user defined maximum value of site categories (default: 25)
          */

          if(tr->rateHetModel == PLL_CAT)
            {        
              rateCategories = pr->partitionData[model]->perSiteRates;
              categories = pr->partitionData[model]->numberOfCategories;
            }
          else  /* GAMMA */
            {        
              rateCategories = pr->partitionData[model]->gammaRates;
              categories = 4;
            }
          
          /* set this pointer to the memory area where space has been reserved a priori for storing the 
             P matrix at the root */
          
          diagptable = pr->partitionData[model]->left;
          
          /* figure out if we need to address tip vectors (a char array that indexes into a precomputed tip likelihood 
             value array) or if we need to address inner vectors */
          
          /* either node p or node q is a tip */
          
          if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
            {                       
              /* q is a tip */
              
              if(isTip(qNumber, tr->mxtips))
                {       
                  /* get the start address of the inner likelihood vector x2 for partition model,
                     note that inner nodes are enumerated/indexed starting at 0 to save allocating some 
                     space for additional pointers */

                  x2_start = pr->partitionData[model]->xVector[p_slot];
                  
                  /* get the corresponding tip vector */
                  
                  tip      = pr->partitionData[model]->yVector[qNumber];

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
                  if (tr->threadID == 0 && pr->partitionData[model]->ascBias)
#else
                  if (pr->partitionData[model]->ascBias)
#endif
                   {
                     x2_start_asc  = &pr->partitionData[model]->ascVector[(pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                     ex2_asc       = &pr->partitionData[model]->ascExpVector[(pNumber - tr->mxtips - 1) * ascWidth];
                   }

                  
                  /* memory saving stuff, let's deal with this later or ask Fernando ;-) */
                  
#if (defined(__SSE3) || defined(__AVX))
                  if(tr->saveMemory)
                    {
                      x2_gap         = &(pr->partitionData[model]->gapVector[pNumber * pr->partitionData[model]->gapVectorLength]);
                      x2_gapColumn   = &(pr->partitionData[model]->gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet]);
                    }
#endif
                  /* per site likelihood scaling */

                  if(!fastScaling)                  
                    ex2 = pr->partitionData[model]->expVector[p_slot];              
                }           
              else
                {       
                  /* p is a tip, same as above */
                  
                  x2_start = pr->partitionData[model]->xVector[q_slot];
                  tip = pr->partitionData[model]->yVector[pNumber];


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
                  if (tr->threadID == 0 && pr->partitionData[model]->ascBias)
#else
                  if (pr->partitionData[model]->ascBias)
#endif
                   {
                     x2_start_asc  = &pr->partitionData[model]->ascVector[(qNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                     ex2_asc       = &pr->partitionData[model]->ascExpVector[(qNumber - tr->mxtips - 1) * ascWidth];
                   }
                  
#if (defined(__SSE3) || defined(__AVX))
                  if(tr->saveMemory)
                    {
                      x2_gap         = &(pr->partitionData[model]->gapVector[qNumber * pr->partitionData[model]->gapVectorLength]);
                      x2_gapColumn   = &(pr->partitionData[model]->gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet]);
                    }
#endif

                  /* per site likelihood scaling */

                  if(!fastScaling)                  
                    ex2 = pr->partitionData[model]->expVector[q_slot];             
                }
            }
          else
            {  
              
              assert(p_slot != q_slot);
              /* neither p nor q are tips, hence we need to get the addresses of two inner vectors */
              
              x1_start = pr->partitionData[model]->xVector[p_slot];
              x2_start = pr->partitionData[model]->xVector[q_slot];
              
              /* memory saving option */
              
#if (defined(__SSE3) || defined(__AVX))
              if(tr->saveMemory)
                {
                  x1_gap = &(pr->partitionData[model]->gapVector[pNumber * pr->partitionData[model]->gapVectorLength]);
                  x2_gap = &(pr->partitionData[model]->gapVector[qNumber * pr->partitionData[model]->gapVectorLength]);
                  x1_gapColumn   = &pr->partitionData[model]->gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
                  x2_gapColumn   = &pr->partitionData[model]->gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
                }
#endif
                      
              /* per site likelihood scaling */

              if(!fastScaling)
                {
                  ex1      = pr->partitionData[model]->expVector[p_slot];
                  ex2      = pr->partitionData[model]->expVector[q_slot];     
                }
              
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
              if (tr->threadID == 0 && pr->partitionData[model]->ascBias)
#else
              if (pr->partitionData[model]->ascBias)
#endif
               {
                 x1_start_asc  = &pr->partitionData[model]->ascVector[(pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                 x2_start_asc  = &pr->partitionData[model]->ascVector[(qNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];

                 ex1_asc       = &pr->partitionData[model]->ascExpVector[(pNumber - tr->mxtips - 1) * ascWidth];
                 ex2_asc       = &pr->partitionData[model]->ascExpVector[(qNumber - tr->mxtips - 1) * ascWidth];
               }



            }
          
          
          /* if we are using a per-partition branch length estimate, the branch has an index, otherwise, for a joint branch length
             estimate over all partitions we just use the branch length value with index 0 */
          
          if(pr->perGeneBranchLengths)
            z = pz[model];
          else
            z = pz[0];
          
          /* calc P-Matrix at root for branch z connecting nodes p and q */
          
          if(pr->partitionData[model]->dataType == PLL_AA_DATA && (pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X))
            calcDiagptableFlex_LG4(z, 4, pr->partitionData[model]->gammaRates, pr->partitionData[model]->EIGN_LG4, diagptable, 20);
          else
            calcDiagptable(z, states, categories, rateCategories, pr->partitionData[model]->EIGN, diagptable);
          
#if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
          
          /* generic slow functions, memory saving option is not implemented for these */
          
          assert(!tr->saveMemory);
          
          /* decide wheter CAT or GAMMA is used and compute log like */
          if(tr->rateHetModel == PLL_CAT)
            partitionLikelihood = evaluateCAT_FLEX(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt, 
                                                x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, getPerSiteLikelihoods);
          else
            partitionLikelihood = evaluateGAMMA_FLEX(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, getPerSiteLikelihoods);
#else
   
          /* if we want to compute the per-site likelihoods, we use the generic evaluate function implementations 
             for this, because the slowdown is not that dramatic */

          if(getPerSiteLikelihoods)
            {         
#ifdef __MIC_NATIVE
                          // not supported on MIC!
                          assert(0 && "Per-site LH calculations is not implemented on Intel MIC");
#else
               if(tr->rateHetModel == PLL_CAT)
                {
                   if(tr->saveMemory)
                     partitionLikelihood = evaluateCAT_FLEX_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                 x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                 tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, PLL_TRUE,
                                                                 x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                   else
                     partitionLikelihood = evaluateCAT_FLEX(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                            x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                            tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, PLL_TRUE);
                }
              else
                {
                  if(tr->saveMemory)
                    partitionLikelihood = evaluateGAMMA_FLEX_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                                  x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                  tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, PLL_TRUE, 
                                                                  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);              
                  else
                    partitionLikelihood = evaluateGAMMA_FLEX(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                             x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                             tip, width, diagptable, states, pr->partitionData[model]->perSiteLikelihoods, PLL_TRUE);
                }
#endif
            }
          else
            {
              /* for the optimized functions we have a dedicated, optimized function implementation 
                 for each rate heterogeneity and data type combination, we switch over the number of states 
                 and the rate heterogeneity model */
              
              switch(states)
                {         
                case 2: /* binary */
                  assert (!tr->saveMemory);
                  if (tr->rateHetModel == PLL_CAT)
                   {
                     partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                  x1_start, x2_start, pr->partitionData[model]->tipVector, 
                                                                  tip, width, diagptable, fastScaling);
                   }
                  else
                   {
                     partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2, pr->partitionData[model]->wgt,
                                                                   x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                   tip, width, diagptable, fastScaling);                 
                   }
                  break;
                case 4: /* DNA */
                  {

#ifdef __MIC_NATIVE

                  /* CAT & memory saving are not supported on MIC */

                  assert(!tr->saveMemory);
                  assert(tr->rateHetModel == PLL_GAMMA);

                  partitionLikelihood =  evaluateGTRGAMMA_MIC(ex1, ex2, pr->partitionData[model]->wgt,
                                              x1_start, x2_start, pr->partitionData[model]->tipVector,
                                              tip, width, diagptable, fastScaling);
#else
                    if(tr->rateHetModel == PLL_CAT)
                      {                           
                        if(tr->saveMemory)
                          partitionLikelihood =  evaluateGTRCAT_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                     x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                     tip, width, diagptable, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                        else
                          partitionLikelihood =  evaluateGTRCAT(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                tip, width, diagptable);
                      }
                    else
                      {         
                        if(tr->saveMemory)                 
                          partitionLikelihood =  evaluateGTRGAMMA_GAPPED_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                                              x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                              tip, width, diagptable,
                                                                              x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);                  
                        else
                          partitionLikelihood =  evaluateGTRGAMMA(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                                  x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                  tip, width, diagptable);                                
                      }
#endif
                  }
                  break;                                   
                case 20: /* proteins */
                  {

#ifdef __MIC_NATIVE

                  /* CAT & memory saving are not supported on MIC */

                  assert(!tr->saveMemory);
                  assert(tr->rateHetModel == PLL_GAMMA);

                  if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                    partitionLikelihood =  evaluateGTRGAMMAPROT_LG4_MIC(pr->partitionData[model]->wgt,
                                                                    x1_start, x2_start, pr->partitionData[model]->tipVector_LG4,
                                                                    tip, width, diagptable, pr->partitionData[model]->lg4x_weights);
                  else
                        partitionLikelihood =  evaluateGTRGAMMAPROT_MIC(ex1, ex2, pr->partitionData[model]->wgt,
                                              x1_start, x2_start, pr->partitionData[model]->tipVector,
                                              tip, width, diagptable, fastScaling);

//                  printf("tip: %p, width: %d,  lh: %f\n", tip, width, partitionLikelihood);
//                  int g;
//                  if (x1_start)
//                                        for (g = 0; g < 20; ++g)
//                                                printf("%f \t", x1_start[g]);
//                  printf("\n");
//                  if (x2_start)
//                                        for (g = 0; g < 20; ++g)
//                                                printf("%f \t", x2_start[g]);
#else

                      if(tr->rateHetModel == PLL_CAT)
                      {                           
                        if(tr->saveMemory)
                          partitionLikelihood = evaluateGTRCATPROT_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                        x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                        tip, width, diagptable,  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                        else
                          partitionLikelihood = evaluateGTRCATPROT(fastScaling, ex1, ex2, pr->partitionData[model]->rateCategory, pr->partitionData[model]->wgt,
                                                                   x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                   tip, width, diagptable);               
                      }
                    else
                      {                                               
                        if(tr->saveMemory)
                          partitionLikelihood = evaluateGTRGAMMAPROT_GAPPED_SAVE(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                                                 x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                                 tip, width, diagptable,
                                                                                 x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
                        else
                      {
                        if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                          partitionLikelihood =  evaluateGTRGAMMAPROT_LG4((int *)NULL, (int *)NULL, pr->partitionData[model]->wgt,
                                                                          x1_start, x2_start, pr->partitionData[model]->tipVector_LG4,
                                                                          tip, width, diagptable, PLL_TRUE, pr->partitionData[model]->lg4x_weights);
                        else
                          partitionLikelihood = evaluateGTRGAMMAPROT(fastScaling, ex1, ex2, pr->partitionData[model]->wgt,
                                                                     x1_start, x2_start, pr->partitionData[model]->tipVector,
                                                                     tip, width, diagptable);           
                      }
                      }
#endif
                  }
                  break;                            
                default:
                  assert(0);        
                }
            }
#endif
              
          /* check that there was no major numerical screw-up, the log likelihood should be < 0.0 always */
          
          assert(partitionLikelihood < 0.0);
          
          /* now here is a nasty part, for each partition and each node we maintain an integer counter to count how often 
             how many entries per node were scaled by a constant factor. Here we use this information generated during Felsenstein's 
             pruning algorithm by the newview() functions to undo the preceding scaling multiplications at the root, for mathematical details 
             you should actually read:
             
             A. Stamatakis: "Orchestrating the Phylogenetic Likelihood Function on Emerging Parallel Architectures". 
             In B. Schmidt, editor, Bioinformatics: High Performance Parallel Computer Architectures, 85-115, CRC Press, Taylor & Francis, 2010.
             
             There's a copy of this book in my office 
          */
          
          if(fastScaling)
            partitionLikelihood += (pr->partitionData[model]->globalScaler[pNumber] + pr->partitionData[model]->globalScaler[qNumber]) * log(PLL_MINLIKELIHOOD);
          
          /* now we have the correct log likelihood for the current partition after undoing scaling multiplications */           
          
          /* finally, we also store the per partition log likelihood which is important for optimizing the alpha parameter 
             of this partition for example */

          /* asc bias stuff */

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
          if (tr->threadID == 0 && pr->partitionData[model]->ascBias)
#else
          if (pr->partitionData[model]->ascBias)
#endif
           {
             size_t
               i;
             
             int        
               w = 0;
             
             double                                
               correction;

             switch(tr->rateHetModel)
               {
               case PLL_CAT:
                 {
                   double 
                     rates = 1.0;
                   
                   //need to re-calculate P-matrix for the correction here assuming a rate of 1.0 
                   calcDiagptable(z, states, 1, &rates, pr->partitionData[model]->EIGN, diagptable);
                   
                   
                   correction = evaluateCatAsc(ex1_asc, ex2_asc, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector,
                                               tip, ascWidth, diagptable, ascWidth);
                 }
                 break;
               case PLL_GAMMA:                       
                 correction = evaluateGammaAsc(ex1_asc, ex2_asc, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector,
                                               tip, ascWidth, diagptable, ascWidth);
                 break;
               default:
                 assert(0);
               }
             
             
             
             for(i = (size_t)pr->partitionData[model]->lower; i < (size_t)pr->partitionData[model]->upper; i++)
               w += tr->aliaswgt[i];

             partitionLikelihood = partitionLikelihood - (double)w * log(1.0 - correction);                  
              
           }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
          if(!(pr->partitionData[model]->ascBias && tr->threadID == 0))
           {
#endif
             if(partitionLikelihood >= 0.0)
               {
                 printf("positive log like: %f for partition %d\n", partitionLikelihood, model);
                 assert(0);
               }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
           }
#endif

          
          pr->partitionData[model]->partitionLH = partitionLikelihood;
        }
      else
        {
          /* if the current thread does not have a single site of this partition
             it is important to set the per partition log like to 0.0 because 
             of the reduction operation that will take place later-on.
             That is, the values of tr->perPartitionLH across all threads 
             need to be in a consistent state, always !
          */
          
          if(width == 0)            
            pr->partitionData[model]->partitionLH = 0.0;
        }
    }


#ifdef DEBUG_PERSITE_LNL
  /* per persite-stuff */
  {
    int model = 0; 
    for(model = 0; model < pr->numberOfPartitions ; ++model)
      {
        int j= 0; 
        pInfo *partition  =  pr->partitionData[model]; 
        for(j = 0;  j < partition->width; ++j)
          printf("[%d] lnl[%d]=%f\n", tr->threadID, j, partition->perSiteLikelihoods[j]); 

      }
  }

#endif
}



/** @ingroup evaluateLikelihoodGroup
    @brief Evaluate the log likelihood of the tree topology

    Evaluate the log likelihood of the tree topology of instance \a tr by
    assuming a virtual root between nodes \a p and \a p->back. If
    \a fullTraversal is set to \b PLL_TRUE then the log likelihood vectors for
    each node are recomputed from scratch.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Specifies the virtual root, which is assumed to be a (virtual node) connecting \a p and \a p->back

    @param fullTraversal
      If set to \b PLL_TRUE, then the likelihood vectors at all nodes are recomputed, otherwise only the
      necessary vectors (those that are not oriented in the right direction) are recomputed.

    @param getPerSiteLikelihoods
      Also compute and store (in \a tr->lhs) the log likelihood of each site of the (compressed) alignment

    @note
      If \a getPerSiteLikelihoods is set to \b PLL_TRUE, then make sure that \a tr->fastScaling is set to
      \b PLL_FALSE, otherwise an assertion will fail.
*/
void pllEvaluateLikelihood (pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean fullTraversal, pllBoolean getPerSiteLikelihoods)
{
  /* now this may be the entry point of the library to compute 
     the log like at a branch defined by p and p->back == q */

  volatile double 
    result = 0.0;

  nodeptr 
    q = p->back; 
  

  pllBoolean
        p_recom = PLL_FALSE, /* if one of was missing, we will need to force recomputation */
        q_recom = PLL_FALSE;

  int
    i,
    model,
    numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions : 1;

  /* if evaluate shall return the per-site log likelihoods 
     fastScaling needs to be disabled, otherwise this will 
     not work */

  if(getPerSiteLikelihoods)          
    assert(!(tr->fastScaling)); 

  /* set the first entry of the traversal descriptor to contain the indices
     of nodes p and q */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          

  /* copy the branch lengths of the tree into the first entry of the traversal descriptor.
     if -M is not used tr->numBranches must be 1 */

  for(i = 0; i < numBranches; i++)
    tr->td[0].ti[0].qz[i] =  q->z[i];

  /* recom part */
  if(tr->useRecom)
  {
    int slot = -1;
    if(!isTip(q->number, tr->mxtips))
    {
      q_recom = getxVector(tr->rvec, q->number, &slot, tr->mxtips);
      tr->td[0].ti[0].slot_q = slot;
    }
    if(!isTip(p->number, tr->mxtips))
    {
      p_recom = getxVector(tr->rvec, p->number, &slot, tr->mxtips);
      tr->td[0].ti[0].slot_p = slot;
    }
    if(!isTip(p->number, tr->mxtips) &&  !isTip(q->number, tr->mxtips))
      assert(tr->td[0].ti[0].slot_q != tr->td[0].ti[0].slot_p);
  }


  /* now compute how many conditionals must be re-computed/re-oriented by newview
     to be able to calculate the likelihood at the root defined by p and q.
     */

  /* one entry in the traversal descriptor is already used, hence set the tarversal length counter to 1 */
  tr->td[0].count = 1;

  if(fullTraversal)
  {
    assert(isTip(q->back->number, tr->mxtips));
    computeTraversal(tr, q, PLL_FALSE, numBranches);
  }
  else
  {
    if(p_recom || needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips))
      computeTraversal(tr, p, PLL_TRUE, numBranches);

    if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
      computeTraversal(tr, q, PLL_TRUE, numBranches);
  }


  /* now we copy this partition execute mask into the traversal descriptor which must come from the 
     calling program, the logic of this should not form part of the library */

  storeExecuteMaskInTraversalDescriptor(tr, pr);

  /* also store in the traversal descriptor that something has changed i.e., in the parallel case that the 
     traversal descriptor list of nodes needs to be broadcast once again */

  tr->td[0].traversalHasChanged = PLL_TRUE;
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

  /* now here we enter the fork-join region for Pthreads */


  /* start the parallel region and tell all threads to compute the log likelihood for 
     their fraction of the data. This call is implemented in the case switch of execFunction in axml.c
     */
  if(getPerSiteLikelihoods)
    {
      memset(tr->lhs, 0, sizeof(double) * tr->originalCrunchedLength); 
      pllMasterBarrier(tr, pr, PLL_THREAD_EVALUATE_PER_SITE_LIKES);
    }
  else
    pllMasterBarrier (tr, pr, PLL_THREAD_EVALUATE);

  /* and now here we explicitly do the reduction operation , that is add over the 
     per-thread and per-partition log likelihoods to obtain the overall log like 
     over all sites and partitions */

 
  /* 
     for unpartitioned data that's easy, we just sum over the log likes computed 
     by each thread, thread 0 stores his results in reductionBuffer[0] thread 1 in 
     reductionBuffer[1] and so on 
     */

  /* This reduction for the partitioned case is more complicated because each thread 
     needs to store the partial log like of each partition and we then need to collect 
     and add everything */

#else
  /* and here is just the sequential case, we directly call pllEvaluateIterative() above 
     without having to tell the threads/processes that they need to compute this function now */

  pllEvaluateIterative(tr, pr, getPerSiteLikelihoods); //PLL_TRUE

  /*
    if we want to obtain per-site rates they have initially been stored 
     in arrays that are associated to the partition, now we 
     copy them into the vector tr->lhs[].
     We may also chose that the user needs to rpovide an array, but this can be decided later-on.
  */

  if(getPerSiteLikelihoods) //PLL_TRUE
    {
      for(model = 0; model < pr->numberOfPartitions; model++)
        memcpy(&(tr->lhs[pr->partitionData[model]->lower]), pr->partitionData[model]->perSiteLikelihoods, pr->partitionData[model]->width  * sizeof(double));
    }

#endif

  for(model = 0; model < pr->numberOfPartitions; model++)
    result += pr->partitionData[model]->partitionLH;

  /* set the tree data structure likelihood value to the total likelihood */

  tr->likelihood = result;    

  /* the code below is mainly for testing if the per-site log 
     likelihoods we have stored in tr->lhs yield the same 
     likelihood as the likelihood we computed. 
     For numerical reasons we need to make a dirt PLL_ABS(difference) < epsilon
     comparison */
     
  if(getPerSiteLikelihoods) //PLL_TRUE
    {
      double 
        likelihood = 0;
      int i; 

      /* note that in tr->lhs, we just store the likelihood of 
         one representative of a potentially compressed pattern,
         hence, we need to multiply the elemnts with the pattern 
         weight vector */


      for(i = 0; i < tr->originalCrunchedLength; i++)
        {
//          printf("lhs[%d]=%f * %d\n", i, tr->lhs[i], tr->aliaswgt[i]); 
          likelihood += (tr->lhs[i]   * tr->aliaswgt[i] );
        }
         
      if( PLL_ABS(tr->likelihood - likelihood) > 0.00001)
        {
  //        printf("likelihood was %f\t summed/weighted per-site-lnl was %f\n", tr->likelihood, likelihood); 
        }

        assert(PLL_ABS(tr->likelihood - likelihood) < 0.00001);
    }


  if(tr->useRecom)
  {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
  }

  /* do some bookkeeping to have traversalHasChanged in a consistent state */

  tr->td[0].traversalHasChanged = PLL_FALSE;
}


void perSiteLogLikelihoods(pllInstance *tr, partitionList *pr, double *logLikelihoods)
{
#if (!defined(_USE_PTHREADS) && !defined(_FINE_GRAIN_MPI))
  double 
    //likelihood,
    accumulatedPerSiteLikelihood = 0.0;

  size_t
    localCount,
    i,
    //globalCounter,
    lower,
    upper;
  int model;
#endif
  /* compute the likelihood of the tree with the standard function to:
     1. obtain the current score for error checking
     2. store a full tree traversal in the traversal descriptor that 
     will then be used for calculating per-site log likelihoods 
     for each site individually and independently */

  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

  //likelihood = tr->likelihood;

  /* now compute per-site log likelihoods using the respective functions */

#if (defined( _USE_PTHREADS ) || defined(_FINE_GRAIN_MPI))
  /* here we need a barrier to invoke a parallel region that calls 
     function 
     perSiteLogLikelihoodsPthreads(tree *tr, partitionList *pr, double *lhs, int n, int tid)
     defined above and subsequently collects the per-site log likelihoods 
     computed by the threads and stored in local per-thread memory 
     and stores them in buffer tr->lhs.
     This corresponds to a gather operation in MPI.
     */

  pllMasterBarrier (tr, pr, PLL_THREAD_PER_SITE_LIKELIHOODS);

  /* 
     when the parallel region has terminated, the per-site log likelihoods 
     are stored in array tr->lhs of the master thread which we copy to the result buffer
  */
  
  memcpy(logLikelihoods, tr->lhs, sizeof(double) * tr->originalCrunchedLength);


#else

  /* sequential case: just loop over all partitions and compute per site log likelihoods */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    lower = pr->partitionData[model]->lower;
    upper = pr->partitionData[model]->upper;

    for(i = lower, localCount = 0; i < upper; i++, localCount++)
    {
        double l = 0;
      /* 
         we need to switch of rate heterogeneity implementations here.
         when we have PSR we actually need to provide the per-site rate 
         to the function evaluatePartialGeneric() that computes the 
         per-site log likelihood.
         Under GAMMA, the rate will just be ignored, here we just set it to 1.0
         */

      switch(tr->rateHetModel)
      {
        case PLL_CAT:
          l = evaluatePartialGeneric (tr, pr, i, pr->partitionData[model]->perSiteRates[pr->partitionData[model]->rateCategory[localCount]], model);
          break;
        case PLL_GAMMA:
          l = evaluatePartialGeneric (tr, pr, i, 1.0, model);
          break;
        default:
          assert(0);
      }

      /* store value in result array and add the likelihood of this site to the overall likelihood */

      logLikelihoods[i] = l;
      accumulatedPerSiteLikelihood += l;
    } 
  }


  /* error checking. We need a dirt PLL_ABS() < epsilon here, because the implementations 
     (standard versus per-site) are pretty different and hence slight numerical 
     deviations are expected */

  assert(PLL_ABS(tr->likelihood - accumulatedPerSiteLikelihood) < 0.00001);
  
#endif
  


}

#if (defined(__SSE3) || defined(__AVX))
static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
                                     double *x1_start, double *x2_start, double *tipVector,                   
                                     unsigned char *tipX1, int n, double *diagptable_start, const pllBoolean fastScaling)
{
  double  sum = 0.0, term;       
  int     i;
#if (!defined(__SSE3) && !defined(__AVX))
  int j;  
#endif
  double  *diagptable, *x1, *x2;                            
 
  if(tipX1)
    {          
      for (i = 0; i < n; i++) 
        {
#if (defined(__SSE3) || defined(__AVX))
          PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
#endif
          x1 = &(tipVector[2 * tipX1[i]]);
          x2 = &(x2_start[2 * i]);
          
          diagptable = &(diagptable_start[2 * cptr[i]]);                          
        
#if (defined(__SSE3) || defined(__AVX))
          _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
          
          if(fastScaling)
            term = log(fabs(t[0] + t[1]));
          else
            term = log(fabs(t[0] + t[1])) + (ex2[i] * log(PLL_MINLIKELIHOOD));                           
#else               
          for(j = 0, term = 0.0; j < 2; j++)                         
            term += x1[j] * x2[j] * diagptable[j];            
                 
          if(fastScaling)
            term = log(fabs(term));
          else
            term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));                                                      
#endif    

          sum += wptr[i] * term;
        }       
    }               
  else
    {
      for (i = 0; i < n; i++) 
        {       
#if (defined(__SSE3) || defined(__AVX))
		  PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
#endif                  
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
          
          diagptable = &diagptable_start[2 * cptr[i]];            
#if (defined(__SSE3) || defined(__AVX))
          _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
          
          if(fastScaling)
            term = log(fabs(t[0] + t[1]));
          else
            term = log(fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));                        
#else     
          for(j = 0, term = 0.0; j < 2; j++)
            term += x1[j] * x2[j] * diagptable[j];   
          
          if(fastScaling)
            term = log(fabs(term));
          else
            term = log(fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
#endif
          
          sum += wptr[i] * term;
        }          
    }
       
  return  sum;         
} 


static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
                                      double *x1_start, double *x2_start, 
                                      double *tipVector, 
                                      unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling)
{
  double   sum = 0.0, term;    
  int     i, j;
#if (!defined(__SSE3) && !defined(__AVX))
  int k;
#endif 
  double  *x1, *x2;             

  if(tipX1)
    {          
      for (i = 0; i < n; i++)
        {
#if (defined(__SSE3) || defined(__AVX))
		  PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
          __m128d termv, x1v, x2v, dv;
#endif
          x1 = &(tipVector[2 * tipX1[i]]);       
          x2 = &x2_start[8 * i];                                
#if (defined(__SSE3) || defined(__AVX))
          termv = _mm_set1_pd(0.0);                
          
          for(j = 0; j < 4; j++)
            {
              x1v = _mm_load_pd(&x1[0]);
              x2v = _mm_load_pd(&x2[j * 2]);
              dv   = _mm_load_pd(&diagptable[j * 2]);
              
              x1v = _mm_mul_pd(x1v, x2v);
              x1v = _mm_mul_pd(x1v, dv);
              
              termv = _mm_add_pd(termv, x1v);                 
            }
          
          _mm_store_pd(t, termv);               
          
          if(fastScaling)
            term = log(0.25 * (fabs(t[0] + t[1])));
          else
            term = log(0.25 * (fabs(t[0] + t[1]))) + (ex2[i] * log(PLL_MINLIKELIHOOD));       
#else
          for(j = 0, term = 0.0; j < 4; j++)
            for(k = 0; k < 2; k++)
              term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];                                                
          
          if(fastScaling)
            term = log(0.25 * fabs(term));
          else
            term = log(0.25 * fabs(term)) + ex2[i] * log(PLL_MINLIKELIHOOD);
#endif   
          
          sum += wptr[i] * term;
        }         
    }
  else
    {         
      for (i = 0; i < n; i++) 
        {
#if (defined(__SSE3) || defined(__AVX))
		  PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
          __m128d termv, x1v, x2v, dv;
#endif                            
          x1 = &x1_start[8 * i];
          x2 = &x2_start[8 * i];
                  
#if (defined(__SSE3) || defined(__AVX))
          termv = _mm_set1_pd(0.0);                
          
          for(j = 0; j < 4; j++)
            {
              x1v = _mm_load_pd(&x1[j * 2]);
              x2v = _mm_load_pd(&x2[j * 2]);
              dv   = _mm_load_pd(&diagptable[j * 2]);
              
              x1v = _mm_mul_pd(x1v, x2v);
              x1v = _mm_mul_pd(x1v, dv);
              
              termv = _mm_add_pd(termv, x1v);                 
            }
          
          _mm_store_pd(t, termv);
          
          
          if(fastScaling)
            term = log(0.25 * (fabs(t[0] + t[1])));
          else
            term = log(0.25 * (fabs(t[0] + t[1]))) + ((ex1[i] +ex2[i]) * log(PLL_MINLIKELIHOOD));     
#else     
          for(j = 0, term = 0.0; j < 4; j++)
            for(k = 0; k < 2; k++)
              term += x1[j * 2 + k] * x2[j * 2 + k] * diagptable[j * 2 + k];                                          

          if(fastScaling)
            term = log(0.25 * fabs(term));
          else
            term = log(0.25 * fabs(term)) + (ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD);
#endif

          sum += wptr[i] * term;
        }                       
    }

  return sum;
} 
#endif



/* below are the optimized function versions with geeky intrinsics */

/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree under the GAMMA model of rate heterogeneity and LG4 model of evolution
    
    This is the same as ::evaluateGAMMA_FLEX but for the LG4 model. It contains two implementations,
    one which is the generic, and one that is optimized with SSE3 instructions. The two implementations
    are separated by preprocessor macros.
    The difference from ::evaluateGAMMA_FLEX is that we have 4 different tipVectors computed from the 4 different
    Q matrix decompositions.
    Please check ::evaluateGAMMA_FLEX for more information and a description of the common
    input parameters.
*/
static double evaluateGTRGAMMAPROT_LG4(int *ex1, int *ex2, int *wptr,
                                       double *x1, double *x2,  
                                       double *tipVector[4], 
                                       unsigned char *tipX1, int n, double *diagptable, const pllBoolean fastScaling,
                                       double * lg4_weights)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
        {
#if (defined(__SSE3) || defined(__AVX))
          __m128d tv = _mm_setzero_pd();
                                  
          for(j = 0, term = 0.0; j < 4; j++)
            {
              double *d = &diagptable[j * 20];

              __m128d
              	  t = _mm_setzero_pd(),
              	  w = _mm_set1_pd(lg4_weights[j]);

              left = &(tipVector[j][20 * tipX1[i]]);
              right = &(x2[80 * i + 20 * j]);
              for(l = 0; l < 20; l+=2)
                {
                  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
                  t = _mm_add_pd(t, _mm_mul_pd(mul, _mm_load_pd(&d[l])));
                }
              tv = _mm_add_pd(tv, _mm_mul_pd(t, w));
            }

          tv = _mm_hadd_pd(tv, tv);
          _mm_storel_pd(&term, tv);
          

#else                             
          for(j = 0, term = 0.0; j < 4; j++)
            {
        	  double t = 0.0;

              left = &(tipVector[j][20 * tipX1[i]]);
              right = &(x2[80 * i + 20 * j]);

              for(l = 0; l < 20; l++)
                t += left[l] * right[l] * diagptable[j * 20 + l];

              term += lg4_weights[j] * t;
            }     
#endif
          
          if(fastScaling)
            term = log(fabs(term));
          else
            term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));

          sum += wptr[i] * term;

        }               
    }              
  else
    {
      for (i = 0; i < n; i++) 
        {                                    
#if (defined(__SSE3) || defined(__AVX))
          __m128d tv = _mm_setzero_pd();                          
              
          for(j = 0, term = 0.0; j < 4; j++)
            {
              double *d = &diagptable[j * 20];

              __m128d
              t = _mm_setzero_pd(),
              w = _mm_set1_pd(lg4_weights[j]);

              left  = &(x1[80 * i + 20 * j]);
              right = &(x2[80 * i + 20 * j]);
              
              for(l = 0; l < 20; l+=2)
                {
                  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
                  t = _mm_add_pd(t, _mm_mul_pd(mul, _mm_load_pd(&d[l])));
                }
              tv = _mm_add_pd(tv, _mm_mul_pd(t, w));
            }
          tv = _mm_hadd_pd(tv, tv);
          _mm_storel_pd(&term, tv);       
#else
          for(j = 0, term = 0.0; j < 4; j++)
            {
        	  double t = 0.0;

              left  = &(x1[80 * i + 20 * j]);
              right = &(x2[80 * i + 20 * j]);       
              
              for(l = 0; l < 20; l++)
                t += left[l] * right[l] * diagptable[j * 20 + l];

              term += lg4_weights[j] * t;
            }
#endif
          
          if(fastScaling)
            term = log(fabs(term));
          else
            term = log(fabs(term)) + ((ex1[i] + ex2[i])*log(PLL_MINLIKELIHOOD));
          
          sum += wptr[i] * term;
        }         
    }

  return  sum;
}

#if (defined(__SSE3) || defined(__AVX))
/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b GAMMA model of rate heterogeneity 
    and the memory saving technique (Optimized SSE3 version for AA data)
 
    This is the SSE3 optimized version of ::evaluateGAMMA_FLEX_SAVE for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateGAMMA_FLEX_SAVE for more information and
    a description of the input parameters
*/
static double evaluateGTRGAMMAPROT_GAPPED_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                                double *x1, double *x2,  
                                                double *tipVector, 
                                                unsigned char *tipX1, int n, double *diagptable, 
                                                double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)                                    
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  
    *left, 
    *right,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x1v,
    *x2v;              
  __m128d tv;

  if(tipX1)
  {               
    for (i = 0; i < n; i++) 
    {
      if(x2_gap[i / 32] & mask32[i % 32])
        x2v = x2_gapColumn;
      else
      {
        x2v = x2_ptr;
        x2_ptr += 80;
      }

	  //TUNG: Standard C does not allow declaration after executable statement
	  tv = _mm_setzero_pd();
      //__m128d tv = _mm_setzero_pd();
      left = &(tipVector[20 * tipX1[i]]);                 

      for(j = 0, term = 0.0; j < 4; j++)
      {
        double *d = &diagptable[j * 20];
        right = &(x2v[20 * j]);
        for(l = 0; l < 20; l+=2)
        {
          __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
          tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));                
        }                               
      }

      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);


      if(!fastScaling)
        term = log(0.25 * fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));    

      sum += wptr[i] * term;
    }                   
  }              
  else
  {
    for (i = 0; i < n; i++) 
    {
      if(x1_gap[i / 32] & mask32[i % 32])
        x1v = x1_gapColumn;
      else
      {
        x1v = x1_ptr;
        x1_ptr += 80;
      }

      if(x2_gap[i / 32] & mask32[i % 32])
        x2v = x2_gapColumn;
      else
      {
        x2v = x2_ptr;
        x2_ptr += 80;
      }

      //__m128d tv = _mm_setzero_pd(); 
	  tv = _mm_setzero_pd();

      for(j = 0, term = 0.0; j < 4; j++)
      {
        double *d = &diagptable[j * 20];
        left  = &(x1v[20 * j]);
        right = &(x2v[20 * j]);

        for(l = 0; l < 20; l+=2)
        {
          __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
          tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));                
        }                               
      }
      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);   


       if(!fastScaling)
        term = log(0.25 * fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));


      sum += wptr[i] * term;
    }         
  }

  return  sum;
}



/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b GAMMA model of rate heterogeneity 
    (Optimized SSE3 version for AA data)
 
    This is the SSE3 optimized version of ::evaluateGAMMA_FLEX for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateGAMMA_FLEX for more information and
    a description of the common input parameters
*/
static double evaluateGTRGAMMAPROT (const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                    double *x1, double *x2,  
                                    double *tipVector, 
                                    unsigned char *tipX1, int n, double *diagptable)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              

  if(tipX1)
  {               
    for (i = 0; i < n; i++) 
    {

      __m128d tv = _mm_setzero_pd();
      left = &(tipVector[20 * tipX1[i]]);                 

      for(j = 0, term = 0.0; j < 4; j++)
      {
        double *d = &diagptable[j * 20];
        right = &(x2[80 * i + 20 * j]);
        for(l = 0; l < 20; l+=2)
        {
          __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
          tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));                
        }                               
      }
      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);


      if(!fastScaling)
        term = log(0.25 * fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));


      sum += wptr[i] * term;
    }                   
  }              
  else
  {
    for (i = 0; i < n; i++) 
    {                                
      __m128d tv = _mm_setzero_pd();                      

      for(j = 0, term = 0.0; j < 4; j++)
      {
        double *d = &diagptable[j * 20];
        left  = &(x1[80 * i + 20 * j]);
        right = &(x2[80 * i + 20 * j]);

        for(l = 0; l < 20; l+=2)
        {
          __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
          tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));                
        }                               
      }
      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);   


       if(!fastScaling)
        term = log(0.25 * fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(term));


      sum += wptr[i] * term;
    }
  }

  return  sum;
}


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b CAT model of rate heterogeneity 
    (Optimized SSE3 version for AA data)
 
    This is the SSE3 optimized version of ::evaluateCAT_FLEX for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateCAT_FLEX for more information and
    a description of the common input parameters
*/
static double evaluateGTRCATPROT (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                  double *x1, double *x2, double *tipVector,
                                  unsigned char *tipX1, int n, double *diagptable_start)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  __m128d tv;

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {           
      left = &(tipVector[20 * tipX1[i]]);
      right = &(x2[20 * i]);

      diagptable = &diagptable_start[20 * cptr[i]];                      

	  //TUNG: Standard C does not allow declaration after executable statement
	  tv = _mm_setzero_pd();
      //__m128d tv = _mm_setzero_pd();        

      for(l = 0; l < 20; l+=2)
      {
        __m128d lv = _mm_load_pd(&left[l]);
        __m128d rv = _mm_load_pd(&right[l]);
        __m128d mul = _mm_mul_pd(lv, rv);
        __m128d dv = _mm_load_pd(&diagptable[l]);

        tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));                  
      }                         

      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);

      if(!fastScaling)
        term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));

      sum += wptr[i] * term;
    }      
  }    
  else
  {

    for (i = 0; i < n; i++) 
    {                                 
      left  = &x1[20 * i];
      right = &x2[20 * i];

      diagptable = &diagptable_start[20 * cptr[i]];             

      __m128d tv = _mm_setzero_pd();        

      for(l = 0; l < 20; l+=2)
      {
        __m128d lv = _mm_load_pd(&left[l]);
        __m128d rv = _mm_load_pd(&right[l]);
        __m128d mul = _mm_mul_pd(lv, rv);
        __m128d dv = _mm_load_pd(&diagptable[l]);

        tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));                  
      }                         

      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);

      if(!fastScaling)
        term = log(fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));  

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b CAT model of rate heterogeneity with memory saving 
    (Optimized SSE3 version for AA data)
 
    This is the SSE3 optimized version of ::evaluateCAT_FLEX_SAVE for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateCAT_FLEX_SAVE for more information and
    a description of the common input parameters
*/
static double evaluateGTRCATPROT_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                       double *x1, double *x2, double *tipVector,
                                       unsigned char *tipX1, int n, double *diagptable_start, 
                                       double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
        term,
        *diagptable,  
        *left, 
        *right,
        *left_ptr = x1,
        *right_ptr = x2;

  int     
    i, 
    l;                           

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {           
      left = &(tipVector[20 * tipX1[i]]);

      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
      {
        right = right_ptr;
        right_ptr += 20;
      }          

      diagptable = &diagptable_start[20 * cptr[i]];                      

      __m128d tv = _mm_setzero_pd();        

      for(l = 0; l < 20; l+=2)
      {
        __m128d lv = _mm_load_pd(&left[l]);
        __m128d rv = _mm_load_pd(&right[l]);
        __m128d mul = _mm_mul_pd(lv, rv);
        __m128d dv = _mm_load_pd(&diagptable[l]);

        tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));                  
      }                         

      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);

      if(!fastScaling)
        term = log(fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));

      sum += wptr[i] * term;
    }      
  }    
  else
  {

    for (i = 0; i < n; i++) 
    {                                     
      if(isGap(x1_gap, i))
        left = x1_gapColumn;
      else
      {
        left = left_ptr;
        left_ptr += 20;
      }

      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
      {
        right = right_ptr;
        right_ptr += 20;
      }

      diagptable = &diagptable_start[20 * cptr[i]];             

      __m128d tv = _mm_setzero_pd();        

      for(l = 0; l < 20; l+=2)
      {
        __m128d lv = _mm_load_pd(&left[l]);
        __m128d rv = _mm_load_pd(&right[l]);
        __m128d mul = _mm_mul_pd(lv, rv);
        __m128d dv = _mm_load_pd(&diagptable[l]);

        tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));                  
      }                         

      tv = _mm_hadd_pd(tv, tv);
      _mm_storel_pd(&term, tv);

      if(!fastScaling)
        term = log(fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(term));  

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b CAT model of rate heterogeneity with memory saving 
    (Optimized SSE3 version for DNA data)
 
    This is the SSE3 optimized version of ::evaluateCAT_FLEX_SAVE for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateCAT_FLEX_SAVE for more information and
    a description of the common input parameters
*/
static double evaluateGTRCAT_SAVE (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                                   double *x1_start, double *x2_start, double *tipVector,                     
                                   unsigned char *tipX1, int n, double *diagptable_start,
                                   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double  sum = 0.0, term;       
  int     i;

  double  *diagptable, 
          *x1, 
          *x2,
          *x1_ptr = x1_start,
          *x2_ptr = x2_start;

  if(tipX1)
  {           
    for (i = 0; i < n; i++) 
    {   
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      x1 = &(tipVector[4 * tipX1[i]]);

      if(isGap(x2_gap, i))
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;
        x2_ptr += 4;
      }

      diagptable = &diagptable_start[4 * cptr[i]];

      x1v1 =  _mm_load_pd(&x1[0]);
      x1v2 =  _mm_load_pd(&x1[2]);
      x2v1 =  _mm_load_pd(&x2[0]);
      x2v2 =  _mm_load_pd(&x2[2]);
      dv1  =  _mm_load_pd(&diagptable[0]);
      dv2  =  _mm_load_pd(&diagptable[2]);

      x1v1 = _mm_mul_pd(x1v1, x2v1);
      x1v1 = _mm_mul_pd(x1v1, dv1);

      x1v2 = _mm_mul_pd(x1v2, x2v2);
      x1v2 = _mm_mul_pd(x1v2, dv2);

      x1v1 = _mm_add_pd(x1v1, x1v2);

      _mm_store_pd(t, x1v1);

      if(!fastScaling)
        term = log(fabs(t[0] + t[1])) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(t[0] + t[1]));



      sum += wptr[i] * term;
    }   
  }               
  else
  {
    for (i = 0; i < n; i++) 
    { 
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      if(isGap(x1_gap, i))
        x1 = x1_gapColumn;
      else
      {
        x1 = x1_ptr;
        x1_ptr += 4;
      }

      if(isGap(x2_gap, i))
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;
        x2_ptr += 4;
      }

      diagptable = &diagptable_start[4 * cptr[i]];      

      x1v1 =  _mm_load_pd(&x1[0]);
      x1v2 =  _mm_load_pd(&x1[2]);
      x2v1 =  _mm_load_pd(&x2[0]);
      x2v2 =  _mm_load_pd(&x2[2]);
      dv1  =  _mm_load_pd(&diagptable[0]);
      dv2  =  _mm_load_pd(&diagptable[2]);

      x1v1 = _mm_mul_pd(x1v1, x2v1);
      x1v1 = _mm_mul_pd(x1v1, dv1);

      x1v2 = _mm_mul_pd(x1v2, x2v2);
      x1v2 = _mm_mul_pd(x1v2, dv2);

      x1v1 = _mm_add_pd(x1v1, x1v2);

      _mm_store_pd(t, x1v1);


       if(!fastScaling)
        term = log(fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(t[0] + t[1]));

      sum += wptr[i] * term;
    }    
  }

  return  sum;         
} 


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b GAMMA model of rate heterogeneity with memory saving 
    (Optimized SSE3 version for DNA data)
 
    This is the SSE3 optimized version of ::evaluateGAMMA_FLEX_SAVE for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateGAMMA_FLEX_SAVE for more information and
    a description of the common input parameters
*/
static double evaluateGTRGAMMA_GAPPED_SAVE(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                                           double *x1_start, double *x2_start, 
                                           double *tipVector, 
                                           unsigned char *tipX1, const int n, double *diagptable,
                                           double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   sum = 0.0, term;    
  int     i, j;
  double  
    *x1, 
    *x2,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;



  if(tipX1)
  {        


    for (i = 0; i < n; i++)
    {
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d termv, x1v, x2v, dv;

      x1 = &(tipVector[4 * tipX1[i]]);   
      if(x2_gap[i / 32] & mask32[i % 32])
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;     
        x2_ptr += 16;
      }


      termv = _mm_set1_pd(0.0);            

      for(j = 0; j < 4; j++)
      {
        x1v = _mm_load_pd(&x1[0]);
        x2v = _mm_load_pd(&x2[j * 4]);
        dv   = _mm_load_pd(&diagptable[j * 4]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);

        x1v = _mm_load_pd(&x1[2]);
        x2v = _mm_load_pd(&x2[j * 4 + 2]);
        dv   = _mm_load_pd(&diagptable[j * 4 + 2]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);
      }

      _mm_store_pd(t, termv);            

       if(!fastScaling)
        term = log(0.25 * fabs(t[0] + t[1])) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(t[0] + t[1]));


      sum += wptr[i] * term;
    }     
  }
  else
  {        

    for (i = 0; i < n; i++) 
    {
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d termv, x1v, x2v, dv;

      if(x1_gap[i / 32] & mask32[i % 32])
        x1 = x1_gapColumn;
      else
      {
        x1 = x1_ptr;              
        x1_ptr += 16;
      }

      if(x2_gap[i / 32] & mask32[i % 32])
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;
        x2_ptr += 16;
      }

      termv = _mm_set1_pd(0.0);          

      for(j = 0; j < 4; j++)
      {
        x1v = _mm_load_pd(&x1[j * 4]);
        x2v = _mm_load_pd(&x2[j * 4]);
        dv   = _mm_load_pd(&diagptable[j * 4]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);

        x1v = _mm_load_pd(&x1[j * 4 + 2]);
        x2v = _mm_load_pd(&x2[j * 4 + 2]);
        dv   = _mm_load_pd(&diagptable[j * 4 + 2]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);
      }

      _mm_store_pd(t, termv);

      if(!fastScaling)
        term = log(0.25 * fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(t[0] + t[1]));


      sum += wptr[i] * term;
    }                           
  }

  return sum;
} 


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b GAMMA model of rate heterogeneity (Optimized SSE3 version for DNA data)
 
    This is the SSE3 optimized version of ::evaluateGAMMA_FLEX for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateGAMMA_FLEX for more information and
    a description of the common input parameters
*/
static double evaluateGTRGAMMA(const pllBoolean fastScaling, int *ex1, int *ex2, int *wptr,
                               double *x1_start, double *x2_start, 
                               double *tipVector, 
                               unsigned char *tipX1, const int n, double *diagptable)
{
  double   sum = 0.0, term;    
  int     i, j;

  double  *x1, *x2;             



  if(tipX1)
  {             
    for (i = 0; i < n; i++)
    {
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d termv, x1v, x2v, dv;

      x1 = &(tipVector[4 * tipX1[i]]);   
      x2 = &x2_start[16 * i];    


      termv = _mm_set1_pd(0.0);            

      for(j = 0; j < 4; j++)
      {
        x1v = _mm_load_pd(&x1[0]);
        x2v = _mm_load_pd(&x2[j * 4]);
        dv   = _mm_load_pd(&diagptable[j * 4]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);

        x1v = _mm_load_pd(&x1[2]);
        x2v = _mm_load_pd(&x2[j * 4 + 2]);
        dv   = _mm_load_pd(&diagptable[j * 4 + 2]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);
      }

      _mm_store_pd(t, termv);


       if(!fastScaling)
        term = log(0.25 * fabs(t[0] + t[1])) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(t[0] + t[1]));



      sum += wptr[i] * term;
    }     
  }
  else
  {        
    for (i = 0; i < n; i++) 
    {
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d termv, x1v, x2v, dv;


      x1 = &x1_start[16 * i];
      x2 = &x2_start[16 * i];             


      termv = _mm_set1_pd(0.0);          

      for(j = 0; j < 4; j++)
      {
        x1v = _mm_load_pd(&x1[j * 4]);
        x2v = _mm_load_pd(&x2[j * 4]);
        dv   = _mm_load_pd(&diagptable[j * 4]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);

        x1v = _mm_load_pd(&x1[j * 4 + 2]);
        x2v = _mm_load_pd(&x2[j * 4 + 2]);
        dv   = _mm_load_pd(&diagptable[j * 4 + 2]);

        x1v = _mm_mul_pd(x1v, x2v);
        x1v = _mm_mul_pd(x1v, dv);

        termv = _mm_add_pd(termv, x1v);
      }

      _mm_store_pd(t, termv);

      if(!fastScaling)
        term = log(0.25 * fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(0.25 * fabs(t[0] + t[1]));



      sum += wptr[i] * term;
    }                           
  }

  return sum;
} 


/** @ingroup evaluateLikelihoodGroup
    @brief Evaluation of log likelihood of a tree using the \b CAT model of rate heterogeneity (Optimized SSE3 version for DNA data)
 
    This is the SSE3 optimized version of ::evaluateCAT_FLEX for evaluating the log
    likelihood at some edge whose two end-points (nodes) have the conditional likelihood
    vectors \a x1 and \a x2. Please check ::evaluateCAT_FLEX for more information and
    a description of the common input parameters
*/
static double evaluateGTRCAT (const pllBoolean fastScaling, int *ex1, int *ex2, int *cptr, int *wptr,
                              double *x1_start, double *x2_start, double *tipVector,                  
                              unsigned char *tipX1, int n, double *diagptable_start)
{
  double  sum = 0.0, term;       
  int     i;

  double  *diagptable, *x1, *x2;                            

  if(tipX1)
  {           
    for (i = 0; i < n; i++) 
    {   
    	PLL_ALIGN_BEGIN	double t[2] PLL_ALIGN_END;
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      x1 = &(tipVector[4 * tipX1[i]]);
      x2 = &x2_start[4 * i];

      diagptable = &diagptable_start[4 * cptr[i]];


      x1v1 =  _mm_load_pd(&x1[0]);
      x1v2 =  _mm_load_pd(&x1[2]);
      x2v1 =  _mm_load_pd(&x2[0]);
      x2v2 =  _mm_load_pd(&x2[2]);
      dv1  =  _mm_load_pd(&diagptable[0]);
      dv2  =  _mm_load_pd(&diagptable[2]);

      x1v1 = _mm_mul_pd(x1v1, x2v1);
      x1v1 = _mm_mul_pd(x1v1, dv1);

      x1v2 = _mm_mul_pd(x1v2, x2v2);
      x1v2 = _mm_mul_pd(x1v2, dv2);

      x1v1 = _mm_add_pd(x1v1, x1v2);

      _mm_store_pd(t, x1v1);

       if(!fastScaling)
        term = log(fabs(t[0] + t[1])) + (ex2[i] * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(t[0] + t[1]));


      sum += wptr[i] * term;
    }   
  }               
  else
  {
    for (i = 0; i < n; i++) 
    { 
    	PLL_ALIGN_BEGIN double t[2] PLL_ALIGN_END;
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      x1 = &x1_start[4 * i];
      x2 = &x2_start[4 * i];

      diagptable = &diagptable_start[4 * cptr[i]];      


      x1v1 =  _mm_load_pd(&x1[0]);
      x1v2 =  _mm_load_pd(&x1[2]);
      x2v1 =  _mm_load_pd(&x2[0]);
      x2v2 =  _mm_load_pd(&x2[2]);
      dv1  =  _mm_load_pd(&diagptable[0]);
      dv2  =  _mm_load_pd(&diagptable[2]);

      x1v1 = _mm_mul_pd(x1v1, x2v1);
      x1v1 = _mm_mul_pd(x1v1, dv1);

      x1v2 = _mm_mul_pd(x1v2, x2v2);
      x1v2 = _mm_mul_pd(x1v2, dv2);

      x1v1 = _mm_add_pd(x1v1, x1v2);

      _mm_store_pd(t, x1v1);

      if(!fastScaling)
        term = log(fabs(t[0] + t[1])) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
      else
        term = log(fabs(t[0] + t[1]));


      sum += wptr[i] * term;
    }    
  }

  return  sum;         
} 





#endif
