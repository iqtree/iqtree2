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
 * @file bipartitionList.c
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

#ifdef __SSE3
#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
/*#include <tmmintrin.h>*/
#endif
#endif

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif


/** @file makenewzGenericSpecial.c
 *  
 *  @brief Branch length optimization
 */



/* pointers to reduction buffers for storing and gathering the first and second derivative 
   of the likelihood in Pthreads and MPI */

#ifdef IS_PARALLEL
void branchLength_parallelReduce(pllInstance *tr, double *dlnLdlz,  double *d2lnLdlz2, int numBranches ) ;
//extern double *globalResult;
#endif


extern const unsigned int mask32[32];

#if (defined(__SSE3) || defined(__AVX))
static void sumGAMMA_BINARY(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
                            unsigned char *tipX1, unsigned char *tipX2, int n);
static void coreGTRGAMMA_BINARY(const int upper, double *sumtable,
                                volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr);
static void coreGTRCAT_BINARY(int upper, int numberOfCategories, double *sum,
                              volatile double *d1, volatile double *d2, 
                              double *rptr, double *EIGN, int *cptr, double lz, int *wgt);
static void sumCAT_BINARY(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
                          unsigned char *tipX1, unsigned char *tipX2, int n);
#endif

/*******************/


/* generic function to get the required pointers to the data associated with the left and right node that define a branch */

static void getVects(pllInstance *tr, 
                     partitionList *pr, 
                     unsigned char **tipX1, unsigned char **tipX2, 
                     double **x1_start, double **x2_start, 
                     int *tipCase, 
                     int model, 
                     double **x1_gapColumn, double **x2_gapColumn, 
                     unsigned int **x1_gap, unsigned int **x2_gap,
                     double ** x1_start_asc,
                     double ** x2_start_asc)
{
  int    
    rateHet = (int)discreteRateCategories(tr->rateHetModel),
            states = pr->partitionData[model]->states,
            pNumber, 
            qNumber; 

  /* get the left and right node number of the nodes defining the branch we want to optimize */

  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;

  /* get the index where the ancestral vector is expected to be found */
  int p_slot, q_slot;
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
   

  /* initialize to NULL */

  *x1_start = (double*)NULL,
  *x2_start = (double*)NULL;
  
  *tipX1 = (unsigned char*)NULL,
  *tipX2 = (unsigned char*)NULL;

  *x1_start_asc = NULL;
  *x2_start_asc = NULL;

  /* switch over the different tip cases again here */

  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
  {      
    if(!( isTip(pNumber, tr->mxtips) && isTip(qNumber, tr->mxtips)) )
    {
      *tipCase = PLL_TIP_INNER;
      if(isTip(qNumber, tr->mxtips))
      {
        *tipX1 = pr->partitionData[model]->yVector[qNumber];
        *x2_start = pr->partitionData[model]->xVector[p_slot];

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
        if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
          if(pr->partitionData[model]->ascBias)
#endif
          {
            *x2_start_asc = &pr->partitionData[model]->ascVector[(pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
          }

        if(tr->saveMemory)
        {
          *x2_gap = &(pr->partitionData[model]->gapVector[pNumber * pr->partitionData[model]->gapVectorLength]);
          *x2_gapColumn   = &pr->partitionData[model]->gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
        }
      }
      else
      {
        *tipX1 = pr->partitionData[model]->yVector[pNumber];
        *x2_start = pr->partitionData[model]->xVector[q_slot];

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
        if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
          if(pr->partitionData[model]->ascBias)
#endif  
          {
            *x2_start_asc = &pr->partitionData[model]->ascVector[(qNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
          }

        if(tr->saveMemory)
        {
          *x2_gap = &(pr->partitionData[model]->gapVector[qNumber * pr->partitionData[model]->gapVectorLength]);
          *x2_gapColumn   = &pr->partitionData[model]->gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
        }
      }
    }
    else
    {
      /* note that tip tip should normally not occur since this means that we are trying to optimize 
         a branch in a two-taxon tree. However, this has been inherited be some RAxML function 
         that optimized pair-wise distances between all taxa in a tree */

      *tipCase = PLL_TIP_TIP;
      *tipX1 = pr->partitionData[model]->yVector[pNumber];
      *tipX2 = pr->partitionData[model]->yVector[qNumber];
    }
  }
  else
  {
    *tipCase = PLL_INNER_INNER;

    *x1_start = pr->partitionData[model]->xVector[p_slot];
    *x2_start = pr->partitionData[model]->xVector[q_slot];

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
        if(pr->partitionData[model]->ascBias)
#endif
        {
          *x1_start_asc = &pr->partitionData[model]->ascVector[(pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
          *x2_start_asc = &pr->partitionData[model]->ascVector[(qNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
        }           
    if(tr->saveMemory)
    {
      *x1_gap = &(pr->partitionData[model]->gapVector[pNumber * pr->partitionData[model]->gapVectorLength]);
      *x1_gapColumn   = &pr->partitionData[model]->gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];

      *x2_gap = &(pr->partitionData[model]->gapVector[qNumber * pr->partitionData[model]->gapVectorLength]);
      *x2_gapColumn   = &pr->partitionData[model]->gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
    }
  }

}


/* this is actually a pre-computation and storage of values that remain constant while we change the value of the branch length 
   we want to adapt. the target pointer sumtable is a single pre-allocated array that has the same 
   size as a conditional likelihood vector at an inner node.

   So if we want to do a Newton-Raphson optimization we only execute this function once in the beginning for each new branch we are considering !
   */

#if (!defined(__SSE3) && !defined(__AVX))
static void sumCAT_FLEX(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, const int states)
{
  int 
    i, 
    l;

  double 
    *sum, 
    *left, 
    *right;

  switch(tipCase)
  {

    /* switch over possible configurations of the nodes p and q defining the branch */

    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
      {
        left  = &(tipVector[states * tipX1[i]]);
        right = &(tipVector[states * tipX2[i]]);
        sum = &sumtable[states * i];

        /* just multiply the values with each other for each site, note the similarity with evaluate() 
           we precompute the product which will remain constant and then just multiply this pre-computed 
           product with the changing P matrix exponentaions that depend on the branch lengths */

        for(l = 0; l < states; l++)
          sum[l] = left[l] * right[l];
      }
      break;
    case PLL_TIP_INNER:

      /* same as for PLL_TIP_TIP only that 
         we now access on tip vector and one 
         inner vector. 

         You may also observe that we do not consider using scaling vectors anywhere here.

         This is because we are interested in the first and second derivatives of the likelihood and 
         hence the addition of the log() of the scaling factor times the number of scaling events
         becomes obsolete through the derivative */

      for (i = 0; i < n; i++)
      {
        left = &(tipVector[states * tipX1[i]]);
        right = &x2[states * i];
        sum = &sumtable[states * i];

        for(l = 0; l < states; l++)
          sum[l] = left[l] * right[l];
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        left  = &x1[states * i];
        right = &x2[states * i];
        sum = &sumtable[states * i];

        for(l = 0; l < states; l++)
          sum[l] = left[l] * right[l];
      }
      break;
    default:
      assert(0);
  }
}
#endif



#if (!defined(__SSE3) && !defined(__AVX))

/* same thing for GAMMA models. The only noteworthy thing here is that we have an additional inner loop over the 
   number of discrete gamma rates. The data access pattern is also different since for tip vector accesses through our 
   lookup table, we do not distnguish between rates 

   Note the different access pattern in PLL_TIP_INNER:

   left = &(tipVector[states * tipX1[i]]);        
   right = &(x2[span * i + l * states]);

*/

static void sumGAMMA_FLEX(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, const int states)
{
  int 
    i, 
    l, 
    k;

  const int 
    span = 4 * states;

  double 
    *left, 
    *right, 
    *sum;




  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
      {
        left  = &(tipVector[states * tipX1[i]]);
        right = &(tipVector[states * tipX2[i]]);

        for(l = 0; l < 4; l++)
        {
          sum = &sumtable[i * span + l * states];

          for(k = 0; k < states; k++)
            sum[k] = left[k] * right[k];

        }
      }
      break;
    case PLL_TIP_INNER:
      //reorder_back( x2, n, span );
      for(i = 0; i < n; i++)
      {
        left = &(tipVector[states * tipX1[i]]);

        for(l = 0; l < 4; l++)
        {
          right = &(x2[span * i + l * states]);
          sum = &sumtable[i * span + l * states];

          for(k = 0; k < states; k++)
            sum[k] = left[k] * right[k];

        }
      }
      //reorder( x2, n, span );
      break;
    case PLL_INNER_INNER:
      //reorder_back( x1, n, span );
      //reorder_back( x2, n, span );
      for(i = 0; i < n; i++)
      {
        for(l = 0; l < 4; l++)
        {
          left  = &(x1[span * i + l * states]);
          right = &(x2[span * i + l * states]);
          sum   = &(sumtable[i * span + l * states]);


          for(k = 0; k < states; k++)
            sum[k] = left[k] * right[k];
        }
      }
      //reorder( x1, n, span );
      //reorder( x2, n, span );
      break;
    default:
      assert(0);
  }
}
#endif

/* optimized functions for branch length optimization */


#if (defined(__SSE3) || defined(__AVX))

static void sumCAT_SAVE(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static void sumGAMMA_GAPPED_SAVE(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGAMMAPROT_GAPPED_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static void sumGAMMAPROT_LG4(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector[4],
                             unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGTRCATPROT_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static void coreGTRGAMMAPROT_LG4(double *gammaRates, double *EIGN[4], double *sumtable, int upper, int *wrptr,
                                 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz,
                                 double * lg4_weights);

static void coreGTRGAMMA(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);

static void coreGTRCAT(int upper, int numberOfCategories, double *sum,
    volatile double *d1, volatile double *d2, int *wgt, 
    double *rptr, double *EIGN, int *cptr, double lz);


static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz);

static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
    int *wgt, volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable);

#endif


/* now this is the core function of the newton-Raphson based branch length optimization that actually computes 
   the first and second derivative of the likelihood given a new proposed branch length lz */

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

static double coreCatAsc(double *EIGN, double *sumtable, int upper,
			 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, const int numStates,
			 double *ascScaler)
{
  double  
    diagptable[1024], 
    lh = 0.0,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr;

  int     
    i,     
    l;  

 
  ki = 1.0;
  kisqr = 1.0;

  for(l = 1; l < numStates; l++)
    {
      diagptable[l * 4]     = exp(EIGN[l-1] * ki * lz);
      diagptable[l * 4 + 1] = EIGN[l-1] * ki;
      diagptable[l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      double
	*sum = &sumtable[i * numStates],
	tmp,
	inv_Li   = 0.0,
	dlnLidlz = 0.0,
	d2lnLidlz2 = 0.0;

    
      inv_Li += sum[0];

      for(l = 1; l < numStates; l++)
	{
	  inv_Li     += (tmp = diagptable[l * 4] * sum[l]);
	  dlnLidlz   += tmp * diagptable[l * 4 + 1];
	  d2lnLidlz2 += tmp * diagptable[l * 4 + 2];
	}	            
            
      inv_Li = fabs(inv_Li);             
       
      lh        += inv_Li * ascScaler[i];
      dlnLdlz   += dlnLidlz * ascScaler[i];
      d2lnLdlz2 += d2lnLidlz2 * ascScaler[i];
    } 

  *ext_dlnLdlz   = (dlnLdlz / (lh - 1.0));
  *ext_d2lnLdlz2 = (((lh - 1.0) * (d2lnLdlz2) - (dlnLdlz * dlnLdlz)) / ((lh - 1.0) * (lh - 1.0)));  

  return lh;
}


static double coreGammaAsc(double *gammaRates, double *EIGN, double *sumtable, int upper,
			   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, const int numStates,
			   double *ascScaler)
{
  double  
    diagptable[1024], 
    lh = 0.0,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr;

  int     
    i, 
    j, 
    l;  

  const int 
    gammaStates = 4 * numStates;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < numStates; l++)
	{
	  diagptable[i * gammaStates + l * 4]     = exp(EIGN[l-1] * ki * lz);
	  diagptable[i * gammaStates + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * gammaStates + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      double
	*sum = &sumtable[i * gammaStates],
	tmp,
	inv_Li   = 0.0,
	dlnLidlz = 0.0,
	d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * numStates];

	  for(l = 1; l < numStates; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * gammaStates + l * 4] * sum[j * numStates + l]);
	      dlnLidlz   += tmp * diagptable[j * gammaStates + l * 4 + 1];
	      d2lnLidlz2 += tmp * diagptable[j * gammaStates + l * 4 + 2];
	    }	  
	}    
            
      inv_Li = 0.25 * fabs(inv_Li);         
      dlnLidlz *= 0.25;
      d2lnLidlz2 *= 0.25;
       
      lh        += inv_Li * ascScaler[i];
      dlnLdlz   += dlnLidlz * ascScaler[i];
      d2lnLdlz2 += d2lnLidlz2 * ascScaler[i];
    } 

  *ext_dlnLdlz   = (dlnLdlz / (lh - 1.0));
  *ext_d2lnLdlz2 = (((lh - 1.0) * (d2lnLdlz2) - (dlnLdlz * dlnLdlz)) / ((lh - 1.0) * (lh - 1.0)));  

  return lh;
}

static void sumCatAsc(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			int n, const int numStates)
{
  int i, k;
  double *left, *right, *sum;

  unsigned char 
    tip[32];

  ascertainmentBiasSequence(tip, numStates);

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[numStates * tip[i]]);
	  right = &(tipVector[numStates * tip[i]]);

	  
	  sum = &sumtable[i * numStates];
	  
	  for(k = 0; k < numStates; k++)
	    sum[k] = left[k] * right[k];	  
	}
      break;
    case PLL_TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[numStates * tip[i]]);

	  
	  right = &(x2[i * numStates]);
	  sum = &sumtable[i * numStates];

	  for(k = 0; k < numStates; k++)
	    sum[k] = left[k] * right[k];	 
	}
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  left  = &(x1[i * numStates]);
	  right = &(x2[i * numStates]);
	  sum   = &(sumtable[i * numStates]);

	  for(k = 0; k < numStates; k++)
	    sum[k] = left[k] * right[k];	 
	}
      break;
    default:
      assert(0);
    }
}

static void sumGammaAsc(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			int n, const int numStates)
{
  int i, l, k;
  double *left, *right, *sum;

  const int gammaStates = numStates * 4;

  unsigned char 
    tip[32];

  ascertainmentBiasSequence(tip, numStates);

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[numStates * tip[i]]);
	  right = &(tipVector[numStates * tip[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * gammaStates + l * numStates];
	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case PLL_TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[numStates * tip[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum = &sumtable[i * gammaStates + l * numStates];

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[gammaStates * i + l * numStates]);
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum   = &(sumtable[i * gammaStates + l * numStates]);

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}




#if (!defined(__AVX) && !defined(__SSE3))
static void coreCAT_FLEX(int upper, int numberOfCategories, double *sum,
    volatile double *d1, volatile double *d2, int *wgt,
    double *rptr, double *EIGN, int *cptr, double lz, const int states)
    /* rptr perSiteRates pointer, cptr rateCategory pointer */
{
  int 
    i, 
    l;

  double 
    *d, 

    /* arrays to store stuff we can pre-compute */
    *d_start = NULL,
    *e = NULL,
    *s = NULL,
    *dd = NULL,
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  rax_posix_memalign ((void **) &d_start, PLL_BYTE_ALIGNMENT, numberOfCategories * states * sizeof(double));
  rax_posix_memalign ((void **) &e,       PLL_BYTE_ALIGNMENT, (states * sizeof(double)));
  rax_posix_memalign ((void **) &s,       PLL_BYTE_ALIGNMENT, states * sizeof(double));
  rax_posix_memalign ((void **) &dd,      PLL_BYTE_ALIGNMENT, states * sizeof(double)),
  d = d_start;

  e[0] = 0.0;
  s[0] = 0.0; 
  dd[0] = 0.0;


  /* we are pre-computing values for computing the first and second derivative of P(lz)
     since this requires an exponetial that the only thing we really have to derive here */

  for(l = 1; l < states; l++)
  { 
    s[l]  = EIGN[l];
    e[l]  = EIGN[l] * EIGN[l];     
    dd[l] = s[l] * lz;
  }

  /* compute the P matrices and their derivatives for 
     all per-site rate categories */

  for(i = 0; i < numberOfCategories; i++)
  {      
    d[states * i] = 1.0;
    for(l = 1; l < states; l++)
      d[states * i + l] = exp(dd[l] * rptr[i]);
  }


  /* now loop over the sites in this partition to obtain the per-site 1st and 2nd derivatives */

  for (i = 0; i < upper; i++)
  {    
    /* get the correct p matrix for the rate at the current site i */

    d = &d_start[states * cptr[i]];      

    /* this is the likelihood at site i, NOT the log likelihood, we don't need the log 
       likelihood to compute derivatives ! */

    inv_Li     = sum[states * i]; 

    /* those are for storing the first and second derivative of the Likelihood at site i */

    dlnLidlz   = 0.0;
    d2lnLidlz2 = 0.0;

    /* now multiply the likelihood and the first and second derivative with the 
       appropriate derivatives of P(lz) */

    for(l = 1; l < states; l++)
    {
      double
        tmpv = d[l] * sum[states * i + l];

      inv_Li     += tmpv;                 
      dlnLidlz   += tmpv * s[l];       
      d2lnLidlz2 += tmpv * e[l];
    }     

    /* below we are implementing the other mathematical operations that are required 
       to obtain the deirivatives */

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    /* compute the accumulated first and second derivatives of this site */

    dlnLdlz  += wgt[i] * rptr[cptr[i]] * dlnLidlz;
    d2lnLdlz2 += wgt[i] * rptr[cptr[i]] * rptr[cptr[i]] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  /* 
     set the result values, i.e., the sum of the per-site first and second derivatives of the likelihood function 
     for this partition. 
     */

  *d1  = dlnLdlz;
  *d2 = d2lnLdlz2;

  /* free the temporary arrays */

  rax_free(d_start);
  rax_free(e);
  rax_free(s);
  rax_free(dd);
}

static void coreGAMMA_FLEX(int upper, double *sumtable, volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, 
    double *EIGN, double *gammaRates, double lz, int *wrptr, const int states)
{
  double  
    *sum, 
    diagptable[1024], /* TODO make this dynamic */
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr,
    tmp,
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2;

  int     
    i, 
    j, 
    l;  

  const int 
    gammaStates = 4 * states;

  /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

  for(i = 0; i < 4; i++)
  {
    ki = gammaRates[i];
    kisqr = ki * ki;

    for(l = 1; l < states; l++)
    {
      diagptable[i * gammaStates + l * 4]     = exp(EIGN[l] * ki * lz);
      diagptable[i * gammaStates + l * 4 + 1] = EIGN[l] * ki;
      diagptable[i * gammaStates + l * 4 + 2] = EIGN[l] * EIGN[l] * kisqr;
    }
  }

  /* loop over sites in this partition */

  for (i = 0; i < upper; i++)
  {
    /* access the array with pre-computed values */
    sum = &sumtable[i * gammaStates];

    /* initial per-site likelihood and 1st and 2nd derivatives */

    inv_Li   = 0.0;
    dlnLidlz = 0.0;
    d2lnLidlz2 = 0.0;

    /* loop over discrete GAMMA rates */

    for(j = 0; j < 4; j++)
    {
      inv_Li += sum[j * states];

      for(l = 1; l < states; l++)
      {
        inv_Li     += (tmp = diagptable[j * gammaStates + l * 4] * sum[j * states + l]);
        dlnLidlz   +=  tmp * diagptable[j * gammaStates + l * 4 + 1];
        d2lnLidlz2 +=  tmp * diagptable[j * gammaStates + l * 4 + 2];
      }
    }

    /* finalize derivative computation */
    /* note that wrptr[] here unlike in CAT above is the 
       integer weight vector of the current site 

       The operations:

       EIGN[l] * ki;
       EIGN[l] * EIGN[l] * kisqr;

       that are hidden in CAT in wrptr (at least the * ki and * ki *ki part of them 
       are done explicitely here 

*/

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz   += wrptr[i] * dlnLidlz;
    d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

}
#endif

//void sumGAMMA_FLEX_reorder(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
//    unsigned char *tipX1, unsigned char *tipX2, int n, const int states);

/** @brief Precompute values (sumtable) from the 2 likelihood vectors of a given branch
 *
 * @warning These precomputations are stored in \a tr->partitionData[model].sumBuffer, which is used by function \a execCore
 *
 * @param tr
 *   Library instance
 *
 * @warning the given branch is implicitly defined in \a tr by these nodes:
 * pNumber = tr->td[0].ti[0].pNumber;
 * qNumber = tr->td[0].ti[0].qNumber;
 *
 *
 * @note This function should be called only once at the very beginning of each Newton-Raphson procedure for optimizing barnch lengths. It initially invokes an iterative newview call to get a consistent pair of vectors at the left and the right end of the branch and thereafter invokes the one-time only precomputation of values (sumtable) that can be re-used in each Newton-Raphson iteration. Once this function has been called we can execute the actual NR procedure
 *
 *
 */
void makenewzIterative(pllInstance *tr, partitionList * pr)
{
  int 
    model, 
    tipCase;

  double
    *x1_start     = NULL,
    *x2_start     = NULL,
    *x1_start_asc = NULL,
    *x2_start_asc = NULL;


  unsigned char
    *tipX1,
    *tipX2;

  double
    *x1_gapColumn = (double*)NULL,
    *x2_gapColumn = (double*)NULL;

  unsigned int
    *x1_gap = (unsigned int*)NULL,
    *x2_gap = (unsigned int*)NULL;                            

  /* call newvieIterative to get the likelihood arrays to the left and right of the branch */

  pllNewviewIterative(tr, pr, 1);


  /* 
     loop over all partoitions to do the precomputation of the sumTable buffer 
     This is analogous to the pllNewviewIterative() and pllEvaluateIterative() 
     implementations.
     */

  for(model = 0; model < pr->numberOfPartitions; model++)
  { 
    int 
      width = pr->partitionData[model]->width;

    if(tr->td[0].executeModel[model] && width > 0)
    {
      int          
        states = pr->partitionData[model]->states;


      getVects(tr, pr, &tipX1, &tipX2, &x1_start, &x2_start, &tipCase, model, &x1_gapColumn, &x2_gapColumn, &x1_gap, &x2_gap, &x1_start_asc, &x2_start_asc);

#if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
      assert(!tr->saveMemory);
      if(tr->rateHetModel == PLL_CAT)
        sumCAT_FLEX(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
            width, states);
      else
        //sumGAMMA_FLEX_reorder(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
          sumGAMMA_FLEX(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
            width, states);
#else
      switch(states)
      {
      case 2: /* BINARY */
          assert(!tr->saveMemory);
          if (tr->rateHetModel == PLL_CAT)
            sumCAT_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                          width);
          else
            sumGAMMA_BINARY(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                            width);
          break;
      case 4: /* DNA */
#ifdef __MIC_NATIVE
      assert(!tr->saveMemory);
      assert(tr->rateHetModel == PLL_GAMMA);

      sumGTRGAMMA_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
          width);
#else
          if(tr->rateHetModel == PLL_CAT)
          {
            if(tr->saveMemory)
              sumCAT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumCAT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width);
          }
          else
          {
            if(tr->saveMemory)
              sumGAMMA_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumGAMMA(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width);
          }
#endif
          break;                
        case 20: /* proteins */
#ifdef __MIC_NATIVE
          assert(!tr->saveMemory);
          assert(tr->rateHetModel == PLL_GAMMA);

              if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                          sumGTRGAMMAPROT_LG4_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4, tipX1, tipX2,
                                  width);
              else
                          sumGTRGAMMAPROT_MIC(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                                  width);
#else

            if(tr->rateHetModel == PLL_CAT)
          {
            if(tr->saveMemory)
              sumGTRCATPROT_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else                      
              sumGTRCATPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width);
          }
          else
          {

            if(tr->saveMemory)
              sumGAMMAPROT_GAPPED_SAVE(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
              else
                    {
                      if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                        sumGAMMAPROT_LG4(tipCase,  pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector_LG4,
                                         tipX1, tipX2, width);
            else
              sumGAMMAPROT(tipCase, pr->partitionData[model]->sumBuffer, x1_start, x2_start, pr->partitionData[model]->tipVector,
                  tipX1, tipX2, width);
                    }
          }
#endif
          break;                
        default:
          assert(0);
      }
#endif
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if (pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
      if (pr->partitionData[model]->ascBias)
#endif
       {
            int pNumber = tr->td[0].ti[0].pNumber, qNumber =
                    tr->td[0].ti[0].qNumber, i, *ex1_asc =
                    &pr->partitionData[model]->ascExpVector[(pNumber
                            - tr->mxtips - 1) * states], *ex2_asc =
                    &pr->partitionData[model]->ascExpVector[(qNumber
                            - tr->mxtips - 1) * states];
            switch (tipCase)
            {
            case PLL_TIP_TIP:
                assert(0);
                break;
            case PLL_TIP_INNER:
                if (isTip(pNumber, tr->mxtips))
                {
                    for (i = 0; i < states; i++)
                        pr->partitionData[model]->ascScaler[i] = pow(
                                PLL_MINLIKELIHOOD, (double) ex2_asc[i]);
                }
                else
                {
                    for (i = 0; i < states; i++)
                        pr->partitionData[model]->ascScaler[i] = pow(
                                PLL_MINLIKELIHOOD, (double) ex1_asc[i]);
                }
                break;
            case PLL_INNER_INNER:
                for (i = 0; i < states; i++)
                    pr->partitionData[model]->ascScaler[i] = pow(
                            PLL_MINLIKELIHOOD,
                            (double) (ex1_asc[i] + ex2_asc[i]));
                break;
            default:
                assert(0);
            }
         if (tr->rateHetModel == PLL_CAT)
           sumCatAsc  (tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
         else
           sumGammaAsc(tipCase, pr->partitionData[model]->ascSumBuffer, x1_start_asc, x2_start_asc, pr->partitionData[model]->tipVector, states, states);
       }
    }
  }
}


/** @brief Compute first and second derivatives of the likelihood with respect to a given branch length 
 *
 * @param tr
 *   library instance
 *
 * @param _dlnLdlz 
 *   First derivative dl/dlz
 *
 * @param _d2lnLdlz2
 *   Second derivative d(dl/dlz)/dlz
 *
 * @warning \a makenewzIterative should have been called to precompute \a tr->partitionData[model].sumBuffer at the given branch
 *
 * @note  this function actually computes the first and second derivatives of the likelihood for a given branch stored in tr->coreLZ[model] Note that in the parallel case coreLZ must always be broadcasted together with the traversal descriptor, at least for optimizing branch lengths 
 *
 */
void execCore(pllInstance *tr, partitionList *pr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2)
{
  int model, branchIndex;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  double lz;

  _dlnLdlz[0]   = 0.0;
  _d2lnLdlz2[0] = 0.0;

  /* loop over partitions */

  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      width = pr->partitionData[model]->width;

    /* check if we (the present thread for instance) needs to compute something at 
       all for the present partition */

    if(tr->td[0].executeModel[model] && width > 0)
    {
      int           
        states = pr->partitionData[model]->states;

      double 
        *sumBuffer       = (double*)NULL;


      volatile double
        dlnLdlz   = 0.0,
                  d2lnLdlz2 = 0.0;

      /* set a pointer to the part of the pre-computed sumBuffer we are going to access */

      sumBuffer = pr->partitionData[model]->sumBuffer;

      /* figure out if we are optimizing branch lengths individually per partition or jointly across 
         all partitions. If we do this on a per partition basis, we also need to compute and store 
         the per-partition derivatives of the likelihood separately, otherwise not */

      if(numBranches > 1)
      {
        branchIndex = model;          
        lz = tr->td[0].parameterValues[model];
        _dlnLdlz[model]   = 0.0;
        _d2lnLdlz2[model] = 0.0;
      }
      else
      {
        branchIndex = 0;              
        lz = tr->td[0].parameterValues[0];
      }

#if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
      /* compute first and second derivatives with the slow generic functions */

      if(tr->rateHetModel == PLL_CAT)
        coreCAT_FLEX(width, pr->partitionData[model]->numberOfCategories, sumBuffer,
            &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->wgt,
            pr->partitionData[model]->perSiteRates, pr->partitionData[model]->EIGN,  pr->partitionData[model]->rateCategory, lz, states);
      else
        coreGAMMA_FLEX(width, sumBuffer,
            &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->EIGN, pr->partitionData[model]->gammaRates, lz,
            pr->partitionData[model]->wgt, states);
#else
      switch(states)
       {    
         case 2: /* BINARY */
           if (tr->rateHetModel == PLL_CAT)
              coreGTRCAT_BINARY(width, 
                                pr->partitionData[model]->numberOfCategories, 
                                sumBuffer,
                                &dlnLdlz, 
                                &d2lnLdlz2, 
                                pr->partitionData[model]->perSiteRates, 
                                pr->partitionData[model]->EIGN,  
                                pr->partitionData[model]->rateCategory, 
                                lz, 
                                pr->partitionData[model]->wgt);
           else
              coreGTRGAMMA_BINARY(width, 
                                   sumBuffer,
                                   &dlnLdlz, 
                                   &d2lnLdlz2, 
                                   pr->partitionData[model]->EIGN,
                                   pr->partitionData[model]->gammaRates, 
                                   lz,
                                   pr->partitionData[model]->wgt);
           break;
         case 4: /* DNA */
#ifdef __MIC_NATIVE
           assert(tr->rateHetModel == PLL_GAMMA);

           coreGTRGAMMA_MIC(width, 
                            sumBuffer,
                            &dlnLdlz, 
                            &d2lnLdlz2, 
                            pr->partitionData[model]->EIGN, 
                            pr->partitionData[model]->gammaRates, 
                            lz,
                            pr->partitionData[model]->wgt);
#else
          if(tr->rateHetModel == PLL_CAT)
            coreGTRCAT(width, pr->partitionData[model]->numberOfCategories, sumBuffer,
                &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->wgt,
                pr->partitionData[model]->perSiteRates, pr->partitionData[model]->EIGN,  pr->partitionData[model]->rateCategory, lz);
          else 
            coreGTRGAMMA(width, sumBuffer,
                &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->EIGN, pr->partitionData[model]->gammaRates, lz,
                pr->partitionData[model]->wgt);

#endif
          break;                    
        case 20: /* proteins */

#ifdef __MIC_NATIVE
      assert(tr->rateHetModel == PLL_GAMMA);

          if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                  coreGTRGAMMAPROT_LG4_MIC(width, sumBuffer,
                          &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->EIGN_LG4, pr->partitionData[model]->gammaRates, lz,
                          pr->partitionData[model]->wgt, pr->partitionData[model]->lg4x_weights);
          else
                  coreGTRGAMMAPROT_MIC(width, sumBuffer,
                          &dlnLdlz, &d2lnLdlz2, pr->partitionData[model]->EIGN, pr->partitionData[model]->gammaRates, lz,
                          pr->partitionData[model]->wgt);
#else

          if(tr->rateHetModel == PLL_CAT)
            coreGTRCATPROT(pr->partitionData[model]->EIGN, lz, pr->partitionData[model]->numberOfCategories,  pr->partitionData[model]->perSiteRates,
                pr->partitionData[model]->rateCategory, width,
                pr->partitionData[model]->wgt,
                &dlnLdlz, &d2lnLdlz2,
                sumBuffer);
            else
                { 
                  if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                    coreGTRGAMMAPROT_LG4(pr->partitionData[model]->gammaRates, pr->partitionData[model]->EIGN_LG4,
                                         sumBuffer, width, pr->partitionData[model]->wgt,
                                         &dlnLdlz, &d2lnLdlz2, lz, pr->partitionData[model]->lg4x_weights);
          else

            coreGTRGAMMAPROT(pr->partitionData[model]->gammaRates, pr->partitionData[model]->EIGN,
                sumBuffer, width, pr->partitionData[model]->wgt,
                &dlnLdlz, &d2lnLdlz2, lz);
            
                }
#endif
          break;                   
        default:
          assert(0);
      }
#endif

      /* store first and second derivative */
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
     if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
     if(pr->partitionData[model]->ascBias)
#endif  
       {
         double correction = 0;
         int             w = 0;
         
         volatile double 
           d1 = 0.0,
           d2 = 0.0;                   
         
         for(size_t i = (size_t)pr->partitionData[model]->lower; i < (size_t)pr->partitionData[model]->upper; i++)
           w += tr->aliaswgt[i];     
         
          switch(tr->rateHetModel)
            {
            case PLL_CAT:
              correction = coreCatAsc(pr->partitionData[model]->EIGN, pr->partitionData[model]->ascSumBuffer, states,
                                        &d1,  &d2, lz, states, pr->partitionData[model]->ascScaler);
              break;
            case PLL_GAMMA:
              correction = coreGammaAsc(pr->partitionData[model]->gammaRates, pr->partitionData[model]->EIGN, pr->partitionData[model]->ascSumBuffer, states,
                                        &d1,  &d2, lz, states, pr->partitionData[model]->ascScaler);
              break;
            default:
              assert(0);
            }        
         correction = 1.0 - correction; //Never used!
     
         /* Lewis correction */
         _dlnLdlz[branchIndex]   =  _dlnLdlz[branchIndex] + dlnLdlz - (double)w * d1;
         _d2lnLdlz2[branchIndex] =  _d2lnLdlz2[branchIndex] + d2lnLdlz2-  (double)w * d2;
       }  
      else
       {
         _dlnLdlz[branchIndex]   = _dlnLdlz[branchIndex]   + dlnLdlz;
         _d2lnLdlz2[branchIndex] = _d2lnLdlz2[branchIndex] + d2lnLdlz2;
       }
    }
    else if(width == 0 && (numBranches > 1)) {
          /* set to 0 to make the reduction operation consistent */
          _dlnLdlz[model]   = 0.0;
        _d2lnLdlz2[model] = 0.0;
      }                                    
  }
}


/* the function below actually implements the iterative Newton-Raphson procedure.
   It is particularly messy and hard to read because for the case of per-partition branch length 
   estimates it needs to keep track of whetehr the Newton Raphson procedure has 
   converged for each partition individually. 

   The rational efor doing it like this is also provided in:


   A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009,

*/

static void topLevelMakenewz(pllInstance *tr, partitionList * pr, double *z0, int _maxiter, double *result)
{
  double   z[PLL_NUM_BRANCHES], zprev[PLL_NUM_BRANCHES], zstep[PLL_NUM_BRANCHES];
  volatile double  dlnLdlz[PLL_NUM_BRANCHES], d2lnLdlz2[PLL_NUM_BRANCHES];
  int i, maxiter[PLL_NUM_BRANCHES], model;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;
  pllBoolean firstIteration = PLL_TRUE;
  pllBoolean outerConverged[PLL_NUM_BRANCHES];
  pllBoolean loopConverged;


  /* figure out if this is on a per partition basis or jointly across all partitions */



  /* initialize loop convergence variables etc. 
     maxiter is the maximum number of NR iterations we are going to do before giving up */

  for(i = 0; i < numBranches; i++)
  {
    z[i] = z0[i];
    maxiter[i] = _maxiter;
    outerConverged[i] = PLL_FALSE;
    tr->curvatOK[i]       = PLL_TRUE;
  }


  /* nested do while loops of Newton-Raphson */

  do
  {

    /* check if we ar done for partition i or if we need to adapt the branch length again */

    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_TRUE)
      {
        tr->curvatOK[i] = PLL_FALSE;

        zprev[i] = z[i];

        zstep[i] = (1.0 - PLL_ZMAX) * z[i] + PLL_ZMIN;
      }
    }

    for(i = 0; i < numBranches; i++)
    {
      /* other case, the outer loop hasn't converged but we are trying to approach 
         the maximum from the wrong side */

      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
      {
        double lz;

        if (z[i] < PLL_ZMIN) z[i] = PLL_ZMIN;
        else if (z[i] > PLL_ZMAX) z[i] = PLL_ZMAX;
        lz    = log(z[i]);

        tr->coreLZ[i] = lz;
      }
    }


    /* set the execution mask */

    if(numBranches > 1)
    {
      for(model = 0; model < pr->numberOfPartitions; model++)
      {
        if(pr->partitionData[model]->executeModel)
          pr->partitionData[model]->executeModel = !tr->curvatOK[model];

      }
    }
    else
    {
      for(model = 0; model < pr->numberOfPartitions; model++)
        pr->partitionData[model]->executeModel = !tr->curvatOK[0];
    }


    /* store it in traversal descriptor */

    storeExecuteMaskInTraversalDescriptor(tr, pr);

    /* store the new branch length values to be tested in traversal descriptor */

    storeValuesInTraversalDescriptor(tr, pr, &(tr->coreLZ[0]));

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

    /* if this is the first iteration of NR we will need to first do this one-time call 
       of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
       subsequent calls to pllMasterBarrier(PLL_THREAD_MAKENEWZ, tr); will not require this
       */

    if(firstIteration)
      {
        tr->td[0].traversalHasChanged = PLL_TRUE; 
        pllMasterBarrier (tr, pr, PLL_THREAD_MAKENEWZ_FIRST);
        firstIteration = PLL_FALSE; 
        tr->td[0].traversalHasChanged = PLL_FALSE; 
      }
    else 
      pllMasterBarrier(tr, pr, PLL_THREAD_MAKENEWZ);
    branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2, numBranches);
#else 
    /* sequential part, if this is the first newton-raphson implementation,
       do the precomputations as well, otherwise just execute the computation
       of the derivatives */
    if(firstIteration)
      {
        makenewzIterative(tr, pr);
        firstIteration = PLL_FALSE;
      }
    execCore(tr, pr, dlnLdlz, d2lnLdlz2);
#endif

    /* do a NR step, if we are on the correct side of the maximum that's okay, otherwise 
       shorten branch */

    for(i = 0; i < numBranches; i++)
    {
      if(outerConverged[i] == PLL_FALSE && tr->curvatOK[i] == PLL_FALSE)
      {
        if ((d2lnLdlz2[i] >= 0.0) && (z[i] < PLL_ZMAX))
          zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
        else
          tr->curvatOK[i] = PLL_TRUE;
      }
    }

    /* do the standard NR step to obrain the next value, depending on the state for eahc partition */

    for(i = 0; i < numBranches; i++)
    {
      if(tr->curvatOK[i] == PLL_TRUE && outerConverged[i] == PLL_FALSE)
      {
        if (d2lnLdlz2[i] < 0.0)
        {
          double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
          if (tantmp < 100)
          {
            z[i] *= exp(tantmp);
            if (z[i] < PLL_ZMIN)
              z[i] = PLL_ZMIN;

            if (z[i] > 0.25 * zprev[i] + 0.75)
              z[i] = 0.25 * zprev[i] + 0.75;
          }
          else
            z[i] = 0.25 * zprev[i] + 0.75;
        }
        if (z[i] > PLL_ZMAX) z[i] = PLL_ZMAX;

        /* decrement the maximum number of itarations */

        maxiter[i] = maxiter[i] - 1;

        /* check if the outer loop has converged */

        //old code below commented out, integrated new PRELIMINARY BUG FIX !
        //this needs further work at some point!

        /*
        if(maxiter[i] > 0 && (PLL_ABS(z[i] - zprev[i]) > zstep[i]))
          outerConverged[i] = PLL_FALSE;
        else
          outerConverged[i] = PLL_TRUE;
        */

        if((PLL_ABS(z[i] - zprev[i]) > zstep[i]))
         {
           /* We should make a more informed decision here,
              based on the log like improvement */

           if(maxiter[i] < -20)
            {
              z[i] = z0[i];
              outerConverged[i] = PLL_TRUE;
            }
           else
             outerConverged[i] = PLL_FALSE;
         }
        else
          outerConverged[i] = PLL_TRUE;
      }
    }

    /* check if the loop has converged for all partitions */

    loopConverged = PLL_TRUE;
    for(i = 0; i < numBranches; i++)
      loopConverged = loopConverged && outerConverged[i];
  }
  while (!loopConverged);


  /* reset  partition execution mask */

  for(model = 0; model < pr->numberOfPartitions; model++)
    pr->partitionData[model]->executeModel = PLL_TRUE;

  /* copy the new branches in the result array of branches.
     if we don't do a per partition estimate of 
     branches this will only set result[0]
     */

  for(i = 0; i < numBranches; i++)
    result[i] = z[i];
}


/** @brief Optimize branch length value(s) of a given branch with the Newton-Raphtson procedure 
 *
 * @warning A given branch may have one or several branch length values (up to PLL_NUM_BRANCHES), usually the later refers to partition-specific branch length values. Thus z0 and result represent collections rather than double values. The number of branch length values is given by \a tr->numBranches 
 *
 * @param tr
 *   Library instance
 *
 * @param p
 *   One node that defines the branch (p->z)
 *
 * @param q
 *   The other node side of the branch (usually p->back), but the branch length can be estimated even if p and q are
 *   not connected, e.g. before the insertion of a subtree.
 *
 * @param z0 
 *   Initial branch length value(s) for the given branch \a p->z 
 *
 * @param maxiter 
 *   Maximum number of iterations in the Newton-Raphson procedure 
 *
 * @param result 
 *   Resulting branch length value(s) for the given branch \a p->z 
 *
 * @param mask 
 *   Specifies if a mask to track partition convergence (\a tr->partitionConverged) is being used.
 *
 * @sa typical values for \a maxiter are constants \a iterations and \a PLL_NEWZPERCYCLE
 * @note Requirement: q->z == p->z
 */
void makenewzGeneric(pllInstance *tr, partitionList * pr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, pllBoolean mask)
{
  int i;
  //boolean originalExecute[PLL_NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  pllBoolean 
    p_recom = PLL_FALSE, /* if one of was missing, we will need to force recomputation */
    q_recom = PLL_FALSE;

  /* the first entry of the traversal descriptor stores the node pair that defines 
     the branch */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < numBranches; i++)
  {
    //originalExecute[i] =  pr->partitionData[i]->executeModel;
    tr->td[0].ti[0].qz[i] =  z0[i];
    if(mask)
    {
      if (tr->partitionConverged[i])
        pr->partitionData[i]->executeModel = PLL_FALSE;
      else
        pr->partitionData[i]->executeModel = PLL_TRUE;
    }
  }
  if (tr->useRecom)
  {
    int
      slot = -1;
      //count = 0;

    /* Ensure p and q get a unpinnable slot in physical memory */
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
  }


  /* compute the traversal descriptor of the likelihood vectors  that need to be re-computed 
     first in makenewzIterative */

  tr->td[0].count = 1;

  if(p_recom || needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips))
    computeTraversal(tr, p, PLL_TRUE, numBranches);

  if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
    computeTraversal(tr, q, PLL_TRUE, numBranches);

  /* call the Newton-Raphson procedure */

  topLevelMakenewz(tr, pr, z0, maxiter, result);

  /* Mark node as unpinnable */
  if(tr->useRecom)
  {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
  }

  /* fix eceuteModel this seems to be a bit redundant with topLevelMakenewz */ 

  for(i = 0; i < numBranches; i++)
    pr->partitionData[i]->executeModel = PLL_TRUE;
}


/* below are, once again the optimized functions */

#if (defined(__SSE3) || defined(__AVX))


static void sumCAT_BINARY(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
                          unsigned char *tipX1, unsigned char *tipX2, int n)

{
  int i;
  
#if (!defined(__SSE3) && !defined(__AVX))
  int j;
#endif
  double *x1, *x2;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
        {
          x1 = &(tipVector[2 * tipX1[i]]);
          x2 = &(tipVector[2 * tipX2[i]]);

#if (!defined(__SSE3) && !defined(__AVX))
          for(j = 0; j < 2; j++)
            sum[i * 2 + j]     = x1[j] * x2[j];
#else
          _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));
#endif
        }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &(tipVector[2 * tipX1[i]]);
          x2 = &x2_start[2 * i];

#if (!defined(__SSE3) && !defined(__AVX))
          for(j = 0; j < 2; j++)
            sum[i * 2 + j]     = x1[j] * x2[j];
#else
          _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));  
#endif
        }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
#if (!defined(__SSE3) && !defined(__AVX))
          for(j = 0; j < 2; j++)
            sum[i * 2 + j]     = x1[j] * x2[j];
#else
          _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));   
#endif
        }
      break;
    default:
      assert(0);
    }
}


static void sumCAT_SAVE(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int i;
  double 
    *x1, 
    *x2,    
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        if(isGap(x2_gap, i))
          x2 = x2_gapColumn;
        else
        {
          x2 = x2_ptr;
          x2_ptr += 4;
        }

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
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

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));

      }    
      break;
    default:
      assert(0);
  }
}

static void sumGAMMA_BINARY(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
                            unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j;
#if (!defined(_USE_PTHREADS) && !defined(_FINE_GRAIN_MPI))
  int k;
#endif

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++)
        {
          x1 = &(tipVector[2 * tipX1[i]]);
          x2 = &(tipVector[2 * tipX2[i]]);
          sum = &sumtable[i * 8];
#if (!defined(_USE_PTHREADS) && !defined(_FINE_GRAIN_MPI))
          for(j = 0; j < 4; j++)
            for(k = 0; k < 2; k++)
              sum[j * 2 + k] = x1[k] * x2[k];
#else
          for(j = 0; j < 4; j++)
            _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));         
#endif
        }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
        {
          x1  = &(tipVector[2 * tipX1[i]]);
          x2  = &x2_start[8 * i];
          sum = &sumtable[8 * i];

#if (!defined(_USE_PTHREADS) && !defined(_FINE_GRAIN_MPI))
          for(j = 0; j < 4; j++)
            for(k = 0; k < 2; k++)
              sum[j * 2 + k] = x1[k] * x2[j * 2 + k];
#else
          for(j = 0; j < 4; j++)
            _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[j * 2] )));
#endif
        }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1  = &x1_start[8 * i];
          x2  = &x2_start[8 * i];
          sum = &sumtable[8 * i];
#if (!defined(_USE_PTHREADS) && !defined(_FINE_GRAIN_MPI))
          for(j = 0; j < 4; j++)
            for(k = 0; k < 2; k++)
              sum[j * 2 + k] = x1[j * 2 + k] * x2[j * 2 + k];
#else
          for(j = 0; j < 4; j++)
            _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[j * 2] ), _mm_load_pd( &x2[j * 2] )));
#endif
        }
      break;
    default:
      assert(0);
    }
}


static void sumGAMMA_GAPPED_SAVE(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double 
    *x1, 
    *x2, 
    *sum,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

  int i, j, k; 

  switch(tipCase)
  {
    case PLL_TIP_TIP:     
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);
        sum = &sumtable[i * 16];

        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[k] )));
      }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1  = &(tipVector[4 * tipX1[i]]);

        if(x2_gap[i / 32] & mask32[i % 32])
          x2 = x2_gapColumn;
        else
        {
          x2  = x2_ptr;
          x2_ptr += 16;
        }

        sum = &sumtable[16 * i];

        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[j * 4 + k] )));
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        if(x1_gap[i / 32] & mask32[i % 32])
          x1 = x1_gapColumn;
        else
        {
          x1  = x1_ptr;
          x1_ptr += 16;
        }

        if(x2_gap[i / 32] & mask32[i % 32])
          x2 = x2_gapColumn;
        else
        {
          x2  = x2_ptr;
          x2_ptr += 16;
        }

        sum = &sumtable[16 * i];


        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[j * 4 + k] ), _mm_load_pd( &x2[j * 4 + k] )));
      }
      break;
    default:
      assert(0);
  }
}




static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j, k;

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);
        sum = &sumtable[i * 16];

        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[k] )));
      }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1  = &(tipVector[4 * tipX1[i]]);
        x2  = &x2_start[16 * i];
        sum = &sumtable[16 * i];

        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[j * 4 + k] )));
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        x1  = &x1_start[16 * i];
        x2  = &x2_start[16 * i];
        sum = &sumtable[16 * i];

        for(j = 0; j < 4; j++)      
          for(k = 0; k < 4; k+=2)
            _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[j * 4 + k] ), _mm_load_pd( &x2[j * 4 + k] )));
      }
      break;
    default:
      assert(0);
  }
}


static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i;
  double 
    *x1, 
    *x2;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &x2_start[4 * i];

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &x1_start[4 * i];
        x2 = &x2_start[4 * i];

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));

      }    
      break;
    default:
      assert(0);
  }
}
static void sumGAMMAPROT_GAPPED_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int i, l, k;
  double 
    *left, 
    *right, 
    *sum,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x1v,
    *x2v;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);

        for(l = 0; l < 4; l++)
        {
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }

        }
      }
      break;
    case PLL_TIP_INNER:
      for(i = 0; i < n; i++)
      {
        left = &(tipVector[20 * tipX1[i]]);

        if(x2_gap[i / 32] & mask32[i % 32])
          x2v = x2_gapColumn;
        else
        {
          x2v = x2_ptr;
          x2_ptr += 80;
        }

        for(l = 0; l < 4; l++)
        {
          right = &(x2v[l * 20]);
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }
        }
      }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
      {
        if(x1_gap[i / 32] & mask32[i % 32])
          x1v = x1_gapColumn;
        else
        {
          x1v  = x1_ptr;
          x1_ptr += 80;
        }

        if(x2_gap[i / 32] & mask32[i % 32])
          x2v = x2_gapColumn;
        else
        {
          x2v  = x2_ptr;
          x2_ptr += 80;
        }

        for(l = 0; l < 4; l++)
        {
          left  = &(x1v[l * 20]);
          right = &(x2v[l * 20]);
          sum   = &(sumtable[i * 80 + l * 20]);

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }
        }
      }
      break;
    default:
      assert(0);
  }
}


static void sumGAMMAPROT_LG4(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector[4],
                             unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
        {         
          for(l = 0; l < 4; l++)
            {
              left  = &(tipVector[l][20 * tipX1[i]]);
              right = &(tipVector[l][20 * tipX2[i]]);

              sum = &sumtable[i * 80 + l * 20];
#ifdef __SSE3
              for(k = 0; k < 20; k+=2)
                {
                  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
                  
                  _mm_store_pd(&sum[k], sumv);           
                }
#else
              for(k = 0; k < 20; k++)
                sum[k] = left[k] * right[k];
#endif
            }
        }
      break;
    case PLL_TIP_INNER:
      for(i = 0; i < n; i++)
        {
         

          for(l = 0; l < 4; l++)
            { 
              left = &(tipVector[l][20 * tipX1[i]]);
              right = &(x2[80 * i + l * 20]);
              sum = &sumtable[i * 80 + l * 20];
#ifdef __SSE3
              for(k = 0; k < 20; k+=2)
                {
                  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
                  
                  _mm_store_pd(&sum[k], sumv);           
                }
#else
              for(k = 0; k < 20; k++)
                sum[k] = left[k] * right[k];
#endif
            }
        }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
        {
          for(l = 0; l < 4; l++)
            {
              left  = &(x1[80 * i + l * 20]);
              right = &(x2[80 * i + l * 20]);
              sum   = &(sumtable[i * 80 + l * 20]);

#ifdef __SSE3
              for(k = 0; k < 20; k+=2)
                {
                  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
                  
                  _mm_store_pd(&sum[k], sumv);           
                }
#else
              for(k = 0; k < 20; k++)
                sum[k] = left[k] * right[k];
#endif
            }
        }
      break;
    default:
      assert(0);
    }
}


static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for(i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);

        for(l = 0; l < 4; l++)
        {
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }

        }
      }
      break;
    case PLL_TIP_INNER:
      for(i = 0; i < n; i++)
      {
        left = &(tipVector[20 * tipX1[i]]);

        for(l = 0; l < 4; l++)
        {
          right = &(x2[80 * i + l * 20]);
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }

        }
      }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
      {
        for(l = 0; l < 4; l++)
        {
          left  = &(x1[80 * i + l * 20]);
          right = &(x2[80 * i + l * 20]);
          sum   = &(sumtable[i * 80 + l * 20]);


          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);                 
          }
        }
      }
      break;
    default:
      assert(0);
  }
}


static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);
        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }

      }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
      {
        left = &(tipVector[20 * tipX1[i]]);
        right = &x2[20 * i];
        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }

      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        left  = &x1[20 * i];
        right = &x2[20 * i];
        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }

      }
      break;
    default:
      assert(0);
  }
}


static void sumGTRCATPROT_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int 
    i, 
    l;

  double 
    *sum, 
    *left, 
    *right,
    *left_ptr = x1,
    *right_ptr = x2;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);
        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }

      }
      break;
    case PLL_TIP_INNER:
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

        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }

      }
      break;
    case PLL_INNER_INNER:
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

        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);           
        }
      }
      break;
    default:
      assert(0);
  }
}

static void coreGTRGAMMA(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  double 
    dlnLdlz = 0.0,
            d2lnLdlz2 = 0.0,
            ki, 
            kisqr,  
            inv_Li, 
            dlnLidlz, 
            d2lnLidlz2,  
		*sum;
	PLL_ALIGN_BEGIN double
            diagptable0[16] PLL_ALIGN_END,
            diagptable1[16] PLL_ALIGN_END,
            diagptable2[16] PLL_ALIGN_END;

  int     
    i, 
    j, 
    l;

  for(i = 0; i < 4; i++)
  {
    ki = gammaRates[i];
    kisqr = ki * ki;

    diagptable0[i * 4] = 1.0;
    diagptable1[i * 4] = 0.0;
    diagptable2[i * 4] = 0.0;

    for(l = 1; l < 4; l++)
    {
      diagptable0[i * 4 + l] = exp(EIGN[l] * ki * lz);
      diagptable1[i * 4 + l] = EIGN[l] * ki;
      diagptable2[i * 4 + l] = EIGN[l] * EIGN[l] * kisqr;
    }
  }

  for (i = 0; i < upper; i++)
  { 
    __m128d a0 = _mm_setzero_pd();
    __m128d a1 = _mm_setzero_pd();
    __m128d a2 = _mm_setzero_pd();

    sum = &sumtable[i * 16];         

    for(j = 0; j < 4; j++)
    {                   
      double       
        *d0 = &diagptable0[j * 4],
        *d1 = &diagptable1[j * 4],
        *d2 = &diagptable2[j * 4];

      for(l = 0; l < 4; l+=2)
      {
        __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 4 + l]));
        a0 = _mm_add_pd(a0, tmpv);
        a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
        a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
      }           
    }

    a0 = _mm_hadd_pd(a0, a0);
    a1 = _mm_hadd_pd(a1, a1);
    a2 = _mm_hadd_pd(a2, a2);

    _mm_storel_pd(&inv_Li, a0);     
    _mm_storel_pd(&dlnLidlz, a1);
    _mm_storel_pd(&d2lnLidlz2, a2); 

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;     

    dlnLdlz   += wrptr[i] * dlnLidlz;
    d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }


  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}

static void coreGTRCAT_BINARY(int upper, int numberOfCategories, double *sum,
                              volatile double *d1, volatile double *d2, 
                              double *rptr, double *EIGN, int *cptr, double lz, int *wgt)
{
  int i;
  double
    *d, *d_start = NULL,
    tmp_0, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[2];
  double dd1;

  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];


  d = d_start = (double *)rax_malloc(numberOfCategories * sizeof(double));

  dd1 = e[0] * lz;

  for(i = 0; i < numberOfCategories; i++)
    d[i] = exp(dd1 * rptr[i]);

  for (i = 0; i < upper; i++)
    {
      double
        r = rptr[cptr[i]],
        wr1 = r * wgt[i],
        wr2 = r * r * wgt[i];
      
      d = &d_start[cptr[i]];

      inv_Li = sum[2 * i];
      inv_Li += (tmp_0 = d[0] * sum[2 * i + 1]);

      inv_Li = 1.0/fabs(inv_Li);

      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(d_start);
}


static void coreGTRCAT(int upper, int numberOfCategories, double *sum,
    volatile double *d1, volatile double *d2, int *wgt,
    double *rptr, double *EIGN, int *cptr, double lz)
{
  int i;
  double
    *d, *d_start = NULL,
    inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  PLL_ALIGN_BEGIN double e1[4] PLL_ALIGN_END;
  PLL_ALIGN_BEGIN double e2[4] PLL_ALIGN_END;
  double dd1, dd2, dd3;

  __m128d
    e1v[2],
    e2v[2];

  e1[0] = 0.0;
  e2[0] = 0.0;
  e1[1] = EIGN[1];
  e2[1] = EIGN[1] * EIGN[1];
  e1[2] = EIGN[2];
  e2[2] = EIGN[2] * EIGN[2];
  e1[3] = EIGN[3];
  e2[3] = EIGN[3] * EIGN[3];

  e1v[0]= _mm_load_pd(&e1[0]);
  e1v[1]= _mm_load_pd(&e1[2]);

  e2v[0]= _mm_load_pd(&e2[0]);
  e2v[1]= _mm_load_pd(&e2[2]);

  rax_posix_memalign ((void **) &d_start, PLL_BYTE_ALIGNMENT, numberOfCategories * 4 * sizeof(double));
  d = d_start;

  dd1 = EIGN[1] * lz;
  dd2 = EIGN[2] * lz;
  dd3 = EIGN[3] * lz;

  for(i = 0; i < numberOfCategories; i++)
  {
    d[i * 4 + 0] = 1.0;
    d[i * 4 + 1] = exp(dd1 * rptr[i]);
    d[i * 4 + 2] = exp(dd2 * rptr[i]);
    d[i * 4 + 3] = exp(dd3 * rptr[i]);
  }

  for (i = 0; i < upper; i++)
  {
    double *s = &sum[4 * i];
    d = &d_start[4 * cptr[i]];  

    __m128d tmp_0v =_mm_mul_pd(_mm_load_pd(&d[0]),_mm_load_pd(&s[0]));
    __m128d tmp_1v =_mm_mul_pd(_mm_load_pd(&d[2]),_mm_load_pd(&s[2]));

    __m128d inv_Liv    = _mm_add_pd(tmp_0v, tmp_1v);      

    __m128d dlnLidlzv   = _mm_add_pd(_mm_mul_pd(tmp_0v, e1v[0]), _mm_mul_pd(tmp_1v, e1v[1]));     
    __m128d d2lnLidlz2v = _mm_add_pd(_mm_mul_pd(tmp_0v, e2v[0]), _mm_mul_pd(tmp_1v, e2v[1]));


    inv_Liv   = _mm_hadd_pd(inv_Liv, inv_Liv);
    dlnLidlzv = _mm_hadd_pd(dlnLidlzv, dlnLidlzv);
    d2lnLidlz2v = _mm_hadd_pd(d2lnLidlz2v, d2lnLidlz2v);                 

    _mm_storel_pd(&inv_Li, inv_Liv);     
    _mm_storel_pd(&dlnLidlz, dlnLidlzv);                 
    _mm_storel_pd(&d2lnLidlz2, d2lnLidlz2v);      

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz  += wgt[i] * rptr[cptr[i]] * dlnLidlz;
    d2lnLdlz2 += wgt[i] * rptr[cptr[i]] * rptr[cptr[i]] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(d_start);
}

#if (!defined(__SSE3) && !defined(__AVX))
static void coreGTRGAMMA_BINARY(const int upper, double *sumtable,
                                volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i, j;
  double
    *diagptable, *diagptable_start, *sum,
    tmp_1, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  diagptable = diagptable_start = (double *)rax_malloc(sizeof(double) * 12);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      diagptable[i * 3]     = exp (EIGN[1] * ki * lz);
      diagptable[i * 3 + 1] = EIGN[1] * ki;
      diagptable[i * 3 + 2] = EIGN[1] * EIGN[1] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      diagptable = diagptable_start;
      sum = &(sumtable[i * 8]);

      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;

      for(j = 0; j < 4; j++)
        {
          inv_Li += sum[2 * j];

          tmp_1      =  diagptable[3 * j] * sum[2 * j + 1];
          inv_Li     += tmp_1;
          dlnLidlz   += tmp_1 * diagptable[3 * j + 1];
          d2lnLidlz2 += tmp_1 * diagptable[3 * j + 2];
        }

      inv_Li = 1.0 / fabs(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;


      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(diagptable_start);
}
#else
static void coreGTRGAMMA_BINARY(const int upper, double *sumtable,
                                volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
	double
		dlnLdlz = 0.0,
		d2lnLdlz2 = 0.0,
		ki,
		kisqr,
		inv_Li,
		dlnLidlz,
		d2lnLidlz2,
		*sum;
	PLL_ALIGN_BEGIN double
		diagptable0[8] PLL_ALIGN_END,
		diagptable1[8] PLL_ALIGN_END,
		diagptable2[8] PLL_ALIGN_END;
    
  int     
    i, 
    j;
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 2] = 1.0;
      diagptable1[i * 2] = 0.0;
      diagptable2[i * 2] = 0.0;
     
      diagptable0[i * 2 + 1] = exp(EIGN[0] * ki * lz);
      diagptable1[i * 2 + 1] = EIGN[0] * ki;
      diagptable2[i * 2 + 1] = EIGN[0] * EIGN[0] * kisqr;    
    }

  for (i = 0; i < upper; i++)
    { 
      __m128d a0 = _mm_setzero_pd();
      __m128d a1 = _mm_setzero_pd();
      __m128d a2 = _mm_setzero_pd();

      sum = &sumtable[i * 8];         

      for(j = 0; j < 4; j++)
        {                       
          double           
            *d0 = &diagptable0[j * 2],
            *d1 = &diagptable1[j * 2],
            *d2 = &diagptable2[j * 2];
                         
          __m128d tmpv = _mm_mul_pd(_mm_load_pd(d0), _mm_load_pd(&sum[j * 2]));
          a0 = _mm_add_pd(a0, tmpv);
          a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(d1)));
          a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(d2)));
                          
        }

      a0 = _mm_hadd_pd(a0, a0);
      a1 = _mm_hadd_pd(a1, a1);
      a2 = _mm_hadd_pd(a2, a2);

      _mm_storel_pd(&inv_Li, a0);     
      _mm_storel_pd(&dlnLidlz, a1);
      _mm_storel_pd(&d2lnLidlz2, a2); 

      inv_Li = 1.0 / fabs(inv_Li);
     
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;     

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

 
  *d1   = dlnLdlz;
  *d2 = d2lnLdlz2; 
}


#endif

static void coreGTRGAMMAPROT_LG4(double *gammaRates, double *EIGN[4], double *sumtable, int upper, int *wrptr,
                                 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz,
                                 double * lg4_weights)
{
	double  *sum;
	PLL_ALIGN_BEGIN double
    diagptable0[80] PLL_ALIGN_END,
    diagptable1[80] PLL_ALIGN_END,
    diagptable2[80] PLL_ALIGN_END;    
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 20] = 1.0;
      diagptable1[i * 20] = 0.0;
      diagptable2[i * 20] = 0.0;

      for(l = 1; l < 20; l++)
        {
          diagptable0[i * 20 + l] = exp(EIGN[i][l] * ki * lz);
          diagptable1[i * 20 + l] = EIGN[i][l] * ki;
          diagptable2[i * 20 + l] = EIGN[i][l] * EIGN[i][l] * kisqr;
        }
    }

  for (i = 0; i < upper; i++)
    { 

      double
      	  inv_Li = 0.0,
      	  dlnLidlz = 0.0,
      	  d2lnLidlz2 = 0.0;

      sum = &sumtable[i * 80];         

      for(j = 0; j < 4; j++)
        {                       
          double
          	l0,
          	l1,
          	l2,
            *d0 = &diagptable0[j * 20],
            *d1 = &diagptable1[j * 20],
            *d2 = &diagptable2[j * 20];
                 
          __m128d a0 = _mm_setzero_pd();
          __m128d a1 = _mm_setzero_pd();
          __m128d a2 = _mm_setzero_pd();

          for(l = 0; l < 20; l+=2)
            {
              __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 20 +l]));
              a0 = _mm_add_pd(a0, tmpv);
              a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
              a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
            }             

          a0 = _mm_hadd_pd(a0, a0);
      	  a1 = _mm_hadd_pd(a1, a1);
      	  a2 = _mm_hadd_pd(a2, a2);

      	 _mm_storel_pd(&l0, a0);
      	 _mm_storel_pd(&l1, a1);
      	 _mm_storel_pd(&l2, a2);

      	 inv_Li     += lg4_weights[j] * l0;
      	 dlnLidlz   += lg4_weights[j] * l1;
     	 d2lnLidlz2 += lg4_weights[j] * l2;
      }

      inv_Li = 1.0 / fabs (inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}



static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
	double  *sum;
	PLL_ALIGN_BEGIN double
		diagptable0[80] PLL_ALIGN_END,
		diagptable1[80] PLL_ALIGN_END,
		diagptable2[80] PLL_ALIGN_END;

  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
  {
    ki = gammaRates[i];
    kisqr = ki * ki;

    diagptable0[i * 20] = 1.0;
    diagptable1[i * 20] = 0.0;
    diagptable2[i * 20] = 0.0;

    for(l = 1; l < 20; l++)
    {
      diagptable0[i * 20 + l] = exp(EIGN[l] * ki * lz);
      diagptable1[i * 20 + l] = EIGN[l] * ki;
      diagptable2[i * 20 + l] = EIGN[l] * EIGN[l] * kisqr;
    }
  }

  for (i = 0; i < upper; i++)
  { 
    __m128d a0 = _mm_setzero_pd();
    __m128d a1 = _mm_setzero_pd();
    __m128d a2 = _mm_setzero_pd();

    sum = &sumtable[i * 80];         

    for(j = 0; j < 4; j++)
    {                   
      double       
        *d0 = &diagptable0[j * 20],
        *d1 = &diagptable1[j * 20],
        *d2 = &diagptable2[j * 20];

      for(l = 0; l < 20; l+=2)
      {
        __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 20 +l]));
        a0 = _mm_add_pd(a0, tmpv);
        a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
        a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
      }           
    }

    a0 = _mm_hadd_pd(a0, a0);
    a1 = _mm_hadd_pd(a1, a1);
    a2 = _mm_hadd_pd(a2, a2);

    _mm_storel_pd(&inv_Li, a0);
    _mm_storel_pd(&dlnLidlz, a1);
    _mm_storel_pd(&d2lnLidlz2, a2);

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz   += wrptr[i] * dlnLidlz;
    d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}



static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
    int *wgt, volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start = NULL, *sum;
  PLL_ALIGN_BEGIN double 
    e[20] PLL_ALIGN_END, 
    s[20] PLL_ALIGN_END, 
    dd[20] PLL_ALIGN_END;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  rax_posix_memalign ((void **)&d_start, PLL_BYTE_ALIGNMENT, numberOfCategories * 20 * sizeof(double));
  d1 = d_start; 

  e[0] = 0.0;
  s[0] = 0.0; 

  for(l = 1; l < 20; l++)
  {
    e[l]  = EIGN[l] * EIGN[l];
    s[l]  = EIGN[l];
    dd[l] = s[l] * lz;
  }

  for(i = 0; i < numberOfCategories; i++)
  {      
    d1[20 * i] = 1.0;
    for(l = 1; l < 20; l++)
      d1[20 * i + l] = exp(dd[l] * rptr[i]);
  }

  for (i = 0; i < upper; i++)
  {
    __m128d a0 = _mm_setzero_pd();
    __m128d a1 = _mm_setzero_pd();
    __m128d a2 = _mm_setzero_pd();

    d1 = &d_start[20 * cptr[i]];
    sum = &sumtable[20 * i];

    for(l = 0; l < 20; l+=2)
    {     
      __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d1[l]), _mm_load_pd(&sum[l]));

      a0 = _mm_add_pd(a0, tmpv);
      __m128d sv = _mm_load_pd(&s[l]);    

      a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, sv));
      __m128d ev = _mm_load_pd(&e[l]);    

      a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, ev));
    }

    a0 = _mm_hadd_pd(a0, a0);
    a1 = _mm_hadd_pd(a1, a1);
    a2 = _mm_hadd_pd(a2, a2);

    _mm_storel_pd(&inv_Li, a0);     
    _mm_storel_pd(&dlnLidlz, a1);                 
    _mm_storel_pd(&d2lnLidlz2, a2);

    inv_Li = 1.0 / fabs (inv_Li);

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz  += wgt[i] * rptr[cptr[i]] * dlnLidlz;
    d2lnLdlz2 += wgt[i] * rptr[cptr[i]] * rptr[cptr[i]] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}




#endif



