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
 * @file newviewGenericSpecial.c
 *  
 * @brief Functions that deal (mostly) with conditional likelihood (re)computation
 */

#include "mem_alloc.h"
#include "systypes.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

#ifdef __MIC_NATIVE
#include "mic_native.h"
#endif


#ifdef __SSE3
#include <stdint.h>
#ifdef __SSE3
#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif
#endif
#include "cycle.h"

static void computeTraversalInfo(nodeptr, traversalInfo *, int *, int, int, pllBoolean, recompVectors *, pllBoolean);
static void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, pllBoolean saveMem, int maxCat, const int states);
#if (defined(__SSE3) && !defined(__AVX))
static void newviewGTRGAMMAPROT_LG4(int tipCase,
                                    double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
                                    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                    int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);

static void newviewGTRGAMMA_GAPPED_SAVE(int tipCase,
                                        double *x1_start, double *x2_start, double *x3_start,
                                        double *EV, double *tipVector,
                                        int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                        const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                        unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
                                        double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn);

static void newviewGTRGAMMA(int tipCase,
                            double *x1_start, double *x2_start, double *x3_start,
                            double *EV, double *tipVector,
                            int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                            const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling
                            );

static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
                           double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                           int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                           int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling);


static void newviewGTRCAT_SAVE( int tipCase,  double *EV,  int *cptr,
                                double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                                int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

static void newviewGTRGAMMAPROT_GAPPED_SAVE(int tipCase,
                                            double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                                            int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                            int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                            unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,  
                                            double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
                                            );

static void newviewGTRGAMMAPROT(int tipCase,
                                double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                                int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling);

static void newviewGTRCATPROT(int tipCase, double *extEV,
                              int *cptr,
                              double *x1, double *x2, double *x3, double *tipVector,
                              int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                              int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling);

static void newviewGTRCATPROT_SAVE(int tipCase, double *extEV,
                                   int *cptr,
                                   double *x1, double *x2, double *x3, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

#endif
#if (defined(__AVX) || defined(__SSE3))
static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);
static void newviewGTRGAMMA_BINARY(int tipCase,
                                   double *x1_start, double *x2_start, double *x3_start,
                                   double *EV, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling);
#endif

/* required to compute the absolute values of double precision numbers with SSE3 */

PLL_ALIGN_BEGIN const union PLL_ALIGN_END
{
  uint64_t i[2];
  __m128d m;
} absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};



#endif

static int pllGetTransitionMatrixNormal (pllInstance * tr, partitionList * pr, nodeptr p, int model, int rate, double * outBuffer);
static int pllGetTransitionMatrixLG4 (partitionList * pr, nodeptr p, int model, double * outBuffer);

extern const char binaryStateNames[2];  /**< @brief Alphabet of binary states */
extern const char dnaStateNames[4];     /**< @brief DNA alphabet  */
extern const char protStateNames[20];   /**< @brief Amino-acid alphabet */
extern const unsigned int mask32[32];   /**< @brief Contains the first 32 powers of 2, i.e. 2^0 upto 2^31 */

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

static void newviewAscCat(int tipCase,
			  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			  int *ex3, 
			  const int n, double *left, double *right, 			    
			  const int numStates)
{
  double
    *le, *ri, *v, *vl, *vr,
    ump_x1, ump_x2, x1px2;
  
  int 
    i, l, j, scale;

 
  unsigned char 
    tip[32];

  ascertainmentBiasSequence(tip, numStates);
  
  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[0];
	    ri = &right[0];

	    vl = &(tipVector[numStates * tip[i]]);
	    vr = &(tipVector[numStates * tip[i]]);
	    v  = &x3[numStates * i];

	    for(l = 0; l < numStates; l++)
	      v[l] = 0.0;

	    for(l = 0; l < numStates; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < numStates; j++)
		  {
		    ump_x1 += vl[j] * le[l * numStates + j];
		    ump_x2 += vr[j] * ri[l * numStates + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < numStates; j++)
		  v[j] += x1px2 * extEV[l * numStates + j];
	      }	    
	  }
      }
      break;
    case PLL_TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[0];
	    ri = &right[0];

	    vl = &(tipVector[numStates * tip[i]]);
	    vr = &x2[numStates * i];
	    v  = &x3[numStates * i];

	    for(l = 0; l < numStates; l++)
	      v[l] = 0.0;

	    for(l = 0; l < numStates; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < numStates; j++)
		  {
		    ump_x1 += vl[j] * le[l * numStates + j];
		    ump_x2 += vr[j] * ri[l * numStates + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < numStates; j++)
		  v[j] += x1px2 * extEV[l * numStates + j];
	      }

	    scale = 1;
	    for(l = 0; scale && (l < numStates); l++)
	      scale = ((v[l] < PLL_MINLIKELIHOOD) && (v[l] > PLL_MINUSMINLIKELIHOOD));	    

	    if(scale)
	      {
		for(l = 0; l < numStates; l++)
		  v[l] *= PLL_TWOTOTHE256;
			
		ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[0];
	  ri = &right[0];

	  vl = &x1[numStates * i];
	  vr = &x2[numStates * i];
	  v = &x3[numStates * i];

	  for(l = 0; l < numStates; l++)
	    v[l] = 0.0;

	  for(l = 0; l < numStates; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < numStates; j++)
		{
		  ump_x1 += vl[j] * le[l * numStates + j];
		  ump_x2 += vr[j] * ri[l * numStates + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < numStates; j++)
		v[j] += x1px2 * extEV[l * numStates + j];
	    }

	   scale = 1;
	   for(l = 0; scale && (l < numStates); l++)
	     scale = ((v[l] < PLL_MINLIKELIHOOD) && (v[l] > PLL_MINUSMINLIKELIHOOD));
	  
	   if(scale)
	     {
	       for(l = 0; l < numStates; l++)
		 v[l] *= PLL_TWOTOTHE256;
	      
	       ex3[i]  += 1;	     
	     }
	}
      break;
    default:
      assert(0);
    }
  
 

}


static void newviewAscGamma(int tipCase,
			    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			    int *ex3, 
			    const int n, double *left, double *right, 			    
			    const int numStates)
{
  
  int  
    i, j, l, k, scale;
  
  const int 
    statesSquare = numStates * numStates,
    gammaStates = 4 * numStates;

  double 
    *vl, *vr, al, ar, *v, x1px2;

  unsigned char 
    tip[32];

  ascertainmentBiasSequence(tip, numStates);
  
  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tip[i]]);
		vr = &(tipVector[numStates * tip[i]]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }	    
	  }
      }
      break;
    case PLL_TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tip[i]]);
		vr = &(x2[gammaStates * i + numStates * k]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }
	   
	    v = &x3[gammaStates * i];
	    scale = 1;
	    for(l = 0; scale && (l < gammaStates); l++)
	      scale = (PLL_ABS(v[l]) < PLL_MINLIKELIHOOD);

	    if(scale)
	      {		
		for(l = 0; l < gammaStates; l++)
		  v[l] *= PLL_TWOTOTHE256;
		
		ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[gammaStates * i + numStates * k]);
	     vr = &(x2[gammaStates * i + numStates * k]);
	     v =  &(x3[gammaStates * i + numStates * k]);

	     for(l = 0; l < numStates; l++)
	       v[l] = 0;

	     for(l = 0; l < numStates; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < numStates; j++)
		   {
		     al += vl[j] * left[k * statesSquare + l * numStates + j];
		     ar += vr[j] * right[k * statesSquare + l * numStates + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < numStates; j++)
		   v[j] += x1px2 * extEV[numStates * l + j];
	       }
	   }
	 
	 v = &(x3[gammaStates * i]);
	 scale = 1;
	 for(l = 0; scale && (l < gammaStates); l++)
	   scale = ((PLL_ABS(v[l]) <  PLL_MINLIKELIHOOD));

	 if(scale)
	   {	    
	     for(l = 0; l < gammaStates; l++)
	       v[l] *= PLL_TWOTOTHE256;
	     
	     ex3[i]  += 1;	    
	   }
       }
      break;
    default:
      assert(0);
    }  
}


/* generic function for computing the P matrices, for computing the conditional likelihood at a node p, given child nodes q and r 
   we compute P(z1) and P(z2) here */

/** @brief Computes two P matrices for two edges.

    Generic function for computing the P matrices of two nodes based on their edges. This is used to 
    (later) compute the the conditional likelihood at a node p which has two descendants \a q and \r, 
    which in turn have the edges \a z1 and \a z2 that connect them with \a p. Given those edges, we
    compute two P matrices \a P(z1) and \a P(z2) which are stored in the arrays \a left and \a right.
 
    The following value is computed here: 
    \f[
     EI\cdot exp( EIGN \cdot z)
     \f]
     to fill up the P matrix.
     
    @param z1    Branch length leading to left descendant node (let's call it \a q)
    @param z2    Branch length leading to right descendant node (let's call it \a r)
    @param rptr  Array of values for rate categories
    @param EI    Inverse eigenvectors of Q-matrix
    @param EIGN  Eigenvalues of Q-matrix
    @param numberOfCategories How many rate heterogeneity categories we have, depending on GAMMA and CAT
    @param left  Where to store the left P matrix (for node \a q)
    @param right Where to store the right P matrix (for node \a r)
    @param saveMem If set to \b PLL_TRUE, memory saving technique is enabled
    @param maxCat Maximum number of rate categories
    @param states Number of states for the particular data (4 for DNA or 20 for AA)
*/
static void 
makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, pllBoolean saveMem, int maxCat, const int states)
{
  int  i, j, k, statesSquare = states * states;

  /* assign some space for pre-computing and later re-using functions */

  double 
    *lz1 = (double*)rax_malloc(sizeof(double) * states),
    *lz2 = (double*)rax_malloc(sizeof(double) * states),
    *d1 = (double*)rax_malloc(sizeof(double) * states),
    *d2 = (double*)rax_malloc(sizeof(double) * states);

  /* multiply branch lengths with eigenvalues */

  for(i = 1; i < states; i++)
  {
    lz1[i] = EIGN[i] * z1;
    lz2[i] = EIGN[i] * z2;
  }


  /* loop over the number of rate categories, this will be 4 for the GAMMA model and 
     variable for the CAT model */

  for(i = 0; i < numberOfCategories; i++)
  {
    /* exponentiate the rate multiplied by the branch */

    for(j = 1; j < states; j++)
    {
      d1[j] = exp(rptr[i] * lz1[j]);
      d2[j] = exp(rptr[i] * lz2[j]);

    }

    /* now fill the P matrices for the two branch length values */

    for(j = 0; j < states; j++)
    {
      /* left and right are pre-allocated arrays */

      left[statesSquare * i  + states * j] = 1.0;
      right[statesSquare * i + states * j] = 1.0;         

      for(k = 1; k < states; k++)
      {
        left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
        right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
      }
    }
  }


  /* if memory saving is enabled and we are using CAT we need to do one additional P matrix 
     calculation for a rate of 1.0 to compute the entries of a column/tree site comprising only gaps */


  if(saveMem)
  {
    i = maxCat;

    for(j = 1; j < states; j++)
    {
      d1[j] = exp (lz1[j]);
      d2[j] = exp (lz2[j]);
    }

    for(j = 0; j < states; j++)
    {
      left[statesSquare * i  + states * j] = 1.0;
      right[statesSquare * i + states * j] = 1.0;

      for(k = 1; k < states; k++)
      {
        left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
        right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
      }
    }
  }

  /* free the temporary buffers */

  rax_free(lz1);
  rax_free(lz2);
  rax_free(d1);
  rax_free(d2);
}


/** Compute the transition probability matrix for a given branch

    Computes the transition probability matrix for the branch \a p->z and partition \a model given the
    PLL instance \a tr and list of partitions \a pr. The result is stored in \a outBuffer which must
    be of sufficient size, i.e states * states * (numberOfRateCategories + 1) * sizeof(double);

    @param tr  PLL instance
    @param pr  List of partitions
    @param model  Partition index for which to take the branch length
    @param p  Adjacent node to the edge we want to compute the trans. prob. matrix
    @param outBuffer Output buffer where to store the transition probability matrix

*/
int pllGetTransitionMatrix (pllInstance * tr, partitionList * pr, nodeptr p, int model, int rate, double * outBuffer)
{
  if (tr->rateHetModel == PLL_CAT)
   {
     if (rate >= pr->partitionData[model]->numberOfCategories) return (PLL_FALSE);
   }
  else
   {
     if (rate >= 4) return (PLL_FALSE);
   }

  if (pr->partitionData[model]->dataType == PLL_AA_DATA &&
		  (pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X))
    return (pllGetTransitionMatrixLG4 (pr, p, model, outBuffer));
    
    
  return (pllGetTransitionMatrixNormal (tr, pr, p, model, rate, outBuffer));
}


/* TODO: Fix this function according to pllGetTransitionMatrixNormal */
static int pllGetTransitionMatrixLG4 (partitionList * pr, nodeptr p, int model, double * outBuffer)
{
  int
    i, j, k,
    states = pr->partitionData[model]->states,
    numberOfCategories = 4;
  double
    d[64],
    *  rptr = pr->partitionData[model]->gammaRates,
    ** EI   = pr->partitionData[model]->EI_LG4,
    ** EIGN = pr->partitionData[model]->EIGN_LG4;

  assert (states == 20);

  for (i = 0; i < numberOfCategories; ++i)
   {
     for (j = 1; j < states; ++j)
      {
        d[j] = exp(rptr[i] * EIGN[i][j] * p->z[model]);
      }
     for (j = 0; j < states; ++ j)
      {
        outBuffer[states * states * i + states * j] = 1.0;
        for (k = 1; k < states; ++k) 
         {
           outBuffer[states * states * i + states * j + k] = d[k] * EI[i][states * j + k];
         }
      }
   }
  return (PLL_TRUE);
}

static int pllGetTransitionMatrixNormal (pllInstance * tr, partitionList * pr, nodeptr p, int model, int rate, double * outBuffer)
{
  int 
    i, j, k,
    /* numberOfCategories, */
    states = pr->partitionData[model]->states;
  double
    * d = (double *)rax_malloc(sizeof(double) * states),
    * rptr,
    * EI   = pr->partitionData[model]->EI,
    * EIGN = pr->partitionData[model]->EIGN,
    * EV = pr->partitionData[model]->EV;
  
  double lz = (p->z[model] > PLL_ZMIN) ? log(p->z[model]) : log(PLL_ZMIN);                        

  if (tr->rateHetModel == PLL_CAT)
   {
     rptr               = pr->partitionData[model]->perSiteRates;
     /* numberOfCategories = pr->partitionData[model]->numberOfCategories; */
   }
  else
   {
     rptr               = pr->partitionData[model]->gammaRates;
     /* numberOfCategories = 4; */
   }

  for (i = 0; i < states * states; ++ i) outBuffer[i] = 0;

  d[0] = 1.0;
  for (j = 1; j < states; ++ j)
   {
     d[j] = exp(rptr[rate] * EIGN[j] * lz);
   }

  for (i = 0; i < states; ++ i)
   {
     for (j = 0; j < states; ++ j)
      {
        for (k = 0; k < states; ++ k)
         {
           outBuffer[states * i + j] += (d[k] * EI[states * i + k] * EV[states * j + k]);
         }
      }
   }

  assert (!tr->saveMemory);
  // TODO: Fix the following snippet
  //if (tr->saveMemory)
  // {
  //   i = tr->maxCategories;
  //   
  //   for (j = 1; j < states; ++j)
  //    {
  //      d[j] = EXP(EIGN[j] * p->z[model]);
  //    }

  //   for (j = 0; j < states; ++j)
  //    {
  //      outBuffer[states * states * i + states * j] = 1.0;
  //      for (k = 1; k < states; ++k)
  //       {
  //         outBuffer[states * states * i + states * j + k] = d[k] * EI[states * j + k];
  //       }
  //    }
  // }

  rax_free(d);

  return (PLL_TRUE);
}


/** @brief Compute two P matrices for two edges for the LG4 model
    
    Computing the P matrices of two nodes based on their edges for the LG4 model. This is used to 
    (later) compute the the conditional likelihood at a node p which has two descendants \a q and \r, 
    which in turn have the edges \a z1 and \a z2 that connect them with \a p. Given those edges, we
    compute two P matrices \a P(z1) and \a P(z2) which are stored in the arrays \a left and \a right.

    @param z1
      Branch length leading to left descendant node (let's call it \a q)
     
    @param z2
      Branch length leading to right descendant node (let's call it \a r)

    @param rptr
      Array of values for rate categories

    @param EI
      Inverse eigenvectors of 4 Q-matrices
     
    @param EIGN
      Eigenvalues of 4 Q-matrix

    @param numberOfCategories
      How many rate heterogeneity categories we have, depending on GAMMA and CAT
     
    @param left
      Where to store the left P matrix (for node \a q)
     
    @param right
      Where to store the right P matrix (for node \a r)

    @param numStates
      Number of states for the particular data (4 for DNA or 20 for AA)

    @todo
      Present the maths here as in ::makeP

*/
static void makeP_FlexLG4(double z1, double z2, double *rptr, double *EI[4],  double *EIGN[4], int numberOfCategories, double *left, double *right, const int numStates)
{
  int 
    i,
    j,
    k;
  
  const int
    statesSquare = numStates * numStates;

  double    
    d1[64],  
    d2[64];

  assert(numStates <= 64);
       
  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 1; j < numStates; j++)
        {
          d1[j] = exp (rptr[i] * EIGN[i][j] * z1);
          d2[j] = exp (rptr[i] * EIGN[i][j] * z2);
        }

      for(j = 0; j < numStates; j++)
        {
          left[statesSquare * i  + numStates * j] = 1.0;
          right[statesSquare * i + numStates * j] = 1.0;

          for(k = 1; k < numStates; k++)
            {
              left[statesSquare * i + numStates * j + k]  = d1[k] * EI[i][numStates * j + k];
              right[statesSquare * i + numStates * j + k] = d2[k] * EI[i][numStates * j + k];
            }
        }
    }  
}

#if (!defined(__AVX) && !defined(__SSE3))

/** @brief Computation of conditional likelihood arrays for CAT
 
    This is a generic, slow but readable function implementation for computing the 
     conditional likelihood arrays at p, given child nodes q and r using the CAT
     mode of rate heterogeneity. Depending whether \a q, resp. \r, are tips or internal
     nodes (indicated by \a tipCase) the conditional likelihoods are computed based on
     \a x1 if \a q is an inner node or \a tipX1 if it is a tip, resp. \a x2 if \a r
     is an inner node or \a tipX2 if it is a tip. Output array \a ex3 stores the
     number of times the likelihood of each site for each internal node has been scaled.
     The conditional likelihood vectors for any possible base-pair (which is useful when
     \a q or \a r are tips) has been already precomputed from the eigenvalues of the Q
     matrix in the array \a tipVector. In case the conditional likelihood for a particular
     site is very small in terms of a floating point number, then it is multiplied by a
     very large number (scaling), and then number of times it has been scaled (per node) is
     stored in the array \a ex3, if \a fastScaling is set to \b PLL_FALSE. Otherwise, the
     total number of scalings for all sites and all nodes is stored in a single variable
     \a scalerIncrement.

    @param tipCase
      Can be either \b PLL_TIP_TIP, or \b PLL_TIP_INNER or \b PLL_INNER_INNER, and describes the
      descendants of the node for which we currently compute the condition likelihood
      vector, i.e. whether they are both tips (leaves), or one is tip and the other
      an inner node, or both are inner nodes.

    @param extEV
      Eigenvectors of Q matrix
      
    @param cptr
      Array where the rate for each site in the compressed partition alignment is stored

    @param x1
      Conditional likelihood vectors of the first child node, in case it is an internal node

    @param x2
      Conditional likelihood vectors of the second child node, in case it is an internal node

    @param x3
      Pointer to where the computed conditional likelihood vector of node \a p will be stored

    @param tipVector
      Vector contining sums of left eigenvectors for likelihood computation at tips.

    @param ex3
      Pointer to an array of whose elements correspond to the number of times the likelihood of
      a particular site of a particular internal nodeis scaled. Those elements are incremented
      at every scaling operation and only if \a fastScaling flag is set to \b PLL_FALSE. This 
      array will be used later when evaluating the likelihood of the whole tree.

    @param tipX1
      Pointer to the alignment data (sequence) of first child node, in case it is a tip

    @param tipX2
      Pointer to the alignment data (sequence) of second child node, in case it is a tip

    @param n
      Number of sites for which we are doing the evaluation. For the single-thread version this is the number of sites in the
      current partition, for multi-threads this is the number of sites assigned to the running thread from the current partition.

    @param left
      Pointer to the P matrix of the left child

    @param right
      Pointer to the P matrix of the right child

    @param wgt
      Array of weights for each site

    @param scalerIncrement
      Where to store the number of scalings carried out in case \a fastScaling is set to \b PLL_TRUE.

    @param fastScaling
      If set to \b PLL_TRUE, only the total number of scalings for all sites of the partition will be
      stored in \a scalerIncrement, otherwise per-site scalings are stored in the array \a ex3. 

    @param states
      Number of states for the particular data (4 for DNA or 20 for AA)
 */
static void newviewCAT_FLEX(int tipCase, double *extEV,
                            int *cptr,
                            double *x1, double *x2, double *x3, double *tipVector,
                            int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                            int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling, const int states)
{
  double
    *le, 
    *ri, 
    *v, 
    *vl, 
    *vr,
    ump_x1, 
    ump_x2, 
    x1px2;

  int 
    i, 
    l, 
    j, 
    scale, 
    addScale = 0;

  const int 
    statesSquare = states * states;


  /* here we switch over the different cases for efficiency, but also because 
     each case accesses different data types.

     We consider three cases: either q and r are both tips, q or r are tips, and q and r are inner 
     nodes.
     */


  switch(tipCase)
  {

    /* both child nodes of p weher we want to update the conditional likelihood are tips */
    case PLL_TIP_TIP:     
      /* loop over sites */
      for (i = 0; i < n; i++)
      {
        /* set a pointer to the P-Matrices for the rate category of this site */
        le = &left[cptr[i] * statesSquare];
        ri = &right[cptr[i] * statesSquare];

        /* pointers to the likelihood entries of the tips q (vl) and r (vr) 
           We will do reading accesses to these values only.
           */
        vl = &(tipVector[states * tipX1[i]]);
        vr = &(tipVector[states * tipX2[i]]);

        /* address of the conditional likelihood array entres at site i. This is 
           a writing access to v */
        v  = &x3[states * i];

        /* initialize v */
        for(l = 0; l < states; l++)
          v[l] = 0.0;

        /* loop over states to compute the cond likelihoods at p (v) */

        for(l = 0; l < states; l++)
        {             
          ump_x1 = 0.0;
          ump_x2 = 0.0;

          /* le and ri are the P-matrices */

          for(j = 0; j < states; j++)
          {
            ump_x1 += vl[j] * le[l * states + j];
            ump_x2 += vr[j] * ri[l * states + j];
          }

          x1px2 = ump_x1 * ump_x2;

          /* multiply with matrix of eigenvectors extEV */

          for(j = 0; j < states; j++)
            v[j] += x1px2 * extEV[l * states + j];
        }          
      }    
      break;
    case PLL_TIP_INNER:      

      /* same as above, only that now vl is a tip and vr is the conditional probability vector 
         at an inner node. Note that, if we have the case that either q or r is a tip, the 
         nodes will be flipped to ensure that tipX1 always points to the sequence at the tip.
         */

      for (i = 0; i < n; i++)
      {
        le = &left[cptr[i] * statesSquare];
        ri = &right[cptr[i] * statesSquare];

        /* access tip vector lookup table */
        vl = &(tipVector[states * tipX1[i]]);

        /* access conditional likelihoo arrays */
        /* again, vl and vr are reading accesses, while v is a writing access */
        vr = &x2[states * i];
        v  = &x3[states * i];

        /* same as in the loop above */

        for(l = 0; l < states; l++)
          v[l] = 0.0;

        for(l = 0; l < states; l++)
        {
          ump_x1 = 0.0;
          ump_x2 = 0.0;

          for(j = 0; j < states; j++)
          {
            ump_x1 += vl[j] * le[l * states + j];
            ump_x2 += vr[j] * ri[l * states + j];
          }

          x1px2 = ump_x1 * ump_x2;

          for(j = 0; j < states; j++)
            v[j] += x1px2 * extEV[l * states + j];
        }

        /* now let's check for numerical scaling. 
           The maths in RAxML are a bit non-standard to avoid/economize on arithmetic operations 
           at the virtual root and for branch length optimization and hence values stored 
           in the conditional likelihood vectors can become negative.
           Below we check if all absolute values stored at position i of v are smaller 
           than a pre-defined value in pll.h. If they are all smaller we can then safely 
           multiply them by a large, constant number PLL_TWOTOTHE256 (without numerical overflow) 
           that is also speced in pll.h */

        scale = 1;
        for(l = 0; scale && (l < states); l++)
          scale = ((v[l] < PLL_MINLIKELIHOOD) && (v[l] > PLL_MINUSMINLIKELIHOOD));         

        if(scale)
        {
          for(l = 0; l < states; l++)
            v[l] *= PLL_TWOTOTHE256;

          /* if we have scaled the entries to prevent underflow, we need to keep track of how many scaling 
             multiplications we did per node such as to undo them at the virtual root, e.g., in 
             evaluateGeneric() 
             Note here, that, if we scaled the site we need to increment the scaling counter by the wieght, i.e., 
             the number of sites this potentially compressed pattern represents ! */ 

          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];   
          
        }
      }   
      break;
    case PLL_INNER_INNER:

      /* same as above, only that the two child nodes q and r are now inner nodes */

      for(i = 0; i < n; i++)
      {
        le = &left[cptr[i] * statesSquare];
        ri = &right[cptr[i] * statesSquare];

        /* index conditional likelihood vectors of inner nodes */

        vl = &x1[states * i];
        vr = &x2[states * i];
        v = &x3[states * i];

        for(l = 0; l < states; l++)
          v[l] = 0.0;

        for(l = 0; l < states; l++)
        {
          ump_x1 = 0.0;
          ump_x2 = 0.0;

          for(j = 0; j < states; j++)
          {
            ump_x1 += vl[j] * le[l * states + j];
            ump_x2 += vr[j] * ri[l * states + j];
          }

          x1px2 =  ump_x1 * ump_x2;

          for(j = 0; j < states; j++)
            v[j] += x1px2 * extEV[l * states + j];            
        }

        scale = 1;
        for(l = 0; scale && (l < states); l++)
          scale = ((v[l] < PLL_MINLIKELIHOOD) && (v[l] > PLL_MINUSMINLIKELIHOOD));

        if(scale)
        {
          for(l = 0; l < states; l++)
            v[l] *= PLL_TWOTOTHE256;
          
          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];    
        }
      }
      break;
    default:
      assert(0);
  }

  /* increment the scaling counter by the additional scalings done at node p */

  if(fastScaling)
    *scalerIncrement = addScale;
}

/** @brief Computation of conditional likelihood arrays for \b GAMMA
 
    This is a generic, slow but readable function implementation for computing the 
     conditional likelihood arrays at \a p, given child nodes \a q and \a r using the \b GAMMA
     model of rate heterogeneity. Depending whether \a q, resp. \r, are tips or internal
     nodes (indicated by \a tipCase) the conditional likelihoods are computed based on
     \a x1 if \a q is an inner node or \a tipX1 if it is a tip, resp. \a x2 if \a r
     is an inner node or \a tipX2 if it is a tip. Output array \a ex3 stores the
     number of times the likelihood of each site for each internal node has been scaled.
     The conditional likelihood vectors for any possible base-pair (which is useful when
     \a q or \a r are tips) has been already precomputed from the eigenvalues of the Q
     matrix in the array \a tipVector. In case the conditional likelihood for a particular
     site is very small in terms of a floating point number, then it is multiplied by a
     very large number (scaling), and then number of times it has been scaled (per node) is
     stored in the array \a ex3, if \a fastScaling is set to \b PLL_FALSE. Otherwise, the
     total number of scalings for all sites and all nodes is stored in a single variable
     \a scalerIncrement.

    @param tipCase
      Can be either \b PLL_TIP_TIP, or \b PLL_TIP_INNER or \b PLL_INNER_INNER, and describes the
      descendants of the node for which we currently compute the condition likelihood
      vector, i.e. whether they are both tips (leaves), or one is tip and the other
      an inner node, or both are inner nodes.

    @param x1
      Conditional likelihood vectors of the first child node, in case it is an internal node

    @param x2
      Conditional likelihood vectors of the second child node, in case it is an internal node

    @param x3
      Pointer to where the computed conditional likelihood vector of node \a p will be stored

    @param extEV
      Eigenvectors of Q matrix

    @param tipVector
      Vector contining sums of left eigenvectors for likelihood computation at tips.

    @param ex3
      Pointer to an array of whose elements correspond to the number of times the likelihood of
      a particular site of a particular internal nodeis scaled. Those elements are incremented
      at every scaling operation and only if \a fastScaling flag is set to \b PLL_FALSE. This 
      array will be used later when evaluating the likelihood of the whole tree.

    @param tipX1
      Pointer to the alignment data (sequence) of first child node, in case it is a tip

    @param tipX2
      Pointer to the alignment data (sequence) of second child node, in case it is a tip

    @param n
      Number of sites to be processed

    @param left
      Pointer to the P matrix of the left child

    @param right
      Pointer to the P matrix of the right child

    @param wgt
      Array of weights for each site

    @param scalerIncrement
      Where to store the number of scalings carried out in case \a fastScaling is set to \b PLL_TRUE.

    @param fastScaling
      If set to \b PLL_TRUE, only the total number of scalings for all sites of the partition will be
      stored in \a scalerIncrement, otherwise per-site scalings are stored in the array \a ex3. 

    @param states
      Number of states for the particular data (4 for DNA or 20 for AA)

    @param maxStateValue
      Number of all possible base-pairs including degenerate characters, i.e. 16 for  DNA and 23 for AA
 */
static void newviewGAMMA_FLEX(int tipCase,
                              double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                              int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                              int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling, const int states, const int maxStateValue)
{
  double  
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr, 
    al, 
    ar;

  int  
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

  const int     
    statesSquare = states * states,
                 span = states * 4,
                 /* this is required for doing some pre-computations that help to save 
                    numerical operations. What we are actually computing here are additional lookup tables 
                    for each possible state a certain data-type can assume.
                    for DNA with ambuguity coding this is 15, for proteins this is 22 or 23, since there 
                    also exist one or two amibguity codes for protein data.
                    Essentially this is very similar to the tip vectors which we also use as lookup tables */
                 precomputeLength = maxStateValue * span;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        /* allocate pre-compute memory space */

        double 
          *umpX1 = (double*)rax_malloc(sizeof(double) * precomputeLength),
          *umpX2 = (double*)rax_malloc(sizeof(double) * precomputeLength);

        /* multiply all possible tip state vectors with the respective P-matrices 
        */

        for(i = 0; i < maxStateValue; i++)
        {
          v = &(tipVector[states * i]);

          for(k = 0; k < span; k++)
          {

            umpX1[span * i + k] = 0.0;
            umpX2[span * i + k] = 0.0;

            for(l = 0; l < states; l++)
            {
              umpX1[span * i + k] +=  v[l] *  left[k * states + l];
              umpX2[span * i + k] +=  v[l] * right[k * states + l];
            }

          }
        }

        for(i = 0; i < n; i++)
        {
          /* access the precomputed arrays (pre-computed multiplication of conditional with the tip state) 
          */

          uX1 = &umpX1[span * tipX1[i]];
          uX2 = &umpX2[span * tipX2[i]];

          /* loop over discrete GAMMA rates */

          for(j = 0; j < 4; j++)
          {
            /* the rest is the same as for CAT */
            v = &x3[i * span + j * states];

            for(k = 0; k < states; k++)
              v[k] = 0.0;

            for(k = 0; k < states; k++)
            {              
              x1px2 = uX1[j * states + k] * uX2[j * states + k];

              for(l = 0; l < states; l++)                                                       
                v[l] += x1px2 * extEV[states * k + l];               
            }

          }        
        }

        /* free precomputed vectors */

        rax_free(umpX1);
        rax_free(umpX2);
      }
      break;
    case PLL_TIP_INNER:
      {
        /* we do analogous pre-computations as above, with the only difference that we now do them 
           only for one tip vector */

        double 
          *umpX1 = (double*)rax_malloc(sizeof(double) * precomputeLength),
          *ump_x2 = (double*)rax_malloc(sizeof(double) * states);

        /* precompute P and left tip vector product */

        for(i = 0; i < maxStateValue; i++)
        {
          v = &(tipVector[states * i]);

          for(k = 0; k < span; k++)
          {

            umpX1[span * i + k] = 0.0;

            for(l = 0; l < states; l++)
              umpX1[span * i + k] +=  v[l] * left[k * states + l];


          }
        }

        for (i = 0; i < n; i++)
        {
          /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */

          uX1 = &umpX1[span * tipX1[i]];

          /* loop over discrete GAMMA rates */

          for(k = 0; k < 4; k++)
          {
            v = &(x2[span * i + k * states]);

            for(l = 0; l < states; l++)
            {
              ump_x2[l] = 0.0;

              for(j = 0; j < states; j++)
                ump_x2[l] += v[j] * right[k * statesSquare + l * states + j];
            }

            v = &(x3[span * i + states * k]);

            for(l = 0; l < states; l++)
              v[l] = 0;

            for(l = 0; l < states; l++)
            {
              x1px2 = uX1[k * states + l]  * ump_x2[l];
              for(j = 0; j < states; j++)
                v[j] += x1px2 * extEV[l * states  + j];
            }
          }

          /* also do numerical scaling as above. Note that here we need to scale 
             4 * 4 values for DNA or 4 * 20 values for protein data.
             If they are ALL smaller than our threshold, we scale. Note that,
             this can cause numerical problems with GAMMA, if the values generated 
             by the four discrete GAMMA rates are too different.

             For details, see: 

             F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees"

*/


          v = &x3[span * i];
          scale = 1;
          for(l = 0; scale && (l < span); l++)
            scale = (PLL_ABS(v[l]) <  PLL_MINLIKELIHOOD);


          if (scale)
          {
            for(l = 0; l < span; l++)
              v[l] *= PLL_TWOTOTHE256;
            
            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];                   
          }
        }

        rax_free(umpX1);
        rax_free(ump_x2);
      }
      break;
    case PLL_INNER_INNER:

      /* same as above, without pre-computations */

      for (i = 0; i < n; i++)
      {
        for(k = 0; k < 4; k++)
        {
          vl = &(x1[span * i + states * k]);
          vr = &(x2[span * i + states * k]);
          v =  &(x3[span * i + states * k]);


          for(l = 0; l < states; l++)
            v[l] = 0;


          for(l = 0; l < states; l++)
          {              

            al = 0.0;
            ar = 0.0;

            for(j = 0; j < states; j++)
            {
              al += vl[j] * left[k * statesSquare + l * states + j];
              ar += vr[j] * right[k * statesSquare + l * states + j];
            }

            x1px2 = al * ar;

            for(j = 0; j < states; j++)
              v[j] += x1px2 * extEV[states * l + j];

          }
        }

        v = &(x3[span * i]);
        scale = 1;
        for(l = 0; scale && (l < span); l++)
          scale = ((PLL_ABS(v[l]) <  PLL_MINLIKELIHOOD));

        if(scale)
        {  
          for(l = 0; l < span; l++)
            v[l] *= PLL_TWOTOTHE256;
          
          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];           
        }
      }
      break;
    default:
      assert(0);
  }

  /* as above, increment the global counter that counts scaling multiplications by the scaling multiplications 
     carried out for computing the likelihood array at node p */

  if(fastScaling)
    *scalerIncrement = addScale;
}


/* Candidate for deletion */
/*
static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
                           double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                           int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                           int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  double
    ump_x1, ump_x2, x1px2[4];
  int i, j, k, scale, addScale = 0;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[4 * tipX1[i]]);
            x2 = &(tipVector[4 * tipX2[i]]);
            x3 = &x3_start[4 * i];

            le =  &left[cptr[i] * 16];
            ri =  &right[cptr[i] * 16];

            for(j = 0; j < 4; j++)
              {
                ump_x1 = 0.0;
                ump_x2 = 0.0;
                for(k = 0; k < 4; k++)
                  {
                    ump_x1 += x1[k] * le[j * 4 + k];
                    ump_x2 += x2[k] * ri[j * 4 + k];
                  }
                x1px2[j] = ump_x1 * ump_x2;
              }

            for(j = 0; j < 4; j++)
              x3[j] = 0.0;

            for(j = 0; j < 4; j++)
              for(k = 0; k < 4; k++)
                x3[k] += x1px2[j] * EV[j * 4 + k];          
          }
      }
      break;
    case PLL_TIP_INNER:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[4 * tipX1[i]]);
            x2 = &x2_start[4 * i];
            x3 = &x3_start[4 * i];

            le =  &left[cptr[i] * 16];
            ri =  &right[cptr[i] * 16];

            for(j = 0; j < 4; j++)
              {
                ump_x1 = 0.0;
                ump_x2 = 0.0;
                for(k = 0; k < 4; k++)
                  {
                    ump_x1 += x1[k] * le[j * 4 + k];
                    ump_x2 += x2[k] * ri[j * 4 + k];
                  }
                x1px2[j] = ump_x1 * ump_x2;
              }

            for(j = 0; j < 4; j++)
              x3[j] = 0.0;

            for(j = 0; j < 4; j++)
              for(k = 0; k < 4; k++)
                x3[k] +=  x1px2[j] *  EV[4 * j + k];       

            scale = 1;
            for(j = 0; j < 4 && scale; j++)
              scale = (x3[j] < PLL_MINLIKELIHOOD && x3[j] > PLL_MINUSMINLIKELIHOOD);               
                    
            if(scale)
              {             
                for(j = 0; j < 4; j++)
                  x3[j] *= PLL_TWOTOTHE256;
                
                if(useFastScaling)
                  addScale += wgt[i];
                else
                  ex3[i]  += 1;         
              }      
          }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &x1_start[4 * i];
          x2 = &x2_start[4 * i];
          x3 = &x3_start[4 * i];

          le = &left[cptr[i] * 16];
          ri = &right[cptr[i] * 16];

          for(j = 0; j < 4; j++)
            {
              ump_x1 = 0.0;
              ump_x2 = 0.0;
              for(k = 0; k < 4; k++)
                {
                  ump_x1 += x1[k] * le[j * 4 + k];
                  ump_x2 += x2[k] * ri[j * 4 + k];
                }
              x1px2[j] = ump_x1 * ump_x2;
            }

          for(j = 0; j < 4; j++)
            x3[j] = 0.0;

          for(j = 0; j < 4; j++)
            for(k = 0; k < 4; k++)
              x3[k] +=  x1px2[j] *  EV[4 * j + k];
        
          scale = 1;
          for(j = 0; j < 4 && scale; j++)
            scale = (x3[j] < PLL_MINLIKELIHOOD && x3[j] > PLL_MINUSMINLIKELIHOOD);

          if(scale)
            {               
              for(j = 0; j < 4; j++)
                x3[j] *= PLL_TWOTOTHE256;
              
              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;           
            }     
        }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}
*/
#if 0
static void newviewGTRGAMMA_BINARY(int tipCase,
                                   double *x1_start, double *x2_start, double *x3_start,
                                   double *EV, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling
                                   )
{
  double
    *x1, *x2, *x3;
  double
    ump_x1,
    ump_x2,
    x1px2[4];
  int i, j, k, l, scale, addScale = 0;


  /* C-OPT figure out if we are at an inner node who has two tips/leaves
     as descendants TIP_TIP, a tip and another inner node as descendant
     TIP_INNER, or two inner nodes as descendants INNER_INNER */

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &(tipVector[2 * tipX2[i]]);
            x3 = &x3_start[i * 8];

            for(j = 0; j < 8; j++)
              x3[j] = 0.0;

            for (j = 0; j < 4; j++)
              {
                for (k = 0; k < 2; k++)
                  {
                    ump_x1 = 0.0;
                    ump_x2 = 0.0;

                    for (l=0; l < 2; l++)
                      {
                        ump_x1 += x1[l] * left[ j*4 + k*2 + l];
                        ump_x2 += x2[l] * right[j*4 + k*2 + l];
                      }

                    x1px2[k] = ump_x1 * ump_x2;
                  }

                for(k = 0; k < 2; k++)
                  for (l = 0; l < 2; l++)
                    x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

              }    
          }
      }
      break;
    case PLL_TIP_INNER:
      {
         for (i = 0; i < n; i++)
           {
             x1 = &(tipVector[2 * tipX1[i]]);
             x2 = &x2_start[i * 8];
             x3 = &x3_start[i * 8];

             for(j = 0; j < 8; j++)
               x3[j] = 0.0;

             for (j = 0; j < 4; j++)
               {
                 for (k = 0; k < 2; k++)
                   {
                     ump_x1 = 0.0;
                     ump_x2 = 0.0;

                     for (l=0; l < 2; l++)
                       {
                         ump_x1 += x1[l] * left[ j*4 + k*2 + l];
                         ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
                       }

                     x1px2[k] = ump_x1 * ump_x2;
                   }

                 for(k = 0; k < 2; k++)
                   for (l = 0; l < 2; l++)
                     x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

               }            

             scale = 1;
             for(l = 0; scale && (l < 8); l++)
               scale = (PLL_ABS(x3[l]) <  PLL_MINLIKELIHOOD);

             if(scale)
               {
                 for (l=0; l < 8; l++)
                   x3[l] *= PLL_TWOTOTHE256;
                 
                 if(useFastScaling)
                   addScale += wgt[i];
                 else
                   ex3[i]  += 1;               
               }

           }
      }
      break;
    case PLL_INNER_INNER:

      /* C-OPT here we don't do any pre-computations
         This should be the most compute intensive loop of the three
         cases here. If we have one or two tips as descendants
         we can take a couple of shortcuts */


     for (i = 0; i < n; i++)
       {
         x1 = &x1_start[i * 8];
         x2 = &x2_start[i * 8];
         x3 = &x3_start[i * 8];

         for(j = 0; j < 8; j++)
           x3[j] = 0.0;

         for (j = 0; j < 4; j++)
           {
             for (k = 0; k < 2; k++)
               {
                 ump_x1 = 0.0;
                 ump_x2 = 0.0;

                 for (l=0; l < 2; l++)
                   {
                     ump_x1 += x1[j*2 + l] * left[ j*4 + k*2 + l];
                     ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
                   }

                 x1px2[k] = ump_x1 * ump_x2;
               }

             for(k = 0; k < 2; k++)
               for (l = 0; l < 2; l++)
                 x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

           }
         
         scale = 1;
         for(l = 0; scale && (l < 8); l++)
           scale = (PLL_ABS(x3[l]) <  PLL_MINLIKELIHOOD);


         if(scale)
           {
             for (l=0; l<8; l++)
               x3[l] *= PLL_TWOTOTHE256;

             if(useFastScaling)
               addScale += wgt[i];
             else
               ex3[i]  += 1;      
           }
       }
     break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
				  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
				  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  double
    ump_x1, ump_x2, x1px2[2];
  int i, j, k, scale, addScale = 0;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &(tipVector[2 * tipX2[i]]);
	    x3 = &x3_start[2 * i];	    

	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 2; j++)
	      for(k = 0; k < 2; k++)
		x3[k] += x1px2[j] * EV[j * 2 + k];	   
	  }
      }
      break;
    case PLL_TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[2 * i];
	    x3 = &x3_start[2 * i];
	    
	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 2; j++)
	      for(k = 0; k < 2; k++)
		x3[k] +=  x1px2[j] *  EV[2 * j + k];	   

	    scale = 1;
	    for(j = 0; j < 2 && scale; j++)
	      scale = (x3[j] < PLL_MINLIKELIHOOD && x3[j] > PLL_MINUSMINLIKELIHOOD);

	    if(scale)
	      {
		for(j = 0; j < 2; j++)
		  x3[j] *= PLL_TWOTOTHE256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	       
	      }
	  }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  x3 = &x3_start[2 * i];

	  le = &left[cptr[i] * 4];
	  ri = &right[cptr[i] * 4];

	  for(j = 0; j < 2; j++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      for(k = 0; k < 2; k++)
		{
		  ump_x1 += x1[k] * le[j * 2 + k];
		  ump_x2 += x2[k] * ri[j * 2 + k];
		}
	      x1px2[j] = ump_x1 * ump_x2;
	    }

	  for(j = 0; j < 2; j++)
	    x3[j] = 0.0;

	  for(j = 0; j < 2; j++)
	    for(k = 0; k < 2; k++)
	      x3[k] +=  x1px2[j] *  EV[2 * j + k];	  

	  scale = 1;
	  for(j = 0; j < 2 && scale; j++)
	    scale = (x3[j] < PLL_MINLIKELIHOOD && x3[j] > PLL_MINUSMINLIKELIHOOD);

	  if(scale)
	    {
	      for(j = 0; j < 2; j++)
		x3[j] *= PLL_TWOTOTHE256;

	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	   
	    }
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}
#endif    /* end if 0 */
#endif

#if (defined(__AVX) || defined(__SSE3))
static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
                                  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
                                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  int i, l, scale, addScale = 0;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
        for(i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &(tipVector[2 * tipX2[i]]);
            x3 = &x3_start[2 * i];         

            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }            
          }
      }
      break;
    case PLL_TIP_INNER:
      {
        for (i = 0; i < n; i++)
          {
            x1 = &(tipVector[2 * tipX1[i]]);
            x2 = &x2_start[2 * i];
            x3 = &x3_start[2 * i];
            
            le =  &left[cptr[i] * 4];
            ri =  &right[cptr[i] * 4];

            _mm_store_pd(x3, _mm_setzero_pd());     
                     
            for(l = 0; l < 2; l++)
              {                                                                                                                          
                __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
                __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
                
                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);
                
                al = _mm_mul_pd(al, ar);
                
                __m128d vv  = _mm_load_pd(x3);
                __m128d EVV = _mm_load_pd(&EV[2 * l]);
                
                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
                
                _mm_store_pd(x3, vv);                                                     
              }  
            
            __m128d minlikelihood_sse = _mm_set1_pd(PLL_MINLIKELIHOOD);
         
            scale = 1;
            
            __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;                         
            
            if(scale)
              {
                __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
                
                __m128d ex3v = _mm_load_pd(x3);           
                _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                                 
                
                if(useFastScaling)
                  addScale += wgt[i];
                else
                  ex3[i]  += 1;   
              }                    
          }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
        {
          x1 = &x1_start[2 * i];
          x2 = &x2_start[2 * i];
          x3 = &x3_start[2 * i];

          le = &left[cptr[i] * 4];
          ri = &right[cptr[i] * 4];

          _mm_store_pd(x3, _mm_setzero_pd());       
          
          for(l = 0; l < 2; l++)
            {                                                                                                                            
              __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
              __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
              
              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);
              
              al = _mm_mul_pd(al, ar);
              
              __m128d vv  = _mm_load_pd(x3);
              __m128d EVV = _mm_load_pd(&EV[2 * l]);
              
              vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
              
              _mm_store_pd(x3, vv);                                                       
            }                             

          __m128d minlikelihood_sse = _mm_set1_pd(PLL_MINLIKELIHOOD);
         
          scale = 1;
                  
          __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;                   
         
          if(scale)
            {
              __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
                    
              __m128d ex3v = _mm_load_pd(x3);             
              _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));                                           
             
              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;     
           }             
        }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling
				   )
{
  double
    *x1, *x2, *x3;
 
  int i, k, l, scale, addScale = 0; 

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 x2  = &(tipVector[2 * tipX2[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     	    
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
       }
      break;
    case PLL_TIP_INNER:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }      
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
       {	 
	 for(k = 0; k < 4; k++)
	   {	     
	     x1 = &(x1_start[8 * i + 2 * k]);
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }
      break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


#endif




/* The function below computes partial traversals only down to the point/node in the tree where the 
   conditional likelihhod vector summarizing a subtree is already oriented in the correct direction */


/** @brief Compute a partial or full traversal descriptor for a subtree of the topology

   Unless the \a partialTraversal is set to \b PLL_TRUE, compute a partial traversal descriptor down 
   to the point/node in the tree where the conditional likelihood vector representing a subtree is
   already oriented in the correct direction. The elements of the traversal descriptor are stored in
   \a ti and a \a counter keeps track of the number of elements.

   @param p
     Root of the  subtree for which we want to compute the traversal descriptor. The two descendents are \a p->next->back and \a p->next->next->back

   @param ti
i    Traversal descriptor element structure

   @param counter
     Number of elements in the traversal descriptor. Updated when an element is added

   @param maxTips
     Number of tips in the tree structure

   @param numBranches
     Number of branches
   
   @param partialTraversal
     If \b PLL_TRUE, a partial traversal descriptor is computed, otherwise a full

   @param rvec
     Parameter concerning ancestral state recomputation. Please document

   @param useRecom
     If \b PLL_TRUE, then ancestral state recomputation is enabled.
   
   @todo Fill in the ancestral recomputation parameter information 
 */
static void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, pllBoolean partialTraversal, recompVectors *rvec, pllBoolean useRecom)
{
  /* if it's a tip we don't do anything */

  if(isTip(p->number, maxTips))
    return;

  {
    int 
      i;

    /* recom default values */
    int slot = -1,
        unpin1 = -1, 
        unpin2 = -1;
    /* get the left and right descendants */

    nodeptr 
      q = p->next->back,
        r = p->next->next->back;   

    /* if the left and right children are tips there is not that much to do */
    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
    {
      /* fix the orientation of p->x */

      if (! p->x)
        getxnode(p);    
      
      assert(p->x);

      /* add the current node triplet p,q,r to the traversal descriptor */
      ti[*counter].tipCase = PLL_TIP_TIP;
      ti[*counter].pNumber = p->number;
      ti[*counter].qNumber = q->number;
      ti[*counter].rNumber = r->number;


      /* copy branches to traversal descriptor */
      for(i = 0; i < numBranches; i++)
      {     
        ti[*counter].qz[i] = q->z[i];
        ti[*counter].rz[i] = r->z[i];
      }

      /* recom - add the slot to the traversal descriptor */
      if(useRecom)
      {
        getxVector(rvec, p->number, &slot, maxTips);
        ti[*counter].slot_p = slot;
        ti[*counter].slot_q = -1;
        ti[*counter].slot_r = -1;
      }

      /* increment length counter */

      *counter = *counter + 1;
    }
    else
    {
      /* if either r or q are tips, flip them to make sure that the tip data is stored 
         for q */
      if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
      {     
        if(isTip(r->number, maxTips))
        {
          nodeptr 
            tmp = r;
          r = q;
          q = tmp;
        }


        /* if the orientation of the liklihood vector at r is not correct we need to re-compute it 
           and descend into its subtree to figure out if there are more vrctors in there to re-compute and 
           re-orient */

        if(needsRecomp(useRecom, rvec, r, maxTips) || !partialTraversal) 
          computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
        else
          {
            if(useRecom)
              /* the node is available,  now make sure it will not be unpinned until it is read */
              protectNode(rvec, r->number, maxTips);
          }
        /* Now that r is oriented, we can safely set the orientation of p */
        if(! p->x)
          getxnode(p);   

        /* make sure that everything is consistent now */

        assert(p->x && r->x);

        /* store data for p, q, r in the traversal descriptor */

        ti[*counter].tipCase = PLL_TIP_INNER;
        ti[*counter].pNumber = p->number;
        ti[*counter].qNumber = q->number;
        ti[*counter].rNumber = r->number;

        for(i = 0; i < numBranches; i++)
        {       
          ti[*counter].qz[i] = q->z[i];
          ti[*counter].rz[i] = r->z[i];
        }

        if(useRecom)
        {
          getxVector(rvec, r->number, &slot, maxTips);
          ti[*counter].slot_r = slot;

          getxVector(rvec, p->number, &slot, maxTips);
          ti[*counter].slot_p = slot;

          ti[*counter].slot_q = -1;

          unpin2 = r->number; /* when PLL_TIP_INNER finishes, the INNER input vector r can be unpinned*/
        }

        *counter = *counter + 1;
      }
      else
      {
        /* same as above, only now q and r are inner nodes. Hence if they are not 
           oriented correctly they will need to be recomputed and we need to descend into the 
           respective subtrees to check if everything is consistent in there, potentially expanding 
           the traversal descriptor */
        if(( useRecom && (!partialTraversal) ) || 
            ( useRecom && needsRecomp(useRecom, rvec, q, maxTips) && needsRecomp(useRecom, rvec, r, maxTips) ))
        {
          /* PLL_INNER_INNER and recomputation implies that the order we descend q and r matters, 
           * if we are in a partial traversal, this is only relevant if both require recomputation
           * see TODOFER add ref. */

          int q_stlen = rvec->stlen[q->number - maxTips - 1],
              r_stlen = rvec->stlen[q->number - maxTips - 1];
          assert(q_stlen >= 2 && q_stlen <= maxTips - 1);
          assert(r_stlen >= 2 && r_stlen <= maxTips - 1);

          if(q_stlen > r_stlen)
          {
            computeTraversalInfo(q, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
            computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
          }
          else
          {
            computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
            computeTraversalInfo(q, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
          }
        }
        else
        {
          /* Now the order does not matter */
          /* If we are in a recomputation and partial, only either q or r will be descended */

          if(!partialTraversal || needsRecomp(useRecom, rvec, q, maxTips))
            computeTraversalInfo(q, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
          else
          {
            if(useRecom)
              /* the node is available,  now make sure it will not be unpinned until it is read */
              protectNode(rvec, q->number, maxTips);
          }

          if(!partialTraversal || needsRecomp(useRecom, rvec, r, maxTips))
            computeTraversalInfo(r, ti, counter, maxTips, numBranches, partialTraversal, rvec, useRecom);
          else
          {
            if(useRecom)
              protectNode(rvec, r->number, maxTips);
          }
        }


        if(! p->x)
          getxnode(p);

        /* check that the vector orientations are consistent now */

        assert(p->x && r->x && q->x);

        ti[*counter].tipCase = PLL_INNER_INNER;
        ti[*counter].pNumber = p->number;
        ti[*counter].qNumber = q->number;
        ti[*counter].rNumber = r->number;

        if(useRecom)
        {
          /* We check that the strategy cannot re-use slots */
          getxVector(rvec, q->number, &slot, maxTips);
          ti[*counter].slot_q = slot;

          getxVector(rvec, r->number, &slot, maxTips);
          ti[*counter].slot_r = slot;
          assert(slot != ti[*counter].slot_q);

          getxVector(rvec, p->number, &slot, maxTips);
          ti[*counter].slot_p = slot;
          assert(slot != ti[*counter].slot_q);
          assert(slot != ti[*counter].slot_r);

          /* And at these point both input INNER can be marked as unpinned */
          unpin2 = r->number;
          unpin1 = q->number;
        }

        for(i = 0; i < numBranches; i++)
        {       
          ti[*counter].qz[i] = q->z[i];
          ti[*counter].rz[i] = r->z[i];
        }

        *counter = *counter + 1;
      }
    }
    if(useRecom)
    {
      /* Mark the nodes as unpinnable(will be unpinned while executing the replacement strategy only if required)*/
      unpinNode(rvec, unpin1, maxTips);
      unpinNode(rvec, unpin2, maxTips);
    }
  }
}

/* below are the optimized unrolled, and vectorized versions of the above generi cfunctions 
   for computing the conditional likelihood at p given child nodes q and r. The actual implementation is located at the end/bottom of this 
   file.
   */
/* now this is the function that just iterates over the length of the traversal descriptor and 
   just computes the conditional likelihhod arrays in the order given by the descriptor.
   So in a sense, this function has no clue that there is any tree-like structure 
   in the traversal descriptor, it just operates on an array of structs of given length */ 


/** @brief Compute the conditional likelihood for each entry (node) of the traversal descriptor

    Computes the conditional likelihood vectors for each entry (node) in the already computed
    traversal descriptor, starting from the \a startIndex entry.
     
    @param tr
      PLL instance

    @param pr
      List of partitions

    @param startIndex
      From which node to start computing the conditional likelihood vectors in the traversal
      descriptor
     
    @note This function just iterates over the length of the traversal descriptor and 
      computes the conditional likelihhod arrays in the order given by the descriptor.
      So in a sense, this function has no clue that there is any tree-like structure 
      in the traversal descriptor, it just operates on an array of structs of given length.
 */
void pllNewviewIterative (pllInstance *tr, partitionList *pr, int startIndex)
{
  traversalInfo 
    *ti   = tr->td[0].ti;

  int 
    i, 
    model;

  int 
    p_slot = -1, 
    q_slot = -1, 
    r_slot = -1;

#ifdef _DEBUG_RECOMPUTATION
  /* recom */
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
#else
  countTraversal(tr);
#endif
  /* E recom */
#endif

  /* loop over traversal descriptor length. Note that on average we only re-compute the conditionals on 3 -4 
     nodes in RAxML */

  for(i = startIndex; i < tr->td[0].count; i++)
  {

    traversalInfo 
      *tInfo = &ti[i];
    
    /* Note that the slots refer to different things if recomputation is applied */
    if(tr->useRecom)
      {
        /* a slot has been assigned while computing the traversal descriptor  */
        p_slot = tInfo->slot_p;
        q_slot = tInfo->slot_q;
        r_slot = tInfo->slot_r;
      }
    else
      {
        /* a fixed slot is always given for each inner node, we only need an offset to get the right index */
        p_slot = tInfo->pNumber - tr->mxtips - 1;
        q_slot = tInfo->qNumber - tr->mxtips - 1;
        r_slot = tInfo->rNumber - tr->mxtips - 1;
      }

    /* now loop over all partitions for nodes p, q, and r of the current traversal vector entry */

    for(model = 0; model < pr->numberOfPartitions; model++)
    {
      /* number of sites in this partition */
      size_t            
        width  = (size_t)pr->partitionData[model]->width;

      /* this conditional statement is exactly identical to what we do in pllEvaluateIterative */

      if(tr->td[0].executeModel[model] && width > 0)
      {       
        double
          *x1_start = (double*)NULL,
          *x2_start = (double*)NULL,
          *x3_start = pr->partitionData[model]->xVector[p_slot],
          *left     = (double*)NULL,
          *right    = (double*)NULL,            
#if (defined(__SSE3) || defined(__AVX))
          *x1_gapColumn = (double*)NULL,
          *x2_gapColumn = (double*)NULL,
          *x3_gapColumn = (double*)NULL,
#endif
          *rateCategories = (double*)NULL,
          *x1_ascColumn = NULL,
          *x2_ascColumn = NULL,
          *x3_ascColumn = NULL;

        int
          categories,
          scalerIncrement = 0,

          /* integer wieght vector with pattern compression weights */

          *wgt = pr->partitionData[model]->wgt;

        /* pointers for per-site scaling array at node p */
        
        int      
          *ex3     = NULL,
          *ex3_asc = NULL;

        /* select fastScaling or per-site scaling of conidtional likelihood entries */

        pllBoolean
          fastScaling = tr->fastScaling;

#if (defined(__SSE3) || defined(__AVX))
        unsigned int
          *x1_gap = (unsigned int*)NULL,
          *x2_gap = (unsigned int*)NULL,
          *x3_gap = (unsigned int*)NULL;
#endif

        unsigned char
          *tipX1 = (unsigned char *)NULL,
          *tipX2 = (unsigned char *)NULL;

        double 
          qz, 
          rz;        

        size_t
#if (defined(__SSE3) || defined(__AVX))
          gapOffset = 0,
#endif
          rateHet = discreteRateCategories(tr->rateHetModel),
          ascWidth = (size_t)pr->partitionData[model]->states,

          /* get the number of states in the data stored in partition model */
          
          states = (size_t)pr->partitionData[model]->states,
          
          /* get the length of the current likelihood array stored at node p. This is 
             important mainly for the SEV-based memory saving option described in here:
             
             F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees".
             
             So pr->partitionData[model]->xSpaceVector[i] provides the length of the allocated conditional array of partition model
             and node i 
          */
          
          availableLength = pr->partitionData[model]->xSpaceVector[p_slot],
          requiredLength = 0;        
        
        /* figure out what kind of rate heterogeneity approach we are using */

        if(tr->rateHetModel == PLL_CAT)
          {              
            rateCategories = pr->partitionData[model]->perSiteRates;
            categories = pr->partitionData[model]->numberOfCategories;
          }
        else
          {                              
            rateCategories = pr->partitionData[model]->gammaRates;
            categories = 4;
          }

        /* memory saving stuff, not important right now, but if you are interested ask Fernando */

#if (defined(__SSE3) || defined(__AVX))
        if(tr->saveMemory)
          {
            size_t
              j,
              setBits = 0;                
            
            gapOffset = states * (size_t)getUndetermined(pr->partitionData[model]->dataType);
            
            x1_gap = &(pr->partitionData[model]->gapVector[tInfo->qNumber * pr->partitionData[model]->gapVectorLength]);
            x2_gap = &(pr->partitionData[model]->gapVector[tInfo->rNumber * pr->partitionData[model]->gapVectorLength]);
            x3_gap = &(pr->partitionData[model]->gapVector[tInfo->pNumber * pr->partitionData[model]->gapVectorLength]);
            
            for(j = 0; j < (size_t)pr->partitionData[model]->gapVectorLength; j++)
              {              
                x3_gap[j] = x1_gap[j] & x2_gap[j];
                setBits += (size_t)(bitcount_32_bit(x3_gap[j])); 
              }
            
            requiredLength = (width - setBits)  * rateHet * states * sizeof(double);            
          }
        else
#endif
          {
            /* if we are not trying to save memory the space required to store an inner likelihood array 
               is the number of sites in the partition times the number of states of the data type in the partition 
               times the number of discrete GAMMA rates (1 for CAT essentially) times 8 bytes */
            requiredLength  =  virtual_width( width ) * rateHet * states * sizeof(double);
            
            //                   printf( "req: %d %d %d %d\n", requiredLength, width, virtual_width(width), model );
          }
        
        /* Initially, even when not using memory saving no space is allocated for inner likelihood arrats hence 
           availableLength will be zero at the very first time we traverse the tree.
           Hence we need to allocate something here */

        if(requiredLength != availableLength)
          {               
            /* if there is a vector of incorrect length assigned here i.e., x3 != NULL we must free 
               it first */
            if(x3_start)
              rax_free(x3_start);
            
            /* allocate memory: note that here we use a byte-boundary aligned malloc, because we need the vectors
               to be aligned at 16 BYTE (SSE3) or 32 BYTE (AVX) boundaries! */
            
            rax_posix_memalign ((void **)&x3_start, PLL_BYTE_ALIGNMENT, requiredLength);              
            
            /* update the data structures for consistent bookkeeping */
            pr->partitionData[model]->xVector[p_slot]      = x3_start;
            pr->partitionData[model]->xSpaceVector[p_slot] = requiredLength;
          }
        

        /* 
           if we are not using fast scaling, we need to assign memory for storing 
           integer vectors at each inner node that are as long as the sites of the 
           partition. IMPORTANT: while this looks as if this might be a memory saving trick 
           it is not. The ex3 vectors will be allocated once during the very first tree 
           traversal and then never again because they will always have the required length!
        */

        if(!fastScaling)
          {
            size_t
              availableExpLength = pr->partitionData[model]->expSpaceVector[p_slot],
              requiredExpLength  = width * sizeof(int);
            
            ex3 = pr->partitionData[model]->expVector[p_slot];
            
            if(requiredExpLength != availableExpLength)
              {
                if(ex3)
                  rax_free(ex3);
                
                rax_posix_memalign ((void **)&ex3, PLL_BYTE_ALIGNMENT, requiredExpLength);               
                
                pr->partitionData[model]->expVector[p_slot] = ex3;
                
                pr->partitionData[model]->expSpaceVector[p_slot] = requiredExpLength;
              }
          }

        /* now just set the pointers for data accesses in the newview() implementations above to the corresponding values 
           according to the tip case */
        
        switch(tInfo->tipCase)
          {
          case PLL_TIP_TIP:           
            tipX1    = pr->partitionData[model]->yVector[tInfo->qNumber];
            tipX2    = pr->partitionData[model]->yVector[tInfo->rNumber];

#if (defined(__SSE3) || defined(__AVX))
            if(tr->saveMemory)
              {
                x1_gapColumn   = &(pr->partitionData[model]->tipVector[gapOffset]);
                x2_gapColumn   = &(pr->partitionData[model]->tipVector[gapOffset]);
                x3_gapColumn   = &(pr->partitionData[model]->gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet]);
              }
#endif            

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
            if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
            if(pr->partitionData[model]->ascBias)
#endif
             {
              size_t
                k;
              
              x3_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
              ex3_asc      = &pr->partitionData[model]->ascExpVector[(tInfo->pNumber - tr->mxtips - 1) * ascWidth];

              for(k = 0; k < ascWidth; k++)
                ex3_asc[k] = 0;               
             }
            /* if we do per-site log likelihood scaling, and both child nodes are tips,
               just initialize the vector with zeros, i.e., no scaling events */

            if(!fastScaling)
              {
                size_t
                  k;                                 

                for(k = 0; k < width; k++)
                  ex3[k] = 0;
              }
            break;
          case PLL_TIP_INNER:                
            tipX1    =  pr->partitionData[model]->yVector[tInfo->qNumber];
            x2_start = pr->partitionData[model]->xVector[r_slot];
            assert(r_slot != p_slot);
            
#if (defined(__SSE3) || defined(__AVX))
            if(tr->saveMemory)
              { 
                x1_gapColumn   = &(pr->partitionData[model]->tipVector[gapOffset]);
                x2_gapColumn   = &pr->partitionData[model]->gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
                x3_gapColumn   = &pr->partitionData[model]->gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
              }
#endif

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
            if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
              if(pr->partitionData[model]->ascBias)
#endif      
              {   
                size_t
                  k;

                int 
                  *ex2_asc;
                
                x2_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->rNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                x3_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                
                ex2_asc = &pr->partitionData[model]->ascExpVector[(tInfo->rNumber - tr->mxtips - 1) * ascWidth];
                ex3_asc = &pr->partitionData[model]->ascExpVector[(tInfo->pNumber - tr->mxtips - 1) * ascWidth];

                for(k = 0; k < ascWidth; k++)
                  ex3_asc[k] = ex2_asc[k];
              }
            
            /* if one child node is not a tip, just copy the values from there, coudl also be done with memcpy of course 
               the elements of ex3[] will then potentially be further incremented in the actual newview() if scaling events 
               take place */

            if(!fastScaling)
              {
                size_t 
                  k;
                int
                  *ex2 = pr->partitionData[model]->expVector[r_slot];                
                      
                for(k = 0; k < width; k++)
                  ex3[k] = ex2[k];
              }
            break;
          case PLL_INNER_INNER:                              
            x1_start       = pr->partitionData[model]->xVector[q_slot];
            x2_start       = pr->partitionData[model]->xVector[r_slot];
            assert(r_slot != p_slot);
            assert(q_slot != p_slot);
            assert(q_slot != r_slot);
            
#if (defined(__SSE3) || defined(__AVX))
            if(tr->saveMemory)
              {
                x1_gapColumn   = &pr->partitionData[model]->gapColumn[(tInfo->qNumber - tr->mxtips - 1) * states * rateHet];
                x2_gapColumn   = &pr->partitionData[model]->gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
                x3_gapColumn   = &pr->partitionData[model]->gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
              }
#endif

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
              if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
              if(pr->partitionData[model]->ascBias)
#endif          
               {                
                 size_t
                   k;

                 int 
                   *ex1_asc,
                   *ex2_asc;
                 
                 x1_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->qNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                 x2_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->rNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                 x3_ascColumn = &pr->partitionData[model]->ascVector[(tInfo->pNumber - tr->mxtips - 1) * pr->partitionData[model]->ascOffset];
                 
                 ex1_asc = &pr->partitionData[model]->ascExpVector[(tInfo->qNumber - tr->mxtips - 1) * ascWidth];
                 ex2_asc = &pr->partitionData[model]->ascExpVector[(tInfo->rNumber - tr->mxtips - 1) * ascWidth];
                 ex3_asc = &pr->partitionData[model]->ascExpVector[(tInfo->pNumber - tr->mxtips - 1) * ascWidth];

                 for(k = 0; k < ascWidth; k++)
                   ex3_asc[k] = ex1_asc[k] + ex2_asc[k];
               }
            /* both child nodes are inner nodes, thus the initial value of the scaling vector 
               ex3 is the sum of the scaling values of the left and right child node */

            if(!fastScaling)
              {
                size_t
                  k;
                      
                int            
                  *ex1      = pr->partitionData[model]->expVector[q_slot],
                  *ex2      = pr->partitionData[model]->expVector[r_slot];                    
                      
                  for(k = 0; k < width; k++)
                    ex3[k] = ex1[k] + ex2[k];
              }
            break;
          default:
            assert(0);
          }

        /* set the pointers to the left and right P matrices to the pre-allocated memory space for storing them */

        left  = pr->partitionData[model]->left;
        right = pr->partitionData[model]->right;

        /* if we use per-partition branch length optimization 
           get the branch length of partition model and take the log otherwise 
           use the joint branch length among all partitions that is always stored 
           at index [0] */

        if(pr->perGeneBranchLengths)
        {
          qz = tInfo->qz[model];                                    
          rz = tInfo->rz[model];                  
        }
        else
        {
          qz = tInfo->qz[0];
          rz = tInfo->rz[0];
        }

        qz = (qz > PLL_ZMIN) ? log(qz) : log(PLL_ZMIN);                        
        rz = (rz > PLL_ZMIN) ? log(rz) : log(PLL_ZMIN);                       

        /* compute the left and right P matrices */

        if(pr->partitionData[model]->dataType == PLL_AA_DATA &&
        		(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X))
                makeP_FlexLG4(qz, rz, pr->partitionData[model]->gammaRates,
                              pr->partitionData[model]->EI_LG4,
                              pr->partitionData[model]->EIGN_LG4,
                              4, left, right, 20);
        else
        makeP(qz, rz, rateCategories,   pr->partitionData[model]->EI,
              pr->partitionData[model]->EIGN, categories,
              left, right, tr->saveMemory, tr->maxCategories, states);


#if (!defined(__SSE3) && !defined(__AVX) && !defined(__MIC_NATIVE))
        assert(!tr->saveMemory);

        /* figure out if we need to compute the CAT or GAMMA model of rate heterogeneity */

        if(tr->rateHetModel == PLL_CAT)
         {

           newviewCAT_FLEX(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                           x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                           ex3, tipX1, tipX2,
                           width, left, right, wgt, &scalerIncrement, fastScaling, states);
         }
        else 
         {
            newviewGAMMA_FLEX(tInfo->tipCase,
                 x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                 0, tipX1, tipX2,
                 width, left, right, wgt, &scalerIncrement, fastScaling, states, getUndetermined(pr->partitionData[model]->dataType) + 1);
         }
#else
        /* dedicated highly optimized functions. Analogously to the functions in evaluateGeneric() 
           we also siwtch over the state number */

        switch(states)
        {               
        case 2:
          assert (!tr->saveMemory);
          if (tr->rateHetModel == PLL_CAT)
           {
             newviewGTRCAT_BINARY(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                  x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                  ex3, tipX1, tipX2,
                                  width, left, right, wgt, &scalerIncrement, fastScaling);
           }
          else
           {
             newviewGTRGAMMA_BINARY(tInfo->tipCase,
                                    x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                    ex3, tipX1, tipX2,
                                    width, left, right, wgt, &scalerIncrement, fastScaling);                  
           }
          break;

        case 4: /* DNA */
#ifdef __MIC_NATIVE

              /* CAT & memory saving are not supported on MIC */

              assert(!tr->saveMemory);
              assert(tr->rateHetModel == PLL_GAMMA);

              newviewGTRGAMMA_MIC(tInfo->tipCase,
                                x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                ex3, tipX1, tipX2,
                                width, left, right, wgt, &scalerIncrement, fastScaling);
#else
          if(tr->rateHetModel == PLL_CAT)
            {                                
              
              if(tr->saveMemory)
#ifdef __AVX
                newviewGTRCAT_AVX_GAPPED_SAVE(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                              x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                              ex3, tipX1, tipX2,
                                              width, left, right, wgt, &scalerIncrement, fastScaling, x1_gap, x2_gap, x3_gap,
                                              x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#else
                newviewGTRCAT_SAVE(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                   x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                   ex3, tipX1, tipX2,
                                   width, left, right, wgt, &scalerIncrement, fastScaling, x1_gap, x2_gap, x3_gap,
                                   x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
              else
#ifdef __AVX
                newviewGTRCAT_AVX(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                  x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                  ex3, tipX1, tipX2,
                                  width, left, right, wgt, &scalerIncrement, fastScaling);
#else
              newviewGTRCAT(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                            x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                            ex3, tipX1, tipX2,
                            width, left, right, wgt, &scalerIncrement, fastScaling);
#endif
            }
          else
            {
              
              if(tr->saveMemory)
#ifdef __AVX
                newviewGTRGAMMA_AVX_GAPPED_SAVE(tInfo->tipCase,
                                                x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                                ex3, tipX1, tipX2,
                                                width, left, right, wgt, &scalerIncrement, fastScaling,
                                                x1_gap, x2_gap, x3_gap, 
                                                x1_gapColumn, x2_gapColumn, x3_gapColumn);

#else
              newviewGTRGAMMA_GAPPED_SAVE(tInfo->tipCase,
                                          x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                          ex3, tipX1, tipX2,
                                          width, left, right, wgt, &scalerIncrement, fastScaling,
                                          x1_gap, x2_gap, x3_gap, 
                                          x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
              else
#ifdef __AVX
                newviewGTRGAMMA_AVX(tInfo->tipCase,
                                    x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                    ex3, tipX1, tipX2,
                                    width, left, right, wgt, &scalerIncrement, fastScaling);
#else
              newviewGTRGAMMA(tInfo->tipCase,
                              x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                              ex3,tipX1, tipX2,
                              width, left, right, wgt, &scalerIncrement, fastScaling);
#endif
            }
#endif

            break;                  
          case 20: /* proteins */

#ifdef __MIC_NATIVE

                        /* CAT & memory saving are not supported on MIC */

                        assert(!tr->saveMemory);
                        assert(tr->rateHetModel == PLL_GAMMA);

                        if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                        {
                                  newviewGTRGAMMAPROT_LG4_MIC(tInfo->tipCase,
                            x1_start, x2_start, x3_start, pr->partitionData[model]->EV_LG4, pr->partitionData[model]->tipVector_LG4,
                            tipX1, tipX2,
                            width, left, right, wgt, &scalerIncrement);
                        }
                        else
                        {
                                  newviewGTRGAMMAPROT_MIC(tInfo->tipCase,
                                                x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                                ex3, tipX1, tipX2,
                                                width, left, right, wgt, &scalerIncrement, fastScaling);
                        }
#else

            if(tr->rateHetModel == PLL_CAT)
            {


              if(tr->saveMemory)
#ifdef __AVX
                newviewGTRCATPROT_AVX_GAPPED_SAVE(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                                  x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                                  ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, fastScaling, 
                                                  x1_gap, x2_gap, x3_gap,
                                                  x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#else
              newviewGTRCATPROT_SAVE(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                     x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                     ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, fastScaling, x1_gap, x2_gap, x3_gap,
                                     x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
#endif
              else
#ifdef __AVX
                newviewGTRCATPROT_AVX(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                      x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                      ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, fastScaling);
#else
              newviewGTRCATPROT(tInfo->tipCase,  pr->partitionData[model]->EV, pr->partitionData[model]->rateCategory,
                                x1_start, x2_start, x3_start, pr->partitionData[model]->tipVector,
                                ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, fastScaling);                     
#endif
            }
            else
            {

              

              if(tr->saveMemory)
#ifdef __AVX
                newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(tInfo->tipCase,
                                                    x1_start, x2_start, x3_start,
                                                    pr->partitionData[model]->EV,
                                                    pr->partitionData[model]->tipVector,
                                                    ex3, tipX1, tipX2,
                                                    width, left, right, wgt, &scalerIncrement, fastScaling,
                                                    x1_gap, x2_gap, x3_gap,
                                                    x1_gapColumn, x2_gapColumn, x3_gapColumn);
#else
                newviewGTRGAMMAPROT_GAPPED_SAVE(tInfo->tipCase,
                                                x1_start, x2_start, x3_start,
                                                pr->partitionData[model]->EV,
                                                pr->partitionData[model]->tipVector,
                                                ex3, tipX1, tipX2,
                                                width, left, right, wgt, &scalerIncrement, fastScaling,
                                                x1_gap, x2_gap, x3_gap,
                                                x1_gapColumn, x2_gapColumn, x3_gapColumn);
#endif
            
             else
                        {
                          if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
                            {
#ifdef __AVX 
                              newviewGTRGAMMAPROT_AVX_LG4(tInfo->tipCase,
                                                          x1_start, x2_start, x3_start,
                                                          pr->partitionData[model]->EV_LG4,
                                                          pr->partitionData[model]->tipVector_LG4,
                                                          (int*)NULL, tipX1, tipX2,
                                                          width, left, right, wgt, &scalerIncrement, PLL_TRUE);
#else
                              newviewGTRGAMMAPROT_LG4(tInfo->tipCase,
                                                      x1_start, x2_start, x3_start,
                                                      pr->partitionData[model]->EV_LG4,
                                                      pr->partitionData[model]->tipVector_LG4,
                                                      (int*)NULL, tipX1, tipX2,
                                                      width, left, right, 
                                                      wgt, &scalerIncrement, PLL_TRUE);
#endif                      
                            }
              else
#ifdef __AVX
                newviewGTRGAMMAPROT_AVX(tInfo->tipCase,
                                        x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                        ex3, tipX1, tipX2,
                                        width, left, right, wgt, &scalerIncrement, fastScaling);
#else
              newviewGTRGAMMAPROT(tInfo->tipCase,
                                  x1_start, x2_start, x3_start, pr->partitionData[model]->EV, pr->partitionData[model]->tipVector,
                                  ex3, tipX1, tipX2,
                                  width, left, right, wgt, &scalerIncrement, fastScaling);
#endif                 
            }   
        }
#endif
            
            break;      
          default:
            assert(0);
        }
#endif


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
       if(pr->partitionData[model]->ascBias && tr->threadID == 0)
#else
       if(pr->partitionData[model]->ascBias)
#endif         
         {
           switch(tr->rateHetModel)
             {
             case PLL_CAT:
               {
                 double 
                   rates = 1.0;
                 
                 //need to re-calculate transition probabilities assuming a rate of 1.0 
                 makeP(qz, rz, 
                       &rates,  
                       pr->partitionData[model]->EI,
                       pr->partitionData[model]->EIGN,
                       1, 
                       left, right, 
                       tr->saveMemory,
                       tr->maxCategories,
                       states);
                 
                 newviewAscCat(tInfo->tipCase,
                               x1_ascColumn, x2_ascColumn, x3_ascColumn,
                               pr->partitionData[model]->EV,
                               pr->partitionData[model]->tipVector,
                               ex3_asc,
                               states, left, right, states);
               }
               break;
             case PLL_GAMMA:
               newviewAscGamma(tInfo->tipCase,
                               x1_ascColumn, x2_ascColumn, x3_ascColumn,
                               pr->partitionData[model]->EV,
                               pr->partitionData[model]->tipVector,
                               ex3_asc,
                               states, left, right, states);                        
               break;
             default:
               assert(0);
             }
         }


        /* important step, here we essentiallt recursively compute the number of scaling multiplications 
           at node p: it's the sum of the number of scaling multiplications already conducted 
           for computing nodes q and r plus the scaling multiplications done at node p */

        if(fastScaling)
          {
            pr->partitionData[model]->globalScaler[tInfo->pNumber] =
              pr->partitionData[model]->globalScaler[tInfo->qNumber] +
              pr->partitionData[model]->globalScaler[tInfo->rNumber] +
              (unsigned int)scalerIncrement;
            
            /* check that we are not getting an integer overflow ! */

            assert(pr->partitionData[model]->globalScaler[tInfo->pNumber] < INT_MAX);
          }
        
        /* show the output vector */
      } 
    }
  }
}

/** @brief Compute the traversal descriptor of the subtree rooted at \a p.
    
    Computes the traversal descriptor of the subtree with root \a p. By traversal
    descriptory we essentially mean a preorder traversal of the unrooted topology
    by rooting it at a node \a p.
    If \a partialTraversal is set to \b PLL_TRUE then subtrees which are oriented
    correctly (i.e. if root node \a r of a subtree has \a r->x == 1) are not
    included in the traversal descriptor.

    @param tr
      PLL instance

    @param p
      Node assumed to be the root

    @param partialTraversal
      If set to \b PLL_TRUE, then a partial traversal descriptor is computed.

    @param numBranches
      Number of branches (either per-partition branch or joint branch estimate)
*/
void computeTraversal(pllInstance *tr, nodeptr p, pllBoolean partialTraversal, int numBranches)
{
  /* Only if we apply recomputations we need the additional step of updating the subtree lengths */
  if(tr->useRecom)
  {
    int traversal_counter = 0;
    if(partialTraversal)
      computeTraversalInfoStlen(p, tr->mxtips, tr->rvec, &traversal_counter);
    else
      computeFullTraversalInfoStlen(p, tr->mxtips, tr->rvec);
  }
  computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, numBranches, partialTraversal, tr->rvec, tr->useRecom);
}


/** @brief Computes the conditional likelihood vectors of all nodes in the subtree rooted at \a p
  
    Compute the conditional likelihood vectors of all nodes in the subtree rooted at node \a p. The
    conditional likelihood vector at node \a p is recomputed regardless of whether the orientation (i.e. \a p->x) 
    is correct or not, and, recursuvely, the likelihoods at each node in the subtree as needed and if necessary.
    In case \a masked is set to \b PLL_TRUE, the computation will not take place at partitions for which the 
    conditional likelihood has converged (for example as a reult of previous branch length optimization).
    
    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Root of the subtree for which we want to recompute the conditional likelihood vectors

    @param masked
      If set to \b PLL_TRUE, then likelihood vectors of partitions that are converged are
      not recomputed.
 */
void pllUpdatePartials (pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean masked)
{  
  /* if it's a tip there is nothing to do */

  if(isTip(p->number, tr->mxtips))
    return;

  /* the first entry of the traversal descriptor is always reserved for evaluate or branch length optimization calls,
     hence we start filling the array at the second entry with index one. This is not very nice and should be fixed 
     at some point */

  tr->td[0].count = 0;

  /* compute the traversal descriptor, which will include nodes-that-need-update descending the subtree  p */
  computeTraversal(tr, p, PLL_TRUE, pr->perGeneBranchLengths?pr->numberOfPartitions : 1);

  /* the traversal descriptor has been recomputed -> not sure if it really always changes, something to 
     optimize in the future */
  tr->td[0].traversalHasChanged = PLL_TRUE;

  /* We do a masked newview, i.e., do not execute newvies for each partition, when for example 
     doing a branch length optimization on the entire tree when branches are estimated on a per partition basis.

     you may imagine that for partition 5 the branch length optimization has already converged whereas 
     for partition 6 we still need to go over the tree again.

     This is explained in more detail in:

     A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009

     The external pllBoolean array tr->partitionConverged[] contains exactly that information and is copied
     to executeModel and subsequently to the executeMask of the traversal descriptor 

*/


  if(masked)
  {
    int model;

    for(model = 0; model < pr->numberOfPartitions; model++)
    {
      if(tr->partitionConverged[model])
        pr->partitionData[model]->executeModel = PLL_FALSE;
      else
        pr->partitionData[model]->executeModel = PLL_TRUE;
    }
  }

  /* if there is something to re-compute */

  if(tr->td[0].count > 0)
  {
    /* store execute mask in traversal descriptor */

    storeExecuteMaskInTraversalDescriptor(tr, pr);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
    /* do the parallel for join for pthreads
       not that we do not need a reduction operation here, but just a barrier to make 
       sure that all threads are done with their partition */

    pllMasterBarrier(tr, pr, PLL_THREAD_NEWVIEW);
#else
    /* in the sequential case we now simply call pllNewviewIterative() */

    pllNewviewIterative(tr, pr, 0);
#endif

  }

  /* clean up */

  if(masked)
  {
    int model;

    for(model = 0; model < pr->numberOfPartitions; model++)
      pr->partitionData[model]->executeModel = PLL_TRUE;
  }

  tr->td[0].traversalHasChanged = PLL_FALSE;
}

/* function to compute the marginal ancestral probability vector at a node p for CAT/PSR model */

/** @brief Compute the marginal ancestral probability vector for CAT/PSR model
    
    Computes the marginal ancestral probability vector for CAT/PSR model, given the conditional likelihood
    vector \a x3 of some node, and a zero branch length P matrix \a diagptable.

    @param x3
      Conditional likelihood of the node for which we are computing the ancestral vector

    @param ancestralBuffer
      Buffer where to store the marginal ancestral probability vector

    @param diagptable
      A zero branch length P matrix

    @param n
      Number of sites in the partition to process (in the case of MPI/PTHREADS, the number of sites in the partition assigned to the current thread/process)

    @param numStates
      Number of states

    @param cptr
      Array where the rate for each site in the compressed partition alignment is stored
      
 */
static void ancestralCat(double *x3, double *ancestralBuffer, double *diagptable, const int n, const int numStates, int *cptr)
{ 
  double 
    *term = (double*)rax_malloc(sizeof(double) * numStates);

  int 
    i;

  const int
    statesSquare = numStates * numStates;
  
  for(i = 0; i < n; i++)
    {
      double 
        sum = 0.0,
        *v = &x3[numStates * i],
        *ancestral = &ancestralBuffer[numStates * i],
        *d = &diagptable[cptr[i] * statesSquare];            

      int 
        l,
        j;

      for(l = 0; l < numStates; l++)
        {
          double 
            ump_x1 = 0.0;
      
          for(j = 0; j < numStates; j++)        
            ump_x1 += v[j] * d[l * numStates + j];

          sum += ump_x1;
          term[l] = ump_x1;      
        }
                
      for(l = 0; l < numStates; l++)          
        ancestral[l] = term[l] / sum;   
    }
   
  rax_free(term);
}


/* compute marginal ancestral states for GAMMA models,
   for the euqation to obtain marginal ancestral states 
   see Ziheng Yang's book */

/** @brief Compute the marginal ancestral probability vector for GAMMA model
    
    Computes the marginal ancestral probability vector for the GAMMA model, given the conditional likelihood
    vector \a x3 of some node, and a zero branch length P matrix \a diagptable.

    @param x3
      Conditional likelihood of the node for which we are computing the ancestral vector

    @param ancestralBuffer
      Buffer where to store the marginal ancestral probability vector

    @param diagptable
      A zero branch length P matrix

    @param n
      Number of sites in the partition to process (in the case of MPI/PTHREADS, the number of sites in the partition assigned to the current thread/process)

    @param numStates
      Number of states

    @param gammaStates
      Number of GAMMA categories times number of states
      
 */
static void ancestralGamma(double *x3, double *ancestralBuffer, double *diagptable, const int n, const int numStates, const int gammaStates)
{
  int 
    i;

  const int
    statesSquare = numStates * numStates;

  double    
    *term = (double*)rax_malloc(sizeof(double) * numStates);                  
  
  for(i = 0; i < n; i++)
    {
      double 
        sum = 0.0,
        *_v = &x3[gammaStates * i],
        *ancestral = &ancestralBuffer[numStates * i];  
      
      int
        k,
        j,
        l;
      
      for(l = 0; l < numStates; l++)
        term[l] = 0.0;

      for(k = 0; k < 4; k++)
        {
          double 
            *v =  &(_v[numStates * k]);

          for(l = 0; l < numStates; l++)
            {
              double
                al = 0.0;
              
              for(j = 0; j < numStates; j++)        
                al += v[j] * diagptable[k * statesSquare + l * numStates + j];
          
              term[l] += al;
              sum += al;
            }
        }
  
      for(l = 0; l < numStates; l++)        
        ancestral[l] = term[l] / sum;       
    }
   
  rax_free(term);
}

/* compute dedicated zero branch length P matrix */
/** @brief Compute a dedicated zero branch length P matrix
   
    Computes a P matrix by assuming a branch length of zero. This is used
    for the marginal ancestral probabilities recomputation.

    @param rptr
      Array of values for rate categories

    @param EI
      Inverse eigenvector of Q matrix

    @param EIGN
      Eigenvalues of Q matrix

    @param numberOfCategories
      Number of rate categories

    @param left
      Where to store the resulting P matrix

    @param numStates
      Number of states
 */
static void calc_diagp_Ancestral(double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, const int numStates)
{
  int 
    i,
    j,
    k;
  
  const int   
    statesSquare = numStates * numStates;

  double 
    z1 = 0.0,
    lz1[64],
    d1[64];

  assert(numStates <= 64);
     
  for(i = 0; i < numStates; i++)    
    lz1[i] = EIGN[i] * z1;
     

  for(i = 0; i < numberOfCategories; i++)
    {
      d1[0] = 1.0;

      for(j = 1; j < numStates; j++)    
        d1[j] = exp(rptr[i] * lz1[j]);
         
      for(j = 0; j < numStates; j++)
        {
          left[statesSquare * i  + numStates * j] = 1.0;         

          for(k = 1; k < numStates; k++)            
            left[statesSquare * i + numStates * j + k]  = d1[k] * EI[numStates * j + k];             
        }
    }  
}

/** @brief A very simple iterative function, we only access the conditional likelihood vector at node \a p
 *
 *
 */
void newviewAncestralIterative(pllInstance *tr, partitionList *pr)
{
  traversalInfo 
    *ti    = tr->td[0].ti,
    *tInfo = &ti[0];

  int    
    model,
    p_slot = -1;

  /* make sure that the traversal descriptor has length 1 */

  assert(tr->td[0].count == 1);
  assert(!tr->saveMemory);

  /* get the index to the conditional likelihood vector depending on whether recomputation is used or not */

  if(tr->useRecom)    
    p_slot = tInfo->slot_p;         
  else    
    p_slot = tInfo->pNumber - tr->mxtips - 1;         

  /* now loop over all partitions for nodes p of the current traversal vector entry */

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      /* number of sites in this partition */
      size_t            
        width  = (size_t)pr->partitionData[model]->width;

      /* this conditional statement is exactly identical to what we do in pllEvaluateIterative */

      if(tr->td[0].executeModel[model] && width > 0)
        {             
          double         
            *x3_start = pr->partitionData[model]->xVector[p_slot],
//          *left     = (double*)NULL,
//          *right    = (double*)NULL,                 
            *rateCategories = (double*)NULL,
            *diagptable = (double*)NULL;

          int
            categories;
        
          size_t                  
            states = (size_t)pr->partitionData[model]->states,
            availableLength = pr->partitionData[model]->xSpaceVector[p_slot],
            requiredLength = 0,
            rateHet = discreteRateCategories(tr->rateHetModel);   

        /* figure out what kind of rate heterogeneity approach we are using */

          if(tr->rateHetModel == PLL_CAT)
            {            
              rateCategories = pr->partitionData[model]->perSiteRates;
              categories     = pr->partitionData[model]->numberOfCategories;
            }
          else
            {                            
              rateCategories = pr->partitionData[model]->gammaRates;
              categories     = 4;
            }
          
          /* allocate some space for a special P matrix with a branch length of 0 into which we mingle 
             the eignevalues. This will allow us to obtain real probabilites from the internal RAxML 
             representation */

          rax_posix_memalign ((void **)&diagptable, PLL_BYTE_ALIGNMENT, categories * states * states * sizeof(double));
          
          requiredLength  =  virtual_width( width ) * rateHet * states * sizeof(double);
          
          /* make sure that this vector had already been allocated. This must be PLL_TRUE since we first invoked a standard newview() on this */

          assert(requiredLength == availableLength);                                     

          /* now compute the special P matrix */

          calc_diagp_Ancestral(rateCategories, pr->partitionData[model]->EI,  pr->partitionData[model]->EIGN, categories, diagptable, states);
          
          /* switch over the rate heterogeneity model 
             and call generic functions that compute the marginal ancestral states and 
             store them in pr->partitionData[model]->ancestralBuffer
          */

          if(tr->rateHetModel == PLL_CAT)       
            ancestralCat(x3_start, pr->partitionData[model]->ancestralBuffer, diagptable, width, states, pr->partitionData[model]->rateCategory);
          else
            ancestralGamma(x3_start, pr->partitionData[model]->ancestralBuffer, diagptable, width, states, categories * states);
          
          rax_free(diagptable);                   
        }       
    }
}

/** @brief Computes the Conditional Likelihood Vector (CLV) for each rate of some internal node.

    Computes the conditional likelihood vectors of node \a p for each rate, given the partition
    index \a partition. The result is placed in the array \a outProbs, which must be pre-allocated
    by the caller, and must be of size \a sites * categories * states * sizeof(double). The structure of
    the resulting array is the following:
    For each site we have \a categories * states cells of size \a double. Those cells are divided per rate
    category, i.e. first \a states cells are the probabilities for the states of rate 1 (ordered alphabetically
    by base name), next \a states cells for rate 2 and so on.

    @param tr   PLL instance
    @param pr     List of partitions
    @param p Node for which we want to compute the CLV
    @param partition   Index of the partition for which to compute the CLV
    @param outProbs    Pre-allocated array where the result will be stored

    @returns Returns \b PLL_TRUE on success, \b PLL_FALSE on failure

    @todo       Fix to work with CAT
*/
int pllGetCLV (pllInstance * tr, partitionList * pr, nodeptr p, int partition, double * outProbs)
{
  size_t i, j, k, l;

  if (tr->rateHetModel != PLL_GAMMA) return (PLL_FALSE);

  int p_slot;
  size_t states = (size_t)pr->partitionData[partition]->states;

  double
    *term = (double*)rax_malloc(sizeof(double) * states);

  if(tr->useRecom)
    p_slot = p->number;
  else
    p_slot = p->number - tr->mxtips - 1;

  size_t width = (size_t) pr->partitionData[partition]->width;
  double * diagptable = NULL;
  double * rateCategories = pr->partitionData[partition]->gammaRates;
  double * x3 = pr->partitionData[partition]->xVector[p_slot];
  size_t categories = 4;

  rax_posix_memalign ((void **)&diagptable, PLL_BYTE_ALIGNMENT, categories * states * states * sizeof (double));

  calc_diagp_Ancestral(rateCategories, pr->partitionData[partition]->EI,  pr->partitionData[partition]->EIGN, categories, diagptable, states);

  for (i = 0; i < width; ++ i)
   {
     double
       *_v  = &x3[categories * states * i],
       *clv = &outProbs[categories * states * i];

     for (k = 0; k < categories; ++ k)
      {
        double
         sum = 0.0,
         *v = &(_v[states * k]);

        for (l = 0; l < states; ++ l)
         {
           double al = 0.0;

           for (j = 0; j < states; ++ j)
             al += v[j] * diagptable[k * states * states + l * states + j];

           term[l] = al;
           sum += al;
         }
        for (l = 0; l < states; ++ l)
           clv[k * categories + l] = term[l] / sum;
      }
   }

  rax_free(term);
  rax_free(diagptable);

  return (PLL_TRUE);
}

/* this is very similar to pllUpdatePartials, except that it also computes the marginal ancestral probabilities 
   at node p. To simplify the code I am re-using newview() here to first get the likelihood vector p->x at p
   and then I deploy newviewAncestralIterative(tr); that should always only have a traversal descriptor of lenth 1,
   to do some mathematical transformations that are required to obtain the marginal ancestral probabilities from 
   the conditional likelihood array at p.

   Note that the marginal ancestral probability vector summarizes the subtree rooted at p! */

/** @brief Computes the conditional likelihood vectors of all nodes in the subtree rooted at \a p
    and the marginal ancestral probabilities at node \a p

    Compute the conditional likelihood vectors of all nodes in the subtree rooted at node \a p. The
    conditional likelihood vector at node \a p is recomputed regardless of whether the orientation (i.e. \a p->x)
    is correct or not, and, recursively, the likelihoods at each node in the subtree as needed and if necessary.
    In addition, the marginal ancestral probability vector for node \a p is also computed.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Node for which we want to compute the ancestral vector

    @note
      This function is not implemented with the saveMemory technique. 
*/
void pllUpdatePartialsAncestral(pllInstance *tr, partitionList *pr, nodeptr p)
{
  /* error check, we don't need to compute anything for tips */
  
  if(isTip(p->number, tr->mxtips))
    {
      printf("You are trying to compute the ancestral states on a tip node of the tree\n");
      assert(0);
    }

  /* doesn't work yet in conjunction with SEVs, can be implemented though at some point 
     if urgently required */

  if(tr->saveMemory)
    {
      printf("ancestral state implementation will not work with memory saving (SEVs) enabled!\n");
      printf("returning without computing anything ... \n");
      return;
    }

  /* first call pllUpdatePartials() with mask set to PLL_FALSE such that the likelihood vector is there ! */

  pllUpdatePartials(tr, pr, p, PLL_FALSE);

  /* now let's compute the ancestral states using this vector ! */
  
  /* to make things easy and reduce code size, let's re-compute a standard traversal descriptor for node p,
     hence we need to set the count to 0 */

  tr->td[0].count = 0;

  computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, pr->perGeneBranchLengths?pr->numberOfPartitions : 1, PLL_TRUE, tr->rvec, tr->useRecom);

  tr->td[0].traversalHasChanged = PLL_TRUE;

  /* here we actually assert, that the traversal descriptor only contains one node triplet p, p->next->back, p->next->next->back
     this must be PLL_TRUE because we have alread invoked the standard pllUpdatePartials() on p.
  */ 

  assert(tr->td[0].count == 1);  
  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* use the pthreads barrier to invoke newviewAncestralIterative() on a per-thread basis */

  pllMasterBarrier (tr, pr, PLL_THREAD_NEWVIEW_ANCESTRAL);
#else
  /* now call the dedicated function that does the mathematical transformation of the 
     conditional likelihood vector at p to obtain the marginal ancestral states */

  newviewAncestralIterative(tr, pr);
#endif

  tr->td[0].traversalHasChanged = PLL_FALSE;

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* invoke another parallel region to gather the marginal ancestral probabilities 
     from the threads/MPI processes */

  pllMasterBarrier (tr, pr, PLL_THREAD_GATHER_ANCESTRAL);
#endif

  
}

/* returns the character representation of an enumerated DNA or AA state */

/** @brief Get the character representation of an enumerated DNA or AA state
    
    Returns the character representation of the enumarates DNA or AA state,
    from the constant arrays \a dnaStateNames (for DNA) or \a protStateNames (for proteins).

    @param dataType
      Type of data, i.e. \b PLL_DNA_DATA or \b PLL_AA_DATA

    @param state
      The number which we want to decode to a letter

    @return
      Returns the decoded character
 */
static char getStateCharacter(int dataType, int state)
{
  char 
    result;

  switch(dataType)
    {    
    case PLL_BINARY_DATA:
       result = binaryStateNames[state];
       break;
    case PLL_DNA_DATA:
       result = dnaStateNames[state];
      break;
    case PLL_AA_DATA:
      result =  protStateNames[state];
      break;    
    default:
      assert(0);
      result = '\x0';
      break;
  }

  return  result;
}

/** @brief Prints the ancestral state information for a node \a p to the terminal 
 
    Prints the ancestral state information for a node \a p to the terminal. 
    The ancestral state sequence, resp. marginal ancestral state probabilities, is printed
    depending on whether \a \a printStates, resp. \a printProbs, is set to \b PLL_TRUE.

    @param p
      The node for which to print the ancestral state sequence

    @param printStates
      If set to \b PLL_TRUE then the ancestral state sequence is printed

    @param printProbs
      If set to \b PLL_TRUE then the marginal ancestral state probabilities are printed

    @param tr
      PLL instance

    @param pr
      List of partitions
 
    @note  Here one can see how to store the ancestral probabilities in a dedicated data structure
 */
void printAncestralState(nodeptr p, pllBoolean printStates, pllBoolean printProbs, pllInstance *tr, partitionList *pr)
{
#ifdef _USE_PTHREADS
  size_t 
    accumulatedOffset = 0;
#endif

  int
    j,
    k,
    model,
    globalIndex = 0;
  
  /* allocate an array of structs for storing ancestral prob vector info/data */

  ancestralState 
    *a = (ancestralState *)rax_malloc(sizeof(ancestralState) * tr->originalCrunchedLength);   

  /* loop over partitions */

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      int            
        i,
        width = pr->partitionData[model]->upper - pr->partitionData[model]->lower,
        states = pr->partitionData[model]->states;
      
      /* set pointer to ancestral probability vector */

#ifdef _USE_PTHREADS
      double
        *ancestral = &tr->ancestralVector[accumulatedOffset];
#else
      double 
        *ancestral = pr->partitionData[model]->ancestralBuffer;
#endif        
      
      /* loop over the sites of the partition */

      for(i = 0; i < width; i++, globalIndex++)
        {
          double
            equal = 1.0 / (double)states,
            max = -1.0;
            
          pllBoolean
            approximatelyEqual = PLL_TRUE;

          int
            max_l = -1,
            l;
          
          char 
            c;

          /* stiore number of states for this site */

          a[globalIndex].states = states;

          /* alloc space for storing marginal ancestral probabilities */

          a[globalIndex].probs = (double *)rax_malloc(sizeof(double) * states);
          
          /* loop over states to store probabilities and find the maximum */

          for(l = 0; l < states; l++)
            {
              double 
                value = ancestral[states * i + l];

              if(value > max)
                {
                  max = value;
                  max_l = l;
                }
              
              /* this is used for discretizing the ancestral state sequence, if all marginal ancestral 
                 probabilities are approximately equal we output a ? */

              approximatelyEqual = approximatelyEqual && (PLL_ABS(equal - value) < 0.000001);
              
              a[globalIndex].probs[l] = value;                
            }

          
          /* figure out the discrete ancestral nucleotide */

          if(approximatelyEqual)
            c = '?';      
          else
            c = getStateCharacter(pr->partitionData[model]->dataType, max_l);
          
          a[globalIndex].c = c;   
        }

#ifdef _USE_PTHREADS
      accumulatedOffset += width * states;
#endif            
    }

  /* print marginal ancestral probs to terminal */

  if(printProbs)
    {
      printf("%d\n", p->number);
      
      for(k = 0; k < tr->originalCrunchedLength; k++)
        {
          for(j = 0; j < a[k].states; j++)
            printf("%f ", a[k].probs[j]);
          printf("\n");      
        }
      
      printf("\n");
    }
 
  /* print discrete state ancestrakl sequence to terminal */

  if(printStates)
    {
      printf("%d ", p->number);

      for(k = 0; k < tr->originalCrunchedLength; k++)          
        printf("%c", a[k].c);   
  
      printf("\n");
    }
  
  /* free the ancestral state data structure */
          
  for(j = 0; j < tr->originalCrunchedLength; j++)
    rax_free(a[j].probs);  

  rax_free(a);
}

void pllGetAncestralState(pllInstance *tr, partitionList *pr, nodeptr p, double * outProbs, char * outSequence)
{
#ifdef _USE_PTHREADS
  size_t 
    accumulatedOffset = 0;
#endif

  int
    j,
    k,
    model,
    globalIndex = 0;
     
  pllUpdatePartialsAncestral(tr, pr, p);
  
  /* allocate an array of structs for storing ancestral prob vector info/data */

  ancestralState 
    *a = (ancestralState *)rax_malloc(sizeof(ancestralState) * tr->originalCrunchedLength);   

  /* loop over partitions */

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      int            
        i,
        width = pr->partitionData[model]->upper - pr->partitionData[model]->lower,
        states = pr->partitionData[model]->states;
      
      /* set pointer to ancestral probability vector */

#ifdef _USE_PTHREADS
      double
        *ancestral = &tr->ancestralVector[accumulatedOffset];
#else
      double 
        *ancestral = pr->partitionData[model]->ancestralBuffer;
#endif        
      
      /* loop over the sites of the partition */

      for(i = 0; i < width; i++, globalIndex++)
        {
          double
            equal = 1.0 / (double)states,
            max = -1.0;
            
          pllBoolean
            approximatelyEqual = PLL_TRUE;

          int
            max_l = -1,
            l;
          
          char 
            c;

          /* stiore number of states for this site */

          a[globalIndex].states = states;

          /* alloc space for storing marginal ancestral probabilities */

          a[globalIndex].probs = (double *)rax_malloc(sizeof(double) * states);
          
          /* loop over states to store probabilities and find the maximum */

          for(l = 0; l < states; l++)
            {
              double 
                value = ancestral[states * i + l];

              if(value > max)
                {
                  max = value;
                  max_l = l;
                }
              
              /* this is used for discretizing the ancestral state sequence, if all marginal ancestral 
                 probabilities are approximately equal we output a ? */

              approximatelyEqual = approximatelyEqual && (PLL_ABS(equal - value) < 0.000001);
              
              a[globalIndex].probs[l] = value;                
            }

          
          /* figure out the discrete ancestral nucleotide */

          if(approximatelyEqual)
            c = '?';      
          else
            c = getStateCharacter(pr->partitionData[model]->dataType, max_l);
          
          a[globalIndex].c = c;   
        }

#ifdef _USE_PTHREADS
      accumulatedOffset += width * states;
#endif            
    }

  /* print marginal ancestral probs to terminal */

  for(k = 0; k < tr->originalCrunchedLength; k++)
    {
      for(j = 0; j < a[k].states; j++)
        outProbs[k * a[k].states + j] = a[k].probs[j];
    }
 
  /* print discrete state ancestrakl sequence to terminal */

  for(k = 0; k < tr->originalCrunchedLength; k++)          
      outSequence[k] = a[k].c;
  outSequence[tr->originalCrunchedLength] = 0;
  
  /* free the ancestral state data structure */
          
  for(j = 0; j < tr->originalCrunchedLength; j++)
    rax_free(a[j].probs);  

  rax_free(a);
}
/* optimized function implementations */


/**
 *  @defgroup group1 Optimized functions
 *  This is the optimized functions group
 */

#if (!defined(__AVX) && defined(__SSE3))

/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR GAMMA with memory saving (Optimized SSE3 version for DNA data)

    This is the SSE3 optimized version of ::newviewGAMMA_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b GAMMA
    model of rate heterogeneity. The memory saving technique is incorporated.

    @note
    For more details and function argument description check the function ::newviewGAMMA_FLEX
*/
static void newviewGTRGAMMA_GAPPED_SAVE(int tipCase,
                                        double *x1_start, double *x2_start, double *x3_start,
                                        double *EV, double *tipVector,
                                        int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                        const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                        unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
                                        double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn)
{
  int     
    i, 
    j, 
    k, 
    l,
    addScale = 0, 
    scaleGap = 0;

  double
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start,       
    max;
  PLL_ALIGN_BEGIN double
    maxima[2] PLL_ALIGN_END,
    EV_t[16] PLL_ALIGN_END;

  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);      



  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        double *uX1, *uX2;
        PLL_ALIGN_BEGIN double umpX1[256] PLL_ALIGN_END, umpX2[256] PLL_ALIGN_END;


        for (i = 1; i < 16; i++)
        {           
          __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
          __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));       

          for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {                            
              __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
            }

          for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {
              __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);

            }
        }                 

        uX1 = &umpX1[240];
        uX2 = &umpX2[240];                          

        for (j = 0; j < 4; j++)
        {                                                                                  
          __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
          __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );

          __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
          __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );

          __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
          __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );                                                 

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6]; 
          __m128d EV_t_l3_k2 = EVV[7];

          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

          _mm_store_pd( &x3_gapColumn[j * 4 + 0], EV_t_l0_k0 );
          _mm_store_pd( &x3_gapColumn[j * 4 + 2], EV_t_l2_k0 );    
        }  


        x3 = x3_start;

        for (i = 0; i < n; i++)
        {           
          if(!(x3_gap[i / 32] & mask32[i % 32]))             
          {
            uX1 = &umpX1[16 * tipX1[i]];
            uX2 = &umpX2[16 * tipX2[i]];                                        

            for (j = 0; j < 4; j++)
            {                                                                              
              __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
              __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );


              __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
              __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );


              //
              // multiply left * right
              //

              __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
              __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );


              //
              // multiply with EV matrix (!?)
              //

              __m128d EV_t_l0_k0 = EVV[0];
              __m128d EV_t_l0_k2 = EVV[1];
              __m128d EV_t_l1_k0 = EVV[2];
              __m128d EV_t_l1_k2 = EVV[3];
              __m128d EV_t_l2_k0 = EVV[4];
              __m128d EV_t_l2_k2 = EVV[5];
              __m128d EV_t_l3_k0 = EVV[6]; 
              __m128d EV_t_l3_k2 = EVV[7];

              EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
              EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
              EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

              EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
              EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

              EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
              EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

              EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
              EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
              EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

              EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
              EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
              EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

              EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

              _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
              _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
            }

            x3 += 16;
          }
        }
      }
      break;
    case PLL_TIP_INNER:
      { 
        double 
          *uX1;
        PLL_ALIGN_BEGIN double
          umpX1[256] PLL_ALIGN_END;

        for (i = 1; i < 16; i++)
        {
          __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
          __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));       

          for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {            
              __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);                
            }
        }

        {
          __m128d maxv =_mm_setzero_pd();

          scaleGap = 0;

          x2 = x2_gapColumn;                     
          x3 = x3_gapColumn;

          uX1 = &umpX1[240];         

          for (j = 0; j < 4; j++)
          {                                
            double *x2_p = &x2[j*4];
            double *right_k0_p = &right[j*16];
            double *right_k1_p = &right[j*16 + 1*4];
            double *right_k2_p = &right[j*16 + 2*4];
            double *right_k3_p = &right[j*16 + 3*4];
            __m128d x2_0 = _mm_load_pd( &x2_p[0] );
            __m128d x2_2 = _mm_load_pd( &x2_p[2] );

            __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
            __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
            __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
            __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
            __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
            __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
            __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
            __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );

            right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
            right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

            right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
            right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
            right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

            right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
            right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

            right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
            right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
            right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

            __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
            __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );

            __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
            __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );

            __m128d EV_t_l0_k0 = EVV[0];
            __m128d EV_t_l0_k2 = EVV[1];
            __m128d EV_t_l1_k0 = EVV[2];
            __m128d EV_t_l1_k2 = EVV[3];
            __m128d EV_t_l2_k0 = EVV[4];
            __m128d EV_t_l2_k2 = EVV[5];
            __m128d EV_t_l3_k0 = EVV[6]; 
            __m128d EV_t_l3_k2 = EVV[7];

            EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
            EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

            EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
            EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

            EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

            EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
            EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

            EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
            EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
            EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

            values[j * 2]     = EV_t_l0_k0;
            values[j * 2 + 1] = EV_t_l2_k0;                                

            maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
            maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));                                    
          }


          _mm_store_pd(maxima, maxv);

          max = PLL_MAX(maxima[0], maxima[1]);

          if(max < PLL_MINLIKELIHOOD)
          {
            scaleGap = 1;

            __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

            _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));       
            _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
            _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
            _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
            _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));       
            _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
            _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
            _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));                        
          }
          else
          {
            _mm_store_pd(&x3[0], values[0]);       
            _mm_store_pd(&x3[2], values[1]);
            _mm_store_pd(&x3[4], values[2]);
            _mm_store_pd(&x3[6], values[3]);
            _mm_store_pd(&x3[8], values[4]);       
            _mm_store_pd(&x3[10], values[5]);
            _mm_store_pd(&x3[12], values[6]);
            _mm_store_pd(&x3[14], values[7]);
          }
        }                       

        x3 = x3_start;

        for (i = 0; i < n; i++)
        {
          if((x3_gap[i / 32] & mask32[i % 32]))
          {            
            if(scaleGap)
            {   
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];                  
            }
          }
          else
          {                              
            __m128d maxv =_mm_setzero_pd();              

            if(x2_gap[i / 32] & mask32[i % 32])
              x2 = x2_gapColumn;
            else
            {
              x2 = x2_ptr;
              x2_ptr += 16;
            }

            uX1 = &umpX1[16 * tipX1[i]];             


            for (j = 0; j < 4; j++)
            {                              
              double *x2_p = &x2[j*4];
              double *right_k0_p = &right[j*16];
              double *right_k1_p = &right[j*16 + 1*4];
              double *right_k2_p = &right[j*16 + 2*4];
              double *right_k3_p = &right[j*16 + 3*4];
              __m128d x2_0 = _mm_load_pd( &x2_p[0] );
              __m128d x2_2 = _mm_load_pd( &x2_p[2] );

              __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
              __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
              __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
              __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
              __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
              __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
              __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
              __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );


              right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
              right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

              right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
              right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

              right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
              right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
              right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);


              right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
              right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


              right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
              right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

              right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
              right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
              right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

              {
                //
                // load left side from tip vector
                //

                __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
                __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );


                //
                // multiply left * right
                //

                __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
                __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );


                //
                // multiply with EV matrix (!?)
                //                                

                __m128d EV_t_l0_k0 = EVV[0];
                __m128d EV_t_l0_k2 = EVV[1];
                __m128d EV_t_l1_k0 = EVV[2];
                __m128d EV_t_l1_k2 = EVV[3];
                __m128d EV_t_l2_k0 = EVV[4];
                __m128d EV_t_l2_k2 = EVV[5];
                __m128d EV_t_l3_k0 = EVV[6]; 
                __m128d EV_t_l3_k2 = EVV[7];


                EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
                EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
                EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

                EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
                EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

                EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
                EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

                EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
                EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
                EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

                EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
                EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
                EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

                EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

                values[j * 2]     = EV_t_l0_k0;
                values[j * 2 + 1] = EV_t_l2_k0;                            

                maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
                maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));                
              }            
            }


            _mm_store_pd(maxima, maxv);

            max = PLL_MAX(maxima[0], maxima[1]);

            if(max < PLL_MINLIKELIHOOD)
            {
              __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

              _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));     
              _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
              _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
              _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
              _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));     
              _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
              _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
              _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));      

              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];

            }
            else
            {
              _mm_store_pd(&x3[0], values[0]);     
              _mm_store_pd(&x3[2], values[1]);
              _mm_store_pd(&x3[4], values[2]);
              _mm_store_pd(&x3[6], values[3]);
              _mm_store_pd(&x3[8], values[4]);     
              _mm_store_pd(&x3[10], values[5]);
              _mm_store_pd(&x3[12], values[6]);
              _mm_store_pd(&x3[14], values[7]);
            }            

            x3 += 16;
          }
        }
      }
      break;
    case PLL_INNER_INNER:         
      {
        __m128d maxv =_mm_setzero_pd();

        scaleGap = 0;

        x1 = x1_gapColumn;                  
        x2 = x2_gapColumn;          
        x3 = x3_gapColumn;

        for (j = 0; j < 4; j++)
        {

          double *x1_p = &x1[j*4];
          double *left_k0_p = &left[j*16];
          double *left_k1_p = &left[j*16 + 1*4];
          double *left_k2_p = &left[j*16 + 2*4];
          double *left_k3_p = &left[j*16 + 3*4];

          __m128d x1_0 = _mm_load_pd( &x1_p[0] );
          __m128d x1_2 = _mm_load_pd( &x1_p[2] );

          __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
          __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
          __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
          __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
          __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
          __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
          __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
          __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);


          double *x2_p = &x2[j*4];
          double *right_k0_p = &right[j*16];
          double *right_k1_p = &right[j*16 + 1*4];
          double *right_k2_p = &right[j*16 + 2*4];
          double *right_k3_p = &right[j*16 + 3*4];
          __m128d x2_0 = _mm_load_pd( &x2_p[0] );
          __m128d x2_2 = _mm_load_pd( &x2_p[2] );

          __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
          __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
          __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
          __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
          __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
          __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
          __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
          __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);                                    

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );                                          

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6]; 
          __m128d EV_t_l3_k2 = EVV[7];

          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );


          values[j * 2] = EV_t_l0_k0;
          values[j * 2 + 1] = EV_t_l2_k0;                           

          maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
          maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
        }

        _mm_store_pd(maxima, maxv);

        max = PLL_MAX(maxima[0], maxima[1]);

        if(max < PLL_MINLIKELIHOOD)
        {
          __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

          scaleGap = 1;

          _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));         
          _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
          _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
          _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
          _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));         
          _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
          _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
          _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));                      
        }
        else
        {
          _mm_store_pd(&x3[0], values[0]);         
          _mm_store_pd(&x3[2], values[1]);
          _mm_store_pd(&x3[4], values[2]);
          _mm_store_pd(&x3[6], values[3]);
          _mm_store_pd(&x3[8], values[4]);         
          _mm_store_pd(&x3[10], values[5]);
          _mm_store_pd(&x3[12], values[6]);
          _mm_store_pd(&x3[14], values[7]);
        }
      }


      x3 = x3_start;

      for (i = 0; i < n; i++)
      { 
        if(x3_gap[i / 32] & mask32[i % 32])
        {            
          if(scaleGap)
          {     
            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];                              
          }
        }
        else
        {
          __m128d maxv =_mm_setzero_pd();                   

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


          for (j = 0; j < 4; j++)
          {

            double *x1_p = &x1[j*4];
            double *left_k0_p = &left[j*16];
            double *left_k1_p = &left[j*16 + 1*4];
            double *left_k2_p = &left[j*16 + 2*4];
            double *left_k3_p = &left[j*16 + 3*4];

            __m128d x1_0 = _mm_load_pd( &x1_p[0] );
            __m128d x1_2 = _mm_load_pd( &x1_p[2] );

            __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
            __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
            __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
            __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
            __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
            __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
            __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
            __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );

            left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
            left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

            left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
            left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

            left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
            left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
            left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

            left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
            left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

            left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
            left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

            left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
            left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
            left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);


            //
            // multiply/add right side
            //
            double *x2_p = &x2[j*4];
            double *right_k0_p = &right[j*16];
            double *right_k1_p = &right[j*16 + 1*4];
            double *right_k2_p = &right[j*16 + 2*4];
            double *right_k3_p = &right[j*16 + 3*4];
            __m128d x2_0 = _mm_load_pd( &x2_p[0] );
            __m128d x2_2 = _mm_load_pd( &x2_p[2] );

            __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
            __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
            __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
            __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
            __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
            __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
            __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
            __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );

            right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
            right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

            right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
            right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
            right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

            right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
            right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


            right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
            right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
            right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);     

            //
            // multiply left * right
            //

            __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
            __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );


            //
            // multiply with EV matrix (!?)
            //       

            __m128d EV_t_l0_k0 = EVV[0];
            __m128d EV_t_l0_k2 = EVV[1];
            __m128d EV_t_l1_k0 = EVV[2];
            __m128d EV_t_l1_k2 = EVV[3];
            __m128d EV_t_l2_k0 = EVV[4];
            __m128d EV_t_l2_k2 = EVV[5];
            __m128d EV_t_l3_k0 = EVV[6]; 
            __m128d EV_t_l3_k2 = EVV[7];


            EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
            EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

            EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
            EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

            EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

            EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
            EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );


            EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
            EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
            EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );


            values[j * 2] = EV_t_l0_k0;
            values[j * 2 + 1] = EV_t_l2_k0;                         

            maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
            maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
          }


          _mm_store_pd(maxima, maxv);

          max = PLL_MAX(maxima[0], maxima[1]);

          if(max < PLL_MINLIKELIHOOD)
          {
            __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

            _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));       
            _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
            _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
            _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
            _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));       
            _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
            _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
            _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));        

            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];

          }
          else
          {
            _mm_store_pd(&x3[0], values[0]);       
            _mm_store_pd(&x3[2], values[1]);
            _mm_store_pd(&x3[4], values[2]);
            _mm_store_pd(&x3[6], values[3]);
            _mm_store_pd(&x3[8], values[4]);       
            _mm_store_pd(&x3[10], values[5]);
            _mm_store_pd(&x3[12], values[6]);
            _mm_store_pd(&x3[14], values[7]);
          }      



          x3 += 16;

        }
      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}



/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR GAMMA (Optimized SSE3 version for DNA data)

    This is the SSE3 optimized version of ::newviewGAMMA_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b GAMMA
    model of rate heterogeneity.

    @note
    For more details and function argument description check the function ::newviewGAMMA_FLEX
*/
static void newviewGTRGAMMA(int tipCase,
                            double *x1_start, double *x2_start, double *x3_start,
                            double *EV, double *tipVector,
                            int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                            const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling
                            )
{
  int 
    i, 
    j, 
    k, 
    l,
    addScale = 0;

  //int scaling = 0;

  double
    *x1,
    *x2,
    *x3,
    max;
  PLL_ALIGN_BEGIN double
    maxima[2] PLL_ALIGN_END,
    EV_t[16] PLL_ALIGN_END;

  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        double *uX1, *uX2;
        PLL_ALIGN_BEGIN double umpX1[256] PLL_ALIGN_END, umpX2[256] PLL_ALIGN_END;


        for (i = 1; i < 16; i++)
        {
          __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
          __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));       

          for (j = 0; j < 4; j++)

            for (k = 0; k < 4; k++) {
              __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
            }

          for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {
              __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);

            }
        }       

        for (i = 0; i < n; i++)
        {
          x3 = &x3_start[i * 16];


          uX1 = &umpX1[16 * tipX1[i]];
          uX2 = &umpX2[16 * tipX2[i]];                      

          for (j = 0; j < 4; j++)
          {                                                                                
            __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
            __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );


            __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
            __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );


            //
            // multiply left * right
            //

            __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
            __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );


            //
            // multiply with EV matrix (!?)
            //

            __m128d EV_t_l0_k0 = EVV[0];
            __m128d EV_t_l0_k2 = EVV[1];
            __m128d EV_t_l1_k0 = EVV[2];
            __m128d EV_t_l1_k2 = EVV[3];
            __m128d EV_t_l2_k0 = EVV[4];
            __m128d EV_t_l2_k2 = EVV[5];
            __m128d EV_t_l3_k0 = EVV[6]; 
            __m128d EV_t_l3_k2 = EVV[7];

            EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
            EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

            EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
            EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

            EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
            EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

            EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
            EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

            EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
            EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
            EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

            _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
            _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
          }
        }
      }
      break;
    case PLL_TIP_INNER:
      { 
        double *uX1;
        PLL_ALIGN_BEGIN double umpX1[256] PLL_ALIGN_END;


        for (i = 1; i < 16; i++)
        {
          __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
          __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));       

          for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {            
              __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
              __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);

              __m128d acc = _mm_setzero_pd();

              acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
              acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));

              acc = _mm_hadd_pd(acc, acc);
              _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);                
            }
        }

        for (i = 0; i < n; i++)
        {
          __m128d maxv =_mm_setzero_pd();

          x2 = &x2_start[i * 16];
          x3 = &x3_start[i * 16];

          uX1 = &umpX1[16 * tipX1[i]];       

          for (j = 0; j < 4; j++)
          {

            //
            // multiply/add right side
            //
            double *x2_p = &x2[j*4];
            double *right_k0_p = &right[j*16];
            double *right_k1_p = &right[j*16 + 1*4];
            double *right_k2_p = &right[j*16 + 2*4];
            double *right_k3_p = &right[j*16 + 3*4];
            __m128d x2_0 = _mm_load_pd( &x2_p[0] );
            __m128d x2_2 = _mm_load_pd( &x2_p[2] );

            __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
            __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
            __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
            __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
            __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
            __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
            __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
            __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );



            right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
            right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

            right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
            right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
            right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
            right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);


            right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
            right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


            right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
            right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
            right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
            right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

            {
              //
              // load left side from tip vector
              //

              __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
              __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );


              //
              // multiply left * right
              //

              __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
              __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );


              //
              // multiply with EV matrix (!?)
              //                                  

              __m128d EV_t_l0_k0 = EVV[0];
              __m128d EV_t_l0_k2 = EVV[1];
              __m128d EV_t_l1_k0 = EVV[2];
              __m128d EV_t_l1_k2 = EVV[3];
              __m128d EV_t_l2_k0 = EVV[4];
              __m128d EV_t_l2_k2 = EVV[5];
              __m128d EV_t_l3_k0 = EVV[6]; 
              __m128d EV_t_l3_k2 = EVV[7];


              EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
              EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
              EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

              EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
              EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

              EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
              EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

              EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
              EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
              EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

              EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
              EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
              EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

              EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

              values[j * 2]     = EV_t_l0_k0;
              values[j * 2 + 1] = EV_t_l2_k0;                              

              maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
              maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));                  
            }
          }


          _mm_store_pd(maxima, maxv);

          max = PLL_MAX(maxima[0], maxima[1]);

          if(max < PLL_MINLIKELIHOOD)
          {
            __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

            _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));       
            _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
            _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
            _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
            _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));       
            _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
            _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
            _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));        

             if(!fastScaling)
               ex3[i] += 1;
             else
               addScale += wgt[i];

          }
          else
          {
            _mm_store_pd(&x3[0], values[0]);       
            _mm_store_pd(&x3[2], values[1]);
            _mm_store_pd(&x3[4], values[2]);
            _mm_store_pd(&x3[6], values[3]);
            _mm_store_pd(&x3[8], values[4]);       
            _mm_store_pd(&x3[10], values[5]);
            _mm_store_pd(&x3[12], values[6]);
            _mm_store_pd(&x3[14], values[7]);
          }
        }
      }
      break;
    case PLL_INNER_INNER:

      for (i = 0; i < n; i++)
      {
        __m128d maxv =_mm_setzero_pd();


        x1 = &x1_start[i * 16];
        x2 = &x2_start[i * 16];
        x3 = &x3_start[i * 16];

        for (j = 0; j < 4; j++)
        {

          double *x1_p = &x1[j*4];
          double *left_k0_p = &left[j*16];
          double *left_k1_p = &left[j*16 + 1*4];
          double *left_k2_p = &left[j*16 + 2*4];
          double *left_k3_p = &left[j*16 + 3*4];

          __m128d x1_0 = _mm_load_pd( &x1_p[0] );
          __m128d x1_2 = _mm_load_pd( &x1_p[2] );

          __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
          __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
          __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
          __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
          __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
          __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
          __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
          __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);


          //
          // multiply/add right side
          //
          double *x2_p = &x2[j*4];
          double *right_k0_p = &right[j*16];
          double *right_k1_p = &right[j*16 + 1*4];
          double *right_k2_p = &right[j*16 + 2*4];
          double *right_k3_p = &right[j*16 + 3*4];
          __m128d x2_0 = _mm_load_pd( &x2_p[0] );
          __m128d x2_2 = _mm_load_pd( &x2_p[2] );

          __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
          __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
          __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
          __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
          __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
          __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
          __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
          __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);       

          //
          // multiply left * right
          //

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );


          //
          // multiply with EV matrix (!?)
          //         

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6]; 
          __m128d EV_t_l3_k2 = EVV[7];


          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );


          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );


          values[j * 2] = EV_t_l0_k0;
          values[j * 2 + 1] = EV_t_l2_k0;                           

          maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
          maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
        }


        _mm_store_pd(maxima, maxv);

        max = PLL_MAX(maxima[0], maxima[1]);

        if(max < PLL_MINLIKELIHOOD)
        {
          __m128d sv = _mm_set1_pd(PLL_TWOTOTHE256);

          _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));         
          _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
          _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
          _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
          _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));         
          _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
          _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
          _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));          

           if(!fastScaling)
             ex3[i] += 1;
           else
             addScale += wgt[i];        
        }
        else
        {
          _mm_store_pd(&x3[0], values[0]);         
          _mm_store_pd(&x3[2], values[1]);
          _mm_store_pd(&x3[4], values[2]);
          _mm_store_pd(&x3[6], values[3]);
          _mm_store_pd(&x3[8], values[4]);         
          _mm_store_pd(&x3[10], values[5]);
          _mm_store_pd(&x3[12], values[6]);
          _mm_store_pd(&x3[14], values[7]);
        }        
      }

      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}


/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR CAT (Optimized SSE3 version for DNA data)

    This is the SSE3 optimized version of ::newviewCAT_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b CAT
    model of rate heterogeneity.

    @note
    For more details and function argument description check the function ::newviewCAT_FLEX
*/
static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
                           double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                           int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                           int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling)
{
  double
    *le,
    *ri,
    *x1,
    *x2, 
    *x3;
  PLL_ALIGN_BEGIN double
    EV_t[16] PLL_ALIGN_END;

  int 
    i, 
    j, 
    scale, 
    addScale = 0;

  __m128d
    minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD ),
                      sc = _mm_set1_pd(PLL_TWOTOTHE256),
                      EVV[8];  

  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];

  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);

  switch(tipCase)
  {
    case PLL_TIP_TIP:      
      for (i = 0; i < n; i++)
      {  
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        x3 = &x3_start[i * 4];

        le =  &left[cptr[i] * 16];
        ri =  &right[cptr[i] * 16];

        __m128d x1_0 = _mm_load_pd( &x1[0] );
        __m128d x1_2 = _mm_load_pd( &x1[2] );

        __m128d left_k0_0 = _mm_load_pd( &le[0] );
        __m128d left_k0_2 = _mm_load_pd( &le[2] );
        __m128d left_k1_0 = _mm_load_pd( &le[4] );
        __m128d left_k1_2 = _mm_load_pd( &le[6] );
        __m128d left_k2_0 = _mm_load_pd( &le[8] );
        __m128d left_k2_2 = _mm_load_pd( &le[10] );
        __m128d left_k3_0 = _mm_load_pd( &le[12] );
        __m128d left_k3_2 = _mm_load_pd( &le[14] );

        left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
        left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

        left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
        left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
        left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

        left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
        left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

        left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
        left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
        left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

        __m128d x2_0 = _mm_load_pd( &x2[0] );
        __m128d x2_2 = _mm_load_pd( &x2[2] );

        __m128d right_k0_0 = _mm_load_pd( &ri[0] );
        __m128d right_k0_2 = _mm_load_pd( &ri[2] );
        __m128d right_k1_0 = _mm_load_pd( &ri[4] );
        __m128d right_k1_2 = _mm_load_pd( &ri[6] );
        __m128d right_k2_0 = _mm_load_pd( &ri[8] );
        __m128d right_k2_2 = _mm_load_pd( &ri[10] );
        __m128d right_k3_0 = _mm_load_pd( &ri[12] );
        __m128d right_k3_2 = _mm_load_pd( &ri[14] );

        right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
        right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

        right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
        right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
        right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

        right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
        right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

        right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
        right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
        right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);         

        __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
        __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );           

        __m128d EV_t_l0_k0 = EVV[0];
        __m128d EV_t_l0_k2 = EVV[1];
        __m128d EV_t_l1_k0 = EVV[2];
        __m128d EV_t_l1_k2 = EVV[3];
        __m128d EV_t_l2_k0 = EVV[4];
        __m128d EV_t_l2_k2 = EVV[5];
        __m128d EV_t_l3_k0 = EVV[6];
        __m128d EV_t_l3_k2 = EVV[7];

        EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
        EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

        EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
        EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

        EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

        EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
        EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

        EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
        EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
        EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );      

        _mm_store_pd(x3, EV_t_l0_k0);
        _mm_store_pd(&x3[2], EV_t_l2_k0);                                   
      }
      break;
    case PLL_TIP_INNER:      
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &x2_start[4 * i];
        x3 = &x3_start[4 * i];

        le =  &left[cptr[i] * 16];
        ri =  &right[cptr[i] * 16];

        __m128d x1_0 = _mm_load_pd( &x1[0] );
        __m128d x1_2 = _mm_load_pd( &x1[2] );

        __m128d left_k0_0 = _mm_load_pd( &le[0] );
        __m128d left_k0_2 = _mm_load_pd( &le[2] );
        __m128d left_k1_0 = _mm_load_pd( &le[4] );
        __m128d left_k1_2 = _mm_load_pd( &le[6] );
        __m128d left_k2_0 = _mm_load_pd( &le[8] );
        __m128d left_k2_2 = _mm_load_pd( &le[10] );
        __m128d left_k3_0 = _mm_load_pd( &le[12] );
        __m128d left_k3_2 = _mm_load_pd( &le[14] );

        left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
        left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

        left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
        left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
        left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

        left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
        left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

        left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
        left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
        left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

        __m128d x2_0 = _mm_load_pd( &x2[0] );
        __m128d x2_2 = _mm_load_pd( &x2[2] );

        __m128d right_k0_0 = _mm_load_pd( &ri[0] );
        __m128d right_k0_2 = _mm_load_pd( &ri[2] );
        __m128d right_k1_0 = _mm_load_pd( &ri[4] );
        __m128d right_k1_2 = _mm_load_pd( &ri[6] );
        __m128d right_k2_0 = _mm_load_pd( &ri[8] );
        __m128d right_k2_2 = _mm_load_pd( &ri[10] );
        __m128d right_k3_0 = _mm_load_pd( &ri[12] );
        __m128d right_k3_2 = _mm_load_pd( &ri[14] );

        right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
        right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

        right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
        right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
        right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

        right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
        right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

        right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
        right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
        right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);         

        __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
        __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

        __m128d EV_t_l0_k0 = EVV[0];
        __m128d EV_t_l0_k2 = EVV[1];
        __m128d EV_t_l1_k0 = EVV[2];
        __m128d EV_t_l1_k2 = EVV[3];
        __m128d EV_t_l2_k0 = EVV[4];
        __m128d EV_t_l2_k2 = EVV[5];
        __m128d EV_t_l3_k0 = EVV[6];
        __m128d EV_t_l3_k2 = EVV[7];


        EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
        EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

        EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
        EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

        EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

        EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
        EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

        EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
        EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
        EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );                                       

        scale = 1;

        __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
        else
        {
          v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
        }

        if(scale)
        {                     
          _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
          _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));                   

           if(!fastScaling)
             ex3[i] += 1;
           else
             addScale += wgt[i];          
        }       
        else
        {
          _mm_store_pd(x3, EV_t_l0_k0);
          _mm_store_pd(&x3[2], EV_t_l2_k0);
        }


      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &x1_start[4 * i];
        x2 = &x2_start[4 * i];
        x3 = &x3_start[4 * i];

        le =  &left[cptr[i] * 16];
        ri =  &right[cptr[i] * 16];

        __m128d x1_0 = _mm_load_pd( &x1[0] );
        __m128d x1_2 = _mm_load_pd( &x1[2] );

        __m128d left_k0_0 = _mm_load_pd( &le[0] );
        __m128d left_k0_2 = _mm_load_pd( &le[2] );
        __m128d left_k1_0 = _mm_load_pd( &le[4] );
        __m128d left_k1_2 = _mm_load_pd( &le[6] );
        __m128d left_k2_0 = _mm_load_pd( &le[8] );
        __m128d left_k2_2 = _mm_load_pd( &le[10] );
        __m128d left_k3_0 = _mm_load_pd( &le[12] );
        __m128d left_k3_2 = _mm_load_pd( &le[14] );

        left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
        left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

        left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
        left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
        left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
        left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

        left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
        left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

        left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
        left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
        left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
        left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

        __m128d x2_0 = _mm_load_pd( &x2[0] );
        __m128d x2_2 = _mm_load_pd( &x2[2] );

        __m128d right_k0_0 = _mm_load_pd( &ri[0] );
        __m128d right_k0_2 = _mm_load_pd( &ri[2] );
        __m128d right_k1_0 = _mm_load_pd( &ri[4] );
        __m128d right_k1_2 = _mm_load_pd( &ri[6] );
        __m128d right_k2_0 = _mm_load_pd( &ri[8] );
        __m128d right_k2_2 = _mm_load_pd( &ri[10] );
        __m128d right_k3_0 = _mm_load_pd( &ri[12] );
        __m128d right_k3_2 = _mm_load_pd( &ri[14] );

        right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
        right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

        right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
        right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
        right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
        right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

        right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
        right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

        right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
        right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
        right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
        right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);         

        __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
        __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

        __m128d EV_t_l0_k0 = EVV[0];
        __m128d EV_t_l0_k2 = EVV[1];
        __m128d EV_t_l1_k0 = EVV[2];
        __m128d EV_t_l1_k2 = EVV[3];
        __m128d EV_t_l2_k0 = EVV[4];
        __m128d EV_t_l2_k2 = EVV[5];
        __m128d EV_t_l3_k0 = EVV[6];
        __m128d EV_t_l3_k2 = EVV[7];


        EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
        EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

        EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
        EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

        EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
        EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

        EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
        EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

        EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
        EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
        EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

        EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );                                              

        scale = 1;

        __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
        else
        {
          v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
        }

        if(scale)
        {                     
          _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
          _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));                   

          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];   
        }       
        else
        {
          _mm_store_pd(x3, EV_t_l0_k0);
          _mm_store_pd(&x3[2], EV_t_l2_k0);
        }

      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}
#endif

/** @brief Check whether the position \a pos in bitvector \a x is a gap
    
    @param x
      A bitvector represented by unsigned integers

    @param pos
      Position to check in \a x if it is set (i.e. it is a gap) 

    @return
      Returns the value of the bit vector (\b 1 if set, \b 0 if not)
*/
//#ifndef __clang__
//__inline
//#endif
pllBoolean isGap(unsigned int *x, int pos)
{
  return (x[pos / 32] & mask32[pos % 32]);
}

/** @brief Check whether the position \a pos in bitvector \a x is \b NOT a gap
    
    @param x
      A bitvector represented by unsigned integers

    @param pos
      Position to check in \a x if it is \b NOT set (i.e. it is \b NOT a gap) 

    @return
      Returns the value of the bit vector (\b 1 if set, \b 0 if not)
*/
//#ifndef __clang__
//__inline
//#endif
pllBoolean noGap(unsigned int *x, int pos)
{
  return (!(x[pos / 32] & mask32[pos % 32]));
}

#if (!defined(__AVX) && defined(__SSE3))
/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR CAT with memory saving (Optimized SSE3 version for DNA data)

    This is the SSE3 optimized version of ::newviewCAT_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b CAT
    model of rate heterogeneity. The memory saving technique is incorporated.

    @note
    For more details and function argument description check the function ::newviewCAT_FLEX
*/
static void newviewGTRCAT_SAVE( int tipCase,  double *EV,  int *cptr,
                                double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
                                int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le,
    *ri,
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start, 
    *x3_ptr = x3_start;
  PLL_ALIGN_BEGIN double
    EV_t[16] PLL_ALIGN_END;

  int 
    i, 
    j, 
    scale, 
    scaleGap = 0,
    addScale = 0;

  __m128d
    minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD ),
                      sc = _mm_set1_pd(PLL_TWOTOTHE256),
                      EVV[8];  

  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];

  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);

  {
    x1 = x1_gapColumn;        
    x2 = x2_gapColumn;
    x3 = x3_gapColumn;

    le =  &left[maxCats * 16];           
    ri =  &right[maxCats * 16];                                                  

    __m128d x1_0 = _mm_load_pd( &x1[0] );
    __m128d x1_2 = _mm_load_pd( &x1[2] );

    __m128d left_k0_0 = _mm_load_pd( &le[0] );
    __m128d left_k0_2 = _mm_load_pd( &le[2] );
    __m128d left_k1_0 = _mm_load_pd( &le[4] );
    __m128d left_k1_2 = _mm_load_pd( &le[6] );
    __m128d left_k2_0 = _mm_load_pd( &le[8] );
    __m128d left_k2_2 = _mm_load_pd( &le[10] );
    __m128d left_k3_0 = _mm_load_pd( &le[12] );
    __m128d left_k3_2 = _mm_load_pd( &le[14] );

    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

    __m128d x2_0 = _mm_load_pd( &x2[0] );
    __m128d x2_2 = _mm_load_pd( &x2[2] );

    __m128d right_k0_0 = _mm_load_pd( &ri[0] );
    __m128d right_k0_2 = _mm_load_pd( &ri[2] );
    __m128d right_k1_0 = _mm_load_pd( &ri[4] );
    __m128d right_k1_2 = _mm_load_pd( &ri[6] );
    __m128d right_k2_0 = _mm_load_pd( &ri[8] );
    __m128d right_k2_2 = _mm_load_pd( &ri[10] );
    __m128d right_k3_0 = _mm_load_pd( &ri[12] );
    __m128d right_k3_2 = _mm_load_pd( &ri[14] );

    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);     

    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

    __m128d EV_t_l0_k0 = EVV[0];
    __m128d EV_t_l0_k2 = EVV[1];
    __m128d EV_t_l1_k0 = EVV[2];
    __m128d EV_t_l1_k2 = EVV[3];
    __m128d EV_t_l2_k0 = EVV[4];
    __m128d EV_t_l2_k2 = EVV[5];
    __m128d EV_t_l3_k0 = EVV[6];
    __m128d EV_t_l3_k2 = EVV[7];

    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );                                   

    if(tipCase != PLL_TIP_TIP)
    {    
      scale = 1;

      __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
      if(_mm_movemask_pd( v1 ) != 3)
        scale = 0;
      else
      {
        v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
      }

      if(scale)
      {               
        _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
        _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));                     

        scaleGap = PLL_TRUE;       
      } 
      else
      {
        _mm_store_pd(x3, EV_t_l0_k0);
        _mm_store_pd(&x3[2], EV_t_l2_k0);
      }
    }
    else
    {
      _mm_store_pd(x3, EV_t_l0_k0);
      _mm_store_pd(&x3[2], EV_t_l2_k0);
    }
  }


  switch(tipCase)
  {
    case PLL_TIP_TIP:      
      for (i = 0; i < n; i++)
      {
        if(noGap(x3_gap, i))
        {
          x1 = &(tipVector[4 * tipX1[i]]);
          x2 = &(tipVector[4 * tipX2[i]]);

          x3 = x3_ptr;

          if(isGap(x1_gap, i))
            le =  &left[maxCats * 16];
          else            
            le =  &left[cptr[i] * 16];    

          if(isGap(x2_gap, i))
            ri =  &right[maxCats * 16];
          else            
            ri =  &right[cptr[i] * 16];

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);       

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );                 

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];

          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );    

          _mm_store_pd(x3, EV_t_l0_k0);
          _mm_store_pd(&x3[2], EV_t_l2_k0);                                 

          x3_ptr += 4;
        }
      }
      break;
    case PLL_TIP_INNER:      
      for (i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)
            {
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];
            }
        }
        else
        {             
          x1 = &(tipVector[4 * tipX1[i]]);

          x2 = x2_ptr;
          x3 = x3_ptr;

          if(isGap(x1_gap, i))
            le =  &left[maxCats * 16];
          else
            le =  &left[cptr[i] * 16];

          if(isGap(x2_gap, i))
          {              
            ri =  &right[maxCats * 16];
            x2 = x2_gapColumn;
          }
          else
          {
            ri =  &right[cptr[i] * 16];
            x2 = x2_ptr;
            x2_ptr += 4;
          }                               

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);       

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];


          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );                                     

          scale = 1;

          __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
          else
          {
            v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }

          if(scale)
          {                   
            _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
            _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));                 
            
            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];         
          }     
          else
          {
            _mm_store_pd(x3, EV_t_l0_k0);
            _mm_store_pd(&x3[2], EV_t_l2_k0);
          }

          x3_ptr += 4;
        }

      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)
            {
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];
            }
        }
        else
        {            
          x3 = x3_ptr;

          if(isGap(x1_gap, i))
          {
            x1 = x1_gapColumn;
            le =  &left[maxCats * 16];
          }
          else
          {
            le =  &left[cptr[i] * 16];
            x1 = x1_ptr;
            x1_ptr += 4;
          }

          if(isGap(x2_gap, i))  
          {
            x2 = x2_gapColumn;
            ri =  &right[maxCats * 16];     
          }
          else
          {
            ri =  &right[cptr[i] * 16];
            x2 = x2_ptr;
            x2_ptr += 4;
          }                               

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);       

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];


          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );                                            

          scale = 1;

          __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
          else
          {
            v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }

          if(scale)
          {                   
            _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
            _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));                 

            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];         
          }     
          else
          {
            _mm_store_pd(x3, EV_t_l0_k0);
            _mm_store_pd(&x3[2], EV_t_l2_k0);
          }

          x3_ptr += 4;
        }
      }
      break;
    default:
      assert(0);
  }


  if(fastScaling)
    *scalerIncrement = addScale;
}

/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR GAMMA with memory saving (Optimized SSE3 version for AA data)

    This is the SSE3 optimized version of ::newviewGAMMA_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b GAMMA
    model of rate heterogeneity. The memory saving technique is incorporated.

    @note
    For more details and function argument description check the function ::newviewGAMMA_FLEX
*/
static void newviewGTRGAMMAPROT_GAPPED_SAVE(int tipCase,
                                            double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                                            int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                            int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                            unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,  
                                            double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
                                            )
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0,   
       gapScaling = 0;
  double 
    *vl, *vr, *x1v, *x2v,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x3_ptr = x3;



  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        double umpX1[1840], umpX2[1840];

        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];
            double *rr =  &right[k * 20];

            __m128d umpX1v = _mm_setzero_pd();
            __m128d umpX2v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
              umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));                                 
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
            umpX2v = _mm_hadd_pd(umpX2v, umpX2v);

            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);
            _mm_storel_pd(&umpX2[80 * i + k], umpX2v);
          }
        }

        {
          uX1 = &umpX1[1760];
          uX2 = &umpX2[1760];

          for(j = 0; j < 4; j++)
          {
            v = &x3_gapColumn[j * 20];

            __m128d zero =  _mm_setzero_pd();
            for(k = 0; k < 20; k+=2)                                
              _mm_store_pd(&v[k], zero);

            for(k = 0; k < 20; k++)
            { 
              double *eev = &extEV[k * 20];
              x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(l = 0; l < 20; l+=2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d ee = _mm_load_pd(&eev[l]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[l], vv);
              }
            }
          }        
        }       

        for(i = 0; i < n; i++)
        {
          if(!(x3_gap[i / 32] & mask32[i % 32]))
          {
            uX1 = &umpX1[80 * tipX1[i]];
            uX2 = &umpX2[80 * tipX2[i]];

            for(j = 0; j < 4; j++)
            {
              v = &x3_ptr[j * 20];


              __m128d zero =  _mm_setzero_pd();
              for(k = 0; k < 20; k+=2)                              
                _mm_store_pd(&v[k], zero);

              for(k = 0; k < 20; k++)
              { 
                double *eev = &extEV[k * 20];
                x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
                __m128d x1px2v = _mm_set1_pd(x1px2);

                for(l = 0; l < 20; l+=2)
                {
                  __m128d vv = _mm_load_pd(&v[l]);
                  __m128d ee = _mm_load_pd(&eev[l]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                  _mm_store_pd(&v[l], vv);
                }
              }
            }      
            x3_ptr += 80;
          }
        }
      }
      break;
    case PLL_TIP_INNER:
      {
        double umpX1[1840], ump_x2[20];


        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];

            __m128d umpX1v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));                                                 
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);                               
            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);          

          }
        }

        {
          uX1 = &umpX1[1760];

          for(k = 0; k < 4; k++)
          {
            v = &(x2_gapColumn[k * 20]);

            for(l = 0; l < 20; l++)
            {              
              double *r =  &right[k * 400 + l * 20];
              __m128d ump_x2v = _mm_setzero_pd();           

              for(j = 0; j < 20; j+= 2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d rr = _mm_load_pd(&r[j]);
                ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
              }

              ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

              _mm_storel_pd(&ump_x2[l], ump_x2v);                                    
            }

            v = &(x3_gapColumn[20 * k]);

            __m128d zero =  _mm_setzero_pd();
            for(l = 0; l < 20; l+=2)                                
              _mm_store_pd(&v[l], zero);

            for(l = 0; l < 20; l++)
            {
              double *eev = &extEV[l * 20];
              x1px2 = uX1[k * 20 + l]  * ump_x2[l];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d ee = _mm_load_pd(&eev[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[j], vv);
              }                             
            }                   

          }

          { 
            v = x3_gapColumn;
            __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

            scale = 1;
            for(l = 0; scale && (l < 80); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }             
          }


          if (scale)
          {
            gapScaling = 1;
            __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

            for(l = 0; l < 80; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);                  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));      
            }                                                          
          }
        }

        for (i = 0; i < n; i++)
        {           
          if((x3_gap[i / 32] & mask32[i % 32]))
          {            
            if(gapScaling)
            {   
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];                  
            }
          }
          else
          {
            uX1 = &umpX1[80 * tipX1[i]];

            if(x2_gap[i / 32] & mask32[i % 32])
              x2v = x2_gapColumn;
            else
            {
              x2v = x2_ptr;
              x2_ptr += 80;
            }

            for(k = 0; k < 4; k++)
            {
              v = &(x2v[k * 20]);

              for(l = 0; l < 20; l++)
              {            
                double *r =  &right[k * 400 + l * 20];
                __m128d ump_x2v = _mm_setzero_pd();         

                for(j = 0; j < 20; j+= 2)
                {
                  __m128d vv = _mm_load_pd(&v[j]);
                  __m128d rr = _mm_load_pd(&r[j]);
                  ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
                }

                ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

                _mm_storel_pd(&ump_x2[l], ump_x2v);                                  
              }

              v = &x3_ptr[20 * k];

              __m128d zero =  _mm_setzero_pd();
              for(l = 0; l < 20; l+=2)                              
                _mm_store_pd(&v[l], zero);

              for(l = 0; l < 20; l++)
              {
                double *eev = &extEV[l * 20];
                x1px2 = uX1[k * 20 + l]  * ump_x2[l];
                __m128d x1px2v = _mm_set1_pd(x1px2);

                for(j = 0; j < 20; j+=2)
                {
                  __m128d vv = _mm_load_pd(&v[j]);
                  __m128d ee = _mm_load_pd(&eev[j]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                  _mm_store_pd(&v[j], vv);
                }                                   
              }                 

            }


            { 
              v = x3_ptr;
              __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

              scale = 1;
              for(l = 0; scale && (l < 80); l += 2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d v1 = _mm_and_pd(vv, absMask.m);
                v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
                if(_mm_movemask_pd( v1 ) != 3)
                  scale = 0;
              }           
            }


            if (scale)
            {
              __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

              for(l = 0; l < 80; l+=2)
              {
                __m128d ex3v = _mm_load_pd(&v[l]);                
                _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));    
              }                           
              
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];                   
            }

            x3_ptr += 80;
          }
        }
      }
      break;
    case PLL_INNER_INNER:
      {
        for(k = 0; k < 4; k++)
        {
          vl = &(x1_gapColumn[20 * k]);
          vr = &(x2_gapColumn[20 * k]);
          v =  &(x3_gapColumn[20 * k]);

          __m128d zero =  _mm_setzero_pd();
          for(l = 0; l < 20; l+=2)                                  
            _mm_store_pd(&v[l], zero);

          for(l = 0; l < 20; l++)
          {              
            {
              __m128d al = _mm_setzero_pd();
              __m128d ar = _mm_setzero_pd();

              double *ll   = &left[k * 400 + l * 20];
              double *rr   = &right[k * 400 + l * 20];
              double *EVEV = &extEV[20 * l];

              for(j = 0; j < 20; j+=2)
              {
                __m128d lv  = _mm_load_pd(&ll[j]);
                __m128d rv  = _mm_load_pd(&rr[j]);
                __m128d vll = _mm_load_pd(&vl[j]);
                __m128d vrr = _mm_load_pd(&vr[j]);

                al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
              }                  

              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);

              al = _mm_mul_pd(al, ar);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv  = _mm_load_pd(&v[j]);
                __m128d EVV = _mm_load_pd(&EVEV[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                _mm_store_pd(&v[j], vv);
              }                                           
            }            

          }
        }


        { 
          v = x3_gapColumn;
          __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

          scale = 1;
          for(l = 0; scale && (l < 80); l += 2)
          {
            __m128d vv = _mm_load_pd(&v[l]);
            __m128d v1 = _mm_and_pd(vv, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }               
        }

        if (scale)
        {
          gapScaling = 1;
          __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

          for(l = 0; l < 80; l+=2)
          {
            __m128d ex3v = _mm_load_pd(&v[l]);            
            _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));        
          }                               


        }
      }

      for (i = 0; i < n; i++)
      {
        if(x3_gap[i / 32] & mask32[i % 32])
        {            
          if(gapScaling)
          {     
            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];                              
          }
        }
        else
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

          for(k = 0; k < 4; k++)
          {
            vl = &(x1v[20 * k]);
            vr = &(x2v[20 * k]);
            v =  &x3_ptr[20 * k];

            __m128d zero =  _mm_setzero_pd();
            for(l = 0; l < 20; l+=2)                                
              _mm_store_pd(&v[l], zero);

            for(l = 0; l < 20; l++)
            {            
              {
                __m128d al = _mm_setzero_pd();
                __m128d ar = _mm_setzero_pd();

                double *ll   = &left[k * 400 + l * 20];
                double *rr   = &right[k * 400 + l * 20];
                double *EVEV = &extEV[20 * l];

                for(j = 0; j < 20; j+=2)
                {
                  __m128d lv  = _mm_load_pd(&ll[j]);
                  __m128d rv  = _mm_load_pd(&rr[j]);
                  __m128d vll = _mm_load_pd(&vl[j]);
                  __m128d vrr = _mm_load_pd(&vr[j]);

                  al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                  ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
                }                

                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);

                al = _mm_mul_pd(al, ar);

                for(j = 0; j < 20; j+=2)
                {
                  __m128d vv  = _mm_load_pd(&v[j]);
                  __m128d EVV = _mm_load_pd(&EVEV[j]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                  _mm_store_pd(&v[j], vv);
                }                                                 
              }          

            }
          }



          { 
            v = x3_ptr;
            __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

            scale = 1;
            for(l = 0; scale && (l < 80); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }             
          }


          if (scale)
          {
            __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

            for(l = 0; l < 80; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);                  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));      
            }                             

            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];                         
          }
          x3_ptr += 80;
        }
      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;  
}



/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR GAMMA (Optimized SSE3 version for AA data)

    This is the SSE3 optimized version of ::newviewGAMMA_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b GAMMA
    model of rate heterogeneity.

    @note
    For more details and function argument description check the function ::newviewGAMMA_FLEX
*/
static void newviewGTRGAMMAPROT(int tipCase,
                                double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                                int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling)
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr;



  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        double umpX1[1840], umpX2[1840];

        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];
            double *rr =  &right[k * 20];

            __m128d umpX1v = _mm_setzero_pd();
            __m128d umpX2v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
              umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));                                 
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
            umpX2v = _mm_hadd_pd(umpX2v, umpX2v);

            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);
            _mm_storel_pd(&umpX2[80 * i + k], umpX2v);

          }
        }

        for(i = 0; i < n; i++)
        {
          uX1 = &umpX1[80 * tipX1[i]];
          uX2 = &umpX2[80 * tipX2[i]];

          for(j = 0; j < 4; j++)
          {
            v = &x3[i * 80 + j * 20];


            __m128d zero =  _mm_setzero_pd();
            for(k = 0; k < 20; k+=2)                                
              _mm_store_pd(&v[k], zero);

            for(k = 0; k < 20; k++)
            { 
              double *eev = &extEV[k * 20];
              x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(l = 0; l < 20; l+=2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d ee = _mm_load_pd(&eev[l]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[l], vv);
              }
            }


          }        
        }
      }
      break;
    case PLL_TIP_INNER:
      {
        double umpX1[1840], ump_x2[20];


        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];

            __m128d umpX1v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));                                                 
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);                               
            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);          


          }
        }

        for (i = 0; i < n; i++)
        {
          uX1 = &umpX1[80 * tipX1[i]];

          for(k = 0; k < 4; k++)
          {
            v = &(x2[80 * i + k * 20]);

            for(l = 0; l < 20; l++)
            {              
              double *r =  &right[k * 400 + l * 20];
              __m128d ump_x2v = _mm_setzero_pd();           

              for(j = 0; j < 20; j+= 2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d rr = _mm_load_pd(&r[j]);
                ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
              }

              ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

              _mm_storel_pd(&ump_x2[l], ump_x2v);                                    
            }

            v = &(x3[80 * i + 20 * k]);

            __m128d zero =  _mm_setzero_pd();
            for(l = 0; l < 20; l+=2)                                
              _mm_store_pd(&v[l], zero);

            for(l = 0; l < 20; l++)
            {
              double *eev = &extEV[l * 20];
              x1px2 = uX1[k * 20 + l]  * ump_x2[l];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d ee = _mm_load_pd(&eev[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[j], vv);
              }                             
            }                   

          }


          { 
            v = &(x3[80 * i]);
            __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

            scale = 1;
            for(l = 0; scale && (l < 80); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }             
          }


          if (scale)
          {

            __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

            for(l = 0; l < 80; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);                  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));      
            }                             


            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];

          }
        }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
      {
        for(k = 0; k < 4; k++)
        {
          vl = &(x1[80 * i + 20 * k]);
          vr = &(x2[80 * i + 20 * k]);
          v =  &(x3[80 * i + 20 * k]);


          __m128d zero =  _mm_setzero_pd();
          for(l = 0; l < 20; l+=2)                                  
            _mm_store_pd(&v[l], zero);


          for(l = 0; l < 20; l++)
          {              

            {
              __m128d al = _mm_setzero_pd();
              __m128d ar = _mm_setzero_pd();

              double *ll   = &left[k * 400 + l * 20];
              double *rr   = &right[k * 400 + l * 20];
              double *EVEV = &extEV[20 * l];

              for(j = 0; j < 20; j+=2)
              {
                __m128d lv  = _mm_load_pd(&ll[j]);
                __m128d rv  = _mm_load_pd(&rr[j]);
                __m128d vll = _mm_load_pd(&vl[j]);
                __m128d vrr = _mm_load_pd(&vr[j]);

                al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
              }                  

              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);

              al = _mm_mul_pd(al, ar);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv  = _mm_load_pd(&v[j]);
                __m128d EVV = _mm_load_pd(&EVEV[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                _mm_store_pd(&v[j], vv);
              }                                           
            }            

          }
        }



        { 
          v = &(x3[80 * i]);
          __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

          scale = 1;
          for(l = 0; scale && (l < 80); l += 2)
          {
            __m128d vv = _mm_load_pd(&v[l]);
            __m128d v1 = _mm_and_pd(vv, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }               
        }


        if (scale)
        {

          __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

          for(l = 0; l < 80; l+=2)
          {
            __m128d ex3v = _mm_load_pd(&v[l]);            
            _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));        
          }                               


          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];
        }
      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}



/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR CAT (Optimized SSE3 version for AA data)

    This is the SSE3 optimized version of ::newviewCAT_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b CAT
    model of rate heterogeneity.

    @note
    For more details and function argument description check the function ::newviewCAT_FLEX
*/
static void newviewGTRCATPROT(int tipCase, double *extEV,
                              int *cptr,
                              double *x1, double *x2, double *x3, double *tipVector,
                              int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                              int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;

  int i, l, j, scale, addScale = 0;

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        for (i = 0; i < n; i++)
        {
          le = &left[cptr[i] * 400];
          ri = &right[cptr[i] * 400];

          vl = &(tipVector[20 * tipX1[i]]);
          vr = &(tipVector[20 * tipX2[i]]);
          v  = &x3[20 * i];

          for(l = 0; l < 20; l+=2)
            _mm_store_pd(&v[l], _mm_setzero_pd());                      


          for(l = 0; l < 20; l++)
          {
            __m128d x1v = _mm_setzero_pd();
            __m128d x2v = _mm_setzero_pd();      
            double 
              *ev = &extEV[l * 20],
              *lv = &le[l * 20],
              *rv = &ri[l * 20];

            for(j = 0; j < 20; j+=2)
            {
              x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                  
              x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
            }

            x1v = _mm_hadd_pd(x1v, x1v);
            x2v = _mm_hadd_pd(x2v, x2v);

            x1v = _mm_mul_pd(x1v, x2v);

            for(j = 0; j < 20; j+=2)
            {
              __m128d vv = _mm_load_pd(&v[j]);
              vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
              _mm_store_pd(&v[j], vv);
            }               

          }        
        }
      }
      break;
    case PLL_TIP_INNER:
      {
        for (i = 0; i < n; i++)
        {
          le = &left[cptr[i] * 400];
          ri = &right[cptr[i] * 400];

          vl = &(tipVector[20 * tipX1[i]]);
          vr = &x2[20 * i];
          v  = &x3[20 * i];

          for(l = 0; l < 20; l+=2)
            _mm_store_pd(&v[l], _mm_setzero_pd());                      



          for(l = 0; l < 20; l++)
          {

            __m128d x1v = _mm_setzero_pd();
            __m128d x2v = _mm_setzero_pd();     
            double 
              *ev = &extEV[l * 20],
              *lv = &le[l * 20],
              *rv = &ri[l * 20];

            for(j = 0; j < 20; j+=2)
            {
              x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                  
              x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
            }

            x1v = _mm_hadd_pd(x1v, x1v);
            x2v = _mm_hadd_pd(x2v, x2v);

            x1v = _mm_mul_pd(x1v, x2v);

            for(j = 0; j < 20; j+=2)
            {
              __m128d vv = _mm_load_pd(&v[j]);
              vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
              _mm_store_pd(&v[j], vv);
            }               

          }

          {         
            __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

            scale = 1;
            for(l = 0; scale && (l < 20); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }             
          }


          if(scale)
          {

            __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

            for(l = 0; l < 20; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));                  
            }

            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];         
          }
        }
      }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
      {
        le = &left[cptr[i] * 400];
        ri = &right[cptr[i] * 400];

        vl = &x1[20 * i];
        vr = &x2[20 * i];
        v = &x3[20 * i];


        for(l = 0; l < 20; l+=2)
          _mm_store_pd(&v[l], _mm_setzero_pd());                        


        for(l = 0; l < 20; l++)
        {

          __m128d x1v = _mm_setzero_pd();
          __m128d x2v = _mm_setzero_pd();
          double 
            *ev = &extEV[l * 20],
            *lv = &le[l * 20],
            *rv = &ri[l * 20];


          for(j = 0; j < 20; j+=2)
          {
            x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                    
            x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
          }

          x1v = _mm_hadd_pd(x1v, x1v);
          x2v = _mm_hadd_pd(x2v, x2v);

          x1v = _mm_mul_pd(x1v, x2v);

          for(j = 0; j < 20; j+=2)
          {
            __m128d vv = _mm_load_pd(&v[j]);
            vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
            _mm_store_pd(&v[j], vv);
          }                 

        }

        {           
          __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

          scale = 1;
          for(l = 0; scale && (l < 20); l += 2)
          {
            __m128d vv = _mm_load_pd(&v[l]);
            __m128d v1 = _mm_and_pd(vv, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }               
        }


        if(scale)
        {

          __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

          for(l = 0; l < 20; l+=2)
          {
            __m128d ex3v = _mm_load_pd(&v[l]);            
            _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));        
          }                               


          if(!fastScaling)
            ex3[i] += 1;
          else
            addScale += wgt[i];    
        }
      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}

/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for GTR CAT with memory saving (Optimized SSE3 version for AA data)

    This is the SSE3 optimized version of ::newviewCAT_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b CAT
    model of rate heterogeneity.

    @note
    For more details and function argument description check the function ::newviewCAT_FLEX
*/
static void newviewGTRCATPROT_SAVE(int tipCase, double *extEV,
                                   int *cptr,
                                   double *x1, double *x2, double *x3, double *tipVector,
                                   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                   int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling,
                                   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
                                   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le, 
    *ri, 
    *v, 
    *vl, 
    *vr,
    *x1_ptr = x1,
    *x2_ptr = x2, 
    *x3_ptr = x3;

  int 
    i, 
    l, 
    j, 
    scale, 
    scaleGap = 0,
    addScale = 0;

  {
    vl = x1_gapColumn;        
    vr = x2_gapColumn;
    v = x3_gapColumn;

    le = &left[maxCats * 400];
    ri = &right[maxCats * 400];   

    for(l = 0; l < 20; l+=2)
      _mm_store_pd(&v[l], _mm_setzero_pd());                    

    for(l = 0; l < 20; l++)
    {
      __m128d x1v = _mm_setzero_pd();
      __m128d x2v = _mm_setzero_pd();
      double 
        *ev = &extEV[l * 20],
        *lv = &le[l * 20],
        *rv = &ri[l * 20];


      for(j = 0; j < 20; j+=2)
      {
        x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                
        x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
      }

      x1v = _mm_hadd_pd(x1v, x1v);
      x2v = _mm_hadd_pd(x2v, x2v);

      x1v = _mm_mul_pd(x1v, x2v);

      for(j = 0; j < 20; j+=2)
      {
        __m128d vv = _mm_load_pd(&v[j]);
        vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
        _mm_store_pd(&v[j], vv);
      }                 
    }

    if(tipCase != PLL_TIP_TIP)
    {       
      __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

      scale = 1;
      for(l = 0; scale && (l < 20); l += 2)
      {
        __m128d vv = _mm_load_pd(&v[l]);
        __m128d v1 = _mm_and_pd(vv, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
      }                 

      if(scale)
      {
        __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

        for(l = 0; l < 20; l+=2)
        {
          __m128d ex3v = _mm_load_pd(&v[l]);              
          _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));  
        }                                 

        scaleGap = PLL_TRUE;       
      }
    }
  }

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        for (i = 0; i < n; i++)
        {
          if(noGap(x3_gap, i))
          {             
            vl = &(tipVector[20 * tipX1[i]]);
            vr = &(tipVector[20 * tipX2[i]]);
            v  = x3_ptr;

            if(isGap(x1_gap, i))
              le =  &left[maxCats * 400];
            else                  
              le =  &left[cptr[i] * 400];         

            if(isGap(x2_gap, i))
              ri =  &right[maxCats * 400];
            else                  
              ri =  &right[cptr[i] * 400];

            for(l = 0; l < 20; l+=2)
              _mm_store_pd(&v[l], _mm_setzero_pd());                    

            for(l = 0; l < 20; l++)
            {
              __m128d x1v = _mm_setzero_pd();
              __m128d x2v = _mm_setzero_pd();    
              double 
                *ev = &extEV[l * 20],
                *lv = &le[l * 20],
                *rv = &ri[l * 20];

              for(j = 0; j < 20; j+=2)
              {
                x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                
                x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
              }

              x1v = _mm_hadd_pd(x1v, x1v);
              x2v = _mm_hadd_pd(x2v, x2v);

              x1v = _mm_mul_pd(x1v, x2v);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
                _mm_store_pd(&v[j], vv);
              }            
            }

            x3_ptr += 20;

          }   
        }
      }
      break;
    case PLL_TIP_INNER:
      {
        for (i = 0; i < n; i++)
        {
          if(isGap(x3_gap, i))
          {
            if(scaleGap)
              {
                if(!fastScaling)
                  ex3[i] += 1;
                else
                  addScale += wgt[i];
              }
          }
          else
          {      
            vl = &(tipVector[20 * tipX1[i]]);

            vr = x2_ptr;
            v = x3_ptr;

            if(isGap(x1_gap, i))
              le =  &left[maxCats * 400];
            else
              le =  &left[cptr[i] * 400];

            if(isGap(x2_gap, i))
            {            
              ri =  &right[maxCats * 400];
              vr = x2_gapColumn;
            }
            else
            {
              ri =  &right[cptr[i] * 400];
              vr = x2_ptr;
              x2_ptr += 20;
            }                                             

            for(l = 0; l < 20; l+=2)
              _mm_store_pd(&v[l], _mm_setzero_pd());                               

            for(l = 0; l < 20; l++)
            {
              __m128d x1v = _mm_setzero_pd();
              __m128d x2v = _mm_setzero_pd();   
              double 
                *ev = &extEV[l * 20],
                *lv = &le[l * 20],
                *rv = &ri[l * 20];

              for(j = 0; j < 20; j+=2)
              {
                x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                
                x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
              }

              x1v = _mm_hadd_pd(x1v, x1v);
              x2v = _mm_hadd_pd(x2v, x2v);

              x1v = _mm_mul_pd(x1v, x2v);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
                _mm_store_pd(&v[j], vv);
              }             
            }

            {       
              __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

              scale = 1;
              for(l = 0; scale && (l < 20); l += 2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d v1 = _mm_and_pd(vv, absMask.m);
                v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
                if(_mm_movemask_pd( v1 ) != 3)
                  scale = 0;
              }           
            }


            if(scale)
            {
              __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

              for(l = 0; l < 20; l+=2)
              {
                __m128d ex3v = _mm_load_pd(&v[l]);
                _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));                
              }
              
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];       
            }
            x3_ptr += 20;
          }
        }
      }
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)
            {
              if(!fastScaling)
                ex3[i] += 1;
              else
                addScale += wgt[i];
            }
        }
        else
        {                    
          v = x3_ptr;

          if(isGap(x1_gap, i))
          {
            vl = x1_gapColumn;
            le =  &left[maxCats * 400];
          }
          else
          {
            le =  &left[cptr[i] * 400];
            vl = x1_ptr;
            x1_ptr += 20;
          }

          if(isGap(x2_gap, i))  
          {
            vr = x2_gapColumn;
            ri =  &right[maxCats * 400];            
          }
          else
          {
            ri =  &right[cptr[i] * 400];
            vr = x2_ptr;
            x2_ptr += 20;
          }                               

          for(l = 0; l < 20; l+=2)
            _mm_store_pd(&v[l], _mm_setzero_pd());                      

          for(l = 0; l < 20; l++)
          {
            __m128d x1v = _mm_setzero_pd();
            __m128d x2v = _mm_setzero_pd();
            double 
              *ev = &extEV[l * 20],
              *lv = &le[l * 20],
              *rv = &ri[l * 20];

            for(j = 0; j < 20; j+=2)
            {
              x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));                  
              x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
            }

            x1v = _mm_hadd_pd(x1v, x1v);
            x2v = _mm_hadd_pd(x2v, x2v);

            x1v = _mm_mul_pd(x1v, x2v);

            for(j = 0; j < 20; j+=2)
            {
              __m128d vv = _mm_load_pd(&v[j]);
              vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
              _mm_store_pd(&v[j], vv);
            }               

          }

          {         
            __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );

            scale = 1;
            for(l = 0; scale && (l < 20); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }             
          }

          if(scale)
          {
            __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);

            for(l = 0; l < 20; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);                  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));      
            }                             

            if(!fastScaling)
              ex3[i] += 1;
            else
              addScale += wgt[i];          
          }
          x3_ptr += 20;
        }
      }
      break;
    default:
      assert(0);
  }

  if(fastScaling)
    *scalerIncrement = addScale;
}


/** @ingroup group1
 *  @brief Computation of conditional likelihood arrray for the GTR GAMMA and for the LG4 model (Optimized SSE3 version for AA data)

    This is the SSE3 optimized version of ::newviewGAMMA_FLEX for computing the conditional
    likelihood arrays at some node \a p, given child nodes \a q and \a r using the \b GAMMA
    model of rate heterogeneity and the LG4 model of evolution. Note that the original unoptimized
    function does not incorporate the LG4 model.

    @note
    For more details and function argument description check the function ::newviewGAMMA_FLEX
*/
static void newviewGTRGAMMAPROT_LG4(int tipCase,
                                    double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
                                    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                                    int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr;
#ifndef __SSE3
  double al, ar;
#endif



  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
        double umpX1[1840], umpX2[1840];

        for(i = 0; i < 23; i++)
          {
           

            for(k = 0; k < 80; k++)
              {
                
                v = &(tipVector[k / 20][20 * i]);
#ifdef __SSE3
                double *ll =  &left[k * 20];
                double *rr =  &right[k * 20];
                
                __m128d umpX1v = _mm_setzero_pd();
                __m128d umpX2v = _mm_setzero_pd();

                for(l = 0; l < 20; l+=2)
                  {
                    __m128d vv = _mm_load_pd(&v[l]);
                    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
                    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));                                   
                  }
                
                umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
                umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
                
                _mm_storel_pd(&umpX1[80 * i + k], umpX1v);
                _mm_storel_pd(&umpX2[80 * i + k], umpX2v);
#else
                umpX1[80 * i + k] = 0.0;
                umpX2[80 * i + k] = 0.0;

                for(l = 0; l < 20; l++)
                  {
                    umpX1[80 * i + k] +=  v[l] *  left[k * 20 + l];
                    umpX2[80 * i + k] +=  v[l] * right[k * 20 + l];
                  }
#endif
              }
          }

        for(i = 0; i < n; i++)
          {
            uX1 = &umpX1[80 * tipX1[i]];
            uX2 = &umpX2[80 * tipX2[i]];

            for(j = 0; j < 4; j++)
              {
                v = &x3[i * 80 + j * 20];

#ifdef __SSE3
                __m128d zero =  _mm_setzero_pd();
                for(k = 0; k < 20; k+=2)                                    
                  _mm_store_pd(&v[k], zero);

                for(k = 0; k < 20; k++)
                  { 
                    double *eev = &extEV[j][k * 20];
                    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
                    __m128d x1px2v = _mm_set1_pd(x1px2);

                    for(l = 0; l < 20; l+=2)
                      {
                        __m128d vv = _mm_load_pd(&v[l]);
                        __m128d ee = _mm_load_pd(&eev[l]);

                        vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
                        
                        _mm_store_pd(&v[l], vv);
                      }
                  }

#else

                for(k = 0; k < 20; k++)
                  v[k] = 0.0;

                for(k = 0; k < 20; k++)
                  {                
                    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
                   
                    for(l = 0; l < 20; l++)                                                     
                      v[l] += x1px2 * extEV[j][20 * k + l];                  
                  }
#endif
              }    
          }
      }
      break;
    case PLL_TIP_INNER:
      {
        double umpX1[1840], ump_x2[20];


        for(i = 0; i < 23; i++)
          {
           

            for(k = 0; k < 80; k++)
              { 
                v = &(tipVector[k / 20][20 * i]);
#ifdef __SSE3
                double *ll =  &left[k * 20];
                                
                __m128d umpX1v = _mm_setzero_pd();
                
                for(l = 0; l < 20; l+=2)
                  {
                    __m128d vv = _mm_load_pd(&v[l]);
                    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));                                                   
                  }
                
                umpX1v = _mm_hadd_pd(umpX1v, umpX1v);                           
                _mm_storel_pd(&umpX1[80 * i + k], umpX1v);              
#else       
                umpX1[80 * i + k] = 0.0;

                for(l = 0; l < 20; l++)
                  umpX1[80 * i + k] +=  v[l] * left[k * 20 + l];
#endif

              }
          }

        for (i = 0; i < n; i++)
          {
            uX1 = &umpX1[80 * tipX1[i]];

            for(k = 0; k < 4; k++)
              {
                v = &(x2[80 * i + k * 20]);
#ifdef __SSE3              
                for(l = 0; l < 20; l++)
                  {                
                    double *r =  &right[k * 400 + l * 20];
                    __m128d ump_x2v = _mm_setzero_pd();     
                    
                    for(j = 0; j < 20; j+= 2)
                      {
                        __m128d vv = _mm_load_pd(&v[j]);
                        __m128d rr = _mm_load_pd(&r[j]);
                        ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
                      }
                     
                    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
                    
                    _mm_storel_pd(&ump_x2[l], ump_x2v);                              
                  }

                v = &(x3[80 * i + 20 * k]);

                __m128d zero =  _mm_setzero_pd();
                for(l = 0; l < 20; l+=2)                                    
                  _mm_store_pd(&v[l], zero);
                  
                for(l = 0; l < 20; l++)
                  {
                    double *eev = &extEV[k][l * 20];
                    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
                    __m128d x1px2v = _mm_set1_pd(x1px2);
                  
                    for(j = 0; j < 20; j+=2)
                      {
                        __m128d vv = _mm_load_pd(&v[j]);
                        __m128d ee = _mm_load_pd(&eev[j]);
                        
                        vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
                        
                        _mm_store_pd(&v[j], vv);
                      }                             
                  }                     
#else
                for(l = 0; l < 20; l++)
                  {
                    ump_x2[l] = 0.0;

                    for(j = 0; j < 20; j++)
                      ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
                  }

                v = &(x3[80 * i + 20 * k]);

                for(l = 0; l < 20; l++)
                  v[l] = 0;

                for(l = 0; l < 20; l++)
                  {
                    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
                    for(j = 0; j < 20; j++)
                      v[j] += x1px2 * extEV[k][l * 20  + j];
                  }
#endif
              }
           
#ifdef __SSE3
            { 
              v = &(x3[80 * i]);
              __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );
              
              scale = 1;
              for(l = 0; scale && (l < 80); l += 2)
                {
                  __m128d vv = _mm_load_pd(&v[l]);
                  __m128d v1 = _mm_and_pd(vv, absMask.m);
                  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
                  if(_mm_movemask_pd( v1 ) != 3)
                    scale = 0;
                }                 
            }
#else
            v = &x3[80 * i];
            scale = 1;
            for(l = 0; scale && (l < 80); l++)
              scale = (PLL_ABS(v[l]) <  PLL_MINLIKELIHOOD );
#endif

            if (scale)
              {
#ifdef __SSE3
               __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
               
               for(l = 0; l < 80; l+=2)
                 {
                   __m128d ex3v = _mm_load_pd(&v[l]);             
                   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto)); 
                 }                                
#else
                for(l = 0; l < 80; l++)
                  v[l] *= PLL_TWOTOTHE256;
#endif

                if(useFastScaling)
                  addScale += wgt[i];
                else
                  ex3[i]  += 1;        
              }
          }
      }
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
       {
         for(k = 0; k < 4; k++)
           {
             vl = &(x1[80 * i + 20 * k]);
             vr = &(x2[80 * i + 20 * k]);
             v =  &(x3[80 * i + 20 * k]);

#ifdef __SSE3
             __m128d zero =  _mm_setzero_pd();
             for(l = 0; l < 20; l+=2)                               
               _mm_store_pd(&v[l], zero);
#else
             for(l = 0; l < 20; l++)
               v[l] = 0;
#endif

             for(l = 0; l < 20; l++)
               {                 
#ifdef __SSE3
                 {
                   __m128d al = _mm_setzero_pd();
                   __m128d ar = _mm_setzero_pd();

                   double *ll   = &left[k * 400 + l * 20];
                   double *rr   = &right[k * 400 + l * 20];
                   double *EVEV = &extEV[k][20 * l];
                   
                   for(j = 0; j < 20; j+=2)
                     {
                       __m128d lv  = _mm_load_pd(&ll[j]);
                       __m128d rv  = _mm_load_pd(&rr[j]);
                       __m128d vll = _mm_load_pd(&vl[j]);
                       __m128d vrr = _mm_load_pd(&vr[j]);
                       
                       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
                     }                   
                       
                   al = _mm_hadd_pd(al, al);
                   ar = _mm_hadd_pd(ar, ar);
                   
                   al = _mm_mul_pd(al, ar);

                   for(j = 0; j < 20; j+=2)
                     {
                       __m128d vv  = _mm_load_pd(&v[j]);
                       __m128d EVV = _mm_load_pd(&EVEV[j]);

                       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                       _mm_store_pd(&v[j], vv);
                     }                                            
                 }               
#else
                 al = 0.0;
                 ar = 0.0;

                 for(j = 0; j < 20; j++)
                   {
                     al += vl[j] * left[k * 400 + l * 20 + j];
                     ar += vr[j] * right[k * 400 + l * 20 + j];
                   }

                 x1px2 = al * ar;

                 for(j = 0; j < 20; j++)
                   v[j] += x1px2 * extEV[k][20 * l + j];
#endif
               }
           }
         

#ifdef __SSE3
         { 
           v = &(x3[80 * i]);
           __m128d minlikelihood_sse = _mm_set1_pd( PLL_MINLIKELIHOOD );
           
           scale = 1;
           for(l = 0; scale && (l < 80); l += 2)
             {
               __m128d vv = _mm_load_pd(&v[l]);
               __m128d v1 = _mm_and_pd(vv, absMask.m);
               v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
               if(_mm_movemask_pd( v1 ) != 3)
                 scale = 0;
             }            
         }
#else
         v = &(x3[80 * i]);
         scale = 1;
         for(l = 0; scale && (l < 80); l++)
           scale = ((PLL_ABS(v[l]) <  PLL_MINLIKELIHOOD ));
#endif

         if (scale)
           {
#ifdef __SSE3
               __m128d twoto = _mm_set_pd(PLL_TWOTOTHE256, PLL_TWOTOTHE256);
               
               for(l = 0; l < 80; l+=2)
                 {
                   __m128d ex3v = _mm_load_pd(&v[l]);             
                   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto)); 
                 }                                
#else        
             for(l = 0; l < 80; l++)
               v[l] *= PLL_TWOTOTHE256;
#endif

             if(useFastScaling)
               addScale += wgt[i];
             else
               ex3[i]  += 1;      
           }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}
#endif


