/*  pAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with
 *  thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <unistd.h>
#endif



#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

#ifdef __SIM_SSE3
#include <xmmintrin.h>
#include <pmmintrin.h>
/*#include <tmmintrin.h>*/
#endif

/* pointers to reduction buffers for storing and gathering the first and second derivative 
   of the likelihood in Pthreads and MPI */

#if IS_PARALLEL
void branchLength_parallelReduce(tree *tr, double *dlnLdlz,  double *d2lnLdlz2 ) ; 
extern double *globalResult;
#endif



extern const unsigned int mask32[32];

/*******************/


/* generic function to get the required pointers to the data associated with the left and right node that define a branch */

static void getVects(tree *tr, unsigned char **tipX1, unsigned char **tipX2, double **x1_start, double **x2_start, int *tipCase, int model,
    double **x1_gapColumn, double **x2_gapColumn, unsigned int **x1_gap, unsigned int **x2_gap)
{
  int    
    rateHet = (int)discreteRateCategories(tr->rateHetModel),
            states = tr->partitionData[model].states,
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

  /* switch over the different tip cases again here */

  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
  {      
    if(!( isTip(pNumber, tr->mxtips) && isTip(qNumber, tr->mxtips)) )
    {
      *tipCase = TIP_INNER;
      if(isTip(qNumber, tr->mxtips))
      {
        *tipX1 = tr->partitionData[model].yVector[qNumber];
        *x2_start = tr->partitionData[model].xVector[p_slot];

        if(tr->saveMemory)
        {
          *x2_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
          *x2_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];  
        }
      }
      else
      {
        *tipX1 = tr->partitionData[model].yVector[pNumber];
        *x2_start = tr->partitionData[model].xVector[q_slot];

        if(tr->saveMemory)
        {
          *x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
          *x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
        }
      }
    }
    else
    {
      /* note that tip tip should normally not occur since this means that we are trying to optimize 
         a branch in a two-taxon tree. However, this has been inherited be some RAxML function 
         that optimized pair-wise distances between all taxa in a tree */

      *tipCase = TIP_TIP;
      *tipX1 = tr->partitionData[model].yVector[pNumber];
      *tipX2 = tr->partitionData[model].yVector[qNumber];
    }
  }
  else
  {
    *tipCase = INNER_INNER;

    *x1_start = tr->partitionData[model].xVector[p_slot];
    *x2_start = tr->partitionData[model].xVector[q_slot];

    if(tr->saveMemory)
    {
      *x1_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
      *x1_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet]; 

      *x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
      *x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet]; 
    }
  }

}


/* this is actually a pre-computation and storage of values that remain constant while we change the value of the branch length 
   we want to adapt. the target pointer sumtable is a single pre-allocated array that has the same 
   size as a conditional likelihood vector at an inner node.

   So if we want to do a Newton-Rpahson optimization we only execute this function once in the beginning for each new branch we are considering !
   */

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

    case TIP_TIP:
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
    case TIP_INNER:

      /* same as for TIP_TIP only that 
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
    case INNER_INNER:
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



#ifndef _OPTIMIZED_FUNCTIONS

/* same thing for GAMMA models. The only noteworthy thing here is that we have an additional inner loop over the 
   number of discrete gamma rates. The data access pattern is also different since for tip vector accesses through our 
   lookup table, we do not distnguish between rates 

   Note the different access pattern in TIP_INNER:

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
    case TIP_TIP:
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
    case TIP_INNER:
      reorder_back( x2, n, span );
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
      reorder( x2, n, span );
      break;
    case INNER_INNER:
      reorder_back( x1, n, span );
      reorder_back( x2, n, span );
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
      reorder( x1, n, span );
      reorder( x2, n, span );
      break;
    default:
      assert(0);
  }
}
#endif

/* optimized functions for branch length optimization */


#ifdef _OPTIMIZED_FUNCTIONS

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

static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

static void sumGTRCATPROT_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

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

    *d_start = (double *)malloc_aligned(numberOfCategories * states * sizeof(double)),
    *e =(double *)malloc_aligned(states * sizeof(double)),
    *s = (double *)malloc_aligned(states * sizeof(double)),
    *dd = (double *)malloc_aligned(states * sizeof(double)),
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

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
      d[states * i + l] = EXP(dd[l] * rptr[i]);
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

    inv_Li = 1.0/ inv_Li;

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

  free(d_start);
  free(e);
  free(s);
  free(dd);
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
      diagptable[i * gammaStates + l * 4]     = EXP(EIGN[l] * ki * lz);
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

    inv_Li = 1.0 / inv_Li;

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz   += wrptr[i] * dlnLidlz;
    d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

}

void sumGAMMA_FLEX_reorder(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, const int states);
/* the function below is called only once at the very beginning of each Newton-Raphson procedure for optimizing barnch lengths.
   It initially invokes an iterative newview call to get a consistent pair of vectors at the left and the right end of the 
   branch and thereafter invokes the one-time only precomputation of values (sumtable) that can be re-used in each Newton-Raphson 
   iteration. Once this function has been called we can execute the actual NR procedure */

void makenewzIterative(tree *tr)
{
  int 
    model, 
    tipCase;

  double
    *x1_start = (double*)NULL,
    *x2_start = (double*)NULL;


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

  newviewIterative(tr, 1);


  /* 
     loop over all partoitions to do the precomputation of the sumTable buffer 
     This is analogous to the newviewIterative() and evaluateIterative() 
     implementations.
     */

  for(model = 0; model < tr->NumberOfModels; model++)
  { 
    int 
      width = tr->partitionData[model].width;

    if(tr->td[0].executeModel[model] && width > 0)
    {
      int 	   
        states = tr->partitionData[model].states;


      getVects(tr, &tipX1, &tipX2, &x1_start, &x2_start, &tipCase, model, &x1_gapColumn, &x2_gapColumn, &x1_gap, &x2_gap);

#ifndef _OPTIMIZED_FUNCTIONS
      assert(!tr->saveMemory);
      if(tr->rateHetModel == CAT)
        sumCAT_FLEX(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
            width, states);
      else
        sumGAMMA_FLEX_reorder(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
            width, states);
#else
      switch(states)
      {

        case 4: /* DNA */
          if(tr->rateHetModel == CAT)
          {
            if(tr->saveMemory)
              sumCAT_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumCAT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
                  width);
          }
          else
          {
            if(tr->saveMemory)
              sumGAMMA_GAPPED_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumGAMMA(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
                  width);
          }
          break;		
        case 20: /* proteins */
          if(tr->rateHetModel == CAT)
          {
            if(tr->saveMemory)
              sumGTRCATPROT_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
                  tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else	      	      
              sumGTRCATPROT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
                  tipX1, tipX2, width);
          }
          else
          {

            if(tr->saveMemory)
              sumGAMMAPROT_GAPPED_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
                  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
            else
              sumGAMMAPROT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
                  tipX1, tipX2, width);

          }
          break;		
        default:
          assert(0);
      }
#endif
    }
  }
}



/* this function actually computes the first and second derivatives of the likelihood for a given branch stored 
   in tr->coreLZ[model] Note that in the parallel case coreLZ must always be broadcasted together with the 
   traversal descriptor, at least for optimizing branch lengths */

void execCore(tree *tr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2)
{
  int model, branchIndex;
  double lz;

  _dlnLdlz[0]   = 0.0;
  _d2lnLdlz2[0] = 0.0;

  /* loop over partitions */

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    int 
      width = tr->partitionData[model].width;

    /* check if we (the present thread for instance) needs to compute something at 
       all for the present partition */

    if(tr->td[0].executeModel[model] && width > 0)
    {
      int 	    
        states = tr->partitionData[model].states;

      double 
        *sumBuffer       = (double*)NULL;


      volatile double
        dlnLdlz   = 0.0,
                  d2lnLdlz2 = 0.0;

      /* set a pointer to the part of the pre-computed sumBuffer we are going to access */

      sumBuffer = tr->partitionData[model].sumBuffer;

      /* figure out if we are optimizing branch lengths individually per partition or jointly across 
         all partitions. If we do this on a per partition basis, we also need to compute and store 
         the per-partition derivatives of the likelihood separately, otherwise not */

      if(tr->numBranches > 1)
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

#ifndef _OPTIMIZED_FUNCTIONS

      /* compute first and second derivatives with the slow generic functions */

      if(tr->rateHetModel == CAT)
        coreCAT_FLEX(width, tr->partitionData[model].numberOfCategories, sumBuffer,
            &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].wgt,
            tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz, states);
      else
        coreGAMMA_FLEX(width, sumBuffer,
            &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
            tr->partitionData[model].wgt, states);
#else
      switch(states)
      {	   
        case 4: /* DNA */
          if(tr->rateHetModel == CAT)
            coreGTRCAT(width, tr->partitionData[model].numberOfCategories, sumBuffer,
                &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].wgt, 
                tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz);
          else 
            coreGTRGAMMA(width, sumBuffer,
                &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
                tr->partitionData[model].wgt);

          break;		    
        case 20: /* proteins */
          if(tr->rateHetModel == CAT)
            coreGTRCATPROT(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
                tr->partitionData[model].rateCategory, width,
                tr->partitionData[model].wgt, 
                &dlnLdlz, &d2lnLdlz2,
                sumBuffer);
          else

            coreGTRGAMMAPROT(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
                sumBuffer, width, tr->partitionData[model].wgt,
                &dlnLdlz, &d2lnLdlz2, lz);

          break;		   
        default:
          assert(0);
      }
#endif

      /* store first and second derivative */

      _dlnLdlz[branchIndex]   = _dlnLdlz[branchIndex]   + dlnLdlz;
      _d2lnLdlz2[branchIndex] = _d2lnLdlz2[branchIndex] + d2lnLdlz2;
    }
    else
    {
      /* set to 0 to make the reduction operation consistent */

      if(width == 0 && (tr->numBranches > 1))
      {
        _dlnLdlz[model]   = 0.0;
        _d2lnLdlz2[model] = 0.0;
      }			       	    	   
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

static void topLevelMakenewz(tree *tr, double *z0, int _maxiter, double *result)
{
  double   z[NUM_BRANCHES], zprev[NUM_BRANCHES], zstep[NUM_BRANCHES];
  volatile double  dlnLdlz[NUM_BRANCHES], d2lnLdlz2[NUM_BRANCHES];
  int i, maxiter[NUM_BRANCHES], model;
  boolean firstIteration = TRUE;
  boolean outerConverged[NUM_BRANCHES];
  boolean loopConverged;


  /* figure out if this is on a per partition basis or jointly across all partitions */



  /* initialize loop convergence variables etc. 
     maxiter is the maximum number of NR iterations we are going to do before giving up */

  for(i = 0; i < tr->numBranches; i++)
  {
    z[i] = z0[i];
    maxiter[i] = _maxiter;
    outerConverged[i] = FALSE;
    tr->curvatOK[i]       = TRUE;
  }


  /* nested do while loops of Newton-Raphson */

  do
  {

    /* check if we ar done for partition i or if we need to adapt the branch length again */

    for(i = 0; i < tr->numBranches; i++)
    {
      if(outerConverged[i] == FALSE && tr->curvatOK[i] == TRUE)
      {
        tr->curvatOK[i] = FALSE;

        zprev[i] = z[i];

        zstep[i] = (1.0 - zmax) * z[i] + zmin;
      }
    }

    for(i = 0; i < tr->numBranches; i++)
    {
      /* other case, the outer loop hasn't converged but we are trying to approach 
         the maximum from the wrong side */

      if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
      {
        double lz;

        if (z[i] < zmin) z[i] = zmin;
        else if (z[i] > zmax) z[i] = zmax;
        lz    = log(z[i]);

        tr->coreLZ[i] = lz;
      }
    }


    /* set the execution mask */

    if(tr->numBranches > 1)
    {
      for(model = 0; model < tr->NumberOfModels; model++)
      {
        if(tr->executeModel[model])
          tr->executeModel[model] = !tr->curvatOK[model];

      }
    }
    else
    {
      for(model = 0; model < tr->NumberOfModels; model++)
        tr->executeModel[model] = !tr->curvatOK[0];
    }


    /* store it in traversal descriptor */

    storeExecuteMaskInTraversalDescriptor(tr); 

    /* store the new branch length values to be tested in traversal descriptor */

    storeValuesInTraversalDescriptor(tr, &(tr->coreLZ[0]));

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

    /* if this is the first iteration of NR we will need to first do this one-time call 
       of maknewzIterative() Note that, only this call requires broadcasting the traversal descriptor,
       subsequent calls to masterBarrier(THREAD_MAKENEWZ, tr); will not require this
       */

    if(firstIteration)
      {
	tr->td[0].traversalHasChanged = TRUE; 
	masterBarrier(THREAD_MAKENEWZ_FIRST, tr); 
	firstIteration = FALSE; 
	tr->td[0].traversalHasChanged = FALSE; 
      }
    else 
      masterBarrier(THREAD_MAKENEWZ, tr); 
    branchLength_parallelReduce(tr, (double*)dlnLdlz, (double*)d2lnLdlz2); 
#else 
    /* sequential part, if this is the first newton-raphson implementation,
       do the precomputations as well, otherwise just execute the computation
       of the derivatives */
    if(firstIteration)
      {
	makenewzIterative(tr);
	firstIteration = FALSE;
      }
    execCore(tr, dlnLdlz, d2lnLdlz2);
#endif

    /* do a NR step, if we are on the correct side of the maximum that's okay, otherwise 
       shorten branch */

    for(i = 0; i < tr->numBranches; i++)
    {
      if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
      {
        if ((d2lnLdlz2[i] >= 0.0) && (z[i] < zmax))
          zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
        else
          tr->curvatOK[i] = TRUE;
      }
    }

    /* do the standard NR step to obrain the next value, depending on the state for eahc partition */

    for(i = 0; i < tr->numBranches; i++)
    {
      if(tr->curvatOK[i] == TRUE && outerConverged[i] == FALSE)
      {
        if (d2lnLdlz2[i] < 0.0)
        {
          double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
          if (tantmp < 100)
          {
            z[i] *= EXP(tantmp);
            if (z[i] < zmin)
              z[i] = zmin;

            if (z[i] > 0.25 * zprev[i] + 0.75)
              z[i] = 0.25 * zprev[i] + 0.75;
          }
          else
            z[i] = 0.25 * zprev[i] + 0.75;
        }
        if (z[i] > zmax) z[i] = zmax;

        /* decrement the maximum number of itarations */

        maxiter[i] = maxiter[i] - 1;

        /* check if the outer loop has converged */
        if(maxiter[i] > 0 && (ABS(z[i] - zprev[i]) > zstep[i]))
          outerConverged[i] = FALSE;
        else
          outerConverged[i] = TRUE;
      }
    }

    /* check if the loop has converged for all partitions */

    loopConverged = TRUE;
    for(i = 0; i < tr->numBranches; i++)
      loopConverged = loopConverged && outerConverged[i];
  }
  while (!loopConverged);


  /* reset  partition execution mask */

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->executeModel[model] = TRUE;

  /* copy the new branches in the result array of branches.
     if we don't do a per partition estimate of 
     branches this will only set result[0]
     */

  for(i = 0; i < tr->numBranches; i++)
    result[i] = z[i];
}

/* function called from RAxML to optimize a given branch with current branch lengths z0 
   between nodes p and q.
   The new branch lengths will be stored in result */

void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask)
{
  int i;
  boolean originalExecute[NUM_BRANCHES];

  boolean 
    p_recom = FALSE, /* if one of was missing, we will need to force recomputation */
    q_recom = FALSE;

  /* the first entry of the traversal descriptor stores the node pair that defines 
     the branch */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)
  {
    originalExecute[i] =  tr->executeModel[i];
    tr->td[0].ti[0].qz[i] =  z0[i];
    if(mask)
    {
      if(tr->partitionConverged[i])
        tr->executeModel[i] = FALSE;
      else
        tr->executeModel[i] = TRUE;
    }
  }
  if(tr->useRecom)
  {
    int
      slot = -1,
      count = 0;

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
    computeTraversal(tr, p, TRUE);

  if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
    computeTraversal(tr, q, TRUE);

  /* call the Newton-Raphson procedure */

  topLevelMakenewz(tr, z0, maxiter, result);

  /* Mark node as unpinnable */
  if(tr->useRecom)
  {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
  }

  /* fix eceuteModel this seems to be a bit redundant with topLevelMakenewz */ 

  for(i = 0; i < tr->numBranches; i++)
    tr->executeModel[i] = TRUE;
}


/* below are, once again the optimized functions */

#ifdef _OPTIMIZED_FUNCTIONS
static inline boolean isGap(unsigned int *x, int pos)
{
  return (x[pos / 32] & mask32[pos % 32]);
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
    case TIP_TIP:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case TIP_INNER:
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
    case INNER_INNER:
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
    case TIP_TIP:     
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
    case TIP_INNER:
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
    case INNER_INNER:
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
    case TIP_TIP:
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
    case TIP_INNER:
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
    case INNER_INNER:
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
    case TIP_TIP:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &x2_start[4 * i];

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case INNER_INNER:
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
    case TIP_TIP:
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
    case TIP_INNER:
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
    case INNER_INNER:
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


static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
  {
    case TIP_TIP:
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
    case TIP_INNER:
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
    case INNER_INNER:
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
    case TIP_TIP:
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
    case TIP_INNER:
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
    case INNER_INNER:
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
    case TIP_TIP:
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
    case TIP_INNER:
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
    case INNER_INNER:
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
            *sum, 
            diagptable0[16] __attribute__ ((aligned (BYTE_ALIGNMENT))),
            diagptable1[16] __attribute__ ((aligned (BYTE_ALIGNMENT))),
            diagptable2[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));    

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
      diagptable0[i * 4 + l] = EXP(EIGN[l] * ki * lz);
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

    inv_Li = 1.0 / inv_Li;

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;     

    dlnLdlz   += wrptr[i] * dlnLidlz;
    d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }


  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}



static void coreGTRCAT(int upper, int numberOfCategories, double *sum,
    volatile double *d1, volatile double *d2, int *wgt,
    double *rptr, double *EIGN, int *cptr, double lz)
{
  int i;
  double
    *d, *d_start,
    inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e1[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double e2[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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

  d = d_start = (double *)malloc_aligned(numberOfCategories * 4 * sizeof(double));

  dd1 = EIGN[1] * lz;
  dd2 = EIGN[2] * lz;
  dd3 = EIGN[3] * lz;

  for(i = 0; i < numberOfCategories; i++)
  {
    d[i * 4 + 0] = 1.0;
    d[i * 4 + 1] = EXP(dd1 * rptr[i]);
    d[i * 4 + 2] = EXP(dd2 * rptr[i]);
    d[i * 4 + 3] = EXP(dd3 * rptr[i]);
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

    inv_Li = 1.0/inv_Li;

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz  += wgt[i] * rptr[cptr[i]] * dlnLidlz;
    d2lnLdlz2 += wgt[i] * rptr[cptr[i]] * rptr[cptr[i]] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(d_start);
}




static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, 
          diagptable0[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
          diagptable1[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
          diagptable2[80] __attribute__ ((aligned (BYTE_ALIGNMENT)));    
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
      diagptable0[i * 20 + l] = EXP(EIGN[l] * ki * lz);
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

    inv_Li = 1.0 / inv_Li;

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
  double *d1, *d_start, *sum;
  double 
    e[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
    s[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
    dd[20] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)malloc_aligned(numberOfCategories * 20 * sizeof(double));

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
      d1[20 * i + l] = EXP(dd[l] * rptr[i]);
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

    inv_Li = 1.0/inv_Li;

    dlnLidlz   *= inv_Li;
    d2lnLidlz2 *= inv_Li;

    dlnLdlz  += wgt[i] * rptr[cptr[i]] * dlnLidlz;
    d2lnLdlz2 += wgt[i] * rptr[cptr[i]] * rptr[cptr[i]] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
  }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}




#endif



