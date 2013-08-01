/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 
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
#endif


/* optimized implementation for computing per-site log likelihoods under CAT and GAMMA for DNA and protein data */

#ifdef _OPTIMIZED_FUNCTIONS
static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   unsigned  char **yVector, int mxtips);

static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, unsigned char **yVector, 
					int branchReference, int mxtips);

static inline void computeVectorGTRGAMMAPROT(double *lVector, int *eVector, double *gammaRates, int i, double qz, double rz,
					     traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					     unsigned  char **yVector, int mxtips);

static double evaluatePartialGTRGAMMAPROT(int i, int counter,  traversalInfo *ti, double qz,
					  int w, double *EIGN, double *EI, double *EV,
					  double *tipVector, unsigned char **yVector, 
					  double *gammaRates,
					  int branchReference, int mxtips);

static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips);

static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips);

static double evaluatePartialGTRGAMMA(int i, int counter,  traversalInfo *ti, double qz,
				      int w, double *EIGN, double *EI, double *EV,
				      double *tipVector, unsigned char **yVector, 
				      double *gammaRates,
				      int branchReference, int mxtips);
#endif

/* the next two functions are generic non-optimized versions of the per-site log likelihood calculations,
   but only under the CAT model. There are no generic implementations available for GAMMA yet, since 
   these functions were not needed in RAxML. However there exist optimized functions for GAMMA further below.
   The only use of the CAT functions was to optimize per-site rates based on their likelihood for the CAT 
   model of rate heterogeneity. */


static inline void computeVectorCAT_FLEX(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					 traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					 unsigned char **yVector, int mxtips, const int states)
{      
  /* allocate some space we need */
 
  double  
    *d1 =    (double *)malloc(sizeof(double) * states), 
    *d2 =    (double *)malloc(sizeof(double) * states),  
    *x1px2 = (double *)malloc(sizeof(double) * states), 
    ump_x1, 
    ump_x2,    
    lz1, 
    lz2,
    *x1, 
    *x2, 
    *x3;
  
  int 
    scale,
    j, 
    k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  /* 
     lVector holds the space for computing ancestral probablities on a single column of the tree 
     hence under CAT we index the current space required to store the parent ancestral probability vector 
     by multiplying the number of states with the offset in the array given by the inner node number
   */

  x3  = &lVector[states * (pNumber  - mxtips)];  
 
  /* do a case switch to figure out how to index the child nodes x1 and x2,
     analogous to the standard newview implementation.
     Note the index i that we use to index the specific tip poistion/index 
     for which we want to compute the per-site log likelihood */

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[states * yVector[qNumber][i]]);
      x2 = &(tipVector[states * yVector[rNumber][i]]);    
      break;
    case TIP_INNER:     
      x1 = &(tipVector[states * yVector[qNumber][i]]);
      x2 = &(lVector[states * (rNumber - mxtips)]);           
      break;
    case INNER_INNER:            
      x1 = &(lVector[states * (qNumber - mxtips)]);
      x2 = &(lVector[states * (rNumber - mxtips)]);     
      break;
    default:
      assert(0);
    }
     
  /* multiply the branch lengths with the evolutionary rate */

  lz1 = qz * ki;  
  lz2 = rz * ki;
  

  /* exponentiate the branch lengths using the eigenvalues */

  d1[0] = x1[0];
  d2[0] = x2[0];


  for(j = 1; j < states; j++)
    {
      d1[j] = x1[j] * EXP(EIGN[j] * lz1);
      d2[j] = x2[j] * EXP(EIGN[j] * lz2);	    
    }
 
 
  /* now loop over all states */

  for(j = 0; j < states; j++)
    {         
      ump_x1 = 0.0;
      ump_x2 = 0.0;

      for(k = 0; k < states; k++)
	{
	  ump_x1 += d1[k] * EI[j * states + k];
	  ump_x2 += d2[k] * EI[j * states + k];
	}
      
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < states; j++)
    x3[j] = 0.0;

  /* multiply the result of looping over all states with the eigenvector matrix EV */

  for(j = 0; j < states; j++)          
    for(k = 0; k < states; k++)	
      x3[k] +=  x1px2[j] *  EV[states * j + k];	   
      
  /* now determine if we need to scale the #states entries in x[3] to avoid 
     numerical underflow. */
     

  scale = 1;
  for(j = 0; scale && (j < states); j++)
    scale = ((x3[j] < minlikelihood) && (x3[j] > minusminlikelihood));
  
  /* if we need to scale, we multiply all probabilities of the site with 2^256 
     and increment the scaling counter by 1. 
     The counter eVector is used for tracking/counting the number of scaling events 
     at the site i for which we are computing the per-site log likelihood such that 
     we can "undo" the scaling multiplications when we compute the log likelihood of the site 
     at the virtual root */
  
  if(scale)
    {
      for(j = 0; j < states; j++)
	x3[j] *= twotothe256;       
      *eVector = *eVector + 1;
    }	              

  free(d1);
  free(d2);
  free(x1px2);
       
  return;
}


/* the following function computes the per-site log likelihood of a given site i at the virtual root of the tree.
   as input it takes the indeix i, of the site, the evolutionary rate ki (for computing Q^(rt) where r = ki) 
   the traversalDescriptor defining the full tree traversal (felsenstein pruning algo) 
   the branch length at the root qz, the weigth of the site pattern w, i.e., how many identical sites have been compressed 
   into the current site pattern, the eigenvalues etc (EIGN, EI, EV) associated to the Eigenvector/Eigenvalue decomposition 
   of the given instataneous substitution matrix Q, the tipVector lookup table for obtaining tip probability vectors, 
   a pointer to the raw sequence data at the tips, a branch index (to get the correct branch length/index into the correct branch 
   if -M is used, i.e., a per-partition branch length estimate is deployed, and finally the maximum number of tips in the comprehensive tree 
   as well as the number of states in the current model. */

static double evaluatePartialCAT_FLEX(int i, double ki, int counter,  traversalInfo *ti, double qz,
				      int w, double *EIGN, double *EI, double *EV,
				      double *tipVector, unsigned  char **yVector, 
				      int branchReference, int mxtips, const int states)
{
  int 
    scale = 0, 
    k;
  
  double 
    /* lVector is a temporary buffer to store the ancestral probability vactors of 
       a single site, thus we allocate states * mxtips space for storing probability values.
       Essentially  only (states * (mxtips - 2)) space would be required, but I was to lazy 
       to think if it has to be -1 or -2 here */
       
    *lVector = (double *)malloc_aligned(sizeof(double) * states * mxtips),
    *d = (double *)malloc_aligned(sizeof(double) * states),
    lz, 
    term, 
    *x1, 
    *x2; 

  

  traversalInfo 
    *trav = &ti[0];
 
  /* make sure that at one end of the branch into which we have placed the virtual root 
     there actually is a tip!*/

  assert(isTip(trav->pNumber, mxtips));
     
  /* for the tip we alread have the data, so just set the left probability vector to the 
     corresponding address in the pre-computed tipVector[] lookup table */

  x1 = &(tipVector[states *  yVector[trav->pNumber][i]]);   

  /* now iterate over the traversal descriptor that contains the nodes of the tree in the order required 
     by the Felsenstein pruning algorithm */

  for(k = 1; k < counter; k++)    
    {
      /* obtain the branch lengths and take the logarithms */
      
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      /* invoke essentially a newview() for one site on the entry k of the traversal descriptor.
	 counter should always correspond to the number of inner nodes in the tree for which we need
	 to compute ancestral probability values */

      computeVectorCAT_FLEX(lVector, &scale, ki, i, qz, rz, &ti[k], 
			    EIGN, EI, EV, 
			    tipVector, yVector, mxtips, states);       
    }
   
  /* now the ancestral probability values for site i at the node to the right of the virtual root 
     are available and correctly computed, such that we can set the pointer to the right vector x2
     to the corresponding entry */

  x2 = &lVector[states * (trav->qNumber - mxtips)]; 

  /* a paranoic assertion */

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
 
  /* now just compute the log likelihood score of this site */
      
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = 1.0; 

  for(k = 1; k < states; k++)
    d[k] = EXP (EIGN[k] * lz);
  
  term = 0.0;

  for(k = 0; k < states; k++) 
    term += x1[k] * x2[k] * d[k];       

  /* note the "scale * LOG(minlikelihood)" term here which we use to undo/revert the scaling multiplications 
     such that we obtain a correct log likelihood score. The integer variable scale, contains the number of times 
     we had to scale (multiply by 2^256) for site i only during a full tree traversal using Felsenstein's algorithm */

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  /* multiply with the site pattern weight (site pattern compression factor */

  term = term * w;

  /* free the memory space used for likelihood computations on this site */

  free(lVector);  
  free(d);

  return  term;
}

/* this is the top-level function that can be called from other parts of the code.
   As input it takes the tree data structure, the site index, the evolutionary rate ki, 
   and the model index (partition index. It will return the 
   log likelihood of site i. 
   An important pre-condition is that the tree traversal descriptor must contain 
   a full tree traversal starting at a tip !

   Note that, if you wamt to obtain per-site log likes for other altered model parameters such 
   as the Q matrix, you will have do re-invoke the eigenvalue/eigenvector decomposition prior 
   to calling the function below.
*/

double evaluatePartialGeneric (tree *tr, int i, double ki, int _model)
{
  double 
    result;
  
  
  int     
    branchReference,

    /* number of states of the data type in this partition */
    states = tr->partitionData[_model].states;
    
  /* SOS ATTENTION: note the different indexing used for the parallel and sequential versions ! */

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  int index = i; 
#else
  int index = i - tr->partitionData[_model].lower;
#endif
  
  /* here we figure out if all partitions are linked via the same branch length, that is,
     if we are conducting a joint branch length estimate or a per-partition branch length estimate */

  if(tr->numBranches > 1)
    branchReference = _model;
  else
    branchReference = 0;

  /* for the generic function implementation we only offer the CAT implementation for computing/optimizing per-site evolutionary rates */

#ifndef _OPTIMIZED_FUNCTIONS
  if(tr->rateHetModel == CAT)
    result = evaluatePartialCAT_FLEX(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
				     tr->partitionData[_model].wgt[index],
				     tr->partitionData[_model].EIGN, 
				     tr->partitionData[_model].EI, 
				     tr->partitionData[_model].EV,
				     tr->partitionData[_model].tipVector,
				     tr->partitionData[_model].yVector, branchReference, tr->mxtips, states);
  else
    /* 
       the per-site site likelihood function should only be called for the CAT model
       under the GAMMA model this is required only for estimating per-site protein models 
       which has however been removed in this version of the code
    */
    assert(0); 
  
 
#else
  /* switch over the number of states of the data in the current model/partition */
  switch(states)
    {
   
    case 4:   /* DNA */
      /* switch over CAT versus GAMMA and pass all model parameters for the respective partition to the respective functions */
      if(tr->rateHetModel == CAT)      
	result = evaluatePartialGTRCAT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
				       tr->partitionData[_model].wgt[index],
				       tr->partitionData[_model].EIGN, 
				       tr->partitionData[_model].EI, 
				       tr->partitionData[_model].EV,
				       tr->partitionData[_model].tipVector,
				       tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      else	
	result = evaluatePartialGTRGAMMA(index, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					 tr->partitionData[_model].wgt[index],
					 tr->partitionData[_model].EIGN, 
					 tr->partitionData[_model].EI, 
					 tr->partitionData[_model].EV,
					 tr->partitionData[_model].tipVector, 
					 tr->partitionData[_model].yVector, 
					 tr->partitionData[_model].gammaRates,
					 branchReference, tr->mxtips);	
	
      break;
    case 20: /* proteins */     
      if(tr->rateHetModel == CAT)
	result = evaluatePartialGTRCATPROT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					   tr->partitionData[_model].wgt[index],
					   tr->partitionData[_model].EIGN, 
					   tr->partitionData[_model].EI, 
					   tr->partitionData[_model].EV,
					   tr->partitionData[_model].tipVector, 
					   tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      else
	result =  evaluatePartialGTRGAMMAPROT(index, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					      tr->partitionData[_model].wgt[index],
					      tr->partitionData[_model].EIGN, 
					      tr->partitionData[_model].EI, 
					      tr->partitionData[_model].EV,
					      tr->partitionData[_model].tipVector, 
					      tr->partitionData[_model].yVector, 
					      tr->partitionData[_model].gammaRates,
					      branchReference, tr->mxtips);
      break;   
    default:
      assert(0);
    }
  #endif
 

  return result;
}

#ifdef _OPTIMIZED_FUNCTIONS

/* optimized function implementations for computing per-site log likelihoods under CAT and GAMMA for protein and 
   DNA data. 
   The structure is analoguous as above with some data- and model-specific optimizations and vectorizations.
*/

static inline void computeVectorGTRGAMMAPROT(double *lVector, int *eVector, double *gammaRates, int i, double qz, double rz,
					     traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					     unsigned  char **yVector, int mxtips)
{       
  double   
    *x1, 
    *x2, 
    *x3;  
  
  int
    s,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber,
    index1[4],
    index2[4];
  
 
  x3  = &(lVector[80 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);     
      for(s = 0; s < 4; s++)
	{
	  index1[s] = 0;
	  index2[s] = 0;
	}
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[80 * (rNumber - mxtips)]);   
      for(s = 0; s < 4; s++)       
	index1[s] = 0;
      for(s = 0; s < 4; s++)     
	index2[s] = s;                     
      break;
    case INNER_INNER:            
      x1 = &(lVector[80 * (qNumber - mxtips)]);
      x2 = &(lVector[80 * (rNumber - mxtips)]); 
      for(s = 0; s < 4; s++)
	{
	  index1[s] = s;
	  index2[s] = s;
	}                
      break;    
    default:
      assert(0);
    }
     
  {
    double  
      e1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      e2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      d1[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      d2[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      lz1, lz2;  
    
    int 
      l, 
      k, 
      scale, 
      j;
     
    for(j = 0; j < 4; j++)
      {
	lz1 = qz * gammaRates[j];            
	lz2 = rz * gammaRates[j];        

	e1[0] = 1.0;
	e2[0] = 1.0;
    
	for(l = 1; l < 20; l++)
	  {
	    e1[l] = EXP(EIGN[l] * lz1);
	    e2[l] = EXP(EIGN[l] * lz2);
	  }

	for(l = 0; l < 20; l+=2)
	  {
	    __m128d d1v = _mm_mul_pd(_mm_load_pd(&x1[20 * index1[j] + l]), _mm_load_pd(&e1[l]));
	    __m128d d2v = _mm_mul_pd(_mm_load_pd(&x2[20 * index2[j] + l]), _mm_load_pd(&e2[l]));
	    
	    _mm_store_pd(&d1[l], d1v);
	    _mm_store_pd(&d2[l], d2v);	
	  }

	__m128d zero = _mm_setzero_pd();

	for(l = 0; l < 20; l+=2)
	  _mm_store_pd(&x3[j * 20 + l], zero);
                
	for(l = 0; l < 20; l++)
	  { 	      
	    double *ev = &EV[l * 20];
	    __m128d ump_x1v = _mm_setzero_pd();
	    __m128d ump_x2v = _mm_setzero_pd();
	    __m128d x1px2v;
	    
	    for(k = 0; k < 20; k+=2)
	      {       
		__m128d eiv = _mm_load_pd(&EI[20 * l + k]);
		__m128d d1v = _mm_load_pd(&d1[k]);
		__m128d d2v = _mm_load_pd(&d2[k]);
		ump_x1v = _mm_add_pd(ump_x1v, _mm_mul_pd(d1v, eiv));
		ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(d2v, eiv));	  
	      }

	    ump_x1v = _mm_hadd_pd(ump_x1v, ump_x1v);
	    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

	    x1px2v = _mm_mul_pd(ump_x1v, ump_x2v);

	    for(k = 0; k < 20; k+=2)
	      {
		__m128d ex3v = _mm_load_pd(&x3[j * 20 + k]);
		__m128d EVV  = _mm_load_pd(&ev[k]);
		ex3v = _mm_add_pd(ex3v, _mm_mul_pd(x1px2v, EVV));
		
		_mm_store_pd(&x3[j * 20 + k], ex3v);	   	   
	      }
	  }        
      }
    
    scale = 1;
    for(l = 0; scale && (l < 80); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	      
	__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

	for(l = 0; l < 80; l+=2)
	  {
	    __m128d ex3v = _mm_mul_pd(_mm_load_pd(&x3[l]),twoto);
	    _mm_store_pd(&x3[l], ex3v);	
	  }

	*eVector = *eVector + 1;
      }
    
    return;      
  }
}

static  void computeVectorGTRGAMMA(double *lVector, int *eVector, double *gammaRates, int i, double qz, double rz,
					 traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					 unsigned  char **yVector, int mxtips)
{       
  double   
    *x1, 
    *x2, 
    *x3;   

  int
    s,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber,
    index1[4],
    index2[4];
  
 
  x3  = &(lVector[16 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:          
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(tipVector[4 * yVector[rNumber][i]]);     
      
      for(s = 0; s < 4; s++)
	{
	  index1[s] = 0;
	  index2[s] = 0;
	}
      break;
    case TIP_INNER:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(lVector[16 * (rNumber - mxtips)]);   
      for(s = 0; s < 4; s++)       
	{
	  index1[s] = 0;      
	  index2[s] = s;  
	}
      break;
    case INNER_INNER:            
      x1 = &(lVector[16 * (qNumber - mxtips)]);
      x2 = &(lVector[16 * (rNumber - mxtips)]);       
      for(s = 0; s < 4; s++)
	{
	  index1[s] = s;
	  index2[s] = s;
	}                
      break;    
    default:
      assert(0);
    }
     
  {
    double  
      e1[4] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      e2[4] __attribute__ ((aligned (BYTE_ALIGNMENT))),
      d1[4] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      d2[4] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
      lz1, lz2;  
    
    int 
      l, 
      k, 
      scale, 
      j;
     
    for(j = 0; j < 4; j++)
      {
	lz1 = qz * gammaRates[j];            
	lz2 = rz * gammaRates[j];        

	e1[0] = 1.0;
	e2[0] = 1.0;
    
	for(l = 1; l < 4; l++)
	  {
	    e1[l] = EXP(EIGN[l] * lz1);
	    e2[l] = EXP(EIGN[l] * lz2);
	  }

	for(l = 0; l < 4; l+=2)
	  {
	    __m128d d1v = _mm_mul_pd(_mm_load_pd(&x1[4 * index1[j] + l]), _mm_load_pd(&e1[l]));
	    __m128d d2v = _mm_mul_pd(_mm_load_pd(&x2[4 * index2[j] + l]), _mm_load_pd(&e2[l]));
	    
	    _mm_store_pd(&d1[l], d1v);
	    _mm_store_pd(&d2[l], d2v);	
	  }

	__m128d zero = _mm_setzero_pd();

	for(l = 0; l < 4; l+=2)
	  _mm_store_pd(&x3[j * 4 + l], zero);
                
	for(l = 0; l < 4; l++)
	  { 	      
	    double *ev = &EV[l * 4];
	    __m128d ump_x1v = _mm_setzero_pd();
	    __m128d ump_x2v = _mm_setzero_pd();
	    __m128d x1px2v;
	    
	    for(k = 0; k < 4; k+=2)
	      {       
		__m128d eiv = _mm_load_pd(&EI[4 * l + k]);
		__m128d d1v = _mm_load_pd(&d1[k]);
		__m128d d2v = _mm_load_pd(&d2[k]);
		ump_x1v = _mm_add_pd(ump_x1v, _mm_mul_pd(d1v, eiv));
		ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(d2v, eiv));	  
	      }

	    ump_x1v = _mm_hadd_pd(ump_x1v, ump_x1v);
	    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

	    x1px2v = _mm_mul_pd(ump_x1v, ump_x2v);

	    for(k = 0; k < 4; k+=2)
	      {
		__m128d ex3v = _mm_load_pd(&x3[j * 4 + k]);
		__m128d EVV  = _mm_load_pd(&ev[k]);
		ex3v = _mm_add_pd(ex3v, _mm_mul_pd(x1px2v, EVV));
		
		_mm_store_pd(&x3[j * 4 + k], ex3v);	   	   
	      }
	  }        
      }
    
  
    scale = 1;
    for(l = 0; scale && (l < 16); l++)
      scale = (ABS(x3[l]) < minlikelihood);	       	      	      	       	       
    
    if(scale)
      {	      
	__m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	
	for(l = 0; l < 16; l+=2)
	  {
	    __m128d ex3v = _mm_mul_pd(_mm_load_pd(&x3[l]),twoto);
	    _mm_store_pd(&x3[l], ex3v);	
	  }
	
	*eVector = *eVector + 1;
      }  
    
    return;      
  }
}


static double evaluatePartialGTRGAMMAPROT(int i, int counter,  traversalInfo *ti, double qz,
					  int w, double *EIGN, double *EI, double *EV,
					  double *tipVector, unsigned char **yVector, 
					  double *gammaRates,
					  int branchReference, int mxtips)
{
  double lz, term;       
  double  d[80];
  double   *x1, *x2; 
  int scale = 0, k, l, j;
  double 
    *lVector = (double *)malloc_aligned(sizeof(double) * 80 * mxtips),
    myEI[400]  __attribute__ ((aligned (BYTE_ALIGNMENT)));

  traversalInfo 
    *trav = &ti[0];

  for(k = 0; k < 20; k++)
    {         
      for(l = 0; l < 20; l++)
	myEI[k * 20 + l] = EI[k * 20 + l];
    }

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRGAMMAPROT(lVector, &scale, gammaRates, i, qz, rz, 
				&ti[k], EIGN, myEI, EV, 
				tipVector, yVector, mxtips);
    }
   
  x2 = &lVector[80 * (trav->qNumber - mxtips)];       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  lz = qz;

  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz);
  
  
  
  for(j = 0; j < 4; j++)
    {
      d[20 * j] = 1.0;
      for(l = 1; l < 20; l++)
	d[20 * j + l] = EXP(EIGN[l] * lz * gammaRates[j]);
    }

 
  for(j = 0, term = 0.0; j < 4; j++)
    {
      for(l = 0; l < 20; l++)
	term += x1[l] * x2[20 * j + l] * d[j * 20 + l];	      
    }
  
  term = LOG(0.25 * FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
 
  return  term;
}

static double evaluatePartialGTRGAMMA(int i, int counter,  traversalInfo *ti, double qz,
				      int w, double *EIGN, double *EI, double *EV,
				      double *tipVector, unsigned char **yVector, 
				      double *gammaRates,
				      int branchReference, int mxtips)
{
  double lz, term;       
  double  d[16];
  double   *x1, *x2; 
  int scale = 0, k, l, j;
  double 
    *lVector = (double *)malloc_aligned(sizeof(double) * 16 * mxtips),
    myEI[16]  __attribute__ ((aligned (BYTE_ALIGNMENT)));

  traversalInfo 
    *trav = &ti[0];

  for(k = 0; k < 4; k++)
    {           
      for(l = 0; l < 4; l++)
	myEI[k * 4 + l] = EI[k * 4 + l];
    }

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[4 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRGAMMA(lVector, &scale, gammaRates, i, qz, rz, 
				&ti[k], EIGN, myEI, EV, 
				tipVector, yVector, mxtips);
    }
   
  x2 = &lVector[16 * (trav->qNumber - mxtips)];       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  
  for(j = 0; j < 4; j++)
    {
      d[4 * j] = 1.0;
      for(l = 1; l < 4; l++)
	d[4 * j + l] = EXP(EIGN[l] * lz * gammaRates[j]);
    }

 
  for(j = 0, term = 0.0; j < 4; j++)
    {
      for(l = 0; l < 4; l++)
	term += x1[l] * x2[4 * j + l] * d[j * 4 + l];	      
    }

  term = LOG(0.25 * FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
  
  return  term;
}




static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips)
{       
  double  d1[3], d2[3],  ump_x1, ump_x2, x1px2[4], lz1, lz2; 
  double *x1, *x2, *x3;
  int j, k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[4 * (pNumber  - mxtips)];  
 

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(tipVector[4 * yVector[rNumber][i]]);    
      break;
    case TIP_INNER:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &lVector[4 * (rNumber - mxtips)];           
      break;
    case INNER_INNER:            
      x1 = &lVector[4 * (qNumber - mxtips)];
      x2 = &lVector[4 * (rNumber - mxtips)];     
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
  for(j = 0; j < 3; j++)
    {
      d1[j] = 
	x1[j + 1] * 
	EXP(EIGN[j + 1] * lz1);
      d2[j] = x2[j + 1] * EXP(EIGN[j + 1] * lz2);	    
    }
 
 
  for(j = 0; j < 4; j++)
    {     
      ump_x1 = x1[0];
      ump_x2 = x2[0];
      for(k = 0; k < 3; k++)
	{
	  ump_x1 += d1[k] * EI[j * 4 + k + 1];
	  ump_x2 += d2[k] * EI[j * 4 + k + 1];
	}
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 4; j++)
    x3[j] = 0.0;

  for(j = 0; j < 4; j++)          
    for(k = 0; k < 4; k++)	
      x3[k] +=  x1px2[j] *  EV[4 * j + k];	   
      
  
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
      x3[3] < minlikelihood && x3[3] > minusminlikelihood)
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;
      x3[2]   *= twotothe256;     
      x3[3]   *= twotothe256;     
      *eVector = *eVector + 1;
    }	              

  return;
}








static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips)
{
  double lz, term;       
  double  d[3];
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc_aligned(sizeof(double) * 4 * mxtips);    

  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[4 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)    
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRCAT(lVector, &scale, ki, i, qz, rz, &ti[k], 
			  EIGN, EI, EV, 
			  tipVector, yVector, mxtips);       
    }
   
  x2 = &lVector[4 * (trav->qNumber - mxtips)]; 

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = EXP (EIGN[1] * lz);
  d[1] = EXP (EIGN[2] * lz);
  d[2] = EXP (EIGN[3] * lz);       	   
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d[0];
  term += x1[2] * x2[2] * d[1];
  term += x1[3] * x2[3] * d[2];     

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);  

  return  term;
}

/**********************************************************************************/

static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips)
{       
  double  d1[20], d2[20],  ump_x1, ump_x2, x1px2[20], lz1, lz2; 
  double *x1, *x2, *x3;
  int j, k,
    scale = 1,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[20 * (pNumber  - mxtips)];  
 

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);    
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &lVector[20 * (rNumber - mxtips)];           
      break;
    case INNER_INNER:            
      x1 = &lVector[20 * (qNumber - mxtips)];
      x2 = &lVector[20 * (rNumber - mxtips)];     
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
   d1[0] = x1[0];
   d2[0] = x2[0];

  for(j = 1; j < 20; j++)
    {
      d1[j] = x1[j] * EXP(EIGN[j] * lz1);
      d2[j] = x2[j] * EXP(EIGN[j] * lz2);	    
    }
 
 
  for(j = 0; j < 20; j++)
    {        
      ump_x1 = 0;
      ump_x2 = 0;

      for(k = 0; k < 20; k++)
	{
	  ump_x1 += d1[k] * EI[j * 20 + k];
	  ump_x2 += d2[k] * EI[j * 20 + k];
	}
      
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 20; j++)
    x3[j] = 0.0;

  for(j = 0; j < 20; j++)          
    for(k = 0; k < 20; k++)	
      x3[k] +=  x1px2[j] *  EV[20 * j + k];	   
      
  scale = 1;
  for(k = 0; (k < 20) && scale; k++)    
    scale = ((x3[k] < minlikelihood) && (x3[k] > minusminlikelihood));    

  if(scale)
    {	        

      for(k = 0; k < 20; k++)
	x3[k]   *= twotothe256;
         
      *eVector = *eVector + 1;
    }	              

  return;
}








static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips)
{
  double lz, term;       
  double  d[20];
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc_aligned(sizeof(double) * 20 * mxtips);    

  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)    
    {
      double 
	qz = ti[k].qz[branchReference],
	rz = ti[k].rz[branchReference];
      
      qz = (qz > zmin) ? log(qz) : log(zmin);
      rz = (rz > zmin) ? log(rz) : log(zmin);

      computeVectorGTRCATPROT(lVector, &scale, ki, i, qz, rz, &ti[k], 
			  EIGN, EI, EV, 
			  tipVector, yVector, mxtips);       
    }
   
  x2 = &lVector[20 * (trav->qNumber - mxtips)]; 

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = 1.0;
  
  for(k = 1; k < 20; k++)
    d[k] =  EXP (EIGN[k] * lz);

        	   
  term =  0.0;
  for(k = 0; k < 20; k++)
    term += x1[k] * x2[k] * d[k];     

  term = LOG(FABS(term)) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);  

  return  term;
}

/******************************************/



#endif
