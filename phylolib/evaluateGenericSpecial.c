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

/* the set of functions in here computes the log likelihood at a given branch (the virtual root of a tree) */

/* includes for using SSE3 intrinsics */

#ifdef __SIM_SSE3
#include <xmmintrin.h>
#include <pmmintrin.h>
/*#include <tmmintrin.h>*/
#endif

#if (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
void perSiteLogLikelihoodsPthreads(tree *tr, double *lhs, int n, int tid); 
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


static void calcDiagptable(const double z, const int states, const int numberOfCategories, const double *rptr, const double *EIGN, double *diagptable)
{
  int 
    i, 
    l;

  double 
    lz,
    *lza = (double *)malloc(sizeof(double) * states);

  /* transform the root branch length to the log and check if it is not too small */

  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  /* do some pre-computations to avoid redundant computations further below */

  for(i = 0; i < states; i++)      
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
      diagptable[i * states + l] = EXP(rptr[i] * lza[l]);
  }

  free(lza);
}

#ifndef _OPTIMIZED_FUNCTIONS
/* below a a slow generic implementation of the likelihood computation at the root under the GAMMA model */

static double evaluateGAMMA_FLEX(int *wptr,
    double *x1_start, double *x2_start, 
    double *tipVector, 
    unsigned char *tipX1, const int n, double *diagptable, const int states)
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


  int vn = virtual_width(n);  

  printf( "n: %d %d\n", n, vn );

  if( tipX1 == 0 ) {
    reorder_back( x1_start, vn, span );
  }
  reorder_back( x2_start, vn, span );

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

      term = LOG(0.25 * FABS(term));

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

      term = LOG(0.25 * FABS(term));

      sum += wptr[i] * term;
    }                      	
  }
  if( tipX1 == 0 ) {
    reorder( x1_start, vn, span );
  }
  reorder( x2_start, vn, span );


  return sum;
} 

/* a generic and slow implementation of the CAT model of rate heterogeneity */

static double evaluateCAT_FLEX (int *cptr, int *wptr,
    double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, int n, double *diagptable_start, const int states)
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

      term = LOG(FABS(term));

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

      term = LOG(FABS(term));	 

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 
#endif

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



#ifdef _OPTIMIZED_FUNCTIONS
/* GAMMA for proteins with memory saving */

static double evaluateGTRGAMMAPROT_GAPPED_SAVE (int *wptr,
    double *x1, double *x2,  
    double *tipVector, 
    unsigned char *tipX1, int n, double *diagptable, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);


/* GAMMA for proteins */

static double evaluateGTRGAMMAPROT (int *wptr,
    double *x1, double *x2,  
    double *tipVector, 
    unsigned char *tipX1, int n, double *diagptable);

/* CAT for proteins */

static double evaluateGTRCATPROT (int *cptr, int *wptr,
    double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, int n, double *diagptable_start);


/* CAT for proteins with memory saving */

static double evaluateGTRCATPROT_SAVE (int *cptr, int *wptr,
    double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, int n, double *diagptable_start, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

/* analogous DNA fuctions */

static double evaluateGTRCAT_SAVE (int *cptr, int *wptr,
    double *x1_start, double *x2_start, double *tipVector, 		      
    unsigned char *tipX1, int n, double *diagptable_start,
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA_GAPPED_SAVE(int *wptr,
    double *x1_start, double *x2_start, 
    double *tipVector, 
    unsigned char *tipX1, const int n, double *diagptable,
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap);

static double evaluateGTRGAMMA(int *wptr,
    double *x1_start, double *x2_start, 
    double *tipVector, 
    unsigned char *tipX1, const int n, double *diagptable);


static double evaluateGTRCAT (int *cptr, int *wptr,
    double *x1_start, double *x2_start, double *tipVector, 		      
    unsigned char *tipX1, int n, double *diagptable_start);


#endif


/* This is the core function for computing the log likelihood at a branch */

void evaluateIterative(tree *tr)
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

  newviewIterative(tr, 1);  

  /* after the above call we are sure that we have properly and consistently computed the 
     conditionals to the right and left of the virtual root and we can now invoke the 
     the log likelihood computation */

  /* we need to loop over all partitions. Note that we may have a mix of DNA, protein binary data etc partitions */

  for(model = 0; model < tr->NumberOfModels; model++)
  {    
    /* whats' the number of sites of this partition (at the current thread) */
    int 	    
      width = tr->partitionData[model].width;

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
        rateHet = (int)discreteRateCategories(tr->rateHetModel),
                categories,

                /* get the number of states in the partition, e.g.: 4 = DNA, 20 = Protein */

                states = tr->partitionData[model].states;

      double 
        *rateCategories = (double*)NULL,
        z, 
        partitionLikelihood = 0.0, 	   
        *x1_start   = (double*)NULL, 
        *x2_start   = (double*)NULL,
        *diagptable = (double*)NULL,  
        *x1_gapColumn = (double*)NULL,
        *x2_gapColumn = (double*)NULL;

      unsigned int
        *x1_gap = (unsigned int*)NULL,
        *x2_gap = (unsigned int*)NULL;	 

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

      if(tr->rateHetModel == CAT)
      {	     
        rateCategories = tr->partitionData[model].perSiteRates;
        categories = tr->partitionData[model].numberOfCategories;
      }
      else
      {	     
        rateCategories = tr->partitionData[model].gammaRates;
        categories = 4;
      }

      /* set this pointer to the memory area where space has been reserved a priori for storing the 
         P matrix at the root */

      diagptable = tr->partitionData[model].left;

      /* figure out if we need to address tip vectors (a char array that indexes into a precomputed tip likelihood 
         value array or if we need to address inner vectors */

      /* either node p or node q is a tip */

      if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
      {	        	    
        /* q is a tip */

        if(isTip(qNumber, tr->mxtips))
        {	
          /* get the start address of the inner likelihood vector x2 for partition model,
             note that inner nodes are enumerated/indexed starting at 0 to save allocating some 
             space for additional pointers */

          x2_start = tr->partitionData[model].xVector[p_slot];		  

          /* get the corresponding tip vector */

          tip      = tr->partitionData[model].yVector[qNumber];	 

          /* memory saving stuff, let's deal with this later or ask Fernando ;-) */

          if(tr->saveMemory)
          {
            x2_gap         = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
            x2_gapColumn   = &(tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet]);
          }
        }           
        else
        {	
          /* p is a tip, same as above */

          x2_start = tr->partitionData[model].xVector[q_slot];		  		  
          tip = tr->partitionData[model].yVector[pNumber];

          if(tr->saveMemory)
          {
            x2_gap         = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
            x2_gapColumn   = &(tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet]);
          }

        }
      }
      else
      {  

        assert(p_slot != q_slot);
        /* neither p nor q are tips, hence we need to get the addresses of two inner vectors */

        x1_start = tr->partitionData[model].xVector[p_slot];
        x2_start = tr->partitionData[model].xVector[q_slot];

        /* memory saving option */

        if(tr->saveMemory)
        {
          x1_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
          x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
          x1_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
          x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
        }

      }


      /* if we are using a per-partition branch length estimate, the branch has an index, otherwise, for a joint branch length
         estimate over all partitions we just use the branch length value with index 0 */

      if(tr->numBranches > 1)
        z = pz[model];
      else
        z = pz[0];

      /* calc P-Matrix at root for branch z connecting nodes p and q */

      calcDiagptable(z, states, categories, rateCategories, tr->partitionData[model].EIGN, diagptable);	 

#ifndef _OPTIMIZED_FUNCTIONS

      /* generic slow functions, memory saving option is not implemented for these */

      assert(!tr->saveMemory);

      /* decide wheter CAT or GAMMA is used and compute log like */

      if(tr->rateHetModel == CAT)
        partitionLikelihood = evaluateCAT_FLEX(tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
            x1_start, x2_start, tr->partitionData[model].tipVector, 
            tip, width, diagptable, states);
      else
        partitionLikelihood = evaluateGAMMA_FLEX(tr->partitionData[model].wgt,
            x1_start, x2_start, tr->partitionData[model].tipVector,
            tip, width, diagptable, states);
#else

      /* for the optimized functions we have a dedicated, optimized function implementation 
         for each rate heterogeneity and data type combination, we switch over the number of states 
         and the rate heterogeneity model */

      switch(states)
      { 	  
        case 4: /* DNA */
          {
            if(tr->rateHetModel == CAT)
            {		  		  
              if(tr->saveMemory)
                partitionLikelihood =  evaluateGTRCAT_SAVE(tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector, 
                    tip, width, diagptable, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
              else
                partitionLikelihood =  evaluateGTRCAT(tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector, 
                    tip, width, diagptable);
            }
            else
            {		
              if(tr->saveMemory)		   
                partitionLikelihood =  evaluateGTRGAMMA_GAPPED_SAVE(tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable,
                    x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);		    
              else
                partitionLikelihood =  evaluateGTRGAMMA(tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable); 		    		  
            }
          }
          break;	  	   		   
        case 20: /* proteins */
          {
            if(tr->rateHetModel == CAT)
            {		   		  
              if(tr->saveMemory)
                partitionLikelihood = evaluateGTRCATPROT_SAVE(tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable,  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
              else
                partitionLikelihood = evaluateGTRCATPROT(tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable);		  
            }
            else
            {		    		    		      
              if(tr->saveMemory)
                partitionLikelihood = evaluateGTRGAMMAPROT_GAPPED_SAVE(tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable,
                    x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);

              else
                partitionLikelihood = evaluateGTRGAMMAPROT(tr->partitionData[model].wgt,
                    x1_start, x2_start, tr->partitionData[model].tipVector,
                    tip, width, diagptable);		         
            }
          }
          break;	      		    
        default:
          assert(0);	    
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


      partitionLikelihood += (tr->partitionData[model].globalScaler[pNumber] + tr->partitionData[model].globalScaler[qNumber]) * LOG(minlikelihood);	  

      /* now we have the correct log likelihood for the current partition after undoing scaling multiplications */	  	 

      /* finally, we also store the per partition log likelihood which is important for optimizing the alpha parameter 
         of this partition for example */

      tr->perPartitionLH[model] = partitionLikelihood; 	  
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
        tr->perPartitionLH[model] = 0.0;	   
    }
  }
}




void evaluateGeneric (tree *tr, nodeptr p, boolean fullTraversal)
{
  /* now this may be the entry point of the library to compute 
     the log like at a branch defined by p and p->back == q */

  volatile double 
    result = 0.0;

  nodeptr 
    q = p->back; 
  int 
    i,
    model;

  boolean
        p_recom = FALSE, /* if one of was missing, we will need to force recomputation */
        q_recom = FALSE;


  /* set the first entry of the traversal descriptor to contain the indices
     of nodes p and q */

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          

  /* copy the branch lengths of the tree into the first entry of the traversal descriptor.
     if -M is not used tr->numBranches must be 1 */

  for(i = 0; i < tr->numBranches; i++)    
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
    computeTraversal(tr, q, FALSE);
  }
  else
  {
    if(p_recom || needsRecomp(tr->useRecom, tr->rvec, p, tr->mxtips))
      computeTraversal(tr, p, TRUE);

    if(q_recom || needsRecomp(tr->useRecom, tr->rvec, q, tr->mxtips))
      computeTraversal(tr, q, TRUE);
  }


  /* now we copy this partition execute mask into the traversal descriptor which must come from the 
     calling program, the logic of this should not form part of the library */

  storeExecuteMaskInTraversalDescriptor(tr);  

  /* also store in the traversal descriptor that something has changed i.e., in the parallel case that the 
     traversal descriptor list of nodes needs to be broadcast once again */

  tr->td[0].traversalHasChanged = TRUE;
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

  /* now here we enter the fork-join region for Pthreads */


  /* start the parallel region and tell all threads to compute the log likelihood for 
     their fraction of the data. This call is implemented in the case switch of execFunction in axml.c
     */

  masterBarrier(THREAD_EVALUATE, tr); 

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
  /* and here is just the sequential case, we directly call evaluateIterative() above 
     without having to tell the threads/processes that they need to compute this function now */

  evaluateIterative(tr);  
#endif

  for(model = 0; model < tr->NumberOfModels; model++)
    result += tr->perPartitionLH[model];
  /* set the tree data structure likelihood value to the total likelihood */

  tr->likelihood = result;    

  if(tr->useRecom)
  {
    unpinNode(tr->rvec, p->number, tr->mxtips);
    unpinNode(tr->rvec, q->number, tr->mxtips);
  }

  /* do some bookkeeping to have traversalHasChanged in a consistent state */

  tr->td[0].traversalHasChanged = FALSE;
}


void perSiteLogLikelihoods(tree *tr, double *logLikelihoods)
{
  double 
    likelihood,
    accumulatedPerSiteLikelihood = 0.0;

  size_t
    localCount,
    i,
    globalCounter,
    model,
    lower,
    upper;

  /* compute the likelihood of the tree with the standard function to:
     1. obtain the current score for error checking
     2. store a full tree traversal in the traversal descriptor that 
     will then be used for calculating per-site log likelihoods 
     for each site individually and independently */

  evaluateGeneric (tr, tr->start, TRUE);

  likelihood = tr->likelihood;

  /* now compute per-site log likelihoods using the respective functions */

#if (defined( _USE_PTHREADS ) || defined(_FINE_GRAIN_MPI))
  /* here we need a barrier to invoke a parallel region that calls 
     function 
     perSiteLogLikelihoodsPthreads(tree *tr, double *lhs, int n, int tid)
     defined above and subsequently collects the per-site log likelihoods 
     computed by the threads and stored in local per-thread memory 
     and stores them in buffer tr->lhs.
     This corresponds to a gather operation in MPI.
     */

  masterBarrier(THREAD_PER_SITE_LIKELIHOODS, tr);

  /* 
     when the parallel region has terminated, the per-site log likelihoods 
     are stored in array tr->lhs of the master thread which we copy to the result buffer
  */
  
  memcpy(logLikelihoods, tr->lhs, sizeof(double) * tr->originalCrunchedLength);


#else

  /* sequential case: just loop over all partitions and compute per site log likelihoods */

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    lower = tr->partitionData[model].lower;
    upper = tr->partitionData[model].upper;

    for(i = lower, localCount = 0; i < upper; i++, localCount++)
    {
      double 
        l;

      /* 
         we need to switch of rate heterogeneity implementations here.
         when we have PSR we actually need to provide the per-site rate 
         to the function evaluatePartialGeneric() that computes the 
         per-site log likelihood.
         Under GAMMA, the rate will just be ignored, here we just set it to 1.0
         */

      switch(tr->rateHetModel)
      {
        case CAT:
          l = evaluatePartialGeneric (tr, i, tr->partitionData[model].perSiteRates[tr->partitionData[model].rateCategory[localCount]], model);
          break;
        case GAMMA:
          l = evaluatePartialGeneric (tr, i, 1.0, model);
          break;
        default:
          assert(0);
      }

      /* store value in result array and add the likelihood of this site to the overall likelihood */

      logLikelihoods[i] = l;
      accumulatedPerSiteLikelihood += l;
    } 
  }


  /* error checking. We need a dirt ABS() < epsilon here, because the implementations 
     (standard versus per-site) are pretty different and hence slight numerical 
     deviations are expected */

  assert(ABS(tr->likelihood - accumulatedPerSiteLikelihood) < 0.00001);
  
#endif
  


}




/* below are the optimized function versions with geeky intrinsics */

#ifdef _OPTIMIZED_FUNCTIONS

static inline boolean isGap(unsigned int *x, int pos)
{
  return (x[pos / 32] & mask32[pos % 32]);
}

static double evaluateGTRGAMMAPROT_GAPPED_SAVE (int *wptr,
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

      __m128d tv = _mm_setzero_pd();
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



      term = LOG(0.25 * FABS(term));	  

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

      __m128d tv = _mm_setzero_pd();	 	  	  

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


      term = LOG(0.25 * FABS(term));


      sum += wptr[i] * term;
    }         
  }

  return  sum;
}



static double evaluateGTRGAMMAPROT (int *wptr,
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



      term = LOG(0.25 * FABS(term));


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


      term = LOG(0.25 * FABS(term));


      sum += wptr[i] * term;
    }
  }

  return  sum;
}


static double evaluateGTRCATPROT (int *cptr, int *wptr,
    double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, int n, double *diagptable_start)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {	       	
      left = &(tipVector[20 * tipX1[i]]);
      right = &(x2[20 * i]);

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


      term = LOG(FABS(term));

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

      term = LOG(FABS(term));	 

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 


static double evaluateGTRCATPROT_SAVE (int *cptr, int *wptr,
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


      term = LOG(FABS(term));

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

      term = LOG(FABS(term));	 

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 


static double evaluateGTRCAT_SAVE (int *cptr, int *wptr,
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
      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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

      term = LOG(FABS(t[0] + t[1]));



      sum += wptr[i] * term;
    }	
  }               
  else
  {
    for (i = 0; i < n; i++) 
    { 
      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(FABS(t[0] + t[1]));

      sum += wptr[i] * term;
    }    
  }

  return  sum;         
} 


static double evaluateGTRGAMMA_GAPPED_SAVE(int *wptr,
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
      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(0.25 * FABS(t[0] + t[1]));


      sum += wptr[i] * term;
    }     
  }
  else
  {        

    for (i = 0; i < n; i++) 
    {

      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(0.25 * FABS(t[0] + t[1]));


      sum += wptr[i] * term;
    }                      	
  }

  return sum;
} 


static double evaluateGTRGAMMA(int *wptr,
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

      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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



      term = LOG(0.25 * FABS(t[0] + t[1]));



      sum += wptr[i] * term;
    }     
  }
  else
  {        
    for (i = 0; i < n; i++) 
    {

      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(0.25 * FABS(t[0] + t[1]));



      sum += wptr[i] * term;
    }                      	
  }

  return sum;
} 


static double evaluateGTRCAT (int *cptr, int *wptr,
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

      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(FABS(t[0] + t[1]));


      sum += wptr[i] * term;
    }	
  }               
  else
  {
    for (i = 0; i < n; i++) 
    { 

      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


      term = LOG(FABS(t[0] + t[1]));


      sum += wptr[i] * term;
    }    
  }

  return  sum;         
} 





#endif
