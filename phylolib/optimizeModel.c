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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
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

static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[20];




/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


/* the following function is used to set rates in the Q matrix 
   the data structure called symmetryVector is used to 
   define the symmetries between rates as they are specified 
   in some of the secondary structure substitution models that 
   generally don't use GTR matrices but more restricted forms thereof */

static void setRateModel(tree *tr, int model, double rate, int position)
{
  int
    states   = tr->partitionData[model].states,
    numRates = (states * states - states) / 2;

  if(tr->partitionData[model].dataType == DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(tr->partitionData[model].dataType != BINARY_DATA); 

  assert(rate >= RATE_MIN && rate <= RATE_MAX);

  if(tr->partitionData[model].nonGTR)
    {    
      int 
	i, 
	k = tr->partitionData[model].symmetryVector[position];

      assert(tr->partitionData[model].dataType == SECONDARY_DATA ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_7);

      if(k == -1)
	tr->partitionData[model].substRates[position] = 0.0;
      else
	{
	  if(k == tr->partitionData[model].symmetryVector[numRates - 1])
	    {
	      for(i = 0; i < numRates - 1; i++)
		if(tr->partitionData[model].symmetryVector[i] == k)
		  tr->partitionData[model].substRates[position] = 1.0;
	    }
	  else
	    {
	      for(i = 0; i < numRates - 1; i++)
		{
		  if(tr->partitionData[model].symmetryVector[i] == k)
		    tr->partitionData[model].substRates[i] = rate; 
		}	      	     
	    }
	}
    }
  else
    tr->partitionData[model].substRates[position] = rate;
}


/* 
   the following three functions are used to link/unlink parameters 
   between partitions. This should work in a generic way, however 
   this is so far mainly used for linking unlinking GTR matrix parameter 
   estimates across different protein data partitions.
   Generally this mechanism can also be used for linking/inlinking alpha paremeters 
   between partitions and the like.
   However, all alpha parameter estimates for all partitions and GTR estimates for 
   DNA partitions are unlinked by default. This is actually hard-coded 
   in here. 
*/

/* initializwe a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
*/

static linkageList* initLinkageList(int *linkList, tree *tr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  
  linkageList 
    *ll = (linkageList*)malloc(sizeof(linkageList));
    
  /* figure out how many distinct parameters we need to estimate 
     in total, if all parameters are linked the result will be 1 if all 
     are unlinked the result will be tr->NumberOfModels */
  
  for(i = 0; i < tr->NumberOfModels; i++)
    {
      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  /* allocate the linkage list data structure that containes information which parameters of which partition are 
     linked with each other.

     Note that we need a separate invocation of initLinkageList() and a separate linkage list 
     for each parameter type */

  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)malloc(sizeof(linkageData) * numberOfModels);

  /* noe loop over the number of free parameters and assign the corresponding partitions to each parameter */

  for(i = 0; i < numberOfModels; i++)
    {
      /* 
	 the valid flag is used for distinguishing between DNA and protein data partitions.
	 This can be used to enable/disable parameter optimization for the paremeter 
	 associated to the corresponding partitions. This deature is used in optRatesGeneric 
	 to first optimize all DNA GTR rate matrices and then all PROT GTR rate matrices */

      ll->ld[i].valid = TRUE;
      partitions = 0;

      /* now figure out how many partitions share this joint parameter */

      for(k = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  partitions++;	    

      /* assign a list to store the partitions that share the parameter */

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)malloc(sizeof(int) * partitions);
      
      /* now store the respective partition indices in this list */
      
      for(k = 0, pos = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  /* return the linkage list for the parameter */

  return ll;
}

/* dedicated helper function to initialize the linkage list, that is, essentiaylly compute 
   the integer vector int *linkList used above for linking GTR models.
   
   Once again, this is hard-coded in RAxML, because users can not influence the linking.

*/
   

static linkageList* initLinkageListGTR(tree *tr)
{
  int
    i,
    *links = (int*)malloc(sizeof(int) * tr->NumberOfModels),
    firstAA = tr->NumberOfModels + 2,
    countGTR = 0,
    countOtherModel = 0;
  
  linkageList
    * ll;

  /* here we only want to figure out if either all prot data partitions 
     are supposed to use a joint GTR prot subst matrix or not 

    We either allow ALL prot partitions to use a shared/joint estimate of the GTR matrix or not,
    things like having one prot partition evolving under WAG and the others under a joint GTR estimate are 
    not allowed.
  */

  for(i = 0; i < tr->NumberOfModels; i++)
    {     
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  if(tr->partitionData[i].protModels == GTR)
	    {
	      if(i < firstAA)
		firstAA = i;
	      countGTR++;
	    }
	  else
	    countOtherModel++;
	}
    }
  
  assert((countGTR > 0 && countOtherModel == 0) || (countGTR == 0 && countOtherModel > 0) ||  (countGTR == 0 && countOtherModel == 0));

  /* if there is no joint GTR matrix optimization for protein data partitions we can unlink rate matrix calculations for all partitions */

  if(countGTR == 0)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	links[i] = i;
    }
  else
    {
      /* otherwise we let all partitions, except for the protein partitions use 
	 unlinked rate matrices while we link the GTR rate matrices of all 
	 protein data partitions */
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  switch(tr->partitionData[i].dataType)
	    {	   
	    case DNA_DATA:
	    case BINARY_DATA:
	    case GENERIC_32:
	    case GENERIC_64:
	    case SECONDARY_DATA:
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7: 
	      links[i] = i;
	      break;
	    case AA_DATA:	  
	      links[i] = firstAA;
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  

  /* we can now pass an appropriate integer vector to the linkage list initialization function :-) */

  ll = initLinkageList(links, tr); 

  free(links);
  
  return ll;
}

/* free linkage list data structure */

static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    free(ll->ld[i].partitionList);         

  free(ll->ld);
  free(ll);   
}

#define ALPHA_F 0
#define INVAR_F 1
#define RATE_F  2


/* function that evaluates the change to a parameter */

static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll)
{ 
  int i, k, pos;

  switch(whichFunction)
    {    
    case RATE_F:
      /* loop over linkage list entries for Q matrix rates */

      for(i = 0, pos = 0; i < ll->entries; i++)
	{
	  /* if valid, i.e. we want to assess the change for a certain partition data type, e.g., DNA */
	  if(ll->ld[i].valid)
	    {
	      /* if the iterative numerical procedure for the partitions sharing the parameter 
		 has already converged don't do anything.
		 This stuff is required for parallel load balanicing as described in:
		 
		 A. Stamatakis, M. Ott: "Load Balance in the Phylogenetic Likelihood Kernel". Proceedings of ICPP 2009, Vienna, Austria, September 2009.
	      */
	      if(converged[pos])
		{		 
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
		}
	      else
		{
		  /* loop over partitions that share the same parameter */

		  for(k = 0; k < ll->ld[i].partitions; k++)
		    {
		      int index = ll->ld[i].partitionList[k];		  
		      
		      /* update them to new, proposed value for which we want to obtain the likelihood */
	      
		      setRateModel(tr, index, value[pos], rateNumber);  

		      /* in the case of rates, i.e., Q matrix, re-compute the corresponding eigenvector eigenvalue decomposition */

		      initReversibleGTR(tr, index);		 
		    }
		}
	      pos++;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;	     
	    }
	 
	}

      assert(pos == numberOfModels);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))      
      masterBarrier(THREAD_OPT_RATE, tr);
#else
      /* and compute the likelihood by doing a full tree traversal :-) */
      evaluateGeneric(tr, tr->start, TRUE);      
#endif     
      
      /* update likelihoods and the sum of per-partition likelihoods for those partitions that share the parameter.
	 that's the likelihood we actually want to optimize here */

      for(i = 0, pos = 0; i < ll->entries; i++)	
	{
	  if(ll->ld[i].valid)
	    {
	      result[pos] = 0.0;
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];

		  assert(tr->perPartitionLH[index] <= 0.0);

		  result[pos] -= tr->perPartitionLH[index];
		  
		}
	      pos++;
	    }
	   for(k = 0; k < ll->ld[i].partitions; k++)
	     {
	       int index = ll->ld[i].partitionList[k];
	       tr->executeModel[index] = TRUE;
	     }	  
	}

      assert(pos == numberOfModels);
      break;
    case ALPHA_F:

      /* analoguos to the rate stuff above */

      for(i = 0; i < ll->entries; i++)
	{
	  if(converged[i])
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
	    }
	  else
	    {
	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int index = ll->ld[i].partitionList[k];
		  tr->executeModel[index] = TRUE;
		  tr->partitionData[index].alpha = value[i];

		  /* re-compute the discrete gamma function approximation for the new alpha parameter */
		  makeGammaCats(tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useMedian);
		}
	    }
	}
#if (defined( _USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
      masterBarrier(THREAD_OPT_ALPHA, tr);
#else  
      evaluateGeneric(tr, tr->start, TRUE);
#endif
            
      for(i = 0; i < ll->entries; i++)	
	{	  
	  result[i] = 0.0;
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {
	      int index = ll->ld[i].partitionList[k];
	      	      
	      assert(tr->perPartitionLH[index] <= 0.0);		
	      
	      result[i] -= tr->perPartitionLH[index];	            
	      tr->executeModel[index] = TRUE;
	    }
	}
      break;
    default:
      assert(0);	
    }

}

/* generic implementation of Brent's algorithm for one-dimensional parameter optimization */

static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double lim_inf, double lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)malloc(sizeof(double) * numberOfModels),
    *b     = (double *)malloc(sizeof(double) * numberOfModels),
    *d     = (double *)malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)malloc(sizeof(double) * numberOfModels),
    *p     = (double *)malloc(sizeof(double) * numberOfModels),
    *q     = (double *)malloc(sizeof(double) * numberOfModels),
    *r     = (double *)malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)malloc(sizeof(double) * numberOfModels),
    *u     = (double *)malloc(sizeof(double) * numberOfModels),
    *v     = (double *)malloc(sizeof(double) * numberOfModels),
    *w     = (double *)malloc(sizeof(double) * numberOfModels),
    *x     = (double *)malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)malloc(sizeof(double) * numberOfModels),
    *e     = (double *)malloc(sizeof(double) * numberOfModels);
  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf && a[i] <= lim_sup);
      assert(b[i] >= lim_inf && b[i] <= lim_sup);
      assert(x[i] >= lim_inf && x[i] <= lim_sup);
      assert(v[i] >= lim_inf && v[i] <= lim_sup);
      assert(w[i] >= lim_inf && w[i] <= lim_sup);
    }
  
  

  for(iter = 1; iter <= ITMAX; iter++)
    {
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  free(converged);
	  free(a);
	  free(b);
	  free(d);
	  free(etemp);
	  free(fu);
	  free(fv);
	  free(fw);
	  free(fx);
	  free(p);
	  free(q);
	  free(r);
	  free(tol1);
	  free(tol2);
	  free(u);
	  free(v);
	  free(w);
	  free(x);
	  free(xm);
	  free(e);
	  return;
	}     

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	     	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
  
	      xm[i] = 0.5 * (a[i] + b[i]);
	      tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
	  
	      if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
		{		 
		  result[i] =  -fx[i];
		  xmin[i]   = x[i];
		  converged[i] = TRUE;		  
		}
	      else
		{
		  if(fabs(e[i]) > tol1[i])
		    {		     
		      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
		      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
		      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
		      q[i] = 2.0 * (q[i] - r[i]);
		      if(q[i] > 0.0)
			p[i] = -p[i];
		      q[i] = fabs(q[i]);
		      etemp[i] = e[i];
		      e[i] = d[i];
		      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
			d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
		      else
			{
			  d[i] = p[i] / q[i];
			  u[i] = x[i] + d[i];
			  if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
			    d[i] = SIGN(tol1[i], xm[i] - x[i]);
			}
		    }
		  else
		    {		     
		      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
		    }
		  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]): (x[i] +SIGN(tol1[i], d[i])));
		}

	      if(!converged[i])
		assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
                 
      evaluateChange(tr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      if(fu[i] <= fx[i])
		{
		  if(u[i] >= x[i])
		    a[i] = x[i];
		  else
		    b[i] = x[i];
		  
		  SHFT(v[i],w[i],x[i],u[i]);
		  SHFT(fv[i],fw[i],fx[i],fu[i]);
		}
	      else
		{
		  if(u[i] < x[i])
		    a[i] = u[i];
		  else
		    b[i] = u[i];
		  
		  if(fu[i] <= fw[i] || w[i] == x[i])
		    {
		      v[i] = w[i];
		      w[i] = u[i];
		      fv[i] = fw[i];
		      fw[i] = fu[i];
		    }
		  else
		    {
		      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
			{
			  v[i] = u[i];
			  fv[i] = fu[i];
			}
		    }	    
		}
	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
	      assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
    }

  free(converged);
  free(a);
  free(b);
  free(d);
  free(etemp);
  free(fu);
  free(fv);
  free(fw);
  free(fx);
  free(p);
  free(q);
  free(r);
  free(tol1);
  free(tol2);
  free(u);
  free(v);
  free(w);
  free(x);
  free(xm);
  free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}

/* generic bracketing function required for Brent's algorithm. For details please see the corresponding chapter in the book Numerical Recipees in C */

static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double lim_inf, double lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll)
{
  double 
    *ulim = (double *)malloc(sizeof(double) * numberOfModels),
    *u    = (double *)malloc(sizeof(double) * numberOfModels),
    *r    = (double *)malloc(sizeof(double) * numberOfModels),
    *q    = (double *)malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)malloc(sizeof(int) * numberOfModels),
    *endState = (int *)malloc(sizeof(int) * numberOfModels);

  boolean *converged = (boolean *)malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup) 	
	param[i] = ax[i] = lim_sup;
      
      if(param[i] < lim_inf) 
	param[i] = ax[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
   
  
  evaluateChange(tr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup) 
	param[i] = bx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = bx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
  evaluateChange(tr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup) 
	param[i] = cx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = cx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
 
  evaluateChange(tr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll);

   while(1) 
     {       
       allConverged = TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
	 allConverged = allConverged && converged[i];

       if(allConverged)
	 {
	   for(i = 0; i < numberOfModels; i++)
	     {	       
	       if(ax[i] > lim_sup) 
		 ax[i] = lim_sup;
	       if(ax[i] < lim_inf) 
		 ax[i] = lim_inf;

	       if(bx[i] > lim_sup) 
		 bx[i] = lim_sup;
	       if(bx[i] < lim_inf) 
		 bx[i] = lim_inf;
	       
	       if(cx[i] > lim_sup) 
		 cx[i] = lim_sup;
	       if(cx[i] < lim_inf) 
		 cx[i] = lim_inf;
	     }

	   free(converged);
	   free(ulim);
	   free(u);
	   free(r);
	   free(q);
	   free(fu);
	   free(dum); 
	   free(temp);
	   free(state);   
	   free(endState);
	   return 0;
	   
	 }

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {
	       switch(state[i])
		 {
		 case 0:
		   endState[i] = 0;
		   if(!(fb[i] > fc[i]))		         
		     converged[i] = TRUE;		       		     
		   else
		     {
		   
		       if(ax[i] > lim_sup) 
			 ax[i] = lim_sup;
		       if(ax[i] < lim_inf) 
			 ax[i] = lim_inf;
		       if(bx[i] > lim_sup) 
			 bx[i] = lim_sup;
		       if(bx[i] < lim_inf) 
			 bx[i] = lim_inf;
		       if(cx[i] > lim_sup) 
			 cx[i] = lim_sup;
		       if(cx[i] < lim_inf) 
			 cx[i] = lim_inf;
		       
		       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			 (2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		       if(u[i] > lim_sup) 
			 u[i] = lim_sup;
		       if(u[i] < lim_inf) 
			 u[i] = lim_inf;
		       if(ulim[i] > lim_sup) 
			 ulim[i] = lim_sup;
		       if(ulim[i] < lim_inf) 
			 ulim[i] = lim_inf;
		       
		       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			 {
			   param[i] = u[i];
			   if(param[i] > lim_sup) 			     
			     param[i] = u[i] = lim_sup;
			   if(param[i] < lim_inf)
			     param[i] = u[i] = lim_inf;
			   endState[i] = 1;
			 }
		       else 
			 {
			   if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			     {
			       param[i] = u[i];
			       if(param[i] > lim_sup) 
				 param[i] = u[i] = lim_sup;
			       if(param[i] < lim_inf) 
				 param[i] = u[i] = lim_inf;
			       endState[i] = 2;
			     }		  	       
			   else
			     {
			       if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				 {
				   u[i] = ulim[i];
				   param[i] = u[i];	
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = ulim[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = ulim[i] = lim_inf;
				   endState[i] = 0;
				 }		  		
			       else 
				 {		  
				   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				   param[i] = u[i];
				   endState[i] = 0;
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = lim_inf;
				 }
			     }	  
			 }
		     }
		   break;
		 case 1:
		   endState[i] = 0;
		   break;
		 case 2:
		   endState[i] = 3;
		   break;
		 default:
		   assert(0);
		 }
	       assert(param[i] >= lim_inf && param[i] <= lim_sup);
	     }
	 }
             
       evaluateChange(tr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll);

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {	       
	       switch(endState[i])
		 {
		 case 0:
		   fu[i] = temp[i];
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 case 1:
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {
		       ax[i]=(bx[i]);
		       bx[i]=u[i];
		       fa[i]=(fb[i]);
		       fb[i]=fu[i]; 
		       converged[i] = TRUE;		      
		     } 
		   else 
		     {
		       if (fu[i] > fb[i]) 
			 {
			   assert(u[i] >= lim_inf && u[i] <= lim_sup);
			   cx[i]=u[i];
			   fc[i]=fu[i];
			   converged[i] = TRUE;			  
			 }
		       else
			 {		   
			   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			   param[i] = u[i];
			   if(param[i] > lim_sup) {param[i] = u[i] = lim_sup;}
			   if(param[i] < lim_inf) {param[i] = u[i] = lim_inf;}	  
			   state[i] = 1;		 
			 }		  
		     }
		   break;
		 case 2: 
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {		     
		       SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
		       state[i] = 2;
		     }	   
		   else
		     {
		       state[i] = 0;
		       SHFT(ax[i],bx[i],cx[i],u[i]);
		       SHFT(fa[i],fb[i],fc[i],fu[i]);
		     }
		   break;	   
		 case 3:		  
		   SHFT(fb[i],fc[i],fu[i], temp[i]);
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 default:
		   assert(0);
		 }
	     }
	 }
    }
   

   assert(0);
   free(converged);
   free(ulim);
   free(u);
   free(r);
   free(q);
   free(fu);
   free(dum); 
   free(temp);
   free(state);   
   free(endState);

  

   return(0);
}


/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/




/* function for optimizing alpha parameter */


static void optAlpha(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i, 
    k,
    numberOfModels = ll->entries;
  
  double 
    lim_inf     = ALPHA_MIN,
    lim_sup     = ALPHA_MAX;
  double
    *startLH    = (double *)malloc(sizeof(double) * numberOfModels),
    *startAlpha = (double *)malloc(sizeof(double) * numberOfModels),
    *endAlpha   = (double *)malloc(sizeof(double) * numberOfModels),
    *_a         = (double *)malloc(sizeof(double) * numberOfModels),
    *_b         = (double *)malloc(sizeof(double) * numberOfModels),
    *_c         = (double *)malloc(sizeof(double) * numberOfModels),
    *_fa        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fb        = (double *)malloc(sizeof(double) * numberOfModels),
    *_fc        = (double *)malloc(sizeof(double) * numberOfModels),
    *_param     = (double *)malloc(sizeof(double) * numberOfModels),
    *result     = (double *)malloc(sizeof(double) * numberOfModels),
    *_x         = (double *)malloc(sizeof(double) * numberOfModels);   

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
   int revertModel = 0;
#endif   

   evaluateGeneric(tr, tr->start, TRUE);
   
   /* 
     at this point here every worker has the traversal data it needs for the 
     search, so we won't re-distribute it he he :-)
  */

  for(i = 0; i < numberOfModels; i++)
    {
      /* make sure that we want to optimize alpha here */

      assert(ll->ld[i].valid);
      
      /* get the representative alpha value for all partitions that share it.
	 In RAxML we have hard-coded all alpha parameters to be estimated separately.
	 hence linking alpha parameters has not been tested for a long long time, be cautious 
      */
      
      startAlpha[i] = tr->partitionData[ll->ld[i].partitionList[0]].alpha;
      _a[i] = startAlpha[i] + 0.1;
      _b[i] = startAlpha[i] - 0.1;      
      if(_b[i] < lim_inf) 
	_b[i] = lim_inf;

      startLH[i] = 0.0;
      
      for(k = 0; k < ll->ld[i].partitions; k++)	
	{
	  /* sum over all per-partition log likelihood scores that share the same alpha parameter */
	  startLH[i] += tr->perPartitionLH[ll->ld[i].partitionList[k]];
	  
	  /* make sure that all copies of the linked alpha parameter have been updated appropriately */
	  assert(tr->partitionData[ll->ld[i].partitionList[0]].alpha ==  tr->partitionData[ll->ld[i].partitionList[k]].alpha);
	}
    }					  

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, -1, ALPHA_F, tr, ll);       
  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, ALPHA_F, -1, tr, ll, lim_inf, lim_sup);

  for(i = 0; i < numberOfModels; i++)
    endAlpha[i] = result[i];
  
  /* if for some reason we couldn't improve the likelihood further, restore the old alpha value for all partitions sharing it */

  for(i = 0; i < numberOfModels; i++)
    {
      if(startLH[i] > endAlpha[i])
	{    	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {	      
	      tr->partitionData[ll->ld[i].partitionList[k]].alpha = startAlpha[i];

	      makeGammaCats(tr->partitionData[ll->ld[i].partitionList[k]].alpha, tr->partitionData[ll->ld[i].partitionList[k]].gammaRates, 4, tr->useMedian); 
	    }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  revertModel++;
#endif
	}  
    }

  /* broadcast new alpha value to all parallel threads/processes */

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  if(revertModel > 0)
    masterBarrier(THREAD_COPY_ALPHA, tr);
#endif
  
  free(startLH);
  free(startAlpha);
  free(endAlpha);
  free(result);
  free(_a);
  free(_b);
  free(_c);
  free(_fa);
  free(_fb);
  free(_fc);
  free(_param);
  free(_x);  

}

/* optimize rates in the Q matrix */

static void optRates(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int 
    i, 
    k, 
    j, 
    pos,
    numberOfRates = ((states * states - states) / 2) - 1;
    
  double lim_inf = RATE_MIN;
  double lim_sup = RATE_MAX;
  double 
    *startRates;
  double 
    *startLH= (double *)malloc(sizeof(double) * numberOfModels),
    *endLH  = (double *)malloc(sizeof(double) * numberOfModels),
    *_a     = (double *)malloc(sizeof(double) * numberOfModels),
    *_b     = (double *)malloc(sizeof(double) * numberOfModels),
    *_c     = (double *)malloc(sizeof(double) * numberOfModels),
    *_fa    = (double *)malloc(sizeof(double) * numberOfModels),
    *_fb    = (double *)malloc(sizeof(double) * numberOfModels),
    *_fc    = (double *)malloc(sizeof(double) * numberOfModels),
    *_param = (double *)malloc(sizeof(double) * numberOfModels),
    *result = (double *)malloc(sizeof(double) * numberOfModels),
    *_x     = (double *)malloc(sizeof(double) * numberOfModels); 

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
   int revertModel = 0;
#endif

   assert(states != -1);

  startRates = (double *)malloc(sizeof(double) * numberOfRates * numberOfModels);

  evaluateGeneric(tr, tr->start, TRUE);
  
  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  /* get the initial rates */

  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  /* loop over all partitions that are linked via the corresponding GTR matrx, i.e., that share the same GTR matrix */

	  for(j = 0; j < ll->ld[i].partitions; j++)
	    {
	      int 
		index = ll->ld[i].partitionList[j];
	      
	      startLH[pos] += tr->perPartitionLH[index];
	      
	      for(k = 0; k < numberOfRates; k++)
		startRates[pos * numberOfRates + k] = tr->partitionData[index].substRates[k];      
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
  
  /* 
     Now loop over all rates in the matrix, e.g., 5 if it's DNA, and 189 if it's protein data.
     Note that we do this on a rate by rate basis, i.e., first try to optimize a->c, a->g, etc.
   */

  for(i = 0; i < numberOfRates; i++)
    {     
      for(k = 0, pos = 0; k < ll->entries; k++)
	{
	  if(ll->ld[k].valid)
	    {
	      int index = ll->ld[k].partitionList[0];
	      _a[pos] = tr->partitionData[index].substRates[i] + 0.1;
	      _b[pos] = tr->partitionData[index].substRates[i] - 0.1;
	      
	      if(_a[pos] < lim_inf) _a[pos] = lim_inf;
	      if(_a[pos] > lim_sup) _a[pos] = lim_sup;
	      
	      if(_b[pos] < lim_inf) _b[pos] = lim_inf;
	      if(_b[pos] > lim_sup) _b[pos] = lim_sup;    
	      pos++;
	    }
	}                       	     

      assert(pos == numberOfModels);

      brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, i, RATE_F, tr, ll);
      
      for(k = 0; k < numberOfModels; k++)
	{
	  assert(_a[k] >= lim_inf && _a[k] <= lim_sup);
	  assert(_b[k] >= lim_inf && _b[k] <= lim_sup);	  
	  assert(_c[k] >= lim_inf && _c[k] <= lim_sup);	    
	}      

      brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, RATE_F, i, tr,  ll, lim_inf, lim_sup);
	
      for(k = 0; k < numberOfModels; k++)
	endLH[k] = result[k];	           
      
      /* if the proposed new rate does not improve the likelihood undo the change */

      for(k = 0, pos = 0; k < ll->entries; k++)
	{
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  revertModel = 0;
#endif
	  if(ll->ld[k].valid)
	    { 
	      if(startLH[pos] > endLH[pos])
		{		  
		  for(j = 0; j < ll->ld[k].partitions; j++)
		    {
		      int index = ll->ld[k].partitionList[j];
		      tr->partitionData[index].substRates[i] = startRates[pos * numberOfRates + i];

		      initReversibleGTR(tr, index);
		    }
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
		  revertModel++;
#endif
		}
	      pos++;
	    }
	}

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if(revertModel > 0)
	masterBarrier(THREAD_COPY_RATES, tr);
#endif
      
      assert(pos == numberOfModels);
    }

 
  free(startLH);
  free(endLH);
  free(result);
  free(_a);
  free(_b);
  free(_c);
  free(_fa);
  free(_fb);
  free(_fc);
  free(_param);
  free(_x);  
  free(startRates);
}


/* figure out if all AA models have been assigned a joint GTR matrix */

static boolean AAisGTR(tree *tr)
{
  int i, count = 0;

  for(i = 0; i < tr->NumberOfModels; i++)   
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  count++;
	  if(tr->partitionData[i].protModels != GTR)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}


/* generic substitiution matrix (Q matrix) optimization */

static void optRatesGeneric(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitions  = 0,
    secondaryPartitions = 0,
    secondaryModel = -1,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* 
     first optimize all rates in DNA data partition matrices. That's where we use the valid field in the 
     linkage list data structure. 
   */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  ll->ld[i].valid = TRUE;
	  dnaPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  /* if we have dna partitions in our dataset, let's optimize all 5 rates in their substitution matrices */

  if(dnaPartitions > 0)
    optRates(tr, modelEpsilon, ll, dnaPartitions, states);
  

  /* now we do the same thing for secondary structure models  */

   for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	  /* 
	     we only have one type of secondary data models in one analysis, RAxML allows for only one secondary data 
	     partition !
	   */
	case SECONDARY_DATA_6:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_6;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA_7: 
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_7;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:	
	case DNA_DATA:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }

   /* if there are secondary struct data partitions (and there can only be at most one), let's optimize 
      their rates, tehre are models with 6, 7 and 16 states */
   
   if(secondaryPartitions > 0)
     {
       assert(secondaryPartitions == 1);

       switch(secondaryModel)
	 {
	 case SECONDARY_DATA:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_6:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_7:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break; 
	 default:
	   assert(0);
	 }
     }

  /* then AA for GTR */

   /* now if all AA partitions share a joint GTR subst matrix, let's do a joint estimate 
      of the 189 rates across all of them. Otherwise we don't need to optimize anything since 
      we will be using one of the fixed models like WAG, JTT, etc */

  if(AAisGTR(tr))
    {
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case AA_DATA:
	      states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitions == 1);     
      
      optRates(tr, modelEpsilon, ll, aaPartitions, states);
    }
  
  /* then multi-state data partitions

     now here we have to be careful, because every multi-state partition can actually 
     have a distinct number of states, so we will go to every multi-state partition separately,
     parallel efficiency for this will suck, but what the hell .....
  */

  if(tr->multiStateModel == GTR_MULTI_STATE)
    {     
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case GENERIC_32:
	      {
		int k;
		
		states = tr->partitionData[ll->ld[i].partitionList[0]].states;			      

		ll->ld[i].valid = TRUE;
		
		for(k = 0; k < ll->entries; k++)
		  if(k != i)
		    ll->ld[k].valid = FALSE;
		
		optRates(tr, modelEpsilon, ll, 1, states);
	      }
	      break;
	    case AA_DATA:	    
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	    case GENERIC_64:
	      break;
	    default:
	      assert(0);
	    }	 
	}           
    }

  /* done with all partitions, so we can set all entries in the linkage list to valid again :-) */

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}





/*********************FUNCTIONS FOR PSR/CAT model of rate heterogeneity ***************************************/






static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorizePartition(tree *tr, rateCategorize *rc, int model, int lower, int upper)
{
  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
      {
	double
	  temp = tr->patrat[i];

	int
	  found = 0;
	
	for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;
		tr->rateCategory[i] = k; 
		break;
	      }
	  }
	
	if(!found)
	  {
	    min = fabs(temp - rc[0].rate);
	    tr->rateCategory[i] = 0;

	    for(k = 1; k < tr->partitionData[model].numberOfCategories; k++)
	      {
		diff = fabs(temp - rc[k].rate);

		if(diff < min)
		  {
		    min = diff;
		    tr->rateCategory[i] = k;
		  }
	      }
	  }
      }

  for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
    tr->partitionData[model].perSiteRates[k] = rc[k].rate; 
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model, 
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      int 
	localIndex = 0;

      boolean 
	execute = ((tr->manyPartitions && isThisMyPartition(tr, tid, model)) || (!tr->manyPartitions));

      if(execute)
	for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	  {
	    if(tr->manyPartitions || (i % n == tid))
	      {
	      
		double initialRate, initialLikelihood, 
		  leftLH, rightLH, leftRate, rightRate, v;
		const double epsilon = 0.00001;
		int k;	      
		
		tr->patrat[i] = tr->patratStored[i];     
		initialRate = tr->patrat[i];
		
		initialLikelihood = evaluatePartialGeneric(tr, localIndex, initialRate, model); /* i is real i ??? */
		
		
		leftLH = rightLH = initialLikelihood;
		leftRate = rightRate = initialRate;
		
		k = 1;
		
		while((initialRate - k * lower_spacing > 0.0001) && 
		      ((v = evaluatePartialGeneric(tr, localIndex, initialRate - k * lower_spacing, model)) 
		       > leftLH) && 
		      (fabs(leftLH - v) > epsilon))  
		  {	  
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif
		    
		    leftLH = v;
		    leftRate = initialRate - k * lower_spacing;
		    k++;	  
		  }      
		
		k = 1;
		
		while(((v = evaluatePartialGeneric(tr, localIndex, initialRate + k * upper_spacing, model)) > rightLH) &&
		      (fabs(rightLH - v) > epsilon))    	
		  {
#ifndef WIN32
		    if(isnan(v))
		      assert(0);
#endif     
		    rightLH = v;
		    rightRate = initialRate + k * upper_spacing;	 
		    k++;
		  }           
		
		if(rightLH > initialLikelihood || leftLH > initialLikelihood)
		  {
		    if(rightLH > leftLH)	    
		      {	     
			tr->patrat[i] = rightRate;
			lhs[i] = rightLH;
		      }
		    else
		      {	      
			tr->patrat[i] = leftRate;
			lhs[i] = leftLH;
		      }
		  }
		else
		  lhs[i] = initialLikelihood;
		
		tr->patratStored[i] = tr->patrat[i];
		localIndex++;
	      }
	  }
      assert(localIndex == tr->partitionData[model].width);    
    }
}



#else


static void optRateCatModel(tree *tr, int model, double lower_spacing, double upper_spacing, double *lhs)
{
  int lower = tr->partitionData[model].lower;
  int upper = tr->partitionData[model].upper;
  int i;
  for(i = lower; i < upper; i++)
    {
      double initialRate, initialLikelihood, 
	leftLH, rightLH, leftRate, rightRate, v;
      const double epsilon = 0.00001;
      int k;
      
      tr->patrat[i] = tr->patratStored[i];     
      initialRate = tr->patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(tr, i, initialRate, model); 
      
      
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
	    ((v = evaluatePartialGeneric(tr, i, initialRate - k * lower_spacing, model)) 
	     > leftLH) && 
	    (fabs(leftLH - v) > epsilon))  
	{	  
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif
	  
	  leftLH = v;
	  leftRate = initialRate - k * lower_spacing;
	  k++;	  
	}      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(tr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
	    (fabs(rightLH - v) > epsilon))    	
	{
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif     
	  rightLH = v;
	  rightRate = initialRate + k * upper_spacing;	 
	  k++;
	}           
  
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	{
	  if(rightLH > leftLH)	    
	    {	     
	      tr->patrat[i] = rightRate;
	      lhs[i] = rightLH;
	    }
	  else
	    {	      
	      tr->patrat[i] = leftRate;
	      lhs[i] = leftLH;
	    }
	}
      else
	lhs[i] = initialLikelihood;
      
      tr->patratStored[i] = tr->patrat[i];
    }

}


#endif



/* 
   set scaleRates to FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

void updatePerSiteRates(tree *tr, boolean scaleRates)
{
  int 
    i,
    model;

  if(tr->numBranches > 1)
    {            
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 	       
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;
	  
	  if(scaleRates)
	    {
	      double 
		scaler = 0.0,       
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	  
	      scaler = 1.0 / ((double)accRat);
	  	  
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] *= scaler;	    

	      accRat = 0.0;	 
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}	         

	      accRat /= ((double)accWgt);	  

	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }
	  else
	    {
	      double 		   
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	      
	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }

	  
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  {
	    int 
	      localCount = 0;
	    
	    for(i = lower, localCount = 0; i < upper; i++, localCount++)
	      {	    	      
		tr->partitionData[model].rateCategory[localCount] = tr->rateCategory[i];
	      }
	  }
#endif
	}
    }
  else
    {
      int
	accWgt = 0;

      double 
	scaler = 0.0,       
	accRat = 0.0; 

      if(scaleRates)
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}
	    }
	  
	  accRat /= ((double)accWgt);
	  
	  scaler = 1.0 / ((double)accRat);
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] *= scaler;
	    }

	  for(model = 0, accRat = 0.0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}
	    }           

	  accRat /= ((double)accWgt);	  

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
      else
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}
	    }
	  
	  accRat /=  (double)accWgt;

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
         
       for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 
	    localCount = 0,
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;

	}         
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      for(model = 0; model < tr->NumberOfModels; model++)	
	{   	  	  	 
	  int 
	    localCount,
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;	  
	  
	  for(i = lower, localCount = 0; i < upper; i++, localCount++)
	      tr->partitionData[model].rateCategory[localCount] = tr->rateCategory[i];
	}
#endif
    }
  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif               
}

static void optimizeRateCategories(tree *tr, int _maxCategories)
{
  assert(_maxCategories > 0);

  if(_maxCategories > 1)
    {
      double  
	temp,  
	lower_spacing, 
	upper_spacing,
	initialLH = tr->likelihood,	
	*ratStored = (double *)malloc(sizeof(double) * tr->originalCrunchedLength),
	/**lhs =       (double *)malloc(sizeof(double) * tr->originalCrunchedLength),*/
	**oldCategorizedRates = (double **)malloc(sizeof(double *) * tr->NumberOfModels);

      int  
	i,
	k,
	maxCategories = _maxCategories,
	*oldCategory =  (int *)malloc(sizeof(int) * tr->originalCrunchedLength),
	model,
	*oldNumbers = (int *)malloc(sizeof(int) * tr->NumberOfModels);
  
      assert(isTip(tr->start->number, tr->mxtips));         
      
      evaluateGeneric(tr, tr->start, TRUE);

      if(tr->optimizeRateCategoryInvocations == 1)
	{
	  lower_spacing = 0.5 / ((double)(tr->optimizeRateCategoryInvocations));
	  upper_spacing = 1.0 / ((double)(tr->optimizeRateCategoryInvocations));
	}
      else
	{
	  lower_spacing = 0.05 / ((double)(tr->optimizeRateCategoryInvocations));
	  upper_spacing = 0.1 / ((double)(tr->optimizeRateCategoryInvocations));
	}
      
      if(lower_spacing < 0.001)
	lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
	upper_spacing = 0.001;
      
      tr->optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations + 1;

      memcpy(oldCategory, tr->rateCategory, sizeof(int) * tr->originalCrunchedLength);	     
      memcpy(ratStored,   tr->patratStored, sizeof(double) * tr->originalCrunchedLength);

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldNumbers[model]          = tr->partitionData[model].numberOfCategories;

	  oldCategorizedRates[model] = (double *)malloc(sizeof(double) * tr->maxCategories);
	  
	  memcpy(oldCategorizedRates[model], tr->partitionData[model].perSiteRates, tr->maxCategories * sizeof(double));	  	 	  
	}      
      
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      /*tr->lhs = lhs;*/
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      masterBarrier(THREAD_RATE_CATS, tr);      
#else      
      for(model = 0; model < tr->NumberOfModels; model++)      
	optRateCatModel(tr, model, lower_spacing, upper_spacing, tr->lhs);
#endif     

      for(model = 0; model < tr->NumberOfModels; model++)
	{     
	  int 
	    where = 1,
	    found = 0,
	    width = tr->partitionData[model].upper -  tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper,
	    lower = tr->partitionData[model].lower;
	    
	  rateCategorize 
	    *rc = (rateCategorize *)malloc(sizeof(rateCategorize) * width);		 
	
	  for (i = 0; i < width; i++)
	    {
	      rc[i].accumulatedSiteLikelihood = 0.0;
	      rc[i].rate = 0.0;
	    }  
	
	  rc[0].accumulatedSiteLikelihood = tr->lhs[lower];
	  rc[0].rate = tr->patrat[lower];
	
	  tr->rateCategory[lower] = 0;
	
	  for (i = lower + 1; i < upper; i++) 
	    {
	      temp = tr->patrat[i];
	      found = 0;
	    
	      for(k = 0; k < where; k++)
		{
		  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
		    {
		      found = 1;						
		      rc[k].accumulatedSiteLikelihood += tr->lhs[i];	
		      break;
		    }
		}
	    
	      if(!found)
		{	    
		  rc[where].rate = temp;	    
		  rc[where].accumulatedSiteLikelihood += tr->lhs[i];	    
		  where++;
		}
	    }
	
	  qsort(rc, where, sizeof(rateCategorize), catCompare);
	
	  if(where < maxCategories)
	    {
	      tr->partitionData[model].numberOfCategories = where;
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	  else
	    {
	      tr->partitionData[model].numberOfCategories = maxCategories;	
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	
	  free(rc);
	}
        	
      updatePerSiteRates(tr, TRUE);	

      evaluateGeneric(tr, tr->start, TRUE);
      
      if(tr->likelihood < initialLH)
	{	 		  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      tr->partitionData[model].numberOfCategories = oldNumbers[model];
	      memcpy(tr->partitionData[model].perSiteRates, oldCategorizedRates[model], tr->maxCategories * sizeof(double));
	    }	      
	  
	  memcpy(tr->patratStored, ratStored, sizeof(double) * tr->originalCrunchedLength);
	  memcpy(tr->rateCategory, oldCategory, sizeof(int) * tr->originalCrunchedLength);	     
	  
	  updatePerSiteRates(tr, FALSE);
	  
	  evaluateGeneric(tr, tr->start, TRUE);

	  /* printf("REVERT: %1.40f %1.40f\n", initialLH, tr->likelihood); */

	  assert(initialLH == tr->likelihood);
	}
          
      for(model = 0; model < tr->NumberOfModels; model++)
	free(oldCategorizedRates[model]);
                   
      free(oldCategorizedRates);
      free(oldCategory);
      free(ratStored);       
      /*      free(lhs); */
      free(oldNumbers);
    }
}
  

/************************* end of functions for CAT model of rate heterogeneity */




/*****************************************************************************************************/

/* reset all branche lengths in tree to default values */

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < tr->numBranches; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < tr->numBranches; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}

/* print the protein GTR substitution matrix to file. This is only printed if we are actually estimating a 
   GTR model for all protein data partitions in the data */

static void printAAmatrix(tree *tr, double epsilon)
{
  if(AAisGTR(tr))
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].dataType == AA_DATA) 
	    {
	      char gtrFileName[1024];
	      char epsilonStr[1024];
	      FILE *gtrFile;
	      double *rates = tr->partitionData[model].substRates;
	      double *f     = tr->partitionData[model].frequencies;
	      double q[20][20];
	      int    r = 0;
	      int i, j;

	      assert(tr->partitionData[model].protModels == GTR);

	      sprintf(epsilonStr, "%f", epsilon);

	      strcpy(gtrFileName, workdir);
	      strcat(gtrFileName, "RAxML_proteinGTRmodel.");
	      strcat(gtrFileName, run_id);
	      strcat(gtrFileName, "_");
	      strcat(gtrFileName, epsilonStr);

	      gtrFile = myfopen(gtrFileName, "wb");

	      for(i = 0; i < 20; i++)
		for(j = 0; j < 20; j++)
		  q[i][j] = 0.0;

	      for(i = 0; i < 19; i++)
		for(j = i + 1; j < 20; j++)
		  q[i][j] = rates[r++];

	      for(i = 0; i < 20; i++)
		for(j = 0; j <= i; j++)
		  {
		    if(i == j)
		      q[i][j] = 0.0;
		    else
		      q[i][j] = q[j][i];
		  }
	   
	      for(i = 0; i < 20; i++)
		{
		  for(j = 0; j < 20; j++)		
		    fprintf(gtrFile, "%1.80f ", q[i][j]);
		
		  fprintf(gtrFile, "\n");
		}
	      for(i = 0; i < 20; i++)
		fprintf(gtrFile, "%1.80f ", f[i]);
	      fprintf(gtrFile, "\n");

	      fclose(gtrFile);

	      printBothOpen("\nPrinted intermediate AA substitution matrix to file %s\n\n", gtrFileName);
	      
	      break;
	    }

	}	  
    }
}

/* 
   automatically compute the best protein substitution model for the dataset at hand.
 */

static void autoProtein(tree *tr)
{
  int 
    countAutos = 0,
    i,
    model;
    
   for(model = 0; model < tr->NumberOfModels; model++)	      
     if(tr->partitionData[model].protModels == AUTO)
       countAutos++;
  
  if(countAutos > 0)
    {
      int 
	numProteinModels = AUTO,
	*bestIndex = (int*)malloc(sizeof(int) * tr->NumberOfModels);

      double
	*bestScores = (double*)malloc(sizeof(double) * tr->NumberOfModels);

      /*printf("Entry: %f\n", tr->likelihood);*/
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  bestIndex[model] = -1;
	  bestScores[model] = unlikely;
	}

      for(i = 0; i < numProteinModels; i++)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	   
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = i;

		  initReversibleGTR(tr, model);  
		}
	    }
	  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	  masterBarrier(THREAD_COPY_RATES, tr);     
#endif
	  
	  resetBranches(tr);
	  evaluateGeneric(tr, tr->start, TRUE);  
	  treeEvaluate(tr, 16);// 0.5 * 32 = 16.0   

	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  if(tr->perPartitionLH[model] > bestScores[model])
		    {
		      bestScores[model] = tr->perPartitionLH[model];
		      bestIndex[model] = i;		      
		    }
		}

	    }       
	}

      printBothOpen("\n\n");
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    {
	      tr->partitionData[model].autoProtModels = bestIndex[model];

	      initReversibleGTR(tr, model);  

	      printBothOpen("Partition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[tr->partitionData[model].autoProtModels], bestScores[model]);
	    }	 
	}
      
      printBothOpen("\n\n");
            
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      masterBarrier(THREAD_COPY_RATES, tr);     
#endif

      resetBranches(tr);
      evaluateGeneric(tr, tr->start, TRUE); 
      treeEvaluate(tr, 16); // 0.5 * 32 = 16
      
      /*printf("Exit: %f\n", tr->likelihood);*/
      
      free(bestIndex);
      free(bestScores);
    }
}

/* iterative procedure for optimizing all model parameters */

void modOpt(tree *tr, double likelihoodEpsilon)
{ 
  int i, catOpt = 0; 
  double 
    currentLikelihood,
    modelEpsilon = 0.0001;

  /* linkage lists for alpha, p-invar has actually been ommitted in this version of the code 
     and the GTR subst matrices */

  linkageList *alphaList;
  linkageList *invarList;
  linkageList *rateList; 

  int *unlinked = (int *)malloc(sizeof(int) * tr->NumberOfModels);

  modelEpsilon = 0.0001;

  /* set the integer vector for linking parameters to all parameters being unlinked, 
     i.e. a separate/independent estimate of parameters will be conducted for eahc partition */

  for(i = 0; i < tr->NumberOfModels; i++)
    unlinked[i] = i;

  /* alpha parameters and p-invar parameters are unlinked.
     this is the point where I actually hard-coded this in RAxML */

  alphaList = initLinkageList(unlinked, tr);
  invarList = initLinkageList(unlinked, tr);

  /* call the dedicated function for linking the GTR matrix across all AA data partitions 
     If we have only DNA data all GTR matrix estimates will be unlinked.
     */

  rateList  = initLinkageListGTR(tr);

  tr->start = tr->nodep[1];

  do
  {           
    //printBothOpen("cur LH: %f\n", tr->likelihood);
    currentLikelihood = tr->likelihood;     

    optRatesGeneric(tr, modelEpsilon, rateList);

    evaluateGeneric(tr, tr->start, TRUE);                                       

    autoProtein(tr);

    treeEvaluate(tr, 2); // 0.0625 * 32 = 2.0    

    switch(tr->rateHetModel)
    {
      case GAMMA:      
        optAlpha(tr, modelEpsilon, alphaList); 
        evaluateGeneric(tr, tr->start, TRUE); 	 	 
        treeEvaluate(tr, 3); // 0.1 * 32 = 3.2  	 
        break;
      case CAT:
        if(catOpt < 3)
        {	      	     	     
          optimizeRateCategories(tr, tr->categories);	      	     	      	      
          catOpt++;
        }
        break;	  
      default:
        assert(0);
    }                   



    printAAmatrix(tr, fabs(currentLikelihood - tr->likelihood));    
  }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  

  free(unlinked);
  freeLinkageList(alphaList);
  freeLinkageList(rateList);
  freeLinkageList(invarList);  
}

