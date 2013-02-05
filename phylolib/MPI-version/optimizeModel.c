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
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[20];


extern int processes;
extern int processID;

/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


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





static linkageList* initLinkageList(int *linkList, tree *tr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  linkageList* ll = (linkageList*)malloc(sizeof(linkageList));
      
  for(i = 0; i < tr->NumberOfModels; i++)
    {
      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)malloc(sizeof(linkageData) * numberOfModels);


  for(i = 0; i < numberOfModels; i++)
    {
      ll->ld[i].valid = TRUE;
      partitions = 0;

      for(k = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  partitions++;	    

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)malloc(sizeof(int) * partitions);
      
      for(k = 0, pos = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  return ll;
}


static linkageList* initLinkageListGTR(tree *tr)
{
  int
    i,
    *links = (int*)malloc(sizeof(int) * tr->NumberOfModels),
    firstAA = tr->NumberOfModels + 2,
    countGTR = 0,
    countOtherModel = 0;
  linkageList* ll;

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

  if(countGTR == 0)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	links[i] = i;
    }
  else
    {
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
  

  ll = initLinkageList(links, tr);

  free(links);
  
  return ll;
}



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



static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll)
{ 
  int i, k, pos;

  switch(whichFunction)
    {    
    case RATE_F:     
      for(i = 0, pos = 0; i < ll->entries; i++)
	{
	  if(ll->ld[i].valid)
	    {
	      if(converged[pos])
		{
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
		}
	      else
		{
		  for(k = 0; k < ll->ld[i].partitions; k++)
		    {
		      int index = ll->ld[i].partitionList[k];		  	      
		      setRateModel(tr, index, value[pos], rateNumber);  

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


      evaluateGeneric(tr, tr->start, TRUE);  
      
   
      
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

		  makeGammaCats(tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useMedian);
		}
	    }
	}

      evaluateGeneric(tr, tr->start, TRUE);
            
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


   evaluateGeneric(tr, tr->start, TRUE);
   /* 
     at this point here every worker has the traversal data it needs for the 
     search, so we won't re-distribute it he he :-)
  */

  for(i = 0; i < numberOfModels; i++)
    {
      assert(ll->ld[i].valid);

      startAlpha[i] = tr->partitionData[ll->ld[i].partitionList[0]].alpha;
      _a[i] = startAlpha[i] + 0.1;
      _b[i] = startAlpha[i] - 0.1;      
      if(_b[i] < lim_inf) 
	_b[i] = lim_inf;

      startLH[i] = 0.0;
      
      for(k = 0; k < ll->ld[i].partitions; k++)	
	{
	  startLH[i] += tr->perPartitionLH[ll->ld[i].partitionList[k]];
	  assert(tr->partitionData[ll->ld[i].partitionList[0]].alpha ==  tr->partitionData[ll->ld[i].partitionList[k]].alpha);
	}
    }					  

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, -1, ALPHA_F, tr, ll);       
  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, result, numberOfModels, ALPHA_F, -1, tr, ll, lim_inf, lim_sup);

  for(i = 0; i < numberOfModels; i++)
    endAlpha[i] = result[i];
  
  for(i = 0; i < numberOfModels; i++)
    {
      if(startLH[i] > endAlpha[i])
	{    	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {	      
	      tr->partitionData[ll->ld[i].partitionList[k]].alpha = startAlpha[i];
#ifndef _LOCAL_DISCRETIZATION
	      makeGammaCats(tr->partitionData[ll->ld[i].partitionList[k]].alpha, tr->partitionData[ll->ld[i].partitionList[k]].gammaRates, 4, tr->useMedian); 	      
#endif		
	    }

	}  
    }



  
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


   assert(states != -1);

  startRates = (double *)malloc(sizeof(double) * numberOfRates * numberOfModels);

  evaluateGeneric(tr, tr->start, TRUE);
  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  for(j = 0; j < ll->ld[i].partitions; j++)
	    {
	      int index = ll->ld[i].partitionList[j];
	      
	      startLH[pos] += tr->perPartitionLH[index];
	      for(k = 0; k < numberOfRates; k++)
		startRates[pos * numberOfRates + k] = tr->partitionData[index].substRates[k];      
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
  
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
      
      for(k = 0, pos = 0; k < ll->entries; k++)
	{

	  if(ll->ld[k].valid)
	    { 
	      if(startLH[pos] > endLH[pos])
		{		  
		  for(j = 0; j < ll->ld[k].partitions; j++)
		    {
		      int index = ll->ld[k].partitionList[j];
		      tr->partitionData[index].substRates[i] = startRates[pos * numberOfRates + i];
#ifndef _LOCAL_DISCRETIZATION 	             	  
		      initReversibleGTR(tr, index);		     
#endif
		    }

		}
	      pos++;
	    }
	}


      
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

  /* first do DNA */

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

  if(dnaPartitions > 0)
    optRates(tr, modelEpsilon, ll, dnaPartitions, states);
  

  /* then SECONDARY */

   for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	  /* we only have one type of secondary data models in one analysis */
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
  
  /* then multi-state */

  /* 
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

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}





/*********************FUNCTIONS FOR APPROXIMATE MODEL OPTIMIZATION ***************************************/






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




static void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model, 
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      int 
	localIndex = 0;
     
      for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	{
	  if(i % n == tid)
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






/* 
   set scaleRates to FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

static size_t calcSendBufferSize(tree *tr, int tid)
{
  size_t 
    length = 0;

  int 
    model;  

  for(model = 0; model < tr->NumberOfModels; model++)
    if(isThisMyPartition(tr, tid, model))      
      length += ((size_t)tr->partitionData[model].upper - (size_t)tr->partitionData[model].lower);            

  return length;
}

static void gatherCatsWorker(tree *tr, int tid)
{  
  size_t 
    model,   
    sendBufferSize = calcSendBufferSize(tr, tid);

  int   
    *catBufSend = (int *)malloc(sendBufferSize * sizeof(int));

  double 
    *rateBufSend = (double *)malloc(sendBufferSize * sizeof(double)),
    *patBufSend = (double *)malloc(sendBufferSize * sizeof(double)),
    *patStoredBufSend =  (double *)malloc(sendBufferSize * sizeof(double)),
    *lhsBufSend = (double *)malloc(sendBufferSize * sizeof(double));

  size_t
    offsets = 0; 

  for(model = 0, offsets = 0; model < (size_t)tr->NumberOfModels; model++)
    {               	
      size_t
	start = (size_t)tr->partitionData[model].lower,
	width = (size_t)tr->partitionData[model].upper - (size_t)tr->partitionData[model].lower;
      
      if(isThisMyPartition(tr, tid, model))
	{	 	  	  
	  memcpy(&catBufSend[offsets],       &tr->rateCategory[start], sizeof(int) * width);
	  memcpy(&rateBufSend[offsets],      tr->partitionData[model].perSiteRates, sizeof(double) * tr->maxCategories);
	  memcpy(&patBufSend[offsets],       &tr->patrat[start],       sizeof(double) * width);
	  memcpy(&patStoredBufSend[offsets], &tr->patratStored[start], sizeof(double) * width);
	  memcpy(&lhsBufSend[offsets],       &tr->lhs[start],          sizeof(double) * width);
	  
	  offsets += width;	  
	}		 		
    }
    
  MPI_Gatherv(catBufSend,       sendBufferSize, MPI_INT,    (int*)NULL, (int*)NULL, (int*)NULL, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Gatherv(rateBufSend,      sendBufferSize, MPI_DOUBLE, (double*)NULL, (int*)NULL, (int*)NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(patBufSend,       sendBufferSize, MPI_DOUBLE, (double*)NULL, (int*)NULL, (int*)NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(patStoredBufSend, sendBufferSize, MPI_DOUBLE, (double*)NULL, (int*)NULL, (int*)NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(lhsBufSend,       sendBufferSize, MPI_DOUBLE, (double*)NULL, (int*)NULL, (int*)NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  free(catBufSend);
  free(rateBufSend);
  free(patBufSend);
  free(patStoredBufSend);
  free(lhsBufSend);
}


static void gatherCatsMaster(tree *tr, int tid, int n)
{
  int  
    offsets = 0,
    *countArray  = (int *)calloc(n, sizeof(int)),
    *offsetArray = (int *)calloc(n, sizeof(int)),
    *modelOffsets = (int *)calloc(tr->NumberOfModels, sizeof(int)), 
    i,
    p;
  
  size_t
    model,
    sendBufferSize = calcSendBufferSize(tr, tid),
    recvBufferSize = (size_t)tr->originalCrunchedLength;

  int   
    *catBufSend = (int *)malloc(sendBufferSize * sizeof(int)),
    *catBufRecv = (int *)malloc(recvBufferSize * sizeof(int));

  double   
    *rateBufSend =       (double *)malloc(sendBufferSize * sizeof(double)),   
    *patBufSend =        (double *)malloc(sendBufferSize * sizeof(double)),
    *patStoredBufSend =  (double *)malloc(sendBufferSize * sizeof(double)),
    *lhsBufSend =        (double *)malloc(sendBufferSize * sizeof(double)),
    
     *rateBufRecv =      (double *)malloc(recvBufferSize * sizeof(double)),
    *patBufRecv =        (double *)malloc(recvBufferSize * sizeof(double)),
    *patStoredBufRecv =  (double *)malloc(recvBufferSize * sizeof(double)),
    *lhsBufRecv =        (double *)malloc(recvBufferSize * sizeof(double));

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    countArray[tr->partitionAssignment[model]] += (int)(tr->partitionData[model].upper - tr->partitionData[model].lower);

  for(i = 0, offsets = 0; i < n; i++)
    {
      offsetArray[i] = offsets;
      offsets += countArray[i];
    }

  for(p = 0; p < n; p++)    
    {
      int
	localOffset = 0,
	globalOffset = offsetArray[p];
      
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	if(tr->partitionAssignment[model] == p)
	  {
	    modelOffsets[model] = globalOffset + localOffset;
	    localOffset += (tr->partitionData[model].upper - tr->partitionData[model].lower);
	  }
    }


  assert((size_t)offsets == recvBufferSize); 

  for(model = 0, offsets = 0; model < (size_t)tr->NumberOfModels; model++)
    {               	        
      size_t       
	start = (size_t)tr->partitionData[model].lower,
	width = (size_t)tr->partitionData[model].upper - (size_t)tr->partitionData[model].lower;
      
      if(isThisMyPartition(tr, tid, model))
	{		  	  
	  memcpy(&catBufSend[offsets],       &tr->rateCategory[start], sizeof(int)    * width);
	  memcpy(&rateBufSend[offsets],      tr->partitionData[model].perSiteRates, sizeof(double) * tr->maxCategories);
	  memcpy(&patBufSend[offsets],       &tr->patrat[start],       sizeof(double) * width);
	  memcpy(&patStoredBufSend[offsets], &tr->patratStored[start], sizeof(double) * width);
	  memcpy(&lhsBufSend[offsets],       &tr->lhs[start],          sizeof(double) * width);	 

	  offsets += width; 	 	 
	}	      
    }
   
  MPI_Gatherv(catBufSend,       sendBufferSize, MPI_INT,    catBufRecv,       countArray, offsetArray, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Gatherv(rateBufSend,      sendBufferSize, MPI_DOUBLE, rateBufRecv,      countArray, offsetArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(patBufSend,       sendBufferSize, MPI_DOUBLE, patBufRecv,       countArray, offsetArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(patStoredBufSend, sendBufferSize, MPI_DOUBLE, patStoredBufRecv, countArray, offsetArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(lhsBufSend,       sendBufferSize, MPI_DOUBLE, lhsBufRecv,       countArray, offsetArray, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {        
      size_t
	start  = (size_t)tr->partitionData[model].lower,
	width  = (size_t)tr->partitionData[model].upper - (size_t)tr->partitionData[model].lower;            
      
      memcpy(&tr->rateCategory[start], &catBufRecv[modelOffsets[model]],       sizeof(int) * width);
      memcpy(tr->partitionData[model].perSiteRates,       &rateBufRecv[modelOffsets[model]],       sizeof(double) * tr->maxCategories);
      memcpy(&tr->patrat[start],       &patBufRecv[modelOffsets[model]],       sizeof(double) * width);
      memcpy(&tr->patratStored[start], &patStoredBufRecv[modelOffsets[model]], sizeof(double) * width);
      memcpy(&tr->lhs[start],          &lhsBufRecv[modelOffsets[model]],       sizeof(double) * width); 
      
      {	
	int 	  
	  *numCAT = (int *)calloc(tr->maxCategories, sizeof(int)),
	  k, 
	  numCats = 0;

	for(k = tr->partitionData[model].lower; k < tr->partitionData[model].upper; k++)
	  numCAT[tr->rateCategory[k]] = 1;

	for(k = 0; k < tr->maxCategories; k++)
	  if(numCAT[k] == 1)
	    numCats++;
	
	if(isThisMyPartition(tr, processID, model))
	  assert(tr->partitionData[model].numberOfCategories == numCats);

	tr->partitionData[model].numberOfCategories = numCats;

	

	free(numCAT);
      }

    }
  
  /* 
     {
    double sum = 0.0;
    
    for(model = 0; model < (size_t)tr->NumberOfModels; model++)
      {        
	int i;
	printf("\nModel %d cats %d\n", model, tr->partitionData[model].numberOfCategories);
	for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)	
	  {
	    printf("%d %f %f %f\n", tr->rateCategory[i], tr->patrat[i], tr->patratStored[i], tr->lhs[i]);     
	    sum += tr->lhs[i];
	    }	
      }
    
    printf("likelihood: %f\n", sum);
    }
  */
  
  free(modelOffsets);
  free(countArray);
  free(offsetArray);
 
  free(catBufSend);
  free(rateBufSend);
  free(patBufSend);
  free(patStoredBufSend);
  free(lhsBufSend);

  free(catBufRecv);
  free(rateBufRecv);
  free(patBufRecv);
  free(patStoredBufRecv);
  free(lhsBufRecv);
}



static void updatePerSiteRatesManyPartitions(tree *tr, boolean scaleRates)
{
  int 
    i,
    model; 

  if(tr->numBranches > 1)
    {            
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(isThisMyPartition(tr, processID, model))
	    {
	      int 	       
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper,	      
		localCount = 0;
	      
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
	      
	      for(i = lower; i < upper; i++)
		{
		  double
		    w = ((double)tr->aliaswgt[i]);	      
		  
		  double
		    wtemp,
		    temp = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  wtemp = temp * w;
		  
		  tr->wr[i]  = wtemp;
		  tr->wr2[i] = temp * wtemp;
		}	            	  	  	      
	
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{	    	      
		  tr->partitionData[model].wr[localCount]  = tr->wr[i];
		  tr->partitionData[model].wr2[localCount] = tr->wr2[i];
		  tr->partitionData[model].rateCategory[localCount] = tr->rateCategory[i];		
		}	      
	    }	 
	}
    }
  else
    {
      int
	accWgt = 0;

      double 
	scaler = 0.0,       
	accRat = 0.0,
	dwgt,
	a[2],
	r[2];
      
      if(scaleRates)
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, processID, model))
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
	    }
	  	  	    	      
	  a[0] = (double)accRat;
	  a[1] = (double)accWgt;

	  MPI_Reduce(a, r, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Bcast(r, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	      
	  accRat = r[0];
	  dwgt   = r[1];
	  
	  accRat /= dwgt;	     	    
	 	  
	  scaler = 1.0 / ((double)accRat);

	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, processID, model))
		{
		  for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		    tr->partitionData[model].perSiteRates[i] *= scaler;
		}
	    }

	  for(model = 0, accRat = 0.0; model < tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, processID, model))
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
	    }           
	
	 
	  a[0] = (double)accRat;
	  a[1] = (double)accWgt;
	  
	  MPI_Reduce(a, r, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Bcast(r, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  
	  accRat = r[0];
	  dwgt   = r[1];
	  
	  accRat /= dwgt;	     	     		  

	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
      else
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, processID, model))
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
	    }
	  
	   

	  a[0] = (double)accRat;
	  a[1] = (double)accWgt;
	  
	  MPI_Reduce(a, r, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Bcast(r, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  
	  accRat = r[0];
	  dwgt   = r[1];
	  
	  accRat /= dwgt;	     
	 
	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
      
       for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(isThisMyPartition(tr, processID, model))
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  double
		    w = ((double)tr->aliaswgt[i]);	      
		  
		  double
		    wtemp,
		    temp = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
		  
		  wtemp = temp * w;
		  
		  tr->wr[i]  = wtemp;
		  tr->wr2[i] = temp * wtemp;
		}
	    }
	}         

      for(model = 0; model < tr->NumberOfModels; model++)
	{  
	  if(isThisMyPartition(tr, processID, model))	  	  	 
	    {
	      int 
		localCount,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;	  
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{	    	      
		  tr->partitionData[model].wr[localCount]  = tr->wr[i];
		  tr->partitionData[model].wr2[localCount] = tr->wr2[i];
		  tr->partitionData[model].rateCategory[localCount] = tr->rateCategory[i];
		}
	    }
	}
    }  

  if(processID == 0)    
    gatherCatsMaster(tr, processID, processes);    
  else    
    gatherCatsWorker(tr, processID);    
}

static void gatherRatesFewPartitions(tree *tr, int tid)
{
  int n = processes;

  if(tid == 0)
    {
       int
          localCounter,
          model,
          i,
          sendBufferSize = (tr->originalCrunchedLength / n) + 1,
          recvBufferSize = sendBufferSize * n;

        double 	  	
          *patBufSend = (double *)malloc(sendBufferSize * sizeof(double)),
          *patStoredBufSend =  (double *)malloc(sendBufferSize * sizeof(double)),
          *lhsBufSend = (double *)malloc(sendBufferSize * sizeof(double)),
          *patBufRecv = (double *)malloc(recvBufferSize * sizeof(double)),
          *patStoredBufRecv =  (double *)malloc(recvBufferSize * sizeof(double)),
          *lhsBufRecv = (double *)malloc(recvBufferSize * sizeof(double));                

        for(model = 0, localCounter = 0; model < tr->NumberOfModels; model++)
	  {               	      	     
	    for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	      if(i % n == tid)
		{
		  patBufSend[localCounter] = tr->patrat[i];
		  patStoredBufSend[localCounter] = tr->patratStored[i];
		  lhsBufSend[localCounter] = tr->lhs[i];
		  localCounter++;
		}	    
	  }

        MPI_Gather(patBufSend,       sendBufferSize, MPI_DOUBLE, patBufRecv,       sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(patStoredBufSend, sendBufferSize, MPI_DOUBLE, patStoredBufRecv, sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(lhsBufSend,       sendBufferSize, MPI_DOUBLE, lhsBufRecv,       sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

        for(model = 0; model < tr->NumberOfModels; model++)
	  {   	      
	    for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	      {
		int 
		  offset = i % n,
		  position = i / n;
		
		tr->patrat[i]       = patBufRecv[offset * sendBufferSize + position];
		tr->patratStored[i] = patStoredBufRecv[offset * sendBufferSize + position];
		tr->lhs[i]                = lhsBufRecv[offset * sendBufferSize + position];		     
	      }	   	   
	  }

        free(patBufSend);
        free(patStoredBufSend);
        free(lhsBufSend);
        free(patBufRecv);
        free(patStoredBufRecv);
        free(lhsBufRecv);
    }
  else
    {
      int 
	model,
	i,
	localCounter,
	sendBufferSize = (tr->originalCrunchedLength / n) + 1;
      
      double 
	*localDummy = (double*)NULL,	      
	*patBufSend = (double *)malloc(sendBufferSize * sizeof(double)),
	*patStoredBufSend =  (double *)malloc(sendBufferSize * sizeof(double)),
	*lhsBufSend = (double *)malloc(sendBufferSize * sizeof(double));
      
      for(model = 0, localCounter = 0; model < tr->NumberOfModels; model++)
	{               		  		  
	  for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	    if(i % n == tid)
	      {
		patBufSend[localCounter] = tr->patrat[i];
		patStoredBufSend[localCounter] = tr->patratStored[i];
		lhsBufSend[localCounter] = tr->lhs[i];
		localCounter++;
	      }		  
	}
      
      MPI_Gather(patBufSend,       sendBufferSize, MPI_DOUBLE, localDummy, sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(patStoredBufSend, sendBufferSize, MPI_DOUBLE, localDummy, sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(lhsBufSend,       sendBufferSize, MPI_DOUBLE, localDummy, sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      free(patBufSend);
      free(patStoredBufSend);
      free(lhsBufSend);
    }
}

static void broadcastRatesFewPartitions(tree *tr, int tid)
{ 
  int 
    model,
    n = processes;
  
  MPI_Bcast(tr->patrat, tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tr->patratStored, tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for(model = 0; model < tr->NumberOfModels; model++)
    { 
      MPI_Bcast(&(tr->partitionData[model].numberOfCategories), 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(tr->partitionData[model].perSiteRates, tr->partitionData[model].numberOfCategories, MPI_DOUBLE, 0, MPI_COMM_WORLD);	   
    }


  MPI_Bcast(tr->rateCategory, tr->originalCrunchedLength, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(tr->wr,                 tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(tr->wr2,                tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int
	i,
	localCounter;

      for(localCounter = 0, i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	{
	  if(i % n == tid)
	    {		 
	      tr->partitionData[model].rateCategory[localCounter] = tr->rateCategory[i];
	      tr->partitionData[model].wr[localCounter]             = tr->wr[i];
	      tr->partitionData[model].wr2[localCounter]            = tr->wr2[i];		 

	      localCounter++;
	    }
	}
    }
 
  MPI_Barrier(MPI_COMM_WORLD);
}

static void updatePerSiteRatesFewPartitions(tree *tr, boolean scaleRates)
{  
  int 
    i,
    model;

  /*gatherRatesFewPartitions(tr, processID);*/
  
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
    }
      
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 	    
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;
      
      for(i = lower; i < upper; i++)
	{
	  double
	    w = ((double)tr->aliaswgt[i]);	      
	  
	  double
	    wtemp,
	    temp = tr->partitionData[model].perSiteRates[tr->rateCategory[i]];
	  
	  wtemp = temp * w;
	  
	  tr->wr[i]  = wtemp;
	  tr->wr2[i] = temp * wtemp;
	}
    }         
  
  broadcastRatesFewPartitions(tr, processID);
}

void updatePerSiteRates(tree *tr, boolean scaleRates)
{
  if(tr->rateHetModel != CAT)
    return;

  if(tr->manyPartitions)
    updatePerSiteRatesManyPartitions(tr, scaleRates);
  else
    updatePerSiteRatesFewPartitions(tr, scaleRates);
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

      if(optimizeRateCategoryInvocations == 1)
	{
	  lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
	}
      else
	{
	  lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
	}
      
      if(lower_spacing < 0.001)
	lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
	upper_spacing = 0.001;
      
      optimizeRateCategoryInvocations++;

      memcpy(oldCategory, tr->rateCategory, sizeof(int) * tr->originalCrunchedLength);	     
      memcpy(ratStored,   tr->patratStored, sizeof(double) * tr->originalCrunchedLength);

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldNumbers[model]          = tr->partitionData[model].numberOfCategories;

	  oldCategorizedRates[model] = (double *)malloc(sizeof(double) * tr->maxCategories);
	  
	  memcpy(oldCategorizedRates[model], tr->partitionData[model].perSiteRates, tr->maxCategories * sizeof(double));	  	 	  
	}            	     

      if(tr->manyPartitions)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)      
	    if(isThisMyPartition(tr, processID, model))
	       optRateCatModel(tr, model, lower_spacing, upper_spacing, tr->lhs);
	}
      else
	{
	  optRateCatPthreads(tr, lower_spacing, upper_spacing, tr->lhs, processes, processID);
	  gatherRatesFewPartitions(tr, processID);
	}     
           
      for(model = 0; model < tr->NumberOfModels; model++)
	{    
	  boolean 
	    execute;
	  
	  if(tr->manyPartitions)
	    execute = isThisMyPartition(tr, processID, model);
	  else
	    execute = (processID == 0);
 
	  if(execute)
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
	}            

      updatePerSiteRates(tr, TRUE);	
     
      /*
	if(processID == 0)
	{	  
	  double sum = 0.0;
	  int model;
	  
	  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	    {        
	      int i;
	      printf("\nModel %d cats %d\n", model, tr->partitionData[model].numberOfCategories);
	      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)	
		{
		  printf("%d %f %f %f\n", tr->rateCategory[i], tr->patrat[i], tr->patratStored[i], tr->lhs[i]); 
		  sum += tr->lhs[i];
		}	
	    }	
	  
	  printf("likelihood: %f\n", sum);
	}
      */
     
      evaluateGeneric(tr, tr->start, TRUE);

      /* printf("%f \n", tr->likelihood); */
      
      /*MPI_Finalize();*/
      
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
  






/*****************************************************************************************************/

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
#ifndef _LOCAL_DISCRETIZATION 
		  initReversibleGTR(tr, model);  
#endif
		}
	    }
	  

	  
	  resetBranches(tr);
	  evaluateGeneric(tr, tr->start, TRUE);  
	  treeEvaluate(tr, 0.5);     

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
#ifndef _LOCAL_DISCRETIZATION 
	      initReversibleGTR(tr, model);  
#endif
	      printBothOpen("Partition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[tr->partitionData[model].autoProtModels], bestScores[model]);
	    }	 
	}
      
      printBothOpen("\n\n");
            


      resetBranches(tr);
      evaluateGeneric(tr, tr->start, TRUE); 
      treeEvaluate(tr, 0.5);
      
      /*printf("Exit: %f\n", tr->likelihood);*/
      
      free(bestIndex);
      free(bestScores);
    }
}



void modOpt(tree *tr, double likelihoodEpsilon)
{ 
  int 
    i, 
    catOpt = 0,
    *unlinked = (int *)malloc(sizeof(int) * tr->NumberOfModels);
  
  double 
    currentLikelihood,
    modelEpsilon = 0.0001;
  
  linkageList 
    *alphaList,
    *invarList,
    *rateList;      

  for(i = 0; i < tr->NumberOfModels; i++)
    unlinked[i] = i;
  
  alphaList = initLinkageList(unlinked, tr);
  invarList = initLinkageList(unlinked, tr);
  rateList  = initLinkageListGTR(tr);
   
  tr->start = tr->nodep[1];
                 
  do
    {           
      currentLikelihood = tr->likelihood;     
     
     

      optRatesGeneric(tr, modelEpsilon, rateList);
      
      
      
      

      evaluateGeneric(tr, tr->start, TRUE);                                       
     
      
     

      autoProtein(tr);

      treeEvaluate(tr, 0.0625);      

      

      switch(tr->rateHetModel)
	{
	case GAMMA:      
	  optAlpha(tr, modelEpsilon, alphaList); 
	  
	  evaluateGeneric(tr, tr->start, TRUE); 	 	 
	  
	  treeEvaluate(tr, 0.1);	  	 
	 
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

