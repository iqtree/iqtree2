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
 * @file treeIO.c
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

extern char *likelihood_key;
extern char *ntaxa_key;
extern char *smoothed_key;
extern int partCount;

int countTips(nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))  
    return 1;    
  {
    nodeptr q;
    int tips = 0;

    q = p->next;
    while(q != p)
      { 
	tips += countTips(q->back, numsp);
	q = q->next;
      } 
    
    return tips;
  }
}


static double getBranchLength(pllInstance *tr, partitionList *pr, int perGene, nodeptr p)
{
  double 
    z = 0.0,
    x = 0.0;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  assert(perGene != PLL_NO_BRANCHES);
	      
  if(numBranches == 1)
    {
      assert(tr->fracchange != -1.0);
      z = p->z[0];
      if (z < PLL_ZMIN) 
	z = PLL_ZMIN;      	 
      
      x = -log(z) * tr->fracchange;           
    }
  else
    {
      if(perGene == PLL_SUMMARIZE_LH)
	{
	  int 
	    i;
	  
	  double 
	    avgX = 0.0;
		      
	  for(i = 0; i < numBranches; i++)
	    {
	      assert(pr->partitionData[i]->partitionContribution != -1.0);
	      assert(pr->partitionData[i]->fracchange != -1.0);
	      z = p->z[i];
	      if(z < PLL_ZMIN) 
		z = PLL_ZMIN;      	 
	      x = -log(z) * pr->partitionData[i]->fracchange;
	      avgX += x * pr->partitionData[i]->partitionContribution;
	    }

	  x = avgX;
	}
      else
	{	
	  assert(pr->partitionData[perGene]->fracchange != -1.0);
	  assert(perGene >= 0 && perGene < numBranches);
	  
	  z = p->z[perGene];
	  
	  if(z < PLL_ZMIN) 
	    z = PLL_ZMIN;      	 
	  
	  x = -log(z) * pr->partitionData[perGene]->fracchange;
	}
    }

  return x;
}

static char *pllTreeToNewickREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean printBranchLengths, pllBoolean printNames,
			    pllBoolean printLikelihood, pllBoolean rellTree, pllBoolean finalPrint, int perGene, pllBoolean branchLabelSupport, pllBoolean printSHSupport)
{
  char  *nameptr;            
      
  if(isTip(p->number, tr->mxtips)) 
    {	       	  
      if(printNames)
	{
	  nameptr = tr->nameList[p->number];     
	  sprintf(treestr, "%s", nameptr);
	}
      else
	sprintf(treestr, "%d", p->number);    
	
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = pllTreeToNewickREC(treestr, tr, pr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      *treestr++ = ',';
      treestr = pllTreeToNewickREC(treestr, tr, pr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = pllTreeToNewickREC(treestr, tr, pr, p->back, printBranchLengths, printNames, printLikelihood, rellTree,
				   finalPrint, perGene, branchLabelSupport, printSHSupport);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {	      	 
      if(printBranchLengths && !rellTree)
	sprintf(treestr, ":0.0;\n");
      else
	sprintf(treestr, ";\n");	 	  	
    }
  else 
    {                   
      if(rellTree || branchLabelSupport || printSHSupport)
	{	 	 
	  if(( !isTip(p->number, tr->mxtips)) && 
	     ( !isTip(p->back->number, tr->mxtips)))
	    {			      
	      assert(p->bInf != (branchInfo *)NULL);
	      
	      if(rellTree)
		sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	      if(branchLabelSupport)
		sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f[%d]", getBranchLength(tr, pr, perGene, p), p->bInf->support);
	      
	    }
	  else		
	    {
	      if(rellTree || branchLabelSupport)
		sprintf(treestr, ":%8.20f", p->z[0]);	
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f", getBranchLength(tr, pr, perGene, p));
	    }
	}
      else
	{
	  if(printBranchLengths)	    
	    sprintf(treestr, ":%8.20f", getBranchLength(tr, pr, perGene, p));
	  else	    
	    sprintf(treestr, "%s", "\0");	    
	}      
    }
  
  while (*treestr) treestr++;
  return  treestr;
}


char *pllTreeToNewick(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean printBranchLengths, pllBoolean printNames, pllBoolean printLikelihood,
		  pllBoolean rellTree, pllBoolean finalPrint, int perGene, pllBoolean branchLabelSupport, pllBoolean printSHSupport)
{ 

  if(rellTree)
    assert(!branchLabelSupport && !printSHSupport);

  if(branchLabelSupport)
    assert(!rellTree && !printSHSupport);

  if(printSHSupport)
    assert(!branchLabelSupport && !rellTree);

 
  pllTreeToNewickREC(treestr, tr, pr, p, printBranchLengths, printNames, printLikelihood, rellTree,
		 finalPrint, perGene, branchLabelSupport, printSHSupport);  
    
  
  while (*treestr) treestr++;
  
  return treestr;
}

