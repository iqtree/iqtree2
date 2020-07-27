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
 * @file searchAlgo.c
 * @brief Collection of routines for performing likelihood computation and branch optimization.
 *
 * Detailed description to appear soon.
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
#include <errno.h>

#include "pll.h"
#include "pllInternal.h"

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;

double treeOptimizeRapid(pllInstance *tr, partitionList *pr, int mintrav, int maxtrav, bestlist *bt, infoList *iList);
nniMove getBestNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, double curLH);
void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, nniMove* nniList, int* cnt, int* cnt_nni, double curLH);


static int cmp_nni(const void* nni1, const void* nni2);
static void pllTraverseUpdate (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav, pllRearrangeList * bestList);
static int pllStoreRearrangement (pllRearrangeList * bestList, pllRearrangeInfo * rearr);
static int pllTestInsertBIG (pllInstance * tr, partitionList * pr, nodeptr p, nodeptr q, pllRearrangeList * bestList);
static int pllTestSPR (pllInstance * tr, partitionList * pr, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList);
static void pllCreateSprInfoRollback (pllInstance * tr, pllRearrangeInfo * rearr, int numBranches);
static void pllCreateNniInfoRollback (pllInstance * tr, pllRearrangeInfo * rearr);
static void pllCreateRollbackInfo (pllInstance * tr, pllRearrangeInfo * rearr, int numBranches);
static void pllRollbackNNI (pllInstance * tr, partitionList * pr, pllRollbackInfo * ri);
static void pllRollbackSPR (partitionList * pr, pllRollbackInfo * ri);

extern partitionLengths pLengths[PLL_MAX_MODEL];

pllBoolean initrav (pllInstance *tr, partitionList *pr, nodeptr p)
{ 
  nodeptr  q;

  if (!isTip(p->number, tr->mxtips)) 
  {      
    q = p->next;

    do 
    {	   
      if (! initrav(tr, pr, q->back))  return PLL_FALSE;
      q = q->next;	
    } 
    while (q != p);  

    pllUpdatePartials(tr, pr, p, PLL_FALSE);
  }

  return PLL_TRUE;
} 


/** @brief Optimize the length of a specific branch

    Optimize the length of the branch connecting \a p and \a p->back
    for each partition (\a tr->numBranches) in library instance \a tr.
 
    @param tr
      The library instance

    @param pr
      Partition list
 
    @param p
      Endpoints of branch to be optimized 
*/
void update(pllInstance *tr, partitionList *pr, nodeptr p)
{       
  nodeptr  q; 
  int i;
  double   z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  #ifdef _DEBUG_UPDATE
    double 
      startLH;
  
    pllEvaluateLikelihood (tr, p);
  
    startLH = tr->likelihood;
  #endif

  q = p->back;   

  for(i = 0; i < numBranches; i++)
    z0[i] = q->z[i];    

  if(numBranches > 1)
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_TRUE);
  else
    makenewzGeneric(tr, pr, p, q, z0, PLL_NEWZPERCYCLE, z, PLL_FALSE);

  for(i = 0; i < numBranches; i++)
  {         
    if(!tr->partitionConverged[i])
    {	  
      if(PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ)  
      {	      
        tr->partitionSmoothed[i] = PLL_FALSE;
      }	 

      p->z[i] = q->z[i] = z[i];	 
    }
  }
 
  #ifdef _DEBUG_UPDATE
    pllEvaluateLikelihood (tr, p);
  
    if(tr->likelihood <= startLH)
      {
        if(fabs(tr->likelihood - startLH) > 0.01)
  	{
  	  printf("%f %f\n", startLH, tr->likelihood);
  	  assert(0);      
  	}
      }
  #endif
}

/** @brief Branch length optimization of subtree

    Optimize the length of branch connected by \a p and \a p->back, and the
    lengths of all branches in the subtrees rooted at \a p->next and \a p->next->next

    @param tr
      The library instance

    @param pr
      Partition list

    @param p
      Endpoint of branches to be optimized
*/
void smooth (pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  update(tr, pr, p);    /*  Adjust branch */

  if (! isTip(p->number, tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      smooth(tr, pr, q->back);
      q = q->next;
    }	

    if(numBranches > 1 && !tr->useRecom)
      pllUpdatePartials(tr, pr,p, PLL_TRUE);
    else
      pllUpdatePartials(tr, pr,p, PLL_FALSE);
  }
} 

/**  @brief Check whether the branches in all partitions have been optimized
 
     Check if all branches in all partitions have reached the threshold for
     optimization. If at least one branch can be optimized further return \b PLL_FALSE.

     @param tr
       The library instance 

     @return
       If at least one branch can be further optimized return \b PLL_FALSE,
       otherwise \b PLL_TRUE.
             
*/
static pllBoolean allSmoothed(pllInstance *tr, int numBranches)
{
  int i;
  pllBoolean result = PLL_TRUE;

  for(i = 0; i < numBranches; i++)
  {
    if(tr->partitionSmoothed[i] == PLL_FALSE)
      result = PLL_FALSE;
    else
      tr->partitionConverged[i] = PLL_TRUE;
  }

  return result;
}


/** @brief Optimize all branch lenghts of a tree
  
    Perform \a maxtimes rounds of branch length optimization by running smooth()
    on all neighbour nodes of node \a tr->start.

    @param tr
      The library instance

    @param maxtimes
      Number of optimization rounds to perform
*/
/* do maxtimes rounds of branch length optimization */
void smoothTree (pllInstance *tr, partitionList *pr, int maxtimes)
{
	nodeptr  p, q;
	int i, count = 0;
    int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

	p = tr->start;
	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = PLL_FALSE;

	while (--maxtimes >= 0)
	{
		for(i = 0; i < numBranches; i++)
			tr->partitionSmoothed[i] = PLL_TRUE;

		smooth(tr, pr, p->back);
		if (!isTip(p->number, tr->mxtips))
		{
			q = p->next;
			while (q != p)
			{
				smooth(tr, pr, q->back);
				q = q->next;
			}
		}
		count++;

		if (allSmoothed(tr, numBranches)) break;
	}

	for(i = 0; i < numBranches; i++)
		tr->partitionConverged[i] = PLL_FALSE;
} 


/** @brief Optimize the branch length of edges around a specific node
    
    Optimize \a maxtimes the branch length of all (3) edges around a given node 
    \a p of the tree of library instance \a tr.

    @param tr
      The library instance

    @param p
      The node around which to optimize the edges

    @param maxtimes
      Number of optimization rounds to perform
*/
void localSmooth (pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
  if (isTip(p->number, tr->mxtips)) return;

  for(i = 0; i < PLL_NUM_BRANCHES; i++)
    tr->partitionConverged[i] = PLL_FALSE;	

  while (--maxtimes >= 0) 
  {     
    for(i = 0; i < PLL_NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = PLL_TRUE;

    q = p;
    do 
    {
      update(tr, pr, q);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < PLL_NUM_BRANCHES; i++)
  {
    tr->partitionSmoothed[i] = PLL_FALSE; 
    tr->partitionConverged[i] = PLL_FALSE;
  }
}




/** @brief Reset an \a infoList

    Resets an \a infoList by setting elements \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a PLL_UNLIKELY, respectively.

    @param iList
      The given \a infoList.
*/
static void resetInfoList(infoList *iList)
{
  int 
    i;

  iList->valid = 0;

  for(i = 0; i < iList->n; i++)    
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = PLL_UNLIKELY;
  }    
}

/** @brief Initialize an \a infoList

    Initialize an \a infoList by creating a \a bestInfo list structure
    of \a n elements and setting the attributes \a node and \a likelihood
    of each element of the \a bestInfo list structure to \b NULL and
    \a PLL_UNLIKELY, respectively.

    @param iList
      The given \a infoList.

    @param n
      Number of elements to be created in the \a bestInfo list.
*/
static void initInfoList(infoList *iList, int n)
{
  int 
    i;

  iList->n = n;
  iList->valid = 0;
  iList->list = (bestInfo *)rax_malloc(sizeof(bestInfo) * (size_t)n);

  for(i = 0; i < n; i++)
  {
    iList->list[i].node = (nodeptr)NULL;
    iList->list[i].likelihood = PLL_UNLIKELY;
  }
}

/** @brief Deallocate the contents of an \a infoList
    
    Deallocate the contents of a given \a infoList by freeing
    the memory used by its \a bestInfo list structure.

    @param iList
      The \a infoList to be used.
*/
static void freeInfoList(infoList *iList)
{ 
  rax_free(iList->list);   
}


/** @brief Insert a record in an \a infoList

    Insert the pair \a likelihood and \node into list \a iList 
    \b only if there already exists a pair in \a iList 
    whose \a likelihood attribute is smaller than the given \a 
    likelihood. The insertion is done by replacing the smallest
    likelihood pair with the new pair.

    @param node
      The given node

    @param likelihood
      The given likelihood

    @param iList
      The given \a infoList where the record will possibly be appended.
*/
static void insertInfoList(nodeptr node, double likelihood, infoList *iList)
{
  int 
    i,
    min = 0;

  double 
    min_l =  iList->list[0].likelihood;

  for(i = 1; i < iList->n; i++)
  {
    if(iList->list[i].likelihood < min_l)
    {
      min = i;
      min_l = iList->list[i].likelihood;
    }
  }

  if(likelihood > min_l)
  {
    iList->list[min].likelihood = likelihood;
    iList->list[min].node = node;
    if(iList->valid < iList->n)
      iList->valid += 1;
  }
}


/** @brief  Optimize branch lengths of region

    Optimize the branch lenghts of only a specific region. The branch optimization starts
    at a node \a p and is carried out in all nodes with distance upto \a region edges from 
    \a p.

    @param tr
      The library instance.
    
    @param p
      Node to start branch optimization from.

    @param region
      The allowed node distance from \p for which to still perform branch optimization.
*/
void smoothRegion (pllInstance *tr, partitionList *pr, nodeptr p, int region)
{ 
  nodeptr  q;

  update(tr, pr, p);   /* Adjust branch */

  if (region > 0)
  {
    if (!isTip(p->number, tr->mxtips)) 
    {                                 
      q = p->next;
      while (q != p) 
      {
        smoothRegion(tr, pr, q->back, --region);
        q = q->next;
      }	

      pllUpdatePartials(tr, pr,p, PLL_FALSE);
    }
  }
}

/** @brief Wrapper function for optimizing the branch length of a region \a maxtimes times

    Optimize the branch lengths of a specific region \a maxtimes times. The branch optimization
    starts at a given node \a p and is carried out in all nodes with distance upto \a region
    from \a p.

    @param tr
      The library instance.

    @param p
      Node to start branch optimization from.

    @param maxtimes
      Number of times to perform branch optimization.

    @param region
      The allwed node distance from \p for which to still perform branch optimization.

    @todo
      In the previous version (before the model-sep merge) the loops were controlled by tr->numBranches,
      and now they are controlled by a constant PLL_NUM_BRANCHES. What is right?
*/
void regionalSmooth (pllInstance *tr, partitionList *pr, nodeptr p, int maxtimes, int region)
{
  nodeptr  q;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if (isTip(p->number, tr->mxtips)) return;            /* Should be an error */

  for(i = 0; i < PLL_NUM_BRANCHES; i++)
    tr->partitionConverged[i] = PLL_FALSE;

  while (--maxtimes >= 0) 
  {	
    for(i = 0; i < PLL_NUM_BRANCHES; i++)
      tr->partitionSmoothed[i] = PLL_TRUE;

    q = p;
    do 
    {
      smoothRegion(tr, pr, q, region);
      q = q->next;
    } 
    while (q != p);

    if (allSmoothed(tr, numBranches))
      break;
  }

  for(i = 0; i < PLL_NUM_BRANCHES; i++) {
    tr->partitionSmoothed[i] = PLL_FALSE;
    tr->partitionConverged[i] = PLL_FALSE;
  }
} 




/** @brief Split the tree into two components and optimize new branch length

   Split the tree into two components. The disconnection point is node \a p.
   First, a branch length is computed for the newly created branch between nodes
   \a p->next->back and \a p->next->next->back and then the two nodes are
   connected (hookup). Disconnection is done by setting \a p->next->next->back
   and \a p->next->back to \b NULL.

   @param tr
     The library instance

   @param p
     The node at which the tree should be decomposed into two components.

   @param numBranches
     Number of branches per partition

   @return
     Node from the disconnected component

   @todo
     Why do we return this node?

   @image html removeBIG.png "The diagram shows in blue color the new edge that is created and in red the edges that are removed"
*/
nodeptr  removeNodeBIG (pllInstance *tr, partitionList *pr, nodeptr p, int numBranches)
{  
//  double   zqr[numBranches], result[numBranches];
    double* zqr = (double*)rax_malloc(numBranches * sizeof(double));
    double* result = (double*)rax_malloc(numBranches * sizeof(double));
  nodeptr  q, r;
  int i;

  q = p->next->back;
  r = p->next->next->back;

  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        

  makenewzGeneric(tr, pr, q, r, zqr, PLL_ITERATIONS, result, PLL_FALSE);

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 

  p->next->next->back = p->next->back = (node *) NULL;

  rax_free(result);
  rax_free(zqr);
  return  q; 
}

/** @brief Split the tree into two components and recompute likelihood

    Split the tree into two component. The disconnection point is node \a p.
    Set the branch length of the new node between \a p->next->back and
    \a p->next->next->back to \a tr->currentZQR and then decompose the tree
    into two components by setting \a p->next->back and \a p->next->next->back
    to \b NULL.

    @param tr
      The library instance

    @param p
      The node at which the tree should be decomposed into two components.

    @return q
      the node after \a p

    @todo
      Why do we return this node? Why do we set to tr->currentZQR and not compute
      new optimized length? What is tr->currentZQR? 
*/
nodeptr  removeNodeRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q, r;

  q = p->next->back;
  r = p->next->next->back;  

  pllUpdatePartials(tr, pr,q, PLL_FALSE);
  pllUpdatePartials(tr, pr,r, PLL_FALSE);

  hookup(q, r, tr->currentZQR, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  p->next->next->back = p->next->back = (node *) NULL;

  return  q;
}

/** @brief Connect two disconnected tree components
   
   Connect two disconnected components by specifying an internal edge from one
   component and a leaf from the other component. The internal edge \a e is the
   edge between \a q and \a q->back. The leaf is specified by \a p.
   Edge \a e is removed and two new edges are created. The first one is an edge
   between \a p->next and \a q, and the second one is between \a p->next->next
   and \a q->back. The new likelihood vector for node \a p is computed.

   @note The function makes use of the \a thoroughInsertion flag

   @todo
     What is tr->lzi ? What is thorough insertion? Why do we optimize branch lengths
     that will be removed? Add explanation

   @image html pll.png "The diagram shows in blue colors the new edges that are created and in red the edge that is removed" 
*/
pllBoolean insertBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  int i;
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  r = q->back;
  s = p->back;

  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];

  if(tr->thoroughInsertion)
  { 
      double* zqr = (double*)rax_malloc(numBranches * sizeof(double));
      double* zqs = (double*)rax_malloc(numBranches * sizeof(double));
      double* zrs = (double*)rax_malloc(numBranches * sizeof(double));
	  double lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax; 
      double *defaultArray = (double*)rax_malloc(numBranches*sizeof(double));
      double* e1 = (double*)rax_malloc(numBranches * sizeof(double));
      double* e2 = (double*)rax_malloc(numBranches * sizeof(double));
      double* e3 = (double*)rax_malloc(numBranches*sizeof(double));
    double *qz;

    qz = q->z;

    for(i = 0; i < numBranches; i++)
      defaultArray[i] = PLL_DEFAULTZ;

    makenewzGeneric(tr, pr, q, r, qz, PLL_ITERATIONS, zqr, PLL_FALSE);
    /* the branch lengths values will be estimated using q, r and s
     * q-s are not connected, but both q and s have a valid LH vector , so we can call makenewzGeneric  to get a value for
     * lzsum, which is then use to generate reasonable starting values e1, e2, e3 for the new branches we create after the       insertion
     */

    makenewzGeneric(tr, pr, q, s, defaultArray, PLL_ITERATIONS, zqs, PLL_FALSE);
    makenewzGeneric(tr, pr, r, s, defaultArray, PLL_ITERATIONS, zrs, PLL_FALSE);


    for(i = 0; i < numBranches; i++)
    {
      lzqr = (zqr[i] > PLL_ZMIN) ? log(zqr[i]) : log(PLL_ZMIN); 
      lzqs = (zqs[i] > PLL_ZMIN) ? log(zqs[i]) : log(PLL_ZMIN);
      lzrs = (zrs[i] > PLL_ZMIN) ? log(zrs[i]) : log(PLL_ZMIN);
      lzsum = 0.5 * (lzqr + lzqs + lzrs);

      lzq = lzsum - lzrs;
      lzr = lzsum - lzqs;
      lzs = lzsum - lzqr;
      lzmax = log(PLL_ZMAX);

      if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
      else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
      else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          

      e1[i] = exp(lzq);
      e2[i] = exp(lzr);
      e3[i] = exp(lzs);
    }
    hookup(p->next,       q, e1, numBranches);
    hookup(p->next->next, r, e2, numBranches);
    hookup(p,             s, e3, numBranches);      		  
	rax_free(e3);
	rax_free(e2);
	rax_free(e1);
	rax_free(defaultArray);
	rax_free(zrs);
	rax_free(zqs);
	rax_free(zqr);

  }
  else
  {       
	  double  *z = (double*) rax_malloc(numBranches*sizeof(double));

    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      

      if(z[i] < PLL_ZMIN) 
        z[i] = PLL_ZMIN;
      if(z[i] > PLL_ZMAX)
        z[i] = PLL_ZMAX;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
	rax_free(z);
  }

  pllUpdatePartials(tr, pr,p, PLL_FALSE);

  if(tr->thoroughInsertion)
  {     
    localSmooth(tr, pr, p, PLL_MAX_LOCAL_SMOOTHING_ITERATIONS);
    for(i = 0; i < numBranches; i++)
    {
      tr->lzq[i] = p->next->z[i];
      tr->lzr[i] = p->next->next->z[i];
      tr->lzs[i] = p->z[i];            
    }
  }           

  return  PLL_TRUE;
}

/** @brief Connect two disconnected tree components without optimizing branch lengths
   
   Connect two disconnected components by specifying an internal edge from one
   component and a leaf from the other component. The internal edge \a e is the
   edge between \a q and \a q->back. The leaf is specified by \a p.
   Edge \a e is removed and two new edges are created. The first one is an edge
   between \a p->next and \a q, and the second one is between \a p->next->next
   and \a q->back. The new likelihood vector for node \a p is computed.

   @note The function makes use of the \a thoroughInsertion flag

   @todo
     What is the difference between this and insertBIG? 
*/
pllBoolean insertRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;

  r = q->back;
  s = p->back;

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  {       
    double  z[PLL_NUM_BRANCHES];
    int i;

    for(i = 0; i < numBranches; i++)
    {
      double zz;
      zz = sqrt(q->z[i]);     
      if(zz < PLL_ZMIN) 
        zz = PLL_ZMIN;
      if(zz > PLL_ZMAX)
        zz = PLL_ZMAX;
      z[i] = zz;
    }

    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }   

  pllUpdatePartials(tr, pr,p, PLL_FALSE);

  return  PLL_TRUE;
}


static void restoreTopologyOnly(pllInstance *tr, bestlist *bt, int numBranches)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[PLL_NUM_BRANCHES], pz[PLL_NUM_BRANCHES], p1z[PLL_NUM_BRANCHES], p2z[PLL_NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;

  p1 = p->next->back;
  p2 = p->next->next->back;

  //memcpy(p1z, p1->z, numBranches*sizeof(double));
  //memcpy(p2z, p2->z, numBranches*sizeof(double));
  //memcpy(qz, q->z, numBranches*sizeof(double));
  //memcpy(pz, p->z, numBranches*sizeof(double));
  for(i = 0; i < numBranches; i++)
  {
    p1z[i] = p1->z[i];
    p2z[i] = p2->z[i];
  }

  hookup(p1, p2, tr->currentZQR, numBranches);

  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }

  r = q->back;
  s = p->back;

  if(tr->thoroughInsertion)
  {                        
    hookup(p->next,       q, tr->currentLZQ, numBranches);
    hookup(p->next->next, r, tr->currentLZR, numBranches);
    hookup(p,             s, tr->currentLZS, numBranches);
  }
  else
  { 	
    double  z[PLL_NUM_BRANCHES];	
    for(i = 0; i < numBranches; i++)
    {
      z[i] = sqrt(q->z[i]);      
      if(z[i] < PLL_ZMIN)
        z[i] = PLL_ZMIN;
      if(z[i] > PLL_ZMAX)
        z[i] = PLL_ZMAX;
    }
    hookup(p->next,       q, z, numBranches);
    hookup(p->next->next, r, z, numBranches);
  }     

  tr->likelihood = tr->bestOfNode;

  saveBestTree(bt, tr, numBranches);

  tr->likelihood = currentLH;

  hookup(q, r, qz, numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)    
    hookup(p, s, pz, numBranches);

  hookup(p->next,       p1, p1z, numBranches);
  hookup(p->next->next, p2, p2z, numBranches);
}

/** @brief Test the 
*/
pllBoolean testInsertBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  double  qz[PLL_NUM_BRANCHES], pz[PLL_NUM_BRANCHES];
  nodeptr  r;
  double startLH = tr->endLH;
  int i;

  r = q->back; 
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }

  if (! insertBIG(tr, pr, p, q))       return PLL_FALSE;

  pllEvaluateLikelihood (tr, pr, p->next->next, PLL_FALSE, PLL_FALSE);

  if(tr->likelihood > tr->bestOfNode)
  {
    tr->bestOfNode = tr->likelihood;
    tr->insertNode = q;
    tr->removeNode = p;   
    for(i = 0; i < numBranches; i++)
    {
      tr->currentZQR[i] = tr->zqr[i];           
      tr->currentLZR[i] = tr->lzr[i];
      tr->currentLZQ[i] = tr->lzq[i];
      tr->currentLZS[i] = tr->lzs[i];      
    }
  }

  if(tr->likelihood > tr->endLH)
  {			  
    tr->insertNode = q;
    tr->removeNode = p;   
    for(i = 0; i < numBranches; i++)
      tr->currentZQR[i] = tr->zqr[i];      
    tr->endLH = tr->likelihood;                      
  }        

  /* reset the topology so that it is the same as it was before calling insertBIG */
  hookup(q, r, qz, numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)
  {
    nodeptr s = p->back;
    hookup(p, s, pz, numBranches);
  } 

  if((tr->doCutoff) && (tr->likelihood < startLH))
  {
    tr->lhAVG += (startLH - tr->likelihood);
    tr->lhDEC++;
    if((startLH - tr->likelihood) >= tr->lhCutoff)
      return PLL_FALSE;	    
    else
      return PLL_TRUE;
  }
  else
    return PLL_TRUE;
}


/** @brief Recursively traverse tree and test insertion

    Recursively traverses the tree structure starting from node \a q and
    tests the insertion of the component specified by leaf \a p at the edge
    between \a q and \a q->back.

    @param tr
      PLL instance

    @param pr
      List of partitions
    @param p
      Leaf node of one tree component

    @param q
      Endpoint node of the edge to test the insertion

    @param mintrav
      Minimum radius around \a q to test the insertion

    @param maxtrav
      Maximum radius around \a q to test the insertion\
*/
void addTraverseBIG(pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
  {              
    if (! testInsertBIG(tr, pr, p, q))  return;

  }

  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
  {    
    addTraverseBIG(tr, pr, p, q->next->back, mintrav, maxtrav);
    addTraverseBIG(tr, pr, p, q->next->next->back, mintrav, maxtrav);
  }
} 




/** @brief  Compute the  best SPR movement

    Compute all SPR moves starting from \a p in the space defined by \a mintrav and
    \a maxtrav and store the best in the \a tr structure.

    @param tr
      PLL instancve

    @param pr
      List of partitions

    @param p
      Node from which to start the SPR moves testing

    @param mintrav
      Minimum distance from \a p where to start testing SPRs

    @param maxtrav
      Maximum distance from \a p where to test SPRs

    @return
       0,1 or \b PLL_BADREAR

    @todo
      fix the return value
*/
int rearrangeBIG(pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav)
{  
  double   p1z[PLL_NUM_BRANCHES], p2z[PLL_NUM_BRANCHES], q1z[PLL_NUM_BRANCHES], q2z[PLL_NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  pllBoolean doP = PLL_TRUE, doQ = PLL_TRUE;
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  if (maxtrav < 1 || mintrav > maxtrav)  return (0);
  q = p->back;




  if (!isTip(p->number, tr->mxtips) && doP) 
  {     
    p1 = p->next->back;
    p2 = p->next->next->back;


    if(!isTip(p1->number, tr->mxtips) || !isTip(p2->number, tr->mxtips))
    {
      for(i = 0; i < numBranches; i++)
      {
        p1z[i] = p1->z[i];
        p2z[i] = p2->z[i];	   	   
      }

      if (! removeNodeBIG(tr, pr, p,  numBranches)) return PLL_BADREAR;

      if (!isTip(p1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p1->next->back,
            mintrav, maxtrav);         

        addTraverseBIG(tr, pr, p, p1->next->next->back,
            mintrav, maxtrav);          
      }

      if (!isTip(p2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, p, p2->next->back,
            mintrav, maxtrav);
        addTraverseBIG(tr, pr, p, p2->next->next->back,
            mintrav, maxtrav);          
      }

      hookup(p->next,       p1, p1z, numBranches);
      hookup(p->next->next, p2, p2z, numBranches);
      pllUpdatePartials(tr, pr,p, PLL_FALSE);
    }
  }  

  if (!isTip(q->number, tr->mxtips) && maxtrav > 0 && doQ) 
  {
    q1 = q->next->back;
    q2 = q->next->next->back;

    /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
      ((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
    if (
        (
         ! isTip(q1->number, tr->mxtips) && 
         (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
        )
        ||
        (
         ! isTip(q2->number, tr->mxtips) && 
         (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
        )
       )
    {

      for(i = 0; i < numBranches; i++)
      {
        q1z[i] = q1->z[i];
        q2z[i] = q2->z[i];
      }

      if (! removeNodeBIG(tr, pr, q, numBranches)) return PLL_BADREAR;

      mintrav2 = mintrav > 2 ? mintrav : 2;

      if (/*! q1->tip*/ !isTip(q1->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q1->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q1->next->next->back,
            mintrav2 , maxtrav);         
      }

      if (/*! q2->tip*/ ! isTip(q2->number, tr->mxtips)) 
      {
        addTraverseBIG(tr, pr, q, q2->next->back,
            mintrav2 , maxtrav);
        addTraverseBIG(tr, pr, q, q2->next->next->back,
            mintrav2 , maxtrav);          
      }	   

      hookup(q->next,       q1, q1z, numBranches);
      hookup(q->next->next, q2, q2z, numBranches);

      pllUpdatePartials(tr, pr,q, PLL_FALSE);
    }
  } 

  return  1;
} 




/** @brief Perform an SPR move?

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param mintrav

    @param maxtrav

    @param adef

    @param bt

    @param iList

*/
double treeOptimizeRapid(pllInstance *tr, partitionList *pr, int mintrav, int maxtrav, bestlist *bt, infoList *iList)
{
  int i, index, *perm = (int*)NULL;   

  nodeRectifier(tr);

  if (maxtrav > tr->mxtips - 3) {
      maxtrav = tr->mxtips - 3;
  }

  resetInfoList(iList);

  resetBestTree(bt);

  tr->startLH = tr->endLH = tr->likelihood;

  if(tr->doCutoff)
  {
    if(tr->bigCutoff)
    {	  
      if(tr->itCount == 0)    
        tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
      else    		 
        tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
    }
    else
    {
      if(tr->itCount == 0)    
        tr->lhCutoff = tr->likelihood / -1000.0;    
      else    		 
        tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
    }    

    tr->itCount = tr->itCount + 1;
    tr->lhAVG = 0;
    tr->lhDEC = 0;
  }

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
  {           
    tr->bestOfNode = PLL_UNLIKELY;
    //James B. Was doing a null de-reference of perm here, 
    //if tr->permuteTreeoptimize was non-zero.  No longer!
    index = i;     

    if(rearrangeBIG(tr, pr, tr->nodep[index], mintrav, maxtrav))
    {    
      if(tr->thoroughInsertion)
      {
        if(tr->endLH > tr->startLH)                 	
        {			   
          /* commit the best SPR found by rearrangeBIG */
          restoreTreeFast(tr, pr);    
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 		  
          if(tr->bestOfNode != PLL_UNLIKELY)
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }	   
      }
      else
      {
        insertInfoList(tr->nodep[index], tr->bestOfNode, iList);	    
        if(tr->endLH > tr->startLH)                 	
        {		      
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
        }	    	  
      }
    }     
  }     

  if(!tr->thoroughInsertion)
  {           
    tr->thoroughInsertion = PLL_TRUE;  

    for(i = 0; i < iList->valid; i++)
    { 	  
      tr->bestOfNode = PLL_UNLIKELY;

      if(rearrangeBIG(tr, pr, iList->list[i].node, mintrav, maxtrav))
      {	  
        if(tr->endLH > tr->startLH)                 	
        {	 	     
          restoreTreeFast(tr, pr);
          tr->startLH = tr->endLH = tr->likelihood;	 
          saveBestTree(bt, tr, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
        }
        else
        { 

          if(tr->bestOfNode != PLL_UNLIKELY)
          {	     
            restoreTopologyOnly(tr, bt, pr->perGeneBranchLengths?pr->numberOfPartitions:1);
          }	
        }      
      }
    }       

    tr->thoroughInsertion = PLL_FALSE;
  }

  if(tr->permuteTreeoptimize)
    rax_free(perm);

  return tr->startLH;     
}




pllBoolean testInsertRestoreBIG (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{    
  if(tr->thoroughInsertion)
  {
    if (! insertBIG(tr, pr, p, q))       return PLL_FALSE;

    pllEvaluateLikelihood (tr, pr, p->next->next, PLL_FALSE, PLL_FALSE);
  }
  else
  {
    if (! insertRestoreBIG(tr, pr, p, q))       return PLL_FALSE;

    {
      nodeptr x, y;
      x = p->next->next;
      y = p->back;

      if(! isTip(x->number, tr->mxtips) && isTip(y->number, tr->mxtips))
      {
        while ((! x->x)) 
        {
          if (! (x->x))
            pllUpdatePartials(tr, pr,x, PLL_FALSE);
        }
      }

      if(isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! y->x)) 
        {		  
          if (! (y->x))
            pllUpdatePartials(tr, pr,y, PLL_FALSE);
        }
      }

      if(!isTip(x->number, tr->mxtips) && !isTip(y->number, tr->mxtips))
      {
        while ((! x->x) || (! y->x)) 
        {
          if (! (x->x))
            pllUpdatePartials(tr, pr,x, PLL_FALSE);
          if (! (y->x))
            pllUpdatePartials(tr, pr,y, PLL_FALSE);
        }
      }				      	

    }

    tr->likelihood = tr->endLH;
  }

  return PLL_TRUE;
} 

void restoreTreeFast(pllInstance *tr, partitionList *pr)
{
  removeNodeRestoreBIG(tr, pr, tr->removeNode);
  testInsertRestoreBIG(tr, pr, tr->removeNode, tr->insertNode);
}

/*
static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written == nmemb);
}

static void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, stream);

  assert(bytes_read == nmemb);
}

static void readTree(pllInstance *tr, partitionList *pr, FILE *f)
{
  int 
    nodeNumber,   
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    startAddress;

  myfread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];


  myfread(&startAddress, sizeof(nodeptr), 1, f);

  myfread(tr->nodeBaseAddress, sizeof(node), x, f);

  {
    int i;    

    size_t         
      offset;

    pllBoolean 
      addIt;

    if(startAddress > tr->nodeBaseAddress)
    {
      addIt = PLL_FALSE;
      offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
    }
    else
    {
      addIt = PLL_TRUE;
      offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
    }       

    for(i = 0; i < x; i++)
    {      	
      if(addIt)
      {	    
        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
      }
      else
      {

        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);	
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);	   
      } 
    }

  }

  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

  printBothOpen("RAxML Restart with likelihood: %1.50f\n", tr->likelihood);
}

static void readCheckpoint(pllInstance *tr, partitionList *pr)
{
  int  
    restartErrors = 0,
                  model; 

  FILE 
    *f = myfopen(binaryCheckpointInputName, "r");
*/
  /* cdta */   
/*
  myfread(&(tr->ckp), sizeof(checkPointState), 1, f);



  if(tr->ckp.searchConvergenceCriterion != tr->searchConvergenceCriterion)
  {
    printf("restart error, you are trying to re-start a run where the ML search criterion was turned %s\n", (tr->ckp.searchConvergenceCriterion)?"ON":"OFF");
    restartErrors++;
  }  

  if(tr->ckp.rateHetModel !=  tr->rateHetModel)
  {
    printf("restart error, you are trying to re-start a run with a different model of rate heterogeneity, the checkpoint was obtained under: %s\n", (tr->ckp.rateHetModel == PLL_GAMMA)?"GAMMA":"PSR");
    restartErrors++;
  }  

  if(tr->ckp.maxCategories !=  tr->maxCategories)
  {
    printf("restart error, you are trying to re-start a run with %d per-site rate categories, the checkpoint was obtained with: %d\n", tr->maxCategories, tr->ckp.maxCategories);
    restartErrors++;
  }

  if(tr->ckp.NumberOfModels != pr->numberOfPartitions)
  {
    printf("restart error, you are trying to re-start a run with %d partitions, the checkpoint was obtained with: %d partitions\n", (int)pr->numberOfPartitions, tr->ckp.NumberOfModels);
    restartErrors++;      
  }

  if(tr->ckp.numBranches != pr->perGeneBranchLengths?pr->numberOfPartitions:1)
  {
    printf("restart error, you are trying to re-start a run where independent per-site branch length estimates were turned %s\n", (tr->ckp.numBranches > 1)?"ON":"OFF");
    restartErrors++;
  }

  if(tr->ckp.originalCrunchedLength != tr->originalCrunchedLength)
  {
    printf("restart error, you are trying to re-start a run with %d site patterns, the checkpoint was obtained with: %d site patterns\n", tr->ckp.originalCrunchedLength, tr->originalCrunchedLength);
    restartErrors++; 
  }

  if(tr->ckp.mxtips != tr->mxtips)
  {
    printf("restart error, you are trying to re-start a run with %d taxa, the checkpoint was obtained with: %d taxa\n", tr->mxtips, tr->ckp.mxtips);
    restartErrors++; 
  }

  if(strcmp(tr->ckp.seq_file, seq_file) != 0)
  {
    printf("restart error, you are trying to re-start from alignemnt file %s, the checkpoint was obtained with file: %s\n", tr->ckp.seq_file, seq_file);
    restartErrors++; 
  }

  printf("REstart errors: %d\n", restartErrors);

  if(restartErrors > 0)
  {
    printf("User induced errors with the restart from checkpoint, exiting ...\n");

    if(restartErrors > 4)
      printf(" ... maybe you should do field work instead of trying to use a computer ...\n");
    if(restartErrors > 6)
      printf(" ... kala eisai telios ilithios;\n");

    exit(-1);
  }

  tr->ntips = tr->mxtips;

  tr->startLH    = tr->ckp.tr_startLH;
  tr->endLH      = tr->ckp.tr_endLH;
  tr->likelihood = tr->ckp.tr_likelihood;
  tr->bestOfNode = tr->ckp.tr_bestOfNode;

  tr->lhCutoff   = tr->ckp.tr_lhCutoff;
  tr->lhAVG      = tr->ckp.tr_lhAVG;
  tr->lhDEC      = tr->ckp.tr_lhDEC;
  tr->itCount    = tr->ckp.tr_itCount;
  tr->thoroughInsertion       = tr->ckp.tr_thoroughInsertion;



  accumulatedTime = tr->ckp.accumulatedTime;
*/
  /* printf("Accumulated time so far: %f\n", accumulatedTime); */
/*
  tr->optimizeRateCategoryInvocations = tr->ckp.tr_optimizeRateCategoryInvocations;


  myfread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  if(tr->searchConvergenceCriterion)
  {
    int bCounter = 0;

    if((tr->ckp.state == PLL_FAST_SPRS && tr->ckp.fastIterations > 0) ||
        (tr->ckp.state == PLL_SLOW_SPRS && tr->ckp.thoroughIterations > 0))
    { 

#ifdef _DEBUG_CHECKPOINTING    
      printf("parsing Tree 0\n");
#endif

      treeReadTopologyString(tr->tree0, tr);   

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 0, PLL_BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }

    bCounter = 0;

    if((tr->ckp.state == PLL_FAST_SPRS && tr->ckp.fastIterations > 1) ||
        (tr->ckp.state == PLL_SLOW_SPRS && tr->ckp.thoroughIterations > 1))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 1\n");
#endif

      treeReadTopologyString(tr->tree1, tr); 

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->h, 1, PLL_BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, PLL_FALSE, PLL_FALSE, tr->threadID);

      assert(bCounter == tr->mxtips - 3);
    }
  }

  myfread(tr->rateCategory, sizeof(int), tr->originalCrunchedLength, f);
  myfread(tr->patrat, sizeof(double), tr->originalCrunchedLength, f);
  myfread(tr->patratStored, sizeof(double), tr->originalCrunchedLength, f);

*/
  /* need to read this as well in checkpoints, otherwise the branch lengths 
     in the output tree files will be wrong, not the internal branch lengths though */
/*
  //TODO: Same problem as writing the checkpoint
  //myfread(tr->fracchanges,  sizeof(double), pr->numberOfPartitions, f);
  myfread(&(tr->fracchange),   sizeof(double), 1, f);
*/
  /* pInfo */
/*
  for(model = 0; model < pr->numberOfPartitions; model++)
  {
    int 
      dataType = pr->partitionData[model]->dataType;

    myfread(&(pr->partitionData[model]->numberOfCategories), sizeof(int), 1, f);
    myfread(pr->partitionData[model]->perSiteRates, sizeof(double), tr->maxCategories, f);
    myfread(pr->partitionData[model]->EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfread(pr->partitionData[model]->EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfread(pr->partitionData[model]->EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfread(pr->partitionData[model]->frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfread(pr->partitionData[model]->tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfread(pr->partitionData[model]->substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfread(&(pr->partitionData[model]->alpha), sizeof(double), 1, f);
    
    if(pr->partitionData[model]->protModels == PLL_LG4M || pr->partitionData[model]->protModels == PLL_LG4X)
	{
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {
	      myfread(pr->partitionData[model]->EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	      myfread(pr->partitionData[model]->EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	      myfread(pr->partitionData[model]->EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	      myfread(pr->partitionData[model]->frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	      myfread(pr->partitionData[model]->tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	      myfread(pr->partitionData[model]->substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	    }
	}

    pllMakeGammaCats(pr->partitionData[model]->alpha, pr->partitionData[model]->gammaRates, 4, tr->useMedian);
  }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier (tr, pr, PLL_THREAD_COPY_INIT_MODEL);
#endif

  updatePerSiteRates(tr, pr, PLL_FALSE);

  readTree(tr, pr, f);

  fclose(f); 

}

void restart(pllInstance *tr, partitionList *pr)
{  
  readCheckpoint(tr, pr);

  switch(tr->ckp.state)
  {
    case PLL_REARR_SETTING:      
      break;
    case PLL_FAST_SPRS:
      break;
    case PLL_SLOW_SPRS:
      break;
    default:
      assert(0);
  }
}
*/

/* The number of maximum smoothing iterations is given explicitely */
/** @brief Optimize branch lenghts and evaluate likelihood of topology
    
    Optimize the branch lengths \a maxSmoothIterations times and evaluate
    the likelihood of tree. The resulting likelihood is placed in
    \a tr->likelihood

    @param tr
      The PLL instance

    @param pr
      List of partitions

    @param maxSmoothIterations
      Number of times to optimize branch lengths
*/
void
pllOptimizeBranchLengths (pllInstance *tr, partitionList *pr, int maxSmoothIterations)       /* Evaluate a user tree */
{
  smoothTree(tr, pr, maxSmoothIterations); /* former (32 * smoothFactor) */

  pllEvaluateLikelihood (tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
}

/** @brief Perform an NNI move

    Modify the tree topology of instance \a tr by performing an NNI (Neighbour Neighbor
    Interchange) move at node \a p. Let \a q be \a p->back. If \a swap is set to \b PLL_NNI_P_NEXT 
    then the subtrees rooted at \a p->next->back and \a q->next->back will be swapped. Otherwise,
    if \a swap is set to \b PLL_NNI_P_NEXTNEXT then the subtrees rooted at \a p->next->next->back and
    \a q->next->back are swapped. For clarity, see the illustration.

    @param tr
      PLL instance

    @param p
      Node to use as origin for performing NNI

    @param swap
      Which node to use for the NNI move. \b PLL_NNI_P_NEXT uses node p->next while \b PLL_NNI_P_NEXTNEXT uses p->next->next

    @return
      In case of success \b PLL_TRUE, otherwise \b PLL_FALSE

    @todo
      Started error checking here. Instead of checking the errors in the specified way, implement a variadic
      function where we pass the results of each check and the error code we want to assign if there is at
      least one negative result

    @image html nni.png "In case \a swap is set to \b PLL_NNI_P_NEXT then the dashed red edge between \a p and \a r is removed and the blue edges are created. If \a swap is set to \b PLL_INIT_P_NEXTNEXT then the dashed red edge between \a p and \a s is removed and the green edges are created. In both cases the black dashed edge is removed"
*/
int pllTopologyPerformNNI(pllInstance * tr, nodeptr p, int swap)
{
  nodeptr       q, r;

  q = p->back;
  if (isTip(q->number, tr->mxtips))
   {
     errno = PLL_NNI_Q_TIP;
     return (PLL_FALSE);
   }
  if (isTip(p->number, tr->mxtips))
   {
     errno = PLL_NNI_P_TIP;
     return (PLL_FALSE);
   }
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == PLL_NNI_P_NEXT)
   {
     r = p->next->back;
     hookupFull(p->next, q->next->back, q->next->z);
     hookupFull(q->next, r,             p->next->z);
   }
  else
   {
     r = p->next->next->back;
     hookupFull(p->next->next, q->next->back, q->next->z);
     hookupFull(q->next,       r,             p->next->next->z);
   }

  return PLL_TRUE;
}

/** @brief Compares 2 NNI moves */
static int cmp_nni(const void* nni1, const void* nni2) {
	nniMove* myNNI1 = (nniMove*) nni1;
	nniMove* myNNI2 = (nniMove*) nni2;
	return (int) (1000000.f * myNNI1->deltaLH - 1000000.f * myNNI2->deltaLH);
}

/** @brief Gets the best NNI move for a branch

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Node to use as origin for performing NNI

    @param curLH
      The current likelihood

    @return
      The best NNI move

*/
nniMove getBestNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p,
		double curLH) {
	nodeptr q = p->back;
	assert( ! isTip(p->number, tr->mxtips));
	assert( ! isTip(q->number, tr->mxtips));
#ifdef _DEBUG_NNI
	pllTreeToNewick(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
#endif

	/* Backup the current branch length */
	double z0[PLL_NUM_BRANCHES];
	int i;
	for (i = 0; i < pr->numberOfPartitions; i++) {
		z0[i] = p->z[i];
	}
#ifdef _DEBUG_NNI
	double lhOld = tr->likelihood;
	printf("lhOld: %f \n", lhOld);
#endif
	double lh0 = curLH;


#ifdef _DEBUG_NNI
	printf("lh0: %f \n", lh0);
#endif
	nniMove nni0; // nni0 means no NNI move is done
	nni0.p = p;
	nni0.nniType = 0;
	nni0.deltaLH = 0;
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni0.z[i] = p->z[i];
	}

	/* Save the scaling factor */
	// Now try to do an NNI move of type 1
	pllTopologyPerformNNI(tr, p, PLL_NNI_P_NEXT);
	double lh1 = tr->likelihood;
	/* Update branch lengths */
	pllUpdatePartials(tr, pr, p, PLL_FALSE);
	pllUpdatePartials(tr, pr, q, PLL_FALSE);
	update(tr, pr, p);
	pllEvaluateLikelihood (tr, pr, p, PLL_FALSE, PLL_FALSE);

	nniMove nni1;
	nni1.p = p;
	nni1.nniType = 1;
	// Store the optimized und unoptimized central branch length
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni1.z[i] = p->z[i];
		nni1.z0[i] = z0[i];
	}
	nni1.likelihood = lh1;
	nni1.deltaLH = lh1 - lh0;
#ifdef _DEBUG_NNI
	printf("Delta likelihood of the 1.NNI move: %f\n", nni1.deltaLH);
#endif

	/* Restore previous NNI move */
	pllTopologyPerformNNI(tr, p, PLL_NNI_P_NEXT);
	/* Restore the old branch length */
	for (i = 0; i < pr->numberOfPartitions; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}

#ifdef _DEBUG_NNI
	printf("Restore topology\n");
	pllTreeToNewick(tr->tree_string, tr, tr->start->back, TRUE, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
	fprintf(stderr, "%s\n", tr->tree_string);
	pllEvaluateLikelihood (tr, tr->start, TRUE);
	printf("Likelihood after restoring from NNI 1: %f\n", tr->likelihood);
#endif
	/* Try to do an NNI move of type 2 */
	pllTopologyPerformNNI(tr, p, 2);
	double lh2 = tr->likelihood;
	/* Update branch lengths */
	pllUpdatePartials(tr, pr, p, PLL_FALSE);
	pllUpdatePartials(tr, pr, q, PLL_FALSE);
	update(tr, pr, p);
	pllEvaluateLikelihood (tr, pr, p, PLL_FALSE, PLL_FALSE);

	// Create the nniMove struct to store this move
	nniMove nni2;
	nni2.p = p;
	nni2.nniType = 2;

	// Store the optimized and unoptimized central branch length
	for (i = 0; i < pr->numberOfPartitions; i++) {
		nni2.z[i] = p->z[i];
		nni2.z0[i] = z0[i];
	}
	nni2.likelihood = lh2;
	nni2.deltaLH = lh2 - lh0;
#ifdef _DEBUG_NNI
	printf("Delta likelihood of the 2.NNI move: %f\n", nni2.deltaLH);
#endif

	/* Restore previous NNI move */
	pllTopologyPerformNNI(tr, p, 2);
	pllUpdatePartials(tr, pr, p, PLL_FALSE);
	pllUpdatePartials(tr, pr, p->back, PLL_FALSE);
	/* Restore the old branch length */
	for (i = 0; i < pr->numberOfPartitions; i++) {
		p->z[i] = z0[i];
		p->back->z[i] = z0[i];
	}
	if (nni1.deltaLH > 0 && nni1.deltaLH >= nni2.deltaLH) {
		return nni1;
	} else if (nni1.deltaLH > 0 && nni1.deltaLH < nni2.deltaLH) {
		return nni2;
	} else if (nni1.deltaLH < 0 && nni2.deltaLH > 0) {
		return nni2;
	} else {
		return nni0;
	}
}

/** @brief ??? Not sure */
void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p,
		nniMove* nniList, int* cnt, int* cnt_nni, double curLH) {
	if (!isTip(p->number, tr->mxtips)) {
		nniList[*cnt] = getBestNNIForBran(tr, pr, p, curLH);
		if (nniList[*cnt].deltaLH != 0.0) {
			*cnt_nni = *cnt_nni + 1;
		}
		*cnt = *cnt + 1;
		nodeptr q = p->next;
		while (q != p) {
			evalNNIForSubtree(tr, pr, q->back, nniList, cnt, cnt_nni, curLH);
			q = q->next;
		}
	}
}

/** @brief Perform an NNI search

    Modify the tree topology of instance and model parameters \a tr by performing a NNI (Neighbour Neighbor
    Interchange) moves \a p.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param estimateModel
      Determine wheter the model parameters should be optimized

    @return
      In case of success \b PLL_TRUE, otherwise \b PLL_FALSE

*/
int pllNniSearch(pllInstance * tr, partitionList *pr, int estimateModel) {

	double curScore = tr->likelihood;

	/* Initialize the NNI list */
	nniMove* nniList = (nniMove*) malloc((tr->mxtips - 3) * sizeof(nniMove));
	int i;
	/* fill up the NNI list */
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	int cnt = 0; // number of visited internal branches during NNI evaluation
	int cnt_nni = 0; // number of positive NNI found
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, nniList, &cnt, &cnt_nni, curScore);
		q = q->next;
	}
    if (cnt_nni == 0) {
        free(nniList); //James B. 23-Jul-2020 (memory leak)
        return 0.0;
    }

	nniMove* impNNIList = (nniMove*) malloc(cnt_nni * sizeof(nniMove));
	int j = 0;
	for (i = 0; i < tr->mxtips - 3; i++) {
		if (nniList[i].deltaLH > 0.0) {
			impNNIList[j] = nniList[i];
			j++;
		}
	}
	// sort impNNIList
	qsort(impNNIList, cnt_nni, sizeof(nniMove), cmp_nni);

	// creating a list of non-conflicting positive NNI
	nniMove* nonConfNNIList = (nniMove*) calloc(cnt_nni, sizeof(nniMove));

	// the best NNI will always be taken
	nonConfNNIList[0] = impNNIList[cnt_nni - 1];

	// Filter out conflicting NNI
	int numNonConflictNNI = 1; // size of the non-conflicting NNI list;
	int k;
	for (k = cnt_nni - 2; k >= 0; k--) {
		int conflict = PLL_FALSE;
		int j;
		for (j = 0; j < numNonConflictNNI; j++) {
			if (impNNIList[k].p->number == nonConfNNIList[j].p->number
					|| impNNIList[k].p->number
							== nonConfNNIList[j].p->back->number) {
				conflict = PLL_TRUE;
				break;
			}
		}
		if (conflict) {
			continue;
		} else {
			nonConfNNIList[numNonConflictNNI] = impNNIList[k];
			numNonConflictNNI++;
		}
	}

	// Applying non-conflicting NNI moves
	double delta = 1.0; // portion of NNI moves to apply
	int notImproved;
	do {
		notImproved = PLL_FALSE;
		int numNNI2Apply = ceil(numNonConflictNNI * delta);
		for (i = 0; i < numNNI2Apply; i++) {
			// Just do the topological change
			pllTopologyPerformNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType);
			pllUpdatePartials(tr, pr, nonConfNNIList[i].p, PLL_FALSE);
			pllUpdatePartials(tr, pr, nonConfNNIList[i].p->back, PLL_FALSE);
			// Apply the store branch length
			int j;
			for (j = 0; j < pr->numberOfPartitions; j++) {
				nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z[j];
				nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z[j];
			}
		}
		// Re-optimize all branches
		smoothTree(tr, pr, 2);
		pllEvaluateLikelihood (tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
		if (estimateModel) {
			modOpt(tr, pr, 0.1);
		}
		pllEvaluateLikelihood (tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
		if (tr->likelihood < curScore) {
#ifdef _DEBUG_NNI
			printf("Tree likelihood gets worse after applying NNI\n");
			printf("curScore = %30.20f\n", curScore);
			printf("newScore = %30.20f\n", tr->likelihood);
			printf("Rolling back the tree\n");
#endif
			for (i = 0; i < numNNI2Apply; i++) {
				pllTopologyPerformNNI(tr, nonConfNNIList[i].p, nonConfNNIList[i].nniType);
				// Restore the branch length
				int j;
				for (j = 0; j < pr->numberOfPartitions; j++) {
					nonConfNNIList[i].p->z[j] = nonConfNNIList[i].z0[j];
					nonConfNNIList[i].p->back->z[j] = nonConfNNIList[i].z0[j];
				}
			}
			pllEvaluateLikelihood (tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
#ifdef _DEBUG_NNI
			printf("Tree likelihood after rolling back = %f \n",
					tr->likelihood);
#endif
			notImproved = PLL_TRUE & (numNNI2Apply > 1);
			delta = delta * 0.5;
		}
	} while (notImproved);
	free(nniList);
	free(impNNIList);
	free(nonConfNNIList);

	return PLL_TRUE;
}


/** @defgroup rearrangementGroup Topological rearrangements
    
    This set of functions handles the rearrangement of the tree topology
*/


/** @ingroup rearrangementGroup
    @brief Create a list for storing topology rearrangements
 
    Allocates space and initializes a structure that will hold information
    of \a max topological rearrangements

    @param max
      Maximum number of elements that the structure should hold
    
    @note This should be called for creating a storage space (list) for
    routines such as ::pllRearrangeSearch which compute the best NNI/PR/TBR rearrangements.
*/
pllRearrangeList * pllCreateRearrangeList (int max)
{
  pllRearrangeList * bl;

  bl = (pllRearrangeList *) malloc (sizeof (pllRearrangeList));

  bl->max_entries = max;
  bl->entries     = 0;
  bl->rearr       = (pllRearrangeInfo *) malloc (max * sizeof (pllRearrangeInfo));

  return bl;
}

/** @ingroup rearrangementGroup
    @brief Deallocator for topology rearrangements list
    
    Call this to destroy (deallocate) the memory taken by the \a bestList which holds
    topological rearrangements

    @param bestList
      Pointer to the list to be deallocated
*/
void pllDestroyRearrangeList (pllRearrangeList ** bestList)
{
  pllRearrangeList * bl;

  bl = *bestList;

  rax_free (bl->rearr);
  rax_free (bl);

  *bestList = NULL;
}


/** @ingroup rearrangementGroup
    @brief Store a rearrangement move to the list of best rearrangement moves

     Checks if the likelihood yielded by the rearrangement move described in \a rearr
     is better than any in the sorted list \a bestList. If it is, or
     if there is still space in \a bestList, the info about the
     move is inserted in the list.

     @param bestList
       The list of information about the best rearrangement moves

     @param rearr
       Info about the current rearrangement move

     @return
       Returns \b PLL_FALSE if the rearrangement move doesn't make it in the list, otherwise \b PLL_TRUE
*/
static int pllStoreRearrangement (pllRearrangeList * bestList, pllRearrangeInfo * rearr)
 {
   /* naive implementation of saving rearrangement moves */
   int i;

   for (i = 0; i < bestList->entries; ++ i)
    {
      /* Does the new rearrangement yield a better likelihood that the current in the list */
      if (rearr->likelihood > bestList->rearr[i].likelihood)
       {
         /* is there enough space in the array ? */
         if (bestList->entries < bestList->max_entries)
          {
            /* slide the entries to the right and overwrite the i-th element with the new item */
            memmove (&(bestList->rearr[i + 1]), &(bestList->rearr[i]), (bestList->entries - i ) * sizeof (pllRearrangeInfo));
            ++ bestList->entries;
          }
         else
          {
            memmove (&(bestList->rearr[i + 1]), &(bestList->rearr[i]), (bestList->entries - i - 1 ) * sizeof (pllRearrangeInfo));
          }
         memcpy (&(bestList->rearr[i]), rearr, sizeof (pllRearrangeInfo));
         return (PLL_TRUE);
       }
    }
   if (bestList->entries < bestList->max_entries)
    {
      memcpy (&(bestList->rearr[bestList->entries]), rearr, sizeof (pllRearrangeInfo));
      ++ bestList->entries;
      return (PLL_TRUE);
    }

   return (PLL_FALSE);
 }

/** @ingroup rearrangementGroup
    @brief Internal function for testing and saving an SPR move
    
    Checks the likelihood of the placement of the pruned subtree specified by \a p
    to node \a q. If the likelihood is better than some in the sorted list 
    \a bestList, or if there is still free space in \a bestList, then the SPR 
    move is recorded (in \a bestList)

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Root of the subtree that is to be pruned

    @param q
      Where to place the pruned subtree (between \a q and \a q->back

    @param bestList
      Where to store the SPR move

    @note Internal function which is not part of the PLL API and therefore should not be
    called by the user

    @return
*/
static int
pllTestInsertBIG (pllInstance * tr, partitionList * pr, nodeptr p, nodeptr q, pllRearrangeList * bestList)
{
  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;
  pllRearrangeInfo rearr;

  double  qz[PLL_NUM_BRANCHES], pz[PLL_NUM_BRANCHES];
  nodeptr  r;
  //double startLH = tr->endLH;
  int i;

  r = q->back; 
  for(i = 0; i < numBranches; i++)
  {
    qz[i] = q->z[i];
    pz[i] = p->z[i];
  }

  if (! insertBIG(tr, pr, p, q))       return PLL_FALSE;

  pllEvaluateLikelihood (tr, pr, p->next->next, PLL_FALSE, PLL_FALSE);
  
  rearr.rearrangeType  = PLL_REARRANGE_SPR;
  rearr.likelihood     = tr->likelihood;
  rearr.SPR.removeNode = p;
  rearr.SPR.insertNode = q;
  for (i = 0; i < numBranches; ++ i)
   {
     rearr.SPR.zqr[i] = tr->zqr[i];
   }

  pllStoreRearrangement (bestList, &rearr);

/*
  if(tr->likelihood > tr->bestOfNode)
  {
    pllStoreRearrangement (bestList, rearr)
    tr->bestOfNode = tr->likelihood;
    tr->insertNode = q;
    tr->removeNode = p;   
    for(i = 0; i < numBranches; i++)
    {
      tr->currentZQR[i] = tr->zqr[i];           
      tr->currentLZR[i] = tr->lzr[i];
      tr->currentLZQ[i] = tr->lzq[i];
      tr->currentLZS[i] = tr->lzs[i];      
    }
  }

  if(tr->likelihood > tr->endLH)
  {			  
    
    tr->insertNode = q;
    tr->removeNode = p;   
    for(i = 0; i < numBranches; i++)
      tr->currentZQR[i] = tr->zqr[i];      
    tr->endLH = tr->likelihood;                      
  }        
*/
  /* reset the topology so that it is the same as it was before calling insertBIG */
  hookup(q, r, qz, numBranches);

  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(tr->thoroughInsertion)
  {
    nodeptr s = p->back;
    hookup(p, s, pz, numBranches);
  } 

/*
  if((tr->doCutoff) && (tr->likelihood < startLH))
  {
    tr->lhAVG += (startLH - tr->likelihood);
    tr->lhDEC++;
    if((startLH - tr->likelihood) >= tr->lhCutoff)
      return PLL_FALSE;	    
    else
      return PLL_TRUE;
  }
  else
    return PLL_TRUE;
  */
  return (PLL_TRUE);
}

/** @ingroup rearrangementGroup
    @brief Internal function for recursively traversing a tree and testing a possible subtree insertion

    Recursively traverses the tree rooted at \a q in the direction of \a q->next->back and \a q->next->next->back
    and at each node tests the placement of the pruned subtree rooted at \a p by calling the function
    \a pllTestInsertBIG, which in turn saves the computed SPR in \a bestList if a) there is still space in
    the \a bestList or b) if the likelihood of the SPR is better than any of the ones in \a bestList.

    @note This function is not part of the API and should not be called by the user.
*/
static void pllTraverseUpdate (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav, pllRearrangeList * bestList)
{  
  if (--mintrav <= 0) 
  {              
    if (! pllTestInsertBIG(tr, pr, p, q, bestList))  return;

  }

  if ((!isTip(q->number, tr->mxtips)) && (--maxtrav > 0)) 
  {    
    pllTraverseUpdate(tr, pr, p, q->next->back, mintrav, maxtrav, bestList);
    pllTraverseUpdate(tr, pr, p, q->next->next->back, mintrav, maxtrav, bestList);
  }
} 


/** @ingroup rearrangementGroup
    @brief Internal function for computing SPR moves

    Compute a list of at most \a max SPR moves that can be performed by pruning
    the subtree rooted at node \a p and testing all possible placements in a
    radius of at least \a mintrav nodes and at most \a maxtrav nodes from \a p.
    Note that \a tr->thoroughInsertion affects the behaviour of the function (see note).

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Node specifying the root of the pruned subtree, i.e. where to prune.

    @param mintrav
      Minimum distance from \a p where to try inserting the pruned subtree

    @param maxtrav
      Maximum distance from \a p where to try inserting the pruned subtree

    @param bestList
      The list of best topological rearrangements

    @note This function is not part of the API and should not be called by the user
    as it is called internally by the API function \a pllComputeSPR. 
    Also, setting \a tr->thoroughInsertion affects this function. For each tested SPR
    the new branch lengths will also be optimized. This computes better likelihoods
    but also slows down the method considerably.
*/
static int pllTestSPR (pllInstance * tr, partitionList * pr, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList)
{
  nodeptr 
    p1, p2, q, q1, q2;
  double 
    p1z[PLL_NUM_BRANCHES], p2z[PLL_NUM_BRANCHES], q1z[PLL_NUM_BRANCHES], q2z[PLL_NUM_BRANCHES];
  int
    mintrav2, i;
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  if (maxtrav < 1 || mintrav > maxtrav) return (PLL_FALSE);
  q = p->back;

  if (!isTip (p->number, tr->mxtips))
   {
     p1 = p->next->back;
     p2 = p->next->next->back;

     if (!isTip (p1->number, tr->mxtips) || !isTip (p2->number, tr->mxtips))
      {
        /* save branch lengths before splitting the tree in two components */
        for (i = 0; i < numBranches; ++ i)
         {
           p1z[i] = p1->z[i];
           p2z[i] = p2->z[i];
         }

        /* split the tree in two components */
        if (! removeNodeBIG (tr, pr, p, numBranches)) return PLL_BADREAR;

        /* recursively traverse and perform SPR on the subtree rooted at p1 */
        if (!isTip (p1->number, tr->mxtips))
         {
           pllTraverseUpdate (tr, pr, p, p1->next->back,       mintrav, maxtrav, bestList);
           pllTraverseUpdate (tr, pr, p, p1->next->next->back, mintrav, maxtrav, bestList);
         }

        /* recursively traverse and perform SPR on the subtree rooted at p2 */
        if (!isTip (p2->number, tr->mxtips))
         {
           pllTraverseUpdate (tr, pr, p, p2->next->back,       mintrav, maxtrav, bestList);
           pllTraverseUpdate (tr, pr, p, p2->next->next->back, mintrav, maxtrav, bestList);
         }

        /* restore the topology as it was before the split */
        hookup (p->next,       p1, p1z, numBranches);
        hookup (p->next->next, p2, p2z, numBranches);
        pllUpdatePartials (tr, pr, p, PLL_FALSE);
      }
   }

  if (!isTip (q->number, tr->mxtips) && maxtrav > 0)
   {
     q1 = q->next->back;
     q2 = q->next->next->back;

    /* why so many conditions? Why is it not analogous to the previous if for node p? */
    if (
        (
         ! isTip(q1->number, tr->mxtips) && 
         (! isTip(q1->next->back->number, tr->mxtips) || ! isTip(q1->next->next->back->number, tr->mxtips))
        )
        ||
        (
         ! isTip(q2->number, tr->mxtips) && 
         (! isTip(q2->next->back->number, tr->mxtips) || ! isTip(q2->next->next->back->number, tr->mxtips))
        )
       )
     {
       for (i = 0; i < numBranches; ++ i)
        {
          q1z[i] = q1->z[i];
          q2z[i] = q2->z[i];
        }

       if (! removeNodeBIG (tr, pr, q, numBranches)) return PLL_BADREAR;

       mintrav2 = mintrav > 2 ? mintrav : 2;

       if (!isTip (q1->number, tr->mxtips))
        {
          pllTraverseUpdate (tr, pr, q, q1->next->back,       mintrav2, maxtrav, bestList);
          pllTraverseUpdate (tr, pr, q, q1->next->next->back, mintrav2, maxtrav, bestList);
        }

       if (!isTip (q2->number, tr->mxtips))
        {
          pllTraverseUpdate (tr, pr, q, q2->next->back,       mintrav2, maxtrav, bestList);
          pllTraverseUpdate (tr, pr, q, q2->next->next->back, mintrav2, maxtrav, bestList);
        }

       hookup (q->next,       q1, q1z, numBranches);
       hookup (q->next->next, q2, q2z, numBranches);
       pllUpdatePartials (tr, pr, q, PLL_FALSE);
     }
   }
  return (PLL_TRUE);
}

/** @ingroup rearrangementGroup
    @brief Compute a list of possible SPR moves
    
    Iteratively tries all possible SPR moves that can be performed by
    pruning the subtree rooted at \a p and testing all possible placements
    in a radius of at least \a mintrav nodea and at most \a maxtrav nodes from
    \a p. Note that \a tr->thoroughInsertion affects the behaviour of the function (see note).

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Node specifying the root of the pruned subtree, i.e. where to prune.

    @param mintrav
      Minimum distance from \a p where to try inserting the pruned subtree

    @param maxtrav
      Maximum distance from \a p where to try inserting the pruned subtree

    @note
      Setting \a tr->thoroughInsertion affects this function. For each tested SPR
      the new branch lengths will also be optimized. This computes better likelihoods
      but also slows down the method considerably.
*/
static void 
pllComputeSPR (pllInstance * tr, partitionList * pr, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList)
{

  tr->startLH = tr->endLH = tr->likelihood;

  /* TODO: Add cutoff code */

  tr->bestOfNode = PLL_UNLIKELY;

  pllTestSPR (tr, pr, p, mintrav, maxtrav, bestList);
}

/** @ingroup rearrangementGroup
    @brief Return the yielded likelihood of an NNI move, without altering the topology

    This function performs the NNI move of type \a swapType at node \a p, optimizes
    the branch with endpoints \a p  and \a p->back and evalutes the resulting likelihood.
    It then restores the topology  to the origin and returns the likelihood that the NNI
    move yielded.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Where to perform the NNI move

    @param swapType
      What type of NNI move to perform

    @return
      The likelihood yielded from the NNI
*/
static double 
pllTestNNILikelihood (pllInstance * tr, partitionList * pr, nodeptr p, int swapType)
{
  double lh;
  double z0[PLL_NUM_BRANCHES];
  int i;

  /* store the origin branch lengths and likelihood. The original branch lengths could
  be passed as a parameter in order to avoid duplicate computations because of the two
  NNI moves */
  for (i = 0; i < pr->numberOfPartitions; ++ i)
   {
     z0[i] = p->z[i];
   }

  /* perform NNI */
  pllTopologyPerformNNI(tr, p, swapType);
  /* recompute the likelihood vectors of the two subtrees rooted at p and p->back,
     optimize the branch lengths and evaluate the likelihood  */
  pllUpdatePartials (tr, pr, p,       PLL_FALSE);
  pllUpdatePartials (tr, pr, p->back, PLL_FALSE);
  update (tr, pr, p);
  pllEvaluateLikelihood (tr, pr, p, PLL_FALSE, PLL_FALSE);
  lh = tr->likelihood;

  /* restore topology */
  pllTopologyPerformNNI(tr, p, swapType);
  pllUpdatePartials (tr, pr, p,       PLL_FALSE);
  pllUpdatePartials (tr, pr, p->back, PLL_FALSE);
  //update (tr, pr, p);
  pllEvaluateLikelihood (tr, pr, p, PLL_FALSE, PLL_FALSE);
  for (i = 0; i < pr->numberOfPartitions; ++ i)
   {
     p->z[i] = p->back->z[i] = z0[i];
   }

  return lh;
}
/** @ingroup rearrangementGroup
    @brief Compares NNI likelihoods at a node and store in the rearrangement list

    Compares the two possible NNI moves that can be performed at node \a p, and
    if the likelihood improves from the one of the original topology, then 
    it picks the one that yields the highest likelihood and tries to insert it in
    the list of best rearrangement moves

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param bestList
      Rearrangement moves list
*/
static void pllTestNNI (pllInstance * tr, partitionList * pr, nodeptr p, pllRearrangeList * bestList)
{
  double lh0, lh1, lh2;
  pllRearrangeInfo rearr;

  /* store the original likelihood */
  lh0 = tr->likelihood;

  lh1 = pllTestNNILikelihood (tr, pr, p, PLL_NNI_P_NEXT);
  lh2 = pllTestNNILikelihood (tr, pr, p, PLL_NNI_P_NEXTNEXT);

  if (lh0 > lh1 && lh0 > lh2) return;

  /* set the arrangement structure */
  rearr.rearrangeType  = PLL_REARRANGE_NNI;
  rearr.likelihood     = PLL_MAX (lh1, lh2);
  rearr.NNI.originNode = p;
  rearr.NNI.swapType   = (lh1 > lh2) ? PLL_NNI_P_NEXT : PLL_NNI_P_NEXTNEXT;

  /* try to store it in the best list */
  pllStoreRearrangement (bestList, &rearr);
}

/** @ingroup rearrangementGroup
    @brief Recursive traversal of the tree structure for testing NNI moves
 
    Recursively traverses the tree structure and tests all allowed NNI
    moves in the area specified by \a mintrav and \a maxtrav. For more
    information and details on the function arguments check ::pllSearchNNI
*/
static void 
pllTraverseNNI (pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList)
{
  if (isTip (p->number, tr->mxtips)) return;

  /* if we are at the right radius then compute the NNIs for nodes p->next and p->next->next */
  if (!mintrav)
   {
     pllTestNNI (tr, pr, p->next, bestList);
     pllTestNNI (tr, pr, p->next->next, bestList);
   }
  
  /* and then avoid computing the NNIs for nodes p->next->back and p->next->next->back as they are
  the same to the ones computed in the above two lines. This way we do not need to resolve conflicts
  later on as in the old code */
  if (maxtrav)
   {
     if (!isTip (p->next->back->number, tr->mxtips))       
       pllTraverseNNI (tr, pr, p->next->back,       mintrav ? mintrav - 1 : 0, maxtrav - 1, bestList);
     if (!isTip (p->next->next->back->number, tr->mxtips)) 
       pllTraverseNNI (tr, pr, p->next->next->back, mintrav ? mintrav - 1 : 0, maxtrav - 1, bestList);
   }
}

/** @ingroup rearrangementGroup
    @brief Compute a list of possible NNI moves
    
    Iteratively tries all possible NNI moves at each node that is at
    least \a mintrav and at most \a maxtrav nodes far from node \a p.
    At each NNI move, the likelihood is tested and if it is higher than
    the likelihood of an element in the sorted (by likelihood) list 
    \a bestList, or if there is still empty space in \a bestList, it is
    inserted at the corresponding position.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param p
      Node specifying the point where the NNI will be performed.

    @param mintrav
      Minimum distance from \a p where the NNI can be tested 

    @param maxtrav
      Maximum distance from \a p where to try NNIs
*/
static void
pllSearchNNI (pllInstance * tr, partitionList * pr, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList)
{
  /* avoid conflicts by precomputing the NNI of the first node */

  if (mintrav == 0) 
  pllTestNNI (tr, pr, p, bestList);
  
  pllTraverseNNI (tr, pr, p, mintrav, maxtrav, bestList);
  if (maxtrav)
    pllTraverseNNI (tr, pr, p->back, mintrav, maxtrav - 1, bestList);

}

/** @ingroup rearrangementGroup
    @brief Create rollback information for an SPR move
    
    Creates a structure of type ::pllRollbackInfo and fills it with rollback
    information about the SPR move described in \a rearr. The rollback info
    is stored in the PLL instance in a LIFO manner.

    @param tr
      PLL instance

    @param rearr
      Description of the SPR move

    @param numBranches
      Number of partitions
*/
static void 
pllCreateSprInfoRollback (pllInstance * tr, pllRearrangeInfo * rearr, int numBranches)
{
  pllRollbackInfo * sprRb;
  nodeptr p, q;
  int i;

  p = rearr->SPR.removeNode;
  q = rearr->SPR.insertNode;

  sprRb = (pllRollbackInfo *) rax_malloc (sizeof (pllRollbackInfo) + 4 * numBranches * sizeof (double));
  sprRb->SPR.zp   = (double *) ((char *)sprRb + sizeof (pllRollbackInfo));
  sprRb->SPR.zpn  = (double *) ((char *)sprRb + sizeof (pllRollbackInfo) + numBranches * sizeof (double));
  sprRb->SPR.zpnn = (double *) ((char *)sprRb + sizeof (pllRollbackInfo) + 2 * numBranches * sizeof (double));
  sprRb->SPR.zqr  = (double *) ((char *)sprRb + sizeof (pllRollbackInfo) + 3 * numBranches * sizeof (double));

  for (i = 0; i < numBranches; ++ i)
   {
     sprRb->SPR.zp[i]   = p->z[i];
     sprRb->SPR.zpn[i]  = p->next->z[i];
     sprRb->SPR.zpnn[i] = p->next->next->z[i];
     sprRb->SPR.zqr[i]  = q->z[i];
   }

  sprRb->SPR.pn  = p->next->back;
  sprRb->SPR.pnn = p->next->next->back;
  sprRb->SPR.r   = q->back;
  sprRb->SPR.q   = q;
  sprRb->SPR.p   = p;

  sprRb->rearrangeType = PLL_REARRANGE_SPR;

  pllStackPush (&(tr->rearrangeHistory), (void *) sprRb);
}

/** @ingroup rearrangementGroup
    @brief Create rollback information for an NNI move

    Creates a structure of type ::pllRollbackInfo and fills it with rollback
    information about the SPR move described in \a rearr. The rollback info
    is stored in the PLL instance in a LIFO manner

    @param tr
      PLL instance

    @param rearr
      Description of the NNI move
*/
static void
pllCreateNniInfoRollback (pllInstance * tr, pllRearrangeInfo * rearr)
{
  /*TODO: add the branches ? */
  pllRollbackInfo * ri;

  ri = (pllRollbackInfo *) rax_malloc (sizeof (pllRollbackInfo));

  ri->rearrangeType = PLL_REARRANGE_NNI;

  ri->NNI.origin   = rearr->NNI.originNode;
  ri->NNI.swapType = rearr->NNI.swapType;

  pllStackPush (&(tr->rearrangeHistory), (void *) ri);
  
}


/** @ingroup rearrangementGroup
    @brief Generic function for creating rollback information

    Creates a structure of type ::pllRollbackInfo and fills it with rollback
    information about the move described in \a rearr. The rollback info
    is stored in the PLL instance in a LIFO manner

    @param tr
      PLL instance

    @param rearr
      Description of the NNI move

    @param numBranches
      Number of partitions
*/
static void
pllCreateRollbackInfo (pllInstance * tr, pllRearrangeInfo * rearr, int numBranches)
{
  switch (rearr->rearrangeType)
   {
     case PLL_REARRANGE_NNI:
       pllCreateNniInfoRollback (tr, rearr);
       break;
     case PLL_REARRANGE_SPR:
       pllCreateSprInfoRollback (tr, rearr, numBranches);
       break;
     default:
       break;
   }

}


/** @ingroup rearrangementGroup
    @brief Rollback an SPR move

    Perform a rollback (undo) on the last SPR move.
    
    @param tr
      PLL instance

    @param pr
      List of partitions

    @param ri
      Rollback information
*/
static void
pllRollbackSPR (partitionList * pr, pllRollbackInfo * ri)
{
  int numBranches;

  numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  hookup (ri->SPR.p->next,       ri->SPR.pn,      ri->SPR.zpn,  numBranches);
  hookup (ri->SPR.p->next->next, ri->SPR.pnn,     ri->SPR.zpnn, numBranches); 
  hookup (ri->SPR.p,             ri->SPR.p->back, ri->SPR.zp,   numBranches);
  hookup (ri->SPR.q,             ri->SPR.r,       ri->SPR.zqr,  numBranches);

  rax_free (ri);
}

/** @ingroup rearrangementGroup
    @brief Rollback an NNI move

    Perform a rollback (undo) on the last NNI move.
    
    @param tr
      PLL instance

    @param pr
      List of partitions

    @param ri
      Rollback information
*/
static void
pllRollbackNNI (pllInstance * tr, partitionList * pr, pllRollbackInfo * ri)
{
  nodeptr p = ri->NNI.origin;

  pllTopologyPerformNNI(tr, p, ri->NNI.swapType);
  pllUpdatePartials (tr, pr, p,       PLL_FALSE);
  pllUpdatePartials (tr, pr, p->back, PLL_FALSE);
  update (tr, pr, p);
  pllEvaluateLikelihood (tr, pr, p, PLL_FALSE, PLL_FALSE);
  
  rax_free (ri);
}

/** @ingroup rearrangementGroup
    @brief Rollback the last committed rearrangement move
    
    Perform a rollback (undo) on the last committed rearrangement move.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @return
      Returns \b PLL_TRUE is the rollback was successful, otherwise \b PLL_FALSE
      (if no rollback was done)
*/
int 
pllRearrangeRollback (pllInstance * tr, partitionList * pr)
{
  pllRollbackInfo * ri;
  
  ri = (pllRollbackInfo *) pllStackPop (&(tr->rearrangeHistory));
  if (!ri) return (PLL_FALSE);

  switch (ri->rearrangeType)
   {
     case PLL_REARRANGE_NNI:
       pllRollbackNNI (tr, pr, ri);
       break;
     case PLL_REARRANGE_SPR:
       pllRollbackSPR (pr, ri);
       break;
     default:
       rax_free (ri);
       return (PLL_FALSE);
   }

  return (PLL_TRUE);
  
}


/** @ingroup rearrangementGroup
    @brief Commit a rearrangement move

    Applies the rearrangement move specified in \a rearr to the tree topology in \a tr. 
    In case of SPR moves, if
    \a tr->thoroughInsertion is set to \b PLL_TRUE, the new branch lengths are also optimized. 
    The function stores rollback information in pllInstance::rearrangeHistory if \a saveRollbackInfo
    is set to \b PLL_TRUE. This way, the rearrangement move can be rolled back (undone) by calling
    ::pllRearrangeRollback

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param rearr
      An element of a \a pllRearrangeInfo structure that contains information about the rearrangement move

    @param saveRollbackInfo
      If set to \b PLL_TRUE, rollback info will be kept for undoing the rearrangement move
*/
void
pllRearrangeCommit (pllInstance * tr, partitionList * pr, pllRearrangeInfo * rearr, int saveRollbackInfo)
{
  int numBranches;

  numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  if (saveRollbackInfo)
    pllCreateRollbackInfo (tr, rearr, numBranches);

  switch (rearr->rearrangeType)
   {
     case PLL_REARRANGE_NNI:
       pllTopologyPerformNNI(tr, rearr->NNI.originNode, rearr->NNI.swapType);
       pllUpdatePartials (tr, pr, rearr->NNI.originNode, PLL_FALSE);
       pllUpdatePartials (tr, pr, rearr->NNI.originNode->back, PLL_FALSE);
       update (tr, pr, rearr->NNI.originNode);
       pllEvaluateLikelihood (tr, pr, rearr->NNI.originNode, PLL_FALSE, PLL_FALSE);
       break;
     case PLL_REARRANGE_SPR:
       removeNodeBIG (tr, pr, rearr->SPR.removeNode, numBranches);
       insertBIG     (tr, pr, rearr->SPR.removeNode, rearr->SPR.insertNode);
       break;
     default:
       break;
   }
}


/******** new rearrangement functions ****************/

/* change this to return the number of new elements in the list */
/** @ingroup rearrangementGroup
    @brief Search for rearrangement topologies
    
    Search for possible rearrangement moves of type \a rearrangeType in the
    annular area defined by the minimal resp. maximal radii \a mintrav resp.
    \a maxtrav. If the resulting likelihood is better than the current, try
    to insert the move specification in \a bestList, which is a sorted list
    that holds the rearrange info of the best moves sorted by likelihood
    (desccending order).

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param rearrangeType
      Type of rearrangement. Can be \b PLL_REARRANGE_SPR or \b PLL_REARRANGE_NNI

    @param p
      Point of origin, i.e. where to start searching from

    @param mintrav
      The minimal radius of the annulus

    @param maxtrav
      The maximal radius of the annulus

    @param bestList
      List that holds the details of the best rearrangement moves found

    @note
      If \a bestList is not empty, the existing entries will not be altered unless
      better rearrangement moves (that means yielding better likelihood) are found
      and the list is full, in which case the entries with the worst likelihood will be
      thrown away.
*/
void
pllRearrangeSearch (pllInstance * tr, partitionList * pr, int rearrangeType, nodeptr p, int mintrav, int maxtrav, pllRearrangeList * bestList)
{
  switch (rearrangeType)
   {
     case PLL_REARRANGE_SPR:
       pllComputeSPR (tr, pr, p, mintrav, maxtrav, bestList);
       break;

     case PLL_REARRANGE_NNI:
       pllSearchNNI (tr, pr, p, mintrav, maxtrav, bestList);
       break;

     case PLL_REARRANGE_TBR:
       break;
     default:
       break;
   }
}


static int
determineRearrangementSetting(pllInstance *tr, partitionList *pr,
    bestlist *bestT, bestlist *bt)
{
  int i, mintrav, maxtrav, bestTrav, impr, index, MaxFast, *perm = (int*) NULL;
  double startLH;
  pllBoolean cutoff;

  MaxFast = 26;

  startLH = tr->likelihood;

  cutoff = tr->doCutoff;
  tr->doCutoff = PLL_FALSE;

  mintrav = 1;

  bestTrav = maxtrav = 5;

  impr = 1;

  resetBestTree(bt);

  if (tr->permuteTreeoptimize)
    {
      int n = tr->mxtips + tr->mxtips - 2;
      perm = (int *) rax_malloc(sizeof(int) * (n + 1));
      makePermutation(perm, n, tr);
    }

  while (impr && maxtrav < MaxFast)
    {
      recallBestTree(bestT, 1, tr, pr);
      nodeRectifier(tr);

      if (maxtrav > tr->ntips - 3)
        maxtrav = tr->ntips - 3;

      tr->startLH = tr->endLH = tr->likelihood;

      for (i = 1; i <= tr->mxtips + tr->mxtips - 2; i++) {
          if (perm != NULL)
              index = perm[i];
          else
              index = i;


          tr->bestOfNode = PLL_UNLIKELY;
          if (rearrangeBIG(tr, pr, tr->nodep[index], mintrav, maxtrav))
            {
              if (tr->endLH > tr->startLH)
                {
                  restoreTreeFast(tr, pr);
                  tr->startLH = tr->endLH = tr->likelihood;
                }
            }
        }

      pllOptimizeBranchLengths(tr, pr, 8);
      saveBestTree(bt, tr,
          pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

      if (tr->likelihood > startLH)
        {
          startLH = tr->likelihood;
          bestTrav = maxtrav;
          impr = 1;
        }
      else
        {
          impr = 0;
        }
      maxtrav += 5;

      if (tr->doCutoff)
        {
          tr->lhCutoff = (tr->lhAVG) / ((double) (tr->lhDEC));

          tr->itCount = tr->itCount + 1;
          tr->lhAVG = 0;
          tr->lhDEC = 0;
        }
    }

  recallBestTree(bt, 1, tr, pr);
  tr->doCutoff = cutoff;

  if (tr->permuteTreeoptimize)
    rax_free(perm);

  return bestTrav;
}


static void hash_dealloc_bipentry (void * entry)
{
  pllBipartitionEntry * e = (pllBipartitionEntry *)entry;

  if(e->bitVector)     rax_free(e->bitVector);
  if(e->treeVector)    rax_free(e->treeVector);
  if(e->supportVector) rax_free(e->supportVector);

}

/** @ingroup rearrangementGroup
    @brief RAxML algorithm for ML search

    RAxML algorithm for searching the Maximum Likelihood tree and model.

    @param tr
      PLL instance

    @param pr
      List of partitions

    @param estimateModel
      If true, model parameters are optimized in a ML framework.

    @note
      For datasets with a large number of taxa, setting tr->searchConvergenceCriterion to
    PLL_TRUE can improve the execution time in up to 50% looking for topology convergence.
*/
int
pllRaxmlSearchAlgorithm(pllInstance * tr, partitionList * pr,
    pllBoolean estimateModel)
{
  pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
  pllOptimizeBranchLengths(tr, pr, 32);

  unsigned int vLength = 0;
  int i, impr, bestTrav, rearrangementsMax = 0, rearrangementsMin = 0,
      thoroughIterations = 0, fastIterations = 0;

  double lh, previousLh, difference, epsilon;
  bestlist *bestT, *bt;
  infoList iList;
  pllOptimizeBranchLengths(tr, pr, 32);

  pllHashTable *h = NULL;
  //hashtable *h = NULL;
  unsigned int **bitVectors = (unsigned int**) NULL;

  /* Security check... These variables might have not been initialized! */
  if (tr->stepwidth == 0) tr->stepwidth = 5;
  if (tr->max_rearrange == 0) tr->max_rearrange = 21;

  if (tr->searchConvergenceCriterion)
    {
      bitVectors = initBitVector(tr->mxtips, &vLength);
      //h = initHashTable(tr->mxtips * 4);
      h = pllHashInit (tr->mxtips * 4);
    }

  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  bt = (bestlist *) rax_malloc(sizeof(bestlist));
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips);

  initInfoList(&iList, 50);

  epsilon = tr->likelihoodEpsilon;

  tr->thoroughInsertion = 0;

  if (estimateModel)
    {
      pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
      pllOptimizeModelParameters(tr, pr, 10.0);
    }
  else
    pllOptimizeBranchLengths(tr, pr, 64);

  saveBestTree(bestT, tr,
      pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

  if (!tr->initialSet)
    bestTrav = tr->bestTrav = determineRearrangementSetting(tr, pr, bestT, bt);
  else
    bestTrav = tr->bestTrav = tr->initial;

  if (estimateModel)
    {
      pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
      pllOptimizeModelParameters(tr, pr, 5.0);
    }
  else
    pllOptimizeBranchLengths(tr, pr, 32);

  saveBestTree(bestT, tr,
      pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
  impr = 1;
  if (tr->doCutoff)
    tr->itCount = 0;

  while (impr)
    {
      recallBestTree(bestT, 1, tr, pr);

      if (tr->searchConvergenceCriterion)
        {
          int bCounter = 0;

          if (fastIterations > 1)
            cleanupHashTable(h, (fastIterations % 2));

          bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips,
              vLength, h, fastIterations % 2, PLL_BIPARTITIONS_RF,
              (branchInfo *) NULL, &bCounter, 1, PLL_FALSE, PLL_FALSE, 0);

          assert(bCounter == tr->mxtips - 3);

          if (fastIterations > 0)
            {
              double rrf = convergenceCriterion(h, tr->mxtips);

              if (rrf <= 0.01) /* 1% cutoff */
                {
                  cleanupHashTable(h, 0);
                  cleanupHashTable(h, 1);
                  goto cleanup_fast;
                }
            }
        }

      fastIterations++;

      pllOptimizeBranchLengths(tr, pr, 32);

      saveBestTree(bestT, tr,
          pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

      lh = previousLh = tr->likelihood;

      treeOptimizeRapid(tr, pr, 1, bestTrav, bt, &iList);

      impr = 0;

      for (i = 1; i <= bt->nvalid; i++)
        {
          recallBestTree(bt, i, tr, pr);

          pllOptimizeBranchLengths(tr, pr, 8);

          difference = (
              (tr->likelihood > previousLh) ?
                  tr->likelihood - previousLh : previousLh - tr->likelihood);
          if (tr->likelihood > lh && difference > epsilon)
            {
              impr = 1;
              lh = tr->likelihood;
              saveBestTree(bestT, tr,
                  pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
            }
        }
    }

  if (tr->searchConvergenceCriterion)
    {
      cleanupHashTable(h, 0);
      cleanupHashTable(h, 1);
    }

  cleanup_fast:

  tr->thoroughInsertion = 1;
  impr = 1;

  recallBestTree(bestT, 1, tr, pr);
  if (estimateModel)
    {
      pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
      pllOptimizeModelParameters(tr, pr, 1.0);
    }
  else
    pllOptimizeBranchLengths(tr, pr, 32);

  while (1)
    {
      recallBestTree(bestT, 1, tr, pr);
      if (impr)
        {
          rearrangementsMin = 1;
          rearrangementsMax = tr->stepwidth;

          if (tr->searchConvergenceCriterion)
            {
              int bCounter = 0;

              if (thoroughIterations > 1)
                cleanupHashTable(h, (thoroughIterations % 2));

              bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back,
                  tr->mxtips, vLength, h, thoroughIterations % 2,
                  PLL_BIPARTITIONS_RF, (branchInfo *) NULL, &bCounter, 1,
                  PLL_FALSE, PLL_FALSE, 0);

              assert(bCounter == tr->mxtips - 3);

              if (thoroughIterations > 0)
                {
                  double rrf = convergenceCriterion(h, tr->mxtips);

                  if (rrf <= 0.01) /* 1% cutoff */
                    {
                      goto cleanup;
                    }
                }
            }

          thoroughIterations++;
        }
      else
        {
          rearrangementsMax += tr->stepwidth;
          rearrangementsMin += tr->stepwidth;
          if (rearrangementsMax > tr->max_rearrange)
            goto cleanup;
        }
      pllOptimizeBranchLengths(tr, pr, 32);

      previousLh = lh = tr->likelihood;
      saveBestTree(bestT, tr,
          pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);

      treeOptimizeRapid(tr, pr, rearrangementsMin, rearrangementsMax, bt,
          &iList);

      impr = 0;

      for (i = 1; i <= bt->nvalid; i++)
        {
          recallBestTree(bt, i, tr, pr);

          pllOptimizeBranchLengths(tr, pr, 8);

          difference = (
              (tr->likelihood > previousLh) ?
                  tr->likelihood - previousLh : previousLh - tr->likelihood);
          if (tr->likelihood > lh && difference > epsilon)
            {
              impr = 1;
              lh = tr->likelihood;
              saveBestTree(bestT, tr,
                  pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
            }
        }

    }

  cleanup:
  if (tr->searchConvergenceCriterion)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      rax_free(bitVectors);
      //freeHashTable(h);
      //rax_free(h);
      pllHashDestroy(&h, hash_dealloc_bipentry);
    }

  freeBestTree(bestT);
  rax_free(bestT);
  freeBestTree(bt);
  rax_free(bt);

  freeInfoList(&iList);

  if (estimateModel) {
      pllOptimizeModelParameters(tr, pr, epsilon);
  }
  pllOptimizeBranchLengths(tr, pr, 64);

  return 0;
}

