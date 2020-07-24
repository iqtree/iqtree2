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
 * @file recom.c
 * @brief Functions used for recomputation of vectors (only a fraction of LH vectors stored in RAM)   
 */
#include "mem_alloc.h"
#include "systypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "pll.h"
#include "pllInternal.h"

/** @brief Locks node \a nodenum to force it remains availably in memory
 *
 * @warning If a node is available we dont need to recompute it, but we neet to make sure it is not unpinned while buildding the rest of the traversal descriptor, i.e. unpinnable must be PLL_FALSE at this point, it will automatically be set to PLL_TRUE, after the counter post-order instructions have been executed 
Omitting this call the traversal will likely still work as long as num_allocated_nodes >> log n, but wrong inner vectors will be used at the wrong moment of pllNewviewIterative, careful! 
 *
 *  @param rvec 
 *    Recomputation info
 *
 *  @param nodenum
 *    Node id that must remain available in memory 
 *
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
void protectNode(recompVectors *rvec, int nodenum, int mxtips)
{

  int slot;
  slot = rvec->iNode[nodenum - mxtips - 1];
  assert(slot != PLL_NODE_UNPINNED);
  assert(rvec->iVector[slot] == nodenum);

  if(rvec->unpinnable[slot])
    rvec->unpinnable[slot] = PLL_FALSE;
}

/** @brief Checks if \a nodenum  is currently pinned (available in RAM)
 *
 *  @note shall we document static functions? 
 * 
 *  @param rvec 
 *    Recomputation info
 *
 *  @param nodenum
 *    Node id to be checked
 *
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
static pllBoolean isNodePinned(recompVectors *rvec, int nodenum, int mxtips)
{
  assert(nodenum > mxtips);

  if(rvec->iNode[nodenum - mxtips - 1] == PLL_NODE_UNPINNED)
    return PLL_FALSE;
  else
    return PLL_TRUE;
}

/** @brief Checks if the likelihood entries at node \a p should be updated
 *
 * A node needs update if one of the following holds:
 *    1. It is not oriented (p->x == 0) 
 *    2. We are applying recomputations and node \a p is not currently available in RAM
 *  
 *  @param recompute 
 *    PLL_TRUE if recomputation is currently applied 
 *
 *  @param p
 *    Node to check whether it is associated with the likelihood vector
 *
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
pllBoolean needsRecomp(pllBoolean recompute, recompVectors *rvec, nodeptr p, int mxtips)
{ 
  if((!p->x) || (recompute && !isNodePinned(rvec, p->number, mxtips)))
    return PLL_TRUE;
  else
    return PLL_FALSE;
}



/** @brief Allocates memory for recomputation structure
 *  
 *  
 *  @todo this should not depend on tr (\a vectorRecomFraction should be a parameter)
 *    PLL_TRUE if recomputation is currently applied 
 *
 */
void allocRecompVectorsInfo(pllInstance *tr)
{
  recompVectors 
    *v = (recompVectors *) rax_malloc(sizeof(recompVectors));

  int 
    num_inner_nodes = tr->mxtips - 2,
                    num_vectors, 
                    i;

  assert(tr->vectorRecomFraction > PLL_MIN_RECOM_FRACTION);
  assert(tr->vectorRecomFraction < PLL_MAX_RECOM_FRACTION);

  num_vectors = (int) (1 + tr->vectorRecomFraction * (float)num_inner_nodes); 

  int theoretical_minimum_of_vectors = 3 + ((int)(log((double)tr->mxtips)/log(2.0)));
  //printBothOpen("Try to use %d ancestral vectors, min required %d\n", num_vectors, theoretical_minimum_of_vectors);

  assert(num_vectors >= theoretical_minimum_of_vectors);
  assert(num_vectors < tr->mxtips);


  v->numVectors = num_vectors; /* use minimum bound theoretical */

  /* init vectors tracking */

  v->iVector         = (int *) rax_malloc((size_t)num_vectors * sizeof(int));
  v->unpinnable      = (pllBoolean *) rax_malloc((size_t)num_vectors * sizeof(pllBoolean));

  for(i = 0; i < num_vectors; i++)
  {
    v->iVector[i]         = PLL_SLOT_UNUSED;
    v->unpinnable[i]      = PLL_FALSE;
  }

  v->iNode      = (int *) rax_malloc((size_t)num_inner_nodes * sizeof(int));
  v->stlen      = (int *) rax_malloc((size_t)num_inner_nodes * sizeof(int));

  for(i = 0; i < num_inner_nodes; i++)
  {
    v->iNode[i] = PLL_NODE_UNPINNED;
    v->stlen[i] = PLL_INNER_NODE_INIT_STLEN;
  }

  v->allSlotsBusy = PLL_FALSE;

  /* init nodes tracking */

  v->maxVectorsUsed = 0;
  tr->rvec = v;
}

/** @brief Find the slot id with the minimum cost to be recomputed.
 *  
 *  The minum cost is defined as the minimum subtree size. In general, the closer a vector is to the tips, 
 *  the less recomputations are required to re-establish its likelihood entries
 *
 *  @todo remove _DEBUG_RECOMPUTATION code
 *  
 *  @param v
 *
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
static int findUnpinnableSlotByCost(recompVectors *v, int mxtips)
{
  int 
    i, 
    slot, 
    cheapest_slot = -1, 
    min_cost = mxtips * 2; /* more expensive than the most expensive*/
#ifdef _DEBUG_RECOMPUTATION 
  double straTime = gettime();
#endif 


  for(i = 0; i < mxtips - 2; i++)
  {
    slot = v->iNode[i];
    if(slot != PLL_NODE_UNPINNED)
    {
      assert(slot >= 0 && slot < v->numVectors);

      if(v->unpinnable[slot])
      {
        assert(v->stlen[i] > 0);

        if(v->stlen[i] < min_cost)
        {
          min_cost = v->stlen[i];
          cheapest_slot = slot;
          /* if the slot costs 2 you can break cause there is nothing cheaper to recompute */
          if(min_cost == 2)
            break;
        }
      }
    }
  }
  assert(min_cost < mxtips * 2 && min_cost >= 2);
  assert(cheapest_slot >= 0);
  return cheapest_slot;
}

static void unpinAtomicSlot(recompVectors *v, int slot, int mxtips)
{
  int 
    nodenum = v->iVector[slot];

  v->iVector[slot] = PLL_SLOT_UNUSED;

  if(nodenum != PLL_SLOT_UNUSED)  
    v->iNode[nodenum - mxtips - 1] = PLL_NODE_UNPINNED; 
}

/** @brief Finds the cheapest slot and unpins it
 *
 */
static int findUnpinnableSlot(recompVectors *v, int mxtips)
{
  int     
    slot_unpinned = findUnpinnableSlotByCost(v, mxtips);

  assert(slot_unpinned >= 0);
  assert(v->unpinnable[slot_unpinned]);

  unpinAtomicSlot(v, slot_unpinned, mxtips);

  return slot_unpinned;
}

/** @brief Finds a free slot 
 * 
 *  If all slots are occupied, it will find the cheapest slot and unpin it
 *
 */
static int findFreeSlot(recompVectors *v, int mxtips)
{
  int 
    slotno = -1, 
           i;

  assert(v->allSlotsBusy == PLL_FALSE);

  for(i = 0; i < v->numVectors; i++)
  {
    if(v->iVector[i] == PLL_SLOT_UNUSED)
    {
      slotno = i;
      break;
    } 
  }

  if(slotno == -1)
  {
    v->allSlotsBusy = PLL_TRUE;
    slotno = findUnpinnableSlot(v, mxtips);
  }

  return slotno;
}


/** @brief Pins node \a nodenum to slot \a slot
 *  
 *  The slot is initialized as non-unpinnable (ensures that the contents of the vector will not be overwritten)
 *
 *  @param nodenum
 *    node id
 *
 *  @param slot
 *    slot id 
 *    
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
static void pinAtomicNode(recompVectors *v, int nodenum, int slot, int mxtips)
{
  v->iVector[slot] = nodenum;
  v->iNode[nodenum - mxtips - 1] = slot;
  v->unpinnable[slot] = PLL_FALSE;
}

static int pinNode(recompVectors *rvec, int nodenum, int mxtips)
{
  int 
    slot;

  assert(!isNodePinned(rvec, nodenum, mxtips));

  if(rvec->allSlotsBusy)
    slot = findUnpinnableSlot(rvec, mxtips);
  else
    slot = findFreeSlot(rvec, mxtips);

  assert(slot >= 0);

  pinAtomicNode(rvec, nodenum, slot, mxtips);

  if(slot > rvec->maxVectorsUsed)
    rvec->maxVectorsUsed = slot;

  assert(slot == rvec->iNode[nodenum - mxtips - 1]);

  return slot;
}

/** @brief Marks node \a nodenum as unpinnable
 *  
 *  The slot holding the node \a nodenum is added to the pool of slot candidates that can be overwritten.
 *
 *  @param v
 *    Recomputation info
 *    
 *  @param nodenum
 *    node id
 *    
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
void unpinNode(recompVectors *v, int nodenum, int mxtips)
{
  if(nodenum <= mxtips)
    return;
  else
  {
    int 
      slot = -1;

    assert(nodenum > mxtips);
    slot = v->iNode[nodenum-mxtips-1];
    assert(slot >= 0 && slot < v->numVectors); 

    if(slot >= 0 && slot < v->numVectors)
      v->unpinnable[slot] = PLL_TRUE;
  }
}


/** @brief Get a pinned slot \a slot that holds the likelihood vector for inner node \a nodenum
 *  
 *  If node \a node nodenum is not pinned to any slot yet, the minimum cost replacement strategy is used.
 *
 *  @param v
 *    Recomputation info
 *    
 *  @param nodenum
 *    node id
 *    
 *  @param slot
 *    slot id
 *
 *  @param mxtips
 *    Number of tips in the tree
 *
 */
pllBoolean getxVector(recompVectors *rvec, int nodenum, int *slot, int mxtips)
{
  pllBoolean 
    slotNeedsRecomp = PLL_FALSE;

  *slot = rvec->iNode[nodenum - mxtips - 1];

  if(*slot == PLL_NODE_UNPINNED)
  {
    *slot = pinNode(rvec, nodenum, mxtips); /* now we will run the replacement strategy */
    slotNeedsRecomp = PLL_TRUE;
  }

  assert(*slot >= 0 && *slot < rvec->numVectors);

  rvec->unpinnable[*slot] = PLL_FALSE;

  return slotNeedsRecomp;
}


#ifdef _DEBUG_RECOMPUTATION

static int subtreeSize(nodeptr p, int maxTips)
{
  if(isTip(p->number, maxTips))
    return 1;
  else   
    return (subtreeSize(p->next->back, maxTips) + subtreeSize(p->next->next->back, maxTips));
}

#endif

/** @brief Annotes unoriented tree nodes \a tr with their subtree size 
 *  
 *  This function recursively updates the subtree size of each inner node.
 *  @note The subtree size of node \a p->number is the number of nodes included in the subtree where node record \a p is the virtual root. 
 *
 *  @param p
 *    Pointer to node 
 *    
 *  @param maxTips
 *    Number of tips in the tree
 *
 *  @param rvec 
 *    Recomputation info
 *    
 *  @param count
 *    Number of visited nodes 
 */
void computeTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec, int *count) 
{
  if(isTip(p->number, maxTips))
    return;
  else
  {          
    nodeptr 
      q = p->next->back,
        r = p->next->next->back;

    *count += 1;
    /* set xnode info at this point */     

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))  
    {
      rvec->stlen[p->number - maxTips - 1] = 2;	

#ifdef _DEBUG_RECOMPUTATION
      assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
    }
    else
    {
      if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
      {	     
        nodeptr 
          tmp;

        if(isTip(r->number, maxTips))
        {
          tmp = r;
          r = q;
          q = tmp;
        }

        if(!r->x)
          computeTraversalInfoStlen(r, maxTips, rvec, count);

        rvec->stlen[p->number - maxTips - 1] = rvec->stlen[r->number - maxTips - 1] + 1;

#ifdef _DEBUG_RECOMPUTATION	      
        assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
      }
      else
      {		 
        if(!r->x)
          computeTraversalInfoStlen(r, maxTips, rvec, count);
        if(!q->x)
          computeTraversalInfoStlen(q, maxTips, rvec, count); 

        rvec->stlen[p->number - maxTips - 1] = rvec->stlen[q->number - maxTips - 1] + rvec->stlen[r->number - maxTips - 1];	

#ifdef _DEBUG_RECOMPUTATION
        assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
      }
    }
  }
}




/* pre-compute the node stlens (this needs to be known prior to running the strategy) */
/** @brief Annotes all tree nodes \a tr with their subtree size 
 *  
 *  Similar to \a computeTraversalInfoStlen, but does a full traversal ignoring orientation.
 *  The minum cost is defined as the minimum subtree size. In general, the closer a vector is to the tips, 
 *  the less recomputations are required to re-establish its likelihood entries
 *
 *  @param p
 *    Pointer to node 
 *    
 *  @param maxTips
 *    Number of tips in the tree
 *
 *  @param rvec 
 *    Recomputation info
 */
void computeFullTraversalInfoStlen(nodeptr p, int maxTips, recompVectors *rvec) 
{
  if(isTip(p->number, maxTips))
    return;
  else
  {    
    nodeptr 
      q = p->next->back,
        r = p->next->next->back;     

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
    {	  
      rvec->stlen[p->number - maxTips - 1] = 2;

#ifdef _DEBUG_RECOMPUTATION
      assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
    }
    else
    {	    
      if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
      {	  	      
        nodeptr 
          tmp;

        if(isTip(r->number, maxTips))
        {
          tmp = r;
          r = q;
          q = tmp;
        }

        computeFullTraversalInfoStlen(r, maxTips, rvec);

        rvec->stlen[p->number - maxTips - 1] = rvec->stlen[r->number - maxTips - 1] + 1;	   

#ifdef _DEBUG_RECOMPUTATION
        assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
      }
      else
      {	    	     	      
        computeFullTraversalInfoStlen(r, maxTips, rvec);
        computeFullTraversalInfoStlen(q, maxTips, rvec); 

        rvec->stlen[p->number - maxTips - 1] = rvec->stlen[q->number - maxTips - 1] + rvec->stlen[r->number - maxTips - 1];
#ifdef _DEBUG_RECOMPUTATION
        assert(rvec->stlen[p->number - maxTips - 1] == subtreeSize(p, maxTips));
#endif
      }
    }
  }
}


#ifdef _DEBUG_RECOMPUTATION

void allocTraversalCounter(pllInstance *tr)
{
  traversalCounter 
    *tc;

  int 
    k;

  tc = (traversalCounter *)rax_malloc(sizeof(traversalCounter));

  tc->travlenFreq = (unsigned int *)rax_malloc(tr->mxtips * sizeof(int));

  for(k = 0; k < tr->mxtips; k++)
    tc->travlenFreq[k] = 0;

  tc->tt = 0;
  tc->ti = 0;
  tc->ii = 0;
  tc->numTraversals = 0;
  tr->travCounter = tc;
}

/* recomp */
/* code to track traversal descriptor stats */

void countTraversal(pllInstance *tr)
{
  traversalInfo 
    *ti   = tr->td[0].ti;
  int i;
  traversalCounter *tc = tr->travCounter; 
  tc->numTraversals += 1;

  /*
  printBothOpen("trav #%d(%d):",tc->numTraversals, tr->td[0].count);
  */

  for(i = 1; i < tr->td[0].count; i++)
  {
    traversalInfo *tInfo = &ti[i];

    /* 
       printBothOpen(" %d q%d r%d |",  tInfo->pNumber, tInfo->qNumber, tInfo->rNumber);
       printBothOpen("%d",  tInfo->pNumber);
       */
    switch(tInfo->tipCase)
    {
      case PLL_TIP_TIP: 
        tc->tt++; 
        /* printBothOpen("T"); */
        break;		  
      case PLL_TIP_INNER: 
        tc->ti++; 
        /* printBothOpen("M"); */
        break;		  

      case PLL_INNER_INNER: 
        tc->ii++; 
        /* printBothOpen("I"); */
        break;		  
      default: 
        assert(0);
    }
    /* printBothOpen(" "); */
  }
  /* printBothOpen(" so far T %d, M %d, I %d \n", tc->tt, tc->ti,tc->ii); */
  tc->travlenFreq[tr->td[0].count] += 1;
}


/*
void printTraversalInfo(pllInstance *tr)
{
  int 
    k, 
    total_steps = 0;

  printBothOpen("Traversals : %d \n", tr->travCounter->numTraversals);
  printBothOpen("Traversals tt: %d \n", tr->travCounter->tt);
  printBothOpen("Traversals ti: %d \n", tr->travCounter->ti);
  printBothOpen("Traversals ii: %d \n", tr->travCounter->ii);
  printBothOpen("all: %d \n", tr->travCounter->tt + tr->travCounter->ii + tr->travCounter->ti);
  printBothOpen("Traversals len freq  : \n");
  
  for(k = 0; k < tr->mxtips; k++)
  {
    total_steps += tr->travCounter->travlenFreq[k] * (k - 1);
    if(tr->travCounter->travlenFreq[k] > 0)
      printBothOpen("len %d : %d\n", k, tr->travCounter->travlenFreq[k]);
  }
  printBothOpen("all steps: %d \n", total_steps);
}
*/
/*end code to track traversal descriptor stats */
/* E recomp */

/*
void printVector(double *vector, int len, char *name)
{ 
  int i;
  printBothOpen("LHVECTOR %s :", name);
  for(i=0; i < len; i++)
  {
    printBothOpen("%.2f ", vector[i]);
    if(i>10)
    {
      printBothOpen("...");
      break; 
    }
  } 
  printBothOpen("\n");
} 
*/

#endif

