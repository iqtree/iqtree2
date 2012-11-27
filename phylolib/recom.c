#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "axml.h"

void protectNode(recompVectors *rvec, int nodenum, int mxtips)
{
  /* If a node is available we dont need to recompute it, but we neet to maker sure it is not unpinned while buildding the rest of the traversal descriptor, i.e. unpinnable must be false at this point, it will automatically be set to true, after the *counter post-order instructions have been executed */
  /* This might be a but in RAxML-light, Omitting this code will likely still work as long as num_allocated_nodes >> log n, but wrong inner vectors will be used at the wrong moment of newviewIterative, careful! */

  int slot;
  slot = rvec->iNode[nodenum - mxtips - 1];
  assert(slot != NODE_UNPINNED);
  assert(rvec->iVector[slot] == nodenum);

  if(rvec->unpinnable[slot])
    rvec->unpinnable[slot] = FALSE;
}

static boolean isNodePinned(recompVectors *rvec, int nodenum, int mxtips)
{
  assert(nodenum > mxtips);

  if(rvec->iNode[nodenum - mxtips - 1] == NODE_UNPINNED)
    return FALSE;
  else
    return TRUE;
}

boolean needsRecomp(boolean recompute, recompVectors *rvec, nodeptr p, int mxtips)
{ 
  /* A node needs reomputation if one of the following holds:
     1. It is not oriented (p->x == 0) 
     2. We are applying recomputations and p is not available on RAM
   */
  if((!p->x) || (recompute && !isNodePinned(rvec, p->number, mxtips)))
    return TRUE;
  else
    return FALSE;
}



void allocRecompVectorsInfo(tree *tr)
{
  recompVectors 
    *v = (recompVectors *) malloc(sizeof(recompVectors));

  int 
    num_inner_nodes = tr->mxtips - 2,
                    num_vectors, 
                    i;

  assert(tr->vectorRecomFraction > MIN_RECOM_FRACTION);
  assert(tr->vectorRecomFraction < MAX_RECOM_FRACTION);

  num_vectors = (int) (1 + tr->vectorRecomFraction * (float)num_inner_nodes); 

  int theoretical_minimum_of_vectors = 3 + ((int)(log((double)tr->mxtips)/log(2.0)));
  printBothOpen("Try to use %d ancestral vectors, min required %d\n", num_vectors, theoretical_minimum_of_vectors);

  assert(num_vectors >= theoretical_minimum_of_vectors);
  assert(num_vectors < tr->mxtips);


  v->numVectors = num_vectors; /* use minimum bound theoretical */

  /* init vectors tracking */

  v->iVector         = (int *) malloc((size_t)num_vectors * sizeof(int));
  v->unpinnable      = (boolean *) malloc((size_t)num_vectors * sizeof(boolean));

  for(i = 0; i < num_vectors; i++)
  {
    v->iVector[i]         = SLOT_UNUSED;
    v->unpinnable[i]      = FALSE;
  }

  v->iNode      = (int *) malloc((size_t)num_inner_nodes * sizeof(int));
  v->stlen      = (int *) malloc((size_t)num_inner_nodes * sizeof(int));

  for(i = 0; i < num_inner_nodes; i++)
  {
    v->iNode[i] = NODE_UNPINNED;
    v->stlen[i] = INNER_NODE_INIT_STLEN;
  }

  v->allSlotsBusy = FALSE;

  /* init nodes tracking */

  v->maxVectorsUsed = 0;
  tr->rvec = v;
}


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

    if(slot != NODE_UNPINNED)
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
#ifdef _DEBUG_RECOMPUTATION 
  v->recomStraTime += (gettime() - straTime);
#endif 

  return cheapest_slot;
}

static void unpinAtomicSlot(recompVectors *v, int slot, int mxtips)
{
  int 
    nodenum = v->iVector[slot];

  v->iVector[slot] = SLOT_UNUSED;

  if(nodenum != SLOT_UNUSED)  
    v->iNode[nodenum - mxtips - 1] = NODE_UNPINNED; 
}

static int findUnpinnableSlot(recompVectors *v, int mxtips)
{
  int     
    slot_unpinned = findUnpinnableSlotByCost(v, mxtips);

  assert(slot_unpinned >= 0);
  assert(v->unpinnable[slot_unpinned]);

  unpinAtomicSlot(v, slot_unpinned, mxtips);

  return slot_unpinned;
}

static int findFreeSlot(recompVectors *v, int mxtips)
{
  int 
    slotno = -1, 
           i;

  assert(v->allSlotsBusy == FALSE);

  for(i = 0; i < v->numVectors; i++)
  {
    if(v->iVector[i] == SLOT_UNUSED)
    {
      slotno = i;
      break;
    } 
  }

  if(slotno == -1)
  {
    v->allSlotsBusy = TRUE;
    slotno = findUnpinnableSlot(v, mxtips);
  }

  return slotno;
}


static void pinAtomicNode(recompVectors *v, int nodenum, int slot, int mxtips)
{
  v->iVector[slot] = nodenum;
  v->iNode[nodenum - mxtips - 1] = slot;
  v->unpinnable[slot] = FALSE;
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
      v->unpinnable[slot] = TRUE;
  }
}

/*
#ifdef _DEBUG_RECOMPUTATION
static void printSlots(recompVectors *rvec, int nodenum)
{
  int i;
  printBothOpen("n %d , current slots:", nodenum);
  for(i=0;i<rvec->numVectors;i++)
    printBothOpen(" %d (%s)", rvec->iVector[i], (rvec->unpinnable[i] ? " ": "x"));
  printBothOpen("\n");
}
#endif
*/


boolean getxVector(recompVectors *rvec, int nodenum, int *slot, int mxtips)
{
#ifdef _DEBUG_RECOMPUTATION
  double 
    tstart = gettime();
#endif

  boolean 
    slotNeedsRecomp = FALSE;

  *slot = rvec->iNode[nodenum - mxtips - 1];

  if(*slot == NODE_UNPINNED)
  {
    *slot = pinNode(rvec, nodenum, mxtips); /* now we will run the replacement strategy */
    slotNeedsRecomp = TRUE;
  }

  assert(*slot >= 0 && *slot < rvec->numVectors);

  rvec->unpinnable[*slot] = FALSE;
#ifdef _DEBUG_RECOMPUTATION
  rvec->pinTime += gettime() - tstart;
#endif

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

void allocTraversalCounter(tree *tr)
{
  traversalCounter 
    *tc;

  int 
    k;

  tc = (traversalCounter *)malloc(sizeof(traversalCounter));

  tc->travlenFreq = (unsigned int *)malloc(tr->mxtips * sizeof(int));

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

void countTraversal(tree *tr)
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
      case TIP_TIP: 
        tc->tt++; 
        /* printBothOpen("T"); */
        break;		  
      case TIP_INNER: 
        tc->ti++; 
        /* printBothOpen("M"); */
        break;		  

      case INNER_INNER: 
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



void printTraversalInfo(tree *tr)
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


