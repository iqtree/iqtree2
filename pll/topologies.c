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
 * @file topologies.c
 * @brief Miscellanous functions working with tree topology
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

static void saveTopolRELLRec(pllInstance *tr, nodeptr p, topolRELL *tpl, int *i, int numsp)
{
  int k;
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;      
      while(q != p)
	{	  
	  tpl->connect[*i].p = q;
	  tpl->connect[*i].q = q->back; 
	  
	  if(tr->grouped ||  tr->constrained)
	    {
	      tpl->connect[*i].cp = tr->constraintVector[q->number];
	      tpl->connect[*i].cq = tr->constraintVector[q->back->number]; 
	    }
	  
	  for(k = 0; k < PLL_NUM_BRANCHES; k++)
	    tpl->connect[*i].z[k] = q->z[k];
	  *i = *i + 1;

	  saveTopolRELLRec(tr, q->back, tpl, i, numsp);
	  q = q->next;
	}
    }
}

static void saveTopolRELL(pllInstance *tr, topolRELL *tpl)
{
  nodeptr p = tr->start;
  int k, i = 0;
      
  tpl->likelihood = tr->likelihood;
  tpl->start      = 1;
      
  tpl->connect[i].p = p;
  tpl->connect[i].q = p->back;
  
  if(tr->grouped ||  tr->constrained)
    {
      tpl->connect[i].cp = tr->constraintVector[p->number];
      tpl->connect[i].cq = tr->constraintVector[p->back->number]; 
    }

  for(k = 0; k < PLL_NUM_BRANCHES; k++)
    tpl->connect[i].z[k] = p->z[k];
  i++;
      
  saveTopolRELLRec(tr, p->back, tpl, &i, tr->mxtips);

  assert(i == 2 * tr->mxtips - 3);
}


static void restoreTopolRELL(pllInstance *tr, topolRELL *tpl, int numBranches)
{
  int i;
  
  for (i = 0; i < 2 * tr->mxtips - 3; i++) 
    {
      hookup(tpl->connect[i].p, tpl->connect[i].q, tpl->connect[i].z,  numBranches);
      tr->constraintVector[tpl->connect[i].p->number] = tpl->connect[i].cp;
      tr->constraintVector[tpl->connect[i].q->number] = tpl->connect[i].cq;
    }
  

  tr->likelihood = tpl->likelihood;
  tr->start      = tr->nodep[tpl->start];
  /* TODO */
}



/** @brief Initializes space as large as the tree
  *
  * @param rl
  *   RELL 
  *
  * @param tr
  *   PLL instance
  *
  * @param n
  *   Number of
  *
  * @todo
  *   Don't know what is this used for. Something with RELL?
  *
  */
void initTL(topolRELL_LIST *rl, pllInstance *tr, int n)
{
  int i;

  rl->max = n; 
  rl->t = (topolRELL **)rax_malloc(sizeof(topolRELL *) * n);

  for(i = 0; i < n; i++)
    {
      rl->t[i] = (topolRELL *)rax_malloc(sizeof(topolRELL));
      rl->t[i]->connect = (connectRELL *)rax_malloc((2 * tr->mxtips - 3) * sizeof(connectRELL));
      rl->t[i]->likelihood = PLL_UNLIKELY;     
    }
}

/** @brief Deallocate the space associated with this structure
  *
  * @paral rl
  *   This structure
  *
  * @todo
  *   fill the description
  */
void freeTL(topolRELL_LIST *rl)
{
  int i;
  for(i = 0; i < rl->max; i++)    
    {
      rax_free(rl->t[i]->connect);          
      rax_free(rl->t[i]);
    }
  rax_free(rl->t);
}


void restoreTL(topolRELL_LIST *rl, pllInstance *tr, int n, int numBranches)
{
  assert(n >= 0 && n < rl->max);    

  restoreTopolRELL(tr, rl->t[n], numBranches);
}



/** @brief Reset this structure
  *
  * Reset the likelihoods in this structure
  *
  * @param rl
  *   This structure
  *
  * @todo
  *   Complete this
  */
void resetTL(topolRELL_LIST *rl)
{
  int i;

  for(i = 0; i < rl->max; i++)    
    rl->t[i]->likelihood = PLL_UNLIKELY;          
}



/** @brief Save 
  *
  * Save this topology?
  *
  * @todo 
  *  Complete this
  */
void saveTL(topolRELL_LIST *rl, pllInstance *tr, int index)
{ 
  assert(index >= 0 && index < rl->max);    
    
  if(tr->likelihood > rl->t[index]->likelihood)        
    saveTopolRELL(tr, rl->t[index]); 
}


static void  *tipValPtr (nodeptr p)
{ 
  return  (void *) & p->number;
}


static int  cmpTipVal (void *v1, void *v2)
{
  int  i1, i2;
  
  i1 = *((int *) v1);
  i2 = *((int *) v2);
  return  (i1 < i2) ? -1 : ((i1 == i2) ? 0 : 1);
}


/*  These are the only routines that need to UNDERSTAND topologies */

/** @brief Allocate and initialize space for a tree topology
    
    Allocate and initialize a \a topol structure for a tree topology of
    \a maxtips tips

    @param
      Number of tips of topology

    @return
      Pointer to the allocated \a topol structure
*/
topol  *setupTopol (int maxtips)
{
  topol   *tpl;

  if (! (tpl = (topol *) rax_malloc(sizeof(topol))) || 
      ! (tpl->links = (connptr) rax_malloc((2*maxtips-3) * sizeof(pllConnect))))
    {
      printf("ERROR: Unable to get topology memory");
      tpl = (topol *) NULL;
    }
  else 
    {
      tpl->likelihood  = PLL_UNLIKELY;
      tpl->start       = (node *) NULL;
      tpl->nextlink    = 0;
      tpl->ntips       = 0;
      tpl->nextnode    = 0;    
      tpl->scrNum      = 0;     /* position in sorted list of scores */
      tpl->tplNum      = 0;     /* position in sorted list of trees */	      
    }
  
  return  tpl;
} 


/** @brief Deallocate the space occupied by a \a topol structure
    
    Deallocate the space occupied by a \a topol structure

    @param tpl
      The \a topol structure that is to be deallocated
*/
void freeTopol (topol *tpl)
{
  rax_free(tpl->links);
  rax_free(tpl);
} 


static int saveSubtree (nodeptr p, topol *tpl, int numsp, int numBranches)  
{
  connptr  r, r0;
  nodeptr  q, s;
  int      t, t0, t1, k;

  r0 = tpl->links;
  r = r0 + (tpl->nextlink)++;
  r->p = p;
  r->q = q = p->back;

  for(k = 0; k < numBranches; k++)
    r->z[k] = p->z[k];

  r->descend = 0;                     /* No children (yet) */

  if (isTip(q->number, numsp)) 
    {
      r->valptr = tipValPtr(q);         /* Assign value */
    }
  else 
    {                              /* Internal node, look at children */
      s = q->next;                      /* First child */
      do 
	{
	  t = saveSubtree(s, tpl, numsp, numBranches);        /* Generate child's subtree */

	  t0 = 0;                         /* Merge child into list */
	  t1 = r->descend;
	  while (t1 && (cmpTipVal(r0[t1].valptr, r0[t].valptr) < 0)) {
	    t0 = t1;
	    t1 = r0[t1].sibling;
          }
	  if (t0) r0[t0].sibling = t;  else  r->descend = t;
	  r0[t].sibling = t1;

	  s = s->next;                    /* Next child */
        } while (s != q);

      r->valptr = r0[r->descend].valptr;   /* Inherit first child's value */
      }                                 /* End of internal node processing */

  return  (r - r0);
}

/** @brief Get the node with the smallest tip value
    
    Recursively finds and returns the tip with the smallest value around a node
    \a p0, or returns \a p0 if it is a tip.

    @param p0
      Node around which to at which the recursion starts

    @param numsp
      Number of species (tips) in the tree

    @todo
      Why do we return p0 immediately if it is a tip? Perhaps one of the two other nodes,
      i.e. p0->next and p0->next->next, is a tip as well with a smaller number than p0.
*/
static nodeptr minSubtreeTip (nodeptr  p0, int numsp)
{ 
  nodeptr  minTip, p, testTip;

  if (isTip(p0->number, numsp)) 
    return p0;

  p = p0->next;

  minTip = minSubtreeTip(p->back, numsp);

  while ((p = p->next) != p0) 
    {
      testTip = minSubtreeTip(p->back, numsp);
      if (cmpTipVal(tipValPtr(testTip), tipValPtr(minTip)) < 0)
        minTip = testTip;
    }
  return minTip;
} 


/** @brief
*/
static nodeptr  minTreeTip (nodeptr  p, int numsp)
{
  nodeptr  minp, minpb;

  minp  = minSubtreeTip(p, numsp);
  minpb = minSubtreeTip(p->back, numsp);
  return (cmpTipVal(tipValPtr(minp), tipValPtr(minpb)) < 0 ? minp : minpb);
}

/** @brief Save the tree topology in a \a topol structure
    
    Save the current tree topology in \a topol structure \a tpl.

*/
void saveTree (pllInstance *tr, topol *tpl, int numBranches)
/*  Save a tree topology in a standard order so that first branches
 *  from a node contain lower value tips than do second branches from
 *  the node.  The root tip should have the lowest value of all.
 */
{
  connptr  r;  
  
  tpl->nextlink = 0;                             /* Reset link pointer */
  r = tpl->links + saveSubtree(minTreeTip(tr->start, tr->mxtips), tpl, tr->mxtips, numBranches);  /* Save tree */
  r->sibling = 0;
  
  tpl->likelihood = tr->likelihood;
  tpl->start      = tr->start;
  tpl->ntips      = tr->ntips;
  tpl->nextnode   = tr->nextnode;    
  
} /* saveTree */


/* @brief Transform tree to a given topology and evaluate likelihood

   Transform our current tree topology to the one stored in \a tpl and
   evaluates the likelihood

   @param tr
     PLL instance

   @param pr
     List of partitions

   @return
     \b PLL_TRUE

   @todo
     Remove the return value, unnecessary

*/
pllBoolean restoreTree (topol *tpl, pllInstance *tr, partitionList *pr)
{ 
  connptr  r;
  nodeptr  p, p0;    
  int  i;

  /* first of all set all backs to NULL so that tips do not point anywhere */
  for (i = 1; i <= 2*(tr->mxtips) - 2; i++) 
    {  
      /* Uses p = p->next at tip */
      p0 = p = tr->nodep[i];
      do 
	{
	  p->back = (nodeptr) NULL;
	  p = p->next;
	} 
      while (p != p0);
    }

  /*  Copy connections from topology */

  /* then connect the nodes together */
  for (r = tpl->links, i = 0; i < tpl->nextlink; r++, i++)     
    hookup(r->p, r->q, r->z, pr->perGeneBranchLengths?pr->numberOfPartitions:1);

  tr->likelihood = tpl->likelihood;
  tr->start      = tpl->start;
  tr->ntips      = tpl->ntips;
  
  tr->nextnode   = tpl->nextnode;    

  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
  return PLL_TRUE;
}



/** @brief Initialize a list of best trees
    
    Initialize a list that will contain the best \a newkeep tree topologies,
    i.e. the ones that yield the best likelihood. Inside the list initialize
    space for \a newkeep + 1 topologies of \a numsp tips. The additional
    topology is the starting one

    @param bt
      Pointer to \a bestlist to be initialized

    @param newkeep
      Number of new topologies to keep

    @param numsp
      Number of species (tips)

    @return
      number of tree topology slots in the list (minus the starting one)

    @todo
      Is there a reason that this function is so complex? Many of the checks
      are unnecessary as the function is called only at two places in the
      code with newkeep=1 and newkeep=20
*/
int initBestTree (bestlist *bt, int newkeep, int numsp)
{ /* initBestTree */
  int  i;

  bt->nkeep = 0;

  if (bt->ninit <= 0) 
    {
      if (! (bt->start = setupTopol(numsp)))  return  0;
      bt->ninit    = -1;
      bt->nvalid   = 0;
      bt->numtrees = 0;
      bt->best     = PLL_UNLIKELY;
      bt->improved = PLL_FALSE;
      bt->byScore  = (topol **) rax_malloc((newkeep + 1) * sizeof(topol *));
      bt->byTopol  = (topol **) rax_malloc((newkeep + 1) * sizeof(topol *));
      if (! bt->byScore || ! bt->byTopol) {
        printf( "initBestTree: malloc failure\n");
        return 0;
      }
    }
  else if (PLL_ABS(newkeep) > bt->ninit) {
    if (newkeep <  0) newkeep = -(bt->ninit);
    else newkeep = bt->ninit;
  }

  if (newkeep < 1) {    /*  Use negative newkeep to clear list  */
    newkeep = -newkeep;
    if (newkeep < 1) newkeep = 1;
    bt->nvalid = 0;
    bt->best = PLL_UNLIKELY;
  }
  
  if (bt->nvalid >= newkeep) {
    bt->nvalid = newkeep;
    bt->worst = bt->byScore[newkeep]->likelihood;
  }
  else 
    {
      bt->worst = PLL_UNLIKELY;
    }
  
  for (i = bt->ninit + 1; i <= newkeep; i++) 
    {    
      if (! (bt->byScore[i] = setupTopol(numsp)))  break;
      bt->byTopol[i] = bt->byScore[i];
      bt->ninit = i;
    }
  
  return  (bt->nkeep = PLL_MIN(newkeep, bt->ninit));
} /* initBestTree */



void resetBestTree (bestlist *bt)
{ /* resetBestTree */
  bt->best     = PLL_UNLIKELY;
  bt->worst    = PLL_UNLIKELY;
  bt->nvalid   = 0;
  bt->improved = PLL_FALSE;
} /* resetBestTree */


pllBoolean  freeBestTree(bestlist *bt)
{ /* freeBestTree */
  while (bt->ninit >= 0)  freeTopol(bt->byScore[(bt->ninit)--]);
    
  /* VALGRIND */

  rax_free(bt->byScore);
  rax_free(bt->byTopol);

  /* VALGRIND END */

  freeTopol(bt->start);
  return PLL_TRUE;
} /* freeBestTree */


/*  Compare two trees, assuming that each is in standard order.  Return
 *  -1 if first preceeds second, 0 if they are identical, or +1 if first
 *  follows second in standard order.  Lower number tips preceed higher
 *  number tips.  A tip preceeds a corresponding internal node.  Internal
 *  nodes are ranked by their lowest number tip.
 */

static int  cmpSubtopol (connptr p10, connptr p1, connptr p20, connptr p2)
{
  connptr  p1d, p2d;
  int  cmp;
  
  if (! p1->descend && ! p2->descend)          /* Two tips */
    return cmpTipVal(p1->valptr, p2->valptr);
  
  if (! p1->descend) return -1;                /* p1 = tip, p2 = node */
  if (! p2->descend) return  1;                /* p2 = tip, p1 = node */
  
  p1d = p10 + p1->descend;
  p2d = p20 + p2->descend;
  while (1) {                                  /* Two nodes */
    if ((cmp = cmpSubtopol(p10, p1d, p20, p2d)))  return cmp; /* Subtrees */
    if (! p1d->sibling && ! p2d->sibling)  return  0; /* Lists done */
    if (! p1d->sibling) return -1;             /* One done, other not */
    if (! p2d->sibling) return  1;             /* One done, other not */
    p1d = p10 + p1d->sibling;                  /* Neither done */
    p2d = p20 + p2d->sibling;
  }
}



static int  cmpTopol (void *tpl1, void *tpl2)
{ 
  connptr  r1, r2;
  int      cmp;    
  
  r1 = ((topol *) tpl1)->links;
  r2 = ((topol *) tpl2)->links;
  cmp = cmpTipVal(tipValPtr(r1->p), tipValPtr(r2->p));
  if (cmp)      	
    return cmp;     
  return  cmpSubtopol(r1, r1, r2, r2);
} 



static int  cmpTplScore (void *tpl1, void *tpl2)
{ 
  double  l1, l2;
  
  l1 = ((topol *) tpl1)->likelihood;
  l2 = ((topol *) tpl2)->likelihood;
  return  (l1 > l2) ? -1 : ((l1 == l2) ? 0 : 1);
}



/*  Find an item in a sorted list of n items.  If the item is in the list,
 *  return its index.  If it is not in the list, return the negative of the
 *  position into which it should be inserted.
 */

static int  findInList (void *item, void *list[], int n, int (* cmpFunc)(void *, void *))
{
  int  mid, hi, lo, cmp = 0;
  
  if (n < 1) return  -1;                    /*  No match; first index  */
  
  lo = 1;
  mid = 0;
  hi = n;
  while (lo < hi) {
    mid = (lo + hi) >> 1;
    cmp = (* cmpFunc)(item, list[mid-1]);
    if (cmp) {
      if (cmp < 0) hi = mid;
      else lo = mid + 1;
    }
    else  return  mid;                        /*  Exact match  */
  }
  
  if (lo != mid) {
    cmp = (* cmpFunc)(item, list[lo-1]);
    if (cmp == 0) return lo;
  }
  if (cmp > 0) lo++;                         /*  Result of step = 0 test  */
  return  -lo;
} 



static int  findTreeInList (bestlist *bt, pllInstance *tr, int numBranches)
{
  topol  *tpl;
  
  tpl = bt->byScore[0];
  saveTree(tr, tpl, numBranches);
  return  findInList((void *) tpl, (void **) (& (bt->byTopol[1])),
		     bt->nvalid, cmpTopol);
} 


/** @brief Save the current tree in the \a bestlist structure
    
    Save the current tree topology in \a bestlist structure \a bt.

    @param tr
      The PLL instance
    
    @param bt
      The \a bestlist structure
    
    @param numBranches
      Number of branches u

    @return
      it is never used

    @todo
      What to do with the return value? Should we simplify the code?
*/
int  saveBestTree (bestlist *bt, pllInstance *tr, int numBranches)
{    
  topol  *tpl, *reuse;
  int  tplNum, scrNum, reuseScrNum, reuseTplNum, i, oldValid, newValid;
  
  tplNum = findTreeInList(bt, tr, numBranches);
  tpl = bt->byScore[0];
  oldValid = newValid = bt->nvalid;
  
  if (tplNum > 0) {                      /* Topology is in list  */
    reuse = bt->byTopol[tplNum];         /* Matching topol  */
    reuseScrNum = reuse->scrNum;
    reuseTplNum = reuse->tplNum;
  }
  /* Good enough to keep? */
  else if (tr->likelihood < bt->worst)  return 0;
  
  else {                                 /* Topology is not in list */
    tplNum = -tplNum;                    /* Add to list (not replace) */
    if (newValid < bt->nkeep) bt->nvalid = ++newValid;
    reuseScrNum = newValid;              /* Take worst tree */
    reuse = bt->byScore[reuseScrNum];
    reuseTplNum = (newValid > oldValid) ? newValid : reuse->tplNum;
    if (tr->likelihood > bt->start->likelihood) bt->improved = PLL_TRUE;
  }
  
  scrNum = findInList((void *) tpl, (void **) (& (bt->byScore[1])),
		      oldValid, cmpTplScore);
  scrNum = PLL_ABS(scrNum);
  
  if (scrNum < reuseScrNum)
    for (i = reuseScrNum; i > scrNum; i--)
      (bt->byScore[i] = bt->byScore[i-1])->scrNum = i;
  
  else if (scrNum > reuseScrNum) {
    scrNum--;
    for (i = reuseScrNum; i < scrNum; i++)
      (bt->byScore[i] = bt->byScore[i+1])->scrNum = i;
  }
  
  if (tplNum < reuseTplNum)
    for (i = reuseTplNum; i > tplNum; i--)
      (bt->byTopol[i] = bt->byTopol[i-1])->tplNum = i;
  
  else if (tplNum > reuseTplNum) {
    tplNum--;
    for (i = reuseTplNum; i < tplNum; i++)
      (bt->byTopol[i] = bt->byTopol[i+1])->tplNum = i;
  }
  
  
  
  tpl->scrNum = scrNum;
  tpl->tplNum = tplNum;
  bt->byTopol[tplNum] = bt->byScore[scrNum] = tpl;
  bt->byScore[0] = reuse;
  
  if (scrNum == 1)  bt->best = tr->likelihood;
  if (newValid == bt->nkeep) bt->worst = bt->byScore[newValid]->likelihood;
  
  return  scrNum;
} 


/** @brief Restore the best tree from \a bestlist structure
    
    Restore the \a rank-th best tree from the \a bestlist structure \a bt.

    @param bt
      The \a bestlist structure containing the stored best trees

    @param rank
      The rank (by score) of the tree we want to retrieve

    @param tr
      PLL instance

    @param pr
      List of partitions

    @return
      Index (rank) of restored topology in \a bestlist
*/
int  recallBestTree (bestlist *bt, int rank, pllInstance *tr, partitionList *pr)
{ 
  if (rank < 1)  rank = 1;
  if (rank > bt->nvalid)  rank = bt->nvalid;
  if (rank > 0)  if (! restoreTree(bt->byScore[rank], tr, pr)) return PLL_FALSE;
  return  rank;
}




