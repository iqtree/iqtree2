#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "axml.h"


static void insertTaxon (nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, numBranches);
  hookupDefault(p->next->next, r, numBranches); 
} 

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTreeRandom (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  
    p, 
    s;
  
  int  
    i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  
  s = buildNewTip(tr, tr->nodep[ir]);
  
  insertTaxon(s, p, tr->numBranches);
}

static int randomInt(int n, tree *tr)
{
  int 
    res = (int)((double)(n) * randum(&tr->randomNumberSeed));

  assert(res >= 0 && res < n);
  
  return res;
}

static void makePermutation(int *perm, int n, tree *tr)
{    
  int  
    i, 
    j, 
    k;    

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {    
      k =  randomInt(n + 1 - i, tr); /*(int)((double)(n + 1 - i) * randum(&tr->randomNumberSeed));*/

      assert(i + k <= n);
      
      j        = perm[i];
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

static int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp)
{
  if(isTip(p->number, numsp))
    return 0;
  else
    {
      branches[*counter] = p->next;
      branches[*counter + 1] = p->next->next;
      
      *counter = *counter + 2;
      
      return ((2 + markBranches(branches, p->next->back, counter, numsp) + 
	       markBranches(branches, p->next->next->back, counter, numsp)));
    }
}



void makeRandomTree(tree *tr)
{  
  nodeptr 
    p, 
    f, 
    randomBranch,
    *branches = (nodeptr *)malloc(sizeof(nodeptr) * (2 * tr->mxtips));    
  
  int 
    nextsp, 
    *perm = (int *)malloc((tr->mxtips + 1) * sizeof(int)), 
    branchCounter;                      
  
  makePermutation(perm, tr->mxtips, tr);              
  
  tr->ntips = 0;       	       
  tr->nextnode = tr->mxtips + 1;    
  
  buildSimpleTreeRandom(tr, perm[1], perm[2], perm[3]);
  
  while(tr->ntips < tr->mxtips) 
    {	             
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];            
      
      buildNewTip(tr, p);  	
      
      f = findAnyTip(tr->start, tr->mxtips);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(branches, f, &branchCounter, tr->mxtips);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[randomInt(branchCounter, tr)];
      
      insertTaxon(p->back, randomBranch, tr->numBranches);      
    }
  
  free(perm);            
  free(branches);
}

