#include "mem_alloc.h"

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "pll.h"


 
  
/** @brief Reorder nodes in PLL tree

    Re-order the internal nodes of the tree of PLL instance \a tr in a preorder
    traversal such that they start from \a p
    
    @param tr
      PLL instance

    @param np
      Array of node pointers

    @param p
      Node from where the preorder traversal should start

    @param count

    @todo
      why not insert a break in the for loop when the node is found?

*/
static void reorderNodes(pllInstance *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if(isTip(p->number, tr->mxtips))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	} 
      
      assert(found != 0);
     
      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}

void nodeRectifier(pllInstance *tr)
{
  nodeptr *np = (nodeptr *)rax_malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = PLL_FALSE;

  /* TODO why is tr->rooted set to PLL_FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 

 
  rax_free(np);
}

nodeptr findAnyTip(nodeptr p, int numsp)
{ 
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
} 

