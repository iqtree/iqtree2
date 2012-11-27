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
#include "axml.h"


 
  
static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
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

void nodeRectifier(tree *tr)
{
  nodeptr *np = (nodeptr *)malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = FALSE;

  /* TODO why is tr->rooted set to FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 

 
  free(np);
}

nodeptr findAnyTip(nodeptr p, int numsp)
{ 
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
} 

