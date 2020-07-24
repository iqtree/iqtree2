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
 * @file bipartitionList.c
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

/*
static pllBoolean treeNeedString(const char *fp, char c1, int *position)
{
  char 
    c2 = fp[(*position)++];
  
  if(c2 == c1)  
    return PLL_TRUE;
  else  
    {   
      int 
	lower = PLL_MAX(0, *position - 20),
	upper = *position + 20;
      
      printf("Tree Parsing ERROR: Expecting '%c', found: '%c'\n", c1, c2); 
      printf("Context: \n");
      
      while(lower < upper && fp[lower])
	printf("%c", fp[lower++]);
      
      printf("\n");

      return PLL_FALSE;
  }
} 


static pllBoolean treeLabelEndString (char ch)
{
  switch(ch) 
    {   
    case '\0':  
    case '\t':  
    case '\n':  
    case '\r': 
    case ' ':
    case ':':  
    case ',':   
    case '(':   
    case ')':  
    case ';':
      return PLL_TRUE;
    default:
      break;
    }
  
  return PLL_FALSE;
} 

static pllBoolean  treeGetLabelString (const char *fp, char *lblPtr, int maxlen, int *position)
{
  char 
    ch;
  
  pllBoolean  
    done, 
    lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *)NULL; 
  else 
    if(lblPtr == NULL) 
      maxlen = 0;

  ch = fp[(*position)++];
  
  done = treeLabelEndString(ch);

  lblfound = !done;  

  while(!done) 
    {      
      if(treeLabelEndString(ch)) 
	break;     

      if(--maxlen >= 0) 
	*lblPtr++ = ch;
      
      ch = fp[(*position)++];      
    }
  
  (*position)--; 

  if (lblPtr != NULL) 
    *lblPtr = '\0';

  return lblfound;
}

static pllBoolean  treeFlushLabelString(const char *fp, int *position)
{ 
  return  treeGetLabelString(fp, (char *) NULL, (int) 0, position);
} 

static pllBoolean treeProcessLengthString (const char *fp, double *dptr, int *position)
{ 
  (*position)++;
  
  if(sscanf(&fp[*position], "%lf", dptr) != 1) 
    {
      printf("ERROR: treeProcessLength: Problem reading branch length\n");     
      assert(0);
    }

  while(fp[*position] != ',' && fp[*position] != ')' && fp[*position] != ';')
    *position = *position + 1;
  
  return  PLL_TRUE;
}

static int treeFlushLenString (const char *fp, int *position)
{
  double  
    dummy;  
  
  char     
    ch;

  ch = fp[(*position)++];
 
  if(ch == ':') 
    {     
      if(!treeProcessLengthString(fp, &dummy, position)) 
	return 0;
      return 1;	  
    }
    
  (*position)--;

  return 1;
} 

static int treeFindTipByLabelString(char  *str, pllInstance *tr)                    
{
  int lookup = lookupWord(str, tr->nameHash);

  if(lookup > 0)
    {
      assert(! tr->nodep[lookup]->back);
      return lookup;
    }
  else
    { 
      printf("ERROR: Cannot find tree species: %s\n", str);
      return  0;
    }
}

static int treeFindTipNameString (const char *fp, pllInstance *tr, int *position)
{
  char    str[PLL_NMLNGTH + 2];
  int      n;

  if (treeGetLabelString (fp, str, PLL_NMLNGTH + 2, position))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   
  return  n;
} 

static pllBoolean addElementLenString(const char *fp, pllInstance *tr, nodeptr p, int *position)
{
  nodeptr  
    q;
  
  int      
    n, 
    fres;

  char 
    ch;
  
  if ((ch = fp[(*position)++]) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return PLL_FALSE;
	    }
	  else 
	    {	   
	      tr->rooted = PLL_TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (!addElementLenString(fp, tr, q->next, position))        
	return PLL_FALSE;
      if (!treeNeedString(fp, ',', position))             
	return PLL_FALSE;
      if (!addElementLenString(fp, tr, q->next->next, position))  
	return PLL_FALSE;
      if (!treeNeedString(fp, ')', position))             
	return PLL_FALSE;
      
     
      treeFlushLabelString(fp, position);
    }
  else 
    {   
      (*position)--;
     
      if ((n = treeFindTipNameString(fp, tr, position)) <= 0)          
	return PLL_FALSE;
      q = tr->nodep[n];
      
      if (tr->start->number > n)  
	tr->start = q;
      (tr->ntips)++;
    }
  
     
  fres = treeFlushLenString(fp, position);
  if(!fres) 
    return PLL_FALSE;
  
  hookupDefault(p, q);

  return PLL_TRUE;          
}



void treeReadTopologyString(char *treeString, pllInstance *tr)
{ 
  char 
    *fp = treeString;

  nodeptr  
    p;
  
  int
    position = 0, 
    i;
  
  char 
    ch;   
    

  for(i = 1; i <= tr->mxtips; i++)    
    tr->nodep[i]->back = (node *)NULL;      
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;           
    }
      
  tr->start       = tr->nodep[1];
  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;    
  tr->rooted      = PLL_FALSE;      
  
  p = tr->nodep[(tr->nextnode)++]; 
   
  assert(fp[position++] == '(');  
    
  if (! addElementLenString(fp, tr, p, &position))                 
    assert(0);
  
  if (! treeNeedString(fp, ',', &position))                
    assert(0);
   
  if (! addElementLenString(fp, tr, p->next, &position))           
    assert(0);

  if(!tr->rooted) 
    {
      if ((ch = fp[position++]) == ',') 
	{ 
	  if (! addElementLenString(fp, tr, p->next->next, &position)) 
	    assert(0);	 
	}
      else 
	assert(0);     
    }
  else
    assert(0);
        
  if (! treeNeedString(fp, ')', &position))                
    assert(0);

  treeFlushLabelString(fp, &position);
  
  if (!treeFlushLenString(fp, &position))                         
    assert(0);
  
  if (!treeNeedString(fp, ';', &position))       
    assert(0);
    
  if(tr->rooted)     
    assert(0);           
  else           
    tr->start = tr->nodep[1];   

  printf("Tree parsed\n");

} 
*/
