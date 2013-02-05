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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"


static boolean treeNeedString(const char *fp, char c1, int *position)
{
  char 
    c2 = fp[(*position)++];
  
  if(c2 == c1)  
    return TRUE;
  else  
    {   
      int 
	lower = MAX(0, *position - 20),
	upper = *position + 20;
      
      printf("Tree Parsing ERROR: Expecting '%c', found: '%c'\n", c1, c2); 
      printf("Context: \n");
      
      while(lower < upper && fp[lower])
	printf("%c", fp[lower++]);
      
      printf("\n");

      return FALSE;
  }
} 



static boolean treeLabelEndString (char ch)
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
      return TRUE;
    default:
      break;
    }
  
  return FALSE;
} 

static boolean  treeGetLabelString (const char *fp, char *lblPtr, int maxlen, int *position)
{
  char 
    ch;
  
  boolean  
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

static boolean  treeFlushLabelString(const char *fp, int *position)
{ 
  return  treeGetLabelString(fp, (char *) NULL, (int) 0, position);
} 


static boolean treeProcessLengthString (const char *fp, double *dptr, int *position)
{ 
  (*position)++;
  
  if(sscanf(&fp[*position], "%lf", dptr) != 1) 
    {
      printf("ERROR: treeProcessLength: Problem reading branch length\n");     
      assert(0);
    }

  while(fp[*position] != ',' && fp[*position] != ')' && fp[*position] != ';')
    *position = *position + 1;
  
  return  TRUE;
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

static int treeFindTipByLabelString(char  *str, tree *tr)                    
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

static int treeFindTipNameString (const char *fp, tree *tr, int *position)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabelString(fp, str, nmlngth+2, position))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   
  return  n;
} 

static boolean addElementLenString(const char *fp, tree *tr, nodeptr p, int *position)
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
	      return FALSE;
	    }
	  else 
	    {	   
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (!addElementLenString(fp, tr, q->next, position))        
	return FALSE;
      if (!treeNeedString(fp, ',', position))             
	return FALSE;
      if (!addElementLenString(fp, tr, q->next->next, position))  
	return FALSE;
      if (!treeNeedString(fp, ')', position))             
	return FALSE;
      
     
      treeFlushLabelString(fp, position);
    }
  else 
    {   
      (*position)--;
     
      if ((n = treeFindTipNameString(fp, tr, position)) <= 0)          
	return FALSE;
      q = tr->nodep[n];
      
      if (tr->start->number > n)  
	tr->start = q;
      (tr->ntips)++;
    }
  
     
  fres = treeFlushLenString(fp, position);
  if(!fres) 
    return FALSE;
  
  hookupDefault(p, q, tr->numBranches);

  return TRUE;          
}




void treeReadTopologyString(char *treeString, tree *tr)
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
  tr->rooted      = FALSE;      
  
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
