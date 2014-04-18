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

#include "mem_alloc.h"

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
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

extern char infoFileName[1024];
extern char tree_file[1024];
extern char *likelihood_key;
extern char *ntaxa_key;
extern char *smoothed_key;
extern int partCount;
extern double masterTime;





stringHashtable *initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)rax_malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)rax_calloc((size_t)tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}


static hashNumberType  hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}

 

void addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)rax_malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)rax_malloc(((size_t)strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}

int lookupWord(char *s, stringHashtable *h)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return p->nodeNumber;	  	
    }

  return -1;
}


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


  
static char *TreeInner2StringREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, boolean printBranchLengths, boolean printNames,
			    boolean printLikelihood, boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport, boolean printInnerNodes)
{
  /* TODOFER simplify this, should be used just to print inner nodes for testing */
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
    treestr = TreeInner2StringREC(treestr, tr, pr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
        finalPrint, perGene, branchLabelSupport, printSHSupport, printInnerNodes);
    *treestr++ = ',';
    treestr = TreeInner2StringREC(treestr, tr, pr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
        finalPrint, perGene, branchLabelSupport, printSHSupport, printInnerNodes);
    if(p == tr->start->back) 
    {
      *treestr++ = ',';
      treestr = TreeInner2StringREC(treestr, tr, pr, p->back, printBranchLengths, printNames, printLikelihood, rellTree,
          finalPrint, perGene, branchLabelSupport, printSHSupport, printInnerNodes);
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
    if(rellTree || branchLabelSupport || printSHSupport || printInnerNodes)
    {	 	 
      if(( !isTip(p->number, tr->mxtips)) && 
          ( !isTip(p->back->number, tr->mxtips)))
      {			      
        //assert(p->bInf != (branchInfo *)NULL);

        if(rellTree)
          sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
        if(branchLabelSupport)
          sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
        if(printSHSupport)
          sprintf(treestr, ":%8.20f[%d]", getBranchLength(tr, pr, perGene, p), p->bInf->support);
        if(printInnerNodes)
          sprintf(treestr, "[%d]", p->number);

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


static char *pllTreeToNewickREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, boolean printBranchLengths, boolean printNames,
			    boolean printLikelihood, boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport)
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


char *pllTreeToNewick(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood,
		  boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport)
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


/*=======================================================================*/
/*                         Read a tree from a file                       */
/*=======================================================================*/


/*  1.0.A  Processing of quotation marks in comment removed
 */

static int treeFinishCom (FILE *fp, char **strp)
{
  int  ch;
  
  while ((ch = getc(fp)) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(fp, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
} /* treeFinishCom */


static int treeGetCh (FILE *fp)         /* get next nonblank, noncomment character */
{ /* treeGetCh */
  int  ch;

  while ((ch = getc(fp)) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[') {                   /* comment; find its end */
      if ((ch = treeFinishCom(fp, (char **) NULL)) == EOF)  break;
    }
    else  break;
  }
  
  return  ch;
} /* treeGetCh */


static boolean treeLabelEnd (int ch)
{
  switch (ch) 
    {
    case EOF:  
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


static boolean  treeGetLabel (FILE *fp, char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = getc(fp);
  done = treeLabelEnd(ch);

  lblfound = ! done;
  quoted = (ch == '\'');
  if (quoted && ! done) 
    {
      ch = getc(fp); 
      done = (ch == EOF);
    }

  while (! done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = getc(fp); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = getc(fp);
      if (ch == EOF) break;
    }

  if (ch != EOF)  (void) ungetc(ch, fp);

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}


static boolean  treeFlushLabel (FILE *fp)
{ 
  return  treeGetLabel(fp, (char *) NULL, (int) 0);
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


static int treeFindTipName(FILE *fp, pllInstance *tr)
{
  char    str[PLL_NMLNGTH + 2];
  int      n;

  if(treeGetLabel(fp, str, PLL_NMLNGTH + 2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
} 



static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = PLL_TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = PLL_TRUE;
    }
    else {
      waswhite = PLL_FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
} /* treeEchoContext */


static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return PLL_FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  PLL_FALSE;
  }
  
  return  PLL_TRUE;
}


static int treeFlushLen (FILE  *fp)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(fp);
  
  if (ch == ':') 
    {
      ch = treeGetCh(fp);
      
      ungetc(ch, fp);
      if(!treeProcessLength(fp, & dummy)) return 0;
      return 1;	  
    }
  
  
  
  if (ch != EOF) (void) ungetc(ch, fp);
  return 1;
} 





static boolean treeNeedCh (FILE *fp, int c1, char *where)
{
  int  c2;
  
  if ((c2 = treeGetCh(fp)) == c1)  return PLL_TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where);
  if (c2 == EOF) 
    {
      printf("End-of-File");
    }
  else 
    {      	
      ungetc(c2, fp);
      treeEchoContext(fp, stdout, 40);
    }
  putchar('\n');

  if(c1 == ':')    
    printf("RAxML may be expecting to read a tree that contains branch lengths\n");

  return PLL_FALSE;
} 



static boolean addElementLen (FILE *fp, pllInstance *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(fp)) == '(') 
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
	      assert(!readNodeLabels);
	      tr->rooted = PLL_TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (! addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return PLL_FALSE;
      if (! treeNeedCh(fp, ',', "in"))             return PLL_FALSE;
      if (! addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return PLL_FALSE;
      if (! treeNeedCh(fp, ')', "in"))             return PLL_FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (fp, label, 10))
	    {	
	      int val = sscanf(label, "%d", &support);
      
	      assert(val == 1);

	      /*printf("LABEL %s Number %d\n", label, support);*/
	      p->support = q->support = support;
	      /*printf("%d %d %d %d\n", p->support, q->support, p->number, q->number);*/
	      assert(p->number > tr->mxtips && q->number > tr->mxtips);
	      *lcount = *lcount + 1;
	    }
	}
      else	
	(void) treeFlushLabel(fp);
    }
  else 
    {   
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)          return PLL_FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (! treeNeedCh(fp, ':', "in"))                 return PLL_FALSE;
      if (! treeProcessLength(fp, &branch))            return PLL_FALSE;
      
      /*printf("Branch %8.20f %d\n", branch, tr->numBranches);*/
      hookupFull(p, q, &branch);
    }
  else
    {
      fres = treeFlushLen(fp);
      if(!fres) return PLL_FALSE;
      
      hookupDefault(p, q);
    }
  return PLL_TRUE;          
} 











static nodeptr uprootTree (pllInstance *tr, nodeptr p, boolean readBranchLengths, boolean readConstraint, int numBranches)
{
  nodeptr  q, r, s, start;
  int      n, i;              

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips - 1; i++)
    assert(i == tr->nodep[i]->number);

  
  if(isTip(p->number, tr->mxtips) || p->back) 
    {
      printf("ERROR: Unable to uproot tree.\n");
      printf("       Inappropriate node marked for removal.\n");
      assert(0);
    }
  
  assert(p->back == (nodeptr)NULL);
  
  tr->nextnode = tr->nextnode - 1;

  assert(tr->nextnode < 2 * tr->mxtips);
  
  n = tr->nextnode;               
  
  assert(tr->nodep[tr->nextnode]);

  if (n != tr->mxtips + tr->ntips - 1) 
    {
      printf("ERROR: Unable to uproot tree.  Inconsistent\n");
      printf("       number of tips and nodes for rooted tree.\n");
      assert(0);
    }

  q = p->next->back;                  /* remove p from tree */
  r = p->next->next->back;
  assert(p->back == (nodeptr)NULL);
    
  if(readBranchLengths)
    {
      double b[PLL_NUM_BRANCHES];
      int i;
      for(i = 0; i < numBranches; i++)
	b[i] = (r->z[i] + q->z[i]);
      hookup (q, r, b, numBranches);
    }
  else    
    hookupDefault(q, r);

  if(readConstraint && tr->grouped)
    {    
      if(tr->constraintVector[p->number] != 0)
	{
	  printf("Root node to remove should have top-level grouping of 0\n");
	  assert(0);
	}
    }  
 
  assert(!(isTip(r->number, tr->mxtips) && isTip(q->number, tr->mxtips))); 

  assert(p->number > tr->mxtips);

  if(tr->ntips > 2 && p->number != n) 
    {    	
      q = tr->nodep[n];            /* transfer last node's conections to p */
      r = q->next;
      s = q->next->next;
      
      if(readConstraint && tr->grouped)	
	tr->constraintVector[p->number] = tr->constraintVector[q->number];       
      
      hookup(p,             q->back, q->z, numBranches);   /* move connections to p */
      hookup(p->next,       r->back, r->z, numBranches);
      hookup(p->next->next, s->back, s->z, numBranches);
      
      q->back = q->next->back = q->next->next->back = (nodeptr) NULL;
    }
  else    
    p->back = p->next->back = p->next->next->back = (nodeptr) NULL;
  
  assert(tr->ntips > 2);
  
  start = findAnyTip(tr->nodep[tr->mxtips + 1], tr->mxtips);
  
  assert(isTip(start->number, tr->mxtips));
  tr->rooted = PLL_FALSE;
  return  start;
}


int treeReadLen (FILE *fp, pllInstance *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly)
{
  nodeptr  
    p;
  
  int      
    i, 
    ch, 
    lcount = 0; 

  for (i = 1; i <= tr->mxtips; i++) 
    {
      tr->nodep[i]->back = (node *) NULL; 
      if(topologyOnly)
	tr->nodep[i]->support = -1;
    }

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;

      if(topologyOnly)
	{
	  tr->nodep[i]->support = -2;
	  tr->nodep[i]->next->support = -2;
	  tr->nodep[i]->next->next->support = -2;
	}
    }

  if(topologyOnly)
    tr->start       = tr->nodep[tr->mxtips];
  else
    tr->start       = tr->nodep[1];

  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;      
 
  for(i = 0; i < PLL_NUM_BRANCHES; i++)
    tr->partitionSmoothed[i] = PLL_FALSE;
  
  tr->rooted      = PLL_FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
      
  if(!topologyOnly)
    assert(readBranches == PLL_FALSE && readNodeLabels == PLL_FALSE);
  
       
  if (! addElementLen(fp, tr, p, readBranches, readNodeLabels, &lcount))                 
    assert(0);
  if (! treeNeedCh(fp, ',', "in"))                
    assert(0);
  if (! addElementLen(fp, tr, p->next, readBranches, readNodeLabels, &lcount))
    assert(0);
  if (! tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{ 
	  if (! addElementLen(fp, tr, p->next->next, readBranches, readNodeLabels, &lcount))
	    assert(0);	    
	}
      else 
	{                                    /*  A rooted format */
	  tr->rooted = PLL_TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}	
    }
  else 
    {      
      p->next->next->back = (nodeptr) NULL;
    }
  if (! treeNeedCh(fp, ')', "in"))                
    assert(0);

  if(topologyOnly)
    assert(!(tr->rooted && readNodeLabels));

  (void) treeFlushLabel(fp);
  
  if (! treeFlushLen(fp))                         
    assert(0);
 
  if (! treeNeedCh(fp, ';', "at end of"))       
    assert(0);
  
  if (tr->rooted) 
    {     
      assert(!readNodeLabels);

      p->next->next->back = (nodeptr) NULL;      
      //DIEGO: CHECK THIS
      tr->start = uprootTree(tr, p->next->next, PLL_FALSE, PLL_FALSE, PLL_NUM_BRANCHES);
      if (! tr->start)                              
	{
	  printf("FATAL ERROR UPROOTING TREE\n");
	  assert(0);
	}    
    }
  else    
    tr->start = findAnyTip(p, tr->mxtips);    
  
  
 
  assert(tr->ntips == tr->mxtips);
  
 
   
  
  return lcount;
}








void getStartingTree(pllInstance *tr)
{
  FILE *treeFile = myfopen(tree_file, "rb");

  tr->likelihood = PLL_UNLIKELY;
   
  treeReadLen(treeFile, tr, PLL_FALSE, PLL_FALSE, PLL_FALSE);
               
  fclose(treeFile);
 
  tr->start = tr->nodep[1];
}



