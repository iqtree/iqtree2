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
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
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

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
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

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((strlen(s) + 1) * sizeof(char));

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


static double getBranchLength(tree *tr, int perGene, nodeptr p)
{
  double 
    z = 0.0,
    x = 0.0;

  assert(perGene != NO_BRANCHES);
	      
  if(tr->numBranches == 1)
    {
      assert(tr->fracchange != -1.0);
    
      z = p->z[0];
      if (z < zmin) 
	z = zmin;      	 
      
      x = -log(z) * tr->fracchange;  
      /* printf("%f %f %f\n", tr->fracchange, x, z);           */
    }
  else
    {
      if(perGene == SUMMARIZE_LH)
	{
	  int 
	    i;
	  
	  double 
	    avgX = 0.0;
		      
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      assert(tr->partitionContributions[i] != -1.0);
	      assert(tr->fracchanges[i] != -1.0);
	      z = p->z[i];
	      if(z < zmin) 
		z = zmin;      	 
	      x = -log(z) * tr->fracchanges[i];
	      avgX += x * tr->partitionContributions[i];
	    }

	  x = avgX;
	}
      else
	{	
	  assert(tr->fracchanges[perGene] != -1.0);
	  assert(perGene >= 0 && perGene < tr->numBranches);
	  
	  z = p->z[perGene];
	  
	  if(z < zmin) 
	    z = zmin;      	 
	  
	  x = -log(z) * tr->fracchanges[perGene];	  
	}
    }

  return x;
}


  


static char *Tree2StringREC(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, 
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
      treestr = Tree2StringREC(treestr, tr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      *treestr++ = ',';
      treestr = Tree2StringREC(treestr, tr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2StringREC(treestr, tr, p->back, printBranchLengths, printNames, printLikelihood, rellTree, 
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
	      assert(0);
		      
	      /*assert(p->bInf != (branchInfo *)NULL);*/
	      
	      /*if(rellTree)
		sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	      if(branchLabelSupport)
		sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f[%d]", getBranchLength(tr, perGene, p), p->bInf->support);
	      */
	      
	    }
	  else		
	    {
	      if(rellTree || branchLabelSupport)
		sprintf(treestr, ":%8.20f", p->z[0]);	
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f", getBranchLength(tr, perGene, p));
	    }
	}
      else
	{
	  if(printBranchLengths)	    
	    sprintf(treestr, ":%8.20f", getBranchLength(tr, perGene, p));	      	   
	  else	    
	    sprintf(treestr, "%s", "\0");	    
	}      
    }
  
  while (*treestr) treestr++;
  return  treestr;
}





    




char *Tree2String(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, 
		  boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport)
{ 

  if(rellTree)
    assert(!branchLabelSupport && !printSHSupport);

  if(branchLabelSupport)
    assert(!rellTree && !printSHSupport);

  if(printSHSupport)
    assert(!branchLabelSupport && !rellTree);

 
  Tree2StringREC(treestr, tr, p, printBranchLengths, printNames, printLikelihood, rellTree, 
		 finalPrint, perGene, branchLabelSupport, printSHSupport);  
    
  
  while (*treestr) treestr++;
  
  return treestr;
}


void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission)
{  
  FILE *treeFile;
  char extendedTreeFileName[1024];
  char buf[16];
  int i;

  assert(adef->perGeneBranchLengths);
     
  for(i = 0; i < tr->numBranches; i++)	
    {
      strcpy(extendedTreeFileName, fileName);
      sprintf(buf,"%d", i);
      strcat(extendedTreeFileName, ".PARTITION.");
      strcat(extendedTreeFileName, buf);
      /*printf("Partitiuon %d file %s\n", i, extendedTreeFileName);*/
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, i, FALSE, FALSE);
      treeFile = myfopen(extendedTreeFileName, permission);
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);
    }  
    
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
      return TRUE;
    default:
      break;
    }
  return FALSE;
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


static int treeFindTipName(FILE *fp, tree *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(fp, str, nmlngth+2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
} 



static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
} /* treeEchoContext */


static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  FALSE;
  }
  
  return  TRUE;
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
  
  if ((c2 = treeGetCh(fp)) == c1)  return TRUE;
  
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

  return FALSE;
} 



static boolean addElementLen (FILE *fp, tree *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
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
	      return FALSE;
	    }
	  else 
	    {
	      assert(!readNodeLabels);
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (! addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return FALSE;
      if (! treeNeedCh(fp, ',', "in"))             return FALSE;
      if (! addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return FALSE;
      if (! treeNeedCh(fp, ')', "in"))             return FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (fp, label, 10))
	    {	
	      int val = sscanf(label, "%d", &support);
      
	      assert(val == 1);

	      /*printf("LABEL %s Number %d\n", label, support);*/
	      /*p->support = q->support = support;*/
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
      if ((n = treeFindTipName(fp, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (! treeNeedCh(fp, ':', "in"))                 return FALSE;
      if (! treeProcessLength(fp, &branch))            return FALSE;
      
      /*printf("Branch %8.20f %d\n", branch, tr->numBranches);*/
      hookup(p, q, &branch, tr->numBranches);
    }
  else
    {
      fres = treeFlushLen(fp);
      if(!fres) return FALSE;
      
      hookupDefault(p, q, tr->numBranches);
    }
  return TRUE;          
} 











static nodeptr uprootTree (tree *tr, nodeptr p, boolean readBranchLengths, boolean readConstraint)
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
      double b[NUM_BRANCHES];
      int i;
      for(i = 0; i < tr->numBranches; i++)
	b[i] = (r->z[i] + q->z[i]);
      hookup (q, r, b, tr->numBranches);
    }
  else    
    hookupDefault(q, r, tr->numBranches);    

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
      
      hookup(p,             q->back, q->z, tr->numBranches);   /* move connections to p */
      hookup(p->next,       r->back, r->z, tr->numBranches);
      hookup(p->next->next, s->back, s->z, tr->numBranches);           
      
      q->back = q->next->back = q->next->next->back = (nodeptr) NULL;
    }
  else    
    p->back = p->next->back = p->next->next->back = (nodeptr) NULL;
  
  assert(tr->ntips > 2);
  
  start = findAnyTip(tr->nodep[tr->mxtips + 1], tr->mxtips);
  
  assert(isTip(start->number, tr->mxtips));
  tr->rooted = FALSE;
  return  start;
}


int treeReadLen (FILE *fp, tree *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly)
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
      /*if(topologyOnly)
	tr->nodep[i]->support = -1;*/
    }

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;

      /*if(topologyOnly)
	{
	  tr->nodep[i]->support = -2;
	  tr->nodep[i]->next->support = -2;
	  tr->nodep[i]->next->next->support = -2;
	  }*/
    }

  if(topologyOnly)
    tr->start       = tr->nodep[tr->mxtips];
  else
    tr->start       = tr->nodep[1];

  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;      
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
      
  if(!topologyOnly)
    assert(readBranches == FALSE && readNodeLabels == FALSE);
  
       
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
	  tr->rooted = TRUE;
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
      tr->start = uprootTree(tr, p->next->next, FALSE, FALSE);      
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








void getStartingTree(tree *tr)
{
  FILE *treeFile = myfopen(tree_file, "rb");

  tr->likelihood = unlikely;
   
  treeReadLen(treeFile, tr, FALSE, FALSE, FALSE);
               
  fclose(treeFile);
 
  tr->start = tr->nodep[1];
}



