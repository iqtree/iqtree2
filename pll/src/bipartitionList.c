/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
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
 *  stamatak@ics.forth.gr
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *  
 *  Alexandros Stamatakis: "An Efficient Program for phylogenetic Inference Using Simulated Annealing". 
 *  Proceedings of IPDPS2005,  Denver, Colorado, April 2005.
 *  
 *  AND
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */
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
#include <stdint.h>
#include <assert.h>
#include "pll.h"



extern const unsigned int mask32[32];


static void getxnodeBips (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->xBips || (s = s->next)->xBips)
    {
      p->xBips = s->xBips;
      s->xBips = 0;
    }

  assert(p->xBips);
}


entry *initEntry(void)
{
  entry *e = (entry*)rax_malloc(sizeof(entry));

  e->bitVector     = (unsigned int*)NULL;
  e->treeVector    = (unsigned int*)NULL;
  e->supportVector = (int*)NULL;
  e->bipNumber  = 0;
  e->bipNumber2 = 0;
  e->supportFromTreeset[0] = 0;
  e->supportFromTreeset[1] = 0;
  e->next       = (entry*)NULL;

  return e;
} 

hashtable *initHashTable(hashNumberType n)
{
  /* 
     init with primes 
    
     static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
  */

  /* init with powers of two */

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  
  hashtable *h = (hashtable*)rax_malloc(sizeof(hashtable));
  
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

 

  h->table = (entry**)rax_calloc((size_t)tableSize, sizeof(entry*));
  h->tableSize = tableSize;  
  h->entryCount = 0;  

  return h;
}




void freeHashTable(hashtable *h)
{
  hashNumberType
    i,
    entryCount = 0;
   

  for(i = 0; i < h->tableSize; i++)
    {
      if(h->table[i] != NULL)
	{
	  entry *e = h->table[i];
	  entry *previous;	 

	  do
	    {
	      previous = e;
	      e = e->next;

	      if(previous->bitVector)
		rax_free(previous->bitVector);

	      if(previous->treeVector)
		rax_free(previous->treeVector);

	      if(previous->supportVector)
		rax_free(previous->supportVector);
	      
	      rax_free(previous);	      
	      entryCount++;
	    }
	  while(e != NULL);	  
	}

    }

  assert(entryCount == h->entryCount);
 
  rax_free(h->table);
}



void cleanupHashTable(hashtable *h, int state)
{
  hashNumberType
    k,
    entryCount = 0,
    removeCount = 0;
 
  assert(state == 1 || state == 0);

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];
	  entry *start     = (entry*)NULL;
	  entry *lastValid = (entry*)NULL;
	  	  
	  do
	    {	   	 	      	
	      if(state == 0)
		{
		  e->treeVector[0] = e->treeVector[0] & 2;	
		  assert(!(e->treeVector[0] & 1));
		}
	      else
		{
		  e->treeVector[0] = e->treeVector[0] & 1;
		  assert(!(e->treeVector[0] & 2));
		}
	      
	      if(e->treeVector[0] != 0)
		{
		  if(!start)
		    start = e;
		  lastValid = e;
		  e = e->next;
		}	  
	      else
		{
		  entry *remove = e;
		  e = e->next;
		  
		  removeCount++;

		  if(lastValid)		    		    
		    lastValid->next = remove->next;

		  if(remove->bitVector)
		    rax_free(remove->bitVector);
		  if(remove->treeVector)
		    rax_free(remove->treeVector);
		  if(remove->supportVector)
		    rax_free(remove->supportVector);
		  rax_free(remove);		 
		}
	      
	      entryCount++;	     	     
	    }
	  while(e != NULL);	 

	  if(!start)
	    {
	      assert(!lastValid);
	      h->table[k] = NULL;
	    }
	  else
	    {
	      h->table[k] = start;
	    }	 	 
	}    
    }

  assert(entryCount ==  h->entryCount);  

  h->entryCount -= removeCount;
}











unsigned int **initBitVector(int mxtips, unsigned int *vectorLength)
{
  unsigned int 
    **bitVectors = (unsigned int **)rax_malloc(sizeof(unsigned int*) * 2 * (size_t)mxtips);
  
  int 
    i;

  if(mxtips % PLL_MASK_LENGTH == 0)
    *vectorLength = mxtips / PLL_MASK_LENGTH;
  else
    *vectorLength = 1 + (mxtips / PLL_MASK_LENGTH); 
  
  for(i = 1; i <= mxtips; i++)
    {
      bitVectors[i] = (unsigned int *)rax_calloc((size_t)(*vectorLength), sizeof(unsigned int));
      assert(bitVectors[i]);
      bitVectors[i][(i - 1) / PLL_MASK_LENGTH] |= mask32[(i - 1) % PLL_MASK_LENGTH];
    }
  
  for(i = mxtips + 1; i < 2 * mxtips; i++) 
    {
      bitVectors[i] = (unsigned int *)rax_malloc(sizeof(unsigned int) * (size_t)(*vectorLength));
      assert(bitVectors[i]);
    }

  return bitVectors;
}

void freeBitVectors(unsigned int **v, int n)
{
  int i;

  for(i = 1; i < n; i++)
    rax_free(v[i]);
}





static void newviewBipartitions(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, int processID)
{
  
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    
    
    
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    unsigned 
      int i;      
    
    assert(processID == 0);
    

    while(!p->xBips)
      {	
	if(!p->xBips)
	  getxnodeBips(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(!r->xBips)
	      {
		if(!r->xBips)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength, processID);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((!r->xBips) || (!q->xBips))
	      {
		if(!q->xBips)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength, processID);
		if(!r->xBips)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength, processID);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	 
	  }

      }     
  }     
}




static void insertHashRF(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber, int treeVectorLength, hashNumberType position, int support, 
			 boolean computeWRF)
{     
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / PLL_MASK_LENGTH] |= mask32[treeNumber % PLL_MASK_LENGTH];
	      if(computeWRF)
		{
		  e->supportVector[treeNumber] = support;
		 
		  assert(0 <= treeNumber && treeNumber < treeVectorLength * PLL_MASK_LENGTH);
		}
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
       
      /*e->bitVector  = (unsigned int*rax_callocvectorLength, sizeof(unsigned int));*/
      e->bitVector = (unsigned int*)rax_malloc_aligned((size_t)vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      e->treeVector = (unsigned int*)rax_calloc((size_t)treeVectorLength, sizeof(unsigned int));
      if(computeWRF)
	e->supportVector = (int*)rax_calloc((size_t)treeVectorLength * PLL_MASK_LENGTH, sizeof(int));

      e->treeVector[treeNumber / PLL_MASK_LENGTH] |= mask32[treeNumber % PLL_MASK_LENGTH];
      if(computeWRF)
	{
	  e->supportVector[treeNumber] = support;
	 
	  assert(0 <= treeNumber && treeNumber < treeVectorLength * PLL_MASK_LENGTH);
	}

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 
       
      /*e->bitVector  = (unsigned int*rax_callocvectorLength, sizeof(unsigned int)); */

      e->bitVector = (unsigned int*)rax_malloc_aligned((size_t)vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      e->treeVector = (unsigned int*)rax_calloc((size_t)treeVectorLength, sizeof(unsigned int));
      if(computeWRF)	
	e->supportVector = (int*)rax_calloc((size_t)treeVectorLength * PLL_MASK_LENGTH, sizeof(int));


      e->treeVector[treeNumber / PLL_MASK_LENGTH] |= mask32[treeNumber % PLL_MASK_LENGTH];
      if(computeWRF)
	{
	  e->supportVector[treeNumber] = support;
	 
	  assert(0 <= treeNumber && treeNumber < treeVectorLength * PLL_MASK_LENGTH);
	}

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}



void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, 
			     int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF, int processID)
{
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr 
	q = p->next;          

      do 
	{
	  bitVectorInitravSpecial(bitVectors, q->back, numsp, vectorLength, h, treeNumber, function, bInf, countBranches, treeVectorLength, traverseOnly, computeWRF, processID);
	  q = q->next;
	}
      while(q != p);
           
      newviewBipartitions(bitVectors, p, numsp, vectorLength, processID);
      
      assert(p->xBips);

      assert(!traverseOnly);     

      if(!(isTip(p->back->number, numsp)))
	{
	  unsigned int 
	    *toInsert  = bitVectors[p->number];
	  
	  hashNumberType 
	    position = p->hash % h->tableSize;
	 
	  assert(!(toInsert[0] & 1));
	  assert(!computeWRF);
	  
	  switch(function)
	    {	     
	    case PLL_BIPARTITIONS_RF:	     
	      insertHashRF(toInsert, h, vectorLength, treeNumber, treeVectorLength, position, 0, computeWRF);
	      *countBranches =  *countBranches + 1;
	      break;
	    default:
	      assert(0);
	    }	  	  
	}
      
    }
}






















double convergenceCriterion(hashtable *h, int mxtips)
{
  int      
    rf = 0; 

  unsigned int 
    collisions = 0,
    k = 0, 
    entryCount = 0;
  
  double    
    rrf;  

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  unsigned int 
	    slotCollisions = 0;

	  do
	    {
	      unsigned int *vector = e->treeVector;	     
	      if(((vector[0] & 1) > 0) + ((vector[0] & 2) > 0) == 1)
		rf++;	     
	      
	      entryCount++;
	      slotCollisions++;
	      e = e->next;
	    }
	  while(e != NULL);

	  collisions += (slotCollisions - 1);
	}     
    }

  assert(entryCount == h->entryCount);  
      
  rrf = (double)rf/((double)(2 * (mxtips - 3)));  

#ifdef _DEBUG_CHECKPOINTING
  printf("Collisions: %u\n", collisions);
#endif

  return rrf;
}




