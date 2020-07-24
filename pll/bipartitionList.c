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
#include "pllInternal.h"


static pllBipartitionEntry *initEntry(void);
static void getxnodeBips (nodeptr p);
static void newviewBipartitions(unsigned int **bitVectors, 
                                nodeptr p, 
                                int numsp, 
                                unsigned int vectorLength, 
                                int processID);

static void insertHashRF(unsigned int *bitVector, 
                         pllHashTable *h, 
                         unsigned int vectorLength, 
                         int treeNumber, 
                         int treeVectorLength, 
                         hashNumberType position, 
                         int support, 
                         pllBoolean computeWRF);

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


static pllBipartitionEntry *initEntry(void)
{
  pllBipartitionEntry * e = (pllBipartitionEntry *)rax_malloc(sizeof(pllBipartitionEntry));

  e->bitVector     = (unsigned int*)NULL;
  e->treeVector    = (unsigned int*)NULL;
  e->supportVector = (int*)NULL;
  e->bipNumber  = 0;
  e->bipNumber2 = 0;
  e->supportFromTreeset[0] = 0;
  e->supportFromTreeset[1] = 0;
  e->next       = (pllBipartitionEntry *)NULL;

  return e;
} 

void cleanupHashTable(pllHashTable *h, int state)
{
  unsigned int
    k,
    entryCount = 0,
    removeCount = 0;
 
  assert(state == 1 || state == 0);

  for(k = 0, entryCount = 0; k < h->size; k++)       
    { 
      pllHashItem * start     = NULL;
      pllHashItem * lastValid = NULL;
      
      pllHashItem * hitem = h->Items[k];
      while (hitem)
       {                           
         pllBipartitionEntry *e = (pllBipartitionEntry *)(hitem->data);
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
               start = hitem;
             lastValid = hitem;
             hitem = hitem->next;
           }         
         else
           {
             pllHashItem *tmp = hitem;
             pllBipartitionEntry *remove = e;
             hitem = hitem->next;
             
             removeCount++;

             if(lastValid) lastValid->next = hitem;

             if(remove->bitVector)     rax_free(remove->bitVector);
             if(remove->treeVector)    rax_free(remove->treeVector);
             if(remove->supportVector) rax_free(remove->supportVector);
             rax_free(remove);              
             rax_free(tmp);
           }
         entryCount++;
       }

      if(!start)
        {
          assert(!lastValid);
          h->Items[k] = NULL;
        }
      else
        {
          h->Items[k] = start;
        }            
    }

  assert(entryCount ==  h->entries);
  h->entries-= removeCount;
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


static void newviewBipartitions(unsigned int **bitVectors, 
                                nodeptr p, 
                                int numsp, 
                                unsigned int vectorLength, 
                                int processID)
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




static void insertHashRF(unsigned int *bitVector, 
                         pllHashTable *h, 
                         unsigned int vectorLength, 
                         int treeNumber, 
                         int treeVectorLength, 
                         hashNumberType position, 
                         int support, 
                         pllBoolean computeWRF)
{
  pllBipartitionEntry * e;
  pllHashItem * hitem;

  if(h->Items[position] != NULL)
    {
      for (hitem = h->Items[position]; hitem; hitem = hitem->next)
        { 
          e = (pllBipartitionEntry *)(hitem->data);
          
          if (!memcmp(bitVector, e->bitVector, vectorLength * sizeof(unsigned int)))
            {
              e->treeVector[treeNumber / PLL_MASK_LENGTH] |= mask32[treeNumber % PLL_MASK_LENGTH];
              if(computeWRF)
                {
                  e->supportVector[treeNumber] = support;
                  assert(0 <= treeNumber && treeNumber < treeVectorLength * PLL_MASK_LENGTH);
                }
              return;
            }
        }
    }
  e = initEntry(); 
       
  rax_posix_memalign ((void **)&(e->bitVector), PLL_BYTE_ALIGNMENT, (size_t)vectorLength * sizeof(unsigned int));
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
  
  pllHashAdd (h, position, NULL, (void *)e);
}



void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, pllHashTable *h, int treeNumber, int function, branchInfo *bInf, 
                             int *countBranches, int treeVectorLength, pllBoolean traverseOnly, pllBoolean computeWRF, int processID)
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
            position = p->hash % h->size;
         
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

double convergenceCriterion(pllHashTable *h, int mxtips)
{
  int      
    rf = 0; 

  unsigned int 
    k = 0, 
    entryCount = 0;
  
  double    
    rrf;  

  pllHashItem * hitem;

  for(k = 0, entryCount = 0; k < h->size; k++)          
    {      
      for (hitem = h->Items[k]; hitem; hitem = hitem->next)
       {
         pllBipartitionEntry* e = (pllBipartitionEntry*)(hitem->data);
         unsigned int *vector = e->treeVector;          

         if(((vector[0] & 1) > 0) + ((vector[0] & 2) > 0) == 1)
           rf++;        
          
         entryCount++;
         e = e->next;
       }
    }

  assert(entryCount == h->entries);  
  rrf = (double)rf/((double)(2 * (mxtips - 3)));  
  return rrf;
}
