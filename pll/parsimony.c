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
 * @file parsimony.c
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

#if defined(__MIC_NATIVE)

#include <immintrin.h>

#define INTS_PER_VECTOR 16
#define LONG_INTS_PER_VECTOR 8
#define INT_TYPE __m512i
#define CAST double*
#define SET_ALL_BITS_ONE _mm512_set1_epi32(0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm512_setzero_epi32()
#define VECTOR_LOAD _mm512_load_epi32
#define VECTOR_STORE  _mm512_store_epi32
#define VECTOR_BIT_AND _mm512_and_epi32
#define VECTOR_BIT_OR  _mm512_or_epi32
#define VECTOR_AND_NOT _mm512_andnot_epi32

#elif defined(__AVX)

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <immintrin.h>
#include <pmmintrin.h>
#endif

#define ULINT_SIZE 64
#define INTS_PER_VECTOR 8
#define LONG_INTS_PER_VECTOR 4
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd

#elif (defined(__SSE3))

#if defined(__ARM_NEON)
#include "sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif
  
#define INTS_PER_VECTOR 4
#ifdef __i386__
#define ULINT_SIZE 32
#define LONG_INTS_PER_VECTOR 4
#else
#define ULINT_SIZE 64
#define LONG_INTS_PER_VECTOR 2
#endif
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128

#endif

#include "pll.h"
#include "pllInternal.h"

extern const unsigned int mask32[32]; 

static __inline unsigned int vectorPopcount(INT_TYPE v)
{
  unsigned long
    counts[LONG_INTS_PER_VECTOR] __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)));

  int    
    i,
    sum = 0;

  VECTOR_STORE((CAST)counts, v);

  for(i = 0; i < LONG_INTS_PER_VECTOR; i++)
     sum += __builtin_popcountl(counts[i]);

  return ((unsigned int)sum);
}

static __inline void storePerSiteScores (partitionList * pr, int model, INT_TYPE v, unsigned int offset)
{
  unsigned long
    counts[LONG_INTS_PER_VECTOR] __attribute__ ((aligned (PLL_BYTE_ALIGNMENT)));
  parsimonyNumber * buf;

  int    
    i,
    j;
  
  VECTOR_STORE((CAST)counts, v);

  for (i = 0; i < LONG_INTS_PER_VECTOR; ++i)
   {
     buf = &(pr->partitionData[model]->perSiteParsScores[offset * PLL_PCF + i * ULINT_SIZE]);
     for (j = 0; j < ULINT_SIZE; ++ j)
        buf[j] += ((counts[i] >> j) & 1);
   }
  
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->xPars || (s = s->next)->xPars)
    {
      p->xPars = s->xPars;
      s->xPars = 0;
    }

  assert(p->next->xPars || p->next->next->xPars || p->xPars);

}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, pllBoolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->xPars)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
         computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
        computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->xPars) 
        computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->xPars) 
        computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}

/* check whether site contains at least 2 different letters, i.e.
   whether it will generate a score */
static pllBoolean isInformative(pllInstance *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
        
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {      
      nucleotide = tr->yVector[j][site];            
      check[nucleotide] = 1;
      assert(bitVector[nucleotide] > 0);                   
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
        informativeCounter++;    
    } 
          
  if(informativeCounter > 1)
    return PLL_TRUE;    

  return PLL_FALSE;          
}

static void compressDNA(pllInstance *tr, partitionList *pr, int *informative, int perSiteScores)
{
  size_t
    totalNodes,
    i,
    model;
   
  totalNodes = 2 * (size_t)tr->mxtips;

 

  for(model = 0; model < (size_t) pr->numberOfPartitions; model++)
    {
      size_t
        k,
        states = (size_t)pr->partitionData[model]->states,
        compressedEntries,
        compressedEntriesPadded,
        entries = 0, 
        lower = pr->partitionData[model]->lower,
        upper = pr->partitionData[model]->upper;

      parsimonyNumber 
        **compressedTips = (parsimonyNumber **)rax_malloc(states * sizeof(parsimonyNumber*)),
        *compressedValues = (parsimonyNumber *)rax_malloc(states * sizeof(parsimonyNumber));
      
      for(i = lower; i < upper; i++)    
        if(informative[i])
          entries += (size_t)tr->aliaswgt[i];     
  
      compressedEntries = entries / PLL_PCF;

      if(entries % PLL_PCF != 0)
        compressedEntries++;

#if (defined(__SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
        compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
        compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      
      rax_posix_memalign ((void **) &(pr->partitionData[model]->parsVect), PLL_BYTE_ALIGNMENT, (size_t)compressedEntriesPadded * states * totalNodes * sizeof(parsimonyNumber));
      if (perSiteScores)
       {
         rax_posix_memalign ((void **) &(pr->partitionData[model]->perSiteParsScores), PLL_BYTE_ALIGNMENT, (size_t)pr->partitionData[model]->width* sizeof (parsimonyNumber));
         for (i = 0; i < (size_t)pr->partitionData[model]->width; ++i) pr->partitionData[model]->perSiteParsScores[i] = 0;
       }

     
      for(i = 0; i < compressedEntriesPadded * states * totalNodes; i++)      
        pr->partitionData[model]->parsVect[i] = 0;

      for(i = 0; i < (size_t)tr->mxtips; i++)
        {
          size_t
            w = 0,
            compressedIndex = 0,
            compressedCounter = 0,
            index = 0;

          for(k = 0; k < states; k++)
            {
              compressedTips[k] = &(pr->partitionData[model]->parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
              compressedValues[k] = 0;
            }                
              
          for(index = lower; index < (size_t)upper; index++)
            {
              if(informative[index])
                {
                  const unsigned int 
                    *bitValue = getBitVector(pr->partitionData[model]->dataType);

                  parsimonyNumber 
                    value = bitValue[tr->yVector[i + 1][index]];          
              
                  for(w = 0; w < (size_t)tr->aliaswgt[index]; w++)
                    {      
                      for(k = 0; k < states; k++)
                        {
                          if(value & mask32[k])
                            compressedValues[k] |= mask32[compressedCounter];
                        }
                     
                      compressedCounter++;
                  
                      if(compressedCounter == PLL_PCF)
                        {
                          for(k = 0; k < states; k++)
                            {
                              compressedTips[k][compressedIndex] = compressedValues[k];
                              compressedValues[k] = 0;
                            }                    
                          
                          compressedCounter = 0;
                          compressedIndex++;
                        }
                    }
                }
            }
          
          for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
            {   
              for(;compressedCounter < PLL_PCF; compressedCounter++)              
                for(k = 0; k < states; k++)
                  compressedValues[k] |= mask32[compressedCounter];               
          
              for(k = 0; k < states; k++)
                {
                  compressedTips[k][compressedIndex] = compressedValues[k];
                  compressedValues[k] = 0;
                }                     
              
              compressedCounter = 0;
            }           
        }
  
      pr->partitionData[model]->parsimonyLength = compressedEntriesPadded;

      rax_free(compressedTips);
      rax_free(compressedValues);
    }
  
  rax_posix_memalign ((void **) &(tr->parsimonyScore), PLL_BYTE_ALIGNMENT, sizeof(unsigned int) * totalNodes);  
          
  for(i = 0; i < totalNodes; i++) 
    tr->parsimonyScore[i] = 0;
}

static void determineUninformativeSites(pllInstance *tr, partitionList *pr, int *informative)
{
  int 
    model,
    number = 0,
    i;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */


  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      for(i = pr->partitionData[model]->lower; i < pr->partitionData[model]->upper; i++)
        {
           if(isInformative(tr, pr->partitionData[model]->dataType, i))
             informative[i] = 1;
           else
             {
               informative[i] = 0;
               number++;
             }  
        }      
    }

  /* printf("Uninformative Patterns: %d\n", number); */
}

void pllInitParsimonyStructures(pllInstance *tr, partitionList *pr, pllBoolean perSiteScores)
{
  int 
    i,
    *informative = (int *)rax_malloc(sizeof(int) * (size_t)tr->originalCrunchedLength);

  for (i = 0; i < pr->numberOfPartitions; ++ i)
     rax_free (pr->partitionData[i]->parsVect);

  rax_free (tr->parsimonyScore);
 
  determineUninformativeSites(tr, pr, informative);

  compressDNA(tr, pr, informative, perSiteScores);

  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 1; i++)
    {
      nodeptr 
        p = tr->nodep[i];

      p->xPars             = 1;
      p->next->xPars       = 0;
      p->next->next->xPars = 0;
    }

  tr->ti = (int*)rax_malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  rax_free(informative); 
}

static void newviewParsimonyIterativeFast(pllInstance *tr, partitionList *pr, pllBoolean perSiteScores)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;

  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
        totalScore = 0;

      size_t
        pNumber = (size_t)ti[index],
        qNumber = (size_t)ti[index + 1],
        rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < pr->numberOfPartitions; model++)
        {
          size_t
            k,
            states = pr->partitionData[model]->states,
            width = pr->partitionData[model]->parsimonyLength;
            
          unsigned int  
            i;      
                 
          switch(states)
            {
            case 2:       
              {
                parsimonyNumber
                  *left[2],
                  *right[2],
                  *this[2];

                for(k = 0; k < 2; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 2 * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * pNumber) + width * k]);
                  }

                for(i = 0; i < width; i += INTS_PER_VECTOR)
                  {               
                    INT_TYPE
                      s_r, s_l, v_N,
                      l_A, l_C,
                      v_A, v_C;          
                    
                    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
                    l_A = VECTOR_BIT_AND(s_l, s_r);
                    v_A = VECTOR_BIT_OR(s_l, s_r);
                    
                    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
                    l_C = VECTOR_BIT_AND(s_l, s_r);
                    v_C = VECTOR_BIT_OR(s_l, s_r);                                                                
                    
                    v_N = VECTOR_BIT_OR(l_A, l_C);
                    
                    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
                    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);            
                    if (perSiteScores)
                       storePerSiteScores (pr, model, v_N, i);
                  }
              }
              break;
            case 4:
              {
                parsimonyNumber
                  *left[4],
                  *right[4],
                  *this[4];

                for(k = 0; k < 4; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 4 * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * pNumber) + width * k]);
                  }
                for(i = 0; i < width; i += INTS_PER_VECTOR)
                  {               
                    INT_TYPE
                      s_r, s_l, v_N,
                      l_A, l_C, l_G, l_T,
                      v_A, v_C, v_G, v_T;                
                    
                    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
                    l_A = VECTOR_BIT_AND(s_l, s_r);
                    v_A = VECTOR_BIT_OR(s_l, s_r);
                    
                    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
                    l_C = VECTOR_BIT_AND(s_l, s_r);
                    v_C = VECTOR_BIT_OR(s_l, s_r);
                    
                    s_l = VECTOR_LOAD((CAST)(&left[2][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[2][i]));
                    l_G = VECTOR_BIT_AND(s_l, s_r);
                    v_G = VECTOR_BIT_OR(s_l, s_r);
                    
                    s_l = VECTOR_LOAD((CAST)(&left[3][i]));
                    s_r = VECTOR_LOAD((CAST)(&right[3][i]));
                    l_T = VECTOR_BIT_AND(s_l, s_r);
                    v_T = VECTOR_BIT_OR(s_l, s_r);
                    
                    v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));                                
                    
                    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
                    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
                    VECTOR_STORE((CAST)(&this[2][i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
                    VECTOR_STORE((CAST)(&this[3][i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);  
                    
                    if (perSiteScores)
                       storePerSiteScores (pr, model, v_N, i);
                  }
              }
              break;
            case 20:
              {
                parsimonyNumber
                  *left[20],
                  *right[20],
                  *this[20];

                for(k = 0; k < 20; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 20 * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * pNumber) + width * k]);
                  }

                for(i = 0; i < width; i += INTS_PER_VECTOR)
                  {               
                    size_t j;
                    
                    INT_TYPE
                      s_r, s_l, 
                      v_N = SET_ALL_BITS_ZERO,
                      l_A[20], 
                      v_A[20];           
                    
                    for(j = 0; j < 20; j++)
                      {
                        s_l = VECTOR_LOAD((CAST)(&left[j][i]));
                        s_r = VECTOR_LOAD((CAST)(&right[j][i]));
                        l_A[j] = VECTOR_BIT_AND(s_l, s_r);
                        v_A[j] = VECTOR_BIT_OR(s_l, s_r);
                        
                        v_N = VECTOR_BIT_OR(v_N, l_A[j]);
                      }
                    
                    for(j = 0; j < 20; j++)                 
                      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);

                    if (perSiteScores)
                       storePerSiteScores (pr, model, v_N, i);
                  }
              }
              break;
            default:
              {
                parsimonyNumber
                  *left[32], 
                  *right[32],
                  *this[32];

                assert(states <= 32);
                
                for(k = 0; k < states; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * states * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * states * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * states * pNumber) + width * k]);
                  }

                for(i = 0; i < width; i += INTS_PER_VECTOR)
                  {               
                    size_t j;
                    
                    INT_TYPE
                      s_r, s_l, 
                      v_N = SET_ALL_BITS_ZERO,
                      l_A[32], 
                      v_A[32];           
                    
                    for(j = 0; j < states; j++)
                      {
                        s_l = VECTOR_LOAD((CAST)(&left[j][i]));
                        s_r = VECTOR_LOAD((CAST)(&right[j][i]));
                        l_A[j] = VECTOR_BIT_AND(s_l, s_r);
                        v_A[j] = VECTOR_BIT_OR(s_l, s_r);
                        
                        v_N = VECTOR_BIT_OR(v_N, l_A[j]);
                      }
                    
                    for(j = 0; j < states; j++)             
                      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);

                    if (perSiteScores)
                       storePerSiteScores (pr, model, v_N, i);
                  }                             
              }
            }            
        }
      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}

static unsigned int evaluateParsimonyIterativeFast(pllInstance *tr, partitionList *pr, pllBoolean perSiteScores)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr, pr, perSiteScores);

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      size_t
        k,
        states = pr->partitionData[model]->states,
        width  = pr->partitionData[model]->parsimonyLength,
        i;

       switch(states)
         {
         case 2:
           {
             parsimonyNumber
               *left[2],
               *right[2];
             
             for(k = 0; k < 2; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * 2 * pNumber) + width * k]);
               }     
             
             for(i = 0; i < width; i += INTS_PER_VECTOR)
               {                                               
                 INT_TYPE      
                   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
                   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),            
                   v_N = VECTOR_BIT_OR(l_A, l_C);
                 
                 v_N = VECTOR_AND_NOT(v_N, allOne);
                 
                 sum += vectorPopcount(v_N);
                  if (perSiteScores)
                    storePerSiteScores (pr, model, v_N, i);
               }
           }
           break;
         case 4:
           {
             parsimonyNumber
               *left[4],
               *right[4];
      
             for(k = 0; k < 4; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * 4 * pNumber) + width * k]);
               }        

             for(i = 0; i < width; i += INTS_PER_VECTOR)
               {                                                
                 INT_TYPE      
                   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
                   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),
                   l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[2][i])), VECTOR_LOAD((CAST)(&right[2][i]))),
                   l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[3][i])), VECTOR_LOAD((CAST)(&right[3][i]))),
                   v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     
                 
                 v_N = VECTOR_AND_NOT(v_N, allOne);
                 
                 sum += vectorPopcount(v_N);
                  if (perSiteScores)
                    storePerSiteScores (pr, model, v_N, i);
               }                 
           }
           break;
         case 20:
           {
             parsimonyNumber
               *left[20],
               *right[20];
             
              for(k = 0; k < 20; k++)
                {
                  left[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * qNumber) + width * k]);
                  right[k] = &(pr->partitionData[model]->parsVect[(width * 20 * pNumber) + width * k]);
                }  
           
              for(i = 0; i < width; i += INTS_PER_VECTOR)
                {                              
                  int 
                    j;
                  
                  INT_TYPE      
                    l_A,
                    v_N = SET_ALL_BITS_ZERO;     
                  
                  for(j = 0; j < 20; j++)
                    {
                      l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
                      v_N = VECTOR_BIT_OR(l_A, v_N);
                    }
                  
                  v_N = VECTOR_AND_NOT(v_N, allOne);
                  
                  sum += vectorPopcount(v_N);          
                  if (perSiteScores)
                    storePerSiteScores (pr, model, v_N, i);
                }
           }
           break;
         default:
           {
             parsimonyNumber
               *left[32],  
               *right[32]; 

             assert(states <= 32);

             for(k = 0; k < states; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * states * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * states * pNumber) + width * k]);
               }  
           
             for(i = 0; i < width; i += INTS_PER_VECTOR)
               {                               
                 size_t
                   j;
                 
                 INT_TYPE      
                   l_A,
                   v_N = SET_ALL_BITS_ZERO;     
                 
                 for(j = 0; j < states; j++)
                   {
                     l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
                     v_N = VECTOR_BIT_OR(l_A, v_N);
                   }
                 
                 v_N = VECTOR_AND_NOT(v_N, allOne);
                 
                 sum += vectorPopcount(v_N);           
                 if (perSiteScores)
                   storePerSiteScores (pr, model, v_N, i);
               }
           }
         }
    }
  
  return sum;
}

unsigned int pllEvaluateParsimony(pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean full, pllBoolean perSiteScores)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
        computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
        computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->xPars)
        computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->xPars)
        computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;

  result = evaluateParsimonyIterativeFast(tr, pr, perSiteScores);

  return result;
}
