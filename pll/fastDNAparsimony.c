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
 * @file fastDNAparsimony.c
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
//#define LONG_INTS_PER_VECTOR 8
#define LONG_INTS_PER_VECTOR (64/sizeof(long))
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

#define INTS_PER_VECTOR 8
//#define LONG_INTS_PER_VECTOR 4
#define LONG_INTS_PER_VECTOR (32/sizeof(long))
#define INT_TYPE __m256d
#define CAST double*
//#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
//#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define SET_ALL_BITS_ONE _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF))
#define SET_ALL_BITS_ZERO _mm256_castsi256_pd(_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000))
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
//#define LONG_INTS_PER_VECTOR 4
#define LONG_INTS_PER_VECTOR (16/sizeof(long))
#else
//#define LONG_INTS_PER_VECTOR 2
#define LONG_INTS_PER_VECTOR (16/sizeof(long))
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

#if defined (_MSC_VER)
#	if defined ( __SSE4_2__ ) || defined (__AVX__)
#if defined(__ARM_NEON)
#        include "sse2neon.h"
#else
#include <nmmintrin.h>
#endif
#		define __builtin_popcount _mm_popcnt_u32
#		define __builtin_popcountl _mm_popcnt_u64
#	else
#       if defined(CLANG_UNDER_VS)
#           define _mm_popcnt_u64 _popcnt64
#       else
#		include <intrin.h>
static __inline uint32_t __builtin_popcount (uint32_t a) {
		// popcnt instruction not available
		uint32_t b = a - ((a >> 1) & 0x55555555);
		uint32_t c = (b & 0x33333333) + ((b >> 2) & 0x33333333);
		uint32_t d = (c + (c >> 4)) & 0x0F0F0F0F;
		uint32_t e = d * 0x01010101;
		return   e >> 24;
	}
//#		define __builtin_popcount __popcnt
#		define __builtin_popcountl __popcnt64
#       endif
#	endif
#endif

static pllBoolean tipHomogeneityCheckerPars(pllInstance *tr, nodeptr p, int grouping);

extern const unsigned int mask32[32]; 
/* vector-specific stuff */


extern double masterTime;

/************************************************ pop count stuff ***********************************************/

 unsigned int bitcount_32_bit(unsigned int i)
{
  return ((unsigned int) __builtin_popcount(i));
}

/* bit count for 64 bit integers */

//__inline unsigned int bitcount_64_bit(uint64_t i)
//{
//  return ((unsigned int) __builtin_popcountl(i));
//}

/* bit count for 128 bit SSE3 and 256 bit AVX registers */

#if (defined(__SSE3) || defined(__AVX))

#if defined(_WIN32) &&!defined(WIN64)
 /* emulate with 32-bit version */
static __inline unsigned int vectorPopcount(INT_TYPE v)
{
PLL_ALIGN_BEGIN unsigned int counts[INTS_PER_VECTOR] PLL_ALIGN_END;

  int
    i,
    sum = 0;

  VECTOR_STORE((CAST)counts, v);

  for(i = 0; i < INTS_PER_VECTOR; i++)
    sum += __builtin_popcount(counts[i]);

  return ((unsigned int)sum);
}
#else

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
#endif

#endif



/********************************DNA FUNCTIONS *****************************************************************/


static int checkerPars(pllInstance *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
        return group;

      group = checkerPars(tr, p->next->back);
      if(group != -9) 
        return group;

      group = checkerPars(tr, p->next->next->back);
      if(group != -9) 
        return group;

      return -9;
    }
}

static pllBoolean tipHomogeneityCheckerPars(pllInstance *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->constraintVector[p->number] != grouping) 
        return PLL_FALSE;
      else 
        return PLL_TRUE;
    }
  else
    {   
      return  (tipHomogeneityCheckerPars(tr, p->next->back, grouping) && tipHomogeneityCheckerPars(tr, p->next->next->back,grouping));      
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







#if (defined(__SSE3) || defined(__AVX))

static void newviewParsimonyIterativeFast(pllInstance *tr, partitionList *pr)
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
                  *here[2];

                for(k = 0; k < 2; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 2 * rNumber) + width * k]);
                    here[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * pNumber) + width * k]);
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
                    
                    VECTOR_STORE((CAST)(&here[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
                    VECTOR_STORE((CAST)(&here[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);            
                  }
              }
              break;
            case 4:
              {
                parsimonyNumber
                  *left[4],
                  *right[4],
                  *here[4];

                for(k = 0; k < 4; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 4 * rNumber) + width * k]);
                    here[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * pNumber) + width * k]);
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
                    
                    VECTOR_STORE((CAST)(&here[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
                    VECTOR_STORE((CAST)(&here[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
                    VECTOR_STORE((CAST)(&here[2][i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
                    VECTOR_STORE((CAST)(&here[3][i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);  
                  }
              }
              break;
            case 20:
              {
                parsimonyNumber
                  *left[20],
                  *right[20],
                  *here[20];

                for(k = 0; k < 20; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 20 * rNumber) + width * k]);
                    here[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * pNumber) + width * k]);
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
                      VECTOR_STORE((CAST)(&here[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);
                  }
              }
              break;
            default:
              {
                parsimonyNumber
                  *left[32], 
                  *right[32],
                  *here[32];

                assert(states <= 32);
                
                for(k = 0; k < states; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * states * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * states * rNumber) + width * k]);
                    here[k]  = &(pr->partitionData[model]->parsVect[(width * states * pNumber) + width * k]);
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
                      VECTOR_STORE((CAST)(&here[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));                                                                    
                    
                    v_N = VECTOR_AND_NOT(v_N, allOne);
                    
                    totalScore += vectorPopcount(v_N);
                  }                             
              }
            }            
        }

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



static unsigned int evaluateParsimonyIterativeFast(pllInstance *tr, partitionList *pr)
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
    newviewParsimonyIterativeFast(tr, pr);

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
                 
                 if(sum >= bestScore)
                   return sum;                         
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
                 
                 if(sum >= bestScore)            
                   return sum;          
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
                  
                  if(sum >= bestScore)      
                    return sum;                        
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
                 
                 if(sum >= bestScore)         
                   return sum;                 
               }
           }
         }
    }
  
  return sum;
}


#else
static void newviewParsimonyIterativeFast(pllInstance *tr, partitionList * pr)
{    
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
                
                parsimonyNumber
                   o_A,
                   o_C,
                   t_A,
                   t_C, 
                   t_N;
                
                for(k = 0; k < 2; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 2 * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * pNumber) + width * k]);
                  }

                for(i = 0; i < width; i++)
                  {               
                    t_A = left[0][i] & right[0][i];
                    t_C = left[1][i] & right[1][i];                

                    o_A = left[0][i] | right[0][i];
                    o_C = left[1][i] | right[1][i];
                  
                    t_N = ~(t_A | t_C);   

                    this[0][i] = t_A | (t_N & o_A);
                    this[1][i] = t_C | (t_N & o_C);                
                    
                    totalScore += ((unsigned int) __builtin_popcount(t_N));
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

                parsimonyNumber
                   o_A,
                   o_C,
                   o_G,
                   o_T,
                   t_A,
                   t_C,
                   t_G,
                   t_T, 
                   t_N;

                for(i = 0; i < width; i++)
                  {               
                    t_A = left[0][i] & right[0][i];
                    t_C = left[1][i] & right[1][i];
                    t_G = left[2][i] & right[2][i];       
                    t_T = left[3][i] & right[3][i];

                    o_A = left[0][i] | right[0][i];
                    o_C = left[1][i] | right[1][i];
                    o_G = left[2][i] | right[2][i];       
                    o_T = left[3][i] | right[3][i];

                    t_N = ~(t_A | t_C | t_G | t_T);       

                    this[0][i] = t_A | (t_N & o_A);
                    this[1][i] = t_C | (t_N & o_C);
                    this[2][i] = t_G | (t_N & o_G);
                    this[3][i] = t_T | (t_N & o_T); 
                    
                    totalScore += ((unsigned int) __builtin_popcount(t_N));
                  }
              }
              break;
            case 20:
              {
                parsimonyNumber
                  *left[20],
                  *right[20],
                  *this[20];

                parsimonyNumber
                  o_A[20],
                  t_A[20],        
                  t_N;

                for(k = 0; k < 20; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * 20 * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * pNumber) + width * k]);
                  }

                for(i = 0; i < width; i++)
                  {               
                    size_t k;
                    
                    t_N = 0;

                    for(k = 0; k < 20; k++)
                      {
                        t_A[k] = left[k][i] & right[k][i];
                        o_A[k] = left[k][i] | right[k][i];
                        t_N = t_N | t_A[k];
                      }
                    
                    t_N = ~t_N;

                    for(k = 0; k < 20; k++)                   
                      this[k][i] = t_A[k] | (t_N & o_A[k]);                
                    
                    totalScore += ((unsigned int) __builtin_popcount(t_N));
                  }
              }
              break;
            default:
              {         
                parsimonyNumber
                  *left[32],
                  *right[32],
                  *this[32];
                
                parsimonyNumber
                  o_A[32],
                  t_A[32],        
                  t_N;
                
                assert(states <= 32);
                
                for(k = 0; k < states; k++)
                  {
                    left[k]  = &(pr->partitionData[model]->parsVect[(width * states * qNumber) + width * k]);
                    right[k] = &(pr->partitionData[model]->parsVect[(width * states * rNumber) + width * k]);
                    this[k]  = &(pr->partitionData[model]->parsVect[(width * states * pNumber) + width * k]);
                  }
                
                for(i = 0; i < width; i++)
                  {               
                    t_N = 0;
                    
                    for(k = 0; k < states; k++)
                      {
                        t_A[k] = left[k][i] & right[k][i];
                        o_A[k] = left[k][i] | right[k][i];
                        t_N = t_N | t_A[k];
                      }
                    
                    t_N = ~t_N;
                    
                    for(k = 0; k < states; k++)               
                      this[k][i] = t_A[k] | (t_N & o_A[k]);                
                    
                    totalScore += ((unsigned int) __builtin_popcount(t_N));
                  }
              }                       
            } 
        }

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}


static unsigned int evaluateParsimonyIterativeFast(pllInstance *tr, partitionList * pr)
{
  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr, pr); 

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
               t_A,
               t_C,           
               t_N,
               *left[2],
               *right[2];
             
             for(k = 0; k < 2; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * 2 * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * 2 * pNumber) + width * k]);
               }     
             
             for(i = 0; i < width; i++)
               {                                               
                 t_A = left[0][i] & right[0][i];
                 t_C = left[1][i] & right[1][i];
                 
                  t_N = ~(t_A | t_C);

                  sum += ((unsigned int) __builtin_popcount(t_N));
                 
                 if(sum >= bestScore)
                   return sum;                         
               }
           }
           break;
         case 4:
           {
             parsimonyNumber
               t_A,
               t_C,
               t_G,
               t_T,
               t_N,
               *left[4],
               *right[4];
      
             for(k = 0; k < 4; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * 4 * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * 4 * pNumber) + width * k]);
               }        

             for(i = 0; i < width; i++)
               {                                                
                  t_A = left[0][i] & right[0][i];
                  t_C = left[1][i] & right[1][i];
                  t_G = left[2][i] & right[2][i];         
                  t_T = left[3][i] & right[3][i];

                  t_N = ~(t_A | t_C | t_G | t_T);

                  sum += ((unsigned int) __builtin_popcount(t_N));
                 
                 if(sum >= bestScore)            
                   return sum;          
               }                 
           }
           break;
         case 20:
           {
             parsimonyNumber
               t_A,
               t_N,
               *left[20],
               *right[20];
             
              for(k = 0; k < 20; k++)
                {
                  left[k]  = &(pr->partitionData[model]->parsVect[(width * 20 * qNumber) + width * k]);
                  right[k] = &(pr->partitionData[model]->parsVect[(width * 20 * pNumber) + width * k]);
                }  
           
              for(i = 0; i < width; i++)
                { 
                  t_N = 0;
                  
                  for(k = 0; k < 20; k++)
                    {
                      t_A = left[k][i] & right[k][i];
                      t_N = t_N | t_A;
                    }
               
                  t_N = ~t_N;

                  sum += ((unsigned int) __builtin_popcount(t_N));
                  
                  if(sum >= bestScore)      
                    return sum;                        
                }
           }
           break;
         default:
           {
             parsimonyNumber
               t_A,
               t_N,
               *left[32], 
               *right[32];  

             assert(states <= 32);

             for(k = 0; k < states; k++)
               {
                 left[k]  = &(pr->partitionData[model]->parsVect[(width * states * qNumber) + width * k]);
                 right[k] = &(pr->partitionData[model]->parsVect[(width * states * pNumber) + width * k]);
               }  
           
             for(i = 0; i < width; i++)
               {                               
                 t_N = 0;
                  
                 for(k = 0; k < states; k++)
                   {
                     t_A = left[k][i] & right[k][i];
                     t_N = t_N | t_A;
                   }
               
                  t_N = ~t_N;

                  sum += ((unsigned int) __builtin_popcount(t_N));
                                                 
                 if(sum >= bestScore)                     
                   return sum;                     
               }                     
           }
         }
    }
  
  return sum;
}

#endif






static unsigned int evaluateParsimony(pllInstance *tr, partitionList *pr, nodeptr p, pllBoolean full)
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

  result = evaluateParsimonyIterativeFast(tr, pr);

  return result;
}


static void newviewParsimony(pllInstance *tr, partitionList *pr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, PLL_FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr, pr);
  }
}





/****************************************************************************************************************************************/

static void insertParsimony (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q);
  hookupDefault(p->next->next, r);
   
  newviewParsimony(tr, pr, p);
} 



static nodeptr buildNewTip (pllInstance *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTree (pllInstance *tr, partitionList *pr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = PLL_MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq]);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, pr, s, p);
}


static void testInsertParsimony (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, pllBoolean saveBranches)
{ 
  unsigned int 
    mp;
 
  nodeptr  
    r = q->back;   

  pllBoolean
    doIt = PLL_TRUE;

  int numBranches = pr->perGeneBranchLengths?pr->numberOfPartitions:1;

  if(tr->grouped)
    {
      int 
        rNumber = tr->constraintVector[r->number],
        qNumber = tr->constraintVector[q->number],
        pNumber = tr->constraintVector[p->number];

      doIt = PLL_FALSE;
     
      if(pNumber == -9)
        pNumber = checkerPars(tr, p->back);
      if(pNumber == -9)
        doIt = PLL_TRUE;
      else
        {
          if(qNumber == -9)
            qNumber = checkerPars(tr, q);

          if(rNumber == -9)
            rNumber = checkerPars(tr, r);

          if(pNumber == rNumber || pNumber == qNumber)
            doIt = PLL_TRUE;       
        }
    }

  if(doIt)
    {
      double* z = (double*)rax_malloc(numBranches*sizeof(double));
      
      if(saveBranches)
        {
          int i;
          
          for(i = 0; i < numBranches; i++)
            z[i] = q->z[i];
        }

      insertParsimony(tr, pr, p, q);
  
      mp = evaluateParsimony(tr, pr, p->next->next, PLL_FALSE);

      if(mp < tr->bestParsimony)
        {
          tr->bestParsimony = mp;
          tr->insertNode = q;
          tr->removeNode = p;
        }
      
      if(saveBranches)
        hookup(q, r, z, numBranches);
      else
        hookupDefault(q, r);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
      rax_free(z);
    }
       
  return;
} 


static void restoreTreeParsimony(pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupDefault(p->next,       q);
  hookupDefault(p->next->next, r);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, PLL_FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr, pr);
}


static void addTraverseParsimony (pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q, int mintrav, int maxtrav, pllBoolean doAll, pllBoolean saveBranches)
{        
  if (doAll || (--mintrav <= 0))               
    testInsertParsimony(tr, pr, p, q, saveBranches);

  if (((q->number > tr->mxtips)) && ((--maxtrav > 0) || doAll))
    {         
      addTraverseParsimony(tr, pr, p, q->next->back, mintrav, maxtrav, doAll, saveBranches);
      addTraverseParsimony(tr, pr, p, q->next->next->back, mintrav, maxtrav, doAll, saveBranches);
    }
}





static void makePermutationFast(int *perm, int n, pllInstance *tr)
{    
  int  
    i, 
    j, 
    k;

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum(&tr->randomNumberSeed);

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

//static nodeptr  removeNodeParsimony (nodeptr p, tree *tr)
static nodeptr  removeNodeParsimony (nodeptr p)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupDefault(q, r);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}

static int rearrangeParsimony(pllInstance *tr, partitionList *pr, nodeptr p, int mintrav, int maxtrav, pllBoolean doAll)
{   
  nodeptr  
    p1, 
    p2, 
    q, 
    q1, 
    q2;
  
  int      
    mintrav2; 

  pllBoolean
    doP = PLL_TRUE,
    doQ = PLL_TRUE;
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  assert(mintrav == 1);

  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(tr->constrained)
    {    
      if(! tipHomogeneityCheckerPars(tr, p->back, 0))
        doP = PLL_FALSE;
        
      if(! tipHomogeneityCheckerPars(tr, q->back, 0))
        doQ = PLL_FALSE;
                        
      if(doQ == PLL_FALSE && doP == PLL_FALSE)
        return 0;
    }  

  if((p->number > tr->mxtips) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if ((p1->number > tr->mxtips) || (p2->number > tr->mxtips)) 
        {                 
          //removeNodeParsimony(p, tr);          
          removeNodeParsimony(p);                

          if ((p1->number > tr->mxtips)) 
            {
              addTraverseParsimony(tr, pr, p, p1->next->back, mintrav, maxtrav, doAll, PLL_FALSE);
              addTraverseParsimony(tr, pr, p, p1->next->next->back, mintrav, maxtrav, doAll, PLL_FALSE);
            }
         
          if ((p2->number > tr->mxtips)) 
            {
              addTraverseParsimony(tr, pr, p, p2->next->back, mintrav, maxtrav, doAll, PLL_FALSE);
              addTraverseParsimony(tr, pr, p, p2->next->next->back, mintrav, maxtrav, doAll, PLL_FALSE);
            }
            
           
          hookupDefault(p->next,       p1);
          hookupDefault(p->next->next, p2);

          newviewParsimony(tr, pr, p);
        }
    }  
       
  if ((q->number > tr->mxtips) && (maxtrav > 0) && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

      if (
          (
           (q1->number > tr->mxtips) && 
           ((q1->next->back->number > tr->mxtips) || (q1->next->next->back->number > tr->mxtips))
           )
          ||
          (
           (q2->number > tr->mxtips) && 
           ((q2->next->back->number > tr->mxtips) || (q2->next->next->back->number > tr->mxtips))
           )
          )
        {          

          //removeNodeParsimony(q, tr);
          removeNodeParsimony(q);
          
          mintrav2 = mintrav > 2 ? mintrav : 2;
          
          if ((q1->number > tr->mxtips)) 
            {
              addTraverseParsimony(tr, pr, q, q1->next->back, mintrav2 , maxtrav, doAll, PLL_FALSE);
              addTraverseParsimony(tr, pr, q, q1->next->next->back, mintrav2 , maxtrav, doAll, PLL_FALSE);
            }
         
          if ((q2->number > tr->mxtips)) 
            {
              addTraverseParsimony(tr, pr, q, q2->next->back, mintrav2 , maxtrav, doAll, PLL_FALSE);
              addTraverseParsimony(tr, pr, q, q2->next->next->back, mintrav2 , maxtrav, doAll, PLL_FALSE);
            }      
           
          hookupDefault(q->next,       q1);
          hookupDefault(q->next->next, q2);
           
          newviewParsimony(tr, pr, q);
        }
    }

  return 1;
} 


static void restoreTreeRearrangeParsimony(pllInstance *tr, partitionList *pr)
{    
  removeNodeParsimony(tr->removeNode);  
  //removeNodeParsimony(tr->removeNode, tr);  
  restoreTreeParsimony(tr, pr, tr->removeNode, tr->insertNode);
}

/*
static pllBoolean isInformative2(pllInstance *tr, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = 15;

  unsigned char
    nucleotide,
    target = 0;
        
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {      
      nucleotide = tr->yVector[j][site];            
      check[nucleotide] =  check[nucleotide] + 1;                  
    }
  
  
  if(check[1] > 1)
    {
      informativeCounter++;    
      target = target | 1;
    }
  if(check[2] > 1)
    {
      informativeCounter++; 
      target = target | 2;
    }
  if(check[4] > 1)
    {
      informativeCounter++; 
      target = target | 4;
    }
  if(check[8] > 1)
    {
      informativeCounter++; 
      target = target | 8;
    }
          
  if(informativeCounter >= 2)
    return PLL_TRUE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
        {
          if(j == 3 || j == 5 || j == 6 || j == 7 || j == 9 || j == 10 || j == 11 || 
             j == 12 || j == 13 || j == 14)
            {
              if(check[j] > 1)
                {
                  if(!(target & j))
                    return PLL_TRUE;
                }
            }
        } 
    }
     
  return PLL_FALSE;          
}
*/

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
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);                   
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
        informativeCounter++;    
    } 
          
  if(informativeCounter <= 1)
    return PLL_FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
        {
          if(check[j] > 1)
            return PLL_TRUE;
        } 
    }
     
  return PLL_FALSE;          
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


static void reorderNodes(pllInstance *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if((p->number <= tr->mxtips))    
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



static void nodeRectifierPars(pllInstance *tr)
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


  
static void compressDNA(pllInstance *tr, partitionList *pr, int *informative)
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



static void stepwiseAddition(pllInstance *tr, partitionList *pr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, PLL_FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;
    
  mp = evaluateParsimonyIterativeFast(tr, pr);
  
  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {         
      stepwiseAddition(tr, pr, p, q->next->back);
      stepwiseAddition(tr, pr, p, q->next->next->back);
    }
}



void allocateParsimonyDataStructures(pllInstance *tr, partitionList *pr)
{
  int 
    i,
    *informative = (int *)rax_malloc(sizeof(int) * (size_t)tr->originalCrunchedLength);
 
  determineUninformativeSites(tr, pr, informative);

  compressDNA(tr, pr, informative);

  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 1; i++)
    {
      nodeptr 
        p = tr->nodep[i];

      p->xPars = 1;
      p->next->xPars = 0;
      p->next->next->xPars = 0;
    }

  tr->ti = (int*)rax_malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  rax_free(informative); 
}

void pllFreeParsimonyDataStructures(pllInstance *tr, partitionList *pr)
{
  size_t 
    model;

  rax_free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) pr->numberOfPartitions; ++model)
    rax_free(pr->partitionData[model]->parsVect);
  
  rax_free(tr->ti);
}


void pllMakeParsimonyTreeFast(pllInstance *tr, partitionList *pr, int sprDist)
{   
  nodeptr  
    p, 
    f;    

  int 
    i, 
    nextsp,
    *perm        = (int *)rax_malloc((size_t)(tr->mxtips + 1) * sizeof(int));  

  unsigned int 
    randomMP, 
    startMP;         
  
  assert(!tr->constrained);

  makePermutationFast(perm, tr->mxtips, tr);
  
  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
  
  buildSimpleTree(tr, pr, perm[1], perm[2], perm[3]);
  
  f = tr->start;       
  
  while(tr->ntips < tr->mxtips) 
    {   
      nodeptr q;
      
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped)
        {
          int 
            number = p->back->number;            

          tr->constraintVector[number] = -9;
        }
          
      stepwiseAddition(tr, pr, q, f->back);
      
      {
        nodeptr   
          r = tr->insertNode->back;
        
        int counter = 4;
        
        hookupDefault(q->next,       tr->insertNode);
        hookupDefault(q->next->next, r);
        
        computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, PLL_FALSE);              
        tr->ti[0] = counter;
        
        newviewParsimonyIterativeFast(tr, pr);
      }
    }    
  
  nodeRectifierPars(tr);
  
  randomMP = tr->bestParsimony;        
  
  do
    {
      startMP = randomMP;
      nodeRectifierPars(tr);
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
        {
          rearrangeParsimony(tr, pr, tr->nodep[i], 1, sprDist, PLL_FALSE);
          if(tr->bestParsimony < randomMP)
            {           
              restoreTreeRearrangeParsimony(tr, pr);
              randomMP = tr->bestParsimony;
            }
        }                          
    }
  while(randomMP < startMP);
  
  rax_free(perm);
} 
