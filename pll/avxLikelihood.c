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
 * @file avxLikelihood.c
 *
 * @brief AVX versions of the likelihood functions
 *
 * AVX versions of the likelihood functions
 */
#include "systypes.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <stdint.h>

#if defined(__ARM_NEON)
#include "utils/sse2neon.h"
#else
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#endif

#include <assert.h>

#ifdef _FMA
#include <x86intrin.h>
#define FMAMACC(a,b,c) _mm256_fmadd_pd(b,c,a)
#endif

#include "pll.h"
#include "pllInternal.h"

extern const unsigned int mask32[32];

PLL_ALIGN_BEGIN const union PLL_ALIGN_END
{
  uint64_t i[4];
  __m256d m;
  
} absMask_AVX = {{0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL}};



static __inline __m256d hadd4(__m256d v, __m256d u)
{ 
  __m256d
    a, b;
  
  v = _mm256_hadd_pd(v, v);
  a = _mm256_permute2f128_pd(v, v, 1);
  v = _mm256_add_pd(a, v);

  u = _mm256_hadd_pd(u, u);
  b = _mm256_permute2f128_pd(u, u, 1);
  u = _mm256_add_pd(b, u);

  v = _mm256_mul_pd(v, u);	
  
  return v;
}

static __inline __m256d hadd3(__m256d v)
{ 
  __m256d
    a;
  
  v = _mm256_hadd_pd(v, v);
  a = _mm256_permute2f128_pd(v, v, 1);
  v = _mm256_add_pd(a, v);
  
  return v;
}


void  newviewGTRGAMMA_AVX(int tipCase,
			 double *x1, double *x2, double *x3,
			 double *extEV, double *tipVector,
			 int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			 const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling
			 )
{
 
  int  
    i, 
    k, 
    scale, 
    addScale = 0;
 
  __m256d 
    minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD),
    twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
 

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	double 
	  *uX1, *uX2;
	PLL_ALIGN_BEGIN double
	  umpX1[1024] PLL_ALIGN_END,
	  umpX2[1024] PLL_ALIGN_END;

	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i * 4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&right[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX2[i * 64 + j * 16 + k * 4], left1);
		}	    
	  }   	
	  

	for(i = 0; i < n; i++)
	  {	    		 	    
	    uX1 = &umpX1[64 * tipX1[i]];
	    uX2 = &umpX2[64 * tipX2[i]];		  
	    
	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		  xv = _mm256_setzero_pd();
	       
		int 
		  l;
		
		for(l = 0; l < 4; l++)
		  {	       	     				      	      																	   
		    __m256d
		      x1v =  _mm256_mul_pd(_mm256_load_pd(&uX1[k * 16 + l * 4]), _mm256_load_pd(&uX2[k * 16 + l * 4]));
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
#ifdef _FMA
		    xv = FMAMACC(xv,x1v,evv);
#else						  
		    xv = _mm256_add_pd(xv, _mm256_mul_pd(x1v, evv));
#endif
		  }
		
		_mm256_store_pd(&x3[16 * i + 4 * k], xv);
	      }	         	   	    
	  }
      }
      break;
    case PLL_TIP_INNER:
      {
	double 
	  *uX1;
	PLL_ALIGN_BEGIN double
	  umpX1[1024] PLL_ALIGN_END;

	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i*4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}	 	   
	  }   	
	
	for(i = 0; i < n; i++)
	  { 
	    __m256d
	      xv[4];	    	   
	    
	    scale = 1;
	    uX1 = &umpX1[64 * tipX1[i]];

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   		 
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k]  = _mm256_setzero_pd();
		  
		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d  
		      x1v = _mm256_load_pd(&uX1[k * 16 + l * 4]),		     
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x2v = hadd3(x2v);
		    x1v = _mm256_mul_pd(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
			
#ifdef _FMA
		    xv[k] = FMAMACC(xv[k],x1v,evv);
#else			  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
#endif
		  }
		    
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }	    

	    if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i] += 1;
	      }

	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    case PLL_INNER_INNER:
      {
	for(i = 0; i < n; i++)
	  {	
	    __m256d
	      xv[4];
	    
	    scale = 1;

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		 
		  xvl = _mm256_load_pd(&(x1[i * 16 + k * 4])),
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k] = _mm256_setzero_pd();

		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d 
		      x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x1v = hadd4(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
						  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		  }
		
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }

	     if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i] += 1;		
	      }
		
	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
  
}

void  newviewGTRGAMMA_AVX_GAPPED_SAVE(int tipCase,
				      double *x1_start, double *x2_start, double *x3_start,
				      double *extEV, double *tipVector,
				      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				      const int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
				      unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
				      double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
				      )
{
 
  int  
    i, 
    k, 
    scale,
    scaleGap,
    addScale = 0;
 
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD ),
    twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
 
  double
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	double 
	  *uX1, *uX2;
	PLL_ALIGN_BEGIN double
	  umpX1[1024] PLL_ALIGN_END,
	  umpX2[1024] PLL_ALIGN_END;

	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i * 4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&right[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX2[i * 64 + j * 16 + k * 4], left1);
		}	    
	  }   	
	  
	x3 = x3_gapColumn;

	{
	  uX1 = &umpX1[960];
	  uX2 = &umpX2[960];		  
	  
	  for(k = 0; k < 4; k++)
	    {
	      __m256d	   
		xv = _mm256_setzero_pd();
	      
	      int 
		l;
	      
	      for(l = 0; l < 4; l++)
		{	       	     				      	      																	   
		  __m256d
		    x1v =  _mm256_mul_pd(_mm256_load_pd(&uX1[k * 16 + l * 4]), _mm256_load_pd(&uX2[k * 16 + l * 4]));
		  
		  __m256d 
		    evv = _mm256_load_pd(&extEV[l * 4]);
#ifdef _FMA
		  xv = FMAMACC(xv,x1v,evv);
#else						  
		  xv = _mm256_add_pd(xv, _mm256_mul_pd(x1v, evv));
#endif
		}
		    
	      _mm256_store_pd(&x3[4 * k], xv);
	    }
	}
	
	x3 = x3_start;

	for(i = 0; i < n; i++)
	  {		    	    	
	    if(!(x3_gap[i / 32] & mask32[i % 32]))	     
	      {
		uX1 = &umpX1[64 * tipX1[i]];
		uX2 = &umpX2[64 * tipX2[i]];		  
	    
		for(k = 0; k < 4; k++)
		  {
		    __m256d	   
		      xv = _mm256_setzero_pd();
	       
		    int 
		      l;
		
		    for(l = 0; l < 4; l++)
		      {	       	     				      	      																	   
			__m256d
			  x1v =  _mm256_mul_pd(_mm256_load_pd(&uX1[k * 16 + l * 4]), _mm256_load_pd(&uX2[k * 16 + l * 4]));
			
			__m256d 
			  evv = _mm256_load_pd(&extEV[l * 4]);
#ifdef _FMA
			xv = FMAMACC(xv,x1v,evv);
#else						  
			xv = _mm256_add_pd(xv, _mm256_mul_pd(x1v, evv));
#endif
		      }
		    
		    _mm256_store_pd(&x3[4 * k], xv);
		  }

		x3 += 16;
	      }
	  }
      }
      break;
    case PLL_TIP_INNER:
      {
	double 
	  *uX1;
	PLL_ALIGN_BEGIN double
	  umpX1[1024] PLL_ALIGN_END;
       
	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i*4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}	 	   
	  }	

	{ 
	  __m256d
	    xv[4];
	  
	  scaleGap = 1;
	  uX1 = &umpX1[960];

	  x2 = x2_gapColumn;			 
	  x3 = x3_gapColumn;

	  for(k = 0; k < 4; k++)
	    {
	      __m256d	   		 
		xvr = _mm256_load_pd(&(x2[k * 4]));

	      int 
		l;

	      xv[k]  = _mm256_setzero_pd();
		  
	      for(l = 0; l < 4; l++)
		{	       	     				      	      															
		  __m256d  
		    x1v = _mm256_load_pd(&uX1[k * 16 + l * 4]),		     
		    x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		  x2v = hadd3(x2v);
		  x1v = _mm256_mul_pd(x1v, x2v);			
		
		  __m256d 
		    evv = _mm256_load_pd(&extEV[l * 4]);
			
#ifdef _FMA
		  xv[k] = FMAMACC(xv[k],x1v,evv);
#else			  
		  xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
#endif
		}
		    
	      if(scaleGap)
		{
		  __m256d 	     
		    v1 = _mm256_and_pd(xv[k], absMask_AVX.m);
		  
		  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		  if(_mm256_movemask_pd( v1 ) != 15)
		    scaleGap = 0;
		}
	    }
	
	  if(scaleGap)
	    {
	      xv[0] = _mm256_mul_pd(xv[0], twoto);
	      xv[1] = _mm256_mul_pd(xv[1], twoto);
	      xv[2] = _mm256_mul_pd(xv[2], twoto);
	      xv[3] = _mm256_mul_pd(xv[3], twoto);	    
	    }

	  _mm256_store_pd(&x3[0],      xv[0]);
	  _mm256_store_pd(&x3[4],  xv[1]);
	  _mm256_store_pd(&x3[8],  xv[2]);
	  _mm256_store_pd(&x3[12], xv[3]);
	}
	
	x3 = x3_start;
	
	for(i = 0; i < n; i++)
	  {
	    if((x3_gap[i / 32] & mask32[i % 32]))
	      {
		if(scaleGap)
		  {
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;
		  }
	      }
	    else
	      {
		if(x2_gap[i / 32] & mask32[i % 32])
		  x2 = x2_gapColumn;
		else
		  {
		    x2 = x2_ptr;
		    x2_ptr += 16;
		  }
		
		__m256d
		  xv[4];	    	   
		
		scale = 1;
		uX1 = &umpX1[64 * tipX1[i]];
		
		for(k = 0; k < 4; k++)
		  {
		    __m256d	   		 
		      xvr = _mm256_load_pd(&(x2[k * 4]));
		    
		    int 
		      l;
		    
		    xv[k]  = _mm256_setzero_pd();
		    
		    for(l = 0; l < 4; l++)
		      {	       	     				      	      															
			__m256d  
			  x1v = _mm256_load_pd(&uX1[k * 16 + l * 4]),		     
			  x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
			x2v = hadd3(x2v);
			x1v = _mm256_mul_pd(x1v, x2v);			
			
			__m256d 
			  evv = _mm256_load_pd(&extEV[l * 4]);
			
#ifdef _FMA
			xv[k] = FMAMACC(xv[k],x1v,evv);
#else			  
			xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
#endif
		      }
		    
		    if(scale)
		      {
			__m256d 	     
			  v1 = _mm256_and_pd(xv[k], absMask_AVX.m);
			
			v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
			
			if(_mm256_movemask_pd( v1 ) != 15)
			  scale = 0;
		      }
		  }	    
	      
		if(scale)
		  {
		    xv[0] = _mm256_mul_pd(xv[0], twoto);
		    xv[1] = _mm256_mul_pd(xv[1], twoto);
		    xv[2] = _mm256_mul_pd(xv[2], twoto);
		    xv[3] = _mm256_mul_pd(xv[3], twoto);

		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i] += 1;		   
		  }
	      
		_mm256_store_pd(&x3[0],      xv[0]);
		_mm256_store_pd(&x3[4],  xv[1]);
		_mm256_store_pd(&x3[8],  xv[2]);
		_mm256_store_pd(&x3[12], xv[3]);
	      
		x3 += 16;
	      }
	  }
      }
      break;
    case PLL_INNER_INNER:
      {          
	{		
	  x1 = x1_gapColumn;	     	    
	  x2 = x2_gapColumn;	    
	  x3 = x3_gapColumn;

	  __m256d
	    xv[4];
	    
	  scaleGap = 1;

	  for(k = 0; k < 4; k++)
	    {
	      __m256d	   
		
		xvl = _mm256_load_pd(&(x1[k * 4])),
		xvr = _mm256_load_pd(&(x2[k * 4]));

	      int 
		l;

	      xv[k] = _mm256_setzero_pd();

	      for(l = 0; l < 4; l++)
		{	       	     				      	      															
		  __m256d 
		    x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		    x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
		  
		  x1v = hadd4(x1v, x2v);			
		  
		  __m256d 
		    evv = _mm256_load_pd(&extEV[l * 4]);
		  
		  xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		}
		
	      if(scaleGap)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scaleGap = 0;
		  }
	    }

	  if(scaleGap)
	    {
	      xv[0] = _mm256_mul_pd(xv[0], twoto);
	      xv[1] = _mm256_mul_pd(xv[1], twoto);
	      xv[2] = _mm256_mul_pd(xv[2], twoto);
	      xv[3] = _mm256_mul_pd(xv[3], twoto);	       
	    }
		
	  _mm256_store_pd(&x3[0],  xv[0]);
	  _mm256_store_pd(&x3[4],  xv[1]);
	  _mm256_store_pd(&x3[8],  xv[2]);
	  _mm256_store_pd(&x3[12], xv[3]);
	}	  
      
	x3 = x3_start;

	for(i = 0; i < n; i++)
	  {
	    if(x3_gap[i / 32] & mask32[i % 32])
	      {	     
		if(scaleGap)
		  {
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1; 	       
		  }
	      }
	    else
	      {	
		if(x1_gap[i / 32] & mask32[i % 32])
		  x1 = x1_gapColumn;
		else
		  {
		    x1 = x1_ptr;
		    x1_ptr += 16;
		  }
	     
		if(x2_gap[i / 32] & mask32[i % 32])
		  x2 = x2_gapColumn;
		else
		  {
		    x2 = x2_ptr;
		    x2_ptr += 16;
		  }

		__m256d
		  xv[4];
	    
		scale = 1;

		for(k = 0; k < 4; k++)
		  {
		    __m256d	   
		      
		      xvl = _mm256_load_pd(&(x1[k * 4])),
		      xvr = _mm256_load_pd(&(x2[k * 4]));
		    
		    int 
		      l;
		    
		    xv[k] = _mm256_setzero_pd();
		    
		    for(l = 0; l < 4; l++)
		      {	       	     				      	      															
			__m256d 
			  x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
			  x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
			x1v = hadd4(x1v, x2v);			
			
			__m256d 
			  evv = _mm256_load_pd(&extEV[l * 4]);
			
			xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		      }
		    
		    if(scale)
		      {
			__m256d 	     
			  v1 = _mm256_and_pd(xv[k], absMask_AVX.m);
			
			v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
			
			if(_mm256_movemask_pd( v1 ) != 15)
			  scale = 0;
		      }
		  }

		if(scale)
		  {
		    xv[0] = _mm256_mul_pd(xv[0], twoto);
		    xv[1] = _mm256_mul_pd(xv[1], twoto);
		    xv[2] = _mm256_mul_pd(xv[2], twoto);
		    xv[3] = _mm256_mul_pd(xv[3], twoto);
		    
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i] += 1;
		  }
		
		_mm256_store_pd(&x3[0],      xv[0]);
		_mm256_store_pd(&x3[4],  xv[1]);
		_mm256_store_pd(&x3[8],  xv[2]);
		_mm256_store_pd(&x3[12], xv[3]);
	      
		x3 += 16;
	      }
	  }
      }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
  
}




void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
			   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1,
    *x2;
    
  int 
    i, 
    addScale = 0;
   
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD ),
    twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
  
  switch(tipCase)
    {
    case PLL_TIP_TIP:      
      for (i = 0; i < n; i++)
	{	 
	  int 
	    l;
	  
	  le = &left[cptr[i] * 16];
	  ri = &right[cptr[i] * 16];

	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  __m256d	   
	    vv = _mm256_setzero_pd();
	   	   	    
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else				
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
	    }	  		  

	  _mm256_store_pd(&x3_start[4 * i], vv);	    	   	    
	}
      break;
    case PLL_TIP_INNER:      
      for (i = 0; i < n; i++)
	{
	  int 
	    l;

	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];	 
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m256d	   
	    vv = _mm256_setzero_pd();
	  
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
				
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else	      
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));
#endif
	    }	  		  
	  
	  
	  __m256d 	     
	    v1 = _mm256_and_pd(vv, absMask_AVX.m);

	  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	    
	  if(_mm256_movemask_pd( v1 ) == 15)
	    {	     	      
	      vv = _mm256_mul_pd(vv, twoto);	      
	      
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i] += 1;	      	     
	    }       
	  
	  _mm256_store_pd(&x3_start[4 * i], vv);	 	  	  
	}
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  int 
	    l;

	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m256d	   
	    vv = _mm256_setzero_pd();
	  
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else						
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
	    }	  		  

	 
	  __m256d 	     
	    v1 = _mm256_and_pd(vv, absMask_AVX.m);

	  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	    
	  if(_mm256_movemask_pd( v1 ) == 15)
	    {	
	      vv = _mm256_mul_pd(vv, twoto);
	      
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i] += 1;	   
	    }	

	  _mm256_store_pd(&x3_start[4 * i], vv);
	  	  
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
}


void newviewGTRCAT_AVX_GAPPED_SAVE(int tipCase,  double *EV,  int *cptr,
				   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le,
    *ri,
    *x1,
    *x2, 
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start, 
    *x3_ptr = x3_start;
  
  int 
    i, 
    scaleGap = 0,
    addScale = 0;
   
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD ),
    twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
  

  {
    int 
      l;

    x1 = x1_gapColumn;	      
    x2 = x2_gapColumn;
    x3 = x3_gapColumn;    	 
	  	  
    le =  &left[maxCats * 16];
    ri =  &right[maxCats * 16];

    __m256d	   
      vv = _mm256_setzero_pd();
	  
    for(l = 0; l < 4; l++)
      {	       	     				      	      															
	__m256d 
	  x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
	  x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
	
	x1v = hadd4(x1v, x2v);			
	
	__m256d 
	  evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
	vv = FMAMACC(vv,x1v,evv);
#else						
	vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
      }	  		  

    if(tipCase != PLL_TIP_TIP)
      {
	__m256d 	     
	  v1 = _mm256_and_pd(vv, absMask_AVX.m);
    
	v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
    
	if(_mm256_movemask_pd( v1 ) == 15)
	  {
	    vv = _mm256_mul_pd(vv, twoto);	      	 
	    scaleGap = 1;
	  }
      }
    
    _mm256_store_pd(x3, vv);    
  }

  switch(tipCase)
    {
    case PLL_TIP_TIP:      
      for (i = 0; i < n; i++)
	{ 
	  if(noGap(x3_gap, i))
	    {	 
	      int 
		l;
	      
	      x1 = &(tipVector[4 * tipX1[i]]);
	      x2 = &(tipVector[4 * tipX2[i]]);

	      x3 = x3_ptr;

	      if(isGap(x1_gap, i))
		le =  &left[maxCats * 16];
	      else	  	  
		le =  &left[cptr[i] * 16];	  
	      
	      if(isGap(x2_gap, i))
		ri =  &right[maxCats * 16];
	      else	 	  
		ri =  &right[cptr[i] * 16];
	  	  
	      __m256d	   
		vv = _mm256_setzero_pd();
	      
	      for(l = 0; l < 4; l++)
		{	       	     				      	      															
		  __m256d 
		    x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		    x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
		  
		  x1v = hadd4(x1v, x2v);			
		  
		  __m256d 
		    evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
		  vv = FMAMACC(vv,x1v,evv);
#else				
		  vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
		}	  		  

	      _mm256_store_pd(x3, vv);	 
	      
	      x3_ptr += 4;
	    }
	}
      break;
    case PLL_TIP_INNER:      
      for (i = 0; i < n; i++)
	{ 
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)
		{
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i] += 1;		   		    
		}	       
	    }
	  else
	    {
	      int 
		l;

	      x1 = &(tipVector[4 * tipX1[i]]);    
	      x3 = x3_ptr;

	      if(isGap(x1_gap, i))
		le =  &left[maxCats * 16];
	      else
		le =  &left[cptr[i] * 16];
	  
	      if(isGap(x2_gap, i))
		{		 
		  ri =  &right[maxCats * 16];
		  x2 = x2_gapColumn;
		}
	      else
		{
		  ri =  &right[cptr[i] * 16];
		  x2 = x2_ptr;
		  x2_ptr += 4;
		}	  	 

	      __m256d	   
		vv = _mm256_setzero_pd();
	      
	      for(l = 0; l < 4; l++)
		{	       	     				      	      															
		  __m256d 
		    x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		    x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
		  
		  x1v = hadd4(x1v, x2v);			
		  
		  __m256d 
		    evv = _mm256_load_pd(&EV[l * 4]);
		  
#ifdef _FMA
		  vv = FMAMACC(vv,x1v,evv);
#else	      
		  vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));
#endif
		}	  		  
	  
	  
	      __m256d 	     
		v1 = _mm256_and_pd(vv, absMask_AVX.m);
	      
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) == 15)
		{	     	      
		  vv = _mm256_mul_pd(vv, twoto);	      
		  
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i] += 1;		 
		}       
	  
	      _mm256_store_pd(x3, vv);	 	  	  

	      x3_ptr += 4;
	    }
	}
      break;
    case PLL_INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)		   		    
		{
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i] += 1;
		}	      
	    }
	  else
	    {
	      int 
		l;
	      
	      x3 = x3_ptr;
	      
	      if(isGap(x1_gap, i))
		{
		  x1 = x1_gapColumn;
		  le =  &left[maxCats * 16];
		}
	      else
		{
		  le =  &left[cptr[i] * 16];
		  x1 = x1_ptr;
		  x1_ptr += 4;
		}

	      if(isGap(x2_gap, i))	
		{
		  x2 = x2_gapColumn;
		  ri =  &right[maxCats * 16];	    
		}
	      else
		{
		  ri =  &right[cptr[i] * 16];
		  x2 = x2_ptr;
		  x2_ptr += 4;
		}	 	  	  	  
	  
	      __m256d	   
		vv = _mm256_setzero_pd();
	      
	      for(l = 0; l < 4; l++)
		{	       	     				      	      															
		  __m256d 
		    x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		    x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
		  
		  x1v = hadd4(x1v, x2v);			
		  
		  __m256d 
		    evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
		  vv = FMAMACC(vv,x1v,evv);
#else						
		  vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
		}	  		  
	      
	      
	      __m256d 	     
		v1 = _mm256_and_pd(vv, absMask_AVX.m);
	      
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) == 15)
		{	
		  vv = _mm256_mul_pd(vv, twoto);	      
		  
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i] += 1;		
		}	
	      
	      _mm256_store_pd(x3, vv);
	      
	      x3_ptr += 4;
	    }	  	  
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
}

void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
			       int *cptr,
			       double *x1, double *x2, double *x3, double *tipVector,
			       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			       int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;

  int i, l, scale, addScale = 0;

#ifdef _FMA
  int k;
#endif

  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {	   
	    le = &left[cptr[i] * 400];
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);
	    vr = &(tipVector[20 * tipX2[i]]);
	    v  = &x3[20 * i];	    	    	   	    

	    __m256d vv[5];
	    
	    vv[0] = _mm256_setzero_pd();
	    vv[1] = _mm256_setzero_pd();
	    vv[2] = _mm256_setzero_pd();
	    vv[3] = _mm256_setzero_pd();
	    vv[4] = _mm256_setzero_pd();	   	    

	    for(l = 0; l < 20; l++)
	      {	       
		__m256d 
		  x1v = _mm256_setzero_pd(),
		  x2v = _mm256_setzero_pd();	
				
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];														

#ifdef _FMA		
		for(k = 0; k < 20; k += 4) 
		  {
		    __m256d vlv = _mm256_load_pd(&vl[k]);
		    __m256d lvv = _mm256_load_pd(&lv[k]);
		    x1v = FMAMACC(x1v,vlv,lvv);
		    __m256d vrv = _mm256_load_pd(&vr[k]);
		    __m256d rvv = _mm256_load_pd(&rv[k]);
		    x2v = FMAMACC(x2v,vrv,rvv);
		  }
#else		
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));

		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));	
#endif

		x1v = hadd4(x1v, x2v);			
#ifdef _FMA
		for(k = 0; k < 5; k++) 
		  {
		    __m256d evv = _mm256_load_pd(&ev[k*4]);
		    vv[k] = FMAMACC(vv[k],x1v,evv);
		  }	  
#else		
		__m256d 
		  evv[5];
	    	
		evv[0] = _mm256_load_pd(&ev[0]);
		evv[1] = _mm256_load_pd(&ev[4]);
		evv[2] = _mm256_load_pd(&ev[8]);
		evv[3] = _mm256_load_pd(&ev[12]);
		evv[4] = _mm256_load_pd(&ev[16]);		
		
		vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
		vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
		vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
		vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
		vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      		      	  
#endif
	      }
	    _mm256_store_pd(&v[0], vv[0]);
	    _mm256_store_pd(&v[4], vv[1]);
	    _mm256_store_pd(&v[8], vv[2]);
	    _mm256_store_pd(&v[12], vv[3]);
	    _mm256_store_pd(&v[16], vv[4]);
	  }
      }
      break;
    case PLL_TIP_INNER:      	
      for (i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 400];
	  ri = &right[cptr[i] * 400];
	  
	  vl = &(tipVector[20 * tipX1[i]]);
	  vr = &x2[20 * i];
	  v  = &x3[20 * i];	   
	  
	  __m256d vv[5];
	  
	  vv[0] = _mm256_setzero_pd();
	  vv[1] = _mm256_setzero_pd();
	  vv[2] = _mm256_setzero_pd();
	  vv[3] = _mm256_setzero_pd();
	  vv[4] = _mm256_setzero_pd();
	  
	 

	  for(l = 0; l < 20; l++)
	    {	       
	      __m256d 
		x1v = _mm256_setzero_pd(),
		x2v = _mm256_setzero_pd();	
	      
	      double 
		*ev = &extEV[l * 20],
		*lv = &le[l * 20],
		*rv = &ri[l * 20];														
#ifdef _FMA
	      for(k = 0; k < 20; k += 4) 
		{
		  __m256d vlv = _mm256_load_pd(&vl[k]);
		  __m256d lvv = _mm256_load_pd(&lv[k]);
		  x1v = FMAMACC(x1v,vlv,lvv);
		  __m256d vrv = _mm256_load_pd(&vr[k]);
		  __m256d rvv = _mm256_load_pd(&rv[k]);
		  x2v = FMAMACC(x2v,vrv,rvv);
		}
#else	      
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
	      
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));
#endif

	      x1v = hadd4(x1v, x2v);			
	      
	      __m256d 
		evv[5];
	      
	      evv[0] = _mm256_load_pd(&ev[0]);
	      evv[1] = _mm256_load_pd(&ev[4]);
	      evv[2] = _mm256_load_pd(&ev[8]);
	      evv[3] = _mm256_load_pd(&ev[12]);
	      evv[4] = _mm256_load_pd(&ev[16]);		

#ifdef _FMA
	      for(k = 0; k < 5; k++)
		vv[k] = FMAMACC(vv[k],x1v,evv[k]);		 
#else	      
	      vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
	      vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
	      vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
	      vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
	      vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
	    }	  

	   	     
	  __m256d minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD );
	  
	  scale = 1;
	  
	  for(l = 0; scale && (l < 20); l += 4)
	    {	       
	      __m256d 
		v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) != 15)
		scale = 0;
	    }	    	  	  
	 

	  if(scale)
	    {
	      __m256d 
		twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
	      
	      for(l = 0; l < 20; l += 4)
		vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
	  
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	      
	    }

	  _mm256_store_pd(&v[0], vv[0]);
	  _mm256_store_pd(&v[4], vv[1]);
	  _mm256_store_pd(&v[8], vv[2]);
	  _mm256_store_pd(&v[12], vv[3]);
	  _mm256_store_pd(&v[16], vv[4]);	       
	}
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 400];
	  ri = &right[cptr[i] * 400];

	  vl = &x1[20 * i];
	  vr = &x2[20 * i];
	  v = &x3[20 * i];

	  __m256d vv[5];
	  
	  vv[0] = _mm256_setzero_pd();
	  vv[1] = _mm256_setzero_pd();
	  vv[2] = _mm256_setzero_pd();
	  vv[3] = _mm256_setzero_pd();
	  vv[4] = _mm256_setzero_pd();
	  
	  for(l = 0; l < 20; l++)
	    {	       
	      __m256d 
		x1v = _mm256_setzero_pd(),
		x2v = _mm256_setzero_pd();	
	      
	      double 
		*ev = &extEV[l * 20],
		*lv = &le[l * 20],
		*rv = &ri[l * 20];														
	      
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
	      
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));

	      x1v = hadd4(x1v, x2v);			
#ifdef _FMA
	       for(k = 0; k < 5; k++) 
		 {
		   __m256d evv = _mm256_load_pd(&ev[k*4]);
		   vv[k] = FMAMACC(vv[k],x1v,evv);
		 }
#else	      
	      __m256d 
		evv[5];
	      
	      evv[0] = _mm256_load_pd(&ev[0]);
	      evv[1] = _mm256_load_pd(&ev[4]);
	      evv[2] = _mm256_load_pd(&ev[8]);
	      evv[3] = _mm256_load_pd(&ev[12]);
	      evv[4] = _mm256_load_pd(&ev[16]);		
	      
	      vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
	      vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
	      vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
	      vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
	      vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
	    }	  

	   	     
	  __m256d minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD );
	  
	  scale = 1;
	  
	  for(l = 0; scale && (l < 20); l += 4)
	    {	       
	      __m256d 
		v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) != 15)
		scale = 0;
	    }	    	  	  

	  if(scale)
	    {
	      __m256d 
		twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
	      
	      for(l = 0; l < 20; l += 4)
		vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
	  
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	      
	    }

	  _mm256_store_pd(&v[0], vv[0]);
	  _mm256_store_pd(&v[4], vv[1]);
	  _mm256_store_pd(&v[8], vv[2]);
	  _mm256_store_pd(&v[12], vv[3]);
	  _mm256_store_pd(&v[16], vv[4]);
	 
	}
      break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;
}

void newviewGTRCATPROT_AVX_GAPPED_SAVE(int tipCase, double *extEV,
				       int *cptr,
				       double *x1, double *x2, double *x3, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
				       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le, 
    *ri, 
    *v, 
    *vl, 
    *vr,
    *x1_ptr = x1,
    *x2_ptr = x2, 
    *x3_ptr = x3;
  
  int 
    i, 
    l, 
    scale, 
    addScale = 0,
    scaleGap = 0;

#ifdef _FMA
  int k;
#endif

  {
    le = &left[maxCats * 400];
    ri = &right[maxCats * 400];
    
    vl = x1_gapColumn;
    vr = x2_gapColumn;
    v  = x3_gapColumn;

    __m256d vv[5];
    
    vv[0] = _mm256_setzero_pd();
    vv[1] = _mm256_setzero_pd();
    vv[2] = _mm256_setzero_pd();
    vv[3] = _mm256_setzero_pd();
    vv[4] = _mm256_setzero_pd();
    
    for(l = 0; l < 20; l++)
      {	       
	__m256d 
	  x1v = _mm256_setzero_pd(),
	  x2v = _mm256_setzero_pd();	
	
	double 
	  *ev = &extEV[l * 20],
	  *lv = &le[l * 20],
	  *rv = &ri[l * 20];														
	
	x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
	x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
	x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
	x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
	x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
	
	x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
	x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
	x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
	x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
	x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));
	
	x1v = hadd4(x1v, x2v);			
#ifdef _FMA
	for(k = 0; k < 5; k++) 
	  {
	    __m256d evv = _mm256_load_pd(&ev[k*4]);
	    vv[k] = FMAMACC(vv[k],x1v,evv);
	  }
#else	      
	__m256d 
	  evv[5];
	
	evv[0] = _mm256_load_pd(&ev[0]);
	evv[1] = _mm256_load_pd(&ev[4]);
	evv[2] = _mm256_load_pd(&ev[8]);
	evv[3] = _mm256_load_pd(&ev[12]);
	evv[4] = _mm256_load_pd(&ev[16]);		
	
	vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
	vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
	vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
	vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
	vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
      }	  


     if(tipCase != PLL_TIP_TIP)
       {
	 __m256d minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD );
	  
	 scale = 1;
	  
	 for(l = 0; scale && (l < 20); l += 4)
	   {	       
	     __m256d 
	       v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
	     v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	     
	     if(_mm256_movemask_pd( v1 ) != 15)
	       scale = 0;
	   }	    	  	  

	 if(scale)
	   {
	      __m256d 
		twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
	      
	      for(l = 0; l < 20; l += 4)
		vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 	      	     	      
	   
	      scaleGap = 1;
	   }
       }

     _mm256_store_pd(&v[0], vv[0]);
     _mm256_store_pd(&v[4], vv[1]);
     _mm256_store_pd(&v[8], vv[2]);
     _mm256_store_pd(&v[12], vv[3]);
     _mm256_store_pd(&v[16], vv[4]);     
  }



  switch(tipCase)
    {
    case PLL_TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    if(noGap(x3_gap, i))	   
	      {	    
		vl = &(tipVector[20 * tipX1[i]]);
		vr = &(tipVector[20 * tipX2[i]]);
		v  = x3_ptr;	    	    	   	    

		if(isGap(x1_gap, i))
		  le =  &left[maxCats * 400];
		else	  	  
		  le =  &left[cptr[i] * 400];	  
		
		if(isGap(x2_gap, i))
		  ri =  &right[maxCats * 400];
		else	 	  
		  ri =  &right[cptr[i] * 400];

		__m256d vv[5];
		
		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();	   	    
		
		for(l = 0; l < 20; l++)
		  {	       
		    __m256d 
		      x1v = _mm256_setzero_pd(),
		      x2v = _mm256_setzero_pd();	
		    
		    double 
		      *ev = &extEV[l * 20],
		      *lv = &le[l * 20],
		      *rv = &ri[l * 20];														
		    
#ifdef _FMA		
		    for(k = 0; k < 20; k += 4) 
		      {
			__m256d vlv = _mm256_load_pd(&vl[k]);
			__m256d lvv = _mm256_load_pd(&lv[k]);
			x1v = FMAMACC(x1v,vlv,lvv);
			__m256d vrv = _mm256_load_pd(&vr[k]);
			__m256d rvv = _mm256_load_pd(&rv[k]);
			x2v = FMAMACC(x2v,vrv,rvv);
		      }
#else		
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
		    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));	
#endif
		    
		    x1v = hadd4(x1v, x2v);			
#ifdef _FMA
		    for(k = 0; k < 5; k++) 
		      {
			__m256d evv = _mm256_load_pd(&ev[k*4]);
			vv[k] = FMAMACC(vv[k],x1v,evv);
		      }	  
#else		
		    __m256d 
		      evv[5];
		    
		    evv[0] = _mm256_load_pd(&ev[0]);
		    evv[1] = _mm256_load_pd(&ev[4]);
		    evv[2] = _mm256_load_pd(&ev[8]);
		    evv[3] = _mm256_load_pd(&ev[12]);
		    evv[4] = _mm256_load_pd(&ev[16]);		
		    
		    vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
		    vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
		    vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
		    vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
		    vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      		      	  
#endif
		  }
		
		_mm256_store_pd(&v[0], vv[0]);
		_mm256_store_pd(&v[4], vv[1]);
		_mm256_store_pd(&v[8], vv[2]);
		_mm256_store_pd(&v[12], vv[3]);
		_mm256_store_pd(&v[16], vv[4]);

		x3_ptr += 20;
	      }
	  }
      }
      break;
    case PLL_TIP_INNER:      	
      for (i = 0; i < n; i++)
	{
	  if(isGap(x3_gap, i))
	    {
	      if(scaleGap)
		{
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i] += 1;		   		    
		}	     
	    }
	  else
	    {
	      vl = &(tipVector[20 * tipX1[i]]);

	      vr = x2_ptr;
	      v = x3_ptr;
	      
	      if(isGap(x1_gap, i))
		le =  &left[maxCats * 400];
	      else
		le =  &left[cptr[i] * 400];
	      
	      if(isGap(x2_gap, i))
		{		 
		  ri =  &right[maxCats * 400];
		  vr = x2_gapColumn;
		}
	      else
		{
		  ri =  &right[cptr[i] * 400];
		  vr = x2_ptr;
		  x2_ptr += 20;
		}	  	  
	  
	      __m256d vv[5];
	      
	      vv[0] = _mm256_setzero_pd();
	      vv[1] = _mm256_setzero_pd();
	      vv[2] = _mm256_setzero_pd();
	      vv[3] = _mm256_setzero_pd();
	      vv[4] = _mm256_setzero_pd();
	      	      	      
	      for(l = 0; l < 20; l++)
		{	       
		  __m256d 
		    x1v = _mm256_setzero_pd(),
		    x2v = _mm256_setzero_pd();	
		  
		  double 
		    *ev = &extEV[l * 20],
		    *lv = &le[l * 20],
		    *rv = &ri[l * 20];														
#ifdef _FMA
		  for(k = 0; k < 20; k += 4) 
		    {
		      __m256d vlv = _mm256_load_pd(&vl[k]);
		      __m256d lvv = _mm256_load_pd(&lv[k]);
		      x1v = FMAMACC(x1v,vlv,lvv);
		      __m256d vrv = _mm256_load_pd(&vr[k]);
		      __m256d rvv = _mm256_load_pd(&rv[k]);
		      x2v = FMAMACC(x2v,vrv,rvv);
		    }
#else	      
		  x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
		  x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
		  x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
		  x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
		  x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
		  
		  x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
		  x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
		  x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
		  x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
		  x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));
#endif
		  
		  x1v = hadd4(x1v, x2v);			
		  
		  __m256d 
		    evv[5];
		  
		  evv[0] = _mm256_load_pd(&ev[0]);
		  evv[1] = _mm256_load_pd(&ev[4]);
		  evv[2] = _mm256_load_pd(&ev[8]);
		  evv[3] = _mm256_load_pd(&ev[12]);
		  evv[4] = _mm256_load_pd(&ev[16]);		
		  
#ifdef _FMA
		  for(k = 0; k < 5; k++)
		    vv[k] = FMAMACC(vv[k],x1v,evv[k]);		 
#else	      
		  vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
		  vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
		  vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
		  vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
		  vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
		}	  

	   	     
	      __m256d minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD );
	  
	      scale = 1;
	      
	      for(l = 0; scale && (l < 20); l += 4)
		{	       
		  __m256d 
		    v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
		  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		  
		  if(_mm256_movemask_pd( v1 ) != 15)
		    scale = 0;
		}	    	  	  
	 
	      if(scale)
		{
		  __m256d 
		    twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
		  
		  for(l = 0; l < 20; l += 4)
		    vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
		  
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i]  += 1;	      
		}

	      _mm256_store_pd(&v[0], vv[0]);
	      _mm256_store_pd(&v[4], vv[1]);
	      _mm256_store_pd(&v[8], vv[2]);
	      _mm256_store_pd(&v[12], vv[3]);
	      _mm256_store_pd(&v[16], vv[4]);	       
	      
	      x3_ptr += 20;
	    }
	}    
      break;
    case PLL_INNER_INNER:
      for(i = 0; i < n; i++)
	{
	   if(isGap(x3_gap, i))
	     {
	       if(scaleGap)		   		    
		 {
		   if(useFastScaling)
		     addScale += wgt[i];
		   else
		     ex3[i] += 1;
		 }		 	       
	     }
	   else
	     {

	        v = x3_ptr;

		if(isGap(x1_gap, i))
		  {
		    vl = x1_gapColumn;
		    le =  &left[maxCats * 400];
		  }
		else
		  {
		    le =  &left[cptr[i] * 400];
		    vl = x1_ptr;
		    x1_ptr += 20;
		  }
		
		if(isGap(x2_gap, i))	
		  {
		    vr = x2_gapColumn;
		    ri =  &right[maxCats * 400];	    
		  }
		else
		  {
		    ri =  &right[cptr[i] * 400];
		    vr = x2_ptr;
		    x2_ptr += 20;
		  }	 	  	 
		
		__m256d vv[5];
		
		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l++)
		  {	       
		    __m256d 
		      x1v = _mm256_setzero_pd(),
		      x2v = _mm256_setzero_pd();	
		    
		    double 
		      *ev = &extEV[l * 20],
		      *lv = &le[l * 20],
		      *rv = &ri[l * 20];														
		    
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
		    x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
		    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
		    x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));
		    
		    x1v = hadd4(x1v, x2v);			
#ifdef _FMA
		    for(k = 0; k < 5; k++) 
		      {
			__m256d evv = _mm256_load_pd(&ev[k*4]);
			vv[k] = FMAMACC(vv[k],x1v,evv);
		      }
#else	      
		    __m256d 
		      evv[5];
		    
		    evv[0] = _mm256_load_pd(&ev[0]);
		    evv[1] = _mm256_load_pd(&ev[4]);
		    evv[2] = _mm256_load_pd(&ev[8]);
		    evv[3] = _mm256_load_pd(&ev[12]);
		    evv[4] = _mm256_load_pd(&ev[16]);		
		    
		    vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
		    vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
		    vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
		    vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
		    vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
		  }	  

	   	     
		__m256d minlikelihood_avx = _mm256_set1_pd( PLL_MINLIKELIHOOD );
		
		scale = 1;
		
		for(l = 0; scale && (l < 20); l += 4)
		  {	       
		    __m256d 
		      v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }	    	  	  
		
		if(scale)
		  {
		    __m256d 
		      twoto = _mm256_set1_pd(PLL_TWOTOTHE256);
		    
		    for(l = 0; l < 20; l += 4)
		      vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
		    
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;	      
		  }

		_mm256_store_pd(&v[0], vv[0]);
		_mm256_store_pd(&v[4], vv[1]);
		_mm256_store_pd(&v[8], vv[2]);
		_mm256_store_pd(&v[12], vv[3]);
		_mm256_store_pd(&v[16], vv[4]);

		 x3_ptr += 20;
	     }
	}   
      break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;
}



void newviewGTRGAMMAPROT_AVX_LG4(int tipCase,
				 double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
				 int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
				 double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling) 
{
  double	
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr;
  
  int	
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

 
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif


#if GCC_VERSION < 40500 && defined (__GNUC__)
   __m256d
    bitmask = _mm256_set_pd(0,0,0,-1);
#else
  __m256i
    bitmask = _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1);
#endif 
  
  switch(tipCase) 
    {
    case PLL_TIP_TIP: 
      {
       
    PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  umpX2[1840] PLL_ALIGN_END;

	
	for(i = 0; i < 23; i++) 
	  {	    	    
	    for(k = 0; k < 80; k++) 
	      {
		double 
		  *ll =  &left[k * 20],
		  *rr =  &right[k * 20];
		
		__m256d 
		  umpX1v = _mm256_setzero_pd(),
		  umpX2v = _mm256_setzero_pd();
		
		v = &(tipVector[k / 20][20 * i]);

		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
#ifdef _FMA
		    __m256d llv = _mm256_load_pd(&ll[l]);
		    umpX1v = FMAMACC(umpX1v,vv,llv);
		    __m256d rrv = _mm256_load_pd(&rr[l]);
		    umpX2v = FMAMACC(umpX2v,vv,rrv);
#else		    
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&ll[l])));
		    umpX2v = _mm256_add_pd(umpX2v,_mm256_mul_pd(vv,_mm256_load_pd(&rr[l])));
#endif
		  }
		
		umpX1v = hadd3(umpX1v);
		umpX2v = hadd3(umpX2v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
		_mm256_maskstore_pd(&umpX2[80 * i + k], bitmask, umpX2v);
	      } 
	  }

	for(i = 0; i < n; i++) 
	  {	    
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];
	   
	    for(j = 0; j < 4; j++) 
	      {     	
		__m256d vv[5];  

		v = &x3[i * 80 + j * 20];
			
		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();

		for(k = 0; k < 20; k++) 
		  {			 
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];

		    __m256d x1px2v = _mm256_set1_pd(x1px2);		    
		    
		    __m256d extEvv = _mm256_load_pd(&extEV[j][20 * k]);
#ifdef _FMA
		    vv[0] = FMAMACC(vv[0],x1px2v,extEvv);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[0],vv[0]);
		    
		    extEvv = _mm256_load_pd(&extEV[j][20 * k + 4]);
#ifdef _FMA
		    vv[1] = FMAMACC(vv[1],x1px2v,extEvv);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

		    extEvv = _mm256_load_pd(&extEV[j][20 * k + 8]);
#ifdef _FMA
		    vv[2] = FMAMACC(vv[2],x1px2v,extEvv);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[8],vv[2]);

		    extEvv = _mm256_load_pd(&extEV[j][20 * k + 12]);
#ifdef _FMA
		    vv[3] = FMAMACC(vv[3],x1px2v,extEvv);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[12],vv[3]);

		    extEvv = _mm256_load_pd(&extEV[j][20 * k + 16]);
#ifdef _FMA
		    vv[4] = FMAMACC(vv[4],x1px2v,extEvv);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[16],vv[4]);
		  } 
	      } 
	  } 
      } 
      break;
    case PLL_TIP_INNER: 
      {

    	  PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  ump_x2[20] PLL_ALIGN_END;

	for(i = 0; i < 23; i++) 
	  {	   
	    for(k = 0; k < 80; k++) 
	      {
		__m256d umpX1v = _mm256_setzero_pd();
		
		 v = &(tipVector[k / 20][20 * i]);

		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    __m256d leftv = _mm256_load_pd(&left[k * 20 + l]);
#ifdef _FMA
		   
		    umpX1v = FMAMACC(umpX1v, vv, leftv);
#else
		    umpX1v = _mm256_add_pd(umpX1v, _mm256_mul_pd(vv, leftv));
#endif
		  }
		umpX1v = hadd3(umpX1v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
	      } 
	  }
	
	for (i = 0; i < n; i++) 
	  {	   
	    uX1 = &umpX1[80 * tipX1[i]];
	   	    
	    for(k = 0; k < 4; k++) 
	      {
		v = &(x2[80 * i + k * 20]);
		
		for(l = 0; l < 20; l++) 
		  {
		    __m256d ump_x2v = _mm256_setzero_pd();
		    		  
		    __m256d vv = _mm256_load_pd(&v[0]);
		    __m256d rightv = _mm256_load_pd(&right[k*400+l*20+0]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    vv = _mm256_load_pd(&v[4]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+4]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[8]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+8]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[12]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+12]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[16]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+16]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    ump_x2v = hadd3(ump_x2v);
		    _mm256_maskstore_pd(&ump_x2[l], bitmask, ump_x2v);
		  }
		
		v = &(x3[80 * i + 20 * k]);
	

		__m256d vv[5]; 

		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l++) 
		  {
		    x1px2 = uX1[k * 20 + l]	* ump_x2[l];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);	
	    		 
#ifdef _FMA
		    __m256d ev = _mm256_load_pd(&extEV[k][l * 20 + 0]);
		    vv[0] = FMAMACC(vv[0],x1px2v, ev);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[k][l * 20 + 0])));
#endif
		    _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[k][l * 20 + 4]);
		    vv[1] = FMAMACC(vv[1],x1px2v, ev);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[k][l * 20 + 4])));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[k][l * 20 + 8]);
		    vv[2] = FMAMACC(vv[2],x1px2v, ev);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[k][l * 20 + 8])));
#endif
		    _mm256_store_pd(&v[8],vv[2]);
		    
#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[k][l * 20 + 12]);
		    vv[3] = FMAMACC(vv[3],x1px2v, ev);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[k][l * 20 + 12])));
#endif
		    _mm256_store_pd(&v[12],vv[3]);


#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[k][l * 20 + 16]);
		    vv[4] = FMAMACC(vv[4],x1px2v, ev);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[k][l * 20 + 16])));
#endif
		    _mm256_store_pd(&v[16],vv[4]);

		  } 
	      }
	   
	    v = &x3[80 * i];
	    __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 4) 
	      {
		__m256d vv = _mm256_load_pd(&v[l]);
		__m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		if(_mm256_movemask_pd(vv_abs) != 15)
		  scale = 0;
	      }
	    
	    if(scale) 
	      {		
		__m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
		for(l = 0; l < 80; l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		  }
		if(useFastScaling)
		  addScale += wgt[i];				
		else
		  ex3[i] += 1;
	      } 
	  } 
      } 
      break;
    case PLL_INNER_INNER:      
      for(i = 0; i < n; i++) 
	{ 
	  scale = 1;
	  
	  for(k = 0; k < 4; k++) 
	    {
	      vl = &(x1[80 * i + 20 * k]);
	      vr = &(x2[80 * i + 20 * k]);
	      v  = &(x3[80 * i + 20 * k]);	      	   

	      __m256d vv[5]; 
	      
	      vv[0] = _mm256_setzero_pd();
	      vv[1] = _mm256_setzero_pd();
	      vv[2] = _mm256_setzero_pd();
	      vv[3] = _mm256_setzero_pd();
	      vv[4] = _mm256_setzero_pd();
	      
	      for(l = 0; l < 20; l++) 
		{		  
		  __m256d al = _mm256_setzero_pd();
		  __m256d ar = _mm256_setzero_pd();
       		  
		  __m256d leftv  = _mm256_load_pd(&left[k * 400 + l * 20 + 0]);
		  __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 0]);
		  __m256d vlv = _mm256_load_pd(&vl[0]);
		  __m256d vrv = _mm256_load_pd(&vr[0]);
		  
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));		  
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 4]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 4]);
		  vlv = _mm256_load_pd(&vl[4]);
		  vrv = _mm256_load_pd(&vr[4]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 8]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 8]);
		  vlv = _mm256_load_pd(&vl[8]);
		  vrv = _mm256_load_pd(&vr[8]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 12]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 12]);
		  vlv = _mm256_load_pd(&vl[12]);
		  vrv = _mm256_load_pd(&vr[12]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 16]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 16]);
		  vlv = _mm256_load_pd(&vl[16]);
		  vrv = _mm256_load_pd(&vr[16]);

#ifdef _FMA		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  /**************************************************************************************************************/

		  al = hadd3(al);
		  ar = hadd3(ar);
		  al = _mm256_mul_pd(ar,al);
		  
		  /************************************************************************************************************/
#ifdef _FMA		    
		  __m256d ev =  _mm256_load_pd(&extEV[k][20 * l + 0]);
		  vv[0] = FMAMACC(vv[0], al, ev);		 
#else
		  vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(al, _mm256_load_pd(&extEV[k][20 * l + 0])));			  		 		  
#endif
		  _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[k][20 * l + 4]);
		  vv[1] = FMAMACC(vv[1], al, ev);		 
#else
		  vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(al, _mm256_load_pd(&extEV[k][20 * l + 4])));		  		 
#endif
		  _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[k][20 * l + 8]);
		  vv[2] = FMAMACC(vv[2], al, ev);		 
#else
		  vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(al, _mm256_load_pd(&extEV[k][20 * l + 8])));		  		 
#endif
		  _mm256_store_pd(&v[8],vv[2]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[k][20 * l + 12]);
		  vv[3] = FMAMACC(vv[3], al, ev);		 
#else
		  vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(al, _mm256_load_pd(&extEV[k][20 * l + 12])));		  		 
#endif
		  _mm256_store_pd(&v[12],vv[3]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[k][20 * l + 16]);
		  vv[4] = FMAMACC(vv[4], al, ev);		 
#else
		  vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(al, _mm256_load_pd(&extEV[k][20 * l + 16])));			 	  
#endif
		  _mm256_store_pd(&v[16],vv[4]);		 
		} 
	    }
	  v = &(x3[80 * i]);
	  scale = 1;
	  __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);	 

	  for(l = 0; scale && (l < 80); l += 4) 
	    {
	      __m256d vv = _mm256_load_pd(&v[l]);
	      __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
	      vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
	      if(_mm256_movemask_pd(vv_abs) != 15)
		scale = 0;	     
	    }

	  if(scale) 
	    {		     	      
	      __m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
	      for(l = 0; l < 80; l += 4) 
		{
		  __m256d vv = _mm256_load_pd(&v[l]);
		  _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		}
	      if(useFastScaling)
		addScale += wgt[i];					
	      else
		ex3[i] += 1;
	    } 
	}
      break;
    default:
      assert(0);
    }
 
  if(useFastScaling)
    *scalerIncrement = addScale;
}
 

void newviewGTRGAMMAPROT_AVX(int tipCase,
			     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			     int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
			     double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling) 
{
  double	
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr;
  
  int	
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

 
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif


#if GCC_VERSION < 40500 && defined(__GNUC__)
   __m256d
    bitmask = _mm256_set_pd(0,0,0,-1);
#else
  __m256i
    bitmask = _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1);
#endif 
  
  switch(tipCase) 
    {
    case PLL_TIP_TIP: 
      {
       
    PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  umpX2[1840] PLL_ALIGN_END;

	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);
	    
	    for(k = 0; k < 80; k++) 
	      {
		double 
		  *ll =  &left[k * 20],
		  *rr =  &right[k * 20];
		
		__m256d 
		  umpX1v = _mm256_setzero_pd(),
		  umpX2v = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
#ifdef _FMA
		    __m256d llv = _mm256_load_pd(&ll[l]);
		    umpX1v = FMAMACC(umpX1v,vv,llv);
		    __m256d rrv = _mm256_load_pd(&rr[l]);
		    umpX2v = FMAMACC(umpX2v,vv,rrv);
#else		    
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&ll[l])));
		    umpX2v = _mm256_add_pd(umpX2v,_mm256_mul_pd(vv,_mm256_load_pd(&rr[l])));
#endif
		  }
		
		umpX1v = hadd3(umpX1v);
		umpX2v = hadd3(umpX2v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
		_mm256_maskstore_pd(&umpX2[80 * i + k], bitmask, umpX2v);
	      } 
	  }

	for(i = 0; i < n; i++) 
	  {	    
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];
	   
	    for(j = 0; j < 4; j++) 
	      {     	
		__m256d vv[5];  

		v = &x3[i * 80 + j * 20];
			
		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();

		for(k = 0; k < 20; k++) 
		  {			 
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];

		    __m256d x1px2v = _mm256_set1_pd(x1px2);		    
		    
		    __m256d extEvv = _mm256_load_pd(&extEV[20 * k]);
#ifdef _FMA
		    vv[0] = FMAMACC(vv[0],x1px2v,extEvv);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[0],vv[0]);
		    
		    extEvv = _mm256_load_pd(&extEV[20 * k + 4]);
#ifdef _FMA
		    vv[1] = FMAMACC(vv[1],x1px2v,extEvv);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 8]);
#ifdef _FMA
		    vv[2] = FMAMACC(vv[2],x1px2v,extEvv);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[8],vv[2]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 12]);
#ifdef _FMA
		    vv[3] = FMAMACC(vv[3],x1px2v,extEvv);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[12],vv[3]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 16]);
#ifdef _FMA
		    vv[4] = FMAMACC(vv[4],x1px2v,extEvv);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[16],vv[4]);
		  } 
	      } 
	  } 
      } 
      break;
    case PLL_TIP_INNER: 
      {

    	  PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  ump_x2[20] PLL_ALIGN_END;

	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++) 
	      {
		__m256d umpX1v = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    __m256d leftv = _mm256_load_pd(&left[k * 20 + l]);
#ifdef _FMA
		   
		    umpX1v = FMAMACC(umpX1v, vv, leftv);
#else
		    umpX1v = _mm256_add_pd(umpX1v, _mm256_mul_pd(vv, leftv));
#endif
		  }
		umpX1v = hadd3(umpX1v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
	      } 
	  }
	
	for (i = 0; i < n; i++) 
	  {	   
	    uX1 = &umpX1[80 * tipX1[i]];
	   	    
	    for(k = 0; k < 4; k++) 
	      {
		v = &(x2[80 * i + k * 20]);
		
		for(l = 0; l < 20; l++) 
		  {
		    __m256d ump_x2v = _mm256_setzero_pd();
		    		  
		    __m256d vv = _mm256_load_pd(&v[0]);
		    __m256d rightv = _mm256_load_pd(&right[k*400+l*20+0]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    vv = _mm256_load_pd(&v[4]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+4]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[8]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+8]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[12]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+12]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[16]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+16]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    ump_x2v = hadd3(ump_x2v);
		    _mm256_maskstore_pd(&ump_x2[l], bitmask, ump_x2v);
		  }
		
		v = &(x3[80 * i + 20 * k]);
	

		__m256d vv[5]; 

		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l++) 
		  {
		    x1px2 = uX1[k * 20 + l]	* ump_x2[l];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);	
	    		 
#ifdef _FMA
		    __m256d ev = _mm256_load_pd(&extEV[l * 20 + 0]);
		    vv[0] = FMAMACC(vv[0],x1px2v, ev);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 0])));
#endif
		    _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 4]);
		    vv[1] = FMAMACC(vv[1],x1px2v, ev);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 4])));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 8]);
		    vv[2] = FMAMACC(vv[2],x1px2v, ev);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 8])));
#endif
		    _mm256_store_pd(&v[8],vv[2]);
		    
#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 12]);
		    vv[3] = FMAMACC(vv[3],x1px2v, ev);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 12])));
#endif
		    _mm256_store_pd(&v[12],vv[3]);


#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 16]);
		    vv[4] = FMAMACC(vv[4],x1px2v, ev);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 16])));
#endif
		    _mm256_store_pd(&v[16],vv[4]);

		  } 
	      }
	   
	    v = &x3[80 * i];
	    __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 4) 
	      {
		__m256d vv = _mm256_load_pd(&v[l]);
		__m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		if(_mm256_movemask_pd(vv_abs) != 15)
		  scale = 0;
	      }
	    
	    if(scale) 
	      {		
		__m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
		for(l = 0; l < 80; l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		  }
		if(useFastScaling)
		  addScale += wgt[i];				
		else
		  ex3[i] += 1;
	      } 
	  } 
      } 
      break;
    case PLL_INNER_INNER:      
      for(i = 0; i < n; i++) 
	{ 
	  scale = 1;
	  
	  for(k = 0; k < 4; k++) 
	    {
	      vl = &(x1[80 * i + 20 * k]);
	      vr = &(x2[80 * i + 20 * k]);
	      v  = &(x3[80 * i + 20 * k]);	      	   

	      __m256d vv[5]; 
	      
	      vv[0] = _mm256_setzero_pd();
	      vv[1] = _mm256_setzero_pd();
	      vv[2] = _mm256_setzero_pd();
	      vv[3] = _mm256_setzero_pd();
	      vv[4] = _mm256_setzero_pd();
	      
	      for(l = 0; l < 20; l++) 
		{		  
		  __m256d al = _mm256_setzero_pd();
		  __m256d ar = _mm256_setzero_pd();
       		  
		  __m256d leftv  = _mm256_load_pd(&left[k * 400 + l * 20 + 0]);
		  __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 0]);
		  __m256d vlv = _mm256_load_pd(&vl[0]);
		  __m256d vrv = _mm256_load_pd(&vr[0]);
		  
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));		  
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 4]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 4]);
		  vlv = _mm256_load_pd(&vl[4]);
		  vrv = _mm256_load_pd(&vr[4]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 8]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 8]);
		  vlv = _mm256_load_pd(&vl[8]);
		  vrv = _mm256_load_pd(&vr[8]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 12]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 12]);
		  vlv = _mm256_load_pd(&vl[12]);
		  vrv = _mm256_load_pd(&vr[12]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 16]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 16]);
		  vlv = _mm256_load_pd(&vl[16]);
		  vrv = _mm256_load_pd(&vr[16]);

#ifdef _FMA		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  /**************************************************************************************************************/

		  al = hadd3(al);
		  ar = hadd3(ar);
		  al = _mm256_mul_pd(ar,al);
		  
		  /************************************************************************************************************/
#ifdef _FMA		    
		  __m256d ev =  _mm256_load_pd(&extEV[20 * l + 0]);
		  vv[0] = FMAMACC(vv[0], al, ev);		 
#else
		  vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 0])));			  		 		  
#endif
		  _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 4]);
		  vv[1] = FMAMACC(vv[1], al, ev);		 
#else
		  vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 4])));		  		 
#endif
		  _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 8]);
		  vv[2] = FMAMACC(vv[2], al, ev);		 
#else
		  vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 8])));		  		 
#endif
		  _mm256_store_pd(&v[8],vv[2]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 12]);
		  vv[3] = FMAMACC(vv[3], al, ev);		 
#else
		  vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 12])));		  		 
#endif
		  _mm256_store_pd(&v[12],vv[3]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 16]);
		  vv[4] = FMAMACC(vv[4], al, ev);		 
#else
		  vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 16])));			 	  
#endif
		  _mm256_store_pd(&v[16],vv[4]);		 
		} 
	    }
	  v = &(x3[80 * i]);
	  scale = 1;
	  __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);	 

	  for(l = 0; scale && (l < 80); l += 4) 
	    {
	      __m256d vv = _mm256_load_pd(&v[l]);
	      __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
	      vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
	      if(_mm256_movemask_pd(vv_abs) != 15)
		scale = 0;	     
	    }

	  if(scale) 
	    {		     	      
	      __m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
	      for(l = 0; l < 80; l += 4) 
		{
		  __m256d vv = _mm256_load_pd(&v[l]);
		  _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		}
	      if(useFastScaling)
		addScale += wgt[i];					
	      else
		ex3[i] += 1;
	    } 
	}
      break;
    default:
      assert(0);
    }
 
  if(useFastScaling)
    *scalerIncrement = addScale;
}



void newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(int tipCase,
					 double *x1_start, double *x2_start, double *x3_start, double *extEV, double *tipVector,
					 int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
					 double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean useFastScaling,
					 unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					 double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn) 
{
  double	
    *x1 = x1_start,
    *x2 = x2_start,
    *x3_ptr = x3_start,
    *x2_ptr = x2_start,
    *x1_ptr = x1_start,
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr;
  
  int	
    i, 
    j, 
    l, 
    k, 
    gapScaling = 0,
    scale, 
    addScale = 0;

 
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif


#if GCC_VERSION < 40500 && defined(__GNUC__)
   __m256d
    bitmask = _mm256_set_pd(0,0,0,-1);
#else
  __m256i
    bitmask = _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1);
#endif 
  
  switch(tipCase) 
    {
    case PLL_TIP_TIP: 
      {       
    	  PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  umpX2[1840] PLL_ALIGN_END;



	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);
	    
	    for(k = 0; k < 80; k++) 
	      {
		double 
		  *ll =  &left[k * 20],
		  *rr =  &right[k * 20];
		
		__m256d 
		  umpX1v = _mm256_setzero_pd(),
		  umpX2v = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
#ifdef _FMA
		    __m256d llv = _mm256_load_pd(&ll[l]);
		    umpX1v = FMAMACC(umpX1v,vv,llv);
		    __m256d rrv = _mm256_load_pd(&rr[l]);
		    umpX2v = FMAMACC(umpX2v,vv,rrv);
#else		    
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&ll[l])));
		    umpX2v = _mm256_add_pd(umpX2v,_mm256_mul_pd(vv,_mm256_load_pd(&rr[l])));
#endif
		  }
		
		umpX1v = hadd3(umpX1v);
		umpX2v = hadd3(umpX2v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
		_mm256_maskstore_pd(&umpX2[80 * i + k], bitmask, umpX2v);
	      } 
	  }

	
	{	    
	  uX1 = &umpX1[1760];
	  uX2 = &umpX2[1760];
	  
	  for(j = 0; j < 4; j++) 
	    {     	
	      __m256d vv[5];  
	      
	      v = &x3_gapColumn[j * 20];
	      
	      vv[0] = _mm256_setzero_pd();
	      vv[1] = _mm256_setzero_pd();
	      vv[2] = _mm256_setzero_pd();
	      vv[3] = _mm256_setzero_pd();
	      vv[4] = _mm256_setzero_pd();
	      
	      for(k = 0; k < 20; k++) 
		{			 
		  x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		  
		  __m256d x1px2v = _mm256_set1_pd(x1px2);		    
		  
		  __m256d extEvv = _mm256_load_pd(&extEV[20 * k]);
#ifdef _FMA
		  vv[0] = FMAMACC(vv[0],x1px2v,extEvv);
#else
		  vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v,extEvv));
#endif
		  _mm256_store_pd(&v[0],vv[0]);
		  
		  extEvv = _mm256_load_pd(&extEV[20 * k + 4]);
#ifdef _FMA
		  vv[1] = FMAMACC(vv[1],x1px2v,extEvv);
#else
		  vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v,extEvv));
#endif
		  _mm256_store_pd(&v[4],vv[1]);
		  
		  extEvv = _mm256_load_pd(&extEV[20 * k + 8]);
#ifdef _FMA
		  vv[2] = FMAMACC(vv[2],x1px2v,extEvv);
#else
		  vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v,extEvv));
#endif
		  _mm256_store_pd(&v[8],vv[2]);
		  
		  extEvv = _mm256_load_pd(&extEV[20 * k + 12]);
#ifdef _FMA
		  vv[3] = FMAMACC(vv[3],x1px2v,extEvv);
#else
		  vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v,extEvv));
#endif
		  _mm256_store_pd(&v[12],vv[3]);
		  
		  extEvv = _mm256_load_pd(&extEV[20 * k + 16]);
#ifdef _FMA
		  vv[4] = FMAMACC(vv[4],x1px2v,extEvv);
#else
		  vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v,extEvv));
#endif
		  _mm256_store_pd(&v[16],vv[4]);
		} 
	    } 
	}

	
	for(i = 0; i < n; i++) 
	  {
	    if(!(x3_gap[i / 32] & mask32[i % 32]))
	      {	    
		uX1 = &umpX1[80 * tipX1[i]];
		uX2 = &umpX2[80 * tipX2[i]];
	   
		for(j = 0; j < 4; j++) 
		  {     	
		    __m256d vv[5];  
		    
		    v = &x3_ptr[j * 20];
			
		    vv[0] = _mm256_setzero_pd();
		    vv[1] = _mm256_setzero_pd();
		    vv[2] = _mm256_setzero_pd();
		    vv[3] = _mm256_setzero_pd();
		    vv[4] = _mm256_setzero_pd();

		    for(k = 0; k < 20; k++) 
		      {			 
			x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
			
			__m256d x1px2v = _mm256_set1_pd(x1px2);		    
			
			__m256d extEvv = _mm256_load_pd(&extEV[20 * k]);
#ifdef _FMA
			vv[0] = FMAMACC(vv[0],x1px2v,extEvv);
#else
			vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[0],vv[0]);
			
			extEvv = _mm256_load_pd(&extEV[20 * k + 4]);
#ifdef _FMA
			vv[1] = FMAMACC(vv[1],x1px2v,extEvv);
#else
			vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[4],vv[1]);
			
			extEvv = _mm256_load_pd(&extEV[20 * k + 8]);
#ifdef _FMA
			vv[2] = FMAMACC(vv[2],x1px2v,extEvv);
#else
			vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[8],vv[2]);
			
			extEvv = _mm256_load_pd(&extEV[20 * k + 12]);
#ifdef _FMA
			vv[3] = FMAMACC(vv[3],x1px2v,extEvv);
#else
			vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[12],vv[3]);
			
			extEvv = _mm256_load_pd(&extEV[20 * k + 16]);
#ifdef _FMA
			vv[4] = FMAMACC(vv[4],x1px2v,extEvv);
#else
			vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[16],vv[4]);
		      } 
		  }
		x3_ptr += 80;		  
	      }
	  }
      }
      break;
    case PLL_TIP_INNER: 
      {
    	  PLL_ALIGN_BEGIN double
	  umpX1[1840] PLL_ALIGN_END,
	  ump_x2[20] PLL_ALIGN_END;



	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++) 
	      {
		__m256d umpX1v = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    __m256d leftv = _mm256_load_pd(&left[k * 20 + l]);
#ifdef _FMA
		   
		    umpX1v = FMAMACC(umpX1v, vv, leftv);
#else
		    umpX1v = _mm256_add_pd(umpX1v, _mm256_mul_pd(vv, leftv));
#endif
		  }
		umpX1v = hadd3(umpX1v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
	      } 
	  }

	{	   
	  uX1 = &umpX1[1760];
	   	    
	  for(k = 0; k < 4; k++) 
	    {
	      v = &(x2_gapColumn[k * 20]);
		
		for(l = 0; l < 20; l++) 
		  {
		    __m256d ump_x2v = _mm256_setzero_pd();
		    		  
		    __m256d vv = _mm256_load_pd(&v[0]);
		    __m256d rightv = _mm256_load_pd(&right[k*400+l*20+0]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    vv = _mm256_load_pd(&v[4]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+4]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[8]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+8]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[12]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+12]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[16]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+16]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    ump_x2v = hadd3(ump_x2v);
		    _mm256_maskstore_pd(&ump_x2[l], bitmask, ump_x2v);
		  }
		
		v = &x3_gapColumn[20 * k];
	
		__m256d vv[5]; 

		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l++) 
		  {
		    x1px2 = uX1[k * 20 + l]	* ump_x2[l];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);	
	    		 
#ifdef _FMA
		    __m256d ev = _mm256_load_pd(&extEV[l * 20 + 0]);
		    vv[0] = FMAMACC(vv[0],x1px2v, ev);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 0])));
#endif
		    _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 4]);
		    vv[1] = FMAMACC(vv[1],x1px2v, ev);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 4])));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 8]);
		    vv[2] = FMAMACC(vv[2],x1px2v, ev);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 8])));
#endif
		    _mm256_store_pd(&v[8],vv[2]);
		    
#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 12]);
		    vv[3] = FMAMACC(vv[3],x1px2v, ev);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 12])));
#endif
		    _mm256_store_pd(&v[12],vv[3]);


#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 16]);
		    vv[4] = FMAMACC(vv[4],x1px2v, ev);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 16])));
#endif
		    _mm256_store_pd(&v[16],vv[4]);

		  } 
	      }
	   
	    v = x3_gapColumn;
	    __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 4) 
	      {
		__m256d vv = _mm256_load_pd(&v[l]);
		__m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		if(_mm256_movemask_pd(vv_abs) != 15)
		  scale = 0;
	      }
	    
	    if(scale) 
	      {		
		__m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
		gapScaling = 1;

		for(l = 0; l < 80; l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		  }	
	      } 
	}       
	
	for (i = 0; i < n; i++) 
	  {	   
	    if((x3_gap[i / 32] & mask32[i % 32]))
	      {	       
		if(gapScaling)
		  {
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;
		  }
	      }
	    else
	      {		
		uX1 = &umpX1[80 * tipX1[i]];
		
		if(x2_gap[i / 32] & mask32[i % 32])
		  x2 = x2_gapColumn;
		else
		  {
		    x2 = x2_ptr;
		    x2_ptr += 80;
		  }	      
	    
		for(k = 0; k < 4; k++) 
		  {
		    v = &(x2[k * 20]);
		    
		    for(l = 0; l < 20; l++) 
		      {
			__m256d ump_x2v = _mm256_setzero_pd();
		    	
			__m256d vv = _mm256_load_pd(&v[0]);
			__m256d rightv = _mm256_load_pd(&right[k*400+l*20+0]);
#ifdef _FMA
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
			
			vv = _mm256_load_pd(&v[4]);
			rightv = _mm256_load_pd(&right[k*400+l*20+4]);
#ifdef _FMA
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
			
			vv = _mm256_load_pd(&v[8]);
			rightv = _mm256_load_pd(&right[k*400+l*20+8]);
#ifdef _FMA
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
			
			vv = _mm256_load_pd(&v[12]);
			rightv = _mm256_load_pd(&right[k*400+l*20+12]);
#ifdef _FMA
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
			
			vv = _mm256_load_pd(&v[16]);
			rightv = _mm256_load_pd(&right[k*400+l*20+16]);
#ifdef _FMA
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
			
			ump_x2v = hadd3(ump_x2v);
			_mm256_maskstore_pd(&ump_x2[l], bitmask, ump_x2v);
		      }
		  
		    
		    v = &x3_ptr[k * 20];
		    
		    __m256d vv[5]; 
		    
		    vv[0] = _mm256_setzero_pd();
		    vv[1] = _mm256_setzero_pd();
		    vv[2] = _mm256_setzero_pd();
		    vv[3] = _mm256_setzero_pd();
		    vv[4] = _mm256_setzero_pd();
		    
		    for(l = 0; l < 20; l++) 
		      {
			x1px2 = uX1[k * 20 + l]	* ump_x2[l];
			__m256d x1px2v = _mm256_set1_pd(x1px2);	
			
#ifdef _FMA
			__m256d ev = _mm256_load_pd(&extEV[l * 20 + 0]);
			vv[0] = FMAMACC(vv[0],x1px2v, ev);
#else
			vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 0])));
#endif
			_mm256_store_pd(&v[0],vv[0]);
			
#ifdef _FMA
			ev = _mm256_load_pd(&extEV[l * 20 + 4]);
			vv[1] = FMAMACC(vv[1],x1px2v, ev);
#else
			vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 4])));
#endif
			_mm256_store_pd(&v[4],vv[1]);
			
#ifdef _FMA
			ev = _mm256_load_pd(&extEV[l * 20 + 8]);
			vv[2] = FMAMACC(vv[2],x1px2v, ev);
#else
			vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 8])));
#endif
			_mm256_store_pd(&v[8],vv[2]);
			
#ifdef _FMA
			ev = _mm256_load_pd(&extEV[l * 20 + 12]);
			vv[3] = FMAMACC(vv[3],x1px2v, ev);
#else
			vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 12])));
#endif
			_mm256_store_pd(&v[12],vv[3]);
			
			
#ifdef _FMA
			ev = _mm256_load_pd(&extEV[l * 20 + 16]);
			vv[4] = FMAMACC(vv[4],x1px2v, ev);
#else
			vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 16])));
#endif
			_mm256_store_pd(&v[16],vv[4]);
			
		      } 
		  }
		
		v = x3_ptr;
		__m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);
		scale = 1;
		for(l = 0; scale && (l < 80); l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		    vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		    if(_mm256_movemask_pd(vv_abs) != 15)
		      scale = 0;
		  }
	    
		if(scale) 
		  {		
		    __m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
		    for(l = 0; l < 80; l += 4) 
		      {
			__m256d vv = _mm256_load_pd(&v[l]);
			_mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		      }
		    if(useFastScaling)
		      addScale += wgt[i];				
		    else
		      ex3[i] += 1;
		  }	      
		x3_ptr += 80;
	      }
	  }
      }
      break;
    case PLL_INNER_INNER:    	  
      for(k = 0; k < 4; k++) 
	{
	  vl = &(x1_gapColumn[20 * k]);
	  vr = &(x2_gapColumn[20 * k]);
	  v  = &(x3_gapColumn[20 * k]);	      	   

	  __m256d vv[5]; 
	  
	  vv[0] = _mm256_setzero_pd();
	  vv[1] = _mm256_setzero_pd();
	  vv[2] = _mm256_setzero_pd();
	  vv[3] = _mm256_setzero_pd();
	  vv[4] = _mm256_setzero_pd();
	  
	  for(l = 0; l < 20; l++) 
	    {		  
	      __m256d al = _mm256_setzero_pd();
	      __m256d ar = _mm256_setzero_pd();
	      
	      __m256d leftv  = _mm256_load_pd(&left[k * 400 + l * 20 + 0]);
	      __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 0]);
	      __m256d vlv = _mm256_load_pd(&vl[0]);
	      __m256d vrv = _mm256_load_pd(&vr[0]);
	      
#ifdef _FMA
	      
	      al = FMAMACC(al, vlv, leftv);
	      ar = FMAMACC(ar, vrv, rightv);
#else
	      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
	      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));		  
#endif
	      
	      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 4]);
	      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 4]);
	      vlv = _mm256_load_pd(&vl[4]);
	      vrv = _mm256_load_pd(&vr[4]);
#ifdef _FMA
	      
	      al = FMAMACC(al, vlv, leftv);
	      ar = FMAMACC(ar, vrv, rightv);
#else
	      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
	      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
	      
	      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 8]);
	      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 8]);
	      vlv = _mm256_load_pd(&vl[8]);
	      vrv = _mm256_load_pd(&vr[8]);
#ifdef _FMA
	      
	      al = FMAMACC(al, vlv, leftv);
	      ar = FMAMACC(ar, vrv, rightv);
#else
	      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
	      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
	      
	      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 12]);
	      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 12]);
	      vlv = _mm256_load_pd(&vl[12]);
	      vrv = _mm256_load_pd(&vr[12]);
#ifdef _FMA
	      
	      al = FMAMACC(al, vlv, leftv);
	      ar = FMAMACC(ar, vrv, rightv);
#else
	      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
	      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
	      
	      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 16]);
	      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 16]);
	      vlv = _mm256_load_pd(&vl[16]);
	      vrv = _mm256_load_pd(&vr[16]);
	      
#ifdef _FMA		    
	      al = FMAMACC(al, vlv, leftv);
	      ar = FMAMACC(ar, vrv, rightv);
#else
	      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
	      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
	      
	      /**************************************************************************************************************/
	      
	      al = hadd3(al);
	      ar = hadd3(ar);
	      al = _mm256_mul_pd(ar,al);
	      
	      /************************************************************************************************************/
#ifdef _FMA		    
	      __m256d ev =  _mm256_load_pd(&extEV[20 * l + 0]);
	      vv[0] = FMAMACC(vv[0], al, ev);		 
#else
	      vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 0])));			  		 		  
#endif
	      _mm256_store_pd(&v[0],vv[0]);
	      
#ifdef _FMA		    
	      ev =  _mm256_load_pd(&extEV[20 * l + 4]);
	      vv[1] = FMAMACC(vv[1], al, ev);		 
#else
	      vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 4])));		  		 
#endif
	      _mm256_store_pd(&v[4],vv[1]);
	      
#ifdef _FMA		    
	      ev =  _mm256_load_pd(&extEV[20 * l + 8]);
	      vv[2] = FMAMACC(vv[2], al, ev);		 
#else
	      vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 8])));		  		 
#endif
	      _mm256_store_pd(&v[8],vv[2]);
	      
#ifdef _FMA		    
	      ev =  _mm256_load_pd(&extEV[20 * l + 12]);
	      vv[3] = FMAMACC(vv[3], al, ev);		 
#else
	      vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 12])));		  		 
#endif
	      _mm256_store_pd(&v[12],vv[3]);
	      
#ifdef _FMA		    
	      ev =  _mm256_load_pd(&extEV[20 * l + 16]);
	      vv[4] = FMAMACC(vv[4], al, ev);		 
#else
	      vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 16])));			 	  
#endif
	      _mm256_store_pd(&v[16],vv[4]);		 
	    } 
	}
	
      v = x3_gapColumn;
      scale = 1;
      __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);	 
      
      for(l = 0; scale && (l < 80); l += 4) 
	{
	  __m256d vv = _mm256_load_pd(&v[l]);
	  __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
	  vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
	  if(_mm256_movemask_pd(vv_abs) != 15)
	    scale = 0;	     
	}

      if(scale) 
	{		     	      
	  __m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
	  gapScaling = 1;

	  for(l = 0; l < 80; l += 4) 
	    {
	      __m256d vv = _mm256_load_pd(&v[l]);
	      _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
	    }
	  
	} 
   
     

      for(i = 0; i < n; i++) 
	{   
	  
	  if(x3_gap[i / 32] & mask32[i % 32])
	    {	     
	      if(gapScaling)
		{
		  if(useFastScaling)
		    addScale += wgt[i];
		  else
		    ex3[i]  += 1; 	       
		}
	    }
	  else
	    {
	      if(x1_gap[i / 32] & mask32[i % 32])
		x1 = x1_gapColumn;
	      else
		{
		  x1 = x1_ptr;
		  x1_ptr += 80;
		}

	      if(x2_gap[i / 32] & mask32[i % 32])
		x2 = x2_gapColumn;
	      else
		{
		  x2 = x2_ptr;
		  x2_ptr += 80;
		}	   
	  
	      for(k = 0; k < 4; k++) 
		{
		  vl = &(x1[20 * k]);
		  vr = &(x2[20 * k]);
		  v  = &(x3_ptr[20 * k]);	      	   
		  
		  __m256d vv[5]; 
		  
		  vv[0] = _mm256_setzero_pd();
		  vv[1] = _mm256_setzero_pd();
		  vv[2] = _mm256_setzero_pd();
		  vv[3] = _mm256_setzero_pd();
		  vv[4] = _mm256_setzero_pd();
		  
		  for(l = 0; l < 20; l++) 
		    {		  
		      __m256d al = _mm256_setzero_pd();
		      __m256d ar = _mm256_setzero_pd();
		      
		      __m256d leftv  = _mm256_load_pd(&left[k * 400 + l * 20 + 0]);
		      __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 0]);
		      __m256d vlv = _mm256_load_pd(&vl[0]);
		      __m256d vrv = _mm256_load_pd(&vr[0]);
		      
#ifdef _FMA
		      
		      al = FMAMACC(al, vlv, leftv);
		      ar = FMAMACC(ar, vrv, rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));		  
#endif
		      
		      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 4]);
		      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 4]);
		      vlv = _mm256_load_pd(&vl[4]);
		      vrv = _mm256_load_pd(&vr[4]);
#ifdef _FMA
		      
		      al = FMAMACC(al, vlv, leftv);
		      ar = FMAMACC(ar, vrv, rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
		      
		      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 8]);
		      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 8]);
		      vlv = _mm256_load_pd(&vl[8]);
		      vrv = _mm256_load_pd(&vr[8]);
#ifdef _FMA
		      
		      al = FMAMACC(al, vlv, leftv);
		      ar = FMAMACC(ar, vrv, rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
		      
		      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 12]);
		      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 12]);
		      vlv = _mm256_load_pd(&vl[12]);
		      vrv = _mm256_load_pd(&vr[12]);
#ifdef _FMA
		      
		      al = FMAMACC(al, vlv, leftv);
		      ar = FMAMACC(ar, vrv, rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
		      
		      leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 16]);
		      rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 16]);
		      vlv = _mm256_load_pd(&vl[16]);
		      vrv = _mm256_load_pd(&vr[16]);
		      
#ifdef _FMA		    
		      al = FMAMACC(al, vlv, leftv);
		      ar = FMAMACC(ar, vrv, rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
		      
		      /**************************************************************************************************************/
		      
		      al = hadd3(al);
		      ar = hadd3(ar);
		      al = _mm256_mul_pd(ar,al);
		      
		      /************************************************************************************************************/
#ifdef _FMA		    
		      __m256d ev =  _mm256_load_pd(&extEV[20 * l + 0]);
		      vv[0] = FMAMACC(vv[0], al, ev);		 
#else
		      vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 0])));			  		 		  
#endif
		      _mm256_store_pd(&v[0],vv[0]);
		      
#ifdef _FMA		    
		      ev =  _mm256_load_pd(&extEV[20 * l + 4]);
		      vv[1] = FMAMACC(vv[1], al, ev);		 
#else
		      vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 4])));		  		 
#endif
		      _mm256_store_pd(&v[4],vv[1]);
		      
#ifdef _FMA		    
		      ev =  _mm256_load_pd(&extEV[20 * l + 8]);
		      vv[2] = FMAMACC(vv[2], al, ev);		 
#else
		      vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 8])));		  		 
#endif
		      _mm256_store_pd(&v[8],vv[2]);
		      
#ifdef _FMA		    
		      ev =  _mm256_load_pd(&extEV[20 * l + 12]);
		      vv[3] = FMAMACC(vv[3], al, ev);		 
#else
		      vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 12])));		  		 
#endif
		      _mm256_store_pd(&v[12],vv[3]);
		      
#ifdef _FMA		    
		      ev =  _mm256_load_pd(&extEV[20 * l + 16]);
		      vv[4] = FMAMACC(vv[4], al, ev);		 
#else
		      vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 16])));			 	  
#endif
		      _mm256_store_pd(&v[16],vv[4]);		 
		    }
		}
	      
	      v = x3_ptr;
	      scale = 1;
	      
	      __m256d minlikelihood_avx = _mm256_set1_pd(PLL_MINLIKELIHOOD);	 
	      
	      for(l = 0; scale && (l < 80); l += 4) 
		{
		  __m256d vv = _mm256_load_pd(&v[l]);
		  __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		  vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		  if(_mm256_movemask_pd(vv_abs) != 15)
		    scale = 0;	     
		}
	      
	      if(scale) 
		{		     	      
		  __m256d PLL_TWOTOTHE256v = _mm256_set_pd(PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256,PLL_TWOTOTHE256);
		  for(l = 0; l < 80; l += 4) 
		    {
		      __m256d vv = _mm256_load_pd(&v[l]);
		      _mm256_store_pd(&v[l],_mm256_mul_pd(vv,PLL_TWOTOTHE256v));
		    }
		  if(useFastScaling)
		    addScale += wgt[i];					
		  else
		    ex3[i] += 1;
		}  
	      x3_ptr += 80;
	    }
	}
      break;
    default:
		assert(0);
		break;
	}
 
  if(useFastScaling)
    *scalerIncrement = addScale;
}
