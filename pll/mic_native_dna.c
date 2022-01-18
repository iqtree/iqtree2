#include <omp.h>
#if !defined(__ARM_NEON)
#include <immintrin.h>
#include <string.h>
#include <math.h>

#include "pll.h"
#include "mic_native.h"

static const int states = 4;
static const int statesSquare = 16;
static const int span = 4 * 4;
static const int maxStateValue = 16;

__inline void mic_broadcast16x64(const double* inv, double* outv)
{
    __mmask8 k1 = _mm512_int2mask(0x0F);
    __mmask8 k2 = _mm512_int2mask(0xF0);
    for(int l = 0; l < 16; l += 2)
    {
        __m512d t = _mm512_setzero_pd();
        t = _mm512_mask_extload_pd(t, k1, &inv[(l%4)*4 + l/4], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        t = _mm512_mask_extload_pd(t, k2, &inv[((l+1)%4)*4 + (l+1)/4], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

        _mm512_store_pd(&outv[l*4], t);
    }
}

void newviewGTRGAMMA_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling)
{
    __m512d minlikelihood_MIC = _mm512_set1_pd(PLL_MINLIKELIHOOD);
    __m512d twotothe256_MIC = _mm512_set1_pd(PLL_TWOTOTHE256);
    __m512i absMask_MIC = _mm512_set1_epi64(0x7fffffffffffffffULL);

	int addScale = 0;

    double aEV[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));

    #pragma ivdep
    for (int l = 0; l < 64; ++l)
    {
        aEV[l] = extEV[(l / 16) * 4 + (l % 4)];
    }

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        /* multiply all possible tip state vectors with the respective P-matrices
        */

            double umpX1[256] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double umpX2[256] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            for(int k = 0; k < 256; ++k)
            {
                umpX1[k] = 0.0;
                umpX2[k] = 0.0;
            }

            for(int i = 0; i < maxStateValue; ++i)
            {
              for(int l = 0; l < states; ++l)
              {
                  #pragma ivdep
                  for(int k = 0; k < span; ++k)
                  {
                      umpX1[16 * i + k] +=  tipVector[i * 4 + l] *  left[k * 4 + l];
                      umpX2[16 * i + k] +=  tipVector[i * 4 + l] * right[k * 4 + l];
                  }
              }
            }

        double auX[64] __attribute__((align(64)));

        for(int i = 0; i < n; ++i)
        {
            _mm_prefetch((const char*) (const char*) &x3[span*(i+8)], _MM_HINT_ET1);
            _mm_prefetch((const char*) &x3[span*(i+8) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char*) &x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char*) &x3[span*(i+1) + 8], _MM_HINT_ET0);

            const double *uX1 = &umpX1[16 * tipX1[i]];
            const double *uX2 = &umpX2[16 * tipX2[i]];

            double uX[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double* v = &x3[i * 16];

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v[l] = 0.;
            }

            mic_broadcast16x64(uX, auX);

            for (int j = 0; j < 4; ++j)
            {
                #pragma ivdep
                #pragma vector aligned
                #pragma vector nontemporal
                for(int k = 0; k < 16; ++k)
                {
                    v[k] += auX[j*16 + k] * aEV[j*16 + k];
                }
            }

            // init scaling counter for the site
            if (!fastScaling)
                ex3[i] = 0;

        } // sites loop

      }
      break;
    case PLL_TIP_INNER:
      {
        /* we do analogous pre-computations as above, with the only difference that we now do them
        only for one tip vector */

          double umpX1[256] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        /* precompute P and left tip vector product */

        for(int k = 0; k < 256; ++k)
        {
            umpX1[k] = 0.0;
        }

        for(int i = 0; i < 16; ++i)
        {
          for(int l = 0; l < 4; ++l)
          {
              #pragma ivdep
              for(int k = 0; k < 16; ++k)
              {
                  umpX1[16 * i + k] +=  tipVector[i * 4 + l] *  left[k * 4 + l];
              }
          }
        }

        // re-arrange right matrix for better memory layout
        double aRight[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int l = 0; l < 16; l++)
            {
                aRight[j*16 + l] = right[l*4 + j];
            }
        }

        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char*) &x2[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2[span*(i+16) + 8], _MM_HINT_T1);
            _mm_prefetch((const char*) &x3[span*(i+16)], _MM_HINT_ET1);
            _mm_prefetch((const char*) &x3[span*(i+16) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char*) &x2[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char*) &x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char*) &x3[span*(i+1) + 8], _MM_HINT_ET0);

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX2[l] = 0.;
            }

            double aV2[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            const double* v2 = &(x2[16 * i]);

            mic_broadcast16x64(v2, aV2);

            for(int j = 0; j < 4; j++)
            {
                #pragma ivdep
                #pragma vector aligned
                for(int l = 0; l < 16; l++)
                {
                    uX2[l] += aV2[j*16 + l] * aRight[j*16 + l];
                }
            }

            double* v3 = &(x3[span * i]);

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            double auX[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            mic_broadcast16x64(uX, auX);

            for (int j = 0; j < 4; ++j)
            {
                #pragma ivdep
                #pragma vector aligned
                for(int k = 0; k < 16; ++k)
                {
                    v3[k] += auX[j*16 + k] * aEV[j*16 + k];
                }
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax1 = _mm512_reduce_gmax_pd(t1);
            __m512d t2 = _mm512_load_pd(&v3[8]);
            t2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t2), absMask_MIC));
            double vmax2 = _mm512_reduce_gmax_pd(t2);

            if(vmax1 < PLL_MINLIKELIHOOD && vmax2 < PLL_MINLIKELIHOOD)
            {
				t1 = _mm512_mul_pd(t1, twotothe256_MIC);
				_mm512_store_pd(&v3[0], t1);
				t2 = _mm512_mul_pd(t2, twotothe256_MIC);
				_mm512_store_pd(&v3[8], t2);

                if(!fastScaling)
                  ex3[i] += 1;
                else
                  addScale += wgt[i];
            }
        } // site loop
      }
      break;
    case PLL_INNER_INNER:
    {
      /* same as above, without pre-computations */

        // re-arrange right matrix for better memory layout
        double aLeft[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double aRight[64] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int l = 0; l < 16; l++)
            {
                aLeft[j*16 + l] = left[l*4 + j];
                aRight[j*16 + l] = right[l*4 + j];
            }
        }

        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char*) &x1[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x1[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char*) &x3[span*(i+8)], _MM_HINT_ET1);
            _mm_prefetch((const char*) &x3[span*(i+8) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char*) &x1[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x1[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char*) &x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char*) &x3[span*(i+1) + 8], _MM_HINT_ET0);

            double uX1[16] __attribute__((align(64)));
            double uX2[16] __attribute__((align(64)));
            double uX[16] __attribute__((align(64)));

            for(int l = 0; l < 16; l++)
            {
              uX1[l] = 0.;
              uX2[l] = 0.;
            }

            double aV1[64] __attribute__((align(64)));
            double aV2[64] __attribute__((align(64)));

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            mic_broadcast16x64(v1, aV1);

            mic_broadcast16x64(v2, aV2);

            for(int j = 0; j < 4; j++)
            {
                #pragma ivdep
                #pragma vector aligned
                for(int l = 0; l < 16; l++)
                {
                    uX1[l] += aV1[j*16 + l] * aLeft[j*16 + l];
                    uX2[l] += aV2[j*16 + l] * aRight[j*16 + l];
                }
            }

            double* v3 =  &(x3[span * i]);

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            double auX[64] __attribute__((align(64)));
            mic_broadcast16x64(uX, auX);

            for(int j = 0; j < 4; ++j)
            {
                #pragma ivdep
                #pragma vector aligned
                for(int k = 0; k < 16; ++k)
                {
                    v3[k] += auX[j*16 + k] * aEV[j*16 + k];
                }
            }


            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax1 = _mm512_reduce_gmax_pd(t1);
            __m512d t2 = _mm512_load_pd(&v3[8]);
            t2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t2), absMask_MIC));
            double vmax2 = _mm512_reduce_gmax_pd(t2);

            if(vmax1 < PLL_MINLIKELIHOOD && vmax2 < PLL_MINLIKELIHOOD)
            {
				t1 = _mm512_mul_pd(t1, twotothe256_MIC);
				_mm512_store_pd(&v3[0], t1);
				t2 = _mm512_mul_pd(t2, twotothe256_MIC);
				_mm512_store_pd(&v3[8], t2);

                if(!fastScaling)
                  ex3[i] += 1;
                else
                  addScale += wgt[i];
            }
        }
    } break;
    default:
//      assert(0);
      break;
  }

  /* as above, increment the global counter that counts scaling multiplications by the scaling multiplications
     carried out for computing the likelihood array at node p */

  if (fastScaling)
  {
      *scalerIncrement = addScale;
  }

}

double evaluateGTRGAMMA_MIC(int *ex1, int *ex2, int *wgt,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling)
{
	double sum = 0.0;

    /* the left node is a tip */
    if(tipX1)
    {

        double aTipVec[256] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int k = 0; k < 16; k++)
        {
            for(int l = 0; l < 4; l++)
            {
                aTipVec[k*16 + l] = aTipVec[k*16 + 4 + l] = aTipVec[k*16 + 8 + l] = aTipVec[k*16 + 12 + l] = tipVector[k*4 + l];
            }
        }

        /* loop over the sites of this partition */
        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char*) &x2_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2_start[span*(i+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char*) &x2_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2_start[span*(i+1) + 8], _MM_HINT_T0);

          /* access pre-computed tip vector values via a lookup table */
          const double *x1 = &(aTipVec[16 * tipX1[i]]);
          /* access the other(inner) node at the other end of the branch */
          const double *x2 = &(x2_start[span * i]);

          double term = 0.;

          #pragma ivdep
          #pragma vector aligned
          for(int j = 0; j < span; j++)
              term += x1[j] * x2[j] * diagptable[j];

          if(!fastScaling)
              term = log(0.25 * term) + (ex2[i] * log(PLL_MINLIKELIHOOD));
          else
              term = log(0.25 * term);

          sum += wgt[i] * term;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char*) &x1_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x1_start[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2_start[span*(i+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char*) &x1_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x1_start[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2_start[span*(i+1) + 8], _MM_HINT_T0);

          const double *x1 = &(x1_start[span * i]);
          const double *x2 = &(x2_start[span * i]);

          double term = 0.;

          #pragma ivdep
          #pragma vector aligned
          for(int j = 0; j < span; j++)
              term += x1[j] * x2[j] * diagptable[j];

          if(!fastScaling)
              term = log(0.25 * fabs(term)) + ((ex1[i] + ex2[i]) * log(PLL_MINLIKELIHOOD));
          else
              term = log(0.25 * term);

          sum += wgt[i] * term;
        }
    }

    return sum;
}

void sumGTRGAMMA_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
	double aTipVec[256] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    for(int k = 0; k < 16; k++)
    {
        for(int l = 0; l < 4; l++)
        {
            aTipVec[k*16 + l] = aTipVec[k*16 + 4 + l] = aTipVec[k*16 + 8 + l] = aTipVec[k*16 + 12 + l] = tipVector[k*4 + l];
        }
    }

    switch(tipCase)
    {
      case PLL_TIP_TIP:
      {
        for(int i = 0; i < n; i++)
        {
            const double *left  = &(aTipVec[16 * tipX1[i]]);
            const double *right = &(aTipVec[16 * tipX2[i]]);
            double* sum = &sumtable[i * span];

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
              sum[l] = left[l] * right[l];
            }
        }
      } break;
      case PLL_TIP_INNER:
      {
        for(int i = 0; i < n; i++)
        {
          _mm_prefetch((const char*) &x2_start[span*(i+32)], _MM_HINT_T1);
          _mm_prefetch((const char*) &x2_start[span*(i+32) + 8], _MM_HINT_T1);

          _mm_prefetch((const char*) &x2_start[span*(i+4)], _MM_HINT_T0);
          _mm_prefetch((const char*) &x2_start[span*(i+4) + 8], _MM_HINT_T0);

          const double *left = &(aTipVec[16 * tipX1[i]]);
          const double *right = &(x2_start[span * i]);
          double* sum = &sumtable[i * span];

          #pragma ivdep
          #pragma vector aligned nontemporal
          for(int l = 0; l < span; l++)
          {
              sum[l] = left[l] * right[l];
          }
        }
      } break;
      case PLL_INNER_INNER:
      {
        for(int i = 0; i < n; i++)
        {
            _mm_prefetch((const char*) &x1_start[span*(i+32)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x1_start[span*(i+32) + 8], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2_start[span*(i+32)], _MM_HINT_T1);
            _mm_prefetch((const char*) &x2_start[span*(i+32) + 8], _MM_HINT_T1);

            _mm_prefetch((const char*) &x1_start[span*(i+4)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x1_start[span*(i+4) + 8], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2_start[span*(i+4)], _MM_HINT_T0);
            _mm_prefetch((const char*) &x2_start[span*(i+4) + 8], _MM_HINT_T0);

            const double *left  = &(x1_start[span * i]);
            const double *right = &(x2_start[span * i]);
            double* sum = &sumtable[i * span];

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
                sum[l] = left[l] * right[l];
            }
        }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMA_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt)
{
	double diagptable0[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable1[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable2[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable01[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable02[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*4] = 1.;
        diagptable1[i*4] = 0.;
        diagptable2[i*4] = 0.;

        for(int l = 1; l < states; l++)
        {
          diagptable0[i * 4 + l]  = exp(EIGN[l] * ki * lz);
          diagptable1[i * 4 + l] = EIGN[l] * ki;
          diagptable2[i * 4 + l] = EIGN[l] * EIGN[l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < 16; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % 8 == 0 ? upper / 8 : upper / 8 + 1;

    double dlnLBuf[8] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double d2lnLBuf[8] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    for (int j = 0; j < 8; ++j)
    {
        dlnLBuf[j] = 0.;
        d2lnLBuf[j] = 0.;
    }

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        _mm_prefetch((const char*) &sumtable[i * span * 8], _MM_HINT_T0);
        _mm_prefetch((const char*) &sumtable[i * span * 8 + 8], _MM_HINT_T0);

        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * 8];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[8] __attribute__((align(64)));
        double d1Buf[8] __attribute__((align(64)));
        double d2Buf[8] __attribute__((align(64)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < 8; j++)
        {
            _mm_prefetch((const char*) &sum[span*(j+8)], _MM_HINT_T1);
            _mm_prefetch((const char*) &sum[span*(j+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char*) &sum[span*(j+1)], _MM_HINT_T0);
            _mm_prefetch((const char*) &sum[span*(j+1) + 8], _MM_HINT_T0);

            __m512d d0_1 = _mm512_load_pd(&diagptable0[0]);
            __m512d d0_2 = _mm512_load_pd(&diagptable0[8]);

            __m512d d01_1 = _mm512_load_pd(&diagptable01[0]);
            __m512d d01_2 = _mm512_load_pd(&diagptable01[8]);

            __m512d d02_1 = _mm512_load_pd(&diagptable02[0]);
            __m512d d02_2 = _mm512_load_pd(&diagptable02[8]);

            __m512d s_1 = _mm512_load_pd(&sum[j*16]);
            __m512d s_2 = _mm512_load_pd(&sum[j*16 + 8]);
            __m512d inv_1 = _mm512_mul_pd(d0_1, s_1);
            __m512d d1_1 = _mm512_mul_pd(d01_1, s_1);
            __m512d d2_1 = _mm512_mul_pd(d02_1, s_1);

            __m512d inv_2 = _mm512_fmadd_pd(d0_2, s_2, inv_1);
            __m512d d1_2 = _mm512_fmadd_pd(d01_2, s_2, d1_1);
            __m512d d2_2 = _mm512_fmadd_pd(d02_2, s_2, d2_1);

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_2 = _mm512_add_pd (inv_2, _mm512_swizzle_pd(inv_2, _MM_SWIZ_REG_CDAB));
            inv_2 = _mm512_add_pd (inv_2, _mm512_swizzle_pd(inv_2, _MM_SWIZ_REG_BADC));
            inv_2 = _mm512_add_pd (inv_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_2), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_2);

            d1_2 = _mm512_add_pd (d1_2, _mm512_swizzle_pd(d1_2, _MM_SWIZ_REG_CDAB));
            d1_2 = _mm512_add_pd (d1_2, _mm512_swizzle_pd(d1_2, _MM_SWIZ_REG_BADC));
            d1_2 = _mm512_add_pd (d1_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_2), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_2);

            d2_2 = _mm512_add_pd (d2_2, _mm512_swizzle_pd(d2_2, _MM_SWIZ_REG_CDAB));
            d2_2 = _mm512_add_pd (d2_2, _mm512_swizzle_pd(d2_2, _MM_SWIZ_REG_BADC));
            d2_2 = _mm512_add_pd (d2_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_2), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_2);
        }

        _mm512_store_pd(&invBuf[0], invVec);
        _mm512_store_pd(&d1Buf[0], d1Vec);
        _mm512_store_pd(&d2Buf[0], d2Vec);

        #pragma ivdep
        #pragma vector aligned
        for (int j = 0; j < 8; ++j)
        {
            const double inv_Li = 1.0 / invBuf[j];

            const double d1 = d1Buf[j] * inv_Li;
            const double d2 = d2Buf[j] * inv_Li;

            dlnLBuf[j] += wgt[i * 8 + j] * d1;
            d2lnLBuf[j] += wgt[i * 8 + j] * (d2 - d1 * d1);
        }
    } // site loop

    double dlnLdlz = 0.;
    double d2lnLdlz2 = 0.;
    for (int j = 0; j < 8; ++j)
    {
        dlnLdlz += dlnLBuf[j];
        d2lnLdlz2 += d2lnLBuf[j];
    }

    *ext_dlnLdlz   = dlnLdlz;
    *ext_d2lnLdlz2 = d2lnLdlz2;
}
