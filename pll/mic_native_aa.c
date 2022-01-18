#include <omp.h>
#if !defined(__ARM_NEON)
#include <immintrin.h>
#endif
#include <string.h>
#include <math.h>

#include "pll.h"
#include "mic_native.h"

static const int states = 20;
static const int statesSquare = 20 * 20;
static const int span = 20 * 4;
static const int maxStateValue = 23;

__inline void mic_fma4x80(const double* inv, double* outv, double* mulv)
{
    __mmask8 k1 = _mm512_int2mask(0x0F);
    __mmask8 k2 = _mm512_int2mask(0xF0);
    for(int l = 0; l < 80; l += 40)
    {
        __m512d t = _mm512_setzero_pd();

        t = _mm512_extload_pd(&inv[l], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d m = _mm512_load_pd(&mulv[l]);
        __m512d acc = _mm512_load_pd(&outv[l]);
        __m512d r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l], r);

        m = _mm512_load_pd(&mulv[l + 8]);
        acc = _mm512_load_pd(&outv[l + 8]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 8], r);

        t = _mm512_mask_extload_pd(t, k1, &inv[l], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        t = _mm512_mask_extload_pd(t, k2, &inv[l+20], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

        m = _mm512_load_pd(&mulv[l + 16]);
        acc = _mm512_load_pd(&outv[l + 16]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 16], r);

        t = _mm512_extload_pd(&inv[l+20], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        m = _mm512_load_pd(&mulv[l + 24]);
        acc = _mm512_load_pd(&outv[l + 24]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 24], r);

        m = _mm512_load_pd(&mulv[l + 32]);
        acc = _mm512_load_pd(&outv[l + 32]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 32], r);
    }
}


void newviewGTRGAMMAPROT_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement, const pllBoolean fastScaling)
{
  __m512d minlikelihood_MIC = _mm512_set1_pd(PLL_MINLIKELIHOOD);
  __m512d twotothe256_MIC = _mm512_set1_pd(PLL_TWOTOTHE256);
  __m512i absMask_MIC = _mm512_set1_epi64(0x7fffffffffffffffULL);

  int addScale = 0;

  double aEV[1600] __attribute__((align(PLL_BYTE_ALIGNMENT)));

  #pragma ivdep
  for (int l = 0; l < 1600; ++l)
  {
      aEV[l] = extEV[(l / span) * states + (l % states)];
  }

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        /* multiply all possible tip state vectors with the respective P-matrices
        */

        double umpX1[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double umpX2[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        for(int i = 0; i < maxStateValue; ++i)
        {
          for(int k = 0; k < span; ++k)
          {
              umpX1[i * span + k] = 0.0;
              umpX2[i * span + k] = 0.0;

              #pragma ivdep
              for(int l = 0; l < states; ++l)
              {
                  umpX1[i * span + k] +=  tipVector[i * states + l] *  left[k * states + l];
                  umpX2[i * span + k] +=  tipVector[i * states + l] * right[k * states + l];
              }
          }
        }

        for (int i = 0; i < n; i++)
        {
            const double *uX1 = &umpX1[span * tipX1[i]];
            const double *uX2 = &umpX2[span * tipX2[i]];

            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double* v3 = &x3[i * span];

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
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

          double umpX1[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        /* precompute P and left tip vector product */

        for(int i = 0; i < maxStateValue; ++i)
        {
          for(int k = 0; k < span; ++k)
          {
              umpX1[i * span + k] = 0.0;

              #pragma ivdep
              for(int l = 0; l < states; ++l)
              {
                  umpX1[i * span + k] +=  tipVector[i * states + l] *  left[k * states + l];
              }
          }
        }

        // re-arrange right matrix for better memory layout
        double aRight[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < states; k++)
            {
                for(int l = 0; l < states; l++)
                {
                    aRight[k * span + j * states + l] = right[j * statesSquare +  l * states + k];
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
            }

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            double mx[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = PLL_MAX(vmax, vmax2);
            }

            if (vmax < PLL_MINLIKELIHOOD)
            {
                #pragma vector aligned nontemporal
                for(int l = 0; l < span; l++)
                  v3[l] *= PLL_TWOTOTHE256;

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
        double aLeft[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double aRight[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < states; k++)
            {
                for(int l = 0; l < states; l++)
                {
                    aLeft[k * span + j * states + l] = left[j * statesSquare + l * states + k];
                    aRight[k * span + j * states + l] = right[j * statesSquare + l * states + k];
                }
            }
        }

        for (int i = 0; i < n; i++)
        {

            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T1);
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
            }


            double uX1[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX1[l] = 0.;
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                    _mm_prefetch((const char *)&aLeft[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v1[k], uX1, &aLeft[k * span]);
                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            double mx[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = PLL_MAX(vmax, vmax2);
            }

            if (vmax < PLL_MINLIKELIHOOD)
            {
                #pragma vector aligned nontemporal
                for(int l = 0; l < span; l++)
                  v3[l] *= PLL_TWOTOTHE256;

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

  *scalerIncrement = addScale;

}



double evaluateGTRGAMMAPROT_MIC(int *ex1, int *ex2, int *wgt, double *x1_start, double *x2_start, double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable, const pllBoolean fastScaling)
{
    double sum = 0.0;

    /* the left node is a tip */
    if(tipX1)
    {
        double aTipVec[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int k = 0; k < maxStateValue; k++)
        {
            for(int l = 0; l < states; l++)
            {
                aTipVec[k*span + l] = aTipVec[k*span + states + l] = aTipVec[k*span + 2*states + l] = aTipVec[k*span + 3*states + l] = tipVector[k*states + l];
            }
        }

        /* loop over the sites of this partition */
        for (int i = 0; i < n; i++)
        {
          /* access pre-computed tip vector values via a lookup table */
          const double *x1 = &(aTipVec[span * tipX1[i]]);
          /* access the other(inner) node at the other end of the branch */
          const double *x2 = &(x2_start[span * i]);

          double term = 0.;

          #pragma ivdep
          #pragma vector aligned
          for(int j = 0; j < span; j++) {
              term += x1[j] * x2[j] * diagptable[j];
          }

          if(!fastScaling)
              term = log(0.25 * fabs(term)) + (ex2[i] * log(PLL_MINLIKELIHOOD));
          else
              term = log(0.25 * fabs(term));

          sum += wgt[i] * term;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *) &x1_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x1_start[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &x1_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x1_start[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+1) + 8], _MM_HINT_T0);

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
              term = log(0.25 * fabs(term));

          sum += wgt[i] * term;
        }
    }

    return sum;
}

void sumGTRGAMMAPROT_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
    double aTipVec[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    for(int k = 0; k < maxStateValue; k++)
    {
        for(int l = 0; l < states; l++)
        {
            aTipVec[k*span + l] = aTipVec[k*span + states + l] = aTipVec[k*span + 2*states + l] = aTipVec[k*span + 3*states + l] = tipVector[k*states + l];
        }
    }

    switch(tipCase)
    {
      case PLL_TIP_TIP:
      {
        for(int i = 0; i < n; i++)
        {
            const double *left  = &(aTipVec[span * tipX1[i]]);
            const double *right = &(aTipVec[span * tipX2[i]]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
      case PLL_TIP_INNER:
      {
        for(int i = 0; i < n; i++)
        {
          _mm_prefetch((const char *) &x2_start[span*(i+16)], _MM_HINT_T1);
          _mm_prefetch((const char *) &x2_start[span*(i+16) + 8], _MM_HINT_T1);

          _mm_prefetch((const char *) &x2_start[span*(i+2)], _MM_HINT_T0);
          _mm_prefetch((const char *) &x2_start[span*(i+2) + 8], _MM_HINT_T0);

          const double *left = &(aTipVec[span * tipX1[i]]);
          const double *right = &(x2_start[span * i]);

          #pragma ivdep
          #pragma vector aligned nontemporal
          for(int l = 0; l < span; l++)
          {
              sumtable[i * span + l] = left[l] * right[l];
          }
        }
      } break;
      case PLL_INNER_INNER:
      {
        for(int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *) &x1_start[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x1_start[span*(i+16) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+16) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &x1_start[span*(i+2)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x1_start[span*(i+2) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+2)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+2) + 8], _MM_HINT_T0);

            const double *left  = &(x1_start[span * i]);
            const double *right = &(x2_start[span * i]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMAPROT_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt)
{
    double diagptable0[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable1[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable01[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable02[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*states] = 1.;
        diagptable1[i*states] = 0.;
        diagptable2[i*states] = 0.;

        for(int l = 1; l < states; l++)
        {
          diagptable0[i * states + l]  = exp(EIGN[l] * ki * lz);
          diagptable1[i * states + l] = EIGN[l] * ki;
          diagptable2[i * states + l] = EIGN[l] * EIGN[l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < span; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % PLL_VECTOR_WIDTH == 0 ? upper / PLL_VECTOR_WIDTH : upper / PLL_VECTOR_WIDTH + 1;

    double dlnLdlz = 0.;
    double d2lnLdlz2 = 0.;

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        _mm_prefetch((const char *) &sumtable[i * span * 8], _MM_HINT_T0);
        _mm_prefetch((const char *) &sumtable[i * span * 8 + 8], _MM_HINT_T0);

        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * PLL_VECTOR_WIDTH];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[PLL_VECTOR_WIDTH] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double d1Buf[PLL_VECTOR_WIDTH] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double d2Buf[PLL_VECTOR_WIDTH] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < PLL_VECTOR_WIDTH; j++)
        {
            _mm_prefetch((const char *) &sum[span*(j+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &sum[span*(j+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &sum[span*(j+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &sum[span*(j+1) + 8], _MM_HINT_T0);

            __m512d inv_1 = _mm512_setzero_pd();
            __m512d d1_1 = _mm512_setzero_pd();
            __m512d d2_1 = _mm512_setzero_pd();

            for (int offset = 0; offset < span; offset += 8)
            {
                __m512d d0_1 = _mm512_load_pd(&diagptable0[offset]);
                __m512d d01_1 = _mm512_load_pd(&diagptable01[offset]);
                __m512d d02_1 = _mm512_load_pd(&diagptable02[offset]);
                __m512d s_1 = _mm512_load_pd(&sum[j*span + offset]);

                inv_1 = _mm512_fmadd_pd(d0_1, s_1, inv_1);
                d1_1 = _mm512_fmadd_pd(d01_1, s_1, d1_1);
                d2_1 = _mm512_fmadd_pd(d02_1, s_1, d2_1);
            }

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_CDAB));
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_BADC));
            inv_1 = _mm512_add_pd (inv_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_1), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_1);

            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_CDAB));
            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_BADC));
            d1_1 = _mm512_add_pd (d1_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_1), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_1);

            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_CDAB));
            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_BADC));
            d2_1 = _mm512_add_pd (d2_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_1), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_1);
        }

        _mm512_store_pd(&invBuf[0], invVec);
        _mm512_store_pd(&d1Buf[0], d1Vec);
        _mm512_store_pd(&d2Buf[0], d2Vec);

        #pragma ivdep
        #pragma vector aligned
        for (int j = 0; j < PLL_VECTOR_WIDTH; ++j)
        {
            const double inv_Li = 1.0 / invBuf[j];

            const double d1 = d1Buf[j] * inv_Li;
            const double d2 = d2Buf[j] * inv_Li;

            dlnLdlz += wgt[i * PLL_VECTOR_WIDTH + j] * d1;
            d2lnLdlz2 += wgt[i * PLL_VECTOR_WIDTH + j] * (d2 - d1 * d1);
        }
    } // site loop

    *ext_dlnLdlz   = dlnLdlz;
    *ext_d2lnLdlz2 = d2lnLdlz2;
}


/****
 *       PROTEIN - LG4
 */

void newviewGTRGAMMAPROT_LG4_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement)
{

  __m512d minlikelihood_MIC = _mm512_set1_pd(PLL_MINLIKELIHOOD);
  __m512d twotothe256_MIC = _mm512_set1_pd(PLL_TWOTOTHE256);
  __m512i absMask_MIC = _mm512_set1_epi64(0x7fffffffffffffffULL);

  int addScale = 0;

  double aEV[1600] __attribute__((align(PLL_BYTE_ALIGNMENT)));

  #pragma ivdep
  for (int l = 0; l < 1600; ++l)
  {
      aEV[l] = extEV[(l % span) / states][(l / span) * states + (l % states)];
  }

  switch(tipCase)
  {
    case PLL_TIP_TIP:
      {
        /* multiply all possible tip state vectors with the respective P-matrices
        */

        double umpX1[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double umpX2[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        for(int i = 0; i < 23; ++i)
        {
          for(int k = 0; k < span; ++k)
          {
              umpX1[i * span + k] = 0.0;
              umpX2[i * span + k] = 0.0;
              double *tipv = &(tipVector[k / states][i * states]);


              #pragma ivdep
              for(int l = 0; l < states; ++l)
              {
                  umpX1[i * span + k] +=  tipv[l] *  left[k * states + l];
                  umpX2[i * span + k] +=  tipv[l] * right[k * states + l];
              }
          }
        }

        for (int i = 0; i < n; i++)
        {
            const double *uX1 = &umpX1[span * tipX1[i]];
            const double *uX2 = &umpX2[span * tipX2[i]];

            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double* v3 = &x3[i * span];

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

        } // sites loop
      }
      break;
    case PLL_TIP_INNER:
      {
        /* we do analogous pre-computations as above, with the only difference that we now do them
        only for one tip vector */

          double umpX1[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        /* precompute P and left tip vector product */

        for(int i = 0; i < 23; ++i)
        {
          for(int k = 0; k < span; ++k)
          {
              umpX1[i * span + k] = 0.0;
              double *tipv = &(tipVector[k / states][i * states]);

              #pragma ivdep
              for(int l = 0; l < states; ++l)
              {
                  umpX1[i * span + k] +=  tipv[l] *  left[k * states + l];
              }
          }
        }

        // re-arrange right matrix for better memory layout
        double aRight[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < states; k++)
            {
                for(int l = 0; l < states; l++)
                {
                    aRight[k * span + j * states + l] = right[j * statesSquare +  l * states + k];
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
            }

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
				#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
				#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }


            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            double mx[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = PLL_MAX(vmax, vmax2);
            }

            if (vmax < PLL_MINLIKELIHOOD)
            {
                #pragma vector aligned nontemporal
                for(int l = 0; l < span; l++)
                  v3[l] *= PLL_TWOTOTHE256;

                addScale += wgt[i];
            }
        } // site loop

      }
      break;
    case PLL_INNER_INNER:
    {
      /* same as above, without pre-computations */

        // re-arrange right matrix for better memory layout
        double aLeft[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double aRight[4 * statesSquare] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < states; k++)
            {
                for(int l = 0; l < states; l++)
                {
                    aLeft[k * span + j * states + l] = left[j * statesSquare + l * states + k];
                    aRight[k * span + j * states + l] = right[j * statesSquare + l * states + k];
                }
            }
        }

        for (int i = 0; i < n; i++)
        {

            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T1);
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
            }


            double uX1[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX1[l] = 0.;
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
				#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                    _mm_prefetch((const char *)&aLeft[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v1[k], uX1, &aLeft[k * span]);
                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
				#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            double mx[16] __attribute__((align(PLL_BYTE_ALIGNMENT)));
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = PLL_MAX(vmax, vmax2);
            }

            if (vmax < PLL_MINLIKELIHOOD)
            {
                #pragma vector aligned nontemporal
                for(int l = 0; l < span; l++)
                  v3[l] *= PLL_TWOTOTHE256;

                addScale += wgt[i];
            }
        }
    } break;
    default:
//      assert(0);
      break;
  }

  *scalerIncrement = addScale;

}



double evaluateGTRGAMMAPROT_LG4_MIC(int *wgt, double *x1_start, double *x2_start, double *tipVector[4],
                 unsigned char *tipX1, const int n, double *diagptable)
{
    double sum = 0.0;

    /* the left node is a tip */
    if(tipX1)
    {
        double aTipVec[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        for(int k = 0; k < 23; k++)
        {
            for(int j = 0; j < 4; j++)
            {
				for(int l = 0; l < states; l++)
				{
					aTipVec[k*span + j*states + l] = tipVector[j][k*states + l];
				}
            }
        }

        /* loop over the sites of this partition */
        for (int i = 0; i < n; i++)
        {
			/* access pre-computed tip vector values via a lookup table */
			const double *x1 = &(aTipVec[span * tipX1[i]]);
			/* access the other(inner) node at the other end of the branch */
			const double *x2 = &(x2_start[span * i]);

			#pragma unroll(10)
			for (int k = 0; k < span; k += 8)
			{
				_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
			}

			double term = 0.;

			#pragma ivdep
			#pragma vector aligned
			#pragma noprefetch x2
			for(int j = 0; j < span; j++) {
			  term += x1[j] * x2[j] * diagptable[j];
			}

			term = log(0.25 * fabs(term));

			sum += wgt[i] * term;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
			#pragma unroll(10)
			for (int k = 0; k < span; k += 8)
			{
				_mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

				_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
			}

			const double *x1 = &(x1_start[span * i]);
			const double *x2 = &(x2_start[span * i]);

			double term = 0.;

			#pragma ivdep
			#pragma vector aligned
			#pragma noprefetch x1 x2
			for(int j = 0; j < span; j++)
			  term += x1[j] * x2[j] * diagptable[j];

			term = log(0.25 * fabs(term));

			sum += wgt[i] * term;
        }
    }

    return sum;
}

void sumGTRGAMMAPROT_LG4_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector[4],
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
    double aTipVec[1840] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    for(int k = 0; k < maxStateValue; k++)
    {
        for(int j = 0; j < 4; j++)
        {
			for(int l = 0; l < states; l++)
			{
				aTipVec[k*span + j*states + l] = tipVector[j][k*states + l];
			}
        }
    }

    switch(tipCase)
    {
      case PLL_TIP_TIP:
      {
        for(int i = 0; i < n; i++)
        {
            const double *left  = &(aTipVec[span * tipX1[i]]);
            const double *right = &(aTipVec[span * tipX2[i]]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
      case PLL_TIP_INNER:
      {
        for(int i = 0; i < n; i++)
        {
			#pragma unroll(10)
			for (int k = 0; k < span; k += 8)
			{
				_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
			}

          const double *left = &(aTipVec[span * tipX1[i]]);
          const double *right = &(x2_start[span * i]);

          #pragma ivdep
          #pragma vector aligned nontemporal
		  #pragma noprefetch right
          for(int l = 0; l < span; l++)
          {
              sumtable[i * span + l] = left[l] * right[l];
          }
        }
      } break;
      case PLL_INNER_INNER:
      {
        for(int i = 0; i < n; i++)
        {
			#pragma unroll(10)
			for (int k = 0; k < span; k += 8)
			{
				_mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

				_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
			}

            const double *left  = &(x1_start[span * i]);
            const double *right = &(x2_start[span * i]);

            #pragma ivdep
            #pragma vector aligned nontemporal
			#pragma noprefetch left right
            for(int l = 0; l < span; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMAPROT_LG4_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN[4], double *gammaRates, double lz, int *wgt)
{
    double diagptable0[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable1[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable2[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable01[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));
    double diagptable02[span] __attribute__((align(PLL_BYTE_ALIGNMENT)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*states] = 1.;
        diagptable1[i*states] = 0.;
        diagptable2[i*states] = 0.;

        for(int l = 1; l < states; l++)
        {
          diagptable0[i * states + l]  = exp(EIGN[i][l] * ki * lz);
          diagptable1[i * states + l] = EIGN[i][l] * ki;
          diagptable2[i * states + l] = EIGN[i][l] * EIGN[i][l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < span; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % 8 == 0 ? upper / 8 : upper / 8 + 1;

    double dlnLdlz = 0.;
    double d2lnLdlz2 = 0.;

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * 8];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[8] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double d1Buf[8] __attribute__((align(PLL_BYTE_ALIGNMENT)));
        double d2Buf[8] __attribute__((align(PLL_BYTE_ALIGNMENT)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < 8; j++)
        {

        	#pragma unroll(10)
			for (int k = 0; k < span; k += 8)
			{
				_mm_prefetch((const char *) &sum[span*(j+2) + k], _MM_HINT_T1);
				_mm_prefetch((const char *) &sum[span*(j+1) + k], _MM_HINT_T0);
			}

            __m512d inv_1 = _mm512_setzero_pd();
            __m512d d1_1 = _mm512_setzero_pd();
            __m512d d2_1 = _mm512_setzero_pd();

            for (int offset = 0; offset < span; offset += 8)
            {
                __m512d d0_1 = _mm512_load_pd(&diagptable0[offset]);
                __m512d d01_1 = _mm512_load_pd(&diagptable01[offset]);
                __m512d d02_1 = _mm512_load_pd(&diagptable02[offset]);
                __m512d s_1 = _mm512_load_pd(&sum[j*span + offset]);

                inv_1 = _mm512_fmadd_pd(d0_1, s_1, inv_1);
                d1_1 = _mm512_fmadd_pd(d01_1, s_1, d1_1);
                d2_1 = _mm512_fmadd_pd(d02_1, s_1, d2_1);
            }

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_CDAB));
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_BADC));
            inv_1 = _mm512_add_pd (inv_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_1), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_1);

            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_CDAB));
            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_BADC));
            d1_1 = _mm512_add_pd (d1_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_1), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_1);

            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_CDAB));
            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_BADC));
            d2_1 = _mm512_add_pd (d2_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_1), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_1);
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

            dlnLdlz += wgt[i * 8 + j] * d1;
            d2lnLdlz2 += wgt[i * 8 + j] * (d2 - d1 * d1);
        }
    } // site loop

    *ext_dlnLdlz   = dlnLdlz;
    *ext_d2lnLdlz2 = d2lnLdlz2;
}

