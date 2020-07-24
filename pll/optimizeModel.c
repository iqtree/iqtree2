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
 * @file optimizeModel.c
 *
 * @brief Model optimization routines
 */ 

#include "mem_alloc.h"
#include "systypes.h"

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

static const double MNBRAK_GOLD =    1.618034;          /**< Golden ratio */
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =       1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern char ratesFileName[1024];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[PLL_NUM_PROT_MODELS];

static void optParamGeneric(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType);
// FLAG for easier debugging of model parameter optimization routines 

//#define _DEBUG_MOD_OPT


/*********************FUNCTIONS FOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


/* the following function is used to set rates in the Q matrix 
   the data structure called symmetryVector is used to 
   define the symmetries between rates as they are specified 
   in some of the secondary structure substitution models that 
   generally don't use GTR matrices but more restricted forms thereof */

/** @brief Set a specific rate in the substitition matrix
  *
  * This function is used to set the \a position-th substitution rate of
  * partition \a index to \a rate.
  *
  * @param pr
  *   List of partitions
  *
  * @param model
  *   Index of partition
  *
  * @param rate
  *   The new value to which to set the specific substition rate
  *
  * @param posititon
  *   The number of the substition rate
  */
static void setRateModel(partitionList *pr, int model, double rate, int position)
{
  int
    states   = pr->partitionData[model]->states,
    numRates = (states * states - states) / 2;

  if(pr->partitionData[model]->dataType == PLL_DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(pr->partitionData[model]->dataType != PLL_BINARY_DATA);

  assert(rate >= PLL_RATE_MIN && rate <= PLL_RATE_MAX);

  if(pr->partitionData[model]->nonGTR)
    {    
      int 
        i, 
        index    = pr->partitionData[model]->symmetryVector[position],
        lastRate = pr->partitionData[model]->symmetryVector[numRates - 1];
           
      for(i = 0; i < numRates; i++)
        {       
          if(pr->partitionData[model]->symmetryVector[i] == index)
            {
              if(index == lastRate)
                pr->partitionData[model]->substRates[i] = 1.0;
              else
                pr->partitionData[model]->substRates[i] = rate;      
            }
          
          //printf("%f ", tr->partitionData[model].substRates[i]);
        }
      //printf("\n");
    }
  else
    pr->partitionData[model]->substRates[position] = rate;
}

//LIBRARY: the only thing that we will need to do here is to 
//replace linkList by a string and also add some error correction 
//code

/* 
   the following three functions are used to link/unlink parameters 
   between partitions. This should work in a generic way, however 
   this is so far mainly used for linking unlinking GTR matrix parameter 
   estimates across different protein data partitions.
   Generally this mechanism can also be used for linking/inlinking alpha paremeters 
   between partitions and the like.
   However, all alpha parameter estimates for all partitions and GTR estimates for 
   DNA partitions are unlinked by default. This is actually hard-coded 
   in here. 
*/

/* initializwe a parameter linkage list for a certain parameter type (can be whatever).
   the input is an integer vector that contaions NumberOfModels (numberOfPartitions) elements.

   if we want to have all alpha parameters unlinked and have say 4 partitions the input 
   vector would look like this: {0, 1, 2, 3}, if we want to link partitions 0 and 3 the vector 
   should look like this: {0, 1, 2, 0} 
*/






/* dedicated helper function to initialize the linkage list, that is, essentiaylly compute 
   the integer vector int *linkList used above for linking GTR models.
   
   Once again, this is hard-coded in RAxML, because users can not influence the linking.

*/
   

/* free linkage list data structure */

#define ALPHA_F    0
#define RATE_F     1
#define FREQ_F     2
#define LXRATE_F   3
#define LXWEIGHT_F 4

static void updateWeights(partitionList *pr, int model, int rate, double value)
{
    int j;
    double w = 0.0;
    assert(rate >= 0 && rate < 4);
    pr->partitionData[model]->lg4x_weightExponents[rate] = value;
    for (j = 0; j < 4; j++)
        w += exp(pr->partitionData[model]->lg4x_weightExponents[j]);
    for (j = 0; j < 4; j++)
        pr->partitionData[model]->lg4x_weights[j] = exp(
                pr->partitionData[model]->lg4x_weightExponents[j]) / w;
}

static void optimizeWeights(pllInstance *tr, partitionList *pr, double modelEpsilon, linkageList *ll,
        int numberOfModels)
{
    int i;
    double initialLH = 0.0, finalLH = 0.0;
    pllEvaluateLikelihood(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
    initialLH = tr->likelihood;
    for (i = 0; i < 4; i++)
        optParamGeneric(tr, pr, modelEpsilon, ll, numberOfModels, i, -1000000.0,
                200.0, LXWEIGHT_F);

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
    pllMasterBarrier(tr, pr, PLL_THREAD_COPY_LG4X_RATES);
#endif

    pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
    finalLH = tr->likelihood;
    if (finalLH < initialLH)
        printf("Final: %f initial: %f\n", finalLH, initialLH);
    assert(finalLH >= initialLH);
}

/** @brief Wrapper function for changing a specific model parameter to the specified value
  *
  * Change the \a rateNumber-th model parameter of the type specified by \a whichParameterType to
  * the value \a value.
  * This routine is usually called by model optimization routines to restore the original
  * model parameter vlaue when optimization leads to worse likelihood than the original, or
  * when optimizing routines and testing the new parameter.
  * In case of changing a frequency or substitution rate the Q matrix is also decomposed (into
  * eigenvalues and eigenvectors)
  *
  * @param index
  *   Index of partition
  *
  * @param rateNumber
  *   The index of the model parameter
  *
  * @param value
  *   The value to which the parameter must be changed
  *
  * @param whichParameterType
  *   Type of model parameter. Can be \b RATE_F, \b ALPHA_F or \b FREQ_F, that is substitution rates,
  *   alpha rates, or base frequencies rates
  */   
static void changeModelParameters(int index, int rateNumber, double value, int whichParameterType, pllInstance *tr, partitionList * pr)
{
  switch(whichParameterType)
    {
    case RATE_F:
      setRateModel(pr, index, value, rateNumber);  
      pllInitReversibleGTR(tr, pr, index);          
      break;
    case ALPHA_F:
      pr->partitionData[index]->alpha = value;
      pllMakeGammaCats(pr->partitionData[index]->alpha, pr->partitionData[index]->gammaRates, 4, tr->useMedian);
      break;
    case FREQ_F:
      {
        int 
          states = pr->partitionData[index]->states,
          j;

        double 
          w = 0.0;

        pr->partitionData[index]->freqExponents[rateNumber] = value;

        for(j = 0; j < states; j++)
          w += exp(pr->partitionData[index]->freqExponents[j]);

        for(j = 0; j < states; j++)              
          pr->partitionData[index]->frequencies[j] = exp(pr->partitionData[index]->freqExponents[j]) / w;
        
        pllInitReversibleGTR(tr, pr, index);
      }
      break;
    case LXRATE_F:
        pr->partitionData[index]->gammaRates[rateNumber] = value;
        break;
    case LXWEIGHT_F:
        updateWeights(pr, index, rateNumber, value);
        break;
    default:
      assert(0);
    }
}

/* function that evaluates the change to a parameter */
/** @brief Evaluate the change of a parameter
 *
 *  Evaluate the likelihood for each entry \a i in the linkage list when changing the
 *  \a rateNumber-th parameter of type \a whichFunction (\b ALPHA_F, \b RATE_F 
 *  or \b FREQ_F) to \a value[i]. The resulting likelihood for each entry list \a i in the
 *  linkage list is then stored in \a result[i]
 *
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param rateNumber
 *    Index of the parameter to optimize 
 *
 *  @param value
 *
 *  @param result
 *    An array where the total likelihood of each entry list \a i in the linkage list \a ll  is stored when evaluating the new \a i-th parameter of array \a value
 *
 *  @param converged
 *
 *  @param whichFunction
 *    Type of the model parameter. Possible values are \b ALPHA_F, \b RATE_F and \b FREQ_F
 *
 *  @param numberOfModels
 *    Number of partitions for which we are optimizing 
 *
 *  @param ll
 *    Linkage list
 *
 *  @param modelEpsilon
 *    Epsilon threshold
 */
static void evaluateChange(pllInstance *tr, partitionList *pr, int rateNumber, double *value, double *result, pllBoolean* converged, int whichFunction, int numberOfModels, linkageList *ll, double modelEpsilon)
{ 
  int 
    i, 
    k, 
    pos;

  pllBoolean
    atLeastOnePartition = PLL_FALSE;

  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
        {
          if(converged[pos])
            {
              for(k = 0; k < ll->ld[i].partitions; k++)
                pr->partitionData[ll->ld[i].partitionList[k]]->executeModel = PLL_FALSE;
            }
          else
            {
              atLeastOnePartition = PLL_TRUE;
              for(k = 0; k < ll->ld[i].partitions; k++)
                {
                  int 
                    index = ll->ld[i].partitionList[k];


                  changeModelParameters(index, rateNumber, value[pos], whichFunction, tr, pr);
                }
            }
          pos++;
        }
      else
        {
          for(k = 0; k < ll->ld[i].partitions; k++)
            pr->partitionData[ll->ld[i].partitionList[k]]->executeModel = PLL_FALSE;
        }      
    }

  assert(pos == numberOfModels);

    //some error checks for individual model parameters
    switch (whichFunction)
    {
    case RATE_F:
        assert(rateNumber != -1);
        break;
    case ALPHA_F:
        break;
    case LXRATE_F:
        assert(rateNumber != -1);
        break;
    case LXWEIGHT_F:
        assert(rateNumber != -1);
        break;
    case FREQ_F:
        break;
    default:
        assert(0);
    }

    switch (whichFunction)
    {
    case RATE_F:
    case ALPHA_F:
    case LXRATE_F:
    case FREQ_F:
        pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
        break;
    case LXWEIGHT_F:
        pllEvaluateLikelihood(tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
        break;
    default:
        assert(0);
    }
    //nested optimization for LX4 model, now optimize the weights!
    if (whichFunction == LXRATE_F && atLeastOnePartition)
    {
        pllBoolean *buffer = (pllBoolean*) malloc(
                pr->numberOfPartitions* sizeof(pllBoolean));

        for (i = 0; i < pr->numberOfPartitions; i++) {
            buffer[i] = pr->partitionData[i]->executeModel;
            pr->partitionData[i]->executeModel = PLL_FALSE;
        }

        for (i = 0; i < ll->entries; i++)
        {
            int index = ll->ld[i].partitionList[0];
            if (ll->ld[i].valid)
                pr->partitionData[index]->executeModel = PLL_TRUE;
        }
        optimizeWeights(tr, pr, modelEpsilon, ll, numberOfModels);

        for (i = 0; i < pr->numberOfPartitions; i++) {
            pr->partitionData[i]->executeModel = buffer[i];
        }

        free(buffer);
    }

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

   switch (whichFunction)
    {
      case RATE_F:
        pllMasterBarrier(tr, pr, PLL_THREAD_OPT_RATE);
        break;
      case ALPHA_F:
        pllMasterBarrier(tr, pr, PLL_THREAD_OPT_ALPHA);
        break;
      case FREQ_F:
        pllMasterBarrier(tr, pr, PLL_THREAD_OPT_RATE);
        break;
      case LXRATE_F:
        pllMasterBarrier(tr, pr, PLL_THREAD_OPT_LG4X_RATE);
        break;
      case LXWEIGHT_F:
        pllMasterBarrier(tr, pr, PLL_THREAD_OPT_LG4X_RATE);
        break;
      default:
        break;
    }
#else
   //commented out evaluate below in the course of the LG4X integration
   //pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
#endif     


  for(i = 0, pos = 0; i < ll->entries; i++)     
    {
      if(ll->ld[i].valid)
        {
          result[pos] = 0.0;
          
          for(k = 0; k < ll->ld[i].partitions; k++)
            {
              int 
                index = ll->ld[i].partitionList[k];

              assert(pr->partitionData[index]->partitionLH <= 0.0);
              result[pos] -= pr->partitionData[index]->partitionLH;
              
            }
          pos++;
        }
      for(k = 0; k < ll->ld[i].partitions; k++)
        {
          int index = ll->ld[i].partitionList[k];
          pr->partitionData[index]->executeModel = PLL_TRUE;
        }         
    }
  
  assert(pos == numberOfModels);   
}

/* generic implementation of Brent's algorithm for one-dimensional parameter optimization */

/** @brief Brent's algorithm
 *
 *  Generic implementation of Brent's algorithm for one-dimensional parameter optimization
 *
 *  @param ax
 *
 *  @param bx
 *
 *  @param cx
 *
 *  @param fb
 *
 *  @param tol
 *
 *  @param xmin
 *
 *  @param result
 *
 *  @param numberOfModels
 *    Number of partitions for which we are optimizing 
 *
 *  @param whichFunction
 *    Type of the model parameter. Possible values are \b ALPHA_F, \b RATE_F and \b FREQ_F
 *
 *  @param rateNumber
 *     Index of the parameter to optimize 
 *   
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param ll
 *    Linkage list
 *
 *  @param lim_inf
 *    Lower bound for the rate assignment
 *
 *  @param lim_sup
 *    Upper bound for the rate assignment
 *
 *  @todo
 *     Fill the rest of the entries. Also, why not preallocate all memory instead of allocating
 *     at every call? We can save a lot of time which is lost due to function calls, finding free
 *     memory blocks by allocation strategy, and also prevent mem fragmentation.
 */
static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
                         int whichFunction, int rateNumber, pllInstance *tr, partitionList *pr, linkageList *ll, double lim_inf, double lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *b     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *d     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *p     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *v     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *w     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *x     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *e     = (double *)rax_malloc(sizeof(double) * numberOfModels);
  pllBoolean *converged = (pllBoolean *)rax_malloc(sizeof(pllBoolean) * numberOfModels);
  pllBoolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = PLL_FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf && a[i] <= lim_sup);
      assert(b[i] >= lim_inf && b[i] <= lim_sup);
      assert(x[i] >= lim_inf && x[i] <= lim_sup);
      assert(v[i] >= lim_inf && v[i] <= lim_sup);
      assert(w[i] >= lim_inf && w[i] <= lim_sup);
    }
  
  

  for(iter = 1; iter <= PLL_ITMAX; iter++)
    {
      allConverged = PLL_TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
        allConverged = allConverged && converged[i];

      if(allConverged)
        {
          rax_free(converged);
          rax_free(a);
          rax_free(b);
          rax_free(d);
          rax_free(etemp);
          rax_free(fu);
          rax_free(fv);
          rax_free(fw);
          rax_free(fx);
          rax_free(p);
          rax_free(q);
          rax_free(r);
          rax_free(tol1);
          rax_free(tol2);
          rax_free(u);
          rax_free(v);
          rax_free(w);
          rax_free(x);
          rax_free(xm);
          rax_free(e);
          return;
        }     

      for(i = 0; i < numberOfModels; i++)
        {
          if(!converged[i])
            {                 
              assert(a[i] >= lim_inf && a[i] <= lim_sup);
              assert(b[i] >= lim_inf && b[i] <= lim_sup);
              assert(x[i] >= lim_inf && x[i] <= lim_sup);
              assert(v[i] >= lim_inf && v[i] <= lim_sup);
              assert(w[i] >= lim_inf && w[i] <= lim_sup);
  
              xm[i] = 0.5 * (a[i] + b[i]);
              tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
          
              if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
                {                
                  result[i] =  -fx[i];
                  xmin[i]   = x[i];
                  converged[i] = PLL_TRUE;                
                }
              else
                {
                  if(fabs(e[i]) > tol1[i])
                    {                
                      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
                      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
                      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
                      q[i] = 2.0 * (q[i] - r[i]);
                      if(q[i] > 0.0)
                        p[i] = -p[i];
                      q[i] = fabs(q[i]);
                      etemp[i] = e[i];
                      e[i] = d[i];
                      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
                        d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
                      else
                        {
                          d[i] = p[i] / q[i];
                          u[i] = x[i] + d[i];
                          if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
                            d[i] = PLL_SIGN(tol1[i], xm[i] - x[i]);
                        }
                    }
                  else
                    {                
                      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
                    }
                  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]) : (x[i] + PLL_SIGN(tol1[i], d[i])));
                }

              if(!converged[i])
                assert(u[i] >= lim_inf && u[i] <= lim_sup);
            }
        }
                 
      evaluateChange(tr, pr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll, tol);

      for(i = 0; i < numberOfModels; i++)
        {
          if(!converged[i])
            {
              if(fu[i] <= fx[i])
                {
                  if(u[i] >= x[i])
                    a[i] = x[i];
                  else
                    b[i] = x[i];
                  
                  PLL_SHFT(v[i],w[i],x[i],u[i]);
                  PLL_SHFT(fv[i],fw[i],fx[i],fu[i]);
                }
              else
                {
                  if(u[i] < x[i])
                    a[i] = u[i];
                  else
                    b[i] = u[i];
                  
                  if(fu[i] <= fw[i] || w[i] == x[i])
                    {
                      v[i] = w[i];
                      w[i] = u[i];
                      fv[i] = fw[i];
                      fw[i] = fu[i];
                    }
                  else
                    {
                      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
                        {
                          v[i] = u[i];
                          fv[i] = fu[i];
                        }
                    }       
                }
              
              assert(a[i] >= lim_inf && a[i] <= lim_sup);
              assert(b[i] >= lim_inf && b[i] <= lim_sup);
              assert(x[i] >= lim_inf && x[i] <= lim_sup);
              assert(v[i] >= lim_inf && v[i] <= lim_sup);
              assert(w[i] >= lim_inf && w[i] <= lim_sup);
              assert(u[i] >= lim_inf && u[i] <= lim_sup);
            }
        }
    }

  rax_free(converged);
  rax_free(a);
  rax_free(b);
  rax_free(d);
  rax_free(etemp);
  rax_free(fu);
  rax_free(fv);
  rax_free(fw);
  rax_free(fx);
  rax_free(p);
  rax_free(q);
  rax_free(r);
  rax_free(tol1);
  rax_free(tol2);
  rax_free(u);
  rax_free(v);
  rax_free(w);
  rax_free(x);
  rax_free(xm);
  rax_free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}

/* generic bracketing function required for Brent's algorithm. For details please see the corresponding chapter in the book Numerical Recipees in C */

/** @brief Bracketing function
 *
 *  Generic bracketing function required for Brent's algorithm.
 *  
 *  @param param
 *
 *  @param ax
 *
 *  @param bx
 *
 *  @param cx
 *
 *  @param fa
 *
 *  @param fb
 *
 *  @param fc
 *
 *  @param lim_inf
 *    Lower bound for the rate assignment
 *
 *  @param lim_sup
 *    Upper bound for the rate assignment
 *
 *  @param numberOfModels
 *    Number of partitions for which we are optimizing 
 *
 *  @param rateNumber
 *     Index of the parameter to optimize 
 *
 *  @param whichFunction
 *    Type of the model parameter. Possible values are \b ALPHA_F, \b RATE_F and \b FREQ_F
 *
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param ll
 *    Linkage list
 *
 *  @param modelEpsilon
 *
 *  @return
 *    Fill this
 *
 *  @todo
 *    Fill remaining details
 */
static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
                       double *fc, double lim_inf, double lim_sup, 
                       int numberOfModels, int rateNumber, int whichFunction, pllInstance *tr, partitionList *pr,
                       linkageList *ll, double modelEpsilon)
{
  double 
    *ulim = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)rax_malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)rax_malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)rax_malloc(sizeof(int) * numberOfModels),
    *endState = (int *)rax_malloc(sizeof(int) * numberOfModels);

  pllBoolean *converged = (pllBoolean *)rax_malloc(sizeof(pllBoolean) * numberOfModels);
  pllBoolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = PLL_FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup)    
        param[i] = ax[i] = lim_sup;
      
      if(param[i] < lim_inf) 
        param[i] = ax[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
   
  
  evaluateChange(tr, pr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll, modelEpsilon);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup) 
        param[i] = bx[i] = lim_sup;
      if(param[i] < lim_inf) 
        param[i] = bx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
  evaluateChange(tr, pr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll, modelEpsilon);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
        {         
          PLL_SHFT(dum[i],ax[i],bx[i],dum[i]);
          PLL_SHFT(dum[i],fa[i],fb[i],dum[i]);
        }
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup) 
        param[i] = cx[i] = lim_sup;
      if(param[i] < lim_inf) 
        param[i] = cx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
 
  evaluateChange(tr, pr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll, modelEpsilon);

   while(1) 
     {       
       allConverged = PLL_TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
         allConverged = allConverged && converged[i];

       if(allConverged)
         {
           for(i = 0; i < numberOfModels; i++)
             {         
               if(ax[i] > lim_sup) 
                 ax[i] = lim_sup;
               if(ax[i] < lim_inf) 
                 ax[i] = lim_inf;

               if(bx[i] > lim_sup) 
                 bx[i] = lim_sup;
               if(bx[i] < lim_inf) 
                 bx[i] = lim_inf;
               
               if(cx[i] > lim_sup) 
                 cx[i] = lim_sup;
               if(cx[i] < lim_inf) 
                 cx[i] = lim_inf;
             }

           rax_free(converged);
           rax_free(ulim);
           rax_free(u);
           rax_free(r);
           rax_free(q);
           rax_free(fu);
           rax_free(dum); 
           rax_free(temp);
           rax_free(state);   
           rax_free(endState);
           return 0;
           
         }

       for(i = 0; i < numberOfModels; i++)
         {
           if(!converged[i])
             {
               switch(state[i])
                 {
                 case 0:
                   endState[i] = 0;
                   if(!(fb[i] > fc[i]))                  
                     converged[i] = PLL_TRUE;                                
                   else
                     {
                   
                       if(ax[i] > lim_sup) 
                         ax[i] = lim_sup;
                       if(ax[i] < lim_inf) 
                         ax[i] = lim_inf;
                       if(bx[i] > lim_sup) 
                         bx[i] = lim_sup;
                       if(bx[i] < lim_inf) 
                         bx[i] = lim_inf;
                       if(cx[i] > lim_sup) 
                         cx[i] = lim_sup;
                       if(cx[i] < lim_inf) 
                         cx[i] = lim_inf;
                       
                       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
                       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
                       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
                         (2.0 * PLL_SIGN(PLL_MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
                       
                       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
                       
                       if(u[i] > lim_sup) 
                         u[i] = lim_sup;
                       if(u[i] < lim_inf) 
                         u[i] = lim_inf;
                       if(ulim[i] > lim_sup) 
                         ulim[i] = lim_sup;
                       if(ulim[i] < lim_inf) 
                         ulim[i] = lim_inf;
                       
                       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
                         {
                           param[i] = u[i];
                           if(param[i] > lim_sup)                            
                             param[i] = u[i] = lim_sup;
                           if(param[i] < lim_inf)
                             param[i] = u[i] = lim_inf;
                           endState[i] = 1;
                         }
                       else 
                         {
                           if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
                             {
                               param[i] = u[i];
                               if(param[i] > lim_sup) 
                                 param[i] = u[i] = lim_sup;
                               if(param[i] < lim_inf) 
                                 param[i] = u[i] = lim_inf;
                               endState[i] = 2;
                             }                         
                           else
                             {
                               if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
                                 {
                                   u[i] = ulim[i];
                                   param[i] = u[i];     
                                   if(param[i] > lim_sup) 
                                     param[i] = u[i] = ulim[i] = lim_sup;
                                   if(param[i] < lim_inf) 
                                     param[i] = u[i] = ulim[i] = lim_inf;
                                   endState[i] = 0;
                                 }                              
                               else 
                                 {                
                                   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
                                   param[i] = u[i];
                                   endState[i] = 0;
                                   if(param[i] > lim_sup) 
                                     param[i] = u[i] = lim_sup;
                                   if(param[i] < lim_inf) 
                                     param[i] = u[i] = lim_inf;
                                 }
                             }    
                         }
                     }
                   break;
                 case 1:
                   endState[i] = 0;
                   break;
                 case 2:
                   endState[i] = 3;
                   break;
                 default:
                   assert(0);
                 }
               assert(param[i] >= lim_inf && param[i] <= lim_sup);
             }
         }
             
       evaluateChange(tr, pr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll, modelEpsilon);

       for(i = 0; i < numberOfModels; i++)
         {
           if(!converged[i])
             {         
               switch(endState[i])
                 {
                 case 0:
                   fu[i] = temp[i];
                   PLL_SHFT(ax[i],bx[i],cx[i],u[i]);
                   PLL_SHFT(fa[i],fb[i],fc[i],fu[i]);
                   state[i] = 0;
                   break;
                 case 1:
                   fu[i] = temp[i];
                   if (fu[i] < fc[i]) 
                     {
                       ax[i]=(bx[i]);
                       bx[i]=u[i];
                       fa[i]=(fb[i]);
                       fb[i]=fu[i]; 
                       converged[i] = PLL_TRUE;               
                     } 
                   else 
                     {
                       if (fu[i] > fb[i]) 
                         {
                           assert(u[i] >= lim_inf && u[i] <= lim_sup);
                           cx[i]=u[i];
                           fc[i]=fu[i];
                           converged[i] = PLL_TRUE;                       
                         }
                       else
                         {                 
                           u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
                           param[i] = u[i];
                           if(param[i] > lim_sup) {param[i] = u[i] = lim_sup;}
                           if(param[i] < lim_inf) {param[i] = u[i] = lim_inf;}    
                           state[i] = 1;                 
                         }                
                     }
                   break;
                 case 2: 
                   fu[i] = temp[i];
                   if (fu[i] < fc[i]) 
                     {               
                       PLL_SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
                       state[i] = 2;
                     }     
                   else
                     {
                       state[i] = 0;
                       PLL_SHFT(ax[i],bx[i],cx[i],u[i]);
                       PLL_SHFT(fa[i],fb[i],fc[i],fu[i]);
                     }
                   break;          
                 case 3:                  
                   PLL_SHFT(fb[i],fc[i],fu[i], temp[i]);
                   PLL_SHFT(ax[i],bx[i],cx[i],u[i]);
                   PLL_SHFT(fa[i],fb[i],fc[i],fu[i]);
                   state[i] = 0;
                   break;
                 default:
                   assert(0);
                 }
             }
         }
    }
   

   assert(0);
   rax_free(converged);
   rax_free(ulim);
   rax_free(u);
   rax_free(r);
   rax_free(q);
   rax_free(fu);
   rax_free(dum); 
   rax_free(temp);
   rax_free(state);   
   rax_free(endState);

  

   return(0);
}

/*******************************************************************************************************/
/******** LG4X ***************************************************************************************/

void pllOptLG4X(pllInstance *tr, partitionList * pr, double modelEpsilon,
        linkageList *ll, int numberOfModels)
{
    int i;
    double lg4xScaler, *lg4xScalers = (double *) calloc(pr->numberOfPartitions,
            sizeof(double)), wgtsum = 0.0;
    for (i = 0; i < 4; i++)
        optParamGeneric(tr, pr, modelEpsilon, ll, numberOfModels, i, PLL_LG4X_RATE_MIN,
                PLL_LG4X_RATE_MAX, LXRATE_F);
    for (i = 0; i < pr->numberOfPartitions; i++)
        lg4xScalers[i] = 1.0;
    for (i = 0; i < ll->entries; i++)
    {
        if (ll->ld[i].valid)
        {
            int j, index = ll->ld[i].partitionList[0];
            double averageRate = 0.0;
            assert(ll->ld[i].partitions == 1);
            for (j = 0; j < 4; j++)
                averageRate += pr->partitionData[index]->gammaRates[j];
            averageRate /= 4.0;
            lg4xScalers[index] = averageRate;
        }
    }
    if (pr->numberOfPartitions > 1)
    {
        for (i = 0; i < pr->numberOfPartitions; i++)
            pr->partitionData[i]->fracchange = pr->partitionData[i]->rawFracchange * (1.0 / lg4xScalers[i]);
    }
    for (i = 0; i < pr->numberOfPartitions; i++)
        wgtsum += (double) pr->partitionData[i]->partitionWeight;
    lg4xScaler = 0.0;
    for (i = 0; i < pr->numberOfPartitions; i++)
    {
        double fraction = (double) pr->partitionData[i]->partitionWeight / wgtsum;
        lg4xScaler += (fraction * lg4xScalers[i]);
    }
    tr->fracchange = tr->rawFracchange * (1.0 / lg4xScaler);
    free(lg4xScalers);
}

/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/


//this function is required for implementing the LG4X model later-on 

/** @brief Optimize alpha rates
  *
  * Generic routine for alpha rates optimization
  *
  * @param tr
  *   PLL instance
  *
  * @param pr
  *   List of partitions
  *
  * @param modelEpsilon
  *   Don't know yet
  *
  * @param ll
  *   Linkage list
  *
  * @todo
  *   Implement the LG4X model
  */
void pllOptAlphasGeneric(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    non_LG4X_Partitions = 0,
    LG4X_Partitions  = 0;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do non-LG4X partitions */

  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
        {
        case PLL_DNA_DATA:                          
        case PLL_BINARY_DATA:
        case PLL_SECONDARY_DATA:
        case PLL_SECONDARY_DATA_6:
        case PLL_SECONDARY_DATA_7:
        case PLL_GENERIC_32:
        case PLL_GENERIC_64:
            if (pr->partitionData[ll->ld[i].partitionList[0]]->optimizeAlphaParameter)
            {
                ll->ld[i].valid = PLL_TRUE;
                non_LG4X_Partitions++;
            }
            else
                ll->ld[i].valid = PLL_FALSE;
            break;
        case PLL_AA_DATA:
            if (pr->partitionData[ll->ld[i].partitionList[0]]->optimizeAlphaParameter)
            {
                if (pr->partitionData[ll->ld[i].partitionList[0]]->protModels == PLL_LG4X)
                {
                    LG4X_Partitions++;
                    ll->ld[i].valid = PLL_FALSE;
                }
                else
                {
                    ll->ld[i].valid = PLL_TRUE;
                    non_LG4X_Partitions++;
                }
            }
            else
                ll->ld[i].valid = PLL_FALSE;
            break;
        default:
            assert(0);
        }      
    }   

 

  if(non_LG4X_Partitions > 0)    
    optParamGeneric(tr, pr, modelEpsilon, ll, non_LG4X_Partitions, -1, PLL_ALPHA_MIN, PLL_ALPHA_MAX, ALPHA_F);
  
  /* then LG4x partitions */

  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
        {
        case PLL_DNA_DATA:                          
        case PLL_BINARY_DATA:
        case PLL_SECONDARY_DATA:
        case PLL_SECONDARY_DATA_6:
        case PLL_SECONDARY_DATA_7:
        case PLL_GENERIC_32:
        case PLL_GENERIC_64:
          ll->ld[i].valid = PLL_FALSE;    
          break;
        case PLL_AA_DATA:     
          if(pr->partitionData[ll->ld[i].partitionList[0]]->protModels == PLL_LG4X)
            ll->ld[i].valid = PLL_TRUE;
          else
            ll->ld[i].valid = PLL_FALSE;                    
          break;
        default:
          assert(0);
        }      
    }   
  
  if(LG4X_Partitions > 0)
    pllOptLG4X(tr, pr, modelEpsilon, ll, LG4X_Partitions);

  for(i = 0; ll && i < ll->entries; i++)
    ll->ld[i].valid = PLL_TRUE;
}

/** @brief Optimize model parameters
  *
  * Function for optimizing the \a rateNumber-th model parameter of type \a whichParameterTYpe,
  * i.e. alpha rate, substitution rate, or base frequency rate, in all partitions with the \a
  * valid flag set to \b PLL_TRUE.
  *
  * @param tr
  *   PLL instance
  *
  * @param pr
  *   List of partitions
  *   
  * @param modelEpsilon
  *    A parameter passed for Brent / Brak
  *
  * @param ll
  *   Linkage list
  * 
  * @param numberOfModels
  *   Number of partitions for which we are optimizing 
  *
  * @param rateNumber
  *  Index of the parameter to optimize 
  *
  * @param lim_inf
  *  Lower bound for the rate assignment
  *
  * @param lim_sup
  *  Upper bound for the rate assignment
  *
  * @param whichParameterType
  *  Type of the model parameter. Possible values are \b ALPHA_F, \b RATE_F and \b FREQ_F
  *
  * @todo
  *    Describe the modelEpsilon parameter in detail
  */
static void optParamGeneric(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType)
{
  int
    l,
    k, 
    j, 
    pos;

  double
    *startRates     = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startWeights   = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startExponents = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startValues = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *startLH     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *endLH       = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_a          = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_b          = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_c          = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fa         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fb         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fc         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_param      = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_x          = (double *)rax_malloc(sizeof(double) * numberOfModels);
   
  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
    if (whichParameterType == LXWEIGHT_F)
        pllEvaluateLikelihood (tr, pr, tr->start, PLL_FALSE, PLL_FALSE);
    else
    {
        pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
        if (whichParameterType == LXRATE_F)
        {
            int j;
            for (j = 0; j < pr->numberOfPartitions; j++)
                pr->partitionData[j]->lg4x_weightLikelihood = pr->partitionData[j]->partitionLH;
        }
    }
  
#ifdef  _DEBUG_MOD_OPT
  double
    initialLH = tr->likelihood;
#endif

  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  /* store in startValues the values of the old parameters */
  for(l = 0, pos = 0; ll && l < ll->entries; l++)
    {
      if(ll->ld[l].valid)
        {
          endLH[pos] = PLL_UNLIKELY;
          startLH[pos] = 0.0;

          for(j = 0; j < ll->ld[l].partitions; j++)
            {
              int 
                index = ll->ld[l].partitionList[j];
              
              startLH[pos] += pr->partitionData[index]->partitionLH;
              
              switch(whichParameterType)
                {
                case ALPHA_F:
                  startValues[pos] = pr->partitionData[index]->alpha;
                  break;
                case RATE_F:
                  startValues[pos] = pr->partitionData[index]->substRates[rateNumber];      
                  break;
                case FREQ_F:
                  startValues[pos] = pr->partitionData[index]->freqExponents[rateNumber];
                  break;
                case LXRATE_F:
                    assert(rateNumber >= 0 && rateNumber < 4);
                    startValues[pos] =
                            pr->partitionData[index]->gammaRates[rateNumber];
                    memcpy(&startRates[pos * 4],
                            pr->partitionData[index]->gammaRates,
                            4 * sizeof(double));
                    memcpy(&startExponents[pos * 4],
                            pr->partitionData[index]->lg4x_weightExponents,
                            4 * sizeof(double));
                    memcpy(&startWeights[pos * 4],
                            pr->partitionData[index]->lg4x_weights,
                            4 * sizeof(double));
                    break;
                case LXWEIGHT_F:
                    assert(rateNumber >= 0 && rateNumber < 4);
                    startValues[pos] =
                            pr->partitionData[index]->lg4x_weightExponents[rateNumber];
                    break;
                default:
                  assert(0);
                }
            }
          pos++;
        }
    }  

  assert(pos == numberOfModels);
   
  for(k = 0, pos = 0; ll && k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
        {
          _a[pos] = startValues[pos] + 0.1;
          _b[pos] = startValues[pos] - 0.1;

          if(_a[pos] < lim_inf) 
            _a[pos] = lim_inf;
          
          if(_a[pos] > lim_sup) 
            _a[pos] = lim_sup;
              
          if(_b[pos] < lim_inf) 
            _b[pos] = lim_inf;
          
          if(_b[pos] > lim_sup) 
            _b[pos] = lim_sup;    

          pos++;
        }
    }                                

  assert(pos == numberOfModels);

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, rateNumber, whichParameterType, tr, pr, ll, modelEpsilon);
      
  for(k = 0; k < numberOfModels; k++)
    {
      assert(_a[k] >= lim_inf && _a[k] <= lim_sup);
      assert(_b[k] >= lim_inf && _b[k] <= lim_sup);       
      assert(_c[k] >= lim_inf && _c[k] <= lim_sup);         
    }      

  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, endLH, numberOfModels, whichParameterType, rateNumber, tr,  pr, ll, lim_inf, lim_sup);
        
  for(k = 0, pos = 0; ll && k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
        { 
          if(startLH[pos] > endLH[pos])
            {
              //if the initial likelihood was better than the likelihodo after optimization, we set the values back 
              //to their original values 

              for(j = 0; j < ll->ld[k].partitions; j++)
                {
                  int 
                    index = ll->ld[k].partitionList[j];
                  
                  if (whichParameterType == LXRATE_F)
                    {
                        memcpy(pr->partitionData[index]->lg4x_weights,
                                &startWeights[pos * 4], sizeof(double) * 4);
                        memcpy(pr->partitionData[index]->gammaRates,
                                &startRates[pos * 4], sizeof(double) * 4);
                        memcpy(pr->partitionData[index]->lg4x_weightExponents,
                                &startExponents[pos * 4], 4 * sizeof(double));
                    }

                    changeModelParameters(index, rateNumber, startValues[pos], whichParameterType, tr, pr); 
                }
            }
          else
            {
              //otherwise we set the value to the optimized value 
              //this used to be a bug in standard RAxML, before I fixed it 
              //I was not using _x[pos] as value that needs to be set 

              for(j = 0; j < ll->ld[k].partitions; j++)
                {
                  int 
                    index = ll->ld[k].partitionList[j];
                  
                  changeModelParameters(index, rateNumber, _x[pos], whichParameterType, tr, pr);

                  if (whichParameterType == LXWEIGHT_F)
                    {
                        if (endLH[pos]
                                > pr->partitionData[index]->lg4x_weightLikelihood)
                        {
                            memcpy(pr->partitionData[index]->lg4x_weightsBuffer,
                                    pr->partitionData[index]->lg4x_weights,
                                    sizeof(double) * 4);
                            memcpy(
                                    pr->partitionData[index]->lg4x_weightExponentsBuffer,
                                    pr->partitionData[index]->lg4x_weightExponents,
                                    sizeof(double) * 4);
                            pr->partitionData[index]->lg4x_weightLikelihood =
                                    endLH[pos];
                        }
                    }
                    if (whichParameterType == LXRATE_F)
                    {
                        memcpy(pr->partitionData[index]->lg4x_weights,
                                pr->partitionData[index]->lg4x_weightsBuffer,
                                sizeof(double) * 4);
                        memcpy(pr->partitionData[index]->lg4x_weightExponents,
                                pr->partitionData[index]->lg4x_weightExponentsBuffer,
                                sizeof(double) * 4);
                    }
                }
            }
          pos++;
        }
    }

  #if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      if (whichParameterType == LXRATE_F || whichParameterType == LXWEIGHT_F) {
        pllMasterBarrier(tr, pr, PLL_THREAD_COPY_LG4X_RATES);
      } else {
        pllMasterBarrier(tr, pr, PLL_THREAD_COPY_RATES);
      }

//    switch(whichParameterType)
//      {
//      case FREQ_F:
//      case RATE_F:
//          pllMasterBarrier(tr, pr, PLL_THREAD_COPY_RATES);
//        break;
//      case ALPHA_F:
//          pllMasterBarrier(tr, pr, PLL_THREAD_COPY_ALPHA);
//        break;
//      case LXRATE_F:
//      case LXWEIGHT_F:
//          pllMasterBarrier(tr, pr, PLL_THREAD_COPY_LG4X_RATES);
//        break;
//      default:
//        assert(0);
//      }

  #endif    

    
  assert(pos == numberOfModels);

  rax_free(startLH);
  rax_free(endLH);
  rax_free(_a);
  rax_free(_b);
  rax_free(_c);
  rax_free(_fa);
  rax_free(_fb);
  rax_free(_fc);
  rax_free(_param);
  rax_free(_x);
  rax_free(startValues);
  rax_free(startRates);
  rax_free(startWeights);
  rax_free(startExponents);

#ifdef _DEBUG_MOD_OPT
  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

  if(tr->likelihood < initialLH)
    printf("%f %f\n", tr->likelihood, initialLH);
  assert(tr->likelihood >= initialLH);
#endif
}

//******************** rate optimization functions ***************************************************/

/** @brief Wrapper function for optimizing base frequency rates
  *
  * Wrapper function for optimizing base frequency rates of \a numberOfModels partitions. 
  * The function iteratively calls the function \a optParamGeneric for optimizing each of the \a states
  * parameters
  *
  * @param tr
  *   PLL instance
  *
  * @param pr
  *   List of partitions
  *
  * @param modelEpsilon
  *   Dont know yet
  *
  * @param ll
  *   Linkage list
  *
  * @param numberOfModels
  *   Number of partitions that we are optimizing
  *
  * @param states
  *   Number of states
  */
static void optFreqs(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{ 
  int 
    rateNumber;

  double
    freqMin = -1000000.0,
    freqMax = 200.0;
  
  for(rateNumber = 0; rateNumber < states; rateNumber++)
    optParamGeneric(tr, pr, modelEpsilon, ll, numberOfModels, rateNumber, freqMin, freqMax, FREQ_F);   
}

/** @brief Optimize base frequencies 
 *  
 *  Wrapper function for optimizing base frequencies
 *
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param modelEpsilon
 *    
 *
 *  @param ll
 *    Linkage list
 *
 */
void pllOptBaseFreqs(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    states,
    dnaPartitions = 0,
    aaPartitions  = 0,
    binPartitions = 0;

  /* first do DNA */

  /* Set the valid flag in linkage list to PLL_TRUE for all DNA partitions */
  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
        {
        case PLL_DNA_DATA:  
          states = pr->partitionData[ll->ld[i].partitionList[0]]->states; 
          if(pr->partitionData[ll->ld[i].partitionList[0]]->optimizeBaseFrequencies)
            {
              ll->ld[i].valid = PLL_TRUE;
              dnaPartitions++;              
            }
          else
             ll->ld[i].valid = PLL_FALSE;
          break;       
        case PLL_BINARY_DATA:
        case PLL_AA_DATA:
          ll->ld[i].valid = PLL_FALSE;
          break;
        default:
          assert(0);
        }      
    }   

  /* Optimize the frequency rates of all DNA partitions */
  if(dnaPartitions > 0)
    optFreqs(tr, pr, modelEpsilon, ll, dnaPartitions, states);
  
  /* then AA */

  /* find all partitions that have frequency optimization enabled */ 
  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
        {
        case PLL_AA_DATA:
          states = pr->partitionData[ll->ld[i].partitionList[0]]->states;             
          if(pr->partitionData[ll->ld[i].partitionList[0]]->optimizeBaseFrequencies)
            {
              ll->ld[i].valid = PLL_TRUE;
              aaPartitions++;           
            }
          else
            ll->ld[i].valid = PLL_FALSE; 
          break;
        case PLL_DNA_DATA:      
        case PLL_BINARY_DATA:
          ll->ld[i].valid = PLL_FALSE;
          break;
        default:
          assert(0);
        }        
    }

  if (aaPartitions > 0) {
      optFreqs(tr, pr, modelEpsilon, ll, aaPartitions, states);
  }
  /* then binary */
  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	{
	case PLL_BINARY_DATA:	  
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states; 	      
	  if(pr->partitionData[ll->ld[i].partitionList[0]]->optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = PLL_TRUE;
	      binPartitions++;		
	    }
	  else
	    ll->ld[i].valid = PLL_FALSE; 
	  break;
	case PLL_DNA_DATA:	  
	case PLL_AA_DATA:      
	case PLL_SECONDARY_DATA:
	case PLL_SECONDARY_DATA_6:
	case PLL_SECONDARY_DATA_7:
	case PLL_GENERIC_32:
	case PLL_GENERIC_64:	    
	  ll->ld[i].valid = PLL_FALSE;
	  break;
	default:
	  assert(0);
	}	 
    }

  if(binPartitions > 0)      
    optFreqs(tr, pr, modelEpsilon, ll, binPartitions, states);

  /* done */

  for(i = 0; ll && i < ll->entries; i++)
    ll->ld[i].valid = PLL_TRUE;
}



/* new version for optimizing rates, an external loop that iterates over the rates */
/** @brief Wrapper function for optimizing substitution rates
  *
  * Wrapper function for optimizing substitution rates of \a numberOfModels partitions. 
  * The function determines the  number of free parameters and iteratively calls the 
  * function \a optParamGeneric for optimizing each parameter
  *
  * @param tr
  *   PLL instance
  *
  * @param pr
  *   List of partitions
  *
  * @param modelEpsilon
  *   Dont know yet
  *
  * @param ll
  *   Linkage list
  *
  * @param numberOfModels
  *   Number of partitions that we are optimizing
  *
  * @param states
  *   Number of states
  */
static void optRates(pllInstance *tr, partitionList * pr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, pr, modelEpsilon, ll, numberOfModels, rateNumber, PLL_RATE_MIN, PLL_RATE_MAX, RATE_F);
}


/* figure out if all AA models have been assigned a joint GTR matrix */

/** @brief Check whether all protein partitions have been assigned a joint GTR matrix
  *
  * Check whether there exists at least one protein partition and whether all
  * protein partitions have been assigned a joint GTR matrix.
  *
  * @param pr
  *   List of partitions
  *
  * @return
  *   Return \b PLL_TRUE in case there exists at least one protein partition and all of
  *   protein partitions are assigned a joint GTR matrix. Otherwise return \b PLL_FALSE
  */
static pllBoolean AAisGTR(partitionList *pr)
{
  int i, count = 0;

  for(i = 0; i < pr->numberOfPartitions; i++)
    {
      if(pr->partitionData[i]->dataType == PLL_AA_DATA)
        {
          count++;
          if(pr->partitionData[i]->protModels != PLL_GTR)
            return PLL_FALSE;
        }
    }

  if(count == 0)
    return PLL_FALSE;

  return PLL_TRUE;
}


/* generic substitiution matrix (Q matrix) optimization */

/** @brief Optimize substitution rates
  *
  * Generic routine for substitution matrix (Q matrix) optimization
  *
  * @param tr
  *   PLL instance
  *
  * @param pr
  *   List of partitions
  *
  * @param modelEpsilon
  *   Don't know yet
  *
  * @param ll
  *   Linkage list
  */
void pllOptRatesGeneric(pllInstance *tr, partitionList *pr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitions  = 0,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* 
     first optimize all rates in DNA data partition matrices. That's where we use the valid field in the 
     linkage list data structure. 
   */

  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
        {
          case PLL_DNA_DATA:  
            states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	    if(pr->partitionData[ll->ld[i].partitionList[0]]->optimizeSubstitutionRates)
	      {
		ll->ld[i].valid = PLL_TRUE;
		++ dnaPartitions;  
	      }
	    else	      
	      ll->ld[i].valid = PLL_FALSE;	      
            break;
          case PLL_BINARY_DATA:
          case PLL_AA_DATA:
          case PLL_SECONDARY_DATA:
          case PLL_SECONDARY_DATA_6:
          case PLL_SECONDARY_DATA_7:
          case PLL_GENERIC_32:
          case PLL_GENERIC_64:
            ll->ld[i].valid = PLL_FALSE;
            break;
          default:
            assert(0);
        }      
    }   

  /* if we have dna partitions in our dataset, let's optimize all 5 rates in their substitution matrices */

  if(dnaPartitions > 0)
    optRates(tr, pr, modelEpsilon, ll, dnaPartitions, states);
  
  /* AA partitions evolving under a GTR model do not need to be linked any more, this responsibility now remains 
     with the library user !
   */
  
  for(i = 0; ll && i < ll->entries; i++)
    {
      switch(pr->partitionData[ll->ld[i].partitionList[0]]->dataType)
	{
	case PLL_AA_DATA:
	  states = pr->partitionData[ll->ld[i].partitionList[0]]->states;
	  if(pr->partitionData[ll->ld[i].partitionList[0]]->optimizeSubstitutionRates)
	    {
	      ll->ld[i].valid = PLL_TRUE;
	      aaPartitions++;
	    }
	  else
	    ll->ld[i].valid = PLL_FALSE;
	  break;
	case PLL_DNA_DATA:          
	case PLL_BINARY_DATA:
	case PLL_SECONDARY_DATA:        
	case PLL_SECONDARY_DATA_6:
	case PLL_SECONDARY_DATA_7:
	  ll->ld[i].valid = PLL_FALSE;
	  break;
	default:
	  assert(0);
	}    
    }
  
  if(aaPartitions > 0)
    optRates(tr, pr, modelEpsilon, ll, aaPartitions, states); 

  /* done with all partitions, so we can set all entries in the linkage list to valid again :-) */

  for(i = 0; ll && i < ll->entries; i++)
    ll->ld[i].valid = PLL_TRUE;
}





/*********************FUNCTIONS FOR PSR/CAT model of rate heterogeneity ***************************************/






static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorizePartition(pllInstance *tr, partitionList *pr, rateCategorize *rc, int model, int lower, int upper)
{
  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
      {
        double
          temp = tr->patrat[i];

        int
          found = 0;
        
        for(k = 0; k < pr->partitionData[model]->numberOfCategories; k++)
          {
            if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
              {
                found = 1;
                tr->rateCategory[i] = k; 
                break;
              }
          }
        
        if(!found)
          {
            min = fabs(temp - rc[0].rate);
            tr->rateCategory[i] = 0;

            for(k = 1; k < pr->partitionData[model]->numberOfCategories; k++)
              {
                diff = fabs(temp - rc[k].rate);

                if(diff < min)
                  {
                    min = diff;
                    tr->rateCategory[i] = k;
                  }
              }
          }
      }

  for(k = 0; k < pr->partitionData[model]->numberOfCategories; k++)
    pr->partitionData[model]->perSiteRates[k] = rc[k].rate;
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

void optRateCatPthreads(pllInstance *tr, partitionList *pr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model, 
    i;

  for(model = 0; model < pr->numberOfPartitions; model++)
    {      
      int 
        localIndex = 0;

      pllBoolean 
        execute = ((tr->manyPartitions && isThisMyPartition(pr, tid, model)) || (!tr->manyPartitions));

      if(execute)
        for(i = pr->partitionData[model]->lower;  i < pr->partitionData[model]->upper; i++)
          {
            if(tr->manyPartitions || (i % n == tid))
              {
              
                double initialRate, initialLikelihood, 
                  leftLH, rightLH, leftRate, rightRate, v;
                const double epsilon = 0.00001;
                int k;        
                
                tr->patrat[i] = tr->patratStored[i];     
                initialRate = tr->patrat[i];
                
                initialLikelihood = evaluatePartialGeneric(tr, pr, localIndex, initialRate, model); /* i is real i ??? */
                
                
                leftLH = rightLH = initialLikelihood;
                leftRate = rightRate = initialRate;
                
                k = 1;
                
                while((initialRate - k * lower_spacing > 0.0001) && 
                      ((v = evaluatePartialGeneric(tr, pr, localIndex, initialRate - k * lower_spacing, model))
                       > leftLH) && 
                      (fabs(leftLH - v) > epsilon))  
                  {       
#if !defined(WIN32) && !defined(WIN64)
                    if(isnan(v))
                      assert(0);
#endif
                    
                    leftLH = v;
                    leftRate = initialRate - k * lower_spacing;
                    k++;          
                  }      
                
                k = 1;
                
                while(((v = evaluatePartialGeneric(tr, pr, localIndex, initialRate + k * upper_spacing, model)) > rightLH) &&
                      (fabs(rightLH - v) > epsilon))            
                  {
#if !defined(WIN32) && !defined(WIN64)
                    if(isnan(v))
                      assert(0);
#endif     
                    rightLH = v;
                    rightRate = initialRate + k * upper_spacing;         
                    k++;
                  }           
                
                if(rightLH > initialLikelihood || leftLH > initialLikelihood)
                  {
                    if(rightLH > leftLH)            
                      {      
                        tr->patrat[i] = rightRate;
                        lhs[i] = rightLH;
                      }
                    else
                      {       
                        tr->patrat[i] = leftRate;
                        lhs[i] = leftLH;
                      }
                  }
                else
                  lhs[i] = initialLikelihood;
                
                tr->patratStored[i] = tr->patrat[i];
                localIndex++;
              }
          }
      assert(localIndex == pr->partitionData[model]->width);
    }
}



#else

/** @brief Optimize rates for CAT model
 *
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param model
 *    Partition index
 *
 *  @param lower_specing
 *
 *  @param upper_spacing
 *
 *  @param lhs
 */
static void optRateCatModel(pllInstance *tr, partitionList *pr, int model, double lower_spacing, double upper_spacing, double *lhs)
{
  int lower = pr->partitionData[model]->lower;
  int upper = pr->partitionData[model]->upper;
  int i;
  for(i = lower; i < upper; i++)
    {
      double initialRate, initialLikelihood, 
        leftLH, rightLH, leftRate, rightRate, v;
      const double epsilon = 0.00001;
      int k;
      
      tr->patrat[i] = tr->patratStored[i];     
      initialRate = tr->patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(tr, pr, i, initialRate, model);
      
      
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
            ((v = evaluatePartialGeneric(tr, pr, i, initialRate - k * lower_spacing, model))
             > leftLH) && 
            (fabs(leftLH - v) > epsilon))  
        {         
#if !defined(WIN32) && !defined(WIN64)
          if(isnan(v))
            assert(0);
#endif
          
          leftLH = v;
          leftRate = initialRate - k * lower_spacing;
          k++;    
        }      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(tr, pr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
            (fabs(rightLH - v) > epsilon))      
        {
#if !defined(WIN32) && !defined(WIN64)
          if(isnan(v))
            assert(0);
#endif     
          rightLH = v;
          rightRate = initialRate + k * upper_spacing;   
          k++;
        }           
  
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
        {
          if(rightLH > leftLH)      
            {        
              tr->patrat[i] = rightRate;
              lhs[i] = rightLH;
            }
          else
            {         
              tr->patrat[i] = leftRate;
              lhs[i] = leftLH;
            }
        }
      else
        lhs[i] = initialLikelihood;
      
      tr->patratStored[i] = tr->patrat[i];
    }

}


#endif



/* 
   set scaleRates to PLL_FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

void updatePerSiteRates(pllInstance *tr, partitionList *pr, pllBoolean scaleRates)
{
  int 
    i,
    model;

  if(pr->perGeneBranchLengths && pr->numberOfPartitions > 1)
    {            
      for(model = 0; model < pr->numberOfPartitions; model++)
        {
          int          
            lower = pr->partitionData[model]->lower,
            upper = pr->partitionData[model]->upper;
          
          if(scaleRates)
            {
              double 
                scaler = 0.0,       
                accRat = 0.0; 

              int 
                accWgt     = 0;
              
              for(i = lower; i < upper; i++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
                  
                  accWgt += w;
                  
                  accRat += (w * rate);
                }          
          
              accRat /= ((double)accWgt);
          
              scaler = 1.0 / ((double)accRat);
                  
              for(i = 0; i < pr->partitionData[model]->numberOfCategories; i++)
                pr->partitionData[model]->perSiteRates[i] *= scaler;

              accRat = 0.0;      
              
              for(i = lower; i < upper; i++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);        
                  
                  accRat += (w * rate);
                }                

              accRat /= ((double)accWgt);         

              assert(PLL_ABS(1.0 - accRat) < 1.0E-5);
            }
          else
            {
              double               
                accRat = 0.0; 

              int 
                accWgt     = 0;
              
              for(i = lower; i < upper; i++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
                  
                  accWgt += w;
                  
                  accRat += (w * rate);
                }          
          
              accRat /= ((double)accWgt);
              
              assert(PLL_ABS(1.0 - accRat) < 1.0E-5);
            }

          
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
          {
            int 
              localCount = 0;
            
            for(i = lower, localCount = 0; i < upper; i++, localCount++)
              {               
                pr->partitionData[model]->rateCategory[localCount] = tr->rateCategory[i];
              }
          }
#endif
        }
    }
  else
    {
      int
        accWgt = 0;

      double 
        scaler = 0.0,       
        accRat = 0.0; 

      if(scaleRates)
        {
          for(model = 0, accRat = 0.0, accWgt = 0; model < pr->numberOfPartitions; model++)
            {
              int 
                localCount = 0,
                lower = pr->partitionData[model]->lower,
                upper = pr->partitionData[model]->upper;
              
              for(i = lower, localCount = 0; i < upper; i++, localCount++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
                  
                  accWgt += w;
                  
                  accRat += (w * rate);
                }
            }
          
          accRat /= ((double)accWgt);
          
          scaler = 1.0 / ((double)accRat);
          
          for(model = 0; model < pr->numberOfPartitions; model++)
            {
              for(i = 0; i < pr->partitionData[model]->numberOfCategories; i++)
                pr->partitionData[model]->perSiteRates[i] *= scaler;
            }

          for(model = 0, accRat = 0.0; model < pr->numberOfPartitions; model++)
            {
              int 
                localCount = 0,
                lower = pr->partitionData[model]->lower,
                upper = pr->partitionData[model]->upper;
              
              for(i = lower, localCount = 0; i < upper; i++, localCount++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);        
                  
                  accRat += (w * rate);
                }
            }           

          accRat /= ((double)accWgt);     

          assert(PLL_ABS(1.0 - accRat) < 1.0E-5);
        }
      else
        {
          for(model = 0, accRat = 0.0, accWgt = 0; model < pr->numberOfPartitions; model++)
            {
              int 
                localCount = 0,
                lower = pr->partitionData[model]->lower,
                upper = pr->partitionData[model]->upper;
              
              for(i = lower, localCount = 0; i < upper; i++, localCount++)
                {
                  int 
                    w = tr->aliaswgt[i];
                  
                  double
                    rate = pr->partitionData[model]->perSiteRates[tr->rateCategory[i]];
                  
                  assert(0 <= tr->rateCategory[i] && tr->rateCategory[i] < tr->maxCategories);
                  
                  accWgt += w;
                  
                  accRat += (w * rate);
                }
            }
          
          accRat /=  (double)accWgt;

          assert(PLL_ABS(1.0 - accRat) < 1.0E-5);
        }
         
         /*
       for(model = 0; model < pr->numberOfPartitions; model++)
        {
          int 
            localCount = 0,
            lower = pr->partitionData[model]->lower,
            upper = pr->partitionData[model]->upper;

        }  */       
#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      for(model = 0; model < pr->numberOfPartitions; model++)
        {                        
          int 
            localCount,
            lower = pr->partitionData[model]->lower,
            upper = pr->partitionData[model]->upper;
          
          for(i = lower, localCount = 0; i < upper; i++, localCount++)
              pr->partitionData[model]->rateCategory[localCount] = tr->rateCategory[i];
        }
#endif
    }
  
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  pllMasterBarrier(tr, pr, PLL_THREAD_COPY_RATE_CATS);
#endif               
}

/** @brief Optimize rate categories for CAT model
 *
 *  Optimize rate categories for CAT model
 *
 *  @param tr
 *    PLL instance
 *
 *  @param pr
 *    List of partitions
 *
 *  @param _maxCategories
 *    Number of categories
 */
static void optimizeRateCategories(pllInstance *tr, partitionList *pr, int _maxCategories)
{
  assert(_maxCategories > 0);

  if(_maxCategories > 1)
    {
      double  
        temp,  
        lower_spacing, 
        upper_spacing,
        initialLH = tr->likelihood,     
        *ratStored = (double *)rax_malloc(sizeof(double) * tr->originalCrunchedLength),
        /**lhs =       (double *)malloc(sizeof(double) * tr->originalCrunchedLength),*/
        **oldCategorizedRates = (double **)rax_malloc(sizeof(double *) * pr->numberOfPartitions);

      int  
        i,
        k,
        maxCategories = _maxCategories,
        *oldCategory =  (int *)rax_malloc(sizeof(int) * tr->originalCrunchedLength),
        model,
        *oldNumbers = (int *)rax_malloc(sizeof(int) * pr->numberOfPartitions);
  
      assert(isTip(tr->start->number, tr->mxtips));         
      
      pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

      if(tr->optimizeRateCategoryInvocations == 1)
        {
          lower_spacing = 0.5 / ((double)(tr->optimizeRateCategoryInvocations));
          upper_spacing = 1.0 / ((double)(tr->optimizeRateCategoryInvocations));
        }
      else
        {
          lower_spacing = 0.05 / ((double)(tr->optimizeRateCategoryInvocations));
          upper_spacing = 0.1 / ((double)(tr->optimizeRateCategoryInvocations));
        }
      
      if(lower_spacing < 0.001)
        lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
        upper_spacing = 0.001;
      
      tr->optimizeRateCategoryInvocations = tr->optimizeRateCategoryInvocations + 1;

      memcpy(oldCategory, tr->rateCategory, sizeof(int) * tr->originalCrunchedLength);       
      memcpy(ratStored,   tr->patratStored, sizeof(double) * tr->originalCrunchedLength);

      for(model = 0; model < pr->numberOfPartitions; model++)
        {
          oldNumbers[model]          = pr->partitionData[model]->numberOfCategories;

          oldCategorizedRates[model] = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
          
          memcpy(oldCategorizedRates[model], pr->partitionData[model]->perSiteRates, tr->maxCategories * sizeof(double));
        }      
      
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
      /*tr->lhs = lhs;*/
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      pllMasterBarrier(tr, pr, PLL_THREAD_RATE_CATS);
#else      
      for(model = 0; model < pr->numberOfPartitions; model++)
        optRateCatModel(tr, pr, model, lower_spacing, upper_spacing, tr->lhs);
#endif     

      for(model = 0; model < pr->numberOfPartitions; model++)
        {     
          int 
            where = 1,
            found = 0,
            width = pr->partitionData[model]->upper -  pr->partitionData[model]->lower,
            upper = pr->partitionData[model]->upper,
            lower = pr->partitionData[model]->lower;
            
          rateCategorize 
            *rc = (rateCategorize *)rax_malloc(sizeof(rateCategorize) * width);          
        
          for (i = 0; i < width; i++)
            {
              rc[i].accumulatedSiteLikelihood = 0.0;
              rc[i].rate = 0.0;
            }  
        
          rc[0].accumulatedSiteLikelihood = tr->lhs[lower];
          rc[0].rate = tr->patrat[lower];
        
          tr->rateCategory[lower] = 0;
        
          for (i = lower + 1; i < upper; i++) 
            {
              temp = tr->patrat[i];
              found = 0;
            
              for(k = 0; k < where; k++)
                {
                  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
                    {
                      found = 1;                                                
                      rc[k].accumulatedSiteLikelihood += tr->lhs[i];    
                      break;
                    }
                }
            
              if(!found)
                {           
                  rc[where].rate = temp;            
                  rc[where].accumulatedSiteLikelihood += tr->lhs[i];        
                  where++;
                }
            }
        
          qsort(rc, where, sizeof(rateCategorize), catCompare);
        
          if(where < maxCategories)
            {
              pr->partitionData[model]->numberOfCategories = where;
              categorizePartition(tr, pr, rc, model, lower, upper);
            }
          else
            {
              pr->partitionData[model]->numberOfCategories = maxCategories;
              categorizePartition(tr, pr, rc, model, lower, upper);
            }
        
          rax_free(rc);
        }
                
      updatePerSiteRates(tr, pr, PLL_TRUE);

      pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
      
      if(tr->likelihood < initialLH)
        {                         
          for(model = 0; model < pr->numberOfPartitions; model++)
            {
              pr->partitionData[model]->numberOfCategories = oldNumbers[model];
              memcpy(pr->partitionData[model]->perSiteRates, oldCategorizedRates[model], tr->maxCategories * sizeof(double));
            }         
          
          memcpy(tr->patratStored, ratStored, sizeof(double) * tr->originalCrunchedLength);
          memcpy(tr->rateCategory, oldCategory, sizeof(int) * tr->originalCrunchedLength);           
          
          updatePerSiteRates(tr, pr, PLL_FALSE);
          
          pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

          /* printf("REVERT: %1.40f %1.40f\n", initialLH, tr->likelihood); */

          assert(initialLH == tr->likelihood);
        }
          
      for(model = 0; model < pr->numberOfPartitions; model++)
        rax_free(oldCategorizedRates[model]);
                   
      rax_free(oldCategorizedRates);
      rax_free(oldCategory);
      rax_free(ratStored);       
      /*     rax_free(lhs); */
      rax_free(oldNumbers);
    }
}
  

/************************* end of functions for CAT model of rate heterogeneity */




/*****************************************************************************************************/

/* reset all branche lengths in tree to default values */

/** @brief Reset all branch lengths to default values
  
    Reset all branch lengths in the tree instance to default values (\b PLL_DEFAULTZ)

    @param tr
      PLL instance
  */
void resetBranches(pllInstance *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
	  p->z[0] = PLL_DEFAULTZ;
	  if (tr->perGeneBranchLengths)
        for(i = 1; i < PLL_NUM_BRANCHES; i++)
          p->z[i] = PLL_DEFAULTZ;
        
      q = p->next;
      while(q != p)
        {       
    	  q->z[0] = PLL_DEFAULTZ;
    	  if (tr->perGeneBranchLengths)
            for(i = 1; i < PLL_NUM_BRANCHES; i++)
              q->z[i] = PLL_DEFAULTZ;
          q = q->next;
        }
      p++;
    }
}

/**
 * @brief Adjust frequencies in case some base frequency is close to zero.
 */
static void smoothFrequencies(double *frequencies, int numberOfFrequencies) {
	int countScale = 0, l, loopCounter = 0;

	for (l = 0; l < numberOfFrequencies; l++)
		if (frequencies[l] < PLL_FREQ_MIN)
			countScale++;

	if (countScale > 0) {
		while (countScale > 0) {
			double correction = 0.0, factor = 1.0;

			for (l = 0; l < numberOfFrequencies; l++) {
				if (frequencies[l] == 0.0)
					correction += PLL_FREQ_MIN;
				else if (frequencies[l] < PLL_FREQ_MIN) {
					correction += (PLL_FREQ_MIN - frequencies[l]);
					factor -= (PLL_FREQ_MIN - frequencies[l]);
				}
			}

			countScale = 0;

			for (l = 0; l < numberOfFrequencies; l++) {
				if (frequencies[l] >= PLL_FREQ_MIN)
					frequencies[l] = frequencies[l] - (frequencies[l] * correction * factor);
				else
					frequencies[l] = PLL_FREQ_MIN;

				if (frequencies[l] < PLL_FREQ_MIN)
					countScale++;
			}
			assert(loopCounter < 100);
			loopCounter++;
		}
	}
}

/**
 * @brief Evaluate all possible protein models
 */
static void optimizeProteinModels(pllInstance *tr, partitionList * pr, int *bestIndex, double *bestScores, pllBoolean empiricalFreqs)
{
	int modelIndex, partitionIndex,
	    numProteinModels = PLL_AUTO;

	for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
		bestIndex[partitionIndex] = -1;
		bestScores[partitionIndex] = PLL_UNLIKELY;
	}

	if (empiricalFreqs) {
		double ** freqs = pllBaseFrequenciesInstance(tr, pr);
		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			smoothFrequencies(freqs[partitionIndex], PLL_NUM_AA_STATES);
			memcpy(pr->partitionData[partitionIndex]->empiricalFrequencies, freqs[partitionIndex], PLL_NUM_AA_STATES*sizeof(double));
		}
		free(freqs);
	}

	for (modelIndex = 0; modelIndex < numProteinModels; modelIndex++) {
		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO) {

				pr->partitionData[partitionIndex]->autoProtModels = modelIndex;
				pr->partitionData[partitionIndex]->protUseEmpiricalFreqs =
						empiricalFreqs;

				assert(!pr->partitionData[partitionIndex]->optimizeBaseFrequencies);

				pllInitReversibleGTR(tr, pr, partitionIndex);
			}
		}

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
		pllMasterBarrier (tr, pr, PLL_THREAD_COPY_RATES);
#endif

		/* optimize branch lengths */
		resetBranches(tr);
		pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
		pllOptimizeBranchLengths(tr, pr, 16);

		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO) {
				if (pr->partitionData[partitionIndex]->partitionLH > bestScores[partitionIndex]) {
					/* improved best score */
					bestScores[partitionIndex] = pr->partitionData[partitionIndex]->partitionLH;
					bestIndex[partitionIndex] = modelIndex;
				}
			}
		}
	}
}

/* 
   automatically compute the best protein substitution model for the dataset at hand.
 */

/** @brief Compute the best protein substitution model
  *
  * Automatically compute the best protein substitution model for the dataset
  * at hand
  *
  * @param tr
  *   The PLL instance
  *
  * @param pr
  *   List of partitions
  *
  */
static void autoProtein(pllInstance *tr, partitionList *pr)
{
	int countAutos = 0, partitionIndex;

	/* count the number of partitions with model set to PLL_AUTO */
	for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++)
		if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO)
			countAutos++;

	/* if there are partitions with model set to PLL_AUTO compute the best model */
	if (countAutos > 0) {
		int *bestIndex = (int*) rax_malloc(
				sizeof(int) * pr->numberOfPartitions),
		    *bestIndexEmpFreqs = (int*) rax_malloc(
				sizeof(int) * pr->numberOfPartitions),
		    *oldIndex =
				(int*) rax_malloc(sizeof(int) * pr->numberOfPartitions);

		//pllBoolean *oldFreqs = (pllBoolean*) malloc(  
		//		sizeof(pllBoolean) * pr->numberOfPartitions);
        //  JB 10-Jul-2020 Never used.

		double startLH,
		      *bestScores = (double*) rax_malloc(
				sizeof(double) * pr->numberOfPartitions),
			  *bestScoresEmpFreqs = (double*) rax_malloc(
				sizeof(double) * pr->numberOfPartitions);

		topolRELL_LIST *rl = (topolRELL_LIST *) rax_malloc(
				sizeof(topolRELL_LIST));

		initTL(rl, tr, 1);
		saveTL(rl, tr, 0);

		pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

		/* store the initial likelihood of the tree with the currently assigned protein models */
		startLH = tr->likelihood;

		/* save the currently assigned protein model for each PLL_AUTO partition */
		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			oldIndex[partitionIndex] = pr->partitionData[partitionIndex]->autoProtModels;
			//oldFreqs[partitionIndex] = pr->partitionData[partitionIndex]->protUseEmpiricalFreqs; //JB 10-Jul-2020 Never used
			bestIndex[partitionIndex] = -1;
			bestScores[partitionIndex] = PLL_UNLIKELY;
		}

		/* evaluate all models with fixed base frequencies */
		optimizeProteinModels(tr, pr, bestIndex, bestScores, PLL_FALSE);
		/* evaluate all models with fixed empirical frequencies */
		optimizeProteinModels(tr, pr, bestIndexEmpFreqs, bestScoresEmpFreqs, PLL_TRUE);

		/* model selection */
		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO) {
				int bestIndexFixed = bestIndex[partitionIndex],
				    bestIndexEmp = bestIndexEmpFreqs[partitionIndex];

				double bestLhFixed = bestScores[partitionIndex],
					   bestLhEmp = bestScoresEmpFreqs[partitionIndex],
					   samples = 0.0,
					   freeParamsFixed = 0.0,
					   freeParamsEmp = 0.0;

				samples = pr->partitionData[partitionIndex]->partitionWeight;
				assert(samples > 0.0 && samples >= pr->partitionData[partitionIndex]->width);

				assert(tr->ntips == tr->mxtips);
				freeParamsFixed = freeParamsEmp = (2 * tr->ntips - 3);
				freeParamsEmp += 19.0;

				switch (tr->rateHetModel) {
				case PLL_CAT:
					freeParamsFixed +=
							(double) pr->partitionData[partitionIndex]->numberOfCategories;
					freeParamsEmp +=
							(double) pr->partitionData[partitionIndex]->numberOfCategories;
					break;
				case PLL_GAMMA:
					freeParamsFixed += 1.0;
					freeParamsEmp += 1.0;
					break;
				default:
					assert(0);
				}

				switch (tr->autoProteinSelectionType) {
				case PLL_AUTO_ML:
					if (bestLhFixed > bestLhEmp) {
						pr->partitionData[partitionIndex]->autoProtModels =
								bestIndexFixed;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 0;
					} else {
						pr->partitionData[partitionIndex]->autoProtModels = bestIndexEmp;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 1;
					}
					break;
				case PLL_AUTO_BIC: {
					//BIC: -2 * lnL + k * ln(n)
					double bicFixed = -2.0 * bestLhFixed
							+ freeParamsFixed * log(samples),
						   bicEmp = -2.0
							* bestLhEmp + freeParamsEmp * log(samples);

					if (bicFixed < bicEmp) {
						pr->partitionData[partitionIndex]->autoProtModels =
								bestIndexFixed;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 0;
					} else {
						pr->partitionData[partitionIndex]->autoProtModels = bestIndexEmp;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 1;
					}
				}
					break;
				case PLL_AUTO_AIC: {
					//AIC: 2 * (k - lnL)
					double aicFixed = 2.0 * (freeParamsFixed - bestLhFixed),
							aicEmp = 2.0 * (freeParamsEmp - bestLhEmp);

					if (aicFixed < aicEmp) {
						pr->partitionData[partitionIndex]->autoProtModels =
								bestIndexFixed;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 0;
					} else {
						pr->partitionData[partitionIndex]->autoProtModels = bestIndexEmp;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 1;
					}
				}
					break;
				case PLL_AUTO_AICC: {
					//AICc: AIC + (2 * k * (k + 1))/(n - k - 1)
					double aiccFixed, aiccEmp;

					/*
					 * Even though samples and freeParamsFixed are fp variables, they are actually integers.
					 * That's why we are comparing with a 0.5 threshold.
					 */

					if (fabs(samples - freeParamsFixed - 1.0) < 0.5)
						aiccFixed = 0.0;
					else
						aiccFixed = (2.0 * (freeParamsFixed - bestLhFixed))
								+ ((2.0 * freeParamsFixed
										* (freeParamsFixed + 1.0))
										/ (samples - freeParamsFixed - 1.0));

					if (fabs(samples - freeParamsEmp - 1.0) < 0.5)
						aiccEmp = 0.0;
					else
						aiccEmp = (2.0 * (freeParamsEmp - bestLhEmp))
								+ ((2.0 * freeParamsEmp * (freeParamsEmp + 1.0))
										/ (samples - freeParamsEmp - 1.0));

					if (aiccFixed < aiccEmp) {
						pr->partitionData[partitionIndex]->autoProtModels =
								bestIndexFixed;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 0;
					} else {
						pr->partitionData[partitionIndex]->autoProtModels = bestIndexEmp;
						pr->partitionData[partitionIndex]->protUseEmpiricalFreqs = 1;
					}
				}
					break;
				default:
					assert(0);
				}

				pllInitReversibleGTR(tr, pr, partitionIndex);
			}
		}

		resetBranches(tr);
		pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
		pllOptimizeBranchLengths(tr, pr, 64);

		/* set the protein model of PLL_AUTO partitions to the best computed and reset model parameters */
		for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
			if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO) {
				pr->partitionData[partitionIndex]->autoProtModels = bestIndex[partitionIndex];
				pllInitReversibleGTR(tr, pr, partitionIndex);
			}
		}

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
		pllMasterBarrier(tr, pr, PLL_THREAD_COPY_RATES);
#endif

		/* compute again the likelihood of the tree */
		resetBranches(tr);
		pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
		pllOptimizeBranchLengths(tr, pr, 64);

		/* check if the likelihood of the tree with the new protein models assigned to PLL_AUTO partitions is better than the with the old protein models */
		if (tr->likelihood < startLH) {
			for (partitionIndex = 0; partitionIndex < pr->numberOfPartitions; partitionIndex++) {
				if (pr->partitionData[partitionIndex]->protModels == PLL_AUTO) {
					pr->partitionData[partitionIndex]->autoProtModels = oldIndex[partitionIndex];
					pllInitReversibleGTR(tr, pr, partitionIndex);
				}
			}

			//this barrier needs to be called in the library
			//#ifdef _USE_PTHREADS
			//pllMasterBarrier(tr, pr, PLL_THREAD_COPY_RATES);
			//#endif

			/* Restore the topology. rl holds the topology before the optimization. However,
			 since the topology doesn't change - only the branch lengths do - maybe we
			 could write a new routine that will store only the branch lengths and restore them */
			restoreTL(rl, tr, 0,
					pr->perGeneBranchLengths ? pr->numberOfPartitions : 1);
			pllEvaluateLikelihood(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
		}

		assert(tr->likelihood >= startLH);

		freeTL(rl);
		rax_free(rl);

		rax_free(oldIndex);
		rax_free(bestIndex);
		rax_free(bestIndexEmpFreqs);
		rax_free(bestScores);
		rax_free(bestScoresEmpFreqs);
	}
}


/* iterative procedure for optimizing all model parameters */

/* @brief Optimize all model parameters
 *
 * Iterative procedure for optimizing all model parameters
 *
 * @param tr
 *   PLL instance
 *
 * @param pr
 *   List of partitions
 *
 * @param likelihoodEpsilon
 *   Optimize model parameters until we get a difference of \a likelihoodEpsilon
 *
 * @todo
 *   Describe likelihoodEpsilon. Understand the TODO marked blocks.
 */
void modOpt(pllInstance *tr, partitionList *pr, double likelihoodEpsilon)
{ 
  int catOpt = 0; 
  double 
    inputLikelihood,
    currentLikelihood,
    modelEpsilon = 0.0001;

  /* linkage lists for alpha, p-invar has actually been ommitted in this version of the code 
     and the GTR subst matrices */

  linkageList
    *alphaList = pr->alphaList,
    *rateList  = pr->rateList,
    *freqList  = pr->freqList;

  modelEpsilon = 0.0001;

  // test code for library
  if (0)
   {
     
      //assuming that we have three partitions for testing here 

      //alphaList = initLinkageListString("0,1,2", pr);
      //rateList  = initLinkageListString("0,1,1", pr);
    
      //init_Q_MatrixSymmetries("0,1,2,3,4,5", pr, 0);
      //init_Q_MatrixSymmetries("0,1,2,3,4,4", pr, 1);
      //init_Q_MatrixSymmetries("0,1,1,2,3,4", pr, 2);
      
      //function that checks that partitions that have linked Q matrices as in our example above
      //will not have different configurations of the Q matrix as set by the init_Q_MatrixSymmetries() function
      //e.g., on would have HKY and one would have GTR, while the user claimes that they are linked
      //in our example, the Q matrices of partitions 1 and 2 are linked 
      //but we set different matrix symmetries via 
      // init_Q_MatrixSymmetries("0,1,2,3,4,4", tr, 1);
      // and
      // init_Q_MatrixSymmetries("0,1,1,2,3,4", tr, 2);
      //
      //the function just let's assertions fail for the time being .....

      //checkMatrixSymnmetriesAndLinkage(pr, rateList);

  /* alpha parameters and p-invar parameters are unlinked.
     this is the point where I actually hard-coded this in RAxML */

  /* call the dedicated function for linking the GTR matrix across all AA data partitions 
     If we have only DNA data all GTR matrix estimates will be unlinked.
     */
   }
  else
   {
     //alphaList = initLinkageList(unlinked, pr);
     //freqList  = initLinkageList(unlinked, pr);
     //rateList  = initLinkageListGTR(pr);
   }

  tr->start = tr->nodep[1];

  /* This check is here to make sure that the likelihood 
     computed prior to entering modOpt() is consistent 
     with the likelihood when entering modOpt().
     This allows us to ensure that we didn't forget to update anything prior 
     to entereing this function.
   */
  inputLikelihood = tr->likelihood;
  pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
  assert (inputLikelihood == tr->likelihood);

  do
  {           
    //printBothOpen("cur LH: %f\n", tr->likelihood);
    currentLikelihood = tr->likelihood;     

#ifdef _DEBUG_MOD_OPT
      printf ("start: %f\n", currentLikelihood);
#endif

    pllOptRatesGeneric(tr, pr, modelEpsilon, rateList);

    pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef _DEBUG_MOD_OPT
    printf ("after rates %f\n", tr->likelihood);
#endif

    autoProtein(tr, pr);

    pllOptimizeBranchLengths(tr, pr, 2); // 0.0625 * 32 = 2.0

#ifdef _DEBUG_MOD_OPT
    pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
    printf("after br-len 1 %f\n", tr->likelihood); 
#endif

    pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

    pllOptBaseFreqs(tr, pr, modelEpsilon, freqList);
    
    pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
    
    pllOptimizeBranchLengths(tr, pr, 2); // 0.0625 * 32 = 2.0

#ifdef _DEBUG_MOD_OPT
    pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE); 
    printf("after pllOptBaseFreqs 1 %f\n", tr->likelihood);
#endif 

    switch(tr->rateHetModel)
    {
      case PLL_GAMMA:      
        pllOptAlphasGeneric (tr, pr, modelEpsilon, alphaList);
        pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);

#ifdef _DEBUG_MOD_OPT
          printf("after alphas %f\n", tr->likelihood); 
#endif

        pllOptimizeBranchLengths(tr, pr, 3); // 0.1 * 32 = 3.2

#ifdef _DEBUG_MOD_OPT
          pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);  
          printf("after br-len 2 %f\n", tr->likelihood); 
#endif
        break;
      case PLL_CAT:
        if(catOpt < 3)
        {                            
          pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);  
          optimizeRateCategories(tr, pr, tr->categories);
#ifdef _DEBUG_MOD_OPT
            pllEvaluateLikelihood (tr, pr, tr->start, PLL_TRUE, PLL_FALSE);  
            printf("after cat-opt %f\n", tr->likelihood); 
#endif
          catOpt++;
        }
        break;    
      default:
        assert(0);
    }                   

    if(tr->likelihood < currentLikelihood)
     {
      printf("%.20f %.20f\n", tr->likelihood, currentLikelihood);
      printf("Difference: %.20f\n",tr->likelihood - currentLikelihood);
    }
    assert (tr->likelihood - currentLikelihood > 0.000000000000001);
    //assert(tr->likelihood > currentLikelihood);

  }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  
  /* TODO: Why do we check the computed likelihood with the currentLikelihood which is the likelihood before THIS optimization loop? Why dont we
     rather check it with the initial likelihood (the one before calling modOpt)? Isn't it possible to have a deadlock? */

  
}

