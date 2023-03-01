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
 * @file genericParallelization.c
 */
#include "mem_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef MEASURE_TIME_PARALLEL
#include <time.h>
#endif

#include <assert.h>

#include "genericParallelization.h"
#include "pllInternal.h"
#include "pll.h"
// BQM: this causes compiling error when MPI was not installed
// GMJB: Fixed as per below
#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
  #include <mpi.h>
#endif
/** @file genericParallelization.c
    
    @brief Generic master-worker parallelization with either pthreads or MPI. 
    
    Worker threads/processes mostly work on a local
    tree. Implementationwise, MPI operations are abstracted as good as
    possible via defines (that translate to no-ops or memcpy-calls in
    the pthreads version).

    @todo the code still contains many memory copy operations that
    could be executed more efficiently in-place  
*/



void perSiteLogLikelihoodsPthreads(pllInstance *tr, partitionList *pr, double *lhs, int n, int tid);
void broadcastAfterRateOpt(pllInstance *tr, pllInstance *localTree, partitionList *pr, int n, int tid);
void branchLength_parallelReduce(pllInstance *tr, double *dlnLdlz,  double *d2lnLdlz2, int numBranches );
void pllMasterPostBarrier(pllInstance *tr, partitionList *pr, int jobType);
static void distributeYVectors(pllInstance *localTree, pllInstance *tr, partitionList *localPr);
static void distributeWeights(pllInstance *localTree, pllInstance *tr, partitionList *localPr);
static pllBoolean execFunction(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n);

static void *likelihoodThread(void *tData); 

static void multiprocessorScheduling(pllInstance * tr, partitionList *pr, int tid);

static void computeFraction(partitionList *localPr, int tid, int n);
static void computeFractionMany(partitionList *localPr, int tid);
static void initializePartitionsMaster(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n);

#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
static char* addBytes(char *buf, void *toAdd, size_t numBytes); 
static char* popBytes(char *buf, void *result, size_t numBytes); 
static void defineTraversalInfoMPI(void);
static pllBoolean pllWorkerTrap(pllInstance *tr, partitionList *pr);
#endif

#ifdef _USE_PTHREADS
static pthread_t *threads;
static threadData *tData;
#endif

extern volatile int jobCycle; 
extern volatile int threadJob;          /**< current job to be done by worker threads/processes */
extern pllBoolean treeIsInitialized; 

#ifdef MEASURE_TIME_PARALLEL
extern double masterTimePerPhase; 
double timeBuffer[NUM_PAR_JOBS]; 
double timePerRegion[NUM_PAR_JOBS]; 
#endif

extern const char* getJobName(int tmp); 

//extern double *globalResult; 
extern volatile char *barrierBuffer;


#if defined(_FINE_GRAIN_MPI) || ( defined(_IQTREE_MPI) && defined(CLANG_UNDER_VS) )
extern MPI_Datatype TRAVERSAL_MPI; 

/** @brief Pthreads helper function for adding bytes to communication buffer.

    Copy from \toAdd to \a buf \a numBytes bytes

    @param buf
      Where to place bytes

    @pram toAdd
      Where to copy them from

    @para numBytes
      How many to copy

    @return
      Pointer to the end of placed data in communication buffer (first free slot)
 */ 
static char* addBytes(char *buf, void *toAdd, size_t numBytes)
{
  memcpy(buf, toAdd, numBytes);  
  return buf + numBytes;  
}

/** @brief Pthreads helper function for removing bytes from communication buffer
    
    Copies \a numBytes from communication buffer \a buf to some local buffer \a buf

    @param buf
      Where to store the bytes

    @param result
      Where to copy from

    @param numBytes
      How many to copy
    
    @return
      Pointer to the end of read data in communication buffer (first free slot)
 */ 
static char* popBytes(char *buf, void *result, size_t numBytes)
{
  memcpy(result, buf, numBytes); 
  return buf + numBytes;   
}

/** @brief Lock the MPI slave processes prior allocating partitions

    MPI slave processes are locked and wait until the master process
    has read the number of partitions, which it then broadcasts
    to slaves, effectively unlocking them. The slave processes will
    then allocate their own data structures and be locked in the
    likelihood function.

    @param tr
      PLL instance
    
    @todo
      This function should not be called by the user. It is called
      at \a pllCreateInstance. Probably this function should be removed
      and inline code be placed in \a pllCreateInstance.
*/
void pllLockMPI (pllInstance * tr)
{
  int numberOfPartitions;
  partitionList * pr;

  if (!MASTER_P) 
   {
     //MPI_Bcast (&numberOfPartitions, 1, MPI_INT, MPI_ROOT, MPI_COMM_WORLD);
     MPI_Bcast (&numberOfPartitions, 1, MPI_INT, 0, MPI_COMM_WORLD);
     pr = (partitionList *) rax_calloc (1, sizeof (partitionList));
     pr->numberOfPartitions = numberOfPartitions;

     pllWorkerTrap (tr, pr);
     MPI_Barrier (MPI_COMM_WORLD);
     MPI_Finalize ();
     exit(0);
   }
}

/** Finalize MPI run

    Finalizes MPI run by synchronizing all processes (master + slaves) with a
    barrier so that all free their allocated resources. Then \a MPI_Finalize ()
    is called.

    @todo
      Similarly as with the \a pllLockMPI function, this should not be called
      by the user, but it is called implicitly at the end of \a pllDestroyInstance.
      Probably this function should be removed and inline code be placed in
      \a pllDestroyInstance.
*/
void pllFinalizeMPI (void)
{
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize ();
}

/**
   @brief Sets up the MPI environment.  

   Calls the \a MPI_Init function and makes sure all processes store
   their process ID and the total number of processes, using a barrier.
   
   @note this should be the first call that is executed in your main
   method.
   
   @param argc   
     Address of argc from main
   @param argv   
     Address of argv from main
 */
void pllInitMPI(int * argc, char **argv[])
{  
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);

  /* if(MASTER_P) */
  /*   printf("\nThis is RAxML Process Number: %d (MASTER)\n", processID); */
  MPI_Barrier(MPI_COMM_WORLD);

}


/**
   @brief Traps worker MPI processes.    
   
   @note  This function should be called immediately after initMPI()

   @param tr 
     PLL instance 

   @param pr
     List of partitions

   @return
     Returns /b PLL_FALSE if the callee was the master thread/process, otherwise /b PLL_TRUE
 */ 
static pllBoolean pllWorkerTrap(pllInstance *tr, partitionList *pr)
{
  /// @note for the broadcasting, we need to, if the tree structure has already been initialized 
  treeIsInitialized = PLL_FALSE; 

  if(NOT MASTER_P) 
    {
      threadData tData; 
      tData.tr = tr; 
      tData.threadNumber = processID;
      tData.pr = pr;
      
      likelihoodThread(&tData);

      /* notice: the next call MUST be the return call from the main method */
      return PLL_TRUE; 
    }
  return PLL_FALSE; 
}


#define ELEMS_IN_TRAV_INFO  9
/** @brief Create a datastructure for sending the traversal descriptor.
    
    @note This seems to be a very safe method to define your own mpi
   datatypes (often there are problems with padding). But it is not
   entirely for the weak of heart...
 */ 
static void defineTraversalInfoMPI (void)
{
  MPI_Datatype *result  = &TRAVERSAL_MPI; 

  int i ; 
  MPI_Aint base; 
  int blocklen[ELEMS_IN_TRAV_INFO+1] = {1, 1, 1, 1, PLL_NUM_BRANCHES, PLL_NUM_BRANCHES, 1,1,1,1}; 
  MPI_Aint disp[ELEMS_IN_TRAV_INFO+1];
  MPI_Datatype type[ELEMS_IN_TRAV_INFO+1] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_UB}; 
  traversalInfo desc[2]; 

  MPI_Get_address( desc, disp);
  MPI_Get_address( &(desc[0].pNumber), disp + 1 );
  MPI_Get_address( &(desc[0].qNumber), disp + 2 );  
  MPI_Get_address( &(desc[0].rNumber), disp + 3); 
  MPI_Get_address( desc[0].qz, disp + 4 );
  MPI_Get_address( desc[0].rz, disp + 5 );
  MPI_Get_address( &(desc[0].slot_p), disp + 6);
  MPI_Get_address( &(desc[0].slot_q), disp + 7);
  MPI_Get_address( &(desc[0].slot_r), disp + 8);
  MPI_Get_address( desc + 1, disp + 9);

  base = disp[0]; 
  for(i = 0; i < ELEMS_IN_TRAV_INFO+1; ++i)
    disp[i] -= base;

  MPI_Type_create_struct( ELEMS_IN_TRAV_INFO+1 , blocklen, disp, type, result);
  MPI_Type_commit(result);
}


#endif


/********************/
/* PTHREAD-SPECIFIC */
/********************/
#ifdef _USE_PTHREADS

#ifndef _PORTABLE_PTHREADS
/** @brief Pins a thread to a core (for efficiency). 

    This is a non-portable function that works only on some linux distributions of pthreads.
    It sets the affinity of each thread to a specific core so that the performance is not
    degraded due to threads migration.

    @note 
      It is only called if \a _PORTABLE_PTHREADS is not defined

    @param tid the thread id
 */ 
void pinToCore(int tid)
{
  static int nextCore = 0;

  cpu_set_t cpuset;

  CPU_ZERO(&cpuset);    
  CPU_SET(nextCore++, &cpuset);

  if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
    {
      assert(0);
    }
}
#endif

/**  Start PThreads

     Start JOINABLE threads by executing \a pthread_create. The threads
     are attached to the \a pllLikelihoodThread function

     @param tr
       PLL instance

     @param pr
       List of partitions

     @todo
       This function should never be called by the user. It is called
       implicitly at \a pllInitModel. Perhaps we should add a check
       or inline the code
 */ 
void pllStartPthreads (pllInstance *tr, partitionList *pr)
{
  pthread_attr_t attr;
  int rc, t;
  treeIsInitialized = PLL_FALSE; 

  jobCycle        = 0;
  threadJob       = 0;

  /* printf("\nThis is the RAxML Master Pthread\n");   */

#if (NOT defined(_USE_PTHREADS) && defined( MEASURE_TIME_PARALLEL))
  timeBuffer = rax_calloc(NUM_PAR_JOBS * tr->numberOfThreads, sizeof(double)); 
#endif

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  threads    = (pthread_t *)rax_malloc((size_t)tr->numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)rax_malloc((size_t)tr->numberOfThreads * sizeof(threadData));

  barrierBuffer            = (volatile char *)  rax_malloc(sizeof(volatile char)   *  (size_t)tr->numberOfThreads);

  for(t = 0; t < tr->numberOfThreads; t++)
    barrierBuffer[t] = 0;

  for(t = 1; t < tr->numberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].pr  = pr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
  pthread_attr_destroy (&attr);
}

/** Stop PThread
    
    Stop threads by \a pthread_join

    @param  tr
      PLL instance

    @todo
      This function should never be called by the user. It is implicitly called
      at \a pllPartitionsDestroy. We should inline the code
*/
void pllStopPthreads (pllInstance * tr)
{
  int i;

  for (i = 1; i < tr->numberOfThreads; ++ i)
   {
     pthread_join (threads[i], NULL);
   }
 
  rax_free (threads);
  rax_free (tData);
  rax_free ((void *)barrierBuffer);
  rax_free (globalResult);

}
#endif


/** Compute per-site log likelihoods (PThreads version) 

    Worker threads evaluate the likelihood on their sites

    @param tr 
      Tree instance

    @param lhs
      Likelihood array

    @param n
      Number of threads

    @param tid
      Thread id
 */ 
void perSiteLogLikelihoodsPthreads(pllInstance *tr, partitionList *pr, double *lhs, int n, int tid)
{
  size_t 
    model, 
    i;

  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    {      
      size_t 
	localIndex = 0;

      /* decide if this partition is handled by the thread when -Q is ativated 
	 or when -Q is not activated figure out which sites have been assigned to the 
	 current thread */

      pllBoolean 
	execute = ((tr->manyPartitions && isThisMyPartition(pr, tid, model)) || (!tr->manyPartitions));
      //On Windows, this will result in a compilation error if you don't have an MPI API 
      //downloaded and installed.  
      //I've only built with the Microsoft MPI (MS-MPI), v10.1.2 (November 2019).

      /* if the entire partition has been assigned to this thread (-Q) or if -Q is not activated 
	 we need to compute some per-site log likelihoods with thread tid for this partition */

      if(execute)
	for(i = (size_t)(pr->partitionData[model]->lower);  i < (size_t)(pr->partitionData[model]->upper); i++)
	  {
	    /* if -Q is active we compute all per-site log likelihoods for the partition,
	       othwerise we only compute those that have been assigned to thread tid 
	       using the cyclic distribution scheme */

	    if(tr->manyPartitions || (i % n == (size_t)tid))
	      {
	        double l = 0;

	        /* now compute the per-site log likelihood at the current site */

	        switch(tr->rateHetModel)
		        {
		        case PLL_CAT:
		        l = evaluatePartialGeneric (tr, pr, localIndex, pr->partitionData[model]->perSiteRates[pr->partitionData[model]->rateCategory[localIndex]], model);
		        break;
		        case PLL_GAMMA:
		        l = evaluatePartialGeneric (tr, pr, localIndex, 1.0, model);
		        break;
		        default:
		        assert(0);
		        }

	        /* store it in an array that is local in memory to the current thread,
		        see function collectDouble() in axml.c for understanding how we then collect these 
		        values stored in local arrays from the threads */

	        lhs[i] = l;

	        localIndex++;
	      }
	  }
    }
}

/** @brief Check if a partition is assign to a thread/process.

    Checks whether partition \a model from partition list \a localPr is
    assigned to be processed by process/thread with id \a tid.

    @param localTree
      Local PLL instance

    @param tid 
      Thread/Process id

    @param model
      Partition number
 */ 
pllBoolean isThisMyPartition(partitionList *localPr, int tid, int model)
{ 
  if(localPr->partitionData[model]->partitionAssignment == tid)
    return PLL_TRUE;
  else
    return PLL_FALSE;
}

/** @brief Computes partition size for all partitions (in case full partitions are assigns to workers). 

    @param localPr the local partitions instance
    
    @param tid thread id    
 */ 
static void computeFractionMany(partitionList *localPr, int tid)
{
  int
    sites = 0;

  int   
    model;

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      if(isThisMyPartition(localPr, tid, model))
	{	 
    	  localPr->partitionData[model]->width = localPr->partitionData[model]->upper - localPr->partitionData[model]->lower;
	  sites += localPr->partitionData[model]->width;
	}
      else       	  
    	  localPr->partitionData[model]->width = 0;
    }


}


/** @brief Computes partition size for all partitions (for cyclic distribution of sites)
    
    @param localPr the local partitions instance
    @param tid thread id
    @param n number of workers
 */ 
static void computeFraction(partitionList *localPr, int tid, int n)
{
  int
    i,
    model;

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      int width = 0;

      for(i = localPr->partitionData[model]->lower; i < localPr->partitionData[model]->upper; i++)
	if(i % n == tid)
	  width++;
      localPr->partitionData[model]->width = width;
    }
}



/** @brief Compare partition sizes. 
    @param p1 pointer to a partition
    @param p2 pointer to another partition
 */ 
static int partCompare(const void *p1, const void *p2)
{
  partitionType 
    *rc1 = (partitionType *)p1,
    *rc2 = (partitionType *)p2;

  int 
    i = rc1->partitionLength,
    j = rc2->partitionLength;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


/** @brief Top-level function for the multi processor scheduling
    scheme (assigns full partitions to workers).
    
   tr->manyPartitions is set to PLL_TRUE if the user has indicated via -Q
   that there are substantially more partitions than threads/cores
   available. In that case we do not distribute sites from each
   partition in a cyclic fashion to the cores , but distribute entire
   partitions to cores.  Achieving a good balance of alignment sites
   to cores boils down to the multi-processor scheduling problem known
   from theoretical comp. sci.  which is NP-complete.  We have
   implemented very simple "standard" heuristics for solving the
   multiprocessor scheduling problem that turn out to work very well
   and are cheap to compute.
   
   @param pr 
     List of partitions

   @param tid
     Id of current process/thread 
*/
static void multiprocessorScheduling(pllInstance * tr, partitionList *pr, int tid)
{
  int 
    s,
    model,
    modelStates[2] = {4, 20},
    numberOfPartitions[2] = {0 , 0},
      arrayLength = sizeof(modelStates) / sizeof(int);

      /* check that we have not addedd any new models for data types with a different number of states
	 and forgot to update modelStates */

      for(model = 0; model < pr->numberOfPartitions; model++)
	{        
	  pllBoolean 
	    exists = PLL_FALSE;

	  for(s = 0; s < arrayLength; s++)
	    {
	      exists = exists || (pr->partitionData[model]->states == modelStates[s]);
	      if(pr->partitionData[model]->states == modelStates[s])
		numberOfPartitions[s] += 1;
	    }

	  assert(exists);
	}

      for(s = 0; s < arrayLength; s++)
	{
	  if(numberOfPartitions[s] > 0)
	    {
	      size_t   
		checkSum = 0,
		sum = 0;

	      int    
		i,
		k,
#ifndef _FINE_GRAIN_MPI
		n = tr->numberOfThreads,
#else
		n = processes,
#endif
		p = numberOfPartitions[s],    
		*assignments = (int *)rax_calloc((size_t)n, sizeof(int));  

	      partitionType 
		*pt = (partitionType *)rax_malloc(sizeof(partitionType) * (size_t)p);



	      for(i = 0, k = 0; i < pr->numberOfPartitions; i++)
		{
		  if(pr->partitionData[i]->states == modelStates[s])
		    {
		      pt[k].partitionNumber = i;
		      pt[k].partitionLength = pr->partitionData[i]->upper - pr->partitionData[i]->lower;
		      sum += (size_t)pt[k].partitionLength;
		      k++;
		    }
		}

	      assert(k == p);

	      qsort(pt, p, sizeof(partitionType), partCompare);    

	      for(i = 0; i < p; i++)
		{
		  int 
		    k, 
		    min = INT_MAX,
		    minIndex = -1;

		  for(k = 0; k < n; k++)	
		    if(assignments[k] < min)
		      {
			min = assignments[k];
			minIndex = k;
		      }

		  assert(minIndex >= 0);

		  assignments[minIndex] +=  pt[i].partitionLength;
		  assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < pr->numberOfPartitions);
		  pr->partitionData[pt[i].partitionNumber]->partitionAssignment = minIndex;
		}

              
              /* Process i gets assignments[i] sites for modelStates[s] state model */

	      for(i = 0; i < n; i++)
		checkSum += (size_t)assignments[i];

	      assert(sum == checkSum);

	      rax_free(assignments);
	      rax_free(pt);
	    }
	}
}



/** @brief Reduce the first and second derivative of the likelihood
    function.
    
    We collect the first and second derivatives from the various
    threads and sum them up. It's similar to what we do in
    pllEvaluateGeneric() with the only difference that we have to collect
    two values (firsrt and second derivative) instead of onyly one (the
    log likelihood

   @warning operates on global reduction buffers \a globalResult
   
   @param tr tree 
   @param dlnLdlz first derivative
   @param d2lnLdlz2 second derivative
*/
void branchLength_parallelReduce(pllInstance *tr, double *dlnLdlz,  double *d2lnLdlz2, int numBranches )
{
#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS

  /* only the master executes this  */
  assert(tr->threadID == 0); 
  
  int b; 
  int t; 
  for(b = 0; b < numBranches; ++b)
    {
      dlnLdlz[b] = 0; 
      d2lnLdlz2[b] = 0; 

      for(t = 0; t < tr->numberOfThreads; ++t)
	{
	  dlnLdlz[b] += globalResult[t * numBranches * 2 + b ];
	  d2lnLdlz2[b] += globalResult[t * numBranches * 2 + numBranches + b];
	}
    }
#else 
  memcpy(dlnLdlz, globalResult, sizeof(double) * numBranches);
  memcpy(d2lnLdlz2, globalResult + numBranches, sizeof(double) * numBranches);
#endif
}



/** @brief Read from buffer or writes rates into buffer.  Return
    number of elems written.

    If \a read is set to \b PLL_TRUE, then the contents \a srcTar are
    copied to \a buf. Otherwise, the contents of \a buf are moved to
    \a srcTar.
   
   @param buf 
     Buffer

   @param srcTar 
     Pointer to either source or destination array

   @param tr
     PLL instance

   @param n number of workers

   @param tid process id

   @param read 
     If read-mode then set to \b PLL_TRUE

   @param countOnly
     if \b PLL_TRUE, simply return the number of elements
*/
static int doublesToBuffer(double *buf, double *srcTar, pllInstance *tr, partitionList *pr, int n, int tid, pllBoolean read, pllBoolean countOnly)
{
  int 
    model,
    i;
  double 
    *initPtr = buf; 

  for(model = 0; model < pr->numberOfPartitions; model++)
    {
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(pr, tid, model))
	    for(i = pr->partitionData[model]->lower; i < pr->partitionData[model]->upper; i++)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i]; 
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++;
	      }	  
	}      
      else
	{
	  for(i = pr->partitionData[model]->lower; i < pr->partitionData[model]->upper; i++)
	    if(i % n == tid)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i];
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++; 
	      }
	}
    }
  
  return buf - initPtr; 
}




/** @brief broadcast rates after rate optimization. 
    
    @param tre Library instance
    @param localTree local library instance 
    @param n number of workers 
    @param tid worker id 
    
    @todo mpi_alltoallv/w may be more efficient, but it is a hell to set up
 */ 
void broadcastAfterRateOpt(pllInstance *tr, pllInstance *localTree, partitionList *pr, int n, int tid)
{				  
  int
    num1 = 0,
    num2 = 0,
    num3 = 0, 
    i ; 
    
  for(i = 0; i < n; ++i)
    {
      double
	allBuf[tr->originalCrunchedLength * 3],
	buf1[tr->originalCrunchedLength],
	buf2[tr->originalCrunchedLength], 
	buf3[tr->originalCrunchedLength]; 

#ifdef _USE_PTHREADS
      if(i != tid)
	continue; 
#endif
      int numDouble = 0; 
      
      /* extract doubles  */

      num1 = doublesToBuffer(buf1, localTree->patrat, tr, pr, n,i, PLL_TRUE, i!= tid);
      num2 = doublesToBuffer(buf2, localTree->patratStored, tr, pr, n,i, PLL_TRUE, i!= tid);
      num3 = doublesToBuffer(buf3, localTree->lhs, tr, pr, n,i, PLL_TRUE, i!= tid);

      /* printf("%d + %d + %d\n", num1, num2, num3);  */

      numDouble += num1 + num2 + num3; 

      /* copy doubles  */
      
      memcpy(allBuf, buf1, num1 * sizeof(double)); 
      memcpy(allBuf + num1, buf2, num2 * sizeof(double)); 
      memcpy(allBuf + (num1 + num2) , buf3, num3 * sizeof(double)); 

      BCAST_BUF(allBuf, numDouble, MPI_DOUBLE, i); 

      memcpy(buf1, allBuf, num1 * sizeof(double)); 
      memcpy(buf2, allBuf + num1, num2 * sizeof(double)); 
      memcpy(buf3, allBuf + (num1 + num2), num3 * sizeof(double)); 
      
      /* re-insert doubles  */
      int assertCtr = 0; 
      assertCtr += doublesToBuffer(buf1, tr->patrat, tr, pr, n,i,PLL_FALSE, PLL_FALSE);
      assertCtr += doublesToBuffer(buf2, tr->patratStored, tr, pr, n,i,PLL_FALSE, PLL_FALSE);
      assertCtr += doublesToBuffer(buf3, tr->lhs, tr, pr, n,i,PLL_FALSE, PLL_FALSE);

      assert(assertCtr == numDouble); 
    }
}


/** @brief Collect doubles from workers to master.
 
    

    @param dst destination array
    @param src source array
    @param tr library instance 
    @param n number of workers 
    @param tid worker id 
 */
static void collectDouble(double *dst, double *src, pllInstance *tr, partitionList *pr, int n, int tid)
{
#ifdef _FINE_GRAIN_MPI    
  int
    assertNum = 0,
    i, 
    displacements[tr->numberOfThreads];
  double 
    buf[tr->originalCrunchedLength],
    resultBuf[tr->originalCrunchedLength]; 

  /* NOTE: This was moved here because it was an additional unnecessary move for the PTHREADS version. I didnt
  have time to check the MPI version, have to get back to this and remove it */
  /* gather own persite log likelihood values into local buffer  */
  int numberCollected = doublesToBuffer(buf, src, tr, pr,n,tid,PLL_TRUE, PLL_FALSE);

  /* this communicates all the values to the master */
  
  int numberPerWorker[tr->numberOfThreads];     
  if(MASTER_P)			/* master counts number to receive, receives and writes back */
    {
      for(i = 0; i < n; ++i)
	{
	  numberPerWorker[i] = doublesToBuffer(buf,src,tr,pr,n,i,PLL_FALSE, PLL_TRUE);
	  displacements[i] = i == 0 ? 0 : displacements[i-1] + numberPerWorker[i-1]; 
	}
      
      MPI_Gatherv(buf, numberCollected, MPI_DOUBLE,
		  resultBuf, numberPerWorker, displacements,  MPI_DOUBLE,
		  0, MPI_COMM_WORLD); 

      double *bufPtr = resultBuf; 
      for(i = 0 ; i < n; ++i)
	{
	  int numberWritten = doublesToBuffer(bufPtr, dst,tr,pr,n,i, PLL_FALSE, PLL_FALSE);
	  bufPtr += numberWritten; 
	  assertNum += numberWritten; 
	}    
      
      assert(assertNum == tr->originalCrunchedLength);
    }
  else 				/* workers only send their buffer   */
    MPI_Gatherv(buf, numberCollected, MPI_DOUBLE, resultBuf, numberPerWorker, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
#else 
  /* pthread version only writes to global space  */  

  //assertNum = doublesToBuffer(buf, dst,tr,pr,n,tid, PLL_FALSE, PLL_FALSE);
  doublesToBuffer (dst, src, tr, pr, n, tid, PLL_TRUE, PLL_FALSE);
  //assert(assertNum == numberCollected); 
#endif
}



/** @brief broadcast a new alpha (for the GAMMA model)
    @param localTree local library instance
    @param tr library instance
    @param tid worker id 
 */
static void broadCastAlpha(partitionList *localPr, partitionList *pr)
{
  int  i, 
    model; 

#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
    int bufSize = localPr->numberOfPartitions * 4 * sizeof(double);
  char bufDbl[bufSize]; 
  char *bufPtrDbl = bufDbl;   
#endif

  RECV_BUF(bufDbl, bufSize, MPI_BYTE); 

  for(model = 0; model < localPr->numberOfPartitions; model++)
    for(i = 0; i < 4; ++i)
      ASSIGN_BUF_DBL(localPr->partitionData[model]->gammaRates[i], pr->partitionData[model]->gammaRates[i]);
  
  SEND_BUF(bufDbl, bufSize, MPI_BYTE);  
}

/** @brief broadcast new LG4X weights
    @param localTree local library instance
    @param tr library instance
    @param tid worker id
 */
static void broadCastLg4xWeights(partitionList *localPr, partitionList *pr)
{
  int  i,
    model;

#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
    int bufSize = localPr->numberOfPartitions * 4 * sizeof(double);
  char bufDbl[bufSize];
  char *bufPtrDbl = bufDbl;
#endif

  RECV_BUF(bufDbl, bufSize, MPI_BYTE);

  for(model = 0; model < localPr->numberOfPartitions; model++)
    for(i = 0; i < 4; ++i)
      ASSIGN_BUF_DBL(localPr->partitionData[model]->lg4x_weights[i], pr->partitionData[model]->lg4x_weights[i]);

  SEND_BUF(bufDbl, bufSize, MPI_BYTE);
}

static void copyLG4(partitionList *localPr, partitionList *pr)
{
    int model, i, k;

    /* determine size of buffer needed first */
    int bufSize = 0;

#ifdef _FINE_GRAIN_MPI
    for(model = 0; model < localPr->numberOfPartitions; ++model )
      {
        const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]);
        bufSize += 4*(pl->eignLength + pl->evLength + pl->eiLength + pl->tipVectorLength + pl->substRatesLength + pl->frequenciesLength) * sizeof(double) ;
      }
#endif

    char bufDbl[bufSize];
    char *bufPtrDbl = bufDbl;

    RECV_BUF(bufDbl, bufSize, MPI_BYTE);

    for (model = 0; model < localPr->numberOfPartitions; model++)
    {
        pInfo * localInfo = localPr->partitionData[model];
        pInfo * info = pr->partitionData[model];

        if (info->protModels == PLL_LG4M || info->protModels == PLL_LG4X)
        {
            for (k = 0; k < 4; k++)
            {
                const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]);

                for (i = 0; i < pl->eignLength; ++i)
                    ASSIGN_BUF_DBL(
                            localPr->partitionData[model]->EIGN_LG4[k][i],
                            pr->partitionData[model]->EIGN_LG4[k][i]);
                for (i = 0; i < pl->evLength; ++i)
                    ASSIGN_BUF_DBL(localPr->partitionData[model]->EV_LG4[k][i],
                            pr->partitionData[model]->EV_LG4[k][i]);
                for (i = 0; i < pl->eiLength; ++i)
                    ASSIGN_BUF_DBL(localPr->partitionData[model]->EI_LG4[k][i],
                            pr->partitionData[model]->EI_LG4[k][i]);
                for (i = 0; i < pl->substRatesLength; ++i)
                    ASSIGN_BUF_DBL(
                            localPr->partitionData[model]->substRates_LG4[k][i],
                            pr->partitionData[model]->substRates_LG4[k][i]);
                for (i = 0; i < pl->frequenciesLength; ++i)
                    ASSIGN_BUF_DBL(
                            localPr->partitionData[model]->frequencies_LG4[k][i],
                            pr->partitionData[model]->frequencies_LG4[k][i]);
                for (i = 0; i < pl->tipVectorLength; ++i)
                    ASSIGN_BUF_DBL(
                            localPr->partitionData[model]->tipVector_LG4[k][i],
                            pr->partitionData[model]->tipVector_LG4[k][i]);
            }
        }
    }
    SEND_BUF(bufDbl, bufSize, MPI_BYTE); /*  */
}

/** @brief Master broadcasts rates.
    
    @param localTree local library instance
    @param tr library instance
    @param tid worker id     
 */ 
static void broadCastRates(partitionList *localPr, partitionList *pr)
{
  int 
    model;

  /* determine size of buffer needed first */
  int bufSize = 0;
#ifdef _FINE_GRAIN_MPI
  for(model = 0; model < localPr->numberOfPartitions; ++model )
    {	  
      const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]); /* this is constant, isnt it?  */
      bufSize += (pl->eignLength + pl->evLength + pl->eiLength + pl->tipVectorLength) * sizeof(double) ;
    }
#endif

  char
      bufDbl[bufSize];
    char *bufPtrDbl = bufDbl;

  RECV_BUF(bufDbl, bufSize, MPI_BYTE);
  int i ; 

  for(model = 0; model < localPr->numberOfPartitions; model++)
    {
      const partitionLengths *pl = getPartitionLengths(pr->partitionData[model]); /* this is constant, isnt it?  */

      for(i = 0; i < pl->eignLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EIGN[i], pr->partitionData[model]->EIGN[i]);
      for(i = 0; i < pl->evLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EV[i],pr->partitionData[model]->EV[i]);
      for(i = 0; i  < pl->eiLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->EI[i], pr->partitionData[model]->EI[i]);
      for(i = 0; i < pl->tipVectorLength; ++i)
	ASSIGN_BUF_DBL(localPr->partitionData[model]->tipVector[i],   pr->partitionData[model]->tipVector[i]);
    }
  SEND_BUF(bufDbl, bufSize, MPI_BYTE); /*  */

  copyLG4(localPr, pr);
}

/** @brief Evaluate the likelihood of this topology (PThreads/MPI implementation)

    Evaluate the likelihood of the topology described in the PLL instance. First
    every thread calls \a pllEvaluateIterative where it computes the log likelihoods
    for the  portion of each assigned partition. The results (for all partition) are stored
    as elements of a local buffer array (\a buf). This is done by all threads. Subsequently, 
    an \a MPI_Reduce operation sums the contents of corresponding elements of the local
    buffer arrays into another array (\a targetBuf) which are the log likelihoods of
    each (complete) partition. Finally, the last array is copied to the master thread/process.
    In addition, if \a getPerSiteLikelihoods is enabled the log likelihoods for each site
    in the (compressed) alignment are stored in the array \a tr->lhs.

    @param tr
      PLL instance
    @param tr
      Local (thread/process) PLL instance

    @param pr
      Local (thread/process) list of partitions

    @param tid
      Thread/Process ID

    @param getPerSiteLikelihoods 
      If set to \b PLL_TRUE, compute the log likelihood for each site. 
 */ 
static void reduceEvaluateIterative(pllInstance *tr, pllInstance *localTree, partitionList *localPr, int tid, pllBoolean getPerSiteLikelihoods)
{
  int model;

  pllEvaluateIterative(localTree, localPr, getPerSiteLikelihoods);

  /* when this is done we need to write the per-thread log likelihood to the 
     global reduction buffer. Tid is the thread ID, hence thread 0 will write its 
     results to reductionBuffer[0] thread 1 to reductionBuffer[1] etc.

     the actual sum over the entries in the reduction buffer will then be computed 
     by the master thread which ensures that the sum is determinsitic */

  
  /* if (getPerSiteLikelihoods == PLL_TRUE) store per-site likelihoods in array tr->lhs */
  if(getPerSiteLikelihoods)
    {    
#ifdef _FINE_GRAIN_MPI
      int n = processes; 
#else 
      int n = tr->numberOfThreads; 
#endif

      /* rearrange per site likelihoods into single local array for gathering */
      int i ; 
      for(model = 0; model < localPr->numberOfPartitions; ++model)
	{
	  pInfo *partition = localPr->partitionData[model]; 
	  pllBoolean isMyPartition  = isThisMyPartition(localPr, tid, model);

	  int ctr = 0; 
	  for(i = partition->lower; i < partition->upper; ++i)
	    {
	      if(tr->manyPartitions && isMyPartition)
		localTree->lhs[i] = partition->perSiteLikelihoods[ ctr++]; 
	      else if(NOT tr->manyPartitions && (i % n) == tid)
		localTree->lhs[i] = partition->perSiteLikelihoods[ctr++];
	    }
	}
      
      /* gather all the double into the global array */
      collectDouble(tr->lhs, localTree->lhs, localTree, localPr,  n, tid); 
    }

  /* printf("collecting done\n" ); */
#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
  /* 
     aberer: I implemented this as a mpi_gather operation into this buffer, 
     pthreads version emulates this gather; 
     master takes care of the reduction; 
  */

  double 
    buf[localPr->numberOfPartitions];

  for(model = 0; model < localPr->numberOfPartitions; ++model)
    buf[model] = localPr->partitionData[model]->partitionLH;

  /* either make reproducible or efficient */
  ASSIGN_GATHER(globalResult, buf, localPr->numberOfPartitions, PLL_DOUBLE, tid);

  /* printf("gather worked\n"); */
#else 
  /* the efficient mpi version: a proper reduce  */
  double 
    buf[localPr->numberOfPartitions];
  
  for(model = 0; model < localPr->numberOfPartitions; ++model)
    buf[model] = localPr->partitionData[model]->partitionLH;

  double 
    targetBuf[localPr->numberOfPartitions];
  
  memset(targetBuf, 0, sizeof(double) * localPr->numberOfPartitions);

  MPI_Reduce(buf, targetBuf, localPr->numberOfPartitions, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if(MASTER_P) 
    {
      for(model = 0; model < localPr->numberOfPartitions; ++model) {
	localPr->partitionData[model]->partitionLH = targetBuf[model];
      }
    }
#endif
}



/*@ @brief Broadcast the traversal descriptor to worker threads. 

  The one below is a hack we are re-assigning the local pointer to
  the global one the memcpy version below is just for testing and
  preparing the fine-grained MPI BlueGene version

  @param localTree local library instance
  @param tr library instance
*/
/* TODO: we should reset this at some point, the excplicit copy is just done for testing */
__inline static void broadcastTraversalInfo(pllInstance *localTree, pllInstance *tr, partitionList *localPr)
{
  /* @todo these two regions could be joined */
#ifdef _USE_PTHREADS
  /* memcpy -> memmove (see ticket #43). This function is sometimes called with localTree == tr,
   * in which case some memcpy implementations can corrupt the buffers.
   */
  
  localTree->td[0].functionType =            tr->td[0].functionType;
  localTree->td[0].count =                   tr->td[0].count ;
  localTree->td[0].traversalHasChanged =     tr->td[0].traversalHasChanged;

  memmove(localTree->td[0].executeModel,    tr->td[0].executeModel,    sizeof(pllBoolean) * localPr->numberOfPartitions);
  memmove(localTree->td[0].parameterValues, tr->td[0].parameterValues, sizeof(double) * localPr->numberOfPartitions);
  
  if(localTree->td[0].traversalHasChanged)
    memmove(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));

#else
  /* MPI */
  /* like in raxml-light: first we send a small message, if the
     travesalDescriptor is longer, then resend */
  
  int length = treeIsInitialized ? localPr->numberOfPartitions : 0;
  char broadCastBuffer[messageSize(length)]; 
  char *bufPtr = broadCastBuffer; 
  int i; 

  RECV_BUF(broadCastBuffer, messageSize(length), MPI_BYTE); 

  ASSIGN_BUF(localTree->td[0].functionType, tr->td[0].functionType , int);   
  ASSIGN_BUF(localTree->td[0].count,  tr->td[0].count , int); 
  ASSIGN_BUF(localTree->td[0].traversalHasChanged, tr->td[0].traversalHasChanged , int); 

  if(treeIsInitialized)  
    { 
      for(i = 0; i < localPr->numberOfPartitions; ++i)
	{
	  ASSIGN_BUF(localTree->td[0].executeModel[i],      tr->td[0].executeModel[i], int); 
	  ASSIGN_BUF(localTree->td[0].parameterValues[i],	 tr->td[0].parameterValues[i], double); 
	}      

      for(i = 0; i < TRAVERSAL_LENGTH; ++i )
	ASSIGN_BUF(localTree->td[0].ti[i], tr->td[0].ti[i], traversalInfo); 
    }
    
  SEND_BUF(broadCastBuffer, messageSize(length), MPI_BYTE); 

  /* now we send the second part of the traversal descriptor, if we
     exceed the pre-set number of elements */
  if(treeIsInitialized && localTree->td[0].count > TRAVERSAL_LENGTH) 
    {
      /* lets use the MPI_Datatype for this thing, what I've read it's
	 supposed to be more secure and efficient */
      MPI_Bcast(localTree->td[0].ti + TRAVERSAL_LENGTH, localTree->td[0].count - TRAVERSAL_LENGTH, TRAVERSAL_MPI, 0, MPI_COMM_WORLD );
    }
#endif
}


/** @brief helper that yields a string representation of a parallel region. 
    
    @param type type of parallel region
 */ 
const char* getJobName(int type)
{
  switch(type)  
    {
    case  PLL_THREAD_NEWVIEW:       
      return "PLL_THREAD_NEWVIEW";
    case PLL_THREAD_EVALUATE: 
      return "PLL_THREAD_EVALUATE";
    case PLL_THREAD_MAKENEWZ: 
      return "PLL_THREAD_MAKENEWZ";
    case PLL_THREAD_MAKENEWZ_FIRST: 
      return "PLL_THREAD_MAKENEWZ_FIRST";
    case PLL_THREAD_RATE_CATS: 
      return "PLL_THREAD_RATE_CATS";
    case PLL_THREAD_COPY_RATE_CATS: 
      return "PLL_THREAD_COPY_RATE_CATS";
    case PLL_THREAD_COPY_INIT_MODEL: 
      return "PLL_THREAD_COPY_INIT_MODEL";
    case PLL_THREAD_INIT_PARTITION: 
      return "PLL_THREAD_INIT_PARTITION";
    case PLL_THREAD_OPT_ALPHA: 
      return "PLL_THREAD_OPT_ALPHA";
    case PLL_THREAD_OPT_RATE: 
      return "PLL_THREAD_OPT_RATE";
    case PLL_THREAD_COPY_ALPHA: 
      return "PLL_THREAD_COPY_ALPHA";
    case PLL_THREAD_COPY_RATES: 
      return "PLL_THREAD_COPY_RATES";
    case PLL_THREAD_PER_SITE_LIKELIHOODS: 
      return "PLL_THREAD_PER_SITE_LIKELIHOODS";
    case PLL_THREAD_NEWVIEW_ANCESTRAL: 
      return "PLL_THREAD_NEWVIEW_ANCESTRAL";
    case PLL_THREAD_GATHER_ANCESTRAL: 
      return "PLL_THREAD_GATHER_ANCESTRAL";
    case PLL_THREAD_EXIT_GRACEFULLY: 
      return "PLL_THREAD_EXIT_GRACEFULLY";
    case PLL_THREAD_EVALUATE_PER_SITE_LIKES:
      return "PLL_THREAD_EVALUATE_PER_SITE_LIKES";
    default: assert(0); 
      return "Unrecognized Job Type";
    }
}

/**
   @brief Generic entry point for parallel regions (mostly broadcasts
   traversal descriptor first).

   This function here handles all parallel regions in the Pthreads
   version, when we enter this function pllMasterBarrier() has been called
   by the master thread from within the sequential part of the
   program, tr is the library instance (tree) at the master thread, 
   localTree is the library instance (tree) at the worker threads

   While this is not necessary, adress spaces of threads are indeed
   separated for easier transition to a distributed memory paradigm
   
   @param tr library instance
   @param localTree local library instance 
   @param tid worker id 
   @param n number of workers 
*/
static pllBoolean execFunction(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{
  int
    i,
    model,
    localCounter;

#ifdef MEASURE_TIME_PARALLEL
  double timeForParallelRegion = gettime();
#endif


#ifdef _USE_PTHREADS
  /* some stuff associated with the barrier implementation using Pthreads and busy wait */
  int currentJob = threadJob >> 16;
#endif

  /* here the master sends and all threads/processes receive the traversal descriptor */
  broadcastTraversalInfo(localTree, tr, localPr);

#ifdef _USE_PTHREADS
  /* make sure that nothing is going wrong */
  assert(currentJob == localTree->td[0].functionType);
#else   
  localTree = tr; 
  int currentJob = localTree->td[0].functionType; 
#endif

#ifdef DEBUG_PARALLEL
  printf("[%d] working on %s\n", tid, getJobName(currentJob)); 
#endif  

  switch(currentJob)
    { 
    case PLL_THREAD_NEWVIEW: 
      /* just a newview on the fraction of sites that have been assigned to this thread */

      pllNewviewIterative(localTree, localPr, 0);
      break;     
    case PLL_THREAD_EVALUATE: 
      reduceEvaluateIterative(tr, localTree, localPr, tid, PLL_FALSE);
      break;	
    case PLL_THREAD_MAKENEWZ_FIRST:

      /* this is the first call from within makenewz that requires getting the likelihood vectors to the left and 
         right of the branch via newview and doing some precomputations.
	 
         For details see comments in makenewzGenericSpecial.c 
      */
    case  PLL_THREAD_MAKENEWZ:
      {	
	double
	  dlnLdlz[PLL_NUM_BRANCHES],
	  d2lnLdlz2[PLL_NUM_BRANCHES]; 

	if(localTree->td[0].functionType == PLL_THREAD_MAKENEWZ_FIRST)
	  makenewzIterative(localTree, localPr);
	execCore(localTree, localPr, dlnLdlz, d2lnLdlz2);

	/* gather the first and second derivatives that have been written by each thread */
	/* as for evaluate above, the final sum over the derivatives will be computed by the 
	   master thread in its sequential part of the code */

	int numBranches = localPr->perGeneBranchLengths?localPr->numberOfPartitions:1;

#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
	/* MPI: implemented as a gather again, pthreads: just buffer copying */	
	double buf[ 2 * numBranches];
	memcpy( buf, dlnLdlz, numBranches * sizeof(double) );
	memcpy(buf + numBranches, d2lnLdlz2, numBranches * sizeof(double));

	ASSIGN_GATHER(globalResult, buf,  2 * numBranches, PLL_DOUBLE, tid);
#else 	
	double result[numBranches];
	memset(result,0, numBranches * sizeof(double));
	MPI_Reduce( dlnLdlz , result , numBranches, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MASTER_P)
	  memcpy(globalResult, result, sizeof(double) * numBranches);
	
	memset(result,0,numBranches * sizeof(double));
	MPI_Reduce( d2lnLdlz2 , result , numBranches, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MASTER_P)
	  memcpy(globalResult + numBranches, result, sizeof(double) * numBranches);
#endif
      }

      break;

    case PLL_THREAD_INIT_PARTITION:       

      /* broadcast data and initialize and allocate arrays in partitions */
      
      initializePartitionsMaster(tr, localTree, pr, localPr, tid, n);

      break;          
    case PLL_THREAD_COPY_ALPHA: 
    case PLL_THREAD_OPT_ALPHA:
      /* this is when we have changed the alpha parameter, inducing a change in the discrete gamma rate categories.
	 this is called when we are optimizing or sampling (in the Bayesioan case) alpha parameter values */
      
      /* distribute the new discrete gamma rates to the threads */
      broadCastAlpha(localPr,pr);

      /* compute the likelihood, note that this is always a full tree traversal ! */
      if(localTree->td[0].functionType == PLL_THREAD_OPT_ALPHA)
	reduceEvaluateIterative(tr, localTree, localPr, tid, PLL_FALSE);

      break;
    case PLL_THREAD_OPT_RATE:
    case PLL_THREAD_COPY_RATES:

      /* if we are optimizing the rates in the transition matrix Q this induces recomputing the eigenvector eigenvalue 
	 decomposition and the tipVector as well because of the special numerics in RAxML, the matrix of eigenvectors 
	 is "rotated" into the tip lookup table.

	 Hence if the sequential part of the program that steers the Q matrix rate optimization has changed a rate we
	 need to broadcast all eigenvectors, eigenvalues etc to each thread 
      */

      broadCastRates(localPr, pr);

      /* now evaluate the likelihood of the new Q matrix, this always requires a full tree traversal because the changes need
	 to be propagated throughout the entire tree */

      if(localTree->td[0].functionType == PLL_THREAD_OPT_RATE)
	reduceEvaluateIterative(tr, localTree, localPr, tid, PLL_FALSE);

      break;
    case PLL_THREAD_COPY_LG4X_RATES:

        broadCastLg4xWeights(localPr, pr);
        broadCastAlpha(localPr, pr);

        assert(localPr->partitionData[0]->lg4x_weights[0] == pr->partitionData[0]->lg4x_weights[0]);

        break;
    case PLL_THREAD_OPT_LG4X_RATE:

        broadCastLg4xWeights(localPr, pr);
        broadCastAlpha(localPr, pr);

        assert(localPr->partitionData[0]->lg4x_weights[0] == pr->partitionData[0]->lg4x_weights[0]);

        /* compute the likelihood, note that this is always a full tree traversal ! */
        reduceEvaluateIterative(tr, localTree, localPr, tid, PLL_FALSE);

        break;
    case PLL_THREAD_COPY_INIT_MODEL:
      {

	/* need to be very careful here ! PLL_THREAD_COPY_INIT_MODEL is also used when the program is restarted 
	   it is hence not sufficient to just initialize everything by the default values ! */

	broadCastRates(localPr, pr);
	broadCastAlpha(localPr, pr); /* isnt that only executed when we are on gamma?  */
	broadCastLg4xWeights(localPr, pr);

	/*
	  copy initial model parameters, the Q matrix and alpha are initially, when we start our likelihood search 
	  set to default values. 
	  Hence we need to copy all those values that are required for computing the likelihood 
	  with newview(), evaluate() and makenez() to the private memory of the threads 
	*/


	if( localTree->rateHetModel == PLL_CAT) /* TRICKY originally this should only be executed by workers  */
	  {
#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
	    int bufSize = 2 * localTree->originalCrunchedLength * sizeof(double); 
	    char bufDbl[bufSize], 
	      *bufPtrDbl = bufDbl; 
#endif

	    RECV_BUF(bufDbl, bufSize,MPI_BYTE); 

	    /* this should be local  */
	    for(model = 0; model < localPr->numberOfPartitions; model++)
	      localPr->partitionData[model]->numberOfCategories      = pr->partitionData[model]->numberOfCategories;


	    /* this is only relevant for the PSR model, we can worry about this later */
	    for(i = 0; i < localTree->originalCrunchedLength; ++i)
	      {
		ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]);
		ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	      }

	    SEND_BUF(bufDbl, bufSize, MPI_BYTE); 
	  }
      } 
      break;    
    case PLL_THREAD_RATE_CATS: 
      {
	/* this is for optimizing per-site rate categories under PSR, let's worry about this later */

	ASSIGN_DBL( localTree->lower_spacing,  tr->lower_spacing);
	ASSIGN_DBL( localTree->upper_spacing,  tr->upper_spacing);

	optRateCatPthreads(localTree, localPr, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

	broadcastAfterRateOpt(tr, localTree, localPr, n,  tid);
      }
      break;
    case PLL_THREAD_COPY_RATE_CATS:
      {
	/* 
	   this is invoked when we have changed the per-site rate category assignment
	   In essence it distributes the new per site rates to all threads 

	   The pthread-version here simply assigns everything as ought to
	   be. The MPI-version is configured to write to a buffer instead
	   and SEND (master) or RECV (workers) it.

	*/

	/* 
	   start of communication part 
	*/

	int i, 
	  /* buf[localPr->numberOfPartitions], */
	  /* assertCtr = 0,  */
	  dblBufSize = 0; 

#if defined(FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
	int bufSize = localPr->numberOfPartitions * sizeof(int); 
	char buf[bufSize]; 
	char *bufPtr = buf; 
#endif
     
	RECV_BUF(buf, bufSize, MPI_BYTE);

	for( model = 0; model < localPr->numberOfPartitions; ++model)
	  {
	    ASSIGN_BUF(localPr->partitionData[model]->numberOfCategories, pr->partitionData[model]->numberOfCategories, int);
	    dblBufSize += localPr->partitionData[model]->numberOfCategories * sizeof(double);
	  }

	SEND_BUF(buf, bufSize, MPI_BYTE); 


	dblBufSize += 2 * localTree->originalCrunchedLength * sizeof(double); 

#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
	char bufDbl[dblBufSize],
	  *bufPtrDbl = bufDbl;
#endif

	RECV_BUF(bufDbl, dblBufSize, MPI_BYTE); 

	for(i = 0; i < localTree->originalCrunchedLength; ++i)
	  {	 
	    ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]); 
	    ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	  }

	for( model = 0; model < localPr->numberOfPartitions; ++model)
	  for(i = 0; i < localPr->partitionData[model]->numberOfCategories; i++)
	    ASSIGN_BUF_DBL(localPr->partitionData[model]->perSiteRates[i], pr->partitionData[model]->perSiteRates[i]);

	SEND_BUF(bufDbl, dblBufSize, MPI_BYTE); 


	/* lets test, if it is a good idea to send around the basic categories  */
#ifdef _FINE_GRAIN_MPI
	/* TODO this is inefficient, but is seems to have a small impact on performance */
	MPI_Bcast(tr->rateCategory, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD); 
#endif


	/* 
	   now re-assign values 
	*/
	for(model = 0; model < localPr->numberOfPartitions; model++)
	  {
	    if(localTree->manyPartitions)
	      {
		if(isThisMyPartition(localPr, tid, model))
		  for(localCounter = 0, i = localPr->partitionData[model]->lower;  i < localPr->partitionData[model]->upper; i++, localCounter++)
		    {	     
		      localPr->partitionData[model]->rateCategory[localCounter] = tr->rateCategory[i];
		    } 
	      }
	    else	  
	      {
		for(localCounter = 0, i = localPr->partitionData[model]->lower;  i < localPr->partitionData[model]->upper; i++)
		  {
		    if(i % n == tid)
		      {		 
			localPr->partitionData[model]->rateCategory[localCounter] = tr->rateCategory[i];

			localCounter++;
		      }
		  }
	      }
	  }
      }
      break;
    case PLL_THREAD_PER_SITE_LIKELIHOODS:      
      {

	/* compute per-site log likelihoods for the sites/partitions 
	   that are handled by this thread */
	perSiteLogLikelihoodsPthreads(localTree, localPr, localTree->lhs, n, tid);

	/* do a parallel gather operation, the threads will write their results 
	   into the global buffer tr->lhs that will then contain all per-site log likelihoods
	   in the proper order 
	*/

	collectDouble(tr->lhs,                localTree->lhs,                  localTree, localPr, n, tid);

      }
      break;
      /* check for errors */
    case PLL_THREAD_NEWVIEW_ANCESTRAL:       
      assert(0);
      break; 
    case PLL_THREAD_GATHER_ANCESTRAL:
      assert(0); 
      break; 
    case PLL_THREAD_EXIT_GRACEFULLY: 
      {
	/* cleans up the workers memory */

#ifdef _USE_PTHREADS
	/* TODO destroying the tree does not work yet in a highly
	   generic manner. */

	if(NOT MASTER_P)
	  {
	    pllPartitionsDestroy (localTree, &localPr);
	    /* pllTreeDestroy (localTree); */
	  }
	else 
	  {
	    //pllPartitionsDestroy (tr, &pr);
	    /* pllTreeDestroy (tr); */
	  }

#else 
	//pllPartitionsDestroy (tr, &pr);
	/* pllTreeDestroy (tr); */
	
	//MPI_Finalize();
	//exit(0); 
#endif	
	return PLL_FALSE; 
      }
      break; 
    case PLL_THREAD_EVALUATE_PER_SITE_LIKES: 
      {
	reduceEvaluateIterative(tr, localTree, localPr, tid, PLL_TRUE);
      }
      break;
    default:
      printf("Job %d\n", currentJob);
      assert(0);
    }

  return PLL_TRUE; 
}




/**  Target function where the threads/processes are trapped

     The threads/processes spend all of their time in this function
     running operations on the data (computing likelihoods).

     @param tData
       Structure that contains the vital information for the thread/process, 
       i.e. PLL instance, list of partitions and thread ID

     @note
       The data in \a tData are different for pthreads and MPI. 
       Expand this section.
 */ 
static void *likelihoodThread(void *tData)
{
  threadData *td = (threadData*)tData;
  pllInstance 
    *tr = td->tr;
  partitionList *pr = td->pr;

#ifdef _USE_PTHREADS
  pllInstance *localTree = rax_calloc(1,sizeof(pllInstance )); 
  partitionList *localPr = rax_calloc(1,sizeof(partitionList));

  int
    myCycle = 0,
    localTrap = 1;

  const int 
    n = td->tr->numberOfThreads,
    tid = td->threadNumber;

#ifndef _PORTABLE_PTHREADS
  pinToCore(tid);
#endif

  /* printf("\nThis is RAxML Worker Pthread Number: %d\n", tid); */

  while(localTrap)
    {

      while (myCycle == threadJob);
      myCycle = threadJob;

      if ((threadJob >> 16) != PLL_THREAD_INIT_PARTITION) {
    	  localPr->perGeneBranchLengths = pr->perGeneBranchLengths;
      	  localPr->numberOfPartitions = pr->numberOfPartitions;
      }
      localTrap = execFunction(tr, localTree, pr, localPr, tid, n);

      barrierBuffer[tid] = 1;     
    }
    rax_free (localTree->td[0].executeModel); //localTree->td[0].executeModel = NULL;
    rax_free (localTree->td[0].parameterValues); //localTree->td[0].parameterValues = NULL;
    rax_free (localTree->rateCategory); //localTree->rateCategory = NULL;
    rax_free (localTree->lhs); //localTree->lhs = NULL;
    rax_free (localTree->patrat); //localTree->patrat = NULL;
    rax_free (localTree->patratStored); //localTree->patratStored = NULL;
    rax_free (localTree->td[0].ti); //localTree->td[0].ti = NULL;
    rax_free (localTree);
#else 
  const int
    n = processes, 
    tid = td->threadNumber;
  int i;

  /* printf("\nThis is RAxML Worker Process Number: %d\n", tid); */

  while(execFunction(tr, tr, pr, pr, tid,n));

  rax_free (tr->lhs);
  rax_free (tr->td[0].ti);
  rax_free (tr->td[0].executeModel);
  rax_free (tr->td[0].parameterValues);
  rax_free (tr->patrat);
  rax_free (tr->patratStored);
  rax_free (tr->aliaswgt);
  rax_free (tr->y_ptr);
  for (i = 0; i < pr->numberOfPartitions; ++ i)
    rax_free (pr->partitionData[i]);
  rax_free (pr->partitionData);
  rax_free (pr);
  rax_free (tr);
#endif

  return (void*)NULL;
}


/**
   @brief Cleanup step once the master barrier succeeded. 

   This is master specific code called once the barrier is
   passed. Stuff such as reduction operations.  If we execute this
   here, we can keep the code mostly free from parallel -specific
   code.
   
   @param tr 
     PLL instance

   @param pr
     List of partitions

   @param jobType 
     Job that is to be executed
*/
void pllMasterPostBarrier(pllInstance *tr, partitionList *pr, int jobType)
{
  assert(tr->threadID == 0); 
  
  switch(jobType)
    {
    case PLL_THREAD_EVALUATE: 
    case PLL_THREAD_OPT_RATE: 
    case PLL_THREAD_OPT_ALPHA:
    case PLL_THREAD_OPT_LG4X_RATE:
    case PLL_THREAD_EVALUATE_PER_SITE_LIKES: 
      {
#ifdef _REPRODUCIBLE_MPI_OR_PTHREADS
	int i,j;
	volatile double partitionResult;	

	for(j = 0; j < pr->numberOfPartitions; j++)
	  {
	    for(i = 0, partitionResult = 0.0; i < tr->numberOfThreads; i++) 
	      partitionResult += globalResult[i * pr->numberOfPartitions+ j];

	    pr->partitionData[j]->partitionLH = partitionResult;
	  }
#endif      

	break; 
      } 
    case PLL_THREAD_PER_SITE_LIKELIHOODS:
      {
	int i; 
	/* now just compute the sum over per-site log likelihoods for error checking */      
	double accumulatedPerSiteLikelihood = 0.; 
	for(i = 0; i < tr->originalCrunchedLength; i++)
	  accumulatedPerSiteLikelihood += tr->lhs[i];

	/* printf("RESULT: %f\t%f", tr->likelihood, accumulatedPerSiteLikelihood);  */
	assert(PLL_ABS(tr->likelihood - accumulatedPerSiteLikelihood) < 0.00001);
      }
      break;
    default: 
      ; 			/* dont do anything on default,
				   mostly, we can skip that */
    } 
}

/**
   @brief A generic master barrier for executing parallel parts of the code

   A generic master barrier through which the master thread/process controls
   the work job execution. Through the parameter \a jobType the master instructs
   the slaves of what type of work they must conduct.

   @param tr
     PLL instance

   @param pr
     List of partitions

   @param jobType 
     Type of job to be conducted
 */ 
void pllMasterBarrier(pllInstance *tr, partitionList *pr, int jobType)
{

#ifdef MEASURE_TIME_PARALLEL
  assert(jobType < NUM_PAR_JOBS); 
  timePerRegion[NUM_PAR_JOBS]  += gettime()- masterTimePerPhase ; 
  masterTimePerPhase = gettime();
#endif

#ifdef _USE_PTHREADS
  const int 
    n = tr->numberOfThreads;

  tr->td[0].functionType = jobType;

  jobCycle = !jobCycle;
  threadJob = (jobType << 16) + jobCycle;

  execFunction(tr, tr, pr, pr, 0, n);

  int 
    i, 
    sum;

  do
    {
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);  

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
#else 
  tr->td[0].functionType = jobType; 
  execFunction(tr,tr,pr,pr,0,processes);
#endif

  /* code executed by the master, once the barrier is crossed */
  pllMasterPostBarrier(tr, pr, jobType);

#ifdef MEASURE_TIME_PARALLEL
  timePerRegion[jobType] += gettime() - masterTimePerPhase; 
  masterTimePerPhase = gettime();
#endif
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

/** @brief Initialize structures for slave process/threads
 
    Allocate all memory structures required by slave threads/processes

    @param tr 
      PLL Instance

    @param localTree 
      A local PLL instance for the slave process/thread which is initialized in this function based on \a tr

    @pram pr
      List of partitions

    @param localPr
      A local list of partitions for the slave process/thread which will be initialized based on \a pr 

    @pram tid
      The slave process/thread ID

    @note
      This function should never be called by the master thread, but is called by master process in MPI implementation.
 */ 
static void assignAndInitPart1(pllInstance *localTree, pllInstance *tr, partitionList *localPr, partitionList *pr, int *tid)
{
  size_t
    model; 
  int
    totalLength = 0; 

#ifdef _USE_PTHREADS
  localTree->threadID = *tid; 
  /* printf("my id is %d\n", *tid);  */
  assert(localTree != tr);
  localTree->numberOfThreads = tr->numberOfThreads;
#else  /* => MPI */
  *tid = processID; 
  localTree->threadID = processID; 
  tr->numberOfThreads = processes;

  int bufSize = (9 + pr->numberOfPartitions* 8) * sizeof(int);
  char buf[bufSize], 
    *bufPtr = buf;  
#endif

  RECV_BUF(buf, bufSize, MPI_BYTE); 

  ASSIGN_BUF( localTree->useRecom,                  tr->useRecom, int);
  ASSIGN_BUF( localTree->rateHetModel,              tr->rateHetModel, int);
  ASSIGN_BUF( localTree->useMedian,                 tr->useMedian, int); 
  ASSIGN_BUF( localTree->saveMemory,                tr->saveMemory, int);
  ASSIGN_BUF( localTree->maxCategories,             tr->maxCategories, int);
  ASSIGN_BUF( localTree->originalCrunchedLength,    tr->originalCrunchedLength, int);
  ASSIGN_BUF( localTree->mxtips,                    tr->mxtips, int);
  ASSIGN_BUF( localPr->numberOfPartitions,          pr->numberOfPartitions, int);
  ASSIGN_BUF( localPr->perGeneBranchLengths,        pr->perGeneBranchLengths, pllBoolean);

  localTree->td[0].count = 0; 

  if(NOT MASTER_P)
    {
      localTree->lhs                     = (double*)rax_calloc((size_t)localTree->originalCrunchedLength, sizeof(double));     
      localPr->partitionData           = (pInfo**)rax_calloc(PLL_NUM_BRANCHES,sizeof(pInfo*));
      for(model = 0; model < (size_t)localPr->numberOfPartitions; model++) {
    	localPr->partitionData[model] = (pInfo*)rax_calloc(1,sizeof(pInfo));
      }
      localTree->td[0].ti              = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)localTree->mxtips);
      localTree->td[0].executeModel    = (pllBoolean *)rax_malloc(sizeof(pllBoolean) * PLL_NUM_BRANCHES);
      localTree->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * PLL_NUM_BRANCHES);
      localTree->patrat       = (double*)rax_malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);
      localTree->patratStored = (double*)rax_malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);            
    }
  
  for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
    {
      ASSIGN_BUF(localPr->partitionData[model]->numberOfCategories,     pr->partitionData[model]->numberOfCategories, int);
      ASSIGN_BUF(localPr->partitionData[model]->states,                 pr->partitionData[model]->states, int);
      ASSIGN_BUF(localPr->partitionData[model]->maxTipStates ,          pr->partitionData[model]->maxTipStates, int);
      ASSIGN_BUF(localPr->partitionData[model]->dataType ,              pr->partitionData[model]->dataType, int);
      ASSIGN_BUF(localPr->partitionData[model]->protModels ,            pr->partitionData[model]->protModels, int);
      ASSIGN_BUF(localPr->partitionData[model]->protUseEmpiricalFreqs , pr->partitionData[model]->protUseEmpiricalFreqs, int);
      ASSIGN_BUF(localPr->partitionData[model]->lower ,                 pr->partitionData[model]->lower, int);
      ASSIGN_BUF(localPr->partitionData[model]->upper ,                 pr->partitionData[model]->upper, int);
      ASSIGN_BUF(localPr->partitionData[model]->ascBias,                pr->partitionData[model]->ascBias, pllBoolean);

      localPr->partitionData[model]->partitionLH = 0.0;      

      totalLength += (localPr->partitionData[model]->upper -  localPr->partitionData[model]->lower);
    }

  SEND_BUF(buf, bufSize, MPI_BYTE); 

  assert(totalLength == localTree->originalCrunchedLength);

  ASSIGN_DBL(localTree->vectorRecomFraction, tr->vectorRecomFraction); 
}
#endif


/** @brief Distribute y-vectors during initialization. 

    Distribute the alignment data to the slave process/threads. Each slave
    copies the data (alignment) from its assigned partition to its local 
    partition structure.

    @param tr 
      PLL instance
    
    @param localTree 
      Local library instance for the current thread

    @param localPr
      Local list of partitions structure for the current thread
 */ 
static void distributeYVectors(pllInstance *localTree, pllInstance *tr, partitionList *localPr)
{
  size_t 
    i,
    n = localTree->numberOfThreads,
    globalCounter = 0,
    localCounter = 0,
    model = 0, 
    j; 
  int tid = localTree->threadID; 
  

  /* distribute the y-vectors */
  for(j = 1 ; j <= (size_t)localTree->mxtips; j++)	
    {
#if defined(_FINE_GRAIN_MPI) || defined(_IQTREE_MPI)
      unsigned char yBuf[tr->originalCrunchedLength]; 	  
      if (MASTER_P) {
          memcpy(yBuf, tr->yVector[j], tr->originalCrunchedLength * sizeof(unsigned char));
      }
      MPI_Bcast(  yBuf, tr->originalCrunchedLength, MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD); 
#endif	  

      for(model = 0, globalCounter = 0; model < (size_t)localPr->numberOfPartitions; model++)
	{
	  if(tr->manyPartitions)
	    {
	      if(isThisMyPartition(localPr, tid, model))
		{
		  assert(localPr->partitionData[model]->upper - localPr->partitionData[model]->lower == localPr->partitionData[model]->width);
		  for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, localCounter++, globalCounter++)
#ifdef _USE_PTHREADS
		    localPr->partitionData[model]->yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		  localPr->partitionData[model]->yVector[j][localCounter] = yBuf[globalCounter];
#endif


		}
	      else
		globalCounter += (localPr->partitionData[model]->upper - localPr->partitionData[model]->lower);
	    }
	  else 
	    {
	      for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, globalCounter++)
		{
		  if(i % (size_t)n == (size_t)tid)
		    {
#ifdef _USE_PTHREADS
		      localPr->partitionData[model]->yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		      localPr->partitionData[model]->yVector[j][localCounter] = yBuf[globalCounter];
#endif
		      ++localCounter; 
		    }
		}	   
	    }
	}
    }
}

/** @brief Distribute the weights in the alignment of slave process/threads

    Allocate space in the local tree structure for the alignment weights. Then
    copy the weights vector from the master process/thread to the slaves.

    @param tr 
      PLL instance
    
    @param localTree 
      Local library instance for the current process/thread

    @param localPr
      Local list of partitions for the current process/thread

    @todo
      The alignment weights should go to the partitions structure rather than the tree structure
 */ 
static void distributeWeights(pllInstance *localTree, pllInstance *tr, partitionList *localPr)
{
  int tid = localTree->threadID; 
  int n = localTree->numberOfThreads; 

  size_t     
    globalCounter = 0,
    i,
    localCounter  = 0,
    model; 



  /* distribute the weights  */
#ifdef _FINE_GRAIN_MPI 		/* need to broadcast a few things first */
  if(NOT MASTER_P)
    tr->aliaswgt = rax_malloc(sizeof(int) * tr->originalCrunchedLength); 
  MPI_Bcast(tr->aliaswgt, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD);      
#endif
  for(model = 0, globalCounter = 0; model < (size_t)localPr->numberOfPartitions; model++)
    { 
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(localPr, tid, model))
	    {
	      assert(localPr->partitionData[model]->upper - localPr->partitionData[model]->lower == localPr->partitionData[model]->width);
	      for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, localCounter++, globalCounter++)
		localPr->partitionData[model]->wgt[localCounter]          = tr->aliaswgt[globalCounter];
	    }
	  else
	    globalCounter += (localPr->partitionData[model]->upper - localPr->partitionData[model]->lower);
	}
      else 
	{ 
	  for(localCounter = 0, i = (size_t)localPr->partitionData[model]->lower;  i < (size_t)localPr->partitionData[model]->upper; i++, globalCounter++)
	    {
	      if(i % (size_t)n == (size_t)tid)
		localPr->partitionData[model]->wgt[localCounter++]       = tr->aliaswgt[globalCounter];
	    }	   
	}
    }
}


/** @brief Initialize the partitioning scheme (master function) in parallel environment.
    
    Initialize the partition scheme in all processes/threads. This is a wrapper function
    that calls all necessary functions for allocating the local structures for slave threads
    and for distributing all necessary data from the master threads, such as alignment data,
    and weight vectors.

    @param tr 
      PLL instance

    @param localTree 
      Local PLL instance for the slave process/thread

    @param pr
      List of partitions

    @param localPr
      Local partition structure for the slave process/thread

    @param tid
      Process/thread id

    @param n 
      Number of processes/threads
*/ 
static void initializePartitionsMaster(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{ 
  size_t
    model;

  treeIsInitialized = PLL_TRUE; 

  ASSIGN_INT(localTree->manyPartitions, tr->manyPartitions);
  ASSIGN_INT(localTree->numberOfThreads, tr->numberOfThreads);  
  ASSIGN_INT(localPr->numberOfPartitions, pr->numberOfPartitions);

#ifdef _USE_PTHREADS
  if(MASTER_P)
    globalResult = rax_calloc((size_t) tr->numberOfThreads * (size_t)pr->numberOfPartitions* 2 ,sizeof(double));
  else 
    assignAndInitPart1(localTree, tr, localPr, pr, &tid);
#else 
  globalResult = (double*)rax_calloc((size_t) tr->numberOfThreads * (size_t)pr->numberOfPartitions* 2 ,sizeof(double));
  assignAndInitPart1(localTree, tr, localPr, pr, &tid);
  defineTraversalInfoMPI();
#endif

  for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
    localPr->partitionData[model]->width        = 0;

  if(tr->manyPartitions)    
    {
      multiprocessorScheduling(localTree, localPr, tid);
      computeFractionMany(localPr, tid);
    }
  else
    computeFraction(localPr, tid, n);

  initializePartitionData(localTree, localPr);

  {
    size_t 
      model,  
      i,      
      countOffset,
      myLength = 0;

    for(model = 0; model < (size_t)localPr->numberOfPartitions; model++)
      myLength += localPr->partitionData[model]->width;

    /* assign local memory for storing sequence data */
    
    localTree->y_ptr = (unsigned char *)rax_malloc(myLength * (size_t)(localTree->mxtips) * sizeof(unsigned char));
    assert(localTree->y_ptr != NULL);

    for(i = 0; i < (size_t)localTree->mxtips; i++)
      {
	for(model = 0, countOffset = 0; model < (size_t)localPr->numberOfPartitions; model++)
	  {	    
	    localPr->partitionData[model]->yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
	    countOffset +=  localPr->partitionData[model]->width;
	  }
	assert(countOffset == myLength);
      }

    /* figure in data */

    distributeWeights(localTree, tr, localPr);

    distributeYVectors(localTree, tr, localPr);

  }

  initMemorySavingAndRecom(localTree, localPr);
}
