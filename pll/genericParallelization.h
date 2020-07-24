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
 * @file genericParallelization.h
 */
#ifndef _GENERIC_PARALL_H 
#define _GENERIC_PARALL_H 


extern double *globalResult; 


/**********/
/* CONFIG */
/**********/

/* #define MEASURE_TIME_PARALLEL */
#define _PORTABLE_PTHREADS
/* #define DEBUG_PARALLEL */ 
/* #define DEBUG_MPI_EACH_SEND */
/* #define _REPRODUCIBLE_MPI_OR_PTHREADS */
#ifdef _USE_PTHREADS
#ifndef _PORTABLE_PTHREADS
void pinToCore(int tid);
#endif
#endif


#define NOT ! 
#define IS_PARALLEL (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI)) 



#ifdef MEASURE_TIME_PARALLEL
#define NUM_PAR_JOBS 16
extern double masterTimePerPhase; 
#endif


/******************/
/* MPI SPECIFIC   */
/******************/
#if defined(_FINE_GRAIN_MPI) || (defined(_IQTREE_MPI) && defined(CLANG_UNDER_VS))
#include <mpi.h>
#ifdef DEBUG_MPI_EACH_SEND
#define DEBUG_PRINT(text, elem) printf(text, elem)
#else 
#define DEBUG_PRINT(text, elem) NULL
#endif

/* for the broadcast of traversal descriptor */
#define TRAVERSAL_LENGTH 5
#define traversalSize sizeof(traversalInfo)
#define messageSize(x)   (3 * sizeof(int) +  x * (sizeof(int)+ sizeof(double)) + TRAVERSAL_LENGTH * traversalSize)

#define VOLATILE_PAR 
#define MASTER_P (processID == 0)
#define POP_OR_PUT_BYTES(bufPtr, elem, type) (MASTER_P ? (bufPtr = addBytes((bufPtr), &(elem), sizeof(type))) : (bufPtr = popBytes((bufPtr), &(elem), sizeof(type))))

#define ASSIGN_INT(x,y) (MPI_Bcast((int*)(&y),1,MPI_INT,0,MPI_COMM_WORLD),DEBUG_PRINT("\tSEND/RECV %d\n", y)) 
#define ASSIGN_BUF(x,y,type) (POP_OR_PUT_BYTES(bufPtr, y,type))
#define ASSIGN_BUF_DBL(x,y) (POP_OR_PUT_BYTES(bufPtrDbl,y, double))
#define ASSIGN_DBL(x,y) (MPI_Bcast(&y,1,MPI_DOUBLE, 0, MPI_COMM_WORLD), DEBUG_PRINT("\tSEND/RECV %f\n", y)) 
#define ASSIGN_DBLS(tar,src,length) MPI_Bcast(tar, length, MPI_DOUBLE, 0, MPI_COMM_WORLD)
#define PLL_DOUBLE MPI_DOUBLE
#define ASSIGN_GATHER(tar,src,length,type,tid) MPI_Gather(src,length,type,tar,length,type,0, MPI_COMM_WORLD)
#define SEND_BUF(buf, bufSize,type) if(MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define RECV_BUF(buf, bufSize,type) if(NOT MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define BCAST_BUF(buf, bufSize,type,who)  MPI_Bcast(buf, bufSize, type, who,MPI_COMM_WORLD )



extern int processes; 
extern int processID; 
#endif 

/*********************/
/* PTHREAD SPECIFIC  */
/*********************/
#ifdef _USE_PTHREADS
#ifndef CLANG_UNDER_VS
#if defined (_MSC_VER) 
#include "pthread.h"
#else
#include <pthread.h>
#endif
#endif
#define _REPRODUCIBLE_MPI_OR_PTHREADS
#define VOLATILE_PAR volatile 
#define MASTER_P (tid == 0)
#define ASSIGN_INT(x,y) (x = y)
#define ASSIGN_BUF(x,y,type) (x = y)
#define ASSIGN_BUF_DBL(x,y) (x = y)
#define ASSIGN_DBL(x,y) (x = y)
#define ASSIGN_DBLS(tar,src,length) memmove(tar, src, length * sizeof(double))
#define PLL_DOUBLE double 	/* just rededining that to make the source code less confusing */
#define ASSIGN_GATHER(tar,src,length,type,tid) (memmove((tar) + (tid) * (length) ,src, length * sizeof(type)))
#define SEND_BUF(buf, bufSize, type) 
#define RECV_BUF(buf, bufSize, type) 
#define BCAST_BUF(buf, bufSize,type,who)  
#define TRAVERSAL_LENGTH 5
#define messageSize(x) 0
#endif


#endif	/* end include guard  */
