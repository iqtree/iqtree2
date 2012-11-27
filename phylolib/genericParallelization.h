#ifndef _GENERIC_PARALL_H 
#define _GENERIC_PARALL_H 


extern double *globalResult; 


/**********/
/* CONFIG */
/**********/

/* #define MEASURE_TIME_PARALLEL */
#define _PORTABLE_PTHREADS
/* #define DEBUG_PARALLEL  */
/* #define DEBUG_MPI_EACH_SEND */
/* #define _REPRODUCIBLE_MPI_OR_PTHREADS */



#define NOT ! 
#define IS_PARALLEL (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI)) 
void *likelihoodThread(void *tData); 



#ifdef MEASURE_TIME_PARALLEL
#define NUM_PAR_JOBS 16
extern double masterTimePerPhase; 
#endif


/******************/
/* MPI SPECIFIC   */
/******************/
#ifdef _FINE_GRAIN_MPI
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

#define ASSIGN_INT(x,y) (MPI_Bcast(&y,1,MPI_INT,0,MPI_COMM_WORLD),DEBUG_PRINT("\tSEND/RECV %d\n", y)) 
#define ASSIGN_BUF(x,y,type) (POP_OR_PUT_BYTES(bufPtr, y,type))
#define ASSIGN_BUF_DBL(x,y) (POP_OR_PUT_BYTES(bufPtrDbl,y, double))
#define ASSIGN_DBL(x,y) (MPI_Bcast(&y,1,MPI_DOUBLE, 0, MPI_COMM_WORLD), DEBUG_PRINT("\tSEND/RECV %f\n", y)) 
#define ASSIGN_DBLS(tar,src,length) MPI_Bcast(tar, length, MPI_DOUBLE, 0, MPI_COMM_WORLD)
#define DOUBLE MPI_DOUBLE
#define ASSIGN_GATHER(tar,src,length,type,tid) MPI_Gather(src,length,type,tar,length,type,0, MPI_COMM_WORLD)
#define SEND_BUF(buf, bufSize,type) if(MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define RECV_BUF(buf, bufSize,type) if(NOT MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define BCAST_BUF(buf, bufSize,type,who)  MPI_Bcast(buf, bufSize, type, who,MPI_COMM_WORLD )



#define INIT_LENGTH
typedef struct  _jobDef
{
  int functionType; 
  int count; 
  int traveralHasChanged;
  
  
} InitjobDescr; 

extern int processes; 
extern int processID; 
char* addBytes(char *buf, void *toAdd, int numBytes); 
char* popBytes(char *buf, void *result, int numBytes); 
#endif 

/*********************/
/* PTHREAD SPECIFIC  */
/*********************/
#ifdef _USE_PTHREADS
#include <pthread.h>
#define _REPRODUCIBLE_MPI_OR_PTHREADS
#define VOLATILE_PAR volatile 
#define MASTER_P (tid == 0)
#define ASSIGN_INT(x,y) (x = y)
#define ASSIGN_BUF(x,y,type) (x = y)
#define ASSIGN_BUF_DBL(x,y) (x = y)
#define ASSIGN_DBL(x,y) (x = y)
#define ASSIGN_DBLS(tar,src,length) memmove(tar, src, length * sizeof(double))
#define DOUBLE double 	/* just rededining that to make the source code less confusing */
#define ASSIGN_GATHER(tar,src,length,type,tid) (memmove((tar) + (tid) * (length) ,src, length * sizeof(type)))
#define SEND_BUF(buf, bufSize, type) 
#define RECV_BUF(buf, bufSize, type) 
#define BCAST_BUF(buf, bufSize,type,who)  
#define TRAVERSAL_LENGTH 5
#define messageSize(x) 0
#endif


#endif	/* end include guard  */
