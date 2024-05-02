/*
 *  ParallelContext.h
 *  ParallelContext class is used for Parallelization Purposes
 *  Created on: Feb 15, 2022
 *      Author: Anastasis
 */

#ifndef TERRACE_PARALLELCONTEXT
#define TERRACE_PARALLELCONTEXT

#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include "global_vars.h"
#include <pthread.h>
#include <atomic>
#include <omp.h>

class ParallelContext{

    public:
    int num_threads;
    
    int taskCount;

    bool working[100];
    bool complete;
    bool condition_met;
    
    /* pthread_mutex_t test_terrace_trees_mutex;
    pthread_mutex_t test_intermediate_trees_mutex;
    pthread_mutex_t time_mutex;
    pthread_mutex_t thread_mutex;
    pthread_mutex_t mutexQueue;
    pthread_cond_t condQueue;
    pthread_attr_t detachedThread;
    */
    
    omp_lock_t test_terrace_trees_mutex;
    omp_lock_t test_intermediate_trees_mutex;
    omp_lock_t time_mutex;
    omp_lock_t thread_mutex;
    omp_lock_t taskReader;
    omp_lock_t printLocker;
    
    std::mutex mutexQueue;
    std::condition_variable condQueue;

    /* 
    std::mutex test_terrace_trees_mutex;
    std::mutex test_intermediate_trees_mutex;
    std::mutex time_mutex;
    std::mutex thread_mutex;
    std::mutex taskReader;
     */

    // have to define the condition variable
    
    // single barrier for all threads before they begin
    // pthread_barrier_t begin_barrier;

    ParallelContext(int _num_threads);
    ~ParallelContext();

    bool CheckIfFinished();

};


#endif /* TERRACE_PARALLELCONTEXT */


// Files that I need to modify:
// terrace.hpp, terrace.cpp
// presenceabsencematric.hpp and cpp
// paralllelcontext.cpp and hpp
