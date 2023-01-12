#include "ParallelContext.hpp"
#include <terrace/global_vars.h>
#include <omp.h>

ParallelContext::ParallelContext(int _num_threads){
    num_threads = _num_threads;
    taskCount = 0;
    
    for (int i = 0; i<num_threads; i++){
         working[i] = false;
    }

    complete = false;
    condition_met = false;

    omp_init_lock(&test_terrace_trees_mutex);
    omp_init_lock(&test_intermediate_trees_mutex);
    omp_init_lock(&time_mutex);
    omp_init_lock(&thread_mutex);
    //omp_init_lock(&mutexQueue);
    omp_init_lock(&taskReader);
    omp_init_lock(&printLocker);
    
    //pthread_attr_init(&detachedThread);
    //pthread_attr_setdetachstate(&detachedThread, PTHREAD_CREATE_DETACHED);

    // Initializing mutexes
    
    /* 
    pthread_mutex_init(&test_terrace_trees_mutex, NULL); // Initialization
    pthread_mutex_init(&test_intermediate_trees_mutex, NULL); // Initialization
    pthread_mutex_init(&time_mutex, NULL); // Initialization
    pthread_mutex_init(&thread_mutex, NULL); // Initialization
    pthread_mutex_init(&working_mutex, NULL); // Initialization
    pthread_mutex_init(&mutexQueue, NULL);
    
    pthread_cond_init(&condQueue, NULL);
    pthread_attr_init(&detachedThread);
    pthread_attr_setdetachstate(&detachedThread, PTHREAD_CREATE_DETACHED); */
}

ParallelContext::~ParallelContext(){

    // destroying mutexes
    omp_destroy_lock(&test_terrace_trees_mutex);
    omp_destroy_lock(&test_intermediate_trees_mutex);
    omp_destroy_lock(&time_mutex);
    omp_destroy_lock(&thread_mutex);
    //omp_destroy_lock(&mutexQueue);
    omp_destroy_lock(&taskReader);
    omp_destroy_lock(&printLocker);

    /* 
    pthread_mutex_destroy(&test_terrace_trees_mutex); // Initialization
    pthread_mutex_destroy(&test_intermediate_trees_mutex); // Initialization
    pthread_mutex_destroy(&time_mutex); // Initialization
    pthread_mutex_destroy(&thread_mutex); // Initialization
    pthread_mutex_destroy(&working_mutex); // Initialization

    pthread_mutex_destroy(&mutexQueue);
    pthread_cond_destroy(&condQueue);
    pthread_attr_destroy(&detachedThread);
    */
}

bool ParallelContext::CheckIfFinished(){
    
    for (int i = 0; i<num_threads; i++){
        if (working[i]){
            return false;
        }
    }

    return true;

}
