//
// Created by tamara on 3/8/21.
//

#ifndef IQTREE_NEURALNETWORK_H
#define IQTREE_NEURALNETWORK_H

#include <onnxruntime_cxx_api.h>
#include <alignment/alignment.h>

#ifdef _CUDA
#include "cuda_runtime.h"
#endif

#ifdef _IQTREE_MPI
#include "utils/MPIHelper.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils/timeutil.h"

class NeuralNetwork {
public:
    /** constructor */
    NeuralNetwork(Alignment *alignment);

    /** destructor */
    virtual ~NeuralNetwork();

    double doAlphaInference();

    /** do nnmodelFind inference
     * @param with_mf: this will output the probabilities
     * */

    string doModelInference(StrVector* model_names = nullptr);

    static void getModelsAboveThreshold(StrVector* model_names, DoubleVector model_probabilities);

    Alignment *alignment;


    // for measure time

    void initializeTimer();

    void startTimer();

    void stopTimer();

    static bool time_initialized;


#if defined(_OPENMP) && defined (_IQTREE_MPI)
    static DoubleVector run_time_array; // run time for each thread
    static double cpu_time;
    static double wall_time;
#elif defined(_OPENMP)
    static DoubleVector run_time_array; // run time for each thread
    static DoubleVector cpu_time_array;
    static DoubleVector wall_time_array;
#elif defined(_IQTREE_MPI)
    static double cpu_time;
    static double wall_time;
#else
    static double cpu_time;
    static double wall_time;
#endif


#ifdef _CUDA
    static float gpu_time;
#endif



private:
#ifdef _CUDA
    cudaEvent_t start, stop;
#endif
    double local_cpu_time;
    double local_wall_time;
};



#endif //IQTREE_NEURALNETWORK_H
