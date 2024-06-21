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

    static void getModelsAboveThreshold(StrVector* model_names, float* floatarr, int element_count);

    Alignment *alignment;

#ifdef _CUDA
    static float gpu_time;

private:
    cudaEvent_t start, stop;
#endif

};



#endif //IQTREE_NEURALNETWORK_H
