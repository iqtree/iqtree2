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
    string doModelInference();

    Alignment *alignment;

#ifdef _CUDA
    static float gpu_time;

private:
    cudaEvent_t start, stop;
#endif

};



#endif //IQTREE_NEURALNETWORK_H
