//
// Created by tamara on 3/8/21.
//

#ifndef IQTREE_NEURALNETWORK_H
#define IQTREE_NEURALNETWORK_H

#include <onnxruntime_cxx_api.h>
#include <alignment/alignment.h>

class NeuralNetwork {
public:
    /** constructor */
    NeuralNetwork(Alignment *alignment);

    /** destructor */
    virtual ~NeuralNetwork();

    double doAlphaInference();
    string doModelInference();

    Alignment *alignment;

};


#endif //IQTREE_NEURALNETWORK_H
