/*
 * modelfactorymixlen.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: minh
 */

#include "modelfactorymixlen.h"

ModelFactoryMixlen::ModelFactoryMixlen(Params &params, PhyloTree *tree, ModelsBlock *models_block) :   
    ModelFactory(params, tree, models_block) {
    if (!model->isMixture())
        outError("Model is not mixture");
    if (params.num_mixlen != model->getNMixtures())
        outError("#mixture categories and #mixture branch lengths do not match");
}

double ModelFactoryMixlen::optimizeParameters(bool fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    double score = ModelFactory::optimizeParameters(fixed_len, write_info, logl_epsilon, gradient_epsilon);
    
    return score;
}
