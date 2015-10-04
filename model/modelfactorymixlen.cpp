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
