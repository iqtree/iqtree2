/*
 * modelfactorymixlen.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: minh
 */

#include "phylotreemixlen.h"
#include "timeutil.h"
#include "model/modelfactorymixlen.h"
#include "model/modelgtr.h"
#include "model/modelmixture.h"

ModelFactoryMixlen::ModelFactoryMixlen(Params &params, PhyloTree *tree, ModelsBlock *models_block) :   
    ModelFactory(params, tree, models_block) {
    if (!tree->isMixlen()) {
        cerr << "ERROR: Please add '-mixlen " << site_rate->getNRate() << "' option into the command line" << endl;
        outError("Sorry for the inconvience, please rerun IQ-TREE with option above");
    }
    if (tree->getMixlen() != site_rate->getNRate()) {
        ((PhyloTreeMixlen*)tree)->setMixlen(site_rate->getNRate());
//        outError("#heterotachy classes and #mixture branch lengths do not match");
    }
}

double ModelFactoryMixlen::optimizeParameters(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {

	PhyloTreeMixlen *tree = (PhyloTreeMixlen*)site_rate->getTree();
	assert(tree);
    
    tree->initializeMixlen(logl_epsilon);

    return ModelFactory::optimizeParameters(fixed_len, write_info, logl_epsilon, gradient_epsilon);


}

int ModelFactoryMixlen::getNParameters() {
	int df = ModelFactory::getNParameters();
    df += site_rate->phylo_tree->branchNum * (site_rate->phylo_tree->getMixlen()-1);
	return df;
}
