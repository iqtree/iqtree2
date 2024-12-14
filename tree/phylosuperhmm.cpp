//
//  phylosuperhmm.cpp
//
//  Created by Thomas Wong on 5/12/24.
//

#include "phylosuperhmm.h"

/**
    constructor
*/
PhyloSuperHmm::PhyloSuperHmm() : PhyloSuperTree() {
}

/**
 constructor
 */
PhyloSuperHmm::PhyloSuperHmm(SuperAlignment *alignment, Params &params, int numTree) : PhyloSuperTree(alignment,false,false) {
    vector<Alignment*>::iterator it;
    for (it = alignment->partitions.begin(); it != alignment->partitions.end(); it++) {
        vector<IQTree*> trees;
        for (int i=0; i<numTree; i++) {
            trees.push_back(new IQTree(*it));
        }
        push_back(new IQTreeMixHmm(params, (*it), trees));
    }
}

/**
    constructor
*/
PhyloSuperHmm::PhyloSuperHmm(SuperAlignment *alignment, PhyloSuperHmm *super_hmm) : PhyloSuperTree(alignment, super_hmm) {
}

/**
 destructor
 */
PhyloSuperHmm::~PhyloSuperHmm() {
    model_factory = NULL;
    model = NULL;
    site_rate = NULL;
    for (reverse_iterator it = rbegin(); it != rend(); it++) {
        delete (*it);
    }
    clear();
}

/**
    set minimum branch length
*/
void PhyloSuperHmm::setMinBranchLen(Params& params) {
    
    if (params.min_branch_length <= 0.0) {
        params.min_branch_length = MAST_MIN_BRANCH_LEN;
    }
    cout << setprecision(7) << "Minimum branch length is set to " << params.min_branch_length << endl;
}

void PhyloSuperHmm::initSettings(Params &params) {
    IQTree::initSettings(params);
    setLikelihoodKernel(params.SSE);
    setNumThreads(params.num_threads);
    for (iterator it = begin(); it != end(); it++) {
        ((IQTreeMixHmm*)(*it))->initSettings(params);
    }
}

/*
void PhyloSuperHmm::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    for (iterator it = begin(); it != end(); it++) {
        ((IQTreeMixHmm*)(*it))->initializeModel(params, model_name, models_block);
    }
}
*/

/**
 * Generate the initial tree (usually used for model parameter estimation)
 */
void PhyloSuperHmm::computeInitialTree(LikelihoodKernel kernel, istream* in) {
    for (iterator it = begin(); it != end(); it++) {
        ((IQTreeMixHmm*)(*it))->computeInitialTree(kernel, in);
    }
}

void PhyloSuperHmm::setRootNode(const char *my_root, bool multi_taxa) {
    for (iterator it = begin(); it != end(); it++) {
        ((IQTreeMixHmm*)(*it))->setRootNode(my_root, multi_taxa);
    }
}

// show the assignment of the categories along sites with max likelihood
// cat_assign_method:
//  0 - the categories along sites is assigned according to the path with maximum probability (default)
//  1 - the categories along sites is assigned according to the max posterior probability
void PhyloSuperHmm::printResults(string prefix, string ext, int cat_assign_method) {
    for (iterator it = begin(); it != end(); it++) {
        string filename = prefix + "." + (*it)->aln->name + ext;
        ((IQTreeMixHmm*)(*it))->printResults(filename.c_str(), cat_assign_method);
    }
}
