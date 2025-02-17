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
PhyloSuperHmm::PhyloSuperHmm(SuperAlignment *alignment, Params &params) : PhyloSuperTree(alignment,false,false) {
    vector<Alignment*>::iterator it;
    for (it = alignment->partitions.begin(); it != alignment->partitions.end(); it++) {
        push_back(new IQTreeMixHmm(params, (*it)));
    }
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
    if (size() == 0)
        return;
    
    iterator it = begin();
    int ntree = ((IQTreeMixHmm*)(*it))->size();
    int* numSiteCat = new int[ntree];
    int* numSiteCatTot = new int[ntree];
    int numSiteSum = 0;
    int k;
    memset(numSiteCatTot, 0, sizeof(int)*ntree);
    for (iterator it = begin(); it != end(); it++) {
        string filename = prefix + "." + (*it)->aln->name + ext;
        ((IQTreeMixHmm*)(*it))->printResults(filename.c_str(), cat_assign_method, numSiteCat);
        for (k = 0; k < ntree; k++)
            numSiteCatTot[k] += numSiteCat[k];
    }
    for (k = 0; k < ntree; k++)
        numSiteSum += numSiteCatTot[k];
    
    // print out overall percentage of sites over the trees
    cout << "Overall percentage of sites over the trees:";
    for (k = 0; k < ntree; k++)
        cout << " " << (double) numSiteCatTot[k] / numSiteSum;
    cout << endl;
    
    delete[] numSiteCat;
    delete[] numSiteCatTot;
}
