//
//  phylosuperhmm.h
//
//  Created by Thomas Wong on 5/12/24.
//

#ifndef phylosuperhmm_h
#define phylosuperhmm_h

#include <stdio.h>
#include "phylosupertree.h"
#include "iqtreemixhmm.h"

class PhyloSuperHmm : public PhyloSuperTree
{
public:
    
    /**
     constructor
     */
    PhyloSuperHmm();
    
    /**
     constructor
     */
    PhyloSuperHmm(SuperAlignment *alignment, Params &params);
    
    /**
     destructor
     */
    ~PhyloSuperHmm();
    
    /**
     @return true if this is a mixture of trees, default: false
     */
    virtual bool isTreeMix() { return true; }
    
    /**
     set minimum branch length
     */
    void setMinBranchLen(Params& params);
    
    void initSettings(Params &params);
    
    // void initializeModel(Params &params, string model_name, ModelsBlock *models_block);
    
    /**
     * Generate the initial tree (usually used for model parameter estimation)
     */
    void computeInitialTree(LikelihoodKernel kernel, istream* in);
    
    void setRootNode(const char *my_root, bool multi_taxa);
    
    // show the assignment of the categories along sites with max likelihood
    // cat_assign_method:
    //  0 - the categories along sites is assigned according to the path with maximum probability (default)
    //  1 - the categories along sites is assigned according to the max posterior probability
    void printResults(string prefix, string ext, int cat_assign_method);
};

#endif
