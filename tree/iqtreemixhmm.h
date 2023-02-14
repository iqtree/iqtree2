//
//  iqtreemixhmm.h
//  tree
//
//  Created by Thomas Wong on 19/01/23.
//

#ifndef iqtreemixhmm_h
#define iqtreemixhmm_h

#include <cmath>
#include "iqtreemix.h"
#include "tree/phylohmm.h"
#include "model/modelhmm.h"
#include "model/modelhmmgm.h"

class IQTreeMixHmm : public IQTreeMix, public PhyloHmm {
public:
    
    /**
     default constructor
     */
    IQTreeMixHmm();
    
    IQTreeMixHmm(Params &params, Alignment *aln, vector<IQTree*> &trees);
    
    // initialize the model
    void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    // initialize the transition model
    void initializeTransitModel(Params &params);
    
    // initialize the parameters
    void initializeParams();
    
    // compute the log of dotproduct of the logorithm arrays
    double logDotProd(double* ln_x, double* ln_y, int n);
    
    // obtain the log-likelihoods for every pattern and every tree
    // output ptn_like_cat[i * ntree + j] : log-likelihood of pattern i and tree j
    void computeLogLikelihoodSiteTree(int updateTree);
    
    // compute backward log-likelihood
    virtual double computeLikelihood(double *pattern_lh = NULL, bool save_log_value = true);
    
    /**
     optimize all branch lengths of one tree
     @param iterations number of iterations to loop through all branches
     */
    void optimizeAllBranchesOneTree(int whichtree, int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);
    
    /**
     optimize all branch lengths of all trees
     @param iterations number of iterations to loop through all branches
     @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);
    
    double optimizeAllBranchLensByBFGS(double gradient_epsilon, double logl_epsilon, int maxsteps = 3);
    /**
     @return true if this is a HMM model
     */
    virtual bool isHMM() { return true; }
    
    virtual void startCheckpoint();
    
    virtual string optimizeModelParameters(bool printInfo, double logl_epsilon);
    
    virtual void setNumThreads(int num_threads);
    
    /**
     test the best number of threads
     */
    virtual int testNumThreads();
    
    virtual int getNParameters();
    
    // print out all the results to a file
    void printResults(const char *filename);
    
    // show the values of the parameters
    void showParameters(ostream& out);
    
    // optimize all substitution models
    double optimizeAllSubstModels(double gradient_epsilon);
    
    // optimize all RHAS models
    double optimizeAllRHASModels(double gradient_epsilon, double score = 0.0);
    
private:
    
    // indicate which tree's unlinked parameters is under optimization
    // -1   : parameters affecting all trees (default)
    // >= 0 : parameters affecting a specific tree
    int optimTree;
    
    // indicate which branch group is under optimization
    // -1  : no branch group is under optimization
    // >=0 : optimizing a specific branch group
    int optimBranchGrp;
    
    // branch lengths of all the trees
    vector<DoubleVector> allbranchlens;
    
    // compute the log-likelihoods for a single tree t
    void computeLogLikelihoodSingleTree(int t);

    // get the branch lengths of all trees to the variable allbranchlens
    void getAllBranchLengths();

    // set the branch lengths of all trees from the variable allbranchlens
    void setAllBranchLengths();
    
    // show the branch lengths of all trees
    void showAllBranchLengths();
    
    //--------------------------------------------
    // optimization of branch lengths using BFGS
    //--------------------------------------------

    // the following three functions are for dimension = 1
    double computeFunction(double x);
    
    double setSingleVariable();
    
    void getSingleVariable(double x);

    // the following four functions are for dimension > 1
    virtual double targetFunk(double x[]);
    
    virtual void setVariables(double *variables);
    
    virtual void getVariables(double *variables);
    
    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

    virtual int getNDim();
    
    double optimizeBranchGroup(int branchgrp, double gradient_epsilon);
    
    void showBranchGrp();

    /**
             If there are multiple branches belonging to the same group
             set all the branches of the same group to their average
     */
    void setAvgLenEachBranchGrp();
    
    // update the ptn_freq array according to the marginal probabilities along each site for each tree
    void computeFreqArray(bool need_computeLike = true, int update_which_tree = -1);
    
    // get marginal probabilities along each site for each tree
    void getMarginalProb(bool need_computeLike = true, int update_which_tree = -1);
};

#endif /* iqtreemixhmm_h */
