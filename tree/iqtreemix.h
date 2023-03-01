//
//  iqtreemix.hpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#ifndef iqtreemix_h
#define iqtreemix_h

#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "iqtree.h"
#include "utils/MPIHelper.h"
#include "model/modelmarkov.h"

class IQTreeMix : public IQTree, public vector<IQTree*> {
public:

    /**
            default constructor
     */
    IQTreeMix();

    IQTreeMix(Params &params, Alignment *aln, vector<IQTree*> &trees);

    /**
            destructor
     */
    virtual ~IQTreeMix();

    /**
            initialization
     */
    // void init(Alignment* aln, Params* params, Checkpoint* chkpt);
    
    void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    // compute the overall likelihood value by combining all the existing likelihood values of the trees
    double computeLikelihood_combine(double *pattern_lh = NULL);

    // this function is designed for the situation that
    // only one tree has been updated.
    double computeLikelihood_oneTreeUpdated(int whichTree);

    virtual double computeLikelihood(double *pattern_lh = NULL);

    virtual double computePatternLhCat(SiteLoglType wsl);

    /**
            compute pattern likelihoods only if the accumulated scaling factor is non-zero.
            Otherwise, copy the pattern_lh attribute
            @param pattern_lh (OUT) pattern log-likelihoods,
                            assuming pattern_lh has the size of the number of patterns
            @param cur_logl current log-likelihood (for sanity check)
            @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
     */
    virtual void computePatternLikelihood(double *pattern_lh = NULL, double *cur_logl = NULL,
            double *pattern_lh_cat = NULL, SiteLoglType wsl = WSL_RATECAT);

    virtual void initializeAllPartialLh();

    virtual void deleteAllPartialLh();

    virtual void clearAllPartialLH(bool make_null = false);

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

    /**
            compute the updated tree weights according to the likelihood values along each site
            prerequisite: computeLikelihood() has been invoked

     */
    double optimizeTreeWeightsByEM(double* pattern_mix_lh, double gradient_epsilon, int max_steps, bool& tree_weight_converge);
    double optimizeTreeWeightsByBFGS(double gradient_epsilon);

    double optimizeBranchLensByBFGS(double gradient_epsilon);

    // save branch lengths of all trees
    // node and dad are always NULL
    void getBranchLengths(vector<DoubleVector> &len, Node *node = NULL, Node *dad = NULL);

    // restore branch lengths of all trees
    // node and dad are always NULL
    void setBranchLengths(vector<DoubleVector> &len, Node *node = NULL, Node *dad = NULL);
    
    virtual void showTree();
    
    /** set the root by name
        @param my_root root node name
        @param multi_taxa TRUE if my_root is a comma-separated list of nodes
     */
    virtual void setRootNode(const char *my_root, bool multi_taxa = false);
    
    /**
        @return true if this is a mixture of trees, default: false
    */
    virtual bool isTreeMix() { return true; }

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    void startCheckpoint();

    void saveCheckpoint();
    
    void restoreCheckpoint();
    
    void setMinBranchLen(Params& params);

    /** set pointer of params variable */
    virtual void setParams(Params* params);

    /*
     * Generate the branch IDs
     * Branches of different trees with the same partition share the same ID
     */
    void computeBranchID();

    /**
     * Generate the initial tree (usually used for model parameter estimation)
     */
    void computeInitialTree(LikelihoodKernel kernel, istream* in = NULL);

    /**
     * setup all necessary parameters
     */
    virtual void initSettings(Params& params);

    /**
     * compute the memory size required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired(size_t ncategory = 1, bool full_mem = false);

    /**
     * compute the memory size for top partitions required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequiredThreaded(size_t ncategory = 1, bool full_mem = false);

    /**
        test the best number of threads
    */
    virtual int testNumThreads();

    /**
        Initialize the tree weights using parsimony scores
        Idea:
        1. Check the parsimony score for each tree along all the sites
        2. Select the sites with different parsimony score between the trees.
        3. For each selected site, we check which parsimony score of the tree is minimum, and assign the site to the tree.
        4. The tree weights are estimated according to the proportion of the sites assigned to each tree.
     */
    void initializeTreeWeights();
    void initializeTreeWeights2();

    string optimizeModelParameters(bool printInfo, double logl_epsilon);

    /**
            print tree to .treefile
            @param params program parameters, field root is taken
     */
    virtual void printResultTree(string suffix = "");

    /**
     * Return the tree string contining taxon names and branch lengths
     * @return
     */
    virtual string getTreeString();

    /**
     @return the average of the tree lengths
     @param node the starting node, NULL to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLength(Node *node = NULL, Node *dad = NULL);

    /**
     @return the average length of all internal branches
     @param node the starting node, NULL to start from the root
     @param dad dad of the node, used to direct the search
     */
    virtual double treeLengthInternal(double epsilon, Node *node = NULL, Node *dad = NULL);

    int getNParameters();

    virtual void drawTree(ostream &out, int brtype = WT_BR_SCALE + WT_INT_NODE, double zero_epsilon = 2e-6);

    /**
            print the tree to the output file in newick format
            @param out the output file.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param brtype type of branch to print
            @return ID of the taxon with smallest ID
     */
    virtual int printTree(ostream &out, int brtype, Node *node, Node *dad = NULL);

    /**
     *  Return best tree string from the candidate set
     *
     *  @param numTrees
     *      Number of best trees to return
     *  @return
     *      A string vector of trees
     */
    virtual vector<string> getBestTrees(int numTrees = 0);
    
    /**
            Read the tree saved with Taxon IDs and branch lengths.
            @param tree_string tree string to read from
            @param updatePLL if true, tree is read into PLL
     */
    virtual void readTreeString(const string &tree_string);

    virtual ModelFactory *getModelFactory() {
        return at(0)->getModelFactory();
    }

    /**
            get rate heterogeneity
            @return associated rate heterogeneity class
     */
    virtual RateHeterogeneity *getRate() {
        return at(0)->getRate();
    }

    virtual ModelSubst *getModel() {
        return at(0)->getModel();
    }

    // show the log-likelihoods and posterior probabilties for each tree along the sites
    void showLhProb(ofstream& out);
 
    // compute parsimony scores for each tree along the patterns
    // results are stored in the array patn_parsimony
    void computeParsimony();
    
    // show the log-likelihoods and posterior probabilties for each tree along the patterns
    void showPatternLhProb(ofstream& out);
    
    // show the log-likelihoods and posterior probabilties for each tree along the patterns
    void showOrderedPatternLhProb(ofstream& out);

    /**
            pattern frequencies
     */
    int* patn_freqs;
    
    /**
            whether pattern is constant
     */
    int* patn_isconst;

    /**
            parsimony scores for each tree along the patterns
     */
    int* patn_parsimony;

    /**
            weights of trees
     */
    vector<double> weights;
    vector<double> tmp_weights; // for optimization
    
    /**
            ratios of parsimony informative sites with max posterior probability for each tree
     */
    vector<double> max_posterior_ratio;
    
    /**
            ratios of parsimony informative sites with max high-enought likelihood for each tree
     */
    vector<double> max_like_ratio;

    /**
            pattern likelihoods for all trees
     */
    double* ptn_like_cat;
    
    /**
            models
     */
    vector<ModelSubst*> models;
    
    /**
            site rates
     */
    vector<RateHeterogeneity*> site_rates;

    /**
            trees assigned for the site rates
     */
    vector<PhyloTree*> site_rate_trees;

    /**
            members of each weight group
     */
    vector<vector<int> > weight_group_member;
    
    /**
            does weight group exist
     */
    bool weightGrpExist;

    /**
            branch ID
            branches of different trees with the same partition are assigned to the same ID
     */
    vector<int> branch_id;
    
    /**
            branch group
            groups of branches of different trees with the same partition
     */
    vector<IntVector> branch_group;
    
    /**
            branch lengths (for optimization)
     */
    vector<DoubleVector> branch_len;
    
    /**
            parsimony scores are computed for each tree along the patterns
     */
    bool parsi_computed;

    /**
            whether the submodels are linked
     */
    bool isLinkModel;
    
    /**
            whether the site rates are linked, if exists
     */
    bool isLinkSiteRate;

private:

    // to separate the submodel names and the site rate names from the full model name
    void separateModel(string modelName);
    
    // reset the ptn_freq array to the original frequencies of the patterns
    void resetPtnOrigFreq();

    /**
            update the ptn_freq array according to the posterior probabilities along each site for each tree

     */
    void computeFreqArray(double* pattern_mix_lh, bool need_computeLike, int update_which_tree = -1);

    /**
            get posterior probabilities along each site for each tree
     */
    void getPostProb(double* pattern_mix_lh, bool need_computeLike, int update_which_tree = -1);
    
    /**
        optimize tree k separately
     */
    void optimizeTreeSeparately(int k, bool printInfo, double gradient_epsilon);
    
    /**
            optimize each tree separately
     */
    void optimizeTreesSeparately(bool printInfo, double logl_epsilon);

    /**
             If there are multiple tree weights belonging to the same group
             set all the tree weights of the same group to their average
     */
    void checkWeightGrp();

    /**
             If there are multiple branches belonging to the same group
             set all the branches of the same group to their average
     */
    void checkBranchGrp();

    // -------------------------------------
    // for BFGS optimzation on tree weights
    // -------------------------------------

    double targetFunk(double x[]);

    // read the tree weights and write into "variables"
    void setVariables(double *variables);

    // read the "variables" and write into tree weights
    void getVariables(double *variables);

    // set the bounds
    void setBounds(double *lower_bound, double *upper_bound, bool* bound_check);
    
    // get the dimension of the variables (for tree weights)
    int getNDim();
    
    // optimization on which variable
    // 1 - tree weights
    // 2 - branch lengths
    int optim_type;
    
    /**
            immediate array for pattern likelihoods during computation
     */
    double* _ptn_like_cat;

    /**
            number of optimization steps, default: number of Trees * 2
     */
    int optimize_steps;
    
    /**
            inputted tree-mixture model name
     */
    string treemix_model;
    
    /**
            individual model names (only one item if linked, more than one otherwise)
     */
    vector<string> model_names;
    
    /**
            individual site rate names (only one item if linked, more than one otherwise)
     */
    vector<string> siterate_names;
    
    /**
            whether there is any site rate
     */
    bool anySiteRate;

    /**
            whether it is a edge-len-restricted model in which the edges of different trees having the same partition have the same lengths
     */
    bool isEdgeLenRestrict;

    /**
            number of trees
     */
    size_t ntree;

    /**
            number of tips
     */
    size_t ntip;

    /**
            number of patterns
     */
    size_t nptn;
    
    /**
            number of parsimony informative sites
     */
    size_t ninformsite;
    
    /**
            number of branches of each tree
     */
    size_t nbranch;

    /**
            initial models
     */
    vector<ModelSubst*> initial_models;
    
    /**
            initial site rates
     */
    vector<RateHeterogeneity*> initial_site_rates;
    
    /**
            is the tree weights fixed
     */
    bool isTreeWeightFixed;
};

#endif /* iqtreemix_h */
