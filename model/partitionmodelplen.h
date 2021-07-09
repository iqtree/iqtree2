//
//  partitionmodelplen.h
//  iqtree
//
//  Created by Olga on 04/05/17.
//
//

#ifndef iqtree_partitionmodelplen_h
#define iqtree_partitionmodelplen_h

#include "tree/phylosupertreeplen.h"
#include "model/partitionmodel.h"

class PartitionModelPlen : public PartitionModel
{
public:
    PartitionModelPlen();
    /**
     constructor
     create partition model with possible rate heterogeneity. Create proper class objects
     for two variables: model and site_rate. It takes the following field of params into account:
     model_name, num_rate_cats, freq_type, store_trans_matrix
     @param params program parameters
     @param tree associated phylogenetic super-tree
     @param report_to_tree tree to report to
     */
    PartitionModelPlen(Params &params, PhyloSuperTreePlen *tree,
                       ModelsBlock *models_block,
                       PhyloTree* report_to_tree);
    
    ~PartitionModelPlen();
    
    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
     save object into the checkpoint
     */
    virtual void saveCheckpoint();
    
    /**
     restore object from the checkpoint
     */
    virtual void restoreCheckpoint();
    
    /**
     * @param brlen_type either BRLEN_OPTIMIZE, BRLEN_FIX or BRLEN_SCALE
     * @return #parameters of the model + # branches
     */
    virtual int getNParameters(int brlen_type) const;

    //virtual int getNDim() const;
    
    /**
     write information
     @param out output stream
     */
    virtual void writeInfo(ostream &out);
    
    /**
     optimize model parameters and tree branch lengths
     NOTE 2016-08-20: refactor the semantic of fixed_len
     @param fixed_len 0: optimize branch lengths, 1: fix branch lengths, 2: scale branch lengths
     @param write_info TRUE to write model parameters every optimization step, FALSE to only print at the end
     @param logl_epsilon log-likelihood epsilon to stop
     @param gradient_epsilon gradient (derivative) epsilon to stop
     @param report_to_tree a tree to report progress to
     @return the best likelihood
     */
    virtual double optimizeParameters(int fixed_len = BRLEN_OPTIMIZE, bool write_info = true,
                                      double logl_epsilon = 0.1, double gradient_epsilon = 0.0001,
                                      PhyloTree* report_to_tree = nullptr);
    
        double optimizeSubtreeModelParameters(PhyloSuperTreePlen *tree, 
                                              int iteration,
                                              double gradient_epsilon,
                                              PhyloTree* report_to_tree);
        double optimizeGeneRate(PhyloSuperTreePlen *tree, double tree_lh, 
                                double cur_lh, double gradient_episilon);

        double optimizeBranchLengthsIfRequested(PhyloSuperTreePlen *tree, 
                                                int    fixed_len, 
                                                int    iteration,
                                                double cur_lh, 
                                                double logl_epsilon,
                                                double gradient_epsilon);
        void   logResultOfParameterOptimization(bool   write_info,
                                                double begin_time,
                                                int    iteration,
                                                PhyloTree* report_to_tree);

    /**
     *  optimize model parameters and tree branch lengths for the +I+G model
     *  using restart strategy.
     * 	@param fixed_len TRUE to fix branch lengths, default is false
     *	@return the best likelihood
     */
    virtual double optimizeParametersGammaInvar(int fixed_len = BRLEN_OPTIMIZE,
                                                bool write_info = true,
                                                double logl_epsilon = 0.1,
                                                double gradient_epsilon = 0.0001,
                                                PhyloTree* report_to_tree = nullptr);
    
    double optimizeGeneRate(double tol);
    
    //	virtual double targetFunk(double x[]);
    //	virtual void getVariables(const double *variables);
    //	virtual void setVariables(double *variables);
    
    /** partition ID currently under optimization of of its rate */
    //    int optimizing_part;
    
    /**
     compute the likelihood for a partition under rate optimization (optimizing_rate).
     Used for the ML optimization of gene rate
     @param value x-value of the function
     @return log-likelihood
     */
    //    virtual double computeFunction(double value);
    
    
};


#endif
