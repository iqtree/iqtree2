//
//  alisimulatorheterogeneity.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#ifndef alisimulatorheterogeneity_h
#define alisimulatorheterogeneity_h

#include "alisimulator.h"

class AliSimulatorHeterogeneity : public AliSimulator
{
protected:

    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad, ofstream &out, vector<string> state_mapping, string &output);
    
    /**
        get site-specific rates based on Continuous Gamma Distribution
    */
    virtual void getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length);
    
    /**
        get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
    */
    void getSiteSpecificRatesDiscrete(double *site_specific_rates, int sequence_length);
    
    /**
      estimate the state from accumulated trans_matrices
    */
    virtual int estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int max_num_states, int dad_state);
    
    /**
      estimate the state from an original trans_matrix
    */
    virtual int estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state);
    
    /**
        initialize site specific model index based on its weights in the mixture model
    */
    void intializeSiteSpecificModelIndex();
    
    /**
        initialize state freqs for all model components (of a mixture model)
    */
    virtual void intializeStateFreqsMixtureModel();
    
    /**
        initialize caching accumulated_trans_matrix
    */
    void intializeCachingAccumulatedTransMatrices(double *cache_trans_matrix, int num_models, int num_rate_categories, int max_num_states, DoubleVector branch_lengths, double *trans_matrix, ModelSubst* model);
    
public:
    
    RateHeterogeneity *rate_heterogeneity;
    vector<short int> site_specific_model_index;
    vector<short int> site_specific_rate_index;
    const int RATE_ZERO_INDEX = -1;
    const int RATE_ONE_INDEX = 0;
    
    /**
        constructor
    */
    AliSimulatorHeterogeneity(Params *params);
    
    /**
        constructor
    */
    AliSimulatorHeterogeneity(AliSimulator *alisimulator);
    
    /**
        get site-specific rates
    */
    void getSiteSpecificRates(double *site_specific_rates, int sequence_length);

    /**
    *  simulate sequences for all nodes in the tree
    */
    virtual void simulateSeqsForTree(string output_filepath);
};

#endif /* alisimulatorheterogeneity_h */
