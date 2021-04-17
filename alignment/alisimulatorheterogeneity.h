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
    virtual void simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad);
    
    /**
        get site-specific rates based on Continuous Gamma Distribution
    */
    virtual void getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length);
    
    /**
        get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
    */
    void getSiteSpecificRatesDiscrete(double *site_specific_rates, int sequence_length);
    
    /**
      estimate the state for the current site of a node in the case WITH Rate Heterogeneity (Gamma/FreeRate Model)
    */
    int estimateStateWithRH(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state);
    
    /**
        initialize site specific model index based on its weights in the mixture model
    */
    void intializeSiteSpecificModelIndex();
    
    /**
        initialize state freqs for all model components (of a mixture model)
    */
    void intializeStateFreqsMixtureModel();
    
    /**
        initialize caching accumulated_trans_matrix
    */
    void intializeCachingAccumulatedTransMax(double *cache_trans_matrix, int num_models, int num_rate_categories, int max_num_states, double branch_length, double *trans_matrix, ModelSubst* model);
    
public:
    
    RateHeterogeneity *rate_heterogeneity;
    IntVector site_specific_model_index;
    IntVector site_specific_rate_index;
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
    virtual void simulateSeqsForTree();
};

#endif /* alisimulatorheterogeneity_h */
