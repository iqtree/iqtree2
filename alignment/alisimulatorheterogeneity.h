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
    int estimateStateWithRH(ModelSubst *model, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state);
    
public:
    
    RateHeterogeneity *rate_heterogeneity;
    
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
