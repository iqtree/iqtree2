//
//  alisimulatorheterogeneityinvar.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#ifndef alisimulatorheterogeneityinvar_h
#define alisimulatorheterogeneityinvar_h

#include "alisimulatorheterogeneity.h"

class AliSimulatorHeterogeneityInvar : public AliSimulatorHeterogeneity
{
protected:
    double invariant_proportion;

    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad);
    
    /**
        get site-specific rates based on Continuous Gamma Distribution
    */
    virtual void getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length);
    
public:
    
    /**
        constructor
    */
    AliSimulatorHeterogeneityInvar(Params *params, double invar_prop);
    
    /**
        constructor
    */
    AliSimulatorHeterogeneityInvar(AliSimulator *alisimulator, double invar_prop);
};

#endif /* alisimulatorheterogeneityinvar_h */
