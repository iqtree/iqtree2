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
        get site-specific rates based on Continuous Gamma Distribution
    */
    virtual void getSiteSpecificRatesContinuousGamma(vector<double> &site_specific_rates, int sequence_length);
    
    /**
      estimate the state from accumulated trans_matrices
    */
    virtual int estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int dad_state, int* rstream);
    
    /**
      estimate the state from an original trans_matrix
    */
    virtual int estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, double branch_length, int dad_state, int site_index, int* rstream);
    
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
