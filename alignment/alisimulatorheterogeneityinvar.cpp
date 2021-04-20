//
//  alisimulatorheterogeneityinvar.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulatorheterogeneityinvar.h"

AliSimulatorHeterogeneityInvar::AliSimulatorHeterogeneityInvar(Params *params, double invar_prop) :
AliSimulatorHeterogeneity(params) {
    invariant_proportion = invar_prop;
}

AliSimulatorHeterogeneityInvar::AliSimulatorHeterogeneityInvar(AliSimulator *alisimulator, double invar_prop):AliSimulatorHeterogeneity(alisimulator){
    invariant_proportion = invar_prop;
}

/**
    get site-specific rates based on Continuous Gamma Distribution
*/
void AliSimulatorHeterogeneityInvar::getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length)
{
    RateContinuousGamma *rate_continuous_gamma = new RateContinuousGammaInvar(rate_heterogeneity->getGammaShape(), params->ran_seed, invariant_proportion);;
    
    rate_continuous_gamma->getSiteSpecificRates(site_specific_rates, sequence_length);
    
    // delete rate_continuous_gamma
    delete rate_continuous_gamma;
}

/**
  estimate the state from accumulated trans_matrices
*/
int AliSimulatorHeterogeneityInvar::estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int max_num_states, int dad_state)
{
    // if this site is invariant -> preserve the dad's state
    if (site_specific_rate == 0)
        return dad_state;
    
    // otherwise, randomly select the state, considering it's dad states, and the accumulated trans_matrices
    return AliSimulatorHeterogeneity::estimateStateFromAccumulatedTransMatrices(cache_trans_matrix, site_specific_rate, site_index, num_rate_categories, max_num_states, dad_state);
}

/**
  estimate the state from an original trans_matrix
*/
int AliSimulatorHeterogeneityInvar::estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state)
{
    // if this site is invariant -> preserve the dad's state
    if (rate == 0)
        return dad_state;
    
    // otherwise, select the state, considering it's dad states, and the transition_probability_matrix
    return AliSimulatorHeterogeneity::estimateStateFromOriginalTransMatrix(model, model_component_index, rate, trans_matrix, max_num_states, branch_length, dad_state);
}
