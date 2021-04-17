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
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulatorHeterogeneityInvar::simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // estimate the sequence for the current neighbor
        (*it)->node->sequence.resize(sequence_length);
        
        // check if trans_matrix could be caching (without rate_heterogeneity or the num of rate_categories is lowr than the threshold (5)) or not
        if (tree->getRateName().empty()
            || (!params->alisim_continuous_gamma && rate_heterogeneity && rate_heterogeneity->getNDiscreteRate() <= params->alisim_max_rate_categories_for_applying_caching))
        {
            int num_models = tree->getModel()->isMixture()?tree->getModel()->getNMixtures():1;
            int num_rate_categories  = tree->getRateName().empty()?1:rate_heterogeneity->getNDiscreteRate();
            double *cache_trans_matrix = new double[num_models*num_rate_categories*max_num_states*max_num_states];
            
            intializeCachingAccumulatedTransMax(cache_trans_matrix, num_models, num_rate_categories, max_num_states, (*it)->length, trans_matrix, model);

            // estimate the sequence
            for (int i = 0; i < sequence_length; i++)
            {
                // if this site is invariant -> preserve the dad's state
                if (site_specific_rates[i] == 0)
                    (*it)->node->sequence[i] = node->sequence[i];
                else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
                {
                    int model_index = site_specific_model_index[i];
                    int category_index = site_specific_rate_index[i];
                    int starting_index = model_index*num_rate_categories*max_num_states*max_num_states + category_index*max_num_states*max_num_states + max_num_states*node->sequence[i];
                
                    (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbabilityMatrix(cache_trans_matrix, starting_index, max_num_states);
                }
            }
            
            // delete cache_trans_matrix
            delete [] cache_trans_matrix;
        }
        // otherwise, estimating the sequence without trans_matrix caching
        else
        {
            for (int i = 0; i < sequence_length; i++)
            {
                // if this site is invariant -> preserve the dad's state
                if (site_specific_rates[i] == 0)
                    (*it)->node->sequence[i] = node->sequence[i];
                else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
                {
                    (*it)->node->sequence[i] = estimateStateWithRH(model, site_specific_model_index[i], site_specific_rates[i], trans_matrix, max_num_states, (*it)->length, node->sequence[i]);
                }
            }
        }
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, (*it)->node, node);
    }
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
