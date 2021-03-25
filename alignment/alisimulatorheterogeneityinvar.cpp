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
        
        for (int i = 0; i < sequence_length; i++)
        {
            // if this site is invariant -> preserve the dad's state
            if (site_specific_rates[i] == 0)
                (*it)->node->sequence[i] = node->sequence[i];
            else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
            {
                (*it)->node->sequence[i] = estimateStateWithRH(model, site_specific_rates[i], trans_matrix, max_num_states, (*it)->length, node->sequence[i]);
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

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulatorHeterogeneityInvar::simulateSeqsForTree(){
    // get variables
    int sequence_length = params->alisim_sequence_length;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // initialize site-specific rates
    double *site_specific_rates = new double[sequence_length];
    getSiteSpecificRates(site_specific_rates, sequence_length);
    
    // simulate Sequences
    simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
        
    // delete the site-specific rates
    delete[] site_specific_rates;
    
    // delete trans_matrix array
    delete[] trans_matrix;
}
