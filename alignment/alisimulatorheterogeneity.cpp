//
//  alisimulatorheterogeneity.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulatorheterogeneity.h"

AliSimulatorHeterogeneity::AliSimulatorHeterogeneity(Params *params) :
AliSimulator(params) {
    rate_heterogeneity = tree->getRate();
}

AliSimulatorHeterogeneity::AliSimulatorHeterogeneity(AliSimulator *alisimulator){
    tree = alisimulator->tree;
    params = alisimulator->params;
    rate_heterogeneity = tree->getRate();
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulatorHeterogeneity::simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // estimate the sequence for the current neighbor
        (*it)->node->sequence.resize(sequence_length);
        
        for (int i = 0; i < sequence_length; i++)
        {
           // randomly select the state, considering it's dad states, and the transition_probability_matrix
            (*it)->node->sequence[i] = estimateStateWithRH(model, site_specific_rates[i], trans_matrix, max_num_states, (*it)->length, node->sequence[i]);
        }
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, (*it)->node, node);
    }
}

int AliSimulatorHeterogeneity::estimateStateWithRH(ModelSubst *model, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state)
{
    // compute the transition matrix
    model->computeTransMatrix(branch_length*rate, trans_matrix);
    
    // iteratively select the state, considering it's dad states, and the transition_probability_matrix
    int starting_index = dad_state*max_num_states;
    return getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states);
}

/**
    get site-specific rates based on Continuous Gamma Distribution
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length)
{
    RateContinuousGamma *rate_continuous_gamma = new RateContinuousGamma(rate_heterogeneity->getGammaShape(), params->ran_seed);
    
    rate_continuous_gamma->getSiteSpecificRates(site_specific_rates, sequence_length);
}

/**
    get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesDiscrete(double *site_specific_rates, int sequence_length)
{
    int num_rate_categories = rate_heterogeneity->getNDiscreteRate();
    
    // initialize the probability array of rate categories
    double *category_probability_matrix = new double[num_rate_categories];
    for (int i = 0; i < num_rate_categories; i++)
        category_probability_matrix[i] = rate_heterogeneity->getProp(i);
    
    // convert the probability matrix of rate categories into an accumulated probability matrix of rate categories
    convertProMatrixIntoAccumulatedProMatrix(category_probability_matrix, 1, num_rate_categories);
    
    // initialize the site-specific rates
    for (int i = 0; i < sequence_length; i++)
    {
        // randomly select a rate from the set of rate categories, considering its probability array.
        int rate_category = getRandomItemWithAccumulatedProbabilityMatrix(category_probability_matrix, 0, num_rate_categories);
        
        // if rate_category == -1 <=> this site is invariant -> return dad's state
        if (rate_category == -1)
            site_specific_rates[i] = 0;
        else // otherwise, get the rate of that rate_category
            site_specific_rates[i] = rate_heterogeneity->getRate(rate_category);
    }
    
    // delete the probability array of rate categories
    delete[] category_probability_matrix;
}

/**
    get site-specific rates
*/
void AliSimulatorHeterogeneity::getSiteSpecificRates(double *site_specific_rates, int sequence_length)
{
    string rate_name = tree->getRateName();
    
    // initalize rates based on continuous gamma distribution
    if ((rate_name.find("+G") != std::string::npos) && params->alisim_continuous_gamma)
    {
        getSiteSpecificRatesContinuousGamma(site_specific_rates, sequence_length);
    }
    // initalize rates based on discrete distribution (gamma/freerate)
    else
    {
        getSiteSpecificRatesDiscrete(site_specific_rates, sequence_length);
    }
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulatorHeterogeneity::simulateSeqsForTree(){
    // get variables
    int sequence_length = params->alisim_sequence_length;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // initialize site-specific rates
    double *site_specific_rates = new double[sequence_length];
    getSiteSpecificRates(site_specific_rates, sequence_length);
    
    simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
        
    // delete the site-specific rates
    delete[] site_specific_rates;
    
    // delete trans_matrix array
    delete[] trans_matrix;
}
