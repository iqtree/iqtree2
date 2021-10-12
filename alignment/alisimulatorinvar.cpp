//
//  alisimulatorinvar.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulatorinvar.h"

AliSimulatorInvar::AliSimulatorInvar(Params *params, double invar_prop) :
AliSimulator(params) {
    invariant_proportion = invar_prop;
}

AliSimulatorInvar::AliSimulatorInvar(AliSimulator *alisimulator, double invar_prop){
    tree = alisimulator->tree;
    params = alisimulator->params;
    num_sites_per_state = alisimulator->num_sites_per_state;
    length_ratio = alisimulator->length_ratio;
    expected_num_sites = alisimulator->expected_num_sites;
    partition_rate = alisimulator->partition_rate;
    invariant_proportion = invar_prop;
    max_length_taxa_name = alisimulator->max_length_taxa_name;
    fundi_items = alisimulator->fundi_items;
    STATE_UNKNOWN = alisimulator->STATE_UNKNOWN;
    max_num_states = alisimulator->max_num_states;
}

/**
    simulate a sequence for a node from a specific branch after all variables has been initializing
*/
void AliSimulatorInvar::simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths)
{
    // compute the transition probability matrix
    model->computeTransMatrix(partition_rate * params->alisim_branch_scale * (*it)->length, trans_matrix);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(trans_matrix, max_num_states, max_num_states);
    
    // estimate the sequence for the current neighbor
    (*it)->node->sequence.resize(sequence_length);
    for (int i = 0; i < sequence_length; i++)
    {
        
        // if this site is invariant or the parent's state is a gap -> preserve the dad's state
        if (site_specific_rates[i] == 0 || node->sequence[i] == STATE_UNKNOWN)
            (*it)->node->sequence[i] = node->sequence[i];
        else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
        {
            int starting_index = node->sequence[i]*max_num_states;
            (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(trans_matrix, starting_index, max_num_states, node->sequence[i]);
        }
    }
}

/**
    initialize variables (e.g., site-specific rate)
*/
void AliSimulatorInvar::initVariables(int sequence_length, bool regenerate_root_sequence)
{
    initSiteSpecificRates(site_specific_rates, sequence_length);
}

/**
  initialize site_specific_rate
*/
void AliSimulatorInvar::initSiteSpecificRates(vector<double> &site_specific_rates, int sequence_length)
{
    site_specific_rates.resize(sequence_length, 1);
    for (int i = 0; i < sequence_length; i++)
    {
        // if this site is invariant -> preserve the dad's state
        if (random_double() <= invariant_proportion)
            site_specific_rates[i] = 0;
        else
            site_specific_rates[i] = 1;
    }
}

/**
*  insert a new sequence into the current sequence
*
*/
void AliSimulatorInvar::insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence)
{
    // initialize new_site_specific_rates for new sequence
    vector<double> new_site_specific_rates;
    initSiteSpecificRates(new_site_specific_rates, new_sequence.size());
    
    // insert new_site_specific_rates into site_specific_rates
    site_specific_rates.insert(site_specific_rates.begin()+position, new_site_specific_rates.begin(), new_site_specific_rates.end());
    
    // insert new_sequence into the current sequence
    AliSimulator::insertNewSequenceForInsertionEvent(indel_sequence, position, new_sequence);
}
