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
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulatorInvar::simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // compute the transition probability matrix
        model->computeTransMatrix(partition_rate*(*it)->length, trans_matrix);
        
        // convert the probability matrix into an accumulated probability matrix
        convertProMatrixIntoAccumulatedProMatrix(trans_matrix, max_num_states, max_num_states);
        
        // estimate the sequence for the current neighbor
        (*it)->node->sequence.resize(sequence_length);
        for (int i = 0; i < sequence_length; i++)
        {
            
            // if this site is invariant -> preserve the dad's state
            if (site_specific_rates[i] == 0)
                (*it)->node->sequence[i] = node->sequence[i];
            else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
            {
                int starting_index = node->sequence[i]*max_num_states;
                (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(trans_matrix, starting_index, max_num_states, node->sequence[i]);
            }
        }
        
        // update the num_children_done_simulation
        node->num_children_done_simulation++;
        // remove the sequence of
        if (!node->isLeaf() && node->num_children_done_simulation >= (node->neighbors.size() - 1))
            vector<int>().swap(node->sequence);
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, (*it)->node, node);
    }
}


/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulatorInvar::simulateSeqsForTree()
{
    // get variables
    int sequence_length = expected_num_sites;
    double invariant_proportion = tree->getRate()->getPInvar();
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // simulate Sequences
    // initialize the site-specific rates
    double *site_specific_rates = new double[sequence_length];
    for (int i = 0; i < sequence_length; i++)
    {
        // if this site is invariant -> preserve the dad's state
        if (random_double() <= invariant_proportion)
            site_specific_rates[i] = 0;
        else
            site_specific_rates[i] = 1;
    }
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // simulate sequences with only Invariant sites option
    simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);

    // delete trans_matrix array
    delete[] trans_matrix;
    
    // delete the site-specific rates
    delete[] site_specific_rates;
    
    // removing constant states if it's necessary
    if (length_ratio > 1)
        removeConstantSites();
}
