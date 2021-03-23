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
        model->computeTransMatrix((*it)->length, trans_matrix);
        
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
                (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbabilityMatrix(trans_matrix, starting_index, max_num_states);
            }
            
        }
        
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
    int sequence_length = params->alisim_sequence_length;
    double invariant_proportion = tree->getRate()->getPInvar();
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
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
    
    // simulate sequences with only Invariant sites option
    simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
    
    // delete the site-specific rates
    delete[] site_specific_rates;
    
    // delete trans_matrix array
    delete[] trans_matrix;
}
