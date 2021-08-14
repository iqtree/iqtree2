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
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulatorInvar::simulateSeqsForTree(map<string,string> input_msa, string output_filepath)
{
    // get variables
    int sequence_length = expected_num_sites;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    ostream *out = NULL;
    vector<string> state_mapping;
    
    // initialize the site-specific rates
    vector<double> site_specific_rates(sequence_length);
    initVariables(sequence_length, site_specific_rates);
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // write output to file (if output_filepath is specified)
    if (output_filepath.length() > 0)
    {
        try {
            // add ".phy" or ".fa" to the output_filepath
            if (params->aln_output_format != IN_FASTA)
                output_filepath = output_filepath + ".phy";
            else
                output_filepath = output_filepath + ".fa";
            if (params->do_compression)
                out = new ogzstream(output_filepath.c_str());
            else
                out = new ofstream(output_filepath.c_str());
            out->exceptions(ios::failbit | ios::badbit);

            // write the first line <#taxa> <length_of_sequence> (for PHYLIP output format)
            if (params->aln_output_format != IN_FASTA)
            {
                int num_leaves = tree->leafNum - ((tree->root->isLeaf() && tree->root->name == ROOT_NAME)?1:0);
                *out <<num_leaves<<" "<< round(expected_num_sites/length_ratio)*num_sites_per_state<< endl;
            }

            // initialize state_mapping (mapping from state to characters)
            initializeStateMapping(num_sites_per_state, tree->aln, state_mapping);
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, output_filepath);
        }
    }
    
    // simulate sequences with only Invariant sites option
    simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root, *out, state_mapping, input_msa);
    
    // close the file if neccessary
    if (output_filepath.length() > 0)
    {
        if (params->do_compression)
            ((ogzstream*)out)->close();
        else
            ((ofstream*)out)->close();
        delete out;
        
        // show the output file name
        cout << "An alignment has just been exported to "<<output_filepath<<endl;
    }

    // delete trans_matrix array
    delete[] trans_matrix;
    
    // removing constant states if it's necessary
    if (length_ratio > 1)
        removeConstantSites();
}

/**
    simulate a sequence for a node from a specific branch after all variables has been initializing
*/
void AliSimulatorInvar::simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, vector<double> site_specific_rates, double *trans_matrix, int max_num_states, Node *node, NeighborVec::iterator it, string lengths)
{
    // compute the transition probability matrix
    model->computeTransMatrix(partition_rate*(*it)->length, trans_matrix);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(trans_matrix, max_num_states, max_num_states);
    
    // estimate the sequence for the current neighbor
    (*it)->node->sequence.resize(sequence_length);
    for (int i = 0; i < sequence_length; i++)
    {
        
        // if this site is invariant or the parent's state is a gap -> preserve the dad's state
        if (site_specific_rates[i] == 0 || node->sequence[i] == max_num_states)
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
void AliSimulatorInvar::initVariables(int sequence_length, vector<double> &site_specific_rates)
{
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
void AliSimulatorInvar::insertNewSequenceForInsertionEvent(Node *node, int position, vector<short int> &new_sequence, vector<double> &site_specific_rates)
{
    // initialize new_site_specific_rates for new sequence
    vector<double> new_site_specific_rates(new_sequence.size());
    initVariables(new_sequence.size(), new_site_specific_rates);
    
    // insert new_site_specific_rates into site_specific_rates
    site_specific_rates.insert(site_specific_rates.begin()+position, new_site_specific_rates.begin(), new_site_specific_rates.end());
    
    // insert new_sequence into the current sequence
    AliSimulator::insertNewSequenceForInsertionEvent(node, position, new_sequence, site_specific_rates);
}
