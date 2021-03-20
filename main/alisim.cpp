/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"


void runAliSim(Params params)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    // show parameters
    showParameters(params);
    
    // read input tree from file
    IQTree *tree = initializeIQTreeFromTreeFile(params, params.sequence_type);
    
    // iteratively generate multiple datasets for each tree
    for (int i = 0; i < params.alisim_dataset_num; i++)
    {
        // initialize output_filepath
        std::string output_filepath(params.user_file);
        output_filepath = output_filepath
        +"_"+params.alisim_output_filename
        +"_"+convertIntToString(i)+".phy";
        
        generateSingleDatasetFromSingleTree(params, tree, output_filepath);
    }
    
    cout << "[Alignment Simulator] Done"<<"\n";
}

void showParameters(Params params)
{
    cout << " - Tree filepath: " << params.user_file <<"\n";
    cout << " - Length of output sequences: " << params.alisim_sequence_length <<"\n";
    if (!params.model_name.empty())
        cout << " - Model: " << params.model_name <<"\n";
    cout << " - Number of output datasets: " << params.alisim_dataset_num<<"\n";
    if (params.alisim_ancestral_sequence >= 0)
        cout << " - Ancestral sequence position: " << params.alisim_dataset_num <<"\n";
}

IQTree *initializeIQTreeFromTreeFile(Params params, char* seq_type)
{
    IQTree *tree = new IQTree();
    bool is_rooted = false;
    tree->readTree(params.user_file, is_rooted);
    initializeAlignment(seq_type, tree);
    initializeModel(params, tree);
    return tree;
}

void initializeAlignment(char* seq_type, IQTree *tree)
{
    tree->aln = new Alignment();
    
    // set the seq_type and the maximum number of bases based on the Seq_type
    tree->aln->seq_type = tree->aln->getSeqType(seq_type);
    
    switch (tree->aln->seq_type) {
    case SEQ_BINARY:
        tree->aln->num_states = 2;
        break;
    case SEQ_PROTEIN:
        tree->aln->num_states = 20;
        break;
    case SEQ_MORPH:
        throw "Sorry! SEQ_MORPH is currently not supported";
        break;
    case SEQ_POMO:
        throw "Sorry! SEQ_POMO is currently not supported";
        break;
    default:
        tree->aln->num_states = 4;
        break;
    }
    
    // add all leaf nodes' name into the alignment
    addLeafNamesToAlignment(tree->aln, tree->root, tree->root);
}

void addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        aln->addSeqName(node->name);
    }
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        addLeafNamesToAlignment(aln, (*it)->node, node);
    }
}

void initializeModel(Params params, IQTree *tree)
{
    tree->aln->model_name = params.model_name;
    ModelsBlock *models_block = readModelsDefinition(params);
    tree->setParams(&params);
    
    tree->initializeModel(params, tree->aln->model_name, models_block);
}

void generateSingleDatasetFromSingleTree(Params params, IQTree *tree, string output_filepath)
{
    // get the ancestral sequence from file or generate it randomly
    IntVector ancestral_sequence = getAncestralSequence(params, tree);
    
    // set ancestral sequence to the root node
    tree->MTree::root->sequence = ancestral_sequence;
    
    // simulate the sequence for each node in the tree by DFS
    simulateSeqsForTree(params.alisim_sequence_length, tree);
    
    // write output to file
    writeSequencesToFile(output_filepath, tree, params.alisim_sequence_length);
}

IntVector getAncestralSequence(Params params, IQTree *tree)
{
    // retrieve the ancestral sequence from input file if its position is specified in the input parameter
    if (params.alisim_ancestral_sequence >= 0)
        return retrieveAncestralSequenceFromInputFile(params.alisim_ancestral_sequence, tree);
    
    // otherwise, randomly generate the sequence
    return generateRandomSequence(params.alisim_sequence_length, tree);
}

IntVector retrieveAncestralSequenceFromInputFile(int sequence_position, IQTree *tree)
{
    IntVector sequence;
    
    // FAKE -> fixed the input sequence instead of retrieved it from the input file
    string sequence_str = "GGAGAGTGTCCTGACCTGGAAGGAATACCTGTAAAGGGGGCGCCATTTATAAAACTACATAGATGGCTCAAAACTAGGACCATAATGCCGGTCCTCAAGG";
    
    sequence.resize(sequence_str.length());
    // convert the input sequence into (numerical states) sequence
    for (int i = 0; i < sequence_str.length(); i++)
        sequence[i] = tree->aln->convertState(sequence_str[i]);
        
    return sequence;
}

IntVector generateRandomSequence(int sequence_length, IQTree *tree)
{
    // initialize sequence
    IntVector sequence;
    sequence.resize(sequence_length);
    
    // get max_num_bases
    int max_num_states = tree->aln->getMaxNumStates();
    
    // if the Frequency Type is FREQ_EQUAL -> randomly generate each site in the sequence follows the normal distribution
    if (tree->getModel()->getFreqType() == FREQ_EQUAL)
    {
        for (int i = 0; i < sequence_length; i++)
            sequence[i] =  random_int(max_num_states);
    }
    else // otherwise, randomly generate each site in the sequence follows the base frequencies defined by the user
    {
        // get the base frequencies
        double *state_freq = new double[max_num_states];
        tree->getModel()->getStateFrequency(state_freq);
        
        // randomly generate each site in the sequence follows the base frequencies defined by the user
        for (int i = 0; i < sequence_length; i++)
        sequence[i] =  getRandomItemWithProbabilityMatrix(state_freq, 0, max_num_states);
        
        // delete state_freq
        delete []  state_freq;
    }
    
    return sequence;
}

int getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items)
{
    // generate a random number
    double random_number = random_double();
    
    // select the current state, considering the random_number, and the probability_matrix
    double accummulated_probability = 0;
    for (int i = 0; i < num_items; i++)
    {
        accummulated_probability += probability_maxtrix[starting_index+i];
        if (random_number <= accummulated_probability)
            return i;
    }
}

void simulateSeqsForTree(int sequence_length, IQTree *tree)
{
    string rate_name = tree->getRateName();
    double invariant_proportion = tree->getRate()->getPInvar();
    
    // simulate Sequences
    // case 1: without rate heterogeneity (Gamma/FreeRate model)
    if ((rate_name.find("+G") == std::string::npos) && (rate_name.find("+R") == std::string::npos))
    {
        simulateSeqsWithoutRH(sequence_length, tree->getModel(), tree->aln->getMaxNumStates(), tree->MTree::root, tree->MTree::root, invariant_proportion);
    }
    else // case 2: with rate heterogeneity (Gamma/FreeRate model)
    {
        simulateSeqsWithRH(sequence_length, tree->getModel(), tree->getRate(), tree->aln->getMaxNumStates(), tree->MTree::root, tree->MTree::root, invariant_proportion);
    }
}

void simulateSeqsWithoutRH(int sequence_length, ModelSubst *model, int max_num_states, Node *node, Node *dad, double invariant_proportion)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // compute the transition probability matrix
        double *trans_matrix = new double[max_num_states*max_num_states];
        model->computeTransMatrix((*it)->length, trans_matrix);
        
        // estimate the sequence for the current neighbor
        // (case 1: without rate heterogeneity (Gamma/FreeRate model)
        // case 1.1: with invariant sites
        if (invariant_proportion > 0)
            (*it)->node->sequence = estimateSeqWithoutRHWithIS(trans_matrix, max_num_states, sequence_length, node->sequence, invariant_proportion);
        else // case 1.2: without invariant sites
            (*it)->node->sequence = estimateSeqWithoutRHWithoutIS(trans_matrix, max_num_states, sequence_length, node->sequence);
        
        // delete trans_matrix array
        delete[] trans_matrix;
        
        // browse 1-step deeper to the neighbor node
        simulateSeqsWithoutRH(sequence_length, model, max_num_states, (*it)->node, node, invariant_proportion);
    }
}

IntVector estimateSeqWithoutRHWithoutIS(double *trans_matrix, int max_num_states, int sequence_length, IntVector dad_sequence)
{
    IntVector sequence;
    sequence.resize(sequence_length);
    
    for (int i = 0; i < sequence_length; i++)
    {
        // iteratively select the state for each site for the current node, considering it's dad states, and the transition_probability_matrix
        int starting_index = dad_sequence[i]*max_num_states;
        sequence[i] = getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states);
    }
    return sequence;
}

IntVector estimateSeqWithoutRHWithIS(double *trans_matrix, int max_num_states, int sequence_length, IntVector dad_sequence, double invariant_proportion)
{
    IntVector sequence;
    sequence.resize(sequence_length);
    
    for (int i = 0; i < sequence_length; i++)
    {
        // generate a random number
        double random_number = random_double();
        
        // if this site is invariant -> preserve the dad's state
        if (random_number <= invariant_proportion)
            sequence[i] = dad_sequence[i];
        else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
        {
            int starting_index = dad_sequence[i]*max_num_states;
            sequence[i] = getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states);
        }
        
    }
    return sequence;
}

void simulateSeqsWithRH(int sequence_length, ModelSubst *model, RateHeterogeneity *rate_heterogeneity, int max_num_states, Node *node, Node *dad, double invariant_proportion)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // initialize transition probability matrix
        double *trans_matrix = new double[max_num_states*max_num_states];
        
        // initialize the probability array of rate categories
        double *category_probability_matrix = new double[rate_heterogeneity->getNDiscreteRate()];
        for (int i = 0; i < rate_heterogeneity->getNDiscreteRate(); i++)
        category_probability_matrix[i] = rate_heterogeneity->getProp(i);
        
        // estimate the sequence for the current neighbor
        // (case 2: with rate heterogeneity (Gamma/FreeRate model)
        // case 2.1: with invariant sites
        if (invariant_proportion > 0)
            (*it)->node->sequence = estimateSeqWithRHWithIS(model, rate_heterogeneity, category_probability_matrix, trans_matrix, max_num_states, sequence_length, (*it)->length, node->sequence, invariant_proportion);
        else // case 2.2: without invariant sites
            (*it)->node->sequence = estimateSeqWithRHWithoutIS(model, rate_heterogeneity, category_probability_matrix, trans_matrix, max_num_states, sequence_length, (*it)->length, node->sequence);
        
        // delete trans_matrix array, and the probability array of rate categories
        delete[] trans_matrix;
        delete[] category_probability_matrix;
        
        // browse 1-step deeper to the neighbor node
        simulateSeqsWithoutRH(sequence_length, model, max_num_states, (*it)->node, node, invariant_proportion);
    }
}

IntVector estimateSeqWithRHWithoutIS(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, int sequence_length, double branch_length, IntVector dad_sequence)
{
    IntVector sequence;
    sequence.resize(sequence_length);
    
    for (int i = 0; i < sequence_length; i++)
    {
        sequence[i] = estimateStateWithRH(model, rate_heterogeneity, category_probability_matrix, trans_matrix, max_num_states, branch_length, dad_sequence[i]);
    }
    return sequence;
}

IntVector estimateSeqWithRHWithIS(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, int sequence_length, double branch_length, IntVector dad_sequence, double invariant_proportion)
{
    IntVector sequence;
    sequence.resize(sequence_length);
    
    for (int i = 0; i < sequence_length; i++)
    {
        // generate a random number
        double random_number = random_double();
        
        // if this site is invariant -> preserve the dad's state
        if (random_number <= invariant_proportion)
            sequence[i] = dad_sequence[i];
        else // otherwise, randomly select the state, considering it's dad states, and the transition_probability_matrix
            sequence[i] = estimateStateWithRH(model, rate_heterogeneity, category_probability_matrix, trans_matrix, max_num_states, branch_length, dad_sequence[i]);
    }
    return sequence;
}

int estimateStateWithRH(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, double branch_length, int dad_state)
{
    // randomly select a rate from the set of rate categories, considering its probability array.
    int rate_category = getRandomItemWithProbabilityMatrix(category_probability_matrix, 0, rate_heterogeneity->getNDiscreteRate());
    double rate = rate_heterogeneity->getRate(rate_category);
    
    // compute the transition matrix
    model->computeTransMatrix(branch_length*rate, trans_matrix);
    
    // iteratively select the state, considering it's dad states, and the transition_probability_matrix
    int starting_index = dad_state*max_num_states;
    return getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states);
}

void writeSequencesToFile(string file_path, IQTree *tree, int sequence_length)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file_path.c_str());
        
        // write the first line <#taxa> <length_of_sequence>
        out <<(tree->leafNum) <<" "<<sequence_length << endl;
        
        // write senquences of leaf nodes to file
        writeASequenceToFile(tree->aln, out, tree->root, tree->root);
        
        // close the file
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_path);
    }
}

void writeASequenceToFile(Alignment *aln, ofstream &out, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        out <<node->name <<" "<<convertEncodedSequenceToReadableSequence(aln, node->sequence) << endl;
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFile(aln, out, (*it)->node, node);
    }
}

string convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence)
{
    string output_sequence = "";

    for (int state : sequence)
        output_sequence = output_sequence + aln->convertStateBackStr(state);
        
    return output_sequence;
    
}
