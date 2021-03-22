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
    simulateSeqsForTree(params, tree);
    
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

void initializeDiscreteRates(double *site_specific_rates, RateHeterogeneity *rate_heterogeneity, int sequence_length)
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

void initializeContinuousGammaRates(double *site_specific_rates, default_random_engine generator, gamma_distribution<double> distribution, int sequence_length, double invariant_proportion)
{
    double sum_rate = 0;
    
    for (int i = 0; i < sequence_length; i++)
    {
        // if this site is invariant -> its rate is zero
        if (random_double() <= invariant_proportion)
            site_specific_rates[i] = 0;
        else
            site_specific_rates[i] = distribution(generator);
        
        // update sum_rate
        sum_rate += site_specific_rates[i];
    }
    
    // compute mean_rate
    double mean_rate = sum_rate/sequence_length;
    
    // normalize the rates
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] /= mean_rate;
    }
}

void initializeContinuousGammaRates(double *site_specific_rates, default_random_engine generator, gamma_distribution<double> distribution, int sequence_length)
{
    double sum_rate = 0;
    
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] = distribution(generator);
        
        // update sum_rate
        sum_rate += site_specific_rates[i];
    }
    
    // compute mean_rate
    double mean_rate = sum_rate/sequence_length;
    
    // normalize the rates
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] /= mean_rate;
    }
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
    
    // if not found, return -1
    return -1;
}

void convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns)
{
    for (int r = 0; r < num_rows; r++)
    {
        for (int c = 1; c < num_columns; c++)
        probability_maxtrix[r*num_rows+c] = probability_maxtrix[r*num_rows+c] + probability_maxtrix[r*num_rows+c-1];
    }
            
}

int getRandomItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, int starting_index, int num_columns)
{
    // generate a random number
    double random_number = random_double();
    
    return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, starting_index, starting_index+(num_columns-1), starting_index)-starting_index;
}

int binarysearchItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, double random_number, int start, int end, int first)
{
    // check search range
    if (start > end)
        return -1; // return -1 ~ not found
    
    // compute the center index
    int center = (start + end)/2;
    
    // if item is found at the center index -> return result
    if ((random_number <= accumulated_probability_maxtrix[center])
        && ((center == first)
            || (random_number > accumulated_probability_maxtrix[center - 1])))
        return center;
    
    // otherwise, search in the left/right side.
    if (random_number <= accumulated_probability_maxtrix[center])
        return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, start, center - 1, first);
    else
        return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, center + 1, end, first);
}

void simulateSeqsForTree(Params params, IQTree *tree)
{
    // get variables
    int sequence_length = params.alisim_sequence_length;
    string rate_name = tree->getRateName();
    double invariant_proportion = tree->getRate()->getPInvar();
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // simulate Sequences
    // case 1: without rate heterogeneity
    if (rate_name.empty())
    {
        simulateSeqsWithoutRH(sequence_length, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
    }
    // case 2: with rate heterogeneity
    else if((rate_name.find("+G") != std::string::npos) || (rate_name.find("+R") != std::string::npos))
    {
        RateHeterogeneity *rate_heterogeneity = tree->getRate();
        
        // initialize site-specific rates
        double *site_specific_rates = new double[sequence_length];
        // initalize rates based on continuous gamma distribution
        if ((rate_name.find("+G") != std::string::npos) && params.alisim_continuous_gamma)
        {
            // initialize gamma distribution
            default_random_engine generator;
            generator.seed(params.ran_seed);
            gamma_distribution<double> distribution(rate_heterogeneity->getGammaShape(),rate_heterogeneity->getGammaShape());
            
            if (invariant_proportion > 0)
                initializeContinuousGammaRates(site_specific_rates, generator, distribution, sequence_length, invariant_proportion);
            else initializeContinuousGammaRates(site_specific_rates, generator, distribution, rate_heterogeneity->getGammaShape(), sequence_length);
        }
        // initalize rates based on discrete distribution (gamma/freerate)
        else
        {
            initializeDiscreteRates(site_specific_rates, rate_heterogeneity, sequence_length);
        }
        
        // simulate sequences with rate heterogeneity
        // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
        if (invariant_proportion > 0)
        {
            simulateSeqsWithRateHeterogeneityWithInvariantSites(sequence_length, model, trans_matrix, site_specific_rates, max_num_states, tree->MTree::root, tree->MTree::root);
        }
        // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
        else
        {
            simulateSeqsWithRateHeterogeneityWithoutInvariantSites(sequence_length, model, trans_matrix, site_specific_rates, max_num_states, tree->MTree::root, tree->MTree::root);
        }
            
        // delete the site-specific rates
        delete[] site_specific_rates;
    }
    // case 2.3: without gamma/freerate model with only invariant sites
    else if (rate_name.find("+I") != std::string::npos)
    {
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
        simulateSeqsWithOnlyInvariantSites(sequence_length, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root, site_specific_rates);
        
        // delete the site-specific rates
        delete[] site_specific_rates;
    }
    
    // delete trans_matrix array
    delete[] trans_matrix;
}

// case 1: without rate heterogeneity
void simulateSeqsWithoutRH(int sequence_length, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad)
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
            // iteratively select the state for each site of the child node, considering it's dad states, and the transition_probability_matrix
            int starting_index = node->sequence[i]*max_num_states;
            (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbabilityMatrix(trans_matrix, starting_index, max_num_states);
        }
        
        // browse 1-step deeper to the neighbor node
        simulateSeqsWithoutRH(sequence_length, model, trans_matrix, max_num_states, (*it)->node, node);
    }
}

// case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
void simulateSeqsWithRateHeterogeneityWithInvariantSites(int sequence_length, ModelSubst *model, double *trans_matrix, double *site_specific_rates, int max_num_states, Node *node, Node *dad)
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
        simulateSeqsWithRateHeterogeneityWithInvariantSites(sequence_length, model, trans_matrix, site_specific_rates, max_num_states, (*it)->node, node);
    }
}

// case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
void simulateSeqsWithRateHeterogeneityWithoutInvariantSites(int sequence_length, ModelSubst *model, double *trans_matrix, double *site_specific_rates, int max_num_states, Node *node, Node *dad)
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
        simulateSeqsWithRateHeterogeneityWithoutInvariantSites(sequence_length, model, trans_matrix, site_specific_rates, max_num_states, (*it)->node, node);
    }
}

// case 2.3: with only invariant sites
void simulateSeqsWithOnlyInvariantSites(int sequence_length, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad, double *site_specific_rates)
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
        simulateSeqsWithOnlyInvariantSites(sequence_length, model, trans_matrix, max_num_states, (*it)->node, node, site_specific_rates);
    }
}

int estimateStateWithRH(ModelSubst *model, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state)
{
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
