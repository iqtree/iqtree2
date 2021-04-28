//
//  alisimulator.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulator.h"

AliSimulator::AliSimulator(Params *input_params)
{
    params = input_params;
    AliSimulator::initializeIQTreeFromTreeFile();
}

/**
    constructor
*/
AliSimulator::AliSimulator(Params *input_params, IQTree *iq_tree)
{
    params = input_params;
    tree = iq_tree;
}

AliSimulator::~AliSimulator()
{
    if (!tree || !(tree->aln)) return;
    
    // delete aln
    delete tree->aln;
    
    // delete tree
    delete tree;
}

/**
*  show all input parameters for AliSim
*/
void AliSimulator::showParameters()
{
    cout << " - Tree filepath: " << params->user_file <<"\n";
    cout << " - Length of output sequences: " << params->alisim_sequence_length <<"\n";
    if (!params->model_name.empty())
        cout << " - Model: " << params->model_name <<"\n";
    cout << " - Number of output datasets: " << params->alisim_dataset_num<<"\n";
    if (params->alisim_ancestral_sequence_name.length() > 0)
        cout << " - Ancestral sequence position: " << params->alisim_dataset_num <<"\n";
}

/**
*  initialize an IQTree instance from input file
*/
void AliSimulator::initializeIQTreeFromTreeFile()
{
    tree = new IQTree();
    bool is_rooted = false;
    tree->readTree(params->user_file, is_rooted);
    initializeAlignment();
    initializeModel();

    // if a Heterotachy model is used -> re-read the PhyloTreeMixlen from file
    if (tree->getRate()->isHeterotachy())
    {
        // initialize a new PhyloTreeMixlen
        IQTree* new_tree = new PhyloTreeMixlen(tree->aln, tree->getRate()->getNRate());
        
        // delete the old tree
        delete tree;
        
        // set the new PhyloTreeMixlen to the new tree
        tree = new_tree;
        
        // re-load the tree/branch-lengths from the file
        tree->IQTree::readTree(params->user_file, is_rooted);
        
        // re-initialize the model
        initializeModel();
    }
}


/**
*  initialize an Alignment instance for IQTree
*/
void AliSimulator::initializeAlignment()
{
    tree->aln = new Alignment();
    
    if (!params->sequence_type)
    {
        // set the seq_type and the maximum number of bases based on the Seq_type
        string model_fullname = params->model_name;
        // if a mixture model is used -> extract the name of the first model component for SeqType detection
        string KEYWORD = "MIX";
        string delimiter = ",";
        if ((model_fullname.length() > KEYWORD.length())
            && (!model_fullname.substr(0, KEYWORD.length()).compare(KEYWORD)))
        {
            // only get the model name, removing additional params (+G,+F,*G,*F,etc)
            model_fullname = model_fullname.substr(0, model_fullname.find("+"));
            model_fullname = model_fullname.substr(0, model_fullname.find("*"));
            
            // validate the input
            if ((model_fullname[KEYWORD.length()]!='{')
                ||(model_fullname[model_fullname.length()-1]!='}')
                ||(model_fullname.find(delimiter) == string::npos))
                outError("Use -m MIX{m1,...,mK} to define a mixture model.");
            
            // remove "MIX{"
            model_fullname.erase(0, KEYWORD.length() + 1);
            
            // get the first model name
            model_fullname = model_fullname.substr(0, model_fullname.find(delimiter));
            
            // remove the weight (if any)
            model_fullname = model_fullname.substr(0, model_fullname.find(":"));
        }
        string model_familyname_with_params = model_fullname.substr(0, model_fullname.find("+"));
        model_familyname_with_params = model_familyname_with_params.substr(0, model_fullname.find("*"));
        string model_familyname = model_familyname_with_params.substr(0, model_familyname_with_params.find("{"));
        detectSeqType(model_familyname.c_str(), tree->aln->seq_type);
        string seq_type_name = convertSeqTypeToSeqTypeName(tree->aln->seq_type);
        params->sequence_type = strcpy(new char[seq_type_name.length() + 1], seq_type_name.c_str());
    }
    else
        tree->aln->seq_type = tree->aln->getSeqType(params->sequence_type);
    
    if (tree->aln->seq_type == SEQ_UNKNOWN)
        outError("Could not detect SequenceType from Model Name. Please check your Model Name or specify the SequenceType by --seqtype <SEQ_TYPE_STR> where <SEQ_TYPE_STR> is BIN, DNA, AA, NT2AA, CODON, or MORPH.");
    
    switch (tree->aln->seq_type) {
    case SEQ_BINARY:
        tree->aln->num_states = 2;
        break;
    case SEQ_PROTEIN:
        tree->aln->num_states = 20;
        break;
    case SEQ_MORPH:
        tree->aln->num_states = params->alisim_num_states_morph;
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
    
    // init Codon (if neccessary)
    if (tree->aln->seq_type == SEQ_CODON)
        tree->aln->initCodon(&params->sequence_type[5]);
}

/**
*  iteratively add name of all leaf nodes into the alignment instance
*/
void AliSimulator::addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        aln->addSeqName(node->name);
    }
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        addLeafNamesToAlignment(aln, (*it)->node, node);
    }
}

/**
*  initialize a Model instance for IQTree
*/
void AliSimulator::initializeModel()
{
    tree->aln->model_name = params->model_name;
    tree->aln->computeUnknownState();
    ModelsBlock *models_block = readModelsDefinition(*params);
    tree->setParams(params);
    
    tree->initializeModel(*params, tree->aln->model_name, models_block);
}

/**
*  generate an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void AliSimulator::generateSingleDatasetFromSingleTree(string output_filepath, IntVector ancestral_sequence)
{
    // set ancestral sequence to the root node
    tree->MTree::root->sequence = ancestral_sequence;
    
    // validate the sequence length (in case of codon)
    validataSeqLengthCodon();
    
    // simulate the sequence for each node in the tree by DFS
    simulateSeqsForTree();
    
    // if constant states are needed to be removed -> firstly create the variant state mask
    int num_variant_states = -1;
    int expected_num_variant_states = params->alisim_sequence_length/params->alisim_sites_per_state;
    IntVector variant_state_mask;
    if (params->alisim_abundant_state_rate > 1)
    {
        createVariantStateMask(variant_state_mask, num_variant_states, expected_num_variant_states, tree->root, tree->root);
        // return error if num_variant_states is less than the expected_num_variant_states
        if (num_variant_states < expected_num_variant_states){
            outError("Unfortunately, after removing constant sites, the number of variant sites is less than the expected sequence length. Please use --abundant-site-rate <ABUNDANT_SITE_RATE> to generate more abundant sites and try again. The current <ABUNDANT_SITE_RATE> is "+ convertDoubleToString(params->alisim_abundant_state_rate));
        }
    }
    
    // write output to file
    writeSequencesToFile(output_filepath, variant_state_mask);
}

/**
    create mask for variant sites
*/
void AliSimulator::createVariantStateMask(IntVector &variant_state_mask, int &num_variant_states, int expected_num_variant_states, Node *node, Node *dad){
    // no need to check the further sites if num_variant_states has exceeded the expected_num_variant_states
    if (num_variant_states >= expected_num_variant_states)
        return;
    
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // initialize the mask (all sites are assumed to be constant)
        if (num_variant_states == -1)
        {
            num_variant_states = 0;
            for (int i = 0; i < node->sequence.size(); i++)
                variant_state_mask.push_back(node->sequence[i]);
        }
        // otherwise, check state by state to update the mask
        else
        {
            for (int i = 0; i < node->sequence.size(); i++)
                if (variant_state_mask[i] != -1 && variant_state_mask[i] != node->sequence[i])
                {
                    // if the current state is changed -> increase num_variant_states, and disable that state
                    variant_state_mask[i] = -1;
                    num_variant_states++;
                    
                    // stop checking further states if num_variant_states has exceeded the expected_num_variant_states
                    if (num_variant_states >= expected_num_variant_states)
                        break;
                }
        }
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        createVariantStateMask(variant_state_mask, num_variant_states, expected_num_variant_states, (*it)->node, node);
    }
}

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void AliSimulator::generateMultipleAlignmentsFromSingleTree()
{
    // Get/Generate the ancestral sequence
    IntVector ancestral_sequence;
    // retrieve the ancestral sequence from input file if its position is specified in the input parameter
    if (params->alisim_ancestral_sequence_name.length() > 0)
        ancestral_sequence = retrieveAncestralSequenceFromInputFile(params->alisim_ancestral_sequence_aln_filepath, params->alisim_ancestral_sequence_name);
    
    // iteratively generate multiple datasets for each tree
    for (int i = 0; i < params->alisim_dataset_num; i++)
    {
        // initialize site specific model index based on its weights (in the mixture model)
        intializeSiteSpecificModelIndex();
        
        // if the ancestral sequence is not specified, randomly generate the sequence
        if (params->alisim_ancestral_sequence_name.length() == 0)
            ancestral_sequence = generateRandomSequence((int)(params->alisim_sequence_length/params->alisim_sites_per_state*params->alisim_abundant_state_rate));
        
        // initialize output_filepath
        std::string output_filepath(params->user_file);
        output_filepath = output_filepath.substr(0, output_filepath.find_last_of("/\\") + 1);
        output_filepath = output_filepath
        +params->alisim_output_filename
        +"_"+convertIntToString(i)+".phy";
        
        generateSingleDatasetFromSingleTree(output_filepath, ancestral_sequence);
    }
}

/**
*  retrieve the ancestral sequence for the root node from an input file
*/
IntVector AliSimulator:: retrieveAncestralSequenceFromInputFile(char *aln_filepath, string sequence_name)
{
    IntVector sequence;
    
    // read sequences from the input file
    Alignment *aln = new Alignment();
    StrVector sequences;
    int nseq = 0, nsite = 0;
    aln->extractSequences(aln_filepath, params->sequence_type, sequences, nseq, nsite);
    StrVector seq_names = aln->getSeqNames();
    
    // delete aln
    delete aln;
    
    // overwrite the output sequence_length
    params->alisim_sequence_length = nsite;
    cout<<" Sequence length is now set equally to the length of ancestral sequence."<<endl;
    
    string sequence_str = "";
    for (int i = 0; i < seq_names.size(); i++)
        if (!sequence_name.compare(seq_names[i]))
        {
            sequence_str = sequences[i];
            break;
        }
    if (sequence_str.length() == 0)
        outError("Sequence name could not be found in the input alignment file.");
    
    // get Max number of states
    int max_num_states = tree->aln->getMaxNumStates();
    
    // get the base frequencies
    double *state_freq = new double[max_num_states];
    getStateFrequenciesFromModel(state_freq);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(state_freq, 1, max_num_states);
    
    // convert the input sequence into (numerical states) sequence
    int sequence_length = params->alisim_sequence_length/params->alisim_sites_per_state;
    sequence.resize(sequence_length);
    ostringstream err_str;
    int num_error = 0;
    for (int i = 0; i < sequence_length; i++)
    {
        if (tree->aln->seq_type == SEQ_CODON)
        {
            int site_index = i*params->alisim_sites_per_state;
            sequence[i] = tree->aln->getCodonStateTypeFromSites(tree->aln->convertState(sequence_str[site_index], SEQ_DNA), tree->aln->convertState(sequence_str[site_index+1], SEQ_DNA), tree->aln->convertState(sequence_str[site_index+2], SEQ_DNA), sequence_name, site_index, err_str, num_error);
        }
        else
        {
            sequence[i] = tree->aln->convertState(sequence_str[i]);
        
            // Handle invalid/unknown state
            if (sequence[i] >= max_num_states)
                sequence[i] = getRandomItemWithAccumulatedProbabilityMatrix(state_freq, 0, max_num_states);
        }
    }
    
    // show error from ancestral sequence (if any)
    if (num_error)
        outError(err_str.str());
    
    // randomly generate abundant states (if necessary)
    int num_abundant_states = params->alisim_sequence_length/params->alisim_sites_per_state*params->alisim_abundant_state_rate - sequence_length;
    for (int i = 0; i < num_abundant_states; i++)
    {
        sequence.push_back(getRandomItemWithAccumulatedProbabilityMatrix(state_freq, 0, max_num_states));
    }
    
    // delete state_freq
    delete [] state_freq;
    
    // show warning
    outWarning("Using an ancestral sequence with base frequencies that are not compatible with the specification of the model may lead to unexpected results.");
        
    return sequence;
}

/**
*  randomly generate the ancestral sequence for the root node
*/
IntVector AliSimulator::generateRandomSequence(int sequence_length)
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
        getStateFrequenciesFromModel(state_freq);
        
        // print model's parameters
        tree->getModel()->writeInfo(cout);
        
        // convert the probability matrix into an accumulated probability matrix
        convertProMatrixIntoAccumulatedProMatrix(state_freq, 1, max_num_states);
        
        // randomly generate each site in the sequence follows the base frequencies defined by the user
        for (int i = 0; i < sequence_length; i++)
            sequence[i] =  getRandomItemWithAccumulatedProbabilityMatrix(state_freq, 0, max_num_states);
        
        // delete state_freq
        delete []  state_freq;
    }
    
    return sequence;
}

void AliSimulator::getStateFrequenciesFromModel(double *state_freqs){
    // firstly, initialize state freqs for mixture models (if neccessary)
    intializeStateFreqsMixtureModel();
    
    // if a mixture model is used -> get weighted sum of state_freq across classes
    if (tree->getModel()->isMixture())
    {
        tree->getModel()->getStateFrequency(state_freqs, -1);
    }
    // get user-defined base frequencies (if any)
    else if ((tree->getModel()->getFreqType() == FREQ_USER_DEFINED)
        || (ModelLieMarkov::validModelName(tree->getModel()->getName()))
             || tree->aln->seq_type == SEQ_CODON)
        tree->getModel()->getStateFrequency(state_freqs);
    else // otherwise, randomly generate the base frequencies
    {
        generateRandomBaseFrequencies(state_freqs, tree->aln->getMaxNumStates());
        tree->getModel()->setStateFrequency(state_freqs);
    }
}

/**
*  randomly generate the base frequencies
*/
void AliSimulator::generateRandomBaseFrequencies(double *base_frequencies, int max_num_bases)
{
    double sum = 0;
    
    // randomly generate the frequencies
    for (int i = 0; i < max_num_bases; i++)
    {
        base_frequencies[i] = random_double();
        sum += base_frequencies[i];
    }
    
    // normalize the frequencies so that sum of them is 1
    for (int i = 0; i < max_num_bases; i++)
    base_frequencies[i] /= sum;
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulator::simulateSeqsForTree()
{
    // get variables
    int sequence_length = params->alisim_sequence_length/params->alisim_sites_per_state*params->alisim_abundant_state_rate;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // simulate Sequences
    simulateSeqs(sequence_length, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
    
    // delete trans_matrix array
    delete[] trans_matrix;
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulator::simulateSeqs(int sequence_length, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad)
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
        simulateSeqs(sequence_length, model, trans_matrix, max_num_states, (*it)->node, node);
    }
}

/**
*  get a random item from a set of items with a probability array
*/
int AliSimulator::getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items)
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


/**
*  convert an probability matrix into an accumulated probability matrix
*/
void AliSimulator::convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns)
{
    for (int r = 0; r < num_rows; r++)
    {
        for (int c = 1; c < num_columns; c++)
        probability_maxtrix[r*num_columns+c] = probability_maxtrix[r*num_columns+c] + probability_maxtrix[r*num_columns+c-1];
    }
            
}

/**
*  get a random item from a set of items with an accumulated probability array by binary search
*/
int AliSimulator::getRandomItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, int starting_index, int num_columns)
{
    // generate a random number
    double random_number = random_double();
    
    return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, starting_index, starting_index+(num_columns-1), starting_index)-starting_index;
}

/**
*  binary search an item from a set with accumulated probability array
*/
int AliSimulator::binarysearchItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, double random_number, int start, int end, int first)
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

/**
*  write all sequences of a tree to an output file
*/
void AliSimulator::writeSequencesToFile(string file_path, IntVector variant_state_mask)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file_path.c_str());
        
        // write the first line <#taxa> <length_of_sequence>
        int leaf_num = tree->leafNum - ((tree->root->isLeaf() && tree->root->name == ROOT_NAME)?1:0);
        out <<leaf_num<<" "<<params->alisim_sequence_length<< endl;
        
        // write senquences of leaf nodes to file with/without gaps copied from the input sequence
        if (params->aln_file && !params->alisim_no_copy_gaps)
        {
            // load input sequences (with gaps)
            vector<string> seq_names;
            vector<string> sequences;
            loadSequences(params->aln_file, seq_names, sequences);
            
            // show a warning if the length of input alignment is unequal to that of simulated sequence
            if (sequences.size() > 0  && sequences[0].length() != params->alisim_sequence_length)
                outWarning("The sequence length of the input alignment is unequal to that of that simulated sequences. Thus, only gaps in the first MIN(input_sequence_length, simulated_sequence_length) sites are copied.");
            
            // write simulated sequence with the gaps copied from the input sequence
            writeASequenceToFileWithGaps(tree->aln, seq_names, sequences, out, variant_state_mask, tree->root, tree->root);
        }
        else
        {
            writeASequenceToFile(tree->aln, out, variant_state_mask, tree->root, tree->root);
        }
        
        // close the file
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_path);
    }
}

/**
*  load sequences from input file
*/
void AliSimulator::loadSequences(char *file_path, vector<string> &seq_names, vector<string> &sequences)
{
    // read sequences from the input file
    Alignment *aln = new Alignment();
    int nseq = 0, nsite = 0;
    aln->extractSequences(params->aln_file, params->sequence_type, sequences, nseq, nsite);
    seq_names = aln->getSeqNames();
}

/**
*  write a sequence of a node to an output file
*/
void AliSimulator::writeASequenceToFile(Alignment *aln, ofstream &out, IntVector variant_state_mask, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        out <<node->name <<" "<<convertEncodedSequenceToReadableSequence(aln, node->sequence, variant_state_mask) << endl;
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFile(aln, out, variant_state_mask, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
*/
string AliSimulator::convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence, IntVector variant_state_mask)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // initialize num_variant_states
    int num_variant_states = 0;
    
    // browse state by state
    for (int i = 0; i < sequence.size(); i++)
    {
        if (variant_state_mask.empty()
                || (!variant_state_mask.empty() && variant_state_mask[i] == -1))
            {
                output_sequence = output_sequence + aln->convertStateBackStr(sequence[i]);
                num_variant_states++;
                
                // stop checking abundant states if the first expected variant states are all converted into readable sites
                if (num_variant_states >= params->alisim_sequence_length/params->alisim_sites_per_state)
                    break;
            }
    }
        
    return output_sequence;
    
};

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
void AliSimulator::writeASequenceToFileWithGaps(Alignment *aln, vector<string> seq_names, vector<string> sequences, ofstream &out, IntVector variant_state_mask, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // retrieve the input sequence of the current node
        for (int i = 0; i < sequences.size(); i++)
            if (!seq_names[i].compare(node->name))
            {
                out <<node->name <<" "<<convertEncodedSequenceToReadableSequenceWithGaps(aln, sequences[i], node->sequence, variant_state_mask) << endl;
            }
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFileWithGaps(aln, seq_names, sequences, out, variant_state_mask, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...) with gaps copied from the input sequence
*/
string AliSimulator::convertEncodedSequenceToReadableSequenceWithGaps(Alignment *aln, string input_sequence, IntVector sequence, IntVector variant_state_mask)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // initialize num_variant_states
    int num_variant_states = 0;
    
    // browse state by state
    for (int i = 0; i < sequence.size(); i++)
    {
        int num_sites_per_state = params->alisim_sites_per_state;
        // handle gaps
        if ((i+1)*num_sites_per_state - 1 < input_sequence.length()
            && ((num_sites_per_state == 3
                &&(input_sequence[i*num_sites_per_state] == '-'
                   || input_sequence[i*num_sites_per_state+1] == '-'
                   || input_sequence[i*num_sites_per_state+2] == '-'))
                || input_sequence[i] == '-')){
            // insert gaps in case with protein data
            if (num_sites_per_state == 3)
            {
                output_sequence = output_sequence + input_sequence[i*num_sites_per_state];
                output_sequence = output_sequence + input_sequence[i*num_sites_per_state+1];
                output_sequence = output_sequence + input_sequence[i*num_sites_per_state+2];
                
                // increase the num_variant_states
                num_variant_states++;
            }
            // insert gaps in normal cases
            else
            {
                output_sequence = output_sequence + "-";
                
                // increase the num_variant_states
                num_variant_states++;
            }
        }
        // if it's not a gap
        else if (variant_state_mask.empty()
                || (!variant_state_mask.empty() && variant_state_mask[i] == -1))
            {
                output_sequence = output_sequence + aln->convertStateBackStr(sequence[i]);
                num_variant_states++;
            }
        
        // stop checking abundant states if the first expected variant states are all converted into readable sites
        if (num_variant_states >= params->alisim_sequence_length/params->alisim_sites_per_state)
            break;
    }
        
    return output_sequence;
}

/**
*  validate sequence length of codon
*
*/
void AliSimulator::validataSeqLengthCodon()
{
    if (tree->aln->seq_type == SEQ_CODON && (params->alisim_sequence_length%3))
        outError("Sequence length of Codon must be divisible by 3. Please check & try again!");
}
