//
//  alisimulator.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulator.h"

AliSimulator::AliSimulator(Params *input_params, int expected_number_sites, double new_partition_rate)
{
    params = input_params;
    AliSimulator::initializeIQTreeFromTreeFile();
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    if (expected_number_sites == -1)
        expected_num_sites = params->alisim_sequence_length/num_sites_per_state*params->alisim_length_ratio;
    else
        expected_num_sites = expected_number_sites;
    partition_rate = new_partition_rate;
}

/**
    constructor
*/
AliSimulator::AliSimulator(Params *input_params, IQTree *iq_tree, int expected_number_sites, double new_partition_rate)
{
    params = input_params;
    tree = iq_tree;
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    if (expected_number_sites == -1)
        expected_num_sites = params->alisim_sequence_length/num_sites_per_state*params->alisim_length_ratio;
    else
        expected_num_sites = expected_number_sites;
    partition_rate = new_partition_rate;
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
*  initialize an IQTree instance from input file
*/
void AliSimulator::initializeIQTreeFromTreeFile()
{
    // handle the case with partition models
    if (params->partition_file) {
        // initilize partition alignments
        Alignment *aln;
        if (params->partition_type == TOPO_UNLINKED)
            aln = new SuperAlignmentUnlinked(*params);
        else
            aln = new SuperAlignment(*params);
        
        // initialize a super tree
        if (params->partition_type == TOPO_UNLINKED) {
            tree = new PhyloSuperTreeUnlinked((SuperAlignment*) aln);
        } else if(params->partition_type != BRLEN_OPTIMIZE){
            // initialize supertree - Proportional Edges case
            tree = new PhyloSuperTreePlen((SuperAlignment*) aln, params->partition_type);
        } else {
            // initialize supertree stuff if user specifies partition file with -sp option
            tree = new PhyloSuperTree((SuperAlignment*) aln);
        }
        tree->setParams(params);
        bool is_rooted = false;
        tree->readTree(params->user_file, is_rooted);
        
        // further initialize super_tree/alignments
        for (int i = 0; i < ((PhyloSuperTree*) tree)->size(); i++)
        {
            // load phylotrees
            IQTree *current_tree = (IQTree *) ((PhyloSuperTree*) tree)->at(i);
            bool is_rooted = false;
            current_tree->readTree(params->user_file, is_rooted, i);
            
            // update the alignment for the current partition
            initializeAlignment(current_tree, current_tree->aln->model_name);
            
            // extract num_sites from partition
            IntVector siteIDs;
            current_tree->aln->extractSiteID(current_tree->aln, current_tree->aln->position_spec.c_str(), siteIDs, -1, true);
            current_tree->aln->setExpectedNumSites(siteIDs.size());
            
            // initialize the model for the current partition
            initializeModel(current_tree, current_tree->aln->model_name);
            
            // if a Heterotachy model is used -> re-read the PhyloTreeMixlen from file
            if (current_tree->getRate()->isHeterotachy())
            {
                // initialize a new PhyloTreeMixlen
                IQTree* new_tree = new PhyloTreeMixlen(current_tree->aln, current_tree->getRate()->getNRate());
                
                // delete the old tree
                delete current_tree;
                
                // set the new PhyloTreeMixlen to the new tree
                current_tree = new_tree;
                
                // re-load the tree/branch-lengths from the file
                current_tree->IQTree::readTree(params->user_file, is_rooted, i);
                
                // re-initialize the model
                initializeModel(current_tree, current_tree->aln->model_name);
            }
            
            // set partition rate
            if (params->partition_type == BRLEN_SCALE && params->alisim_partition_rates.size() > i)
                ((PhyloSuperTree*) tree)->part_info[i].part_rate = params->alisim_partition_rates[i];
        }
    }
    // other cases without partition models
    else
    {
        // initialize tree
        tree = new IQTree();
        bool is_rooted = false;
        tree->readTree(params->user_file, is_rooted);
        tree->setParams(params);
        
        // initialize alignment
        tree->aln = new Alignment();
        initializeAlignment(tree, params->model_name);
        
        // inittialize model
        initializeModel(tree, params->model_name);

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
            initializeModel(tree, params->model_name);
        }
    }
}


/**
*  initialize an Alignment instance for IQTree
*/
void AliSimulator::initializeAlignment(IQTree *tree, string model_fullname)
{
    if (!params->sequence_type)
    {
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
void AliSimulator::initializeModel(IQTree *tree, string model_name)
{
    tree->aln->model_name = model_name;
    tree->aln->computeUnknownState();
    ModelsBlock *models_block = readModelsDefinition(*params);
    
    tree->IQTree::initializeModel(*params, tree->aln->model_name, models_block);
}

/**
    remove all constant sites (in case with +ASC)
*/
void AliSimulator::removeConstantSites(){
    // dummy variables
    int num_variant_states = -1;
    IntVector variant_state_mask;
    
    // create a variant state mask
    createVariantStateMask(variant_state_mask, num_variant_states, expected_num_sites, tree->root, tree->root);
    
    // return error if num_variant_states is less than the expected_num_variant_states
    if (num_variant_states < expected_num_sites){
        outError("Unfortunately, after removing constant sites, the number of variant sites is less than the expected sequence length. Please use --length-ratio <LENGTH_RATIO> to generate more abundant sites and try again. The current <LENGTH_RATIO> is "+ convertDoubleToString(params->alisim_length_ratio));
    }
    
    // get only variant sites for leaves
    getOnlyVariantSites(variant_state_mask, tree->root, tree->root);
}

/**
    only get variant sites
*/
void AliSimulator::getOnlyVariantSites(IntVector variant_state_mask, Node *node, Node *dad){
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // dummy sequence
        IntVector variant_sites;
        
        // initialize the number of variant sites
        int num_variant_states = 0;
        
        // browse sites one by one
        for (int i = 0; i < node->sequence.size(); i++)
            // only get variant sites
            if (variant_state_mask[i] == -1)
            {
                // get the variant site
                variant_sites.push_back(node->sequence[i]);
                num_variant_states++;
                
                // stop checking further states if num_variant_states has exceeded the expected_num_variant_states
                if (num_variant_states >= expected_num_sites)
                    break;
            }
        
        // replace the sequence of the Leaf by variant sites
        node->sequence.clear();
        node->sequence = variant_sites;
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        getOnlyVariantSites(variant_state_mask, (*it)->node, node);
    }
}

/**
*  generate the current partition of an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void AliSimulator::generatePartitionAlignment(IntVector ancestral_sequence)
{
    // if the ancestral sequence is not specified, randomly generate the sequence
    if (ancestral_sequence.size() == 0)
    {
        int sequence_length = expected_num_sites*params->alisim_length_ratio;
        tree->MTree::root->sequence = generateRandomSequence(sequence_length);
    }
    // otherwise, using the ancestral sequence + abundant sites
    else
    {
        // set the ancestral sequence to the root node
        tree->MTree::root->sequence = ancestral_sequence;
        
        // add abundant_sites
        int num_abundant_sites = expected_num_sites * params->alisim_length_ratio - ancestral_sequence.size();
        if (num_abundant_sites > 0)
        {
            IntVector abundant_sites = generateRandomSequence(num_abundant_sites);
            for (int site:abundant_sites)
                tree->MTree::root->sequence.push_back(site);
        }
    }
    
    // validate the sequence length (in case of codon)
    validataSeqLengthCodon();
    
    // simulate the sequence for each node in the tree by DFS
    simulateSeqsForTree();
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
             || tree->aln->seq_type == SEQ_CODON
             || (tree->getModel()->getFreqType() == FREQ_EMPIRICAL && params->aln_file))
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
    int sequence_length = expected_num_sites*params->alisim_length_ratio;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize trans_matrix
    double *trans_matrix = new double[max_num_states*max_num_states];
    
    // simulate Sequences
    simulateSeqs(sequence_length, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root);
    
    // delete trans_matrix array
    delete[] trans_matrix;
    
    // removing constant states if it's necessary
    if (params->alisim_length_ratio > 1)
        removeConstantSites();
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
        model->computeTransMatrix(partition_rate*(*it)->length, trans_matrix);
        
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
*  validate sequence length of codon
*
*/
void AliSimulator::validataSeqLengthCodon()
{
    if (tree->aln->seq_type == SEQ_CODON && (params->alisim_sequence_length%3))
        outError("Sequence length of Codon must be divisible by 3. Please check & try again!");
}

/**
*  update the expected_num_sites due to the change of the sequence_length
*/
void AliSimulator::refreshExpectedNumSites(){
    expected_num_sites = params->alisim_sequence_length/num_sites_per_state*params->alisim_length_ratio;
}
