//
//  alisimulator.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulator.h"
#include "alisimulatorheterogeneity.h"
#include "alisimulatorheterogeneityinvar.h"
#include "alisimulatorinvar.h"

AliSimulator::AliSimulator(Params *input_params, int expected_number_sites, double new_partition_rate)
{
    params = input_params;
    AliSimulator::initializeIQTreeFromTreeFile();
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    STATE_UNKNOWN = tree->aln->STATE_UNKNOWN;
    max_num_states = tree->aln->getMaxNumStates();
    
    // estimating the appropriate length_ratio in cases models with +ASC
    estimateLengthRatio();
    
    if (expected_number_sites == -1)
        expected_num_sites = params->alisim_sequence_length/num_sites_per_state*length_ratio;
    else
        expected_num_sites = expected_number_sites*length_ratio;
    partition_rate = new_partition_rate;
    
    // check if base frequencies for DNA models are specified correctly
    checkBaseFrequenciesDNAModels(tree, params->model_name);
    
    // extract max length of taxa names
    extractMaxTaxaNameLength();
    
    // innialize set of selected sites for permutation in FunDi model
    if (params->alisim_fundi_taxon_set.size()>0)
        fundi_items = selectAndPermuteSites(params->alisim_fundi_proportion, round(expected_num_sites));
}

/**
    constructor
*/
AliSimulator::AliSimulator(Params *input_params, IQTree *iq_tree, int expected_number_sites, double new_partition_rate)
{
    params = input_params;
    tree = iq_tree;
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    STATE_UNKNOWN = tree->aln->STATE_UNKNOWN;
    max_num_states = tree->aln->getMaxNumStates();
    
    // estimating the appropriate length_ratio in cases models with +ASC
    estimateLengthRatio();
    
    if (expected_number_sites == -1)
        expected_num_sites = params->alisim_sequence_length/num_sites_per_state*length_ratio;
    else
        expected_num_sites = expected_number_sites*length_ratio;
    partition_rate = new_partition_rate;
    
    // extract max length of taxa names
    extractMaxTaxaNameLength();
    
    // innialize set of selected sites for permutation in FunDi model
    if (params->alisim_fundi_taxon_set.size()>0)
        fundi_items = selectAndPermuteSites(params->alisim_fundi_proportion, round(expected_num_sites));
}

AliSimulator::~AliSimulator()
{
    if (!tree) return;
    
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
        if (!params->user_file)
            outError("Please supply a tree file by -t <TREE_FILEPATH>");
        tree->readTree(params->user_file, is_rooted);
        
        // compute super_tree_length
        double super_tree_length = ((PhyloSuperTree*) tree)->treeLength();
        
        // sum of rate*n_sites and total sites (for rate normalization)
        double sum = 0;
        int num_sites = 0;
        
        // further initialize super_tree/alignments
        // recording start_time
        auto start = getRealTime();
        
        int i;
        
        for (i = 0; i < ((PhyloSuperTree*) tree)->size(); i++)
        {
            // -Q (params->partition_type == BRLEN_OPTIMIZE) -> tree_line_index = i; otherwise (-p, -q), tree_line_index = 0 (only a tree)
            int tree_line_index = 0;
            if (params->partition_type == BRLEN_OPTIMIZE)
            {
                tree_line_index = i+1;
                // show information for the first time
                if (i == 0)
                {
                    cout<<" The super tree (combining all taxa in all partitions) has been loaded from the first line of the input tree file."<<endl;
                    cout<<" Loading partition trees one by one. Each tree should be specified in a single line in the input tree file."<<endl;
                }
            }
            
            // load phylotrees
            IQTree *current_tree = (IQTree *) ((PhyloSuperTree*) tree)->at(i);
            bool is_rooted = false;
            current_tree->readTree(params->user_file, is_rooted, tree_line_index);
            
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
                current_tree->IQTree::readTree(params->user_file, is_rooted, tree_line_index);
                
                // re-initialize the model
                initializeModel(current_tree, current_tree->aln->model_name);
            }
            
            // set partition rate
            if (params->partition_type == BRLEN_SCALE)
            {
                double current_tree_length = current_tree->aln->tree_len;
                if (current_tree_length <= 0)
                    outError("Please specify tree length for each partition in the input NEXUS file.");
                else
                    ((PhyloSuperTree*) tree)->part_info[i].part_rate = current_tree_length/super_tree_length;
                
                // update sum of rate*n_sites and num_sites (for rate normalization)
                sum += ((PhyloSuperTree*) tree)->part_info[i].part_rate * current_tree->aln->getNSite();
                if (current_tree->aln->seq_type == SEQ_CODON && ((PhyloSuperTree*) tree)->rescale_codon_brlen)
                    num_sites += 3*current_tree->aln->getNSite();
                else
                    num_sites += current_tree->aln->getNSite();
            }
        }
        
        // show the reloading tree time
        auto end = getRealTime();
        cout<<" - Time spent on Loading trees: "<<end-start<<endl;
        
        // normalizing the partition rates (if necessary)
        if (params->partition_type == BRLEN_SCALE)
        {
            sum /= num_sites;
            sum = 1.0/sum;
            
            // check whether normalization is necessary or not
            double epsilon = 0.0001;
            if (sum > 1 + epsilon || sum < 1 - epsilon)
            {
                // show warning
                outWarning("Partitions' rates are normalized so that sum of (partition_rate*partition_sequence_length) of all partitions is 1.");
                
                // update partitions' rates
                for (int i = 0; i < ((PhyloSuperTree*) tree)->size(); i++)
                    ((PhyloSuperTree*) tree)->part_info[i].part_rate  *= sum;
            }
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
    // intializing seq_type if it's unknown
    if (tree->aln->seq_type == SEQ_UNKNOWN)
    {
        // firstly, intializing seq_type from sequence_type if it's not empty
        if (tree->aln->sequence_type.length()>0)
            tree->aln->seq_type = tree->aln->getSeqType(tree->aln->sequence_type.c_str());
        // otherwise, intializing seq_type from sequence_type (in params) if it's not empty
        else
        {
            if (params->sequence_type)
                tree->aln->seq_type = tree->aln->getSeqType(params->sequence_type);
            // otherwise, detect seq_type model's name
            else
            {
                // if a mixture model is used -> extract the name of the first model component for SeqType detection
                string KEYWORD = "MIX";
                string delimiter = ",";
                if ((model_fullname.length() > KEYWORD.length())
                    && (!model_fullname.substr(0, KEYWORD.length()).compare(KEYWORD)))
                {
                    // detect the position of the close_bracket in MIX{...}
                    size_t close_bracket_pos = 0;
                    int num_open_brackets = 0;
                    for (close_bracket_pos = KEYWORD.length(); close_bracket_pos < model_fullname.length(); close_bracket_pos++)
                    {
                        if (model_fullname[close_bracket_pos] == '{')
                            num_open_brackets++;
                        else if (model_fullname[close_bracket_pos] == '}')
                        {
                            num_open_brackets--;
                            if (num_open_brackets == 0)
                                break;
                        }
                    }
                    // only get the model name inside MIX{...}
                    model_fullname = model_fullname.substr(0, close_bracket_pos+1);
                    
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
            }
            if (tree->aln->seq_type != SEQ_UNKNOWN)
                tree->aln->sequence_type = tree->aln->getSeqTypeStr(tree->aln->seq_type);
        }
    }
    
    if (tree->aln->seq_type == SEQ_UNKNOWN)
        outError("Could not detect SequenceType from Model Name. Please check your Model Name or specify the SequenceType by --seqtype <SEQ_TYPE_STR> where <SEQ_TYPE_STR> is BIN, DNA, AA, NT2AA, CODON, or MORPH.");
    
    switch (tree->aln->seq_type) {
    case SEQ_BINARY:
        tree->aln->num_states = 2;
        break;
    case SEQ_DNA:
        tree->aln->num_states = 4;
        break;
    case SEQ_PROTEIN:
        tree->aln->num_states = 20;
        break;
    case SEQ_MORPH:
            // only set num_state if it has not yet set (noting that num_states of Morph could be set in partition file)
            if (tree->aln->num_states == 0)
                tree->aln->num_states = params->alisim_num_states_morph;
        break;
    case SEQ_POMO:
        throw "Sorry! SEQ_POMO is currently not supported";
        break;
    default:
        break;
    }
    
    // add all leaf nodes' name into the alignment
    addLeafNamesToAlignment(tree->aln, tree->root, tree->root);
    
    // init Codon (if neccessary)
    if (tree->aln->seq_type == SEQ_CODON)
        tree->aln->initCodon(&tree->aln->sequence_type[5]);
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
    tree->params = params;
    tree->IQTree::initializeModel(*params, tree->aln->model_name, models_block);
}

/**
    remove all constant sites (in case with +ASC)
*/
void AliSimulator::removeConstantSites(){
    // dummy variables
    int num_variant_states = -1;
    vector<short int> variant_state_mask;
    
    // create a variant state mask
    createVariantStateMask(variant_state_mask, num_variant_states, round(expected_num_sites/length_ratio), tree->root, tree->root);
    
    // return error if num_variant_states is less than the expected_num_variant_states
    if (num_variant_states < round(expected_num_sites/length_ratio)){
        outError("Unfortunately, after removing constant sites, the number of variant sites is less than the expected sequence length. Please use --length-ratio <LENGTH_RATIO> to generate more abundant sites and try again. The current <LENGTH_RATIO> is "+ convertDoubleToString(length_ratio));
    }

    // recording start_time
    auto start = getRealTime();
    
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
#endif
    // get only variant sites for leaves
    getOnlyVariantSites(variant_state_mask, tree->root, tree->root);
    
    // show the time spent on copy variant sites
    auto end = getRealTime();
    cout<<" - Time spent on copying only variant sites: "<<end-start<<endl;
}

/**
    only get variant sites
*/
void AliSimulator::getOnlyVariantSites(vector<short int> variant_state_mask, Node *node, Node *dad){
    if (node->isLeaf() && node->name!=ROOT_NAME) {
#ifdef _OPENMP
#pragma omp task firstprivate(node)
#endif
        {
            // dummy sequence
            vector<short int> variant_sites;
            
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
                    // keep checking if Indels is used
                    if (num_variant_states >= round(expected_num_sites/length_ratio) && params->alisim_insertion_ratio == 0)
                        break;
                }
            
            // replace the sequence of the Leaf by variant sites
            node->sequence.clear();
            node->sequence = variant_sites;
        }
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
void AliSimulator::generatePartitionAlignment(vector<short int> ancestral_sequence, map<string,string> input_msa, string output_filepath)
{
    // if the ancestral sequence is not specified, randomly generate the sequence
    if (ancestral_sequence.size() == 0)
        tree->MTree::root->sequence = generateRandomSequence(expected_num_sites);
    // otherwise, using the ancestral sequence + abundant sites
    else
    {
        // set the ancestral sequence to the root node
        tree->MTree::root->sequence = ancestral_sequence;
        
        // add abundant_sites
        int num_abundant_sites = expected_num_sites - ancestral_sequence.size();
        if (num_abundant_sites > 0)
        {
            vector<short int> abundant_sites = generateRandomSequence(num_abundant_sites);
            for (int site:abundant_sites)
                tree->MTree::root->sequence.push_back(site);
        }
    }
    
    // validate the sequence length (in case of codon)
    validataSeqLengthCodon();
    
    // simulate the sequence for each node in the tree by DFS
    simulateSeqsForTree(input_msa, output_filepath);
}

/**
    create mask for variant sites
*/
void AliSimulator::createVariantStateMask(vector<short int> &variant_state_mask, int &num_variant_states, int expected_num_variant_states, Node *node, Node *dad){
    // no need to check the further sites if num_variant_states has exceeded the expected_num_variant_states
    // keep checking if Indels is used
    if (num_variant_states >= expected_num_variant_states && params->alisim_insertion_ratio == 0)
        return;
    
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // initialize the mask (all sites are assumed to be constant)
        if (num_variant_states == -1)
        {
            num_variant_states = 0;
            variant_state_mask = node->sequence;
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
                    // keep checking if Indels is used
                    if (num_variant_states >= expected_num_variant_states && params->alisim_insertion_ratio == 0)
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
*  by default (initial_freqs = true) freqs could be randomly generated if they are not specified
*/
vector<short int> AliSimulator::generateRandomSequence(int sequence_length, bool initial_freqs)
{
    // initialize sequence
    vector<short int> sequence;
    sequence.resize(sequence_length);
    
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
        // initialize state frequencies if they're not specified and initial_freqs = true
        if (initial_freqs)
            getStateFrequenciesFromModel(tree, state_freq);
        // otherwise, only get the state freqs from the model without re-initializing state freqs
        else
            tree->getModel()->getStateFrequency(state_freq);
        
        // finding the max probability position
        int max_prob_pos = 0;
        for (int i = 1; i < max_num_states; i++)
            if (state_freq[i] > state_freq[max_prob_pos])
                max_prob_pos = i;
        
        // randomly generate the sequence based on the state frequencies
        sequence = generateRandomSequenceFromStateFreqs(sequence_length, state_freq, max_prob_pos);
        
        // delete state_freq
        delete []  state_freq;
    }
    
    return sequence;
}

void AliSimulator::getStateFrequenciesFromModel(IQTree* tree, double *state_freqs){
    // firstly, initialize state freqs for mixture models (if neccessary)
    intializeStateFreqsMixtureModel(tree);
    
    // if a mixture model is used -> get weighted sum of state_freq across classes
    if (tree->getModel()->isMixture())
    {
        tree->getModel()->getStateFrequency(state_freqs, -1);
    }
    // get user-defined base frequencies (if any)
    else if ((tree->getModel()->getFreqType() == FREQ_USER_DEFINED)
        || (ModelLieMarkov::validModelName(tree->getModel()->getName()))
             || tree->aln->seq_type == SEQ_CODON
             || (tree->getModel()->getFreqType() == FREQ_EMPIRICAL && params->alisim_inference_mode))
        tree->getModel()->getStateFrequency(state_freqs);
    else // otherwise, randomly generate the base frequencies
    {
        
        // if sequence_type is dna -> randomly generate base frequencies based on empirical distributions
        if (tree->aln->seq_type == SEQ_DNA)
            random_frequencies_from_distributions(state_freqs);
        // otherwise, randomly generate base frequencies based on uniform distribution
        else
            generateRandomBaseFrequencies(state_freqs);
        tree->getModel()->setStateFrequency(state_freqs);
    }
}

/**
*  randomly generate the base frequencies
*/
void AliSimulator::generateRandomBaseFrequencies(double *base_frequencies)
{
    double sum = 0;
    
    // randomly generate the frequencies
    for (int i = 0; i < max_num_states; i++)
    {
        base_frequencies[i] = random_double();
        sum += base_frequencies[i];
    }
    
    // normalize the frequencies so that sum of them is 1
    for (int i = 0; i < max_num_states; i++)
    base_frequencies[i] /= sum;
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulator::simulateSeqsForTree(map<string,string> input_msa, string output_filepath)
{
    // get variables
    int sequence_length = expected_num_sites;
    ModelSubst *model = tree->getModel();
    ostream *out = NULL;
    vector<string> state_mapping;
    
    // initialize variables (site_specific_rates; site_specific_rate_index; site_specific_model_index)
    initVariables(sequence_length, true);
        
    // initialize trans_matrix
    double *trans_matrix = new double[params->num_threads*max_num_states*max_num_states];
    
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
    
    // rooting the tree if it's unrooted
    if (!tree->rooted)
        rootTree();
    
    // initialize sub_rates, J_Matrix from Q_matrix
    int num_mixture_models = model->getNMixtures();
    sub_rates = new double[num_mixture_models*max_num_states];
    Jmatrix = new double[num_mixture_models*max_num_states*max_num_states];
    extractRatesJMatrix(model);
    
    // simulate Sequences
    simulateSeqs(sequence_length, model, trans_matrix, tree->MTree::root, tree->MTree::root, *out, state_mapping, input_msa);
        
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
    
    // delete sub_rates, J_Matrix
    delete[] sub_rates;
    delete[] Jmatrix;
    
    // removing constant states if it's necessary
    if (length_ratio > 1)
        removeConstantSites();
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulator::simulateSeqs(int &sequence_length, ModelSubst *model, double *trans_matrix, Node *node, Node *dad, ostream &out, vector<string> state_mapping, map<string,string> input_msa)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // reset the num_children_done_simulation
        if (node->num_children_done_simulation >= (node->neighbors.size() - 1))
            node->num_children_done_simulation = 0;
        
        // select the appropriate simulation method
        SIMULATION_METHOD simulation_method = RATE_MATRIX;
        if ((*it)->length > params->alisim_simulation_thresh
            || tree->getRate()->isHeterotachy()
            || (*it)->attributes["model"].length()>0)
            simulation_method = TRANS_PROB_MATRIX;
        
        // dummy variables
        vector<short int> indel_sequence;
        vector<int> index_mapping_by_jump_step;
        
        // if branch_length is too short (less than 1 substitution is expected to occur) -> clone the sequence from the parent's node
        if ((*it)->length < 1.0/sequence_length)
            (*it)->node->sequence = node->sequence;
        else
        {
            // handle indels
            if (params->alisim_insertion_ratio + params->alisim_deletion_ratio != 0 || simulation_method == RATE_MATRIX)
                handleIndels(model, sequence_length, node, it, indel_sequence, index_mapping_by_jump_step, simulation_method);
            
            // only simulate new sequence if the simulation type is Transition Probability Matrix approach
            if (simulation_method == TRANS_PROB_MATRIX)
            {
                // if a model is specify for the current branch -> simulate the sequence based on that branch-specific model
                if ((*it)->attributes["model"].length()>0)
                    branchSpecificEvolution(sequence_length, trans_matrix, node, it);
                // otherwise, simulate the sequence based on the common model
                else
                    simulateASequenceFromBranchAfterInitVariables(model, sequence_length, trans_matrix, node, it);
            }
            // otherwise (Rate_matrix is used as the simulation method) -> use the indel_sequence
            else
                (*it)->node->sequence = indel_sequence;
        }
        
        // permuting selected sites for FunDi model
        if (params->alisim_fundi_taxon_set.size()>0)
        {
            if (node->isLeaf())
                permuteSelectedSites(fundi_items, node);
            if ((*it)->node->isLeaf())
                permuteSelectedSites(fundi_items, (*it)->node);
        }
        
        // merge the simulated sequence with the indel_sequence
        if (params->alisim_insertion_ratio + params->alisim_deletion_ratio != 0 && (*it)->length >= 1.0/sequence_length && simulation_method == TRANS_PROB_MATRIX)
            mergeIndelSequence((*it)->node, indel_sequence, index_mapping_by_jump_step);
        
        // writing and deleting simulated sequence immediately if possible
        writeAndDeleteSequenceImmediatelyIfPossible(out, state_mapping, input_msa, it, node);
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(sequence_length, model, trans_matrix, (*it)->node, node, out, state_mapping, input_msa);
    }
}

/**
    writing and deleting simulated sequence immediately if possible
*/
void AliSimulator::writeAndDeleteSequenceImmediatelyIfPossible(ostream &out, vector<string> state_mapping, map<string,string> input_msa, NeighborVec::iterator it, Node* node)
{
    // write sequence of leaf nodes to file if possible
    if (state_mapping.size() > 0)
    {
        if ((*it)->node->isLeaf())
        {
            string input_sequence = input_msa[(*it)->node->name];
            // convert numerical states into readable characters and write output to file
            if (input_sequence.length()>0)
                // write and copying gaps from the input sequences to the output.
                out<<writeASequenceToFileWithGaps((*it)->node, round(expected_num_sites/length_ratio), num_sites_per_state, input_sequence, state_mapping, params->aln_output_format, max_length_taxa_name);
            else
                // write without copying gaps from the input sequences to the output.
                out<< convertNumericalStatesIntoReadableCharacters((*it)->node, round(expected_num_sites/length_ratio), num_sites_per_state, state_mapping, params->aln_output_format, max_length_taxa_name);
            
            // remove the sequence to release the memory after extracting the sequence
            vector<short int>().swap((*it)->node->sequence);
        }
        
        if (node->isLeaf())
        {
            // avoid writing sequence of __root__
            if (node->name!=ROOT_NAME)
            {
                string input_sequence = input_msa[node->name];
                // convert numerical states into readable characters and write output to file
                if (input_sequence.length()>0)
                    // write and copying gaps from the input sequences to the output.
                    out<<writeASequenceToFileWithGaps(node, round(expected_num_sites/length_ratio), num_sites_per_state, input_sequence, state_mapping, params->aln_output_format, max_length_taxa_name);
                else
                    // write without copying gaps from the input sequences to the output.
                    out<< convertNumericalStatesIntoReadableCharacters(node, round(expected_num_sites/length_ratio), num_sites_per_state, state_mapping, params->aln_output_format, max_length_taxa_name);
            }
            
            // remove the sequence to release the memory after extracting the sequence
            vector<short int>().swap(node->sequence);
        }
    }
    
    // update the num_children_done_simulation
    node->num_children_done_simulation++;
    // remove the sequence of the current node to release the memory if it's an internal node && all of its children have done their simulation && the user does't (use indel && want to output the internal sequences)
    if (!node->isLeaf() && node->num_children_done_simulation >= (node->neighbors.size() - 1)
        && !(params->alisim_insertion_ratio + params->alisim_deletion_ratio != 0 && params->alisim_write_internal_sequences))
    {
        // convert numerical states into readable characters and write internal sequences to file if the user want to do so
        if (params->alisim_write_internal_sequences && state_mapping.size() > 0)
            out<< convertNumericalStatesIntoReadableCharacters(node, round(expected_num_sites/length_ratio), num_sites_per_state, state_mapping, params->aln_output_format, max_length_taxa_name);
        
        // release the memory
        vector<short int>().swap(node->sequence);
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
*  get a random item from a set of items with an accumulated probability array by binary search starting at the max probability
*/
int AliSimulator::getRandomItemWithAccumulatedProbMatrixMaxProbFirst(double *accumulated_probability_maxtrix, int starting_index, int num_columns, int max_prob_position){
    // generate a random number
    double random_number = random_double();
    
    // starting at the probability of unchange first
    if (random_number >= (max_prob_position==0?0:accumulated_probability_maxtrix[starting_index+max_prob_position-1]))
    {
        if (random_number <= accumulated_probability_maxtrix[starting_index+max_prob_position])
            return max_prob_position;
        // otherwise, searching on the right part
        else
            return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, starting_index+max_prob_position+1, starting_index+(num_columns-1), starting_index)-starting_index;
    }
    
    // otherwise, searching on the left part
    return binarysearchItemWithAccumulatedProbabilityMatrix(accumulated_probability_maxtrix, random_number, starting_index, starting_index+max_prob_position-1, starting_index)-starting_index;
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
*  binary search an item from a set with accumulated probability array
*/
int AliSimulator::binarysearchItemWithAccumulatedProbabilityMatrix(vector<double> accumulated_probability_maxtrix, double random_number, int start, int end, int first)
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
    if (tree->aln->seq_type == SEQ_CODON && (!params->partition_file && params->alisim_sequence_length%3))
        outError("Sequence length of Codon must be divisible by 3. Please check & try again!");
}

/**
*  update the expected_num_sites due to the change of the sequence_length
*/
void AliSimulator::refreshExpectedNumSites(){
    expected_num_sites = params->alisim_sequence_length/num_sites_per_state*length_ratio;
}

/**
    estimate length_ratio (for models with +ASC)
*/
void AliSimulator::estimateLengthRatio()
{
    // By default (without +ASC), length_ratio is set at 1
    length_ratio = 1;
        
    // Handle the case with +ASC
    if (tree->getModel() && tree->getSubstName().find("+ASC") != std::string::npos)
    {
        // using the length_ratio in params if it's specified by the user
        if (tree->params->original_params.find("--length-ratio") != std::string::npos)
            length_ratio = params->alisim_length_ratio;
        // otherwise, estimating the length_ratio
        else
        {
            // disable ASC for computing likelihood score
            ASCType asc_type = tree->getModelFactory()->getASC();
            tree->getModelFactory()->setASC(ASC_NONE);
            
            // initialize a string concatenating all characters of all states (eg, ACGT for DNA)
            string all_characters;
            all_characters.resize(max_num_states*num_sites_per_state);
            for (int i = 0; i < max_num_states; i++)
            {
                string characters_from_state = tree->aln->convertStateBackStr(i);
                for (int j = 0; j < num_sites_per_state; j++)
                    all_characters[i*num_sites_per_state+j] = characters_from_state[j];
            }
            
            // convert tree to unrooted tree if it's now rooted
            if (tree->rooted)
            {
                outWarning("The input tree is now converting into unrooted tree.");
                tree->PhyloTree::forceConvertingToUnrooted();
            }
           
            // initialize sequences (a dummy alignment with all sequences are set to all_characters)
            StrVector sequences;
            int nseq = tree->getNumTaxa(), nsite = max_num_states;
            sequences.resize(nseq);
            for (int i = 0; i < nseq; i++)
                sequences[i] = all_characters;
            
            // build all constant site patterns
            char *sequence_type = strcpy(new char[tree->aln->sequence_type.length() + 1], tree->aln->sequence_type.c_str());
            tree->aln->buildPattern(sequences, sequence_type, nseq, nsite*num_sites_per_state);
            
            // compute the likelihood scores of all patterns
            double *patterns_llh = new double[tree->aln->getNPattern()];
            tree->setLikelihoodKernel(params->SSE);
            tree->setNumThreads(params->num_threads);
            tree->initializeAllPartialLh();
            tree->computeLikelihood(patterns_llh);
            
            // initialize the estimated_length_ratio
            double estimated_length_ratio = 0;
            
            // take the sum of all probabilities of all constant patterns
            for (int i = 0; i < max_num_states; i++)
                estimated_length_ratio += exp(patterns_llh[i]);
            
            // delete patterns_llh
            delete [] patterns_llh;
            
            // set ASC type to its original value
            tree->getModelFactory()->setASC(asc_type);
            
            // handle the case when estimated_length_ratio is estimated incorrectly
            if (!isfinite(estimated_length_ratio) || estimated_length_ratio > 1)
                estimated_length_ratio = 0.5;
            
            // update the length_ratio with a 10% (0.1) additional length_ratio (for backup purpose)
            length_ratio = 1/(1-estimated_length_ratio) + 0.1;
        }
    }
}

/**
*  initialize state_mapping (mapping from states into characters)
*
*/
void AliSimulator::initializeStateMapping(int num_sites_per_state, Alignment *aln, vector<string> &state_mapping)
{
    ASSERT(aln);
    
    // initialize state_mapping (mapping from state to characters)
    int total_states = aln->STATE_UNKNOWN + 1;
    state_mapping.resize(total_states);
    for (int i = 0; i< total_states; i++)
        state_mapping[i] = aln->convertStateBackStr(i);

    // add an additional state for gap
    if (num_sites_per_state == 3)
        state_mapping[total_states-1] = "---";
}

/**
*  convert numerical states into readable characters
*
*/
string AliSimulator::convertNumericalStatesIntoReadableCharacters(Node *node, int sequence_length, int num_sites_per_state, vector<string> state_mapping, InputType output_format, int max_length_taxa_name)
{
    ASSERT(sequence_length <= node->sequence.size());
    
    // dummy variables
    std::string output (sequence_length * num_sites_per_state+1, ' ');
    int start_index;
    
    // add node's name
    // in PHYLIP format
    if (output_format != IN_FASTA)
    {
        // add padding to node_name
        string name_with_padding = node->name;
        // write node's id if node's name is empty
        if (name_with_padding.length() == 0) name_with_padding = convertIntToString(node->id);
        ASSERT(max_length_taxa_name >= name_with_padding.length());
        std::string padding (max_length_taxa_name - name_with_padding.length() + 1, ' ');
        name_with_padding += padding;
        output = name_with_padding + output;
        start_index = name_with_padding.length();
    }
    // in FASTA format
    else
    {
        string node_name = node->name;
        // write node's id if node's name is empty
        if (node_name.length() == 0) node_name = convertIntToString(node->id);
        output = ">" + node_name + "\n" + output;
        start_index = node_name.length() + 2;
    }
    output[output.length()-1] = '\n';
    
    // convert normal data
    if (num_sites_per_state == 1)
        for (int i = 0; i < sequence_length; i++)
            output[start_index+i*num_sites_per_state] = state_mapping[node->sequence[i]][0];
    // convert CODON
    else
        for (int i = 0; i < sequence_length; i++)
        {
            output[start_index+i*num_sites_per_state] = state_mapping[node->sequence[i]][0];
            output[start_index+i*num_sites_per_state + 1] = state_mapping[node->sequence[i]][1];
            output[start_index+i*num_sites_per_state + 2] = state_mapping[node->sequence[i]][2];
        }
    
    // return output
    return output;
}

/**
    show warning if base frequencies are set/unset correctly (only check DNA models)
*/
void AliSimulator::checkBaseFrequenciesDNAModels(IQTree* tree, string model_name){
    if (tree->aln && tree->aln->seq_type == SEQ_DNA && !params->partition_file && model_name.find("MIX") == std::string::npos) {
        
        // initializing the list of unequal/equal base frequencies models
        vector<string> unequal_base_frequencies_models = vector<string>{"GTR", "F81", "HKY", "HKY85", "TN", "TN93", "K81u", "TPM2u", "TPM3u", "TIM", "TIM2", "TIM3", "TVM"};
        vector<string> equal_base_frequencies_models = vector<string>{"JC", "JC69", "K80", "K2P", "TNe", "K81", "K3P", "TPM2", "TPM3", "TIMe", "TIM2e", "TIM3e", "TVMe", "SYM"};
        
        // check whether base frequencies are not set for unequal base frequenceies models
        for (string model_item: unequal_base_frequencies_models)
            if (model_name.find(model_item) != std::string::npos && model_name.find("+F") == std::string::npos) {
                outWarning(model_item+" must have unequal base frequencies. The base frequencies could be randomly generated if users do not provide them. However, we strongly recommend users specify the base frequencies by using +F{freq1/.../freqN} for better simulation accuracy.");
                break;
            }
        
        // check whether base frequencies are set for equal base frequenceies models
        for (string model_item: equal_base_frequencies_models)
            if (model_name.find(model_item) != std::string::npos && model_name.find("+F") != std::string::npos) {
                outWarning(model_item+" must have equal base frequencies. Unequal base frequencies specified by users could lead to incorrect simulation. We strongly recommend users to not specify the base frequencies for this model by removing +F{freq1/.../freqN}.");
                break;
            }
    }
}

/**
    extract the maximum length of taxa names
*/
short int AliSimulator::extractMaxTaxaNameLength()
{
    if (tree && tree->aln)
    {
        // if it's a super tree -> check each tree one by one
        if (tree->isSuperTree())
        {
            for (int i = 0 ; i < ((PhyloSuperTree*) tree)->size(); i++)
            {
                IQTree *current_tree = (IQTree *) ((PhyloSuperTree*) tree)->at(i);
                vector<string> seq_names = current_tree->aln->getSeqNames();
                for (int i = 0; i < seq_names.size(); i++)
                    if (seq_names[i].length()>max_length_taxa_name)
                        max_length_taxa_name = seq_names[i].length();
            }
        }
        // otherwise, just check the current tree
        else
        {
            vector<string> seq_names = tree->aln->getSeqNames();
            for (int i = 0; i < seq_names.size(); i++)
                if (seq_names[i].length()>max_length_taxa_name)
                    max_length_taxa_name = seq_names[i].length();
        }
    }
}

/**
    selecting & permuting sites (FunDi models)
*/
vector<FunDi_Item> AliSimulator::selectAndPermuteSites(double proportion, int num_sites){
    ASSERT(proportion<1);
    
    // dummy variables
    vector<FunDi_Item> fundi_items;
    IntVector tmp_selected_sites;
    int num_selected_sites = round(proportion*num_sites);
    
    // select random unique sites one by one
    for (int i = 0; i < num_selected_sites; i++)
    {
        // attempt up to 1000 times to select a random site
        for (int j = 0; j < 1000; j++)
        {
            int random_site = random_int(num_sites);
            
            // check if the random_site has been already selected or not
            if (std::find(tmp_selected_sites.begin(), tmp_selected_sites.end(), random_site) != tmp_selected_sites.end())
                // retry if the random_site has already existed in the selected list
                continue;
            else
            {
                // add the random site to the selected list
                tmp_selected_sites.push_back(random_site);
                break;
            }
        }
        
        if (tmp_selected_sites.size() <= i)
            outError("Failed to select random sites for permutations (of FunDi model) after 1000 attempts");
    }
    
    // select a new position for each of the first num_selected_sites - 1 selected sites
    IntVector position_pool(tmp_selected_sites);
    for (int i = 0; i < num_selected_sites - 1; i++)
    {
        // attempt up to 1000 times
        for (int j = 0; j < 1000; j++)
        {
            int rand_num = random_int(position_pool.size());
            int new_position = position_pool[rand_num];
            
            // if new_position == current_position, then retry
            if (new_position == tmp_selected_sites[i])
                continue;
            // otherwise, it is a valid new position
            else
            {
                FunDi_Item tmp_fundi_item = {tmp_selected_sites[i],new_position};
                fundi_items.push_back(tmp_fundi_item);
                // remove the new_position from the position pool
                position_pool.erase(position_pool.begin() + rand_num);
                break;
            }
        }
        
        if (fundi_items.size() <= i) {
            outError("Failed to select a positions to permute the selected sites (of FunDi model) after 1000 attempts");
        }
    }
    // select a new position for the last selected site
    ASSERT(position_pool.size() == 1);
    if (tmp_selected_sites[tmp_selected_sites.size()-1] != position_pool[0])
    {
        FunDi_Item tmp_fundi_item = {tmp_selected_sites[tmp_selected_sites.size()-1], position_pool[0]};
        fundi_items.push_back(tmp_fundi_item);
    }
    else
    {
        FunDi_Item tmp_fundi_item = {tmp_selected_sites[tmp_selected_sites.size()-1], fundi_items[0].new_position};
        fundi_items.push_back(tmp_fundi_item);
        fundi_items[0].new_position = position_pool[0];
    }
    
    return fundi_items;
}

/**
    permuting selected sites (FunDi models)
*/
void AliSimulator::permuteSelectedSites(vector<FunDi_Item> fundi_items, Node* node)
{
    if (std::find(params->alisim_fundi_taxon_set.begin(), params->alisim_fundi_taxon_set.end(), node->name) != params->alisim_fundi_taxon_set.end()) {
            // caching the current states of all selected sites
            map<int, short int> caching_sites;
            for (int i = 0; i < fundi_items.size(); i++)
                caching_sites[fundi_items[i].selected_site] = node->sequence[fundi_items[i].selected_site];
            
            // permuting sites in FunDi model
            for (int i = 0; i < fundi_items.size(); i++)
                node->sequence[fundi_items[i].new_position] = caching_sites[fundi_items[i].selected_site];
        }
}

/**
    initialize state freqs for all model components (of a mixture model)
*/
void AliSimulator::intializeStateFreqsMixtureModel(IQTree* tree)
{
    // get/init variables
    ModelSubst* model = tree->getModel();
    
    // only initialize state freqs if it's a mixture model && the state freqs have not been estimated by an inference process yet
    if (model->isMixture() && !params->alisim_inference_mode && model->getFreqType() == FREQ_EMPIRICAL)
    {
        // initialize state freqs
        double *state_freq = new double[max_num_states];
        
        // get the weights of model components
        for (int i = 0; i < model->getNMixtures(); i++)
            if (model->getMixtureClass(i)->getFreqType() == FREQ_EMPIRICAL)
            {
                // if sequence_type is dna -> randomly generate base frequencies based on empirical distributions
                if (tree->aln->seq_type == SEQ_DNA)
                    random_frequencies_from_distributions(state_freq);
                // otherwise, randomly generate base frequencies based on uniform distribution
                else
                    generateRandomBaseFrequencies(state_freq);
                
                model->getMixtureClass(i)->setStateFrequency(state_freq);
            }
        
        // delete state_freq
        delete [] state_freq;
    }
}

/**
    branch-specific evolution
*/
void AliSimulator::branchSpecificEvolution(int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it)
{
    // initialize a dummy model for this branch
    string model_full_name = (*it)->attributes["model"];
    // convert separator from "/" to ","
    std::replace(model_full_name.begin(), model_full_name.end(), '/', ',');
    IQTree *tmp_tree = new IQTree();
    tmp_tree->copyPhyloTree(tree, true);
    initializeModel(tmp_tree, model_full_name);
    
    // initialize state frequencies
    double *state_freqs = new double[max_num_states];
    getStateFrequenciesFromModel(tmp_tree, state_freqs);
    delete[] state_freqs;
    
    // check if base frequencies for DNA models are specified correctly
    checkBaseFrequenciesDNAModels(tmp_tree, model_full_name);
    
    // handle Heterotachy model in branch-specific models
    string lengths = "";
    if (tmp_tree->getRate()->isHeterotachy())
    {
        // make sure that the user has specified multiple lengths for the current branch
        lengths = (*it)->attributes["lengths"];
        if (lengths.length() == 0)
            outError("To use Heterotachy model, please specify multiple lengths for the current branch by [&model=...,lengths=<length_0>/.../<length_n>]");
    }
    
    // initialize a new dummy alisimulator
    AliSimulator* tmp_alisimulator = new AliSimulator(params, tmp_tree, expected_num_sites, partition_rate);
    
    // convert alisimulator to the correct type of simulator
    // get variables
    string rate_name = tmp_alisimulator->tree->getRateName();
    double invariant_proportion = tmp_alisimulator->tree->getRate()->getPInvar();
    bool is_mixture_model = tmp_alisimulator->tree->getModel()->isMixture();
    // case 1: without rate heterogeneity or mixture model -> using the current alisimulator (don't need to re-initialize it)
    // case 2: with rate heterogeneity or mixture model
    if ((!rate_name.empty()) || is_mixture_model)
    {
        // if user specifies +I without invariant_rate -> set it to 0
        if (rate_name.find("+I") != std::string::npos && isnan(invariant_proportion)){
            tmp_alisimulator->tree->getRate()->setPInvar(0);
            outWarning("Invariant rate is now set to Zero since it has not been specified");
        }
        
        // case 2.3: with only invariant sites (without gamma/freerate model/mixture models)
        if (!rate_name.compare("+I") && !is_mixture_model)
            tmp_alisimulator = new AliSimulatorInvar(tmp_alisimulator, invariant_proportion);
        else
        {
            // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
            if (invariant_proportion > 0)
                tmp_alisimulator = new AliSimulatorHeterogeneityInvar(tmp_alisimulator, invariant_proportion);
            // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
            else
                tmp_alisimulator = new AliSimulatorHeterogeneity(tmp_alisimulator);
        }
    }
    
    // print model's parameters
    cout<<"Simulating a sequence with branch-specific model named "+tmp_tree->getModel()->getName()<<endl;
    tmp_tree->getModel()->writeInfo(cout);
    
    // simulate the sequence for the current node based on the branch-specific model
    tmp_alisimulator->simulateASequenceFromBranch(tmp_tree->getModel(), sequence_length, trans_matrix, node, it, lengths);
    
    // delete the dummy alisimulator
    delete tmp_alisimulator;
}

/**
    simulate a sequence for a node from a specific branch
*/
void AliSimulator::simulateASequenceFromBranch(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths)
{
    // initialize the site-specific rates
    initVariables(sequence_length);
    
    // check to regenerate the root sequence if the user has specified specific frequencies for root
    if (tree->root->id == node->id && ((*it)->attributes["freqs"]).length() > 0)
        regenerateRootSequenceBranchSpecificModel((*it)->attributes["freqs"], sequence_length, node);
    
    // simulate a sequence for a node from a specific branch after all variables has been initializing
    simulateASequenceFromBranchAfterInitVariables(model, sequence_length, trans_matrix, node, it, lengths);
}

/**
    simulate a sequence for a node from a specific branch after all variables has been initializing
*/
void AliSimulator::simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths)
{
    // compute the transition probability matrix
    model->computeTransMatrix(partition_rate * (*it)->length, trans_matrix);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(trans_matrix, max_num_states, max_num_states);
    
    // estimate the sequence for the current neighbor
    (*it)->node->sequence.resize(sequence_length);
    for (int i = 0; i < sequence_length; i++)
    {
        // if the parent's state is a gap -> the children's state should also be a gap
        if (node->sequence[i] == STATE_UNKNOWN)
            (*it)->node->sequence[i] = STATE_UNKNOWN;
        else
        {
            // iteratively select the state for each site of the child node, considering it's dad states, and the transition_probability_matrix
            int starting_index = node->sequence[i]*max_num_states;
            (*it)->node->sequence[i] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(trans_matrix, starting_index, max_num_states, node->sequence[i]);
        }
    }
}

/**
    initialize variables (e.g., site-specific rate)
*/
void AliSimulator::initVariables(int sequence_length, bool regenerate_root_sequence)
{
    // Do nothing, this method will be overrided in AliSimulatorHeterogeneity and AliSimulatorInvar
}

/**
    regenerate the root sequence if the user has specified specific state frequencies in branch-specific model
*/
void AliSimulator::regenerateRootSequenceBranchSpecificModel(string freqs, int sequence_length, Node* root){
    // initizlize state_freqs
    double* state_freqs = new double[max_num_states];
    
    // parse state_freqs
    int i = 0;
    int max_prob_pos = -1;
    double total_freq = 0;
    while (freqs.length() > 0) {
        // split state_freqs by "/"
        size_t pos = freqs.find('/');
        
        // convert frequency from string to double
        state_freqs[i] = convert_double_with_distribution(freqs.substr(0, pos).c_str());
        total_freq += state_freqs[i];
        
        // update the position with the highest frequency
        if (max_prob_pos == -1 || state_freqs[i] > state_freqs[max_prob_pos])
            max_prob_pos = i;
        
        // remove the current frequency from the list freqs
        if (pos != std::string::npos)
            freqs.erase(0, pos + 1);
        else
            freqs = "";
        
        i++;
    }
    
    // make sure that the number of user-defined frequencies is equal to the number of states
    if (i != max_num_states)
        outError("The number of frequencies ("+convertIntToString(i)+") is different from the number of states ("+convertIntToString(max_num_states)+"). Please check and try again!");
    
    // make sure the sum of all frequencies is equal to 1
    if (fabs(total_freq-1.0) > 1e-5)
    {
        outWarning("Normalizing state frequencies so that sum of them is equal to 1.");
        normalize_frequencies(state_freqs, max_num_states, total_freq);
    }
    
    // re-generate a new sequence for the root from the state frequencies
    root->sequence = generateRandomSequenceFromStateFreqs(sequence_length, state_freqs, max_prob_pos);
    
    // release the memory of state_freqs
    delete[] state_freqs;
}

/**
    generate a random sequence by state frequencies
*/
vector<short int> AliSimulator::generateRandomSequenceFromStateFreqs(int sequence_length, double* state_freqs, int max_prob_pos)
{
    vector<short int> sequence;
    sequence.resize(sequence_length);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(state_freqs, 1, max_num_states);
    
    // randomly generate each site in the sequence follows the base frequencies defined by the user
    for (int i = 0; i < sequence_length; i++)
        sequence[i] =  getRandomItemWithAccumulatedProbMatrixMaxProbFirst(state_freqs, 0, max_num_states, max_prob_pos);
    
    return sequence;
}

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
string AliSimulator::writeASequenceToFileWithGaps(Node *node, int sequence_length, int num_sites_per_state, string input_sequence, vector<string> state_mapping, InputType output_format, int max_length_taxa_name)
{
    // initialize the output sequence with all gaps (to handle the cases with missing taxa in partitions)
    string output (sequence_length * num_sites_per_state+1, '-');
    int start_index;
    
    // add node's name
    // in PHYLIP format
    if (output_format != IN_FASTA)
    {
        // add padding to node_name
        string name_with_padding = node->name;
        // write node's id if node's name is empty
        if (name_with_padding.length() == 0) name_with_padding = convertIntToString(node->id);
        ASSERT(max_length_taxa_name >= name_with_padding.length());
        std::string padding (max_length_taxa_name - name_with_padding.length() + 1, ' ');
        name_with_padding += padding;
        output = name_with_padding + output;
        start_index = name_with_padding.length();
    }
    // in FASTA format
    else
    {
        string node_name = node->name;
        // write node's id if node's name is empty
        if (node_name.length() == 0) node_name = convertIntToString(node->id);
        output = ">" + node_name + "\n" + output;
        start_index = node_name.length() + 2;
    }
    output[output.length()-1] = '\n';
    
    // convert non-empty sequence
    if (node->sequence.size() >= sequence_length)
    {
        // convert normal data
        if (num_sites_per_state == 1)
        {
            for (int i = 0; i < sequence_length; i++){
                // handle gaps
                if ((i+1)*num_sites_per_state - 1 < input_sequence.length()
                    && input_sequence[i] == '-')
                {
                    // insert gaps
                    output[start_index+i*num_sites_per_state] = '-';
                }
                // if it's not a gap
                else
                    output[start_index+i*num_sites_per_state] = state_mapping[node->sequence[i]][0];
            }
        }
        // convert CODON
        else {
            for (int i = 0; i < sequence_length; i++){
                // handle gaps
                if ((i+1)*num_sites_per_state - 1 < input_sequence.length()
                    && (input_sequence[i*num_sites_per_state] == '-'
                            || input_sequence[i*num_sites_per_state+1] == '-'
                            || input_sequence[i*num_sites_per_state+2] == '-')){
                    // insert gaps
                    output[start_index+i*num_sites_per_state] =  input_sequence[i*num_sites_per_state];
                    output[start_index+i*num_sites_per_state + 1] =  input_sequence[i*num_sites_per_state+1];
                    output[start_index+i*num_sites_per_state + 2] =  input_sequence[i*num_sites_per_state+2];
                }
                else
                {
                        output[start_index+i*num_sites_per_state] = state_mapping[node->sequence[i]][0];
                        output[start_index+i*num_sites_per_state + 1] = state_mapping[node->sequence[i]][1];
                        output[start_index+i*num_sites_per_state + 2] = state_mapping[node->sequence[i]][2];
                }
            }
        }
    }
    
    return output;
}

/**
    extract array of substitution rates
*/
double AliSimulator::extractRatesJMatrix(ModelSubst *model)
{
    // get num_mixture_models
    int num_mixture_models = model->getNMixtures();
    double* tmp_Q_matrix = new double[max_num_states*max_num_states];
    
    for (int mixture = 0; mixture <num_mixture_models; mixture++)
    {
        // get the Rate (Q) Matrix
        model->getQMatrix(tmp_Q_matrix, mixture);
        
        // extract sub rate for each state from q_matrix
        int starting_index_sub_rates = mixture * max_num_states;
        for (int i = 0; i < max_num_states; i++)
            sub_rates[starting_index_sub_rates + i] = - tmp_Q_matrix[i * (max_num_states + 1)];
        
        // convert Q_Matrix to J_Matrix;
        int starting_index_J = starting_index_sub_rates * max_num_states;
        for (int i = 0; i < max_num_states; i++)
        for (int j = 0; j < max_num_states; j++)
        {
            if (i == j)
                Jmatrix[starting_index_J+i*max_num_states+j] = 0;
            else
                Jmatrix[starting_index_J+i*max_num_states+j] = tmp_Q_matrix[i*max_num_states+j]/sub_rates[starting_index_sub_rates + i];
        }
    }
    
    // delete tmp_Q_matrix
    delete[] tmp_Q_matrix;
    
    // convert Jmatrix to accumulated Jmatrix
    convertProMatrixIntoAccumulatedProMatrix(Jmatrix, num_mixture_models*max_num_states, max_num_states);
}

/**
    compute the total substitution rate
*/
double AliSimulator::computeTotalSubRate(vector<double> site_specific_rates, vector<short int> site_specific_model_index, vector<short int> sequence)
{
    // initialize variables
    double total_sub_rate = 0;
    int sequence_length = sequence.size();
    
    // compute the total sub rate (Notes: this method could be further optimized to reduce the running time)
    for (int i = 0; i < sequence_length; i++)
    {
        // not compute the substitution rate for gaps/deleted sites
        if (sequence[i] != STATE_UNKNOWN)
        {
            // get mixture_index
            int mixture_index = site_specific_model_index.size() == 0?0:site_specific_model_index[i];
            // if site_specific_rates is empty, the relative site rate should be 1
            if (site_specific_rates.size() > 0)
                total_sub_rate += sub_rates[mixture_index * max_num_states + sequence[i]]*site_specific_rates[i];
            else
                total_sub_rate += sub_rates[mixture_index * max_num_states + sequence[i]];
        }
    }
    
    return total_sub_rate;
}

/**
    initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
*/
void AliSimulator::initVariables4RateMatrix(double &total_sub_rate, int &num_gaps, vector<double> &accummulated_rates, vector<short int> sequence)
{
    // initialize variables
    int sequence_length = sequence.size();
    total_sub_rate = 0;
    num_gaps = 0;
    accummulated_rates.resize(sequence_length, 0);
    
    // compute the total sub rate (Notes: this method could be further optimized to reduce the running time)
    for (int i = 0; i < sequence_length; i++)
    {
        // not compute the substitution rate for gaps/deleted sites
        if (sequence[i] != STATE_UNKNOWN)
        {
            // get the mixture model index
            int mixture_index = site_specific_model_index.size() == 0? 0:site_specific_model_index[i];
            // if site_specific_rates is empty, the relative site rate should be 1
            if (site_specific_rates.size() > 0)
                total_sub_rate += sub_rates[mixture_index * max_num_states + sequence[i]]*site_specific_rates[i];
            else
                total_sub_rate += sub_rates[mixture_index * max_num_states + sequence[i]];
        }
        else
            num_gaps++;
        
        // update accummulated_rates
        accummulated_rates[i] = total_sub_rate;
    }
}

/**
    handle indels
*/
void AliSimulator::handleIndels(ModelSubst *model, int &sequence_length, Node *node, NeighborVec::iterator it, vector<short int> &indel_sequence, vector<int> &index_mapping_by_jump_step, SIMULATION_METHOD simulation_method)
{
    int num_gaps = 0;
    double total_sub_rate = 0;
    vector<double> accummulated_rates;
    // If AliSim is using RATE_MATRIX approach -> initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
    if (simulation_method == RATE_MATRIX)
        initVariables4RateMatrix(total_sub_rate, num_gaps, accummulated_rates, node->sequence);
    else // otherwise, TRANS_PROB_MATRIX approach is used -> only count the number of gaps
        num_gaps = count(node->sequence.begin(), node->sequence.end(), STATE_UNKNOWN);
    double total_ins_rate = params->alisim_insertion_ratio*(sequence_length + 1 - num_gaps);
    double total_del_rate = params->alisim_deletion_ratio*(sequence_length - 1 - num_gaps + computeMeanDelSize(sequence_length));
    double total_event_rate = total_sub_rate + total_ins_rate + total_del_rate;
    
    indel_sequence = node->sequence;
    index_mapping_by_jump_step.resize(sequence_length + 1, 0);
    double branch_length = (*it)->length;
    while (branch_length > 0)
    {
        // generate a waiting time s1 by sampling from the exponential distribution with mean 1/total_event_rate
        double waiting_time = random_double_exponential_distribution(1/total_event_rate);
        
        if (waiting_time > branch_length)
            break;
        else
        {
            // update the remaining length of the current branch
            branch_length -=  waiting_time;
            
            // Determine the event type (insertion, deletion, substitution) occurs
            double random_num = random_double()*total_event_rate;
            EVENT_TYPE event_type = SUBSTITUTION;
            if (random_num < total_ins_rate)
                event_type = INSERTION;
            else if (random_num < total_ins_rate+total_del_rate)
                event_type = DELETION;
            
            // process event
            int length_change = 0;
            switch (event_type)
            {
                case INSERTION:
                {
                    length_change = handleInsertion(sequence_length, index_mapping_by_jump_step, indel_sequence, total_sub_rate, accummulated_rates, simulation_method);
                    break;
                }
                case DELETION:
                {
                    length_change = -handleDeletion(sequence_length, indel_sequence, total_sub_rate, accummulated_rates, simulation_method);
                    break;
                }
                case SUBSTITUTION:
                {
                    if (simulation_method == RATE_MATRIX)
                    {
                        handleSubs(sequence_length, total_sub_rate, accummulated_rates, indel_sequence);
                    }
                    break;
                }
                default:
                    break;
            }
            
            // update total_event_rate
            if (length_change != 0)
            {
                total_ins_rate += params->alisim_insertion_ratio * length_change;
                total_del_rate += params->alisim_deletion_ratio * length_change;
            }
            total_event_rate = total_sub_rate + total_ins_rate + total_del_rate;
        }

    }
    
    // insert gaps to other nodes
    if (index_mapping_by_jump_step[index_mapping_by_jump_step.size() - 1]>0)
    {
        bool stop_inserting_gaps = false;
        insertGapsForInsertionEvents(index_mapping_by_jump_step, (*it)->node->id, tree->root, tree->root, stop_inserting_gaps);
    }
}

/**
*  insert a new sequence into the current sequence
*
*/
void AliSimulator::insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence)
{
    indel_sequence.insert(indel_sequence.begin()+position, new_sequence.begin(), new_sequence.end());
}

/**
*  insert gaps into other nodes when processing Insertion Events
*
*/
void AliSimulator::insertGapsForInsertionEvents(vector<int> index_mapping_by_jump_step, int stopping_node_id, Node *node, Node *dad, bool &stop_inserting_gaps)
{
    // check to stop
    if (stop_inserting_gaps)
        return;
    
    // only insert gaps into not-empty node
    if (node->sequence.size() > 0)
    {
        // resize the sequence with gaps
        node->sequence.resize(index_mapping_by_jump_step.size() + index_mapping_by_jump_step[index_mapping_by_jump_step.size() - 1] - 1, STATE_UNKNOWN);
        
        // update sites to new positions
        for (int i = index_mapping_by_jump_step.size() - 1 - 1; i >= 0; i--)
        {
            if (index_mapping_by_jump_step[i] > 0)
            {
                node->sequence[i + index_mapping_by_jump_step[i]] = node->sequence[i];
                
                // if (*it)->node->sequence[i] is not overrided, it should be a gap
                node->sequence[i] = STATE_UNKNOWN;
            }
            else
                break;
        }
        
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // stop traversing the tree if reaching the stopping node
        if ((*it)->node->id == stopping_node_id)
        {
            stop_inserting_gaps = true;
            break;
        }
        
        // browse 1-step deeper to the neighbor node
        insertGapsForInsertionEvents(index_mapping_by_jump_step, stopping_node_id, (*it)->node, node, stop_inserting_gaps);
    }
}

/**
    handle insertion events
*/
int AliSimulator::handleInsertion(int &sequence_length, vector<int> &index_mapping_by_jump_step, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &accummulated_rates, SIMULATION_METHOD simulation_method)
{
    // Randomly select the position/site (from the set of all sites) where the insertion event occurs based on a uniform distribution between 0 and the current length of the sequence
    int position = selectValidPositionForIndels(sequence_length + 1, indel_sequence);
    
    // Randomly generate the length (length_I) of inserted sites from the indel-length distribution (​​geometric distribution (by default) or user-defined distributions).
    int length = -1;
    for (int i = 0; i < 1000; i++)
    {
        length = generateIndelSize(params->alisim_insertion_distribution);
        
        // a valid length must be greater than 0
        if (length > 0)
            break;
    }
    // validate the length
    if (length <= 0)
        outError("Sorry! Could not generate a positive length (for insertion events) based on the insertion-distribution within 1000 attempts.");
    
    // insert new_sequence into the current sequence
    vector<short int> new_sequence = generateRandomSequence(length, false);
    insertNewSequenceForInsertionEvent(indel_sequence, position, new_sequence);
    
    // if RATE_MATRIX approach is used -> update total_sub_rate and accummulated_rates
    if (simulation_method == RATE_MATRIX)
    {
        // extract new_site_specific_rate
        vector<double> new_site_specific_rates;
        if (site_specific_rates.size() > 0)
        {
            for (int i = position; i < position + length; i++)
                new_site_specific_rates.push_back(site_specific_rates[i]);
        }
        
        // extract new_site_specific_model_index
        vector<short int> new_site_specific_model_index;
        if (site_specific_model_index.size() > 0)
        {
            for (int i = position; i < position + length; i++)
                new_site_specific_model_index.push_back(site_specific_model_index[i]);
        }
        
        // update total_sub_rate
        double diff = computeTotalSubRate(new_site_specific_rates, new_site_specific_model_index, new_sequence);
        total_sub_rate += diff;
        
        // update accummulated_rates
        // update accummulated_rates of the inserted sites
        accummulated_rates.insert(accummulated_rates.begin()+position, length, 0);
        for (int i = position; i < position + length; i++)
        {
            int mixture_index = site_specific_model_index.size() == 0? 0:site_specific_model_index[i];
            double site_rate = site_specific_rates.size() > 0?site_specific_rates[i]:1;
            double prior_acc_rate = i > 0?accummulated_rates[i-1]:0;
            accummulated_rates[i] = prior_acc_rate + site_rate*sub_rates[mixture_index*max_num_states + indel_sequence[i]];
        }
        // update accummulated_rates of the remaining sites
        for (int i = position + length; i < accummulated_rates.size(); i++)
            accummulated_rates[i] += diff;
    }
    
    // update the sequence_length
    sequence_length += length;
    
    // update index_mapping_by_jump_step
    for (int i = index_mapping_by_jump_step.size() - 1; i >= 0; i--)
    {
        if (i + index_mapping_by_jump_step[i] >= position)
            index_mapping_by_jump_step[i] += length;
        else
            break;
    }
    
    // return insertion-size
    return length;
}

/**
    handle deletion events
*/
int AliSimulator::handleDeletion(int sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &accummulated_rates, SIMULATION_METHOD simulation_method)
{
    // Randomly generate the length (length_D) of sites (which will be deleted) from the indel-length distribution.
    int length = -1;
    for (int i = 0; i < 1000; i++)
    {
        length = (int) generateIndelSize(params->alisim_deletion_distribution);
        
        // a valid length must be greater than 0
        if (length > 0)
            break;
    }
    // validate the length
    if (length <= 0)
        outError("Sorry! Could not generate a positive length (for deletion events) based on the deletion-distribution within 1000 attempts.");
    
    // Randomly select the position/site (from the set of all sites) where the deletion event occurs based on a uniform distribution between 0 and the current length of the sequence
    int position = 0;
    int upper_bound = sequence_length - length;
    if (upper_bound > 0)
        position = selectValidPositionForIndels(upper_bound, indel_sequence);
    
    // Replace up to length_D sites by gaps from the sequence starting at the selected location
    int real_deleted_length = 0;
    vector<short int> deleted_sites;
    vector<double> deleted_site_specific_rate;
    vector<short int> deleted_site_specific_model_index;
    int i;
    for (i = 0; i < length && (position + i) < indel_sequence.size(); i++)
    {
        // if the current site is not a gap (has not been deleted) -> replacing it by a gap
        if (indel_sequence[position + i ] != STATE_UNKNOWN)
        {
            // if RATE_MATRIX approach is used -> update deleted_sites and deleted_site_specific_rate to update total_sub_rate later
            if (simulation_method == RATE_MATRIX)
            {
                deleted_sites.push_back(indel_sequence[position + i]);
                if (site_specific_rates.size() > 0)
                    deleted_site_specific_rate.push_back(site_specific_rates[position + i]);
                if (site_specific_model_index.size() > 0)
                    deleted_site_specific_model_index.push_back(site_specific_model_index[position + i]);
            }
            
            // delete site
            indel_sequence[position + i ] = STATE_UNKNOWN;
            real_deleted_length++;
        }
        // otherwise, ignore the current site, moving forward to find a valid site (not a gap)
        else
        {
            i--;
            position++;
        }
        
        // if RATE_MATRIX approach is used -> update accummulated_rates
        if (simulation_method == RATE_MATRIX)
            accummulated_rates[position + i] = position + i > 0?accummulated_rates[position + i - 1]:0;
    }
    
    // if RATE_MATRIX approach is used -> update total_sub_rate, and accummulated_rates
    if (simulation_method == RATE_MATRIX)
    {
        // update total_sub_rate
        double diff = -computeTotalSubRate(deleted_site_specific_rate, deleted_site_specific_model_index, deleted_sites);
        total_sub_rate += diff;
        
        // update accummulated_rates
        for (int j = position + i; j < accummulated_rates.size(); j++)
            accummulated_rates[j] += diff;
    }
    
    // return deletion-size
    return real_deleted_length;
}

/**
    handle substitution events
*/
void AliSimulator::handleSubs(int sequence_length, double &total_sub_rate, vector<double> &accummulated_rates, vector<short int> &indel_sequence)
{
    // select a position where the substitution event occurs
    double random_cumulative_rate = random_double() * total_sub_rate;
    int pos = binarysearchItemWithAccumulatedProbabilityMatrix(accummulated_rates, random_cumulative_rate, 0, sequence_length - 1, 0);
    
    // extract the current state
    short int current_state = indel_sequence[pos];
    
    // estimate the new state
    int mixture_index = site_specific_model_index.size() == 0 ? 0:site_specific_model_index[pos];
    int starting_index = mixture_index*max_num_states*max_num_states + max_num_states*current_state;
    indel_sequence[pos] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(Jmatrix, starting_index, max_num_states, max_num_states/2);
    
    // update total_sub_rate
    double current_site_rate = site_specific_rates.size() == 0 ? 1 : site_specific_rates[pos];
    double diff = current_site_rate*(sub_rates[mixture_index*max_num_states + indel_sequence[pos]] - sub_rates[mixture_index*max_num_states + current_state]);
    total_sub_rate += diff;
    
    // update accummulated_rates
    if (diff != 0)
        for (int i = pos; i < accummulated_rates.size(); i++)
            accummulated_rates[i] += diff;
}

/**
*  randomly select a valid position (not a deleted-site) for insertion/deletion event
*
*/
int AliSimulator::selectValidPositionForIndels(int upper_bound, vector<short int> sequence)
{
    int position = -1;
    for (int i = 0; i < upper_bound; i++)
    {
        position = random_int(upper_bound);
        
        // a valid position must not be a deleted site
        if (position == sequence.size() || sequence[position] != STATE_UNKNOWN)
            break;
    }
    // validate the position
    if (position < sequence.size() && sequence[position] == STATE_UNKNOWN)
        outError("Sorry! Could not select a valid position (not a deleted-site) for insertion/deletion events. You may specify a too high deletion rate, thus almost all sites were deleted. Please try again a a smaller deletion ratio!");
    return position;
}

/**
    merge the simulated sequence with indel_sequence
*/
void AliSimulator::mergeIndelSequence(Node* node, vector<short int> indel_sequence, vector<int> index_mapping_by_jump_step)
{
    // mapping state from the current sequence into indel_sequence
    for (int i = 0; i < index_mapping_by_jump_step.size() - 1; i++)
    {
        // preserving deleted site
        if (indel_sequence[i+index_mapping_by_jump_step[i]] != STATE_UNKNOWN)
            indel_sequence[i+index_mapping_by_jump_step[i]] = node->sequence[i+index_mapping_by_jump_step[i]];
    }
    
    // release the memory for the current sequence
    vector<short int>().swap(node->sequence);
    
    // update the new sequence
    node->sequence = indel_sequence;
}

/**
    generate indel-size from its distribution
*/
int AliSimulator::generateIndelSize(IndelDistribution indel_dis)
{
    int random_size = -1;
    switch (indel_dis.indel_dis_type)
    {
        case NEG_BIN:
            random_size = random_int_nebin(indel_dis.param_1, indel_dis.param_2);
            break;
        case ZIPF:
            random_size = random_int_zipf(indel_dis.param_1, indel_dis.param_2);
            break;
        case LAV:
            random_size = random_int_lav(indel_dis.param_1, indel_dis.param_2);
            break;
        case GEO:
            random_size = random_int_geometric(indel_dis.param_1) + 1;
            break;
        default:
            random_size = random_number_from_distribution(indel_dis.user_defined_dis);
            break;
    }
    return random_size;
}

/**
    compute mean of deletion-size
*/
double AliSimulator::computeMeanDelSize(int sequence_length)
{
    // only compute mean of deletion-size if it has not yet computed
    if (params->alisim_mean_deletion_size == -1)
    {
        // dummy variables
        int total = 0;
        int num_success = 0;
        
        // randomly generate sequence_length random deletion-sizes from the deletion-distribution
        for (int i = 0; i < sequence_length; i++)
        {
            int random_size = generateIndelSize(params->alisim_insertion_distribution);
            
            // only process positive random sizes
            if (random_size > 0)
            {
                total += random_size;
                num_success++;
            }
        }
        
        // make sure we could generate positive size
        if (num_success == 0)
            outError("Could not generate positive deletion-sizes from the deletion-distribution. Please check and try again!");
        else
            params->alisim_mean_deletion_size = (double)total/num_success;
    }
    
    return params->alisim_mean_deletion_size;
}

/**
    root tree
*/
void AliSimulator::rootTree()
{
    // dummy variables
    Node* new_root = new Node();
    Node* second_internal_node;
    
    // extract the intermediate node
    if (tree && tree->root && tree->root->neighbors.size() > 0)
        second_internal_node = tree->root->neighbors[0]->node;
    if (second_internal_node)
    {
        // update new_root
        new_root->name = ROOT_NAME;
        new_root->id = tree->leafNum;
        new_root->sequence = tree->root->sequence;
        
        // change the id of the intermediate node if it's equal to the root's id
        if (second_internal_node->id == new_root->id)
            second_internal_node->id = (new_root->id * 10);
        
        // link new_root node with the intermediate node
        new_root->addNeighbor(second_internal_node, 0);
        second_internal_node->addNeighbor(new_root, 0);
        
        // update related info
        tree->root = new_root;
        tree->rooted = true;
        tree->leafNum++;
    }
}
