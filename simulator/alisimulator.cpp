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
    tree->initSequences();
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    STATE_UNKNOWN = tree->aln->STATE_UNKNOWN;
    max_num_states = tree->aln->getMaxNumStates();
    latest_insertion = NULL;
    first_insertion = NULL;
    
    // estimating the appropriate length_ratio in cases models with +ASC
    estimateLengthRatio();
    
    if (expected_number_sites == -1)
        expected_num_sites = (num_sites_per_state == 1 ? round(params->alisim_sequence_length * length_ratio) : round(params->alisim_sequence_length / num_sites_per_state * length_ratio));
    else
        expected_num_sites = round(expected_number_sites * length_ratio);
    partition_rate = new_partition_rate;
    
    // check if base frequencies for DNA models are specified correctly
    checkBaseFrequenciesDNAModels(tree, params->model_name);
    
    // extract max length of taxa names
    extractMaxTaxaNameLength();
    
    // innialize set of selected sites for permutation in FunDi model
    if (params->alisim_fundi_taxon_set.size()>0)
        selectAndPermuteSites(fundi_items, params->alisim_fundi_proportion, round(expected_num_sites));
}

/**
    constructor
*/
AliSimulator::AliSimulator(Params *input_params, IQTree *iq_tree, int expected_number_sites, double new_partition_rate)
{
    params = input_params;
    tree = iq_tree;
    tree->initSequences();
    num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
    STATE_UNKNOWN = tree->aln->STATE_UNKNOWN;
    max_num_states = tree->aln->getMaxNumStates();
    latest_insertion = NULL;
    first_insertion = NULL;
    
    // estimating the appropriate length_ratio in cases models with +ASC
    estimateLengthRatio();
    
    if (expected_number_sites == -1)
        expected_num_sites = (num_sites_per_state == 1 ? round(params->alisim_sequence_length * length_ratio) : round(params->alisim_sequence_length / num_sites_per_state * length_ratio));
    else
        expected_num_sites = round(expected_number_sites * length_ratio);
    partition_rate = new_partition_rate;
    
    // extract max length of taxa names
    extractMaxTaxaNameLength();
    
    // innialize set of selected sites for permutation in FunDi model
    if (params->alisim_fundi_taxon_set.size()>0)
        selectAndPermuteSites(fundi_items, params->alisim_fundi_proportion, round(expected_num_sites));
}

AliSimulator::~AliSimulator()
{
    // delete first_insertion
    if (first_insertion)
    {
        delete first_insertion;
        first_insertion = NULL;
    }
    
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
        
        // extract names of all taxa in the super tree if topology-unlinked partition is being used
        vector<string> super_taxa_names;
        if (params->partition_type == BRLEN_OPTIMIZE)
            tree->getTaxaName(super_taxa_names);
        
        // compute super_tree_length
        double inverse_super_tree_length = 1.0 / ((PhyloSuperTree*) tree)->treeLength();
        
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
                tree_line_index = i;
                // show information for the first time
                if (i == 0)
                {
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
            extractSiteID(current_tree->aln, current_tree->aln->position_spec.c_str(), siteIDs, false, -1, true);
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
                    ((PhyloSuperTree*) tree)->part_info[i].part_rate = current_tree_length * inverse_super_tree_length;
                
                // update sum of rate*n_sites and num_sites (for rate normalization)
                sum += ((PhyloSuperTree*) tree)->part_info[i].part_rate * current_tree->aln->getNSite();
                if (current_tree->aln->seq_type == SEQ_CODON && ((PhyloSuperTree*) tree)->rescale_codon_brlen)
                    num_sites += 3 * current_tree->aln->getNSite();
                else
                    num_sites += current_tree->aln->getNSite();
            }
            
            // add missing taxa from the current partition tree to the super tree if topology-unlink partition is used
            if (params->partition_type == BRLEN_OPTIMIZE && i > 0)
            {
                vector<string> taxa_names;
                current_tree->getTaxaName(taxa_names);
                
                // iteratively check taxa (in the current tree) exist in the super tree -> if not -> adding new taxon into the super tree
                for (string name: taxa_names)
                    if (std::find(super_taxa_names.begin(), super_taxa_names.end(), name) == super_taxa_names.end())
                    {
                        
                        // find a leaf
                        ASSERT(super_taxa_names.size() > 0);
                        int leaf_index = super_taxa_names.size() > 1?1:0;
                        Node* leaf = tree->findLeafName(super_taxa_names[leaf_index]);
                        if (!leaf || leaf->neighbors.size() == 0) continue;
                        
                        // extract leaf's dad
                        Node* dad = leaf->neighbors[0]->node;
                        
                        // init an internal node
                        Node *internal = new Node();
                    
                        // init a new node for the new taxon
                        Node *new_taxon = new Node();
                        new_taxon->name = name;
                        
                        // update neighbor of dad
                        dad->updateNeighbor(leaf, internal, 0);
                        leaf->updateNeighbor(dad, internal, 0);
                        
                        // add neighbors to internal
                        internal->addNeighbor(dad, 0);
                        internal->addNeighbor(leaf, 0);
                        
                        // add connection between the internal and the new taxon
                        internal->addNeighbor(new_taxon, 0);
                        new_taxon->addNeighbor(internal, 0);
                        
                        // update super_taxa_names
                        super_taxa_names.push_back(name);
                    }
                
                // update the super tree if new taxa were added
                if (tree->leafNum != super_taxa_names.size())
                {
                    tree->leafNum = super_taxa_names.size();
                    tree->nodeNum = tree->leafNum;
                    tree->initializeTree();
                }
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
                
                // manually detect AA data from "NONREV", "GTR20", "Poisson" model and detect DNA data from UNREST model
                if (tree->aln->seq_type == SEQ_UNKNOWN)
                {
                    const char* aa_model_names_plus[] = {"NONREV", "GTR20", "Poisson"};
                    std::transform(model_familyname.begin(), model_familyname.end(), model_familyname.begin(), ::toupper);
                    
                    for (int i = 0; i < sizeof(aa_model_names_plus)/sizeof(char*); i++)
                        if (model_familyname == aa_model_names_plus[i]) {
                            tree->aln->seq_type = SEQ_PROTEIN;
                            break;
                        }
                    
                    // manually detect DNA data from UNREST model
                    if (tree->aln->seq_type == SEQ_UNKNOWN && model_familyname == "UNREST")
                        tree->aln->seq_type = SEQ_DNA;
                }
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
            
            // throw error if users dont specify the number of states when simulating morph data
            if (tree->aln->num_states <= 0)
                outError("Please specify the number of states for morphological data by --seqtype MORPH{<NUM_STATES>}");
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
    delete models_block;
}

/**
    remove all constant sites (in case with +ASC)
*/
void AliSimulator::removeConstantSites(){
    // dummy variables
    int num_variant_states = -1;
    vector<short int> variant_state_mask;
    int rounded_length = round(expected_num_sites * inverse_length_ratio);
    
    // create a variant state mask
    createVariantStateMask(variant_state_mask, num_variant_states, rounded_length, tree->root, tree->root);
    
    // return error if num_variant_states is less than the expected_num_variant_states
    if (num_variant_states < rounded_length){
        outError("Unfortunately, after removing constant sites, the number of variant sites is less than the expected sequence length. Please use --length-ratio <LENGTH_RATIO> to generate more abundant sites and try again. The current <LENGTH_RATIO> is "+ convertDoubleToString(length_ratio));
    }
    
    // if using Indels, update seq_length_indels
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
        seq_length_indels = num_variant_states;

    // recording start_time
    auto start = getRealTime();
    
    // get only variant sites for leaves
    getOnlyVariantSites(variant_state_mask, tree->root, tree->root);
    
    // show the time spent on copy variant sites
    auto end = getRealTime();
    cout<<" - Time spent on copying only variant sites: "<<end-start<<endl;
}

/**
    only get variant sites
*/
void AliSimulator::getOnlyVariantSites(vector<short int> &variant_state_mask, Node *node, Node *dad){
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // dummy sequence
        vector<short int> variant_sites(variant_state_mask.size(),0);
        
        // initialize the number of variant sites
        int num_variant_states = 0;
        int rounded_length = round(expected_num_sites * inverse_length_ratio);
        
        // browse sites one by one
        for (int i = 0; i < node->sequence->sequence_chunks[0].size(); i++)
            // only get variant sites
            if (variant_state_mask[i] == -1)
            {
                // get the variant site
                variant_sites[num_variant_states] = node->sequence->sequence_chunks[0][i];
                num_variant_states++;
                
                // stop checking further states if num_variant_states has exceeded the expected_num_variant_states
                // keep checking if Indels is used
                if (num_variant_states >= rounded_length && params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0)
                    break;
            }
        
        // replace the sequence of the Leaf by variant sites
        node->sequence->sequence_chunks[0].clear();
        variant_sites.resize(num_variant_states);
        node->sequence->sequence_chunks[0] = variant_sites;
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
void AliSimulator::generatePartitionAlignment(vector<short int> &ancestral_sequence, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode)
{
    // reset number of chunks of the root sequence to 1
    tree->MTree::root->sequence->sequence_chunks.resize(1);
    
    // if the ancestral sequence is not specified, randomly generate the sequence
    if (ancestral_sequence.size() == 0)
        generateRandomSequence(expected_num_sites, tree->MTree::root->sequence->sequence_chunks[0]);
    // otherwise, using the ancestral sequence + abundant sites
    else
    {
        // set the ancestral sequence to the root node
        tree->MTree::root->sequence->sequence_chunks[0] = ancestral_sequence;
        
        // add abundant_sites
        int num_abundant_sites = expected_num_sites - ancestral_sequence.size();
        if (num_abundant_sites > 0)
        {
            vector<short int> abundant_sites;
            generateRandomSequence(num_abundant_sites, abundant_sites);
            for (int site:abundant_sites)
                tree->MTree::root->sequence->sequence_chunks[0].push_back(site);
        }
    }
    
    // validate the sequence length (in case of codon)
    validataSeqLengthCodon();
    
    // simulate the sequence for each node in the tree by DFS
    simulateSeqsForTree(input_msa, output_filepath, open_mode);
}

/**
    create mask for variant sites
*/
void AliSimulator::createVariantStateMask(vector<short int> &variant_state_mask, int &num_variant_states, int expected_num_variant_states, Node *node, Node *dad){
    // no need to check the further sites if num_variant_states has exceeded the expected_num_variant_states
    // keep checking if Indels is used
    if (num_variant_states >= expected_num_variant_states && params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0)
        return;
    
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // initialize the mask (all sites are assumed to be constant)
        if (num_variant_states == -1)
        {
            num_variant_states = 0;
            variant_state_mask = node->sequence->sequence_chunks[0];
        }
        // otherwise, check state by state to update the mask
        else
        {
            for (int i = 0; i < node->sequence->sequence_chunks[0].size(); i++)
            {
                if (variant_state_mask[i] != -1 && variant_state_mask[i] != node->sequence->sequence_chunks[0][i] && node->sequence->sequence_chunks[0][i] != STATE_UNKNOWN)
                {
                    // if variant_state_mask is a gap -> update it by the current site of the sequence
                    if (variant_state_mask[i] == STATE_UNKNOWN)
                        variant_state_mask[i] = node->sequence->sequence_chunks[0][i];
                    // otherwise, mask the variant_state_mask of the current site as variant
                    else
                    {
                        // if the current state is changed -> increase num_variant_states, and disable that state
                        variant_state_mask[i] = -1;
                        num_variant_states++;
                        
                        // stop checking further states if num_variant_states has exceeded the expected_num_variant_states
                        // keep checking if Indels is used
                        if (num_variant_states >= expected_num_variant_states && params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0)
                            break;
                    }
                }
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
void AliSimulator::generateRandomSequence(int sequence_length, vector<short int> &sequence, bool initial_freqs)
{
    // if the Frequency Type is FREQ_EQUAL -> randomly generate each site in the sequence follows the normal distribution
    if (tree->getModel()->getFreqType() == FREQ_EQUAL)
    {
        // initialize sequence
        sequence.resize(sequence_length);
        
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
        generateRandomSequenceFromStateFreqs(sequence_length, sequence, state_freq, max_prob_pos);
        
        // delete state_freq
        delete []  state_freq;
    }
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
        tree->getModel()->decomposeRateMatrix();
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
    sum = 1.0 / sum;
    for (int i = 0; i < max_num_states; i++)
        base_frequencies[i] *= sum;
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulator::simulateSeqsForTree(map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode)
{
    // get variables
    int sequence_length = expected_num_sites;
    ModelSubst *model = tree->getModel();
    vector<string> state_mapping;
    int default_segment_length;
    int thread_id = 0;
    bool write_sequences_to_tmp_data = false;
    bool store_seq_at_cache = true;
    int max_depth = 0;
    
    // init variables
    initVariables(sequence_length, output_filepath, state_mapping, model, default_segment_length, max_depth, write_sequences_to_tmp_data, store_seq_at_cache);
    
    // execute one of the AliSim-OpenMP algorithms to simulate sequences
    if (params->alisim_openmp_alg == IM)
        executeIM(thread_id, sequence_length, default_segment_length, model, input_msa, output_filepath, open_mode, write_sequences_to_tmp_data, store_seq_at_cache, max_depth, state_mapping);
    else
        executeEM(thread_id, sequence_length, default_segment_length, model, input_msa, output_filepath, open_mode, write_sequences_to_tmp_data, store_seq_at_cache, max_depth, state_mapping);
    
    // process after simulating sequences
    postSimulateSeqs(sequence_length, output_filepath, write_sequences_to_tmp_data);
}

void AliSimulator::executeEM(int thread_id, int &sequence_length, int default_segment_length, ModelSubst *model, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data, bool store_seq_at_cache, int max_depth, vector<string> &state_mapping)
{
    ostream *single_output = NULL;
    ostream *out = NULL;
    int *rstream = NULL;
    vector<vector<short int>> sequence_cache;
    int actual_segment_length = sequence_length;
    
    // simulate Sequences
    #ifdef _OPENMP
    #pragma omp parallel private(rstream, out, thread_id, sequence_cache, actual_segment_length)
    {
        thread_id = omp_get_thread_num();
        // init random generators
        int ran_seed = params->ran_seed + MPIHelper::getInstance().getProcessID() * 1000 + thread_id;
        init_random(ran_seed, false, &rstream);

        actual_segment_length = thread_id < num_simulating_threads - 1 ? default_segment_length : sequence_length - (num_simulating_threads - 1) * default_segment_length;
    #endif
        // init sequence cache
        if (store_seq_at_cache)
        {
            sequence_cache.resize(max_depth + 1);
            for (int i = 1; i < max_depth + 1; i++)
                sequence_cache[i].resize(actual_segment_length);

            // cache sequence at root
            sequence_cache[0] = tree->root->sequence->sequence_chunks[thread_id];
        }
        
        // init the output stream
        initOutputFile(out, thread_id, actual_segment_length, output_filepath, open_mode, write_sequences_to_tmp_data);
        
        // initialize trans_matrix
        double *trans_matrix = new double[max_num_states * max_num_states];
        simulateSeqs(thread_id, thread_id * default_segment_length, actual_segment_length, sequence_length, model, trans_matrix, sequence_cache, store_seq_at_cache, tree->MTree::root, tree->MTree::root, *out, state_mapping, input_msa, rstream);
        
        // delete trans_matrix array
        delete[] trans_matrix;
        
        // close the output stream
        if (output_filepath.length() > 0 || write_sequences_to_tmp_data)
            closeOutputStream(out, num_threads > 1);
        
        // release sequence cache
        if (store_seq_at_cache)
            vector<vector<short int>>().swap(sequence_cache);
        
        // release mem for rstream
        finish_random(rstream);
        
        // merge output files into a single file if using multithreading
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        mergeOutputFiles(single_output, thread_id, output_filepath, open_mode, write_sequences_to_tmp_data);
            
    #ifdef _OPENMP
    }
    #endif
}

/**
    merge output files
*/
void AliSimulator::mergeOutputFiles(ostream *&single_output, int thread_id, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data)
{
    // skip merging if users want to do so
    if (params->no_merge) return;
    
    if (output_filepath.length() > 0 && !write_sequences_to_tmp_data)
    {
        // merge output files into a single file
        if (num_threads > 1)
        {
            // open single_output stream
            #ifdef _OPENMP
            #pragma omp single
            #endif
            {
                string single_output_filepath;
                if (params->aln_output_format != IN_FASTA)
                    single_output_filepath = output_filepath + ".phy";
                else
                    single_output_filepath = output_filepath + ".fa";
                openOutputStream(single_output, single_output_filepath, open_mode);
                
                // output the first line
                string first_line = "";
                if (params->aln_output_format != IN_FASTA)
                {
                    int num_leaves = tree->leafNum - ((tree->root->isLeaf() && tree->root->name == ROOT_NAME)?1:0);
                    first_line = convertIntToString(num_leaves) + " " + convertIntToString(num_sites_per_state == 1 ? round(expected_num_sites * inverse_length_ratio) : (round(expected_num_sites * inverse_length_ratio) * num_sites_per_state)) + "\n";
                    *single_output << first_line;
                }
                if (!params->do_compression)
                    starting_pos = single_output->tellp();
                else
                    starting_pos = first_line.length();
            }
            
            // dummy variables
            vector<ifstream> input_streams(num_threads);
            int actual_num_lines_per_thread = 0;
            int average_num_lines_per_thread = 0;
            int line_num = 0;
            string line;
            uint64_t pos = 0;
            double inverse_num_threads = 1.0 / num_threads;
            
            // open all files
            for (int i = 0; i < input_streams.size(); i++)
            {
                // add ".phy" or ".fa" to the output_filepath
                string tmp_output_filepath;
                if (params->aln_output_format != IN_FASTA)
                    tmp_output_filepath = output_filepath + "_" + convertIntToString(i + 1) + ".phy";
                else
                    tmp_output_filepath = output_filepath + "_" + convertIntToString(i + 1) + ".fa";
                
                // open an input file
                input_streams[i].open(tmp_output_filepath.c_str(), std::ifstream::binary);
                
                // count num of lines, set start position
                safeGetline(input_streams[i], line);
                int line_length = line.length() + 1;
                // when using FunDi, the entire sequences are written into one of output files -> need to check average_num_lines_per_thread == 0
                if (i == 0 || average_num_lines_per_thread == 0)
                {
                    input_streams[i].seekg(0, input_streams[i].end);
                    uint64_t length = input_streams[i].tellg();
                    int total_num_lines = length/line_length;
                    average_num_lines_per_thread = total_num_lines * inverse_num_threads;
                    
                    // compute the starting line number
                    pos = thread_id * average_num_lines_per_thread;
                    
                    // compute actual_num_lines_per_thread
                    actual_num_lines_per_thread = thread_id < (num_threads - 1) ? average_num_lines_per_thread : total_num_lines - pos;
                    
                    
                }
                // set start position
                input_streams[i].seekg(pos * line_length);
            }
            
            // compute pos in the output file
            pos = starting_pos + (pos - 1) * output_line_length;
            
            // merge chunks of each sequence and write to the single output file
            for (line_num = 0; line_num < actual_num_lines_per_thread; line_num++)
            {
                string output = "";
                for (int j = 0; j < input_streams.size(); j++)
                {
                    safeGetline(input_streams[j], line);
                    output += line;
                }
                
                if (output.length() > 0)
                {
                    // if using FASTA -> break sequence name into new line
                    if (params->aln_output_format == IN_FASTA)
                    {
                        // insert '>' into sequence name
                        output = ">" + output;
                        
                        // insert break line
                        for (int i = 1; i < output.length(); i++)
                            if (output[i-1] == ' ' && output[i] != ' ')
                            {
                                output[i-1] = '\n';
                                break;
                            }
                    }
                    
                    // add the break line for the current output
                    output += "\n";
                    
                    // update pos for the current output line
                    pos += output_line_length;
                
                    // write the concatenated sequence into file
                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                        // jump to the correct position before writing if users want to keep the sequence order
                        if (params->keep_seq_order)
                            single_output->seekp(pos);
                        (*single_output) << output;
                    }
                }
            }
            
            // close all files
            for (int i = 0; i < input_streams.size(); i++)
                input_streams[i].close();
            
            // close single_output stream
            #ifdef _OPENMP
            #pragma omp barrier
            #pragma omp single
            #endif
            closeOutputStream(single_output);
            
            // delete all intermidate files
            // add ".phy" or ".fa" to the output_filepath
            string tmp_output_filepath;
            if (params->aln_output_format != IN_FASTA)
                tmp_output_filepath = output_filepath + "_" + convertIntToString(thread_id + 1) + ".phy";
            else
                tmp_output_filepath = output_filepath + "_" + convertIntToString(thread_id + 1) + ".fa";
            // delete file
            remove(tmp_output_filepath.c_str());
        }
        
        // show the output file name
        #ifdef _OPENMP
        #pragma omp single
        #endif
        {
            string single_output_filepath;
            if (params->aln_output_format != IN_FASTA)
                single_output_filepath = output_filepath + ".phy";
            else
                single_output_filepath = output_filepath + ".fa";
            cout << "An alignment has just been exported to " << single_output_filepath << endl;
        }
    }
}

void AliSimulator::executeIM(int thread_id, int &sequence_length, int default_segment_length, ModelSubst *model, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data, bool store_seq_at_cache, int max_depth, vector<string> &state_mapping)
{
    int actual_segment_length = sequence_length;
    ostream *out = NULL;
    int *rstream = NULL;
    vector<vector<short int>> sequence_cache;
    
    // init the output stream
    initOutputFile(out, thread_id, actual_segment_length, output_filepath, open_mode, write_sequences_to_tmp_data);
    
    // simulate Sequences
    #ifdef _OPENMP
    #pragma omp parallel private(rstream, thread_id, sequence_cache, actual_segment_length)
    {
        thread_id = omp_get_thread_num();
        // init random generators
        int ran_seed = params->ran_seed + MPIHelper::getInstance().getProcessID() * 1000 + thread_id;
        init_random(ran_seed, false, &rstream);
            
        actual_segment_length = thread_id < num_simulating_threads - 1 ? default_segment_length : sequence_length - (num_simulating_threads - 1) * default_segment_length;
    #endif
        // init sequence cache
        if (store_seq_at_cache)
        {
            // don't need to init sequence_cache for the writing thread
            if (!(num_threads > 1 && thread_id == num_threads - 1))
            {
                sequence_cache.resize(max_depth + 1);
                for (int i = 1; i < max_depth + 1; i++)
                    sequence_cache[i].resize(actual_segment_length);
                
                // cache sequence at root
                sequence_cache[0] = tree->root->sequence->sequence_chunks[thread_id];
            }
            
            // init common cache of writing queue
            #ifdef _OPENMP
            #pragma omp single
            #endif
            {
                if (num_threads > 1)
                {
                    // default cache_size_per_thread = num_simulating_threads * 2;
                    cache_size_per_thread = tree->params->mem_limit_factor == 0 ? num_simulating_threads * 2 : ceil(tree->leafNum * tree->params->mem_limit_factor);
                    seq_str_cache.resize(num_simulating_threads * cache_size_per_thread);
                    cache_start_indexes.resize(num_simulating_threads);
                    cache_start_indexes[0] = 0;
                    for (int thread_id = 1; thread_id < num_simulating_threads; thread_id++)
                        cache_start_indexes[thread_id] = cache_start_indexes[thread_id - 1] + cache_size_per_thread;
                }
            }
        }
        
        // writing thread
        if (num_threads > 1 && thread_id == num_threads - 1 && store_seq_at_cache)
        {
            while (num_thread_done < num_simulating_threads)
            {
                // check and write chunks in the cache to files
                writeSeqChunkFromCache(out);
            }
            
            // final round to check and write all chunks in the cache to files
            #ifdef _OPENMP
            #pragma omp flush
            #endif
            writeAllSeqChunkFromCache(out);
            
        }
        // simulating thread
        else
        {
            // initialize trans_matrix
            double *trans_matrix = new double[max_num_states * max_num_states];
            simulateSeqs(thread_id, thread_id * default_segment_length, actual_segment_length, sequence_length, model, trans_matrix, sequence_cache, store_seq_at_cache, tree->MTree::root, tree->MTree::root, *out, state_mapping, input_msa, rstream);
            
            // delete trans_matrix array
            delete[] trans_matrix;
            
            // report simulation done for the current thread
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            num_thread_done++;
        }
        
        // release sequence cache
        if (store_seq_at_cache)
            vector<vector<short int>>().swap(sequence_cache);
        
        // release mem for rstream
        finish_random(rstream);
            
    #ifdef _OPENMP
    }
    #endif
    
    // close the output stream
    if (output_filepath.length() > 0 || write_sequences_to_tmp_data)
        closeOutputStream(out);
}

void AliSimulator::writeSeqChunkFromCache(ostream *&out)
{
    for (int thread_id = 0; thread_id < num_simulating_threads; thread_id++)
    {
        for (int i = cache_start_indexes[thread_id]; i < cache_start_indexes[thread_id] + cache_size_per_thread; i++)
        {
            if (seq_str_cache[i].chunk_status == OCCUPIED)
            {
                #ifdef _OPENMP
                #pragma omp flush
                #endif
                out->seekp(seq_str_cache[i].pos);
                (*out) << seq_str_cache[i].chunk_str;
                
                // update status of the selected slot
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                seq_str_cache[i].chunk_status = EMPTY;
                
                // round robin -> move to slots of other threads
                break;
            }
        }
    }
}

void AliSimulator::writeAllSeqChunkFromCache(ostream *&out)
{
    for (int i = 0; i < seq_str_cache.size(); i++)
    {
        if (seq_str_cache[i].chunk_status == OCCUPIED)
        {
            #ifdef _OPENMP
            #pragma omp flush
            #endif
            out->seekp(seq_str_cache[i].pos);
            (*out) << seq_str_cache[i].chunk_str;
            
            // update status of the selected slot
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            seq_str_cache[i].chunk_status = EMPTY;
        }
    }
}

/**
    process after simulating sequences
*/
void AliSimulator::postSimulateSeqs(int sequence_length, string output_filepath, bool write_sequences_to_tmp_data)
{
    // delete sub_rates, J_Matrix
    delete[] sub_rates;
    delete[] Jmatrix;
    
    // merge chunks into a single sequence if using multiple threads and sequences have not been outputted
    if (num_threads > 1 && (output_filepath.length() == 0 || write_sequences_to_tmp_data))
        mergeChunksAllNodes();
    
    // record the actual (final) seq_length due to Indels
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
        seq_length_indels = sequence_length;
    
    // process delayed Fundi if it is delayed due to Insertion events
    if (params->alisim_fundi_taxon_set.size()>0 && params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
    {
        // update new genome at tips from original genome and the genome tree
        updateNewGenomeIndels(seq_length_indels);
        
        processDelayedFundi(tree->root, tree->root);
    }
    
    // removing constant states if it's necessary
    if (length_ratio > 1)
    {
        // if using Indels, update new genome at tips from original genome and the genome tree
        // skip updating if using Fundi as it must be already updated by Fundi
        if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0 && params->alisim_fundi_taxon_set.size() == 0)
            updateNewGenomeIndels(seq_length_indels);
        
        removeConstantSites();
    }
}

/**
    initialize variables
*/
void AliSimulator::initVariables(int sequence_length, string output_filepath, vector<string> &state_mapping, ModelSubst *model, int &default_segment_length, int &max_depth, bool &write_sequences_to_tmp_data, bool &store_seq_at_cache)
{
    // check if we can store sequences at a fixed cache instead of at nodes
    store_seq_at_cache = params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0 && (output_filepath.length() > 0 || write_sequences_to_tmp_data) && params->alisim_fundi_taxon_set.size() == 0;
    
    // get the number of threads
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp single
    {
        num_threads = omp_get_num_threads();
        num_simulating_threads = num_threads;
        if (params->alisim_openmp_alg == IM && num_threads > 1 && store_seq_at_cache)
            num_simulating_threads = num_threads - 1;
    }
    #endif
    
    default_segment_length = sequence_length / num_simulating_threads;
    seq_name_length = max_length_taxa_name + (params->aln_output_format == IN_FASTA ? 1 : 0);
    output_line_length = round(expected_num_sites * inverse_length_ratio);
    output_line_length = num_sites_per_state == 1 ? (output_line_length + 1 + seq_name_length) : (output_line_length * num_sites_per_state + 1 + seq_name_length);
    
    // reset num_thread_done
    num_thread_done = 0;
    
    // check to use Posterior Mean Rates
    if (tree->params->alisim_rate_heterogeneity!=UNSPECIFIED)
        applyPosRateHeterogeneity = canApplyPosteriorRateHeterogeneity();
    
    // init site to pattern id if the user supplies an input alignment
    if (tree->params->alisim_inference_mode)
        initSite2PatternID(sequence_length);
    
    // initialize variables (site_specific_rates; site_specific_rate_index; site_specific_model_index)
    initVariablesRateHeterogeneity(sequence_length, true);
    
    // check whether we could temporarily write sequences at tips to tmp_data file => a special case: with Indels without FunDi/ASC/Partitions
    write_sequences_to_tmp_data = params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0 && params->alisim_fundi_taxon_set.size() == 0 && length_ratio <= 1 && !params->partition_file;
    
    // initialize state_mapping (mapping from state to characters)
    if (output_filepath.length() > 0 || write_sequences_to_tmp_data)
        initializeStateMapping(num_sites_per_state, tree->aln, state_mapping);
    
    // rooting the tree if it's unrooted
    if (!tree->rooted)
        rootTree();
    
    // compute the switching param to switch between Rate matrix and Probability matrix
    computeSwitchingParam(expected_num_sites);
    
    // initialize sub_rates, J_Matrix from Q_matrix
    int num_mixture_models = model->getNMixtures();
    int num_mixtures_times_num_states = num_mixture_models * max_num_states;
    sub_rates = new double[num_mixtures_times_num_states];
    Jmatrix = new double[num_mixtures_times_num_states * max_num_states];
    extractRatesJMatrix(model);
    
    // init genome_tree, and the initial empty insertion for root if using Indels
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
    {
        // init an empty insertion event
        if (first_insertion)
            delete first_insertion;
        first_insertion = new Insertion();
        latest_insertion = first_insertion;
        
        // init the insertion position for root if it is a leaf (a rooted tree)
        if (tree->root->isLeaf())
            tree->root->sequence->insertion_pos = latest_insertion;
    }
    
    // count the number of gaps at root if Indels is used
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
        tree->root->sequence->num_gaps = count(tree->root->sequence->sequence_chunks[0].begin(), tree->root->sequence->sequence_chunks[0].end(), STATE_UNKNOWN);
    
    // reset variables at nodes (essential when simulating multiple alignments)
    resetTree(max_depth, store_seq_at_cache);
    
    // if using AliSim-OpenMP-EM algorithm, update whether we need to output temporary files in PHYLIP format
    force_output_PHYLIP = params->alisim_openmp_alg == EM && num_threads > 1 && !params->no_merge;
}

/**
    init the output file
*/
void AliSimulator::initOutputFile(ostream *&out, int thread_id, int actual_segment_length, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data)
{
    if (output_filepath.length() > 0 || write_sequences_to_tmp_data)
    {
        // init an output_filepath to temporarily output the sequences (when simulating Indels)
        if (write_sequences_to_tmp_data)
        {
            output_filepath = params->alisim_output_filename + "_" + params->tmp_data_filename + "_" + convertIntToString(MPIHelper::getInstance().getProcessID());
            
            // open the output stream (create new file)
            openOutputStream(out, output_filepath, std::ios_base::out);
        }
        // otherwise, just add ".phy" or ".fa" to the output_filepath
        else
        {
            // only add thread_id to filename if using multiple threads with AliSim-OpenMP-EM
            string thread_id_str = "";
            if (params->alisim_openmp_alg == EM && num_threads > 1)
                thread_id_str = "_" + convertIntToString(thread_id + 1);
            
            // add ".phy" or ".fa" to the output_filepath
            if (params->aln_output_format != IN_FASTA)
                output_filepath = output_filepath + thread_id_str + ".phy";
            else
                output_filepath = output_filepath + thread_id_str + ".fa";
            
            // open the output stream (create new or append an existing file)
            if (params->alisim_openmp_alg == EM && num_threads > 1)
                openOutputStream(out, output_filepath, std::ios_base::out, true);
            else
                openOutputStream(out, output_filepath, open_mode);
        }
        
        // write the first line <#taxa> <length_of_sequence> (for PHYLIP output format)
        string first_line = "";
        if (params->aln_output_format == IN_PHYLIP)
        {
            int num_leaves = tree->leafNum - ((tree->root->isLeaf() && tree->root->name == ROOT_NAME)?1:0);
            
            // if using AliSim-OpenMP-EM algorithm
            if (params->alisim_openmp_alg == EM)
            {
                // only output the first line in the singlethreading mode;
                // in multithreading mode:
                // -> without merging intermediate output files -> also output the first line
                // -> with merging step -> the first line will be output later when merging output files
                if (num_threads == 1
                    || (num_threads > 1 && params->no_merge))
                    *out << num_leaves << " " << round(actual_segment_length * inverse_length_ratio) * num_sites_per_state << endl;
            }
            // if using AliSim-OpenMP-IM algorithm
            else
            {
                first_line = convertIntToString(num_leaves) + " " + convertIntToString(round(expected_num_sites * inverse_length_ratio) * num_sites_per_state) + "\n";
                *out << first_line;
            }
        }
        
        // get the starting position for writing
        if (params->alisim_openmp_alg == IM)
        {
            if (!params->do_compression)
                starting_pos = out->tellp();
            else
                starting_pos = first_line.length();
        }
    }
}

/**
    open an output stream
*/
void AliSimulator::openOutputStream(ostream *&out, string output_filepath, std::ios_base::openmode open_mode, bool force_uncompression)
{
    try {
        if (params->do_compression && !force_uncompression)
            out = new ogzstream(output_filepath.c_str(), open_mode);
        else
            out = new ofstream(output_filepath.c_str(), open_mode);
        out->exceptions(ios::failbit | ios::badbit);
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, output_filepath);
    }
}

/**
    close an output stream
*/
void AliSimulator::closeOutputStream(ostream *&out, bool force_uncompression)
{
    if (params->do_compression && !force_uncompression)
        ((ogzstream*)out)->close();
    else
        ((ofstream*)out)->close();
    delete out;
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulator::simulateSeqs(int thread_id, int segment_start, int &segment_length, int &sequence_length, ModelSubst *model, double *trans_matrix, vector<vector<short int>> &sequence_cache, bool store_seq_at_cache, Node *node, Node *dad, ostream &out, vector<string> &state_mapping, map<string,string> input_msa, int* rstream)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        //  clone the number of gaps from the ancestral sequence if using Indels
        if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
            (*it)->node->sequence->num_gaps = node->sequence->num_gaps;
        
        // get dad_seq_chunk and node_seq_chunk
        vector<short int> *dad_seq_chunk, *node_seq_chunk;
        if (store_seq_at_cache)
        {
            int dad_depth = node->sequence->depth;
            dad_seq_chunk = &sequence_cache[dad_depth];
            node_seq_chunk = &sequence_cache[dad_depth + 1];
        }
        else
        {
            dad_seq_chunk = &node->sequence->sequence_chunks[thread_id];
            node_seq_chunk = &(*it)->node->sequence->sequence_chunks[thread_id];
        }
        
        // select the appropriate simulation method
        SIMULATION_METHOD simulation_method = RATE_MATRIX;
        if (((*it)->length * params->alisim_branch_scale > params->alisim_simulation_thresh && !(model->isMixture() && params->alisim_mixture_at_sub_level))
            || tree->getRate()->isHeterotachy()
            || (*it)->attributes.find("model") != (*it)->attributes.end())
            simulation_method = TRANS_PROB_MATRIX;
        
        // if branch_length is too short (less than 1 substitution is expected to occur) -> clone the sequence from the parent's node
        if ((*it)->length == 0)
        {
            // clone the sequence chunk from the parent node
            (*node_seq_chunk) = (*dad_seq_chunk);
        }
        else
        {
            // only simulate new sequence if the simulation type is Transition Probability Matrix approach
            if (simulation_method == TRANS_PROB_MATRIX)
            {
                // resize the sequence chunk
                if (!store_seq_at_cache)
                    (*node_seq_chunk).resize(segment_length);
                
                // if a model is specify for the current branch -> simulate the sequence based on that branch-specific model
                if ((*it)->attributes.find("model") != (*it)->attributes.end())
                {
                    branchSpecificEvolution(thread_id, sequence_length, *dad_seq_chunk, *node_seq_chunk, store_seq_at_cache, trans_matrix, node, it);
                }
                // otherwise, simulate the sequence based on the common model
                else
                {
                    // simulate the sequence chunk
                    simulateASequenceFromBranchAfterInitVariables(segment_start, model, trans_matrix, *dad_seq_chunk, *node_seq_chunk , node, it, rstream);
                }
                
                // handle indels
                if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
                    simulateSeqByGillespie(segment_start, segment_length, model, *node_seq_chunk, sequence_length, it, simulation_method, rstream);
            }
            // otherwise (Rate_matrix is used as the simulation method) + also handle Indels (if any).
            else
            {
                // clone the sequence chunk from the parent node
                (*node_seq_chunk) = (*dad_seq_chunk);
                
                // Each thread simulate a chunk of sequence using the Gillespie algorithm
                simulateSeqByGillespie(segment_start, segment_length, model, *node_seq_chunk, sequence_length, it, simulation_method, rstream);
            }
        }
        
        // set insertion position for of this node in the list of insertions if using Indels
        if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0 && (*it)->node->isLeaf())
        {
            (*it)->node->sequence->insertion_pos = latest_insertion;
            latest_insertion->phylo_nodes.push_back((*it)->node);
        }
        
        // handle dna error model
        if (model->containDNAerror())
        {
            // only handle DNA error on leaf
            if ((*it)->node->isLeaf())
            {
                // handle all model components one by one
                if (model->isMixture())
                {
                    for (int i = 0; i < model->getNMixtures(); i++)
                    handleDNAerr(segment_start, model->getDNAErrProb(i), *node_seq_chunk, rstream, i);
                }
                // otherwise, handle the DNA model
                else
                    handleDNAerr(segment_start, model->getDNAErrProb(), *node_seq_chunk, rstream);
            }
        }
        
        // write and delete the current chunk of sequence if possible
        writeAndDeleteSequenceChunkIfPossible(thread_id, segment_start, segment_length, *dad_seq_chunk, *node_seq_chunk, store_seq_at_cache, out, state_mapping, input_msa, it, node);
        
        // merge and write sequence in simulations with Indels or FunDi model
        mergeAndWriteSeqIndelFunDi(thread_id, out, sequence_length, state_mapping, input_msa, it, node);
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(thread_id, segment_start, segment_length, sequence_length, model, trans_matrix, sequence_cache, store_seq_at_cache, (*it)->node, node, out, state_mapping, input_msa, rstream);
    }
}

/**
    branch-specific evolution by multi threads
*/
void AliSimulator::branchSpecificEvolution(int thread_id, int sequence_length, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, bool store_seq_at_cache, double *trans_matrix, Node *node, NeighborVec::iterator it)
{
    unsigned short int num_threads_reach_barrier = 0;
    
    // clone sequence chunks from cache
    if (store_seq_at_cache)
    {
        #ifdef _OPENMP
        #pragma omp single nowait
        #endif
        {
            node->sequence->sequence_chunks.resize(num_simulating_threads);
            (*it)->node->sequence->sequence_chunks.resize(num_simulating_threads);
        }
        
        // manual implementation of barrier
        waitAtBarrier(1, (*it)->node);
        
        // clone sequence chunk from cache
        node->sequence->sequence_chunks[thread_id] = dad_seq_chunk;
        (*it)->node->sequence->sequence_chunks[thread_id].resize(dad_seq_chunk.size());
    }
    
    // manual implementation of barrier
    waitAtBarrier(2, (*it)->node);
            
    // only the first thread simulate the sequence
    if (thread_id == 0)
        branchSpecificEvolutionMasterThread(sequence_length, trans_matrix, node, it);
    
    // manual implementation of barrier
    waitAtBarrier(3, (*it)->node);
    
    // cache sequence chunks from node
    if (store_seq_at_cache)
    {
        // cache sequence chunk from node
        node_seq_chunk = (*it)->node->sequence->sequence_chunks[thread_id];
        
        // release memory allocated to node and dad node
        vector<short int>().swap((*it)->node->sequence->sequence_chunks[thread_id]);
        vector<short int>().swap(node->sequence->sequence_chunks[thread_id]);
        
        // manual implementation of barrier
        waitAtBarrier(4, (*it)->node);
        
        // wait to release memory allocated to the dad node
        #ifdef _OPENMP
        #pragma omp single nowait
        #endif
        {
            // release memory allocated to the node
            vector<vector<short int>>().swap((*it)->node->sequence->sequence_chunks);
            vector<vector<short int>>().swap(node->sequence->sequence_chunks);
        }
    }
}

void AliSimulator::waitAtBarrier(const unsigned short int barrier_count, Node* node)
{
    // manual implementation of barrier
    #ifdef _OPENMP
    #pragma omp atomic update
    #endif
    node->sequence->num_threads_reach_barrier++;
    // wait for all threads to reach this step
    unsigned short int num_threads_reach_barrier = node->sequence->num_threads_reach_barrier;
    while (num_threads_reach_barrier < barrier_count * num_simulating_threads)
    {
        #ifdef _OPENMP
        #pragma omp atomic read
        #endif
        num_threads_reach_barrier = node->sequence->num_threads_reach_barrier;
    }
}

/**
    temporarily write internal states to file (when using Indels)
*/
void AliSimulator::writeInternalStatesIndels(Node* node, ostream &out)
{
    out << node->name<<"@"<<node->sequence->sequence_chunks[0].size()<<"@";
    for (int i = 0; i < node->sequence->sequence_chunks[0].size(); i++)
        out << node->sequence->sequence_chunks[0][i]<<" ";
    out<<endl;
    
    // release the memory
    vector<short int>().swap(node->sequence->sequence_chunks[0]);
    
    map_seqname_node[node->name] = node;
}

/**
    merge and write sequence in simulations with Indels or FunDi model
*/
void AliSimulator::mergeAndWriteSeqIndelFunDi(int thread_id, ostream &out, int sequence_length, vector<string> &state_mapping, map<string,string> input_msa, NeighborVec::iterator it, Node* node)
{
    // only handle simulations with Indel or Fundi model in this function
    if (params->alisim_fundi_taxon_set.size() > 0 || params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
    {
        // increase num_threads_done_simulation
        bool this_thread_write_output = false;
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            (*it)->node->sequence->num_threads_done_simulation++;
            
            // make sure only one thread is selected to write the output
            if ((*it)->node->sequence->num_threads_done_simulation == omp_get_num_threads())
                this_thread_write_output = true;
        }
        
        // only one thread is selected to write the output
        if (this_thread_write_output)
        {
            // permuting selected sites for FunDi model. Notes: Delay permuting selected sites if Insertion (in Indels) is used
            if (params->alisim_fundi_taxon_set.size()>0 && params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0)
            {
                if (node->isLeaf())
                    permuteSelectedSites(fundi_items, node);
                if ((*it)->node->isLeaf())
                    permuteSelectedSites(fundi_items, (*it)->node);
            }
            
            // write sequence of leaf nodes to file if possible
            if (state_mapping.size() > 0)
            {
                if ((*it)->node->isLeaf())
                {
                    // if using Indels -> temporarily write out internal states
                    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
                        writeInternalStatesIndels((*it)->node, out);
                    else
                    {
                        // export pre_output string (containing taxon name and ">" or "space" based on the output format)
                        string pre_output = exportPreOutputString((*it)->node, params->aln_output_format, max_length_taxa_name, force_output_PHYLIP);
                        string output(num_sites_per_state == 1 ? sequence_length : (sequence_length * num_sites_per_state), '-');
                        
                        // convert numerical states into readable characters
                        string input_sequence = input_msa[(*it)->node->name];
                        if (input_sequence.length()>0)
                            // extract sequence and copying gaps from the input sequences to the output.
                            exportSequenceWithGaps((*it)->node->sequence->sequence_chunks[0], output, sequence_length, num_sites_per_state, input_sequence, state_mapping);
                        else
                            // extract sequence without copying gaps from the input sequences to the output.
                            convertNumericalStatesIntoReadableCharacters((*it)->node->sequence->sequence_chunks[0], output, sequence_length, num_sites_per_state, state_mapping);
                        
                        // release memory allocated to the sequence chunk
                        vector<short int>().swap((*it)->node->sequence->sequence_chunks[0]);
                        
                        // write output to file
                        #ifdef _OPENMP
                        #pragma omp critical
                        #endif
                            out << pre_output << output << "\n";
                    }
                }
                
                // avoid writing sequence of __root__
                if (node->isLeaf() && (node->name!=ROOT_NAME || params->alisim_write_internal_sequences))
                {
                    // if using Indels -> temporarily write out internal states
                    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
                        writeInternalStatesIndels(node, out);
                    else
                    {
                        // export pre_output string (containing taxon name and ">" or "space" based on the output format)
                        string pre_output = exportPreOutputString(node, params->aln_output_format, max_length_taxa_name, force_output_PHYLIP);
                        string output(num_sites_per_state == 1 ? sequence_length : (sequence_length * num_sites_per_state), '-');
                        
                        // convert numerical states into readable characters
                        string input_sequence = input_msa[node->name];
                        if (input_sequence.length()>0)
                            // extract sequence and copying gaps from the input sequences to the output.
                            exportSequenceWithGaps(node->sequence->sequence_chunks[0], output, sequence_length, num_sites_per_state, input_sequence, state_mapping);
                        else
                            // extract sequence without copying gaps from the input sequences to the output.
                            convertNumericalStatesIntoReadableCharacters(node->sequence->sequence_chunks[0], output, sequence_length, num_sites_per_state, state_mapping);
                        
                        // release memory allocated to the sequence chunk
                        vector<short int>().swap(node->sequence->sequence_chunks[0]);
                        
                        // write output to file
                        #ifdef _OPENMP
                        #pragma omp critical
                        #endif
                            out << pre_output << output << "\n";
                    }
                }
            }
        }
    }
}

/**
    write and delete the current chunk of sequence if possible
*/
void AliSimulator::writeAndDeleteSequenceChunkIfPossible(int thread_id, int segment_start, int segment_length, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, bool store_seq_at_cache, ostream &out, vector<string> &state_mapping, map<string,string> input_msa, NeighborVec::iterator it, Node* node)
{
    // we can only write and delete sequence chunk in normal simulations: without Indels, FunDi, ASC, etc
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio == 0 && state_mapping.size() > 0)
    {
        if (params->alisim_fundi_taxon_set.size() == 0)
        {
            // convert the sequence at leaf
            if ((*it)->node->isLeaf())
            {
                // init a default sequence str
                int sequence_length = round(expected_num_sites * inverse_length_ratio);
                string output(num_sites_per_state == 1 ? segment_length : (segment_length * num_sites_per_state), '-');
                
                // convert numerical states into readable characters
                string input_sequence = input_msa[(*it)->node->name];
                if (input_sequence.length()>0)
                    // extract sequence and copying gaps from the input sequences to the output.
                    exportSequenceWithGaps(node_seq_chunk, output, sequence_length, num_sites_per_state, input_sequence, state_mapping, segment_start, segment_length);
                else
                    // extract sequence without copying gaps from the input sequences to the output.
                    convertNumericalStatesIntoReadableCharacters(node_seq_chunk, output, sequence_length, num_sites_per_state, state_mapping, segment_length);
                
                // output a sequence (chunk) to file (if using AliSim-OpenMP-EM) or store it to common cache (if using AliSim-OpenMP-IM)
                outputOneSequence((*it)->node, output, thread_id, segment_start, out);
            }
                
            // avoid writing sequence of __root__
            if (node->isLeaf() && (node->name!=ROOT_NAME || params->alisim_write_internal_sequences))
            {
                // init a default sequence str
                int sequence_length = round(expected_num_sites * inverse_length_ratio);
                string output(num_sites_per_state == 1 ? segment_length : (segment_length * num_sites_per_state), '-');
                
                // convert numerical states into readable characters
                string input_sequence = input_msa[node->name];
                if (input_sequence.length()>0)
                    // extract sequence and copying gaps from the input sequences to the output.
                    exportSequenceWithGaps(dad_seq_chunk, output, sequence_length, num_sites_per_state, input_sequence, state_mapping, segment_start, segment_length);
                else
                    // extract sequence without copying gaps from the input sequences to the output.
                    convertNumericalStatesIntoReadableCharacters(dad_seq_chunk, output, sequence_length, num_sites_per_state, state_mapping, segment_length);
                
                // output a sequence (chunk) to file (if using AliSim-OpenMP-EM) or store it to common cache (if using AliSim-OpenMP-IM)
                outputOneSequence(node, output, thread_id, segment_start, out);
            }
        }
        
        // convert internal sequence if it's an internal node and the user want to output internal sequences
        if (!(*it)->node->isLeaf() && params->alisim_write_internal_sequences)
        {
            // init a default sequence str
            int sequence_length = round(expected_num_sites * inverse_length_ratio);
            string output(num_sites_per_state == 1 ? segment_length : (segment_length * num_sites_per_state), '-');
            
            // convert numerical states into readable characters
            convertNumericalStatesIntoReadableCharacters(node_seq_chunk, output, sequence_length, num_sites_per_state, state_mapping, segment_length);
            
            // the memory allocated to the current sequence chunk of INTERNAL nodes will be release later
            
            // output a sequence (chunk) to file (if using AliSim-OpenMP-EM) or store it to common cache (if using AliSim-OpenMP-IM)
            outputOneSequence((*it)->node, output, thread_id, segment_start, out);
        }
    }
    
    if (!store_seq_at_cache)
    {
        // update the num_children_done_simulation
        node->sequence->nums_children_done_simulation[thread_id]++;
        
        // delete sequence at internal node if all sequences of its children are simulated
        if (!node->isLeaf() && node->sequence->nums_children_done_simulation[thread_id] >= (node->neighbors.size() - 1) && !(params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0 && params->alisim_write_internal_sequences))
            vector<short int>().swap(dad_seq_chunk);
    }
}

void AliSimulator::outputOneSequence(Node* node, string &output, int thread_id, int segment_start, ostream &out)
{
    // output a sequence with AliSim-OpenMP-EM
    if (params->alisim_openmp_alg == EM)
    {
        // only write sequence name in the first thread
        if (thread_id == 0)
        {
            string pre_output = exportPreOutputString(node, params->aln_output_format, max_length_taxa_name, force_output_PHYLIP);
            out << pre_output << output << "\n";
        }
        // write sequence chunk in other threads
        else
        {
            // if users skip concatenating sequence chunks from intermediate files -> write sequence name
            if (params->no_merge)
            {
                string pre_output = exportPreOutputString(node, params->aln_output_format, max_length_taxa_name);
                out << pre_output << output << "\n";
            }
            // don't write the sequence name in intermediate files
            else
                out << output << "\n";
        }
    }
    // output a sequence with AliSim-OpenMP-IM
    else
    {
        // only write sequence name in the first thread
        if (thread_id == 0)
            output = exportPreOutputString(node, params->aln_output_format, max_length_taxa_name) + output;
        
        // add break-line in the last simulating thread
        if (thread_id == num_simulating_threads - 1)
            output = output + "\n";
        
        //  cache output into the writing queue
        if (num_threads > 1)
        {
            int64_t pos = ((int64_t)node->id) * ((int64_t)output_line_length);
            pos += starting_pos + (num_sites_per_state == 1 ? segment_start : (segment_start * num_sites_per_state)) + (thread_id == 0 ? 0 : seq_name_length);
            cacheSeqChunkStr(pos, output, thread_id);
        }
        // write output to file
        else
            out << output;
    }
}

void AliSimulator::cacheSeqChunkStr(int64_t pos, string seq_chunk_str, int thread_id)
{
    // seek an empty slot in the cache
    int slot_id = -1;
    int start_slot = cache_start_indexes[thread_id];
    int end_slot = start_slot + cache_size_per_thread;
    int i = start_slot;
    while (slot_id == -1)
    {
        // try to find an empty slot
        if (seq_str_cache[i].chunk_status == EMPTY)
            slot_id = i;
        // otherwise, move to next slot
        else
            i >= end_slot - 1 ? i = start_slot : ++i;
    }
    
    #ifdef _OPENMP
    #pragma omp flush
    #endif
    // store the current chunk to the selected slot
    seq_str_cache[slot_id].chunk_str = seq_chunk_str;
    seq_str_cache[slot_id].pos = pos;
    #ifdef _OPENMP
    #pragma omp flush
    #endif
    
    // update the status of the selected slot
    #ifdef _OPENMP
    #pragma omp atomic write
    #endif
    seq_str_cache[slot_id].chunk_status = OCCUPIED;
}

/**
*  export pre_output string (contains taxon name and ">" or "space" based on the output format
*
*/
string AliSimulator::exportPreOutputString(Node *node, InputType output_format, int max_length_taxa_name, bool force_PHYLIP)
{
    string pre_output = "";
    
    // add node's name
    pre_output = node->name;
    // write node's id if node's name is empty
    if (pre_output.length() == 0) pre_output = convertIntToString(node->id);
    pre_output.resize(max_length_taxa_name, ' ');
    
    // in FASTA format
    if (output_format == IN_FASTA && !force_PHYLIP)
    {
        pre_output = ">" + pre_output;
        pre_output[pre_output.length() - 1] = '\n';
    }
    
    return pre_output;
}

/**
*  get a random item from a set of items with a probability array
*/
int AliSimulator::getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items, int* rstream)
{
    // generate a random number
    double random_number = random_double(rstream);
    
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
void AliSimulator::convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns, bool force_round_1)
{
    double* probability_maxtrix_pointer = probability_maxtrix;
    for (int r = 0; r < num_rows; r++, probability_maxtrix_pointer += num_columns)
    {
        for (int c = 1; c < num_columns; c++)
            probability_maxtrix_pointer[c] += probability_maxtrix_pointer[c - 1];
    }
    
    // force rounding the last entry of each row to one to avoid potential bug due to numerical precision, e.g, 0.99999999 < 1
    if (force_round_1)
    {
        probability_maxtrix_pointer = probability_maxtrix;
        int last_entry_index = num_columns - 1;
        for (int r = 0; r < num_rows; r++, probability_maxtrix_pointer += num_columns)
            probability_maxtrix_pointer[last_entry_index] = 1;
    }
}

/**
*  get a random item from a set of items with an accumulated probability array by binary search starting at the max probability
*/
int AliSimulator::getRandomItemWithAccumulatedProbMatrixMaxProbFirst(double *accumulated_probability_maxtrix, int starting_index, int num_columns, int max_prob_position, int* rstream){
    // generate a random number
    double random_number = random_double(rstream);
    
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
    int center = (start + end) * 0.5;
    
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
int AliSimulator::binarysearchItemWithAccumulatedProbabilityMatrix(vector<double> &accumulated_probability_maxtrix, double random_number, int start, int end, int first)
{
    // check search range
    if (start > end)
        return -1; // return -1 ~ not found
    
    // compute the center index
    int center = (start + end) * 0.5;
    
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
    {
        if (params->aln_file || params->alisim_ancestral_sequence_aln_filepath || params->original_params.find("--length") != std::string::npos)
            outError("Sequence length of Codon must be divisible by 3. Please check & try again!");
        else
            params->alisim_sequence_length = 999;
    }
}

/**
*  update the expected_num_sites due to the change of the sequence_length
*/
void AliSimulator::refreshExpectedNumSites(){
    expected_num_sites = round(params->alisim_sequence_length / num_sites_per_state * length_ratio);
}

/**
    estimate length_ratio (for models with +ASC)
*/
void AliSimulator::estimateLengthRatio()
{
    // By default (without +ASC), length_ratio is set at 1
    length_ratio = 1;
    inverse_length_ratio = 1;
        
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
            all_characters.resize(max_num_states * num_sites_per_state);
            int index = 0;
            for (int i = 0; i < max_num_states; i++, index += num_sites_per_state)
            {
                string characters_from_state = tree->aln->convertStateBackStr(i);
                // NHANLT: potential improvement
                // change to use pointer to avoid accessing []
                for (int j = 0; j < num_sites_per_state; j++)
                    all_characters[index + j] = characters_from_state[j];
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
            tree->aln->buildPattern(sequences, sequence_type, nseq, nsite * num_sites_per_state);
            
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
            length_ratio = 1.0 / (1 - estimated_length_ratio) + 0.1;
        }
        
        // compute inverse_length_ratio
        inverse_length_ratio = 1.0 / length_ratio;
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
void AliSimulator::convertNumericalStatesIntoReadableCharacters(vector<short int> &sequence_chunk, string &output, int sequence_length, int num_sites_per_state, vector<string> &state_mapping, int segment_length)
{
    segment_length = segment_length == -1 ? sequence_length : segment_length;
    ASSERT(segment_length <= sequence_chunk.size());
    
    // convert normal data
    if (num_sites_per_state == 1)
    {
        for (int i = 0; i < segment_length; i++)
            output[i] = state_mapping[sequence_chunk[i]][0];
    }
    // convert CODON
    else
    {
        int index = 0;
        for (int i = 0; i < segment_length; i++, index += num_sites_per_state)
        {
            string codon_str = state_mapping[sequence_chunk[i]];
            output[index] = codon_str[0];
            output[index + 1] = codon_str[1];
            output[index + 2] = codon_str[2];
        }
    }
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
void AliSimulator::extractMaxTaxaNameLength()
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
                    if (seq_names[i].length() >= max_length_taxa_name)
                        max_length_taxa_name = seq_names[i].length() + 1;
            }
        }
        // otherwise, just check the current tree
        else
        {
            vector<string> seq_names = tree->aln->getSeqNames();
            for (int i = 0; i < seq_names.size(); i++)
                if (seq_names[i].length() >= max_length_taxa_name)
                    max_length_taxa_name = seq_names[i].length() + 1;
        }
    }
}

/**
    selecting & permuting sites (FunDi models)
*/
void AliSimulator::selectAndPermuteSites(vector<FunDi_Item> &fundi_items, double proportion, int num_sites){
    ASSERT(proportion<1);
    
    // dummy variables
    IntVector tmp_selected_sites;
    int num_selected_sites = round(proportion * num_sites);
    
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
}

/**
    permuting selected sites (FunDi models)
*/
void AliSimulator::permuteSelectedSites(vector<FunDi_Item> &fundi_items, Node* node)
{
    mergeChunks(node);
    
    if (std::find(params->alisim_fundi_taxon_set.begin(), params->alisim_fundi_taxon_set.end(), node->name) != params->alisim_fundi_taxon_set.end())
    {
        // permute selected sites in the sequence
        // caching the current states of all selected sites
        map<int, short int> caching_sites;
        for (int i = 0; i < fundi_items.size(); i++)
            caching_sites[fundi_items[i].selected_site] = node->sequence->sequence_chunks[0][fundi_items[i].selected_site];
        
        // permuting sites in FunDi model
        for (int i = 0; i < fundi_items.size(); i++)
            node->sequence->sequence_chunks[0][fundi_items[i].new_position] = caching_sites[fundi_items[i].selected_site];
    }
}

/**
    process delayed Fundi if it is delayed due to Insertion events
*/
void AliSimulator::processDelayedFundi(Node *node, Node *dad)
{
    // permute selected sites of the current node
    if (node->isLeaf())
        permuteSelectedSites(fundi_items, node);
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        processDelayedFundi((*it)->node, node);
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
    branch-specific evolution by the master thread
*/
void AliSimulator::branchSpecificEvolutionMasterThread(int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it)
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
    tmp_alisimulator->num_threads = num_threads;
    
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
    initVariablesRateHeterogeneity(sequence_length);
    
    // check to regenerate the root sequence if the user has specified specific frequencies for root
    if (tree->root->id == node->id && ((*it)->attributes["freqs"]).length() > 0)
        regenerateRootSequenceBranchSpecificModel((*it)->attributes["freqs"], sequence_length, node);
    
    // simulate chunks of sequence one by one
    int segment_start = 0;
    for (int i = 0; i < node->sequence->sequence_chunks.size(); i++)
    {
        simulateASequenceFromBranchAfterInitVariables(segment_start, model, trans_matrix, node->sequence->sequence_chunks[i], (*it)->node->sequence->sequence_chunks[i], node, it, NULL, lengths);
        segment_start += node->sequence->sequence_chunks[i].size();
    }
}

/**
    simulate a sequence for a node from a specific branch after all variables has been initializing
*/
void AliSimulator::simulateASequenceFromBranchAfterInitVariables(int segment_start, ModelSubst *model, double *trans_matrix, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, Node *node, NeighborVec::iterator it, int* rstream, string lengths)
{
    // compute the transition probability matrix
    model->computeTransMatrix(partition_rate * params->alisim_branch_scale * (*it)->length, trans_matrix);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(trans_matrix, max_num_states, max_num_states);
    
    // estimate the sequence for the current neighbor
    for (int i = 0; i < node_seq_chunk.size(); i++)
    {
        // if the parent's state is a gap -> the children's state should also be a gap
        if (dad_seq_chunk[i] == STATE_UNKNOWN)
            node_seq_chunk[i] = STATE_UNKNOWN;
        else
        {
            // iteratively select the state for each site of the child node, considering it's dad states, and the transition_probability_matrix
            int parent_state = dad_seq_chunk[i];
            node_seq_chunk[i] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(trans_matrix, parent_state * max_num_states, max_num_states, parent_state, rstream);
        }
    }
}

/**
    initialize variables (e.g., site-specific rate)
*/
void AliSimulator::initVariablesRateHeterogeneity(int sequence_length, bool regenerate_root_sequence)
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
    if (fabs(total_freq-1.0) >= 1e-7)
    {
        outWarning("Normalizing state frequencies so that sum of them equals to 1.");
        normalize_frequencies(state_freqs, max_num_states, total_freq);
    }
    
    // re-generate a new sequence for the root from the state frequencies
    generateRandomSequenceFromStateFreqs(sequence_length, root->sequence->sequence_chunks[0], state_freqs, max_prob_pos);
    
    // release the memory of state_freqs
    delete[] state_freqs;
}

/**
    generate a random sequence by state frequencies
*/
void AliSimulator::generateRandomSequenceFromStateFreqs(int sequence_length, vector<short int> &sequence, double* state_freqs, int max_prob_pos)
{
    sequence.resize(sequence_length);
    
    // convert the probability matrix into an accumulated probability matrix
    convertProMatrixIntoAccumulatedProMatrix(state_freqs, 1, max_num_states);
    
    // randomly generate each site in the sequence follows the base frequencies defined by the user
    for (int i = 0; i < sequence_length; i++)
        sequence[i] =  getRandomItemWithAccumulatedProbMatrixMaxProbFirst(state_freqs, 0, max_num_states, max_prob_pos, NULL);
}

/**
*  export a sequence with gaps copied from the input sequence
*/
void AliSimulator::exportSequenceWithGaps(vector<short int> &sequence_chunk, string &output, int sequence_length, int num_sites_per_state, string input_sequence, vector<string> &state_mapping, int segment_start, int segment_length)
{
    segment_length = segment_length == -1 ? sequence_length : segment_length;
    
    // convert non-empty sequence
    if (sequence_chunk.size() >= segment_length)
    {
        // convert normal data
        if (num_sites_per_state == 1)
        {
            for (int i = 0; i < segment_length; i++){
                // handle gaps
                if (segment_start + (i + 1) - 1 < input_sequence.length()
                    && input_sequence[segment_start + i] == '-')
                {
                    // insert gaps
                    output[i] = '-';
                }
                // if it's not a gap
                else
                    output[i] = state_mapping[sequence_chunk[i]][0];
            }
        }
        // convert CODON
        else
        {
            int index = 0;
            int segment_start_plus_index = segment_start;
            for (int i = 0; i < segment_length; i++, index += num_sites_per_state, segment_start_plus_index += num_sites_per_state)
            {
                // handle gaps
                if (segment_start_plus_index + 2 < input_sequence.length()
                    && (input_sequence[segment_start_plus_index] == '-'
                            || input_sequence[segment_start_plus_index + 1] == '-'
                            || input_sequence[segment_start_plus_index + 2] == '-')){
                    // insert gaps
                    output[index] =  input_sequence[segment_start_plus_index];
                    output[index + 1] =  input_sequence[segment_start_plus_index + 1];
                    output[index + 2] =  input_sequence[segment_start_plus_index + 2];
                }
                else
                {
                    string codon_str = state_mapping[sequence_chunk[i]];
                    output[index] = codon_str[0];
                    output[index + 1] = codon_str[1];
                    output[index + 2] = codon_str[2];
                }
            }
        }
    }
}

/**
    extract array of substitution rates
*/
void AliSimulator::extractRatesJMatrix(ModelSubst *model)
{
    // get num_mixture_models
    int num_mixture_models = model->getNMixtures();
    int max_num_states_square = max_num_states * max_num_states;
    double* tmp_Q_matrix = new double[max_num_states_square];
    int starting_index_sub_rates = 0;
    int starting_index_J = 0;
    int max_num_states_plus_one = max_num_states + 1;
    
    for (int mixture = 0; mixture <num_mixture_models; mixture++, starting_index_sub_rates += max_num_states, starting_index_J += max_num_states_square)
    {
        // get the Rate (Q) Matrix
        model->getQMatrix(tmp_Q_matrix, mixture);
        
        // extract sub rate for each state from q_matrix
        int tmp_Q_matrix_index = 0;
        for (int i = 0; i < max_num_states; i++, tmp_Q_matrix_index += max_num_states_plus_one)
            sub_rates[starting_index_sub_rates + i] = - tmp_Q_matrix[tmp_Q_matrix_index];
        
        // convert Q_Matrix to J_Matrix;
        int i_times_max_num_states = 0;
        for (int i = 0; i < max_num_states; i++, i_times_max_num_states += max_num_states)
        {
            double denominator = 1.0 / sub_rates[starting_index_sub_rates + i];
            int tmp_index = starting_index_J + i_times_max_num_states;
            for (int j = 0; j < max_num_states; j++)
            {
                if (i == j)
                    Jmatrix[tmp_index + j] = 0;
                else
                    Jmatrix[tmp_index + j] = tmp_Q_matrix[i_times_max_num_states + j] * denominator;
            }
        }
    }
    
    // delete tmp_Q_matrix
    delete[] tmp_Q_matrix;
    
    // convert Jmatrix to accumulated Jmatrix
    convertProMatrixIntoAccumulatedProMatrix(Jmatrix, num_mixture_models * max_num_states, max_num_states);
}

/**
    initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
*/
void AliSimulator::initVariables4RateMatrix(int segment_start, double &total_sub_rate, int &num_gaps, vector<double> &sub_rate_by_site, vector<short int> &sequence)
{
    // initialize variables
    total_sub_rate = 0;
    num_gaps = 0;
    sub_rate_by_site.resize(sequence.size(), 0);
    
    vector<int> sub_rate_count(max_num_states, 0);
    // compute sub_rate_by_site
    for (int i = 0; i < sequence.size(); i++)
    {
        // not compute the substitution rate for gaps/deleted sites
        if (sequence[i] != STATE_UNKNOWN && (site_specific_rates.size() == 0 || site_specific_rates[segment_start + i] != 0))
        {
            int index = sequence[i];
            sub_rate_count[index]++;
            sub_rate_by_site[i] = sub_rates[index];
        }
        else
        {
            sub_rate_by_site[i] = 0;
            
            if (sequence[i] == STATE_UNKNOWN)
                num_gaps++;
        }
    }
    
    // update total_sub_rate
    for (int i = 0; i < max_num_states; i++)
        total_sub_rate += sub_rate_count[i] * sub_rates[i];
}

/**
    handle indels
*/
void AliSimulator::simulateSeqByGillespie(int segment_start, int &segment_length, ModelSubst *model, vector<short int> &node_seq_chunk, int &sequence_length, NeighborVec::iterator it, SIMULATION_METHOD simulation_method, int *rstream)
{
    int num_gaps = 0;
    double total_sub_rate = 0;
    vector<double> sub_rate_by_site;
    // If AliSim is using RATE_MATRIX approach -> initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
    if (simulation_method == RATE_MATRIX || params->indel_rate_variation)
    {
        initVariables4RateMatrix(segment_start, total_sub_rate, num_gaps, sub_rate_by_site, node_seq_chunk);
        
        // handle cases when total_sub_rate == NaN due to extreme freqs
        if (total_sub_rate != total_sub_rate)
            total_sub_rate = 0;
    }
    else // otherwise, TRANS_PROB_MATRIX approach is used -> only count the number of gaps
        num_gaps = (*it)->node->sequence->num_gaps;
    
    double total_ins_rate = 0;
    double total_del_rate = 0;
    if (params->alisim_insertion_ratio + params->alisim_deletion_ratio > 0)
    {
        // compute total_ins_rate and total_del_rate
        // constant indel-rates
        if (!params->indel_rate_variation)
        {
            int tmp_num = sequence_length - num_gaps;
            total_ins_rate = params->alisim_insertion_ratio * (tmp_num + 1);
            total_del_rate = params->alisim_deletion_ratio * (tmp_num - 1 + computeMeanDelSize(sequence_length));
        }
        // indel-rate variation
        else
        {
            total_ins_rate = params->alisim_insertion_ratio * total_sub_rate;
            total_del_rate = params->alisim_deletion_ratio * total_sub_rate;
        }
    }
    
    // compute total_event_rate
    double total_event_rate = total_ins_rate + total_del_rate;
    // only take into account the substitution rate if using Rate matrix simulation approach
    if (simulation_method == RATE_MATRIX)
        total_event_rate += total_sub_rate;
    
    // dummy variables
    int ori_seq_length = node_seq_chunk.size();
    Insertion* insertion_before_simulation = latest_insertion;
    
    double branch_length = (*it)->length * params->alisim_branch_scale;
    while (branch_length > 0)
    {
        // generate a waiting time s1 by sampling from the exponential distribution with mean 1/total_event_rate
        double waiting_time = random_double_exponential_distribution(1.0 / total_event_rate, rstream);
        
        if (waiting_time > branch_length)
            break;
        else
        {
            // update the remaining length of the current branch
            branch_length -=  waiting_time;
            
            // Determine the event type (insertion, deletion, substitution) occurs
            EVENT_TYPE event_type = SUBSTITUTION;
            if (total_ins_rate > 0 || total_del_rate > 0)
            {
                double random_num = random_double()*total_event_rate;
                if (random_num < total_ins_rate)
                    event_type = INSERTION;
                else if (random_num < total_ins_rate+total_del_rate)
                    event_type = DELETION;
            }
            
            // process event
            int length_change = 0;
            switch (event_type)
            {
                case INSERTION:
                {
                    length_change = handleInsertion(sequence_length, node_seq_chunk, total_sub_rate, sub_rate_by_site, simulation_method);
                    segment_length = sequence_length;
                    break;
                }
                case DELETION:
                {
                    int deletion_length = handleDeletion(sequence_length, node_seq_chunk, total_sub_rate, sub_rate_by_site, simulation_method);
                    length_change = -deletion_length;
                    (*it)->node->sequence->num_gaps += deletion_length;
                    break;
                }
                case SUBSTITUTION:
                {
                    if (simulation_method == RATE_MATRIX)
                    {
                        handleSubs(segment_start, total_sub_rate, sub_rate_by_site, node_seq_chunk, model->getNMixtures(), rstream);
                    }
                    break;
                }
                default:
                    break;
            }
            
            // update total_event_rate
            if (length_change != 0)
            {
                // constant indel-rate
                if (!params->indel_rate_variation)
                {
                    total_ins_rate += params->alisim_insertion_ratio * length_change;
                    total_del_rate += params->alisim_deletion_ratio * length_change;
                }
                // indel-rate variation
                else
                {
                    total_ins_rate = params->alisim_insertion_ratio * total_sub_rate;
                    total_del_rate = params->alisim_deletion_ratio * total_sub_rate;
                }
            }
            
            // re-compute total_event_rate
            total_event_rate = total_ins_rate + total_del_rate;
            // only take into account the substitution rate if using Rate matrix simulation approach
            if (simulation_method == RATE_MATRIX)
                total_event_rate += total_sub_rate;
        }

    }
    
    // if insertion events occur -> insert gaps to other nodes
    if (insertion_before_simulation && insertion_before_simulation->next)
    {
        // init a genome_tree to update new genomes for internal nodes
        GenomeTree* genome_tree = new GenomeTree();
        genome_tree->buildGenomeTree(insertion_before_simulation, ori_seq_length);
        
        // update non-empty internal sequences due to insertions
        updateInternalSeqsIndels(genome_tree, sequence_length, (*it)->node);
        
        // delete genome_tree
        delete genome_tree;
        
        // re-compute the switching param to switch between Rate matrix and Probability matrix
        computeSwitchingParam(sequence_length);
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
*  update internal sequences due to Indels
*
*/
void AliSimulator::updateInternalSeqsIndels(GenomeTree* genome_tree, int seq_length, Node *node)
{
    // if we need to output all internal sequences -> traverse tree from root to the current node to update all internal sequences
    if (params->alisim_write_internal_sequences)
    {
        bool stop_inserting_gaps = false;
        updateInternalSeqsFromRootToNode(genome_tree, seq_length, node->id, tree->root, tree->root, stop_inserting_gaps);
    }
    // otherwise, only need to update sequences on the path from the current node to root
    else
        updateInternalSeqsFromNodeToRoot(genome_tree, seq_length, node);
}

/**
*  update all simulated internal seqs from root to the current node due to insertions
*
*/
void AliSimulator::updateInternalSeqsFromRootToNode(GenomeTree* genome_tree, int seq_length, int stopping_node_id, Node *node, Node* dad, bool &stop_inserting_gaps)
{
    // check to stop
    if (stop_inserting_gaps)
        return;
    
    // if it is a non-empty internal node -> update the current genome by the genome_tree
    if ((!node->isLeaf() || node->name == ROOT_NAME) && node->sequence->sequence_chunks[0].size() > 0)
    {
        node->sequence->num_gaps += seq_length - node->sequence->sequence_chunks[0].size();
        node->sequence->sequence_chunks[0] = genome_tree->exportNewGenome(node->sequence->sequence_chunks[0], seq_length, tree->aln->STATE_UNKNOWN);
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
        updateInternalSeqsFromRootToNode(genome_tree, seq_length, stopping_node_id, (*it)->node, node, stop_inserting_gaps);
    }
}

/**
*  update internal seqs on the path from the current phylonode to root due to insertions
*
*/
void AliSimulator::updateInternalSeqsFromNodeToRoot(GenomeTree* genome_tree, int seq_length, Node *node)
{
    // get parent node
    Node* internal_node = node->sequence->parent;
    
    for (; internal_node;)
    {
        // only update new genome at non-empty internal nodes
        if (!(internal_node->isLeaf()) && internal_node->sequence->sequence_chunks[0].size() > 0)
        {
            internal_node->sequence->num_gaps += seq_length - internal_node->sequence->sequence_chunks[0].size();
            internal_node->sequence->sequence_chunks[0] = genome_tree->exportNewGenome(internal_node->sequence->sequence_chunks[0], seq_length, tree->aln->STATE_UNKNOWN);
        }
        
        // move to the next parent
        internal_node = internal_node->sequence->parent;
    }
}

/**
    handle insertion events
*/
int AliSimulator::handleInsertion(int &sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method)
{
    // Randomly select the position/site (from the set of all sites) where the insertion event occurs
    int position;
    // with constant indel-rate -> based on a uniform distribution between 0 and the current length of the sequence
    if (!params->indel_rate_variation)
        position = selectValidPositionForIndels(sequence_length + 1, indel_sequence);
    // with indel-rate variation -> based on the sub_rate_by_site
    else
    {
        discrete_distribution<> random_discrete_dis(sub_rate_by_site.begin(), sub_rate_by_site.end());
        position = random_discrete_dis(params->generator);
    }
    
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
    vector<short int> new_sequence;
    generateRandomSequence(length, new_sequence, false);
    insertNewSequenceForInsertionEvent(indel_sequence, position, new_sequence);
    
    // if RATE_MATRIX approach is used -> update total_sub_rate and sub_rate_by_site
    if (simulation_method == RATE_MATRIX || params->indel_rate_variation)
    {
        // update sub_rate_by_site of the inserted sites
        double sub_rate_change = 0;
        sub_rate_by_site.insert(sub_rate_by_site.begin()+position, length, 0);
        for (int i = position; i < position + length; i++)
        {
            // NHANLT: potential improvement
            // cache site_specific_model_index[i] * max_num_states
            double sub_rate_from_model = site_specific_model_index.size() == 0 ? sub_rates[indel_sequence[i]] : sub_rates[site_specific_model_index[i] * max_num_states + indel_sequence[i]];
            sub_rate_by_site[i] = site_specific_rates.size() > 0 ? (site_specific_rates[i] * sub_rate_from_model) : sub_rate_from_model;
            sub_rate_change += sub_rate_by_site[i];
        }
        
        // update total_sub_rate
        total_sub_rate += sub_rate_change;
    }
    
    // record the insertion event
    Insertion* new_insertion = new Insertion(position, length, position == sequence_length);
    latest_insertion->next = new_insertion;
    latest_insertion = new_insertion;
    
    // update the sequence_length
    sequence_length += length;
    
    // return insertion-size
    return length;
}

/**
    handle deletion events
*/
int AliSimulator::handleDeletion(int sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method)
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
    
    // Randomly select the position/site (from the set of all sites) where the deletion event occurs
    int position = 0;
    // with constant indel-rate -> based on a uniform distribution between 0 and the current length of the sequence
    if (!params->indel_rate_variation)
    {
        int upper_bound = sequence_length - length;
        if (upper_bound > 0)
            position = selectValidPositionForIndels(upper_bound, indel_sequence);
    }
    // with indel-rate variation -> based on the sub_rate_by_site
    else
    {
        discrete_distribution<> random_discrete_dis(sub_rate_by_site.begin(), sub_rate_by_site.end());
        position = random_discrete_dis(params->generator);
    }
    
    // Replace up to length_D sites by gaps from the sequence starting at the selected location
    int real_deleted_length = 0;
    double sub_rate_change = 0;
    for (int i = 0; i < length && (position + i) < indel_sequence.size(); i++)
    {
        // if the current site is not a gap (has not been deleted) -> replacing it by a gap
        if (indel_sequence[position + i ] != STATE_UNKNOWN)
        {
            indel_sequence[position + i ] = STATE_UNKNOWN;
            real_deleted_length++;
        }
        // otherwise, ignore the current site, moving forward to find a valid site (not a gap)
        else
        {
            i--;
            position++;
        }
        
        // if RATE_MATRIX approach is used -> update sub_rate_by_site
        if (simulation_method == RATE_MATRIX || params->indel_rate_variation)
        {
            sub_rate_change -= sub_rate_by_site[position + i];
            sub_rate_by_site[position + i] = 0;
        }
    }
    
    // if RATE_MATRIX approach is used -> update total_sub_rate
    if (simulation_method == RATE_MATRIX || params->indel_rate_variation)
        total_sub_rate += sub_rate_change;
    
    // return deletion-size
    return real_deleted_length;
}

/**
    handle substitution events
*/
void AliSimulator::handleSubs(int segment_start, double &total_sub_rate, vector<double> &sub_rate_by_site, vector<short int> &indel_sequence, int num_mixture_models, int* rstream)
{
    // select a position where the substitution event occurs
    discrete_distribution<> random_discrete_dis(sub_rate_by_site.begin(), sub_rate_by_site.end());
    int pos = random_discrete_dis(params->generator);
    
    // extract the current state
    short int current_state = indel_sequence[pos];
    
    // estimate the new state
    int mixture_index = 0;
    // randomly select a model component if mixture model at substitution level is used
    if (site_specific_model_index.size() > segment_start + pos)
    {
        if (params->alisim_mixture_at_sub_level)
            mixture_index = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(mixture_accumulated_weight, 0, num_mixture_models, mixture_max_weight_pos, rstream);
        else
            mixture_index = site_specific_model_index[segment_start + pos];
    }
    
    int mixture_index_times_num_states = (mixture_index == 0 ? 0 : (mixture_index * max_num_states));
    int starting_index = (mixture_index_times_num_states + current_state) * max_num_states;
    indel_sequence[pos] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(Jmatrix, starting_index, max_num_states, max_num_states * 0.5, rstream);
    
    // update total_sub_rate
    double sub_rate_change = sub_rates[mixture_index_times_num_states + indel_sequence[pos]] - sub_rates[mixture_index_times_num_states + current_state];
    sub_rate_change = (site_specific_rates.size() == 0 ? sub_rate_change : (sub_rate_change * site_specific_rates[segment_start + pos]));
    total_sub_rate += sub_rate_change;
    
    // update sub_rate_by_site
    sub_rate_by_site[pos] += sub_rate_change;
}

/**
*  randomly select a valid position (not a deleted-site) for insertion/deletion event
*
*/
int AliSimulator::selectValidPositionForIndels(int upper_bound, vector<short int> &sequence)
{
    int position = -1;
    for (int i = 0; i < upper_bound; i++)
    {
        position = random_int(upper_bound);
        
        // try to move to the following site if the selected site is a gap
        if (position < sequence.size() && sequence[position] == STATE_UNKNOWN)
            for (; position < upper_bound; position++)
                if (position == sequence.size() || sequence[position] != STATE_UNKNOWN)
                    break;
        
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
            random_size = random_int_geometric(indel_dis.param_1);
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
            int random_size = generateIndelSize(params->alisim_deletion_distribution);
            
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
        tree->root->sequence = new Sequence();
        
        // change the id of the intermediate node if it's equal to the root's id
        if (second_internal_node->id == new_root->id)
            second_internal_node->id = (tree->nodeNum);
        
        // link new_root node with the intermediate node
        new_root->addNeighbor(second_internal_node, 0);
        second_internal_node->addNeighbor(new_root, 0);
        
        // update related info
        tree->root = new_root;
        tree->rooted = true;
        tree->leafNum++;
        tree->nodeNum++;
    }
}

/**
    compute the switching param to switch between Rate matrix and Probability matrix
*/
void AliSimulator::computeSwitchingParam(int seq_length)
{
    // don't re-set the switching param if the user has specified it
    if (params->original_params.find("--simulation-thresh") == std::string::npos) {
        double a = 1;
        // init 'a' for simulations with discrete rate variation or without rate heterogeneity
        if (!tree->getModelFactory()->is_continuous_gamma)
        {
            if (seq_length >= 1000000)
                a = 1;
            else if (seq_length >= 500000)
                a = 1.1;
            else if (seq_length >= 100000)
                a = 1.4;
            else
                a = 2.226224503;
        }
        // update 'a' for simulations with continuous rate variation
        else
        {
            if (seq_length >= 1000000)
                a = 6;
            else if (seq_length >= 500000)
                a = 7;
            else if (seq_length >= 100000)
                a = 9.1;
            else
                a = 13.3073605;
        }
        
        // compute the switching param
        params->alisim_simulation_thresh = a / seq_length;
    }
}

/**
    change state of sites due to Error model
*/
void AliSimulator::changeSitesErrorModel(vector<int> sites, vector<short int> &sequence_chunk, double error_prop, int* rstream)
{
    // estimate the total of sites need to change
    int num_changes = round(error_prop*sites.size());
    
    // randomly select a site to change
    for (int i = 0; i < num_changes; i++)
    {
        // throw error if the number of available sites is not sufficient
        if (num_changes - i > sites.size())
            outError("Cannot select a site for changing state (to simulate Sequencing Error Model). The proportion of error seems to be too high. You should try again with a smaller proportion of error!");
        
        // randomly select a site
        int selected_index = random_int(sites.size(), rstream);
        int selected_site = sites[selected_index];
        
        // remove the selected_site from the set of available sites
        sites.erase(sites.begin()+selected_index);
        
        // if the selected_site is a gap -> ignore and retry
        if (sequence_chunk[selected_site] == STATE_UNKNOWN)
            i--;
        // otherwise, randomly select a new state
        else
        {
            // select a new state
            short int new_state = random_int(max_num_states, rstream);
            while (new_state == sequence_chunk[selected_site] && max_num_states > 1)
            {
                new_state = random_int(max_num_states, rstream);
            }
            
            // update the selected site
            sequence_chunk[selected_site] = new_state;
        }
    }
}

/**
    handle DNA error
*/
void AliSimulator::handleDNAerr(int segment_start, double error_prop, vector<short int> &sequence_chunk, int* rstream, int model_index)
{
    // dummy variables
    vector<int> sites;
    
    // init vector of available sites
    // extract available sites from site_specific_model if a mixture model is used
    if (model_index >= 0 && site_specific_model_index.size()>0)
    {
        for (int i = 0; i < sequence_chunk.size(); i++)
        if (site_specific_model_index[segment_start + i] == model_index)
            sites.push_back(i);
    }
    // otherwise get all sites
    else
    {
        sites.resize(sequence_chunk.size());
        std::iota(sites.begin(),sites.end(), 0);
    }
    
    // change state of sites due to Error model
    changeSitesErrorModel(sites, sequence_chunk, error_prop, rstream);
}

/**
    TRUE if posterior mean rate can be used
*/
bool AliSimulator::canApplyPosteriorRateHeterogeneity()
{
    bool show_warning_msg = tree->params->original_params.find("--rate-heterogeneity") != std::string::npos;
    
    // without an input alignment
    if (!tree->params->alisim_inference_mode)
    {
        if (show_warning_msg)
            outWarning("Skipping Posterior Mean Rates (or sampling rates from Posterior Distribution) as they can only be used if users supply an input alignment.");
        return false;
    }
    
    // fused mixture models
    if (tree->getModel()->isMixture() && tree->getModel()->isFused())
    {
        if (show_warning_msg)
            outWarning("Skipping Posterior Mean Rates (or sampling rates from Posterior Distribution) as they cannot be used with Fused mixture models.");
        return false;
    }
    
    // without discrete rate heterogeneity
    if (tree->getRateName().find("+G") == std::string::npos && tree->getRateName().find("+R") == std::string::npos)
    {
        if (show_warning_msg)
            outWarning("Skipping Posterior Mean Rates (or sampling rates from Posterior Distribution) as they can be used with only rate heterogeneity based on a discrete Gamma/Free-rate distribution.");
        return false;
    }
    
    // continuous gamma distribution
    if ((tree->getRateName().find("+G") != std::string::npos) && tree->getModelFactory()->is_continuous_gamma)
    {
        if (show_warning_msg)
            outWarning("Skipping Posterior Mean Rates (or sampling rates from Posterior Distribution) as they cannot be used with rate heterogeneity based on a continuous Gamma distribution.");
        return false;
    }

    return true;
}

/**
    init Site to PatternID
*/
void AliSimulator::initSite2PatternID(int length)
{
    ASSERT(tree->params->alisim_inference_mode);
    
    // extract site to pattern id from the input alignment
    tree->aln->getSitePatternIndex(site_to_patternID);
    
    // resize if the input sequence length is different from the output sequence length
    int input_length = site_to_patternID.size();
    if (input_length != length)
    {
        site_to_patternID.resize(length);
        
        // randomly assign a pattern to each additional sites
        int site_id;
        for (int i = input_length; i < length; i++)
        {
            site_id = random_int(input_length);
            site_to_patternID[i] = site_to_patternID[site_id];
        }
        
    }
}

void AliSimulator::updateNewGenomeIndels(int seq_length)
{
    // dummy variables
    int rebuild_indel_his_step = params->rebuild_indel_history_param * tree->leafNum;
    int rebuild_indel_his_thresh = rebuild_indel_his_step;
    int tips_count = 0;
    
    // find the first tip that completed the simulation
    Insertion* insertion = first_insertion;
    for (; insertion && insertion->phylo_nodes.size() == 0; )
        insertion = insertion->next;
    
    ASSERT(insertion && insertion->phylo_nodes.size() > 0);
    
    // build a genome tree from the list of insertions
    GenomeTree* genome_tree = new GenomeTree();
    genome_tree->buildGenomeTree(insertion, insertion->phylo_nodes[0]->sequence->sequence_chunks[0].size(), true);
    
    // export new sequence for the first tip
    for (int i = 0; i < insertion->phylo_nodes.size(); i++)
    {
        tips_count++;
        
        insertion->phylo_nodes[i]->sequence->sequence_chunks[0] = genome_tree->exportNewGenome(insertion->phylo_nodes[i]->sequence->sequence_chunks[0], seq_length, tree->aln->STATE_UNKNOWN);
    
        // delete the insertion_pos of this node as we updated its sequence.
        insertion->phylo_nodes[i]->sequence->insertion_pos = NULL;
    }
    
    // keep track of previous insertion
    Insertion* previous_insertion = insertion;
    
    // move to next insertion
    insertion = insertion->next;
    
    // find the next tip and update it sequence
    for (; insertion; )
    {
        // if we found a tip -> update the genome_tree and export the new sequence for that tip
        if (insertion->phylo_nodes.size() > 0)
        {
            // if it is not the last tip -> update the genome tree
            if (insertion->next)
            {
                // rebuild the indel his if the number of tips (line_num) >= current threshold
                if (tips_count >= rebuild_indel_his_thresh)
                {
                    // detach the insertion and genome nodes
                    for (Insertion* tmp_insertion = insertion; tmp_insertion; )
                    {
                        // detach insertion and genome_nodes
                        tmp_insertion->genome_nodes.clear();
                        
                        // move to the next insertion
                        tmp_insertion = tmp_insertion->next;
                    }
                    
                    // delete and rebuild genome tree
                    delete genome_tree;
                    genome_tree = new GenomeTree();
                    genome_tree->buildGenomeTree(insertion, insertion->phylo_nodes[0]->sequence->sequence_chunks[0].size(), true);
                    
                    // update the next threshold to rebuild the indel his
                    rebuild_indel_his_thresh += rebuild_indel_his_step;
                }
                // otherwise, just update indel his
                else
                    genome_tree->updateGenomeTree(previous_insertion, insertion);
                
                // keep track of previous insertion
                previous_insertion = insertion;
            }
            // otherwise, it is the last tip -> the current sequence is already the latest sequence since there no more insertion occurs
            else
            {
                delete genome_tree;
                genome_tree = new GenomeTree(seq_length);
            }
            
            // export the new sequence for the current tip
            for (int i = 0; i < insertion->phylo_nodes.size(); i++)
            {
                tips_count++;
                
                insertion->phylo_nodes[i]->sequence->sequence_chunks[0] = genome_tree->exportNewGenome(insertion->phylo_nodes[i]->sequence->sequence_chunks[0], seq_length, tree->aln->STATE_UNKNOWN);
            
                // delete the insertion_pos of this node as we updated its sequence.
                insertion->phylo_nodes[i]->sequence->insertion_pos = NULL;
            }
        }
        
        // move to next insertion
        insertion = insertion->next;
    }
    
    // delete genome_tree
    delete genome_tree;
}

/**
*  reset tree (by setting the parents of all node to NULL) -> only using when simulating MSAs with openMP
*
*/
void AliSimulator::resetTree(int &max_depth, bool store_seq_at_cache, Node *node, Node *dad)
{
    // set starting node
    if (!node && !dad)
    {
        node = tree->root;
        dad = tree->root;
        
        // init depth at root
        node->sequence->depth = 0;
        max_depth = 0;
        
        // separate root sequence into chunks
        separateSeqIntoChunks(node);
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // update parent node of the current node
        (*it)->node->sequence->parent = node;
        (*it)->node->sequence->depth = node->sequence->depth + 1;
        if ((*it)->node->sequence->depth > max_depth)
            max_depth = (*it)->node->sequence->depth;
        (*it)->node->sequence->num_threads_done_simulation = 0;
        (*it)->node->sequence->num_threads_reach_barrier = 0;
        if (!store_seq_at_cache)
            (*it)->node->sequence->sequence_chunks.resize(num_threads);
        node->sequence->nums_children_done_simulation.resize(num_threads);
        for (int i = 0; i < num_threads; i++)
            node->sequence->nums_children_done_simulation[i] = 0;
        
        // browse 1-step deeper to the neighbor node
        resetTree(max_depth, store_seq_at_cache, (*it)->node, node);
    }
}

/**
    separate root sequence into chunks
*/
void AliSimulator::separateSeqIntoChunks(Node* node)
{
    if (node->sequence->sequence_chunks.size() != num_simulating_threads)
    {
        vector<short int> root_seq = node->sequence->sequence_chunks[0];
        assert(root_seq.size() == expected_num_sites);
        node->sequence->sequence_chunks.resize(num_simulating_threads);
        int default_segment_length = expected_num_sites / num_simulating_threads;
        
        // resize the first chunk from the root sequence
        node->sequence->sequence_chunks[0].resize(default_segment_length);
        
        // clone other chunks
        int starting_index = default_segment_length;
        for (int i = 1; i < num_simulating_threads; i++, starting_index += default_segment_length)
        {
            int actual_segment_length = i < num_simulating_threads - 1 ? default_segment_length : expected_num_sites - (num_simulating_threads - 1) * default_segment_length;
            node->sequence->sequence_chunks[i].resize(actual_segment_length);
            
            for (int j = 0; j < actual_segment_length; j++)
                node->sequence->sequence_chunks[i][j] = root_seq[starting_index + j];
        }
    }
}

void AliSimulator::mergeChunks(Node* node)
{
    // ignore merging if it was already merged
    if (node->sequence->sequence_chunks.size() == 1) return;
    
    // count num of sites of all chunks
    int total_sites = 0;
    for (int i = 0; i < node->sequence->sequence_chunks.size(); i++)
        total_sites += node->sequence->sequence_chunks[i].size();
    
    // merge other chunks into the first chunk
    int site_pos = node->sequence->sequence_chunks[0].size();
    node->sequence->sequence_chunks[0].resize(total_sites);
    for (int i = 1; i < node->sequence->sequence_chunks.size(); i++)
    {
        for (int j = 0; j < node->sequence->sequence_chunks[i].size(); j++)
        {
            node->sequence->sequence_chunks[0][site_pos] = node->sequence->sequence_chunks[i][j];
            site_pos++;
        }
    }
    
    // resize the sequence_chunks
    node->sequence->sequence_chunks.resize(1);
}

/**
    merge chunks into a single sequence for all nodes in tree
*/
void AliSimulator::mergeChunksAllNodes(Node* node, Node* dad)
{
    // set the starting node
    if (!node && !dad)
    {
        node = tree->root;
        dad = tree->root;
    }
    
    // merge chunks into a single sequence for the current node
    if (node->isLeaf() || params->alisim_write_internal_sequences)
        mergeChunks(node);
    
    // traverse the tree to all other nodes
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // browse 1-step deeper to the neighbor node
        mergeChunksAllNodes((*it)->node, node);
    }
}
