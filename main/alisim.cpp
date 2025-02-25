/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"

void runAliSim(Params &params, Checkpoint *checkpoint)
{
    MPIHelper::getInstance().barrier();
    auto start = getRealTime();
    auto start_cpu = getCPUTime();
    
    // Init variables
    IQTree *tree;
    Alignment *aln;
    bool inference_mode = false;
    
    // check if inference_mode is active
    // case 1: inference_mode = true if an alignment file is supplied
    if (params.aln_file)
        inference_mode = true;
    // case 2: inference_mode = true if alignment files are specified in a partition_file
    else if (params.partition_file) {
        // initilize partition alignments
        SuperAlignment *partition_aln;
        if (params.partition_type == TOPO_UNLINKED)
            partition_aln = new SuperAlignmentUnlinked(params);
        else
            partition_aln = new SuperAlignment(params);
        
        // inference_mode = true if all alignments for all partitions are supplied
        inference_mode = true;
        for (Alignment *par_aln:partition_aln->partitions)
            if (par_aln->aln_file.size() == 0)
                inference_mode = false;
        
        // delete partition_aln
        delete partition_aln;
    }
    params.alisim_inference_mode = inference_mode;
    
    // generate a random tree if neccessary
    if (params.tree_gen != NONE && MPIHelper::getInstance().isMaster())
    {
        // draw a random num_taxa from a user-defined list
        if (!params.alisim_num_taxa_list.empty())
        {
            int rand_index = random_int(params.alisim_num_taxa_list.size());
            params.sub_size = params.alisim_num_taxa_list.at(rand_index);
        }
        // draw a random num_taxa from uniform distribution
        else if (params.alisim_num_taxa_uniform_start > 3)
        {
            int range = params.alisim_num_taxa_uniform_end - params.alisim_num_taxa_uniform_start + 1;
            params.sub_size = params.alisim_num_taxa_uniform_start + random_int(range);
        }
            
        generateRandomTree(params);
        
        // reset flags to make sure IQTree will not re-generate new random starting_tree
        params.start_tree = STT_PLL_PARSIMONY;
        params.tree_gen = NONE;
    }
    // make sure the tree is generated before other MPI processes do further steps
    MPIHelper::getInstance().barrier();
    
    // inferring input parameters if inference mode is active
    if (inference_mode)
    {
        inferInputParameters(params, checkpoint, tree, aln);

        if (params.include_pre_mutations)
        {
            outWarning("Ignore predefined mutations in the input tree since it is not supported in simulations to mimick an input alignment.");
            params.include_pre_mutations = false;
        }
    }
    
    // execute AliSim Simulation
    executeSimulation(params, tree);
    
    // aln and tree are deleted in distructor of AliSimSimulator
    MPIHelper::getInstance().barrier();
    auto end = getRealTime();
    auto end_cpu = getCPUTime();
    cout << "Simulation CPU time: " << fixed << end_cpu - start_cpu << " sec (" <<
        convert_time(end_cpu-start_cpu) << ")" << endl;
    cout << "Simulation wall-clock time: " << fixed << end - start << " sec (" <<
        convert_time(end-start) << ")" << endl;
    cout << endl;
}

/**
*  inferring input parameters for AliSim
*/
void inferInputParameters(Params &params, Checkpoint *checkpoint, IQTree *&tree, Alignment *&aln)
{
    // if user has not specified model_name -> set model_name = "" (to override the default model JC)
    if ((params.original_params.find("-m ") == std::string::npos)
        || (params.original_params.find("-m TEST") != std::string::npos)
        || (params.original_params.find("-m MF") != std::string::npos))
    {
        params.model_name = "";
    }
    
    // infer model's parameters [and a tree]
    runPhyloAnalysis(Params::getInstance(), checkpoint, tree, aln);
    ASSERT(tree && tree->getModel() && tree->aln);
    
    // update model_name
    params.model_name = tree->getModel()->getNameParams();
    
    // initialize the tree_file if it has not been provided by the user
    if (!params.user_file)
    {
        // initialize the tree_file if it has not been provided by the user
        char *pre_fix = params.aln_file?params.aln_file:params.partition_file;
        params.user_file = new char[strlen(pre_fix) + 15];
        strcpy(params.user_file, pre_fix);
        if (params.partition_file && params.partition_type == BRLEN_OPTIMIZE)
            strcat(params.user_file,".parttrees");
        else
            strcat(params.user_file,".treefile");
    }
    // NhanLT: we shouldn't reload the tree
    /*else
    {
        // do not reload the tree from file in case with heterotachy
        if (!tree->getRate()->isHeterotachy())
        {
            bool is_rooted = false;
            tree->readTree(params.user_file, is_rooted);
        }
        
        // handle the case to reload super tree
        if (tree->isSuperTree()){
            // recording start_time
            auto start = getRealTime();
            
            int i;
            for (i = 0; i < ((PhyloSuperTree*) tree)->size(); i++){
                // -Q (params->partition_type == BRLEN_OPTIMIZE) -> tree_line_index = i; otherwise (-p, -q), tree_line_index = 0 (only a tree)
                int tree_line_index = 0;
                if (params.partition_type == BRLEN_OPTIMIZE)
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
                current_tree->readTree(params.user_file, is_rooted, tree_line_index);
                
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
                    current_tree->IQTree::readTree(params.user_file, is_rooted, tree_line_index);
                }
            }
            
            // show the reloading tree time
            auto end = getRealTime();
            cout<<" - Time spent on reloading trees: "<<end-start<<endl;
        }
    }*/
    
    // update sequence_length
    if (params.original_params.find("--length") == std::string::npos)
    {
        // if a partition model is used -> using the total length of all partitions
        if (tree->isSuperTree())
        {
            params.alisim_sequence_length = computeTotalSequenceLengthAllPartitions((PhyloSuperTree*) tree);
        }
        // in normal case (without partitions) -> using the sequence length in the tree
        else
        {
            params.alisim_sequence_length = (tree->aln->seq_type == SEQ_CODON ? (tree->aln->getNSite() * 3) : tree->aln->getNSite());
        }
    }
}

/**
*  compute the total sequence length of all partitions
*/
int computeTotalSequenceLengthAllPartitions(PhyloSuperTree *super_tree)
{
    int total_length = 0;
    // browse partitions one by one
    for (int i = 0; i < super_tree->size(); i++)
    {
        Alignment *aln = super_tree->at(i)->aln;
        total_length += (aln->seq_type == SEQ_CODON ? (aln->getNSite() * 3) : aln->getNSite());
    }
    return total_length;
}

/**
*  generate a random tree
*/
void generateRandomTree(Params &params)
{
    if (params.sub_size < 3 && !params.aln_file) {
        outError(ERR_FEW_TAXA);
    }

    if (!params.user_file) {
        outError("Please specify an output tree file name");
    }
    ////cout << "Random number seed: " << params.ran_seed << endl << endl;

    SplitGraph sg;

    try {

        if (params.tree_gen == YULE_HARDING || params.tree_gen == CATERPILLAR ||
            params.tree_gen == BALANCED || params.tree_gen == UNIFORM || params.tree_gen == STAR_TREE || params.tree_gen == BIRTH_DEATH) {
            if (MPIHelper::getInstance().isMaster() && fileExists(params.user_file) && !params.ignore_checkpoint)
            {
                string tmp_str(params.user_file);
                outError(tmp_str + " exists. Use `-redo` option if you want to overwrite it.");
            }
            ofstream out;
            out.open(params.user_file);
            MTree itree;

            if (params.second_tree) {
                cout << "Generating random branch lengths on tree " << params.second_tree << " ..." << endl;
                itree.readTree(params.second_tree, params.is_rooted);
            } else
            switch (params.tree_gen) {
            case YULE_HARDING:
                cout << "Generating random Yule-Harding tree..." << endl;
                break;
            case UNIFORM:
                cout << "Generating random uniform tree..." << endl;
                break;
            case CATERPILLAR:
                cout << "Generating random caterpillar tree..." << endl;
                break;
            case BALANCED:
                cout << "Generating random balanced tree..." << endl;
                break;
            case STAR_TREE:
                cout << "Generating star tree with random external branch lengths..." << endl;
                break;
            case BIRTH_DEATH:
                cout << "Generating random Birth-death tree..." << endl;
                break;
            default: break;
            }
            ofstream out2;
            if (params.num_zero_len) {
                cout << "Setting " << params.num_zero_len << " internal branches to zero length..." << endl;
                string str = params.user_file;
                str += ".collapsed";
                out2.open(str.c_str());
            }
            for (int i = 0; i < params.repeated_time; i++) {
                MExtTree mtree;
                if (itree.root) {
                    mtree.copyTree(&itree);
                    mtree.generateRandomBranchLengths(params);
                } else {
                    mtree.generateRandomTree(params.tree_gen, params);
                }
                if (params.num_zero_len) {
                    mtree.setZeroInternalBranches(params.num_zero_len);
                    MExtTree collapsed_tree;
                    collapsed_tree.copyTree(&mtree);
                    collapsed_tree.collapseZeroBranches();
                    collapsed_tree.printTree(out2);
                    out2 << endl;
                }
                mtree.printTree(out);
                out << endl;
            }
            out.close();
            cout << params.repeated_time << " tree(s) printed to " << params.user_file << endl;
            if (params.num_zero_len) {
                out2.close();
                cout << params.repeated_time << " collapsed tree(s) printed to " << params.user_file << ".collapsed" << endl;
            }
        }
        // Generate random trees if optioned
        else if (params.tree_gen == CIRCULAR_SPLIT_GRAPH) {
            cout << "Generating random circular split network..." << endl;
            if (MPIHelper::getInstance().isMaster() && fileExists(params.user_file) && !params.ignore_checkpoint)
            {
                string tmp_str(params.user_file);
                outError(tmp_str + " exists. Use `-redo` option if you want to overwrite it.");
            }
            sg.generateCircular(params);
        } else if (params.tree_gen == TAXA_SET) {
            sg.init(params);
            cout << "Generating random taxa set of size " << params.sub_size <<
                " overlap " << params.overlap << " with " << params.repeated_time << " times..." << endl;
            if (MPIHelper::getInstance().isMaster() && fileExists(params.pdtaxa_file) && !params.ignore_checkpoint)
            {
                string tmp_str(params.pdtaxa_file);
                outError(tmp_str + " exists. Use `-redo` option if you want to overwrite it.");
            }
            sg.generateTaxaSet(params.pdtaxa_file, params.sub_size, params.overlap, params.repeated_time);
        }
    } catch (bad_alloc) {
        outError(ERR_NO_MEMORY);
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, params.user_file);
    }

    // calculate the distance
    if (params.run_mode == RunMode::CALC_DIST) {
        if (params.tree_gen == CIRCULAR_SPLIT_GRAPH) {
            cout << "Calculating distance matrix..." << endl;
            sg.calcDistance(params.dist_file);
            cout << "Distances printed to " << params.dist_file << endl;
        }// else {
            //mtree.calcDist(params.dist_file);
        //}
    }

}

std::vector<std::pair<std::string,std::string>> readMutations(const std::string& mutation_file)
{
    // init a vector of node mutations
    std::vector<std::pair<std::string,std::string>> node_mutations;

    // open file, read line by line
    try
    {
        // open file
        if (!fileExists(mutation_file))
            outError("File not found: " + mutation_file);
        std::ifstream file(mutation_file);
        std::string line;

        // read line by line
        while (std::getline(file, line))
        {
            // ignore empty line
            if (!line.length()) continue;

            // remove spaces at the begining of the line
            auto pos = line.find_first_not_of(" ");
            // find the first non-space character
            if (pos != string::npos)
            {
                // if that character is not at the begining of the line -> remove all the spaces at the begining
                if (pos > 0)
                    line = line.substr(pos, line.length() - pos);
            }
            // if not found -> all characters are spaces -> invalid
            else
                outError("Invalid format '" + line + "'. Expected format should be <node_name><spaces><list_of_mutations>");

            // remove spaces at the ending of the line
            pos = line.find_last_not_of(" ");
            // find the last non-space character
            if (pos != string::npos)
            {
                // if that character is not at the end of the line -> remove all the spaces at the end
                if (pos < line.length() - 1)
                    line = line.substr(0, pos + 1);
            }
            // if not found -> all characters are spaces -> invalid
            else
                outError("Invalid format '" + line + "'. Expected format should be <node_name><spaces><list_of_mutations>");

            // validate the input format <node_name><spaces><list_of_mutations>
            if (line.length() < 5)
                outError("Invalid format '" + line + "'. Expected format should be <node_name><spaces><list_of_mutations>");

            // replace \t by a space
            line = regex_replace(line, std::regex("\t"), " ");

            // get the node name
            std::string node_name = "";
            pos = line.find(" ");
            if (pos != std::string::npos)
            {
                node_name = line.substr(0, pos);

                // convert node_name to uppercase
                transform(node_name.begin(), node_name.end(), node_name.begin(), ::toupper);
            }
            else
                outError("Invalid format '" + line + "'. Expected format should be <node_name><spaces><list_of_mutations>");

            // get mutation list at that node
            std::string mutation_str = "";
            // find the first non-space character after the position
            for (; pos < line.length(); ++pos)
                if (line[pos] != ' ')
                    break;
            if (line.length() - pos)
                mutation_str = line.substr(pos, line.length() - pos);
            else
                outError("Invalid format '" + line + "'. Expected format should be <node_name><spaces><list_of_mutations>");

            // record the mutations at that node
            node_mutations.push_back(std::pair<std::string,std::string>(node_name, mutation_str));
        }

        // close file
        file.close();
    }
    catch(std::exception e)
    {
        outError("File not found or invalid format " + mutation_file);
    }

    // return result
    return node_mutations;
}

/**
*  execute AliSim Simulation
*/
void executeSimulation(Params& params, IQTree *&tree)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // disable posterior mean rate (or sampling rate from posterior distribution) if users don't supply input alignment
    if (params.alisim_rate_heterogeneity!=UNSPECIFIED && !params.alisim_inference_mode)
    {
        params.alisim_rate_heterogeneity = UNSPECIFIED;
        
        // show a warning if the user specifies pos_mean rate (or sampling rate from posterior distribution) but doesn't supply the input alignment
        if (params.original_params.find("--rate-heterogeneity") != std::string::npos) {
            outWarning("Skipping --rate-heterogeneity option as it can only be used if users supply an input alignment.");
        }
    }
    
    // disable posterior mean site freqs (or sampling site freqs from posterior distribution) if users don't supply input alignment
    if (params.alisim_stationarity_heterogeneity!=UNSPECIFIED && !params.alisim_inference_mode)
    {
        params.alisim_stationarity_heterogeneity = UNSPECIFIED;
        
        // show a warning if the user specifies posterior mean site freqs (or sampling site freqs from posterior distribution) but doesn't supply the input alignment
        if (params.original_params.find("--state-freqs") != std::string::npos) {
            outWarning("Skipping --state-freqs option as it can only be used if users supply an input alignment.");
        }
    }
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator;
    if (tree && params.alisim_inference_mode)
        alisimulator = new AliSimulator(&params, tree);
    else
        alisimulator = new AliSimulator(&params);
    
    // only unroot tree and stop if the user just wants to do so
    if (alisimulator->params->alisim_only_unroot_tree)
    {
        unrootTree(alisimulator);
        
        // stop the program after unrooting the tree
        return;
    }
    
    // show parameters
    showParameters(params, alisimulator->tree->isSuperTree());
    
    // export tree with new blengths if users want to do so
    if (params.branch_distribution && params.user_file && !params.alisim_inference_mode)
    {
        string tree_path(params.user_file);
        tree_path += ".new_blength";
        std::cout << "Tree with randomly generated branch lengths is outputted at " << tree_path << std::endl;
        ofstream out = ofstream(tree_path.c_str());
        alisimulator->tree->printTree(out);
        if (alisimulator->tree->isSuperTree() && params.partition_type == BRLEN_OPTIMIZE)
        {
            for (int i = 1; i < ((PhyloSuperTree*) alisimulator->tree)->size(); i++)
            {
                out << std::endl;
                ((PhyloSuperTree*) alisimulator->tree)->at(i)->printTree(out);
            }
        }
        out.close();
    }
    
    // load input MSA if any
    map<string,string> input_msa = loadInputMSA(alisimulator);
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    generateMultipleAlignmentsFromSingleTree(alisimulator, input_msa);
    
    // show log file
    if (!(params.suppress_output_flags & OUT_LOG))
        cout << "Screen log file: " << params.out_prefix << ".log" << endl;
    
    // delete alisimulator
    if (alisimulator->tree) delete alisimulator->tree;
    if (alisimulator->first_insertion) delete alisimulator->first_insertion;
    delete alisimulator;
    
    cout << "[Alignment Simulator] Done"<<"\n";
}

/**
*  show all input parameters for AliSim
*/
void showParameters(Params &params, bool is_partition_model)
{
    cout << " - Tree filepath: " << params.user_file <<"\n";
    cout << " - Length of output sequences: " << params.alisim_sequence_length <<"\n";
    if (!params.model_name.empty())
    {
        if (is_partition_model)
            cout << " - Model: " << "Partition model" <<"\n";
        else
            cout << " - Model: " << params.model_name <<"\n";
    }
    cout << " - Number of output datasets: " << params.alisim_dataset_num<<"\n";
}

/**
*  retrieve the ancestral sequence for the root node from an input file
*/
void retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator, vector<short int> &sequence)
{
    // get variables
    char* aln_filepath = new char[super_alisimulator->params->root_ref_seq_aln.length() + 1];
    strcpy(aln_filepath, super_alisimulator->params->root_ref_seq_aln.c_str());
    string sequence_name = super_alisimulator->params->root_ref_seq_name;
    
    // in normal case (without partition) -> using the current tree to load the ancestral sequence
    IQTree *src_tree = super_alisimulator->tree;
    // in case with partitions -> using the first phylotree to load the ancestral sequence
    if (src_tree->isSuperTree())
    {
        // make sure all partitions are using the same sequence_type
        for (int i = 1; i < ((PhyloSuperTree*) src_tree)->size(); i++)
            if (((PhyloSuperTree*) src_tree)->at(i)->aln->seq_type != ((PhyloSuperTree*) src_tree)->at(0)->aln->seq_type)
                outError("To load ancestral sequence from a file, all partitions must use the same sequence_type.");
        
        // using the first phylotree to load the ancestral sequence
        src_tree = ((IQTree*)((PhyloSuperTree*) super_alisimulator->tree)->at(0));
    }

    // read sequences from the input file
    Alignment *aln = new Alignment();
    StrVector sequences;
    int nseq = 0, nsite = 0;
    char *sequence_type = strcpy(new char[src_tree->aln->sequence_type.length() + 1], src_tree->aln->sequence_type.c_str());
    aln->extractSequences(aln_filepath, sequence_type, sequences, nseq, nsite);
    StrVector seq_names = aln->getSeqNames();
    delete[] aln_filepath;

    // delete aln
    delete aln;

    string sequence_str = "";
    for (int i = 0; i < seq_names.size(); i++)
        if (!sequence_name.compare(seq_names[i]))
        {
            sequence_str = sequences[i];
            break;
        }
    if (sequence_str.length() == 0)
        outError("Sequence name could not be found in the input alignment file.");
    
    // overwrite the output sequence_length
    if (super_alisimulator->params->alisim_sequence_length != sequence_str.length())
    {
        super_alisimulator->params->alisim_sequence_length = sequence_str.length();
        outWarning("Sequence length is now set equally to the length of ancestral sequence.");
        super_alisimulator->refreshExpectedNumSites();
    }

    // get Max number of states
    int max_num_states = src_tree->aln->getMaxNumStates();
    
    // convert the input sequence into (numerical states) sequence
    int num_sites_per_state = src_tree->aln->seq_type == SEQ_CODON ? 3 : 1;
    int sequence_length = (src_tree->aln->seq_type == SEQ_CODON ? (super_alisimulator->params->alisim_sequence_length / 3) : super_alisimulator->params->alisim_sequence_length);
    
    // make sure the length of the ancestral sequence must be equal to the total length of all partitions
    if (super_alisimulator->tree->isSuperTree() && sequence_length != super_alisimulator->tree->getAlnNSite())
        outError("The length of the ancestral sequence must be equal to the total length of all partitions");
    
    sequence.resize(sequence_length);
    ostringstream err_str;
    int num_error = 0;
    if (src_tree->aln->seq_type == SEQ_CODON)
    {
        int site_index = 0;
        for (int i = 0; i < sequence_length; i++, site_index += num_sites_per_state)
        {
            // NHANLT: potential improvement
            // change to use pointer of sequence_str to avoid accessing sequence_str[]
            sequence[i] = src_tree->aln->getCodonStateTypeFromSites(src_tree->aln->convertState(sequence_str[site_index], SEQ_DNA), src_tree->aln->convertState(sequence_str[site_index+1], SEQ_DNA), src_tree->aln->convertState(sequence_str[site_index+2], SEQ_DNA), sequence_name, site_index, err_str, num_error);
        }
    }
    else
    {
        for (int i = 0; i < sequence_length; i++)
        {
            sequence[i] = src_tree->aln->convertState(sequence_str[i]);
        
            // Handle invalid/unknown state
            if (sequence[i] >= max_num_states)
                sequence[i] = random_int(max_num_states);
        }
    }
    
    // show error from ancestral sequence (if any)
    if (num_error)
        outError(err_str.str());
    
    // show warning
    outWarning("Using an ancestral sequence with base frequencies that are not compatible with the specification of the model may lead to unexpected results.");
}

void getLockedSites(Node* const node, Node* const dad, std::vector<bool>* const site_locked_vec, Alignment* const aln)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // get locked sites from the predefined mutations at this branch (if any)
        auto atb_it = (*it)->attributes.find(MTree::ANTT_MUT);
        if (atb_it != (*it)->attributes.end())
        {
            // sequence length
            const int seq_length = site_locked_vec->size();

            // parse a list of mutations
            Substitutions pre_mutations = Substitutions(atb_it->second, aln, seq_length);

            // mark those sites locked
            for (auto mut_it = pre_mutations.begin(); mut_it != pre_mutations.end(); ++mut_it)
            {
                // extract position
                const int pos = mut_it->getPosition();

                // vailidate position
                if (pos >= seq_length)
                {
                    outWarning("Ignore a predefined mutation " + aln->convertStateBackStr(mut_it->getOldState()) + convertIntToString((aln->seq_type == SEQ_CODON ? pos * 3 : pos) + Params::getInstance().site_starting_index) + aln->convertStateBackStr(mut_it->getNewState()) + ". Position exceeds the sequence length " + convertIntToString(aln->seq_type == SEQ_CODON ? seq_length * 3 : seq_length));
                }
                // mark the site as locked
                else
                    site_locked_vec->at(pos) = true;
            }
        }

        // browse 1-step deeper to the neighbor node
        getLockedSites((*it)->node, node, site_locked_vec, aln);
    }
}

void createNodeMapping(std::map<std::string, std::pair<Node*, Node*>>& node_mapping, Node* const node, Node* const dad)
{
    // add the current node to the mapping
    if (node->name.length())
    {
        std::string node_name = node->name;
        // convert to uppercase
        transform(node_name.begin(), node_name.end(), node_name.begin(), ::toupper);

        node_mapping.insert(std::pair<std::string, std::pair<Node*, Node*>>(node_name, std::pair<Node*, Node*>(dad, node)));
    }
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        createNodeMapping(node_mapping, (*it)->node, node);
    }
}

void addMutations2Tree(const std::vector<std::pair<std::string, std::string>>& node_mutations, IQTree* const tree)
{
    // validate input
    ASSERT(tree && tree->root);

    // create a mapping between each node name and a pair of pointers <dad_node, node>
    std::map<std::string, std::pair<Node*, Node*>> node_mapping;
    createNodeMapping(node_mapping, tree->root, NULL);

    // browse the list of node_mutations to add mutations of each node to the corresponding node in the tree
    for (const std::pair<std::string, std::string>& mutations : node_mutations)
    {
        // seek the corresponding node
        auto it = node_mapping.find(mutations.first);
        if( it != node_mapping.end())
        {
            // Extract the corresponding pair of nodes
            Node* dad = (it->second).first;
            Node* node = (it->second).second;

            // add attribute (mutations,<mutations_list> to the corresponding branch
            dad->findNeighbor(node)->putAttr(MTree::ANTT_MUT, "{" + mutations.second + "}");
            node->findNeighbor(dad)->putAttr(MTree::ANTT_MUT, "{" + mutations.second + "}");
        }
        // node_name is not found
        else
            outWarning("Parsing predefined mutations. Node " + mutations.first + " is not found in the tree.");
    }
}

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator, map<string,string> input_msa)
{
    // Load ancestral sequence from the input file if user has specified it
    vector<short int> ancestral_sequence;
    if (super_alisimulator->params->root_ref_seq_name.length() > 0)
        retrieveAncestralSequenceFromInputFile(super_alisimulator, ancestral_sequence);
    
    // terminate if users employ more MPI processes than the number of alignments
    if (MPIHelper::getInstance().getNumProcesses() > super_alisimulator->params->alisim_dataset_num)
        outError("You are employing more MPI processes (" + convertIntToString(MPIHelper::getInstance().getNumProcesses()) + ") than the number of alignments (" + convertIntToString(super_alisimulator->params->alisim_dataset_num) + "). Please reduce the number of MPI processes to save the computational resources and try again!");
    
    // BUG FIXED: set auto num_threads to the max #cores
#ifdef _OPENMP
    // num_threads == 0 <=> auto
    if (!super_alisimulator->params->num_threads)
    {
        super_alisimulator->params->num_threads = countPhysicalCPUCores();
        
        // manually set number of threads
        omp_set_num_threads(super_alisimulator->params->num_threads);
    }
    
    // show info
    if (MPIHelper::getInstance().getNumProcesses() == 1)
        cout << " - Number of threads: " << super_alisimulator->params->num_threads << endl;
    else {
        cout << " - Number of threads per MPI process: " << super_alisimulator->params->num_threads << endl;
        cout << " - Number of threads for MPI processes: " << super_alisimulator->params->num_threads * MPIHelper::getInstance().getNumProcesses() << endl;
    }
#endif

    // reset number of OpenMP threads to 1 in simulations with Indels
    if (super_alisimulator->params->num_threads != 1 && super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0)
    {
        outWarning("OpenMP has not yet been supported in simulations with Indels. AliSim is now using a single thread for this simulation.");
        super_alisimulator->params->num_threads = 1;
#ifdef _OPENMP
        omp_set_num_threads(super_alisimulator->params->num_threads);
#endif
    }
    
    // do not support compression when outputting multiple data sets into a same file
    if (Params::getInstance().do_compression && (Params::getInstance().alisim_single_output || super_alisimulator->params->num_threads != 1))
    {
        outWarning("Compression is not supported when either outputting multiple alignments into a single output file or using multithreading. AliSim will output file in normal format.");

        Params::getInstance().do_compression = false;
        super_alisimulator->params->do_compression = false;
    }
    
    // cannot skip concatenating sequence chunks from intermediate files in simulations with FunDi, Partitions, or +ASC models
    if (Params::getInstance().num_threads != 1 && Params::getInstance().no_merge)
    {
        // ignore --no-merge option if using AliSim-OpenMP-IM
        if (Params::getInstance().alisim_openmp_alg == IM)
            outWarning("Ignore --no-merge option as it is only appliable for AliSim-OpenMP-EM algorithm.");
        // otherwise, if using AliSim-OpenMP-EM -> show a warning in cases that we cannot skip merging
        else if (super_alisimulator->tree->isSuperTree()
                  || super_alisimulator->params->alisim_fundi_taxon_set.size() > 0
                  || (super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
                  || super_alisimulator->params->aln_output_format == IN_MAPLE)
        {
            outWarning("Cannot skip merging sequence chunks in simulations with FunDi, Partitions, +ASC models, or when outputting alignment in MAPLE format. AliSim will concatenate sequence chunks from intermediate files into a single output file.");
            
            Params::getInstance().no_merge = false;
            super_alisimulator->params->no_merge = false;
        }
    }
    
    // show a warning if the user wants to write internal sequences in not-supported cases
    if (super_alisimulator->params->alisim_write_internal_sequences
        &&((super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
           || super_alisimulator->tree->isSuperTree() || (super_alisimulator->params->alisim_fundi_taxon_set.size() > 0 && Params::getInstance().num_threads != 1)))
    {
        outWarning("Could not write out the internal sequences when using partition, ASC models, or FunDi model with multithreading. Only sequences at tips will be written to the output file.");
        Params::getInstance().alisim_write_internal_sequences = false;
        super_alisimulator->params->alisim_write_internal_sequences = false;
    }
    
    // check to output a single file
    if (super_alisimulator->params->alisim_single_output && super_alisimulator->params->alisim_dataset_num == 1)
            super_alisimulator->params->alisim_single_output = false;
    
    // don't allow --no-merge and --single-output
    if (super_alisimulator->params->alisim_openmp_alg == EM && super_alisimulator->params->alisim_single_output && super_alisimulator->params->no_merge)
    {
        outWarning("Ignore --single-output option since it is not supported if using with --no-merge option.");
        super_alisimulator->params->alisim_single_output = false;
    }
    
    // don't allow --single-output when outputting MAPLE format
    if (super_alisimulator->params->alisim_single_output && super_alisimulator->params->aln_output_format == IN_MAPLE)
    {
        outWarning("Ignore --single-output option since it is not supported when outputting MAPLE format.");
        super_alisimulator->params->alisim_single_output = false;
    }
        
    // ignore --single-output in version with MPI
#ifdef _IQTREE_MPI
    if (super_alisimulator->params->alisim_single_output)
    {
        outWarning("Ignore --single-output option since it is not supported in IQ-TREE version with MPI. Alignments will be outputted in separated files.");
        super_alisimulator->params->alisim_single_output = false;
    }
#endif
    
    // ignore predefined mutations if using the following models: Partitions, Indels
    // a vector of site status, which denotes whether a site is locked or not. Locked means that site only evolves according to the predefined mutations
    std::vector<bool>* site_locked_vec = nullptr;
    if (super_alisimulator->params->include_pre_mutations)
    {
        if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0 || // Indels
            super_alisimulator->tree->isSuperTree()) // Partitions
        {
            outWarning("Ignore predefined mutations in the input tree since it is not supported in simulations with Indels or Partition models.");
            super_alisimulator->params->include_pre_mutations = false;
        }
        else
        {
            ASSERT(super_alisimulator->tree && super_alisimulator->tree->root && super_alisimulator->tree->aln);

            // show info
            std::cout << "Predefined mutations detected" << std::endl;

            // load predefined mutations from a file (if specified)
            if (super_alisimulator->params->mutation_file.length())
            {
                std::vector<std::pair<std::string, std::string>> node_mutations = readMutations(super_alisimulator->params->mutation_file);

                // add mutations to tree
                if (node_mutations.size())
                    addMutations2Tree(node_mutations, super_alisimulator->tree);
            }

            // init site_locked_vec
            site_locked_vec = new std::vector<bool>(super_alisimulator->expected_num_sites, false);

            // browse the tree to mark all locked sites
            getLockedSites(super_alisimulator->tree->root, NULL, site_locked_vec, super_alisimulator->tree->aln);
        }
    }

    // the output format of the simulated alignment
    InputType actual_output_format = super_alisimulator->params->aln_output_format;
    vector<SeqType> seqtypes;
    vector<std::string> aln_names;
    // If users want to output Maple format -> temporarily output PHYLIP first
    if (actual_output_format == IN_MAPLE)
    {
        super_alisimulator->params->aln_output_format = IN_PHYLIP;
        Params::getInstance().aln_output_format = IN_PHYLIP;
    }
    
    // iteratively generate multiple datasets for each tree
    int proc_ID = MPIHelper::getInstance().getProcessID();
    int nprocs  = MPIHelper::getInstance().getNumProcesses();
    for (int i = proc_ID; i < super_alisimulator->params->alisim_dataset_num; i+=nprocs)
    {
        // parallelize over MPI ranks statically
        //if (i%nprocs != proc_ID) continue;
        
        // If users want to output Maple format -> clear seqtypes and aln_names
        if (actual_output_format == IN_MAPLE)
        {
            seqtypes.clear();
            aln_names.clear();
        }
        
        // record the alignment_id to generate different random seed when simulating different alignment
        super_alisimulator->params->alignment_id = i;
        
        // output the simulated aln at the current execution localtion
        string output_filepath = super_alisimulator->params->alisim_output_filename;
        
        // only add alignment id if users want to generate multiple alignments
        if (super_alisimulator->params->alisim_dataset_num > 1 && !super_alisimulator->params->alisim_single_output)
            output_filepath = output_filepath+"_"+convertIntToString(i+1);
        
        // check whether we should write output to a new file or append it into an existing one
        std::ios_base::openmode open_mode = std::ios_base::out;
        if (i > 0 && super_alisimulator->params->alisim_single_output)
            open_mode = std::ios_base::in|std::ios_base::out|std::ios_base::ate;
        
        // generate multiple alignments one by one
        if (super_alisimulator->tree->isSuperTree())
        {
            PhyloSuperTree* super_tree = ((PhyloSuperTree*) super_alisimulator->tree);
            int total_expected_num_states = super_tree->getAlnNSite();
            
            // override sequence_length by the total length of all partitions
            int total_length = computeTotalSequenceLengthAllPartitions(super_tree);
            if (super_alisimulator->tree->params->alisim_sequence_length != total_length)
            {
                super_alisimulator->tree->params->alisim_sequence_length = total_length;
                outWarning("The sequence_length is now set equally to the total length of all partitions");
                super_alisimulator->refreshExpectedNumSites();
            }
            
            for (int partition_index = 0; partition_index < super_tree->size(); partition_index++)
            {
                // update the alignment_id, taking into account the partition index, so that we use different random seed for each partition in each alignment
                super_alisimulator->params->alignment_id = (i + 1) * 1000000 + partition_index;
                // get variables
                IQTree *current_tree = (IQTree*) super_tree->at(partition_index);
                int expected_num_states_current_tree = current_tree->aln->getNSite();
                int num_sites_per_state = super_tree->at(partition_index)->aln->seq_type == SEQ_CODON?3:1;
                
                // create position_spec in case aln_files are specified in a directory
                if (super_alisimulator->params->partition_file && isDirectory(super_alisimulator->params->partition_file))
                {
                    int total_num_states = super_tree->at(partition_index)->aln->getNSite();
                    if (num_sites_per_state != 1)
                        total_num_states *= num_sites_per_state;
                    ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec = "1-" + convertIntToString(total_num_states);
                }
                
                // extract site_ids of the partition
                string info_spec_str = ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec;
                // convert position_spec from "*" to "start-end"
                if (!info_spec_str.compare("*") && super_tree->at(partition_index)->aln->getNSite() > 0)
                {
                    int total_num_states = super_tree->at(partition_index)->aln->getNSite();
                    if (num_sites_per_state != 1)
                        total_num_states *= num_sites_per_state;
                    
                    info_spec_str = "1-" + convertIntToString(total_num_states);
                    ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec = info_spec_str;
                }
                
                const char* info_spec = info_spec_str.c_str();
                IntVector site_ids;
                extractSiteID(current_tree->aln, info_spec, site_ids, false, total_expected_num_states);

                // extract the ancestral sequence for the current partition from the full ancestral_sequence
                vector<short int> ancestral_sequence_current_tree;
                if (ancestral_sequence.size() > 0)
                {
                    ASSERT(site_ids.size() == expected_num_states_current_tree);
                    ancestral_sequence_current_tree.resize(expected_num_states_current_tree);
                    
                    // extract sites one by one from the full ancestral_sequence
                    for (int j = 0; j < ancestral_sequence_current_tree.size(); j++)
                        ancestral_sequence_current_tree[j] = ancestral_sequence[site_ids[j]];
                }
                
                // stree->part_info[part].part_rate
                double partition_rate = super_tree->params->partition_type == BRLEN_SCALE ? super_tree->part_info[partition_index].part_rate:1;
                // generate alignment for the current tree/partition
                AliSimulator* partition_simulator = new AliSimulator(super_tree->params, current_tree, expected_num_states_current_tree, partition_rate);
                generatePartitionAlignmentFromSingleSimulator(partition_simulator, ancestral_sequence_current_tree, input_msa, site_locked_vec);
                
                // update new genome at tips from the original genome and the genome tree
                // skip updating if using +ASC or Fundi model as they must be already updated
                if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0 && !(partition_simulator->tree->getModelFactory() && partition_simulator->tree->getModelFactory()->getASC() != ASC_NONE) && (partition_simulator->params->alisim_fundi_taxon_set.size() == 0))
                    partition_simulator->updateNewGenomeIndels(partition_simulator->seq_length_indels);
                
                // delete partition_simulator
                if (partition_simulator->first_insertion) delete partition_simulator->first_insertion;
                delete partition_simulator;
            }
        }
        else
        {
            // record the seqtype and alignment names, which will be used later to convert the simulated alignment into Maple format
            if (actual_output_format == IN_MAPLE)
            {
                seqtypes.push_back(super_alisimulator->tree->aln->seq_type);
                aln_names.push_back(output_filepath);
            }
            
            // check whether we could write the output to file immediately after simulating it
            if (super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() == ASC_NONE && super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio == 0)
                generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence, input_msa, site_locked_vec, output_filepath, open_mode);
            // otherwise, writing output to file after completing the simulation
            else
                generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence, input_msa, site_locked_vec);
        }
        
        // merge & write alignments to files if they have not yet been written
        if ((super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
            || super_alisimulator->tree->isSuperTree()
            || super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0)
            mergeAndWriteSequencesToFiles(output_filepath, super_alisimulator, seqtypes, aln_names, open_mode);
        
        // only report model params when simulating the first MSA
        if (i == 0)
        {
            // report model's parameters
            reportSubstitutionProcess(cout, *(super_alisimulator->params), *(super_alisimulator->tree));
            // show omega/kappa/kappa2 when using codon models
            if (super_alisimulator->tree->aln->seq_type == SEQ_CODON)
                super_alisimulator->tree->getModel()->writeInfo(cout);
        }
        
        // remove tmp_data if using Indels
        if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0)
            remove((super_alisimulator->params->alisim_output_filename + "_" + super_alisimulator->params->tmp_data_filename + "_" + convertIntToString(MPIHelper::getInstance().getProcessID())).c_str());
        
        // if users want to output Maple format -> convert PHY into MAPLE and delete PHY
        if (actual_output_format == IN_MAPLE)
        {
            for (auto aln_id = 0 ; aln_id < aln_names.size(); ++ aln_id)
            {
                // initialize a dummy alignment to make sure we'll not change the main alignment when converting the simulated alignment files into Maple format
                Alignment aln;
                aln.seq_type = seqtypes[aln_id];
                
                // convert the simulated alignment files into Maple format
                aln.extractMapleFile(aln_names[aln_id], IN_PHYLIP);
                
                // remove the simulated alignment files (in PHYLIP format)
                remove(getOutputNameWithExt(IN_PHYLIP, aln_names[aln_id]).c_str());
                
                // show the output file name
                cout << "An alignment written to " << getOutputNameWithExt(IN_MAPLE, aln_names[aln_id]) << endl;
            }
        }
        // otherwise print the output file name
        else
        {
            string output_filename = output_filepath;
            
            if (super_alisimulator->params->num_threads != 1 && super_alisimulator->params->alisim_openmp_alg == EM && super_alisimulator->params->no_merge)
            {
                cout << "An alignment has been written to files: "
                << getOutputNameWithExt(super_alisimulator->params->aln_output_format, output_filename + "_1")
                << " - "
                << getOutputNameWithExt(super_alisimulator->params->aln_output_format, output_filename + "_" + convertIntToString(super_alisimulator->params->num_threads)) << endl << endl;
            }
            // each simulated alignment is outputted into a single file
            else if (!super_alisimulator->params->alisim_single_output)
            {
                cout << "An alignment written to "
                << getOutputNameWithExt(super_alisimulator->params->aln_output_format, output_filename) << endl;
                
                // if using indels and outputting unaligned sequences
                if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0
                    && !super_alisimulator->params->alisim_no_export_sequence_wo_gaps)
                {
                    cout << "Unaligned sequences written to "
                    << getOutputNameWithExt(IN_FASTA, output_filename + ".unaligned") << endl;
                }
                
                // add an empty line
                //cout << endl;
            } else if (super_alisimulator->params->alisim_single_output
                       && i == super_alisimulator->params->alisim_dataset_num - 1)
            {
                cout << super_alisimulator->params->alisim_dataset_num
                << " alignments written to "
                << getOutputNameWithExt(super_alisimulator->params->aln_output_format, output_filename) << endl;
                
                // if using indels and outputting unaligned sequences
                if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0
                    && !super_alisimulator->params->alisim_no_export_sequence_wo_gaps)
                {
                    cout << "Unaligned sequences written to "
                    << getOutputNameWithExt(IN_FASTA, output_filename + ".unaligned") << endl;
                }
            }
        }
        
        // delete output alignments (for testing only)
        if (super_alisimulator->params->delete_output)
        {
            string output_filename = output_filepath;
            
            for (int thread_id = 0; thread_id < super_alisimulator->params->num_threads; thread_id++)
            {
                if (super_alisimulator->params->num_threads != 1 && super_alisimulator->params->alisim_openmp_alg == EM && super_alisimulator->params->no_merge)
                    output_filename = output_filepath + "_" + convertIntToString(thread_id + 1);
                
                // add file extension
                output_filename = getOutputNameWithExt(super_alisimulator->params->aln_output_format, output_filename);
                
                // delete the output file
                remove((output_filename).c_str());
                
                // stop deleting if output file was merging
                if (super_alisimulator->params->alisim_openmp_alg == IM || !super_alisimulator->params->no_merge)
                    break;
            }
        }
    }
    
    // output full tree (with internal node names) if outputting internal sequences
    if (super_alisimulator->params->alisim_write_internal_sequences)
        outputTreeWithInternalNames(super_alisimulator);

    // delete site_locked_vec (if necessary)
    if (site_locked_vec)
        delete site_locked_vec;
}

/**
    copy sequences of leaves from a partition tree to super_tree
*/
void copySequencesToSuperTree(IntVector &site_ids, int expected_num_states_super_tree, IQTree *current_tree, int initial_state, Node *node, Node *dad){
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // find the corresponding node in the current_tree
        Node *current_node = current_tree->findLeafName(node->name);

        // initialize sequence of the super_node
        if (node->sequence->sequence_chunks[0].size() != expected_num_states_super_tree)
        {
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            if (node->sequence->sequence_chunks[0].size() != expected_num_states_super_tree)
                node->sequence->sequence_chunks[0].resize(expected_num_states_super_tree, initial_state);
        }
        
        // copy sequence from the current_node to the super_node (if the current node is found)
        if (current_node)
        {
            // copy sites one by one from the current sequence to its position in the sequence of the super_node
            for (int i = 0; i < site_ids.size(); i++)
                node->sequence->sequence_chunks[0][site_ids[i]] = current_node->sequence->sequence_chunks[0][i];
        }
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        copySequencesToSuperTree(site_ids, expected_num_states_super_tree, current_tree, initial_state, (*it)->node, node);
    }
}

/**
*  generate a partition alignment from a single simulator
*/
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *&alisimulator, vector<short int> &ancestral_sequence, map<string,string> input_msa, std::vector<bool>* const site_locked_vec, string output_filepath, std::ios_base::openmode open_mode)
{
    // show an error if continuous gamma is used in inference mode.
    if (alisimulator->params->alisim_inference_mode && alisimulator->tree->getModelFactory() && alisimulator->tree->getModelFactory()->is_continuous_gamma)
        outError("Unfortunately, IQ-Tree has not yet supported Continuous Gamma in phylogeny inference. Therefore, users can only use Continuous Gamma in AliSim without Inference Mode.");
    
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    bool is_mixture_model = alisimulator->tree->getModel()->isMixture();
    
    // case 1: without rate heterogeneity or mixture model -> using the current alisimulator (don't need to re-initialize it)
    AliSimulator *tmp_alisimulator = alisimulator;
    
    // case 2: with rate heterogeneity or mixture model
    if ((!rate_name.empty()) || is_mixture_model)
    {
        // if user specifies +I without invariant_rate -> set it to 0
        if (rate_name.find("+I") != std::string::npos && isnan(invariant_proportion)) {
            tmp_alisimulator->tree->getRate()->setPInvar(0);
            outWarning("Invariant rate is now set to Zero since it has not been specified");
        }
        
        // case 2.3: with only invariant sites (without gamma/freerate model/mixture models)
        if (!rate_name.compare("+I") && !is_mixture_model)
        {
            tmp_alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
        else
        {
            // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
            if (invariant_proportion > 0)
            {
                tmp_alisimulator = new AliSimulatorHeterogeneityInvar(alisimulator, invariant_proportion);
            }
            // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
            else
            {
                tmp_alisimulator = new AliSimulatorHeterogeneity(alisimulator);
            }
        }
    }
    
    tmp_alisimulator->generatePartitionAlignment(ancestral_sequence, input_msa, site_locked_vec, output_filepath, open_mode);
    
    // clone indel data before deleting tmp_alisimulator
    if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0)
    {
        alisimulator->seq_length_indels = tmp_alisimulator->seq_length_indels;
        
        // tmp_alisimulator != alisimulator
        if ((!rate_name.empty()) || is_mixture_model)
        {
            // Bug fixed - only move map_seqname_node if tmp_alisimulator != alisimulator
            alisimulator->map_seqname_node = std::move(tmp_alisimulator->map_seqname_node);
            if (alisimulator->first_insertion) delete alisimulator->first_insertion;
            alisimulator->first_insertion = tmp_alisimulator->first_insertion;
        }
    }
    
    // delete tmp_alisimulator
    if ((!rate_name.empty()) || is_mixture_model)
    {
        delete tmp_alisimulator;

        // bug fixes: avoid accessing to deallocated pointer
        if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0)
            alisimulator->first_insertion = nullptr;
    }

}

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, Alignment *aln, int sequence_length, int num_leaves, AliSimulator *alisimulator, std::ios_base::openmode open_mode)
{
    try {
            // init output_stream for Indels to output aln without gaps
            ostream *out_indels = NULL;
            bool write_indels_output = false;
            if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0
                && !alisimulator->params->alisim_no_export_sequence_wo_gaps)
            {
                write_indels_output = true;
                if (alisimulator->params->do_compression)
                    out_indels = new ogzstream((file_path+".unaligned.fa").c_str(), open_mode);
                else
                    out_indels = new ofstream((file_path+".unaligned.fa").c_str(), open_mode);
            }
        
            // add ".phy" or ".fa" to the output_filepath
            file_path = getOutputNameWithExt(alisimulator->params->aln_output_format, file_path);
            ostream *out;
            if (alisimulator->params->do_compression)
                out = new ogzstream(file_path.c_str(), open_mode);
            else
                out = new ofstream(file_path.c_str(), open_mode);
            out->exceptions(ios::failbit | ios::badbit);

            // write the first line <#taxa> <length_of_sequence> (for PHYLIP output format)
            int seq_length_times_num_sites_per_state = (aln->seq_type == SEQ_CODON ? (sequence_length * 3) : sequence_length);
            string first_line = "";
            uint64_t start_pos = 0;
            if (alisimulator->params->aln_output_format == IN_PHYLIP)
            {
                first_line = convertIntToString(num_leaves) + " " + convertIntToString(seq_length_times_num_sites_per_state) + "\n";
                *out << first_line;
                
                // get the position to write output
                start_pos = first_line.length();
            }

            if (!alisimulator->params->do_compression)
                start_pos = out->tellp();

            // for Windows only, the line break is \r\n instead of only \n
            #if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
            ++start_pos;
            #endif

            uint64_t output_line_length = seq_length_times_num_sites_per_state + 1 + alisimulator->max_length_taxa_name + (alisimulator->params->aln_output_format == IN_FASTA ? 1 : 0);
        
            // for Windows only, the line break is \r\n instead of only \n
            #if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
            output_line_length += alisimulator->params->aln_output_format == IN_FASTA ? 2 : 1;
            #endif
        
            // initialize state_mapping (mapping from state to characters)
            vector<string> state_mapping;
            AliSimulator::initializeStateMapping(alisimulator->num_sites_per_state, aln, state_mapping);
        
            // write sequences at tips to output file from a tmp_data and genome trees => a special case: with Indels without FunDi/ASC/Partitions
            bool write_sequences_from_tmp_data = alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0 && alisimulator->params->alisim_fundi_taxon_set.size() == 0 && !(alisimulator->tree->getModelFactory() && alisimulator->tree->getModelFactory()->getASC() != ASC_NONE) && !alisimulator->tree->isSuperTree();
            if (write_sequences_from_tmp_data)
                writeSeqsFromTmpDataAndGenomeTreesIndels(alisimulator, sequence_length, *out, *out_indels, write_indels_output, state_mapping, alisimulator->params->aln_output_format, alisimulator->max_length_taxa_name);
        
            int num_threads = 1;
            #ifdef _OPENMP
            #pragma omp parallel
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
            #endif
                // browsing all sequences, converting each sequence & caching & writing output string to file
                writeASequenceToFile(aln, sequence_length, num_threads, alisimulator->params->keep_seq_order, start_pos, output_line_length, *out, *out_indels, write_indels_output, state_mapping, alisimulator->params->aln_output_format, alisimulator->max_length_taxa_name, write_sequences_from_tmp_data, alisimulator->tree->root, alisimulator->tree->root);
            #ifdef _OPENMP
            }
            #endif

            // close the output file for Indels
            if (write_indels_output)
            {
                // close the file
                if (alisimulator->params->do_compression)
                    ((ogzstream*)out_indels)->close();
                else
                    ((ofstream*)out_indels)->close();
                delete out_indels;
            }
            
            // close the file
            if (alisimulator->params->do_compression)
                ((ogzstream*)out)->close();
            else
                ((ofstream*)out)->close();
            delete out;
        
            // show actual output sequence length in simulations with Indels
//            if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0)
//                cout << "Output sequence length of " << file_path << ": " << convertIntToString(sequence_length) << endl;
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, file_path);
        }
}

/**
*  merge and write all sequences to output files
*/
void mergeAndWriteSequencesToFiles(string file_path, AliSimulator *alisimulator, vector<SeqType>& seqtypes, vector<std::string>& aln_names, std::ios_base::openmode open_mode){
    // in case with partitions -> merge & write sequences to a single/multiple files
    if (alisimulator->tree->isSuperTree())
    {
        PhyloSuperTree *super_tree = ((PhyloSuperTree *)alisimulator->tree);
        int total_expected_num_states = super_tree->getAlnNSite();
        
        // merge phylotrees using the same alignment file and the same sequence_type and the same num_states (morph) with each other -> write sequences to the corresponding output files
        for (int i = 0; i < super_tree->size(); i++)
        {
            // ignore the current phylotree if it has already merged
            bool already_merged = false;
            for (int j = 0; j < i; j++)
                if (!super_tree->at(i)->aln->aln_file.compare(super_tree->at(j)->aln->aln_file)
                    && super_tree->at(i)->aln->seq_type == super_tree->at(j)->aln->seq_type && super_tree->at(i)->aln->num_states == super_tree->at(j)->aln->num_states)
                {
                    already_merged = true;
                    break;
                }
            
            // merge the sequences of the current phylotree and other phylotrees which use the same alignment file to the super_tree
            if (!already_merged)
            {
                // dummy variables
                auto start = getRealTime();
                string partition_list = "";
                int partition_count = 0;
                
                // the maximum index of sites in all partitions which use the same alignment file.
                int max_site_index = 0;
                int j;
                
                // clear out all sequences in the current super_tree
                clearoutSequencesSuperTree(super_tree->root, super_tree->root);
                #ifdef _OPENMP
                #pragma omp parallel for shared(partition_list, partition_count, max_site_index)
                #endif
                for (j = i; j < super_tree->size(); j++)
                {
                    IQTree *current_tree = (IQTree*) super_tree->at(j);
                    if (!super_tree->at(i)->aln->aln_file.compare(current_tree->aln->aln_file) && super_tree->at(i)->aln->seq_type == super_tree->at(j)->aln->seq_type && super_tree->at(i)->aln->num_states == super_tree->at(j)->aln->num_states)
                    {
                        // extract site_ids of the partition
                        const char* info_spec = ((SuperAlignment*) super_tree->aln)->partitions[j]->CharSet::position_spec.c_str();
                        IntVector site_ids;
                        extractSiteID(current_tree->aln, info_spec, site_ids, false, total_expected_num_states);
                        
                        // copy alignment from the current tree to the super_tree
                        copySequencesToSuperTree(site_ids, total_expected_num_states, current_tree, current_tree->aln->STATE_UNKNOWN, super_tree->root, super_tree->root);
                        
                        // determine the max_site_index for the current partition
                        for (int site_index:site_ids)
                            if (current_tree->max_site_id_mapping <site_index)
                                current_tree->max_site_id_mapping = site_index;
                        #ifdef _OPENMP
                        #pragma omp critical
                        #endif
                        {
                            // update the overall max_site_index
                            if (max_site_index < current_tree->max_site_id_mapping)
                                max_site_index = current_tree->max_site_id_mapping;
                            
                            // update partition_list
                            partition_list = partition_list + "_" + super_tree->at(j)->aln->name;
                            partition_count++;
                        }
                    }
                }
                
                // insert redundant sites (inserted sites due to Indels) to the sequences
                if (super_tree->params->alisim_insertion_ratio + super_tree->params->alisim_deletion_ratio > 0)
                {
                    cout << endl << "Actual sequence length of partitions (due to Indels):" << endl;
                    vector<short int> site_index_step_mapping(max_site_index+1, 0);
                    for (j = i; j < super_tree->size(); j++)
                    {
                        IQTree *current_tree = (IQTree*) super_tree->at(j);
                        if (!super_tree->at(i)->aln->aln_file.compare(current_tree->aln->aln_file) && super_tree->at(i)->aln->seq_type == super_tree->at(j)->aln->seq_type && super_tree->at(i)->aln->num_states == super_tree->at(j)->aln->num_states)
                        {
                            // determine the expected number of sites, the real number of sites, then compute the number of inserted_sites
                            int expected_num_sites = current_tree->aln->getNSite();
                            bool stop = false;
                            int real_sequence_length = 0;
                            determineSequenceLength(current_tree->root, current_tree->root, stop, real_sequence_length);
                            int num_inserted_sites = real_sequence_length - expected_num_sites;
                            
                            // output actual sequence length of the current partition
                            string partition_name = ((SuperAlignment*) super_tree->aln)->partitions[j]->CharSet::name;
                            cout << partition_name << ": " << (current_tree->aln->seq_type == SEQ_CODON ? (real_sequence_length * 3) : real_sequence_length) << endl;
                            
                            if (num_inserted_sites > 0)
                            {
                                // insert Indels sites
                                insertIndelSites(current_tree->max_site_id_mapping + site_index_step_mapping[current_tree->max_site_id_mapping]+1, expected_num_sites, num_inserted_sites, current_tree, super_tree->root, super_tree->root);
                                
                                // update site_index_step_mapping
                                for (int k = current_tree->max_site_id_mapping + 1; k < site_index_step_mapping.size(); k++)
                                    site_index_step_mapping[k] += num_inserted_sites;
                                
                                // update the max_site_index
                                max_site_index += num_inserted_sites;
                            }
                        }
                    }
                    
                }
                
                // show time spent on merging partitions
                auto end = getRealTime();
                cout<<" - Time spent on merging partitions: "<<end-start<<endl;
                
                // initialize output file name
                if (partition_count == super_tree->size())
                    partition_list = "";
                else
                {
                    // instead of using partition_list, use data_type for output file names
                    partition_list = "_"+super_tree->at(i)->aln->sequence_type;
                    if (super_tree->at(i)->aln->seq_type == SEQ_MORPH)
                        partition_list += convertIntToString(super_tree->at(i)->aln->getMaxNumStates());
                }
                
                
                //  get the num_nodes
                int num_nodes = super_tree->leafNum;
                if (alisimulator->params->alisim_write_internal_sequences)
                    num_nodes = super_tree->nodeNum;
                // don't count the fake root
                num_nodes -= ((super_tree->root->isLeaf() && super_tree->root->name == ROOT_NAME)?1:0);
                
                // record the seqtype and alignment names, which will be used later to convert the simulated alignment into Maple format
                seqtypes.push_back(super_tree->at(i)->aln->seq_type);
                aln_names.push_back(file_path + partition_list);
                
                // write the merged sequences to the output file for the current cluster of partitions
                writeSequencesToFile(file_path + partition_list, super_tree->at(i)->aln, max_site_index+1, num_nodes, alisimulator, open_mode);
            }
        }
    }
    // other cases (without partitions), just write sequences to a single file
    else
    {
        int sequence_length = round(alisimulator->expected_num_sites * alisimulator->inverse_length_ratio);
        
        // determine the real sequence_length if Indels is used
        if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0)
            sequence_length = alisimulator->seq_length_indels;
        
        //  get the num_nodes
        int num_nodes = alisimulator->tree->leafNum;
        if (alisimulator->params->alisim_write_internal_sequences)
            num_nodes = alisimulator->tree->nodeNum;
        // don't count the fake root
        num_nodes -= ((alisimulator->tree->root->isLeaf() && alisimulator->tree->root->name == ROOT_NAME)?1:0);
        
        writeSequencesToFile(file_path, alisimulator->tree->aln, sequence_length, num_nodes, alisimulator, open_mode);
    }
}

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(Alignment *aln, int sequence_length, int num_threads, bool keep_seq_order, uint64_t start_pos, uint64_t output_line_length, ostream &out, ostream &out_indels, bool write_indels_output, vector<string> &state_mapping, InputType output_format, int max_length_taxa_name, bool write_sequences_from_tmp_data, Node *node, Node *dad)
{
    // if write_sequences_from_tmp_data and this node is a leaf -> skip this node as its sequence was already written to the output file
    if ((!(node->isLeaf() && write_sequences_from_tmp_data))
        &&((node->isLeaf() && node->name!=ROOT_NAME) || (Params::getInstance().alisim_write_internal_sequences && Params::getInstance().alisim_insertion_ratio + Params::getInstance().alisim_deletion_ratio > 0))) {
        #ifdef _OPENMP
        #pragma omp task firstprivate(node, num_threads, keep_seq_order) shared(out, out_indels, state_mapping)
        #endif
        {
            int num_sites_per_state = aln->seq_type == SEQ_CODON?3:1;
            // initialize the output sequence with all gaps (to handle the cases with missing taxa in partitions)
            string pre_output = AliSimulator::exportPreOutputString(node, output_format, max_length_taxa_name);
            string output(aln->seq_type == SEQ_CODON ? (3 * sequence_length) : sequence_length, '-');
            uint64_t output_pos = start_pos + node->id * output_line_length;
            
            // convert non-empty sequence
            AliSimulator::convertNumericalStatesIntoReadableCharacters(node->sequence->sequence_chunks[0], output, sequence_length, num_sites_per_state, state_mapping);
            
            // preparing output (without gaps) for indels
            string output_indels = "";
            if (write_indels_output)
            {
                // add node's name
                string node_name = node->name;
                // write node's id if node's name is empty
                if (node_name.length() == 0) node_name = convertIntToString(node->id);
                
                output_indels = ">" + node_name + "\n" + output;
                output_indels.erase(remove(output_indels.begin(), output_indels.end(), '-'), output_indels.end());
            }
            
            // concat pre_output & output
            output = pre_output + output;
            
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                // jump to the correct position if user want to keep sequence order
                if (num_threads != 1 && keep_seq_order)
                    out.seekp(output_pos);
                
                // write output to file
                out << output << "\n";
                
                // write aln without gaps for Indels
                if (write_indels_output)
                    out_indels << output_indels << "\n";
            }
        }
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFile(aln, sequence_length, num_threads, keep_seq_order, start_pos, output_line_length, out, out_indels, write_indels_output, state_mapping, output_format, max_length_taxa_name, write_sequences_from_tmp_data, (*it)->node, node);
    }
}

/**
*  clear out all sequences in the super_tree
*
*/
void clearoutSequencesSuperTree(Node *node, Node *dad){
    #ifdef _OPENMP
    #pragma omp task firstprivate(node)
    #endif
    if (node->isLeaf())
        node->sequence->sequence_chunks[0].clear();

    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        clearoutSequencesSuperTree((*it)->node, node);
     }
}

map<string,string> loadInputMSA(AliSimulator *alisimulator)
{
    map<string,string> input_msa;
    // don't load Input MSA if either partitions or ASC model is being used
    if (alisimulator->params->alisim_inference_mode &&
        ((alisimulator->tree->getModelFactory() && alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
        || alisimulator->tree->isSuperTree()
        || alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0))
    {
        outWarning("AliSim will not copy gaps from the input alignment into the output alignments in simulations with Indels/Partitions/+ASC models.");
        return input_msa;
    }
    
    // only load Input MSA if the user has specified an alignment file and wants to copy gaps from the input MSA.
    if (alisimulator->params->aln_file && !alisimulator->params->alisim_no_copy_gaps)
    {
        // initialize dummy variables
        vector<string> seq_names;
        vector<string> sequences;
        int nseq = 0, nsite = 0;
        Alignment *aln = alisimulator->tree->aln;
        char *sequence_type_char = strcpy(new char[aln->sequence_type.length() + 1], aln->sequence_type.c_str());
        
        // read sequences from the input file
        Alignment *tmp_aln = new Alignment();
        tmp_aln->extractSequences(alisimulator->params->aln_file, sequence_type_char, sequences, nseq, nsite);
        seq_names = tmp_aln->getSeqNames();

        // show a warning if the length of input alignment is unequal to that of simulated sequence
        int sequence_length = round(alisimulator->expected_num_sites * alisimulator->inverse_length_ratio);
        if (sequences.size() > 0 && sequences[0].length() != (alisimulator->num_sites_per_state == 1 ? sequence_length : (sequence_length * alisimulator->num_sites_per_state)))
            outWarning("The sequence length of the input alignment is unequal to that of that simulated sequences. Thus, only gaps in the first MIN(input_sequence_length, simulated_sequence_length) sites are copied.");
        
        // return InputMSA;
        for (int i = 0; i < seq_names.size(); i++)
            input_msa.insert(pair<string,string>(seq_names[i], sequences[i]));
        return input_msa;
    }
    return input_msa;
}

/**
*  only unroot tree and stop if the user wants to do so
*
*/
void unrootTree(AliSimulator *alisimulator)
{
    // initialize output_filepath
    string output_filepath(alisimulator->params->user_file);
    output_filepath = output_filepath.substr(0, output_filepath.find_last_of(".") + 1);
    output_filepath = output_filepath + "unrooted.treefile";
    
    if (alisimulator->tree->rooted)
    {
        cout<<"Unrooting the input tree"<<endl;
        alisimulator->tree->PhyloTree::forceConvertingToUnrooted();
    
        cout<<"Outputting the unrooted tree to "+output_filepath<<endl;
    }
    else
        outWarning("The input tree is unrooted, thus, not needing to unroot it.");
    
    // write output to file
    ofstream *out = new ofstream(output_filepath.c_str());
    alisimulator->tree->PhyloTree::printTree(*out);
    ((ofstream*)out)->close();
    delete out;
}

/**
*  determine real sequence length (for Indels)
*
*/
void determineSequenceLength(Node *node, Node *dad, bool &stop, int &sequence_length)
{
    // check to stop
    if (stop)
        return;
    
    // determine the real sequence_length
    if (node->name!=ROOT_NAME && node->sequence->sequence_chunks[0].size() > 0)
    {
        sequence_length = node->sequence->sequence_chunks[0].size();
        stop = true;
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        determineSequenceLength((*it)->node, node, stop, sequence_length);
    }
}


/**
*  insert redundant sites (inserted sites due to Indels) to the sequences of the super tree
*/
void insertIndelSites(int position, int starting_index, int num_inserted_sites, IQTree *current_tree, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // find the corresponding node in the current_tree
        Node *current_node = current_tree->findLeafName(node->name);

        // if current_node is found, inserting sites normally
        if (current_node)
            node->sequence->sequence_chunks[0].insert(node->sequence->sequence_chunks[0].begin()+position, current_node->sequence->sequence_chunks[0].begin()+starting_index, current_node->sequence->sequence_chunks[0].end());
        // otherwise, insert gaps
        else
            node->sequence->sequence_chunks[0].insert(node->sequence->sequence_chunks[0].begin()+position, num_inserted_sites, current_tree->aln->STATE_UNKNOWN);
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        insertIndelSites(position, starting_index, num_inserted_sites, current_tree, (*it)->node, node);
    }
}

/**
*  write sequences to output file from a tmp_data and genome trees => a special case: with Indels without FunDi/ASC/Partitions
*/
void writeSeqsFromTmpDataAndGenomeTreesIndels(AliSimulator* alisimulator, int sequence_length, ostream &out, ostream &out_indels, bool write_indels_output, vector<string> &state_mapping, InputType output_format, int max_length_taxa_name)
{
    // read tmp_data line by line
    igzstream in;
    int line_num = 1;
    string line;
    in.open((Params::getInstance().alisim_output_filename + "_" + Params::getInstance().tmp_data_filename + "_" + convertIntToString(MPIHelper::getInstance().getProcessID())).c_str());
    
    // dummy variables
    GenomeTree* genome_tree = NULL;
    Insertion* previous_insertion = NULL;
    int num_sites_per_state = alisimulator->tree->aln->seq_type == SEQ_CODON ? 3 : 1;
    int seq_length_times_num_sites_per_state = alisimulator->tree->aln->seq_type == SEQ_CODON ? (sequence_length * 3) : sequence_length;
    int rebuild_indel_his_step = alisimulator->params->rebuild_indel_history_param * alisimulator->tree->leafNum;
    int rebuild_indel_his_thresh = rebuild_indel_his_step;

    for (; !in.eof(); line_num++)
    {
        safeGetline(in, line);
        line = line.substr(0, line.find_first_of("\n\r"));
        if (line == "" || (alisimulator->params->aln_output_format == IN_PHYLIP && line_num == 1)) continue;
        
        // extract seq_name
        int index_of_first_at = line.find_first_of("@");
        int index_of_second_at = line.find_first_of("@", index_of_first_at + 1);
        string seq_name = line.substr(0, index_of_first_at);
        
        // retrieve Node from the seq_name
        Node* node = alisimulator->map_seqname_node[seq_name];
        if (!node)
            outError("Oops! Couldn't find the node with name " + seq_name+" . There is something wrong!");
        
        // extract the length of the original sequence
        int seq_length_ori = convert_int(line.substr(index_of_first_at + 1, index_of_second_at - index_of_first_at - 1).c_str());
        
        // extract original sequences
        vector<short int> seq_ori(seq_length_ori, 0);
        string internal_states = line.substr(index_of_second_at + 1, line.length() - index_of_second_at - 1);
        istringstream seq_in(internal_states);
        for (int i = 0; i < seq_length_ori; i++)
            seq_in >> seq_ori[i];
        
        // initialize the output sequence with all gaps (to handle the cases with missing taxa in partitions)
        string pre_output = AliSimulator::exportPreOutputString(node, output_format, max_length_taxa_name);
        string output(seq_length_times_num_sites_per_state, '-');
        
        // build a new genome tree from the list of insertions if the genome tree has not been initialized (~NULL)
        if (!genome_tree)
        {
            genome_tree = new GenomeTree();
            genome_tree->buildGenomeTree(node->sequence->insertion_pos, seq_length_ori, true);
        }
        // otherwise, update the tree by accepted gaps (inserted by previous insertions) as normal characters
        else
        {
            // if it is not the last tip -> rebuild/update the genome tree
            if (node->sequence->insertion_pos->next)
            {
                // rebuild the indel his if the number of tips (line_num) >= current threshold
                if (line_num >= rebuild_indel_his_thresh)
                {
                    // detach the insertion and genome nodes
                    for (Insertion* insertion = node->sequence->insertion_pos; insertion; )
                    {
                        // detach insertion and genome_nodes
                        insertion->genome_nodes.clear();
                        
                        // move to the next insertion
                        insertion = insertion->next;
                    }
                    
                    // delete and rebuild genome tree
                    delete genome_tree;
                    genome_tree = new GenomeTree();
                    genome_tree->buildGenomeTree(node->sequence->insertion_pos, seq_length_ori, true);
                    
                    // update the next threshold to rebuild the indel his
                    rebuild_indel_his_thresh += rebuild_indel_his_step;
                }
                // otherwise, just update indel his
                else
                    genome_tree->updateGenomeTree(previous_insertion, node->sequence->insertion_pos);
            }
            // otherwise, it is the last tip -> the current sequence is already the latest sequence since there no more insertion occurs
            else
            {
                delete genome_tree;
                genome_tree = new GenomeTree(seq_length_ori);
            }
        }
        
        // keep track of previous insertion
        previous_insertion = node->sequence->insertion_pos;
        
        // delete the insertion_pos of this node as we updated its sequence.
        node->sequence->insertion_pos = NULL;
        
        // export sequence of a leaf node from original sequence and genome_tree if using Indels
        genome_tree->exportReadableCharacters(seq_ori, num_sites_per_state, state_mapping, output);
        
        // preparing output (without gaps) for indels
        string output_indels = "";
        if (write_indels_output)
        {
            // add node's name
            string node_name = node->name;
            // write node's id if node's name is empty
            if (node_name.length() == 0) node_name = convertIntToString(node->id);
            
            output_indels = ">" + node_name + "\n" + output;
            output_indels.erase(remove(output_indels.begin(), output_indels.end(), '-'), output_indels.end());
        }
        
        // concat pre_output & output
        output = pre_output + output;
        
        // write output to file
        out << output << "\n";
        
        // write aln without gaps for Indels
        if (write_indels_output)
            out_indels << output_indels << "\n";
    }
    
    // delete the genome tree
    delete genome_tree;
    
    // close the tmp_data file
    in.close();
}

void outputTreeWithInternalNames(AliSimulator* alisimulator)
{
    // don't need to handle supertree here (at current stage) as AliSim doesn't support outputting internal sequences in simulations with partitions
    // set names for internal nodes
    updateInternalNodeName(alisimulator->tree->root);
    
    // output the treefile
    string output_filepath = alisimulator->params->alisim_output_filename + ".full.treefile";
    std::ofstream treefile(output_filepath, std::ios::out);
    alisimulator->tree->printTree(treefile);
    treefile.close();
    
    // show message
    std::cout << "A tree (with internal node names) has been outputted to " << output_filepath << std::endl;
}

void updateInternalNodeName(Node *node, Node *dad)
{
    // if node is an internal and has an empty name -> set its name as its id
    if (!node->isLeaf() && node->name == "")
        node->name = convertIntToString(node->id);
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        updateInternalNodeName((*it)->node, node);
    }
}
