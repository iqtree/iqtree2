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
    }
    
    // execute AliSim Simulation
    executeSimulation(params, tree);
    
    // aln and tree are deleted in distructor of AliSimSimulator
    MPIHelper::getInstance().barrier();
    auto end = getRealTime();
    cout << "Simulation time: " << fixed << end-start << "s" << endl;
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
    if (params.run_mode == CALC_DIST) {
        if (params.tree_gen == CIRCULAR_SPLIT_GRAPH) {
            cout << "Calculating distance matrix..." << endl;
            sg.calcDistance(params.dist_file);
            cout << "Distances printed to " << params.dist_file << endl;
        }// else {
            //mtree.calcDist(params.dist_file);
        //}
    }

}

/**
*  execute AliSim Simulation
*/
void executeSimulation(Params params, IQTree *&tree)
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
    
    // load input MSA if any
    map<string,string> input_msa = loadInputMSA(alisimulator);
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    generateMultipleAlignmentsFromSingleTree(alisimulator, input_msa);
    
    // delete alisimulator
    delete alisimulator;
    
    cout << "[Alignment Simulator] Done"<<"\n";
}

/**
*  show all input parameters for AliSim
*/
void showParameters(Params params, bool is_partition_model)
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
    if (params.alisim_ancestral_sequence_name.length() > 0)
        cout << " - Ancestral sequence position: " << params.alisim_dataset_num <<"\n";
}

/**
*  retrieve the ancestral sequence for the root node from an input file
*/
void retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator, vector<short int> &sequence)
{
    // get variables
    char *aln_filepath = super_alisimulator->params->alisim_ancestral_sequence_aln_filepath;
    string sequence_name = super_alisimulator->params->alisim_ancestral_sequence_name;
    
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
    
    // delete aln
    delete aln;
    
    // overwrite the output sequence_length
    if (super_alisimulator->params->alisim_sequence_length != nsite)
    {
        super_alisimulator->params->alisim_sequence_length = nsite;
        outWarning("Sequence length is now set equally to the length of ancestral sequence.");
        super_alisimulator->refreshExpectedNumSites();
    }
    
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

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator, map<string,string> input_msa)
{
    // Load ancestral sequence from the input file if user has specified it
    vector<short int> ancestral_sequence;
    if (super_alisimulator->params->alisim_ancestral_sequence_name.length() > 0)
        retrieveAncestralSequenceFromInputFile(super_alisimulator, ancestral_sequence);
    
    // terminate if users employ more MPI processes than the number of alignments
    if (MPIHelper::getInstance().getNumProcesses() > super_alisimulator->params->alisim_dataset_num)
        outError("You are employing more MPI processes (" + convertIntToString(MPIHelper::getInstance().getNumProcesses()) + ") than the number of alignments (" + convertIntToString(super_alisimulator->params->alisim_dataset_num) + "). Please reduce the number of MPI processes to save the computational resources and try again!");
    
    // reset number of OpenMP threads to 1 in simulations with Indels
    if (super_alisimulator->params->num_threads > 1 && super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0)
        outError("OpenMP has not yet been supported in simulations with Indels. Please use a single thread for this simulation.");
    
    // do not support compression when outputting multiple data sets into a same file
    if (Params::getInstance().do_compression && (Params::getInstance().alisim_single_output || Params::getInstance().keep_seq_order))
    {
        outWarning("Compression is not supported when either outputting multiple alignments into a single output file or keeping the order of output sequences. AliSim will output file in normal format.");

        Params::getInstance().do_compression = false;
        super_alisimulator->params->do_compression = false;
    }
    
    // cannot skip concatenating sequence chunks from intermediate files in simulations with FunDi, Partitions, or +ASC models
    if (Params::getInstance().num_threads > 1 && Params::getInstance().no_merge)
    {
        // ignore --no-merge option if using AliSim-OpenMP-IM
        if (Params::getInstance().alisim_openmp_alg == IM)
            outWarning("Ignore --no-merge option as it is only appliable for AliSim-OpenMP-EM algorithm.");
        // otherwise, if using AliSim-OpenMP-EM -> show a warning in cases that we cannot skip merging
        else if (super_alisimulator->tree->isSuperTree()
                  || super_alisimulator->params->alisim_fundi_taxon_set.size() > 0
                  || (super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE))
        {
            outWarning("Cannot skip merging sequence chunks in simulations with FunDi, Partitions, or +ASC models. AliSim will concatenate sequence chunks from intermediate files into a single output file.");
            
            Params::getInstance().no_merge = false;
            super_alisimulator->params->no_merge = false;
        }
    }
    
    // show a warning if the user wants to write internal sequences in not-supported cases
    if (super_alisimulator->params->alisim_write_internal_sequences
        &&((super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
           || super_alisimulator->tree->isSuperTree() || (super_alisimulator->params->alisim_fundi_taxon_set.size() > 0 && Params::getInstance().num_threads > 1)))
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
        
    // ignore --single-output in version with MPI
#ifdef _IQTREE_MPI
    if (super_alisimulator->params->alisim_single_output)
    {
        outWarning("Ignore --single-output option since it is not supported in IQ-TREE version with MPI. Alignments will be outputted in separated files.");
        super_alisimulator->params->alisim_single_output = false;
    }
#endif
    
    // iteratively generate multiple datasets for each tree
    for (int i = 0; i < super_alisimulator->params->alisim_dataset_num; i++)
    {
        // parallelize over MPI ranks statically
        int proc_ID = MPIHelper::getInstance().getProcessID();
        int nprocs  = MPIHelper::getInstance().getNumProcesses();
        if (i%nprocs != proc_ID) continue; 
        
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
                generatePartitionAlignmentFromSingleSimulator(partition_simulator, ancestral_sequence_current_tree, input_msa);
                
                // update new genome at tips from the original genome and the genome tree
                // skip updating if using +ASC or Fundi model as they must be already updated
                if (super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0 && !(partition_simulator->tree->getModelFactory() && partition_simulator->tree->getModelFactory()->getASC() != ASC_NONE) && (partition_simulator->params->alisim_fundi_taxon_set.size() == 0))
                    partition_simulator->updateNewGenomeIndels(partition_simulator->seq_length_indels);
            }
        }
        else
        {
            // check whether we could write the output to file immediately after simulating it
            if (super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() == ASC_NONE && super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio == 0)
                generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence, input_msa, output_filepath, open_mode);
            // otherwise, writing output to file after completing the simulation
            else
                generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence, input_msa);
        }
        
        // merge & write alignments to files if they have not yet been written
        if ((super_alisimulator->tree->getModelFactory() && super_alisimulator->tree->getModelFactory()->getASC() != ASC_NONE)
            || super_alisimulator->tree->isSuperTree()
            || super_alisimulator->params->alisim_insertion_ratio + super_alisimulator->params->alisim_deletion_ratio > 0)
            mergeAndWriteSequencesToFiles(output_filepath, super_alisimulator, open_mode);
        
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
        
        // delete output alignments (for testing only)
        if (super_alisimulator->params->delete_output)
        {
            string output_filename = output_filepath;
            
            for (int thread_id = 0; thread_id < super_alisimulator->params->num_threads; thread_id++)
            {
                if (super_alisimulator->params->num_threads > 1 && super_alisimulator->params->alisim_openmp_alg == EM && super_alisimulator->params->no_merge)
                    output_filename = output_filepath + "_" + convertIntToString(thread_id + 1);
                
                // add file extension
                if (super_alisimulator->params->aln_output_format == IN_PHYLIP)
                    output_filename += ".phy";
                else
                    output_filename += ".fa";
                
                // delete the output file
                remove((output_filename).c_str());
                
                // stop deleting if output file was merging
                if (super_alisimulator->params->alisim_openmp_alg == IM || !super_alisimulator->params->no_merge)
                    break;
            }
        }
    }
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
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *&alisimulator, vector<short int> &ancestral_sequence, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode)
{
    // show an error if continuous gamma is used in inference mode.
    if (alisimulator->params->alisim_inference_mode && alisimulator->tree->getModelFactory() && alisimulator->tree->getModelFactory()->is_continuous_gamma)
        outError("Unfortunately, IQ-Tree has not yet supported Continuous Gamma in phylogeny inference. Therefore, users can only use Continuous Gamma in AliSim without Inference Mode.");
    
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    bool is_mixture_model = alisimulator->tree->getModel()->isMixture();
    
    // case 1: without rate heterogeneity or mixture model -> using the current alisimulator (don't need to re-initialize it)
    
    // case 2: with rate heterogeneity or mixture model
    if ((!rate_name.empty()) || is_mixture_model)
    {
        // if user specifies +I without invariant_rate -> set it to 0
        if (rate_name.find("+I") != std::string::npos && isnan(invariant_proportion)) {
            alisimulator->tree->getRate()->setPInvar(0);
            outWarning("Invariant rate is now set to Zero since it has not been specified");
        }
        
        // case 2.3: with only invariant sites (without gamma/freerate model/mixture models)
        if (!rate_name.compare("+I") && !is_mixture_model)
        {
            alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
        else
        {
            // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
            if (invariant_proportion > 0)
            {
                alisimulator = new AliSimulatorHeterogeneityInvar(alisimulator, invariant_proportion);
            }
            // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
            else
            {
                alisimulator = new AliSimulatorHeterogeneity(alisimulator);
            }
        }
    }
    
    alisimulator->generatePartitionAlignment(ancestral_sequence, input_msa, output_filepath, open_mode);
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
            if (alisimulator->params->aln_output_format != IN_FASTA)
                file_path = file_path + ".phy";
            else
                file_path = file_path + ".fa";
            ostream *out;
            if (alisimulator->params->do_compression)
                out = new ogzstream(file_path.c_str(), open_mode);
            else
                out = new ofstream(file_path.c_str(), open_mode);
            out->exceptions(ios::failbit | ios::badbit);

            // write the first line <#taxa> <length_of_sequence> (for PHYLIP output format)
            int seq_length_times_num_sites_per_state = (aln->seq_type == SEQ_CODON ? (sequence_length * 3) : sequence_length);
            string first_line = "";
            if (alisimulator->params->aln_output_format != IN_FASTA)
            {
                first_line = convertIntToString(num_leaves) + " " + convertIntToString(seq_length_times_num_sites_per_state) + "\n";
                *out << first_line;
            }
            
            // get the position to write output
            uint64_t start_pos = first_line.length();
            if (!alisimulator->params->do_compression)
                start_pos = out->tellp();
            uint64_t output_line_length = seq_length_times_num_sites_per_state + 1 + alisimulator->max_length_taxa_name + (alisimulator->params->aln_output_format == IN_FASTA ? 1 : 0);
        
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
        
            // show the output file name
            if (!(MPIHelper::getInstance().getNumProcesses() > 1 && alisimulator->params->alisim_dataset_num > 1))
                cout << "An alignment has just been exported to "<<file_path<<endl;
        
            // show actual output sequence length in simulations with Indels
            if (alisimulator->params->alisim_insertion_ratio + alisimulator->params->alisim_deletion_ratio > 0)
                cout << "Output sequence length: " << convertIntToString(sequence_length) << endl;
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, file_path);
        }
}

/**
*  merge and write all sequences to output files
*/
void mergeAndWriteSequencesToFiles(string file_path, AliSimulator *alisimulator, std::ios_base::openmode open_mode){
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
                
                
                //  get the num_leaves
                int num_leaves = super_tree->leafNum - ((super_tree->root->isLeaf() && super_tree->root->name == ROOT_NAME)?1:0);
                
                // write the merged sequences to the output file for the current cluster of partitions
                writeSequencesToFile(file_path + partition_list, super_tree->at(i)->aln, max_site_index+1, num_leaves, alisimulator, open_mode);
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
        
        //  get the num_leaves
        int num_leaves = alisimulator->tree->leafNum - ((alisimulator->tree->root->isLeaf() && alisimulator->tree->root->name == ROOT_NAME)?1:0);
        writeSequencesToFile(file_path, alisimulator->tree->aln, sequence_length, num_leaves, alisimulator, open_mode);
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
        #pragma omp task firstprivate(node) shared(out, out_indels, state_mapping)
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
                if (num_threads > 1 && keep_seq_order)
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
