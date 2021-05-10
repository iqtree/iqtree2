/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"

void runAliSim(Params &params, Checkpoint *checkpoint)
{
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
    
    // generate a random tree if neccessary
    if (params.tree_gen != NONE)
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
    
    // inferring input parameters if inference mode is active
    if (inference_mode)
    {
        inferInputParameters(params, checkpoint, tree, aln);
    }
    
    // execute AliSim Simulation
    executeSimulation(params, tree, inference_mode);
    
    // aln and tree are deleted in distructor of AliSimSimulator
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
        params.user_file = new char[strlen(pre_fix) + 1];
        strcpy(params.user_file, pre_fix);
        if (params.partition_file && params.partition_type == BRLEN_OPTIMIZE)
            strcat(params.user_file,".parttrees");
        else
            strcat(params.user_file,".treefile");
    }
    // reload tree from tree_file
    else
    {
        // do not reload the tree from file in case with heterotachy
        if (!tree->getRate()->isHeterotachy())
        {
            bool is_rooted = false;
            tree->readTree(params.user_file, is_rooted);
        }
        
        // handle the case to reload super tree
        if (tree->isSuperTree())
        for (int i = 0; i < ((PhyloSuperTree*) tree)->size(); i++)
        {
            // -Q (params->partition_type == BRLEN_OPTIMIZE) -> tree_line_index = i; otherwise (-p, -q), tree_line_index = 0 (only a tree)
            int tree_line_index = params.partition_type == BRLEN_OPTIMIZE?i:0;
            
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
    }
    
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
            int num_sites_per_state = tree->aln->seq_type == SEQ_CODON?3:1;
            params.alisim_sequence_length = tree->aln->getNSite() * num_sites_per_state;
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
        int num_sites_per_state = aln->seq_type == SEQ_CODON?3:1;
        total_length += aln->getNSite() * num_sites_per_state;
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
            if (!overwriteFile(params.user_file)) return;
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
            if (!overwriteFile(params.user_file)) return;
            sg.generateCircular(params);
        } else if (params.tree_gen == TAXA_SET) {
            sg.init(params);
            cout << "Generating random taxa set of size " << params.sub_size <<
                " overlap " << params.overlap << " with " << params.repeated_time << " times..." << endl;
            if (!overwriteFile(params.pdtaxa_file)) return;
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
void executeSimulation(Params params, IQTree *&tree, bool inference_mode)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator;
    if (tree && inference_mode)
        alisimulator = new AliSimulator(&params, tree);
    else
        alisimulator = new AliSimulator(&params);
    
    // show parameters
    showParameters(params, alisimulator->tree->isSuperTree());
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    generateMultipleAlignmentsFromSingleTree(alisimulator, inference_mode);
    
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
IntVector retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator)
{
    // get variables
    IntVector sequence;
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
    int num_sites_per_state = src_tree->aln->seq_type == SEQ_CODON?3:1;
    int sequence_length = super_alisimulator->params->alisim_sequence_length/num_sites_per_state;
    
    // make sure the length of the ancestral sequence must be equal to the total length of all partitions
    if (super_alisimulator->tree->isSuperTree() && sequence_length != super_alisimulator->tree->getAlnNSite())
        outError("The length of the ancestral sequence must be equal to the total length of all partitions");
    
    sequence.resize(sequence_length);
    ostringstream err_str;
    int num_error = 0;
    for (int i = 0; i < sequence_length; i++)
    {
        if (src_tree->aln->seq_type == SEQ_CODON)
        {
            int site_index = i*num_sites_per_state;
            sequence[i] = src_tree->aln->getCodonStateTypeFromSites(src_tree->aln->convertState(sequence_str[site_index], SEQ_DNA), src_tree->aln->convertState(sequence_str[site_index+1], SEQ_DNA), src_tree->aln->convertState(sequence_str[site_index+2], SEQ_DNA), sequence_name, site_index, err_str, num_error);
        }
        else
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
        
    return sequence;
}

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator, bool inference_mode)
{
    // Load ancestral sequence from the input file if user has specified it
    IntVector ancestral_sequence;
    if (super_alisimulator->params->alisim_ancestral_sequence_name.length() > 0)
        ancestral_sequence = retrieveAncestralSequenceFromInputFile(super_alisimulator);
    
    // iteratively generate multiple datasets for each tree
    for (int i = 0; i < super_alisimulator->params->alisim_dataset_num; i++)
    {
        // initialize output_filepath
        std::string output_filepath(super_alisimulator->params->user_file);
        output_filepath = output_filepath.substr(0, output_filepath.find_last_of("/\\") + 1);
        output_filepath = output_filepath
        +super_alisimulator->params->alisim_output_filename
        +"_"+convertIntToString(i);
        
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
                
                // extract site_ids of the partition
                string info_spec_str = ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec;
                // convert position_spec from "*" to "start-end"
                if (!info_spec_str.compare("*") && super_tree->at(partition_index)->aln->getNSite() > 0)
                {
                    int num_sites_per_state = super_tree->at(partition_index)->aln->seq_type == SEQ_CODON?3:1;
                    info_spec_str = "1-" + convertIntToString(super_tree->at(partition_index)->aln->getNSite()*num_sites_per_state);
                    ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec = info_spec_str;
                }
                
                const char* info_spec = info_spec_str.c_str();
                IntVector site_ids;
                current_tree->aln->extractSiteID(current_tree->aln, info_spec, site_ids, total_expected_num_states);

                // extract the ancestral sequence for the current partition from the full ancestral_sequence
                IntVector ancestral_sequence_current_tree;
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
                generatePartitionAlignmentFromSingleSimulator(partition_simulator, ancestral_sequence_current_tree);
            }
        }
        else
            generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence);
        
        // merge & write alignments to files
        mergeAndWriteSequencesToFiles(output_filepath, super_alisimulator, inference_mode);
    }
}

/**
    copy sequences of leaves from a partition tree to super_tree
*/
void copySequencesToSuperTree(IntVector site_ids, int expected_num_states_super_tree, IQTree *super_tree, Node *node, Node *dad){
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // find the corresponding node (super_node) in the super_tree
        Node *super_node = super_tree->findLeafName(node->name);

        // make sure super_node is found
        if (super_node)
        {
            ASSERT(site_ids.size() == node->sequence.size());

            // initialize sequence of the super_node
            if (super_node->sequence.size() != expected_num_states_super_tree)
                super_node->sequence.resize(expected_num_states_super_tree);

            // copy sites one by one from the current sequence to its position in the sequence of the super_node
            for (int i = 0; i < node->sequence.size(); i++)
                super_node->sequence[site_ids[i]] = node->sequence[i];
        }
    }
    
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        // browse 1-step deeper to the neighbor node
        copySequencesToSuperTree(site_ids, expected_num_states_super_tree, super_tree, (*it)->node, node);
    }
}

/**
*  generate a partition alignment from a single simulator
*/
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *alisimulator, IntVector ancestral_sequence)
{
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    bool is_mixture_model = alisimulator->tree->getModel()->isMixture();
    
    // case 1: without rate heterogeneity or mixture model -> using the current alisimulator (don't need to re-initialize it)
    
    // case 2: with rate heterogeneity or mixture model
    if ((!rate_name.empty()) || is_mixture_model)
    {
        // case 2.3: with only invariant sites (without gamma/freerate model/mixture models)
        if (!rate_name.compare("+I"))
        {
            alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
        else
        {
            // if user specifies +I without invariant_rate -> set it to 0
            if (rate_name.find("+I") != std::string::npos && isnan(invariant_proportion)) {
                alisimulator->tree->getRate()->setPInvar(0);
                outWarning("Invariant rate is now set to Zero since it has not been specified");
            }
            
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
    
    alisimulator->generatePartitionAlignment(ancestral_sequence);
}

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, Alignment *aln, int sequence_length, AliSimulator *alisimulator, bool inference_mode)
{
    try {
            // add ".phy" to the file_path
            file_path = file_path + ".phy";
            ofstream out;
            out.exceptions(ios::failbit | ios::badbit);
            out.open(file_path.c_str());

            // write the first line <#taxa> <length_of_sequence>
            int leaf_num = alisimulator->tree->leafNum - ((alisimulator->tree->root->isLeaf() && alisimulator->tree->root->name == ROOT_NAME)?1:0);
            int num_sites_per_state = aln->seq_type == SEQ_CODON?3:1;
            out <<leaf_num<<" "<<sequence_length*num_sites_per_state<< endl;

            // write senquences of leaf nodes to file with/without gaps copied from the input sequence
            if (inference_mode && !alisimulator->params->alisim_no_copy_gaps)
            {
                
                // otherwise, copying gaps from the input alignment to the output sequences
                // load input sequences (with gaps)
                vector<string> seq_names;
                vector<string> sequences;
                // read sequences from the input file
                Alignment *tmp_aln = new Alignment();
                int nseq = 0, nsite = 0;
                char *sequence_type_char = strcpy(new char[aln->sequence_type.length() + 1], aln->sequence_type.c_str());
                string aln_file_str = aln->aln_file.length()?aln->aln_file:alisimulator->params->aln_file;
                char *aln_file_char = strcpy(new char[aln_file_str.length() + 1], aln_file_str.c_str());
                tmp_aln->extractSequences(aln_file_char, sequence_type_char, sequences, nseq, nsite);
                seq_names = tmp_aln->getSeqNames();

                // show a warning if the length of input alignment is unequal to that of simulated sequence
                if (sequences.size() > 0 && sequences[0].length() != sequence_length)
                    outWarning("The sequence length of the input alignment is unequal to that of that simulated sequences. Thus, only gaps in the first MIN(input_sequence_length, simulated_sequence_length) sites are copied.");

                // write simulated sequence with the gaps copied from the input sequence
                writeASequenceToFileWithGaps(aln, sequence_length, seq_names, sequences, out, alisimulator->tree->root, alisimulator->tree->root);
            }
            else
            // write the sequences without copying gaps
                writeASequenceToFile(aln, sequence_length, out, alisimulator->tree->root, alisimulator->tree->root);

            // close the file
            out.close();
        
            // show the output file name
            cout << "An alignment has just been exported to "<<file_path<<endl;
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, file_path);
        }
}

/**
*  merge and write all sequences to output files
*/
void mergeAndWriteSequencesToFiles(string file_path, AliSimulator *alisimulator, bool inference_mode){
    // in case with partitions -> merge & write sequences to a single/multiple files
    if (alisimulator->tree->isSuperTree())
    {
        PhyloSuperTree *super_tree = ((PhyloSuperTree *)alisimulator->tree);
        
        // merge phylotrees using the same alignment file and the same sequence_type with each other -> write sequences to the corresponding output files
        for (int i = 0; i < super_tree->size(); i++)
        {
            // ignore the current phylotree if it has already merged
            bool already_merged = false;
            for (int j = 0; j < i; j++)
                if (!super_tree->at(i)->aln->aln_file.compare(super_tree->at(j)->aln->aln_file)
                    && super_tree->at(i)->aln->seq_type == super_tree->at(j)->aln->seq_type)
                {
                    already_merged = true;
                    break;
                }
            
            // merge the sequences of the current phylotree and other phylotrees which use the same alignment file to the super_tree
            if (!already_merged)
            {
                // dummy variables
                string partition_list = "";
                int partition_count = 0;
                
                // the total number of states (sites) of all partitions which use the same alignment file.
                int total_num_states = 0;
                for (int j = i; j < super_tree->size(); j++)
                {
                    IQTree *current_tree = (IQTree*) super_tree->at(j);
                    if (!super_tree->at(i)->aln->aln_file.compare(current_tree->aln->aln_file) && super_tree->at(i)->aln->seq_type == super_tree->at(j)->aln->seq_type)
                    {
                        // update partition_list
                        partition_list = partition_list + "_" + super_tree->at(j)->aln->name;
                        partition_count++;
                        
                        // extract site_ids of the partition
                        const char* info_spec = ((SuperAlignment*) super_tree->aln)->partitions[j]->CharSet::position_spec.c_str();
                        IntVector site_ids;
                        current_tree->aln->extractSiteID(current_tree->aln, info_spec, site_ids, super_tree->getAlnNSite());
                        
                        // copy alignment from the current tree to the super_tree
                        copySequencesToSuperTree(site_ids, super_tree->getAlnNSite(), super_tree, current_tree->root, current_tree->root);
                        
                        // update total_num_states
                        total_num_states += current_tree->aln->getNSite();
                    }
                }
                
                // initialize output file name
                if (partition_count == super_tree->size())
                    partition_list = "_full";
                
                // write the merged sequences to the output file for the current cluster of partitions
                writeSequencesToFile(file_path + partition_list, super_tree->at(i)->aln, total_num_states, alisimulator, inference_mode);
            }
        }
    }
    // other cases (without partitions), just write sequences to a single file
    else
    {
        int sequence_length = alisimulator->expected_num_sites;
        writeSequencesToFile(file_path, alisimulator->tree->aln, sequence_length, alisimulator, inference_mode);
    }
}

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(Alignment *aln, int sequence_length, ofstream &out, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        out<<node->name<<" "<<convertEncodedSequenceToReadableSequence(aln, sequence_length, node->sequence)<<endl;
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFile(aln, sequence_length, out, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
*/
string convertEncodedSequenceToReadableSequence(Alignment *aln, int sequence_length, IntVector sequence)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // convert states one by one
    ASSERT(sequence_length <= sequence.size());
    for (int i = 0; i < sequence_length; i++)
        output_sequence = output_sequence + aln->convertStateBackStr(sequence[i]);

    return output_sequence;
}

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
void writeASequenceToFileWithGaps(Alignment *aln, int sequence_length, vector<string> seq_names, vector<string> sequences, ofstream &out, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // retrieve the input sequence of the current node
        for (int i = 0; i < sequences.size(); i++)
            if (!seq_names[i].compare(node->name)){
                out <<node->name <<" "<<convertEncodedSequenceToReadableSequenceWithGaps(aln, sequence_length, sequences[i], node->sequence) << endl;
            }
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFileWithGaps(aln, sequence_length, seq_names, sequences, out, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...) with gaps copied from the input sequence
*/
string convertEncodedSequenceToReadableSequenceWithGaps(Alignment *aln, int sequence_length, string input_sequence, IntVector sequence)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // convert states one by one
    ASSERT(sequence_length <= sequence.size());
    for (int i = 0; i < sequence_length; i++)
    {
        // get the number of sites per each state
        int num_sites_per_state = aln->seq_type == SEQ_CODON?3:1;
        
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
            }
            // insert gaps in normal cases
            else
                output_sequence = output_sequence + "-";
        }
        // if it's not a gap
        else
        {
            output_sequence = output_sequence + aln->convertStateBackStr(sequence[i]);
        }
    }

    return output_sequence;
}
