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
    if (params.aln_file)
    {
        inferInputParameters(params, checkpoint, tree, aln);
    }
    
    // run AliSim without inference
    runAliSimWithoutInference(params, tree);
    
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
        params.user_file = new char[strlen(params.aln_file) + 10];
        strcpy(params.user_file, params.aln_file);
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
*  execute AliSim without inference
*/
void runAliSimWithoutInference(Params params, IQTree *&tree)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator;
    if (tree && params.aln_file)
        alisimulator = new AliSimulator(&params, tree);
    else
        alisimulator = new AliSimulator(&params);
    
    // show parameters
    showParameters(params, alisimulator->tree->isSuperTree());
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    generateMultipleAlignmentsFromSingleTree(alisimulator);
    
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
    super_alisimulator->params->alisim_sequence_length = nsite;
    outWarning("Sequence length is now set equally to the length of ancestral sequence.");
    super_alisimulator->refreshExpectedNumSites();
    
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
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator)
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
        +"_"+convertIntToString(i)+".phy";
        
        // mapping from site_id to its algnmennt id
        IntVector aln_ids;
        
        // generate multiple alignments one by one
        if (super_alisimulator->tree->isSuperTree())
        {
            PhyloSuperTree* super_tree = ((PhyloSuperTree*) super_alisimulator->tree);
            int total_expected_num_states = super_tree->getAlnNSite();
            
            // override sequence_length by the total length of all partitions
            super_alisimulator->tree->params->alisim_sequence_length = computeTotalSequenceLengthAllPartitions(super_tree);
            outWarning("The sequence_length is now set at equal to the total length of all partitions");
            super_alisimulator->refreshExpectedNumSites();
            
            // resize aln_ids
            aln_ids.resize(total_expected_num_states);
            
            for (int partition_index = 0; partition_index < super_tree->size(); partition_index++)
            {
                // get variables
                IQTree *current_tree = (IQTree*) super_tree->at(partition_index);
                int expected_num_states_current_tree = current_tree->aln->getNSite();
                
                // extract site_ids of the partition
                const char* info_spec = ((SuperAlignment*) super_tree->aln)->partitions[partition_index]->CharSet::position_spec.c_str();
                IntVector site_ids;
                current_tree->aln->extractSiteID(current_tree->aln, info_spec, site_ids, total_expected_num_states);
                
                // mapping from site_id to its algnmennt id
                for (int site_id: site_ids)
                    aln_ids[site_id] = partition_index;

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
                
                // copy alignment from the current tree to the super_tree
                copySequencesToSuperTree(site_ids, total_expected_num_states, super_tree, current_tree->root, current_tree->root);
            }
        }
        else
            generatePartitionAlignmentFromSingleSimulator(super_alisimulator, ancestral_sequence);
        
        // write alignment to file
        writeSequencesToFile(output_filepath, aln_ids, super_alisimulator);
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
void writeSequencesToFile(string file_path, IntVector aln_ids, AliSimulator *alisimulator)
{
try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(file_path.c_str());

        // write the first line <#taxa> <length_of_sequence>
        int leaf_num = alisimulator->tree->leafNum - ((alisimulator->tree->root->isLeaf() && alisimulator->tree->root->name == ROOT_NAME)?1:0);
        out <<leaf_num<<" "<<alisimulator->params->alisim_sequence_length<< endl;

        bool write_seq_without_gaps = true;
        // write senquences of leaf nodes to file with/without gaps copied from the input sequence
        if (alisimulator->params->aln_file && !alisimulator->params->alisim_no_copy_gaps)
        {
            write_seq_without_gaps = false;
            
            // in normal case (without partition) -> using the sequence type of the current tree to load the input sequence
            string sequence_type = alisimulator->tree->aln->sequence_type;
            // in case with partitions -> using the sequence type of the first phylotree to load the input sequence
            if (alisimulator->tree->isSuperTree())
            {
                // using the first phylotree to load the ancestral sequence
                sequence_type = ((PhyloSuperTree*) alisimulator->tree)->at(0)->aln->sequence_type;
                
                // make sure all partitions are using the same sequence_type
                for (int i = 1; i < ((PhyloSuperTree*) alisimulator->tree)->size(); i++)
                    if (((PhyloSuperTree*) alisimulator->tree)->at(i)->aln->seq_type != ((PhyloSuperTree*) alisimulator->tree)->at(0)->aln->seq_type)
                        write_seq_without_gaps = true;
            }
            
            // if gaps could not be copied due to the differences of sequence types of partitions -> show a warning and just write the sequences without gaps
            if (write_seq_without_gaps)
            {
                // show a warning gaps are not copied due to the differences of sequence types of partitions
                outWarning("Gaps are not copied to the output sequences since the sequence types of all partitions are not the same as each other.");
            }
            // otherwise, copying gaps from the input alignment to the output sequences
            else
            {
                // load input sequences (with gaps)
                vector<string> seq_names;
                vector<string> sequences;
                // read sequences from the input file
                Alignment *aln = new Alignment();
                int nseq = 0, nsite = 0;
                char *sequence_type_char = strcpy(new char[sequence_type.length() + 1], sequence_type.c_str());
                aln->extractSequences(alisimulator->params->aln_file, sequence_type_char, sequences, nseq, nsite);
                seq_names = aln->getSeqNames();

                // show a warning if the length of input alignment is unequal to that of simulated sequence
                if (sequences.size() > 0 && sequences[0].length() != alisimulator->params->alisim_sequence_length)
                    outWarning("The sequence length of the input alignment is unequal to that of that simulated sequences. Thus, only gaps in the first MIN(input_sequence_length, simulated_sequence_length) sites are copied.");

                // write simulated sequence with the gaps copied from the input sequence
                writeASequenceToFileWithGaps(alisimulator->tree, aln_ids, seq_names, sequences, out, alisimulator->tree->root, alisimulator->tree->root);
            }
        }

        // write the sequences without copying gaps
        if (write_seq_without_gaps)
            writeASequenceToFile(alisimulator->tree, aln_ids, out, alisimulator->tree->root, alisimulator->tree->root);

        // close the file
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_path);
    }
}

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(IQTree *tree, IntVector aln_ids, ofstream &out, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        out<<node->name<<" "<<convertEncodedSequenceToReadableSequence(tree, aln_ids, node->sequence)<<endl;
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFile(tree, aln_ids, out, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
*/
string convertEncodedSequenceToReadableSequence(IQTree *tree, IntVector aln_ids, IntVector sequence)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // if a partition model is used -> using the corresponding alignment (mapping by aln_ids) to convert each state
    if (tree->isSuperTree())
    {
        // make sure the size of aln_ids is equal to that of sequence
        ASSERT(aln_ids.size() == sequence.size());
        
        // convert states one by one
        for (int i = 0; i < sequence.size(); i++)
        {
            // get the corresponding alignment for the current state
            Alignment *aln = ((PhyloSuperTree *) tree)->at(aln_ids[i])->aln;
            output_sequence = output_sequence + (aln->convertStateBackStr(sequence[i]));
        }
    }
    else {
        // get the alignment from the tree
        Alignment *aln = tree->aln;
        
        // convert states one by one
        for (int i = 0; i < sequence.size(); i++)
            output_sequence = output_sequence + aln->convertStateBackStr(sequence[i]);
    }
    
    return output_sequence;
}

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
void writeASequenceToFileWithGaps(IQTree *tree, IntVector aln_ids, vector<string> seq_names, vector<string> sequences, ofstream &out, Node *node, Node *dad)
{
    if (node->isLeaf() && node->name!=ROOT_NAME) {
        // retrieve the input sequence of the current node
        for (int i = 0; i < sequences.size(); i++)
        if (!seq_names[i].compare(node->name)){
            out <<node->name <<" "<<convertEncodedSequenceToReadableSequenceWithGaps(tree, aln_ids, sequences[i], node->sequence) << endl;
        }
    }
    
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        writeASequenceToFileWithGaps(tree, aln_ids, seq_names, sequences, out, (*it)->node, node);
    }
}

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...) with gaps copied from the input sequence
*/
string convertEncodedSequenceToReadableSequenceWithGaps(IQTree *tree, IntVector aln_ids, string input_sequence, IntVector sequence)
{
    // initialize the output_sequence
    string output_sequence = "";
    
    // if a partition model is used -> make sure the size of aln_ids is equal to that of sequence
    if (tree->isSuperTree())
        ASSERT(aln_ids.size() == sequence.size());
    
    // convert states one by one
    for (int i = 0; i < sequence.size(); i++)
    {
        // get the corresponding alignment for the current state to convert it
        Alignment *aln = tree->isSuperTree()?((PhyloSuperTree *) tree)->at(aln_ids[i])->aln:tree->aln;
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
