//
//  alisimulator.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#ifndef alisimulator_h
#define alisimulator_h

#include "utils/tools.h"
#include "tree/node.h"
#include "tree/genometree.h"
#include "alignment/superalignment.h"
#include "alignment/superalignmentunlinked.h"
#include "tree/phylotreemixlen.h"
#include "tree/phylosupertree.h"
#include "tree/phylosupertreeplen.h"
#include "tree/phylosupertreeunlinked.h"
#include "model/modelmarkov.h"
#include "model/modelliemarkov.h"
#include "model/ratecontinuousgammainvar.h"
#include "tree/iqtree.h"
#include "main/phylotesting.h"
#include <random>
#include "utils/gzstream.h"
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "utils/MPIHelper.h"
#include "alignment/sequencechunkstr.h"

struct FunDi_Item {
  int selected_site;
  int new_position;
} ;

/**
 *  Specify 3 event types.
 */
enum EVENT_TYPE {
    INSERTION,
    DELETION,
    SUBSTITUTION
};

class AliSimulator{
protected:
    
    /**
    *  initialize an IQTree instance from input file
    */
    void initializeIQTreeFromTreeFile();

    /**
    *  initialize an Alignment instance for IQTree
    */
    void initializeAlignment(IQTree *tree, string model_fullname);

    /**
    *  iteratively add name of all leaf nodes into the alignment instance
    */
    void addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad);

    /**
    *  initialize a Model instance for IQTree
    */
    void initializeModel(IQTree *tree, string model_name);
    
    /**
    *  get state frequencies from model
    */
    void getStateFrequenciesFromModel(IQTree* tree, double *state_freqs);

    /**
    *  randomly generate the ancestral sequence for the root node
    *  by default (initial_freqs = true) freqs could be randomly generated if they are not specified
    */
    void generateRandomSequence(int sequence_length, vector<short int> &sequence, bool initial_freqs = true);
    
    /**
    *  randomly generate the base frequencies
    */
    void generateRandomBaseFrequencies(double *base_frequencies);
    
    /**
    *  get a random item from a set of items with a probability array
    */
    int getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items, int* rstream);
    
    /**
    *  get a random item from a set of items with an accumulated probability array by binary search starting at the max probability
    */
    int getRandomItemWithAccumulatedProbMatrixMaxProbFirst(double *accumulated_probability_maxtrix, int starting_index, int num_columns, int max_prob_position, int* rstream);

    /**
    *  convert an probability matrix into an accumulated probability matrix
    */
    void convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns, bool force_round_1 = true);

    /**
    *  binary search an item from a set with accumulated probability array
    */
    int binarysearchItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, double random_number, int start, int end, int first);
    
    /**
    *  binary search an item from a set with accumulated probability array
    */
    int binarysearchItemWithAccumulatedProbabilityMatrix(vector<double> &accumulated_probability_maxtrix, double random_number, int start, int end, int first);
    
    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    void simulateSeqs(int thread_id, int segment_start, int &segment_length, int &sequence_length, ModelSubst *model, double *trans_matrix, vector<vector<short int>> &sequence_cache, bool store_seq_at_cache, Node *node, Node *dad, ostream &out, vector<string> &state_mapping, map<string, string> input_msa, int* rstream);
    
    /**
    *  reset tree (by reset some variables of nodes)
    *
    */
    void resetTree(int &max_depth, bool store_seq_at_cache, Node *node = NULL, Node *dad = NULL);
    
    /**
    *  validate sequence length of codon
    *
    */
    void validataSeqLengthCodon();
    
    /**
        initialize state freqs for all model components (of a mixture model)
    */
    void intializeStateFreqsMixtureModel(IQTree* tree);
    
    /**
        create mask for variant states
    */
    void createVariantStateMask(vector<short int> &variant_state_mask, int &num_variant_states, int expected_num_variant_states, Node *node, Node *dad);
    
    /**
        remove all constant sites (in case with +ASC)
    */
    void removeConstantSites();
    
    /**
        only get variant sites
    */
    void getOnlyVariantSites(vector<short int> &variant_state_mask, Node *node, Node *dad);
    
    /**
        estimate length_ratio (for models with +ASC)
    */
    void estimateLengthRatio();
    
    /**
        show warning if base frequencies are set/unset correctly (only check DNA models)
    */
    void checkBaseFrequenciesDNAModels(IQTree* tree, string model_name);
    
    /**
        extract the maximum length of taxa names
    */
    void extractMaxTaxaNameLength();
    
    /**
        selecting & permuting sites (FunDi models)
    */
    void selectAndPermuteSites(vector<FunDi_Item> &fundi_items, double proportion, int num_sites);
    
    /**
        permuting selected sites (FunDi models)
    */
    void permuteSelectedSites(vector<FunDi_Item> &fundi_items, Node* node);
    
    /**
        process delayed Fundi if it is delayed due to Insertion events
    */
    void processDelayedFundi(Node *node, Node *dad);
    
    /**
        merge and write sequence in simulations with Indels or FunDi model
    */
    void mergeAndWriteSeqIndelFunDi(int thread_id, ostream &out, int sequence_length, vector<string> &state_mapping, map<string,string> input_msa, NeighborVec::iterator it, Node* node);
    
    /**
        write and delete the current chunk of sequence if possible
    */
    void writeAndDeleteSequenceChunkIfPossible(int thread_id, int segment_start, int segment_length, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, bool store_seq_at_cache, ostream &out, vector<string> &state_mapping, map<string, string> input_msa, NeighborVec::iterator it, Node* node);
    
    /**
        branch-specific evolution by multi threads
    */
    void branchSpecificEvolution(int thread_id, int sequence_length, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, bool store_seq_at_cache, double *trans_matrix, Node *node, NeighborVec::iterator it);
    
    
    /**
        branch-specific evolution by the master thread
    */
    void branchSpecificEvolutionMasterThread(int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it);
    
    /**
        simulate a sequence for a node from a specific branch
    */
    void simulateASequenceFromBranch(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths = "");
    
    /**
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    virtual void simulateASequenceFromBranchAfterInitVariables(int segment_start, ModelSubst *model, double *trans_matrix, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, Node *node, NeighborVec::iterator it, int* rstream, string lengths = "");
    
    /**
        initialize variables
    */
    void initVariables(int sequence_length, string output_filepath, vector<string> &state_mapping, ModelSubst *model, int &default_segment_length, int &max_depth, bool &write_sequences_to_tmp_data, bool &store_seq_at_cache);
    
    /**
        process after simulating sequences
    */
    void postSimulateSeqs(int sequence_length, string output_filepath, bool write_sequences_to_tmp_data);
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    virtual void initVariablesRateHeterogeneity(int sequence_length, bool regenerate_root_sequence = false);
    
    /**
        regenerate the root sequence if the user has specified specific state frequencies in branch-specific model
    */
    void regenerateRootSequenceBranchSpecificModel(string freqs, int sequence_length, Node* root);
    
    /**
        generate a random sequence by state frequencies
    */
    void generateRandomSequenceFromStateFreqs(int sequence_length, vector<short int> &sequence, double* state_freqs, int max_prob_pos);
    
    /**
    *  export a sequence with gaps copied from the input sequence
    */
    void exportSequenceWithGaps(vector<short int> &sequence_chunk, string &output, int sequence_length, int num_sites_per_state, string input_sequence, vector<string> &state_mapping, int segment_start = 0, int segment_length = -1);
    
    /**
        handle indels
    */
    void simulateSeqByGillespie(int segment_start, int &segment_length, ModelSubst *model, vector<short int> &node_seq_chunk, int &sequence_length, NeighborVec::iterator it, SIMULATION_METHOD simulation_method, int *rstream);
    
    /**
        handle substitution events
    */
    void handleSubs(int segment_start, double &total_sub_rate, vector<double> &sub_rate_by_site, vector<short int> &indel_sequence, int num_mixture_models, int* rstream);
    
    /**
        handle insertion events, return the insertion-size
    */
    int handleInsertion(int &sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method);
    
    /**
        handle deletion events, return the deletion-size
    */
    int handleDeletion(int sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method);
    
    /**
        extract array of substitution rates and Jmatrix
    */
    void extractRatesJMatrix(ModelSubst *model);
    
    /**
        initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
    */
    virtual void initVariables4RateMatrix(int segment_start, double &total_sub_rate, int &num_gaps, vector<double> &sub_rate_by_site, vector<short int> &sequence);
    
    /**
    *  insert a new sequence into the current sequence
    *
    */
    virtual void insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence);
    
    /**
    *  update internal sequences due to Indels
    *
    */
    void updateInternalSeqsIndels(GenomeTree* genome_tree, int seq_length, Node *node);
    
    /**
    *  update all simulated internal seqs from root to the current node due to insertions
    *
    */
    void updateInternalSeqsFromRootToNode(GenomeTree* genome_tree, int seq_length, int stopping_node_id, Node *node, Node* dad, bool &stop_inserting_gaps);
    
    /**
    *  update internal seqs on the path from the current phylonode to root due to insertions
    *
    */
    void updateInternalSeqsFromNodeToRoot(GenomeTree* genome_tree, int seq_length, Node *node);
    
    /**
    *  randomly select a valid position (not a deleted-site) for insertion/deletion event
    *
    */
    int selectValidPositionForIndels(int upper_bound, vector<short int> &sequence);
    
    /**
        generate indel-size from its distribution
    */
    int generateIndelSize(IndelDistribution indel_dis);
    
    /**
        compute mean of deletion-size
    */
    double computeMeanDelSize(int sequence_length);
    
    /**
        root Tree
    */
    void rootTree();
    
    /**
        compute the switching param to switch between Rate matrix and Probability matrix
    */
    void computeSwitchingParam(int seq_length);
    
    /**
        change state of sites due to Error model
    */
    void changeSitesErrorModel(vector<int> sites, vector<short int> &sequence, double error_prop, int* rstream);
    
    /**
        handle DNA error
    */
    void handleDNAerr(int segment_start, double error_prop, vector<short int> &sequence, int* rstream, int model_index = -1);
    
    /**
        TRUE if posterior mean rate can be used
    */
    bool canApplyPosteriorRateHeterogeneity();
    
    /**
        init Site to PatternID
    */
    void initSite2PatternID(int length);
    
    /**
        temporarily write internal states to file (when using Indels)
    */
    void writeInternalStatesIndels(Node* node, ostream &out);
    
    /**
        separate root sequence into chunks
    */
    void separateSeqIntoChunks(Node* node);
    
    /**
        merge chunks into a single sequence
    */
    void mergeChunks(Node* node);
    
    /**
        merge chunks into a single sequence for all nodes in tree
    */
    void mergeChunksAllNodes(Node* node = NULL, Node* dad = NULL);
    
    /**
        init the output file
    */
    void initOutputFile(ostream *&out, int thread_id, int actual_segment_length, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data);
    
    /**
        open an output stream
    */
    void openOutputStream(ostream *&out, string output_filepath, std::ios_base::openmode open_mode, bool force_uncompression = false);
    
    /**
        close an output stream
    */
    void closeOutputStream(ostream *&out, bool force_uncompression = false);
    
    /**
        visit cache of each thread in round robin then write sequence chunks one by one
    */
    void writeSeqChunkFromCache(ostream *&output);
    
    /**
        write all remaining chunks from cache
    */
    void writeAllSeqChunkFromCache(ostream *&output);
    
    /**
        cache a sequence chunk (in readable string) into the cache (writing queue)
    */
    void cacheSeqChunkStr(int64_t pos, string seq_chunk_str, int thread_id);
    
    /**
        wait for all threads to reach the manually-implemented-barrier
    */
    void waitAtBarrier(const unsigned short int barrier_count, Node* node);
    
    /**
    *  simulate sequences with AliSim-OpenMP-IM algorithm
    */
    void executeIM(int thread_id, int &sequence_length, int default_segment_length, ModelSubst *model, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data, bool store_seq_at_cache, int max_depth, vector<string> &state_mapping);
    
    /**
    *  simulate sequences with AliSim-OpenMP-EM algorithm
    */
    void executeEM(int thread_id, int &sequence_length, int default_segment_length, ModelSubst *model, map<string,string> input_msa, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data, bool store_seq_at_cache, int max_depth, vector<string> &state_mapping);
    
    /**
        merge output files when using multiple threads
    */
    void mergeOutputFiles(ostream *&single_output, int thread_id, string output_filepath, std::ios_base::openmode open_mode, bool write_sequences_to_tmp_data);
    
    /**
        output a sequence to file (if using AliSim-OpenMP-EM) or store it to common cache (if using AliSim-OpenMP-IM)
    */
    void outputOneSequence(Node* node, string &output, int thread_id, int segment_start, ostream &out);
    
public:
    
    IQTree *tree;
    Params *params;
    int max_num_states;
    int num_sites_per_state;
    int expected_num_sites;
    double partition_rate;
    double length_ratio = 1;
    double inverse_length_ratio = 1;
    short int max_length_taxa_name = 10;
    vector<FunDi_Item> fundi_items;
    short int STATE_UNKNOWN;
    vector<short int> site_specific_model_index;
    vector<short int> site_specific_rate_index;
    vector<double> site_specific_rates;
    const int RATE_ZERO_INDEX = -1;
    const int RATE_ONE_INDEX = 0;
    double* sub_rates;
    double* Jmatrix;
    double* mixture_accumulated_weight = NULL;
    int mixture_max_weight_pos = 0;
    int seq_length_indels = 0; // final seq_length due to indels
    map<string, Node*> map_seqname_node; // mapping sequence name to Node (using when temporarily write sequences at tips to tmp_data file when simulating Indels)
    Insertion* latest_insertion = NULL;
    Insertion* first_insertion = NULL;
    
    // variables to output sequences with multiple threads
    uint64_t starting_pos = 0;
    uint64_t output_line_length = 0;
    uint64_t seq_name_length = 0;
    int num_threads = 1;
    int num_simulating_threads = 1;
    int num_thread_done = 0;
    vector<SequenceChunkStr> seq_str_cache;
    vector<int> cache_start_indexes;
    int cache_size_per_thread;
    bool force_output_PHYLIP = false;
    
    // variables using for posterior mean rates/state frequencies
    bool applyPosRateHeterogeneity = false;
    double* ptn_state_freq = NULL;
    double* ptn_accumulated_state_freq = NULL;
    double* ptn_model_dis = NULL;
    double* ptn_accumulated_rate_dis = NULL;
    DoubleVector pattern_rates;
    IntVector site_to_patternID;
    
    /**
        constructor
    */
    AliSimulator(){};
    
    /**
        constructor
    */
    AliSimulator(Params *params, int expected_number_sites = -1, double new_partition_rate = 1);
    
    /**
        constructor
    */
    AliSimulator(Params *params, IQTree *tree, int expected_number_sites = -1, double new_partition_rate = 1);
    
    /**
        deconstructor
    */
    ~AliSimulator();
    
    /**
    *  simulate sequences for all nodes in the tree
    */
    virtual void simulateSeqsForTree(map<string,string> input_msa, string output_filepath = "", std::ios_base::openmode open_mode = std::ios_base::out);
    
    /**
    *  generate the current partition of an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
    */
    void generatePartitionAlignment(vector<short int> &ancestral_sequence, map<string,string> input_msa, string output_filepath = "", std::ios_base::openmode open_mode = std::ios_base::out);
    
    /**
    *  update the expected_num_sites due to the change of the sequence_length
    */
    void refreshExpectedNumSites();
    
    /**
    *  initialize state_mapping (mapping from states into characters)
    *
    */
    static void initializeStateMapping(int num_sites_per_state, Alignment *aln, vector<string> &state_mapping);
    
    /**
    *  convert numerical states into readable characters
    *
    */
    static void convertNumericalStatesIntoReadableCharacters(vector<short int> &sequence_chunk, string &output, int sequence_length, int num_sites_per_state, vector<string> &state_mapping, int segment_length = -1);
    
    /**
    *  export pre_output string (containing taxon name and ">" or "space" based on the output format)
    *
    */
    static string exportPreOutputString(Node *node, InputType output_format, int max_length_taxa_name, bool force_PHYLIP = false);

    /**
    *  update new genome from original genome and the genome tree for each tips (due to Indels)
    */
    void updateNewGenomeIndels(int seq_length);
};

#endif /* alisimulator_h */
