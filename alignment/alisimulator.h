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
#include "superalignment.h"
#include "superalignmentunlinked.h"
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
    vector<short int> generateRandomSequence(int sequence_length, bool initial_freqs = true);
    
    /**
    *  randomly generate the base frequencies
    */
    void generateRandomBaseFrequencies(double *base_frequencies);
    
    /**
    *  get a random item from a set of items with a probability array
    */
    int getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items);
    
    /**
    *  get a random item from a set of items with an accumulated probability array by binary search starting at the max probability
    */
    int getRandomItemWithAccumulatedProbMatrixMaxProbFirst(double *accumulated_probability_maxtrix, int starting_index, int num_columns, int max_prob_position);

    /**
    *  convert an probability matrix into an accumulated probability matrix
    */
    void convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns);

    /**
    *  binary search an item from a set with accumulated probability array
    */
    int binarysearchItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, double random_number, int start, int end, int first);
    
    /**
    *  binary search an item from a set with accumulated probability array
    */
    int binarysearchItemWithAccumulatedProbabilityMatrix(vector<double> accumulated_probability_maxtrix, double random_number, int start, int end, int first);
    
    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int &sequence_length, ModelSubst *model, double *trans_matrix, Node *node, Node *dad, ostream &out, vector<string> state_mapping, map<string, string> input_msa);
    
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
    void getOnlyVariantSites(vector<short int> variant_state_mask, Node *node, Node *dad);
    
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
    short int extractMaxTaxaNameLength();
    
    /**
        selecting & permuting sites (FunDi models)
    */
    vector<FunDi_Item> selectAndPermuteSites(double proportion, int num_sites);
    
    /**
        permuting selected sites (FunDi models)
    */
    void permuteSelectedSites(vector<FunDi_Item> fundi_items, Node* node);
    
    /**
        process delayed Fundi if it is delayed due to Insertion events
    */
    void processDelayedFundi(Node *node, Node *dad);
    
    /**
        writing and deleting simulated sequence immediately if possible
    */
    void writeAndDeleteSequenceImmediatelyIfPossible(ostream &out, vector<string> state_mapping, map<string,string> input_msa, NeighborVec::iterator it, Node* node);
    
    /**
        branch-specific evolution
    */
    void branchSpecificEvolution(int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it);
    
    /**
        simulate a sequence for a node from a specific branch
    */
    void simulateASequenceFromBranch(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths = "");
    
    /**
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    virtual void simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths = "");
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    virtual void initVariables(int sequence_length, bool regenerate_root_sequence = false);
    
    /**
        regenerate the root sequence if the user has specified specific state frequencies in branch-specific model
    */
    void regenerateRootSequenceBranchSpecificModel(string freqs, int sequence_length, Node* root);
    
    /**
        generate a random sequence by state frequencies
    */
    vector<short int> generateRandomSequenceFromStateFreqs(int sequence_length, double* state_freqs, int max_prob_pos);
    
    /**
    *  export a sequence with gaps copied from the input sequence
    */
    string exportSequenceWithGaps(Node *node, int sequence_length, int num_sites_per_state, string input_sequence, vector<string> state_mapping);
    
    /**
        handle indels
    */
    void handleIndels(ModelSubst *model, int &sequence_length, Node *node, NeighborVec::iterator it, vector<short int> &indel_sequence, vector<int> &index_mapping_by_jump_step, SIMULATION_METHOD simulation_method);
    
    /**
        handle substitution events
    */
    void handleSubs(int sequence_length, double &total_sub_rate, vector<double> &sub_rate_by_site, vector<short int> &indel_sequence, int num_mixture_models);
    
    /**
        handle insertion events, return the insertion-size
    */
    int handleInsertion(int &sequence_length, vector<int> &index_mapping_by_jump_step, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method);
    
    /**
        handle deletion events, return the deletion-size
    */
    int handleDeletion(int sequence_length, vector<short int> &indel_sequence, double &total_sub_rate, vector<double> &sub_rate_by_site, SIMULATION_METHOD simulation_method);
    
    /**
        extract array of substitution rates and Jmatrix
    */
    double extractRatesJMatrix(ModelSubst *model);
    
    /**
        initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
    */
    virtual void initVariables4RateMatrix(double &total_sub_rate, int &num_gaps, vector<double> &sub_rate_by_site, vector<short int> sequence);
    
    /**
    *  insert a new sequence into the current sequence
    *
    */
    virtual void insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence);
    
    /**
    *  insert gaps into other nodes when processing Insertion Events
    *
    */
    void insertGapsForInsertionEvents(vector<int> index_mapping_by_jump_step, int stopping_node_id, Node *node, Node *dad, bool &stop_inserting_gaps);
    
    /**
    *  randomly select a valid position (not a deleted-site) for insertion/deletion event
    *
    */
    int selectValidPositionForIndels(int upper_bound, vector<short int> sequence);
    
    /**
        merge the simulated sequence with indel_sequence
    */
    virtual void mergeIndelSequence(Node* node, vector<short int> indel_sequence, vector<int> index_mapping_by_jump_step);
    
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
        compute the simulation threshold to switch between Rate matrix and Probability matrix
    */
    void computeSimThresh(int seq_length);
    
    /**
        change state of sites due to Error model
    */
    void changeSitesErrorModel(vector<int> sites, vector<short int> &sequence, double error_prop);
    
    /**
        handle DNA error
    */
    void handleDNAerr(double error_prop, vector<short int> &sequence, int model_index = -1);
    
public:
    
    IQTree *tree;
    Params *params;
    int max_num_states;
    int num_sites_per_state;
    int expected_num_sites;
    double partition_rate;
    double length_ratio = 1;
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
    
    // variables using for posterior mean rates/state frequencies
    bool applyPosMeanRate = false;
    double* ptn_state_freq = NULL;
    DoubleVector pattern_rates;
    
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
    virtual void simulateSeqsForTree(map<string,string> input_msa, string output_filepath = "");
    
    /**
    *  generate the current partition of an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
    */
    void generatePartitionAlignment(vector<short int> ancestral_sequence, map<string,string> input_msa, string output_filepath = "");
    
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
    static string convertNumericalStatesIntoReadableCharacters(Node *node, int sequence_length, int num_sites_per_state, vector<string> state_mapping);
    
    /**
    *  export pre_output string (containing taxon name and ">" or "space" based on the output format)
    *
    */
    static string exportPreOutputString(Node *node, InputType output_format, int max_length_taxa_name);
};

#endif /* alisimulator_h */
