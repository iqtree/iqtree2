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
#include "tree/phylotreemixlen.h"
#include "model/modelmarkov.h"
#include "model/modelliemarkov.h"
#include "model/ratecontinuousgammainvar.h"
#include "tree/iqtree.h"
#include "main/phylotesting.h"
#include <random>

class AliSimulator{
protected:
    
    /**
    *  initialize an IQTree instance from input file
    */
    void initializeIQTreeFromTreeFile();

    /**
    *  initialize an Alignment instance for IQTree
    */
    void initializeAlignment();

    /**
    *  iteratively add name of all leaf nodes into the alignment instance
    */
    void addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad);

    /**
    *  initialize a Model instance for IQTree
    */
    void initializeModel();
    
    /**
    *  get state frequencies from model
    */
    void getStateFrequenciesFromModel(double *state_freqs);

    /**
    *  randomly generate the ancestral sequence for the root node
    */
    IntVector generateRandomSequence(int sequence_length);
    
    /**
    *  randomly generate the base frequencies
    */
    void generateRandomBaseFrequencies(double *base_frequencies, int max_num_bases);
    
    /**
    *  get a random item from a set of items with a probability array
    */
    int getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items);

    /**
    *  get a random item from a set of items with an accumulated probability array by binary search
    */
    int getRandomItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, int starting_index, int num_columns);

    /**
    *  convert an probability matrix into an accumulated probability matrix
    */
    void convertProMatrixIntoAccumulatedProMatrix(double *probability_maxtrix, int num_rows, int num_columns);

    /**
    *  binary search an item from a set with accumulated probability array
    */
    int binarysearchItemWithAccumulatedProbabilityMatrix(double *accumulated_probability_maxtrix, double random_number, int start, int end, int first);
    
    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int sequence_length, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad);
    
    /**
    *  validate sequence length of codon
    *
    */
    void validataSeqLengthCodon();
    
    /**
        initialize state freqs for all model components (of a mixture model)
    */
    virtual void intializeStateFreqsMixtureModel(){
        // will be override in alisimulatorheterogeneity
    };
    
    /**
        create mask for variant states
    */
    void createVariantStateMask(IntVector &variant_state_mask, int &num_variant_states, int expected_num_variant_states, Node *node, Node *dad);
    
    /**
        remove all constant sites (in case with +ASC)
    */
    void removeConstantSites();
    
    /**
        only get variant sites
    */
    void getOnlyVariantSites(IntVector variant_state_mask, Node *node, Node *dad);
    
public:
    
    IQTree *tree;
    Params *params;
    int num_sites_per_state;
    int expected_num_sites;
    double partition_rate;
    
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
    virtual void simulateSeqsForTree();
    
    /**
    *  generate the current partition of an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
    */
    void generatePartitionAlignment(IntVector ancestral_sequence);
    
    /**
    *  update the expected_num_sites due to the change of the sequence_length
    */
    void refreshExpectedNumSites();
};

#endif /* alisimulator_h */
