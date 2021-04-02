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
#include "model/modelmarkov.h"
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
    *  generate an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
    */
    void generateSingleDatasetFromSingleTree(string output_filepath, IntVector ancestral_sequence);

    /**
    *  retrieve the ancestral sequence for the root node from an input file
    */
    IntVector retrieveAncestralSequenceFromInputFile(char *aln_filepath, string sequence_name);
    
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
    *  write all sequences of a tree to an output file
    */
    void writeSequencesToFile(string file_path);
    
    /**
    *  load sequences from input file
    */
    void loadSequences(char *file_path, vector<string> &seq_names, vector<string> &sequences);
    
    /**
    *  write a sequence of a node to an output file
    */
    void writeASequenceToFile(Alignment *aln, ofstream &out, Node *node, Node *dad);
    
    /**
    *  write a sequence of a node to an output file with gaps copied from the input sequence
    */
    void writeASequenceToFileWithGaps(Alignment *aln, vector<string> seq_names, vector<string> sequences, ofstream &out, Node *node, Node *dad);

    /**
    *  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
    */
    string convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence);
    
    /**
    *  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...) with gaps copied from the input sequence
    */
    string convertEncodedSequenceToReadableSequenceWithGaps(Alignment *aln, string input_sequence, IntVector sequence);
    
    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int sequence_length, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad);
    
public:
    
    IQTree *tree;
    Params *params;
    
    /**
        constructor
    */
    AliSimulator(){};
    
    /**
        constructor
    */
    AliSimulator(Params *params);
    
    /**
        deconstructor
    */
    ~AliSimulator();

    /**
    *  show all input parameters for AliSim
    */
    virtual void showParameters();
    
    /**
    *  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
    */
    void generateMultipleAlignmentsFromSingleTree();

    /**
    *  simulate sequences for all nodes in the tree
    */
    virtual void simulateSeqsForTree();
    
};

#endif /* alisimulator_h */
