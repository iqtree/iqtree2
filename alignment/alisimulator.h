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
    IntVector retrieveAncestralSequenceFromInputFile(int sequence_position);

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
    *  write a sequence of a node to an output file
    */
    void writeASequenceToFile(Alignment *aln, ofstream &out, Node *node, Node *dad);

    /**
    *  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
    */
    string convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence);
    
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
