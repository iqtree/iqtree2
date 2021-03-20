/*
 *  alisim.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#ifndef alisim_h
#define alisim_h
#include "utils/tools.h"
#include "tree/node.h"
#include "model/modelmarkov.h"
#include "tree/iqtree.h"

using namespace std;

/**
*  execute Alignment Simulator (AliSim)
*/
void runAliSim(Params params);

/**
*  show all input parameters for AliSim
*/
void showParameters(Params params);

/**
*  initialize an IQTree instance from input file
*/
IQTree *initializeIQTreeFromTreeFile(Params params, char* seq_type);

/**
*  initialize an Alignment instance for IQTree
*/
void initializeAlignment(char* seq_type, IQTree *tree);

/**
*  iteratively add name of all leaf nodes into the alignment instance
*/
void addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad);

/**
*  initialize a Model instance for IQTree
*/
void initializeModel(Params params, IQTree *tree);

/**
*  generate an alignment from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateSingleDatasetFromSingleTree(Params params, IQTree *tree, string output_filepath);

/**
*  get the ancestral sequence for the root node (from an input file or randomly generated)
*/
IntVector getAncestralSequence(Params params, IQTree *tree);

/**
*  retrieve the ancestral sequence for the root node from an input file
*/
IntVector retrieveAncestralSequenceFromInputFile(int sequence_position, IQTree *tree);

/**
*  randomly generate the ancestral sequence for the root node
*/
IntVector generateRandomSequence(int sequence_length, IQTree *tree);

/**
*  simulate sequences for all nodes in the tree
*/
void simulateSeqsForTree(int sequence_length, IQTree *tree);

/**
*  simulate sequences for all nodes in the tree in the case WITHOUT Rate Heterogeneity (Gamma/FreeRate Model)
*/
void simulateSeqsWithoutRH(int sequence_length, ModelSubst *model, int max_num_states, Node *node, Node *dad, double invariant_proportion);

/**
*  simulate the sequence for the current node in the case WITHOUT Rate Heterogeneity (Gamma/FreeRate Model) WITHOUT Invariant sites
*/
IntVector estimateSeqWithoutRHWithoutIS(double *trans_matrix, int max_num_states, int sequence_length, IntVector dad_sequence);

/**
*  simulate the sequence for the current node in the case WITHOUT Rate Heterogeneity (Gamma/FreeRate Model) WITH Invariant sites
*/
IntVector estimateSeqWithoutRHWithIS(double *trans_matrix, int max_num_states, int sequence_length, IntVector dad_sequence, double invariant_proportion);

/**
*  simulate sequences for all nodes in the tree in the case WITH Rate Heterogeneity (Gamma/FreeRate Model)
*/
void simulateSeqsWithRH(int sequence_length, ModelSubst *model, RateHeterogeneity *rate_heterogeneity, int max_num_states, Node *node, Node *dad, double invariant_proportion);

/**
*  simulate the sequence for the current node in the case WITH Rate Heterogeneity (Gamma/FreeRate Model) WITHOUT Invariant sites
*/
IntVector estimateSeqWithRHWithoutIS(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, int sequence_length, double branch_length, IntVector dad_sequence);

/**
*  simulate the sequence for the current node in the case WITH Rate Heterogeneity (Gamma/FreeRate Model) WITH Invariant sites
*/
IntVector estimateSeqWithRHWithIS(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, int sequence_length, double branch_length, IntVector dad_sequence, double invariant_proportion);

/**
*  estimate the state for the current site of a node in the case WITH Rate Heterogeneity (Gamma/FreeRate Model)
*/
int estimateStateWithRH(ModelSubst *model, RateHeterogeneity *rate_heterogeneity, double *category_probability_matrix, double *trans_matrix, int max_num_states, double branch_length, int dad_state);

/**
*  get a random item from a set of items with a probability array
*/
int getRandomItemWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int num_items);

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, IQTree *tree, int sequence_length);

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(Alignment *aln, ofstream &out, Node *node, Node *dad);

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
*/
string convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence);

#endif /* alisim_h */
