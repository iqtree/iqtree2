/*
 *  alisim.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#ifndef alisim_h
#define alisim_h

#include "utils/tools.h"
#include "alignment/alisimulatorinvar.h"
#include "alignment/alisimulatorheterogeneityinvar.h"
#include "phyloanalysis.h"
#include "tree/phylosupertree.h"
#include "utils/gzstream.h"
#include <regex>
#include <string.h>

/**
*  execute Alignment Simulator (AliSim)
*/
void runAliSim(Params &params, Checkpoint *checkpoint);

/**
*  execute AliSim without inference
*/
void runAliSimWithoutInference(Params params, IQTree *&tree);

/**
*  inferring input parameters for AliSim
*/
void inferInputParameters(Params &params, Checkpoint *checkpoint, IQTree *&tree, Alignment *&aln);

/**
*  generate a random tree
*/
void generateRandomTree(Params &params);

/**
*  show all input parameters for AliSim
*/
void showParameters(Params params, bool is_partition_model);

/**
*  retrieve the ancestral sequence for the root node from an input file
*/
IntVector retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator);

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator);

/**
*  generate a partition alignment from a single simulator
*/
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *alisimulator, IntVector ancestral_sequence);

/**
*  compute the total sequence length of all partitions
*/
int computeTotalSequenceLengthAllPartitions(PhyloSuperTree *super_tree);

/**
*  copy sequences of leaves from a partition tree to super_tree
*/
void copySequencesToSuperTree(IntVector site_ids, int expected_num_states_super_tree, IQTree *super_tree, Node *node, Node *dad);

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, IntVector aln_ids, AliSimulator *alisimulator);

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(IQTree *tree, IntVector aln_ids, ofstream &out, Node *node, Node *dad);

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
void writeASequenceToFileWithGaps(IQTree *tree, IntVector aln_ids, vector<string> seq_names, vector<string> sequences, ofstream &out, Node *node, Node *dad);

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...)
*/
string convertEncodedSequenceToReadableSequence(IQTree *tree, IntVector aln_ids, IntVector sequence);

/**
*  convert an encoded sequence (with integer numbers) to a readable sequence (with ACGT...) with gaps copied from the input sequence
*/
string convertEncodedSequenceToReadableSequenceWithGaps(IQTree *tree, IntVector aln_ids, string input_sequence, IntVector sequence);

#endif /* alisim_h */
