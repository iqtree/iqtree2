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
*  execute AliSim Simulation
*/
void executeSimulation(Params params, IQTree *&tree, bool inference_mode);

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
vector<short int> retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator);

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator, bool inference_mode);

/**
*  generate a partition alignment from a single simulator
*/
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *alisimulator, vector<short int> ancestral_sequence, string output_filepath = "");

/**
*  compute the total sequence length of all partitions
*/
int computeTotalSequenceLengthAllPartitions(PhyloSuperTree *super_tree);

/**
*  copy sequences of leaves from a partition tree to super_tree
*/
void copySequencesToSuperTree(IntVector site_ids, int expected_num_states_super_tree, IQTree *super_tree, int initial_state, Node *node, Node *dad);

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, Alignment *aln, int sequence_length, int num_leaves, AliSimulator *alisimulator, bool inference_mode);

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(Alignment *aln, int sequence_length, string &output_str, int max_str_length, ostream &out, vector<string> state_mapping, InputType output_format, int max_length_taxa_name, Node *node, Node *dad);

/**
*  write a sequence of a node to an output file with gaps copied from the input sequence
*/
void writeASequenceToFileWithGaps(Alignment *aln, int sequence_length, vector<string> seq_names, vector<string> sequences, string &output_str, int max_str_length, ostream &out, vector<string> state_mapping, InputType output_format, int max_length_taxa_name, Node *node, Node *dad);

/**
*  merge and write all sequences to output files
*/
void mergeAndWriteSequencesToFiles(string file_path, AliSimulator *alisimulator, bool inference_mode);

/**
*  clear out all sequences in the super_tree
*
*/
void clearoutSequencesSuperTree(Node *node, Node *dad);

#endif /* alisim_h */
