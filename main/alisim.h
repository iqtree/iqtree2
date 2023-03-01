/*
 *  alisim.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#ifndef alisim_h
#define alisim_h

#include "utils/tools.h"
#include "utils/MPIHelper.h"
#include "simulator/alisimulatorinvar.h"
#include "simulator/alisimulatorheterogeneityinvar.h"
#include "phyloanalysis.h"
#include "tree/phylosupertree.h"
#include "utils/gzstream.h"
#include <regex>
#include <string.h>
#include "phyloanalysis.h"

/**
*  execute Alignment Simulator (AliSim)
*/
void runAliSim(Params &params, Checkpoint *checkpoint);

/**
*  execute AliSim Simulation
*/
void executeSimulation(Params params, IQTree *&tree);

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
void retrieveAncestralSequenceFromInputFile(AliSimulator *super_alisimulator, vector<short int> &sequence);

/**
*  generate mutiple alignments from a tree (model, alignment instances are supplied via the IQTree instance)
*/
void generateMultipleAlignmentsFromSingleTree(AliSimulator *super_alisimulator, map<string,string> input_msa);

/**
*  generate a partition alignment from a single simulator
*/
void generatePartitionAlignmentFromSingleSimulator(AliSimulator *&alisimulator, vector<short int> &ancestral_sequence, map<string,string> input_msa, string output_filepath = "", std::ios_base::openmode open_mode = std::ios_base::out);

/**
*  compute the total sequence length of all partitions
*/
int computeTotalSequenceLengthAllPartitions(PhyloSuperTree *super_tree);

/**
*  copy sequences of leaves from a partition tree to super_tree
*/
void copySequencesToSuperTree(IntVector &site_ids, int expected_num_states_super_tree, IQTree *current_tree, int initial_state, Node *node, Node *dad);

/**
*  write all sequences of a tree to an output file
*/
void writeSequencesToFile(string file_path, Alignment *aln, int sequence_length, int num_leaves, AliSimulator *alisimulator);

/**
*  write a sequence of a node to an output file
*/
void writeASequenceToFile(Alignment *aln, int sequence_length, int num_threads, bool keep_seq_order, uint64_t start_pos, uint64_t output_line_length,ostream &out, ostream &out_indels, bool write_indels_output, vector<string> &state_mapping, InputType output_format, int max_length_taxa_name, bool write_sequences_from_tmp_data, Node *node, Node *dad);

/**
*  merge and write all sequences to output files
*/
void mergeAndWriteSequencesToFiles(string file_path, AliSimulator *alisimulator, std::ios_base::openmode open_mode = std::ios_base::out);

/**
*  clear out all sequences in the super_tree
*
*/
void clearoutSequencesSuperTree(Node *node, Node *dad);

/**
*  load input MSA if the user wants to copy gaps from the input MSA
*
*/
map<string,string> loadInputMSA(AliSimulator *alisimulator);

/**
*  only unroot tree and stop if the user wants to do so
*
*/
void unrootTree(AliSimulator *alisimulator);

/**
*  determine real sequence length (for Indels)
*
*/
void determineSequenceLength(Node *node, Node *dad, bool &stop, int &sequence_length);

/**
*  insert redundant sites (inserted sites due to Indels) to the sequences of the super tree
*/
void insertIndelSites(int position, int starting_index, int num_inserted_sites, IQTree *current_tree, Node *node, Node *dad);

/**
*  write sequences to output file from a tmp_data and genome trees => a special case: with Indels without FunDi/ASC/Partitions
*/
void writeSeqsFromTmpDataAndGenomeTreesIndels(AliSimulator* alisimulator, int sequence_length, ostream &out, ostream &out_indels, bool write_indels_output, vector<string> &state_mapping, InputType output_format, int max_length_taxa_name);

#endif /* alisim_h */
