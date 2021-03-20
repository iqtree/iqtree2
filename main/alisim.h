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

void runAliSim(Params params);
void showParameters(Params params);
IQTree *initializeIQTreeFromTreeFile(Params params, char* seq_type);
void initializeAlignment(char* seq_type, IQTree *tree);
void addLeafNamesToAlignment(Alignment *aln, Node *node, Node *dad);
void initializeModel(Params params, IQTree *tree);
void generateSingleDatasetFromSingleTree(Params params, IQTree *tree, string output_filepath);
IntVector getAncestralSequence(Params params, IQTree *tree);
IntVector retrieveAncestralSequenceFromInputFile(int sequence_position, IQTree *tree);
IntVector generateRandomSequence(int sequence_length, IQTree *tree);
void simulateSequencesForATreeByDFS(int sequence_length, ModelSubst *model, int max_num_states, Node *node, Node *dad);
IntVector estimateSequenceOfANode(double *trans_matrix, int max_num_states, int sequence_length, IntVector dad_sites);
int getRandomStateWithProbabilityMatrix(double *probability_maxtrix, int starting_index, int max_num_states);
void writeSequencesToFile(string file_path, IQTree *tree, int sequence_length);
void writeASequenceToFile(Alignment *aln, ofstream &out, Node *node, Node *dad);
string convertEncodedSequenceToReadableSequence(Alignment *aln, IntVector sequence);


#endif /* alisim_h */
