/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PHYLOANALYSIS_H
#define PHYLOANALYSIS_H

#include "utils/tools.h"
#include "tree/mexttree.h"
#include "phylotesting.h"
#include "treetesting.h"
#include "tree/upperbounds.h" // Olga: functions for Upper Bounds analysis
#include "utils/pllnni.h"

class PhyloTree;
class IQTree;

/**
	main function to carry out phylogenetic inference
	@param params program parameters
*/
void runPhyloAnalysis(Params &params, Checkpoint *checkpoint);

/**
    carry out phylogenetic inference without deleting IQTree instance
    @param params program parameters
*/
void runPhyloAnalysis(Params &params, Checkpoint *checkpoint, IQTree *&tree, Alignment *&aln);

/**
    Perform separate tree inference across partitions
 */
void runUnlinkedPhyloAnalysis(Params &params, Checkpoint *checkpoint);

void startTreeReconstruction(Params &params, IQTree* &iqtree,
        ModelCheckpoint &model_info);

void runTreeReconstruction(Params &params, IQTree* &tree);

/**
	take the collection of trees from input_trees, it assign support values to target_tree
	and print resulting tree to output_tree. 
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param target_tree tree to assign support value
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting will be written to this file. If NULL,
		output_tree will be named target_tree appended with ".suptree"
*/
void assignBootstrapSupport(const char *input_trees, int burnin, int max_count, const char *target_tree, 
	bool rooted, const char *output_tree, const char *out_prefix, MExtTree &mytree, 
	const char* tree_weight_file, Params *params);

/**
 * assign branch supports from params.user_tree trees file to params.second_tree
 * @param params program parameters
 */
void assignBranchSupportNew(Params &params);

/**
	Compute the consensus tree from the collection of trees from input_trees
	and print resulting tree to output_tree. 
	@param phylo_tree used to optimize branch lengths of the consensus tree. Can be NULL
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param cutoff only incorporate those splits that have support values more than cutoff
	@param weight_threshold minimum weight cutoff
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting consensus tree will be written to this file. If NULL,
		output_tree will be named input_trees appended with ".contree"
*/
void computeConsensusTree(const char *input_trees, int burnin, int max_count, double cutoff, double weight_threshold,
	const char *output_tree, const char *out_prefix, const char* tree_weight_file, Params *params);

/**
	Compute the consensus network from the collection of trees in input_trees.
	print consensus network to output_tree
	@param input_trees collection of input trees to infer split supports
	@param burnin the number trees at the beginning of input_trees to be discarded
	@param max_count max number of trees to load
	@param cutoff only incorporate those splits that have support values more than cutoff
	@param weight_threshold minimum weight cutoff
	@param output_tree (OUT, OVERWRITE IF EXIST) Resulting consensus tree will be written to this file. If NULL,
		output_tree will be named input_trees appended with ".connetwork"
*/
void computeConsensusNetwork(const char *input_trees, int burnin, int max_count, double cutoff,
		int weight_summary, double weight_threshold,
	const char *output_tree, const char *out_prefix, const char* tree_weight_file);

void reportRate(ostream &out, PhyloTree &tree);

void reportSubstitutionProcess(ostream &out, Params &params, IQTree &tree);

void exportAliSimCMD(Params &params, IQTree &tree, ostream &out);

/** compute rootstrap for a user defined tree from a set of trees */
void runRootstrap(Params &params);

#endif
