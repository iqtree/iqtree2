/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
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
#include "phylosupertree.h"
#include "superalignment.h"
#include "superalignmentpairwise.h"

PhyloSuperTree::PhyloSuperTree()
 : IQPTree()
{
}

PhyloSuperTree::PhyloSuperTree(Params &params) :  IQPTree() {
	cout << "Reading partition model file " << params.partition_file << " ..." << endl;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(params.partition_file);
		in.exceptions(ios::badbit);

		Params origin_params = params;
		//memcpy(&origin_params, &params, sizeof(params));
		PartitionInfo info;

		while (!in.eof()) {
			getline(in, info.name, ',');
			if (in.eof()) break;
			getline(in, info.model_name, ',');
			if (info.model_name == "") info.model_name = params.model_name;
			getline(in, info.aln_file, ',');
			if (info.aln_file == "" && params.aln_file) info.aln_file = params.aln_file;
			getline(in, info.sequence_type, ',');
			if (info.sequence_type=="" && params.sequence_type) info.sequence_type = params.sequence_type;
			getline(in, info.position_spec);
			cout << endl << "Reading partition " << info.name << " (model=" << info.model_name << ", aln=" << 
				info.aln_file << ", seq=" << info.sequence_type << ", pos=" << info.position_spec << ") ..." << endl;
			part_info.push_back(info);
			Alignment *part_aln = new Alignment((char*)info.aln_file.c_str(), (char*)info.sequence_type.c_str(), params.intype);
			PhyloTree *tree = new PhyloTree(part_aln);
			push_back(tree);
			params = origin_params;
		}

		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (string str) {
		outError(str);
	}
	aln = new SuperAlignment(this);
	string str = params.out_prefix;
	str += ".part";
	aln->printPhylip(str.c_str());
	str = params.out_prefix;
	str += ".concat";
	((SuperAlignment*)aln)->printCombinedAlignment(str.c_str());
	cout << endl;

}

Node* PhyloSuperTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

Node* PhyloSuperTree::newNode(int node_id, int node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

double PhyloSuperTree::computeDist(int seq1, int seq2, double initial_dist) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0)
        initial_dist = aln->computeDist(seq1, seq2);
    if (initial_dist == MAX_GENETIC_DIST) return initial_dist;
    if (!model_factory || !site_rate) return initial_dist;

    // now optimize the distance based on the model and site rate
    SuperAlignmentPairwise aln_pair(this, seq1, seq2);
    return aln_pair.optimizeDist(initial_dist);
}

void PhyloSuperTree::linkTree(int part, NodeVector &part_taxa, SuperNode *node, SuperNode *dad) {
	if (!node) {
		if (!root->isLeaf()) 
			node = (SuperNode*) root; 
		else 
			node = (SuperNode*)root->neighbors[0]->node;
		if (node->isLeaf()) // two-taxa tree
			dad = (SuperNode*)node->neighbors[0]->node;
	}
	SuperAlignment *aln = (SuperAlignment*) (this->aln);
	if (node->isLeaf() && dad) {
		PhyloNode *node_part = (PhyloNode*)part_taxa[node->id];
		assert(node_part->isLeaf());
		SuperNeighbor *nei = (SuperNeighbor*)(node->neighbors[0]);
		if (nei->link_neighbors.empty()) nei->link_neighbors.resize(size());
	}
	FOR_NEIGHBOR_IT(node, dad, it) {
		linkTree(part, part_taxa, (SuperNode*) (*it)->node, (SuperNode*) node);
	}
}

void PhyloSuperTree::mapTrees() {
	assert(root);
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		string taxa_set = ((SuperAlignment*)aln)->getPattern(part);
		(*it)->copyTree(this, taxa_set);
		//(*it)->drawTree(cout);
		(*it)->initializeAllPartialLh();
	}
}

double PhyloSuperTree::computeLikelihood(double *pattern_lh) {
	double tree_lh = 0.0;
	for (iterator it = begin(); it != end(); it++)
		tree_lh += (*it)->computeLikelihood();
	return tree_lh;
}

double PhyloSuperTree::optimizeAllBranches(int iterations) {
	double tree_lh = 0.0;
	for (iterator it = begin(); it != end(); it++)
		tree_lh += (*it)->optimizeAllBranches(iterations);
	return tree_lh;
}

PhyloSuperTree::~PhyloSuperTree()
{
	for (reverse_iterator it = rbegin(); it != rend(); it++)
		delete (*it);
	clear();
}


