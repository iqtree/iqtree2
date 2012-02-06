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
	SuperNeighbor *nei;
	SuperNeighbor *dad_nei;
	if (dad) {
		nei = (SuperNeighbor*)node->findNeighbor(dad);
		dad_nei = (SuperNeighbor*)dad->findNeighbor(node);
		if (nei->link_neighbors.empty()) nei->link_neighbors.resize(size());
		if (dad_nei->link_neighbors.empty()) dad_nei->link_neighbors.resize(size());
		nei->link_neighbors[part] = NULL;
		dad_nei->link_neighbors[part] = NULL;
	}
	if (node->isLeaf()) {
		assert(dad);
		PhyloNode *node_part = (PhyloNode*)part_taxa[node->id];
		PhyloNode *dad_part = (PhyloNode*)node_part->neighbors[0]->node;
		
		if (node_part) {
			assert(node_part->isLeaf());
			nei->link_neighbors[part] = (PhyloNeighbor*) node_part->neighbors[0];
			dad_nei->link_neighbors[part] = (PhyloNeighbor*)dad_part->findNeighbor(node_part);
		} 
		return;
	}

	vector<PhyloNeighbor*> part_vec;
	vector<PhyloNeighbor*> child_part_vec;
	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		linkTree(part, part_taxa, (SuperNode*) (*it)->node, (SuperNode*) node);
		if (((SuperNeighbor*)*it)->link_neighbors[part]) {
			part_vec.push_back(((SuperNeighbor*)*it)->link_neighbors[part]);
			child_part_vec.push_back(((SuperNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]);
			assert(child_part_vec.back()->node == child_part_vec.front()->node);
		}
	}
	if (!dad || part_vec.empty()) return;

	if (part_vec.size() == 1) {
		nei->link_neighbors[part] = child_part_vec[0];
		dad_nei->link_neighbors[part] = part_vec[0];
		return;
	}
	PhyloNode *node_part = (PhyloNode*) child_part_vec[0]->node;
	PhyloNode *dad_part = NULL;
	FOR_NEIGHBOR(node_part, NULL, it) {
		bool appear = false;
		for (vector<PhyloNeighbor*>::iterator it2 = part_vec.begin(); it2 != part_vec.end(); it2++)
			if ((*it2) == (*it)) { appear = true; break; }
		if (!appear) {
			assert(!dad_part);
			dad_part = (PhyloNode*)(*it)->node;
		}
	}
	nei->link_neighbors[part] = (PhyloNeighbor*)node_part->findNeighbor(dad_part);
	dad_nei->link_neighbors[part] = (PhyloNeighbor*)dad_part->findNeighbor(node_part);
}

void PhyloSuperTree::mapTrees() {
	assert(root);
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		string taxa_set = ((SuperAlignment*)aln)->getPattern(part);
		(*it)->copyTree(this, taxa_set);
		//(*it)->drawTree(cout);
		(*it)->initializeAllPartialLh();
		NodeVector my_taxa;
		(*it)->getTaxa(my_taxa);
		NodeVector part_taxa;
		part_taxa.resize(leafNum, NULL);
		for (int i = 0; i < leafNum; i++)
			part_taxa[i] = my_taxa[((SuperAlignment*)aln)->taxa_index[i][part]];
		linkTree(part, part_taxa);
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


