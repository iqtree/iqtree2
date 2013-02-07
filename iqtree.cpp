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
#include "iqtree.h"
#include "phylosupertree.h"
#include "mexttree.h"
#include "timeutil.h"

extern double t_begin;

IQTree::IQTree() :
		PhyloTree() {
	init();
}

void IQTree::init() {
	k_represent = 0;
	k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
	dist_matrix = NULL;
	nni_count_est = 0.0;
	nni_delta_est = 0;
	curScore = 0.0; // Current score of the tree
	bestScore = 0.0; // Best score found sofar
	cur_pars_score = -1;
	enable_parsimony = false;
	enableHeuris = true; // This is set true when the heuristic started (after N iterations)
	linRegModel = NULL;
	nni_round = 1;
	estimate_nni_cutoff = false;
	nni_cutoff = -1e6;
	nni_sort = false;
	testNNI = false;
	save_all_trees = 0;
	print_tree_lh = false;
	write_intermediate_trees = 0;
	max_candidate_trees = 0;
	logl_cutoff = 0.0;
	len_scale = 10000;
	save_all_br_lens = false;
	duplication_counter = 0;
	raxmlTree = NULL;
	//boot_splits = new SplitGraph;
	//initLeafFrequency();
	nnicut.num_delta = 0;
	nnicut.delta_min = DBL_MAX;
//	if (params->tabu){
//		initLeafFrequency();
//	}

}

IQTree::IQTree(Alignment *aln) :
		PhyloTree(aln) {
	init();
}

double IQTree::computeParsimonyTreePhylolib() {
	allocateParsimonyDataStructures(raxmlTree);
	makeParsimonyTreeFast(raxmlTree);
	freeParsimonyDataStructures(raxmlTree);

	smoothTree(raxmlTree, 32);
	evaluateGeneric(raxmlTree, raxmlTree->start, FALSE);

	return raxmlTree->likelihood;

	/*
	 int printBranchLengths = TRUE;
	 Tree2String(iqtree.raxmlTree->tree_string, iqtree.raxmlTree,
	 iqtree.raxmlTree->start->back, printBranchLengths, TRUE, 0, 0, 0,
	 SUMMARIZE_LH, 0, 0);
	 string parsimony_tree_file = params.out_prefix;
	 parsimony_tree_file += ".parsimonyTree";
	 printString2File(string(iqtree.raxmlTree->tree_string),
	 parsimony_tree_file);
	 */
}

void IQTree::setParams(Params &params) {
	if (params.min_iterations == -1) {
		if (!params.gbo_replicates) {
			if (aln->getNSeq() < 100)
				params.min_iterations = 200;
			else
				params.min_iterations = aln->getNSeq() * 2;
			if (params.iteration_multiple > 1)
				params.min_iterations = aln->getNSeq()
						* params.iteration_multiple;
		} else {
			params.min_iterations = 100;
			params.max_iterations = 1000;
		}
	}
	k_represent = params.k_representative;

	if (params.p_delete == 0.0) {
		if (aln->getNSeq() < 51)
			params.p_delete = 0.5;
		else if (aln->getNSeq() < 100)
			params.p_delete = 0.3;
		else if (aln->getNSeq() < 200)
			params.p_delete = 0.2;
		else
			params.p_delete = 0.1;
	}
	//tree.setProbDelete(params.p_delete);
	if (params.p_delete != 0.0) {
		k_delete = k_delete_min = k_delete_max = round(
				params.p_delete * leafNum);
	} else {
		k_delete = k_delete_min = 10;
		k_delete_max = leafNum / 2;
		if (k_delete_max > 100)
			k_delete_max = 100;
		if (k_delete_max < 20)
			k_delete_max = 20;
		k_delete_stay = ceil(leafNum / k_delete);
	}

	//tree.setIQPIterations(params.stop_condition, params.stop_confidence, params.min_iterations, params.max_iterations);

	stop_rule.setStopCondition(params.stop_condition);
	stop_rule.setConfidenceValue(params.stop_confidence);
	stop_rule.setIterationNum(params.min_iterations, params.max_iterations);

	//tree.setIQPAssessQuartet(params.iqp_assess_quartet);
	iqp_assess_quartet = params.iqp_assess_quartet;
	estimate_nni_cutoff = params.estimate_nni_cutoff;
	nni_cutoff = params.nni_cutoff;
	nni_sort = params.nni_sort;
	testNNI = params.testNNI;

	this->params = &params;

	write_intermediate_trees = params.write_intermediate_trees;

	if (write_intermediate_trees > 2 || params.gbo_replicates > 0) {
		save_all_trees = 1;
	}
	if (params.gbo_replicates > 0) {
		if (params.iqp_assess_quartet != IQP_BOOTSTRAP) {
			save_all_trees = 2;
		}
	}
	if (params.gbo_replicates > 0 && params.do_compression)
		save_all_br_lens = true;
	print_tree_lh = params.print_tree_lh;
	max_candidate_trees = params.max_candidate_trees;
	if (max_candidate_trees == 0)
		max_candidate_trees = aln->getNSeq() * stop_rule.getNumIterations();
	setRootNode(params.root);

	if (params.online_bootstrap && params.gbo_replicates > 0) {
		cout << "Generating " << params.gbo_replicates
				<< " samples for ultra-fast bootstrap..." << endl;
		boot_samples.resize(params.gbo_replicates);
		boot_logl.resize(params.gbo_replicates, -DBL_MAX);
		boot_trees.resize(params.gbo_replicates, -1);
		for (int i = 0; i < params.gbo_replicates; i++) {
			aln->createBootstrapAlignment(boot_samples[i]);
		}
		cout << "Max candidate trees (tau): " << max_candidate_trees << endl;
	}

	nnicut.doNNICut = params.estimate_nni_cutoff;
	if (params.root_state) {
		if (strlen(params.root_state) != 1) 
			outError("Root state must have exactly 1 character");
		root_state = aln->convertState(params.root_state[0]);
		if (root_state < 0 || root_state >= aln->num_states) 
			outError("Invalid root state");
	}
}

IQTree::~IQTree() {
	//if (bonus_values)
	//delete bonus_values;
	//bonus_values = NULL;
	if (dist_matrix)
		delete[] dist_matrix;
	dist_matrix = NULL;

	for (vector<double*>::reverse_iterator it = treels_ptnlh.rbegin();
			it != treels_ptnlh.rend(); it++)
		delete[] (*it);
	treels_ptnlh.clear();
	for (vector<SplitGraph*>::reverse_iterator it2 = boot_splits.rbegin();
			it2 != boot_splits.rend(); it2++)
		delete (*it2);
	//if (boot_splits) delete boot_splits;
	if (raxmlTree)
		delete raxmlTree;
}

double IQTree::getProbDelete() {
	return (double) k_delete / leafNum;
}

void IQTree::resetKDelete() {
	k_delete = k_delete_min;
	k_delete_stay = ceil(leafNum / k_delete);
}

void IQTree::increaseKDelete() {
	if (k_delete >= k_delete_max)
		return;
	k_delete_stay--;
	if (k_delete_stay > 0)
		return;
	k_delete++;
	k_delete_stay = ceil(leafNum / k_delete);
	if (verbose_mode >= VB_MED)
		cout << "Increase k_delete to " << k_delete << endl;
}

void IQTree::setIQPIterations(STOP_CONDITION stop_condition,
		double stop_confidence, int min_iterations, int max_iterations) {
	stop_rule.setStopCondition(stop_condition);
	stop_rule.setConfidenceValue(stop_confidence);
	stop_rule.setIterationNum(min_iterations, max_iterations);
}

RepresentLeafSet* IQTree::findRepresentLeaves(
		vector<RepresentLeafSet*> &leaves_vec, int nei_id, PhyloNode *dad) {
	PhyloNode *node = (PhyloNode*) (dad->neighbors[nei_id]->node);
	int set_id = dad->id * 3 + nei_id;
	if (leaves_vec[set_id])
		return leaves_vec[set_id];
	RepresentLeafSet *leaves = new RepresentLeafSet;
	RepresentLeafSet * leaves_it[2] = { NULL, NULL };
	leaves_vec[set_id] = leaves;
	RepresentLeafSet::iterator last;
	RepresentLeafSet::iterator cur_it;
	int i, j;
	//double admit_height = 1000000;

	leaves->clear();
	if (node->isLeaf()) {
		// set the depth to zero
		//node->height = 0.0;
		leaves->insert(new RepLeaf(node, 0));
	} else {
		for (i = 0, j = 0; i < node->neighbors.size(); i++)
			if (node->neighbors[i]->node != dad) {
				leaves_it[j++] = findRepresentLeaves(leaves_vec, i, node);
			}
		assert(j == 2 && leaves_it[0] && leaves_it[1]);
		if (leaves_it[0]->empty() && leaves_it[1]->empty()) {
			cout << "wrong";
		}
		RepresentLeafSet::iterator lit[2] = { leaves_it[0]->begin(),
				leaves_it[1]->begin() };
		while (leaves->size() < k_represent) {
			int id = -1;
			if (lit[0] != leaves_it[0]->end()
					&& lit[1] != leaves_it[1]->end()) {
				if ((*lit[0])->height < (*lit[1])->height)
					id = 0;
				else if ((*lit[0])->height > (*lit[1])->height)
					id = 1;
				else { // tie, choose at random
					id = random_int(2);
				}
			} else if (lit[0] != leaves_it[0]->end())
				id = 0;
			else if (lit[1] != leaves_it[1]->end())
				id = 1;
			else
				break;
			assert(id < 2 && id >= 0);
			leaves->insert(
					new RepLeaf((*lit[id])->leaf, (*lit[id])->height + 1));
			lit[id]++;
		}
	}
	assert(!leaves->empty());
	/*
	 if (verbose_mode >= VB_MAX) {
	 for (cur_it = leaves->begin(); cur_it != leaves->end(); cur_it++)
	 cout << (*cur_it)->leaf->name << " ";
	 cout << endl;
	 }*/
	return leaves;
}

/*
 void IQPTree::clearRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, Node *node, Node *dad) {
 int nei_id;
 for (nei_id = 0; nei_id < node->neighbors.size(); nei_id++)
 if (node->neighbors[nei_id]->node == dad) break;
 assert(nei_id < node->neighbors.size());
 int set_id = node->id * 3 + nei_id;
 if (leaves_vec[set_id]) {
 for (RepresentLeafSet::iterator rlit = leaves_vec[set_id]->begin(); rlit != leaves_vec[set_id]->end(); rlit++)
 delete (*rlit);
 delete leaves_vec[set_id];
 leaves_vec[set_id] = NULL;
 }
 FOR_NEIGHBOR_IT(node, dad, it) {
 clearRepresentLeaves(leaves_vec, (*it)->node, node);
 }
 }*/

void IQTree::initLeafFrequency(PhyloNode *node, PhyloNode *dad) {
    if (!node) node = (PhyloNode*) root;
    if (node->isLeaf()) {
    	LeafFreq leaf_freq;
    	leaf_freq.leaf = node;
    	leaf_freq.freq = 0;
    	leaf_freqs.push_back(leaf_freq);
    }
	for (NeighborVec::iterator it = node->neighbors.begin();
			it != node->neighbors.end(); it++)
		if ((*it)->node != dad) {
			initLeafFrequency((PhyloNode*) (*it)->node, node);
		}
}

void IQTree::deleteNonTabuLeaves(PhyloNodeVector &del_leaves) {
	// sort node frequency
	sort(leaf_freqs.begin(), leaf_freqs.end());
	// shuffle the sorted list first
	int cur_freq = leaf_freqs[0].freq;
	int startIndex = 0;
	int endIndex = 0;
	int leafCnt = 0;
	while (true) {
		while (leaf_freqs[endIndex].freq == cur_freq) {
			endIndex++;
		}
		random_shuffle ( leaf_freqs.begin() + startIndex, leaf_freqs.begin() + endIndex - 1 );
		do {
			del_leaves.push_back( leaf_freqs[leafCnt].leaf );
			leaf_freqs[leafCnt].freq++;
			leafCnt++;
		} while (leafCnt < k_delete);
		if (leafCnt == k_delete) {
			break;
		}
		startIndex = endIndex;
		cur_freq = leaf_freqs[endIndex].freq;
	}
}

void IQTree::deleteLeaves(PhyloNodeVector &del_leaves) {
	NodeVector taxa;
	// get the vector of taxa
	getTaxa(taxa);
	root = NULL;
	//int num_delete = floor(p_delete * taxa.size());
	int num_delete = k_delete;
	int i;
	if (num_delete > taxa.size() - 4)
		num_delete = taxa.size() - 4;
	// now try to randomly delete some taxa of the probability of p_delete
	for (i = 0; i < num_delete;) {
		int id = random_int(taxa.size());
		if (!taxa[id])
			continue;
		else
			i++;
		PhyloNode *taxon = (PhyloNode*) taxa[id];
		del_leaves.push_back(taxon);
		deleteLeaf(taxon);
		taxa[id] = NULL;
	}
	// set root to the first taxon which was not deleted
	for (i = 0; i < taxa.size(); i++)
		if (taxa[i]) {
			root = taxa[i];
			break;
		}
}

int IQTree::assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2,
		Node *del_leaf) {
	assert(dist_matrix);
	int nseq = aln->getNSeq();
	//int id0 = leaf0->id, id1 = leaf1->id, id2 = leaf2->id;
	double dist0 = dist_matrix[leaf0->id * nseq + del_leaf->id]
			+ dist_matrix[leaf1->id * nseq + leaf2->id];
	double dist1 = dist_matrix[leaf1->id * nseq + del_leaf->id]
			+ dist_matrix[leaf0->id * nseq + leaf2->id];
	double dist2 = dist_matrix[leaf2->id * nseq + del_leaf->id]
			+ dist_matrix[leaf0->id * nseq + leaf1->id];
	if (dist0 < dist1 && dist0 < dist2)
		return 0;
	if (dist1 < dist2)
		return 1;
	return 2;
}

int IQTree::assessQuartetParsimony(Node *leaf0, Node *leaf1, Node *leaf2,
		Node *del_leaf) {
	int score[3] = { 0, 0, 0 };
	for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
		char ch0 = (*it)[leaf0->id];
		char ch1 = (*it)[leaf1->id];
		char ch2 = (*it)[leaf2->id];
		char chd = (*it)[del_leaf->id];
		if (ch0 >= aln->num_states || ch1 >= aln->num_states
				|| ch2 >= aln->num_states || chd >= aln->num_states)
			continue;
		if (chd == ch0 && ch1 == ch2)
			score[0] += (*it).frequency;
		if (chd == ch1 && ch0 == ch2)
			score[1] += (*it).frequency;
		if (chd == ch2 && ch0 == ch1)
			score[2] += (*it).frequency;
	}
	if (score[0] == score[1] && score[0] == score[2]) {
		int id = random_int(3);
		return id;
	}
	if (score[0] > score[1] && score[0] > score[2])
		return 0;
	if (score[1] < score[2])
		return 2;
	return 1;
}

void IQTree::initializeBonus(PhyloNode *node, PhyloNode *dad) {
	if (!node)
		node = (PhyloNode*) root;
	if (dad) {
		PhyloNeighbor *node_nei = (PhyloNeighbor*) node->findNeighbor(dad);
		PhyloNeighbor *dad_nei = (PhyloNeighbor*) dad->findNeighbor(node);
		node_nei->lh_scale_factor = 0.0;
		node_nei->partial_lh_computed = 0;
		dad_nei->lh_scale_factor = 0.0;
		dad_nei->partial_lh_computed = 0;
	}

	FOR_NEIGHBOR_IT(node, dad, it){
	initializeBonus((PhyloNode*) ((*it)->node), node);
}
}

void IQTree::raiseBonus(Neighbor *nei, Node *dad, double bonus) {
	((PhyloNeighbor*) nei)->lh_scale_factor += bonus;
	if (verbose_mode >= VB_DEBUG)
		cout << dad->id << " - " << nei->node->id << " : " << bonus << endl;

	//  FOR_NEIGHBOR_IT(nei->node, dad, it)
	//	raiseBonus((*it), nei->node, bonus);
}

double IQTree::computePartialBonus(Node *node, Node* dad) {
	PhyloNeighbor *node_nei = (PhyloNeighbor*) node->findNeighbor(dad);
	if (node_nei->partial_lh_computed)
		return node_nei->lh_scale_factor;

	FOR_NEIGHBOR_IT(node, dad, it){
	node_nei->lh_scale_factor += computePartialBonus((*it)->node, node);
}
	node_nei->partial_lh_computed = 1;
	return node_nei->lh_scale_factor;
}

double IQTree::doVNS() {
	curScore = optimizeNNI(false);
	return curScore;
}

bool IQTree::containPosNNI(vector<NNIMove> posNNIs) {
	for (vector<NNIMove>::iterator iter = posNNIs.begin();
			iter != posNNIs.end(); iter++) {
		if (iter->delta > 0)
			return true;
	}
	return false;
}

void IQTree::findBestBonus(double &best_score, NodeVector &best_nodes,
		NodeVector &best_dads, Node *node, Node *dad) {
	double score;
	if (!node)
		node = root;
	if (!dad) {
		best_score = 0;
	} else {
		score = computePartialBonus(node, dad) + computePartialBonus(dad, node);
		if (score >= best_score) {
			if (score > best_score) {
				best_score = score;
				best_nodes.clear();
				best_dads.clear();
			}
			best_nodes.push_back(node);
			best_dads.push_back(dad);
		}
		//cout << node->id << " - " << dad->id << " : " << best_score << endl;
	}

	FOR_NEIGHBOR_IT(node, dad, it){
	findBestBonus(best_score, best_nodes, best_dads, (*it)->node, node);
}
}

void IQTree::assessQuartets(vector<RepresentLeafSet*> &leaves_vec,
		PhyloNode *cur_root, PhyloNode *del_leaf) {
	const int MAX_DEGREE = 3;
	RepresentLeafSet * leaves[MAX_DEGREE];
	double bonus[MAX_DEGREE];
	memset(bonus, 0, MAX_DEGREE * sizeof(double));
	int cnt = 0;

	// only work for birfucating tree
	assert(cur_root->degree() == MAX_DEGREE);

	// find the representative leaf set for three subtrees

	FOR_NEIGHBOR_IT(cur_root, NULL, it){
	leaves[cnt] = findRepresentLeaves(leaves_vec, cnt, cur_root);
	cnt++;
}
	for (RepresentLeafSet::iterator i0 = leaves[0]->begin();
			i0 != leaves[0]->end(); i0++)
		for (RepresentLeafSet::iterator i1 = leaves[1]->begin();
				i1 != leaves[1]->end(); i1++)
			for (RepresentLeafSet::iterator i2 = leaves[2]->begin();
					i2 != leaves[2]->end(); i2++) {
				int best_id;
				if (iqp_assess_quartet == IQP_DISTANCE)
					best_id = assessQuartet((*i0)->leaf, (*i1)->leaf,
							(*i2)->leaf, del_leaf);
				else
					best_id = assessQuartetParsimony((*i0)->leaf, (*i1)->leaf,
							(*i2)->leaf, del_leaf);
				bonus[best_id] += 1.0;
			}
	for (cnt = 0; cnt < MAX_DEGREE; cnt++)
		if (bonus[cnt] > 0.0)
			raiseBonus(cur_root->neighbors[cnt], cur_root, bonus[cnt]);

}

void IQTree::reinsertLeavesByParsimony(PhyloNodeVector &del_leaves) {
	PhyloNodeVector::iterator it_leaf;
	assert(root->isLeaf());
	initializeAllPartialPars();
	clearAllPartialLH();
	for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
		initializeAllPartialPars();
		clearAllPartialLH();
        Node *target_node = NULL;
        Node *target_dad = NULL;
	}

}


void IQTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
	PhyloNodeVector::iterator it_leaf;

	//int num_del_leaves = del_leaves.size();
	assert(root->isLeaf());

	for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
		if (verbose_mode >= VB_DEBUG)
			cout << "Reinserting " << (*it_leaf)->name << " (" << (*it_leaf)->id
					<< ")" << endl;
		vector<RepresentLeafSet*> leaves_vec;
		leaves_vec.resize(nodeNum * 3, NULL);
		initializeBonus();
		NodeVector nodes;
		getInternalNodes(nodes);
		if (verbose_mode >= VB_DEBUG)
			drawTree(cout,
					WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE
							| WT_BR_ID);
		//printTree(cout, WT_BR_LEN | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
		for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++) {
			assessQuartets(leaves_vec, (PhyloNode*) (*it), (*it_leaf));
		}
		NodeVector best_nodes, best_dads;
		double best_bonus;
		findBestBonus(best_bonus, best_nodes, best_dads);
		if (verbose_mode >= VB_DEBUG)
			cout << "Best bonus " << best_bonus << " " << best_nodes[0]->id
					<< " " << best_dads[0]->id << endl;
		assert(best_nodes.size() == best_dads.size());
		int node_id = random_int(best_nodes.size());
		if (best_nodes.size() > 1 && verbose_mode >= VB_DEBUG)
			cout << best_nodes.size()
					<< " branches show the same best bonus, branch nr. "
					<< node_id << " is chosen" << endl;

		reinsertLeaf(*it_leaf, best_nodes[node_id], best_dads[node_id]);
		//clearRepresentLeaves(leaves_vec, *it_node, *it_leaf);
		/*if (verbose_mode >= VB_DEBUG) {
		 printTree(cout);
		 cout << endl;
		 }*/
		for (vector<RepresentLeafSet*>::iterator rit = leaves_vec.begin();
				rit != leaves_vec.end(); rit++)
			if ((*rit)) {
				RepresentLeafSet *tit = (*rit);
				for (RepresentLeafSet::iterator rlit = tit->begin();
						rlit != tit->end(); rlit++)
					delete (*rlit);
				delete (*rit);
			}
	}
	initializeTree(); // BQM: re-index nodes and branches s.t. ping-pong neighbors have the same ID

	if (verbose_mode >= VB_DEBUG)
		drawTree(cout,
				WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE
						| WT_BR_ID);
}

double IQTree::doParsimonyReinsertion() {
	PhyloNodeVector del_leaves;
	if (params->tabu) {
		deleteNonTabuLeaves(del_leaves);
	} else {
		deleteLeaves(del_leaves);
	}
	reinsertLeavesByParsimony(del_leaves);
}


double IQTree::doIQP() {
	if (verbose_mode >= VB_DEBUG)
		drawTree(cout,
				WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE
						| WT_BR_ID);
	//double time_begin = getCPUTime();
	PhyloNodeVector del_leaves;

	if (params->tabu) {
		deleteNonTabuLeaves(del_leaves);
	} else {
		deleteLeaves(del_leaves);
	}
	/*
	for (int i = 0; i < del_leaves.size(); i++)
		cout << del_leaves[i]->id << " ";
	cout << endl;
	*/
	reinsertLeaves(del_leaves);

	if (params->raxmllib) {
		// otherwise branch lengths will be optimized using phylolib's kernel
		curScore = 0.00;
	} else {
		// just to make sure IQP does it right
		setAlignment(aln);
		clearAllPartialLH();
		if (params->gbo_replicates)
			curScore = optimizeAllBranches(3, 1.0);
		else {

			// optimize branches at the reinsertion point
			/*
			 for (PhyloNodeVector::iterator dit = del_leaves.begin(); dit != del_leaves.end(); dit++) {
			 PhyloNode *adj_node = (PhyloNode*)(*dit)->neighbors[0]->node;
			 FOR_NEIGHBOR_IT(adj_node, (*dit), it)
			 curScore = optimizeOneBranch(adj_node, (PhyloNode*)(*it)->node);
			 //curScore = optimizeOneBranch(adj_node, (PhyloNode*)(*dit));
			 }
			 */
			curScore = optimizeAllBranches(1);
		}

	}

	//cout << "IQP log-likelihood: " << curScore << endl;
	//cout << "computeLikelihood(): " << computeLikelihood() << endl;

	if (enable_parsimony) {
		cur_pars_score = computeParsimony();
		if (verbose_mode >= VB_MAX) {
			cout << "IQP Likelihood = " << curScore << "  Parsimony = "
					<< cur_pars_score << endl;
		}
	}
	//double time_end = getCPUTime();
	//cout << "IQP Time: " << (time_end - time_begin) << endl;
	return curScore;
}

double IQTree::swapTaxa(PhyloNode *node1, PhyloNode *node2) {
	assert( node1->isLeaf());
	assert( node2->isLeaf());

	PhyloNeighbor *node1nei = (PhyloNeighbor*) *(node1->neighbors.begin());
	PhyloNeighbor *node2nei = (PhyloNeighbor*) *(node2->neighbors.begin());

	node2nei->node->updateNeighbor(node2, node1);
	node1nei->node->updateNeighbor(node1, node2);

	// Update the new neightbors of the 2 nodes
	node1->updateNeighbor(node1->neighbors.begin(), node2nei);
	node2->updateNeighbor(node2->neighbors.begin(), node1nei);

	PhyloNeighbor *node1NewNei = (PhyloNeighbor*) *(node1->neighbors.begin());
	PhyloNeighbor *node2NewNei = (PhyloNeighbor*) *(node2->neighbors.begin());

	// Reoptimize the branch lengths
	optimizeOneBranch(node1, (PhyloNode*) node1NewNei->node);
	this->curScore = optimizeOneBranch(node2, (PhyloNode*) node2NewNei->node);
	//drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
	return this->curScore;
}

double IQTree::perturb(int times) {
	while (times > 0) {
		NodeVector taxa;
		// get the vector of taxa
		getTaxa(taxa);
		int taxonid1 = random_int(taxa.size());
		PhyloNode *taxon1 = (PhyloNode*) taxa[taxonid1];
		PhyloNode *taxon2;
		int *dists = new int[taxa.size()];
		int minDist = 1000000;
		for (int i = 0; i < taxa.size(); i++) {
			if (i == taxonid1)
				continue;
			taxon2 = (PhyloNode*) taxa[i];
			int dist = taxon1->calDist(taxon2);
			dists[i] = dist;
			if (dist >= 7 && dist < minDist)
				minDist = dist;
		}

		int taxonid2 = -1;
		for (int i = 0; i < taxa.size(); i++) {
			if (dists[i] == minDist)
				taxonid2 = i;
		}

		taxon2 = (PhyloNode*) taxa[taxonid2];

		cout << "Swapping node " << taxon1->id << " and node " << taxon2->id
				<< endl;
		cout << "Distance " << minDist << endl;
		curScore = swapTaxa(taxon1, taxon2);
		//taxa.erase( taxa.begin() + taxaID1 );
		//taxa.erase( taxa.begin() + taxaID2 -1 );

		times--;
		delete[] dists;
	}
	curScore = optimizeAllBranches(1);
	return curScore;
}

//double IQPTree::doILS(Params &params, int perturbLevel) {
//
//    string tree_file_name = params.aln_file;
//    tree_file_name += ".treefile";
//    // keep the best tree into a string
//    stringstream best_tree_string;
//    printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
//    bestScore = curScore;PhyloNode* node1 = vec_nonconf_nni.at(i).node1;
//
//    int numIter = params.min_iterations;
//    for (int i=1 ; i <= numIter; i++) {
//
//        if (i > speedUpFromIter) {
//            enableHeuris = true;
//            nbNNI95 = estimateNumNNI();
//            deltaNNI95 = estimateDeltaNNI();
//        }
//
//        cout.precision(10);
//        perturb(perturbLevel);
//
//        optimizeNNI();
//
//        cout.precision(15);
//        cout << "Iteration " << i << " / Log-Likelihood: "
//                << curScore << endl;
//        if (curScore > bestScore + TOL_LIKELIHOOD) {
//            //nni_score = optimizeNNI(true);
//            //curScore = optimizeAllBranches();
//            cout << "BETTER TREE FOUND: " << curScore << endl;
//            bestScore = curScore;
//            best_tree_string.seekp(0);
//            printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
//            printTree(tree_file_name.c_str());
//        } else {
//            /* take back the current best tree */
//            best_tree_string.seekg(0);
//            freeNode();
//            readTree(best_tree_string, rooted);
//            assignLeafNames();
//            initializeAllPartialLh();
//        }
//    }
//
//    return bestScore;
//}


double IQTree::doRandomRestart() {
	int cur_iteration;
	bool maxTimeReached = false;
	setRootNode(params->root);
	// keep the best tree into a string
	stringstream best_tree_string;
	printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
	stringstream best_tree_topo_ss;
	printTree(best_tree_topo_ss, WT_TAXON_ID + WT_SORT_TAXA);
	string best_tree_topo = best_tree_topo_ss.str();
	for (cur_iteration = 2; cur_iteration < params->min_iterations; cur_iteration++) {
		double min_elapsed = (getCPUTime() - params->startTime)/60;
		if (min_elapsed > params->maxtime) {
			maxTimeReached = true;
			break;
		}
		double parsimonyLH = computeParsimonyTreePhylolib();
		cout << "Parsimony score: " << parsimonyLH << endl;
		curScore = optimizeNNIRax();
		if (curScore > bestScore) {
			stringstream cur_tree_topo_ss;
			printTree(cur_tree_topo_ss, WT_TAXON_ID | WT_SORT_TAXA);
			if (cur_tree_topo_ss.str() != best_tree_topo) {
				cout << "BETTER TREE FOUND: " << curScore << endl;
				bestScore = curScore;
				best_tree_string.seekp(0, ios::beg);
				printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
				printResultTree();
			} else {
				// higher likelihood but the same tree topology
				bestScore = curScore;
				cout << "UPDATE BEST LOG-LIKELIHOOD: " << bestScore << endl;
			}
		}
		double cputime_secs = getCPUTime() - params->startTime;
		cout << "Iteration " << cur_iteration << " / LogL: " << curScore << " / CPU time elapsed: "
				<< (int) round(cputime_secs) << "s" << endl;;

	}

	// read in new tree
	int printBranchLengths = TRUE;
	Tree2String(raxmlTree->tree_string, raxmlTree,
			raxmlTree->start->back, printBranchLengths, TRUE, 0, 0,
			0, SUMMARIZE_LH, 0, 0);
	stringstream mytree;
	mytree << raxmlTree->tree_string;
	freeNode();
	readTree(mytree, rooted);
	setAlignment(aln);
}

double IQTree::doIQPNNI() {
	//double bestIQPScore = -DBL_MAX + 100;
	if (testNNI) {
		string str = params->out_prefix;
		str += ".nni";
		outNNI.open(str.c_str());
		outNNI
				<< "cur_lh\tzero_lh\tnni_lh1\tnni_lh2\topt_len\tnnilen1\tnnilen2\tnni_round"
				<< endl;
	}
	time_t begin_time, cur_time;
	time(&begin_time);
	string tree_file_name = params->out_prefix;
	tree_file_name += ".treefile";
	bestScore = curScore;
	//printResultTree(params);
	string treels_name = params->out_prefix;
	treels_name += ".treels";
	string out_lh_file = params->out_prefix;
	out_lh_file += ".treelh";
	string site_lh_file = params->out_prefix;
	site_lh_file += ".sitelh";

	if (params->print_tree_lh) {
		out_treelh.open(out_lh_file.c_str());
		out_sitelh.open(site_lh_file.c_str());
	}

	if (params->write_intermediate_trees)
		out_treels.open(treels_name.c_str());

	if (params->write_intermediate_trees && save_all_trees == 1) {
		printIntermediateTree(
				WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
	}

	//printTree(treels_name.c_str(), WT_NEWLINE | WT_BR_LEN);

	setRootNode(params->root);
	// keep the best tree into a string
	stringstream best_tree_string;
	printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);

	stringstream best_tree_topo_ss;
	printTree(best_tree_topo_ss, WT_TAXON_ID + WT_SORT_TAXA);
	string best_tree_topo = best_tree_topo_ss.str();

	// write tree's loglikelihood to a file (if nni_lh option is enabled)
	ofstream lh_file;
	if (params->nni_lh) {
		// Remove the .treefile ending and add iq-tree.lh ending to the file name
		string aln_file_name;
		aln_file_name.assign(tree_file_name).erase(tree_file_name.size() - 9);
		string lh_file_name = aln_file_name + ".iq-tree.lh";

		lh_file.open((lh_file_name).c_str());
		if (lh_file.is_open()) {
			lh_file.precision(15);
			lh_file << 1;
			lh_file << "\t";
			lh_file << bestScore;
			lh_file << endl;
		} else {
			cout << "Cannot open file " + lh_file_name;
		}
	}
	stop_rule.addImprovedIteration(1);
	int cur_iteration;
	bool speedupMsg = false;
	for (cur_iteration = 2; !stop_rule.meetStopCondition(cur_iteration);
			cur_iteration++) {
		double min_elapsed = (getCPUTime() - params->startTime)/60;
		if (min_elapsed > params->maxtime) {
		  cout << "Maximal running time of " << params->maxtime << " minutes reached" << endl;
			break;
		}
		// estimate logl_cutoff
		if (params->avoid_duplicated_trees && max_candidate_trees > 0
				&& treels_logl.size() > 1000) {
			int num_entries = floor(
					max_candidate_trees
							* ((double) cur_iteration
									/ stop_rule.getNumIterations()));
			if (num_entries < treels_logl.size() * 0.9) {
				DoubleVector logl = treels_logl;
				nth_element(logl.begin(),
						logl.begin() + (treels_logl.size() - num_entries),
						logl.end());
				logl_cutoff = logl[treels_logl.size() - num_entries] - 1.0;
			} else
				logl_cutoff = 0.0;

			if (cur_iteration % 10 == 0) {
				cout << treels.size() << " trees, " << treels_logl.size()
						<< " logls, logl_cutoff= " << logl_cutoff;
				if (params->store_candidate_trees)
					cout << " duplicates= " << duplication_counter << " ("
							<< (int) round(
									100
											* ((double) duplication_counter
													/ treels_logl.size()))
							<< "%)" << endl;
				else
					cout << endl;
			}
		}

		if (estimate_nni_cutoff && nni_info.size() >= 500) {
			estimate_nni_cutoff = false;
			estimateNNICutoff(params);
		}
		if (verbose_mode >= VB_DEBUG)
			cout << "Performing IQP in iteration " << cur_iteration << endl;
		double iqp_score;
		Alignment *saved_aln = aln;

		// randomize the neighbor orders for all nodes
		randomizeNeighbors();

		if (iqp_assess_quartet == IQP_BOOTSTRAP) {
			// create bootstrap sample
			enableHeuris = false;
			Alignment* bootstrap_alignment;
			if (aln->isSuperAlignment())
				bootstrap_alignment = new SuperAlignment;
			else
				bootstrap_alignment = new Alignment;
			bootstrap_alignment->createBootstrapAlignment(aln);
			setAlignment(bootstrap_alignment);
			initializeAllPartialLh();
			clearAllPartialLH();
			curScore = iqp_score = optimizeAllBranches();
		} else {
			if (params->reinsert_par) {
				curScore = doParsimonyReinsertion();
			} else {
				if (!params->raxmllib) {
					curScore = doIQP();
				} else {
					doIQP();
					//printTree( (string(params->out_prefix) + ".iqp_before." + convertIntToString(cur_iteration)).c_str() );
					stringstream iqp_tree_string;

					// TODO: only works for 1 partition
					transformBranchLenRAX(raxmlTree->fracchange);

					printTree(iqp_tree_string);
					treeReadLenString(iqp_tree_string.str().c_str(), raxmlTree, TRUE, FALSE, TRUE);
					smoothTree(raxmlTree, 1);
					evaluateGeneric(raxmlTree, raxmlTree->start, FALSE);
					cout << "IQP log-likelihood = " << raxmlTree->likelihood << endl;
					curScore = raxmlTree->likelihood;
				}
			}

		}

		setRootNode(params->root);

		 //cout.precision(15);
		if (verbose_mode >= VB_DEBUG) {
			string iqp_tree = tree_file_name + "IQP" + convertIntToString(
					cur_iteration);
			printTree(iqp_tree.c_str());
		}

		int skipped = 0;
		int nni_count = 0;
		if (enableHeuris) {
				if (cur_iteration > params->speedup_iter) {
					if (!speedupMsg) {
						speedupMsg = true;
						cout << "SPEED UP HEURISTIC ENABLED!" << endl;
					}
					if (!params->new_heuristic) {
						nni_count_est = estN95();
						nni_delta_est = estDelta95();
					} else {
						double nni_count_est95 = estN95();
						double nni_delta_est95 = estDelta95();
						double nni_count_estMedian = estNMedian();
						double nni_delta_estMedian = estDeltaMedian();
						double maxScore1 = nni_delta_est95 * nni_count_estMedian;
						double maxScore2 = nni_delta_estMedian * nni_count_est95;
						if (maxScore2 > maxScore1) {
							nni_count_est = nni_count_est95;
							nni_delta_est = nni_delta_estMedian;
						} else {
							nni_count_est = nni_count_estMedian;
							nni_delta_est = nni_delta_est95;
						}
						//nni_count_est = estN95();
						//nni_delta_est = estDelta95();
					}
					if (cur_iteration % 10 == 0)
						cout << "Estimated number of NNIs: " << nni_count_est
								<< ", delta-logl per NNI: " << nni_delta_est
								<< endl;
					if (params->raxmllib) {
						curScore = optimizeNNIRax(true, &skipped, &nni_count);
					} else {
						curScore = optimizeNNI(true, &skipped, &nni_count);
					}
				} else {
					if (params->raxmllib){
						curScore = optimizeNNIRax();
					} else {
						curScore = optimizeNNI();
					}
				}
		} else {
			if (params->raxmllib) {
				// Start NNI search using Phylolib kernel
				curScore = optimizeNNIRax();
			} else {
				curScore = optimizeNNI();
			}
		}

		if (iqp_assess_quartet == IQP_BOOTSTRAP) {
			// restore alignment
			delete aln;
			setAlignment(saved_aln);
			initializeAllPartialLh();
			clearAllPartialLH();
			curScore = optimizeAllBranches();
		}

		if (params->nni_lh && lh_file.is_open()) {
			lh_file << cur_iteration;
			lh_file << "\t";
			lh_file << iqp_score;
			lh_file << endl;

			lh_file << cur_iteration;
			lh_file << "\t";
			lh_file << curScore;
			lh_file << endl;
		}

		time(&cur_time);
		double cputime_secs = getCPUTime() - params->startTime;
		double cputime_remaining =
				(stop_rule.getNumIterations() - cur_iteration) * cputime_secs
						/ (cur_iteration - 1);
		/*double remaining_secs = (stop_rule.getNumIterations() - cur_iteration) *
		 elapsed_secs / (cur_iteration - 1);*/
		cout.setf(ios::fixed, ios::floatfield);

		// NNI search was skipped according to the speed up heuristics
		if (!skipped) {
			cout
					<< ((iqp_assess_quartet == IQP_BOOTSTRAP) ?
							"Bootstrap " : "Iteration ") << cur_iteration
					<< " / LogL: " << curScore << " / CPU time elapsed: "
					<< (int) round(cputime_secs) << "s";
			if (cur_iteration > 10 && cputime_secs > 10)
				cout << " (" << (int) round(cputime_remaining) << "s left)";
			cout << endl;
		} else {
			cout
					<< ((iqp_assess_quartet == IQP_BOOTSTRAP) ?
							"Bootstrap " : "Iteration ") << cur_iteration
					<< " skipped after " << nni_count << " NNI / LogL: "
					<< curScore << " / CPU time elapsed: "
					<< (int) round(cputime_secs) << "s";
			if (cur_iteration > 10 && cputime_secs > 10)
				cout << " (" << (int) round(cputime_remaining) << "s left)";
			cout << endl;
		}

		/*
		 if (verbose_mode >= VB_DEBUG) {
		 if (abs(curScore - bestScore) <= 0.0001) {
		 cout << "Found tree with the same score as best score" << endl;
		 if (!copyFile(
		 tree_file_name.c_str(),
		 (tree_file_name + ".bestTree" + convertIntToString(
		 cur_iteration)).c_str()))
		 cout << "Tree file could not be copied successfully";
		 printTree(
		 (tree_file_name + ".sameScoreBestTree"
		 + convertIntToString(cur_iteration)).c_str());
		 }
		 }
		 */

		if (params->write_intermediate_trees && save_all_trees == 1) {
			printIntermediateTree(
					WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
		}

		//printTree(treels_name.c_str(), WT_NEWLINE | WT_APPEND | WT_BR_LEN);

		if (curScore > bestScore) {
			if (params->raxmllib) {
				// read in new tree
				int printBranchLengths = TRUE;
				Tree2String(raxmlTree->tree_string, raxmlTree,
						raxmlTree->start->back, printBranchLengths, TRUE, 0, 0,
						0, SUMMARIZE_LH, 0, 0);
				stringstream mytree;
				mytree << raxmlTree->tree_string;
				//cout << mytree.str() << endl;
				freeNode();
				readTree(mytree, rooted);
				setAlignment(aln);
			}

			//curScore = optimizeAllBranches(100, 0.0001);
			stringstream cur_tree_topo_ss;
			printTree(cur_tree_topo_ss, WT_TAXON_ID | WT_SORT_TAXA);
			if (cur_tree_topo_ss.str() != best_tree_topo) {
				best_tree_topo = cur_tree_topo_ss.str();
				cout << "BETTER TREE FOUND: " << curScore << endl;
				//cout << "Saving new better tree ..." << endl;
				bestScore = curScore;
				best_tree_string.seekp(0, ios::beg);
				printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
				if (params->write_best_trees) {
					ostringstream iter_string;
					iter_string << cur_iteration;
					printResultTree(iter_string.str());
				}
				printResultTree();
				stop_rule.addImprovedIteration(cur_iteration);

				// Variable Neighborhood search idea, reset k_delete if tree is better
				//resetKDelete();
			} else {
				// higher likelihood but the same tree topology
				bestScore = curScore;
				cout << "UPDATE BEST LOG-LIKELIHOOD: " << bestScore << endl;
			}
		} else {
			/* take back the current best tree */
			best_tree_string.seekg(0, ios::beg);
			freeNode();
			readTree(best_tree_string, rooted);
			assignLeafNames();
			if (!params->raxmllib) {
				initializeAllPartialLh();
				clearAllPartialLH();
			}
			//cout << "Rollback tree topology ..." << endl;
			//cout << best_tree_string.str() << endl;
			//printTree(cout);
			//cout << endl;
			//cout << "Recompute best score: " << optimizeAllBranches() << endl;
			// Variable Neighborhood search idea, increase k_delete if tree is not better
			//increaseKDelete();
			//if (curScore > bestScore - 1e-4) // if goes back, increase k_delete once more
			//increaseKDelete();
		}
		if ((cur_iteration) % (params->step_iterations / 2) == 0
				&& params->gbo_replicates) {
			SplitGraph *sg = new SplitGraph;
			summarizeBootstrap(*sg);
			boot_splits.push_back(sg);
			max_candidate_trees = treels_logl.size()
					* (stop_rule.getNumIterations()) / cur_iteration;
			cout << "Setting tau = " << max_candidate_trees << endl;
		}
		if (cur_iteration == stop_rule.getNumIterations()
				&& params->gbo_replicates && !boot_splits.empty()
				&& stop_rule.getNumIterations() + params->step_iterations
						<= params->max_iterations) {
			//SplitGraph *sg = new SplitGraph;
			//summarizeBootstrap(*sg);
			if (!checkBootstrapStopping()) {
				max_candidate_trees = treels_logl.size()
						* (stop_rule.getNumIterations()
								+ params->step_iterations)
						/ stop_rule.getNumIterations();
				stop_rule.setIterationNum(
						stop_rule.getNumIterations() + params->step_iterations,
						params->max_iterations);
				cout << "INFO: Increase number of iterations to "
						<< stop_rule.getNumIterations() << " tau = "
						<< max_candidate_trees << endl;
				//delete boot_splits;
				//boot_splits = sg;
			} //else delete sg;
		}
	}

	int predicted_iteration = stop_rule.getPredictedIteration();
	//cout.unsetf(ios::fixed);

	if (predicted_iteration > cur_iteration) {
		cout << endl << "WARNING: " << predicted_iteration
				<< " iterations are needed to ensure that with a "
				<< floor(params->stop_confidence * 100) << "% confidence"
				<< endl
				<< "         the IQPNNI search will not find a better tree"
				<< endl;
	}

	if (testNNI)
		outNNI.close();
	if (params->write_intermediate_trees)
		out_treels.close();
	if (params->print_tree_lh) {
		out_treelh.close();
		out_sitelh.close();
	}

	return bestScore;
}

/****************************************************************************
 Fast Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/
double IQTree::optimizeNNI(bool beginHeu, int *skipped, int *nni_count_ret) {
	double oldCurScore = curScore;
	bool resetLamda = true;
	bool foundBetterTree = true;
	curLambda = startLambda;
	if (skipped)
		*skipped = 0;
	nni_round = 0; // number of NNI round
	int nni_count = 0;
	if (nni_count_ret)
		*nni_count_ret = nni_count;

	int nni2apply = 0; // number of nni to be applied in this round
	int nonconf_nni = 0; // number of non-conflicting NNIs found in this round
	do {
		if (resetLamda) { // tree get improved, lamda reset
			if (save_all_trees == 2) {
				saveCurrentTree(curScore); // BQM: for new bootstrap
			}
			++nni_round;
			if (verbose_mode >= VB_DEBUG) {
				cout << "Doing NNI round " << nni_round << endl;
				if (verbose_mode >= VB_DEBUG) {
					drawTree(cout,
							WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID
									| WT_NEWLINE);
				}
				if (isSuperTree()) {
					((PhyloSuperTree*) this)->printMapInfo();
				}
			}
			if (beginHeu) {
				double maxScore = curScore
						+ nni_delta_est * (nni_count_est - nni_count);
				if (maxScore < curScore)
					maxScore = curScore;
				if (maxScore <= bestScore) {
					if (skipped)
						*skipped = 1;
					return curScore;
				}
			}
			curLambda = startLambda;
			vec_nonconf_nni.clear();
			mapOptBranLens.clear();
			savedBranLens.clear();
			posNNIs.clear();
			saveBranLens(); // save all current branch lengths
			initPartitionInfo(); // for super tree
			if (!nni_sort) {
				genNNIMoves(params->approximate_nni); // generate all positive NNI moves
			} else {
				genNNIMovesSort(params->approximate_nni);
			}
			if (!containPosNNI(posNNIs)) { // no more positive NNI moves
				if (nni_round == 1)
					/* the local search fails to improve the tree */
					foundBetterTree = false;
				break;
			}
			/* sort all positive NNI moves (descending) */
			sort(posNNIs.begin(), posNNIs.end());
			/* remove conflicting NNIs */
			genNonconfNNIs();
			nonconf_nni = vec_nonconf_nni.size();
		}
		nni2apply = ceil(nonconf_nni * curLambda);
		if (nni2apply == 1)
			curLambda = 0.0;
		changeAllBranches();
		if (nni_round == 1)
			applyNNIs(nni2apply, true);
		else
			applyNNIs(nni2apply, false);
		double newScore = optimizeAllBranches(1);
		if (newScore > curScore) {
			if (enableHeuris) {
				if ( vecImpProNNI.size() < 1000 ) {
					vecImpProNNI.push_back((newScore - curScore) / nni2apply);
				} else {
					vecImpProNNI.erase( vecImpProNNI.begin() );
					vecImpProNNI.push_back((newScore - curScore) / nni2apply);
				}
			}
			nni_count += nni2apply;
			if (nni_count_ret)
				*nni_count_ret = nni_count;
			curScore = newScore; // Update the current score
			resetLamda = true;
		} else {
			curLambda = curLambda / 2;
			if (verbose_mode >= VB_MED) {
				cout << " Tree score gets worse at NNI round: " << nni_round
						<< endl;
				cout << " Number of NNI moves applied: " << nni2apply << endl;
				cout << " Rollback tree with new lambda = " << curLambda
						<< endl;
			}

			/* tree cannot be worse if only 1 NNI is applied */
			if (nni2apply == 1 && newScore < curScore) {
				cout << "THIS IS A BUG !!!" << endl;
				cout
						<< "The tree likelihood is supposed to be greater or equal than "
						<< vec_nonconf_nni.at(0).score << endl;
				cout << "Tree likelihood before the swap is " << curScore
						<< endl;
				cout << "Obtained tree likelihood is " << newScore << endl;
				if (curScore <= oldCurScore)
					foundBetterTree = false;
				// restore the tree by reverting all NNIs
				if (nni_round == 1)
					applyNNIs(nni2apply, true);
				else
					applyNNIs(nni2apply, false);

				// restore the branch lengths
				restoreAllBranLen();
				break;
			}

			// restore the tree by reverting all NNIs
			if (nni_round == 1)
				applyNNIs(nni2apply, true);
			else
				applyNNIs(nni2apply, false);

			// restore the branch lengths
			restoreAllBranLen();
			resetLamda = false;
		}
	} while (true);

	if (foundBetterTree) {
		curScore = optimizeAllBranches(1);
		if (enableHeuris) {
			if ( vecNumNNI.size() < 1000 ) {
				vecNumNNI.push_back(nni_count);
			} else {
				vecNumNNI.erase( vecNumNNI.begin() );
				vecNumNNI.push_back(nni_count);
			}
		}
	} else {
		if (verbose_mode >= VB_MED)
			cout << "Local search could not find any better tree !!!" << endl;
	}

	if (save_all_trees == 2 && params->nni_opt_5branches) {
		curScore = optimizeAllBranches();
		saveCurrentTree(curScore); // BQM: for new bootstrap
		saveNNITrees(); // optimize 5 branches around NNI, this makes program slower
	}
	return curScore;
}

double IQTree::optimizeNNIRax(bool beginHeu, int *skipped, int *nni_count_ret) {
	if (nnicut.num_delta == MAX_NUM_DELTA && nnicut.delta_min == DBL_MAX) {
		estDeltaMin();
		cout << "delta_min = " << nnicut.delta_min << endl;
	}
	int nniRound = 1;
	double curLH = raxmlTree->likelihood;
	//cout << "LH IQP Tree = " << curLH << endl;
	int nniApplied = 0;
	while (true) {
		if (beginHeu) {
			double maxScore = curLH
					+ nni_delta_est * (nni_count_est - nniApplied);
			if (maxScore < curLH) {
				beginHeu = false;
			} else {
				if (maxScore <= bestScore) {
					if (skipped)
						*skipped = 1;
					//cout << "nni_delta_est = " << nni_delta_est << endl;
					//cout << "nni_count_est = " << nni_count_est << endl;
					//cout << "nniApplied = " << nniApplied << endl;
					return curLH;
				}
			}
		}
		int nni_count = 0;
		double deltaNNI = 0.0;
		double newLH = doNNISearch(raxmlTree, &nni_count, &deltaNNI, &nnicut);
		if (newLH == 0.0) {
			break;
		} else {
			//cout << "NNI round " << nniRound << "  LH : " << newLH << endl;
			//cout << "NNI round " << nniRound << "  improvement : " << newLH -curLH << endl;
			curLH = newLH;
			nniRound++;
			//cout << "deltaNNI = " << deltaNNI << endl;
			//cout << "nni_count = " << nni_count << endl;
			if (enableHeuris) {
				if ( vecNumNNI.size() < 1000 ) {
					vecImpProNNI.push_back(deltaNNI);
				} else {
					vecImpProNNI.erase( vecImpProNNI.begin() );
					vecImpProNNI.push_back(deltaNNI);
				}
			}
			nniApplied += nni_count;
			if (nni_count_ret) {
				*nni_count_ret = nniApplied;
			}
		}
	}
	//cout << "nniApplied = " << nniApplied << endl;
	if (enableHeuris) {
		if ( vecNumNNI.size() < 1000 ) {
			vecNumNNI.push_back(nniApplied);
		} else {
			vecNumNNI.erase( vecNumNNI.begin() );
			vecNumNNI.push_back(nniApplied);
		}
	}

	// Re-optimize all branches
	smoothTree(raxmlTree, 1);
	evaluateGeneric(raxmlTree, raxmlTree->start, FALSE);
	return raxmlTree->likelihood;
}

void IQTree::applyNNIs(int nni2apply, bool neg_nni) {
	for (int i = 0; i < nni2apply; i++) {
		if (neg_nni) {
			doNNI(vec_nonconf_nni.at(i));
		} else {
			if (vec_nonconf_nni.at(i).delta <= 0)
				continue;
			doNNI(vec_nonconf_nni.at(i));
		}
		PhyloNode* node1 = vec_nonconf_nni.at(i).node1;
		PhyloNode* node2 = vec_nonconf_nni.at(i).node2;
		string key = bran2string(node1, node2);
		BranLenMap::iterator bran_it = mapOptBranLens.find(key);
		if (bran_it != mapOptBranLens.end()) {
			changeBranLen(node1, node2, bran_it->second);
		}
	}
}

void IQTree::genNonconfNNIs() {
	for (vector<NNIMove>::iterator iterMove = posNNIs.begin();
			iterMove != posNNIs.end(); iterMove++) {
		bool choosen = true;
		for (vector<NNIMove>::iterator iterNextMove = vec_nonconf_nni.begin();
				iterNextMove != vec_nonconf_nni.end(); iterNextMove++) {
			if ((*iterMove).node1 == (*(iterNextMove)).node1
					|| (*iterMove).node2 == (*(iterNextMove)).node1
					|| (*iterMove).node1 == (*(iterNextMove)).node2
					|| (*iterMove).node2 == (*(iterNextMove)).node2) {
				choosen = false;
				break;
			}
		}
		if (choosen) {
			vec_nonconf_nni.push_back(*iterMove);
		}
	}
}

double IQTree::estN95() {
	if (vecNumNNI.size() == 0) {
		return 0;
	} else {
		sort(vecNumNNI.begin(), vecNumNNI.end());
		int index = floor(vecNumNNI.size() * speed_conf);
		return vecNumNNI[index];
	}
}

double IQTree::estNMedian() {
	if (vecNumNNI.size() == 0) {
		return 0;
	} else {
		double median;
		size_t size = vecNumNNI.size();
		sort(vecNumNNI.begin(), vecNumNNI.end());
		if (size % 2 == 0) {
			median = (vecNumNNI[size / 2 + 1] + vecNumNNI[size / 2]) / 2;
		} else {
			median = vecNumNNI[size / 2];
		}
		return median;
	}
}

double IQTree::estDeltaMedian() {
	if (vecImpProNNI.size() == 0) {
		return 0;
	} else {
		double median;
		size_t size = vecImpProNNI.size();
		sort(vecImpProNNI.begin(), vecImpProNNI.end());
		if (size % 2 == 0) {
			median = (vecImpProNNI[size / 2 + 1] + vecImpProNNI[size / 2]) / 2;
		} else {
			median = vecImpProNNI[size / 2];
		}
		return median;
	}
}

inline double IQTree::estDelta95() {
	if (vecImpProNNI.size() == 0) {
		return 0;
	} else {
		sort(vecImpProNNI.begin(), vecImpProNNI.end());
		int index = floor(vecImpProNNI.size() * speed_conf);
		return vecImpProNNI[index];
	}
}

void IQTree::estDeltaMin() {
	sort(nnicut.delta, nnicut.delta + MAX_NUM_DELTA);
	int index = floor(MAX_NUM_DELTA * 0.95);
	nnicut.delta_min = nnicut.delta[index];
}

void IQTree::changeBranLen(PhyloNode *node1, PhyloNode *node2, double newlen) {
	current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
	current_it_back = (PhyloNeighbor*) node2->findNeighbor(node1);
	current_it->length = newlen;
	current_it_back->length = newlen;
	node1->clearReversePartialLh(node2);
	node2->clearReversePartialLh(node1);
}

inline string IQTree::bran2string(PhyloNode* node1, PhyloNode* node2) {
	string key("");
	if (node1->id < node2->id) {
		key += convertIntToString(node1->id) + "->"
				+ convertIntToString(node2->id);
	} else {
		key += convertIntToString(node2->id) + "->"
				+ convertIntToString(node1->id);
	}
	return key;
}

double IQTree::getBranLen(PhyloNode *node1, PhyloNode *node2) {
	current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
	assert(current_it);
	return current_it->length;
}

void IQTree::saveBranLens(PhyloNode *node, PhyloNode *dad) {
	if (!node) {
		node = (PhyloNode*) root;
	}
	if (dad) {
		double len = getBranLen(node, dad);
		string key = bran2string(node, dad);
		savedBranLens.insert(BranLenMap::value_type(key, len));
	}
	FOR_NEIGHBOR_IT(node, dad, it){
	saveBranLens((PhyloNode*) (*it)->node, node);
}
}

void IQTree::restoreAllBranLen(PhyloNode *node, PhyloNode *dad) {
	if (!node) {
		node = (PhyloNode*) root;
	}
	if (dad) {
		string key = bran2string(node, dad);
		current_it = (PhyloNeighbor*) node->findNeighbor(dad);
		assert(current_it);
		current_it_back = (PhyloNeighbor*) dad->findNeighbor(node);
		assert(current_it_back);
		current_it->length = savedBranLens[key];
		current_it_back->length = savedBranLens[key];
	}
	FOR_NEIGHBOR_IT(node, dad, it){
	restoreAllBranLen((PhyloNode*) (*it)->node, node);
}
}

inline double IQTree::getCurScore() {
	return curScore;
}

void IQTree::changeAllBranches(PhyloNode *node, PhyloNode *dad) {
	if (!node) {
		node = (PhyloNode*) root;
	}
	FOR_NEIGHBOR_IT(node, dad, it){
	string key = bran2string ((PhyloNode*) (*it)->node, (PhyloNode*) node);
	BranLenMap::iterator bran_it = mapOptBranLens.find(key);
	if ( bran_it != mapOptBranLens.end() ) {
		double curlen = (*it)->length;
		changeBranLen((PhyloNode*) (*it)->node, (PhyloNode*) node, curlen + curLambda*(bran_it->second-curlen));
	}
	changeAllBranches((PhyloNode*) (*it)->node, (PhyloNode*) node);
}

}

void IQTree::genNNIMoves(bool approx_nni, PhyloNode *node, PhyloNode *dad) {
	if (!node) {
		node = (PhyloNode*) root;
	}
	// internal Branch
	if (!node->isLeaf() && dad && !dad->isLeaf()) {
		NNIMove myMove = getBestNNIForBran(node, dad, approx_nni);
		if (myMove.score != 0) {
			addPositiveNNIMove(myMove);
		}
	}
	FOR_NEIGHBOR_IT(node, dad, it){
	genNNIMoves(approx_nni, (PhyloNode*) (*it)->node, node);
}

}

void IQTree::genNNIMovesSort(bool approx_nni) {
	NodeVector nodes1, nodes2;
	int i;
	double cur_lh = curScore;
	vector<IntBranchInfo> int_branches;

	getInternalBranches(nodes1, nodes2);
	assert(nodes1.size() == leafNum-3 && nodes2.size() == leafNum-3);

	for (i = 0; i < leafNum - 3; i++) {
		IntBranchInfo int_branch;
		PhyloNeighbor *node12_it = (PhyloNeighbor*) nodes1[i]->findNeighbor(
				nodes2[i]);
		//PhyloNeighbor *node21_it = (PhyloNeighbor*) nodes2[i]->findNeighbor(nodes1[i]);
		int_branch.lh_contribution = cur_lh
				- computeLikelihoodZeroBranch(node12_it,
						(PhyloNode*) nodes1[i]);
		if (int_branch.lh_contribution < 0.0)
			int_branch.lh_contribution = 0.0;
		if (int_branch.lh_contribution < fabs(nni_cutoff)) {
			int_branch.node1 = (PhyloNode*) nodes1[i];
			int_branch.node2 = (PhyloNode*) nodes2[i];
			int_branches.push_back(int_branch);
		}
	}
	std::sort(int_branches.begin(), int_branches.end(), int_branch_cmp);
	for (vector<IntBranchInfo>::iterator it = int_branches.begin();
			it != int_branches.end(); it++)
		if (it->lh_contribution >= 0.0) // evaluate NNI if branch contribution is big enough
				{
			NNIMove myMove = getBestNNIForBran(it->node1, it->node2, approx_nni,
					it->lh_contribution);
			if (myMove.score != 0) {
				addPositiveNNIMove(myMove);
				if (!estimate_nni_cutoff)
					for (vector<IntBranchInfo>::iterator it2 = it + 1;
							it2 != int_branches.end(); it2++) {
						if (it2->node1 == it->node1 || it2->node2 == it->node1
								|| it2->node1 == it->node2
								|| it2->node2 == it->node2)
							it2->lh_contribution = -1.0; // do not evaluate this branch later on
					}
			}
		} else { // otherwise, only optimize the branch length
			PhyloNode *node1 = it->node1;
			PhyloNode *node2 = it->node2;
			PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(
					node2);
			PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(
					node1);
			double stored_len = node12_it->length;
			optimizeOneBranch(node1, node2, false);
			string key("");
			if (node1->id < node2->id) {
				key += convertIntToString(node1->id) + "->"
						+ convertIntToString(node2->id);
			} else {
				key += convertIntToString(node2->id) + "->"
						+ convertIntToString(node1->id);
			}

			mapOptBranLens.insert(
					BranLenMap::value_type(key, node12_it->length));
			node12_it->length = stored_len;
			node21_it->length = stored_len;
		}
}

void IQTree::estimateNNICutoff(Params* params) {
	double *delta = new double[nni_info.size()];
	int i;
	for (i = 0; i < nni_info.size(); i++) {
		double lh_score[4];
		memmove(lh_score, nni_info[i].lh_score, 4 * sizeof(double));
		std::sort(lh_score + 1, lh_score + 4); // sort in ascending order
		delta[i] = lh_score[0] - lh_score[2];
		if (verbose_mode >= VB_MED)
		cout << i << ": " << lh_score[0] << " " << lh_score[1] << " " << lh_score[2] << " " << lh_score[3] << endl;
	}
	std::sort(delta, delta + nni_info.size());
	nni_cutoff = delta[nni_info.size() / 20];
	cout << endl << "Estimated NNI cutoff: " << nni_cutoff << endl;
	string file_name = params->out_prefix;
	file_name += ".nnidelta";
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name.c_str());
		for (i = 0; i < nni_info.size(); i++) {
			out << delta[i] << endl;
		}
		out.close();
		cout << "NNI delta printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}
	delete[] delta;
}

NNIMove IQTree::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2,
		bool approx_nni, double lh_contribution) {
	assert(node1->degree() == 3 && node2->degree() == 3);
	NNIMove noMove;
	NNIMove nnimoves[2];
	noMove.score = 0;

	/*
	 *  array contain log-likelihood of 3 trees
	 *  treelhs[0]: current tree
	 *  treelhs[1]: 1.NNI tree
	 *  treelhs[2]: 2.NNI tree
	 */
	double treelhs[3];
	//double approxTreeLH[3];
	PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
	PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);

	/*
	 *  Array to store the branch lengths
	 *  0.Element: the original length of the branch
	 *  1.Element: the optimized length
	 *  2. Element: the optimized length after the 1.NNI
	 *  3.Element: the optimized length after the 2.NNI
	 */

	double node12_len[4];
	node12_len[0] = node12_it->length;
	double bestScore;

	if (approx_nni) {
		double aprroxBran1 = estimateBranchLength(node21_it, node2);
		node12_it->length = aprroxBran1;
		node21_it->length = aprroxBran1;
		bestScore = computeLikelihoodBranch(node21_it, node2);
	} else {
		//bestScore = optimizeOneBranch(node1, node2, false);
		bestScore = curScore;
	}

	// This is done because it could happen that even after
	// optimizing branch length bestScore < curScore
	// Why is it, could it be a BUG?
	/*
	if (bestScore < curScore - 1.0E-6) {
		 cout.precision(15);
		 cout << "Likelihood of the tree reduced after optimizing branch length"
		 << endl;
		 cout << "curScore = " << curScore << endl;
		 cout << "LH after branOPT = " << bestScore << endl;
		bestScore = curScore;
		node12_it->length = node12_len[0];
		node21_it->length = node12_len[0];
	}
	*/
	treelhs[0] = bestScore;
	node12_len[1] = node12_it->length;

	// TEST BQM
	double lh_zero_branch;
	if (lh_contribution < 0.0) {
		// compute likelikelihood if branch length collapse to zero
		lh_zero_branch = computeLikelihoodZeroBranch(node12_it,
				(PhyloNode*) node1);
	} else {
		lh_zero_branch = bestScore - lh_contribution;
	}
	if (lh_zero_branch - bestScore < nni_cutoff) {
		string key = bran2string(node1, node2);
		mapOptBranLens.insert(BranLenMap::value_type(key, node12_len[1]));
		return noMove;
	}

	NNIInfo nni;
	nni.nni_round = nni_round;
	nni.lh_score[0] = lh_zero_branch;
	nni.lh_score[1] = bestScore;
	nni.br_len[1] = node12_it->length;

	if (testNNI) {
		outNNI.precision(8);
		outNNI << bestScore << "\t" << lh_zero_branch;
	}

	//restore the branch length before doing NNI
	node12_it->length = node12_len[0];
	node21_it->length = node12_len[0];

	// save the likelihood vector at the two ends of node1-node2
	double *node1_lh_save = node12_it->partial_lh;
	double *node2_lh_save = node21_it->partial_lh;
	//save scaling vector
	UBYTE *node1_scale_save = node12_it->scale_num;
	UBYTE *node2_scale_save = node21_it->scale_num;
	double node1_lh_scale = node12_it->lh_scale_factor;
	double node2_lh_scale = node21_it->lh_scale_factor;
	// save parsimony vector
	UINT *node1_pars_save = node12_it->partial_pars;
	UINT *node2_pars_save = node21_it->partial_pars;
	node12_it->partial_pars = newBitsBlock();
	node21_it->partial_pars = newBitsBlock();
	// save the first found neighbor of node 1 (excluding node2) in node1_it
	NeighborVec::iterator node1_it;
	for (node1_it = (node1)->neighbors.begin();
			node1_it != (node1)->neighbors.end(); node1_it++)
		if ((*node1_it)->node != (node2))
			break;
	Neighbor *node1_nei = *node1_it;
	double node1_len = node1_nei->length;
	int nniNr = 1;
	int chosenSwap = 1;

	// make alignment 16
	if (((intptr_t) tmp_partial_lh1) % 16 == 0)
		node12_it->partial_lh = tmp_partial_lh1;
	else
		node12_it->partial_lh = tmp_partial_lh1 + 1;

	if (((intptr_t) tmp_partial_lh2) % 16 == 0)
		node21_it->partial_lh = tmp_partial_lh2;
	else
		node21_it->partial_lh = tmp_partial_lh2 + 1;

	node12_it->scale_num = tmp_scale_num1;
	node21_it->scale_num = tmp_scale_num2;

	double bestDelta = 0;
	int lhCount = 0;

	/************************ Swap NNI ********************************/
	for (NeighborVec::iterator node2_it = (node2)->neighbors.begin();
			node2_it != (node2)->neighbors.end(); node2_it++)
		if ((*node2_it)->node != (node1)) {
			nniNr++;
			Neighbor *node2_nei = *node2_it;
			double node2_len = node2_nei->length;
			node1->updateNeighbor(node1_it, node2_nei);
			node2_nei->node->updateNeighbor(node2, node1);
			node2->updateNeighbor(node2_it, node1_nei);
			node1_nei->node->updateNeighbor(node1, node2);
			// clear partial likelihood vector
			node12_it->clearPartialLh();
			node21_it->clearPartialLh();
			lhCount++;
			double newScore;
			if (verbose_mode >= VB_MAX)
				cout << "current branch length: " << node12_len[0] << endl;

			/* compute the score of the swapped topology */
			if (approx_nni) {
				double estBran2 = estimateBranchLength(node21_it, node2);
				node12_it->length = estBran2;
				node21_it->length = estBran2;
				if (verbose_mode >= VB_MAX)
					cout << "approx branch length: " << estBran2 << endl;

				//approxTreeLH[lhCount] = computeLikelihoodBranch(node21_it, node2);
				newScore = computeLikelihoodBranch(node21_it, node2);
				//if (newScore <= treelhs[0])
				//    newScore = optimizeOneBranch(node1, node2, false);
			} else {

				node12_it->length = node12_len[0];
				node21_it->length = node12_len[0];

				newScore = optimizeOneBranch(node1, node2, false);
			}
			if (verbose_mode >= VB_MAX)
				cout << "optimal branch length: " << node12_it->length << endl;
			/*
			 if (approxTreeLH[lhCount] > approxTreeLH[0]) {
			 if (newScore < treelhs[0]) {
			 cout << "FALSE POSITIVE" << endl;
			 cout << "approx LH of current tree = " << approxTreeLH[0] << endl;
			 cout << "approx LH of NNI tree = " << approxTreeLH[lhCount] << endl;
			 cout << "optimal LH of current tree = " << treelhs[0] << endl;
			 cout << "optimal LH of NNI tree = " << newScore << endl;
			 }

			 }*/

			treelhs[nniNr - 1] = newScore;
			double delta = newScore - treelhs[0];

			if (save_all_trees == 2) {
				saveCurrentTree(newScore); // BQM: for new bootstrap
			}
			nni.lh_score[nniNr] = newScore;
			nni.br_len[nniNr] = node12_it->length;

			if (testNNI)
				outNNI << "\t" << newScore;

			// Save the branch length of the NNI nniNr
			node12_len[nniNr] = node12_it->length;

			nnimoves[nniNr - 2].node1Nei_it = node1_it;
			nnimoves[nniNr - 2].node2Nei_it = node2_it;
			nnimoves[nniNr - 2].score = newScore;
			nnimoves[nniNr - 2].node1 = node1;
			nnimoves[nniNr - 2].node2 = node2;
			nnimoves[nniNr - 2].delta = delta;

			if (nnimoves[nniNr - 2].delta > bestDelta + 1e-6) {
				bestDelta = nnimoves[nniNr - 2].delta;
				chosenSwap = nniNr;
			}

			// swap back and recover the branch lengths
			node1->updateNeighbor(node1_it, node1_nei, node1_len);
			node1_nei->node->updateNeighbor(node2, node1, node1_len);
			node2->updateNeighbor(node2_it, node2_nei, node2_len);
			node2_nei->node->updateNeighbor(node1, node2, node2_len);
			node12_it->length = node12_len[0];
			node21_it->length = node12_len[0];
		}

	if (testNNI) {
		outNNI.precision(3);
		outNNI << "\t" << node12_len[1] << "\t" << node12_len[2] << "\t"
				<< node12_len[3] << "\t" << nni_round << endl;
	}

	if (enable_parsimony) {
		delete[] node21_it->partial_pars;
		delete[] node12_it->partial_pars;
		node12_it->partial_pars = node1_pars_save;
		node21_it->partial_pars = node2_pars_save;
	}
	// restore the partial likelihood vector
	node12_it->partial_lh = node1_lh_save;
	node21_it->partial_lh = node2_lh_save;
	node12_it->unclearPartialLh();
	node21_it->unclearPartialLh();
	node12_it->partial_pars = node1_pars_save;
	node21_it->partial_pars = node2_pars_save;
	node12_it->scale_num = node1_scale_save;
	node21_it->scale_num = node2_scale_save;
	node12_it->lh_scale_factor = node1_lh_scale;
	node21_it->lh_scale_factor = node2_lh_scale;
	string key = bran2string(node1, node2);
	mapOptBranLens.insert(BranLenMap::value_type(key, node12_len[chosenSwap]));

	if (estimate_nni_cutoff)
		nni_info.push_back(nni);
	// if a positive NNI is found
	if (chosenSwap > 1) {
		return (nnimoves[chosenSwap - 2]);
	} else {
		return (noMove);
		// if no possitive NNI is found

		// accept an NNI with probability
		/*
		 double prob1 = 	1 / (1 + exp(treelhs[0] - treelhs[1]) + exp(treelhs[2] - treelhs[1]));
		 double prob2 = 	1 / (1 + exp(treelhs[0] - treelhs[2]) + exp(treelhs[1] - treelhs[2]));
		 double prob0 = 1 - prob1 - prob2;
		 double accept_prob = random_double();
		 if ( accept_prob <= prob0 ) {
		 return noMove;
		 } else if ( accept_prob <= prob0 + prob1 ) {
		 if (verbose_mode >= VB_MED)
		 cout << "Accepting a worse NNI with prob = " << prob1 << endl;
		 return nnimoves[0];
		 } else {
		 if (verbose_mode >= VB_MED)
		 cout << "Accepting a worse NNI with prob = " << prob2 << endl;
		 return nnimoves[1];
		 }
		 */
	}
	return noMove;
}

void IQTree::saveCurrentTree(double cur_logl) {
	ostringstream ostr;
	string tree_str;
	StringIntMap::iterator it = treels.end();
	if (params->store_candidate_trees) {
		printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
		tree_str = ostr.str();
		it = treels.find(tree_str);
	}
	int tree_index = -1;
	if (it != treels.end()) { // already in treels
		duplication_counter++;
		tree_index = it->second;
		if (cur_logl <= treels_logl[it->second] + 1e-4) {
			if (cur_logl < treels_logl[it->second] - 5.0)
				if (verbose_mode >= VB_MED)
					cout << "Current lh " << cur_logl
							<< " is much worse than expected "
							<< treels_logl[it->second] << endl;
			return;
		}
		if (verbose_mode >= VB_MAX)
			cout << "Updated logl " << treels_logl[it->second] << " to "
					<< cur_logl << endl;
		treels_logl[it->second] = cur_logl;
		if (save_all_br_lens) {
			ostr.seekp(ios::beg);
			printTree(ostr,
					WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE
							| WT_BR_LEN_ROUNDING);
			treels_newick[it->second] = ostr.str();
		}
		if (boot_samples.empty()) {
			computePatternLikelihood(treels_ptnlh[it->second], &cur_logl);
			return;
		}
		if (verbose_mode >= VB_MED)
			cout << "Update treels_logl[" << tree_index << "] := " << cur_logl
					<< endl;
	} else {
		if (logl_cutoff != 0.0 && cur_logl <= logl_cutoff + 1e-4)
			return;
		tree_index = treels_logl.size();
		if (params->store_candidate_trees)
			treels[tree_str] = tree_index;
		treels_logl.push_back(cur_logl);
		if (verbose_mode >= VB_MED)
			cout << "Add    treels_logl[" << tree_index << "] := " << cur_logl
					<< endl;
	}

	if (write_intermediate_trees)
		printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

	double *pattern_lh = new double[getAlnNPattern()];
	computePatternLikelihood(pattern_lh, &cur_logl);

	if (boot_samples.empty()) {
		// for runGuidedBootstrap
		treels_ptnlh.push_back(pattern_lh);
	} else {
		// online bootstrap
		int nptn = getAlnNPattern();
		int updated = 0;
		for (int sample = 0; sample < boot_samples.size(); sample++) {
			double rell = 0.0;
			for (int ptn = 0; ptn < nptn; ptn++)
				rell += pattern_lh[ptn] * boot_samples[sample][ptn];
			if (rell > boot_logl[sample]) {
				if (tree_str == "") {
					printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
					tree_str = ostr.str();
					it = treels.find(tree_str);
					if (it != treels.end()) {
						tree_index = it->second;
					} else {
						tree_index = treels.size();
						treels[tree_str] = tree_index;
					}
				}
				boot_logl[sample] = rell;
				boot_trees[sample] = tree_index;
				updated++;
			}
		}
		if (updated && verbose_mode >= VB_MED)
			cout << updated << " boot trees updated" << endl;
		/*
		 if (tree_index >= max_candidate_trees/2 && boot_splits->empty()) {
		 // summarize split support half way for stopping criterion
		 cout << "Summarizing current bootstrap supports..." << endl;
		 summarizeBootstrap(*boot_splits);
		 }*/
	}
	if (save_all_br_lens) {
		ostr.seekp(ios::beg);
		printTree(ostr,
				WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE
						| WT_BR_LEN_ROUNDING);
		treels_newick.push_back(ostr.str());
	}
	if (print_tree_lh) {
		out_treelh << cur_logl;
		double prob;
		aln->multinomialProb(pattern_lh, prob);
		out_treelh << "\t" << prob << endl;

		IntVector pattern_index;
		aln->getSitePatternIndex(pattern_index);
		out_sitelh << "Site_Lh   ";
		for (int i = 0; i < getAlnNSite(); i++)
			out_sitelh << " " << pattern_lh[pattern_index[i]];
		out_sitelh << endl;
	}
	if (!boot_samples.empty())
		delete[] pattern_lh;
}

void IQTree::saveNNITrees(PhyloNode *node, PhyloNode *dad) {
	if (!node) {
		node = (PhyloNode*) root;
	}
	if (dad && !node->isLeaf() && !dad->isLeaf()) {
		double *pat_lh1 = new double[aln->getNPattern()];
		double *pat_lh2 = new double[aln->getNPattern()];
		double lh1, lh2;
		computeNNIPatternLh(curScore, lh1, pat_lh1, lh2, pat_lh2, node, dad);
		delete[] pat_lh2;
		delete[] pat_lh1;
	}
	FOR_NEIGHBOR_IT(node, dad, it)saveNNITrees((PhyloNode*) (*it)->node, node);
}

void IQTree::summarizeBootstrap(Params &params, MTreeSet &trees) {
	int sum_weights = trees.sumTreeWeights();
	int i;
	if (verbose_mode >= VB_MED) {
		for (i = 0; i < trees.size(); i++)
			if (trees.tree_weights[i] > 0)
				cout << "Tree " << i + 1 << " weight= "
						<< (double) trees.tree_weights[i] * 100 / sum_weights
						<< endl;
	}
	int max_tree_id = max_element(trees.tree_weights.begin(),
			trees.tree_weights.end()) - trees.tree_weights.begin();
	cout << "max_tree_id = " << max_tree_id + 1 << "   max_weight = "
			<< trees.tree_weights[max_tree_id];
	cout << " (" << (double) trees.tree_weights[max_tree_id] * 100 / sum_weights
			<< "%)" << endl;
	// assign bootstrap support
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	vector<string> taxname;
	taxname.resize(leafNum);
	if (boot_splits.empty()) {
		getTaxaName(taxname);
	} else {
		boot_splits.back()->getTaxaName(taxname);
	}
	/*if (!tree.save_all_trees)
	 trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
	 else
	 trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false);
	 */
	trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false); // do not sort taxa

	cout << sg.size() << " splits found" << endl;

	if (!boot_splits.empty()) {
		// check the stopping criterion for ultra-fast bootstrap
		if (!checkBootstrapStopping())
			cout
					<< "WARNING: bootstrap analysis did not converge, rerun with higher number of iterations"
					<< endl;

	}

	sg.scaleWeight(1.0 / trees.sumTreeWeights(), false, 4);
	string out_file;
	out_file = params.out_prefix;
	out_file += ".splits";
	sg.saveFile(out_file.c_str(), IN_OTHER, true);
	cout << "Split supports printed to star-dot file " << out_file << endl;

	// compute the percentage of appearance
	sg.scaleWeight(100.0, true);
	//	printSplitSet(sg, hash_ss);
	//sg.report(cout);
	cout << "Creating bootstrap support values..." << endl;
	stringstream tree_stream;
	printTree(tree_stream, WT_TAXON_ID | WT_BR_LEN);
	MExtTree mytree;
	mytree.readTree(tree_stream, rooted);
	mytree.assignLeafID();
	mytree.createBootstrapSupport(taxname, trees, sg, hash_ss);

	// now write resulting tree with supports
	tree_stream.seekp(0, ios::beg);
	mytree.printTree(tree_stream);

	// now read resulting tree
	tree_stream.seekg(0, ios::beg);
	freeNode();
	readTree(tree_stream, rooted);
	assignLeafNames();
	initializeAllPartialLh();
	clearAllPartialLH();

	if (!save_all_trees) {
		out_file = params.out_prefix;
		out_file += ".suptree";

		printTree(out_file.c_str());
		cout << "Tree with assigned bootstrap support written to " << out_file
				<< endl;
	}

	out_file = params.out_prefix;
	out_file += ".splits.nex";
	sg.saveFile(out_file.c_str(), IN_NEXUS, true);
	cout << "Split supports printed to NEXUS file " << out_file << endl;

	out_file = params.out_prefix;
	out_file += ".supval";
	writeInternalNodeNames(out_file);

	cout << "Support values written to " << out_file << endl;
}

void IQTree::summarizeBootstrap(Params &params) {
	cout << "Summarizing from " << treels.size() << " candidate trees..."
			<< endl;
	MTreeSet trees;
	IntVector tree_weights;
	tree_weights.resize(treels_logl.size(), 0);
	for (int sample = 0; sample < boot_trees.size(); sample++)
		tree_weights[boot_trees[sample]]++;
	trees.init(treels, rooted, tree_weights);
	summarizeBootstrap(params, trees);
}

void IQTree::summarizeBootstrap(SplitGraph &sg) {
	MTreeSet trees;
	IntVector tree_weights;
	tree_weights.resize(treels_logl.size(), 0);
	for (int sample = 0; sample < boot_trees.size(); sample++)
		tree_weights[boot_trees[sample]]++;
	trees.init(treels, rooted, tree_weights);
	//SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	vector<string> taxname;
	taxname.resize(leafNum);
	getTaxaName(taxname);

	/*if (!tree.save_all_trees)
	 trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
	 else
	 trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false);
	 */
	trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false); // do not sort taxa
}

double computeCorrelation(IntVector &ix, IntVector &iy) {

	assert(ix.size() == iy.size());
	DoubleVector x;
	DoubleVector y;

	double mx = 0.0, my = 0.0; // mean value
	int i;
	x.resize(ix.size());
	y.resize(iy.size());
	for (i = 0; i < x.size(); i++) {
		x[i] = ix[i];
		y[i] = iy[i];
		mx += x[i];
		my += y[i];
	}
	mx /= x.size();
	my /= y.size();
	for (i = 0; i < x.size(); i++) {
		x[i] = x[i] / mx - 1.0;
		y[i] = y[i] / my - 1.0;
	}

	double f1 = 0.0, f2 = 0.0, f3 = 0.0;
	for (i = 0; i < x.size(); i++) {
		f1 += (x[i]) * (y[i]);
		f2 += (x[i]) * (x[i]);
		f3 += (y[i]) * (y[i]);
	}
	if (f2 == 0.0 || f3 == 0.0)
		return 1.0;
	return f1 / (sqrt(f2) * sqrt(f3));
}

bool IQTree::checkBootstrapStopping() {
	if (boot_splits.size() < 2)
		return false;
	IntVector split_supports;
	SplitIntMap split_map;
	int i;
	// collect split supports
	SplitGraph *sg = boot_splits.back();
	SplitGraph *half = boot_splits[(boot_splits.size() - 1) / 2];
	for (i = 0; i < half->size(); i++)
		if (half->at(i)->trivial() == -1) {
			split_map.insertSplit(half->at(i), split_supports.size());
			split_supports.push_back((int) (half->at(i)->getWeight()));
		}

	// collect split supports for new tree collection
	IntVector split_supports_new;
	split_supports_new.resize(split_supports.size(), 0);
	for (i = 0; i < sg->size(); i++)
		if ((*sg)[i]->trivial() == -1) {
			int index;
			Split *sp = split_map.findSplit((*sg)[i], index);
			if (sp) {
				// split found
				split_supports_new[index] = (int) ((*sg)[i]->getWeight());
			} else {
				// new split
				split_supports_new.push_back((int) ((*sg)[i]->getWeight()));
			}
		}
	cout << split_supports_new.size() - split_supports.size()
			<< " new splits compared to old boot_splits" << endl;
	if (split_supports_new.size() > split_supports.size())
		split_supports.resize(split_supports_new.size(), 0);

	// now compute correlation coefficient
	double corr = computeCorrelation(split_supports, split_supports_new);
	cout << "Correlation coefficient: " << corr << endl;
	// printing supports into file
	string outfile = params->out_prefix;
	outfile += ".splitsup";
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(outfile.c_str());
		out << "tau=" << max_candidate_trees / 2 << "\ttau="
				<< treels_logl.size() << endl;
		for (int i = 0; i < split_supports.size(); i++)
			out << split_supports[i] << "\t" << split_supports_new[i] << endl;
		out.close();
		cout << "Split support values printed to " << outfile << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, outfile);
	}
	return (corr >= params->min_correlation);
}

void IQTree::addPositiveNNIMove(NNIMove myMove) {
	posNNIs.push_back(myMove);
}

void IQTree::setRootNode(char *my_root) {
	string root_name;
	if (my_root)
		root_name = my_root;
	else
		root_name = aln->getSeqName(0);
	root = findNodeName(root_name);
	assert(root);
}

void IQTree::printResultTree(string suffix) {
	setRootNode(params->root);
	string tree_file_name = params->out_prefix;
	tree_file_name += ".treefile";
	if (suffix.compare("") != 0) {
		string iter_tree_name = tree_file_name + "." + suffix;
		printTree(iter_tree_name.c_str(),
				WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
	} else {
		printTree(tree_file_name.c_str(),
				WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
	}
	//printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH);
}

void IQTree::printResultTree(ostream &out) {
	setRootNode(params->root);
	printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
}

void IQTree::printIntermediateTree(int brtype) {
	setRootNode(params->root);
	bool duplicated_tree = false;
	double *pattern_lh = NULL;
	double logl = curScore;
	if (params->avoid_duplicated_trees) {
		// estimate logl_cutoff
		stringstream ostr;
		printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
		string tree_str = ostr.str();
		StringIntMap::iterator it = treels.find(tree_str);
		if (it != treels.end()) { // already in treels
			duplicated_tree = true;
			if (curScore > treels_logl[it->second] + 1e-4) {
				if (verbose_mode >= VB_MAX)
					cout << "Updated logl " << treels_logl[it->second] << " to "
							<< curScore << endl;
				treels_logl[it->second] = curScore;
				computeLikelihood(treels_ptnlh[it->second]);
				if (save_all_br_lens) {
					ostr.seekp(ios::beg);
					printTree(ostr,
							WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE
									| WT_BR_LEN_ROUNDING);
					treels_newick[it->second] = ostr.str();
				}
			}
			//pattern_lh = treels_ptnlh[treels[tree_str]];
		} else {
			//cout << __func__ << ": new tree" << endl;
			if (logl_cutoff != 0.0 && curScore <= logl_cutoff + 1e-4)
				duplicated_tree = true;
			else {
				treels[tree_str] = treels_ptnlh.size();
				pattern_lh = new double[aln->getNPattern()];
				computePatternLikelihood(pattern_lh, &logl);
				treels_ptnlh.push_back(pattern_lh);
				treels_logl.push_back(logl);
				if (save_all_br_lens) {
					ostr.seekp(ios::beg);
					printTree(ostr,
							WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE
									| WT_BR_LEN_ROUNDING);
					treels_newick.push_back(ostr.str());
				}
			}
		}
		//cout << tree_str << endl;
	} else {
		if (params->print_tree_lh) {
			pattern_lh = new double[aln->getNPattern()];
			computePatternLikelihood(pattern_lh, &logl);
		}
	}

	if (!duplicated_tree) {
		if (write_intermediate_trees)
			printTree(out_treels, brtype);
		if (params->print_tree_lh) {
			out_treelh.precision(10);
			out_treelh << logl;
			double prob;
			aln->multinomialProb(pattern_lh, prob);
			out_treelh << "\t" << prob << endl;
			if (!(brtype & WT_APPEND))
				out_sitelh << aln->getNSite() << endl;
			out_sitelh << "Site_Lh   ";
			for (int i = 0; i < aln->getNSite(); i++)
				out_sitelh << "\t" << pattern_lh[aln->getPatternID(i)];
			out_sitelh << endl;
			if (!params->avoid_duplicated_trees)
				delete[] pattern_lh;
		}
	}
	if (params->write_intermediate_trees == 1 && save_all_trees != 1) {
		return;
	}
	int x = save_all_trees;
	save_all_trees = 2;
	genNNIMoves(params->approximate_nni);
	save_all_trees = x;
}
