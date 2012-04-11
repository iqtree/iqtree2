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
#include "iqptree.h"
#include "phylosupertree.h"


//TODO Only to test
int cntBranches = 0;
extern clock_t t_begin;

IQPTree::IQPTree() :
PhyloTree() {
    k_represent = 0;
    //p_delete = 0.0;
    k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
    dist_matrix = NULL;
    //bonus_values = NULL;
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
}

IQPTree::IQPTree(Alignment *aln) : PhyloTree(aln) 
{
    k_represent = 0;
    //p_delete = 0.0;
    k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
    dist_matrix = NULL;
    //bonus_values = NULL;    
    nni_count_est = 0.0;
    nni_delta_est = 0;
    curScore = 0.0; // Current score of the tree
    bestScore = 0.0; // Best score found so far
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
}

void IQPTree::setParams(Params &params) {
    if ( params.min_iterations == -1 ) {
        if ( aln->getNSeq() < 100 )
            params.min_iterations = 200;
        else
            params.min_iterations = aln->getNSeq() * 2;
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
		k_delete = k_delete_min = k_delete_max = round(params.p_delete * leafNum);
    } else {
		k_delete = k_delete_min = 10;
		k_delete_max = leafNum / 2;
		if (k_delete_max > 100) k_delete_max = 100;
		if (k_delete_max < 20) k_delete_max = 20;
		k_delete_stay = ceil(leafNum/k_delete);
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
	testNNI	= params.testNNI;
	
	write_intermediate_trees = params.write_intermediate_trees;
	if (write_intermediate_trees > 2 || params.gbo_replicates > 0) 
		save_all_trees = 1;
	if (params.gbo_replicates > 0)
		if (params.iqp_assess_quartet != IQP_BOOTSTRAP) { save_all_trees = 2;  }
	if (params.gbo_replicates > 0 && params.do_compression) save_all_br_lens = true;
	print_tree_lh = params.print_tree_lh;
	max_candidate_trees = params.max_candidate_trees;
	setRootNode(params.root);
}

IQPTree::~IQPTree() {
    //if (bonus_values)
    //delete bonus_values;
    //bonus_values = NULL;
    if (dist_matrix)
        delete[] dist_matrix;
    dist_matrix = NULL;

	for (vector<double* >::reverse_iterator it = treels_ptnlh.rbegin(); it != treels_ptnlh.rend(); it++)
		delete [] (*it);
	treels_ptnlh.clear();
}

double IQPTree::getProbDelete() {
	return (double)k_delete / leafNum;
}

void IQPTree::resetKDelete() {
	k_delete = k_delete_min;
	k_delete_stay = ceil(leafNum / k_delete);
}

void IQPTree::increaseKDelete() {
	if (k_delete >= k_delete_max)
		return;
	k_delete_stay--;
	if (k_delete_stay > 0) return;
	k_delete++;
	k_delete_stay = ceil(leafNum / k_delete);
	if (verbose_mode >= VB_MED)
		cout << "Increase k_delete to " << k_delete << endl;
}

void IQPTree::setIQPIterations(STOP_CONDITION stop_condition,
        double stop_confidence, int min_iterations, int max_iterations) {
    stop_rule.setStopCondition(stop_condition);
    stop_rule.setConfidenceValue(stop_confidence);
    stop_rule.setIterationNum(min_iterations, max_iterations);
}

RepresentLeafSet* IQPTree::findRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, int nei_id,
        PhyloNode *dad) {
    PhyloNode *node = (PhyloNode*) (dad->neighbors[nei_id]->node);
    int set_id = dad->id * 3 + nei_id;
    if (leaves_vec[set_id]) return leaves_vec[set_id];
    RepresentLeafSet *leaves = new RepresentLeafSet;
    RepresentLeafSet * leaves_it[2] = {NULL, NULL};
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
        RepresentLeafSet::iterator lit[2] = {leaves_it[0]->begin(), leaves_it[1]->begin()};
        while (leaves->size() < k_represent) {
            int id = -1;
            if (lit[0] != leaves_it[0]->end() && lit[1] != leaves_it[1]->end()) {
                if ((*lit[0])->height < (*lit[1])->height)
                    id = 0;
                else if ((*lit[0])->height > (*lit[1])->height)
                    id = 1;
                else { // tie, choose at random
                    id = floor((((double)rand()) / RAND_MAX)*2);
					if (id > 1) id = 1;
                }
            } else
                if (lit[0] != leaves_it[0]->end())
                id = 0;
            else if (lit[1] != leaves_it[1]->end())
                id = 1;
            else break;
            assert(id < 2 && id >= 0);
            leaves->insert(new RepLeaf((*lit[id])->leaf, (*lit[id])->height + 1));
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

void IQPTree::deleteLeaves(PhyloNodeVector &del_leaves) {
    NodeVector taxa;
    // get the vector of taxa
    getTaxa(taxa);
    root = NULL;
    //int num_delete = floor(p_delete * taxa.size());
    int num_delete = k_delete;
    int i;
    if (num_delete > taxa.size() - 4) num_delete = taxa.size() - 4;
    // now try to randomly delete some taxa of the probability of p_delete
    for (i = 0; i < num_delete;) {
        int id = floor((((double) rand()) / RAND_MAX) * taxa.size());
        if (id >= taxa.size()) id = taxa.size()-1;
        if (!taxa[id]) continue;
        else i++;
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

int IQPTree::assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2,
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

int IQPTree::assessQuartetParsimony(Node *leaf0, Node *leaf1, Node *leaf2,
        Node *del_leaf) {
    int score[3] = {0, 0, 0};
    for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
        char ch0 = (*it)[leaf0->id];
        char ch1 = (*it)[leaf1->id];
        char ch2 = (*it)[leaf2->id];
        char chd = (*it)[del_leaf->id];
        if (ch0 >= aln->num_states || ch1 >= aln->num_states || ch2
                >= aln->num_states || chd >= aln->num_states)
            continue;
        if (chd == ch0 && ch1 == ch2)
            score[0] += (*it).frequency;
        if (chd == ch1 && ch0 == ch2)
            score[1] += (*it).frequency;
        if (chd == ch2 && ch0 == ch1)
            score[2] += (*it).frequency;
    }
    if (score[0] == score[1] && score[0] == score[2]) {
    	int id = floor(((double) (rand()) / RAND_MAX) * 3);
    	if (id > 2) id = 2;
        return id;
    }
    if (score[0] > score[1] && score[0] > score[2])
        return 0;
    if (score[1] < score[2])
        return 2;
    return 1;
}

void IQPTree::initializeBonus(PhyloNode *node, PhyloNode *dad) {
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

    FOR_NEIGHBOR_IT(node, dad, it) {
        initializeBonus((PhyloNode*) ((*it)->node), node);
    }
}

void IQPTree::raiseBonus(Neighbor *nei, Node *dad, double bonus) {
    ((PhyloNeighbor*) nei)->lh_scale_factor += bonus;
    if (verbose_mode >= VB_DEBUG)
        cout << dad->id << " - " << nei->node->id << " : " << bonus << endl;

    //  FOR_NEIGHBOR_IT(nei->node, dad, it)
    //	raiseBonus((*it), nei->node, bonus);
}

double IQPTree::computePartialBonus(Node *node, Node* dad) {
    PhyloNeighbor *node_nei = (PhyloNeighbor*) node->findNeighbor(dad);
    if (node_nei->partial_lh_computed) return node_nei->lh_scale_factor;

    FOR_NEIGHBOR_IT(node, dad, it) {
        node_nei->lh_scale_factor += computePartialBonus((*it)->node, node);
    }
    node_nei->partial_lh_computed = 1;
    return node_nei->lh_scale_factor;
}

void IQPTree::findBestBonus(double &best_score, NodeVector &best_nodes, NodeVector &best_dads, Node *node, Node *dad) {
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

    FOR_NEIGHBOR_IT(node, dad, it) {
        findBestBonus(best_score, best_nodes, best_dads, (*it)->node, node);
    }
}

void IQPTree::assessQuartets(vector<RepresentLeafSet*> &leaves_vec, PhyloNode *cur_root, PhyloNode *del_leaf) {
    const int MAX_DEGREE = 3;
    RepresentLeafSet * leaves[MAX_DEGREE];
    double bonus[MAX_DEGREE];
    memset(bonus, 0, MAX_DEGREE * sizeof (double));
    int cnt = 0;

    // only work for birfucating tree
    assert(cur_root->degree() == MAX_DEGREE);

    // find the representative leaf set for three subtrees

    FOR_NEIGHBOR_IT(cur_root, NULL, it) {
        leaves[cnt] = findRepresentLeaves(leaves_vec, cnt, cur_root);
        cnt++;
    }
    for (RepresentLeafSet::iterator i0 = leaves[0]->begin(); i0
            != leaves[0]->end(); i0++)
        for (RepresentLeafSet::iterator i1 = leaves[1]->begin(); i1
                != leaves[1]->end(); i1++)
            for (RepresentLeafSet::iterator i2 = leaves[2]->begin(); i2
                    != leaves[2]->end(); i2++) {
                int best_id;
                if (iqp_assess_quartet == IQP_DISTANCE)
                    best_id = assessQuartet((*i0)->leaf, (*i1)->leaf, (*i2)->leaf, del_leaf);
                else
                    best_id = assessQuartetParsimony((*i0)->leaf, (*i1)->leaf, (*i2)->leaf,
                        del_leaf);
                bonus[best_id] += 1.0;
            }
    for (cnt = 0; cnt < MAX_DEGREE; cnt++)
        if (bonus[cnt] > 0.0)
            raiseBonus(cur_root->neighbors[cnt], cur_root, bonus[cnt]);

}

void IQPTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
    PhyloNodeVector::iterator it_leaf;

    //int num_del_leaves = del_leaves.size();
    assert(root->isLeaf());

    for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
        if (verbose_mode >= VB_DEBUG)
            cout << "Reinserting " << (*it_leaf)->name << " (" << (*it_leaf)->id << ")" << endl;
        vector<RepresentLeafSet*> leaves_vec;
        leaves_vec.resize(nodeNum * 3, NULL);
        initializeBonus();
        NodeVector nodes;
        getInternalNodes(nodes);
        if (verbose_mode >= VB_DEBUG)
            drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
        //printTree(cout, WT_BR_LEN | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
        for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++) {
            assessQuartets(leaves_vec, (PhyloNode*) (*it), (*it_leaf));
        }
        NodeVector best_nodes, best_dads;
        double best_bonus;
        findBestBonus(best_bonus, best_nodes, best_dads);
        if (verbose_mode >= VB_DEBUG)
            cout << "Best bonus " << best_bonus << " " << best_nodes[0]->id << " " << best_dads[0]->id << endl;
        assert(best_nodes.size() == best_dads.size());
        int node_id = floor((((double) rand()) / RAND_MAX) * best_nodes.size());
        if (node_id >= best_nodes.size()) node_id = best_nodes.size()-1;
        if (best_nodes.size() > 1 && verbose_mode >= VB_DEBUG)
            cout << best_nodes.size()
            << " branches show the same best bonus, branch nr. "
            << node_id << " is chosen" << endl;

        reinsertLeaf(*it_leaf, best_nodes[node_id],
                best_dads[node_id]);
        //clearRepresentLeaves(leaves_vec, *it_node, *it_leaf);
        /*if (verbose_mode >= VB_DEBUG) {
         printTree(cout);
         cout << endl;
         }*/
        for (vector<RepresentLeafSet*>::iterator rit = leaves_vec.begin(); rit != leaves_vec.end(); rit++)
            if ((*rit)) {
                RepresentLeafSet *tit = (*rit);
                for (RepresentLeafSet::iterator rlit = tit->begin(); rlit != tit->end(); rlit++)
                    delete (*rlit);
                delete (*rit);
            }
    }
	initializeTree(); // BQM: re-index nodes and branches s.t. ping-pong neighbors have the same ID

	if (verbose_mode >= VB_DEBUG)
		drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
}

double IQPTree::doIQP() {
    if (verbose_mode >= VB_DEBUG)
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
    clock_t time_begin = clock();
    PhyloNodeVector del_leaves;
    deleteLeaves(del_leaves);
    reinsertLeaves(del_leaves);
    clock_t time_end = clock();
    if (verbose_mode >= VB_MAX) {
        cout << "IQP Time = " << (double) (time_end - time_begin) / CLOCKS_PER_SEC << endl;
    }
    // just to make sure IQP does it right
    setAlignment(aln);
    clearAllPartialLH();

	// optimize branches at the reinsertion point
/*
	for (PhyloNodeVector::iterator dit = del_leaves.begin(); dit != del_leaves.end(); dit++) {
		PhyloNode *adj_node = (PhyloNode*)(*dit)->neighbors[0]->node;
		FOR_NEIGHBOR_IT(adj_node, (*dit), it)
			curScore = optimizeOneBranch(adj_node, (PhyloNode*)(*it)->node);
		//curScore = optimizeOneBranch(adj_node, (PhyloNode*)(*dit));
	}
	double score = curScore;
*/
	/** BQM: please check careful! */
    //curScore = optimizeAllBranches(1);
    curScore = optimizeAllBranches(5, 1.0);
    //clearAllPartialLH();
    //double curScore2 = computeLikelihood();
    //cout << " diff = " << curScore - score << endl;
    if (enable_parsimony)
    	cur_pars_score = computeParsimony();

    if (verbose_mode >= VB_MAX) {
        cout << "IQP Likelihood = " << curScore << "  Parsimony = " << cur_pars_score << endl;
        //printTree(cout);
        //cout << endl;
    }
    return curScore;
}

/*void get2RandNumb(const int size, int &first, int &second) {
    // pick a random element
    first = floor((((double) rand()) / RAND_MAX) * size);
    // pick a random element from what's left (there is one fewer to choose from)...
    second = floor((((double) rand()) / RAND_MAX) * (size-1));
    // ...and adjust second choice to take into account the first choice
    if (second >= first) {
        ++second;
    }
}*/

double IQPTree::swapTaxa(PhyloNode *node1, PhyloNode *node2) {
    assert( node1->isLeaf() );
    assert( node2->isLeaf() );

    PhyloNeighbor *node1nei = (PhyloNeighbor*) *( node1->neighbors.begin() );
    PhyloNeighbor *node2nei = (PhyloNeighbor*) *( node2->neighbors.begin() );

    node2nei->node->updateNeighbor(node2, node1);
    node1nei->node->updateNeighbor(node1, node2);

    // Update the new neightbors of the 2 nodes
    node1->updateNeighbor(node1->neighbors.begin(), node2nei);
    node2->updateNeighbor(node2->neighbors.begin(), node1nei);

    PhyloNeighbor *node1NewNei = (PhyloNeighbor*) *( node1->neighbors.begin() );
    PhyloNeighbor *node2NewNei = (PhyloNeighbor*) *( node2->neighbors.begin() );

    // Reoptimize the branch lengths
    optimizeOneBranch(node1, (PhyloNode*) node1NewNei->node );
    this->curScore = optimizeOneBranch(node2, (PhyloNode*) node2NewNei->node );
    //drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
    return this->curScore;
}

double IQPTree::perturb(int times) {
    while (times > 0) {
        NodeVector taxa;
        // get the vector of taxa
        getTaxa(taxa);
        int taxonid1 = floor((((double) rand()) / RAND_MAX) * taxa.size());
        if (taxonid1 >= taxa.size()) taxonid1--;
        PhyloNode *taxon1 = (PhyloNode*) taxa[taxonid1];
        PhyloNode *taxon2;
        int dists[taxa.size()];
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

        cout << "Swapping node " << taxon1->id << " and node " << taxon2->id << endl;
        cout << "Distance " << minDist << endl;
        curScore = swapTaxa(taxon1, taxon2);
        //taxa.erase( taxa.begin() + taxaID1 );
        //taxa.erase( taxa.begin() + taxaID2 -1 );

        times--;
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
//    bestScore = curScore;
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
//        clock_t startClock = clock();
//        perturb(perturbLevel);
//        clock_t endClock = clock();
//        cout << "Perturbing Time = " << (double) (endClock - startClock) / CLOCKS_PER_SEC << endl;
//
//        startClock = clock();
//        optimizeNNI();
//        endClock = clock();
//        cout << "NNI Time = " << (double) (endClock - startClock) / CLOCKS_PER_SEC << endl;
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

double IQPTree::doIQPNNI(Params &params) {

	if (testNNI) { 
		string str = params.out_prefix;
		str += ".nni";
		outNNI.open(str.c_str());
		outNNI << "cur_lh\tzero_lh\tnni_lh1\tnni_lh2\topt_len\tnnilen1\tnnilen2\tnni_round" << endl;
	}

    time_t begin_time, cur_time;
    time(&begin_time);
    string tree_file_name = params.out_prefix;
    tree_file_name += ".treefile";
    bestScore = curScore;
    //printResultTree(params);
    string treels_name = params.out_prefix;
    treels_name += ".treels";
	string out_lh_file = params.out_prefix;
	out_lh_file += ".treelh";
	string site_lh_file = params.out_prefix;
	site_lh_file += ".sitelh";

	if (params.print_tree_lh) {
		out_treelh.open(out_lh_file.c_str());
		out_sitelh.open(site_lh_file.c_str());
	}

    if (params.write_intermediate_trees) 
    	out_treels.open(treels_name.c_str());

    printIntermediateTree(WT_SORT_TAXA | WT_NEWLINE | WT_BR_LEN, params);
    //printTree(treels_name.c_str(), WT_NEWLINE | WT_BR_LEN);

    // keep the best tree into a string
    stringstream best_tree_string;
    printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);

    // write tree's loglikelihood to a file (if nni_lh option is enabled)
    ofstream lh_file;
    if (params.nni_lh) {
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
    for (cur_iteration = 2; !stop_rule.meetStopCondition(cur_iteration); cur_iteration++) {

		// estimate logl_cutoff
		if (params.avoid_duplicated_trees && params.max_candidate_trees > 0 && treels_logl.size() > 1000) {
			int num_entries = params.max_candidate_trees * cur_iteration / stop_rule.getNumIterations();
			if (num_entries < treels_logl.size() * 0.9) {
				DoubleVector logl = treels_logl;
				nth_element(logl.begin(), logl.begin() + (treels_logl.size()-num_entries), logl.end());
				logl_cutoff = logl[treels_logl.size()-num_entries] - 5.0;
			} else logl_cutoff = 0.0;

			cout << treels_logl.size() << " entries and logl_cutoff = " << logl_cutoff << endl;
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
        	iqp_score = doIQP();
        }

		setRootNode(params.root);

/*        cout.precision(15);
		if (verbose_mode >= VB_MAX) {
			cout << "IQP score : " << iqp_score << endl;
			printf("Total time used for IQP : %8.6f seconds. \n",
					(double) (-startClock + endClock) / CLOCKS_PER_SEC);
		}

		if (verbose_mode >= VB_DEBUG) {
			string iqp_tree = tree_file_name + "IQP" + convertIntToString(
					cur_iteration);
			printTree(iqp_tree.c_str());
		}
*/

        int skipped = 0;
		int nni_count = 0;
        if (enableHeuris) {
            if (cur_iteration > SPEED_UP_AFTER_ITER_NUM) {
                nni_count_est = estimateNNICount();
                nni_delta_est = estimateNNIDelta();
				cout << "Estimated number of NNIs : "<< nni_count_est << endl;
				cout << "Estimated lh improvement per NNI : " << nni_delta_est << endl;
                optimizeNNI(true, &skipped, &nni_count);
            } else {
                optimizeNNI();
            }
        } else {
            optimizeNNI();
        }

		if (iqp_assess_quartet == IQP_BOOTSTRAP) {
			// restore alignment
			delete aln;
            setAlignment(saved_aln);
            initializeAllPartialLh();
            clearAllPartialLH();
            curScore = optimizeAllBranches();
		}

        if (params.nni_lh && lh_file.is_open()) {
            lh_file << cur_iteration;
            lh_file << "\t";
            lh_file << iqp_score;
            lh_file << endl;

            lh_file << cur_iteration;
            lh_file << "\t";
            lh_file << curScore;
            lh_file << endl;
        }

        cout.precision(6);
        time(&cur_time);
        double elapsed_secs = difftime(cur_time, begin_time);
        double remaining_secs = (stop_rule.getNumIterations() - cur_iteration) *
                elapsed_secs / (cur_iteration - 1);
        cout.setf(ios::fixed, ios::floatfield);
        
        if (!skipped) {
            cout << ((iqp_assess_quartet == IQP_BOOTSTRAP) ? "Bootstrap ":"Iteration ") << cur_iteration << " / LogL: " << curScore 
            	<< " / Time elapsed: " << convert_time(elapsed_secs) << "s";
           if (cur_iteration > 10 && elapsed_secs > 10)
           		cout <<	" (" << convert_time(remaining_secs) << "s left)";
           cout << endl;
        } else {
            cout << ((iqp_assess_quartet == IQP_BOOTSTRAP) ? "Bootstrap ":"Iteration ") << cur_iteration << " skipped after " << nni_count << " NNI / LogL: " << curScore 
            	<< " / Time elapsed: " << convert_time(elapsed_secs) << "s";
           if (cur_iteration > 10 && elapsed_secs > 10)
           		cout <<	" (" << convert_time(remaining_secs) << "s left)";
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

        printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN, params);
        //printTree(treels_name.c_str(), WT_NEWLINE | WT_APPEND | WT_BR_LEN);


        if (curScore > bestScore + TOL_LIKELIHOOD) {
            cout << "BETTER TREE FOUND: " << curScore << endl;
            bestScore = curScore;
            best_tree_string.seekp(0, ios::beg);
            printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);
            //cout << best_tree_string.str() << endl;
            printResultTree(params);
            stop_rule.addImprovedIteration(cur_iteration);

            // Variable Neighborhood search idea, reset k_delete if tree is better
            resetKDelete();
        } else {
            /* take back the current best tree */
            best_tree_string.seekg(0, ios::beg);
            freeNode();
            readTree(best_tree_string, rooted);
            assignLeafNames();
            initializeAllPartialLh();
            clearAllPartialLH();

            // Variable Neighborhood search idea, increase k_delete if tree is not better
            increaseKDelete();

            if (curScore > bestScore - 1e-4) // if goes back, increase k_delete once more
                increaseKDelete();
        }
    }

    int predicted_iteration = stop_rule.getPredictedIteration();
    cout.unsetf(ios::fixed);

    if (predicted_iteration > cur_iteration) {
        cout << endl << "WARNING: " << predicted_iteration
                << " iterations are needed to ensure that with a "
                << floor(params.stop_confidence * 100) << "% confidence" << endl
                << "         the IQPNNI search will not find a better tree" << endl;
    }

	if (testNNI) outNNI.close();
    if (params.write_intermediate_trees) out_treels.close();
	if (params.print_tree_lh) {
		out_treelh.close();
		out_sitelh.close();
	}

    return bestScore;
}

/****************************************************************************
 Fast Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/

double IQPTree::optimizeNNI(bool beginHeu, int *skipped, int *nni_count) {
    int nniIteration = 0; // Number of NNI steps
    bool resetLamda = true;
    int numbNNI = 0; // Total number of NNIs applied in this iteration
	if (nni_count) *nni_count = numbNNI;
    int nniTotal = 0; // Number of non-conflicting NNIs found in a NNI-step    
    bool foundBetterTree = true;
    double curLamb = startLambda;
	if (skipped) *skipped = 0;
    do {
		if (verbose_mode >= VB_MAX) {
        	cout << "Doing NNI round " << nniIteration << endl;
        	if (verbose_mode >= VB_DEBUG) 
        		drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
        	if (isSuperTree()) ((PhyloSuperTree*)this)->printMapInfo();
        }
        nbNNI = 0;
        // Tree get improved, lamda reset
        if (resetLamda) {
            if (save_all_trees == 2) saveCurrentTree(curScore); // BQM: for new bootstrap

            ++nniIteration;
            nni_round = nniIteration;
            //N IQPNNI iterations have been done
            if (beginHeu) {
                double maxScore = curScore + nni_delta_est * (nni_count_est - numbNNI);
                if (maxScore <= bestScore) {
                    if (skipped) *skipped = 1;
                    return curScore;
                }				
            }
            curLamb = startLambda;            
            indeNNIs.clear();
            mapOptBranLens.clear();
            savedBranLens.clear();
            posNNIs.clear();
            //Save all the current branch lengths
            saveBranLens();                
            //Generate all possible NNI moves
			if (!nni_sort)
				genNNIMoves();
			else 
				genNNIMovesSort();
            if (posNNIs.size() == 0) {
                if (nniIteration == 1)
                    foundBetterTree = false;
                break;
            }
            // Sort all the possible moves (this is descending)
            sort(posNNIs.begin(), posNNIs.end());
            // Remove conflicting NNIs
            for (vector<NNIMove>::iterator iterMove = posNNIs.begin(); iterMove
                    != posNNIs.end(); iterMove++) {
                bool choosen = true;
                for (vector<NNIMove>::iterator iterNextMove =
                        indeNNIs.begin(); iterNextMove
                        != indeNNIs.end(); iterNextMove++) {
                    if ((*iterMove).node1 == (*(iterNextMove)).node1
                            || (*iterMove).node2 == (*(iterNextMove)).node1
                            || (*iterMove).node1 == (*(iterNextMove)).node2
                            || (*iterMove).node2 == (*(iterNextMove)).node2) {
                        choosen = false;
                        break;
                    }
                }
                if (choosen) {
                    indeNNIs.push_back(*iterMove);
                }
            }
            nniTotal = indeNNIs.size();
        }
        nbNNI = (int) nniTotal * curLamb;
        if (nbNNI < 1)
            nbNNI = 1;
        //Applying all non-conflicting NNIs
        for (int i = 0; i < nbNNI; i++) {
            // Apply the calculated optimal branch length for the center branch
            applyBranchLengthChange(indeNNIs.at(i).node1, indeNNIs.at(i).node2, false);
            //double scoreBeforeNNI = computeLikelihood();
            doNNI(indeNNIs.at(i));
        }
        double newScore = optimizeAllBranches(1);
        if (newScore > curScore + TOL_LIKELIHOOD) {
            if (enableHeuris)
                vecImpProNNI.push_back((newScore - curScore) / nbNNI);
            numbNNI += nbNNI;                
			if (nni_count) *nni_count = numbNNI;
            curScore = newScore; // Update the current score
            resetLamda = true;            
        } else {
            curLamb = curLamb / 2;
            cout << "The tree didn't improve at NNI step " << nniIteration
                    << " (applied NNIs = " << nbNNI
                    << ") -> Trying new lamda = " << curLamb << endl;
            //Tree cannot be worse if only 1 NNI is applied
            if ( nbNNI == 1) {
                cout << "THIS IS A BUG !!!" << endl;
                cout << "The tree likelihood is supposed to be greater or equal than " << indeNNIs.at(0).score << endl;
                cout << "Tree likelihood before the swap is " << curScore << endl;
                cout << "Obtained tree likelihood is " << newScore << endl;
                exit(1);
            } 
            //Restore the tree by reverting all NNIs
            for (int i = (nbNNI - 1); i >= 0; i--) {
                doNNI(indeNNIs.at(i));
            }
            //Restore the branch lengths
            restoreAllBranLen();
            clearAllPartialLH();
            cout << "Likelihood of the tree is " << computeLikelihood() << endl;
            resetLamda = false;
        }
    } while (true);

    if (foundBetterTree) {
        curScore = optimizeAllBranches();
        if (save_all_trees == 2) saveCurrentTree(curScore); // BQM: for new bootstrap
        if (enableHeuris) {
            vecNumNNI.push_back(numbNNI);                                    
        }
    } else {
        //cout << "NNI search could not find any better tree !!!" << endl;
    }
    return curScore;
}

int IQPTree::estimateNNICount() {
    if ( vecNumNNI.size() == 0 ) {
        return 0;
    }
    else {
        sort(vecNumNNI.begin(), vecNumNNI.end());
        int index = floor ( vecNumNNI.size() * speed_conf );
        return vecNumNNI[index];
    }
}

double IQPTree::estimateNNIDelta() {
    if ( vecImpProNNI.size() == 0 )
        return 0;
    else {
        sort(vecImpProNNI.begin(), vecImpProNNI.end());
        int index = floor( vecImpProNNI.size() * speed_conf );
        return vecImpProNNI[index];
    }
}

void IQPTree::applyAllBranchLengthChanges(PhyloNode *node, PhyloNode *dad) {
    applyChildBranchChanges(node, dad);

    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->node->isLeaf()) {
            applyAllBranchLengthChanges((PhyloNode*) (*it)->node, node);
        }
    }
}

double IQPTree::applyBranchLengthChange(PhyloNode *node1, PhyloNode *node2, bool nonNNIBranch) {
    current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    assert(current_it);
    current_it_back = (PhyloNeighbor*) node2->findNeighbor(node1);
    assert(current_it_back);
    double current_len = current_it->length;
    string key("");
    if (node1->id < node2->id) {
        key += convertIntToString(node1->id) + "->" + convertIntToString(
                node2->id);
    } else {
        key += convertIntToString(node2->id) + "->" + convertIntToString(
                node1->id);
    }
    //double optLen = mapOptBranLens[key]; // 
	// BQM: what happens if key does not exist in the map?
    BranLenMap::iterator it = mapOptBranLens.find(key);
    double optLen = ((it != mapOptBranLens.end()) ? it->second : current_len);
    double newLen;
    if (nonNNIBranch) {
        newLen = current_len + startLambda * (optLen - current_len);
    } else {
        newLen = optLen;
    }
    current_it->length = newLen;
    current_it_back->length = newLen;
    node1->clearReversePartialLh(node2);
    node2->clearReversePartialLh(node1);
    return newLen;
}

double IQPTree::getBranLen(PhyloNode *node1, PhyloNode *node2) {
    current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    assert(current_it);
    return current_it->length;
}

void IQPTree::saveBranLens(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad) {
        double len = getBranLen(node, dad);
        string key("");
        if (node->id < dad->id) {
            key += convertIntToString(node->id) + "->" + convertIntToString(dad->id);
        } else {
            key += convertIntToString(dad->id) + "->" + convertIntToString(node->id);
        }
        savedBranLens.insert(BranLenMap::value_type(key, len));
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        saveBranLens((PhyloNode*) (*it)->node, node);
    }
}

void IQPTree::restoreAllBranLen(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad) {
        string key("");
        if (node->id < dad->id) {
            key += convertIntToString(node->id) + "->" + convertIntToString(dad->id);
        } else {
            key += convertIntToString(dad->id) + "->" + convertIntToString(node->id);
        }
        current_it = (PhyloNeighbor*) node->findNeighbor(dad);
        assert(current_it);
        current_it_back = (PhyloNeighbor*) dad->findNeighbor(node);
        assert(current_it_back);
        current_it->length = savedBranLens[key];
        current_it_back->length = savedBranLens[key];
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        restoreAllBranLen((PhyloNode*) (*it)->node, node);
    }
}

double IQPTree::getCurScore(void) {
    return curScore;
}

void IQPTree::applyChildBranchChanges(PhyloNode *node, PhyloNode *dad) {

    FOR_NEIGHBOR_IT(node, dad, it) {
        bool branchUsed = false;
        for (int i = 0; i < nbNNI; i++) {
            if ((node->id == indeNNIs.at(i).node1->id
                    && (*it)->node->id == indeNNIs.at(i).node2->id)
                    || (node->id == indeNNIs.at(i).node2->id
                    && (*it)->node->id
                    == indeNNIs.at(i).node1->id)) {
                branchUsed = true;
                break;
            }
        }
        if (branchUsed) {
            continue;
        }
        applyBranchLengthChange((PhyloNode*) node,
                (PhyloNode*) (*it)->node, true);
        cntBranches++;
    }

}


void IQPTree::genNNIMoves(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    //Internal Branch
    if (!node->isLeaf() && dad && !dad->isLeaf()) {
        NNIMove myMove = getBestNNIForBran(node, dad);
        if (myMove.score != 0) {
            addPossibleNNIMove(myMove);
        }
    }//External branch
    else if (dad) {
//        double optBranchLen = optimizeOneBranch(node, dad);
//        string key("");
//        if (node->id < dad->id) {
//            key += convertIntToString(node->id) + "->" + convertIntToString(
//                    dad->id);
//        } else {
//            key += convertIntToString(dad->id) + "->" + convertIntToString(
//                    node->id);
//        }
//        mapOptBranLens.insert(MapBranchLength::value_type(key, optBranchLen));
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        genNNIMoves((PhyloNode*) (*it)->node, node);
    }

}

struct IntBranchInfo {
	PhyloNode *node1;
	PhyloNode *node2;
	double lh_contribution; // log-likelihood contribution of this branch: L(T)-L(T|e=0)
};

inline int int_branch_cmp (const IntBranchInfo a, const IntBranchInfo b)
{
	return (a.lh_contribution < b.lh_contribution);
}

void IQPTree::genNNIMovesSort() {
	NodeVector nodes1, nodes2;
	int i;
	double cur_lh = curScore;
	vector<IntBranchInfo> int_branches;

	getInternalBranches(nodes1, nodes2);
	assert(nodes1.size() == leafNum-3 && nodes2.size() == leafNum-3);

	for (i = 0; i < leafNum-3; i++) {
		IntBranchInfo int_branch;
		PhyloNeighbor *node12_it = (PhyloNeighbor*) nodes1[i]->findNeighbor(nodes2[i]);
		PhyloNeighbor *node21_it = (PhyloNeighbor*) nodes2[i]->findNeighbor(nodes1[i]);
		double stored_len = node12_it->length;
		node12_it->length = 0.0;
		node21_it->length = 0.0;
		int_branch.lh_contribution = cur_lh - computeLikelihoodBranch(node12_it, (PhyloNode*) nodes1[i]);
		if (int_branch.lh_contribution < 0.0) int_branch.lh_contribution = 0.0;
		// restore branch length
		node12_it->length = stored_len;
		node21_it->length = stored_len;
		if (int_branch.lh_contribution < fabs(nni_cutoff)) {
			int_branch.node1 = (PhyloNode*) nodes1[i];
			int_branch.node2 = (PhyloNode*) nodes2[i];
			int_branches.push_back(int_branch);
		}
	}
	std::sort(int_branches.begin(), int_branches.end(), int_branch_cmp);
	for (vector<IntBranchInfo>::iterator it = int_branches.begin(); it != int_branches.end(); it++) 
	if (it->lh_contribution >= 0.0) // evaluate NNI if branch contribution is big enough
	{ 
		NNIMove myMove = getBestNNIForBran(it->node1, it->node2, it->lh_contribution);
		if (myMove.score != 0) {
            addPossibleNNIMove(myMove);
			if (!estimate_nni_cutoff)
			for (vector<IntBranchInfo>::iterator it2 = it+1; it2 != int_branches.end(); it2++) {
				if (it2->node1 == it->node1 || it2->node2 == it->node1 ||it2->node1 == it->node2 || it2->node2 == it->node2)
					it2->lh_contribution = -1.0; // do not evaluate this branch later on
			}
		}
	} else { // otherwise, only optimize the branch length
		PhyloNode *node1 = it->node1;
		PhyloNode *node2 = it->node2;
		PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
		PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
		double stored_len = node12_it->length;
		optimizeOneBranch(node1, node2, false);
		string key("");
		if (node1->id < node2->id) {
			key += convertIntToString(node1->id) + "->" + convertIntToString(node2->id);
		} else {
			key += convertIntToString(node2->id) + "->" + convertIntToString(node1->id);
		}    
	
		mapOptBranLens.insert(BranLenMap::value_type(key, node12_it->length));
		// restore branch length
		node12_it->length = stored_len;
		node21_it->length = stored_len;
	} 
}

void IQPTree::estimateNNICutoff(Params &params) {
	double delta[nni_info.size()];
	int i;
	for (i = 0; i < nni_info.size(); i++) {
		double lh_score[4];
		memmove(lh_score,nni_info[i].lh_score, 4*sizeof(double));
		std::sort(lh_score+1,lh_score+4); // sort in ascending order
		delta[i] = lh_score[0] - lh_score[2];
	}
	std::sort(delta,delta+nni_info.size());
	nni_cutoff = delta[nni_info.size()/10];
	cout << endl << "Estimated NNI cutoff: " << nni_cutoff << endl;
	string file_name = params.out_prefix;
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
}

NNIMove IQPTree::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, double lh_contribution) {
    assert(node1->degree() == 3 && node2->degree() == 3);
    NNIMove myMove;
    myMove.score = 0;
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

    //double myLH = computeLikelihood();
    //double lh_branch=computeLikelihoodBranch((PhyloNeighbor*) node1->findNeighbor(node2), (PhyloNode*) node1); 
    double bestScore = optimizeOneBranch(node1, node2, false);    	    
    node12_len[1] = node12_it->length;	

	// TEST BQM
    double lh_zero_branch;
    if (lh_contribution < 0.0) {
		// compute likelikelihood if branch length collapse to zero
		node12_it->length = 0.0;
		node21_it->length = 0.0;
		lh_zero_branch = computeLikelihoodBranch(node12_it, (PhyloNode*) node1); 
		// restore branch length
		node12_it->length = node12_len[1];
		node21_it->length = node12_len[1];
	} else {
		lh_zero_branch = bestScore - lh_contribution;
	} 
	if (lh_zero_branch - bestScore < nni_cutoff) {
		string key("");
		if (node1->id < node2->id) {
			key += convertIntToString(node1->id) + "->" + convertIntToString(node2->id);
		} else {
			key += convertIntToString(node2->id) + "->" + convertIntToString(node1->id);
		}    
	
		mapOptBranLens.insert(BranLenMap::value_type(key, node12_len[1]));

		return myMove;
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



//    if (bestScore < curScore - 0.001) {
//        cout << "This is a BUG. Please report it to the authors !" << endl;
//        cout << "gTree score is supposed to be improved after optimizing branch length !!!" << endl;
//        cout << "ComputeLikelihood = " << curScore << endl;
//        //cout << "myLH = " << myLH << endl;
//        cout << "optimizeOneBranch = " << bestScore << endl; 
//        //cout << "computeLikelihoodBranch = " << lh_branch << endl;
//        cout << "The current length is " << node12_len[0] << endl;
//        cout << "The new length is " << node12_it->length << endl;
//        //cout << "New likelihood is" << computeLikelihood() << endl;
//        cout << "Node 1: " << node1->id << endl;
//        cout << "Node 2: " << node2->id << endl;
//        exit(1);                
//    }
            
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

    // save the first found neighbor of node 1 (excluding node2) in node1_it
    FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)
        break;
    Neighbor *node1_nei = *node1_it;
    double node1_len = node1_nei->length;
    int nniNr = 1;
    int chosenSwap = 1;
    node12_it->partial_lh = tmp_partial_lh1;
    node21_it->partial_lh = tmp_partial_lh2;

    node12_it->scale_num = tmp_scale_num1;
    node21_it->scale_num = tmp_scale_num2;

    FOR_NEIGHBOR_IT(node2, node1, node2_it) {
        nniNr = nniNr + 1;      
        Neighbor *node2_nei = *node2_it;
        double node2_len = node2_nei->length;
        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);
        // clear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();        
        //        double lh_prediction = 100.0;
        // compute score with parsimony, accept topology if parsimony score is not so bad
        int pars_score = -10;
        if (enable_parsimony) {
            // replace partial_pars with a new vector
            node12_it->partial_pars = newBitsBlock();
            node21_it->partial_pars = newBitsBlock();
            pars_score = computeParsimonyBranch(node12_it, node1);
            //            if (linRegModel != NULL)
            //                lh_prediction = linRegModel->getValue(pars_score);
            //            else {
            //                for (int i = 0; i < 3000; i++) {
            //                    if (pars_scores[i] == 0) {
            //                        pars_scores[i] = pars_score;
            //                        newScore = optimizeOneBranch(node1, node2, false);
            //                        lh_scores[i] = newScore;
            //                        break;
            //                    }
            //                }
            //                if (pars_scores[2999] != 0) {
            //                    linRegModel = new Linear(3000, pars_scores, lh_scores);
            //                }
            //            }
            // If enough data points is collected, start linear regression
        }
        //if (lh_prediction > bestScore || pars_score < cur_pars_score)
        if (pars_score < cur_pars_score) {
            // compute the score of the swapped topology            
            double newScore = optimizeOneBranch(node1, node2, false);

            if (save_all_trees == 2) saveCurrentTree(newScore); // BQM: for new bootstrap

			nni.lh_score[nniNr] = newScore;
			nni.br_len[nniNr] = node12_it->length;

			if (testNNI) outNNI << "\t" << newScore;
                 // Save the branch length of the NNI nniNr
                node12_len[nniNr] = node12_it->length;                
            // If score is better, save the NNI move
            if (newScore > bestScore + TOL_LIKELIHOOD) {
//                cout << "Current score of the tree before this NNI round is :" << curScore << endl;
//                cout << "Likelihood of the tree :" << myLH << endl;
//                cout << "Best score found until now : " << bestScore << endl;
//                cout << "Found better score : " << newScore << endl;
                bestScore = newScore;
                chosenSwap = nniNr;
                myMove.node1Nei_it = node1_it;
                myMove.node2Nei_it = node2_it;
                myMove.score = bestScore;
                myMove.node1 = node1;
                myMove.node2 = node2;
            }
        } else {
            //cout << "pars filtered" << endl;
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
		outNNI << "\t" << node12_len[1] << "\t" << node12_len[2] << "\t" << node12_len[3] << "\t" << nni_round << endl;
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
    node12_it->scale_num = node1_scale_save;
    node21_it->scale_num = node2_scale_save;
    node12_it->lh_scale_factor = node1_lh_scale;
    node21_it->lh_scale_factor = node2_lh_scale;
    string key("");
    if (node1->id < node2->id) {
        key += convertIntToString(node1->id) + "->" + convertIntToString(node2->id);
    } else {
        key += convertIntToString(node2->id) + "->" + convertIntToString(node1->id);
    }

    mapOptBranLens.insert(BranLenMap::value_type(key, node12_len[chosenSwap]));

	if (estimate_nni_cutoff) nni_info.push_back(nni);

    //double newLH = computeLikelihood();
//    if (fabs(newLH - curScore) > TOL_LIKELIHOOD) {
//        cout << "Problem !! Tree was not restored correctly " << endl;
//        cout << "newLH = " << newLH << endl;
//        cout << "myLH = " << curScore << endl;
//        exit(0);
//   }
   //clearAllPartialLH();
    return myMove;
}

void IQPTree::saveCurrentTree(double cur_logl) {
	ostringstream ostr;
	printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	string tree_str = ostr.str();
	StringIntMap::iterator it = treels.find(tree_str);
	if (it != treels.end()) {// already in treels
		if (cur_logl <= treels_logl[it->second]+1e-4) {
			if (cur_logl < treels_logl[it->second]-5.0)
			cout << "Current lh " << cur_logl << " is much worse than expected " << treels_logl[it->second] << endl;
			return;
		}
		if (verbose_mode >= VB_MAX) 
			cout << "Updated logl " <<treels_logl[it->second] <<
				" to " << cur_logl << endl;
		treels_logl[it->second] = cur_logl;
		if (save_all_br_lens) {
			ostr.seekp(ios::beg);
			printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
			treels_newick[it->second] = ostr.str();
		}
		computePatternLikelihood(treels_ptnlh[it->second], &cur_logl);
		return;
	}
	else {
		if (logl_cutoff != 0.0 && cur_logl <= logl_cutoff + 1e-4) return;
		treels[tree_str] = treels_ptnlh.size();
	}

	if (write_intermediate_trees) printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

	double *pattern_lh = new double[aln->getNPattern()];
	computePatternLikelihood(pattern_lh, &cur_logl);

	treels_ptnlh.push_back(pattern_lh);
	treels_logl.push_back(cur_logl);
	if (save_all_br_lens) {
		ostr.seekp(ios::beg);
		printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
		treels_newick.push_back(ostr.str());
	}
	if (print_tree_lh) {
		out_treelh << cur_logl;
		double prob;
		aln->multinomialProb(pattern_lh, prob);
		out_treelh << "\t" << prob << endl;

		out_sitelh << "Site_Lh   ";
		for (int i = 0; i < aln->getNSite(); i++)
			out_sitelh << "\t" << pattern_lh[aln->getPatternID(i)];
		out_sitelh << endl;
	}
}

void IQPTree::addPossibleNNIMove(NNIMove myMove) {
    posNNIs.push_back(myMove);
}

void IQPTree::setRootNode(char *my_root) {
    string root_name;
    if (my_root) root_name = my_root;
    else root_name = aln->getSeqName(0);
    root = findNodeName(root_name);
    assert(root);
}

void IQPTree::printResultTree(Params &params) {
    setRootNode(params.root);
    string tree_file_name = params.out_prefix;
    tree_file_name += ".treefile";
    printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
    //printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH);
}

void IQPTree::printResultTree(Params &params, ostream &out) {
    setRootNode(params.root);
    printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
}

void IQPTree::printIntermediateTree(int brtype, Params &params) {
	if (!params.write_intermediate_trees && save_all_trees != 1) return;
	setRootNode(params.root);
	bool duplicated_tree = false;
	double *pattern_lh = NULL;
	double logl = curScore;
	if (params.avoid_duplicated_trees) {
		// estimate logl_cutoff
		stringstream ostr;
		printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
		string tree_str = ostr.str();
		StringIntMap::iterator it = treels.find(tree_str);
		if (it != treels.end()) { // already in treels
			duplicated_tree = true;
			if (curScore > treels_logl[it->second] + 1e-4) {
				if (verbose_mode >= VB_MAX) 
					cout << "Updated logl " <<treels_logl[it->second] <<
						" to " << curScore << endl;
				treels_logl[it->second] = curScore;
				computeLikelihood(treels_ptnlh[it->second]);
				if (save_all_br_lens) {
					ostr.seekp(ios::beg);
					printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
					treels_newick[it->second] = ostr.str();
				}
			}
			//pattern_lh = treels_ptnlh[treels[tree_str]];
		}
		else {
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
					printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
					treels_newick.push_back(ostr.str());
				}
			}
		}
		//cout << tree_str << endl;
	} else {
		if (params.print_tree_lh) {
			pattern_lh = new double[aln->getNPattern()];
			computePatternLikelihood(pattern_lh, &logl);
		}
	}

	if (!duplicated_tree) {
		if (write_intermediate_trees) printTree(out_treels, brtype);
		if (params.print_tree_lh) {
			out_treelh.precision(10);
			out_treelh << logl;
			double prob;
			aln->multinomialProb(pattern_lh, prob);
			out_treelh << "\t" << prob << endl;
			if (!(brtype & WT_APPEND)) out_sitelh << aln->getNSite() << endl;
			out_sitelh << "Site_Lh   ";
			for (int i = 0; i < aln->getNSite(); i++)
				out_sitelh << "\t" << pattern_lh[aln->getPatternID(i)];
			out_sitelh << endl;
			if (!params.avoid_duplicated_trees) delete [] pattern_lh;
		}
	}
	if (params.write_intermediate_trees == 1 && save_all_trees != 1) {
		return;
	}
	int x = save_all_trees;
	save_all_trees = 2;
	genNNIMoves();
	save_all_trees = x;
/*
	ofstream *out = NULL;
	if (params.write_intermediate_trees) out = &out_treels;

	if (params.print_tree_lh)
		if (params.avoid_duplicated_trees) 
			PhyloTree::optimizeNNI(0.0, NULL, NULL, out, WT_NEWLINE | WT_BR_LEN, &out_treelh, &out_sitelh, 
			&treels, &treels_ptnlh, &treels_logl, &max_candidate_trees, &logl_cutoff);
		else
			PhyloTree::optimizeNNI(0.0, NULL, NULL, out, WT_NEWLINE | WT_BR_LEN, &out_treelh, &out_sitelh);
	else
		if (params.avoid_duplicated_trees) 
			PhyloTree::optimizeNNI(0.0, NULL, NULL, out, WT_NEWLINE | WT_BR_LEN, NULL, NULL, 
			&treels, &treels_ptnlh, &treels_logl, &max_candidate_trees, &logl_cutoff);
		else
			PhyloTree::optimizeNNI(0.0, NULL, NULL, out, WT_NEWLINE | WT_BR_LEN);*/
	// now print 1-NNI-away trees
}
