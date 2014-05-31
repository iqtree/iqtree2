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
#include "phylosupertreeplen.h"
#include "mexttree.h"
#include "timeutil.h"
#include "gtrmodel.h"
#include "rategamma.h"
#include <numeric>
#include "pll/pllInternal.h"
#include "nnisearch.h"

Params *globalParam;

IQTree::IQTree() :
        PhyloTree() {
    init();
}

void IQTree::init() {
    k_represent = 0;
    k_delete = k_delete_min = k_delete_max = k_delete_stay = 0;
    dist_matrix = NULL;
    var_matrix = NULL;
    nni_count_est = 0.0;
    nni_delta_est = 0;
    curScore = 0.0; // Current score of the tree
    bestScore = -DBL_MAX; // Best score found so far
    curIteration = 1;
    cur_pars_score = -1;
    enable_parsimony = false;
    estimate_nni_cutoff = false;
    nni_cutoff = -1e6;
    nni_sort = false;
    testNNI = false;
    print_tree_lh = false;
    write_intermediate_trees = 0;
    max_candidate_trees = 0;
    logl_cutoff = 0.0;
    len_scale = 10000;
    save_all_br_lens = false;
    duplication_counter = 0;
    pllInst = NULL;
    pllAlignment = NULL;
    pllPartitions = NULL;
    //boot_splits = new SplitGraph;
    diversification = false;
}

IQTree::IQTree(Alignment *aln) :
        PhyloTree(aln) {
    init();
}
void IQTree::setParams(Params &params) {
    searchinfo.speednni = params.speednni;
    searchinfo.nni_type = params.nni_type;
    optimize_by_newton = params.optimize_by_newton;
    sse = params.SSE;
    setStartLambda(params.lambda);
    if (params.maxtime != 1000000) {
        params.autostop = false;
    }
    if (params.min_iterations == -1) {
        if (!params.gbo_replicates) {
            if (params.autostop) {
                params.min_iterations = aln->getNSeq() * 100;
            } else if (aln->getNSeq() < 100) {
                params.min_iterations = 200;
            } else {
                params.min_iterations = aln->getNSeq() * 2;
            }
            if (params.iteration_multiple > 1)
                params.min_iterations = aln->getNSeq() * params.iteration_multiple;
        } else {
            params.min_iterations = 100;
        }
    }
    if (params.gbo_replicates)
        params.max_iterations = max(params.max_iterations, max(params.min_iterations, 1000));

    k_represent = params.k_representative;

    if (params.p_delete == -1.0) {
        if (aln->getNSeq() < 4)
            params.p_delete = 0.0; // delete nothing
        else if (aln->getNSeq() == 4)
            params.p_delete = 0.25; // just delete 1 leaf
        else if (aln->getNSeq() == 5)
            params.p_delete = 0.4; // just delete 2 leaves
        else if (aln->getNSeq() < 51)
            params.p_delete = 0.5;
        else if (aln->getNSeq() < 100)
            params.p_delete = 0.3;
        else if (aln->getNSeq() < 200)
            params.p_delete = 0.2;
        else if (aln->getNSeq() < 400)
            params.p_delete = 0.1;
        else
            params.p_delete = 0.05;
    }
    //tree.setProbDelete(params.p_delete);
    if (params.p_delete != -1.0) {
        k_delete = k_delete_min = k_delete_max = ceil(params.p_delete * leafNum);
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
    globalParam = &params;

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

    string bootaln_name = params.out_prefix;
    bootaln_name += ".bootaln";
    if (params.print_bootaln) {
        ofstream bootalnout;
    	bootalnout.open(bootaln_name.c_str());
    	bootalnout.close();
    }

    if (params.online_bootstrap && params.gbo_replicates > 0) {
        cout << "Generating " << params.gbo_replicates << " samples for ultrafast bootstrap..." << endl;
        boot_samples.resize(params.gbo_replicates);
        boot_logl.resize(params.gbo_replicates, -DBL_MAX);
        boot_trees.resize(params.gbo_replicates, -1);
        boot_counts.resize(params.gbo_replicates, 0);
        VerboseMode saved_mode = verbose_mode;
        verbose_mode = VB_QUIET;
        for (int i = 0; i < params.gbo_replicates; i++) {
        	if (params.print_bootaln) {
    			Alignment* bootstrap_alignment;
    			if (aln->isSuperAlignment())
    				bootstrap_alignment = new SuperAlignment;
    			else
    				bootstrap_alignment = new Alignment;
    			bootstrap_alignment->createBootstrapAlignment(aln, &(boot_samples[i]), params.bootstrap_spec);
				bootstrap_alignment->printPhylip(bootaln_name.c_str(), true);
				delete bootstrap_alignment;
        	} else
        		aln->createBootstrapAlignment(boot_samples[i], params.bootstrap_spec);
        }
        verbose_mode = saved_mode;
        if (params.print_bootaln) {
        	cout << "Bootstrap alignments printed to " << bootaln_name << endl;
        }

        cout << "Max candidate trees (tau): " << max_candidate_trees << endl;
    }

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

    if (var_matrix)
        delete[] var_matrix;
    var_matrix = NULL;

    for (vector<double*>::reverse_iterator it = treels_ptnlh.rbegin(); it != treels_ptnlh.rend(); it++)
        delete[] (*it);
    treels_ptnlh.clear();
    for (vector<SplitGraph*>::reverse_iterator it2 = boot_splits.rbegin(); it2 != boot_splits.rend(); it2++)
        delete (*it2);
    //if (boot_splits) delete boot_splits;
    if (pllInst)
        pllDestroyInstance(pllInst);
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

void IQTree::setIQPIterations(STOP_CONDITION stop_condition, double stop_confidence, int min_iterations,
        int max_iterations) {
    stop_rule.setStopCondition(stop_condition);
    stop_rule.setConfidenceValue(stop_confidence);
    stop_rule.setIterationNum(min_iterations, max_iterations);
}

RepresentLeafSet* IQTree::findRepresentLeaves(vector<RepresentLeafSet*> &leaves_vec, int nei_id, PhyloNode *dad) {
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
        RepresentLeafSet::iterator lit[2] = { leaves_it[0]->begin(), leaves_it[1]->begin() };
        while (leaves->size() < k_represent) {
            int id = -1;
            if (lit[0] != leaves_it[0]->end() && lit[1] != leaves_it[1]->end()) {
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

void IQTree::initLeafFrequency(PhyloNode *node, PhyloNode *dad) {
    if (!node)
        node = (PhyloNode*) root;
    if (node->isLeaf()) {
        LeafFreq leaf_freq;
        leaf_freq.leaf_id = node->id;
        leaf_freq.freq = 0;
        leaf_freqs.push_back(leaf_freq);
    }
    for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
        if ((*it)->node != dad) {
            initLeafFrequency((PhyloNode*) (*it)->node, node);
        }
}

void IQTree::clearLeafFrequency() {
	for (vector<LeafFreq>::iterator it = leaf_freqs.begin(); it != leaf_freqs.end(); it++) {
		(*it).freq = 0;
	}
}

void IQTree::deleteNonCherryLeaves(PhyloNodeVector &del_leaves) {
    NodeVector cherry_taxa;
    NodeVector noncherry_taxa;
    // get the vector of non cherry taxa
    getNonCherryLeaves(noncherry_taxa, cherry_taxa);
    root = NULL;
    int num_taxa = aln->getNSeq();
    int num_delete = k_delete;
    if (num_delete > num_taxa - 4)
        num_delete = num_taxa - 4;
    if (verbose_mode >= VB_DEBUG) {
        cout << "Deleting " << num_delete << " leaves" << endl;
    }
    vector<unsigned int> indices_noncherry(noncherry_taxa.size());
    //iota(indices_noncherry.begin(), indices_noncherry.end(), 0);
    unsigned int startValue = 0;
    for (vector<unsigned int>::iterator it = indices_noncherry.begin(); it != indices_noncherry.end(); ++it) {
        (*it) = startValue;
        ++startValue;
    }
    random_shuffle(indices_noncherry.begin(), indices_noncherry.end());
    int i;
    for (i = 0; i < num_delete && i < noncherry_taxa.size(); i++) {
        PhyloNode *taxon = (PhyloNode*) noncherry_taxa[indices_noncherry[i]];
        del_leaves.push_back(taxon);
        deleteLeaf(taxon);
        //cout << taxon->id << ", ";
    }
    int j = 0;
    if (i < num_delete) {
        vector<unsigned int> indices_cherry(cherry_taxa.size());
        //iota(indices_cherry.begin(), indices_cherry.end(), 0);
        startValue = 0;
        for (vector<unsigned int>::iterator it = indices_cherry.begin(); it != indices_cherry.end(); ++it) {
            (*it) = startValue;
            ++startValue;
        }
        random_shuffle(indices_cherry.begin(), indices_cherry.end());
        while (i < num_delete) {
            PhyloNode *taxon = (PhyloNode*) cherry_taxa[indices_cherry[j]];
            del_leaves.push_back(taxon);
            deleteLeaf(taxon);
            i++;
            j++;
        }
    }
    root = cherry_taxa[j];
}

void IQTree::deleteNonTabuLeaves(PhyloNodeVector &del_leaves) {
    // sort node frequency
    sort(leaf_freqs.begin(), leaf_freqs.end());
    NodeVector taxa;
    // get the vector of taxa
    getTaxa(taxa);

    // shuffle the sorted list first so that leaves with the same frequency have the same
    // chance of being removed
    int cur_freq = leaf_freqs[0].freq;

    //cout << "size: " << leaf_freqs.size() << endl;
    //cout << "Before randomize" << endl;
    /*
     for (size_t i = 0; i < leaf_freqs.size(); i++) {
     cout << leaf_freqs[i].leaf_id << ":" << leaf_freqs[i].freq << " ";
     }
     cout << endl;
     */
    int startIndex = 0;
    int endIndex = 0;
    while (true) {
        while (leaf_freqs[endIndex].freq == cur_freq && endIndex < leaf_freqs.size()) {
            endIndex++;
        }
        random_shuffle(leaf_freqs.begin() + startIndex, leaf_freqs.begin() + endIndex);
        if (endIndex == leaf_freqs.size())
            break;
        startIndex = endIndex;
        cur_freq = leaf_freqs[endIndex].freq;
    }
    /*
     cout << "After randomize" << endl;
     for (size_t i = 0; i < leaf_freqs.size(); i++) {
     cout << leaf_freqs[i].leaf_id << ":" << leaf_freqs[i].freq << " ";
     }
     cout << endl;
     */
    int leafCnt = 0;
    while (leafCnt < k_delete) {
        PhyloNode *taxon = (PhyloNode*) taxa[leaf_freqs[leafCnt].leaf_id];
        del_leaves.push_back(taxon);
        leaf_freqs[leafCnt].freq++;
        leafCnt++;
        deleteLeaf(taxon);
    }
    root = taxa[leaf_freqs[leafCnt].leaf_id];
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
    if (verbose_mode >= VB_DEBUG) {
        cout << "Deleting " << num_delete << " leaves" << endl;
    }
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

int IQTree::assessQuartet(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf) {
    assert(dist_matrix);
    int nseq = aln->getNSeq();
    //int id0 = leaf0->id, id1 = leaf1->id, id2 = leaf2->id;
    double dist0 = dist_matrix[leaf0->id * nseq + del_leaf->id] + dist_matrix[leaf1->id * nseq + leaf2->id];
    double dist1 = dist_matrix[leaf1->id * nseq + del_leaf->id] + dist_matrix[leaf0->id * nseq + leaf2->id];
    double dist2 = dist_matrix[leaf2->id * nseq + del_leaf->id] + dist_matrix[leaf0->id * nseq + leaf1->id];
    if (dist0 < dist1 && dist0 < dist2)
        return 0;
    if (dist1 < dist2)
        return 1;
    return 2;
}

int IQTree::assessQuartetParsimony(Node *leaf0, Node *leaf1, Node *leaf2, Node *del_leaf) {
    int score[3] = { 0, 0, 0 };
    for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
        char ch0 = (*it)[leaf0->id];
        char ch1 = (*it)[leaf1->id];
        char ch2 = (*it)[leaf2->id];
        char chd = (*it)[del_leaf->id];
        if (ch0 >= aln->num_states || ch1 >= aln->num_states || ch2 >= aln->num_states || chd >= aln->num_states)
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


/*
 bool IQTree::containPosNNI(vector<NNIMove> posNNIs) {
 for (vector<NNIMove>::iterator iter = posNNIs.begin(); iter != posNNIs.end(); iter++) {
 if (iter->newloglh > iter->oldloglh)
 return true;
 }
 return false;
 }
 */

void IQTree::findBestBonus(double &best_score, NodeVector &best_nodes, NodeVector &best_dads, Node *node, Node *dad) {
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

void IQTree::assessQuartets(vector<RepresentLeafSet*> &leaves_vec, PhyloNode *cur_root, PhyloNode *del_leaf) {
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
    for (RepresentLeafSet::iterator i0 = leaves[0]->begin(); i0 != leaves[0]->end(); i0++)
        for (RepresentLeafSet::iterator i1 = leaves[1]->begin(); i1 != leaves[1]->end(); i1++)
            for (RepresentLeafSet::iterator i2 = leaves[2]->begin(); i2 != leaves[2]->end(); i2++) {
                int best_id;
                if (iqp_assess_quartet == IQP_DISTANCE)
                    best_id = assessQuartet((*i0)->leaf, (*i1)->leaf, (*i2)->leaf, del_leaf);
                else
                    best_id = assessQuartetParsimony((*i0)->leaf, (*i1)->leaf, (*i2)->leaf, del_leaf);
                bonus[best_id] += 1.0;
            }
    for (cnt = 0; cnt < MAX_DEGREE; cnt++)
        if (bonus[cnt] > 0.0)
            raiseBonus(cur_root->neighbors[cnt], cur_root, bonus[cnt]);

}

void IQTree::reinsertLeavesByParsimony(PhyloNodeVector &del_leaves) {
    PhyloNodeVector::iterator it_leaf;
    assert(root->isLeaf());
    for (it_leaf = del_leaves.begin(); it_leaf != del_leaves.end(); it_leaf++) {
        //cout << "Add leaf " << (*it_leaf)->id << " to the tree" << endl;
        initializeAllPartialPars();
        clearAllPartialLH();
        Node *target_node = NULL;
        Node *target_dad = NULL;
        Node *added_node = (*it_leaf)->neighbors[0]->node;
        Node *node1 = NULL;
        Node *node2 = NULL;
        //Node *leaf;
        for (int i = 0; i < 3; i++) {
            if (added_node->neighbors[i]->node->id == (*it_leaf)->id) {
                //leaf = added_node->neighbors[i]->node;
            } else if (!node1) {
                node1 = added_node->neighbors[i]->node;
            } else {
                node2 = added_node->neighbors[i]->node;
            }
        }

        //cout << "(" << node1->id << ", " << node2->id << ")" << "----" << "(" << added_node->id << "," << leaf->id << ")" << endl;
        added_node->updateNeighbor(node1, (Node*) 1);
        added_node->updateNeighbor(node2, (Node*) 2);

        addTaxonMPFast(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        target_node->updateNeighbor(target_dad, added_node, -1.0);
        target_dad->updateNeighbor(target_node, added_node, -1.0);
        added_node->updateNeighbor((Node*) 1, target_node, -1.0);
        added_node->updateNeighbor((Node*) 2, target_dad, -1.0);

    }

}

void IQTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
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
        int node_id = random_int(best_nodes.size());
        if (best_nodes.size() > 1 && verbose_mode >= VB_DEBUG)
            cout << best_nodes.size() << " branches show the same best bonus, branch nr. " << node_id << " is chosen"
                    << endl;

        reinsertLeaf(*it_leaf, best_nodes[node_id], best_dads[node_id]);
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

void IQTree::doParsimonyReinsertion() {
    PhyloNodeVector del_leaves;

    deleteLeaves(del_leaves);

    reinsertLeavesByParsimony(del_leaves);
    fixNegativeBranch(false);
}

void IQTree::setBestTree(string treeString, double treeLogl) {
    bestTreeString = treeString;
    bestScore = treeLogl;
}

bool IQTree::updateRefTreeSet(string treeString, double treeLogl) {
    stringstream backupTree;
    printTree(backupTree);
    bool updated = false;
    readTreeString(treeString);
    setRootNode(params->root);
    stringstream treeTopoSS;
    printTree(treeTopoSS, WT_TAXON_ID + WT_SORT_TAXA);
    string treeTopo = treeTopoSS.str();

    double worstLogl = -DBL_MAX;
    if (!refTreeSetSorted.empty()) {
        worstLogl = refTreeSetSorted.begin()->first;
    }
    if (refTreeSet.size() < params->popSize) {
        if (refTreeSet.find(treeTopo) != refTreeSet.end()) {
			cout << "Tree topology already exists in the reference set" << endl;
            updated = false;
        } else {
            refTreeSet.insert(make_pair(treeTopo, treeLogl));
            refTreeSetSorted.insert(make_pair(treeLogl, treeString));
            updated = true;
        }
    } else if (treeLogl > worstLogl && refTreeSet.find(treeTopo) == refTreeSet.end()) {
        if (refTreeSet.size() == params->popSize) {
            // remove the worst tree
            refTreeSetSorted.erase(refTreeSetSorted.begin());
            for (unordered_map<string, double>::iterator it = refTreeSet.begin(); it != refTreeSet.end(); ++it) {
                if (it->second == worstLogl) {
                    refTreeSet.erase(it);
                    break;
                }
            }
        }
        refTreeSet.insert(make_pair(treeTopo, treeLogl));
        refTreeSetSorted.insert(make_pair(treeLogl, treeString));
        assert(refTreeSet.size() == refTreeSetSorted.size() &&
                refTreeSetSorted.size() == params->popSize);
        updated = true;
    } else if (refTreeSet.find(treeTopo) != refTreeSet.end()) {
        cout << "Tree topology is identical with one of the tree in the candidate set" << endl;
    }
    if (updated) {
        printLoglInTreePop();
    }

    backupTree.seekg(0, ios::beg);
    freeNode();
    readTree(backupTree, rooted);
    setAlignment(aln);
    return updated;
}

void IQTree::doRandomNNIs(int numNNI) {
    map<int, Node*> usedNodes;
    NodeVector nodeList1, nodeList2;
    getInternalBranches(nodeList1, nodeList2);
    int numInBran = nodeList1.size();
    assert(numInBran == aln->getNSeq() - 3);
    for (int i = 0; i < numNNI; i++) {
        int index = random_int(numInBran);
        if (usedNodes.find(nodeList1[index]->id) == usedNodes.end()
                && usedNodes.find(nodeList2[index]->id) == usedNodes.end()) {
            doOneRandomNNI(nodeList1[index], nodeList2[index]);
            usedNodes.insert(map<int, Node*>::value_type(nodeList1[index]->id, nodeList1[index]));
            usedNodes.insert(map<int, Node*>::value_type(nodeList2[index]->id, nodeList2[index]));
        } else {
            usedNodes.clear();
            nodeList1.clear();
            nodeList2.clear();
            getInternalBranches(nodeList1, nodeList2);
            doOneRandomNNI(nodeList1[index], nodeList2[index]);
            usedNodes.insert(map<int, Node*>::value_type(nodeList1[index]->id, nodeList1[index]));
            usedNodes.insert(map<int, Node*>::value_type(nodeList2[index]->id, nodeList2[index]));
        }
    }
}

void IQTree::doIQP() {
    if (verbose_mode >= VB_DEBUG)
        drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
    //double time_begin = getCPUTime();
    PhyloNodeVector del_leaves;
    deleteLeaves(del_leaves);
    reinsertLeaves(del_leaves);

    // just to make sure IQP does it right
    setAlignment(aln);

    if (enable_parsimony) {
        cur_pars_score = computeParsimony();
        if (verbose_mode >= VB_MAX) {
            cout << "IQP Likelihood = " << curScore << "  Parsimony = " << cur_pars_score << endl;
        }
    }
}

double IQTree::inputTree2PLL(string treestring, bool computeLH) {
    double res = 0.0;
    // read in the tree string from IQTree kernel
    pllNewickTree *newick = pllNewickParseString(treestring.c_str());
    pllTreeInitTopologyNewick(pllInst, newick, PLL_FALSE);
    pllNewickParseDestroy(&newick);
    if (computeLH) {
        pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
        res = pllInst->likelihood;
    }
    return res;
}

double* IQTree::getModelRatesFromPLL() {
    assert(aln->num_states == 4);
    int numberOfRates = (pllPartitions->partitionData[0]->states * pllPartitions->partitionData[0]->states
            - pllPartitions->partitionData[0]->states) / 2;
    double* rate_params = new double[numberOfRates];
    for (int i = 0; i < numberOfRates; i++) {
        rate_params[i] = pllPartitions->partitionData[0]->substRates[i];
    }
    return rate_params;
}

void IQTree::printPLLModParams() {
    for (int part = 0; part < pllPartitions->numberOfPartitions; part++) {
        //printf("alpha[%d]: %f \n", part, pllPartitions->partitionData[part]->alpha);
        cout << "Alpha[" << part << "]" << ": " << pllPartitions->partitionData[part]->alpha << endl;
        if (aln->num_states == 4) {
            int states, rates;
            states = pllPartitions->partitionData[part]->states;
            rates = ((states * states - states) / 2);
            //printf(rates, "rates[%d] ac ag at cg ct gt: ", part);
            cout << "Rates[" << part << "]: " << " ac ag at cg ct gt: ";
            for (int i = 0; i < rates; i++) {
                //printf(rates,"%f ", pllPartitions->partitionData[part]->substRates[i]);
                cout << pllPartitions->partitionData[part]->substRates[i] << " ";
            }
            cout << endl;
            cout << "Frequencies: ";
            for (int i = 0; i < 4; i++) {
                //printf("%f ", pllPartitions->partitionData[part]->empiricalFrequencies[i]);
                cout << pllPartitions->partitionData[part]->empiricalFrequencies[i] << " ";
            }
            cout << endl;
        }
    }
}

void IQTree::printLoglInTreePop() {
    cout << "Logl of trees in population" << endl;
    for (map<double, string>::iterator it = refTreeSetSorted.begin(); it != refTreeSetSorted.end(); ++it) {
        cout << it->first << " / ";
    }
    cout << endl;
}

void IQTree::printRefTrees() {
    string filename = string(params->out_prefix) + ".reftrees";
    ofstream file;
    file.open(filename.c_str());
    for (map<double, string>::iterator it = refTreeSetSorted.begin(); it != refTreeSetSorted.end(); ++it) {
        file << it->second << endl;
    }
}

double IQTree::getAlphaFromPLL() {
    return pllPartitions->partitionData[0]->alpha;
}

void IQTree::inputModelParam2PLL() {
    // get the alpha parameter
    double alpha = getRate()->getGammaShape();
    if (alpha == 0.0)
        alpha = PLL_ALPHA_MAX;
    if (aln->num_states == 4) {
        // get the rate parameters
        double *rate_param = new double[6];
        getModel()->getRateMatrix(rate_param);
        // get the state frequencies
        double *state_freqs = new double[aln->num_states];
        getModel()->getStateFrequency(state_freqs);

        /* put them into PLL */
        stringstream linkagePattern;
        int partNr;
        for (partNr = 0; partNr < pllPartitions->numberOfPartitions - 1; partNr++) {
            linkagePattern << partNr << ",";
        }
        linkagePattern << partNr;
        char *pattern = new char[linkagePattern.str().length() + 1];
        strcpy(pattern, linkagePattern.str().c_str());
        pllLinkAlphaParameters(pattern, pllPartitions);
        pllLinkFrequencies(pattern, pllPartitions);
        pllLinkRates(pattern, pllPartitions);
        delete[] pattern;

        for (partNr = 0; partNr < pllPartitions->numberOfPartitions; partNr++) {
            pllSetFixedAlpha(alpha, partNr, pllPartitions, pllInst);
            pllSetFixedBaseFrequencies(state_freqs, 4, partNr, pllPartitions, pllInst);
            pllSetFixedSubstitutionMatrix(rate_param, 6, partNr, pllPartitions, pllInst);
        }
        delete[] rate_param;
        delete[] state_freqs;
    } else if (aln->num_states == 20) {
        double *state_freqs = new double[aln->num_states];
        getModel()->getStateFrequency(state_freqs);
        int partNr;
        for (partNr = 0; partNr < pllPartitions->numberOfPartitions; partNr++) {
            pllSetFixedAlpha(alpha, partNr, pllPartitions, pllInst);
            pllSetFixedBaseFrequencies(state_freqs, 20, partNr, pllPartitions, pllInst);
        }
        delete[] state_freqs;
    } else {
        if (params->pll) {
            outError("Phylogenetic likelihood library current does not support data type other than DNA or Protein");
        }
    }
}
double IQTree::swapTaxa(PhyloNode *node1, PhyloNode *node2) {
    assert(node1->isLeaf());
    assert(node2->isLeaf());

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

        cout << "Swapping node " << taxon1->id << " and node " << taxon2->id << endl;
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


double IQTree::doTreeSearch() {
    time_t begin_time, cur_time;
    time(&begin_time);
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
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

    if (params->write_intermediate_trees && save_all_trees != 2) {
        printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
    }

    //printTree(treels_name.c_str(), WT_NEWLINE | WT_BR_LEN);

    setRootNode(params->root);
    // keep the best tree into a string
    stringstream bestTreeStream;
    stringstream bestTopoStream;
    string perturb_tree_string;
    string imd_tree;
    printTree(bestTreeStream, WT_TAXON_ID + WT_BR_LEN);
    printTree(bestTopoStream, WT_TAXON_ID + WT_SORT_TAXA);
    string best_tree_topo = bestTopoStream.str();

    if (!params->autostop) {
        stop_rule.addImprovedIteration(1);
    }
    searchinfo.curFailedIterNum = 0;
    searchinfo.curPerStrength = params->initPerStrength;
    for (curIteration = 2; !stop_rule.meetStopCondition(curIteration); curIteration++) {
        searchinfo.curIterNum = curIteration;
        if (params->autostop) {
            if (searchinfo.curFailedIterNum == params->stopCond) {
                cout << "No better tree was found in the last " << params->stopCond
                        << " iterations. Tree search was stopped after " << curIteration << " iterations!" << endl;
                break;
            }
            if (params->adaptPert) {
               if (searchinfo.curFailedIterNum >= 50 && searchinfo.curFailedIterNum < 75) {
                    searchinfo.curPerStrength = params->initPerStrength * 2;
               } else if (searchinfo.curFailedIterNum > 75) {
                    searchinfo.curPerStrength = params->initPerStrength * 3;
                }
            }

        }
        double min_elapsed = (getCPUTime() - params->startTime) / 60;
        if (min_elapsed > params->maxtime) {
            cout << endl;
            cout << "Maximum running time of " << params->maxtime << " minutes reached" << endl;
            break;
        }
        // estimate logl_cutoff
        if (params->avoid_duplicated_trees && max_candidate_trees > 0 && treels_logl.size() > 1000) {
            int num_entries = floor(max_candidate_trees * ((double) curIteration / stop_rule.getNumIterations()));
            if (num_entries < treels_logl.size() * 0.9) {
                DoubleVector logl = treels_logl;
                nth_element(logl.begin(), logl.begin() + (treels_logl.size() - num_entries), logl.end());
                logl_cutoff = logl[treels_logl.size() - num_entries] - 1.0;
            } else
                logl_cutoff = 0.0;
            if (verbose_mode >= VB_MED) {
                if (curIteration % 10 == 0) {
                    cout << treels.size() << " trees, " << treels_logl.size() << " logls, logl_cutoff= " << logl_cutoff;
                    if (params->store_candidate_trees)
                        cout << " duplicates= " << duplication_counter << " ("
                                << (int) round(100 * ((double) duplication_counter / treels_logl.size())) << "%)"
                                << endl;
                    else
                        cout << endl;
                }
            }
        }

        if (estimate_nni_cutoff && nni_info.size() >= 500) {
            estimate_nni_cutoff = false;
            estimateNNICutoff(params);
        }

        Alignment *saved_aln = aln;

        double perturbScore;
        if (iqp_assess_quartet == IQP_BOOTSTRAP) {
            // create bootstrap sample
            Alignment* bootstrap_alignment;
            if (aln->isSuperAlignment())
                bootstrap_alignment = new SuperAlignment;
            else
                bootstrap_alignment = new Alignment;
            bootstrap_alignment->createBootstrapAlignment(aln, NULL, params->bootstrap_spec);
            setAlignment(bootstrap_alignment);
            initializeAllPartialLh();
            clearAllPartialLH();
            curScore = optimizeAllBranches();
        } else {
            if (params->snni) {
                int numNNI = floor(searchinfo.curPerStrength * (aln->getNSeq() - 3));
                //cout << numNNI << endl;
                vector<string> trees;
                for (map<double, string>::iterator it = refTreeSetSorted.begin(); it != refTreeSetSorted.end(); ++it) {
                    trees.push_back(it->second);
                }
                int index = random_int(trees.size());
                readTreeString(trees[index]);
                doRandomNNIs(numNNI);
            } else {
                doIQP();
            }
            //setAlignment(aln);
            setRootNode(params->root);
            perturb_tree_string = getTreeString();

            if (params->pll) {
                pllNewickTree *perturbTree = pllNewickParseString(perturb_tree_string.c_str());
                assert(perturbTree != NULL);
                pllTreeInitTopologyNewick(pllInst, perturbTree, PLL_FALSE);
                pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_TRUE, PLL_FALSE);
                pllOptimizeBranchLengths(pllInst, pllPartitions, params->numSmoothTree);
                pllNewickParseDestroy(&perturbTree);
                curScore = pllInst->likelihood;
                perturbScore = curScore;
            } else {
                initializeAllPartialLh();
                clearAllPartialLH();
                if (isSuperTree()) {
                    ((PhyloSuperTree*) this)->mapTrees();
                }
                curScore = optimizeAllBranches(1, TOL_LIKELIHOOD, PLL_NEWZPERCYCLE);
                perturbScore = curScore;
            }
        }

        int nni_count = 0;
        int nni_steps = 0;

        if (params->pll) {
            curScore = pllOptimizeNNI(nni_count, nni_steps, searchinfo);
            pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE,
                    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);
            imd_tree = string(pllInst->tree_string);
            readTreeString(imd_tree);
        } else {
            curScore = optimizeNNI(nni_count, nni_steps);
            imd_tree = getTreeString();
        }

        if (iqp_assess_quartet == IQP_BOOTSTRAP) {
            // restore alignment
            delete aln;
            setAlignment(saved_aln);
            initializeAllPartialLh();
            clearAllPartialLH();
        }

        if (isSuperTree()) {
            ((PhyloSuperTree*) this)->computeBranchLengths();
        }

        time(&cur_time);
        double cputime_secs = getCPUTime() - params->startTime;
        double cputime_remaining;
        if (params->maxtime < 1000000) {
            cputime_remaining = params->maxtime * 60 - cputime_secs;
        } else {
            cputime_remaining = (stop_rule.getNumIterations() - curIteration) * cputime_secs / (curIteration - 1);
        }
        cout.setf(ios::fixed, ios::floatfield);

        cout << ((iqp_assess_quartet == IQP_BOOTSTRAP) ? "Bootstrap " : "Iteration ") << curIteration
                << " / Start LogL: " << perturbScore << " / End LogL: " << curScore << " / NNIs: " << nni_count
                << " / NNI steps: " << nni_steps << " / CPU time: " << (int) round(cputime_secs) << "s";

        if (curIteration > 10 && cputime_secs > 10) {
            if (!params->autostop) {
                cout << " (" << (int) round(cputime_remaining) << "s left)";
            }
        }

        cout << endl;

        if (params->write_intermediate_trees && save_all_trees != 2) {
            printIntermediateTree(WT_NEWLINE | WT_APPEND | WT_SORT_TAXA | WT_BR_LEN);
        }

        if (curScore > bestScore) {
            stringstream cur_tree_topo_ss;
            printTree(cur_tree_topo_ss, WT_TAXON_ID | WT_SORT_TAXA);
            if (cur_tree_topo_ss.str() != best_tree_topo) {
                best_tree_topo = cur_tree_topo_ss.str();
                //cout << "Saving new better tree ..." << endl;
                if (params->snni) {
                    /****************************************** START: Optimizing model parameters ***************************************/
                    if (params->pllModOpt) {
                        assert(params->pll);
                        cout << "Optimizing model parameters by PLL ... ";
                        double stime = getCPUTime();
                        pllEvaluateLikelihood(pllInst, pllPartitions, pllInst->start, PLL_FALSE, PLL_FALSE);
                        pllOptimizeModelParameters(pllInst, pllPartitions, 1.0);
                        curScore = pllInst->likelihood;
                        double etime = getCPUTime();
                        cout << etime - stime << " seconds (logl: " << curScore << ")" << endl;
                        pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE,
                                PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
                        imd_tree = string(pllInst->tree_string);
                    } else {
                        if (params->pll) {
                            readTreeString(imd_tree);
                            initializeAllPartialLh();
                            clearAllPartialLH();
                        }
                        double *rate_param_bk = NULL;
                        if (aln->num_states == 4) {
                            rate_param_bk = new double[6];
                            getModel()->getRateMatrix(rate_param_bk);
                        }
                        double alpha_bk = getRate()->getGammaShape();
                        cout << endl;
                        cout << "Re-estimate model parameters ... " << endl;
                        double modOptScore = getModelFactory()->optimizeParameters(params->fixed_branch_length, true, 0.1);
                        if (modOptScore < curScore) {
                            cout << "  BUG: Tree logl gets worse after model optimization!" << endl;
                            cout << "  Old logl: " << curScore << " / " << "new logl: " << modOptScore << endl;
                            readTreeString(imd_tree);
                            clearAllPartialLH();
                            if (aln->num_states == 4) {
                                assert(rate_param_bk != NULL);
                                ((GTRModel*) getModel())->setRateMatrix(rate_param_bk);
                            }
                            dynamic_cast<RateGamma*>(getRate())->setGammaShape(alpha_bk);
                            getModel()->decomposeRateMatrix();
                            cout << "Reset rate parameters!" << endl;
                        } else {
                            curScore = modOptScore;
                            imd_tree = getTreeString();
                            if (params->pll) {
                                deleteAllPartialLh();
                                inputModelParam2PLL();
                                // recompute the curScore using PLL
                                curScore = inputTree2PLL(imd_tree, true);
                            }
                        }
                    }
                    /****************************************** END: Optimizing model parameters ***************************************/
                }

                if (!params->autostop) {
                    stop_rule.addImprovedIteration(curIteration);
                }
                cout << "BETTER TREE FOUND at iteration " << curIteration << ": " << curScore;
                cout << " / CPU time: " << (int) round(getCPUTime() - params->startTime) << "s" << endl << endl;

                // Only increase the number of remaining iterations if a significant improvement is found
                if (curScore - bestScore > 0.1) {
                    searchinfo.curFailedIterNum = 0;
                    searchinfo.curPerStrength = params->initPerStrength;
                }
            } else {
                cout << "UPDATE BEST LOG-LIKELIHOOD: " << curScore << endl;
                searchinfo.curFailedIterNum++;
            }
            setBestTree(imd_tree, curScore);
            if (params->write_best_trees) {
                ostringstream iter_string;
                iter_string << curIteration;
                printResultTree(iter_string.str());
            }
            printResultTree();
        } else {
            searchinfo.curFailedIterNum++;
        }

        // check whether the tree can be put into the reference set
        if (params->snni) {
            updateRefTreeSet(imd_tree, curScore);
        } else {
            // The IQPNNI algorithm
            readTreeString(bestTreeString);
        }

        if ((curIteration) % (params->step_iterations / 2) == 0 && params->gbo_replicates) {
            SplitGraph *sg = new SplitGraph;
            summarizeBootstrap(*sg);
            boot_splits.push_back(sg);
            if (params->max_candidate_trees == 0)
                max_candidate_trees = treels_logl.size() * (stop_rule.getNumIterations()) / curIteration;
            cout << "Setting tau = " << max_candidate_trees << endl;
        }

        if (curIteration == stop_rule.getNumIterations() && params->gbo_replicates && !boot_splits.empty()
                && stop_rule.getNumIterations() + params->step_iterations <= params->max_iterations) {
            //SplitGraph *sg = new SplitGraph;
            //summarizeBootstrap(*sg);
            if (!checkBootstrapStopping()) {
                if (params->max_candidate_trees == 0) {
                    max_candidate_trees = treels_logl.size() * (stop_rule.getNumIterations() + params->step_iterations)
                            / stop_rule.getNumIterations();
                }
                stop_rule.setIterationNum(stop_rule.getNumIterations() + params->step_iterations,
                        params->max_iterations);
                cout << "INFO: Increase number of iterations to " << stop_rule.getNumIterations() << " tau = "
                        << max_candidate_trees << endl;
                //delete boot_splits;
                //boot_splits = sg;
            } //else delete sg;
        }
        cout << endl;
    }

    readTreeString(bestTreeString);

    if (!params->autostop) {
        int predicted_iteration = stop_rule.getPredictedIteration();
        if (predicted_iteration > curIteration) {
            cout << endl << "WARNING: " << predicted_iteration << " iterations are needed to ensure that with a "
                    << floor(params->stop_confidence * 100) << "% confidence" << endl
                    << "         the IQPNNI search will not find a better tree" << endl;
        }
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
double IQTree::optimizeNNI(int &nni_count, int &nni_steps) {
    bool resetLamda = true; // variable indicates whether lambda should be reset
    curLambda = startLambda;
    nni_count = 0;
    int nni2apply = 0; // number of nni to be applied in each NNI steps
    int nonconf_nni = 0; // number of non-conflicting NNIs found in this round
    int MAXSTEPS = 50;
    for (nni_steps = 1; nni_steps <= MAXSTEPS; nni_steps++) {
        double oldScore = curScore;
        if (resetLamda) { // tree get improved, lamda reset
            if (save_all_trees == 2) {
                saveCurrentTree(curScore); // BQM: for new bootstrap
            }
            if (verbose_mode >= VB_DEBUG) {
                cout << "Doing NNI round " << nni_steps << endl;
                if (isSuperTree()) {
                    ((PhyloSuperTree*) this)->printMapInfo();
                }
            }

            curLambda = startLambda;
            vec_nonconf_nni.clear(); // Vector containing non-conflicting positive NNIs
            mapOptBranLens.clear(); // Vector containing branch length of the positive NNIs
            savedBranLens.clear(); // Vector containing all current branch of the tree
            posNNIs.clear(); // Vector containing all positive NNIs
            saveBranLens(); // save all current branch lengths
            initPartitionInfo(); // for super tree
            if (!nni_sort) {
                genNNIMoves(params->approximate_nni); // generate all positive NNI moves
            } else {
                genNNIMovesSort(params->approximate_nni);
            }

            /* sort all positive NNI moves (descending) */
            sort(posNNIs.begin(), posNNIs.end());
            if (verbose_mode >= VB_DEBUG) {
                cout << "curScore: " << curScore << endl;
                for (int i = 0; i < posNNIs.size(); i++) {
                    cout << "Logl of positive NNI " << i << " : " << posNNIs[i].newloglh << endl;
                }
            }

            if (posNNIs.size() == 0) {
                break;
            }

            /* remove conflicting NNIs */
            genNonconfNNIs();
            nonconf_nni = vec_nonconf_nni.size();
            if (verbose_mode >= VB_DEBUG) {
                for (int i = 0; i < vec_nonconf_nni.size(); i++) {
                    cout << "Log-likelihood of non-conflicting NNI " << i << " : " << vec_nonconf_nni[i].newloglh << endl;
                }
            }
        }
        nni2apply = floor(nonconf_nni * curLambda);
        if (nni2apply == 0)
        	nni2apply = 1;
        applyNNIs(nni2apply);

        curScore = optimizeAllBranches(1, TOL_LIKELIHOOD, PLL_NEWZPERCYCLE);

		if (verbose_mode >= VB_DEBUG) {
			cout << "logl: " << curScore << " / NNIs: " << nni2apply << endl;
		}

        if (curScore > oldScore && curScore >= vec_nonconf_nni.at(0).newloglh ) {
        	if (fabs(curScore - oldScore) < 0.001) {
        		break;
        	}
            nni_count += nni2apply;
            resetLamda = true;
        } else {

            /* tree cannot be worse if only 1 NNI is applied */
            if (nni2apply == 1) {
            	if (curScore < vec_nonconf_nni.at(0).newloglh - 0.001)
            		cout << "Error: logl=" << curScore << " < " << vec_nonconf_nni.at(0).newloglh << endl;

                // restore the tree by reverting all NNIs
                applyNNIs(nni2apply, false);
                // restore the branch lengths
                restoreAllBranLen();
                curScore = oldScore;
                break;
            }

			//if (verbose_mode >= VB_MED) {
				cout << "logl=" << curScore << " after applying " << nni2apply << " NNIs for lambda = " << curLambda
						<< " is worse than logl=" << vec_nonconf_nni.at(0).newloglh
						<< " of the best NNI. Roll back tree ..." << endl;
			//}
            curLambda = curLambda * 0.5;
            // restore the tree by reverting all NNIs
            applyNNIs(nni2apply, false);
            // restore the branch lengths
            restoreAllBranLen();
            resetLamda = false;
            curScore = oldScore;
        }
    };

    if (nni_count == 0) {
        cout << "NNI search could not find any better tree for this iteration!" << endl;
    }

    /*
     if (save_all_trees == 2 && params->nni_opt_5branches) {
     curScore = optimizeAllBranches();
     saveCurrentTree(curScore); // BQM: for new bootstrap
     saveNNITrees(); // optimize 5 branches around NNI, this makes program slower
     }*/
    //curScore = optimizeAllBranches(1);
    return curScore;
}


double IQTree::pllOptimizeNNI(int &totalNNICount, int &nniSteps, SearchInfo &searchinfo) {
    searchinfo.numAppliedNNIs = 0;
    searchinfo.curLogl = curScore;
    //cout << "curLogl: " << searchinfo.curLogl << endl;
    const int MAX_NNI_STEPS = 50;
    totalNNICount = 0;
    for (nniSteps = 1; nniSteps <= MAX_NNI_STEPS; nniSteps++) {
        searchinfo.curNumNNISteps = nniSteps;
        searchinfo.posNNIList.clear();
        double newLH = pllDoNNISearch(pllInst, pllPartitions, searchinfo);
        if (searchinfo.curNumAppliedNNIs == 0) { // no positive NNI was found
            searchinfo.curLogl = newLH;
            break;
        } else {
            searchinfo.curLogl = newLH;
            searchinfo.numAppliedNNIs += searchinfo.curNumAppliedNNIs;
        }
    }

    if (nniSteps == (MAX_NNI_STEPS + 1)) {
        cout << "WARNING: NNI search seems to run unusually too long and thus it was stopped!" << endl;
    }

    totalNNICount = searchinfo.numAppliedNNIs;
    pllInst->likelihood = searchinfo.curLogl;
    return searchinfo.curLogl;
}

void IQTree::applyNNIs(int nni2apply, bool changeBran) {
    for (int i = 0; i < nni2apply; i++) {
        doNNI(vec_nonconf_nni.at(i));
        if (!params->leastSquareNNI && changeBran) {
            // apply new branch lengths
            applyNNIBranches(vec_nonconf_nni.at(i));
        }
    }
}


void IQTree::genNonconfNNIs() {
    for (vector<NNIMove>::iterator iterMove = posNNIs.begin(); iterMove != posNNIs.end(); iterMove++) {
        bool choosen = true;
        for (vector<NNIMove>::iterator iterNextMove = vec_nonconf_nni.begin(); iterNextMove != vec_nonconf_nni.end();
                iterNextMove++) {
            if ((*iterMove).node1 == (*(iterNextMove)).node1 || (*iterMove).node2 == (*(iterNextMove)).node1
                    || (*iterMove).node1 == (*(iterNextMove)).node2 || (*iterMove).node2 == (*(iterNextMove)).node2) {
                choosen = false;
                break;
            }
        }
        if (choosen) {
            vec_nonconf_nni.push_back(*iterMove);
        }
    }
}

//double IQTree::estN95() {
//    if (vecNumNNI.size() == 0) {
//        return 0;
//    } else {
//        sort(vecNumNNI.begin(), vecNumNNI.end());
//        int index = floor(vecNumNNI.size() * speed_conf);
//        return vecNumNNI[index];
//    }
//}

double IQTree::getAvgNumNNI() {
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

//inline double IQTree::estDelta95() {
//    if (vecImpProNNI.size() == 0) {
//        return 0;
//    } else {
//        sort(vecImpProNNI.begin(), vecImpProNNI.end());
//        int index = floor(vecImpProNNI.size() * speed_conf);
//        return vecImpProNNI[index];
//    }
//}

int IQTree::getDelete() const {
    return k_delete;
}

void IQTree::setDelete(int _delete) {
    k_delete = _delete;
}

void IQTree::changeBranLen(PhyloNode *node1, PhyloNode *node2, double newlen) {
    node1->findNeighbor(node2)->length = newlen;
    node2->findNeighbor(node1)->length = newlen;
    node1->clearReversePartialLh(node2);
    node2->clearReversePartialLh(node1);
}

double IQTree::getBranLen(PhyloNode *node1, PhyloNode *node2) {
    return  node1->findNeighbor(node2)->length;
}

void IQTree::saveBranLens(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*) root;
    }
    if (dad) {
        double len = getBranLen(node, dad);
        string key = nodePair2String(node, dad);
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
        string key = nodePair2String(node, dad);
        Neighbor* bran_it = node->findNeighbor(dad);
        assert(bran_it);
        Neighbor* bran_it_back = dad->findNeighbor(node);
        assert(bran_it_back);
        assert(savedBranLens.count(key));
        bran_it->length = savedBranLens[key];
        bran_it_back->length = savedBranLens[key];
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
    string key = nodePair2String((PhyloNode*) (*it)->node, (PhyloNode*) node);
    BranLenMap::iterator bran_it = mapOptBranLens.find(key);
    if (bran_it != mapOptBranLens.end()) {
        double curlen = (*it)->length;
        changeBranLen((PhyloNode*) (*it)->node, (PhyloNode*) node, curlen + curLambda * (bran_it->second - curlen));
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
        NNIMove myMove = getBestNNIForBran(node, dad, NULL, approx_nni, params->leastSquareNNI);
        if (myMove.newloglh > curScore + params->loglh_epsilon) {
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
    assert(nodes1.size() == leafNum - 3 && nodes2.size() == leafNum - 3);

    for (i = 0; i < leafNum - 3; i++) {
        IntBranchInfo int_branch;
        PhyloNeighbor *node12_it = (PhyloNeighbor*) nodes1[i]->findNeighbor(nodes2[i]);
        //PhyloNeighbor *node21_it = (PhyloNeighbor*) nodes2[i]->findNeighbor(nodes1[i]);
        int_branch.lh_contribution = cur_lh - computeLikelihoodZeroBranch(node12_it, (PhyloNode*) nodes1[i]);
        if (int_branch.lh_contribution < 0.0)
            int_branch.lh_contribution = 0.0;
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
            NNIMove myMove = getBestNNIForBran(it->node1, it->node2, NULL, approx_nni, it->lh_contribution);
            if (myMove.newloglh > curScore) {
                addPositiveNNIMove(myMove);
                if (!estimate_nni_cutoff)
                    for (vector<IntBranchInfo>::iterator it2 = it + 1; it2 != int_branches.end(); it2++) {
                        if (it2->node1 == it->node1 || it2->node2 == it->node1 || it2->node1 == it->node2
                                || it2->node2 == it->node2)
                            it2->lh_contribution = -1.0; // do not evaluate this branch later on
                    }
            }
        } else { // otherwise, only optimize the branch length
            PhyloNode *node1 = it->node1;
            PhyloNode *node2 = it->node2;
            PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
            PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
            double stored_len = node12_it->length;
            curScore = optimizeOneBranch(node1, node2, false);
            string key("");
            if (node1->id < node2->id) {
                key += convertIntToString(node1->id) + "->" + convertIntToString(node2->id);
            } else {
                key += convertIntToString(node2->id) + "->" + convertIntToString(node1->id);
            }

            mapOptBranLens.insert(BranLenMap::value_type(key, node12_it->length));
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
                    cout << "Current lh " << cur_logl << " is much worse than expected " << treels_logl[it->second]
                            << endl;
            return;
        }
        if (verbose_mode >= VB_MAX)
            cout << "Updated logl " << treels_logl[it->second] << " to " << cur_logl << endl;
        treels_logl[it->second] = cur_logl;
        if (save_all_br_lens) {
            ostr.seekp(ios::beg);
            printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
            treels_newick[it->second] = ostr.str();
        }
        if (boot_samples.empty()) {
            computePatternLikelihood(treels_ptnlh[it->second], &cur_logl);
            return;
        }
        if (verbose_mode >= VB_MAX)
            cout << "Update treels_logl[" << tree_index << "] := " << cur_logl << endl;
    } else {
        if (logl_cutoff != 0.0 && cur_logl <= logl_cutoff + 1e-4)
            return;
        tree_index = treels_logl.size();
        if (params->store_candidate_trees)
            treels[tree_str] = tree_index;
        treels_logl.push_back(cur_logl);
        if (verbose_mode >= VB_MAX)
            cout << "Add    treels_logl[" << tree_index << "] := " << cur_logl << endl;
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
        int nsamples = boot_samples.size();

        for (int sample = 0; sample < nsamples; sample++) {
            double rell = 0.0;

            // TODO: The following parallel is not very efficient, should wrap the above loop
#ifdef _OPENMP
#pragma omp parallel for reduction(+: rell)
#endif
            for (int ptn = 0; ptn < nptn; ptn++)
                rell += pattern_lh[ptn] * boot_samples[sample][ptn];

            if (rell > boot_logl[sample] + params->ufboot_epsilon
                    || (rell > boot_logl[sample] - params->ufboot_epsilon
                            && random_double() <= 1.0 / (boot_counts[sample] + 1))) {
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
                if (rell <= boot_logl[sample] + params->ufboot_epsilon) {
                    boot_counts[sample]++;
                } else {
                    boot_counts[sample] = 1;
                }
                boot_logl[sample] = max(boot_logl[sample], rell);
                boot_trees[sample] = tree_index;
                updated++;
            } /*else if (verbose_mode >= VB_MED && rell > boot_logl[sample] - 0.01) {
             cout << "Info: multiple RELL score trees detected" << endl;
             }*/
        }
        if (updated && verbose_mode >= VB_MAX)
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
        printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
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
    int i, j;
    if (verbose_mode >= VB_MAX) {
        for (i = 0; i < trees.size(); i++)
            if (trees.tree_weights[i] > 0)
                cout << "Tree " << i + 1 << " weight= " << (double) trees.tree_weights[i] * 100 / sum_weights << endl;
    }
    int max_tree_id = max_element(trees.tree_weights.begin(), trees.tree_weights.end()) - trees.tree_weights.begin();
    cout << "max_tree_id = " << max_tree_id + 1 << "   max_weight = " << trees.tree_weights[max_tree_id];
    cout << " (" << (double) trees.tree_weights[max_tree_id] * 100 / sum_weights << "%)" << endl;
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
            cout << "WARNING: bootstrap analysis did not converge, rerun with higher number of iterations" << endl;

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
        cout << "Tree with assigned bootstrap support written to " << out_file << endl;
    }

    out_file = params.out_prefix;
    out_file += ".splits.nex";
    sg.saveFile(out_file.c_str(), IN_NEXUS, false);
    cout << "Split supports printed to NEXUS file " << out_file << endl;

    /*
     out_file = params.out_prefix;
     out_file += ".supval";
     writeInternalNodeNames(out_file);

     cout << "Support values written to " << out_file << endl;
     */

    if (params.print_ufboot_trees) {
        string filename = params.out_prefix;
        filename += ".ufboot";
        ofstream out(filename.c_str());
        for (i = 0; i < trees.size(); i++) {
            NodeVector taxa;
            // change the taxa name from ID to real name
            trees[i]->getOrderedTaxa(taxa);
            for (j = 0; j < taxa.size(); j++)
                taxa[j]->name = aln->getSeqName(taxa[j]->id);
            // now print to file
            for (j = 0; j < trees.tree_weights[i]; j++)
                trees[i]->printTree(out, WT_NEWLINE);
        }
        out.close();
        cout << "UFBoot trees printed to " << filename << endl;
    }

}

void IQTree::summarizeBootstrap(Params &params) {
    cout << "Summarizing from " << treels.size() << " candidate trees..." << endl;
    MTreeSet trees;
    IntVector tree_weights;
    int sample;
    tree_weights.resize(treels_logl.size(), 0);
    for (sample = 0; sample < boot_trees.size(); sample++)
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
    cout << split_supports_new.size() - split_supports.size() << " new splits compared to old boot_splits" << endl;
    if (split_supports_new.size() > split_supports.size())
        split_supports.resize(split_supports_new.size(), 0);

    // now compute correlation coefficient
    double corr = computeCorrelation(split_supports, split_supports_new);
    cout << "Correlation coefficient: " << corr << endl;
    // printing supports into file
    /*
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
     */
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
        printTree(iter_tree_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    } else {
        printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    }
    //printTree(tree_file_name.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH);
}

void IQTree::printResultTree(ostream &out) {
    setRootNode(params->root);
    printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
}

/*
void IQTree::printPhylolibModelParams(const char* suffix) {
    char phyloliModelFile[1024];
    strcpy(phyloliModelFile, params->out_prefix);
    strcat(phyloliModelFile, suffix);
    ofstream modelfile;
    modelfile.open(phyloliModelFile);
    for (int model = 0; model < pllInst->NumberOfModels; model++) {
        cout << "Rate parameters: ";
        for (int i = 0; i < 6; i++) {
            cout << pllInst->partitionData[model].substRates[i] << " ";
            modelfile << pllInst->partitionData[model].substRates[i] << " ";
        }
        cout << endl;
        modelfile << endl;
        cout << "Base frequencies: ";
        for (int i = 0; i < aln->num_states; i++) {
            cout << pll_tree->partitionData[model].frequencies[i] << " ";
            modelfile << pll_tree->partitionData[model].frequencies[i] << " ";
        }
        cout << endl;
        modelfile << endl;
        cout << "Gamma shape :" << pll_tree->partitionData[model].alpha << endl;
        modelfile << pll_tree->partitionData[model].alpha << endl;
    }
}
*/

void IQTree::printPhylolibTree(const char* suffix) {
    pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions, pllInst->start->back, PLL_TRUE, 1, 0, 0, 0,
            PLL_SUMMARIZE_LH, 0, 0);
    char phylolibTree[1024];
    strcpy(phylolibTree, params->out_prefix);
    strcat(phylolibTree, suffix);
    FILE *phylolib_tree = fopen(phylolibTree, "w");
    fprintf(phylolib_tree, "%s", pllInst->tree_string);
    cout << "Tree optimized by Phylolib was written to " << phylolibTree << endl;
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
                    cout << "Updated logl " << treels_logl[it->second] << " to " << curScore << endl;
                treels_logl[it->second] = curScore;
                computeLikelihood(treels_ptnlh[it->second]);
                if (save_all_br_lens) {
                    ostr.seekp(ios::beg);
                    printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
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
                    printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA | WT_BR_LEN | WT_BR_SCALE | WT_BR_LEN_ROUNDING);
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

