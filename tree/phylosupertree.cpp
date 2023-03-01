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
#include "alignment/superalignment.h"
#include "alignment/superalignmentpairwise.h"
#include "main/phylotesting.h"
#include "model/partitionmodel.h"
#include "utils/MPIHelper.h"

PhyloSuperTree::PhyloSuperTree()
 : IQTree()
{
	totalNNIs = evalNNIs = 0;
    rescale_codon_brlen = false;
	// Initialize the counter for evaluated NNIs on subtrees. FOR THIS CASE IT WON'T BE initialized.
}

PhyloSuperTree::PhyloSuperTree(SuperAlignment *alignment, bool new_iqtree) :  IQTree(alignment) {
    totalNNIs = evalNNIs = 0;

    rescale_codon_brlen = false;
    bool has_codon = false;
    vector<Alignment*>::iterator it;
    for (it = alignment->partitions.begin(); it != alignment->partitions.end(); it++)
        if ((*it)->seq_type != SEQ_CODON) {
            rescale_codon_brlen = true;
        } else
            has_codon = true;
    
    rescale_codon_brlen &= has_codon;

    // Initialize the counter for evaluated NNIs on subtrees
    int part = 0;

    StrVector model_names;
    for (it = alignment->partitions.begin(); it != alignment->partitions.end(); it++, part++) {
        PhyloTree *tree;
        if (new_iqtree)
            tree = new IQTree(*it);
        else
            tree = new PhyloTree(*it);
        push_back(tree);
        PartitionInfo info;
        info.cur_ptnlh = NULL;
        info.nniMoves[0].ptnlh = NULL;
        info.nniMoves[1].ptnlh = NULL;
        info.evalNNIs = 0.0;
        part_info.push_back(info);
    }
    
    aln = alignment;
    
}

PhyloSuperTree::PhyloSuperTree(SuperAlignment *alignment, PhyloSuperTree *super_tree) :  IQTree(alignment) {
	totalNNIs = evalNNIs = 0;
    rescale_codon_brlen = super_tree->rescale_codon_brlen;
	part_info = super_tree->part_info;
	for (vector<Alignment*>::iterator it = alignment->partitions.begin(); it != alignment->partitions.end(); it++) {
		PhyloTree *tree = new PhyloTree((*it));
		push_back(tree);
	}
	// Initialize the counter for evaluated NNIs on subtrees
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].evalNNIs = 0.0;
	}

	aln = alignment;
}

void PhyloSuperTree::setModelFactory(ModelFactory *model_fac) {
    PhyloTree::setModelFactory(model_fac);
    if (model_fac) {
        PhyloSuperTree *tree = (PhyloSuperTree*)model_fac->site_rate->phylo_tree;
        for (int part = 0; part != size(); part++) {
            at(part)->setModelFactory(tree->at(part)->getModelFactory());
        }
    } else {
        for (int part = 0; part != size(); part++) {
            at(part)->setModelFactory(NULL);
        }
    }
}

void PhyloSuperTree::setPartInfo(PhyloSuperTree *tree) {
    part_info = tree->part_info;
    int part = 0;
    for (iterator it = begin(); it != end(); it++, part++) {
        part_info[part].evalNNIs = 0.0;
    }
    if (!params->bootstrap_spec) {
        return;
    }
    // 2018-06-03: for -bsam GENE, number of partitions might be reduced
    if (part_info.size() <= size())
        return;
    part_info.erase(part_info.begin()+size(), part_info.end());
    for (part = 0; part < part_info.size(); part++) {
        bool found = false;
        for (int p = 0; p < tree->size(); p++)
            if (tree->at(p)->aln->name == at(part)->aln->name) {
                part_info[part] = tree->part_info[p];
                part_info[part].evalNNIs = 0.0;
                found = true;
                break;
            }
        ASSERT(found);
    }
}


void PhyloSuperTree::setSuperAlignment(Alignment *alignment) {
    PhyloTree::setAlignment(alignment);

    SuperAlignment *saln = (SuperAlignment*)aln;
    for (int i = 0; i < size(); i++)
        at(i)->setAlignment(saln->partitions.at(i));
}

void PhyloSuperTree::setCheckpoint(Checkpoint *checkpoint) {
	IQTree::setCheckpoint(checkpoint);
	for (iterator it = begin(); it != end(); it++)
		(*it)->setCheckpoint(checkpoint);
}

void PhyloSuperTree::saveCheckpoint() {
//    checkpoint->startStruct("PhyloSuperTree");
//    int part = 0;
//    for (iterator it = begin(); it != end(); it++, part++) {
//    	string key = part_info[part].name + ".tree";
//    	checkpoint->put(key, (*it)->getTreeString());
//    }
//    checkpoint->endStruct();
    IQTree::saveCheckpoint();
}

void PhyloSuperTree::restoreCheckpoint() {
    IQTree::restoreCheckpoint();
    
    // first get the newick string of super tree
//    checkpoint->startStruct("PhyloTree");
//    string newick;
//    CKP_RESTORE(newick);
//    checkpoint->endStruct();
//
//    if (newick.empty()) return;
//    
//    // now get partition tree strings
//    checkpoint->startStruct("PhyloSuperTree");
//    int part = 0;
//    for (iterator it = begin(); it != end(); it++, part++) {
//    	string key = part_info[part].name + ".tree";
//    	string part_tree;
//    	if (!checkpoint->get(key, part_tree))
//    		outError("No tree for partition " + part_info[part].name + " found from checkpoint");
//    	newick += part_tree;
//    }
//
//    checkpoint->endStruct();
//
//    readTreeString(newick);

}

void PhyloSuperTree::setParams(Params* params) {
	IQTree::setParams(params);
	for (iterator it = begin(); it != end(); it++) {
		(*it)->setParams(params);
	}
}

void PhyloSuperTree::initSettings(Params &params) {
    IQTree::initSettings(params);
    setLikelihoodKernel(params.SSE);
    setNumThreads(params.num_threads);
    for (iterator it = begin(); it != end(); it++) {
        (*it)->params = &params;
        (*it)->optimize_by_newton = params.optimize_by_newton;
    }
}

void PhyloSuperTree::initSequences(Node* node, Node* dad)
{
    // init sequences for the primary/super tree first
    PhyloTree::initSequences();
    
    // init sequences for each partition trees
    for (iterator it = begin(); it != end(); it++) {
        (*it)->PhyloTree::initSequences();
    }
}

void PhyloSuperTree::setLikelihoodKernel(LikelihoodKernel lk) {
    PhyloTree::setLikelihoodKernel(lk);
    for (iterator it = begin(); it != end(); it++)
        (*it)->setLikelihoodKernel(lk);
}

void PhyloSuperTree::setParsimonyKernel(LikelihoodKernel lk) {
    PhyloTree::setParsimonyKernel(lk);
    for (iterator it = begin(); it != end(); it++)
        (*it)->setParsimonyKernel(lk);
}

void PhyloSuperTree::changeLikelihoodKernel(LikelihoodKernel lk) {
	PhyloTree::changeLikelihoodKernel(lk);
}

void PhyloSuperTree::setNumThreads(int num_threads) {
    PhyloTree::setNumThreads((size() >= num_threads) ? num_threads : 1);
    for (iterator it = begin(); it != end(); it++)
        (*it)->setNumThreads((size() >= num_threads) ? 1 : num_threads);
}

void PhyloSuperTree::printResultTree(string suffix) {
    if (MPIHelper::getInstance().isWorker()) {
        return;
    }
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    
    IQTree::printResultTree(suffix);
    
    string tree_file_name = params->out_prefix;
    tree_file_name += ".parttrees";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(tree_file_name);
        for (iterator it = begin(); it != end(); it++)
            (*it)->printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, tree_file_name);
    }
    if (verbose_mode >= VB_MED)
        cout << "Partition trees printed to " << tree_file_name << endl;
}

string PhyloSuperTree::getTreeString() {
	stringstream tree_stream;
	printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
	for (iterator it = begin(); it != end(); it++)
		(*it)->printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
	return tree_stream.str();
}

void PhyloSuperTree::readTreeString(const string &tree_string) {
	stringstream str;
	str << tree_string;
	str.seekg(0, ios::beg);
	freeNode();
    // bug fix 2017-11-30: in case taxon name happens to be ID
	MTree::readTree(str, rooted);
    assignLeafNames();
//	setAlignment(aln);
	setRootNode(params->root);
	for (iterator it = begin(); it != end(); it++) {
		(*it)->freeNode();
		(*it)->readTree(str, (*it)->rooted);
        (*it)->assignLeafNames();
//		(*it)->setAlignment((*it)->aln);
	}
	linkTrees();
//	if (isSuperTree()) {
//		((PhyloSuperTree*) this)->mapTrees();
//	}
	if (params->pll) {
		ASSERT(0);
		pllReadNewick(getTreeString());
	}
	resetCurScore();

}


/**
 * save branch lengths into a vector
 */
void PhyloSuperTree::saveBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
    ASSERT(getMixlen() == 1); // supertree and treemixlen not allowed together
	int totalBranchNum = branchNum * getMixlen();
	iterator it;
	for (it = begin(); it != end(); it++) {
		totalBranchNum += (*it)->branchNum * (*it)->getMixlen();
	}
	lenvec.resize(startid + totalBranchNum);

	PhyloTree::saveBranchLengths(lenvec, startid);
	startid += branchNum * getMixlen();
	for (iterator it = begin(); it != end(); it++) {
		(*it)->saveBranchLengths(lenvec, startid);
		startid += (*it)->branchNum * (*it)->getMixlen();
	}
}
/**
 * restore branch lengths from a vector previously called with saveBranchLengths
 */
void PhyloSuperTree::restoreBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
	PhyloTree::restoreBranchLengths(lenvec, startid);
	startid += branchNum * getMixlen();
	for (iterator it = begin(); it != end(); it++) {
		(*it)->restoreBranchLengths(lenvec, startid);
		startid += (*it)->branchNum * (*it)->getMixlen();
	}
}

int PhyloSuperTree::collapseInternalBranches(Node *node, Node *dad, double threshold) {
	if (!node) node = root;
    int count = 0;
	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		count += collapseInternalBranches((*it)->node, node, threshold);
	}
	if (node->isLeaf()) {
		return count;
	}
	NeighborVec nei_vec;
	nei_vec.insert(nei_vec.begin(), node->neighbors.begin(), node->neighbors.end());
	for (it = nei_vec.begin(); it != nei_vec.end(); it++) {
		if ((*it)->node != dad && !(*it)->node->isLeaf() && (*it)->length <= threshold) {

			//delete branch of gene-trees
			SuperNeighbor* snei = (SuperNeighbor*)(*it);
			int part = 0;
			for (part = 0; part != size(); part++) {
				if (snei->link_neighbors[part]) {
					SuperNeighbor* snei_back = (SuperNeighbor*)(*it)->node->findNeighbor(node);
					at(part)->removeNode(snei_back->link_neighbors[part]->node, snei->link_neighbors[part]->node);
				}
			}
			// delete the child node
			removeNode(node, (*it)->node);
			count++;
		}
	}
    return count;
}


Node* PhyloSuperTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

Node* PhyloSuperTree::newNode(int node_id, int node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

size_t PhyloSuperTree::getAlnNPattern() {
	size_t num = 0;
	for (iterator it = begin(); it != end(); it++)
		num += (*it)->getAlnNPattern();
	return num;
}

size_t PhyloSuperTree::getAlnNSite() {
	size_t num = 0;
	for (iterator it = begin(); it != end(); it++)
		num += (*it)->getAlnNSite();
	return num;
}

double PhyloSuperTree::computeDist(int seq1, int seq2, double initial_dist, double &var) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0) {
    	if (params->compute_obs_dist)
            initial_dist = aln->computeObsDist(seq1, seq2);
    	else
    		initial_dist = aln->computeDist(seq1, seq2);
    }
    if (initial_dist == MAX_GENETIC_DIST) return initial_dist; // MANUEL: here no d2l is return
    if (!model_factory || !site_rate) return initial_dist; // MANUEL: here no d2l is return

    // now optimize the distance based on the model and site rate
    SuperAlignmentPairwise aln_pair(this, seq1, seq2);
    return aln_pair.optimizeDist(initial_dist, var);
}

void PhyloSuperTree::linkBranch(int part, SuperNeighbor *nei, SuperNeighbor *dad_nei) {
	SuperNode *node = (SuperNode*)dad_nei->node;
	SuperNode *dad = (SuperNode*)nei->node;
	nei->link_neighbors[part] = NULL;
	dad_nei->link_neighbors[part] = NULL;
	vector<PhyloNeighbor*> part_vec;
	vector<PhyloNeighbor*> child_part_vec;

	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		if (((SuperNeighbor*)*it)->link_neighbors[part]) {
			part_vec.push_back(((SuperNeighbor*)*it)->link_neighbors[part]);
			child_part_vec.push_back(((SuperNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]);
			ASSERT(child_part_vec.back()->node == child_part_vec.front()->node || child_part_vec.back()->id == child_part_vec.front()->id);
		}
	}

	if (part_vec.empty())
		return;
	if (part_vec.size() == 1) {
		nei->link_neighbors[part] = child_part_vec[0];
		dad_nei->link_neighbors[part] = part_vec[0];
		return;
	}
	if (part_vec[0] == child_part_vec[1]) {
		// ping-pong, out of sub-tree
		ASSERT(part_vec[1] == child_part_vec[0]);
		return;
	}
	PhyloNode *node_part = (PhyloNode*) child_part_vec[0]->node;
	PhyloNode *dad_part = NULL;
	FOR_NEIGHBOR(node_part, NULL, it) {
		bool appear = false;
		for (vector<PhyloNeighbor*>::iterator it2 = part_vec.begin(); it2 != part_vec.end(); it2++){
			if ((*it2) == (*it)) {
				appear = true; break;
			}
		}
		if (!appear) {
			ASSERT(!dad_part);
			dad_part = (PhyloNode*)(*it)->node;
		}
	}
	nei->link_neighbors[part] = (PhyloNeighbor*)node_part->findNeighbor(dad_part);
	dad_nei->link_neighbors[part] = (PhyloNeighbor*)dad_part->findNeighbor(node_part);
}

void PhyloSuperTree::linkTree(int part, NodeVector &part_taxa, SuperNode *node, SuperNode *dad) {
	if (!node) {
		if (!root->isLeaf())
			node = (SuperNode*) root;
		else
			node = (SuperNode*)root->neighbors[0]->node;
		ASSERT(node);
		if (node->isLeaf()) // two-taxa tree
			dad = (SuperNode*)node->neighbors[0]->node;
	}
	SuperNeighbor *nei = NULL;
	SuperNeighbor *dad_nei = NULL;
	if (dad) {
		nei = (SuperNeighbor*)node->findNeighbor(dad);
		dad_nei = (SuperNeighbor*)dad->findNeighbor(node);
		if (nei->link_neighbors.empty()) nei->link_neighbors.resize(size());
		if (dad_nei->link_neighbors.empty()) dad_nei->link_neighbors.resize(size());
		nei->link_neighbors[part] = NULL;
		dad_nei->link_neighbors[part] = NULL;
	}
	if (node->isLeaf()) {
		ASSERT(dad);
		PhyloNode *node_part = (PhyloNode*)part_taxa[node->id];
		if (node_part) {
			PhyloNode *dad_part = (PhyloNode*)node_part->neighbors[0]->node;
			ASSERT(node_part->isLeaf());
			nei->link_neighbors[part] = (PhyloNeighbor*) node_part->neighbors[0];
			dad_nei->link_neighbors[part] = (PhyloNeighbor*)dad_part->findNeighbor(node_part);
		}
		return;
	}

	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		linkTree(part, part_taxa, (SuperNode*) (*it)->node, (SuperNode*) node);
	}
	if (!dad) return;
	linkBranch(part, nei, dad_nei);
}

void PhyloSuperTree::syncRooting() {
    if (empty() || !front()->root)
        return;
    // check if all trees are similarly rooted or unrooted
    bool same_rooting = true;
    iterator tree;
    for (tree = begin(); tree != end(); tree++) {
        if ((*tree)->rooted != front()->rooted) {
            same_rooting = false;
            break;
        }
    }
    
    if (same_rooting) {
        if (!empty() && rooted != front()->rooted) {
            if (rooted)
                convertToUnrooted();
            else
                convertToRooted();
        }
    } else if (!rooted) {
        // not same rooting between trees
        convertToRooted();
    }
}

void PhyloSuperTree::printMapInfo() {
	NodeVector nodes1, nodes2;
	getBranches(nodes1, nodes2);
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		cout << "Subtree for partition " << part << endl;
		(*it)->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
		for (int i = 0; i < nodes1.size(); i++) {
			PhyloNeighbor *nei1 = ((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
			PhyloNeighbor *nei2 = ((SuperNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
			cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
			if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
			cout << ",";
			if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
			cout << " -> ";
			if (nei2) {
				cout << nei2->id << ":";
				if (nei2->node->isLeaf())
					cout << nei2->node->name;
				else cout << nei2->node->id;
			}
			else cout << -1;
			cout << ",";
			if (nei1)
				if (nei1->node->isLeaf())
					cout << nei1->node->name;
				else cout << nei1->node->id;
			else cout << -1;
			cout << endl;
		}
	}
}


void PhyloSuperTree::mapTrees() {
	ASSERT(root);
    syncRooting();
	int part = 0, i;
	if (verbose_mode >= VB_DEBUG)
		drawTree(cout,  WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
	for (iterator it = begin(); it != end(); it++, part++) {
		string taxa_set;
        Pattern taxa_pat = aln->getPattern(part);
        taxa_set.insert(taxa_set.begin(), taxa_pat.begin(), taxa_pat.end());
		(*it)->copyTree(this, taxa_set);
        if ((*it)->getModel()) {
			(*it)->initializeAllPartialLh();
        }
        (*it)->resetCurScore();
		NodeVector my_taxa, part_taxa;
		(*it)->getOrderedTaxa(my_taxa);
		part_taxa.resize(leafNum, NULL);
		for (i = 0; i < leafNum; i++) {
            int id;
            if (i < aln->getNSeq())
                id = ((SuperAlignment*)aln)->taxa_index[i][part];
            else {
                ASSERT(rooted);
                if ((*it)->rooted)
                    id = (*it)->leafNum-1;
                else
                    id = -1;
            }
			if (id >=0) part_taxa[i] = my_taxa[id];
		}
		if (verbose_mode >= VB_DEBUG) {
			cout << "Subtree for partition " << part << endl;
			(*it)->drawTree(cout,  WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
		}
		linkTree(part, part_taxa);
	}

	if (verbose_mode >= VB_DEBUG) printMapInfo();
}

void PhyloSuperTree::linkTrees() {
	int part = 0;
	iterator it;
	for (it = begin(), part = 0; it != end(); it++, part++) {
		(*it)->initializeTree();
		(*it)->setAlignment((*it)->aln);
        if ((*it)->getModel()) {
			(*it)->initializeAllPartialLh();
        }
        (*it)->resetCurScore();
		NodeVector my_taxa, part_taxa;
		(*it)->getOrderedTaxa(my_taxa);
		part_taxa.resize(leafNum, NULL);
		int i;
		for (i = 0; i < leafNum; i++) {
            int id;
            if (i < aln->getNSeq())
                id = ((SuperAlignment*)aln)->taxa_index[i][part];
            else if ((*it)->rooted)
                id = (*it)->leafNum-1;
            else
                id = -1;
			if (id >=0) part_taxa[i] = my_taxa[id];
		}
		linkTree(part, part_taxa);
	}
}

void PhyloSuperTree::initializeAllPartialLh() {
	for (iterator it = begin(); it != end(); it++) {
		(*it)->initializeAllPartialLh();
	}
}


void PhyloSuperTree::deleteAllPartialLh() {
	for (iterator it = begin(); it != end(); it++) {
		(*it)->deleteAllPartialLh();
	}
}

void PhyloSuperTree::clearAllPartialLH(bool make_null) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->clearAllPartialLH(make_null);
    }
}

int PhyloSuperTree::computeParsimonyBranchObsolete(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst) {
    int score = 0, part = 0;
    SuperNeighbor *dad_nei = (SuperNeighbor*)dad_branch;
    SuperNeighbor *node_nei = (SuperNeighbor*)(dad_branch->node->findNeighbor(dad));
        
    if (branch_subst)
        branch_subst = 0;
    for (iterator it = begin(); it != end(); it++, part++) {
        int this_subst = 0;
        if (dad_nei->link_neighbors[part]) {
            if (branch_subst)
                score += (*it)->computeParsimonyBranch(dad_nei->link_neighbors[part], (PhyloNode*)node_nei->link_neighbors[part]->node, &this_subst);
            else
                score += (*it)->computeParsimonyBranch(dad_nei->link_neighbors[part], (PhyloNode*)node_nei->link_neighbors[part]->node);
        } else
            score += (*it)->computeParsimony();
        if (branch_subst)
            branch_subst += this_subst;
    }
    return score;
}

void PhyloSuperTree::computePartitionOrder() {
    if (!part_order.empty())
        return;
    int i, ntrees = size();
    part_order.resize(ntrees);
    part_order_by_nptn.resize(ntrees);
#ifdef _OPENMP
    int *id = new int[ntrees];
    double *cost = new double[ntrees];
    
    for (i = 0; i < ntrees; i++) {
        Alignment *part_aln = at(i)->aln;
        cost[i] = -((double)part_aln->getNSeq())*part_aln->getNPattern()*part_aln->num_states;
        id[i] = i;
    }
    quicksort(cost, 0, ntrees-1, id);
    for (i = 0; i < ntrees; i++) 
        part_order[i] = id[i];
        
    // compute part_order by number of patterns
    for (i = 0; i < ntrees; i++) {
        Alignment *part_aln = at(i)->aln;
        cost[i] = -((double)part_aln->getNPattern())*part_aln->num_states;
        id[i] = i;
    }
    quicksort(cost, 0, ntrees-1, id);
    for (i = 0; i < ntrees; i++) 
        part_order_by_nptn[i] = id[i];
        
    delete [] cost;
    delete [] id;
    
    if (verbose_mode >= VB_MED) {
        cout << "Partitions ordered by computation costs:" << endl;
        cout << "#nexus" << endl << "begin sets;" << endl;
        for (i = 0; i < ntrees; i++)
            cout << "  charset " << at(part_order[i])->aln->name << " = " << at(part_order[i])->aln->position_spec << ";" << endl;
        cout << "end;" << endl;
    }
#else
    for (i = 0; i < ntrees; i++) {
        part_order[i] = i;
        part_order_by_nptn[i] = i;
    }
#endif // OPENMP
}

double PhyloSuperTree::computeLikelihood(double *pattern_lh) {
	double tree_lh = 0.0;
	int ntrees = size();
	if (pattern_lh) {
		//#ifdef _OPENMP
		//#pragma omp parallel for reduction(+: tree_lh)
		//#endif
		for (int i = 0; i < ntrees; i++) {
			part_info[i].cur_score = at(i)->computeLikelihood(pattern_lh);
			tree_lh += part_info[i].cur_score;
			pattern_lh += at(i)->getAlnNPattern();
		}
	} else {
        if (part_order.empty()) computePartitionOrder();
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(num_threads > 1)
		#endif
		for (int j = 0; j < ntrees; j++) {
            int i = part_order[j];
			part_info[i].cur_score = at(i)->computeLikelihood();
			tree_lh += part_info[i].cur_score;
		}
	}
	return tree_lh;
}

int PhyloSuperTree::getNumLhCat(SiteLoglType wsl) {
    int ncat = 0, it_cat;
    for (iterator it = begin(); it != end(); it++)
        if ((it_cat = (*it)->getNumLhCat(wsl)) > ncat)
            ncat = it_cat;
    return ncat;
}


void PhyloSuperTree::computePatternLikelihood(double *pattern_lh, double *cur_logl, double *ptn_lh_cat, SiteLoglType wsl) {
	size_t offset = 0, offset_lh_cat = 0;
	iterator it;
	for (it = begin(); it != end(); it++) {
		if (ptn_lh_cat)
			(*it)->computePatternLikelihood(pattern_lh + offset, NULL, ptn_lh_cat + offset_lh_cat, wsl);
		else
			(*it)->computePatternLikelihood(pattern_lh + offset);
		offset += (*it)->aln->getNPattern();
        offset_lh_cat += (*it)->aln->getNPattern() * (*it)->getNumLhCat(wsl);
	}
	if (cur_logl) { // sanity check
		double sum_logl = 0;
		offset = 0;
		for (it = begin(); it != end(); it++) {
			int nptn = (*it)->aln->getNPattern();
			for (int j = 0; j < nptn; j++)
				sum_logl += pattern_lh[offset + j] * (*it)->aln->at(j).frequency;
			offset += (*it)->aln->getNPattern();
		}
		if (fabs(sum_logl - *cur_logl) > 0.001) {
            cout << *cur_logl << " " << sum_logl << endl;
//            outError("Wrong PhyloSuperTree::", __func__);
		}
        ASSERT(fabs(sum_logl - *cur_logl) < 0.001);
	}
}

void PhyloSuperTree::computePatternProbabilityCategory(double *ptn_prob_cat, SiteLoglType wsl) {
	size_t offset = 0;
	for (iterator it = begin(); it != end(); it++) {
        (*it)->computePatternProbabilityCategory(ptn_prob_cat + offset, wsl);
        offset += (*it)->aln->getNPattern() * (*it)->getNumLhCat(wsl);
	}
}

double PhyloSuperTree::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
	double tree_lh = 0.0;
	int ntrees = size();
    if (part_order.empty()) computePartitionOrder();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(num_threads > 1)
	#endif
	for (int j = 0; j < ntrees; j++) {
        int i = part_order[j];
		part_info[i].cur_score = at(i)->optimizeAllBranches(my_iterations, tolerance/min(ntrees,10), maxNRStep);
		tree_lh += part_info[i].cur_score;
		if (verbose_mode >= VB_MAX)
			at(i)->printTree(cout, WT_BR_LEN + WT_NEWLINE);
	}

	if (my_iterations >= 100) computeBranchLengths();
	return tree_lh;
}

PhyloSuperTree::~PhyloSuperTree()
{
	for (vector<PartitionInfo>::reverse_iterator pit = part_info.rbegin(); pit != part_info.rend(); pit++) {
		if (pit->nniMoves[1].ptnlh)
			delete [] pit->nniMoves[1].ptnlh;
		pit->nniMoves[1].ptnlh = NULL;
		if (pit->nniMoves[0].ptnlh)
			delete [] pit->nniMoves[0].ptnlh;
		pit->nniMoves[0].ptnlh = NULL;
		if (pit->cur_ptnlh)
			delete [] pit->cur_ptnlh;
		pit->cur_ptnlh = NULL;
	}
	part_info.clear();

	for (reverse_iterator it = rbegin(); it != rend(); it++)
		delete (*it);
	clear();
}


void PhyloSuperTree::initPartitionInfo() {
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].cur_score = 0.0;

		part_info[part].cur_brlen.resize((*it)->branchNum);
		if (params->nni5) {
			part_info[part].nni1_brlen.resize((*it)->branchNum * 5);
			part_info[part].nni2_brlen.resize((*it)->branchNum * 5);
		} else {
			part_info[part].nni1_brlen.resize((*it)->branchNum);
			part_info[part].nni2_brlen.resize((*it)->branchNum);
		}

		(*it)->getBranchLengths(part_info[part].cur_brlen);

		if (save_all_trees == 2 || params->write_intermediate_trees >= 2) {
			// initialize ptnlh for ultrafast bootstrap
			int nptn = (*it)->getAlnNPattern();
			if (!part_info[part].cur_ptnlh)
				part_info[part].cur_ptnlh = new double[nptn];
			if (!part_info[part].nniMoves[0].ptnlh)
				part_info[part].nniMoves[0].ptnlh = new double [nptn];
			if (!part_info[part].nniMoves[1].ptnlh)
				part_info[part].nniMoves[1].ptnlh = new double [nptn];
		} else {
            part_info[part].cur_ptnlh = NULL;
            part_info[part].nniMoves[0].ptnlh = NULL;
            part_info[part].nniMoves[1].ptnlh = NULL;
        }
	}
}

int PhyloSuperTree::getMaxPartNameLength() {
	int namelen = 0;
	for (iterator it = begin(); it != end(); it++)
		namelen = max((int)(*it)->aln->name.length(), namelen);
	return namelen;
}

NNIMove PhyloSuperTree::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves) {
    if (((PhyloNeighbor*)node1->findNeighbor(node2))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        PhyloNode *tmp = node1;
        node1 = node2;
        node2 = tmp;
    }
    NNIMove myMove;
    //myMove.newloglh = 0;
	SuperNeighbor *nei1 = ((SuperNeighbor*)node1->findNeighbor(node2));
	SuperNeighbor *nei2 = ((SuperNeighbor*)node2->findNeighbor(node1));
	ASSERT(nei1 && nei2);
	SuperNeighbor *node1_nei = NULL;
	SuperNeighbor *node2_nei = NULL;
	SuperNeighbor *node2_nei_other = NULL;
	FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)
    if (((PhyloNeighbor*)*node1_it)->direction != TOWARD_ROOT)
    {
		node1_nei = (SuperNeighbor*)(*node1_it);
		break;
	}
	FOR_NEIGHBOR_DECLARE(node2, node1, node2_it) {
		node2_nei = (SuperNeighbor*)(*node2_it);
		break;
	}

	FOR_NEIGHBOR_IT(node2, node1, node2_it_other)
	if ((*node2_it_other) != node2_nei) {
		node2_nei_other = (SuperNeighbor*)(*node2_it_other);
		break;
	}

    // check for compatibility with constraint tree
    bool nni_ok[2] = {true, true};
    int nniid = 0;
	FOR_NEIGHBOR(node2, node1, node2_it) {
        NNIMove nni;
        nni.node1 = node1;
        nni.node2 = node2;
        nni.node1Nei_it = node1->findNeighborIt(node1_nei->node);
        nni.node2Nei_it = node2_it;
        nni_ok[nniid++] = constraintTree.isCompatible(nni);
    }
    ASSERT(nniid == 2);
    myMove.node1 = myMove.node2 = NULL;
    myMove.newloglh = -DBL_MAX;
    // return if both NNIs do not satisfy constraint
    if (!nni_ok[0] && !nni_ok[1]) {
//        ASSERT(!nniMoves);
        if (nniMoves) {
            nniMoves[0].newloglh = nniMoves[1].newloglh = -DBL_MAX;
        }
        return myMove;
    }

	//double bestScore = optimizeOneBranch(node1, node2, false);

	int ntrees = size(), part;
	double nni_score1 = 0.0, nni_score2 = 0.0;
	int local_totalNNIs = 0, local_evalNNIs = 0;

    if (part_order.empty()) computePartitionOrder();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: nni_score1, nni_score2, local_totalNNIs, local_evalNNIs) private(part) schedule(dynamic) if(num_threads>1)
	#endif
	for (int treeid = 0; treeid < ntrees; treeid++) {
        part = part_order_by_nptn[treeid];
		bool is_nni = true;
		local_totalNNIs++;
		FOR_NEIGHBOR_DECLARE(node1, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		FOR_NEIGHBOR(node2, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		if (!is_nni && params->terrace_aware) {
			if (part_info[part].cur_score == 0.0)  {
				part_info[part].cur_score = at(part)->computeLikelihood();
				if (save_all_trees == 2 || nniMoves)
					at(part)->computePatternLikelihood(part_info[part].cur_ptnlh, &part_info[part].cur_score);
			}
			nni_score1 += part_info[part].cur_score;
			nni_score2 += part_info[part].cur_score;
			continue;
		}

		local_evalNNIs++;
		part_info[part].evalNNIs++;

		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];

		int brid = nei1_part->id;

		//NNIMove part_moves[2];
		//part_moves[0].node1Nei_it = NULL;

		// setup subtree NNI correspondingly
		PhyloNode *node1_part = (PhyloNode*)nei2_part->node;
		PhyloNode *node2_part = (PhyloNode*)nei1_part->node;
		part_info[part].nniMoves[0].node1 = part_info[part].nniMoves[1].node1 = node1;
		part_info[part].nniMoves[0].node2 = part_info[part].nniMoves[1].node2 = node2;
		part_info[part].nniMoves[0].node1Nei_it = node1_part->findNeighborIt(node1_nei->link_neighbors[part]->node);
		part_info[part].nniMoves[0].node2Nei_it = node2_part->findNeighborIt(node2_nei->link_neighbors[part]->node);

		part_info[part].nniMoves[1].node1Nei_it = node1_part->findNeighborIt(node1_nei->link_neighbors[part]->node);
		part_info[part].nniMoves[1].node2Nei_it = node2_part->findNeighborIt(node2_nei_other->link_neighbors[part]->node);

		at(part)->getBestNNIForBran((PhyloNode*)nei2_part->node, (PhyloNode*)nei1_part->node, part_info[part].nniMoves);
		// detect the corresponding NNIs and swap if necessary (the swapping refers to the swapping of NNI order)
		if (!((*part_info[part].nniMoves[0].node1Nei_it == node1_nei->link_neighbors[part] &&
				*part_info[part].nniMoves[0].node2Nei_it == node2_nei->link_neighbors[part]) ||
			(*part_info[part].nniMoves[0].node1Nei_it != node1_nei->link_neighbors[part] &&
					*part_info[part].nniMoves[0].node2Nei_it != node2_nei->link_neighbors[part])))
		{
			outError("WRONG");
			NNIMove tmp = part_info[part].nniMoves[0];
			part_info[part].nniMoves[0] = part_info[part].nniMoves[1];
			part_info[part].nniMoves[1] = tmp;
		}
		nni_score1 += part_info[part].nniMoves[0].newloglh;
		nni_score2 += part_info[part].nniMoves[1].newloglh;
		int numlen = 1;
		if (params->nni5) numlen = 5;
		for (int i = 0; i < numlen; i++) {
			part_info[part].nni1_brlen[brid*numlen + i] = part_info[part].nniMoves[0].newLen[i];
			part_info[part].nni2_brlen[brid*numlen + i] = part_info[part].nniMoves[1].newLen[i];
		}

	}
	totalNNIs += local_totalNNIs;
	evalNNIs += local_evalNNIs;
	double nni_scores[2] = {nni_score1, nni_score2};
    
    if (!nni_ok[0]) nni_scores[0] = -DBL_MAX;
    if (!nni_ok[1]) nni_scores[1] = -DBL_MAX;

	myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
	myMove.node1 = node1;
	myMove.node2 = node2;
	if (nni_scores[0] > nni_scores[1]) {
		myMove.swap_id = 1;
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei->node);
		myMove.newloglh = nni_scores[0];
	} else  {
		myMove.swap_id = 2;
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei_other->node);
		myMove.newloglh = nni_scores[1];
	}

	if (save_all_trees != 2 && !nniMoves) return myMove;

	// for bootstrap now
    //now setup pattern likelihoods per partition
	double *save_lh_factor = new double [ntrees];
	double *save_lh_factor_back = new double [ntrees];
	nniid = 0;
	FOR_NEIGHBOR(node2, node1, node2_it) if (nni_ok[nniid]) 
    {

		// do the NNI
		node2_nei = (SuperNeighbor*)(*node2_it);
        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        for (part = 0; part < ntrees; part++) {
			bool is_nni = true;
			FOR_NEIGHBOR_DECLARE(node1, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
			}
			FOR_NEIGHBOR(node2, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
			}
			if (!is_nni)
				memcpy(at(part)->_pattern_lh, part_info[part].cur_ptnlh, at(part)->getAlnNPattern() * sizeof(double));
			else
				memcpy(at(part)->_pattern_lh, part_info[part].nniMoves[nniid].ptnlh, at(part)->getAlnNPattern() * sizeof(double));
    		save_lh_factor[part] = at(part)->current_it->lh_scale_factor;
    		save_lh_factor_back[part] = at(part)->current_it_back->lh_scale_factor;
    		at(part)->current_it->lh_scale_factor = 0.0;
    		at(part)->current_it_back->lh_scale_factor = 0.0;
        }
        if (nniMoves) {
        	nniMoves[nniid].newloglh = nni_scores[nniid];
       		computePatternLikelihood(nniMoves[nniid].ptnlh, &nni_scores[nniid]);
        }
        if (save_all_trees == 2)
        	saveCurrentTree(nni_scores[nniid]);

        // restore information
        for (part = 0; part < ntrees; part++) {
    		at(part)->current_it->lh_scale_factor = save_lh_factor[part];
    		at(part)->current_it_back->lh_scale_factor = save_lh_factor_back[part];
        }

        // swap back to recover the tree
        node1->updateNeighbor(node1_it, node1_nei);
        node1_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node2_nei);
        node2_nei->node->updateNeighbor(node1, node2);
        nniid++;

	}

	delete [] save_lh_factor_back;
	delete [] save_lh_factor;
	return myMove;
}

void PhyloSuperTree::doNNI(NNIMove &move, bool clearLH) {
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	SuperNeighbor *node1_nei = (SuperNeighbor*)*move.node1Nei_it;
	SuperNeighbor *node2_nei = (SuperNeighbor*)*move.node2Nei_it;
	int part = 0;
	iterator it;
	PhyloTree::doNNI(move, clearLH);

	for (it = begin(), part = 0; it != end(); it++, part++) {
		bool is_nni = true;
		FOR_NEIGHBOR_DECLARE(move.node1, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		FOR_NEIGHBOR(move.node2, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		if (!is_nni) {
			// relink the branch if it does not correspond to NNI for partition
			linkBranch(part, nei1, nei2);
			continue;
		}

		NNIMove part_move;
		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
		part_move.node1 = (PhyloNode*)nei2_part->node;
		part_move.node2 = (PhyloNode*)nei1_part->node;
		part_move.node1Nei_it = part_move.node1->findNeighborIt(node1_nei->link_neighbors[part]->node);
		part_move.node2Nei_it = part_move.node2->findNeighborIt(node2_nei->link_neighbors[part]->node);

		(*it)->doNNI(part_move, clearLH);

	}

}

void PhyloSuperTree::changeNNIBrans(NNIMove &move) {
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	iterator it;
	int part;

	for (it = begin(), part = 0; it != end(); it++, part++) {
		bool is_nni = true;
		FOR_NEIGHBOR_DECLARE(move.node1, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		FOR_NEIGHBOR(move.node2, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		if (!is_nni) {
			continue;
		}

		NNIMove part_move;
		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
		int brid = nei1_part->id;
		part_move.node1 = (PhyloNode*)nei2_part->node;
		part_move.node2 = (PhyloNode*)nei1_part->node;
		int numlen = 1;
		if (params->nni5) numlen = 5;
		if (move.swap_id == 1) {
			for (int i = 0; i < numlen; i++)
				part_move.newLen[i] = part_info[part].nni1_brlen[brid*numlen + i];
		} else {
			for (int i = 0; i < numlen; i++)
				part_move.newLen[i] = part_info[part].nni2_brlen[brid*numlen + i];
		}

		(*it)->changeNNIBrans(part_move);

	}

}

//void PhyloSuperTree::restoreAllBrans(PhyloNode *node, PhyloNode *dad) {
//	int part = 0;
//	for (iterator it = begin(); it != end(); it++, part++) {
//		(*it)->setBranchLengths(part_info[part].cur_brlen);
//	}
//}

void PhyloSuperTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
	IQTree::reinsertLeaves(del_leaves);
	mapTrees();
}

void PhyloSuperTree::computeBranchLengths() {
	if (verbose_mode >= VB_DEBUG)
		cout << "Assigning branch lengths for full tree with weighted average..." << endl;
	int part = 0, i;
    iterator it;

	NodeVector nodes1, nodes2;
	getBranches(nodes1, nodes2);
	vector<SuperNeighbor*> neighbors1;
	vector<SuperNeighbor*> neighbors2;
	IntVector occurence;
	occurence.resize(nodes1.size(), 0);
	for (i = 0; i < nodes1.size(); i++) {
		neighbors1.push_back((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]) );
		neighbors2.push_back((SuperNeighbor*)nodes2[i]->findNeighbor(nodes1[i]) );
		neighbors1.back()->length = 0.0;
	}
	for (it = begin(), part = 0; it != end(); it++, part++) {
		IntVector brfreq;
		brfreq.resize((*it)->branchNum, 0);
		for (i = 0; i < nodes1.size(); i++) {
			PhyloNeighbor *nei1 = neighbors1[i]->link_neighbors[part];
			if (!nei1) continue;
			brfreq[nei1->id]++;
		}
		for (i = 0; i < nodes1.size(); i++) {
			PhyloNeighbor *nei1 = neighbors1[i]->link_neighbors[part];
			if (!nei1) continue;
            if ((*it)->aln->seq_type == SEQ_CODON && rescale_codon_brlen) {
                // rescale branch length by 3
                neighbors1[i]->length += (nei1->length) * (*it)->aln->getNSite() / brfreq[nei1->id];
                occurence[i] += (*it)->aln->getNSite()*3;
            } else {
                neighbors1[i]->length += (nei1->length) * (*it)->aln->getNSite() / brfreq[nei1->id];
                occurence[i] += (*it)->aln->getNSite();
            }
			//cout << neighbors1[i]->id << "  " << nodes1[i]->id << nodes1[i]->name <<"," << nodes2[i]->id << nodes2[i]->name <<": " << (nei1->length) / brfreq[nei1->id] << endl;
		}
		//cout << endl;
	}
	for (i = 0; i < nodes1.size(); i++) {
		if (occurence[i])
			neighbors1[i]->length /= occurence[i];
		neighbors2[i]->length = neighbors1[i]->length;
	}
}

string PhyloSuperTree::getModelName() {
	return (string)"Partition model";
}

PhyloTree *PhyloSuperTree::extractSubtree(set<int> &ids) {
	string union_taxa;
	for (auto it = ids.begin(); it != ids.end(); it++) {
		int id = *it;
		if (id < 0 || id >= size())
			outError("Internal error ", __func__);
		string taxa_set;
        Pattern taxa_pat = aln->getPattern(id);
        taxa_set.insert(taxa_set.begin(), taxa_pat.begin(), taxa_pat.end());
		if (it == ids.begin()) union_taxa = taxa_set; else {
			for (int j = 0; j < union_taxa.length(); j++)
				if (taxa_set[j] == 1) union_taxa[j] = 1;
		}
	}
	PhyloTree *tree = new PhyloTree;
	tree->copyTree(this, union_taxa);
	return tree;
}

uint64_t PhyloSuperTree::getMemoryRequired(size_t ncategory, bool full_mem) {
//	uint64_t mem_size = PhyloTree::getMemoryRequired(ncategory);
	// supertree does not need any memory for likelihood vectors!
	uint64_t mem_size = 0;
	for (iterator it = begin(); it != end(); it++)
		mem_size += (*it)->getMemoryRequired(ncategory, full_mem);
	return mem_size;
}

// get memory requirement for ModelFinder
uint64_t PhyloSuperTree::getMemoryRequiredThreaded(size_t ncategory, bool full_mem) {
    // only get the largest k partitions (k=#threads)
    int threads = (params->num_threads != 0) ? params->num_threads : params->num_threads_max;
    threads = min(threads, countPhysicalCPUCores());
    threads = min(threads, (int)size());
    
    // sort partition by computational cost for OpenMP effciency
    uint64_t *part_mem = new uint64_t[size()];
    int i;
    for (i = 0; i < size(); i++) {
        part_mem[i] = at(i)->getMemoryRequired(ncategory, full_mem);
    }
    
    // sort partition memory in increasing order
    quicksort<uint64_t, int>(part_mem, 0, size()-1);
    
    uint64_t mem = 0;
    for (i = size()-threads; i < size(); i++)
        mem += part_mem[i];
    
    delete [] part_mem;
    
    return mem;
}

int PhyloSuperTree::countEmptyBranches(PhyloNode *node, PhyloNode *dad) {
	int count = 0;
    if (!node)
        node = (PhyloNode*)root;

    FOR_NEIGHBOR_IT(node, dad, it) {
    	SuperNeighbor *nei = (SuperNeighbor*)(*it);
    	bool isempty = true;
    	for (PhyloNeighborVec::iterator nit = nei->link_neighbors.begin(); nit != nei->link_neighbors.end(); nit++)
    		if ((*nit)) {
    			isempty = false;
    			break;
    		}
    	if (isempty) count++;
    	count += countEmptyBranches((PhyloNode*)(*it)->node, node);
    }
    return count;
}

/** remove identical sequences from the tree */
void PhyloSuperTree::removeIdenticalSeqs(Params &params) {
	IQTree::removeIdenticalSeqs(params);
	if (removed_seqs.empty()) return;
	// now synchronize aln
	int part = 0;
    SuperAlignment *saln = (SuperAlignment*)aln;
	for (iterator it = begin(); it != end(); it++, part++) {
		if (verbose_mode >= VB_MED) {
			cout << "Partition " << saln->partitions[part]->name << " " << saln->partitions[part]->getNSeq() <<
					" sequences from " << (*it)->aln->getNSeq() << " extracted" << endl;
		}
		(*it)->aln = saln->partitions[part];
	}
	if (verbose_mode >= VB_MED) {
		cout << "Reduced alignment has " << aln->getNSeq() << " sequences with " << getAlnNSite() << " sites and "
				<< getAlnNPattern() << " patterns" << endl;
	}

}

/** reinsert identical sequences into the tree and reset original alignment */
void PhyloSuperTree::reinsertIdenticalSeqs(Alignment *orig_aln) {
	if (removed_seqs.empty()) return;
	IQTree::reinsertIdenticalSeqs(orig_aln);

	// now synchronize aln
	int part = 0;
    for (iterator it = begin(); it != end(); it++, part++) {
//        (*it)->setAlignment(((SuperAlignment*)aln)->partitions[part]);
		(*it)->aln = ((SuperAlignment*)aln)->partitions[part];
    }
	mapTrees();


}

int PhyloSuperTree::fixNegativeBranch(bool force, Node *node, Node *dad) {
	mapTrees();
	int fixed = 0;
	for (iterator it = begin(); it != end(); it++) {
		(*it)->initializeAllPartialPars();
		(*it)->clearAllPartialLH();
		fixed += (*it)->fixNegativeBranch(force);
		(*it)->clearAllPartialLH();
	}
	computeBranchLengths();
	return fixed;
}

/****************************************************************************
        ancestral sequence reconstruction
 ****************************************************************************/

void PhyloSuperTree::initMarginalAncestralState(ostream &out, bool &orig_kernel_nonrev, double* &ptn_ancestral_prob, int* &ptn_ancestral_seq) {
    orig_kernel_nonrev = params->kernel_nonrev;
    if (!orig_kernel_nonrev) {
        // switch to nonrev kernel to compute _pattern_lh_cat_state
        params->kernel_nonrev = true;
        setLikelihoodKernel(sse);
        clearAllPartialLH();
    }

    size_t total_size = 0, total_ptn = 0;

    bool mixed_data = false;

    for (auto it = begin(); it != end(); it++) {
        size_t nptn = (*it)->aln->size();
        size_t nstates = (*it)->model->num_states;
        (*it)->_pattern_lh_cat_state = (*it)->newPartialLh();
        total_size += nptn*nstates;
        total_ptn += nptn;
        if (nstates != front()->model->num_states)
            mixed_data = true;
    }

    ptn_ancestral_prob = aligned_alloc<double>(total_size);
    ptn_ancestral_seq = aligned_alloc<int>(total_ptn);
}

/**
    compute ancestral sequence probability for an internal node by marginal reconstruction
    (Yang, Kumar and Nei 1995)
    @param dad_branch branch leading to an internal node where to obtain ancestral sequence
    @param dad dad of the target internal node
    @param[out] ptn_ancestral_prob pattern ancestral probability vector of dad_branch->node
*/
void PhyloSuperTree::computeMarginalAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad,
    double *ptn_ancestral_prob, int *ptn_ancestral_seq) {

    SuperNeighbor *snei = (SuperNeighbor*)dad_branch;
    SuperNeighbor *snei_back = (SuperNeighbor*)dad_branch->node->findNeighbor(dad);
    int part = 0;
    for (auto it = begin(); it != end(); it++, part++) {
        size_t nptn = (*it)->getAlnNPattern();
        size_t nstates = (*it)->model->num_states;
        if (snei->link_neighbors[part]) {
            (*it)->computeMarginalAncestralState(snei->link_neighbors[part], (PhyloNode*)snei_back->link_neighbors[part]->node,
                ptn_ancestral_prob, ptn_ancestral_seq);
        } else {
            // branch does not exist in partition tree
            double eqprob = 1.0/nstates;
            for (size_t ptn = 0; ptn < nptn; ptn++) {
                for (size_t i = 0; i < nstates; i++)
                    ptn_ancestral_prob[ptn*nstates+i] = eqprob;
                ptn_ancestral_seq[ptn] = (*it)->aln->STATE_UNKNOWN;
            }
        }
        ptn_ancestral_prob += nptn*nstates;
        ptn_ancestral_seq += nptn;
    }
}

void PhyloSuperTree::computeSubtreeAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad,
    double *ptn_ancestral_prob, int *ptn_ancestral_seq)
{
    SuperNeighbor *snei = (SuperNeighbor*)dad_branch;
    SuperNeighbor *snei_back = (SuperNeighbor*)dad_branch->node->findNeighbor(dad);
    int part = 0;
    for (auto it = begin(); it != end(); it++, part++) {
        size_t nptn = (*it)->getAlnNPattern();
        size_t nstates = (*it)->model->num_states;
        if (snei->link_neighbors[part]) {
            (*it)->computeSubtreeAncestralState(snei->link_neighbors[part], (PhyloNode*)snei_back->link_neighbors[part]->node,
                ptn_ancestral_prob, ptn_ancestral_seq);
        } else {
            // branch does not exist in partition tree
            double eqprob = 1.0/nstates;
            for (size_t ptn = 0; ptn < nptn; ptn++) {
                for (size_t i = 0; i < nstates; i++)
                    ptn_ancestral_prob[ptn*nstates+i] = eqprob;
                ptn_ancestral_seq[ptn] = (*it)->aln->STATE_UNKNOWN;
            }
        }
        ptn_ancestral_prob += nptn*nstates;
        ptn_ancestral_seq += nptn;
    }

}

void PhyloSuperTree::writeMarginalAncestralState(ostream &out, PhyloNode *node,
    double *ptn_ancestral_prob, int *ptn_ancestral_seq) {
    int part = 1;
    for (auto it = begin(); it != end(); ++it, ++part) {
        size_t nsites  = (*it)->getAlnNSite();
        int    nstates = (*it)->model->num_states;
        for (size_t site = 0; site < nsites; ++site) {
            int ptn = (*it)->aln->getPatternID(site);
            out << node->name << "\t" << part << "\t" << site+1 << "\t";
            out << (*it)->aln->convertStateBackStr(ptn_ancestral_seq[ptn]);
            double *state_prob = ptn_ancestral_prob + ptn*nstates;
            for (int j = 0; j < nstates; ++j) {
                out << "\t" << state_prob[j];
            }
            out << endl;
        }
        size_t nptn = (*it)->getAlnNPattern();
        ptn_ancestral_prob += nptn*nstates;
        ptn_ancestral_seq += nptn;
    }
}

/**
    end computing ancestral sequence probability for an internal node by marginal reconstruction
*/
void PhyloSuperTree::endMarginalAncestralState(bool orig_kernel_nonrev,
    double* &ptn_ancestral_prob, int* &ptn_ancestral_seq) {
    if (!orig_kernel_nonrev) {
        // switch back to REV kernel
        params->kernel_nonrev = orig_kernel_nonrev;
        setLikelihoodKernel(sse);
        clearAllPartialLH();
    }
    aligned_free(ptn_ancestral_seq);
    aligned_free(ptn_ancestral_prob);

    for (auto it = rbegin(); it != rend(); ++it) {
        aligned_free((*it)->_pattern_lh_cat_state);
        (*it)->_pattern_lh_cat_state = NULL;
    }
}

void PhyloSuperTree::writeSiteLh(ostream &out, SiteLoglType wsl, int partid) {
    int part = 1;
    for (auto it = begin(); it != end(); ++it, ++part)
        (*it)->writeSiteLh(out, wsl, part);
}

void PhyloSuperTree::writeSiteRates(ostream &out, bool bayes, int partid) {

    int part = 1;
    for (iterator it = begin(); it != end(); ++it, ++part) {
        (*it)->writeSiteRates(out, bayes, part);
    }
}

void PhyloSuperTree::writeBranch(ostream &out, Node* node1, Node* node2) {
    SuperNeighbor *nei1 = (SuperNeighbor*)node1->findNeighbor(node2);
    StrVector taxnames;
    if (getNumTaxa(node1, node2) <= leafNum / 2)
        getTaxaName(taxnames, node1, node2);
    else
        getTaxaName(taxnames, node2, node1);
    out << nei1->id+1 << ",";
    bool first = true;
    for (int i = 0; i < taxnames.size(); i++)
        if (!taxnames[i].empty()) {
            if (!first) out << " ";
            out << taxnames[i];
            first = false;
        }
    out << "," << nei1->length;
    for (int part = 0; part != size(); part++) {
        bool present = true;
        FOR_NEIGHBOR_DECLARE(node1, NULL, it) {
            SuperNeighbor *nei = (SuperNeighbor*)(*it);
            if (!nei->link_neighbors[part])
                present = false;
        }
        FOR_NEIGHBOR(node2, NULL, it) {
            SuperNeighbor *nei = (SuperNeighbor*)(*it);
            if (!nei->link_neighbors[part])
                present = false;
        }
        out << ",";
        if (present)
            out << nei1->link_neighbors[part]->length;
    }
}

void PhyloSuperTree::writeBranches(ostream &out) {
    NodeVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int i;
    out << "ID,Taxa,Len";
    for (i = 0; i < size(); i++)
        out << "," << at(i)->aln->name;
    out << endl;
    for (i = 0; i < nodes1.size(); i++) {
        writeBranch(out, nodes1[i], nodes2[i]);
        out << endl;
    }
}

void PhyloSuperTree::printBestPartitionParams(const char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        out << "#nexus" << endl
        << "begin sets;" << endl;
        int part;
        SuperAlignment *saln = (SuperAlignment*)aln;
        for (part = 0; part < size(); part++) {
            string name = saln->partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            out << "  charset " << name << " = ";
            if (!saln->partitions[part]->aln_file.empty()) out << saln->partitions[part]->aln_file << ": ";
            /*if (saln->partitions[part]->seq_type == SEQ_CODON)
                out << "CODON, ";*/
            out << saln->partitions[part]->sequence_type << ", ";
            string pos = saln->partitions[part]->position_spec;
            replace(pos.begin(), pos.end(), ',' , ' ');
            out << pos << ";" << endl;
        }
        out << "  charpartition mymodels =" << endl;
        for (part = 0; part < size(); part++) {
            string name = saln->partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            if (part > 0) out << "," << endl;
            out << "    " << at(part)->getModelNameParams(true) << ": " << name << "{" << at(part)->treeLength() << "}";
        }
        out << ";" << endl;
        out << "end;" << endl;
        out.close();
        cout << "Partition information was printed to " << filename << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}
