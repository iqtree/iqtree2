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
#include "msetsblock.h"
#include "myreader.h"

PhyloSuperTree::PhyloSuperTree()
 : IQTree()
{
	totalNNIs = evalNNIs = 0;
}


PhyloSuperTree::PhyloSuperTree(SuperAlignment *alignment, PhyloSuperTree *super_tree) :  IQTree(alignment) {
	totalNNIs = evalNNIs = 0;
	part_info = super_tree->part_info;
	for (vector<Alignment*>::iterator it = alignment->partitions.begin(); it != alignment->partitions.end(); it++) {
		PhyloTree *tree = new PhyloTree((*it));
		push_back(tree);
	}
	aln = alignment;
}

void PhyloSuperTree::readPartition(Params &params) {
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(params.partition_file);
		in.exceptions(ios::badbit);
		Params origin_params = params;
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

			//info.mem_ptnlh = NULL;
			info.nniMoves[0].ptnlh = NULL;
			info.nniMoves[1].ptnlh = NULL;
			info.cur_ptnlh = NULL;
			part_info.push_back(info);
			Alignment *part_aln = new Alignment((char*)info.aln_file.c_str(), (char*)info.sequence_type.c_str(), params.intype);
			if (!info.position_spec.empty()) {
				Alignment *new_aln = new Alignment();
				new_aln->extractSites(part_aln, info.position_spec.c_str());
				delete part_aln;
				part_aln = new_aln;
			}
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


}

void PhyloSuperTree::readPartitionNexus(Params &params) {
	Params origin_params = params;
	MSetsBlock *sets_block = new MSetsBlock();
    MyReader nexus(params.partition_file);
    nexus.Add(sets_block);
    MyToken token(nexus.inf);
    nexus.Execute(token);

    Alignment *input_aln = NULL;
    if (params.aln_file) {
    	input_aln = new Alignment(params.aln_file, params.sequence_type, params.intype);
    }

    bool empty_partition = true;
    vector<CharSet*>::iterator it;
    for (it = sets_block->charsets.begin(); it != sets_block->charsets.end(); it++)
    	if ((*it)->model_name != "") {
    		empty_partition = false;
    		break;
    	}
    if (empty_partition) {
    	cout << "NOTE: No CharPartition defined, use all CharSets" << endl;
    }

    for (it = sets_block->charsets.begin(); it != sets_block->charsets.end(); it++)
    	if (empty_partition || (*it)->model_name != "") {
			PartitionInfo info;
			info.name = (*it)->name;
			info.model_name = (*it)->model_name;
			info.aln_file = (*it)->aln_file;
			//if (info.aln_file == "" && params.aln_file) info.aln_file = params.aln_file;
			if (info.aln_file == "" && !params.aln_file)
				outError("No input data for partition ", info.name);
			info.sequence_type = (*it)->sequence_type;
			if (info.sequence_type=="" && params.sequence_type) info.sequence_type = params.sequence_type;
			info.position_spec = (*it)->position_spec;
			cout << endl << "Reading partition " << info.name << " (model=" << info.model_name << ", aln=" <<
				info.aln_file << ", seq=" << info.sequence_type << ", pos=" << info.position_spec << ") ..." << endl;
			//info.mem_ptnlh = NULL;
			info.nniMoves[0].ptnlh = NULL;
			info.nniMoves[1].ptnlh = NULL;
			info.cur_ptnlh = NULL;
			part_info.push_back(info);
			Alignment *part_aln;
			if (info.aln_file != "") {
				part_aln = new Alignment((char*)info.aln_file.c_str(), (char*)info.sequence_type.c_str(), params.intype);
			} else {
				part_aln = input_aln;
			}
			if (!info.position_spec.empty() && info.position_spec != "*") {
				Alignment *new_aln = new Alignment();
				new_aln->extractSites(part_aln, info.position_spec.c_str());
				if (part_aln != input_aln) delete part_aln;
				part_aln = new_aln;
			}
			Alignment *new_aln;
			if (params.remove_empty_seq)
				new_aln = part_aln->removeGappySeq();
			else
				new_aln = part_aln;
			if (part_aln != new_aln && part_aln != input_aln) delete part_aln;
			PhyloTree *tree = new PhyloTree(new_aln);
			push_back(tree);
			params = origin_params;
			cout << new_aln->getNSeq() << " sequences and " << new_aln->getNSite() << " sites extracted" << endl;
    	}

    if (input_aln)
    	delete input_aln;
}

void PhyloSuperTree::printPartition(const char *filename) {
	   try {
	        ofstream out;
	        out.exceptions(ios::failbit | ios::badbit);
            out.open(filename);
            out << "#nexus" << endl << "[ partition information for alignment written in .conaln file ]" << endl
            	<< "begin sets;" << endl;
            int part; int start_site;
            for (part = 0, start_site = 1; part < part_info.size(); part++) {
            	string name = part_info[part].name;
            	replace(name.begin(), name.end(), '+', '_');
            	int end_site = start_site + at(part)->getAlnNSite();
            	out << "  charset " << name << " = " << start_site << "-" << end_site-1 << ";" << endl;
            	start_site = end_site;
            }
            out << "  charpartition mymodels =" << endl;
            for (part = 0; part < part_info.size(); part++) {
            	string name = part_info[part].name;
            	replace(name.begin(), name.end(), '+', '_');
            	if (part > 0) out << "," << endl;
            	out << "    " << at(part)->getModelNameParams() << ":" << name;
            }
        	out << ";" << endl;
            out << "end;" << endl;
	        out.close();
	        cout << "Partition information was printed to " << filename << endl;
	    } catch (ios::failure) {
	        outError(ERR_WRITE_OUTPUT, filename);
	    }

}

PhyloSuperTree::PhyloSuperTree(Params &params) :  IQTree() {
	totalNNIs = evalNNIs = 0;

	cout << "Reading partition model file " << params.partition_file << " ..." << endl;
	if (detectInputFile(params.partition_file) == IN_NEXUS)
		readPartitionNexus(params);
	else
		readPartition(params);
	if (part_info.empty())
		outError("No partition found");
	aln = new SuperAlignment(this);
	string str = params.out_prefix;
	//str += ".part";
	//aln->printPhylip(str.c_str());
	str = params.out_prefix;
	str += ".conaln";
	((SuperAlignment*)aln)->printCombinedAlignment(str.c_str());
	cout << "Degree of missing data: " << ((SuperAlignment*)aln)->computeMissingData() << endl;
	cout << endl;

}

void PhyloSuperTree::setParams(Params &params) {
	IQTree::setParams(params);
	for (iterator it = begin(); it != end(); it++) {
		(*it)->params = &params;
		(*it)->sse = params.SSE;
		(*it)->optimize_by_newton = params.optimize_by_newton;
	}

}

Node* PhyloSuperTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

Node* PhyloSuperTree::newNode(int node_id, int node_name) {
    return (Node*) (new SuperNode(node_id, node_name));
}

int PhyloSuperTree::getAlnNPattern() {
	int num = 0;
	for (iterator it = begin(); it != end(); it++)
		num += (*it)->getAlnNPattern();
	return num;
}

int PhyloSuperTree::getAlnNSite() {
	int num = 0;
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
		//cout<<"nei"<<(*it)->node->name<<","<<(*it)->node->id;
		if (((SuperNeighbor*)*it)->link_neighbors[part]) {
			//cout<<" has link_nei "<<((SuperNeighbor*)*it)->link_neighbors[part]->node->id<<endl;
			part_vec.push_back(((SuperNeighbor*)*it)->link_neighbors[part]);
			child_part_vec.push_back(((SuperNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]);
			//cout<<"child"<<((SuperNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]->node->name<<((SuperNeighbor*)(*it)->node->findNeighbor(node))->link_neighbors[part]->node->id<<endl;
			/*if (child_part_vec.size() > 1 && child_part_vec.back()->id == child_part_vec.front()->id)
				cout<<"HERE" << endl;*/
			assert(child_part_vec.back()->node == child_part_vec.front()->node || child_part_vec.back()->id == child_part_vec.front()->id);
		}
		 /* else{
			cout<<" no link_nei"<<endl;
		}*/
	}

	if (part_vec.empty())
		{
		/*cout<<"part_vec is EMPTY"<<endl;
		if (child_part_vec.empty())
			cout<<"child_vec is empty"<<endl;
		else
			cout<<"child_vec is NOT empty"<<endl;
		*/
		return;}
	if (part_vec.size() == 1) {
		nei->link_neighbors[part] = child_part_vec[0];
		dad_nei->link_neighbors[part] = part_vec[0];
		/*cout<<"part_vec.size() = 1 !!!"<<endl;
		if(child_part_vec[0]){
		//cout<<"DAD_part:"<<child_part_vec[0]->node->id<<endl;
		cout<<"NODE_part:"<<part_vec[0]->node->id<<endl;}
		else
			cout<<"NO CHILD VECTOR!!!!!"<<endl;
		*/

		return;
	}
	if (part_vec[0] == child_part_vec[1]) {
		// ping-pong, out of sub-tree
		//cout<<"out of Sub Tree!!!!"<<endl;
		assert(part_vec[1] == child_part_vec[0]);
		return;
	}
	PhyloNode *node_part = (PhyloNode*) child_part_vec[0]->node;
	PhyloNode *dad_part = NULL;
	FOR_NEIGHBOR(node_part, NULL, it) {
		bool appear = false;
		//int i = 0;
		//cout<<"NEI of "<<node_part->id<<" in the scope:"<<(*it)->node->id <<" name "<<(*it)->node->name <<endl;
		for (vector<PhyloNeighbor*>::iterator it2 = part_vec.begin(); it2 != part_vec.end(); it2++){
			//cout<<"part_vec["<<i++<<"]="<<(*it2)->node->id<<" name "<<(*it2)->node->name<<endl;
			if ((*it2) == (*it)) {
				//cout<<"appeared "<<(*it)->node->id<<endl;
				appear = true; break;
			}
		}
		if (!appear) {
			//cout<<"DAD couldn't find this one among part_vec:"<<(*it)->node->id<<endl;
			assert(!dad_part);
			dad_part = (PhyloNode*)(*it)->node;
			//cout<<"DAD:"<<(*it)->node->id<<endl;
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
		assert(node);
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
		assert(dad);
		PhyloNode *node_part = (PhyloNode*)part_taxa[node->id];
		if (node_part) {
			PhyloNode *dad_part = (PhyloNode*)node_part->neighbors[0]->node;
			assert(node_part->isLeaf());
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
	assert(root);
	int part = 0;
	if (verbose_mode >= VB_DEBUG)
		drawTree(cout,  WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
	for (iterator it = begin(); it != end(); it++, part++) {
		string taxa_set = ((SuperAlignment*)aln)->getPattern(part);
		(*it)->copyTree(this, taxa_set);
		if ((*it)->getModel())
			(*it)->initializeAllPartialLh();
		NodeVector my_taxa, part_taxa;
		(*it)->getOrderedTaxa(my_taxa);
		part_taxa.resize(leafNum, NULL);
		int i;
		for (i = 0; i < leafNum; i++) {
			int id = ((SuperAlignment*)aln)->taxa_index[i][part];
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

double PhyloSuperTree::computeLikelihood(double *pattern_lh) {
	double tree_lh = 0.0;
	int ntrees = size();
	if (pattern_lh) {
		//#ifdef _OPENMP
		//#pragma omp parallel for reduction(+: tree_lh)
		//#endif
		for (int i = 0; i < ntrees; i++) {
			tree_lh += at(i)->computeLikelihood(pattern_lh);
			pattern_lh += at(i)->getAlnNPattern();
		}
	} else {
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+: tree_lh)
		#endif
		for (int i = 0; i < ntrees; i++)
			tree_lh += at(i)->computeLikelihood();
	}
	return tree_lh;
}

void PhyloSuperTree::computePatternLikelihood(double *pattern_lh, double *cur_logl, double *ptn_lh_cat) {
	int offset = 0, offset_lh_cat = 0;
	iterator it;
	for (it = begin(); it != end(); it++) {
		if (ptn_lh_cat)
			(*it)->computePatternLikelihood(pattern_lh + offset, NULL, ptn_lh_cat + offset_lh_cat);
		else
			(*it)->computePatternLikelihood(pattern_lh + offset);
		offset += (*it)->aln->getNPattern();
		offset_lh_cat += (*it)->aln->getNPattern() * (*it)->site_rate->getNDiscreteRate();
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
            outError("Wrong PhyloSuperTree::", __func__);
		}
	}
}

double PhyloSuperTree::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
	double tree_lh = 0.0;
	int ntrees = size();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh)
	#endif
	for (int i = 0; i < ntrees; i++) {
		tree_lh += at(i)->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
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
		if (pit->nniMoves[0].ptnlh)
			delete [] pit->nniMoves[0].ptnlh;
		if (pit->cur_ptnlh)
			delete [] pit->cur_ptnlh;
	}
	part_info.clear();

	for (reverse_iterator it = rbegin(); it != rend(); it++)
		delete (*it);
	clear();
}

/*
double PhyloSuperTree::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH) {
	SuperNeighbor *nei1 = ((SuperNeighbor*)node1->findNeighbor(node2));
	SuperNeighbor *nei2 = ((SuperNeighbor*)node2->findNeighbor(node1));
	assert(nei1 && nei2);
	double tree_lh = 0.0;
	int ntrees = size();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh)
	#endif
	for (int part = 0; part < ntrees; part++) {
		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
		double score;
		if (part_info[part].cur_score == 0.0)
			part_info[part].cur_score = at(part)->computeLikelihood();
		if (nei1_part && nei2_part) {
			if (part_info[part].opt_score[nei1_part->id] == 0.0) {
				part_info[part].cur_brlen[nei1_part->id] = nei1_part->length;
				part_info[part].opt_score[nei1_part->id] = at(part)->optimizeOneBranch((PhyloNode*)nei1_part->node, (PhyloNode*)nei2_part->node, clearLH);
				part_info[part].opt_brlen[nei1_part->id] = nei1_part->length;
			}
			score = part_info[part].opt_score[nei1_part->id];
		} else
			score = part_info[part].cur_score;
		tree_lh += score;
	}
	return tree_lh;
}*/

void PhyloSuperTree::initPartitionInfo() {
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].cur_score = 0.0;
		//part_info[part].null_score.clear();
		//part_info[part].null_score.resize((*it)->branchNum, 0.0);
		//part_info[part].opt_score.clear();
		//part_info[part].opt_score.resize((*it)->branchNum, 0.0);
		//part_info[part].nni1_score.clear();
		//part_info[part].nni1_score.resize((*it)->branchNum, 0.0);
		//part_info[part].nni2_score.clear();
		//part_info[part].nni2_score.resize((*it)->branchNum, 0.0);

		part_info[part].cur_brlen.resize((*it)->branchNum, 0.0);
		//part_info[part].opt_brlen.resize((*it)->branchNum, 0.0);
		if (params->nni5) {
			part_info[part].nni1_brlen.resize((*it)->branchNum * 5, 0.0);
			part_info[part].nni2_brlen.resize((*it)->branchNum * 5, 0.0);
		} else {
			part_info[part].nni1_brlen.resize((*it)->branchNum, 0.0);
			part_info[part].nni2_brlen.resize((*it)->branchNum, 0.0);
		}

		(*it)->getBranchLengths(part_info[part].cur_brlen);

		if (save_all_trees == 2) {
			// initialize ptnlh for ultrafast bootstrap
			int nptn = (*it)->getAlnNPattern();
			if (!part_info[part].cur_ptnlh)
				part_info[part].cur_ptnlh = new double[nptn];
			if (!part_info[part].nniMoves[0].ptnlh)
				part_info[part].nniMoves[0].ptnlh = new double [nptn];
			if (!part_info[part].nniMoves[1].ptnlh)
				part_info[part].nniMoves[1].ptnlh = new double [nptn];
		}
		//part_info[part].nni1_ptnlh.resize((*it)->branchNum, NULL);
		//for (dit = part_info[part].nni1_ptnlh.begin(); dit != part_info[part].nni1_ptnlh.end(); dit++, offset += nptn)
		//	(*dit) = offset;
		//part_info[part].nni2_ptnlh.resize((*it)->branchNum, NULL);
		//for (dit = part_info[part].nni2_ptnlh.begin(); dit != part_info[part].nni2_ptnlh.end(); dit++, offset += nptn)
		//	(*dit) = offset;
		/*
		part_info[part].opt_ptnlh.resize((*it)->branchNum, NULL);
		for (dit = part_info[part].opt_ptnlh.begin(); dit != part_info[part].opt_ptnlh.end(); dit++, offset += nptn)
			(*dit) = offset;
		*/
	}
}

int PhyloSuperTree::getMaxPartNameLength() {
	int namelen = 0;
	for (vector<PartitionInfo>::iterator it = part_info.begin(); it != part_info.end(); it++)
		namelen = max((int)it->name.length(), namelen);
	return namelen;
}

NNIMove PhyloSuperTree::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves, bool approx_nni, bool useLS, double lh_contribution) {
    NNIMove myMove;
    //myMove.newloglh = 0;
	SuperNeighbor *nei1 = ((SuperNeighbor*)node1->findNeighbor(node2));
	SuperNeighbor *nei2 = ((SuperNeighbor*)node2->findNeighbor(node1));
	assert(nei1 && nei2);
	SuperNeighbor *node1_nei = NULL;
	SuperNeighbor *node2_nei = NULL;
	SuperNeighbor *node2_nei_other = NULL;
	FOR_NEIGHBOR_DECLARE(node1, node2, node1_it) {
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

	//double bestScore = optimizeOneBranch(node1, node2, false);

	double nni1_score = 0.0, nni2_score = 0.0;
	int ntrees = size(), part;

	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: nni1_score, nni2_score) private(part)
	#endif
	for (part = 0; part < ntrees; part++) {
		bool is_nni = true;
		totalNNIs++;
		FOR_NEIGHBOR_DECLARE(node1, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		FOR_NEIGHBOR(node2, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		if (!is_nni && params->terrace_aware) {
			if (part_info[part].cur_score == 0.0)  {
				part_info[part].cur_score = at(part)->computeLikelihood();
				if (save_all_trees == 2)
					at(part)->computePatternLikelihood(part_info[part].cur_ptnlh, &part_info[part].cur_score);
			}
			nni1_score += part_info[part].cur_score;
			nni2_score += part_info[part].cur_score;
			continue;
		}

		evalNNIs++;

		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];

		int brid = nei1_part->id;

		/*
		if (part_info[part].nni1_score[brid] == 0.0) {
			SwapNNIParam nni_param;
			// nni_param.node1/2_nei tell swapNNIBranch what to swap first
			nni_param.node1_nei = node1_nei->link_neighbors[part];
			nni_param.node2_nei = node2_nei->link_neighbors[part];
			// pass pattern likelihood address to swapNNIBranch
			nni_param.nni1_ptnlh = part_info[part].nni1_ptnlh[brid];
			nni_param.nni2_ptnlh = part_info[part].nni2_ptnlh[brid];

			at(part)->swapNNIBranch(0.0, (PhyloNode*)nei2_part->node, (PhyloNode*)nei1_part->node, &nni_param);

			part_info[part].nni1_score[brid] = nni_param.nni1_score;
			part_info[part].nni2_score[brid] = nni_param.nni2_score;
			part_info[part].nni1_brlen[brid] = nni_param.nni1_brlen;
			part_info[part].nni2_brlen[brid] = nni_param.nni2_brlen;
		}
		nni1_score += part_info[part].nni1_score[brid];
		nni2_score += part_info[part].nni2_score[brid];
		*/

		//NNIMove part_moves[2];
		//part_moves[0].node1Nei_it = NULL;
		at(part)->getBestNNIForBran((PhyloNode*)nei2_part->node, (PhyloNode*)nei1_part->node, part_info[part].nniMoves);
		// detect the corresponding NNIs and swap if necessary
		if (!((*part_info[part].nniMoves[0].node1Nei_it == node1_nei->link_neighbors[part] && *part_info[part].nniMoves[0].node2Nei_it == node2_nei->link_neighbors[part]) ||
			(*part_info[part].nniMoves[0].node1Nei_it != node1_nei->link_neighbors[part] && *part_info[part].nniMoves[0].node2Nei_it != node2_nei->link_neighbors[part])))
		{
			NNIMove tmp = part_info[part].nniMoves[0];
			part_info[part].nniMoves[0] = part_info[part].nniMoves[1];
			part_info[part].nniMoves[1] = tmp;
		}
		nni1_score += part_info[part].nniMoves[0].newloglh;
		nni2_score += part_info[part].nniMoves[1].newloglh;
		int numlen = 1;
		if (params->nni5) numlen = 5;
		for (int i = 0; i < numlen; i++) {
			part_info[part].nni1_brlen[brid*numlen + i] = part_info[part].nniMoves[0].newLen[i];
			part_info[part].nni2_brlen[brid*numlen + i] = part_info[part].nniMoves[1].newLen[i];
		}

	}
	myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
	myMove.node1 = node1;
	myMove.node2 = node2;
	if (nni1_score > nni2_score) {
		myMove.swap_id = 1;
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei->node);
		myMove.newloglh = nni1_score;
	} else  {
		myMove.swap_id = 2;
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei_other->node);
		myMove.newloglh = nni2_score;
	}

	if (save_all_trees != 2) return myMove;

	// for bootstrap now
    //now setup pattern likelihoods per partition
	double *save_lh_factor = new double [ntrees];
	double *save_lh_factor_back = new double [ntrees];
	int nnino = 0;
	FOR_NEIGHBOR(node2, node1, node2_it) {

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
				memcpy(at(part)->_pattern_lh, part_info[part].nniMoves[nnino].ptnlh, at(part)->getAlnNPattern() * sizeof(double));
    		save_lh_factor[part] = at(part)->current_it->lh_scale_factor;
    		save_lh_factor_back[part] = at(part)->current_it_back->lh_scale_factor;
    		at(part)->current_it->lh_scale_factor = 0.0;
    		at(part)->current_it_back->lh_scale_factor = 0.0;
        }
        if (nnino == 0)
        	saveCurrentTree(nni1_score);
        else
        	saveCurrentTree(nni2_score);

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
        nnino++;

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

	//linkTrees();
}

void PhyloSuperTree::applyNNIBranches(NNIMove move) {
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

		(*it)->applyNNIBranches(part_move);

	}

}

void PhyloSuperTree::linkTrees() {
	int part = 0;
	iterator it;
	for (it = begin(), part = 0; it != end(); it++, part++) {

		(*it)->initializeTree();
		NodeVector my_taxa, part_taxa;
		(*it)->getOrderedTaxa(my_taxa);
		part_taxa.resize(leafNum, NULL);
		int i;
		for (i = 0; i < leafNum; i++) {
			int id = ((SuperAlignment*)aln)->taxa_index[i][part];
			if (id >=0) part_taxa[i] = my_taxa[id];
		}
		linkTree(part, part_taxa);
	}
}

void PhyloSuperTree::restoreAllBranLen(PhyloNode *node, PhyloNode *dad) {
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		(*it)->setBranchLengths(part_info[part].cur_brlen);
	}
}

void PhyloSuperTree::reinsertLeaves(PhyloNodeVector &del_leaves) {
	IQTree::reinsertLeaves(del_leaves);
	mapTrees();
}

void PhyloSuperTree::computeBranchLengths() {
	if (verbose_mode >= VB_DEBUG)
		cout << "Assigning branch lengths for full tree with weighted average..." << endl;
	int part = 0, i;
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
	for (iterator it = begin(); it != end(); it++, part++) {
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
			neighbors1[i]->length += (nei1->length) * (*it)->aln->getNSite() / brfreq[nei1->id];
			occurence[i] += (*it)->aln->getNSite();
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

PhyloTree *PhyloSuperTree::extractSubtree(IntVector &ids) {
	string union_taxa;
	int i;
	for (i = 0; i < ids.size(); i++) {
		int id = ids[i];
		if (id < 0 || id >= size())
			outError("Internal error ", __func__);
		string taxa_set = aln->getPattern(id);
		if (i == 0) union_taxa = taxa_set; else {
			for (int j = 0; j < union_taxa.length(); j++)
				if (taxa_set[j] == 1) union_taxa[j] = 1;
		}
	}
	PhyloTree *tree = new PhyloTree;
	tree->copyTree(this, union_taxa);
	return tree;
}

uint64_t PhyloSuperTree::getMemoryRequired() {
	uint64_t mem_size = PhyloTree::getMemoryRequired();
	for (iterator it = begin(); it != end(); it++)
		mem_size += (*it)->getMemoryRequired();
	return mem_size;
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

