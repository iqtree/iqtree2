/*
 * phylosupertreeplen.cpp
 *
 *  Created on: Aug 5, 2013
 *      Author: olga
 */

#include "phylosupertreeplen.h"
#include "alignment/superalignmentpairwiseplen.h"
#include "model/partitionmodelplen.h"
#include <string.h>
#include "utils/timeutil.h"





/**********************************************************
 * class PhyloSuperTreePlen
**********************************************************/


PhyloSuperTreePlen::PhyloSuperTreePlen()
: PhyloSuperTree()
{
	memset(allNNIcases_computed, 0, 5*sizeof(int));
	fixed_rates = false;
}

/*
PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment)
: PhyloSuperTree(alignment)
{
	memset(allNNIcases_computed, 0, 5*sizeof(int));
	fixed_rates = (params->partition_type == BRLEN_FIX) ? true : false;
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].part_rate = 1.0;
		part_info[part].evalNNIs = 0.0;
        if ((*it)->aln->seq_type == SEQ_CODON && rescale_codon_brlen)
            part_info[part].part_rate = 3.0;
	}
}
*/
PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment, int partition_type)
: PhyloSuperTree(alignment)
{
    memset(allNNIcases_computed, 0, 5*sizeof(int));
//    fixed_rates = false;
    fixed_rates = (partition_type == BRLEN_FIX) ? true : false;
    int part = 0;
    bool has_tree_len = false;
    for (iterator it = begin(); it != end(); it++, part++) {
        part_info[part].part_rate = 1.0;
        if (alignment->partitions[part]->tree_len > 0.0) {
            part_info[part].part_rate = alignment->partitions[part]->tree_len;
            has_tree_len = true;
        }
        part_info[part].evalNNIs = 0.0;
        if ((*it)->aln->seq_type == SEQ_CODON && rescale_codon_brlen)
            part_info[part].part_rate *= 3.0;
    }
    
    if (has_tree_len)
        normalizePartRate();
}

PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment, PhyloSuperTree *super_tree)
: PhyloSuperTree(alignment,super_tree)
{
	memset(allNNIcases_computed, 0, 5*sizeof(int));
	fixed_rates = false;
    int part = 0;
    bool has_tree_len = false;
    for (iterator it = begin(); it != end(); it++, part++) {
        part_info[part].part_rate = 1.0;
        if (alignment->partitions[part]->tree_len > 0.0) {
            part_info[part].part_rate = alignment->partitions[part]->tree_len;
            has_tree_len = true;
        }
        part_info[part].evalNNIs = 0.0;
        if ((*it)->aln->seq_type == SEQ_CODON && rescale_codon_brlen)
            part_info[part].part_rate *= 3.0;
    }
    if (has_tree_len)
        normalizePartRate();
}

void PhyloSuperTreePlen::normalizePartRate() {
    double sum = 0.0;
    size_t nsite = 0;
    int i;
    for (i = 0; i < size(); i++) {
        sum += part_info[i].part_rate * at(i)->aln->getNSite();
        if (at(i)->aln->seq_type == SEQ_CODON && rescale_codon_brlen)
            nsite += 3*at(i)->aln->getNSite();
        else
            nsite += at(i)->aln->getNSite();
    }
    sum /= nsite;
    
    //scaleLength(sum);
    sum = 1.0/sum;
    for (i = 0; i < size(); i++)
        part_info[i].part_rate *= sum;

}

void PhyloSuperTreePlen::deleteAllPartialLh() {
	for (iterator it = begin(); it != end(); it++) {
		// reset these pointers so that they are not deleted
		(*it)->central_partial_lh = NULL;
		(*it)->central_scale_num = NULL;
//		(*it)->central_partial_pars = NULL;
		(*it)->_pattern_lh = NULL;
		(*it)->_pattern_lh_cat = NULL;
		(*it)->theta_all = NULL;
        (*it)->buffer_scale_all = NULL;
        (*it)->buffer_partial_lh = NULL;
		(*it)->ptn_freq = NULL;
		(*it)->ptn_freq_computed = false;
        (*it)->ptn_freq_pars = NULL;
		(*it)->ptn_invar = NULL;
        (*it)->nni_partial_lh = NULL;
        (*it)->nni_scale_num = NULL;
	}
    PhyloTree::deleteAllPartialLh();
}

PhyloSuperTreePlen::~PhyloSuperTreePlen()
{
	for (iterator it = begin(); it != end(); it++) {
		// reset these pointers so that they are not deleted
		(*it)->central_partial_lh = NULL;
		(*it)->central_scale_num = NULL;
//		(*it)->central_partial_pars = NULL;
		(*it)->_pattern_lh = NULL;
		(*it)->_pattern_lh_cat = NULL;
		(*it)->theta_all = NULL;
        (*it)->buffer_scale_all = NULL;
        (*it)->buffer_partial_lh = NULL;
		(*it)->ptn_freq = NULL;
		(*it)->ptn_freq_computed = false;
        (*it)->ptn_freq_pars = NULL;
		(*it)->ptn_invar = NULL;
        (*it)->nni_partial_lh = NULL;
        (*it)->nni_scale_num = NULL;
	}
}

void PhyloSuperTreePlen::saveCheckpoint() {
    // bypass PhyloSuperTree
    IQTree::saveCheckpoint();
}

void PhyloSuperTreePlen::restoreCheckpoint() {
    // bypass PhyloSuperTree
    IQTree::restoreCheckpoint();
}

void PhyloSuperTreePlen::printResultTree(string suffix) {
    IQTree::printResultTree(suffix);
}

// -------------------------------------------------------------------------------------------------------------
double PhyloSuperTreePlen::computeDist(int seq1, int seq2, double initial_dist, double &var) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0)
        initial_dist = aln->computeDist(seq1, seq2);
    if (initial_dist == MAX_GENETIC_DIST) return initial_dist;
    if (!model_factory || !site_rate) return initial_dist;

    // now optimize the distance based on the model and site rate
    SuperAlignmentPairwisePlen aln_pair(this, seq1, seq2);
    return aln_pair.optimizeDist(initial_dist, var);
}

void PhyloSuperTreePlen::mapTrees() {
	ASSERT(root);

    syncRooting();
    
    int part = 0;
    // this is important: rescale branch length of codon partitions to be compatible with other partitions.
    // since for codon models, branch lengths = # nucleotide subst per codon site!
    bool noncodon_present = false;
    iterator it;
	for (it = begin(); it != end(); it++) {
		if ((*it)->aln->seq_type != SEQ_CODON) {
			noncodon_present = true;
			break;
		}
	}
	for (it = begin(); it != end(); it++, part++) {
		string taxa_set;
        Pattern taxa_pat = ((SuperAlignment*)aln)->getPattern(part);
        taxa_set.insert(taxa_set.begin(), taxa_pat.begin(), taxa_pat.end());
		(*it)->copyTree(this, taxa_set);

		// the only difference with PhyloSuperTree::mapTrees()
		(*it)->scaleLength(part_info[part].part_rate);

//		if ((*it)->getModel())
//			(*it)->initializeAllPartialLh();
		PhyloNodeVector my_taxa, part_taxa;
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
	if (getModel())
		initializeAllPartialLh();
}

void PhyloSuperTreePlen::linkTrees() {
	mapTrees();


}


double PhyloSuperTreePlen::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
	//initPartitionInfo(); // OLGA: not needed here
	//cout<<"Optimizing all branches"<<endl;
	for(int part = 0; part < size(); part++){
		part_info[part].cur_score = 0.0;
	}

	return PhyloTree::optimizeAllBranches(my_iterations,tolerance, maxNRStep);
}

void PhyloSuperTreePlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH, int maxNRStep) {
	if (rooted && (node1 == root || node2 == root))
	{
		return; // does not optimize virtual branch from root
	}    
	SuperNeighbor *nei1 = (SuperNeighbor*)node1->findNeighbor(node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)node2->findNeighbor(node1);
	int part;

	current_it      = node1->findNeighbor(node2);
    current_it_back = node2->findNeighbor(node1);
	for (part = 0; part < size(); part++) {
		if (((SuperNeighbor*)current_it)->link_neighbors[part]) {
            at(part)->current_it = ((SuperNeighbor*)current_it)->link_neighbors[part];
            at(part)->current_it_back = ((SuperNeighbor*)current_it_back)->link_neighbors[part];
		}
	}
    
	double current_len = current_it->length;
	for (part = 0; part < size(); part++) {
		at(part)->theta_computed = false;
	}

	//this->clearAllPartialLH();
	PhyloTree::optimizeOneBranch(node1, node2, false, maxNRStep);

    if (part_order.empty()) computePartitionOrder();
	// bug fix: assign cur_score into part_info
    #ifdef _OPENMP
    #pragma omp parallel for private(part) schedule(dynamic) if(num_threads > 1)
    #endif    
    for (int partid = 0; partid < size(); partid++) {
        part = part_order_by_nptn[partid];
        if (((SuperNeighbor*)current_it)->link_neighbors[part]) {
            part_info[part].cur_score = at(part)->computeLikelihoodFromBuffer();
        }
    }

	if(clearLH && current_len != current_it->length){
		for (int part = 0; part < size(); part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			if(nei1_part){
				nei1_part->getNode()->clearReversePartialLh(nei2_part->getNode());
				nei2_part->getNode()->clearReversePartialLh(nei1_part->getNode());
			}
		}
	}

//	return tree_lh;
}

double PhyloSuperTreePlen::computeFunction(double value) {

	double tree_lh = 0.0;
	int ntrees = size();

	if (!central_partial_lh) initializeAllPartialLh();

	double lambda = value-current_it->length;
	current_it->length = value;
    current_it_back->length = value;

	SuperNeighbor *nei1 = (SuperNeighbor*)current_it_back->node->findNeighbor(current_it->node);
	SuperNeighbor *nei2 = (SuperNeighbor*)current_it->node->findNeighbor(current_it_back->node);
	ASSERT(nei1 && nei2);

    if (part_order.empty()) computePartitionOrder();
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(num_threads > 1)
    #endif    
	for (int partid = 0; partid < ntrees; partid++) {
            int part = part_order_by_nptn[partid];
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			if (nei1_part && nei2_part) {
				at(part)->current_it = nei1_part;
				at(part)->current_it_back = nei2_part;
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei2_part,nei1_part->getNode());
				tree_lh += part_info[part].cur_score;
			} else {
				if (part_info[part].cur_score == 0.0)
					part_info[part].cur_score = at(part)->computeLikelihood();
				tree_lh += part_info[part].cur_score;
			}
		}
    return -tree_lh;
}

double PhyloSuperTreePlen::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    current_it      = dad_branch;
    current_it_back = dad_branch->getNode()->findNeighbor(dad);
    return -computeFunction(dad_branch->length);
}

double PhyloSuperTreePlen::computeLikelihoodFromBuffer() {
    //return -computeFunction(current_it->length);
	double score = 0.0;
	int part, ntrees = size();
	for (part = 0; part < ntrees; part++) {
//		assert(part_info[part].cur_score != 0.0);
		score += part_info[part].cur_score;
	}
	return score;
}

void PhyloSuperTreePlen::computeFuncDerv(double value, double &df_ret, double &ddf_ret) {
//	double tree_lh = 0.0;
	double df = 0.0;
	double ddf = 0.0;

	int ntrees = size();

	if (!central_partial_lh) initializeAllPartialLh();

	double lambda = value-current_it->length;
	current_it->length = value;
    current_it_back->length = value;

	SuperNeighbor *nei1 = (SuperNeighbor*)current_it_back->node->findNeighbor(current_it->node);
	SuperNeighbor *nei2 = (SuperNeighbor*)current_it->node->findNeighbor(current_it_back->node);
	ASSERT(nei1 && nei2);

    if (part_order.empty()) computePartitionOrder();
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+: df, ddf) schedule(dynamic) if(num_threads > 1)
    #endif    
	for (int partid = 0; partid < ntrees; partid++) {
        int part = part_order_by_nptn[partid];
        double df_aux, ddf_aux;
        PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
        PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
        if (nei1_part && nei2_part) {
            at(part)->current_it = nei1_part;
            at(part)->current_it_back = nei2_part;
            
            nei1_part->length += lambda*part_info[part].part_rate;
            nei2_part->length += lambda*part_info[part].part_rate;
            if(nei1_part->length<-1e-4) {
                cout<<"lambda = "<<lambda<<endl;
                cout<<"NEGATIVE BRANCH len = "<<nei1_part->length<<endl<<" rate = "<<part_info[part].part_rate<<endl;
                ASSERT(0);
                outError("shit!!   ",__func__);
            }
            at(part)->computeLikelihoodDerv(nei2_part, nei1_part->getNode(), &df_aux, &ddf_aux);
            df += part_info[part].part_rate*df_aux;
            ddf += part_info[part].part_rate*part_info[part].part_rate*ddf_aux;
        }
        else {
            if (part_info[part].cur_score == 0.0) {
                part_info[part].cur_score = at(part)->computeLikelihood();
            }
        }
    }
    df_ret = -df;
    ddf_ret = -ddf;
}

NNIMove PhyloSuperTreePlen::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves)
{
    if ((node1->findNeighbor(node2))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        PhyloNode *tmp = node1;
        node1 = node2;
        node2 = tmp;
    }
	ASSERT(node1->degree() == 3 && node2->degree() == 3);

	double backupScore = curScore;


	SwapNNIParam nni_param;
	// nni_param.node1/2_nei tell swapNNIBranch what to swap first

	// ------------------------------------------------------------------
    int cnt;
    NNIMove localNNIMoves[2];
    if (nniMoves==nullptr) {
        nniMoves = localNNIMoves;
        //NNIMove constructor now sets node1 and ptnlh to nullptr (James B. 06-Aug-2020)
    }
    if (nniMoves[0].node1) {
    	// assuming that node1Nei_it and node2Nei_it are defined in nniMoves structure
    	for (cnt = 0; cnt < 2; cnt++) {
    		// sanity check
    		if (!node1->findNeighbor((*nniMoves[cnt].node1Nei_it)->node)) outError(__func__);
    		if (!node2->findNeighbor((*nniMoves[cnt].node2Nei_it)->node)) outError(__func__);
    	}
    } else {
        FOR_NEIGHBOR_IT(node1, node2, node1_it)
        if (((PhyloNeighbor*)*node1_it)->direction != TOWARD_ROOT)
        {
			cnt = 0;
			FOR_NEIGHBOR_IT(node2, node1, node2_it) {
				//   Initialize the 2 NNI moves
				nniMoves[cnt].node1Nei_it = node1_it; // the same neighbor of node1 for cnt = 0 and cnt = 1
				nniMoves[cnt].node2Nei_it = node2_it;
				cnt++;
			}
			break;
        }
    }

    // Initialize node1 and node2 in nniMoves
	nniMoves[0].node1 = nniMoves[1].node1 = node1;
	nniMoves[0].node2 = nniMoves[1].node2 = node2;
    nniMoves[0].newloglh = nniMoves[1].newloglh = -DBL_MAX;

    // check for compatibility with constraint
    // check for consistency with constraint tree
    for (cnt = 0; cnt < 2; cnt++) {
        if (!constraintTree.isCompatible(nniMoves[cnt])) {
            nniMoves[cnt].node1 = nniMoves[cnt].node2 = NULL;
        }
    }

	//--------------------------------------------------------------------------

    if (nniMoves[0].node1 || nniMoves[1].node1)
        this->swapNNIBranch(0.0, node1, node2, &nni_param, nniMoves);


	 // restore curScore
	 curScore = backupScore;

    NNIMove myMove;
    if (nniMoves[0].newloglh > nniMoves[1].newloglh) {
        myMove = nniMoves[0];
        myMove.swap_id = 1;
    } else {
        myMove = nniMoves[1];
        myMove.swap_id = 2;
    }
    return myMove;
}

void PhyloSuperTreePlen::doNNIs(vector<NNIMove> &compatibleNNIs, bool changeBran) {
	IQTree::doNNIs(compatibleNNIs, changeBran);
	mapBranchLen();
	//clearAllPartialLH();
}


void PhyloSuperTreePlen::getNNIType(PhyloNode *node1, PhyloNode *node2, vector<NNIType> &nni_type) {
	int epsilon_cnt, part, ntrees=size();
	nni_type.resize(ntrees, NNI_NO_EPSILON);
	for(part=0; part<ntrees;part++){
		totalNNIs++;
		nni_type[part] = NNI_NO_EPSILON;
		epsilon_cnt = 0;

		FOR_NEIGHBOR_DECLARE(node1,NULL,nit){
			if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { epsilon_cnt++; }
		}
		FOR_NEIGHBOR(node2, node1, nit) {
			if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { epsilon_cnt++; }
		}
		if(epsilon_cnt == 0){
			nni_type[part]=NNI_NO_EPSILON;
		}else if(epsilon_cnt == 1){
			nni_type[part] = NNI_ONE_EPSILON;
		}else if(epsilon_cnt == 2){
			nni_type[part]=NNI_TWO_EPSILON;
		}else if(epsilon_cnt == 3){
			nni_type[part]=NNI_THREE_EPSILON;
		}else {
			nni_type[part] = NNI_MANY_EPSILON;
		}
	}
}

void PhyloSuperTreePlen::doNNI(NNIMove &move, bool clearLH)
{
	//checkBranchLen();
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	SuperNeighbor *node1_nei = (SuperNeighbor*)*move.node1Nei_it;
	SuperNeighbor *node2_nei = (SuperNeighbor*)*move.node2Nei_it;

	int part = 0, ntrees = size();
	iterator it;
	vector<NNIMove> part_move;
	vector<NNIType> is_nni;
	part_move.resize(ntrees);
	getNNIType(move.node1, move.node2, is_nni);

	for (it = begin(), part = 0; it != end(); it++, part++) {
		if(is_nni[part] == NNI_NO_EPSILON){
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			part_move[part].node1 = nei2_part->getNode();
			part_move[part].node2 = nei1_part->getNode();
			part_move[part].node1Nei_it = part_move[part].node1->findNeighborIt(node1_nei->link_neighbors[part]->node);
			part_move[part].node2Nei_it = part_move[part].node2->findNeighborIt(node2_nei->link_neighbors[part]->node);
		}
	}
//	PhyloTree::doNNI(move,clearLH);
	PhyloTree::doNNI(move,false);
	PhyloNode *node1, *node2;

	for (it = begin(), part = 0; it != end(); it++, part++) {
		switch (is_nni[part]) {
		case NNI_NO_EPSILON:
			(*it)->doNNI(part_move[part],clearLH);
			break;
		case NNI_ONE_EPSILON:
			linkBranch(part, nei1, nei2);
			if (clearLH) {
				// clear partial likelihood vector
				node1 = nei2->link_neighbors[part]->getNode();
				node2 = nei1->link_neighbors[part]->getNode();
				nei1->link_neighbors[part]->clearPartialLh();
				nei2->link_neighbors[part]->clearPartialLh();
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_TWO_EPSILON:
			node1 = nei2->link_neighbors[part]->getNode();
			node2 = nei1->link_neighbors[part]->getNode();
			linkBranch(part, nei1, nei2);
			if(clearLH){
				// the check "&& !(PhyloNode*)nei2->link_neighbors[part]" is not needed,
				// since the branch lengths are changed during the optimization
				// and we anyway have to clearReversePartialLh
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_THREE_EPSILON:
			linkBranch(part, nei1, nei2);
			if (clearLH) {
				// clear partial likelihood vector
				node1 = nei2->link_neighbors[part]->getNode();
				node2 = nei1->link_neighbors[part]->getNode();
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_MANY_EPSILON:
			break;
		}
	}
}

double PhyloSuperTreePlen::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, SwapNNIParam *nni_param, NNIMove *nniMoves) {


	int i = 0, id = 0;
	int part, ntrees = size();
    uint64_t total_block_size = 0, total_scale_block_size = 0;
    for (int j = 0; j < ntrees; j++) {
        total_block_size += block_size[j];
        total_scale_block_size += scale_block_size[j];
    }

	/*===========================================================================================
	 * Identify NNIType for partitions
	 *===========================================================================================*/
	vector<NNIType> is_nni;
	getNNIType(node1, node2, is_nni);
	if(verbose_mode >= VB_MED){
		for (part = 0; part < ntrees; part++)
			switch (is_nni[part]) {
			case NNI_NO_EPSILON:
				allNNIcases_computed[0]++;
				break;
			case NNI_ONE_EPSILON:
				allNNIcases_computed[1]++;
				break;
			case NNI_TWO_EPSILON:
				allNNIcases_computed[2]++;
				break;
			case NNI_THREE_EPSILON:
				allNNIcases_computed[3]++;
				break;
			case NNI_MANY_EPSILON:
				allNNIcases_computed[4]++;
				break;
			}
	}
	//==================================================================================================
	// SuperTREE: saving Neighbors and allocating new ones; assign which nodes/neighbors to be swapped.
	//==================================================================================================
	double old_brlen = node1->findNeighbor(node2)->length; // length of the branch between node1 and node2 on SuperTree before NNI
	int IT_NUM = (params->nni5) ? 6 : 2;
	NeighborVec::iterator it, saved_it[6], node_nei_it[4];
	Node* neighbor_nodes[4];

	saved_it[id++] = node1->findNeighborIt(node2);
	saved_it[id++] = node2->findNeighborIt(node1);

	//if (params->nni5) {
		FOR_NEIGHBOR(node1, node2, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node1);
			node_nei_it[i++] = it;
			neighbor_nodes[i-1] = (*it)->node;
		}
		FOR_NEIGHBOR(node2, node1, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node2);
			node_nei_it[i++] = it;
			neighbor_nodes[i-1] = (*it)->node;
		}
	//}

//		cout<<"------NODE_id check-----------------------------------"<<endl;
//		for(part=0; part<ntrees; part++){
//			cout<<"PART = "<<part<<endl;
//			for(id=2; id<6; id++){
//				if(node1->isNeighbor(neighbor_nodes[id-2])){
//					if(((SuperNeighbor*)(node1->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]){
//						cout<<"node1: "<<"id = "<<id<<"; node_id = "<<
//								((SuperNeighbor*)(node1->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]->node->id<<";"<<endl;
//					} else {
//						cout<<"node1: "<<"id = "<<id<<"; no neighbor;"<<endl;
//					}
//				} else if(node2->isNeighbor(neighbor_nodes[id-2])){
//					if(((SuperNeighbor*)(node2->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]){
//						cout<<"node2: "<<"id = "<<id<<"; node_id = "<<
//								((SuperNeighbor*)(node2->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]->node->id<<";"<<endl;
//					} else {
//						cout<<"node2: "<<"id = "<<id<<"; no neighbor;"<<endl;
//					}
//				}
//			}
//			cout<<"------"<<endl;
//			for(id=2; id<6; id++){
//				if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
//					cout<<"id = "<<id<<"; node_id = "<<((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node->id<<";"<<endl;
//				}
//			}
//		}
//		cout<<"------------------------------------------------------"<<endl;

    // reorient partial_lh in case of nni1
    if (!params->nni5) {
        reorientPartialLh(node1->findNeighbor(node2), node1);
        reorientPartialLh(node2->findNeighbor(node1), node2);
    }

	/*------------------------------------------------------------------------------------
	 * Saving original neighbors:
	 * saved_nei[0] - node2 as a neighbor of node1
	 * saved_nei[1] - node1 as a neighbor of node2
	 * IF(nni5Branches)
	 * 		saved_nei[2(3)] - node1 as a neighbor of its nei1(nei2) different from node2
	 * 		saved_nei[4(5)] - node2 as a neighbor of its nei1(nei2) different from node1
	 *------------------------------------------------------------------------------------*/

	SuperNeighbor *saved_nei[6];

	// allocate new Super Neighbor pointers
	for (id = 0; id < IT_NUM; id++) {
		saved_nei[id] = (SuperNeighbor*)(*saved_it[id]);
		*saved_it[id] = saved_nei[id]->newNeighbor();
		for(part = 0; part < ntrees; part++)
			((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
	}

	// Getting NEW Neighbors: get the Neighbors again since they were saved for restoring purpose and replaced by new
	SuperNeighbor *nei1_new = (SuperNeighbor*) node1->findNeighbor(node2);
	SuperNeighbor *nei2_new = (SuperNeighbor*) node2->findNeighbor(node1);

//	/* -------------------------------------------------------------------------------------------
//	 *  NNI details: assigning nodes to be swapped on SuperTree
//	 * -------------------------------------------------------------------------------------------*/
//
//	// node1_nei - one of the node1 neighbors, which is not node2
//	NeighborVec::iterator node1_it = node1->findNeighborIt(nni_param->node1_nei->node);
//	Neighbor *node1_nei = *node1_it;
//
//	// *node2_its[0] - one of the node2 neighbors, which is not node1
//	// *node2_its[1] - second neighbor of node2,   which is not node1
//	vector<NeighborVec::iterator> node2_its;
//	node2_its.push_back(node2->findNeighborIt(nni_param->node2_nei->node));
//
//	FOR_NEIGHBOR_DECLARE(node2, node1, node2_it){
//		FOR_NEIGHBOR_DECLARE(node2,(*node2_it)->node,node2_it2)
//			node2_its.push_back(node2_it2);
//	}
//	assert(node2_its.size() == 2);

	/* =================================================================================================
	 * SubTREEs: saving Neighbors and allocating new ones.
	 * =================================================================================================*/

	/*------------------------------------------------------------------------------------
	 * Variables to be used for saving/restoring purposes on SubTrees
	 *------------------------------------------------------------------------------------*/

	vector<PhyloNeighbor*> sub_saved_nei1,sub_saved_nei2,sub_saved_nei;
	vector<NeighborVec::iterator> sub_saved_it;

	// Saving linked neighbor of node1->findNei(node2)
	sub_saved_nei1.resize(ntrees);
	// Saving linked neighbor of node2->findNei(node1)
	sub_saved_nei2.resize(ntrees);

	sub_saved_nei.resize(6*ntrees);
	sub_saved_it.resize(6*ntrees);

	/*---------------------------------------------------------
	 * For Restoring: saving branch lengths on SubTrees
	 *---------------------------------------------------------*/
	/* NO_EPS:  one/five branches need to be restored in nni1/nni5 cases respectively,
	 * 			but these branches are saved in saved_nei->link_neighbors, we won't store them again
	 * ONE_EPS: three branches need to be restored (stick to ids: 0,...,5)
	 * TWO_EPS: the image of central branch needs to be restored (the id for restoring [6*part+0])
	 * THREE_EPS: one branch needs to be restored: which the central is relinked to after NNI (the id for restoring [6*part+0])
	 * MANY_EPS: nothing to be restored
	 */
	double *sub_saved_branch = new double[6*ntrees];

	/* ---------------------------------------------------------
	 * For Restoring: saving current likelihoods for SubTree
	 * ---------------------------------------------------------*/
	double *saved_cur_score = new double[part_info.size()];
	for (i = 0; i < part_info.size(); i++)
		saved_cur_score[i] = part_info[i].cur_score;

	/* -------------------------------------------------------------------------------------------------------------------
	 * Allocate new PhyloNeighbors:
	 * NO_EPS:  2 or 6 for nni1 and nni5 respectively; update link_neighbors for corresponding SuperNeighbors.
	 * ONE_EPS: 1 or 3 for nni1 and nni5 respectively; update link_neighbors for corresponding SuperNeighbors LATER
	 * 			(since it depends on particular NNI).
	 * -------------------------------------------------------------------------------------------------------------------*/

	// Auxiliary variables: we allocate new PhyloNeighbor for [node_link->findNeighbor(nei_link)]
	Node *node_link, *nei_link;
	SuperNeighbor *nei;

	// For ONE_epsilon case: saves "id" of the neighbors that have an empty image
	int id_eps[part];
    uint64_t lh_addr = 0, scale_addr = 0;
	for(int partid = 0; partid < ntrees; partid++){
        part = part_order[partid];
		if(is_nni[part]==NNI_NO_EPSILON){
			//evalNNIs++;
			//part_info[part].evalNNIs++;

            int mem_id = 0;

			// one branch optimization ------------------------------------------------------------------
			for(id = 0; id < 2; id++){
				/*
					for id=0, nei_link  = node1->find(node2)->link->node = node2_link
				 	for id=0, node_link = node2->find(node1)->link->node = node1_link
					for id=0, saving iterator of neighbor of node1_link, that is node2_link;
						 	  then on this place we'll create a new PhyloNei
				*/

				nei_link  = saved_nei[id]->link_neighbors[part]->node;
				node_link = saved_nei[1-id]->link_neighbors[part]->node;
				sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);

				// Create a new PhyloNeighbor, with new partial lhs, scale number and set the branch id as before
				*sub_saved_it[part*6 + id] = saved_nei[id]->link_neighbors[part]->newNeighbor();
                if (saved_nei[id]->link_neighbors[part]->partial_lh) {
                    ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = nni_partial_lh + (mem_id*total_block_size + lh_addr);
                    ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = nni_scale_num + (mem_id*total_scale_block_size + scale_addr);
                    mem_id++;
                }

				// update link_neighbor[part]: for New SuperNeighbor we set the corresponding new PhyloNeighbor on partition part
				((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
			}

			// optimization on 5 branches ------------------------------------------------------------------
			if(params->nni5){
				for(id = 2; id < 6; id ++){
					nei_link = saved_nei[id]->link_neighbors[part]->node;
					node_link = ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node;
					sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
					*sub_saved_it[part*6 + id] = saved_nei[id]->link_neighbors[part]->newNeighbor();
                    if (saved_nei[id]->link_neighbors[part]->partial_lh) {
                        ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = nni_partial_lh + (mem_id*total_block_size + lh_addr);
                        ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = nni_scale_num + (mem_id*total_scale_block_size + scale_addr);
                        mem_id++;
                    }

					// update link_neighbor[part]
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
				}
                ASSERT(mem_id == 2);
			}


		} else if(is_nni[part]==NNI_ONE_EPSILON){

            int mem_id = 0;
			// Make sure to update all the necessary link_neighbors and take care of branch lengths
			// (increase/decrease by central branch where necessary).

			nei1_new->link_neighbors[part] = saved_nei[0]->link_neighbors[part];
			nei2_new->link_neighbors[part] = saved_nei[1]->link_neighbors[part];

			// Change the length of branch, (node1,node2) WAS linked to (-=)
			nei1_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
			nei2_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
            // Rooted branch length can go below 0 due to numerical precision
            ASSERT(nei1_new->link_neighbors[part]->length >= -1e-6);
            if (nei1_new->link_neighbors[part]->length < 0) {
                nei1_new->link_neighbors[part]->length = 0;
                nei2_new->link_neighbors[part]->length = 0;
            }

			// Allocate three new PhyloNeighbors.
			// For nni1 only one of it will be actually used and which one depends on the NNI.

			// We have this if condition, since saved_nei will be newly allocated neis in nni5 case,
			// while saved_it are the actual neighbors and we don't want to mess them up
			for(id = 2; id < 6; id++){
				if(params->nni5){
					nei = saved_nei[id];
				} else {
					nei = (SuperNeighbor*)(*saved_it[id]);
				}
				if(nei->link_neighbors[part]){
					// nei_link is either node1 or node2 on SubTrees
					nei_link = nei->link_neighbors[part]->node;
					// node_link are nodes neighbors of node1 and node2 on SubTrees
					node_link = ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node;
					sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);

					// Saving branch lengths
					sub_saved_branch[6*part + id] = nei->link_neighbors[part]->length;

					*sub_saved_it[part*6 + id] = nei->link_neighbors[part]->newNeighbor();
                    if (nei->link_neighbors[part]->partial_lh) {
                        ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = nni_partial_lh + (mem_id*total_block_size + lh_addr);
                        ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = nni_scale_num + (mem_id*total_scale_block_size + scale_addr);
                        mem_id++;
                    }

					// If nni5 we update the link neighbors already here, otherwise
					// they will be updated for each NNI within the loop.
					if(params->nni5){
						((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
					}
					//cout<<"saved_it["<<id<<"]; neighbor->node->id = "<<(*sub_saved_it[part*6 + id])->node->id<<endl;
					//cout<<"saved_it["<<id<<"];           node->id = "<<node_link->id<<endl;
				} else {
					id_eps[part] = id;
				}
			}
            ASSERT(mem_id == 1);
		}else if(is_nni[part]==NNI_THREE_EPSILON && params->nni5){
			// you fill out link neighbors vector for newly allocated SuperNeighbors
			for(id = 2; id < 6; id++){
				if(saved_nei[id]->link_neighbors[part]){
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = saved_nei[id]->link_neighbors[part];
				}
			}
		}else if(is_nni[part]==NNI_TWO_EPSILON && params->nni5){
			// you fill out link neighbors vector for newly allocated SuperNeighbors
			for(id = 2; id < 6; id++){
				if(saved_nei[id]->link_neighbors[part]){
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = saved_nei[id]->link_neighbors[part];
				}
			}
		}
        lh_addr += block_size[part];
        scale_addr += scale_block_size[part];
	}

	/* -------------------------------------------------------------------
	 * Variables to store the information about which nodes/neighbors
	 * to be swapped on SubTrees for the corresponding NNI on SuperTree
	 *
	 * node1 -> node1_link[part]
	 * node2 -> node2_link[part]
	 * node1_nei -> node1_link_nei[part]
	 * node2_nei -> node2_link_nei[part]
	 * -------------------------------------------------------------------*/
	vector<PhyloNode*> node1_link,node2_link;
	vector<PhyloNeighbor*> node1_link_nei,node2_link_nei;
	vector<NeighborVec::iterator> node1_link_it, node2_link_it;

	// Nodes which correspond to node1 and node2 on partitions
	node1_link.resize(ntrees);
	node2_link.resize(ntrees);
	// Neighbors of node1_link and node2_link to be swapped during NNI
	node1_link_nei.resize(ntrees);
	node2_link_nei.resize(ntrees);
	// iterators for the neighbors of node1_link and node2_link to be swapped
	node1_link_it.resize(ntrees);
	node2_link_it.resize(ntrees);

	/*===========================================================================================
	 * 	MAIN:
	 * 	- do the NNI swap on SuperTree and perform the corresponding actions on SubTrees;
	 *	- compute the likelihood of swapped topology;
	 *  - swap back;
	 *	- restore if necessary.
	 *===========================================================================================*/
	int cnt;
	for (cnt = 0; cnt < 2; cnt++) if (nniMoves[cnt].node1) // only if nniMove satisfy constraint 
    {
		//cout<<"NNI Loop-----------------------------NNI."<<cnt<<endl;

    	NeighborVec::iterator node1_it = nniMoves[cnt].node1Nei_it;
    	NeighborVec::iterator node2_it = nniMoves[cnt].node2Nei_it;
        Neighbor *node1_nei = *node1_it;
        Neighbor *node2_nei = *node2_it;

		//node2_it = node2_its[cnt];
		//Neighbor *node2_nei = *node2_it;

		// Define which nodes/neighbors to be swapped on SubTree ----------------------------
		for(part=0; part<ntrees; part++)
			if(is_nni[part]==NNI_NO_EPSILON){
				node1_link[part]     = nei2_new->link_neighbors[part]->getNode();
				node2_link[part]     = nei1_new->link_neighbors[part]->getNode();
				node1_link_nei[part] = ((SuperNeighbor*)node1_nei)->link_neighbors[part];
				node1_link_it[part] = node1_link[part]->findNeighborIt(node1_link_nei[part]->node);
				node2_link_nei[part] = ((SuperNeighbor*)node2_nei)->link_neighbors[part];
				node2_link_it[part] = node2_link[part]->findNeighborIt(node2_link_nei[part]->node);
			}

		// Do the NNI swap on SuperTrees ----------------------------------------------------
		node1->updateNeighbor(node1_it, node2_nei);
		node2_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node1_nei);
		node1_nei->node->updateNeighbor(node1, node2);

		// Perform actions in accordance with the type of NNI for a given partition ---------
		for(part = 0; part < ntrees; part++){
			//cout<<"Partition: "<<part<<endl;

			if(is_nni[part]==NNI_NO_EPSILON){
				//cout<<part<<"- NO_EPS: do NNI swap"<<endl;
				//allNNIcases_computed[0] += 1;

                // reorient partial_lh before swap
                at(part)->reorientPartialLh(node1_link[part]->findNeighbor(node2_link[part]), node1_link[part]);
                at(part)->reorientPartialLh(node2_link[part]->findNeighbor(node1_link[part]), node2_link[part]);

				// Do NNI swap on partition
				node1_link[part]->updateNeighbor(node1_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				for(id=0; id<IT_NUM; id++){
					((PhyloNeighbor*)(*sub_saved_it[part*6+id]))->clearPartialLh();
				}
				//checkBranchLen();
			} else if(is_nni[part]==NNI_MANY_EPSILON){
				//cout<<part<<"- MANY_EPS: do nothing"<<endl;
				// the NNI on SuperTree does not change anything on SubTree

			} else if(is_nni[part]==NNI_THREE_EPSILON){
				//cout<<part<<"- THREE_EPS: relink"<<endl;

				// The central branch had no image before the NNI.
				// Relink the central branch and take care of branch lengths.
				// In the end restore one branch (valid for both nni1 and nni5).

				linkBranch(part, nei1_new, nei2_new);
				ASSERT(nei1_new->link_neighbors[part]);

				// Save the branch length
				if(cnt == 0)
					sub_saved_branch[6*part] = nei1_new->link_neighbors[part]->length;

				nei1_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;

				// since the branch length was changed we have to recompute the likelihood of the branch
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei1_new->link_neighbors[part],
						nei2_new->link_neighbors[part]->getNode());

			}else if(is_nni[part]==NNI_TWO_EPSILON){
				//cout<<part<<"- TWO_EPS: relink"<<endl;

				/* In fact, before relinking the image of central branch is NULL (because we allocated
				 * new SuperNeighbor and filled the link_neighbors with NULL for all partitions).
				 * After relinking it can be either NULL or it should relink to the same branch as before.
				 * In the end restore one branch (valid for both nni1 and nni5).*/

				// Save the branch length
				if(cnt == 0)
					sub_saved_branch[6*part] = saved_nei[0]->link_neighbors[part]->length;

				linkBranch(part, nei1_new, nei2_new);
				if(!nei1_new->link_neighbors[part]){
					saved_nei[0]->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
					saved_nei[1]->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
					part_info[part].cur_score = at(part)->computeLikelihoodBranch(saved_nei[0]->link_neighbors[part],
							saved_nei[1]->link_neighbors[part]->getNode());
				}

			}else if(is_nni[part] == NNI_ONE_EPSILON){
				//cout<<part<<"- ONE_EPS: relink, update the link_neighbors"<<endl;

				/* The crazy case, which absorbs most of the bugs:(
				 * Lets say on SuperTree there are five branches, a,b,c,d and central e, and d has an empty image.
				 * The corresponding SubTree has 3 branches, a',b',c'.
				 * Before NNI central branch, e, has an image. Lets say it maps to a'.
				 * After NNI it will be remapped either to b' or c', depending on which nodes will be swapped.
				 * Update the corresponding link_neighbors. Make sure that link_neighbors for central branch e
				 * and for the one it is now mapped to (b' or c'), are the same.
				 * Decrease a' (done before). Increase b' or c' depending on the NNI. Restore three branches.*/

				linkBranch(part, nei1_new, nei2_new);
				ASSERT(nei1_new->link_neighbors[part]);

				//cout<<"nei1_new->link_nei["<<part<<"]->node->id"<<nei1_new->link_neighbors[part]->node->id<<endl;
				//cout<<"nei2_new->link_nei["<<part<<"]->node->id"<<nei2_new->link_neighbors[part]->node->id<<endl;

				ASSERT(nei1_new->link_neighbors[part]->node->findNeighbor(nei2_new->link_neighbors[part]->node));
				ASSERT(nei2_new->link_neighbors[part]->node->findNeighbor(nei1_new->link_neighbors[part]->node));

				// nni1:
				// - you need to update only one link_neighbor with new PhyloNeighbor
				//	 (either node1->findNeighbor(node2) or node2->findNeighbor(node1))
				// - the second is already linked to some existing PhyloNeighbor after linkBranch().
				for(id=2; id<6; id++){
					if(node2->isNeighbor(neighbor_nodes[id-2])){
						// nei2_new should be updated
						if(((SuperNeighbor*)node2->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]){
							//cout<<"node2: "<<"id = "<<id<<"; node_id = "<<((SuperNeighbor*)(node2->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]->node->id<<";"<<endl;
							if(((SuperNeighbor*)node2->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]->node
									== nei1_new->link_neighbors[part]->node){
								//assert(((SuperNeighbor*)node2->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]->node->id
									//	== (*sub_saved_it[part*6 + id])->node->id);
								nei2_new->link_neighbors[part] = (PhyloNeighbor*)(*sub_saved_it[part*6 + id]);
								//cout<<"   nei id = "<<id<<"; node_id = "<<((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node->id<<";"<<endl;
								//cout<<"   sub "<<"id = "<<id<<"; node_id = "<<(*sub_saved_it[part*6 + id])->node->id<<";"<<endl;
								break;
							}
						}
					} else {
						// nei1_new should be updated
						ASSERT(node1->isNeighbor(neighbor_nodes[id-2]));
						if(((SuperNeighbor*)node1->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]){
							//cout<<"node1: "<<"id = "<<id<<"; node_id = "<<((SuperNeighbor*)(node1->findNeighbor(neighbor_nodes[id-2])))->link_neighbors[part]->node->id<<";"<<endl;
							if(((SuperNeighbor*)node1->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]->node
									== nei2_new->link_neighbors[part]->node){
								//assert(((SuperNeighbor*)node1->findNeighbor(neighbor_nodes[id-2]))->link_neighbors[part]->node->id
									//	== (*sub_saved_it[part*6 + id])->node->id);
								nei1_new->link_neighbors[part] = (PhyloNeighbor*)(*sub_saved_it[part*6 + id]);
								//cout<<"   nei id = "<<id<<"; node_id = "<<((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node->id<<";"<<endl;
								//cout<<"   sub "<<"id = "<<id<<"; node_id = "<<(*sub_saved_it[part*6 + id])->node->id<<";"<<endl;
								break;
							}
						}
					}
				}

				// Clear partial likelihoods for all three neighbors nei1/2->find(node1/2)
				if(params->nni5 && cnt == 1){
					for(id=2; id<6; id++){
						if(id != id_eps[part]){
							((PhyloNeighbor*)(*sub_saved_it[part*6 + id]))->clearPartialLh();
						}
					}
				}
				//cout<<"nei1_new->link_nei["<<part<<"]->node->id"<<nei1_new->link_neighbors[part]->node->id<<endl;
				//cout<<"nei2_new->link_nei["<<part<<"]->node->id"<<nei2_new->link_neighbors[part]->node->id<<endl;

				ASSERT(nei1_new->link_neighbors[part]->node->findNeighbor(nei2_new->link_neighbors[part]->node));
				ASSERT(nei2_new->link_neighbors[part]->node->findNeighbor(nei1_new->link_neighbors[part]->node));

				// Increase the branch to which the central is relinked.
				nei1_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;

			} // end of else ONE_EPS case
		} // end of part loop

/*===============================================================================================================================*
 * 											Compute the score of the swapped topology 				  							 *
 *===============================================================================================================================*/
		//cout<<"Before optimization"<<endl;
		//mapBranchLen();
		//checkBranchLen();

		optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
//		double score = computeLikelihoodFromBuffer();
		node1->findNeighbor(node2)->getLength(nniMoves[cnt].newLen[0]);

//		if (verbose_mode >= VB_MED) {
//			cout << "After_nni1 [" << score << "] ";
//			printTree(cout);
//			cout << endl;
//    		//for(part = 0; part < ntrees; part++)
//    		//	cout << is_nni[part] << " ";
//    		//cout << endl;
//    		//cout<<"NNI count = "<<cnt<<endl;
//		}
		//cout<<"After optimization"<<endl;
		//checkBranchLen();

		// %%%%%%%%%%%%%%%%%%%%%%%%  FIVE BRANCH OPTIMIZATION  %%%%%%%%%%%%%%%%%%%%%%%%
		i=1;
	    if (params->nni5) {

	    	// ------ Optimization of branches incident to node1 ---------------
	    	FOR_EACH_SUPER_NEIGHBOR(node1, node2, it, nei_of_node1){
	    		// Clear the partial likelihood of node1 neighbor: only for NO or ONE epsilon cases
	    		for(part = 0; part < ntrees; part++)
	    			if (nei_of_node1->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || is_nni[part]==NNI_ONE_EPSILON)){
	    				node_link = nei_of_node1->link_neighbors[part]->node;
	    				nei_link  = nei2_new->link_neighbors[part]->node; // this should be node 1 on subtree
	    				// the problem is that for ONE_epsilon case node1 on subtree is equal to its neighbor node on subtree
	    				// in this case we have to set nei_link to node2 on subtree
	    				if(node_link->id == nei_link->id){
	    					nei_link = nei1_new->link_neighbors[part]->node;
	    				}
	    				//cout<<"HERE it is: "<<((SuperNeighbor*)(*it))->link_neighbors[part]->node->id<<endl;
	    				//cout<<nei2_new->link_neighbors[part]->node->id<<endl;
						((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
						//cout<<"CASE:"<<is_nni[part]<<"Cleared partial likelihood"<<endl;
	    			}
	    		// Optimize the branch incident to node1
	    		//cout<<"NNI5 : node1 : Before optimization"<<endl;
	    		//checkBranchLen();
	    		optimizeOneBranch(node1, nei_of_node1->getNode(), false, NNI_MAX_NR_STEP);
				node1->findNeighbor(nei_of_node1->node)->getLength(nniMoves[cnt].newLen[i]);
				i++;


	    		//cout<<"NNI5 : node1 : After optimization"<<endl;
	    		//checkBranchLen();
	    	}

	    	// ------ Clear the partial likelihood on the central branch -------
	    	for(part = 0; part < ntrees; part++)
	    		if(((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || is_nni[part]==NNI_ONE_EPSILON)){
	    			((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part]->clearPartialLh();
	    		}

	    	// ------ Optimization of branches incident to node2 ---------------
	    	FOR_NEIGHBOR(node2, node1, it){
	    		// Clear the partial likelihood of node2 neighbor: only for NO or ONE epsilon cases
	    		for(part = 0; part < ntrees; part++){
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || is_nni[part]==NNI_ONE_EPSILON)){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei1_new->link_neighbors[part]->node;
	    				if(node_link->id == nei_link->id){
	    					nei_link = nei2_new->link_neighbors[part]->node;
	    				}
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    				//cout<<"CASE:"<<is_nni[part]<<"Cleared partial likelihood"<<endl;
	    			}
	    		}
	    		// Optimize the branch incident to node2
	    		optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false, NNI_MAX_NR_STEP);
				node2->findNeighbor((*it)->node)->getLength(nniMoves[cnt].newLen[i]);
				i++;
	    	}
	    }

		double score = computeLikelihoodFromBuffer();
		if (verbose_mode >= VB_DEBUG)
			cout << "Log-likelihood: " << score << endl;

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%  END of nni5branch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

		nniMoves[cnt].newloglh = score;

//	    if (verbose_mode >= VB_MED) {
//			//this->clearAllPartialLH();
//			//for(part = 0; part<ntrees; part++){
//			//	at(part)->clearAllPartialLH();
//			//}
//			//cout << "[" << this->computeLikelihood() << "] ";
//			cout << "After_nni5 " << score << " ";
//			printTree(cout);
//			cout << endl;
//			for(part = 0; part < ntrees; part++)
//				cout << is_nni[part] << " ";
//			cout << endl;
//		}

		// FOR SH-aLRT test
		if (nniMoves[cnt].ptnlh)
			computePatternLikelihood(nniMoves[cnt].ptnlh, &score);

	    // Save current tree for ufboot analysis
	    if (save_all_trees == 2) {
	    		saveCurrentTree(score);
	    }

//	    // *************************** STORE INFO ABOUT NNI ***************************
//
//	    // Store information about this NNI for NNImove for SuperTree
//		if (nni_param) {
////			if (verbose_mode >= VB_MAX)
////				printTree(cout, WT_BR_LEN + WT_NEWLINE);
//			if (cnt == 0) {
//				nni_param->nni1_score = score;
//				nni_param->nni1_brlen = nei1_new->length;
//			} else {
//				nni_param->nni2_score = score;
//				nni_param->nni2_brlen = nei1_new->length;
//			}
//		}
//		// ***************************************************************************

		// =============================== RESTORE INFO ==============================
		// Restore the cur_score for partitions
		for (i = 0; i < part_info.size(); i++)
			part_info[i].cur_score = saved_cur_score[i];

		//Restoring branch length on Super Tree
		nei1_new->length = old_brlen;
		nei2_new->length = old_brlen;

// Swap back on SuperTree --------------------------------------------------------------------------------------------------------------
		node1->updateNeighbor(node1_it, node1_nei);
		node1_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node2_nei);
		node2_nei->node->updateNeighbor(node1, node2);

		//Restoring 4 branches around central?
		if(params->nni5){
			for(id=2;id<6;id++){
				(*saved_it[id])->length = saved_nei[id]->length;
				(*node_nei_it[id-2])->length = saved_nei[id]->length;
			}
		}

// Swap back or relink back on SubTrees------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){

			if(is_nni[part]==NNI_NO_EPSILON){

                // reorient partial_lh before swap
                at(part)->reorientPartialLh((PhyloNeighbor*)node1_link[part]->findNeighbor(node2_link[part]), node1_link[part]);
                at(part)->reorientPartialLh((PhyloNeighbor*)node2_link[part]->findNeighbor(node1_link[part]), node2_link[part]);

				node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				//Restoring the branch length on the SubTree
				node1_link[part]->findNeighbor(node2_link[part])->length = saved_nei[0]->link_neighbors[part]->length;
				node2_link[part]->findNeighbor(node1_link[part])->length = saved_nei[0]->link_neighbors[part]->length;

				if(params->nni5){
					for(id = 2; id < 6; id++){
						((SuperNeighbor*)(*saved_it[id]))->link_neighbors[part]->length = saved_nei[id]->link_neighbors[part]->length;
						((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->length = saved_nei[id]->link_neighbors[part]->length;
					}
				}
				//mapBranchLen();

			} else if(is_nni[part]==NNI_ONE_EPSILON){
				//linkCheckRe(part,node1,node2,sub_saved_nei2[part],sub_saved_nei1[part]);
				//linkCheckRe(part,node2,node1,sub_saved_nei1[part],sub_saved_nei2[part]);

				// Relink back
				linkBranch(part, nei1_new, nei2_new);
				ASSERT(nei1_new->link_neighbors[part]->node == saved_nei[0]->link_neighbors[part]->node);
				ASSERT(nei2_new->link_neighbors[part]->node == saved_nei[1]->link_neighbors[part]->node);

				// Restore three branches
				for(id=2; id<6; id++){
					if(((SuperNeighbor*)*saved_it[id])->link_neighbors[part]){
						(*sub_saved_it[part*6+id])->length = sub_saved_branch[6*part + id];
						((SuperNeighbor*)*saved_it[id])->link_neighbors[part]->length = sub_saved_branch[6*part + id];
						((SuperNeighbor*)*node_nei_it[id-2])->link_neighbors[part]->length = sub_saved_branch[6*part + id];
					}
				}

			} else if(is_nni[part]==NNI_THREE_EPSILON){
				nei1_new->link_neighbors[part]->length = sub_saved_branch[6*part];
				nei2_new->link_neighbors[part]->length = sub_saved_branch[6*part];
				//linkBranch(part, nei1_new, nei2_new);
			} else if(is_nni[part]==NNI_TWO_EPSILON){
				//linkBranch(part, nei1_new, nei2_new);
				saved_nei[0]->link_neighbors[part]->length = sub_saved_branch[6*part];
				saved_nei[1]->link_neighbors[part]->length = sub_saved_branch[6*part];
				nei1_new->link_neighbors[part] = NULL;
				nei2_new->link_neighbors[part] = NULL;
			} else if(is_nni[part]==NNI_MANY_EPSILON){
				// There is no need to restore anything
			}
		}

		//cout<<"in NNI1end ---- logL = "<<this->computeLikelihood()<<endl;
	} // end of for(cnt)


//=============================================================================================================================================================
// 							Restoring after 2 NNIs
//=============================================================================================================================================================
// Restoring information for SuperTree ------------------------------------------------------------------------------------------------------------
	// restore the Neighbors*
	for (id = IT_NUM-1; id >= 0; id--) {
		if (*saved_it[id] == current_it) current_it = (SuperNeighbor*) saved_nei[id];
		if (*saved_it[id] == current_it_back) current_it_back = (SuperNeighbor*) saved_nei[id];

		delete (*saved_it[id]);
		(*saved_it[id]) = saved_nei[id];
	 }
	// restore the length of 4 branches around node1, node2
	// since you have restored the neighbors and by this also the correct branch lengths,
	// now just restore branch lengths of the second corresponding neighbors
	FOR_NEIGHBOR(node1, node2, it)
		(*it)->length = (*it)->node->findNeighbor(node1)->length;
	FOR_NEIGHBOR(node2, node1, it)
		(*it)->length = (*it)->node->findNeighbor(node2)->length;

// Restoring information for SubTrees ------------------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){
			if(is_nni[part] == NNI_NO_EPSILON){
				// restore the Neighbors*
				for (i = IT_NUM-1; i >= 0; i--) {
					if((*sub_saved_it[part*6+i])){
//						aligned_free(((PhyloNeighbor*) *sub_saved_it[part*6+i])->scale_num);
//						aligned_free(((PhyloNeighbor*) *sub_saved_it[part*6+i])->partial_lh);
						if (*sub_saved_it[part*6+i] == at(part)->current_it) at(part)->current_it = saved_nei[i]->link_neighbors[part];
						if (*sub_saved_it[part*6+i] == at(part)->current_it_back) at(part)->current_it_back = saved_nei[i]->link_neighbors[part];

						delete (*sub_saved_it[part*6+i]);
						(*sub_saved_it[part*6+i]) = saved_nei[i]->link_neighbors[part];
					}
				}
				// restore the length of 4 branches around node1_link[part], node2_link[part]
				node1_link[part] = saved_nei[1]->link_neighbors[part]->getNode();
				node2_link[part] = saved_nei[0]->link_neighbors[part]->getNode();
				FOR_NEIGHBOR(node1_link[part], node2_link[part], it)
					(*it)->length = (*it)->node->findNeighbor(node1_link[part])->length;
				FOR_NEIGHBOR(node2_link[part], node1_link[part], it)
					(*it)->length = (*it)->node->findNeighbor(node2_link[part])->length;

			} else if(is_nni[part] == NNI_ONE_EPSILON){

				// Delete the allocated neighbors and restore from saved neighbors
				for (id = 5; id >= 2; id--) {
					//if((*sub_saved_it[part*6+id])){
					if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
//						aligned_free(((PhyloNeighbor*) *sub_saved_it[part*6+id])->scale_num);
//						aligned_free(((PhyloNeighbor*) *sub_saved_it[part*6+id])->partial_lh);

						// It was commented, not sure why.. Just keep in mind------------------
						if (*sub_saved_it[part*6+id] == at(part)->current_it)
							at(part)->current_it = saved_nei[id]->link_neighbors[part];
						if (*sub_saved_it[part*6+id] == at(part)->current_it_back)
							at(part)->current_it_back = saved_nei[id]->link_neighbors[part];
						//---------------------------------------------------------------------

						delete (*sub_saved_it[part*6+id]);
						(*sub_saved_it[part*6+id]) = ((SuperNeighbor*)(*saved_it[id]))->link_neighbors[part];
					}
				}
				// Increase the central branch, since the length that was saved, was decreased
				saved_nei[0]->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				saved_nei[1]->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
			}
		}
		//mapBranchLen();
		//cout<<"In the end of swap NNI"<<endl;
		//checkBranchLen();
//------------------------------------------------------------------------------------------------------------------------------------------------
	//if(score_mine != this->computeLikelihood())
	//	cout<<"Something weird happens during NNI evaluation..." << score_mine << " " << computeLikelihood() <<endl;

	delete [] saved_cur_score;
	delete [] sub_saved_branch;
	return cur_score;
}

void PhyloSuperTreePlen::linkCheck(int part,Node* node, Node* dad, PhyloNeighbor* saved_link_dad_nei){
	NeighborVec::iterator it;
	SuperNeighbor *dad_nei = (SuperNeighbor*)dad->findNeighbor(node);
	SuperNeighbor *node_nei = (SuperNeighbor*)node->findNeighbor(dad);
	FOR_NEIGHBOR(node, dad, it){
		if(((SuperNeighbor*)(*it))->link_neighbors[part] == saved_link_dad_nei){
			((SuperNeighbor*)(*it))->link_neighbors[part] = dad_nei->link_neighbors[part];
			((SuperNeighbor*)((*it)->node->findNeighbor(node)))->link_neighbors[part] = node_nei->link_neighbors[part];
			linkCheck(part, (*it)->node, node, saved_link_dad_nei);
		}
	}
}

void PhyloSuperTreePlen::linkCheckRe(int part,Node* node, Node* dad, PhyloNeighbor* saved_link_dad_nei,PhyloNeighbor* saved_link_node_nei){
	NeighborVec::iterator it;
	FOR_NEIGHBOR(node, dad, it){
		if(((SuperNeighbor*)(*it))->link_neighbors[part] == ((SuperNeighbor*)dad->findNeighbor(node))->link_neighbors[part]){
			linkCheckRe(part, (*it)->node, node, saved_link_dad_nei, saved_link_node_nei);
			((SuperNeighbor*)(*it))->link_neighbors[part] = saved_link_dad_nei;
			((SuperNeighbor*)((*it)->node->findNeighbor(node)))->link_neighbors[part] = saved_link_node_nei;
		}
	}
}
//void PhyloSuperTreePlen::restoreAllBrans(PhyloNode *node, PhyloNode *dad) {
//	IQTree::restoreAllBrans(node,dad);
//	mapTrees();
//}

bool PhyloSuperTreePlen::checkBranchLen(){

//	NodeVector nodes1,nodes2;
//	int i;
//	getBranches(nodes1, nodes2);
//	double *checkVAL = new double[branchNum];
//	for(int part = 0; part < size(); part++){
//		memset(checkVAL, 0, at(part)->branchNum*sizeof(double));
//		for (i = 0; i < nodes1.size(); i++){
//			if(((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part])
//				checkVAL[((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part]->id] += nodes1[i]->findNeighbor(nodes2[i])->length * part_info[part].part_rate;
//		}
//		NodeVector nodes1_sub, nodes2_sub;
//		at(part)->getBranches(nodes1_sub, nodes2_sub);
//		for(int j = 0; j<nodes1_sub.size();j++)
//			if(fabs(nodes1_sub[j]->findNeighbor(nodes2_sub[j])->length-checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id])>0.0001){
//				//drawTree(cout, WT_BR_SCALE + WT_INT_NODE + WT_BR_LEN);
//				printMapInfo();
//				cout<<endl;
//				cout<<"Partition = "<<part<<", Branch id = "<<nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id<<endl;
//				outError("Branches on SuperTree and SubTree do not match!!",__func__);
//			}
//
//	}
//	delete [] checkVAL;

	return true;
}

void PhyloSuperTreePlen::mapBranchLen()
{
	NodeVector nodes1,nodes2;
	int i;
	getBranches(nodes1, nodes2);
	double *checkVAL = new double[branchNum];
	for(int part = 0; part < size(); part++){
		memset(checkVAL,0,at(part)->branchNum*sizeof(double));
		for (i = 0; i < nodes1.size(); i++){
			if(((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part])
				checkVAL[((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part]->id] +=
						nodes1[i]->findNeighbor(nodes2[i])->length * part_info[part].part_rate;
		}
		NodeVector nodes1_sub, nodes2_sub;
		at(part)->getBranches(nodes1_sub, nodes2_sub);
		for(int j = 0; j<nodes1_sub.size();j++){
			nodes1_sub[j]->findNeighbor(nodes2_sub[j])->length = checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id];
			nodes2_sub[j]->findNeighbor(nodes1_sub[j])->length = checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id];
		}
	}
	delete [] checkVAL;
}

void PhyloSuperTreePlen::mapBranchLen(int part)
{
	NodeVector nodes1,nodes2;
	int i;
	getBranches(nodes1, nodes2);
	double *checkVAL = new double[branchNum];
    memset(checkVAL,0,at(part)->branchNum*sizeof(double));
    for (i = 0; i < nodes1.size(); i++){
        if(((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part])
            checkVAL[((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part]->id] +=
                    nodes1[i]->findNeighbor(nodes2[i])->length * part_info[part].part_rate;
    }
    NodeVector nodes1_sub, nodes2_sub;
    at(part)->getBranches(nodes1_sub, nodes2_sub);
    for(int j = 0; j<nodes1_sub.size();j++){
        nodes1_sub[j]->findNeighbor(nodes2_sub[j])->length = checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id];
        nodes2_sub[j]->findNeighbor(nodes1_sub[j])->length = checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id];
    }
	delete [] checkVAL;
}

void PhyloSuperTreePlen::printMapInfo() {
	NodeVector nodes1, nodes2;
	getBranches(nodes1, nodes2);
	int part = 0;
	drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_LEN);
	for (iterator it = begin(); it != end(); it++, part++) {
		cout << "Subtree for partition " << part << endl;
		(*it)->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_LEN);
		for (int i = 0; i < nodes1.size(); i++) {
			PhyloNeighbor *nei1 = ((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part];
			PhyloNeighbor *nei2 = ((SuperNeighbor*)nodes2[i]->findNeighbor(nodes1[i]))->link_neighbors[part];
			cout << nodes1[i]->findNeighbor(nodes2[i])->id << ":";
			if (nodes1[i]->isLeaf()) cout << nodes1[i]->name; else cout << nodes1[i]->id;
			cout << ",";
			if (nodes2[i]->isLeaf()) cout << nodes2[i]->name; else cout << nodes2[i]->id;
			cout <<"("<<nodes1[i]->findNeighbor(nodes2[i])->length<<")"<< " -> ";
			if (nei2) {
				cout << nei2->id << ":";
				if (nei2->node->isLeaf())
					cout << nei2->node->name;
				else cout << nei2->node->id;
			}
			else cout << -1;
			cout << ",";
			if (nei1){
				if (nei1->node->isLeaf())
					cout << nei1->node->name;
				else cout << nei1->node->id;
				cout <<"("<<nei1->length<<")";
			}
			else cout << -1;
			cout << endl;
		}
	}
}

void PhyloSuperTreePlen::initPartitionInfo() {
	//PhyloSuperTree::initPartitionInfo();
	for (int part = 0; part < size(); part++){
		if(part_info[part].part_rate == 0.0) { part_info[part].part_rate = 1.0; }
		part_info[part].cur_score = 0.0;
	}
}

void PhyloSuperTreePlen::printNNIcasesNUM(){
	cout<<"For each \"NNI case\" on subtree the number of times it appeared during NNI evaluation:"<<endl;
	cout<<"Case 1: NO_EPS    = "<<allNNIcases_computed[0]<<endl;
	cout<<"Case 2: ONE_EPS   = "<<allNNIcases_computed[1]<<endl;
	cout<<"Case 3: TWO_EPS   = "<<allNNIcases_computed[2]<<endl;
	cout<<"Case 4: THREE_EPS = "<<allNNIcases_computed[3]<<endl;
	cout<<"Case 5: MANY_EPS  = "<<allNNIcases_computed[4]<<endl;
}

void PhyloSuperTreePlen::computeBranchLengths()
{
	}

int PhyloSuperTreePlen::fixNegativeBranch(bool force, PhyloNode *node, PhyloNode *dad) {

	mapTrees();
	int fixed = 0;
	for (iterator it = begin(); it != end(); it++) {
		(*it)->initializeAllPartialPars();
		(*it)->clearAllPartialLH();
		fixed += (*it)->fixNegativeBranch(force);
		(*it)->clearAllPartialLH();
	}
    // FOR OLGA: because this check is not performed, branch lengths of user tree will change even with -fixbr command line
    if (fixed) {
        PhyloSuperTree::computeBranchLengths();

        // it is necessary to map the branch lengths from supertree into gene trees!
        mapTrees();
    }

	return fixed;
}

void PhyloSuperTreePlen::changeNNIBrans(NNIMove &nnimove) {

	PhyloTree::changeNNIBrans(nnimove);
	//mapBranchLen();

}


/**
        initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
 */
void PhyloSuperTreePlen::initializeAllPartialLh() {
	iterator it;
	int part, partid;
	int ntrees = size();

	block_size.resize(ntrees);
	scale_block_size.resize(ntrees);

	vector<uint64_t> mem_size, lh_cat_size, buffer_size;
	mem_size.resize(ntrees);
	lh_cat_size.resize(ntrees);
    buffer_size.resize(ntrees);

	uint64_t
        total_mem_size = 0,
        total_block_size = 0,
        total_scale_block_size = 0,
        total_lh_cat_size = 0,
        total_buffer_size = 0;

	if (part_order.empty())
		computePartitionOrder();

	for (partid = 0; partid < ntrees; partid++) {
		part = part_order[partid];
        it = begin() + part;
        // extra #numStates for ascertainment bias correction
		mem_size[part] = get_safe_upper_limit((*it)->getAlnNPattern()) + get_safe_upper_limit((*it)->aln->num_states);
        size_t mem_cat_size = mem_size[part] * (*it)->getRate()->getNRate() *
				(((*it)->model_factory->fused_mix_rate)? 1 : (*it)->getModel()->getNMixtures());

		block_size[part] = mem_cat_size * (*it)->aln->num_states;
		scale_block_size[part] = mem_cat_size;

		lh_cat_size[part] = mem_size[part] * (*it)->getRate()->getNDiscreteRate() *
				(((*it)->model_factory->fused_mix_rate)? 1 : (*it)->getModel()->getNMixtures());
		total_mem_size += mem_size[part];
		total_block_size += block_size[part];
        total_scale_block_size += scale_block_size[part];
		total_lh_cat_size += lh_cat_size[part];
        total_buffer_size += (buffer_size[part] = (*it)->getBufferPartialLhSize());
	}

    if (!_pattern_lh)
        _pattern_lh = aligned_alloc<double>(total_mem_size);
    at(part_order[0])->_pattern_lh = _pattern_lh;
    if (!_pattern_lh_cat)
        _pattern_lh_cat = aligned_alloc<double>(total_lh_cat_size);
    at(part_order[0])->_pattern_lh_cat = _pattern_lh_cat;
    if (!theta_all)
        theta_all = aligned_alloc<double>(total_block_size);
    if (!buffer_scale_all)
        buffer_scale_all = aligned_alloc<double>(total_mem_size);
    if (!buffer_partial_lh)
        buffer_partial_lh = aligned_alloc<double>(total_buffer_size);
    at(part_order[0])->theta_all = theta_all;
    at(part_order[0])->buffer_scale_all = buffer_scale_all;
    at(part_order[0])->buffer_partial_lh = buffer_partial_lh;
    if (!ptn_freq) {
        ptn_freq = aligned_alloc<double>(total_mem_size);
        ptn_freq_computed = false;
    }
    if (ptn_freq_pars == nullptr) {
        ptn_freq_pars = aligned_alloc<UINT>(total_mem_size);
    }
    at(part_order[0])->ptn_freq = ptn_freq;
    at(part_order[0])->ptn_freq_pars = ptn_freq_pars;
    at(part_order[0])->ptn_freq_computed = false;
    if (!ptn_invar)
        ptn_invar = aligned_alloc<double>(total_mem_size);
    at(part_order[0])->ptn_invar = ptn_invar;

//    size_t IT_NUM = (params->nni5) ? 6 : 2;
    size_t IT_NUM = 2;
    if (!nni_partial_lh) {
        nni_partial_lh = aligned_alloc<double>(IT_NUM*total_block_size);
    }
    at(part_order[0])->nni_partial_lh = nni_partial_lh;
    
    if (!nni_scale_num) {
        nni_scale_num = aligned_alloc<UBYTE>(IT_NUM*total_scale_block_size);
    }
    at(part_order[0])->nni_scale_num = nni_scale_num;

	for (partid = 1; partid < ntrees; partid++) {
        part = part_order[partid-1];
        it = begin() + part_order[partid];
        iterator prev_it = begin()+part_order[partid-1];
		(*it)->_pattern_lh = (*prev_it)->_pattern_lh + mem_size[part];
		(*it)->_pattern_lh_cat = (*prev_it)->_pattern_lh_cat + lh_cat_size[part];
		(*it)->theta_all = (*prev_it)->theta_all + block_size[part];
        (*it)->buffer_scale_all = (*prev_it)->buffer_scale_all + mem_size[part];
        (*it)->buffer_partial_lh = (*prev_it)->buffer_partial_lh + buffer_size[part];
		(*it)->ptn_freq = (*prev_it)->ptn_freq + mem_size[part];
        (*it)->ptn_freq_pars = (*prev_it)->ptn_freq_pars + mem_size[part];
		(*it)->ptn_freq_computed = false;
		(*it)->ptn_invar = (*prev_it)->ptn_invar + mem_size[part];
        (*it)->nni_partial_lh = (*prev_it)->nni_partial_lh + IT_NUM*block_size[part];
        (*it)->nni_scale_num = (*prev_it)->nni_scale_num + IT_NUM*scale_block_size[part];
	}

	// compute total memory for all partitions
	uint64_t total_partial_lh_entries = 0, total_scale_num_entries = 0, total_partial_pars_entries = 0;
	partial_lh_entries.resize(ntrees);
	scale_num_entries.resize(ntrees);
	partial_pars_entries.resize(ntrees);
	for (it = begin(), part = 0; it != end(); it++, part++) {
		(*it)->getMemoryRequired(partial_lh_entries[part], scale_num_entries[part], partial_pars_entries[part]);
		total_partial_lh_entries += partial_lh_entries[part];
		total_scale_num_entries += scale_num_entries[part];
		total_partial_pars_entries += partial_pars_entries[part];
	}

	// allocate central memory for all partitions
	if (!central_partial_lh) {
        try {
        	central_partial_lh = aligned_alloc<double>(total_partial_lh_entries);
        	central_scale_num = aligned_alloc<UBYTE>(total_scale_num_entries);
        } catch (std::bad_alloc &ba) {
        	outError("Not enough memory for partial likelihood vectors (bad_alloc)");
        }
	}
//    if (!central_partial_pars) {
//        try {
//        	central_partial_pars = aligned_alloc<UINT>(total_partial_pars_entries);
//        } catch (std::bad_alloc &ba) {
//        	outError("Not enough memory for partial parsimony vectors (bad_alloc)");
//        }
//    }

    // assign individual chunk just to prevent reallocation of memory, they will not be used
	for (it = begin(); it != end(); it++) {
		(*it)->central_partial_lh = central_partial_lh;
		(*it)->central_scale_num = central_scale_num;
//		(*it)->central_partial_pars = central_partial_pars;
	}

	double *lh_addr = central_partial_lh;
	UBYTE *scale_addr = central_scale_num;
	UINT *pars_addr = central_partial_pars;
	clearAllPartialLH(true);

	initializeAllPartialLh(lh_addr, scale_addr, pars_addr);
    ASSERT((lh_addr - central_partial_lh) < total_partial_lh_entries*sizeof(double) && lh_addr > central_partial_lh);
    tip_partial_lh = NULL;
    tip_partial_pars = NULL;
    for (it = begin(), part = 0; it != end(); it++, part++) {
        (*it)->tip_partial_lh = lh_addr;
        (*it)->tip_partial_pars = pars_addr;
        uint64_t tip_partial_lh_size = (*it)->aln->num_states * ((*it)->aln->STATE_UNKNOWN+1) * (*it)->model->getNMixtures();
        uint64_t tip_partial_pars_size = (*it)->aln->num_states * ((*it)->aln->STATE_UNKNOWN+1);
        //tip_partial_lh_size = ((tip_partial_lh_size+3)/4)*4;
        lh_addr += get_safe_upper_limit(tip_partial_lh_size);
        pars_addr += get_safe_upper_limit_float(tip_partial_pars_size);
    }

    // 2016-09-29: redirect partial_lh when root does not occur in partition tree
    SuperNeighbor *root_nei = (SuperNeighbor*)root->neighbors[0];
    for (it = begin(), part = 0; it != end(); it++, part++) {
        if (root_nei->link_neighbors[part])
            continue;
        NodeVector nodes;
        (*it)->getInternalNodes(nodes);
        for (NodeVector::iterator nit = nodes.begin(); nit != nodes.end(); nit++) {
            bool has_partial_lh = false;
            FOR_NEIGHBOR_IT(*nit, NULL, neiit)
                if ( ((PhyloNeighbor*)(*neiit)->node->findNeighbor(*nit))->partial_lh) {
                    has_partial_lh = true;
                    break;
                }
            if (has_partial_lh)
                continue;
            // add partial_lh
            PhyloNeighbor *back_nei = (PhyloNeighbor*)(*nit)->neighbors[0]->node->findNeighbor(*nit);
            back_nei->partial_lh = lh_addr;
            back_nei->scale_num = scale_addr;
            lh_addr = lh_addr + block_size[part];
            scale_addr = scale_addr + scale_block_size[part];
        }
    }

}

void PhyloSuperTreePlen::initializeAllPartialLh(double* &lh_addr, UBYTE* &scale_addr, UINT* &pars_addr, PhyloNode *node, PhyloNode *dad) {
    if (!node)
        node = getRoot();
    if (dad) {
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        SuperNeighbor *nei = (SuperNeighbor*) node->findNeighbor(dad);
		SuperNeighbor *nei_back = (SuperNeighbor*) dad->findNeighbor(node);
        for (int partid = 0; partid < size(); partid++) {
            int part = part_order[partid];
        	PhyloNeighbor *nei_part = nei->link_neighbors[part];
        	if (!nei_part) continue;
        	PhyloNeighbor *nei_part_back = nei_back->link_neighbors[part];

            if (params->lh_mem_save == LM_PER_NODE) {
                if (!nei_part_back->node->isLeaf()) {
                    if (!nei_part_back->partial_lh) {
                        nei_part_back->partial_lh = lh_addr;
                        nei_part_back->scale_num = scale_addr;
                        lh_addr = lh_addr + block_size[part];
                        scale_addr = scale_addr + scale_block_size[part];
                    }
                } else {
                    nei_part_back->partial_lh = NULL;
                    nei_part_back->scale_num = NULL;
                }
//                nei_part->partial_lh = NULL;
//                nei_part->scale_num = NULL;
            } else {
                if (nei_part->node->isLeaf()) {
                    nei_part->partial_lh = NULL; // do not allocate memory for tip, use tip_partial_lh instead
                    nei_part->scale_num = NULL;
                } else if (!nei_part->partial_lh) {
                    nei_part->partial_lh = lh_addr;
                    nei_part->scale_num = scale_addr;
                    lh_addr = lh_addr + block_size[part];
                    scale_addr = scale_addr + scale_block_size[part];
                }
    //			nei_part->partial_pars = pars_addr;
    //			pars_addr += partial_pars_entries[part];

                nei_part = nei_back->link_neighbors[part];
                if (nei_part->node->isLeaf()) {
                    nei_part->partial_lh = NULL; // do not allocate memory for tip, use tip_partial_lh instead
                    nei_part->scale_num = NULL;
                } else if (!nei_part->partial_lh) {
                    nei_part->partial_lh = lh_addr;
                    nei_part->scale_num = scale_addr;
                    lh_addr = lh_addr + block_size[part];
                    scale_addr = scale_addr + scale_block_size[part];
                }
    //			nei_part->partial_pars = pars_addr;
    //			pars_addr += partial_pars_entries[part];
            }
        }
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, adj) {
        initializeAllPartialLh(lh_addr, scale_addr, pars_addr, adj, node);
    }
}

void PhyloSuperTreePlen::initializeAllPartialLh(int &index, int &indexlh,
                                                bool fullOn, PhyloNode *node, PhyloNode *dad) {
	// this function should not be used, assertion raised if accidentally called
	ASSERT(0);
}

void PhyloSuperTreePlen::reorientPartialLh(PhyloNeighbor* dad_branch, Node *dad) {
    SuperNeighbor *sdad_branch = (SuperNeighbor*) dad_branch;
    SuperNeighbor *snode_branch = (SuperNeighbor*) dad_branch->node->findNeighbor(dad);
    for (int part = 0; part < size(); part++) {
        if (sdad_branch->link_neighbors[part]) {
            at(part)->reorientPartialLh(sdad_branch->link_neighbors[part], snode_branch->link_neighbors[part]->node);
        }
    }
}


string PhyloSuperTreePlen::getTreeString() {
    return PhyloTree::getTreeString();
}

void PhyloSuperTreePlen::readTreeString(const string &tree_string) {
    PhyloTree::readTreeString(tree_string);

}
