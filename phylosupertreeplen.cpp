/*
 * phylosupertreeplen.cpp
 *
 *  Created on: Aug 5, 2013
 *      Author: olga
 */

#include "phylosupertreeplen.h"
#include "superalignmentpairwise.h"
#include <string.h>

/**********************************************************
 * class SuperAlignmentPairwisePlen
**********************************************************/


SuperAlignmentPairwisePlen::SuperAlignmentPairwisePlen(PhyloSuperTreePlen *atree, int seq1, int seq2)
 : SuperAlignmentPairwise((PhyloSuperTree*) atree, seq1, seq2)
{
	part_info = &(atree->part_info);
}

double SuperAlignmentPairwisePlen::computeFunction(double value) {
	int part = 0;
	double lh = 0.0;
	for (vector<AlignmentPairwise*>::iterator it = partitions.begin(); it != partitions.end(); it++, part++) {
		lh += (*it)->computeFunction(part_info->at(part).part_rate*value);
	}
	return lh;
}

double SuperAlignmentPairwisePlen::computeFuncDerv(double value, double &df, double &ddf) {
	int part = 0;
	double lh = 0.0;
	df = 0.0;
	ddf = 0.0;
	for (vector<AlignmentPairwise*>::iterator it = partitions.begin(); it != partitions.end(); it++, part++) {
		double d1, d2;
		lh += (*it)->computeFuncDerv(part_info->at(part).part_rate*value, d1, d2);
		df += part_info->at(part).part_rate*d1;
		ddf += part_info->at(part).part_rate*part_info->at(part).part_rate*d2;
	}
	return lh;
}

SuperAlignmentPairwisePlen::~SuperAlignmentPairwisePlen()
{}

/**********************************************************
 * class PartitionModelPlen
**********************************************************/


PartitionModelPlen::PartitionModelPlen()
        : PartitionModel()
{
	}

PartitionModelPlen::PartitionModelPlen(Params &params, PhyloSuperTreePlen *tree)
        : PartitionModel(params, tree)
{
	}

PartitionModelPlen::~PartitionModelPlen()
{
	}

double PartitionModelPlen::optimizeParameters(bool fixed_len, bool write_info, double epsilon) {
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    double tree_lh = 0.0, cur_lh = 0.0;
    int ntrees = tree->size();
    double tol = 0.2;

	/*#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh)
	#endif*/
    //tree->printMapInfo();
   /* for (int part = 0; part < ntrees; part++) {
    	cout << "LIKELIHOOD | Partition "<<part<<" | "<<tree->at(part)->computeLikelihood()<<endl;
    }*/
    tree_lh = tree->computeLikelihood();
	cout<<"Initial likelihood: "<<tree_lh<<endl;
    for(int i = 1; i < 100; i++){
    	tol = max(tol/2, epsilon);
    	cur_lh = 0.0;
    	for (int part = 0; part < ntrees; part++) {

    		// Subtree model parameters optimization
        	double model_lh = tree->at(part)->getModelFactory()->optimizeParameters(true,false,tol);
        	cur_lh += model_lh;
        	//cout <<"Partition "<<part<<" MODEL:"<<tree->at(part)->getModelName() <<endl;

    	}

    	// Optimizing gene rate
    	//tree->fixed_rates = true;
    	if(!tree->fixed_rates){
    		cur_lh = optimizeGeneRate(tol);
    		/*for (int part = 0; part < ntrees; part++){
    			cout<<"Partition "<<part<<" rate: "<<tree->part_info[part].part_rate<<endl;
    			//tree->at(part)->printTree(cout);
    			//cout<<endl;
    		}*/
    	}

    	// Optimizing branch lengths
    	int my_iter = min(5,i);

    	if(!fixed_len){
    		cur_lh = tree->optimizeAllBranches(my_iter,tol);
    	}
    	cout<<"Current likelihood at step "<<i<<": "<<cur_lh<<endl;
    	if(fabs(cur_lh-tree_lh) < epsilon)
    		break;

    	tree_lh = cur_lh;
    }
    cout <<"OPTIMIZE MODEL has finished"<< endl;
    for(int part = 0; part < ntrees; part++){
    	cout<<"PART RATE "<<part<<" = "<<tree->part_info[part].part_rate<<endl;
    }
    return tree_lh;
}

double PartitionModelPlen::optimizeGeneRate(double tol)
{
	PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
	int ndim = tree->size()-1;

	double *variables   = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool   *bound_check = new bool[ndim+1];
	int i;
	double score;

	// gene rates are optimized by BFGS algorithm

	setVariables(variables);

	for (i = 1; i <= ndim; i++) {
		//cout << variables[i] << endl;
		lower_bound[i] = 1e-4;
		upper_bound[i] = tree->size();
		bound_check[i] = false;
	}

	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, tol);

	getVariables(variables);

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}

double PartitionModelPlen::targetFunk(double x[]) {
	PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();

	double sum = 0.0;
	int part;
	for( part = 0; part < tree->size()-1; part ++){
		sum += x[part+1];
	}
	if (tree->size() - sum < 1e-4) return 1.0e+12;

	for( part = 0, sum = 0.0; part < tree->size(); part ++){
		double rate;
		if (part < tree->size() - 1)
			rate = x[part+1];
		else
			rate = tree->size() - sum;
		sum += rate;
		if(tree->part_info[part].part_rate != rate){
			tree->at(part)->clearAllPartialLH();
			tree->at(part)->scaleLength(rate/tree->part_info[part].part_rate);
			tree->part_info[part].part_rate = rate;
		}
	}
	//getVariables(x);

	return -tree->computeLikelihood();
}

void PartitionModelPlen::getVariables(double *variables) {
	PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
	int ntrees = tree->size()-1;
	double sum = 0.0;
	for(int part = 0; part < ntrees; part++){
		tree->part_info[part].part_rate = variables[part+1];
		sum += variables[part+1];
	}
	tree->part_info[ntrees].part_rate = tree->size() - sum;
}

void PartitionModelPlen::setVariables(double *variables) {
	PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
	int ntrees = tree->size()-1;
	for(int part = 0; part < ntrees; part++){
		variables[part+1] = tree->part_info[part].part_rate;
	}
}

int PartitionModelPlen::getNParameters() {
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
	int df = 0;
    for (PhyloSuperTreePlen::iterator it = tree->begin(); it != tree->end(); it++) {
    	df += (*it)->getModelFactory()->model->getNDim()+(*it)->getModelFactory()->site_rate->getNDim();
		if ( (*it)->getModelFactory()->model->freq_type == FREQ_EMPIRICAL) df +=  (*it)->getModelFactory()->model->num_states-1;
    }
    df += tree->branchNum;
    df += tree->size()-1;
    return df;
}

int PartitionModelPlen::getNDim(){
	PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
	int ndim = tree->size() -1;
	return ndim;
}


/**********************************************************
 * class PhyloSuperTreePlen
**********************************************************/


PhyloSuperTreePlen::PhyloSuperTreePlen()
: PhyloSuperTree()
{
	fixed_rates = false;
}

PhyloSuperTreePlen::PhyloSuperTreePlen(Params &params)
: PhyloSuperTree(params)
{
	fixed_rates = params.partition_fixed_rates;
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].part_rate = 1.0;
	}
}

PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment, PhyloSuperTree *super_tree)
: PhyloSuperTree(alignment,super_tree)
{
	fixed_rates = false;
}

PhyloSuperTreePlen::~PhyloSuperTreePlen()
{}


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
	assert(root);
	int part = 0;
	if (verbose_mode >= VB_DEBUG)
		drawTree(cout,  WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE | WT_BR_ID);
	for (iterator it = begin(); it != end(); it++, part++) {
		string taxa_set = ((SuperAlignment*)aln)->getPattern(part);
		(*it)->copyTree(this, taxa_set);

		// the only difference with PhyloSuperTree::mapTrees()
		(*it)->scaleLength(part_info[part].part_rate);

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

double PhyloSuperTreePlen::optimizeAllBranches(int my_iterations, double tolerance) {
	return PhyloTree::optimizeAllBranches(my_iterations,tolerance);
	cout<<"finished optimizing All"<<endl;
}

double PhyloSuperTreePlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH) {

	SuperNeighbor *nei1 = (SuperNeighbor*)node1->findNeighbor(node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)node2->findNeighbor(node1);

	for (int part = 0; part < size(); part++) {
		at(part)->theta_computed = false;
	}

	double tree_lh = PhyloTree::optimizeOneBranch(node1,node2,clearLH);

	if(clearLH){
		for (int part = 0; part < size(); part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			((PhyloNode*)nei1_part->node)->clearReversePartialLh(((PhyloNode*)nei2_part->node));
			((PhyloNode*)nei2_part->node)->clearReversePartialLh(((PhyloNode*)nei1_part->node));
		}
	}

	return tree_lh;
}

double PhyloSuperTreePlen::computeFunction(double value) {

	double tree_lh = 0.0;
	int ntrees = size();

	double lambda = value-current_it->length;
	current_it->length = value;
    current_it_back->length = value;

	SuperNeighbor *nei1 = (SuperNeighbor*)current_it_back->node->findNeighbor(current_it->node);
	SuperNeighbor *nei2 = (SuperNeighbor*)current_it->node->findNeighbor(current_it_back->node);
	assert(nei1 && nei2);

	for (int part = 0; part < ntrees; part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			if (part_info[part].cur_score == 0.0)
				part_info[part].cur_score = at(part)->computeLikelihood();
			if (nei1_part && nei2_part) {
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei2_part,(PhyloNode*)nei1_part->node);
				tree_lh += part_info[part].cur_score;
			} else {
				tree_lh += part_info[part].cur_score;
			}
		}
    return -tree_lh;
}

double PhyloSuperTreePlen::computeFuncDerv(double value, double &df, double &ddf) {
	double tree_lh = 0.0;
	double df_aux, ddf_aux;
	df = 0.0;
	ddf = 0.0;

	int ntrees = size();

	double lambda = value-current_it->length;
	current_it->length = value;
    current_it_back->length = value;

	SuperNeighbor *nei1 = (SuperNeighbor*)current_it_back->node->findNeighbor(current_it->node);
	SuperNeighbor *nei2 = (SuperNeighbor*)current_it->node->findNeighbor(current_it_back->node);
	assert(nei1 && nei2);

	for (int part = 0; part < ntrees; part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			if (part_info[part].cur_score == 0.0)
				part_info[part].cur_score = at(part)->computeLikelihood();
			if (nei1_part && nei2_part) {
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;

				part_info[part].cur_score = at(part)->computeLikelihoodDerv(nei2_part,(PhyloNode*)nei1_part->node, df_aux, ddf_aux);
				tree_lh += part_info[part].cur_score;
				df -= part_info[part].part_rate*df_aux;
				ddf -= part_info[part].part_rate*part_info[part].part_rate*ddf_aux;
			} else {
				tree_lh += part_info[part].cur_score;
			}
		}
    return -tree_lh;
}

NNIMove PhyloSuperTreePlen::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, bool approx_nni, bool useLS, double lh_contribution)
{
    NNIMove myMove;
    myMove.loglh = 0;

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

	double bestScore = optimizeOneBranch(node1, node2, false);
	double nonNNIScore = bestScore;

	int ntrees = size(), part;

/*	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: nni1_score, nni2_score) private(part)
	#endif
*/
	SwapNNIParam nni_param;
	// nni_param.node1/2_nei tell swapNNIBranch what to swap first
	nni_param.node1_nei = node1_nei;
	nni_param.node2_nei = node2_nei;

	this->swapNNIBranch(0.0, (PhyloNode*)nei2->node, (PhyloNode*)nei1->node, &nni_param);

	// Choose NNI move for SuperTree===========================================
	if (nni_param.nni1_score > bestScore+TOL_LIKELIHOOD) {
		myMove.delta = nni_param.nni1_score - nonNNIScore;
		bestScore = nni_param.nni1_score;
		myMove.swap_id = 1;
		myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei->node);
		myMove.loglh = bestScore;
		myMove.node1 = node1;
		myMove.node2 = node2;
	}

	if (nni_param.nni2_score > bestScore+TOL_LIKELIHOOD) {
		myMove.delta = nni_param.nni2_score - nonNNIScore;
		bestScore = nni_param.nni2_score;
		myMove.swap_id = 2;
		myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei_other->node);
		myMove.loglh = bestScore;
		myMove.node1 = node1;
		myMove.node2 = node2;
	}
	// ========================================================================

	if (save_all_trees != 2) return myMove;

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
    		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];

    		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
    		if (!nei1_part || !nei2_part) {
    			memcpy(at(part)->_pattern_lh, part_info[part].cur_ptnlh, at(part)->getAlnNPattern() * sizeof(double));
    		} else {
				int brid = nei1_part->id;
				bool is_nni = true;
				FOR_NEIGHBOR_DECLARE(node1, node2, nit) {
					if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
				}
				FOR_NEIGHBOR(node2, node1, nit) {
					if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
				}
				if (!is_nni)
					memcpy(at(part)->_pattern_lh, part_info[part].opt_ptnlh[brid], at(part)->getAlnNPattern() * sizeof(double));
				else if (nnino == 0)
					memcpy(at(part)->_pattern_lh, part_info[part].nni1_ptnlh[brid], at(part)->getAlnNPattern() * sizeof(double));
				else
					memcpy(at(part)->_pattern_lh, part_info[part].nni2_ptnlh[brid], at(part)->getAlnNPattern() * sizeof(double));
    		}
    		save_lh_factor[part] = at(part)->current_it->lh_scale_factor;
    		save_lh_factor_back[part] = at(part)->current_it_back->lh_scale_factor;
    		at(part)->current_it->lh_scale_factor = 0.0;
    		at(part)->current_it_back->lh_scale_factor = 0.0;
        }
        if (nnino == 0)
        	saveCurrentTree(nni_param.nni1_score);
        else
        	saveCurrentTree(nni_param.nni2_score);

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

void PhyloSuperTreePlen::doNNI(NNIMove &move)
{
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	SuperNeighbor *node1_nei = (SuperNeighbor*)*move.node1Nei_it;
	SuperNeighbor *node2_nei = (SuperNeighbor*)*move.node2Nei_it;
	int part = 0;
	iterator it;
	double old_brlen = nei1->length;
	PhyloTree::doNNI(move);

	for (it = begin(), part = 0; it != end(); it++, part++) {
		bool is_nni = true;
		FOR_NEIGHBOR_DECLARE(move.node1, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		FOR_NEIGHBOR(move.node2, NULL, nit) {
			if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni = false; break; }
		}
		if (!is_nni) {
			// In the subtree change the length of branch, (node1,node2) WAS linked to
			PhyloNeighbor* nei1_part_a = nei1->link_neighbors[part];
			PhyloNeighbor* nei2_part_a = nei2->link_neighbors[part];
			nei1_part_a->length -= old_brlen * part_info[part].part_rate;
			nei2_part_a->length -= old_brlen * part_info[part].part_rate;

			// Relink the branch if it does not correspond to NNI for partition
			linkBranch(part, nei1, nei2);

			continue;
		}

		NNIMove part_move;
		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
		int brid = nei1_part->id;
		part_move.node1 = (PhyloNode*)nei2_part->node;
		part_move.node2 = (PhyloNode*)nei1_part->node;
		part_move.node1Nei_it = part_move.node1->findNeighborIt(node1_nei->link_neighbors[part]->node);
		part_move.node2Nei_it = part_move.node2->findNeighborIt(node2_nei->link_neighbors[part]->node);

		if (move.swap_id == 1)
			nei1_part->length = nei2_part->length = part_info[part].nni1_brlen[brid];
		else
			nei1_part->length = nei2_part->length = part_info[part].nni2_brlen[brid];

		(*it)->doNNI(part_move);
	}

}

double PhyloSuperTreePlen::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, SwapNNIParam *nni_param) {


	assert(node1->degree() == 3 && node2->degree() == 3);

	//===========================================================================================
	// Prepare the details for NNI: SuperTree
	//===========================================================================================

	//-------------------------------------------------------------------------------------------
	// For restoring purposes -------------------------------------------------------------------

	int i = 0, id = 0;
	int IT_NUM = (params->nni5Branches) ? 6 : 2;
	NeighborVec::iterator it, saved_it[6], node_nei_it[4];
	saved_it[id++] = node1->findNeighborIt(node2);
	saved_it[id++] = node2->findNeighborIt(node1);

	if (params->nni5Branches) {
		FOR_NEIGHBOR(node1, node2, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node1);
			node_nei_it[i++] = it;
		}
		FOR_NEIGHBOR(node2, node1, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node2);
			node_nei_it[i++] = it;
		}
	}
	assert(id == IT_NUM);

	SuperNeighbor *saved_nei[6];
	// save Neighbor and allocate new Neighbor pointer
	for (id = 0; id < IT_NUM; id++) {
		saved_nei[id] = (SuperNeighbor*)(*saved_it[id]);
		*saved_it[id] = new SuperNeighbor(saved_nei[id]->node, saved_nei[id]->length);
	}

	// Getting NEW Neighbors: get the Neighbors again since they were saved for restoring purpose and replaced by new
	SuperNeighbor *nei1_new = (SuperNeighbor*) node1->findNeighbor(node2); // node2|(1-2)
	SuperNeighbor *nei2_new = (SuperNeighbor*) node2->findNeighbor(node1); // node1|(1-2)

	//-------------------------------------------------------------------------------------------
	// NNI details: assigning nodes to be swapped -----------------------------------------------

	// node1_nei - one of the node1 neighbors, which is not node2
	FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)
       break;
	if (nni_param)
		node1_it = node1->findNeighborIt(nni_param->node1_nei->node);
	Neighbor *node1_nei = *node1_it;

	// *node2_its[0] - one of the node2 neighbors, which is not node1
	// *node2_its[1] - second neighbor of node2,   which is not node1
	vector<NeighborVec::iterator> node2_its;
	FOR_NEIGHBOR_DECLARE(node2, node1, node2_it)
		node2_its.push_back(node2_it);
	assert(node2_its.size() == 2);
	if (nni_param && nni_param->node2_nei != (*node2_its[0])) {
		node2_it = node2_its[0];
		node2_its[0] = node2_its[1];
		node2_its[1] = node2_it;
	}

	//===========================================================================================
	// Prepare the details for NNI: SubTrees
	//===========================================================================================

	int part, ntrees = size();
	double old_brlen = saved_nei[0]->length; // length of the branch between node1 and node2 on SuperTree before NNI
	vector<bool> is_nni;
	IntVector brid;
	brid.resize(ntrees);
	is_nni.resize(ntrees);
	/*------------------------------------------------------------------------------------------*
	 * Synchronization:																			*
	 * node1 -> node1_link[part]																*
	 * node2 -> node2_link[part]																*
	 * node1_nei -> node1_link_nei[part]														*
	 * node2_nei -> node2_link_nei[part] - synchronized later, since differs for nni1 and nni2	*
	 *------------------------------------------------------------------------------------------*/
	//PhyloNode *node1_link[ntrees], 	   *node2_link[ntrees];
	//PhyloNeighbor *nei1_new_part[ntrees],  *nei2_new_part[ntrees];
	//PhyloNeighbor *node1_link_nei[ntrees], *node2_link_nei[ntrees];
	vector<PhyloNode*> node1_link,node2_link;
	vector<PhyloNeighbor*> nei1_new_part, nei2_new_part, node1_link_nei, node2_link_nei;
	vector<NeighborVec::iterator> node1_link_it, node2_link_it, sub_saved_it;

	nei1_new_part.resize(ntrees);
	nei2_new_part.resize(ntrees);
	node1_link.resize(ntrees);
	node2_link.resize(ntrees);
	node1_link_nei.resize(ntrees);
	node2_link_nei.resize(ntrees);
	node1_link_it.resize(ntrees);
	node2_link_it.resize(ntrees);
	sub_saved_it.resize(6*ntrees);

	/* to know which nodes to restore:
	 * 		3 spots for node1 as a neighbor of [node2, node1_nei1, node1_nei2] and
	 * 		3 spots for node2 as a neighbor of [node1, node2_nei1, node2_nei2]
	 */


	/* Saved original neighbors:
	 * saved_nei[0] - node2 as a neighbor of node1
	 * saved_nei[1] - node1 as a neighbor of node2
	 * IF(nni5Branches)
	 * 		saved_nei[2(3)] - node1 as a neighbor of its nei1(nei2) different from node2
	 * 		saved_nei[4(5)] - node2 as a neighbor of its nei1(nei2) different from node1
	 */

	Node *node_link, *nei_link;
	for(part = 0; part < ntrees; part++){
		is_nni[part] = true;

//--------------------------------------------------------------------------------------------------------------------
		// Filling out the vector of link_neighbors of new SuperNeighbors on SuperTree with the new PhyloNeighbors on SubTrees
		for(id = 0; id < 2; id++){
			if(saved_nei[id]->link_neighbors[part]){
				nei_link  = saved_nei[id]->link_neighbors[part]->node;
				node_link = saved_nei[1-id]->link_neighbors[part]->node;
				sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
				*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, saved_nei[id]->link_neighbors[part]->length);
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();

				((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back((PhyloNeighbor*)*sub_saved_it[part*6 + id]);

			} else {
				((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
			}
		}
		if(params->nni5Branches){
			for(id = 2; id < 6; id ++){
				if(saved_nei[id]->link_neighbors[part]){
					nei_link = saved_nei[id]->link_neighbors[part]->node;
					node_link = ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node;
					sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
					*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, saved_nei[id]->link_neighbors[part]->length);
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();

					((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back((PhyloNeighbor*)*sub_saved_it[part*6 + id]);
				} else {
					((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
				}
			}
		}
//--------------------------------------------------------------------------------------------------------------------
		if(nei1_new->link_neighbors[part]){
			FOR_NEIGHBOR_DECLARE(node1, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}
			FOR_NEIGHBOR(node2, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}

			nei1_new_part[part] = nei1_new->link_neighbors[part];
			nei2_new_part[part] = nei2_new->link_neighbors[part];

			// partial_lh clear partial likelihood vector
			//nei1_new_part[part]->clearPartialLh();
			//nei2_new_part[part]->clearPartialLh();
		}else{
			is_nni[part]=false;
		}

		if(is_nni[part]){
			// Synchronize node1 and node2 with nodes on SubTree
			node1_link[part] = (PhyloNode*) nei2_new->link_neighbors[part]->node;
			node2_link[part] = (PhyloNode*) nei1_new->link_neighbors[part]->node;

			// Synchronize node1_nei with node on SubTree
			node1_link_nei[part] = ((SuperNeighbor*)node1_nei)->link_neighbors[part];
			node1_link_it[part] = node1_link[part]->findNeighborIt(node1_link_nei[part]->node);
		}
		if(verbose_mode == VB_MAX){
			cout<<"Partition "<<part<<endl;
			at(part)->printTree(cout);
			cout<<endl;
		}

	}

	//===========================================================================================
	// Do the NNI swap and compute the likelihood of swapped topology
	//===========================================================================================
	int cnt;
	for (cnt = 0; cnt < node2_its.size(); cnt++) {
		node2_it = node2_its[cnt];

		// Do the NNI swap on SuperTrees ----------------------------------------------------
		Neighbor *node2_nei = *node2_it;

		node1->updateNeighbor(node1_it, node2_nei);
		node2_nei->node->updateNeighbor(node2, node1);

		node2->updateNeighbor(node2_it, node1_nei);
		node1_nei->node->updateNeighbor(node1, node2);

		// Do the NNI swap on SubTrees or relink (node1,node2) ------------------------------
		for(part = 0; part < ntrees; part++){
			if(!nei1_new->link_neighbors[part]){
				if(part_info[part].cur_score == 0.0){
					part_info[part].cur_score = at(part)->computeLikelihood();
					at(part)->computePatternLikelihood(part_info[part].cur_ptnlh, &part_info[part].cur_score);
				}
			} else {
				brid[part] = nei1_new->link_neighbors[part]->id;
				if(is_nni[part]){
					/*------------------------------------------------------------------*
				 	 * DO NNI by swapping the following neighbors:						*
				 	 * node1_link_nei <-> node2_link_nei								*
				 	 *------------------------------------------------------------------*/
					node2_link_nei[part] = ((SuperNeighbor*)node2_nei)->link_neighbors[part];
					node2_link_it[part] = node2_link[part]->findNeighborIt(node2_link_nei[part]->node);

					// NNI swap
					node1_link[part]->updateNeighbor(node1_link_it[part], node2_link_nei[part]);
					node2_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);

					node2_link[part]->updateNeighbor(node2_link_it[part], node1_link_nei[part]);
					node1_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				}else{
					/*--------------------------------------------------------------*
				 	 * Relink														*
				 	 * the branch (node1, node2) and 								*
				 	 * change the length of old -= and new += branches accordingly	*
				 	 *--------------------------------------------------------------*/

					// Change the length of branch, (node1,node2) WAS linked to (-=)
					nei1_new_part[part]->length -= old_brlen * part_info[part].part_rate;
					nei2_new_part[part]->length -= old_brlen * part_info[part].part_rate;

					if(nei1_new_part[part]->length == 0.0){
						nei1_new_part[part]->length = MIN_BRANCH_LEN;
						nei2_new_part[part]->length = MIN_BRANCH_LEN;
					}
					// Relink the branch if it does not correspond to NNI for partition
					linkBranch(part, nei1_new, nei2_new);

					// Change the length of branch, (node1,node2) is NOW linked to (+=)
					nei1_new_part[part] = nei1_new->link_neighbors[part];
					nei2_new_part[part] = nei2_new->link_neighbors[part];
					nei1_new_part[part]->length += old_brlen * part_info[part].part_rate;
					nei2_new_part[part]->length += old_brlen * part_info[part].part_rate;
				}

				// partial_lh clear partial likelihood vector
				nei1_new_part[part]->clearPartialLh();
				nei2_new_part[part]->clearPartialLh();
			}
		}

		// Compute the score of the swapped topology
		double score = optimizeOneBranch(node1, node2, false);

	    for(part = 0; part < ntrees; part++){
	    	if(nei1_new->link_neighbors[part] && !is_nni[part]){
	    		part_info[part].opt_brlen[brid[part]] = nei1_new_part[part]->length;
	    		at(part)->computePatternLikelihood(part_info[part].opt_ptnlh[brid[part]], &part_info[part].opt_score[brid[part]]);
	    	}
	    }

	    if (params->nni5Branches) {
	    	if (verbose_mode >= VB_DEBUG)
	    		cout << "Log-likelihood: " << score << endl;

	    	FOR_NEIGHBOR(node1, node2, it){
	    		for(part = 0; part < ntrees; part++){
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei2_new->link_neighbors[part]->node;
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		}
	    		score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false);
	    		if (verbose_mode >= VB_DEBUG)
	    			cout << "Log-likelihood: " << score << endl;
	    	}

	    	for(part = 0; part > ntrees; part++){
	    		if(nei2_new->link_neighbors[part])
	    			nei2_new_part[part]->clearPartialLh();
	    	}

	    	FOR_NEIGHBOR(node2, node1, it)
	    	{
	    		for(part = 0; part < ntrees; part++){
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei1_new->link_neighbors[part]->node;
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		}
	    		score = optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false);
	    		if (verbose_mode >= VB_DEBUG)
	    			cout << "Log-likelihood: " << score << endl;
	    	}

	    }

	 /*   for(part = 0; part < ntrees; part++){
	    	if(nei1_new->link_neighbors[part] && is_nni[part]){
		    	if(cnt == 0){
		    		computePatternLikelihood(part_info[part].nni1_ptnlh[brid[part]], &part_info[part].nni1_score[brid[part]]);
		    	}else{
		    		computePatternLikelihood(part_info[part].nni2_ptnlh[brid[part]], &part_info[part].nni2_score[brid[part]]);
		    	}
	    	}
	    }*/

		if (nni_param) {
			if (verbose_mode >= VB_MAX)
				printTree(cout, WT_BR_LEN + WT_NEWLINE);
			if (cnt == 0) {
				nni_param->nni1_score = score;
				nni_param->nni1_brlen = nei1_new->length;
			} else {
				nni_param->nni2_score = score;
				nni_param->nni2_brlen = nei1_new->length;
			}
		}

	    // Swap back on SuperTree -------------------------------------------------------
		node1->updateNeighbor(node1_it, node1_nei);
		node1_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node2_nei);
		node2_nei->node->updateNeighbor(node1, node2);

		// Swap back or relink back on SubTrees -----------------------------------------
		for(part = 0; part < ntrees; part++){
			if(nei1_new->link_neighbors[part]){
				if(is_nni[part]){
					/*------------------------------------------------------------------*
					 * Swap BACK by swapping the following neighbors:					*
					 * node1_link_nei <-> node2_link_nei								*
					 *------------------------------------------------------------------*/
					node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
					node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);

					node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
					node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);
				}else{
					/*--------------------------------------------------------------*
					 * Relink BACK													*
					 * the branch (node1, node2) and 								*
					 * change the length of new -= and old += branches accordingly	*
					 *--------------------------------------------------------------*/

					// Change the length of branch, (node1,node2) WAS relinked to (-=)
					nei1_new_part[part]->length -= old_brlen * part_info[part].part_rate;
					nei2_new_part[part]->length -= old_brlen * part_info[part].part_rate;

					// Relink BACK the branch if it does not correspond to NNI for partition
					linkBranch(part, nei1_new, nei2_new);

					// Change the length of branch, (node1,node2) WAS linked to before (+=)
					nei1_new_part[part] = nei1_new->link_neighbors[part];
					nei2_new_part[part] = nei2_new->link_neighbors[part];
					nei1_new_part[part]->length += old_brlen * part_info[part].part_rate;
					nei2_new_part[part]->length += old_brlen * part_info[part].part_rate;
				}
			}
		}


	}

	//===========================================================================================
	// Restoring
	//===========================================================================================

	// Restoring information for SuperTree ------------------------------------------------------
	// restore the Neighbors*
	for (id = IT_NUM-1; id >= 0; id--) {
		if (*saved_it[id] == current_it) current_it = (SuperNeighbor*) saved_nei[id];
		if (*saved_it[id] == current_it_back) current_it_back = (SuperNeighbor*) saved_nei[id];

		delete (*saved_it[id]);
		(*saved_it[id]) = saved_nei[id];
	 }

	// restore the length of 4 branches around node1, node2
	FOR_NEIGHBOR(node1, node2, it)
		(*it)->length = (*it)->node->findNeighbor(node1)->length;
	FOR_NEIGHBOR(node2, node1, it)
		(*it)->length = (*it)->node->findNeighbor(node2)->length;

	// Restoring information for SubTrees ------------------------------------------------------
	for(part = 0; part < ntrees; part++){
		// restore the Neighbors*
		for (i = IT_NUM-1; i >= 0; i--) {
			if((*sub_saved_it[part*6+i])){
				delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->scale_num;
				delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->partial_lh;
				if (*sub_saved_it[part*6+i] == at(part)->current_it) at(part)->current_it = saved_nei[i]->link_neighbors[part];
				if (*sub_saved_it[part*6+i] == at(part)->current_it_back) at(part)->current_it_back = saved_nei[i]->link_neighbors[part];

				delete (*sub_saved_it[part*6+i]);
				(*sub_saved_it[part*6+i]) = saved_nei[i]->link_neighbors[part];
			}
		}
		// restore the length of 4 branches around node1_link[part], node2_link[part]
		node1_link[part] = (PhyloNode*) ((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part]->node;
		node2_link[part] = (PhyloNode*) ((SuperNeighbor*)node1->findNeighbor(node2))->link_neighbors[part]->node;
		FOR_NEIGHBOR(node1_link[part], node2_link[part], it)
			(*it)->length = (*it)->node->findNeighbor(node1_link[part])->length;
		FOR_NEIGHBOR(node2_link[part], node1_link[part], it)
			(*it)->length = (*it)->node->findNeighbor(node2_link[part])->length;
		}
	return cur_score;
}

void PhyloSuperTreePlen::computeBranchLengths()
{
	}
