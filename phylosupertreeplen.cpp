/*
 * phylosupertreeplen.cpp
 *
 *  Created on: Aug 5, 2013
 *      Author: olga
 */

#include "phylosupertreeplen.h"
#include "superalignmentpairwise.h"
#include <string.h>
#include "timeutil.h"

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

    //tree->initPartitionInfo(); // FOR OLGA: needed here

    //tree_lh = tree->computeLikelihood();
	if (fixed_len) {
    	for(int part = 0; part < ntrees; part++){
    		tree->part_info[part].cur_score = 0.0;
    	}
		tree_lh = tree->computeLikelihood();
	} else {
		tree_lh = tree->optimizeAllBranches(1);
	}

    cout<<"Initial log-likelihood: "<<tree_lh<<endl;
	double begin_time = getCPUTime();
	int i;
    for(i = 1; i < 100; i++){
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
    	int my_iter = min(5,i+1);

    	//cout<<"BEFORE BRANCH OPTIMIZATION"<<endl;
        //tree->printTree(cout);
    	if(!fixed_len){
    		cur_lh = tree->optimizeAllBranches(my_iter,tol);
    	}
    	cout<<"Current log-likelihood at step "<<i<<": "<<cur_lh<<endl;
    	if(fabs(cur_lh-tree_lh) < epsilon)
    		break;

    	tree_lh = cur_lh;
    }
    //tree->printTree(cout);
	//cout<<endl;
    cout <<"OPTIMIZE MODEL has finished"<< endl;
    for(int part = 0; part < ntrees; part++){
    	cout<<"PART RATE "<<part<<" = "<<tree->part_info[part].part_rate<<endl;
    	//tree->at(part)->printTree(cout);
    	//cout<<endl;
    }
	cout << "Parameters optimization took " << i-1 << " rounds (" << getCPUTime()-begin_time << " sec)" << endl << endl;

    //exit(2);
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
			//tree->at(part)->scaleLength(rate/tree->part_info[part].part_rate);
			tree->part_info[part].part_rate = rate;
		}
	}
	tree->mapBranchLen();
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
    if(!tree->fixed_rates)
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
	memset(allNNIcases_computed, 0, 5*sizeof(int));
	fixed_rates = false;
}

PhyloSuperTreePlen::PhyloSuperTreePlen(Params &params)
: PhyloSuperTree(params)
{
	memset(allNNIcases_computed, 0, 5*sizeof(int));
	fixed_rates = (params.partition_type == 'j') ? true : false;
	int part = 0;
	for (iterator it = begin(); it != end(); it++, part++) {
		part_info[part].part_rate = 1.0;
		part_info[part].evalNNIs = 0.0;
	}
}

PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment, PhyloSuperTree *super_tree)
: PhyloSuperTree(alignment,super_tree)
{
	memset(allNNIcases_computed, 0, 5*sizeof(int));
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

double PhyloSuperTreePlen::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
	//initPartitionInfo(); // OLGA: not needed here
	//cout<<"Optimizing all branches"<<endl;
	for(int part = 0; part < size(); part++){
		part_info[part].cur_score = 0.0;
	}
	return PhyloTree::optimizeAllBranches(my_iterations,tolerance, maxNRStep);
}

double PhyloSuperTreePlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH, int maxNRStep) {

	SuperNeighbor *nei1 = (SuperNeighbor*)node1->findNeighbor(node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)node2->findNeighbor(node1);

	for (int part = 0; part < size(); part++) {
		at(part)->theta_computed = false;
	}

	double tree_lh = PhyloTree::optimizeOneBranch(node1,node2,clearLH, maxNRStep);



	if(clearLH){
		for (int part = 0; part < size(); part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			if(nei1_part){
				((PhyloNode*)nei1_part->node)->clearReversePartialLh(((PhyloNode*)nei2_part->node));
				((PhyloNode*)nei2_part->node)->clearReversePartialLh(((PhyloNode*)nei1_part->node));
			}
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
			if (nei1_part && nei2_part) {
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei2_part,(PhyloNode*)nei1_part->node);
				tree_lh += part_info[part].cur_score;
			} else {
				if (part_info[part].cur_score == 0.0)
					part_info[part].cur_score = at(part)->computeLikelihood();
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
			if (nei1_part && nei2_part) {
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;
				if(nei1_part->length<-1e-4){
					cout<<"lambda = "<<lambda<<endl;
					cout<<"NEGATIVE BRANCH len = "<<nei1_part->length<<endl<<" rate = "<<part_info[part].part_rate<<endl;
					outError("shit!!   ",__func__);
				}
				part_info[part].cur_score = at(part)->computeLikelihoodDerv(nei2_part,(PhyloNode*)nei1_part->node, df_aux, ddf_aux);
				tree_lh += part_info[part].cur_score;
				df -= part_info[part].part_rate*df_aux;
				ddf -= part_info[part].part_rate*part_info[part].part_rate*ddf_aux;
			} else {
				if (part_info[part].cur_score == 0.0)
					part_info[part].cur_score = at(part)->computeLikelihood();
				tree_lh += part_info[part].cur_score;
			}
		}
    return -tree_lh;
}

NNIMove PhyloSuperTreePlen::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves, bool approx_nni, bool useLS, double lh_contribution)
{
    NNIMove myMove;
    myMove.newloglh = 0;

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
	//double oldLEN = node1->findNeighbor(node2)->length;

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
	if (nni_param.nni1_score > nni_param.nni2_score) {
		myMove.swap_id = 1;
		myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei->node);
		myMove.newloglh = nni_param.nni1_score;
		myMove.node1 = node1;
		myMove.node2 = node2;
		myMove.newLen[0] = nni_param.nni1_brlen;
		//myMove.oldLen[0] = oldLEN;
	} else {
		myMove.swap_id = 2;
		myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei_other->node);
		myMove.newloglh = nni_param.nni2_score;
		myMove.node1 = node1;
		myMove.node2 = node2;
		myMove.newLen[0] = nni_param.nni2_brlen;
		//myMove.oldLen[0] = oldLEN;
	}
	// ========================================================================
	return myMove;
}

void PhyloSuperTreePlen::applyNNIs(int nni2apply, bool changeBran) {
	IQTree::applyNNIs(nni2apply, changeBran);
	mapBranchLen();
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
	double old_brlen = nei1->length;
	vector<NNIMove> part_move;
	vector<bool> is_nni;
	is_nni.resize(ntrees);
	part_move.resize(ntrees);

	for (it = begin(), part = 0; it != end(); it++, part++) {
		PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
		PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
		is_nni[part] = true;
		if(nei1_part){
			FOR_NEIGHBOR_DECLARE(move.node1, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}
			FOR_NEIGHBOR(move.node2, NULL, nit) {
				if (! ((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}
		} else { is_nni[part] = false;}

		if(is_nni[part]){
			part_move[part].node1 = (PhyloNode*)nei2_part->node;
			part_move[part].node2 = (PhyloNode*)nei1_part->node;

			part_move[part].node1Nei_it = part_move[part].node1->findNeighborIt(node1_nei->link_neighbors[part]->node);
			part_move[part].node2Nei_it = part_move[part].node2->findNeighborIt(node2_nei->link_neighbors[part]->node);
		}else{
			// In the subtree change the length of branch, (node1,node2) WAS linked to
			if(nei1_part){
				nei1_part->length -= old_brlen * part_info[part].part_rate;
				nei2_part->length -= old_brlen * part_info[part].part_rate;
			}
		}
	}
	PhyloTree::doNNI(move,clearLH);
	nei1->length = move.newLen[0];
	nei2->length = move.newLen[0];

	for (it = begin(), part = 0; it != end(); it++, part++) {
		if (!is_nni[part]) {
			// Relink the branch if it does not correspond to NNI for partition
			linkBranch(part, nei1, nei2);
		} else { (*it)->doNNI(part_move[part],clearLH); }

//		PhyloNeighbor* nei1_part_new = nei1->link_neighbors[part];
//		PhyloNeighbor* nei2_part_new = nei2->link_neighbors[part];
//		if(nei1_part_new){
//			int brid = nei1_part_new->id;
//			if (move.swap_id == 1)
//				nei1_part_new->length = nei2_part_new->length = part_info[part].nni1_brlen[brid];
//			else
//				nei1_part_new->length = nei2_part_new->length = part_info[part].nni2_brlen[brid];
//		}
	}
	//mapBranchLen();

//	printTree(cout);
//	cout<<endl;
//checkBranchLen();

}

/**
 * this is to classify the cases which happen on the subtree
 *
 *  NNI_NONE_EPSILON: all 5 branches have images on subtree, this corresponds to change in subtree topology
 * 					  2 partial_lh vectors for -nni1 or 6 partial_lh vectors for -nni5 options
 *  NNI_ONE_EPSILON:  only one of the 5 branches has no image on subtree, this does not change subtree topology, but changes branch length of subtrees
 * 					  we need to allocate partial likelihood memory (1 partial_lh vectors for -nni1 option or 3 partial_lh for -nni5 option)
 * 	NNI_TWO_EPSILON:  two branches (on different sides of central branch) have no images, here after the NNI swap,
 * 					  the image of central branch either does not change or is equal to epsilon (then we decrease the branch length)
 * 					  and no allocation of partial_lh is needed
 * 	NNI_THREE_EPSILON: central and two adjacent edges have no images: after the NNI swap, central branch will have image and we need to relink it
 * 					no allocation of partial_lh is needed
 *  NNI_MANY_EPSILON: more than 3 branches have no images on subtree: nothing changes in subtree and no recomputation of partial likelihood are required
 */
enum NNIType {NNI_NO_EPSILON, NNI_ONE_EPSILON, NNI_TWO_EPSILON, NNI_THREE_EPSILON, NNI_MANY_EPSILON};


double PhyloSuperTreePlen::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, SwapNNIParam *nni_param) {

	assert(node1->degree() == 3 && node2->degree() == 3);
	checkBranchLen();

	int i = 0, id = 0;
	int part, ntrees = size();

	/*===========================================================================================
	 * Identify NNIType for partitions
	 *===========================================================================================*/
	int epsilon_cnt;
	vector<NNIType> is_nni;
	is_nni.resize(ntrees, NNI_NO_EPSILON);
	for(part=0; part<ntrees;part++){
		totalNNIs++;
		is_nni[part] = NNI_NO_EPSILON;
		epsilon_cnt = 0;

		FOR_NEIGHBOR_DECLARE(node1,NULL,nit){
			if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { epsilon_cnt++; }
		}
		FOR_NEIGHBOR(node2, node1, nit) {
			if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { epsilon_cnt++; }
		}

		//cout<<"Partition "<<part<<" : Epsilon = "<<epsilon_cnt<<endl;
		if(epsilon_cnt == 0){
			is_nni[part]=NNI_NO_EPSILON;
		}else if(epsilon_cnt == 1){
			is_nni[part] = NNI_ONE_EPSILON;
		}else if(epsilon_cnt == 2){
			is_nni[part]=NNI_TWO_EPSILON;
		}else if(epsilon_cnt == 3){
			is_nni[part]=NNI_THREE_EPSILON;
		}else {
			is_nni[part] = NNI_MANY_EPSILON;
		}
	}

	//===========================================================================================
	// Prepare the details for NNI: SuperTree
	//===========================================================================================
	double old_brlen = node1->findNeighbor(node2)->length; // length of the branch between node1 and node2 on SuperTree before NNI
	int IT_NUM = (params->nni5) ? 6 : 2;
	NeighborVec::iterator it, saved_it[6], node_nei_it[4];
	saved_it[id++] = node1->findNeighborIt(node2);
	saved_it[id++] = node2->findNeighborIt(node1);

	//if (params->nni5) {
		FOR_NEIGHBOR(node1, node2, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node1);
			node_nei_it[i++] = it;
		}
		FOR_NEIGHBOR(node2, node1, it){
			saved_it[id++] = (*it)->node->findNeighborIt(node2);
			node_nei_it[i++] = it;
		}
	//}
	assert(id == IT_NUM);

	// Save SuperNeighbors
	SuperNeighbor *saved_nei[6];

	// allocate new Super Neighbor pointers
	for (id = 0; id < IT_NUM; id++) {
		saved_nei[id] = (SuperNeighbor*)(*saved_it[id]);
		*saved_it[id] = new SuperNeighbor(saved_nei[id]->node, saved_nei[id]->length);
		(*saved_it[id])->id = saved_nei[id]->id;
		for(part = 0; part < ntrees; part++)
			((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
	}
	/*---------------------------------------------------------
	 * Save branch lengths on SubTrees for restoring purposes
	 *---------------------------------------------------------*/
	/*
	 * NO_EPS:  one/five branches need to be restored in nni1/nni5 cases respectively,
	 * 			but these branches are saved in saved_nei->link_neighbors, we won't store them again
	 *
	 * ONE_EPS: three branches need to be restored (stick to ids: 0,...,5)
	 * 			in nni1 case only the saved_nei for central branch is saved, keep this in mind, when restoring
	 * TWO_EPS: the image of central branch needs to be restored (the id for restoring [6*part+0])
	 * THREE_EPS: one branch needs to be restored: which the central is relinked to after NNI (the id for restoring [6*part+0])
	 *
	 * MANY_EPS: nothing to be restored
	 */
	vector<double*> sub_saved_branch;
	sub_saved_branch.resize(6*ntrees);

	/* Minh bug fix: save the cur_score */
	double *saved_cur_score = new double[part_info.size()];
	for (i = 0; i < part_info.size(); i++)
		saved_cur_score[i] = part_info[i].cur_score;

	/*------------------------------------------------------------------------------------
	 * Saved original neighbors:
	 * saved_nei[0] - node2 as a neighbor of node1
	 * saved_nei[1] - node1 as a neighbor of node2
	 * IF(nni5Branches)
	 * 		saved_nei[2(3)] - node1 as a neighbor of its nei1(nei2) different from node2
	 * 		saved_nei[4(5)] - node2 as a neighbor of its nei1(nei2) different from node1
	 *------------------------------------------------------------------------------------*/

	// Getting NEW Neighbors: get the Neighbors again since they were saved for restoring purpose and replaced by new
	SuperNeighbor *nei1_new = (SuperNeighbor*) node1->findNeighbor(node2);
	SuperNeighbor *nei2_new = (SuperNeighbor*) node2->findNeighbor(node1);

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
	/*------------------------------------------------------------------------------------------*
	 * Synchronization:																			*
	 * node1 -> node1_link[part]																*
	 * node2 -> node2_link[part]																*
	 * node1_nei -> node1_link_nei[part]														*
	 * node2_nei -> node2_link_nei[part] - synchronized later, since differs for nni1 and nni2	*
	 *------------------------------------------------------------------------------------------*/

	Node *node_link, *nei_link;
	vector<PhyloNode*> node1_link,node2_link;
	vector<PhyloNeighbor*> node1_link_nei,node2_link_nei,sub_saved_nei1,sub_saved_nei2,sub_saved_nei;
	vector<NeighborVec::iterator> node1_link_it, node2_link_it, sub_saved_it;

	// Nodes which correspond to node1 and node2 on partitions
	node1_link.resize(ntrees);
	node2_link.resize(ntrees);

	// Neighbors of node1_link and node2_link to be swapped during NNI
	node1_link_nei.resize(ntrees);
	node2_link_nei.resize(ntrees);

	// iterators for the neighbors of node1_link and node2_link to be swapped
	node1_link_it.resize(ntrees);
	node2_link_it.resize(ntrees);


	/*------------------------------------------------------------------------------------
	 * Saving neighbors and their iterators on partitions
	 *------------------------------------------------------------------------------------*/

	// Saving linked neighbor of node1->findNei(node2)
	sub_saved_nei1.resize(ntrees);

	// Saving linked neighbor of node2->findNei(node1)
	sub_saved_nei2.resize(ntrees);

	// This vector of Phylo Neighbors will be used only for ONE_EPS case.
	// Here, we'll have to save and restore 1/3 neis for nni1/nni5. For convenience I'm still gonna use 6*ntrees elements vector.
	sub_saved_nei.resize(6*ntrees);
	sub_saved_it.resize(6*ntrees);

	/*------------------------------------------------------------------------------------
	 * At this point we work only with partitions satisfying is_nni condition:
	 * 		- check if is_nni[part]
	 * 			YES:
	 * 			- allocate new PhyloNeighbor
	 * 			- check if some other branch links to the same branch on SubTree
	 * 			- update link_neighbors[part]
	 * 			- synchronize node1_link[part], node2_link[part], node1_link_nei[part]
	 *------------------------------------------------------------------------------------*/

	for(part = 0; part < ntrees; part++){
		if(is_nni[part]==NNI_NO_EPSILON){
			evalNNIs++;
			part_info[part].evalNNIs++;

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
				*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, saved_nei[id]->link_neighbors[part]->length);
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
				(*sub_saved_it[part*6 + id])->id = saved_nei[id]->link_neighbors[part]->id;

				// update link_neighbor[part]: for New SuperNeighbor we set the corresponding new PhyloNeighbor on partition part
				((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
			}

			// optimization on 5 branches ------------------------------------------------------------------
			if(params->nni5){
				for(id = 2; id < 6; id ++){
					nei_link = saved_nei[id]->link_neighbors[part]->node;
					node_link = ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node;
					sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
					*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, saved_nei[id]->link_neighbors[part]->length);
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
					(*sub_saved_it[part*6 + id])->id = saved_nei[id]->link_neighbors[part]->id;

					// update link_neighbor[part]
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
				}
			}
			// Synchronization: node1, node2, node1_nei; node2_nei varies for two NNIs therefore synchronized later
			node1_link[part] = (PhyloNode*) nei2_new->link_neighbors[part]->node;
			node2_link[part] = (PhyloNode*) nei1_new->link_neighbors[part]->node;
			node1_link_nei[part] = ((SuperNeighbor*)node1_nei)->link_neighbors[part];
			node1_link_it[part] = node1_link[part]->findNeighborIt(node1_link_nei[part]->node);

		} else if(is_nni[part]==NNI_ONE_EPSILON){
			// Allocate three new Phylo Neighbors here. For nni1 only one will be used and for nni5 all three.
			// Make sure to update all the necessary link_neighbors and take care of branch lengths
			// (increase/decrease by central branch where necessary).


		}
	}

	//===========================================================================================
	// Do the NNI swap and compute the likelihood of swapped topology
	//===========================================================================================
	int cnt;
	for (cnt = 0; cnt < node2_its.size(); cnt++) {
		node2_it = node2_its[cnt];
		Neighbor *node2_nei = *node2_it;

		// Synchronize node2_link_nei for NNI on SubTree
		for(part=0; part<ntrees; part++)
			if(is_nni[part]==NNI_NO_EPSILON){
				node2_link_nei[part] = ((SuperNeighbor*)node2_nei)->link_neighbors[part];
				node2_link_it[part] = node2_link[part]->findNeighborIt(node2_link_nei[part]->node);
			}

		// Do the NNI swap on SuperTrees ----------------------------------------------------
		node1->updateNeighbor(node1_it, node2_nei);
		node2_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node1_nei);
		node1_nei->node->updateNeighbor(node1, node2);

// Do the NNI swap on SubTrees or Proceed with partitions requiring relink of (node1,node2) ------------------------------
		for(part = 0; part < ntrees; part++){

			if(is_nni[part]==NNI_NO_EPSILON){
				//allNNIcases_computed[0] += 1;
				// Do NNI swap on partition
				node1_link[part]->updateNeighbor(node1_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				for(id=0; id<IT_NUM; id++){
					((PhyloNeighbor*)(*sub_saved_it[part*6+id]))->clearPartialLh();
				}
			} else if(is_nni[part]==NNI_MANY_EPSILON){
				//no need to change or allocate anything
			} else if(is_nni[part]==NNI_THREE_EPSILON){
				linkBranch(part, nei1_new, nei2_new);
				assert(nei1_new->link_neighbors[part]);

				// Save the branch length
				if(cnt == 0)
					sub_saved_branch[6*part] = new double(nei1_new->link_neighbors[part]->length);

				nei1_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;

				// since the branch length was changed we have to recompute the likelihood of the branch
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei1_new->link_neighbors[part],
						(PhyloNode*)nei2_new->link_neighbors[part]->node);

			}else if(is_nni[part]==NNI_TWO_EPSILON){
				/* Olga: bugfix
				 	 	 In fact, before relinking the image of central branch is NULL
				 	 	 after relinking it can be either NULL or it should relink to the same branch as before.
				*/

				// Save the branch length
				if(cnt == 0)
					sub_saved_branch[6*part] = new double(nei1_new->link_neighbors[part]->length);

				linkBranch(part, nei1_new, nei2_new);
				if(!nei1_new->link_neighbors[part]){
					saved_nei[0]->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
					saved_nei[1]->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
					part_info[part].cur_score = at(part)->computeLikelihoodBranch(saved_nei[0]->link_neighbors[part],
							(PhyloNode*)saved_nei[1]->link_neighbors[part]->node);
				}
			}else if(is_nni[part] == NNI_ONE_EPSILON){

				nei1_new->link_neighbors[part] = saved_nei[0]->link_neighbors[part];
				nei2_new->link_neighbors[part] = saved_nei[1]->link_neighbors[part];

				// Change the length of branch, (node1,node2) WAS linked to (-=)
				nei1_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
				assert(nei1_new->link_neighbors[part]->length >= 0.0);

				// Relink the branch if it does not correspond to NNI for partition
				linkBranch(part, nei1_new, nei2_new);

				if(!params->nni5){
					sub_saved_branch[6*part] = new double (nei1_new->link_neighbors[part]->length);
				}

				// Minh sanity check that old and newly relinked branches MUST be adjacent in subtree!
				int common_nodes = 0;
				if (nei1_new->link_neighbors[part] && saved_nei[0]->link_neighbors[part]) {
					if (nei1_new->link_neighbors[part]->node == saved_nei[0]->link_neighbors[part]->node) common_nodes++;
					if (nei1_new->link_neighbors[part]->node == saved_nei[1]->link_neighbors[part]->node) common_nodes++;
					if (nei2_new->link_neighbors[part]->node == saved_nei[0]->link_neighbors[part]->node) common_nodes++;
					if (nei2_new->link_neighbors[part]->node == saved_nei[1]->link_neighbors[part]->node ) common_nodes++;
					if (common_nodes == 0) {
						printMapInfo();
						cout << "common_nodes = " << common_nodes << endl;
						cout << nei1_new->link_neighbors[part]->node->id << " " << nei2_new->link_neighbors[part]->node->id << " "
							 << saved_nei[0]->link_neighbors[part]->node->id << " " << saved_nei[1]->link_neighbors[part]->node->id << endl;
						outError("linkBranch failed in ", __func__);
					}
				}

				assert(nei1_new->link_neighbors[part]);
				//sub_saved_nei1[part] = nei1_new->link_neighbors[part];
				//sub_saved_nei2[part] = nei2_new->link_neighbors[part];

				// Allocate new PhyloNeighbors

				// one branch optimization
				for(id = 0; id < 2; id++){
					nei_link  =  (id == 0) ? nei1_new->link_neighbors[part]->node : nei2_new->link_neighbors[part]->node;
					node_link =  (id == 1) ? nei1_new->link_neighbors[part]->node : nei2_new->link_neighbors[part]->node;
					sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);

					// Saving neighbor for restoring purposes
					sub_saved_nei[6*part + id] = node_link->findNeighbor(nei_link);

					// the central branch length of the subtree is increased!
					*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, nei1_new->link_neighbors[part]->length + old_brlen * part_info[part].part_rate);
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
					(*sub_saved_it[part*6 + id])->id = nei1_new->link_neighbors[part]->id;
					// update link_neighbor[part]
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];

					// For nni5 we need to update also the link_nei for the corresponding nei,
					// where the central branch was relinked to
					if(params->nni5){
						for(int idh = 2; idh < 6; idh++){
							if(node_link == ((SuperNeighbor*)(*node_nei_it[idh-2]))->link_neighbors[part]->node){
								((SuperNeighbor*)*saved_it[idh])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
							}
						}
					}
				}

				// optimization on 5 branches!!!!!
				if(params->nni5){
					for(id = 2; id < 6; id++){
						if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
							if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node
									!= nei1_new->link_neighbors[part]->node
								&& ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node
								!= nei2_new->link_neighbors[part]->node)
							{
								nei_link = saved_nei[id]->link_neighbors[part]->node;
								node_link = ((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node;
								sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);

								// Save neighbor for restoring purposes
								sub_saved_nei[6*part+id] = node_link->findNeighbor(nei_link);

								*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, ((SuperNeighbor*)(*saved_it[id]))->link_neighbors[part]->length);
								((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
								((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
								(*sub_saved_it[part*6 + id])->id = saved_nei[id]->link_neighbors[part]->id;

								// update link_neighbor[part]
								((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
							}
						}
					}
				}// end of nni5


			} // end of else ONE_EPS case
		} // end of part loop

/*===============================================================================================================================*
 * 											Compute the score of the swapped topology 				  							 *
 *===============================================================================================================================*/
		checkBranchLen();
		double score = optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
		if (verbose_mode >= VB_MED) {
			cout << "[" << score << "] ";
			printTree(cout);
			cout << endl;
		}
		checkBranchLen();

		// %%%%%%%%%%%%%%%%%%%%%%%%  FIVE BRANCH OPTIMIZATION  %%%%%%%%%%%%%%%%%%%%%%%%
	    if (params->nni5) {
	    	FOR_NEIGHBOR(node1, node2, it){
	    		for(part = 0; part < ntrees; part++)
	    			//if(((SuperNeighbor*)(*it))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || NNI_ONE_EPSILON)){
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei2_new->link_neighbors[part]->node;
						((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false, NNI_MAX_NR_STEP);
	    	}
	    	for(part = 0; part < ntrees; part++)
	    		//if(((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || NNI_ONE_EPSILON))
	    		if(((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part])
	    			((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part]->clearPartialLh();
	    	FOR_NEIGHBOR(node2, node1, it){
	    		for(part = 0; part < ntrees; part++)
	    			//if(((SuperNeighbor*)(*it))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON || NNI_ONE_EPSILON)){
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei1_new->link_neighbors[part]->node;
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		score = optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false, NNI_MAX_NR_STEP);
	    	}
	    }
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%  END of nni5branch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	    // Save current tree for ufboot analysis
	    if (save_all_trees == 2) {
	    		saveCurrentTree(score);
	    }

	    // *************************** STORE INFO ABOUT NNI ***************************

	    // Store information about this NNI for NNImove for SuperTree
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
		// ***************************************************************************

		// =============================== RESTORE INFO ==============================
		/* Minh bug fix: restore the cur_score */
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

// Swap back or relink back on SubTrees------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){
			if(is_nni[part]==NNI_NO_EPSILON){
				node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);
			} else if(is_nni[part]==NNI_ONE_EPSILON){
				//linkCheckRe(part,node1,node2,sub_saved_nei2[part],sub_saved_nei1[part]);
				//linkCheckRe(part,node2,node1,sub_saved_nei1[part],sub_saved_nei2[part]);
				for (id = IT_NUM-1; id >= 0; id--) {
					if((*sub_saved_it[part*6+id])){
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+id])->scale_num;
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+id])->partial_lh;
						if (*sub_saved_it[part*6+id] == at(part)->current_it)
							at(part)->current_it = saved_nei[id]->link_neighbors[part];
						if (*sub_saved_it[part*6+id] == at(part)->current_it_back)
							at(part)->current_it_back = saved_nei[id]->link_neighbors[part];

						delete (*sub_saved_it[part*6+id]);
						(*sub_saved_it[part*6+id]) = sub_saved_nei[6*part+id];
					}
					// The branch b on subtree to which the central was relinked to must be set to NULL,
					// but we also have to update the linked nei of corresponding branch on parent tree (which is initially mapped to b)
					// since I do not want to check which from the four branches it is, I set everything to NULL
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = NULL;

				}
				// next three lines could be commented out, oder?
				/*linkBranch(part, nei1_new, nei2_new);
				assert(nei1_new->link_neighbors[part]->node == saved_nei[0]->link_neighbors[part]->node);
				assert(nei2_new->link_neighbors[part]->node == saved_nei[1]->link_neighbors[part]->node);
				*/

				// Restore branches
				/*
				 * for neighbors of node1->link, node2->link
				 * (*it)->length = (*it)->findNeighbor(node1->link->node)->length;
				 *
				 * Same for node2->link
				 */


				// Increase the central branch, since the length that was saved, was decreased
				saved_nei[0]->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				saved_nei[1]->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;


			} else if(is_nni[part]==NNI_THREE_EPSILON){
				nei1_new->link_neighbors[part]->length = sub_saved_branch[6*part];
				nei2_new->link_neighbors[part]->length = sub_saved_branch[6*part];
				//linkBranch(part, nei1_new, nei2_new);
			} else if(is_nni[part]==NNI_TWO_EPSILON){
				saved_nei[0]->link_neighbors[part]->length = sub_saved_branch[6*part];
				saved_nei[1]->link_neighbors[part]->length = sub_saved_branch[6*part];
			} else if(is_nni[part]==NNI_MANY_EPSILON){
				// There is no need to restore anything
			}
		}

/*		for(part = 0; part < ntrees; part++){
			if(is_nni[part]){
				node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);
				//Restoring the branch length on the Sub Tree
				node1_link[part]->findNeighbor(node2_link[part])->length = saved_nei[0]->link_neighbors[part]->length;
				node2_link[part]->findNeighbor(node1_link[part])->length = saved_nei[0]->link_neighbors[part]->length;
			}else{
				// Relink back
				linkBranch(part, nei1_new, nei2_new);

				if(nei1_new->link_neighbors[part]){
					nei1_new->link_neighbors[part]->length += old_brlen*part_info[part].part_rate;
					nei2_new->link_neighbors[part]->length += old_brlen*part_info[part].part_rate;
				}
			}
		}*/

		//mapBranchLen();
		//checkBranchLen();
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
			if(is_nni[part] == NNI_NO_EPSILON){ // for other cases everything is restored within for(cnt) loop
				// restore the Neighbors*
				for (i = IT_NUM-1; i >= 0; i--) {
					if((*sub_saved_it[part*6+i])){
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->scale_num;
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->partial_lh;
						if (*sub_saved_it[part*6+i] == at(part)->current_it) at(part)->current_it = saved_nei[i]->link_neighbors[part];
						if (*sub_saved_it[part*6+i] == at(part)->current_it_back) at(part)->current_it_back = saved_nei[i]->link_neighbors[part];

						delete (*sub_saved_it[part*6+i]);
						(*sub_saved_it[part*6+i]) = saved_nei[i]->link_neighbors[part];
//						cout<<"Restored....."<<endl;
//						cout<<"(*sub_saved_it)->node->id"<<(*sub_saved_it[part*6+i])->node->id<<endl;
					}
				}
				// restore the length of 4 branches around node1_link[part], node2_link[part]
				node1_link[part] = (PhyloNode*)(saved_nei[1]->link_neighbors[part]->node);
				node2_link[part] = (PhyloNode*)(saved_nei[0]->link_neighbors[part]->node);
				FOR_NEIGHBOR(node1_link[part], node2_link[part], it)
					(*it)->length = (*it)->node->findNeighbor(node1_link[part])->length;
				FOR_NEIGHBOR(node2_link[part], node1_link[part], it)
					(*it)->length = (*it)->node->findNeighbor(node2_link[part])->length;

				//nei1_new = (SuperNeighbor*) node1->findNeighbor(node2);
				//nei2_new = (SuperNeighbor*) node2->findNeighbor(node1);
			}
		}
		//mapBranchLen();
		checkBranchLen();
//------------------------------------------------------------------------------------------------------------------------------------------------
	// Minh: delete memory
	delete [] saved_cur_score;
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
void PhyloSuperTreePlen::restoreAllBranLen(PhyloNode *node, PhyloNode *dad) {
	IQTree::restoreAllBranLen(node,dad);
	mapTrees();
}

bool PhyloSuperTreePlen::checkBranchLen(){
	/*

	NodeVector nodes1,nodes2;
	int i;
	getBranches(nodes1, nodes2);
	double *checkVAL = new double[branchNum];
	for(int part = 0; part < size(); part++){
		memset(checkVAL, 0, at(part)->branchNum*sizeof(double));
		for (i = 0; i < nodes1.size(); i++){
			if(((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part])
				checkVAL[((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part]->id] += nodes1[i]->findNeighbor(nodes2[i])->length * part_info[part].part_rate;
		}
		NodeVector nodes1_sub, nodes2_sub;
		at(part)->getBranches(nodes1_sub, nodes2_sub);
		for(int j = 0; j<nodes1_sub.size();j++)
			if(fabs(nodes1_sub[j]->findNeighbor(nodes2_sub[j])->length-checkVAL[nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id])>0.0001){
				//drawTree(cout, WT_BR_SCALE + WT_INT_NODE + WT_BR_LEN);
				printMapInfo();
				cout<<endl;
				cout<<"Partition = "<<part<<", Branch id = "<<nodes1_sub[j]->findNeighbor(nodes2_sub[j])->id<<endl;
				outError("Branches on SuperTree and SubTree do not match!!",__func__);
			}

	}
	delete [] checkVAL;*/

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
	cout<<"Case 1: NNI on SuperTree => NNI on SubTree: "<<allNNIcases_computed[0]<<endl;
	cout<<"Relinking..."<<endl;
	cout<<"Case 2: branch -> empty: "<<allNNIcases_computed[1]<<endl;
	cout<<"Case 3: empty  -> empty: "<<allNNIcases_computed[2]<<endl;
	cout<<"Case 4: branch -> new  : "<<allNNIcases_computed[3]<<endl;
	cout<<"Case 5: branch -> old || empty -> branch: "<<allNNIcases_computed[4]<<endl;
}

void PhyloSuperTreePlen::computeBranchLengths()
{
	}
