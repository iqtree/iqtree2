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

NNIMove PhyloSuperTreePlen::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves)
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

void PhyloSuperTreePlen::doNNIs(int nni2apply, bool changeBran) {
	IQTree::doNNIs(nni2apply, changeBran);
	mapBranchLen();
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

		//cout<<"Partition "<<part<<" : Epsilon = "<<epsilon_cnt<<endl;
		if(epsilon_cnt == 0){
			nni_type[part]=NNI_NO_EPSILON;
//			allNNIcases_computed[0]++;
		}else if(epsilon_cnt == 1){
			nni_type[part] = NNI_ONE_EPSILON;
//			allNNIcases_computed[1]++;
		}else if(epsilon_cnt == 2){
			nni_type[part]=NNI_TWO_EPSILON;
//			allNNIcases_computed[2]++;
		}else if(epsilon_cnt == 3){
			nni_type[part]=NNI_THREE_EPSILON;
//			allNNIcases_computed[3]++;
		}else {
			nni_type[part] = NNI_MANY_EPSILON;
//			allNNIcases_computed[4]++;
		}
	}


}

void PhyloSuperTreePlen::doNNI(NNIMove &move, bool clearLH)
{
	checkBranchLen();
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	SuperNeighbor *node1_nei = (SuperNeighbor*)*move.node1Nei_it;
	SuperNeighbor *node2_nei = (SuperNeighbor*)*move.node2Nei_it;

	int part = 0, ntrees = size();
	iterator it;
	//double old_brlen = nei1->length;
	vector<NNIMove> part_move;
	vector<NNIType> is_nni;
	part_move.resize(ntrees);
	getNNIType(move.node1, move.node2, is_nni);


	for (it = begin(), part = 0; it != end(); it++, part++) {

		if(is_nni[part] == NNI_NO_EPSILON){
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			part_move[part].node1 = (PhyloNode*)nei2_part->node;
			part_move[part].node2 = (PhyloNode*)nei1_part->node;
			part_move[part].node1Nei_it = part_move[part].node1->findNeighborIt(node1_nei->link_neighbors[part]->node);
			part_move[part].node2Nei_it = part_move[part].node2->findNeighborIt(node2_nei->link_neighbors[part]->node);
		}
	}
//	PhyloTree::doNNI(move,clearLH);
	PhyloTree::doNNI(move,false);
	nei1->length = move.newLen[0];
	nei2->length = move.newLen[0];
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
				node1 = (PhyloNode*)nei2->link_neighbors[part]->node;
				node2 = (PhyloNode*)nei1->link_neighbors[part]->node;
				nei1->link_neighbors[part]->clearPartialLh();
				nei2->link_neighbors[part]->clearPartialLh();
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_TWO_EPSILON:
			node1 = (PhyloNode*)nei2->link_neighbors[part]->node;
			node2 = (PhyloNode*)nei1->link_neighbors[part]->node;
			linkBranch(part, nei1, nei2);
			if(clearLH && !(PhyloNode*)nei2->link_neighbors[part]){
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_THREE_EPSILON:
			linkBranch(part, nei1, nei2);
			if (clearLH) {
				// clear partial likelihood vector
				node1 = (PhyloNode*)nei2->link_neighbors[part]->node;
				node2 = (PhyloNode*)nei1->link_neighbors[part]->node;
				node2->clearReversePartialLh(node1);
				node1->clearReversePartialLh(node2);
			}
			break;
		case NNI_MANY_EPSILON:
			break;
		}

	}
	//mapBranchLen();

//	printTree(cout);
//	cout<<endl;
//	cout<<"NNI was applied"<<endl;
//	checkBranchLen();

}




double PhyloSuperTreePlen::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, SwapNNIParam *nni_param) {

//	for (iterator it = begin(); it != end(); it++)
//		if ((*it)->sse != LK_EIGEN_SSE)
//			outError("hey!");

	assert(node1->degree() == 3 && node2->degree() == 3);
	//cout<<"starting NNI evaluation"<<endl;
	checkBranchLen();

	int i = 0, id = 0;
	int part, ntrees = size();

	/*===========================================================================================
	 * Identify NNIType for partitions
	 *===========================================================================================*/
	vector<NNIType> is_nni;
	getNNIType(node1, node2, is_nni);
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

	//==================================================================================================
	// SuperTREE: saving Neighbors and allocating new ones; assign which nodes/neighbors to be swapped.
	//==================================================================================================
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
		*saved_it[id] = new SuperNeighbor(saved_nei[id]->node, saved_nei[id]->length);
		(*saved_it[id])->id = saved_nei[id]->id;
		for(part = 0; part < ntrees; part++)
			((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
	}

	// Getting NEW Neighbors: get the Neighbors again since they were saved for restoring purpose and replaced by new
	SuperNeighbor *nei1_new = (SuperNeighbor*) node1->findNeighbor(node2);
	SuperNeighbor *nei2_new = (SuperNeighbor*) node2->findNeighbor(node1);

	/* -------------------------------------------------------------------------------------------
	 *  NNI details: assigning nodes to be swapped on SuperTree
	 * -------------------------------------------------------------------------------------------*/

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

	for(part = 0; part < ntrees; part++){
		if(is_nni[part]==NNI_NO_EPSILON){
			//evalNNIs++;
			//part_info[part].evalNNIs++;

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

		} else if(is_nni[part]==NNI_ONE_EPSILON){

			// Make sure to update all the necessary link_neighbors and take care of branch lengths
			// (increase/decrease by central branch where necessary).

			nei1_new->link_neighbors[part] = saved_nei[0]->link_neighbors[part];
			nei2_new->link_neighbors[part] = saved_nei[1]->link_neighbors[part];

			// Change the length of branch, (node1,node2) WAS linked to (-=)
			nei1_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
			nei2_new->link_neighbors[part]->length -= old_brlen * part_info[part].part_rate;
			assert(nei1_new->link_neighbors[part]->length >= 0.0);

			// Allocate three new PhyloNeighbors. For nni1 only one of it will be actually used and which one depends on the NNI.
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

					*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, nei->link_neighbors[part]->length);
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
					((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
					(*sub_saved_it[part*6 + id])->id = nei->link_neighbors[part]->id;

					if(params->nni5){
						((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
					}
					//cout<<"saved_it["<<id<<"]; neighbor->node->id = "<<(*sub_saved_it[part*6 + id])->node->id<<endl;
					//cout<<"saved_it["<<id<<"];           node->id = "<<node_link->id<<endl;
				}
			}
		}else if(is_nni[part]==NNI_THREE_EPSILON && params->nni5){
			for(id = 2; id < 6; id++){
				if(saved_nei[id]->link_neighbors[part]){
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = saved_nei[id]->link_neighbors[part];
				}
			}
		}else if(is_nni[part]==NNI_TWO_EPSILON && params->nni5){
			for(id = 2; id < 6; id++){
				if(saved_nei[id]->link_neighbors[part]){
					((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = saved_nei[id]->link_neighbors[part];
				}
			}
		}
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
	for (cnt = 0; cnt < node2_its.size(); cnt++) {
		//cout<<"NNI Loop-----------------------------NNI."<<cnt<<endl;

		node2_it = node2_its[cnt];
		Neighbor *node2_nei = *node2_it;

		// Define which nodes/neighbors to be swapped on SubTree ----------------------------
		for(part=0; part<ntrees; part++)
			if(is_nni[part]==NNI_NO_EPSILON){
				node1_link[part] = (PhyloNode*) nei2_new->link_neighbors[part]->node;
				node2_link[part] = (PhyloNode*) nei1_new->link_neighbors[part]->node;
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
				//cout<<"NO_EPS: do NNI swap"<<endl;
				//allNNIcases_computed[0] += 1;

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
				//cout<<"MANY_EPS: do nothing"<<endl;
				// the NNI on SuperTree does not change anything on SubTree

			} else if(is_nni[part]==NNI_THREE_EPSILON){
				//cout<<"THREE_EPS: relink"<<endl;

				// The central branch had no image before the NNI.
				// Relink the central branch and take care of branch lengths.
				// In the end restore one branch (valid for both nni1 and nni5).

				linkBranch(part, nei1_new, nei2_new);
				assert(nei1_new->link_neighbors[part]);

				// Save the branch length
				if(cnt == 0)
					sub_saved_branch[6*part] = nei1_new->link_neighbors[part]->length;

				nei1_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;

				// since the branch length was changed we have to recompute the likelihood of the branch
				part_info[part].cur_score = at(part)->computeLikelihoodBranch(nei1_new->link_neighbors[part],
						(PhyloNode*)nei2_new->link_neighbors[part]->node);

			}else if(is_nni[part]==NNI_TWO_EPSILON){
				//cout<<"TWO_EPS: relink"<<endl;

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
							(PhyloNode*)saved_nei[1]->link_neighbors[part]->node);
				}

			}else if(is_nni[part] == NNI_ONE_EPSILON){
				//cout<<"ONE_EPS: relink, update the link_neighbors"<<endl;

				/* The crazy case, which absorbs most of the bugs:(
				 * Lets say on SuperTree there are five branches, a,b,c,d and central e, and d has an empty image.
				 * The corresponding SubTree has 3 branches, a',b',c'.
				 * Before NNI central branch, e, has an image. Lets say it maps to a'.
				 * After NNI it will be remapped either to b' or c', depending on which nodes will be swapped.
				 * Update the corresponding link_neighbors. Make sure that link_neighbors for central branch e
				 * and for the one it is now mapped to (b' or c'), are the same.
				 * Decrease a' (done before). Increase b' or c' depending on the NNI. Restore three branches.*/

				if(cnt == 1){
					for(id=2; id<6; id++){
						if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
							if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node == nei1_new->link_neighbors[part]->node){
								nei2_new->link_neighbors[part]->clearPartialLh();
								break;
							} else if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node == nei2_new->link_neighbors[part]->node){
								nei1_new->link_neighbors[part]->clearPartialLh();
								break;
							}
						}
					}

				}

				linkBranch(part, nei1_new, nei2_new);
				assert(nei1_new->link_neighbors[part]);

				//cout<<"nei1_new->link_nei["<<part<<"]->node->id"<<nei1_new->link_neighbors[part]->node->id<<endl;
				//cout<<"nei2_new->link_nei["<<part<<"]->node->id"<<nei2_new->link_neighbors[part]->node->id<<endl;

				assert(nei1_new->link_neighbors[part]->node->findNeighbor(nei2_new->link_neighbors[part]->node));
				assert(nei2_new->link_neighbors[part]->node->findNeighbor(nei1_new->link_neighbors[part]->node));

				// nni1:
				// - you need to update only one link_neighbor with new PhyloNeighbor
				//	 (either node1->findNeighbor(node2) or node2->findNeighbor(node1))
				// - the second is already linked to some existing PhyloNeighbor after linkBranch().
				for(id=2; id<6; id++){
					if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
						if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node == nei1_new->link_neighbors[part]->node){
							nei2_new->link_neighbors[part] = (PhyloNeighbor*)(*sub_saved_it[part*6 + id]);
							if(cnt == 1)
								nei2_new->link_neighbors[part]->clearPartialLh();
							//cout<<"nei2_new->link_neighbors[part] is newed"<<endl;
							break;
						} else if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]->node == nei2_new->link_neighbors[part]->node){
							nei1_new->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
							if(cnt == 1)
								nei1_new->link_neighbors[part]->clearPartialLh();
							//cout<<"nei1_new->link_neighbors[part] is newed"<<endl;
							break;
						}
					}
				}

				//cout<<"nei1_new->link_nei["<<part<<"]->node->id"<<nei1_new->link_neighbors[part]->node->id<<endl;
				//cout<<"nei2_new->link_nei["<<part<<"]->node->id"<<nei2_new->link_neighbors[part]->node->id<<endl;

				assert(nei1_new->link_neighbors[part]->node->findNeighbor(nei2_new->link_neighbors[part]->node));
				assert(nei2_new->link_neighbors[part]->node->findNeighbor(nei1_new->link_neighbors[part]->node));

				// Increase the branch to which the central is relinked.
				nei1_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;
				nei2_new->link_neighbors[part]->length += old_brlen * part_info[part].part_rate;

				// This is done before!
//				// nni5: you need to update also the link_neighbors for corresponding neighbors of node1 and node2
//				if(params->nni5){
//					for(id=2; id<6; id++){
//						if(saved_nei[id]->link_neighbors[part]){
//							((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
//						}
//					}
//				}

			} // end of else ONE_EPS case
		} // end of part loop

/*===============================================================================================================================*
 * 											Compute the score of the swapped topology 				  							 *
 *===============================================================================================================================*/
		//cout<<"Before optimization"<<endl;
		//checkBranchLen();
		double score = optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
		if (verbose_mode >= VB_MED) {
			cout << "[" << score << "] ";
			printTree(cout);
			cout << endl;
		}
		//cout<<"After optimization"<<endl;
		//checkBranchLen();

		// %%%%%%%%%%%%%%%%%%%%%%%%  FIVE BRANCH OPTIMIZATION  %%%%%%%%%%%%%%%%%%%%%%%%
	    if (params->nni5) {
	    	FOR_NEIGHBOR(node1, node2, it){
	    		for(part = 0; part < ntrees; part++)
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON or is_nni[part]==NNI_ONE_EPSILON)){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei2_new->link_neighbors[part]->node; // this should be node 1 on subtree
	    				// the problem is that for ONE_epsilon case node1 on subtree is equal to its neighbor node on subtree
	    				// in this case we have to set nei_link to node2 on subtree
	    				if(node_link->id == nei_link->id){
	    					nei_link = nei1_new->link_neighbors[part]->node;
	    				}
	    				//cout<<"HERE it is: "<<((SuperNeighbor*)(*it))->link_neighbors[part]->node->id<<endl;
	    				//cout<<nei2_new->link_neighbors[part]->node->id<<endl;
						((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		//cout<<"NNI5 : node1 : Before optimization"<<endl;
	    		//checkBranchLen();
	    		score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false, NNI_MAX_NR_STEP);
	    		//cout<<"NNI5 : node1 : After optimization"<<endl;
	    		//checkBranchLen();
	    	}
	    	for(part = 0; part < ntrees; part++)
	    		if(((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON or is_nni[part]==NNI_ONE_EPSILON))
	    			((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part]->clearPartialLh();
	    	FOR_NEIGHBOR(node2, node1, it){
	    		for(part = 0; part < ntrees; part++)
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part] && (is_nni[part]==NNI_NO_EPSILON or is_nni[part]==NNI_ONE_EPSILON)){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei1_new->link_neighbors[part]->node;
	    				if(node_link->id == nei_link->id){
	    					nei_link = nei2_new->link_neighbors[part]->node;
	    				}
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
		if(params->nni5)
		for(id=2;id<6;id++){
			(*saved_it[id])->length = saved_nei[id]->length;
			(*node_nei_it[id-2])->length = saved_nei[id]->length;
		}

// Swap back or relink back on SubTrees------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){

			if(is_nni[part]==NNI_NO_EPSILON){
				node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				//Restoring the branch length on the Sub Tree
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
				assert(nei1_new->link_neighbors[part]->node == saved_nei[0]->link_neighbors[part]->node);
				assert(nei2_new->link_neighbors[part]->node == saved_nei[1]->link_neighbors[part]->node);

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
			//mapBranchLen();
		}
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

			} else if(is_nni[part] == NNI_ONE_EPSILON){

				// Delete the allocated neighbors and restore from saved neighbors
				for (id = 5; id >= 2; id--) {
					//if((*sub_saved_it[part*6+id])){
					if(((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+id])->scale_num;
						delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+id])->partial_lh;
						if (*sub_saved_it[part*6+id] == at(part)->current_it)
							at(part)->current_it = saved_nei[id]->link_neighbors[part];
						if (*sub_saved_it[part*6+id] == at(part)->current_it_back)
							at(part)->current_it_back = saved_nei[id]->link_neighbors[part];

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
		checkBranchLen();
//------------------------------------------------------------------------------------------------------------------------------------------------
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
void PhyloSuperTreePlen::restoreAllBrans(PhyloNode *node, PhyloNode *dad) {
	IQTree::restoreAllBrans(node,dad);
	mapTrees();
}

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
