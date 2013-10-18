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
	fixed_rates = false;
}

PhyloSuperTreePlen::PhyloSuperTreePlen(Params &params)
: PhyloSuperTree(params)
{
	fixed_rates = (params.partition_type == 'j') ? true : false;
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
	//initPartitionInfo(); // OLGA: not needed here
	//cout<<"Optimizing all branches"<<endl;
	for(int part = 0; part < size(); part++){
		part_info[part].cur_score = 0.0;
	}
	return PhyloTree::optimizeAllBranches(my_iterations,tolerance);
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

NNIMove PhyloSuperTreePlen::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, bool approx_nni, bool useLS, double lh_contribution)
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
	double oldLEN = node1->findNeighbor(node2)->length;

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
		myMove.oldLen[0] = oldLEN;
	} else {
		myMove.swap_id = 2;
		myMove.node1Nei_it = node1->findNeighborIt(node1_nei->node);
		myMove.node2Nei_it = node2->findNeighborIt(node2_nei_other->node);
		myMove.newloglh = nni_param.nni2_score;
		myMove.node1 = node1;
		myMove.node2 = node2;
		myMove.newLen[0] = nni_param.nni2_brlen;
		myMove.oldLen[0] = oldLEN;
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
				if (!is_nni) {
					//memcpy(at(part)->_pattern_lh, part_info[part].opt_ptnlh[brid], at(part)->getAlnNPattern() * sizeof(double));
				} else if (nnino == 0)
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

void PhyloSuperTreePlen::applyNNIs(int nni2apply, bool changeBran) {
	IQTree::applyNNIs(nni2apply, changeBran);
	mapBranchLen();
}


void PhyloSuperTreePlen::doNNI(NNIMove &move, bool clearLH)
{
	//checkBranchLen();
	//cout<<endl<<"I did NNI on Super Tree!!:)"<<endl;
	SuperNeighbor *nei1 = (SuperNeighbor*)move.node1->findNeighbor(move.node2);
	SuperNeighbor *nei2 = (SuperNeighbor*)move.node2->findNeighbor(move.node1);
	SuperNeighbor *node1_nei = (SuperNeighbor*)*move.node1Nei_it;
	SuperNeighbor *node2_nei = (SuperNeighbor*)*move.node2Nei_it;

//	cout<<endl<<"node1:"<<move.node1->name<<","<<move.node1->id<<endl;
//	cout<<"node2:"<<move.node2->name<<","<<move.node2->id<<endl;
//	cout<<"node1_nei:"<<(*move.node1Nei_it)->node->name<<","<<(*move.node1Nei_it)->node->id<<endl;
//	cout<<"node2_nei:"<<(*move.node2Nei_it)->node->name<<","<<(*move.node2Nei_it)->node->id<<endl;

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
				if(nei1_part->length<0)
					outError("shit!!!",__func__);
			}
		}
	}
//	cout<<endl<<endl<<"BEFORE I DID NNI o.O"<<endl<<endl;
//	printMapInfo();
	PhyloTree::doNNI(move,clearLH);
	nei1->length = move.newLen[0];
	nei2->length = move.newLen[0];

	for (it = begin(), part = 0; it != end(); it++, part++) {
		if (!is_nni[part]) {
			// Relink the branch if it does not correspond to NNI for partition
			linkBranch(part, nei1, nei2);
			//cout<<endl<<endl<<"AFTER RELINK in I DID NNI o.O"<<endl<<endl;
			//printMapInfo();
		} else { (*it)->doNNI(part_move[part],clearLH); }

//		PhyloNeighbor* nei1_part_new = nei1->link_neighbors[part];
//		PhyloNeighbor* nei2_part_new = nei2->link_neighbors[part];
//		if(nei1_part_new){
//			int brid = nei1_part_new->id;
//			if (move.swap_id == 1)
//				nei1_part_new->length = nei2_part_new->length = part_info[part].nni1_brlen[brid];
//			else
//				nei1_part_new->length = nei2_part_new->length = part_info[part].nni2_brlen[brid];
//
//			if(nei1_part_new->length < 0){
//				outError("Negative branch in doNNI!!!",__func__);
//			}
//
//		}
	}
	//mapBranchLen();

//	cout<<endl<<endl<<"AFTER I DID NNI o.O"<<endl<<endl;
//	printMapInfo();
//	printTree(cout);
//	cout<<endl;
	//checkBranchLen();

}

double PhyloSuperTreePlen::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2, SwapNNIParam *nni_param) {

//	cout<<endl<<"****************************************************************************"<<endl;;
//	cout<<endl<<"NEW branch for NNI analysis ------------------------------------------------------------------------------------------"<<endl;
//	printTree(cout);
//	cout<<endl;
//	printMapInfo();
//	cout<<endl;
	assert(node1->degree() == 3 && node2->degree() == 3);
	checkBranchLen();

	//===========================================================================================
	// Prepare the details for NNI: SuperTree
	//===========================================================================================
	int i = 0, id = 0;
	int part, ntrees = size();
	double old_brlen = node1->findNeighbor(node2)->length; // length of the branch between node1 and node2 on SuperTree before NNI
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
		(*saved_it[id])->id = saved_nei[id]->id;
		for(part = 0; part < ntrees; part++)
			((SuperNeighbor*)*saved_it[id])->link_neighbors.push_back(NULL);
	}


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
	vector<bool> is_nni;
	is_nni.resize(ntrees);

	Node *node_link, *nei_link;
	vector<PhyloNode*> node1_link,node2_link;
	vector<PhyloNeighbor*> nei1_new_part, nei2_new_part, node1_link_nei, node2_link_nei,sub_saved_nei1,sub_saved_nei2;
	vector<NeighborVec::iterator> node1_link_it, node2_link_it, sub_saved_it;

	nei1_new_part.resize(ntrees);
	nei2_new_part.resize(ntrees);
	node1_link.resize(ntrees);
	node2_link.resize(ntrees);
	node1_link_nei.resize(ntrees);
	node2_link_nei.resize(ntrees);
	node1_link_it.resize(ntrees);
	node2_link_it.resize(ntrees);
	sub_saved_nei1.resize(ntrees);
	sub_saved_nei2.resize(ntrees);

	sub_saved_it.resize(6*ntrees);
	/*------------------------------------------------------------------------------------
	 *  to know which nodes to restore:
	 * 		3 spots for node1 as a neighbor of [node2, node1_nei1, node1_nei2] and
	 * 		3 spots for node2 as a neighbor of [node1, node2_nei1, node2_nei2]
	 *------------------------------------------------------------------------------------*/

	/*------------------------------------------------------------------------------------
	 * At this point we work only with partitions satisfying is_nni condition:
	 * 		- check if is_nni[part]
	 * 		- allocate new PhyloNeighbor
	 * 		- check if some other branch links to the same branch on SubTree
	 * 		- update link_neighbors[part]
	 * 		- synchronize node1_link[part], node2_link[part], node1_link_nei[part]
	 *------------------------------------------------------------------------------------*/

	for(part = 0; part < ntrees; part++){
		is_nni[part] = true;
		if(saved_nei[0]->link_neighbors[part]){
			FOR_NEIGHBOR_DECLARE(node1,node2,nit){
				if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}
			FOR_NEIGHBOR(node2, node1, nit) {
				if(!((SuperNeighbor*)*nit)->link_neighbors[part]) { is_nni[part] = false; break; }
			}
		} else { is_nni[part] = false; }


		if(is_nni[part]){
			// one branch optimization
			for(id = 0; id < 2; id++){
				nei_link  = saved_nei[id]->link_neighbors[part]->node;
				node_link = saved_nei[1-id]->link_neighbors[part]->node;
				sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
				*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, saved_nei[id]->link_neighbors[part]->length);
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
				((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();
				(*sub_saved_it[part*6 + id])->id = saved_nei[id]->link_neighbors[part]->id;

				// update link_neighbor[part]
				((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];
			}

			// should be changed something.....
			// optimization on 5 branches
			if(params->nni5Branches){
				for(id = 2; id < 6; id ++){
					/*
					 * here check if neighbor of node1 or node2 has link_neighbor filled in already by linkCheck
					 * if no, allocate new PhyloNeighbor, if yes, do not do anything
					 */
					if (!((SuperNeighbor*)(*node_nei_it[id-2]))->link_neighbors[part]){
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
			}
			// Synchronization: node1, node2, node1_nei; node2_nei varies for two NNIs therefore synchronized later
			node1_link[part] = (PhyloNode*) nei2_new->link_neighbors[part]->node;
			node2_link[part] = (PhyloNode*) nei1_new->link_neighbors[part]->node;
			node1_link_nei[part] = ((SuperNeighbor*)node1_nei)->link_neighbors[part];
			node1_link_it[part] = node1_link[part]->findNeighborIt(node1_link_nei[part]->node);
		}
	}
	//===========================================================================================
	// Do the NNI swap and compute the likelihood of swapped topology
	//===========================================================================================
	int cnt;
	for (cnt = 0; cnt < node2_its.size(); cnt++) {
//		cout<<"====================== Checking NNI "<<cnt<<" ======================"<<endl;
		node2_it = node2_its[cnt];
		Neighbor *node2_nei = *node2_it;

		// Synchronize node2_link_nei
		for(part=0; part<ntrees; part++)
			if(is_nni[part]){
				node2_link_nei[part] = ((SuperNeighbor*)node2_nei)->link_neighbors[part];
				node2_link_it[part] = node2_link[part]->findNeighborIt(node2_link_nei[part]->node);
			}
		// Do the NNI swap on SuperTrees ----------------------------------------------------
		node1->updateNeighbor(node1_it, node2_nei);
		node2_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node1_nei);
		node1_nei->node->updateNeighbor(node1, node2);

		if(verbose_mode == VB_MAX){
		//if(cnt == 0){
			cout<<"============================================"<<endl;
			cout<<"NNI on SuperTree:"<<endl;
			cout<<"node1:"<<node1->name<<" id:"<<node1->id<<endl;
			cout<<"node2:"<<node2->name<<" id:"<<node2->id<<endl;
			cout<<"node1_nei: "<<node1_nei->node->name<<" id:"<<node1_nei->node->id<<endl;
			cout<<"node2_nei: "<<node2_nei->node->name<<" id:"<<node2_nei->node->id<<endl;
		//}
		}

// Do the NNI swap on SubTrees or Proceed with partitions requiring relink of (node1,node2) ------------------------------
		for(part = 0; part < ntrees; part++){
			if(is_nni[part]){
				// Do NNI swap on partition
				node1_link[part]->updateNeighbor(node1_link_it[part], node2_link_nei[part]);
				node2_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
				node2_link[part]->updateNeighbor(node2_link_it[part], node1_link_nei[part]);
				node1_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);

				if(verbose_mode == VB_MAX){
					//cout<<"NNI on SubTree: part = "<<part<<" ==============="<<endl;
					cout<<"node1_link:"<<node1_link[part]->name<<" id:"<<node1_link[part]->id<<endl;
					cout<<"node2_link:"<<node2_link[part]->name<<" id:"<<node2_link[part]->id<<endl;
					cout<<"node1_link_nei: "<<node1_link_nei[part]->node->name<<" id:"<<node1_link_nei[part]->node->id<<endl;
					cout<<"node2_link_nei: "<<node2_link_nei[part]->node->name<<" id:"<<node2_link_nei[part]->node->id<<endl;
				}
				nei1_new_part[part] = nei1_new->link_neighbors[part];
				nei2_new_part[part] = nei2_new->link_neighbors[part];
				if(nei1_new_part[part]){
					nei1_new_part[part]->clearPartialLh();
					nei2_new_part[part]->clearPartialLh();
				}

			} else {
				/*
				 * If it is not NNI on SubTree, but branch is linked to some branch on SubTree, RELINK
				 *
				 * 	- first update the link_neighbor[part] of nei1_new and nei2_new,
				 * let it points to the same link_neighbor as in saved[0]->link_nei and saved[1]->link_nei.
				 * By doing so you will change the length of the branch node1, node2 WAS mapped to and
				 * you can now do relink and NEW link_neighbors will be saved in saved[0] and saved[1].
				 *
				 * 	- then allocate new PhyloNeighbor, change the NEW branch length +=
				 * 	- check if some other branches from SuperTree link to the same branch on SubTree
				 * 		YES: update link_neighbors to new PhyloNeighbors
				 * 	- after optimization: delete new, restore phyloNeighbors and link_neighbors, relink, increase the old branch length
				 */

				if(saved_nei[0]->link_neighbors[part]){

					nei1_new->link_neighbors[part] = saved_nei[0]->link_neighbors[part];
					nei2_new->link_neighbors[part] = saved_nei[1]->link_neighbors[part];

					// Change the length of branch, (node1,node2) WAS linked to (-=)
					nei1_new_part[part] = nei1_new->link_neighbors[part];
					nei2_new_part[part] = nei2_new->link_neighbors[part];

//						cout<<endl<<endl<<"BEFORE DECREASE BRANCH LENGTH"<<endl<<endl;
//						at(part)->printTree(cout);
//						cout<<endl;
//						cout<<"nei1_old-link:"<<nei1_new->link_neighbors[part]->length<<endl;
//						cout<<"nei2_old-link:"<<nei2_new->link_neighbors[part]->length<<endl;
//						cout<<"saved_nei[0]->link:"<<saved_nei[0]->link_neighbors[part]->length<<endl;
//						cout<<"saved_nei[1]->link:"<<saved_nei[1]->link_neighbors[part]->length<<endl;

					nei1_new_part[part]->length -= old_brlen * part_info[part].part_rate;
					nei2_new_part[part]->length -= old_brlen * part_info[part].part_rate;

//						cout<<endl<<endl<<"AFTER DECREASE BRANCH LENGTH"<<endl<<endl;
//						at(part)->printTree(cout);
//						cout<<endl;
//						cout<<"nei1_old-link:"<<nei1_new->link_neighbors[part]->length<<endl;
//						cout<<"nei2_old-link:"<<nei2_new->link_neighbors[part]->length<<endl;

					if(nei1_new_part[part]->length == 0.0){
						nei1_new_part[part]->length = MIN_BRANCH_LEN;
						nei2_new_part[part]->length = MIN_BRANCH_LEN;
					}

//						cout<<"BEFORE RELINKING"<<endl;
//						cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;
//						at(part)->printTree(cout);
//						cout<<endl;
//						cout<<"nei1_old-link:"<<nei1_new->link_neighbors[part]->length<<endl;
//						cout<<"nei2_old-link:"<<nei2_new->link_neighbors[part]->length<<endl;
				}

				// Relink the branch if it does not correspond to NNI for partition
				linkBranch(part, nei1_new, nei2_new);
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

				if(nei1_new->link_neighbors[part]){ // otherwise it's NULL
					sub_saved_nei1[part] = nei1_new->link_neighbors[part];
					sub_saved_nei2[part] = nei2_new->link_neighbors[part];

//						cout<<endl<<endl<<"AFTER relink"<<endl<<endl;
//						//printMapInfo();
//						at(part)->printTree(cout);
//						cout<<endl;
//						cout<<"nei1_new-link:"<<nei1_new->link_neighbors[part]->length<<endl;
//						cout<<"nei2_new-link:"<<nei2_new->link_neighbors[part]->length<<endl;
//						cout<<"saved-link:"<<sub_saved_nei1[part]->length<<endl;
//						cout<<"saved-link:"<<sub_saved_nei2[part]->length<<endl;
//						cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne1_new -> "<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<","<<nei2_new->link_neighbors[part]->node->id<<endl;


					// Allocate new PhyloNeighbors
					// one branch optimization
					for(id = 0; id < 2; id++){
//							cout<<"WHEN allocating new ones"<<endl;
//							cout<<"saved-link:"<<sub_saved_nei1[part]->length<<endl;
//							cout<<"saved-link:"<<sub_saved_nei2[part]->length<<endl;
						nei_link  =  (id == 0) ? sub_saved_nei1[part]->node : sub_saved_nei2[part]->node;
						node_link =  (id == 1) ? sub_saved_nei1[part]->node : sub_saved_nei2[part]->node;
						sub_saved_it[part*6 + id] = node_link->findNeighborIt(nei_link);
						*sub_saved_it[part*6 + id] = new PhyloNeighbor(nei_link, sub_saved_nei1[part]->length + old_brlen * part_info[part].part_rate);
						((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->partial_lh = at(part)->newPartialLh();
						((PhyloNeighbor*) (*sub_saved_it[part*6 + id]))->scale_num = at(part)->newScaleNum();

						PhyloNeighbor *savednei = (id == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
						if (savednei->partial_lh_computed) {
							if (saved_nei[0]->link_neighbors[part] && common_nodes == 1) {
								// Minh: for efficiency, if previous branch maps to some branch, reuse partial_lh if it does not look to old branch
								if (nei_link != saved_nei[0]->link_neighbors[part]->node && nei_link != saved_nei[1]->link_neighbors[part]->node) {

									PhyloNeighbor *thisnei = ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]));
									memcpy(thisnei->partial_lh, savednei->partial_lh, at(part)->getPartialLhBytes());
									memcpy(thisnei->scale_num , savednei->scale_num , at(part)->getScaleNumBytes());
									thisnei->lh_scale_factor = savednei->lh_scale_factor;
									thisnei->partial_lh_computed = savednei->partial_lh_computed;
								}
							} else {
								// Minh: for efficency, if previous branch maps to NULL, reuse partial_lh
								PhyloNeighbor *thisnei = ((PhyloNeighbor*) (*sub_saved_it[part*6 + id]));
								memcpy(thisnei->partial_lh, savednei->partial_lh, at(part)->getPartialLhBytes());
								memcpy(thisnei->scale_num , savednei->scale_num , at(part)->getScaleNumBytes());
								thisnei->lh_scale_factor = savednei->lh_scale_factor;
								thisnei->partial_lh_computed = savednei->partial_lh_computed;
								//if (thisnei->partial_lh_computed)
									//cout << "yh" << endl;
							}
						}

						(*sub_saved_it[part*6 + id])->id = sub_saved_nei1[part]->id;

						// update link_neighbor[part]
						((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = (PhyloNeighbor*)*sub_saved_it[part*6 + id];

					}

//						if(nei1_new->link_neighbors[part] == (PhyloNeighbor*)(*sub_saved_it[part*6+0]))	{cout<<"YESSSS"<<endl;}
//						else { cout<<"NEIN!!!!"<<endl;}
//						if(nei2_new->link_neighbors[part] == (PhyloNeighbor*)(*sub_saved_it[part*6+1])) {cout<<"YESSSS"<<endl;}
//						else {cout<<"NEIN!!!"<<endl;}

					//cout<<"================================"<<endl;
					linkCheck(part,node1,node2,sub_saved_nei2[part]);
					//cout<<"================================"<<endl;
					linkCheck(part,node2,node1,sub_saved_nei1[part]);
					//cout<<"================================"<<endl;


					nei1_new_part[part] = nei1_new->link_neighbors[part];
					nei2_new_part[part] = nei2_new->link_neighbors[part];
					/*if(nei1_new_part[part] && saved_nei[0]->link_neighbors[part]){
						nei1_new_part[part]->clearPartialLh();
						nei2_new_part[part]->clearPartialLh();
					}*/


//						cout<<"AFTER RELINKING, NEW"<<endl;
//						cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;

					// optimization on 5 branches!!!!!
				}else{
					if(saved_nei[0]->link_neighbors[part])
						part_info[part].cur_score = at(part)->computeLikelihoodBranch(saved_nei[0]->link_neighbors[part],(PhyloNode*)saved_nei[1]->link_neighbors[part]->node);
//						cout<<"RELINKed to nothing........."<<endl;
//						for(id = 0; id<2; id++)
//							((SuperNeighbor*)*saved_it[id])->link_neighbors[part] = NULL;
//						if(nei1_new->link_neighbors[part] == NULL){cout<<"YESSSS, link_nei is NULL"<<endl;}
//						else{cout<<"NEIN!!!!"<<endl;}
//						if(nei2_new->link_neighbors[part] == NULL){cout<<"YESSSS, link_nei is NULL"<<endl;}
//						else{cout<<"NEIN!!!"<<endl;}
					nei1_new_part[part] = nei1_new->link_neighbors[part];
					nei2_new_part[part] = nei2_new->link_neighbors[part];
				}
			} // end RELINK

			//part_info[part].cur_score = 0.0;

		} // end of part loop

/*===============================================================================================================================*
 * 											Compute the score of the swapped topology 				  							 *
 *===============================================================================================================================*/
//		cout<<endl<<endl<<"BEFORE OPTIMIZATION"<<endl<<endl;
//		printMapInfo();
		checkBranchLen();

		double score = optimizeOneBranch(node1, node2, false);
		if (verbose_mode >= VB_MED) {
			cout << "[" << score << "] ";
			printTree(cout);
			cout << endl;
		}

//		cout<<endl<<endl<<"AFTER OPTIMIZATION"<<endl<<endl;;
//		printMapInfo();
		checkBranchLen();


		/*  for(part = 0; part < ntrees; part++)
	    	if(nei1_new->link_neighbors[part] && !is_nni[part]){
	    		part_info[part].opt_brlen[brid[part]] = nei1_new_part[part]->length;
	    		at(part)->computePatternLikelihood(part_info[part].opt_ptnlh[brid[part]], &part_info[part].opt_score[brid[part]]);
	    	}*/

		// %%%%%%%%%%%%%%%%%%%%%%%%  FIVE BRANCH OPTIMIZATION  %%%%%%%%%%%%%%%%%%%%%%%%
	    if (params->nni5Branches) {
	    	if (verbose_mode >= VB_DEBUG)
	    		cout << "Log-likelihood: " << score << endl;
	    	FOR_NEIGHBOR(node1, node2, it){
	    		for(part = 0; part < ntrees; part++)
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei2_new->link_neighbors[part]->node;
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false);
	    		if (verbose_mode >= VB_DEBUG)
	    			cout << "Log-likelihood: " << score << endl;
	    	}
	    	for(part = 0; part > ntrees; part++)
	    		if(nei2_new->link_neighbors[part])
	    			nei2_new_part[part]->clearPartialLh();
	    	FOR_NEIGHBOR(node2, node1, it){
	    		for(part = 0; part < ntrees; part++)
	    			if(((SuperNeighbor*)(*it))->link_neighbors[part]){
	    				node_link = ((SuperNeighbor*)(*it))->link_neighbors[part]->node;
	    				nei_link  = nei1_new->link_neighbors[part]->node;
	    				((PhyloNeighbor*)node_link->findNeighbor(nei_link))->clearPartialLh();
	    			}
	    		score = optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false);
	    		if (verbose_mode >= VB_DEBUG)
	    			cout << "Log-likelihood: " << score << endl;
	    	}
	    }
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%  END of nni5branch  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

		/* Minh bug fix: restore the cur_score */
		for (i = 0; i < part_info.size(); i++)
			part_info[i].cur_score = saved_cur_score[i];


	    // *************************** STORE INFO ABOUT NNI ***************************
	    // Store information about this NNI for partitions, this is used for bootstrap...
	    for(part = 0; part < ntrees; part++){
	    	if(nei1_new->link_neighbors[part]){
	    		int brid = nei1_new->link_neighbors[part]->id;
		    	if(cnt == 0){
		    		part_info[part].nni1_brlen[brid] = nei1_new->link_neighbors[part]->length;
		    		//computePatternLikelihood(part_info[part].nni1_ptnlh[brid[part]], &part_info[part].nni1_score[brid[part]]);
		    	}else{
		    		part_info[part].nni2_brlen[brid] = nei1_new->link_neighbors[part]->length;
		    		//computePatternLikelihood(part_info[part].nni2_ptnlh[brid[part]], &part_info[part].nni2_score[brid[part]]);
		    	}
	    	}
	    }
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

		// ========================== RESTORE INFO ABOUT NNI =========================

		//Restoring branch length on Super Tree
		nei1_new->length = old_brlen;
		nei2_new->length = old_brlen;




//		cout<<endl<<endl<<"BEFORE SWAP BACK"<<endl<<endl;;
//		printMapInfo();
		for(part = 0; part < ntrees; part++)
				if(!is_nni[part] && nei1_new->link_neighbors[part]){
//					cout<<"================================"<<endl;
//					cout<<endl<<"PART ššššjjjj"<<part<<endl;
//					cout<<"================================"<<endl;
					linkCheckRe(part,node1,node2,sub_saved_nei2[part],sub_saved_nei1[part]);
//					cout<<"================================"<<endl;
					linkCheckRe(part,node2,node1,sub_saved_nei1[part],sub_saved_nei2[part]);
//					cout<<"================================"<<endl;

//						cout<<endl<<endl<<"PART "<<part<<"     BEFORE deleting NEW"<<endl<<endl;
//						printMapInfo();
						// Restore on SubTree
						for (i = 1; i >= 0; i--) {
							if((*sub_saved_it[part*6+i])){
								delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->scale_num;
								delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->partial_lh;
								if (*sub_saved_it[part*6+i] == at(part)->current_it) at(part)->current_it = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
								if (*sub_saved_it[part*6+i] == at(part)->current_it_back) at(part)->current_it_back = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];

								delete (*sub_saved_it[part*6+i]);
								(*sub_saved_it[part*6+i]) = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
							}
						}
						nei1_new->link_neighbors[part] = sub_saved_nei1[part];
						nei2_new->link_neighbors[part] = sub_saved_nei2[part];

//						cout<<endl<<endl<<"AFTER DELETING NEW"<<endl<<endl;
//						printMapInfo();
					}
// Swap back on SuperTree --------------------------------------------------------------------------------------------------------------
		node1->updateNeighbor(node1_it, node1_nei);
		node1_nei->node->updateNeighbor(node2, node1);
		node2->updateNeighbor(node2_it, node2_nei);
		node2_nei->node->updateNeighbor(node1, node2);
// Swap back or relink back on SubTrees ------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){
			//if(saved_nei[0]->link_neighbors[part]){
				if(is_nni[part]){
					node1_link[part]->updateNeighbor(node1_link_it[part], node1_link_nei[part]);
					node1_link_nei[part]->node->updateNeighbor(node2_link[part], node1_link[part]);
					node2_link[part]->updateNeighbor(node2_link_it[part], node2_link_nei[part]);
					node2_link_nei[part]->node->updateNeighbor(node1_link[part], node2_link[part]);
					//Restoring the branch length on the Sub Tree
					node1_link[part]->findNeighbor(node2_link[part])->length = saved_nei[0]->link_neighbors[part]->length;
					node2_link[part]->findNeighbor(node1_link[part])->length = saved_nei[0]->link_neighbors[part]->length;
				//}
			}else{
//					if(nei1_new->link_neighbors[part]){
//						cout<<endl<<endl<<"BEFORE deleting NEW"<<endl<<endl;
//						printMapInfo();
//						/*cout<<"BEFORE deleting NEW"<<endl;
//						cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;
//						*/
//
//						/*
//						 * delete and restore link_neighbors from saved_nei[0]->link_neighbors and saved[1]->link_neighbors[part]
//						 * make sure you restore link_neighbors also for other branches which linked to the same branch on SubTree
//						 */
//
//						// Restore on SubTree
//						for (i = 1; i >= 0; i--) {
//							if((*sub_saved_it[part*6+i])){
//								delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->scale_num;
//								delete[] ((PhyloNeighbor*) *sub_saved_it[part*6+i])->partial_lh;
//								if (*sub_saved_it[part*6+i] == at(part)->current_it) at(part)->current_it = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
//								if (*sub_saved_it[part*6+i] == at(part)->current_it_back) at(part)->current_it_back = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
//
//								delete (*sub_saved_it[part*6+i]);
//								(*sub_saved_it[part*6+i]) = (i == 0) ? sub_saved_nei1[part] : sub_saved_nei2[part];
//							}
//						}
//
//						// 5branch optimization!!!!!!!!!
//						/*
//						// restore the length of 4 branches around node1_link[part], node2_link[part]
//						node1_link[part] = (PhyloNode*) ((SuperNeighbor*)node2->findNeighbor(node1))->link_neighbors[part]->node;
//						node2_link[part] = (PhyloNode*) ((SuperNeighbor*)node1->findNeighbor(node2))->link_neighbors[part]->node;
//						FOR_NEIGHBOR(node1_link[part], node2_link[part], it)
//							(*it)->length = (*it)->node->findNeighbor(node1_link[part])->length;
//						FOR_NEIGHBOR(node2_link[part], node1_link[part], it)
//							(*it)->length = (*it)->node->findNeighbor(node2_link[part])->length;
//						}*/
//
////						if(nei1_new->link_neighbors[part] == sub_saved_nei1[part]){cout<<"blya"<<endl;}
////						else {cout<<"OKie"<<endl;}
////						if(nei2_new->link_neighbors[part] == sub_saved_nei2[part]){cout<<"blya"<<endl;}
////						else{cout<<"OKie"<<endl;}
//
//						nei1_new->link_neighbors[part] = sub_saved_nei1[part];
//						nei2_new->link_neighbors[part] = sub_saved_nei2[part];
//
//
//						cout<<endl<<endl<<"AFTER DELETING NEW"<<endl<<endl;
//						printMapInfo();
//
//					/*	if(nei1_new->link_neighbors[part] == (PhyloNeighbor*)(*sub_saved_it[part*6+0]))
//							cout<<"YESSSS"<<endl;
//						else
//							cout<<"NEIN!!!!"<<endl;
//						if(nei2_new->link_neighbors[part] == (PhyloNeighbor*)(*sub_saved_it[part*6+1]))
//							cout<<"YESSSS"<<endl;
//						else
//							cout<<"NEIN!!!"<<endl;
//						*/
//					}

//					cout<<endl<<"BEFORE RELINKING BACK"<<endl;
//					if(nei1_new->link_neighbors[part]){
//						cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//						cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;
//					} else {
//						cout<<"it was relinked to nothing...."<<endl;
//					}

					// Relink back
	//				cout<<endl<<endl<<"BEFORE RELINKING BACK"<<endl<<endl;
	//				printMapInfo();

					linkBranch(part, nei1_new, nei2_new);
	//				cout<<endl<<endl<<"AFTER RELINKING BACK"<<endl<<endl;
	//				printMapInfo();

/*
					if(nei1_new->link_neighbors[part] == saved_nei[0]->link_neighbors[part])
						cout<<"GREAT! after relinking nei1->link and saved[0]->link match!!!"<<endl;
					else
						cout<<"NICHT GUT!!!!! after relinking nei1->link and saved[0]->link do not match!!!"<<endl;
					if(nei2_new->link_neighbors[part] == saved_nei[1]->link_neighbors[part])
						cout<<"GREAT! after relinking nei2->link and saved[1]->link match!!!"<<endl;
					else
						cout<<"NICHT GUT!!!!! after relinking nei2->link and saved[1]->link do not match!!!"<<endl;
*/

					if(nei1_new->link_neighbors[part]){
						nei1_new->link_neighbors[part]->length += old_brlen*part_info[part].part_rate;
						nei2_new->link_neighbors[part]->length += old_brlen*part_info[part].part_rate;

//					cout<<"AFTER RELINKING BACK"<<endl;
//					cout<<"ne1_new -> "<<nei1_new->link_neighbors[part]->node->name<<","<<nei1_new->link_neighbors[part]->node->id<<endl;
//					cout<<"ne2_new -> "<<nei2_new->link_neighbors[part]->node->name<<","<<nei2_new->link_neighbors[part]->node->id<<endl;
					}
				}
			//}
//			if(part == ntrees-1){
//				cout<<endl<<endl<<"SWAP BACK, RELINK BACK"<<endl<<endl;
//				printMapInfo();
//			}
		}
		//mapBranchLen();
		checkBranchLen();
	} // end of for(cnt)


//=============================================================================================================================================================
// 							Restoring
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
	FOR_NEIGHBOR(node1, node2, it)
		(*it)->length = (*it)->node->findNeighbor(node1)->length;
	FOR_NEIGHBOR(node2, node1, it)
		(*it)->length = (*it)->node->findNeighbor(node2)->length;

	// Restoring information for SubTrees ------------------------------------------------------------------------------------------------------------
		for(part = 0; part < ntrees; part++){
			if(is_nni[part]){ //for RELINK everything is restored in for(cnt)
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

				nei1_new = (SuperNeighbor*) node1->findNeighbor(node2);
				nei2_new = (SuperNeighbor*) node2->findNeighbor(node1);

			}
		}
		//mapBranchLen();
		checkBranchLen();
//------------------------------------------------------------------------------------------------------------------------------------------------
//		cout<<endl<<endl<<"THE END OF SWAPNNI func"<<endl<<endl;
//		printMapInfo();

		//double new_brlen_last = node1->findNeighbor(node2)->length;
		//if(old_brlen == new_brlen_last)
		//	cout<<"SwapNNI:before and after trial branch length of node1-node2 is the same"<<endl;

	// Minh: delete memory
	delete [] saved_cur_score;

	return cur_score;
}

void PhyloSuperTreePlen::linkCheck(int part,Node* node, Node* dad, PhyloNeighbor* saved_link_dad_nei){
	NeighborVec::iterator it;
	SuperNeighbor *dad_nei = (SuperNeighbor*)dad->findNeighbor(node);
	SuperNeighbor *node_nei = (SuperNeighbor*)node->findNeighbor(dad);
	FOR_NEIGHBOR(node, dad, it){
		//if(((SuperNeighbor*)(*it))->link_neighbors[part]){
			if(((SuperNeighbor*)(*it))->link_neighbors[part] == saved_link_dad_nei){
				((SuperNeighbor*)(*it))->link_neighbors[part] = dad_nei->link_neighbors[part];
				((SuperNeighbor*)((*it)->node->findNeighbor(node)))->link_neighbors[part] = node_nei->link_neighbors[part];
				linkCheck(part, (*it)->node, node, saved_link_dad_nei);
			}
		}

}

void PhyloSuperTreePlen::linkCheckRe(int part,Node* node, Node* dad, PhyloNeighbor* saved_link_dad_nei,PhyloNeighbor* saved_link_node_nei){
	NeighborVec::iterator it;
//	cout<<"linkCheckRE:"<<endl;
//	cout<<"node = "<<node->name<<","<<node->id<<endl;
//	cout<<"dad  = "<<dad->name<<","<<dad->id<<endl;
	FOR_NEIGHBOR(node, dad, it){
//		cout<<"nei->node = "<<(*it)->node->name<<","<<(*it)->node->id<<endl;
//		if(((SuperNeighbor*)(*it))->link_neighbors[part]){
//			cout<<" has linked_nei"<<endl;
			if(((SuperNeighbor*)(*it))->link_neighbors[part] == ((SuperNeighbor*)dad->findNeighbor(node))->link_neighbors[part]){
			//	cout<<"---------------------------------------"<<endl;
				linkCheckRe(part, (*it)->node, node, saved_link_dad_nei, saved_link_node_nei);
			//	cout<<"---------------------------------------"<<endl;
			//	cout<<" has been linked to OLD PhyloNei link_nei"<<endl;
				((SuperNeighbor*)(*it))->link_neighbors[part] = saved_link_dad_nei;
				((SuperNeighbor*)((*it)->node->findNeighbor(node)))->link_neighbors[part] = saved_link_node_nei;

			} else {
			//	cout<<"(*it)->link_nei,            br_id = "<<((SuperNeighbor*)(*it))->link_neighbors[part]->id<<endl;
			//	cout<<"(dad->find(node))->link_nei,br_id = "<<((SuperNeighbor*)dad->findNeighbor(node))->link_neighbors[part]->id<<endl;
			}
		//}
	}
}
void PhyloSuperTreePlen::restoreAllBranLen(PhyloNode *node, PhyloNode *dad) {
	IQTree::restoreAllBranLen(node,dad);
	mapTrees();
}

bool PhyloSuperTreePlen::checkBranchLen(){

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
	delete [] checkVAL;

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
				checkVAL[((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]))->link_neighbors[part]->id] += nodes1[i]->findNeighbor(nodes2[i])->length * part_info[part].part_rate;
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

/* Some case with negative branch:
 * seed 291475,10x1000.03
 * seed 981244,10x1000.05
 */
void PhyloSuperTreePlen::computeBranchLengths()
{
	}
