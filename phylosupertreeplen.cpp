/*
 * phylosupertreeplen.cpp
 *
 *  Created on: Aug 5, 2013
 *      Author: olga
 */

#include "phylosupertreeplen.h"
#include "superalignmentpairwise.h"
#include <string.h>

PhyloSuperTreePlen::PhyloSuperTreePlen()
: PhyloSuperTree()
{}

PhyloSuperTreePlen::PhyloSuperTreePlen(Params &params)
: PhyloSuperTree(params)
{}

PhyloSuperTreePlen::PhyloSuperTreePlen(SuperAlignment *alignment, PhyloSuperTree *super_tree)
: PhyloSuperTree(alignment,super_tree)
{}

PhyloSuperTreePlen::~PhyloSuperTreePlen()
:~PhyloSuperTree()
{}

// SuperAlignmentPairwisePlen ----------------------------------------------------------------------------------

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

// PartitionModelPlen ------------------------------------------------------------------------------------------

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
    double tree_lh = 0.0, cur_lh = 0.0, new_lh = 0.0;
    int ntrees = tree->size();
    double tol = 0.2;

	/*#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh)
	#endif*/

    tree_lh = tree->computeLikelihood();

    for(int i = 1; i < 100; i++){
    	tol = max(tol/2, epsilon);

    	for (int part = 0; part < ntrees; part++) {

    		// Subtree model parameters optimization
        	double model_lh = tree->at(part)->getModelFactory()->model->optimizeParameters();
        	double rate_lh = 0.0;
        	rate_lh = tree->at(part)->getModelFactory()->site_rate->optimizeParameters();
    	}

    	// Optimizing gene rate
    	cur_lh = optimizeGeneRate(tol);

    	// Optimizing branch lengths
    	int my_iter = min(5,i);
    	if(!fixed_len){
    		new_lh = tree->optimizeAllBranches(my_iter,tol);
    	}

    	new_lh = (new_lh != 0) ? new_lh : cur_lh;

    	if(fabs(new_lh-tree_lh) < epsilon)
    		break;

    	tree_lh = new_lh;
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
	for(int part = 0; part < tree->size(); part ++){
			if(tree->part_info[part].part_rate != x[part+1]){
				tree->at(part)->clearAllPartialLH();
			}
		}
	getVariables(x);
	if (tree->part_info[tree->size()-1].part_rate < 1e-4) return 1.0e+12;
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
}

double PhyloSuperTreePlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH) {

	current_it = node1->findNeighbor(node2);
	current_it_back = node2->findNeighbor(node1);

    double negative_lh = 0.0;
    double current_len = current_it->length;
    double ferror, optx;
	int ntrees = size();

	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh)
	#endif

    if (optimize_by_newton) { // Newton-Raphson method
    	// You should change tolerance!! start from small and then increase
            optx = minimizeNewton(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, negative_lh);
    } else
        // Brent method
        optx = minimizeOneDimen(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, &negative_lh, &ferror);
    if (current_len == optx) // if nothing changes, return
        return -negative_lh;

    current_it->length = optx;
    current_it_back->length = optx;

    // how to deal with this stuff?

	/*if (part_info[part].opt_score[nei1_part->id] == 0.0) {
	//	part_info[part].cur_brlen[nei1_part->id] = nei1_part->length;
	//	part_info[part].opt_score[nei1_part->id] = at(part)->optimizeOneBranch((PhyloNode*)nei1_part->node, (PhyloNode*)nei2_part->node, clearLH);
	//	part_info[part].opt_brlen[nei1_part->id] = nei1_part->length;
	//}
	//score = part_info[part].opt_score[nei1_part->id];


    //if (clearLH) {
    //    node1->clearReversePartialLh(node2);
    //    node2->clearReversePartialLh(node1);
    //}*/
    return -negative_lh;
}

double PhyloSuperTreePlen::computeFunction(double value) {

	double tree_lh;
	int ntrees = size();

	double lambda = value-current_it->length;
	current_it->length = value;
    current_it_back->length = value;

	SuperNeighbor *nei1 = ((SuperNeighbor*)current_it_back->node->findNeighbor(current_it->node));
	SuperNeighbor *nei2 = ((SuperNeighbor*)current_it->node->findNeighbor(current_it_back->node));
	assert(nei1 && nei2);

	for (int part = 0; part < ntrees; part++) {
			PhyloNeighbor *nei1_part = nei1->link_neighbors[part];
			PhyloNeighbor *nei2_part = nei2->link_neighbors[part];
			double score;
			if (part_info[part].cur_score == 0.0)
				part_info[part].cur_score = at(part)->computeLikelihood();
			if (nei1_part && nei2_part) {
				nei1_part->length += lambda*part_info[part].part_rate;
				nei2_part->length += lambda*part_info[part].part_rate;
				part_info[part].cur_score = at(part)->computeLikelihood();
			} else {
				score = part_info[part].cur_score;
			}
			tree_lh += score;
		}
    return -tree_lh;
}


// this one has to be changed *********************************

double PhyloSuperTreePlen::computeFuncDerv(double value, double &df, double &ddf) {
    current_it->length = value;
    current_it_back->length = value;
    double lh;
    if (params->fast_branch_opt) {
        // Pre-compute Theta vector
        if (!theta_computed) {
            computeTheta(current_it, (PhyloNode*) current_it_back->node);
            theta_computed = true;
        }
        lh = -computeLikelihoodDervFast(current_it, (PhyloNode*) current_it_back->node, df, ddf);
    } else {
        lh = -computeLikelihoodDerv(current_it, (PhyloNode*) current_it_back->node, df, ddf);
    }
    df = -df;
    ddf = -ddf;

    return lh;
}


