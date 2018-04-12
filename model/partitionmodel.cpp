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
#include "partitionmodel.h"
#include "alignment/superalignment.h"
#include "model/rategamma.h"

PartitionModel::PartitionModel()
        : ModelFactory()
{
	linked_alpha = -1.0;
}

PartitionModel::PartitionModel(Params &params, PhyloSuperTree *tree, ModelsBlock *models_block)
        : ModelFactory()
{
	store_trans_matrix = params.store_trans_matrix;
	is_storing = false;
	joint_optimize = params.optimize_model_rate_joint;
	fused_mix_rate = false;
    linked_alpha = -1.0;

	// create dummy model
	model = new ModelSubst(tree->aln->num_states);
	site_rate = new RateHeterogeneity();
	site_rate->setTree(tree);

//    string model_name = params.model_name;
    PhyloSuperTree::iterator it;
    int part;
    if (params.link_alpha) {
        params.gamma_shape = fabs(params.gamma_shape);
        linked_alpha = params.gamma_shape;
    }
    for (it = tree->begin(), part = 0; it != tree->end(); it++, part++) {
        ASSERT(!((*it)->getModelFactory()));
        string model_name = (*it)->aln->model_name;
        if (model_name == "") // if empty, take model name from command option
        	model_name = params.model_name;
        (*it)->setModelFactory(new ModelFactory(params, model_name, (*it), models_block));
        (*it)->setModel((*it)->getModelFactory()->model);
        (*it)->setRate((*it)->getModelFactory()->site_rate);
//        params.model_name = model_name;
        if ((*it)->aln->getNSeq() < tree->aln->getNSeq() && (*it)->getModel()->freq_type == FREQ_EMPIRICAL && (*it)->aln->seq_type != SEQ_CODON) {
        	// modify state_freq to account for empty sequences
        	(*it)->aln->computeStateFreq((*it)->getModel()->state_freq, (*it)->aln->getNSite() * (tree->aln->getNSeq() - (*it)->aln->getNSeq()));
        	(*it)->getModel()->decomposeRateMatrix();
        }
        
        //string taxa_set = ((SuperAlignment*)tree->aln)->getPattern(part);
        //(*it)->copyTree(tree, taxa_set);
        //(*it)->drawTree(cout);
    }
}

void PartitionModel::setCheckpoint(Checkpoint *checkpoint) {
	ModelFactory::setCheckpoint(checkpoint);
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++)
		(*it)->getModelFactory()->setCheckpoint(checkpoint);
}

void PartitionModel::startCheckpoint() {
    checkpoint->startStruct("PartitionModel");
}

void PartitionModel::saveCheckpoint() {
    startCheckpoint();
    CKP_SAVE(linked_alpha);
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    int part = 0;
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++, part++) {
        checkpoint->startStruct((*it)->aln->name);
        (*it)->getModelFactory()->saveCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();

    CheckpointFactory::saveCheckpoint();
}

void PartitionModel::restoreCheckpoint() {
    CheckpointFactory::restoreCheckpoint();
    startCheckpoint();
    CKP_RESTORE(linked_alpha);

    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    int part = 0;
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++, part++) {
        checkpoint->startStruct((*it)->aln->name);
        (*it)->getModelFactory()->restoreCheckpoint();
        checkpoint->endStruct();
    }

    endCheckpoint();
}

int PartitionModel::getNParameters(int brlen_type) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
	int df = 0;
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++) {
    	df += (*it)->getModelFactory()->getNParameters(brlen_type);
    }
    if (linked_alpha > 0)
        df ++;
    return df;
}

double PartitionModel::computeFunction(double shape) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double res = 0.0;
    linked_alpha = shape;
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++) 
        if ((*it)->getRate()->isGammaRate()) {
            res += (*it)->getRate()->computeFunction(shape);
        }
    if (res == 0.0)
        outError("No partition has Gamma rate heterogeneity!");
	return res;
}

double PartitionModel::optimizeLinkedAlpha(bool write_info, double gradient_epsilon) {
    if (write_info)
        cout << "Optimizing linked gamma shape..." << endl;
	double negative_lh;
	double current_shape = linked_alpha;
	double ferror, optx;
	optx = minimizeOneDimen(site_rate->getTree()->params->min_gamma_shape, current_shape, MAX_GAMMA_SHAPE, max(gradient_epsilon, TOL_GAMMA_SHAPE), &negative_lh, &ferror);
    if (write_info)
        cout << "Linked alpha across partitions: " << linked_alpha << endl;
	return site_rate->getTree()->computeLikelihood();
    
}

double PartitionModel::optimizeParameters(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double tree_lh = 0.0;
    int ntrees = tree->size();

    if (tree->part_order.empty()) tree->computePartitionOrder();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(tree->num_threads > 1)
	#endif
    for (int i = 0; i < ntrees; i++) {
        int part = tree->part_order[i];
    	if (write_info)
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
    		cout << "Optimizing " << tree->at(part)->getModelName() <<
        		" parameters for partition " << tree->at(part)->aln->name <<
        		" (" << tree->at(part)->getModelFactory()->getNParameters(fixed_len) << " free parameters)" << endl;
        }
        tree_lh += tree->at(part)->getModelFactory()->optimizeParameters(fixed_len, write_info && verbose_mode >= VB_MED, 
            logl_epsilon/min(ntrees,10), gradient_epsilon/min(ntrees,10));
    }
    //return ModelFactory::optimizeParameters(fixed_len, write_info);

    if (tree->params->link_alpha) {
        tree_lh = optimizeLinkedAlpha(write_info, gradient_epsilon);
    }
	if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << tree_lh << endl;
    return tree_lh;
}


double PartitionModel::optimizeParametersGammaInvar(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double tree_lh = 0.0;
    int ntrees = tree->size();

    if (tree->part_order.empty()) tree->computePartitionOrder();
	#ifdef _OPENMP
	#pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(tree->num_threads > 1)
	#endif
    for (int i = 0; i < ntrees; i++) {
        int part = tree->part_order[i];
    	if (write_info)
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
    		cout << "Optimizing " << tree->at(part)->getModelName() <<
        		" parameters for partition " << tree->at(part)->aln->name <<
        		" (" << tree->at(part)->getModelFactory()->getNParameters(fixed_len) << " free parameters)" << endl;
        }
        tree_lh += tree->at(part)->getModelFactory()->optimizeParametersGammaInvar(fixed_len, write_info && verbose_mode >= VB_MED, 
            logl_epsilon/min(ntrees,10), gradient_epsilon/min(ntrees,10));
    }
    //return ModelFactory::optimizeParameters(fixed_len, write_info);

    if (tree->params->link_alpha) {
        tree_lh = optimizeLinkedAlpha(write_info, gradient_epsilon);
    }
	if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << tree_lh << endl;
    return tree_lh;
}


PartitionModel::~PartitionModel()
{
}

bool PartitionModel::isUnstableParameters() {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();

	for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++)
		if ((*it)->getModelFactory()->isUnstableParameters()) {
			return true;
		}
	return false;
}
