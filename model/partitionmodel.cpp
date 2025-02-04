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
#include "model/modelmarkov.h"

PartitionModel::PartitionModel()
        : ModelFactory()
{
	linked_alpha = -1.0;
    opt_gamma_invar = false;
    partLike = NULL;
}

PartitionModel::PartitionModel(Params &params, PhyloSuperTree *tree, ModelsBlock *models_block)
        : ModelFactory()
{
	store_trans_matrix = params.store_trans_matrix;
	is_storing = false;
	joint_optimize = params.optimize_model_rate_joint;
	fused_mix_rate = false;
    linked_alpha = -1.0;
    opt_gamma_invar = false;

	// create dummy model
	model = new ModelSubst(tree->aln->num_states);
	site_rate = new RateHeterogeneity();
	site_rate->setTree(tree);
    
    // create an array to store the log-likelihood for each partition
    partLike = new double[tree->size()];

//    string model_name = params.model_name;
    PhyloSuperTree::iterator it;
    int part;
    if (params.link_alpha) {
        params.gamma_shape = fabs(params.gamma_shape);
        linked_alpha = params.gamma_shape;
    }
    double init_by_divmat = false;
    if (params.model_name_init && strcmp(params.model_name_init, "DIVMAT") == 0) {
        init_by_divmat = true;
        params.model_name_init = NULL;
    }
    for (it = tree->begin(), part = 0; it != tree->end(); it++, part++) {
        ASSERT(!((*it)->getModelFactory()));
        string model_name = (*it)->aln->model_name;
        if (model_name == "") // if empty, take model name from command option
        	model_name = params.model_name;
        if ((*it)->isTreeMix()) {
            ((IQTreeMixHmm*)(*it))->initializeModel(params, model_name, models_block);
        } else {
            (*it)->setModelFactory(new ModelFactory(params, model_name, (*it), models_block));
            (*it)->setModel((*it)->getModelFactory()->model);
            (*it)->setRate((*it)->getModelFactory()->site_rate);

            // link models between partitions
            if (params.link_model) {
                (*it)->getModel()->fixParameters(true);
                if (linked_models.find((*it)->getModel()->getName()) == linked_models.end()) {
                    linked_models[(*it)->getModel()->getName()] = (*it)->getModel();
                }
            } else if ((*it)->aln->getNSeq() < tree->aln->getNSeq() && params.partition_type != TOPO_UNLINKED &&
                       (*it)->getModel()->freq_type == FREQ_EMPIRICAL && (*it)->aln->seq_type != SEQ_CODON) {
                // modify state_freq to account for empty sequences
                (*it)->aln->computeStateFreq((*it)->getModel()->state_freq, (*it)->aln->getNSite() * (tree->aln->getNSeq() - (*it)->aln->getNSeq()));
                (*it)->getModel()->decomposeRateMatrix();
            }
        }
        
        //string taxa_set = ((SuperAlignment*)tree->aln)->getPattern(part);
        //(*it)->copyTree(tree, taxa_set);
        //(*it)->drawTree(cout);
    }
    if (init_by_divmat) {
        ASSERT(0 && "init_by_div_mat not working");
        int nstates = linked_models.begin()->second->num_states;
        double *pair_freq = new double[nstates * nstates];
        double *state_freq = new double[nstates];
        tree->aln->computeDivergenceMatrix(pair_freq, state_freq);
        /*
        MatrixXd divmat = Map<Matrix<double,Dynamic, Dynamic, RowMajor> > (pair_freq, nstates, nstates);
        cout << "DivMat: " << endl << divmat << endl;
        auto pi = Map<VectorXd>(state_freq, nstates);
        MatrixXd Q = (pi.asDiagonal() * divmat).log();
        cout << "Q: " << endl << Q << endl;
        cout << "rowsum: " << Q.rowwise().sum() << endl;
        Map<Matrix<double,Dynamic, Dynamic, RowMajor> >(pair_freq, nstates, nstates) = Q;
         */
        ((ModelMarkov*)linked_models.begin()->second)->setFullRateMatrix(pair_freq, state_freq);
        ((ModelMarkov*)linked_models.begin()->second)->decomposeRateMatrix();
        delete [] state_freq;
        delete [] pair_freq;

    } else
    for (auto mit = linked_models.begin(); mit != linked_models.end(); mit++) {
        PhyloSuperTree *stree = (PhyloSuperTree*)site_rate->phylo_tree;
        if (mit->second->freq_type != FREQ_ESTIMATE && mit->second->freq_type != FREQ_EMPIRICAL)
            continue;
        // count state occurrences
        size_t *sum_state_counts = NULL;
        int num_parts = 0;
        for (it = stree->begin(); it != stree->end(); it++) {
            if ((*it)->getModel()->getName() == mit->second->getName()) {
                num_parts++;
                if ((*it)->aln->seq_type == SEQ_CODON)
                    outError("Linking codon models not supported");
                if ((*it)->aln->seq_type == SEQ_POMO)
                    outError("Linking POMO models not supported");
                size_t state_counts[(*it)->aln->STATE_UNKNOWN+1];
                size_t unknown_states = 0;
                if( params.partition_type != TOPO_UNLINKED)
                    unknown_states = (*it)->aln->getNSite() * (tree->aln->getNSeq() - (*it)->aln->getNSeq());
                (*it)->aln->countStates(state_counts, unknown_states);
                if (!sum_state_counts) {
                    sum_state_counts = new size_t[(*it)->aln->STATE_UNKNOWN+1];
                    memset(sum_state_counts, 0, sizeof(size_t)*((*it)->aln->STATE_UNKNOWN+1));
                }
                for (int state = 0; state <= (*it)->aln->STATE_UNKNOWN; ++state) {
                    sum_state_counts[state] += state_counts[state];
                }
            }
        }
        cout << "Linking " << mit->first << " model across " << num_parts << " partitions" << endl;
        int nstates = mit->second->num_states;
        double sum_state_freq[nstates];
        // convert counts to frequencies
        for (it = stree->begin(); it != stree->end(); it++) {
            if ((*it)->getModel()->getName() == mit->second->getName()) {
                (*it)->aln->convertCountToFreq(sum_state_counts, sum_state_freq);
                break;
            }
        }

        cout << "Mean state frequencies:";
        int prec = cout.precision(8);
        for (int state = 0; state < mit->second->num_states; state++)
            cout << " " << sum_state_freq[state];
        cout << endl;
        cout.precision(prec);

        for (it = stree->begin(); it != stree->end(); it++)
            if ((*it)->getModel()->getName() == mit->second->getName()) {
                ((ModelMarkov*)(*it)->getModel())->adaptStateFrequency(sum_state_freq);
                (*it)->getModel()->decomposeRateMatrix();
            }
        delete [] sum_state_counts;
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
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        checkpoint->startStruct(it->first);
        bool fixed = it->second->fixParameters(false);
        it->second->saveCheckpoint();
        it->second->fixParameters(fixed);
        checkpoint->endStruct();
    }
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

    // restore linked models
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        checkpoint->startStruct(it->first);
        for (auto tit = tree->begin(); tit != tree->end(); tit++)
            if ((*tit)->getModel()->getName() == it->first) {
                bool fixed = (*tit)->getModel()->fixParameters(false);
                (*tit)->getModel()->restoreCheckpoint();
                (*tit)->getModel()->fixParameters(fixed);
            }
        checkpoint->endStruct();
    }
    
    endCheckpoint();
}

bool PartitionModel::isReversible() {
    // check that all sub-models must be reversible
    PhyloSuperTree *super_tree = (PhyloSuperTree*)site_rate->getTree();
    for (auto tree : *super_tree) {
        if (!tree->getModelFactory()->isReversible())
            return false; // at least one sub-model is non-reversible
    }
    return true;
}

int PartitionModel::getNParameters(int brlen_type) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
	int df = 0;
    for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++) {
    	df += (*it)->getModelFactory()->getNParameters(brlen_type);
    }
    if (linked_alpha > 0)
        df ++;
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        bool fixed = it->second->fixParameters(false);
        df += it->second->getNDim() + it->second->getNDimFreq();
        it->second->fixParameters(fixed);
    }
    return df;
}

double PartitionModel::computeFunction(double shape) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double res = 0.0;
    int ntrees = tree->size();
    linked_alpha = shape;
    if (tree->part_order.empty()) tree->computePartitionOrder();
#ifdef _OPENMP
#pragma omp parallel for reduction(+: res) schedule(dynamic) if(tree->num_threads > 1)
#endif
    for (int j = 0; j < ntrees; j++) {
        int i = tree->part_order[j];
        if (tree->at(i)->getRate()->isGammaRate())
            res += tree->at(i)->getRate()->computeFunction(shape);
    }
    if (res == 0.0) {
        outError("No partition has Gamma rate heterogeneity!");
    }
	return res;
}

double PartitionModel::optimizeLinkedAlpha(bool write_info, double gradient_epsilon) {
    if (write_info) {
        cout << "Optimizing linked gamma shape..." << endl;
    }
	double negative_lh;
	double current_shape = linked_alpha;
	double ferror, optx;
	optx = minimizeOneDimen(site_rate->getTree()->params->min_gamma_shape, current_shape, MAX_GAMMA_SHAPE, max(gradient_epsilon, TOL_GAMMA_SHAPE), &negative_lh, &ferror);
    double tree_lh = site_rate->getTree()->computeLikelihood();
    if (write_info) {
        cout << "Linked alpha across partitions: " << linked_alpha << endl;
        cout << "Linked alpha log-likelihood: " << tree_lh << endl;
    }
	return tree_lh;
    
}

int PartitionModel::getNDim() {
    return model->getNDim();
}

double PartitionModel::targetFunk(double x[]) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    
    double res = 0;
    int ntrees = tree->size();
    if (tree->part_order.empty()) tree->computePartitionOrder();
    string modelname = model->getName();
    
    // reset the array to store the log-likelihoods for each partition
    memset(partLike, 0, sizeof(double) * ntrees);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(tree->num_threads > 1)
    #endif
    for (int j = 0; j < ntrees; j++) {
        
        int i = tree->part_order[j];
        ModelSubst *part_model = tree->at(i)->getModel();

        if (part_model->getName() != modelname)
            continue;

        bool fixed = part_model->fixParameters(false);
        double ans = part_model->targetFunk(x);
        part_model->fixParameters(fixed);
        partLike[j] = ans;
    }
    
    // calculate the sum
    // different order of the computation due to
    // openmp will not affect the value
    res = 0.0;
    for (int j = 0; j < ntrees; j++) {
        res += partLike[j];
    }
    
    if (res == 0.0)
        outError("No partition has model ", model->getName());
    return res;
}

double PartitionModel::computeMixLh(string &warning) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    int ntrees = tree->size();

    // go through the number of sites of each partition to compute the class weights
    vector<double> weight_array, log_weight_array;
    double sum_sites = 0.0;

    // compute "class weights"
    //if (tree->part_order.empty()) tree->computePartitionOrder();
    for (int j = 0; j < ntrees; j++) {
        //int i = tree->part_order[j];
        int part_sites = tree->at(j)->getAlnNSite();
        weight_array.push_back(part_sites);
        sum_sites += part_sites;
    }
    for (int j = 0; j < ntrees; j++) {
        weight_array[j] /= sum_sites;
    }
    for (int j = 0; j < ntrees; j++) {
        log_weight_array.push_back(log(weight_array[j]));
    }

    //get sets of taxa for each partition tree in advance
    vector<StrVector> t_seqs_vec_array;
    vector<set<string> > t_seqs_set_array;
    t_seqs_vec_array.resize(ntrees);
    t_seqs_set_array.resize(ntrees);

    for (int j = 0; j < ntrees; j++) {
        PhyloTree *t = tree->at(j);
        t->getTaxaName(t_seqs_vec_array[j]);
        t_seqs_set_array[j].insert(t_seqs_vec_array[j].begin(), t_seqs_vec_array[j].end());
    }

    // compute the mixture-based log-likelihood
    double mix_lh = 0.0;
    bool too_much_missing = false;

    for (int j = 0; j < ntrees && (!too_much_missing); j++) {
        //int i = tree->part_order[j];
        Alignment *tree1_aln = tree->at(j)->aln;
        int tree1_nsite = tree1_aln->getNSite();
        StrVector tree1_seqs = t_seqs_vec_array[j];

        // get the site-log-likelihood the the partition under each tree and the corresponding model
        double *lh_array = new double [ntrees*tree1_nsite];

#ifdef _OPENMP
#pragma omp parallel for if(tree->num_threads > 1)
#endif
        for (int k = 0; k < ntrees ; k++) {
            PhyloTree *tree2 = tree->at(k);

            // get the intersection of tree1_aln and tree2.
            set<string> inter_seqs;
            IntVector inter_seqs_id, missing_seqs_id;

            set<string> tree2_seqs_set = t_seqs_set_array[k];
            StrVector tree2_seqs = t_seqs_vec_array[k];
            for (string seq_name : tree1_seqs) {
                int seq_id = tree1_aln->getSeqID(seq_name);
                if (tree2_seqs_set.find(seq_name) != tree2_seqs_set.end()) {
                    inter_seqs.insert(seq_name);
                    inter_seqs_id.push_back(seq_id);
                } else {
                    missing_seqs_id.push_back(seq_id);
                }
            }

            // if the subset has less than 3 sequences, don't compute m-log-likelihood
            if (inter_seqs_id.size() < 3) {

#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    too_much_missing = true;
                    string tree1_name = tree1_aln->name;
                    string tree2_name = tree2->aln->name;
                    int ntaxa = tree->getNumTaxa();
                    warning =
                            "NOTE: Mixture-based log-likelihood conversion is skipped due to too much missing data: at least one of the partitions " +
                            tree1_name + " and " + tree2_name + " show missing data in " +
                            to_string(ntaxa - inter_seqs_id.size()) + " sequences.";
                }

            } else {
                // subset tree1_aln
                Alignment *sub_tree1_aln = NULL;
                if (tree1_seqs.size() != inter_seqs_id.size()) {
                    sub_tree1_aln = new Alignment();
                    sub_tree1_aln->extractSubAlignment(tree1_aln, inter_seqs_id, 0);
                } else {
                    sub_tree1_aln = tree1_aln;
                }

                // subset tree2
                PhyloTree *sub_tree2 = NULL;
                string inter_seqs_set (tree2_seqs.size(), 0);
                for (int l = 0; l < tree2_seqs.size(); l++) {
                    if (inter_seqs.find(tree2_seqs[l]) != inter_seqs.end()) {
                        inter_seqs_set[l] = 1;
                    }
                }

                sub_tree2 = new PhyloTree();
                if (tree2_seqs.size() != inter_seqs_id.size()) {
                    sub_tree2->copyTree(tree2, inter_seqs_set);
                } else {
                    sub_tree2->copyTree(tree2);
                }

                // link sub_tree2 and sub_tree1_aln
                sub_tree2->setAlignment(sub_tree1_aln);
                sub_tree2->setRootNode(tree2->params->root);
                sub_tree2->setModelFactory(tree2->getModelFactory());
                //sub_tree2->setModel(tree2->getModel());
                //sub_tree2->setRate(tree2->getRate());
                sub_tree2->setParams(tree2->params);
                sub_tree2->optimize_by_newton = tree2->params->optimize_by_newton;
                sub_tree2->setLikelihoodKernel(tree2->params->SSE);
                sub_tree2->setNumThreads(tree2->num_threads);
                sub_tree2->ensureNumberOfThreadsIsSet(nullptr);

                sub_tree2->initializeAllPartialLh();

                //compute pattern log-likelihood
                int nptn = sub_tree1_aln->getNPattern();
                double *ptn_lh_array = new double [nptn];
                sub_tree2->computeLikelihood(ptn_lh_array);

                //compute log state frequencies
                int n_states = sub_tree2->getModel()->num_states;
                double *state_freq = new double[n_states];
                sub_tree2->getModel()->getStateFrequency(state_freq);

                vector<double> log_state_freq(n_states);
                for (int n = 0; n < n_states; n++) {
                    log_state_freq[n] = log(state_freq[n]);
                }

                //compute site log-likelihood
                if (tree1_seqs.size() != inter_seqs_id.size()) { //for sequences only appears in tree1, calculate the state frequencies based on the model in tree2
                    for (int l = 0; l < tree1_nsite; l++) {
                        int ptn_id = sub_tree1_aln->getPatternID(l);
                        double site_lh = ptn_lh_array[ptn_id];
                        Pattern p = tree1_aln->at(tree1_aln->getPatternID(l));

                        for (int missing_id: missing_seqs_id) {
                            int char_id = p[missing_id];
                            if (char_id  < n_states) {
                                site_lh += log_state_freq[char_id];
                            } else {
                                // compute ambiguous frequencies
                                int cstate = char_id-n_states+1;
                                double amb_freq = 0;
                                for (int m = 0; m < n_states; m++) {
                                    if ((cstate) & (1 << m)) {
                                        amb_freq += state_freq[m];
                                    }
                                }
                                site_lh += log(amb_freq);
                            }
                        }
                        lh_array[tree1_nsite * k + l] = site_lh;
                    }
                } else {
                    for (int l = 0; l < tree1_nsite; l++) {
                        int ptn_id = sub_tree1_aln->getPatternID(l);
                        lh_array[tree1_nsite * k + l] = ptn_lh_array[ptn_id];
                    }
                }

                //release memory
                if (tree1_seqs.size() != inter_seqs_id.size()) {
                    delete sub_tree1_aln;
                }
                sub_tree2->setModelFactory(NULL);
                sub_tree2->aln = NULL;
                delete sub_tree2;
                delete[] ptn_lh_array;
                delete[] state_freq;
            }

        }

        // compute partition log-likelihood from sites
        if (!too_much_missing){
            double mix_lh_partition = 0.0;
            for (int l = 0; l < tree1_nsite; l++) {
                double weighted_lh, max_lh, mix_lh_site;

                //int ptn_freq = tree1_aln->at(l).frequency;

                for (int k = 0; k < ntrees; k++) {
                    weighted_lh = log_weight_array[k]+lh_array[tree1_nsite*k+l];
                    if (k == 0) {
                        max_lh = weighted_lh;
                    } else if (weighted_lh > max_lh) {
                        max_lh = weighted_lh;
                    }
                }

                double mix_lh_site_original = 0.0;
                for (int k = 0; k < ntrees; k++) {
                    mix_lh_site_original += exp(log_weight_array[k]+lh_array[tree1_nsite*k+l]-max_lh);
                }
                mix_lh_site = max_lh + log(mix_lh_site_original);
                mix_lh_partition += mix_lh_site;
            }
            mix_lh += mix_lh_partition;
        }
        delete[] lh_array; //release array memery
    }

    if (too_much_missing){
        return 1.0;
    } else {
        return mix_lh;
    }
}

void PartitionModel::setVariables(double *variables) {
    model->setVariables(variables);
}

bool PartitionModel::getVariables(double *variables) {
    bool changed = false;
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    for (auto it = tree->begin(); it != tree->end(); it++)
        if ((*it)->getModel()->getName() == model->getName())
            changed |= (*it)->getModel()->getVariables(variables);
    return changed;
}

void PartitionModel::scaleStateFreq(bool sum_one) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    for (auto it = tree->begin(); it != tree->end(); it++)
        if ((*it)->getModel()->getName() == model->getName())
            ((ModelMarkov*)(*it)->getModel())->scaleStateFreq(sum_one);
}

double PartitionModel::optimizeLinkedModel(bool write_info, double gradient_epsilon) {
    int ndim = getNDim();
    
    // return if nothing to be optimized
    if (ndim == 0) return 0.0;
    
    if (write_info)
        cout << "Optimizing linked " << model->getName() << " parameters across all partitions (" << ndim << " free parameters)" << endl;
    
    if (verbose_mode >= VB_MAX)
        cout << "Optimizing " << model->name << " model parameters..." << endl;
    
    //if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);
    
    double *variables = new double[ndim+1]; // used for BFGS numerical recipes
    double *variables2 = new double[ndim+1]; // used for L-BFGS-B
    double *upper_bound = new double[ndim+1];
    double *lower_bound = new double[ndim+1];
    bool *bound_check = new bool[ndim+1];
    double score;
    
    
    // by BFGS algorithm
    setVariables(variables);
    setVariables(variables2);
    ((ModelMarkov*)model)->setBounds(lower_bound, upper_bound, bound_check);
    // expand the bound for linked model
//    for (int i = 1; i <= ndim; i++) {
//        lower_bound[i] = MIN_RATE*0.2;
//        upper_bound[i] = MAX_RATE*2.0;
//    }

//    if (Params::getInstance().optimize_alg.find("BFGS-B") == string::npos)
        score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));
//    else
//        score = -L_BFGS_B(ndim, variables+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_RATE));

    bool changed = getVariables(variables);

    /* 2019-09-05: REMOVED due to numerical issue (NAN) with L-BFGS-B
    // 2017-12-06: more robust optimization using 2 different routines
    // when estimates are at boundary
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_RATE));
    bool changed = getVariables(variables);
    
    if (model->isUnstableParameters()) {
        // parameters at boundary, restart with L-BFGS-B with parameters2
        double score2 = -L_BFGS_B(ndim, variables2+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_RATE));
        if (score2 > score+0.1) {
            if (verbose_mode >= VB_MED)
                cout << "NICE: L-BFGS-B found better parameters with LnL=" << score2 << " than BFGS LnL=" << score << endl;
            changed = getVariables(variables2);
            score = score2;
        } else {
            // otherwise, revert what BFGS found
            changed = getVariables(variables);
        }
    }
    */
    
    // BQM 2015-09-07: normalize state_freq
    if (model->isReversible() && model->freq_type == FREQ_ESTIMATE) {
        scaleStateFreq(true);
        changed = true;
    }
    if (changed) {
        PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
        for (auto it = tree->begin(); it != tree->end(); it++)
            if ((*it)->getModel()->getName() == model->getName())
                (*it)->getModel()->decomposeRateMatrix();
        site_rate->phylo_tree->clearAllPartialLH();
        score = site_rate->phylo_tree->computeLikelihood();
    }
    
    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables2;
    delete [] variables;
    
    /*
    if (write_info) {
        cout << "Linked-model log-likelihood: " << score << endl;
    }
    */

    return score;
}

double PartitionModel::optimizeLinkedModels(bool write_info, double gradient_epsilon) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double tree_lh;
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        ModelSubst *saved_model = model;
        model = it->second;
        PhyloSuperTree::iterator part_tree;
        // un-fix model parameters
        for (part_tree = tree->begin(); part_tree != tree->end(); part_tree++)
            if ((*part_tree)->getModel()->getName() == model->getName())
                (*part_tree)->getModel()->fixParameters(false);
        
        // main call to optimize linked model parameters
        tree_lh = optimizeLinkedModel(write_info, gradient_epsilon);
        
        // fix model parameters again
        for (part_tree = tree->begin(); part_tree != tree->end(); part_tree++)
            if ((*part_tree)->getModel()->getName() == model->getName())
                (*part_tree)->getModel()->fixParameters(true);
        
        saveCheckpoint();
        getCheckpoint()->dump();
        model = saved_model;
    }

    return site_rate->phylo_tree->computeLikelihood();
}

void PartitionModel::reportLinkedModel(ostream &out) {
    if (linked_alpha > 0.0)
        out << "Linked alpha across partitions: " << linked_alpha << endl;
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        out << "Linked model " << it->first << ": " << endl;
        it->second->report(out);
    }
}

bool PartitionModel::isLinkedModel() {
    return Params::getInstance().link_alpha || (linked_models.size()>0);
}

double PartitionModel::optimizeParameters(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
    double prev_tree_lh = -DBL_MAX, tree_lh = 0.0;
    int ntrees = tree->size();

    for (int step = 0; step < Params::getInstance().model_opt_steps; step++) {
        tree_lh = 0.0;
        if (tree->part_order.empty()) tree->computePartitionOrder();
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tree_lh) schedule(dynamic) if(tree->num_threads > 1)
#endif
        for (int i = 0; i < ntrees; i++) {
            int part = tree->part_order[i];
            double score;
            if (tree->at(part)->isTreeMix()) {
                ((IQTreeMixHmm*)tree->at(part))->optimizeModelParameters(write_info && verbose_mode >= VB_MED, logl_epsilon);
                score = ((IQTreeMixHmm*)tree->at(part))->getCurScore();
            } else {
                if (opt_gamma_invar)
                    score = tree->at(part)->getModelFactory()->optimizeParametersGammaInvar(fixed_len,
                                                                                            write_info && verbose_mode >= VB_MED,
                                                                                            logl_epsilon/min(ntrees,10), gradient_epsilon/min(ntrees,10));
                else
                    score = tree->at(part)->getModelFactory()->optimizeParameters(fixed_len,
                                                                                  write_info && verbose_mode >= VB_MED,
                                                                                  logl_epsilon/min(ntrees,10), gradient_epsilon/min(ntrees,10));
            }
            tree_lh += score;
            if (write_info)
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                cout << "Partition " << tree->at(part)->aln->name
                << " / Model: " << tree->at(part)->getModelName()
                << " / df: " << tree->at(part)->getModelFactory()->getNParameters(fixed_len)
                << " / LogL: " << score << endl;
            }
        }
        //return ModelFactory::optimizeParameters(fixed_len, write_info);
        
        if (!isLinkedModel())
            break;

        if (verbose_mode >= VB_MED || write_info)
            cout << step+1 << ". Log-likelihood: " << tree_lh << endl;

        // optimize linked alpha
        if (tree->params->link_alpha) {
            tree_lh = optimizeLinkedAlpha(write_info, gradient_epsilon);
        }

        // optimize linked models
        if (!linked_models.empty()) {
            double new_tree_lh = optimizeLinkedModels(write_info, gradient_epsilon);
            ASSERT(new_tree_lh > tree_lh - 0.1);
            tree_lh = new_tree_lh;
        }
        
        if (verbose_mode >= VB_MED || write_info)
            cout << step+1 << ". Log-likelihood: " << tree_lh << endl;

        if (tree_lh-logl_epsilon*10 < prev_tree_lh)
            break;
        prev_tree_lh = tree_lh;
    }
    
    if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << tree_lh << endl;
    // write linked_models
    if (verbose_mode <= VB_MIN && write_info) {
        for (auto it = linked_models.begin(); it != linked_models.end(); it++)
            it->second->writeInfo(cout);
    }
    return tree_lh;
}


double PartitionModel::optimizeParametersGammaInvar(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    opt_gamma_invar = true;
    double tree_lh = optimizeParameters(fixed_len, write_info, logl_epsilon, gradient_epsilon);
    opt_gamma_invar = false;
    return tree_lh;
}


PartitionModel::~PartitionModel()
{
    if (partLike != NULL)
        delete[] partLike;
}

bool PartitionModel::isUnstableParameters() {
    PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();

	for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++)
		if ((*it)->getModelFactory()->isUnstableParameters()) {
			return true;
		}
	return false;
}
