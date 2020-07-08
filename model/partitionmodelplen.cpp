//
//  partitionmodelplen.cpp
//  iqtree
//
//  Created by Olga on 04/05/17.
//
//

#include <stdio.h>
#include "model/partitionmodelplen.h"
#include "utils/timeutil.h"
#include "model/modelmarkov.h"

/**********************************************************
 * class PartitionModelPlen
 **********************************************************/

//const double MIN_GENE_RATE = 0.001;
//const double MAX_GENE_RATE = 1000.0;
//const double TOL_GENE_RATE = 0.0001;

PartitionModelPlen::PartitionModelPlen()
: PartitionModel()
{
    //    optimizing_part = -1;
}

PartitionModelPlen::PartitionModelPlen(Params &params, PhyloSuperTreePlen *tree, ModelsBlock *models_block)
: PartitionModel(params, tree, models_block)
{
    //    optimizing_part = -1;
}

PartitionModelPlen::~PartitionModelPlen()
{
}

void PartitionModelPlen::startCheckpoint() {
    checkpoint->startStruct("PartitionModelPlen");
}

void PartitionModelPlen::saveCheckpoint() {
    startCheckpoint();
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    if (!tree->fixed_rates) {
        int nrates = tree->part_info.size();
        double *part_rates = new double[nrates];
        for (int i = 0; i < nrates; i++)
            part_rates[i] = tree->part_info[i].part_rate;
        CKP_ARRAY_SAVE(nrates, part_rates);
        delete [] part_rates;
    }
    endCheckpoint();
    PartitionModel::saveCheckpoint();
}

void PartitionModelPlen::restoreCheckpoint() {
    startCheckpoint();
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    if (!tree->fixed_rates) {
        int nrates = tree->part_info.size();
        double *part_rates = new double[nrates];
        if (CKP_ARRAY_RESTORE(nrates, part_rates)) {
            for (int i = 0; i < nrates; i++)
                tree->part_info[i].part_rate = part_rates[i];
            tree->mapTrees();
        }
        delete [] part_rates;
    }
    endCheckpoint();
    PartitionModel::restoreCheckpoint();
}


double PartitionModelPlen::optimizeParameters(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    double tree_lh = 0.0, cur_lh = 0.0;
    int ntrees = tree->size();
    
    
    //tree->initPartitionInfo(); // FOR OLGA: needed here

    unordered_map<string, bool> fixed_params;
    unordered_map<string, ModelSubst*>::iterator it;

    for(int part = 0; part < ntrees; part++){
        tree->part_info[part].cur_score = 0.0;
    }
    
//    if (fixed_len == BRLEN_OPTIMIZE) {
//        tree_lh = tree->optimizeAllBranches(1);
//    } else {
//        tree_lh = tree->computeLikelihood();
//    }
    tree_lh = tree->computeLikelihood();
    
    cout<<"Initial log-likelihood: "<<tree_lh<<endl;
    double begin_time = getRealTime();
    int i;
    for(i = 1; i < tree->params->num_param_iterations; i++){
        cur_lh = 0.0;
        if (tree->part_order.empty()) tree->computePartitionOrder();
#ifdef _OPENMP
#pragma omp parallel for reduction(+: cur_lh) schedule(dynamic) if(tree->num_threads > 1)
#endif
        for (int partid = 0; partid < ntrees; partid++) {
            int part = tree->part_order[partid];
            // Subtree model parameters optimization
            tree->part_info[part].cur_score = tree->at(part)->getModelFactory()->
                optimizeParametersOnly(i+1, gradient_epsilon/min(min(i,ntrees),10),
                                       tree->part_info[part].cur_score);
            if (tree->part_info[part].cur_score == 0.0)
                tree->part_info[part].cur_score = tree->at(part)->computeLikelihood();
            cur_lh += tree->part_info[part].cur_score;
            
            
            // normalize rates s.t. branch lengths are #subst per site
            double mean_rate = tree->at(part)->getRate()->rescaleRates();
            if (fabs(mean_rate-1.0) > 1e-6) {
                if (tree->fixed_rates) {
                    outError("Unsupported -spj. Please use proportion edge-linked partition model (-spp)");
                }
                tree->at(part)->scaleLength(mean_rate);
                tree->part_info[part].part_rate *= mean_rate;
            }
            
        }
        if (tree->params->link_alpha) {
            cur_lh = optimizeLinkedAlpha(write_info, gradient_epsilon);
        }

        // optimize linked models
        if (!linked_models.empty()) {
            double new_cur_lh = optimizeLinkedModels(write_info, gradient_epsilon);
            ASSERT(new_cur_lh > cur_lh - 0.1);
            cur_lh = new_cur_lh;
        }

        if (verbose_mode >= VB_MED)
            cout << "LnL after optimizing individual models: " << cur_lh << endl;
        if (cur_lh <= tree_lh - 1.0) {
            // more info for ASSERTION
            writeInfo(cout);
            tree->printTree(cout, WT_BR_LEN+WT_NEWLINE);
        }
        ASSERT(cur_lh > tree_lh - 1.0 && "individual model opt reduces LnL");
        
        tree->clearAllPartialLH();
        // Optimizing gene rate
        if(!tree->fixed_rates){
            cur_lh = optimizeGeneRate(gradient_epsilon);
            if (verbose_mode >= VB_MED) {
                cout << "LnL after optimizing partition-specific rates: " << cur_lh << endl;
                writeInfo(cout);
            }
            ASSERT(cur_lh > tree_lh - 1.0 && "partition rate opt reduces LnL");
        }
        
        // Optimizing branch lengths
        int my_iter = min(5,i+1);
        
        if (fixed_len == BRLEN_OPTIMIZE){
            double new_lh = tree->optimizeAllBranches(my_iter, logl_epsilon);
            ASSERT(new_lh > cur_lh - 1.0);
            cur_lh = new_lh;
        } else if (fixed_len == BRLEN_SCALE) {
            double scaling = 1.0;
            double new_lh = tree->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
            ASSERT(new_lh > cur_lh - 1.0);
            cur_lh = new_lh;
        }
        cout<<"Current log-likelihood at step "<<i<<": "<<cur_lh<<endl;
        if(fabs(cur_lh-tree_lh) < logl_epsilon) {
            tree_lh = cur_lh;
            break;
        }
        // make sure that the new logl is not so bad compared with previous logl
        ASSERT(cur_lh > tree_lh - 1.0 && "branch length opt reduces LnL");
        tree_lh = cur_lh;
    }
    //    cout <<"OPTIMIZE MODEL has finished"<< endl;
    if (write_info)
        writeInfo(cout);

    // write linked_models
    if (verbose_mode <= VB_MIN && write_info) {
        for (auto it = linked_models.begin(); it != linked_models.end(); it++)
            it->second->writeInfo(cout);
    }

    cout << "Parameters optimization took " << i-1 << " rounds (" << getRealTime()-begin_time << " sec)" << endl << endl;
    
    return tree_lh;
}

double PartitionModelPlen::optimizeParametersGammaInvar(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {
    outError("Thorough I+G parameter optimization does not work with edge-linked partition model yet");
    return 0.0;
}

void PartitionModelPlen::writeInfo(ostream &out) {
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    int ntrees = tree->size();
    if (!tree->fixed_rates) {
        out << "Partition-specific rates: ";
        for(int part = 0; part < ntrees; part++){
            out << " " << tree->part_info[part].part_rate;
        }
        out << endl;
    }
}

double PartitionModelPlen::optimizeGeneRate(double gradient_epsilon)
{
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    // BQM 22-05-2015: change to optimize individual rates
    double score = 0.0;
    size_t nsites = tree->getAlnNSite();
    
    vector<DoubleVector> brlen;
    brlen.resize(tree->branchNum);
    tree->getBranchLengths(brlen);
    double max_brlen = 0.0;
    for (size_t i = 0; i < brlen.size(); ++i) {
        for (size_t j = 0; j < brlen[i].size(); ++j) {
            if (brlen[i][j] > max_brlen) {
                max_brlen = brlen[i][j];
            }
        }
    }
    if (tree->part_order.empty()) tree->computePartitionOrder();
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+: score) schedule(dynamic) if(tree->num_threads > 1)
#endif
    for (int j = 0; j < tree->size(); j++) {
        int i = tree->part_order[j];
        double min_scaling = 1.0/tree->at(i)->getAlnNSite();
        double max_scaling = nsites / tree->at(i)->getAlnNSite();
        if (max_scaling < tree->part_info[i].part_rate)
            max_scaling = tree->part_info[i].part_rate;
        if (min_scaling > tree->part_info[i].part_rate)
            min_scaling = tree->part_info[i].part_rate;
        tree->part_info[i].cur_score = tree->at(i)->optimizeTreeLengthScaling(min_scaling, tree->part_info[i].part_rate, max_scaling, gradient_epsilon);
        score += tree->part_info[i].cur_score;
    }
    // now normalize the rates
    double sum = 0.0;
    size_t nsite = 0;
    for (size_t i = 0; i < tree->size(); ++i) {
        sum += tree->part_info[i].part_rate * tree->at(i)->aln->getNSite();
        if (tree->at(i)->aln->seq_type == SEQ_CODON && tree->rescale_codon_brlen)
            nsite += 3*tree->at(i)->aln->getNSite();
        else
            nsite += tree->at(i)->aln->getNSite();
    }
    sum /= nsite;
    
    if (sum > tree->params->max_branch_length / max_brlen) {
        outWarning("Too high (saturated) partition rates for proportional partition model!");
//        outWarning("Please switch to the edge-equal partition model via -q option instead of -spp");
//        exit(EXIT_FAILURE);
    }
    tree->scaleLength(sum);
    sum = 1.0/sum;
    for (size_t i = 0; i < tree->size(); ++i)
        tree->part_info[i].part_rate *= sum;
    return score;
}


int PartitionModelPlen::getNParameters(int brlen_type) {
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    int df = 0;
    for (PhyloSuperTreePlen::iterator it = tree->begin(); it != tree->end(); it++) {
        df += (*it)->getModelFactory()->model->getNDim() +
        (*it)->getModelFactory()->model->getNDimFreq() +
        (*it)->getModelFactory()->site_rate->getNDim();
    }
    df += tree->branchNum;
    if(!tree->fixed_rates)
        df += tree->size()-1;
    if (linked_alpha > 0.0)
        df ++;
    for (auto it = linked_models.begin(); it != linked_models.end(); it++) {
        bool fixed = it->second->fixParameters(false);
        df += it->second->getNDim() + it->second->getNDimFreq();
        it->second->fixParameters(fixed);
    }
    return df;
}

/*
int PartitionModelPlen::getNDim(){
    PhyloSuperTreePlen *tree = (PhyloSuperTreePlen*)site_rate->getTree();
    int ndim = tree->size() -1;
    return ndim;
}
*/
