//
//  iqtreemixhmm.cpp
//  tree
//
//  Created by Thomas Wong on 19/01/23.
//

#include "iqtreemixhmm.h"

IQTreeMixHmm::IQTreeMixHmm() : IQTreeMix(), PhyloHmm() {
}

IQTreeMixHmm::IQTreeMixHmm(Params &params, Alignment *aln, vector<IQTree*> &trees) : IQTreeMix(params, aln, trees), PhyloHmm(getAlnNSite(), trees.size()) {
    optimTree = -1;
    optimBranchGrp = -1;
}

// initialize the model
void IQTreeMixHmm::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    IQTreeMix::initializeModel(params, model_name, models_block);
    size_t i;
    
    // for all the unlinked substitution models,
    // set all trees to this tree
    if (!isLinkModel) {
        for (i=0; i<models.size(); i++) {
            models[i]->setTree(this);
        }
    }
    
    // handle the linked or unlinked site rate(s)
    site_rates.clear();
    site_rate_trees.clear();
    if (anySiteRate) {
        if (isLinkSiteRate) {
            // for linked site rate model, set all site rates to site_rates[0]
            site_rates.push_back(at(0)->getModelFactory()->site_rate);
            for (i=1; i<ntree; i++) {
                at(i)->getModelFactory()->site_rate = site_rates[0];
                at(i)->setRate(site_rates[0]);
            }
        } else {
            for (i=0; i<ntree; i++) {
                site_rates.push_back(at(i)->getModelFactory()->site_rate);
            }
        }
        // set their trees to this tree
        for (i=0; i<site_rates.size(); i++) {
            site_rates[i]->setTree(this);
        }
        for (i=0; i<ntree; i++) {
            site_rate_trees.push_back(this);
        }
    }
    
    // edge-length-restricted model is not appropriate for this HMM model
    if (isEdgeLenRestrict)
        outError("Edge-length-restricted model is inappropriate for HMM model. Use +T instead of +TR.");
    
    // build the branch ID
    computeBranchID();
    if (verbose_mode >= VB_MED) {
        showBranchGrp();
    }
}

// obtain the log-likelihoods for every site and every tree
// output site_like_cat[i * ntree + j] : log-likelihood of site nsite-i-1 and tree j
void IQTreeMixHmm::computeLogLikelihoodSiteTree(int updateTree) {
    
    if (updateTree > -1) {
        // only update a single tree
        computeLogLikelihoodSingleTree(updateTree);
    } else {
        // update all trees
        // compute likelihood for each tree
        for (int t=0; t<ntree; t++) {
            computeLogLikelihoodSingleTree(t);
        }
    }

    // reorganize the array
    // #pragma omp parallel for schedule(static) num_threads(num_threads) if (num_threads > 1)
    for (int j = 0; j < ntree; j++) {
        int k = nsite;
        int l = j;
        double* ptn_lh_arr = _ptn_like_cat + nptn * j;
        for (int i = 0; i < nsite; i++) {
            k--;
            int ptn = aln->getPatternID(k);
            // ptn_like_cat[i * ntree + j] = log-likelihood of site nsite-i-1 and tree j
            site_like_cat[l] = ptn_lh_arr[ptn];
            l += ntree;
        }
    }
}

double IQTreeMixHmm::computeLikelihood(double *pattern_lh, bool save_log_value) {
    computeLogLikelihoodSiteTree(optimTree);
    return computeBackLike();
}

double IQTreeMixHmm::optimizeBranchGroup(int branchgrp, double gradient_epsilon) {
    double score;
    int ndim;
    // for dimension = 1
    double len;
    // for dimension > 1

    optimTree = -1;
    optimBranchGrp = branchgrp;
    if (branchgrp < branch_group.size()) {
        ndim = branch_group[branchgrp].size();
        if (ndim == 1) {
            len = setSingleVariable();
            double negative_lh;
            double ferror,optx;
            if (verbose_mode >= VB_MED) {
                cout << "[IQTreeMixHmm::optimizeBranchGroup before] branchgrp = " << branchgrp << " single-variable = (" << setprecision(10) << len << ") ndim = 1" << endl;
            }
            optx = minimizeOneDimen(params->min_branch_length, len, MAX_LEN, gradient_epsilon, &negative_lh, &ferror);
            getSingleVariable(optx);
            if (verbose_mode >= VB_MED) {
                cout << "[IQTreeMixHmm::optimizeBranchGroup after] branchgrp = " << branchgrp << " single-variable = (" << setprecision(10) << optx << ") ndim = 1" << endl;
            }
            score = computeLikelihood();
        } else if (ndim > 1) {
            double* variables = new double[ndim + 1];
            double* upper_bound = new double[ndim + 1];
            double* lower_bound = new double[ndim + 1];
            bool* bound_check = new bool[ndim + 1];
            setBounds(lower_bound, upper_bound, bound_check);
            setVariables(variables);
            if (verbose_mode >= VB_MED) {
                cout << "[IQTreeMixHmm::optimizeBranchGroup before] branchgrp = " << branchgrp << " variables = (";
                for (int k = 1; k <= ndim; k++) {
                    if (k>1)
                        cout << ",";
                    cout << setprecision(10) << variables[k];
                }
                cout << ") ndim = " << ndim << endl;
            }
            score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
            getVariables(variables);
            if (verbose_mode >= VB_MED) {
                cout << "[IQTreeMixHmm::optimizeBranchGroup after] branchgrp = " << branchgrp << " variables = (";
                for (int k = 1; k <= ndim; k++) {
                    if (k>1)
                        cout << ",";
                    cout << setprecision(10) << variables[k];
                }
                cout << ") ndim = " << ndim << endl;
            }
            delete[] variables;
            delete[] upper_bound;
            delete[] lower_bound;
            delete[] bound_check;
        } else {
            optimBranchGrp = -1;
            score = computeLikelihood();
        }
    }
    optimBranchGrp = -1;
    return score;
}

double IQTreeMixHmm::optimizeAllBranchLensByBFGS(double gradient_epsilon, double logl_epsilon, int maxsteps) {
    size_t i,j;
    double score, pre_score, step;
    
    // collect the branch lengths of the tree
    getAllBranchLengths();
    score = computeLikelihood();
    step = 0;
    do {
        pre_score = score;
        step++;
        for (i = 0; i < branch_group.size(); i++) {
            score = optimizeBranchGroup(i, gradient_epsilon);
            cout << ".. Current HMM log-likelihood: " << score << endl;
        }
    } while (score - pre_score > logl_epsilon && step < maxsteps);
    return score;
}


void IQTreeMixHmm::startCheckpoint() {
    checkpoint->startStruct("IQTreeMixHmm" + convertIntToString(size()));
}

// ------------------------------------------------------------------

string IQTreeMixHmm::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    size_t i, ptn;
    int step, n, m, substep1, nsubstep1, nsubstep1_start, nsubstep1_max, nsubstep2_start, nsubstep2_max, substep2, nsubstep2, substep2_tot;
    double curr_epsilon;
    double prev_score, prev_score1, prev_score2, score, t_score;
    double gradient_epsilon = 0.0001;
    PhyloTree *ptree;
    
    // the edges with the same partition among the trees are initialized as the same length
    setAvgLenEachBranchGrp();

    cout << setprecision(5) << "Estimate model parameters (epsilon = " << logl_epsilon << ")" << endl;
    
    // minimum value of edge length
    if (verbose_mode >= VB_MED) {
        cout << "Minimum value of edge length is set to: " << setprecision(10) << params->min_branch_length << endl;
    }
    
    score = computeLikelihood();
    cout << "1. Initial HMM log-likelihood: " << score << endl;

    prev_score = score;

    for (step = 0; step < optimize_steps; step++) {
        
        // optimize tree branches
        score = optimizeAllBranchLensByBFGS(gradient_epsilon, logl_epsilon);
        if (verbose_mode >= VB_MED) {
            cout << "after optimizing branches, HMM likelihood = " << score << endl;
        }

        // optimize the linked subsitution model
        if (isLinkModel) {
            models[0]->optimizeParameters(gradient_epsilon);
            if (verbose_mode >= VB_MED) {
                score = computeLikelihood();
                cout << "after optimizing linked subsitution model, HMM likelihood = " << score << endl;
            }
        }

        // optimize the unlinked subsitution models one by one
        if (!isLinkModel) {
            for (i=0; i<ntree; i++) {
                models[i]->optimizeParameters(gradient_epsilon);
            }
            if (verbose_mode >= VB_MED) {
                score = computeLikelihood();
                cout << "after optimizing unlinked subsitution model, HMM likelihood = " << score << endl;
            }
        }

        // optimize the linked site rate model
        if (anySiteRate && isLinkSiteRate) {
            site_rates[0]->optimizeParameters(gradient_epsilon);
            if (verbose_mode >= VB_MED) {
                score = computeLikelihood();
                cout << "after optimizing linked site rate model, HMM likelihood = " << score << endl;
            }
        }

        // optimize the unlinked site rate model
        if (anySiteRate && !isLinkSiteRate) {
            for (i=0; i<ntree; i++) {
                site_rates[i]->optimizeParameters(gradient_epsilon);
            }
            if (verbose_mode >= VB_MED) {
                score = computeLikelihood();
                cout << "after optimizing unlinked site-rate model, HMM likelihood = " << score << endl;
            }
        }

        // optimize transition matrix and prob array
        score = PhyloHmm::optimizeParameters(gradient_epsilon);

        cout << step+2 << ". Current HMM log-likelihood: " << score << endl;
        if (score < prev_score + logl_epsilon) {
            // converged
            break;
        }

        if (verbose_mode >= VB_MED) {
            computeMaxPath();
            showSiteCatMaxLike(cout);
        }

        prev_score = score;
    }

    backLogLike = score;
    setCurScore(score);
    stop_rule.setCurIt(step);
    computeMaxPath();

    return getTreeString();
}

void IQTreeMixHmm::setNumThreads(int num_threads) {

    PhyloTree::setNumThreads(num_threads);

    for (size_t i = 0; i < size(); i++)
        at(i)->setNumThreads(num_threads);
}

/**
    test the best number of threads
*/
int IQTreeMixHmm::testNumThreads() {
    int bestNThres = at(0)->testNumThreads();
    setNumThreads(bestNThres);
    return bestNThres;
}

int IQTreeMixHmm::getNParameters() {
    int df = 0;
    int k;
    size_t i;
    
    if (verbose_mode >= VB_MED)
        cout << endl << "Number of parameters:" << endl;
    for (i=0; i<models.size(); i++) {
        k = models[i]->getNDim() + models[i]->getNDimFreq();
        if (verbose_mode >= VB_MED) {
            if (models.size() == 1)
                cout << " linked model : " << k << endl;
            else
                cout << " model " << i+1 << " : " << k << endl;
        }
        df += k;
    }
    for (i=0; i<site_rates.size(); i++) {
        if (verbose_mode >= VB_MED) {
            if (site_rates.size() == 1)
                cout << " linked site rate : " << site_rates[i]->getNDim() << endl;
            else
                cout << " siterate " << i+1 << " : " << site_rates[i]->getNDim() << endl;
        }
        df += site_rates[i]->getNDim();
    }
    // for branch parameters
    if (params->fixed_branch_length != BRLEN_FIX) {
        for (i=0; i<size(); i++) {
            k = at(i)->getNBranchParameters(BRLEN_OPTIMIZE);
            if (verbose_mode >= VB_MED)
                cout << " branches of tree " << i+1 << " : " << k << endl;
            df += k;
        }
    }
    // for transition matrix
    if (verbose_mode >= VB_MED)
        cout << " transition matrix : " << modelHmm->getNParameters() << endl;
    df += modelHmm->getNParameters();
    // for probability array
    if (verbose_mode >= VB_MED)
        cout << " probability array : " << ntree - 1 << endl;
    df += ntree - 1;

    if (verbose_mode >= VB_MED)
        cout << " == Total : " << df << " == " << endl << endl;
    return df;
}

// print out all the results to a file
void IQTreeMixHmm::printResults(const char *filename) {
    
    size_t i, j;
    ofstream out;
    out.open(filename);
    
    // report the estimated HMM parameters
    showParameters(out);
    out << endl;
    
    // show the assignment of the categories along sites with max likelihood
    showSiteCatMaxLike(out);
    
    out.close();
}

// print out the HMM estimated parameters
void IQTreeMixHmm::showParameters(ostream& out) {
    size_t i, j;
    modelHmm->showParameters(out);
    out << endl;
    out << "Estimated HMM probabilities :" << endl;
    for (i = 0; i < ntree; i++) {
        if (i > 0)
            out << "\t";
        out << fixed << setprecision(5) << prob[i];
    }
    out << endl << endl;
    
    out << "BEST HMM SCORE FOUND :" << fixed << setprecision(5) << backLogLike << endl;
}

// compute the log-likelihoods for a single tree t
void IQTreeMixHmm::computeLogLikelihoodSingleTree(int t) {
    double* pattern_lh_tree = _ptn_like_cat + nptn * t;
    // save the site rate's tree
    PhyloTree* ptree = at(t)->getRate()->getTree();
    // set the tree t as the site rate's tree
    // and compute the likelihood values
    at(t)->getRate()->setTree(at(t));
    at(t)->initializeAllPartialLh();
    at(t)->clearAllPartialLH();
    at(t)->computeLikelihood(pattern_lh_tree, true); // get the log-likelihood values
    // set back the previous site rate's tree
    at(t)->getRate()->setTree(ptree);
}

// get the branch lengths of all trees to the variable allbranchlens
void IQTreeMixHmm::getAllBranchLengths() {
    if (allbranchlens.size() < ntree)
        allbranchlens.resize(ntree);
    for (size_t i=0; i<ntree; i++)
        at(i)->saveBranchLengths(allbranchlens[i]);
}

// set the branch lengths of all trees from the variable allbranchlens
void IQTreeMixHmm::setAllBranchLengths() {
    for (size_t i=0; i<ntree; i++)
        at(i)->restoreBranchLengths(allbranchlens[i]);
}

//--------------------------------------
// optimization of branch lengths
//--------------------------------------

// the following three functions are for dimension = 1
double IQTreeMixHmm::computeFunction(double x) {
    getSingleVariable(x);
    return -computeLikelihood();
}

double IQTreeMixHmm::setSingleVariable() {
    // get the value from branch lengths
    size_t i;
    size_t ndim;
    int treeidx, branchidx;
    double x = 0.0;
    if (optimBranchGrp >= 0) {
        // collect the branch lengths of the tree
        getAllBranchLengths();
        ndim = branch_group[optimBranchGrp].size();
        if (ndim > 0) {
            treeidx = branch_group[optimBranchGrp].at(0) / nbranch;
            branchidx = branch_group[optimBranchGrp].at(0) % nbranch;
            if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
                x = allbranchlens[treeidx].at(branchidx);
            else
                cout << "[IQTreeMixHmm::setSingleVariable] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
        } else {
            cout << "[IQTreeMixHmm::setSingleVariable] Error occurs! ndim = " << ndim << endl;
        }
    } else {
        cout << "[IQTreeMixHmm::setSingleVariable] Error occurs! optimBranchGrp = " << optimBranchGrp << endl;
    }
    return x;
}

void IQTreeMixHmm::getSingleVariable(double x) {
    // save the values to branch lengths
    size_t i;
    size_t ndim;
    int treeidx, branchidx;
    if (optimBranchGrp >= 0) {
        // collect the branch lengths of the tree
        getAllBranchLengths();
        ndim = branch_group[optimBranchGrp].size();
        if (ndim > 0) {
            treeidx = branch_group[optimBranchGrp].at(0) / nbranch;
            branchidx = branch_group[optimBranchGrp].at(0) % nbranch;
            if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
                allbranchlens[treeidx].at(branchidx) = x;
            else
                cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
        } else {
            cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! ndim = " << ndim << endl;
        }
        setAllBranchLengths();
    } else {
        cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! optimBranchGrp = " << optimBranchGrp << endl;
    }
}

// the following four functions are for dimension > 1
double IQTreeMixHmm::targetFunk(double x[]) {
    getVariables(x);
    return -computeLikelihood();
}

void IQTreeMixHmm::setVariables(double *variables) {
    // copy the values from branch lengths
    size_t i;
    size_t ndim;
    int treeidx, branchidx;
    if (optimBranchGrp >= 0) {
        // collect the branch lengths of the tree
        getAllBranchLengths();
        ndim = branch_group[optimBranchGrp].size();
        for (i=0; i<ndim; i++) {
            treeidx = branch_group[optimBranchGrp].at(i) / nbranch;
            branchidx = branch_group[optimBranchGrp].at(i) % nbranch;
            if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
                variables[i+1] = allbranchlens[treeidx].at(branchidx);
            else
                cout << "[IQTreeMixHmm::setVariables] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
        }
        if (ndim == 0) {
            cout << "[IQTreeMixHmm::setVariables] Error occurs! ndim = " << ndim << endl;
        }
    } else {
        cout << "[IQTreeMixHmm::setVariables] Error occurs! optimBranchGrp = " << optimBranchGrp << endl;
    }
}


void IQTreeMixHmm::getVariables(double *variables) {
    // save the values to branch lengths
    size_t i;
    size_t ndim;
    int treeidx, branchidx;
    if (optimBranchGrp >= 0) {
        // collect the branch lengths of the tree
        getAllBranchLengths();
        ndim = branch_group[optimBranchGrp].size();
        for (i=0; i<ndim; i++) {
            treeidx = branch_group[optimBranchGrp].at(i) / nbranch;
            branchidx = branch_group[optimBranchGrp].at(i) % nbranch;
            if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
                allbranchlens[treeidx].at(branchidx) = variables[i+1];
            else
                cout << "[IQTreeMixHmm::getVariables] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
        }
        if (ndim == 0) {
            cout << "[IQTreeMixHmm::getVariables] Error occurs! ndim = " << ndim << endl;
        }
        setAllBranchLengths();
    } else {
        cout << "[IQTreeMixHmm::getVariables] Error occurs! optimBranchGrp = " << optimBranchGrp << endl;
    }
}

void IQTreeMixHmm::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    size_t i;
    size_t ndim;
    if (optimBranchGrp >= 0) {
        ndim = branch_group[optimBranchGrp].size();
        if (verbose_mode >= VB_MED) {
            cout << "[IQTreeMixHmm::setBounds] optimBranchGrp = " << optimBranchGrp << ", ndim = " << ndim << endl;
        }
        for (i = 1; i <= ndim; i++) {
            lower_bound[i] = params->min_branch_length;
            upper_bound[i] = MAX_LEN;
            bound_check[i] = false;
        }
        if (ndim == 0) {
            cout << "[IQTreeMixHmm::setBounds] Error occurs! ndim = " << ndim << endl;
        }
    } else {
        cout << "[IQTreeMixHmm::setBounds] Error occurs! optimBranchGrp = " << optimBranchGrp << endl;
    }
}

void IQTreeMixHmm::showBranchGrp() {
    cout << "Branch Group:" << endl;
    for (size_t i=0; i<branch_group.size(); i++) {
        cout << "  Grp " << i << endl;
        for (size_t j=0; j<branch_group[i].size(); j++) {
            cout << "     " << branch_group[i].at(j) << endl;
        }
    }
}

/**
         If there are multiple branches belonging to the same group
         set all the branches of the same group to their average
 */
void IQTreeMixHmm::setAvgLenEachBranchGrp() {
    size_t i,j;
    size_t treeIdx,branchIdx;
    double grp_len;
    int dim;
    
    // collect the branch lengths of the tree
    getAllBranchLengths();
    
    for (i = 0; i < branch_group.size(); i++) {
        grp_len = 0.0;
        dim = branch_group[i].size();
        for (j = 0; j < dim; j++) {
            treeIdx = branch_group[i].at(j) / nbranch;
            branchIdx = branch_group[i].at(j) % nbranch;
            grp_len += allbranchlens[treeIdx].at(branchIdx);
        }
        grp_len = grp_len / (double)dim;
        for (j = 0; j < dim; j++) {
            treeIdx = branch_group[i].at(j) / nbranch;
            branchIdx = branch_group[i].at(j) % nbranch;
            allbranchlens[treeIdx].at(branchIdx) = grp_len;
        }
    }
    
    // save the updated branch lengths of the tree
    setAllBranchLengths();
}
