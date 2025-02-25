//
//  iqtreemixhmm.cpp
//  tree
//
//  Created by Thomas Wong on 19/01/23.
//

#include "iqtreemixhmm.h"

IQTreeMixHmm::IQTreeMixHmm() : IQTreeMix(), PhyloHmm() {
    optimTree = -1;
    optimBranchGrp = -1;
    objFun = 1; // default: MAST model
    objAlgo = new string[2];
    objAlgo[0] = "HMM";
    objAlgo[1] = "MAST";
    isTMixOptimEngine = false;
    siteTypes = NULL;
}

IQTreeMixHmm::IQTreeMixHmm(Params &params, Alignment *aln) : IQTreeMix(params, aln), PhyloHmm(getAlnNSite(), ntree) {
    optimTree = -1;
    optimBranchGrp = -1;
    objFun = 1; // default: MAST model
    objAlgo = new string[2];
    objAlgo[0] = "HMM";
    objAlgo[1] = "MAST";
    isTMixOptimEngine = false;
    siteTypes = NULL;
    setSiteTypes();
}

IQTreeMixHmm::~IQTreeMixHmm() {
    delete[] objAlgo;
}

// initialize the model
void IQTreeMixHmm::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    IQTreeMix::initializeModel(params, model_name, models_block);
    size_t i;
    
    // build the branch ID
    if (branch_group.size() == 0) {
        computeBranchID();
    }
    if (verbose_mode >= VB_MED) {
        showBranchGrp();
    }

    // initialize the transition model
    initializeTransitModel(params);
}

// initialize the transition model (override it for different transition model)
// this function is invoked inside the constructor
void IQTreeMixHmm::initializeTransitModel(Params &params) {
    // by default, it uses the HMM simple transition matrix
    // the transition probabilities between different categories are the same

    if (params.optimize_params_use_hmm_tm) {
        int ntype = 3; // (0) parsimony informatic, (1) invariant, and (2) uninformatic
        string* types = new string[ntype];
        types[0] = "parsimony informatic";
        types[1] = "invariant";
        types[2] = "uninformatic";
        modelHmm = new ModelHmmTm(ncat, ntype, siteTypes, aln->getNSite(), types);
        delete[] types;
    } else if (params.optimize_params_use_hmm_gm) {
        modelHmm = new ModelHmmGm(ncat);
    } else {
        modelHmm = new ModelHmm(ncat);
    }
    modelHmm->initialize_transitLog();
    
    // show the transition model
    // if (params.treemix_optimize_methods != "mast")
    //    cout << "HMM transition model: " << modelHmm->getFullName() << " (" << modelHmm->getName() << ")"<< endl;

    // set the associated PhyloHmm of modelHmm to this
    modelHmm->setPhyloHmm(this);
}

// get the type of all sites for type-dependent HMM model
// 0 - parsimony informatic; 1 - invariant (including constant and e.g. GS--G-GGG (S = G/C));
// 2 - uninformatic but not invariant (e.g. GTTTTTT)
void IQTreeMixHmm::setSiteTypes() {
    int ptn, i, ctype;
    if (aln->getNSite() > 0) {
        siteTypes = new int[aln->getNSite()];
        for (i = 0; i < aln->getNSite(); i++) {
            ptn = aln->getPatternID(i);
            ctype = 2;
            if (aln->at(ptn).isInformative()) {
                ctype = 0;
            } else if (aln->at(ptn).isConst() || aln->at(ptn).isInvariant()) {
                ctype = 1;
            }
            siteTypes[i] = ctype;
        }
    }
}

// set the tree weights according to the marginal probabilities along the sites
void IQTreeMixHmm::setWeightToMarginalProb() {
    bool need_computeLike = true;
    int update_which_tree = -1; // recompute all trees' likelihoods
    double* mar_prob;
    double weight_sum;
    size_t i,j;
    // reset the tree weights
    for (i = 0; i < ntree; i++) {
        weights[i] = 0.0;
    }
    // get marginal probabilities along each site for each tree
    getMarginalProb(need_computeLike, update_which_tree);
    mar_prob = marginal_prob;
    for (j = 0; j < nsite; j++) {
        for (i = 0; i < ntree; i++) {
            weights[i] += mar_prob[i];
        }
        mar_prob += ntree;
    }
    // normalize the tree weights
    weight_sum = 0.0;
    for (i = 0; i < ntree; i++) {
        weight_sum += weights[i];
    }
    for (i = 0; i < ntree; i++) {
        weights[i] = weights[i] / weight_sum;
    }
    // If there are multiple tree weights belonging to the same group
    // set all the tree weights of the same group to their average
    // when weightGrpExist
    checkWeightGrp();
    // update the weight_logs
    for (i = 0; i < ntree; i++) {
        weight_logs[i] = log(weights[i]);
    }
    /*
    // show the tree weights according to the posterior probability
    cout << "According to the posterior probabilities from the HMM model along the sites, the tree weights are :" << endl;
    for (i=0; i<ntree; i++) {
        if (i>0)
            cout << ",";
        cout << weights[i];
    }
    cout << endl;
    */
}

// obtain the log-likelihoods for every site and every tree
// output site_like_cat[i * ntree + j] : log-likelihood of site nsite-i-1 and tree j
void IQTreeMixHmm::computeLogLikelihoodSiteTree(int updateTree) {
    
    if (isLinkSiteRate && updateTree > 0) {
        // Store the RHAS variables of tree 0 to the array rhas_var
        storeTree0RHAS();
        // Replace the RHAS variables of tree 'updateTree' by those in the array rhas_var
        copyRHASfrTree0(updateTree);
    }
    
    if (updateTree > -1) {
        // only update a single tree
        computeLogLikelihoodSingleTree(updateTree);
    } else {
        // update all trees
        // Store the RHAS variables of tree 0 to the array rhas_var for linked site rate
        if (isLinkSiteRate)
            storeTree0RHAS();
        // compute likelihood for each tree
        for (int t=0; t<ntree; t++) {
            // Replace the RHAS variables of tree t by those in the array rhas_var for linked site rate
            if (isLinkSiteRate)
                copyRHASfrTree0(t);
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
    if (objFun == 1)
        return IQTreeMix::computeLikelihood(pattern_lh, save_log_value);
    
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
    ndim = getNDim();
    if (ndim == 1 || isEdgeLenRestrict) {
        len = setSingleVariable();
        double negative_lh;
        
        
        
        double ferror,optx;
        double bvalue;
        if (verbose_mode >= VB_MED) {
            bvalue = len;
        }
        optx = minimizeOneDimen(params->min_branch_length, len, params->max_branch_length, gradient_epsilon, &negative_lh, &ferror);
        getSingleVariable(optx);
        if (verbose_mode >= VB_MED) {
            cout << "branchgrp:" << branchgrp << "; len:" << setprecision(10)<< bvalue << "->" << optx << "; ndim:" << ndim << endl;
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
    optimBranchGrp = -1;
    return score;
}

double IQTreeMixHmm::optimizeAllBranchLensByBFGS(double gradient_epsilon, double logl_epsilon, int maxsteps) {
    double score, pre_score, step;
    
    // collect the branch lengths of the tree
    getAllBranchLengths();
    score = computeLikelihood();
    step = 0;
    do {
        pre_score = score;
        step++;
        for (int i = 0; i < branch_group.size(); i++) {
            score = optimizeBranchGroup(i, gradient_epsilon);
            if (objFun == 1) // MAST
                cout << ".. Current MAST log-likelihood: " << score << endl;
            else // HMM
                cout << ".. Current HMM log-likelihood: " << score << endl;
        }
    } while (score - pre_score > logl_epsilon && step < maxsteps);
    return score;
}

double IQTreeMixHmm::optimizeAllBranches(double* pattern_mix_lh, int my_iterations, double tolerance, int maxNRStep) {
    double score;
    
    // update the ptn_freq array according to the marginal probabilities along each site for each tree
    computeFreqArray(pattern_mix_lh);
    for (int i = 0; i < ntree; i++) {
        IQTreeMix::optimizeAllBranchesOneTree(i, my_iterations, tolerance, maxNRStep);
    }
    score = computeLikelihood();
    return score;
}

double IQTreeMixHmm::optimizeAllSubstModels(double gradient_epsilon, double* pattern_mix_lh) {
    double score;
    if (isLinkModel) {
        // for linked subsitution model
        // use BFGS method
        resetPtnOrigFreq();
        models[0]->optimizeParameters(gradient_epsilon);
    } else {
        // for unlinked subsitution model
        // use EM method
        computeFreqArray(pattern_mix_lh);
        for (int i = 0; i < ntree; i++) {
            models[i]->optimizeParameters(gradient_epsilon);
        }
    }
    score = computeLikelihood();
    return score;
}

double IQTreeMixHmm::optimizeAllRHASModels(double gradient_epsilon, double score, double* pattern_mix_lh) {
    if (anySiteRate) {
        if (isLinkSiteRate) {
            // for linked RHAS model
            // use BFGS method
            resetPtnOrigFreq();
            site_rates[0]->optimizeParameters(gradient_epsilon);
        } else {
            // for unlinked RHAS model
            // use EM method
            computeFreqArray(pattern_mix_lh);
            for (int i = 0; i < ntree; i++) {
                site_rates[i]->optimizeParameters(gradient_epsilon);
            }
        }
        score = computeLikelihood();
    }
    return score;
}

void IQTreeMixHmm::startCheckpoint() {
    checkpoint->startStruct("IQTreeMixHmm" + convertIntToString(size()));
}

// ------------------------------------------------------------------

string IQTreeMixHmm::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    string ans;

    if (params->treemix_optimize_methods == "hmm") {
        objFun = 0;
        ans = optimizeModelParamHMM(printInfo, logl_epsilon);
    }
    
    else if (params->treemix_optimize_methods == "hmm2mast") {
        objFun = 0;
        optimizeModelParamHMM(printInfo, logl_epsilon);
        ans = optimizeModelParamMAST(printInfo, logl_epsilon);
    }
    
    else if (params->treemix_optimize_methods == "mast") {
        isTMixOptimEngine = true;
        objFun = 1;
        ans = IQTreeMix::optimizeModelParameters(printInfo, params->treemix_eps);
        isTMixOptimEngine = false;
    }

    else if (params->treemix_optimize_methods == "mast2hmm") {
        isTMixOptimEngine = true;
        objFun = 1;
        IQTreeMix::optimizeModelParameters(printInfo, params->treemix_eps);
        isTMixOptimEngine = false; 
        objFun = 0;
        params->HMM_no_avg_brlen = true;
        ans = optimizeModelParamHMM(printInfo, logl_epsilon);
    }
    
    else if (params->treemix_optimize_methods == "hmast") {
        objFun = 1;
        ans = optimizeModelParamMAST(printInfo, params->treemix_eps);
    }

    else {
        outError("Error! Unknown option: " + params->treemix_optimize_methods);
    }
    
    return ans;
}


string IQTreeMixHmm::optimizeModelParamHMM(bool printInfo, double logl_epsilon) {
    int step;
    double prev_score, score;
    double gradient_epsilon = 0.0001;
    
    // the edges with the same partition among the trees are initialized as the same length
    if (params->fixed_branch_length != BRLEN_FIX && !params->HMM_no_avg_brlen) {
        setAvgLenEachBranchGrp();
    } else if (verbose_mode >= VB_MED) {
        cout << "No averaging the branch lengths during initialisation" << endl;
    }

    // change the objective function to backlikelihood
    objFun = 0;

    if (printInfo)
        cout << setprecision(5) << "Estimate model parameters (epsilon = " << logl_epsilon << ")" << endl;
    
    // minimum value of edge length
    if (verbose_mode >= VB_MED) {
        cout << "Minimum value of edge length is set to: " << setprecision(10) << params->min_branch_length << endl;
    }

    // minimum value of HMM same-category transition probability
    if (printInfo && Params::getInstance().HMM_min_stran > 1e-10) {
        cout << "Minimum value of HMM same-category transition probability is set to: " << Params::getInstance().HMM_min_stran << endl;
    }

    score = computeLikelihood();

    if (printInfo)
        cout << "1. Initial HMM log-likelihood: " << score << endl;

    // first optimize prob array
    score = optimizeProbEM();
    if (verbose_mode >= VB_MED)
        cout << "after optimizing probability array, HMM likelihood = " << score << endl;
    
    prev_score = score;

    for (step = 0; step < optimize_steps; step++) {

        // optimize tree branches
        // if params->HMM_no_avg_brlen, then no optimization of branch lengths at the first iteration
        // if (params->fixed_branch_length != BRLEN_FIX &&
        //    !(step==0 && params->HMM_no_avg_brlen)) {
        if (params->fixed_branch_length != BRLEN_FIX) {
            if (isEdgeLenRestrict) {
                // +TR model : branches with the same partition information across the trees are restricted the same
                // use BFGS method
                score = optimizeAllBranchLensByBFGS(gradient_epsilon, logl_epsilon);
            } else {
                score = optimizeAllBranches();
            }
            if (verbose_mode >= VB_MED)
                cout << "after optimizing branch lengths, HMM likelihood = " << score << endl;
        }

        // optimize all subsitution models
        score = optimizeAllSubstModels(gradient_epsilon);
        if (verbose_mode >= VB_MED)
            cout << "after optimizing substution models, HMM likelihood = " << score << endl;

        // optimize all site rate models
        score = optimizeAllRHASModels(gradient_epsilon, score);
        if (verbose_mode >= VB_MED)
            cout << "after optimizing RHAS models, HMM likelihood = " << score << endl;

        // optimize transition matrix and prob array
        score = PhyloHmm::optimizeParameters(gradient_epsilon);
        if (verbose_mode >= VB_MED)
            cout << "after optimizing transition matrix and prob array, HMM likelihood = " << score << endl;

        if (printInfo)
            cout << step+2 << ". Current HMM log-likelihood: " << score << endl;

        if (score < prev_score + logl_epsilon)
            // converged
            break;

        /*
        if (verbose_mode >= VB_MED) {
            // computeMaxPath();
            showSiteCatMaxLike(cout, false);
        }
        */

        prev_score = score;
    }

    backLogLike = score;
    // computeMaxPath();

    // set the tree weights according to the marginal probabilities
    if (!isTreeWeightFixed) {
        setWeightToMarginalProb();
    }
    
    // compute the corresponding score according to the tree mixture formula
    objFun = 1;
    score = computeLikelihood();
    
    if (printInfo)
        cout << "Converted tree mixture likelihood = " << score << endl;

    setCurScore(score);
    stop_rule.setCurIt(step);

    return getTreeString();
}

// Optimize the MAST model starting from the current parameter values
string IQTreeMixHmm::optimizeModelParamMAST(bool printInfo, double logl_epsilon) {
    int step;
    double prev_score, score;
    double gradient_epsilon = 0.0001;
    double* pattern_mix_lh;
    int max_steps_tree_weight = 3;
    bool tree_weight_converge;

    // allocate memory
    pattern_mix_lh = new double[ntree * nptn];
    
    // change the objective function to MAST model
    objFun = 1;

    cout << setprecision(5) << "Estimate MAST model parameters (epsilon = " << logl_epsilon << ")" << endl;

    score = computeLikelihood();
    
    cout << "1. Initial MAST log-likelihood: " << score << endl;

    prev_score = score;

    for (step = 0; step < optimize_steps; step++) {
        
        // optimize tree branches
        if (isEdgeLenRestrict) {
            // +TR model : branches with the same partition information across the trees are restricted the same
            score = optimizeAllBranchLensByBFGS(gradient_epsilon, logl_epsilon);
        } else {
            score = optimizeAllBranches(pattern_mix_lh);
        }

        // optimize all subsitution models
        score = optimizeAllSubstModels(gradient_epsilon, pattern_mix_lh);

        // optimize all site rate models
        score = optimizeAllRHASModels(gradient_epsilon, score, pattern_mix_lh);
        
        // optimize tree weights
        IQTreeMix::optimizeTreeWeightsByEM(pattern_mix_lh, logl_epsilon, max_steps_tree_weight, tree_weight_converge);

        cout << step+2 << ". Current MAST log-likelihood: " << score << endl;

        if (score < prev_score + logl_epsilon)
            // converged
            break;

        prev_score = score;
    }

    setCurScore(score);
    stop_rule.setCurIt(step);

    delete[] pattern_mix_lh;

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
        df += site_rates[i]->getNDim();
        if (verbose_mode >= VB_MED) {
            if (isLinkSiteRate) {
                cout << " linked site rate : " << site_rates[i]->getNDim() << endl;
            } else {
                cout << " siterate " << i+1 << " : " << site_rates[i]->getNDim() << endl;
            }
        }
        if (isLinkSiteRate)
            break;
    }
    // for branch parameters
    if (params->fixed_branch_length != BRLEN_FIX) {
        if (isEdgeLenRestrict) {
            if (verbose_mode >= VB_MED)
                cout << " branch groups (for branch-len-restricted) : " << branch_group.size() << endl;
            df += branch_group.size();
        } else {
            for (i=0; i<size(); i++) {
                k = at(i)->getNBranchParameters(BRLEN_OPTIMIZE);
                if (verbose_mode >= VB_MED)
                    cout << " branches of tree " << i+1 << " : " << k << endl;
                df += k;
            }
        }
    }
    if (objFun == 0) {
        // for transition matrix
        if (verbose_mode >= VB_MED)
            cout << " transition matrix : " << modelHmm->getNParameters() << endl;
        df += modelHmm->getNParameters();
        // for probability array
        if (verbose_mode >= VB_MED)
            cout << " probability array : " << ntree - 1 << endl;
        df += ntree - 1;
    } else {
        // for MAST
        // for tree weight
        if (!isTreeWeightFixed) {
            if (weightGrpExist) {
                if (verbose_mode >= VB_MED)
                    cout << " tree weight groups (for weight-restricted) : " << (weight_group_member.size() - 1) << endl;
                df += (weight_group_member.size() - 1);
            } else {
                if (verbose_mode >= VB_MED)
                    cout << " tree weights : " << (size() - 1) << endl;
                df += ntree - 1;
            }
        }
    }

    if (verbose_mode >= VB_MED)
        cout << " == Total : " << df << " == " << endl << endl;
    return df;
}

// print out all the results to a file
// cat_assign_method:
//  0 - the categories along sites is assigned according to the path with maximum probability (default)
//  1 - the categories along sites is assigned according to the max posterior probability
void IQTreeMixHmm::printResults(const char *filename, int cat_assign_method, int* numSiteCat) {
    
    size_t i, j;
    ofstream out;
    out.open(filename);
    
    // report the estimated HMM parameters
    showParameters(out);
    out << endl;
    
    // show the assignment of the categories along the path with max likelihood
    showSiteCatMaxLike(out, true, cat_assign_method, numSiteCat);
    
    out.close();
}

// print out the marginal probabilities to a file
void IQTreeMixHmm::printMarginalProb(const char *filename) {
    
    size_t i, j;
    ofstream out;
    out.open(filename);

    computeLogLikelihoodSiteTree();
    computeMarginalProb(&out);

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

// show the branch lengths of all trees
void IQTreeMixHmm::showAllBranchLengths() {
    getAllBranchLengths();
    for (size_t i=0; i<ntree; i++) {
        cout << "The branch lengths of tree " << i+1 << endl;
        for (size_t j=0; j<allbranchlens[i].size(); j++) {
            if (j>0)
                cout << ", ";
            cout << allbranchlens[i].at(j);
        }
        cout << endl;
    }

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
    int ndim, i;
    int treeidx, branchidx;
    double x = 0.0;
    // collect the branch lengths of the tree
    getAllBranchLengths();
    ndim = getNDim();
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
    return x;
}

void IQTreeMixHmm::getSingleVariable(double x) {
    // save the values to branch lengths
    int ndim, i;
    int treeidx, branchidx;
    // collect the branch lengths of the tree
    getAllBranchLengths();
    ndim = getNDim();
    if (ndim == 0)
        cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! ndim = " << ndim << endl;
    for (i = 0; i < ndim; i++) {
        treeidx = branch_group[optimBranchGrp].at(i) / nbranch;
        branchidx = branch_group[optimBranchGrp].at(i) % nbranch;
        if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
            allbranchlens[treeidx].at(branchidx) = x;
        else
            cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
    }
    /*
    if (ndim > 0) {
        treeidx = branch_group[optimBranchGrp].at(0) / nbranch;
        branchidx = branch_group[optimBranchGrp].at(0) % nbranch;
        if (treeidx < ntree && branchidx < allbranchlens[treeidx].size())
            allbranchlens[treeidx].at(branchidx) = x;
        else
            cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! treeidx = " << treeidx << ", branchidx = " << branchidx << endl;
    } else {
        cout << "[IQTreeMixHmm::getSingleVariable] Error occurs! ndim = " << ndim << endl;
    }*/
    setAllBranchLengths();
}

// the following four functions are for dimension > 1
double IQTreeMixHmm::targetFunk(double x[]) {

    getVariables(x);
    return -computeLikelihood();
}

void IQTreeMixHmm::setVariables(double *variables) {
    // copy the values from branch lengths
    int ndim, i;
    int treeidx, branchidx;
    
    if (isTMixOptimEngine)
        return IQTreeMix::setVariables(variables);

    // collect the branch lengths of the tree
    getAllBranchLengths();
    ndim = getNDim();
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
}


void IQTreeMixHmm::getVariables(double *variables) {
    // save the values to branch lengths
    int ndim, i;
    int treeidx, branchidx;

    if (isTMixOptimEngine)
        return IQTreeMix::getVariables(variables);

    // collect the branch lengths of the tree
    getAllBranchLengths();
    ndim = getNDim();
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
}

void IQTreeMixHmm::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    int ndim, i;

    if (isTMixOptimEngine)
        return IQTreeMix::setBounds(lower_bound, upper_bound, bound_check);

    ndim = getNDim();
    if (verbose_mode >= VB_MED) {
        cout << "[IQTreeMixHmm::setBounds] optimBranchGrp = " << optimBranchGrp << ", ndim = " << ndim << endl;
    }
    for (i = 1; i <= ndim; i++) {
        lower_bound[i] = params->min_branch_length;
        upper_bound[i] = params->max_branch_length;
        bound_check[i] = false;
    }
    if (ndim == 0) {
        cout << "[IQTreeMixHmm::setBounds] Error occurs! ndim = " << ndim << endl;
    }
}


int IQTreeMixHmm::getNDim() {
    if (isTMixOptimEngine)
        return IQTreeMix::getNDim();
    
    if (optimBranchGrp >= 0 && optimBranchGrp < branch_group.size()) {
        return branch_group[optimBranchGrp].size();
    } else {
        return 0;
    }
}

void IQTreeMixHmm::showBranchGrp() {
    cout << "Branch Group:" << endl;
    for (size_t i=0; i<branch_group.size(); i++) {
        cout << "  Grp " << i << endl;
        for (size_t j=0; j<branch_group[i].size(); j++) {
            if (j > 0)
                cout << ", ";
            else
                cout << "    ";
            cout << branch_group[i].at(j);
        }
        cout << endl;
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
            if (allbranchlens[treeIdx].at(branchIdx) < params->min_branch_length)
                grp_len += params->min_branch_length;
            else
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

// update the ptn_freq array according to the marginal probabilities along each site for each tree
void IQTreeMixHmm::computeFreqArray(double* pattern_mix_lh, bool need_computeLike, int update_which_tree) {
    
    if (objFun == 1 || isTMixOptimEngine) // MAST
        return IQTreeMix::computeFreqArray(pattern_mix_lh, need_computeLike, update_which_tree);
    
    double* mar_prob;
    // get marginal probabilities along each site for each tree
    getMarginalProb(need_computeLike, update_which_tree);
    // #pragma omp parallel for schedule(dynamic) num_threads(num_threads) if (num_threads > 1)
    for (size_t i = 0; i < ntree; i++) {
        PhyloTree* tree = at(i);
        // reset the array ptn_freq
        memset(tree->ptn_freq, 0, sizeof(double)*nptn);
        mar_prob = marginal_prob + i;
        for (size_t j = 0; j < nsite; j++) {
            int ptn = aln->getPatternID(j);
            tree->ptn_freq[ptn] += mar_prob[0];
            mar_prob += ntree;
        }
    }
}

// get marginal probabilities along each site for each tree
void IQTreeMixHmm::getMarginalProb(bool need_computeLike, int update_which_tree) {
    if (need_computeLike) {
        computeLogLikelihoodSiteTree(update_which_tree);
    }
    computeMarginalProb();
}
