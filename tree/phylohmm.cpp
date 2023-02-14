//
//  phylohmm.h
//  tree
//
//  Created by Thomas Wong on 02/02/23.
//

#include "phylohmm.h"

PhyloHmm::PhyloHmm() {
    nsite = ncat = 0;
    prob = NULL;
    prob_log = NULL;
    site_like_cat = NULL;
    site_categories = NULL;
    work_arr = NULL;
    next_cat = NULL;
    bwd_array = NULL;
    fwd_array = NULL;
    marginal_prob = NULL;
    marginal_tran = NULL;
}

PhyloHmm::PhyloHmm(int n_site, int n_cat) {
    
    nsite = n_site;
    ncat = n_cat;
    
    // allocate memory for the arrays
    size_t prob_size = get_safe_upper_limit(ncat);
    size_t site_like_cat_size = get_safe_upper_limit(nsite) * ncat;
    size_t site_transition_size = get_safe_upper_limit(nsite-1) * ncat * ncat;
    
    prob = aligned_alloc<double>(prob_size);
    prob_log = aligned_alloc<double>(prob_size);
    site_like_cat = aligned_alloc<double>(site_like_cat_size);
    site_categories = aligned_alloc<int>(get_safe_upper_limit(nsite));
    work_arr = aligned_alloc<double>(prob_size * 2);
    next_cat = aligned_alloc<int>(site_like_cat_size);
    
    // allocate memory for the backward and forward algorithm
    site_like_cat_fwd = aligned_alloc<double>(site_like_cat_size);
    bwd_array = aligned_alloc<double>(site_like_cat_size);
    fwd_array = aligned_alloc<double>(site_like_cat_size);

    // allocate memory for the marginal probabilities
    marginal_prob = aligned_alloc<double>(site_like_cat_size);
    marginal_tran = aligned_alloc<double>(site_transition_size);

    // the probabilities are initialized as equally distributed
    double init_prob_value = 1.0/((double)ncat);
    for (size_t i = 0; i < ncat; i++)
        prob[i] = init_prob_value;
    computeLogProb();
    
    // initialize the transition model
    initializeTransitModel();
}

PhyloHmm::~PhyloHmm() {
    aligned_free(prob);
    aligned_free(prob_log);
    aligned_free(site_like_cat);
    aligned_free(site_categories);
    aligned_free(work_arr);
    aligned_free(next_cat);
    
    aligned_free(site_like_cat_fwd);
    aligned_free(bwd_array);
    aligned_free(fwd_array);
    aligned_free(marginal_prob);
    
    delete(modelHmm);
}

// initialize the transition model (override it for different transition model)
// this function is invoked inside the constructor
void PhyloHmm::initializeTransitModel() {
    // by default, it uses the HMM simple transition matrix
    // the transition probabilities between different categories are the same
    modelHmm = new ModelHmm(ncat);
    // set the associated PhyloHmm of modelHmm to this
    modelHmm->setPhyloHmm(this);
}

// compute backward log-likelihood
// prerequisite: array site_like_cat has been updated (i.e. computeLogLikelihoodSiteTree() has been invoked)
double PhyloHmm::computeBackLike(bool showInterRst) {
    int showlines = 5;
    size_t pre_k = 0;
    size_t k;
    size_t i,j;
    double* pre_work;
    double* work;
    double* site_lh_arr;
    double* transit_arr;
    site_lh_arr = site_like_cat;
    memcpy(work_arr, site_lh_arr, sizeof(double) * ncat);
    pre_work = work_arr;
    for (i = 1; i < nsite; i++) {
        k = (pre_k + 1) % 2;
        work = work_arr + k * ncat;
        site_lh_arr += ncat;
        transit_arr = modelHmm->getTransitLog();
        for (j = 0; j < ncat; j++) {
            work[j] = logDotProd(transit_arr, pre_work, ncat) + site_lh_arr[j];
            transit_arr += ncat;
        }
        pre_k = k;
        pre_work = work;
        // show the intermediate results
        if (showInterRst && i <= showlines) {
            for (j = 0; j < ncat; j++) {
                if (j > 0)
                    cout << "\t";
                cout << work[j];
            }
            cout << endl;
        }
    }
    return logDotProd(prob_log, pre_work, ncat);
}

// path with max log-likelihood
double PhyloHmm::computeMaxPath() {
    size_t pre_k = 0;
    size_t i,j,k,l;
    double* pre_work;
    double* work;
    double* site_lh_arr;
    double* transit_arr;
    int* next_cat_arr;
    double v;
    site_lh_arr = site_like_cat;
    memcpy(work_arr, site_lh_arr, sizeof(double) * ncat);
    pre_work = work_arr;
    
    for (i = 1; i < nsite; i++) {
        k = (pre_k + 1) % 2;
        work = work_arr + k * ncat;
        site_lh_arr += ncat;
        transit_arr = modelHmm->getTransitLog();
        next_cat_arr = next_cat + (nsite - i - 1) * ncat;
        for (j = 0; j < ncat; j++) {
            work[j] = transit_arr[0] + pre_work[0];
            next_cat_arr[j] = 0;
            for (l = 1; l < ncat; l++) {
                v = transit_arr[l] + pre_work[l];
                if (work[j] < v) {
                    work[j] = v;
                    next_cat_arr[j] = (int) l;
                }
            }
            work[j] += site_lh_arr[j];
            transit_arr += ncat;
        }
        pre_k = k;
        pre_work = work;
    }
    
    // get the start with the max likelihood
    double max_log_like = prob_log[0] + pre_work[0];
    int max_cat = 0;
    for (j = 1; j < ncat; j++) {
        v = prob_log[j] + pre_work[j];
        if (max_log_like < v) {
            max_log_like = v;
            max_cat = (int) j;
        }
    }
    pathLogLike = max_log_like;
    
    // show the max log likelihood
    // cout << "The path with maximum log likelihood: " << pathLogLike << endl;

    // get the assignment of the categories along sites with maximum likelihood
    site_categories[0] = max_cat;
    for (i = 0; i < nsite - 1; i++) {
        max_cat = next_cat[i * ncat + max_cat];
        site_categories[i+1] = max_cat;
    }
    
    pathLogLike = max_log_like;
    
    return max_log_like;
}

// optimize probabilities using EM algorithm
double PhyloHmm::optimizeProbEM() {
    size_t pre_k = 0;
    size_t k;
    size_t i,j;
    double* pre_work;
    double* work;
    double* site_lh_arr;
    double* transit_arr;
    site_lh_arr = site_like_cat;
    memcpy(work_arr, site_lh_arr, sizeof(double) * ncat);
    pre_work = work_arr;
    for (i = 1; i < nsite; i++) {
        k = pre_k ^ 1;
        work = work_arr + k * ncat;
        site_lh_arr += ncat;
        transit_arr = modelHmm->getTransitLog();
        for (j = 0; j < ncat; j++) {
            work[j] = logDotProd(transit_arr, pre_work, ncat) + site_lh_arr[j];
            transit_arr += ncat;
        }
        pre_k = k;
        pre_work = work;
    }
    
    k = pre_k ^ 1;
    work = work_arr + k * ncat;
    // compute the max among prob_log[0]+work[0],prob_log[1]+work[1],...
    for (j = 0; j < ncat; j++) {
        work[j] = prob_log[j] + pre_work[j];
    }
    double max = work[0];
    int max_j = 0;
    for (j = 1; j < ncat; j++) {
        if (max < work[j]) {
            max = work[j];
            max_j = (int) j;
        }
    }
    // exp(prob_log[i] + work[i] - max)
    for (j = 0; j < max_j; j++) {
        work[j] = exp(work[j] - max);
    }
    work[max_j] = 1.0;
    for (j = max_j+1; j < ncat; j++) {
        work[j] = exp(work[j] - max);
    }
    // compute the sum of them
    double sum_like = 0.0;
    for (j = 0; j < ncat; j++) {
        sum_like += work[j];
    }
    // set the new value of prob[i] = pre_work[i] / sum_like
    sum_like = 1.0 / sum_like;
    for (j = 0; j < ncat; j++) {
        prob[j] = work[j] * sum_like;
    }
    computeLogProb();
    
    return logDotProd(prob_log, pre_work, ncat);
}

// optimize the parameters, including the transition matrix and the probability array
double PhyloHmm::optimizeParameters(double gradient_epsilon) {
    // optimize the transition matrix
    
    double score;
    score = modelHmm->optimizeParameters(gradient_epsilon);
    if (verbose_mode >= VB_MED) {
        cout << "after optimizing the transition matrix, HMM likelihood = " << score << endl;
        cout << "modelHmm->tranSameCat : " << modelHmm->tranSameCat << endl;
    }
    // optimize the probability array
    score = optimizeProbEM();
    if (verbose_mode >= VB_MED) {
        cout << "after optimizing the probability array, HMM likelihood = " << score << endl;
        cout << "probability array :";
        for (size_t i = 0; i < ncat; i++) {
            cout << " " << prob[i];
        }
        cout << endl;
    }
    return score;
}

// show the assignment of the categories along sites with max likelihood
void PhyloHmm::showSiteCatMaxLike(ostream& out) {
    size_t i;
    int* numSites; // number of sites for each category
    double* rateSites; // ratio of sites for each category
    
    out << "The assignment of categories along sites with maximum likelihood" << endl;
    out << "Sites\tCategory" << endl;
    int pre_max_cat = site_categories[0];
    int pre_site = 0;
    for (i=1; i<nsite; i++) {
        if (site_categories[i] != pre_max_cat) {
            out << "[" << pre_site + 1 << "," << i << "]\t" << pre_max_cat+1 << endl;
            pre_max_cat = site_categories[i];
            pre_site = (int) i;
        }
    }
    out << "[" << pre_site + 1 << "," << i << "]\t" << pre_max_cat+1 << endl;
    
    numSites = new int[nsite];
    memset(numSites, 0, sizeof(int) * nsite);
    rateSites = new double[nsite];
    for (i=0; i<nsite; i++) {
        numSites[site_categories[i]]++;
    }
    for (i=0; i<ncat; i++)
        rateSites[i] = (double) numSites[i] / nsite;
    
    // show the statistics
    out << "Number of sites for each category:";
    for (i=0; i<ncat; i++)
        out << " " << numSites[i];
    out << endl;
    
    out << "Ratio of sites for each category:";
    for (i=0; i<ncat; i++)
        out << " " << fixed << setprecision(5) << rateSites[i];
    out << endl << endl;

    // show the max log likelihood
    out << "The path with maximum log likelihood: " << fixed << setprecision(5) << pathLogLike << endl;

    delete[] numSites;
    delete[] rateSites;
}

// compute the log values of prob
void PhyloHmm::computeLogProb() {
    size_t i;
    for (i = 0; i < ncat; i++) {
        prob_log[i] = log(prob[i]);
    }
}

// compute backward log-likelihood
// prerequisite: array site_like_cat has been updated (i.e. computeLogLikelihoodSiteTree() has been invoked)
// and save all the intermediate results to the bwd_array array
double PhyloHmm::computeBackLikeArray() {
    size_t pre_k = 0;
    size_t k;
    double* pre_work;
    double* work;
    double* site_lh_arr;
    double* transit_arr;
    double score;
    site_lh_arr = site_like_cat;
    work = bwd_array + (nsite - 1) * ncat;
    memcpy(work, site_lh_arr, sizeof(double) * ncat);
    for (int i = nsite - 1; i >= 1; i--) {
        pre_work = work;
        work = bwd_array + (i - 1) * ncat;
        transit_arr = modelHmm->getTransitLog();
        site_lh_arr += ncat;
        for (int j = 0; j < ncat; j++) {
            work[j] = logDotProd(transit_arr, pre_work, ncat) + site_lh_arr[j];
            transit_arr += ncat;
        }
    }
    score = logDotProd(prob_log, work, ncat);
    return score;
}

// compute forward log-likelihood
// and save all the intermediate results to the fwd_array array
double PhyloHmm::computeFwdLikeArray() {
    size_t pre_k = 0;
    size_t k;
    double* pre_work;
    double* work;
    double* site_lh_arr;
    double* transit_arr;
    double score;
    site_lh_arr = site_like_cat + (nsite-1) * ncat;
    work = fwd_array;
    memcpy(work, prob_log, sizeof(double) * ncat);
    for (int i = 1; i < nsite; i++) {
        pre_work = work;
        work += ncat;
        transit_arr = modelHmm->getTransitLog();
        for (int j = 0; j < ncat; j++) {
            work[j] = logDotProd(transit_arr, pre_work, ncat) + site_lh_arr[j];
            transit_arr += ncat;
        }
        site_lh_arr -= ncat;
    }
    score = logDotProd(site_lh_arr, work, ncat);
    return score;
}

// verify the backLikeArray and FwdLikeArray
void PhyloHmm::checkEachSiteBackFwdLikeArray() {
    double score;
    double* f_array = fwd_array;
    double* b_array = bwd_array;
    for (int i=0; i<nsite; i++) {
        score = logDotProd(f_array, b_array, ncat);
        cout << i+1 << "\t" << score << endl;
        f_array += ncat;
        b_array += ncat;
    }
}

// compute the marginal probabilities for each site
void PhyloHmm::computeMarginalProb() {
    double score;
    double* f_array = fwd_array;
    double* b_array = bwd_array;
    double* mprob = marginal_prob;
    
    computeBackLikeArray();
    computeFwdLikeArray();

    // cout << "Marginal probabilities:" << endl;
    for (int i=0; i<nsite; i++) {
        // cout << i+1;
        score = logDotProd(f_array, b_array, ncat);
        for (int j=0; j<ncat; j++) {
            mprob[j] = exp(f_array[j]+b_array[j]-score);
            // cout << "\t" << mprob[j];
        }
        // cout << endl;
        f_array += ncat;
        b_array += ncat;
        mprob += ncat;
    }
}

// compute the marginal probabilities for transitions between every pair of sites
void PhyloHmm::computeMarginalTransitProb() {
    double score, sum;
    double* f_array = fwd_array;
    double* b_array = bwd_array + ncat;
    double* catlike_array = site_like_cat + ncat * (nsite - 1);
    double* t_array;
    double* mprob = marginal_tran;
    double* t1 = new double[ncat * ncat];
    double* t2 = new double[ncat * ncat];
    int sq_ncat = ncat * ncat;
    int k;
    
    computeBackLikeArray();
    computeFwdLikeArray();

    // cout << "Marginal transition probabilities:" << endl;
    for (int i=0; i<nsite-1; i++) {
        // cout << "(" << i+1 << "," << i+2 << ")";
        t_array = modelHmm->getTransitLog();
        k = 0;
        for (int j1=0; j1<ncat; j1++) {
            for (int j2=0; j2<ncat; j2++) {
                t1[k] = f_array[j1] + catlike_array[j1];
                t2[k] = b_array[j2] + t_array[k];
                k++;
            }
        }
        score = logDotProd(t1, t2, sq_ncat);
        // cout << "[" << score << "]";
        sum = 0.0;
        for (k=0; k<sq_ncat; k++) {
            mprob[k] = exp(t1[k] + t2[k] - score);
            sum += mprob[k];
            // cout << " " << mprob[k];
        }
        // cout << " {" << sum << "}" << endl;
        f_array += ncat;
        b_array += ncat;
        catlike_array -= ncat;
        mprob += ncat * ncat;
    }
    delete[] t1;
    delete[] t2;
}

void PhyloHmm::showSiteLikeCat() {
    int k = 0;
    cout << "Array site_like_cat :" << endl;
    for (int i=0; i<nsite; i++) {
        cout << nsite - i;
        for (int j=0; j<ncat; j++) {
            cout << "\t" << site_like_cat[k++];
        }
        cout << endl;
    }
}

void PhyloHmm::showTransiteLog() {
    double* transit_arr = modelHmm->getTransitLog();
    int k = 0;
    cout << "Array TransiteLog :" << endl;
    for (int i=0; i<ncat; i++) {
        for (int j=0; j<ncat; j++) {
            if (j>0)
                cout << "\t";
            cout << transit_arr[k++];
        }
        cout << endl;
    }
}


