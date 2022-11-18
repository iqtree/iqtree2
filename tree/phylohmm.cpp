//
//  phylohmm.cpp
//  iqtree
//
//  Created by Thomas Wong on 9/11/2022.
//

#include "phylohmm.h"
#include "utils/tools.h"

// constructor
PhyloHMM::PhyloHMM(PhyloTree* phylotree) {
    
    this->phylotree = phylotree;
    
    // get the appropriate type of the categories
    if (phylotree->isTreeMix()) {
        // tree mixture model
        wsl = WSL_TMIXTURE;
    } else if (phylotree->getModel()->isMixture()) {
        // mixture model
        wsl = WSL_MIXTURE;
    } else {
        // non-mixture model
        wsl = WSL_RATECAT;
    }
    
    ncat = phylotree->getNumLhCat(wsl);
    nsite = phylotree->getAlnNSite();
    nptn = phylotree->aln->getNPattern();

    ASSERT(ncat > 1);

    // allocate memory for the arrays
    size_t transit_size = get_safe_upper_limit(ncat * ncat);
    size_t prob_size = get_safe_upper_limit(ncat);
    size_t site_like_cat_size = get_safe_upper_limit(nsite) * ncat;
    size_t ptn_like_cat_size = get_safe_upper_limit(nptn) * ncat;
    size_t ptn_like_size = get_safe_upper_limit(nptn);
    transit = aligned_alloc<double>(transit_size);
    transit_normalize = aligned_alloc<double>(transit_size);
    transit_log = aligned_alloc<double>(transit_size);
    prob = aligned_alloc<double>(prob_size);
    prob_log = aligned_alloc<double>(prob_size);
    site_like_cat = aligned_alloc<double>(site_like_cat_size);
    ptn_like_cat = aligned_alloc<double>(ptn_like_cat_size);
    ptn_like = aligned_alloc<double>(ptn_like_size);
    work_arr = aligned_alloc<double>(prob_size * 2);
    next_cat = aligned_alloc<int>(site_like_cat_size);
    site_categories = aligned_alloc<int>(get_safe_upper_limit(nsite));
    
    // initialize the parameters
    initializeParams();
}

// destructor
PhyloHMM::~PhyloHMM() {
    aligned_free(transit);
    aligned_free(transit_normalize);
    aligned_free(transit_log);
    aligned_free(prob);
    aligned_free(prob_log);
    aligned_free(site_like_cat);
    aligned_free(ptn_like_cat);
    aligned_free(ptn_like);
    aligned_free(work_arr);
    aligned_free(next_cat);
    aligned_free(site_categories);
}

// perform HMM analysis
// and report the results to a file
void PhyloHMM::performHMMAnalysis(const char* filename) {
    
    cout << "Perform HMM Analysis" << endl;

    // load the site likelihood values from the tree
    loadSiteLikeCat();

    // optimize parameters
    optimizeParams();

    // backward log-likelihood
    computeBackLike();
    
    // compute the path with maximum log-likelihood
    computeMaxPath();

    // show the values of the estimated parameters
    showParameters(cout);

    // report the results to a file
    printResults(filename);
}

// load the site likelihood values
void PhyloHMM::loadSiteLikeCat() {
    double* site_lh_arr;
    double* ptn_lh_arr;
    size_t i, j, k, ptn;

    phylotree->computePatternLikelihood(ptn_like, NULL, ptn_like_cat, wsl);
    // the array is ordered from site n to site 1
    site_lh_arr = site_like_cat;
    k = nsite;
    for (i = 0; i < nsite; ++i) {
        --k;
        ptn = phylotree->aln->getPatternID(k);
        ptn_lh_arr = ptn_like_cat + ptn*ncat;
        site_lh_arr = site_like_cat + i*ncat;
        for (j = 0; j < ncat; j++) {
            site_lh_arr[j] = ptn_lh_arr[j];
        }
    }
}

// initialize the parameters
void PhyloHMM::initializeParams() {
    // the transition probability going to the same state is initialized to 0.95
    // and the other probabilities are initialized as equally distributed
    size_t i, j;
    double init_same_tran = 0.95;
    double init_other_tran = (1.0 - init_same_tran) / ((double) ncat - 1.0);
    double init_prob_value = 1.0/((double)ncat);
    for (i = 0; i < ncat; i++) {
        for (j = 0; j < i; j++)
            transit[i * ncat + j] = init_other_tran;
        transit[i * ncat + i] = init_same_tran;
        for (j = i+1; j < ncat; j++)
            transit[i * ncat + j] = init_other_tran;
    }
    for (i = 0; i < ncat; i++)
        prob[i] = init_prob_value;
    computeLogTransits();
    computeLogProb();
}

// compute the log values of transition matrix
void PhyloHMM::computeLogTransits() {
    double* tr;
    double* tr_n;
    double sum;
    size_t i, j;
    for (i = 0; i < ncat; i++) {
        sum = 0.0;
        tr = transit + i * ncat;
        for (j = 0; j < ncat; j++) {
            sum += tr[j];
        }
        sum = 1.0 / sum;
        tr_n = transit_normalize + i * ncat;
        for (j = 0; j < ncat; j++) {
            tr_n[j] = tr[j] * sum;
        }
    }
    for (i = 0; i < ncat * ncat; i++) {
        transit_log[i] = log(transit_normalize[i]);
    }
}

// compute the log values of prob
void PhyloHMM::computeLogProb() {
    size_t i;
    for (i = 0; i < ncat; i++) {
        prob_log[i] = log(prob[i]);
    }
}

// compute the log of dotproduct of the logorithm arrays
double PhyloHMM::logDotProd(double* ln_x, double* ln_y, int n) {
    double max;
    double max_i;
    size_t i;
    double* w;
    double ans;
    
    w = new double[n];
    for (i = 0; i < n; i++) {
        w[i] = ln_x[i] + ln_y[i];
    }
    // find the max
    max = w[0];
    max_i = 0;
    for (i = 1; i < n; i++) {
        if (max < w[i]) {
            max = w[i];
            max_i = i;
        }
    }
    // compute the dotproduct
    ans = 0.0;
    for (i = 0; i < max_i; i++) {
        ans += exp(w[i] - max);
    }
    ans += 1.0;
    for (i = max_i+1; i < n; i++) {
        ans += exp(w[i] - max);
    }
    delete[] w;
    return log(ans) + max;
}

// backward log-likelihood
double PhyloHMM::computeBackLike() {
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
        transit_arr = transit_log;
        for (j = 0; j < ncat; j++) {
            work[j] = logDotProd(transit_arr, pre_work, ncat) + site_lh_arr[j];
            transit_arr += ncat;
        }
        pre_k = k;
        pre_work = work;
    }
    return logDotProd(prob_log, pre_work, ncat);
}

// path with max log-likelihood
double PhyloHMM::computeMaxPath() {
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
        transit_arr = transit_log;
        next_cat_arr = next_cat + (nsite - i - 1) * ncat;
        for (j = 0; j < ncat; j++) {
            work[j] = transit_arr[0] + pre_work[0];
            next_cat_arr[j] = 0;
            for (l = 1; l < ncat; l++) {
                v = transit_arr[l] + pre_work[l];
                if (work[j] < v) {
                    work[j] = v;
                    next_cat_arr[j] = l;
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
            max_cat = j;
        }
    }
    backLogLike = max_log_like;
    
    // show the max log likelihood
    cout << "The path with maximum log likelihood: " << backLogLike << endl;

    // get the assignment of the categories along sites with maximum likelihood
    site_categories[0] = max_cat;
    for (i = 0; i < nsite - 1; i++) {
        max_cat = next_cat[i * ncat + max_cat];
        site_categories[i+1] = max_cat;
    }
    
    return max_log_like;
}


// optimize parameters
void PhyloHMM::optimizeParams(double gradient_epsilon) {
    
    // the array is only for optimization of transition matrix
    size_t arr_size = (ncat - 1) * ncat + 1;
    double *variables = new double[arr_size];
    double *upper_bound = new double[arr_size];
    double *lower_bound = new double[arr_size];
    bool *bound_check = new bool[arr_size];
    double pre_score, score;
    size_t i = 1;
    int ndim;
    
    score = computeBackLike();
    cout << i++ << ". Initial HMM log-likelihood: " << score << endl;
    
    do {
        pre_score = score;
        
        // First optimization of probabilities
        score = optimizeProbEM();

        // Then optimization of transition matrix
        ndim = getNDim();
        setVariables(variables);
        setBounds(lower_bound, upper_bound, bound_check);
        score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
        getVariables(variables);
        updateTransitMatrix();

        cout << i << ". Current HMM log-likelihood: " << score << endl;
            i++;
    } while (score - pre_score > gradient_epsilon);
    
    backLogLike = score;
}

// ------------------------------------------------------------------
// For optimization
// ------------------------------------------------------------------

// optimize probabilities using EM algorithm
double PhyloHMM::optimizeProbEM() {
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
        transit_arr = transit_log;
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
            max_j = j;
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

void PhyloHMM::updateTransitMatrix() {
    size_t i;
    for (i = 0; i < ncat * ncat; i++) {
        transit[i] = transit_normalize[i];
    }
}

// ------------------------------------------------------------------
// For BFGS (only for optimization of transition matrix)
// ------------------------------------------------------------------

int PhyloHMM::getNDim() {
    // number of parameters for transition matrix = (ncat-1) * ncat
    return (ncat - 1) * ncat;
}

void PhyloHMM::setVariables(double *variables) {
    double* var;
    double* tr;
    size_t i,k;
    k = 1;
    // copy the values from transition matrix
    for (i = 0; i < ncat; i++) {
        var = variables + k;
        tr = transit + i*ncat;
        memcpy(var, tr, sizeof(double) * (ncat - 1));
        k += (ncat - 1);
    }
}

void PhyloHMM::getVariables(double *variables) {
    double* var;
    double* tr;
    size_t i, k;
    k = 1;
    // copy the values to transition matrix
    for (i = 0; i < ncat; i++) {
        var = variables + k;
        tr = transit + i*ncat;
        memcpy(tr, var, sizeof(double) * (ncat - 1));
        k += (ncat - 1);
    }
    computeLogTransits();
}

double PhyloHMM::targetFunk(double x[]) {
    getVariables(x);
    return -computeBackLike();
}

void PhyloHMM::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {

    size_t ndim = getNDim();
    size_t i;

    for (i = 1; i <= ndim; i++) {
        //cout << variables[i] << endl;
        lower_bound[i] = MIN_VALUE;
        upper_bound[i] = MAX_VALUE;
        bound_check[i] = false;
    }
}

// -------------------------------------------------------
// For outputting the results / parameters
// -------------------------------------------------------

// print out the site likelihood values (for debugging)
void PhyloHMM::printSiteLikeCat(ostream& out) {
    size_t i, j, k;
    k=0;
    for (i = 0; i < nsite; i++) {
        for (j = 0; j < ncat; j++) {
            if (j > 0)
                out << "\t";
            out << site_like_cat[k];
            k++;
        }
        out << endl;
    }
}

// print out the HMM estimated parameters
void PhyloHMM::showParameters(ostream& out) {
    size_t i, j;
    out << "Estimated HMM transition matrix :" << endl;
    for (i = 0; i < ncat; i++) {
        for (j = 0; j < ncat; j++) {
            if (j > 0)
                out << "\t";
            out << setprecision(10) << transit[i * ncat + j];
        }
        out << endl;
    }
    out << endl;
    out << "Estimated HMM probabilities :" << endl;
    for (i = 0; i < ncat; i++) {
        if (i > 0)
            out << "\t";
        out << setprecision(10) << prob[i];
    }
    out << endl << endl;
    
    out << "BEST HMM SCORE FOUND :" << backLogLike << endl;
}

// show the assignment of the categories along sites with max likelihood
void PhyloHMM::showSiteCatMaxLike(ostream& out) {
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
            pre_site = i;
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
        out << " " << setprecision(5) << rateSites[i];
    out << endl << endl;

    // show the max log likelihood
    out << "The path with maximum log likelihood: " << backLogLike << endl;

    delete[] numSites;
    delete[] rateSites;
}

// print out all the results to a file
void PhyloHMM::printResults(const char *filename) {
    
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
