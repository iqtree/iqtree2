//
//  phylohmm.h
//  tree
//
//  Created by Thomas Wong on 02/02/23.
//

#ifndef phylohmm_h
#define phylohmm_h

#include "model/modelhmm.h"
#include "utils/optimization.h"

using namespace std;

// compute the log of dotproduct of the logorithm arrays
inline double logDotProd(double* ln_x, double* ln_y, int n) {
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

class ModelHmm;

/**
 HMM phylo model base class
 */
class PhyloHmm {
public:
    
    // constructor
    PhyloHmm();
    PhyloHmm(int n_site, int n_cat);
    
    //destructor
    ~PhyloHmm();
    
    // initialize the transition model (override it for different transition model)
    // this function is invoked inside the constructor
    virtual void initializeTransitModel();
    
    // prerequisite for the following three functions:
    //     array site_like_cat has to be updated

    // compute backward log-likelihood
    double computeBackLike();

    // path with max log-likelihood
    double computeMaxPath();

    // optimize probabilities using EM algorithm
    double optimizeProbEM();
    
    // optimize the parameters, including the probability array and the transition matrix
    double optimizeParameters(double gradient_epsilon);

    // show the assignment of the categories along sites with max likelihood
    void showSiteCatMaxLike(ostream& out);

protected:

    // number of sites
    int nsite;
    
    // number of categories
    int ncat;
    
    // transition model
    ModelHmm* modelHmm;
    
    // Parameter: cat probability
    // prob[i] : probability of a site with cat i
    double* prob;
    // sum_i(prob[i]) = 1.0

    // Parameter: log values of cat probability
    // prob[i] : probability of a site with cat i
    double* prob_log;

    // log-likelihood values for each cat on every site
    // site_like_cat[i*ncat+j] = likelihood of site (nsite-i) for cat j
    double* site_like_cat;

    // assignment of categories along the sites with maximum likelihood value
    int* site_categories;

    // backward log-likelihood
    double backLogLike;
    
    // path log-likelihood
    double pathLogLike;

    // working array for computation of backward probabilities
    double* work_arr;

    // working array for computation of the assignment of category with maximum likelihood value
    // for each cat on every site, the category for next site with the max likelihood value
    int* next_cat;

private:

    // compute the log values of prob
    void computeLogProb();
};
#endif
