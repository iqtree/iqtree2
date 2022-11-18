//
//  phylohmm.h
//  iqtree
//
//  Created by Thomas Wong on 9/11/2022.
//

#ifndef phylohmm_h
#define phylohmm_h

#include "phylotree.h"
#include "utils/optimization.h"

#define MIN_VALUE 1e-5
#define MAX_VALUE 10000
#define GRADIENT_EPN 0.001

class PhyloHMM : public Optimization {
public:
    
    // the phylotree
    PhyloTree* phylotree;
    
    // number of categories
    int ncat;
    
    // number of sites
    int nsite;
    
    // number of patterns
    int nptn;
    
    // type of categories
    SiteLoglType wsl;
    
    // Parameter: transition matrix
    // dimension: ncat x ncat
    // transit[i*ncat+j] : transition probability for a site with cat i to the next site with cat j
    double* transit;
    double* transit_normalize;
    // sum_j(transit_normalize[i*ncat+j]) = 1.0

    // Parameter: log values of transition matrix
    // dimension: ncat x ncat
    // transit[i*ncat+j] : transition probability for a site with cat i to the next site with cat j
    double* transit_log;

    // Parameter: cat probability
    // prob[i] : probability of a site with cat i
    double* prob;
    // double* prob_normalize; (obsoleted)
    // sum_i(prob_normalize[i]) = 1.0

    // Parameter: log values of cat probability
    // prob[i] : probability of a site with cat i
    double* prob_log;

    // log-likelihood values for each cat on every site
    // site_like_cat[i*ncat+j] = likelihood of site (nsite-i) for cat j
    double* site_like_cat;
    
    // log-likelihood values for each cat on every pattern
    double* ptn_like_cat;

    // log-likelihood values for every pattern
    double* ptn_like;
    
    // assignment of categories along the sites with maximum likelihood value
    int* site_categories;

    // backward log-likelihood
    double backLogLike;
    
    // path log-likelihood
    double pathLogLike;

    // constructor
    PhyloHMM(PhyloTree* phylotree);
    
    // destructor
    ~PhyloHMM();
    
    // perform HMM analysis
    // and report the results to a file
    void performHMMAnalysis(const char* filename);

    // load the site likelihood values
    void loadSiteLikeCat();

    // backward log-likelihood
    double computeBackLike();

    // maximum log-likelihood
    double computeMaxPath();

    // optimize parameters
    void optimizeParams(double gradient_epsilon = GRADIENT_EPN);
    
    // optimize probabilities using EM algorithm
    double optimizeProbEM();
    
    // print out the site likelihood values
    void printSiteLikeCat(ostream& out);

    // show the values of the parameters
    void showParameters(ostream& out);
    
    // show the assignment of the categories along sites with max likelihood
    void showSiteCatMaxLike(ostream& out);

    // print out all the results to a file
    void printResults(const char *filename);

private:
    
    // working array for computation of backward probabilities
    double* work_arr;

    // working array for computation of the assignment of category with maximum likelihood value
    // for each cat on every site, the category for next site with the max likelihood value
    int* next_cat;

    // initialize the parameters
    void initializeParams();
    
    // compute the log values of transition matrix
    void computeLogTransits();
    
    // compute the log values of prob
    void computeLogProb();
    
    // compute the log of dotproduct of exponential
    double logDotProd(double* x, double* y, int n);
    
    // ---------------------------------------------
    // For optimization
    // ---------------------------------------------

    // update the transition matrix from the matrix "transit_normalize"
    void updateTransitMatrix();
    
    // ---------------------------------------------
    // For BFGS
    // ---------------------------------------------

    int getNDim();
    void setVariables(double *variables);
    void getVariables(double *variables);
    double targetFunk(double x[]);
    void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
};

#endif /* phylohmm_h */
