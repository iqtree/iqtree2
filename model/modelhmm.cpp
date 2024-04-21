//
//  modelhmm.h
//  model
//
//  Created by Thomas Wong on 02/02/23.
//

#include "modelhmm.h"

ModelHmm::ModelHmm(int numcat) {
    ncat = numcat;
    sq_ncat = ncat * ncat;
    transitLog = NULL;
}

/**
 destructor
 */
ModelHmm::~ModelHmm() {
    if (transitLog != NULL)
        aligned_free(transitLog);
}

// initialize transitLog array
void ModelHmm::initialize_transitLog() {
    // memory allocation of the array
    size_t transit_size = get_safe_upper_limit(ncat * ncat);
    if (transitLog != NULL)
        aligned_free(transitLog);
    transitLog = aligned_alloc<double>(transit_size);
    tranSameCat = INITIAL_PROB_SAME_CAT;
    // compute the log values of transition matrix
    computeLogTransits();
}

/**
    set the associated PhyloHmm
    @param phyloHmm the associated PhyloHmm
*/
void ModelHmm::setPhyloHmm(PhyloHmm *phyloHmm) {
    phylo_hmm = phyloHmm;
}

/**
 Optimize the model parameters
 @param gradient_epsilon: epsilon for optimization
 @return log-likelihood value
 */
double ModelHmm::optimizeParameters(double gradient_epsilon) {
    
    // changed to EM algorithm
    return optimizeParametersByEM();
    
    // use BFGS algorithm
    double lower_bound;
    if (Params::getInstance().HMM_min_stran > MIN_TRAN_PROB) {
        lower_bound = Params::getInstance().HMM_min_stran;
    } else {
        lower_bound = MIN_TRAN_PROB;
    }
    double upper_bound = 1.0 - MIN_TRAN_PROB;
    
    // optimization of transition matrix
    double negative_lh;
    double ferror,optx;
    optx = minimizeOneDimen(lower_bound, tranSameCat, upper_bound, gradient_epsilon, &negative_lh, &ferror);
    return -computeFunction(optx);
}

/**
 Optimize the transition matrix by EM algorithm
 @return log-likelihood value
 */
double ModelHmm::optimizeParametersByEM() {
    // Optimize the transition matrix by EM algorithm
    double* marginalTran;
    int sq_ncat = ncat * ncat;
    double optx;
    int i,j,k;
    
    // compute the expected value of transition probability between the same categories
    phylo_hmm->computeMarginalTransitProb();
    marginalTran = phylo_hmm->marginal_tran;
    optx = 0.0;
    for (k=0; k<phylo_hmm->nsite-1; k++) {
        // only consider the transition probs between the same categories
        j = 0;
        for (i=0; i<ncat; i++) {
            optx += marginalTran[j];
            j += (ncat + 1);
        }
        marginalTran += sq_ncat;
    }
    optx = optx / (double)(phylo_hmm->nsite-1);
    if (optx < MIN_TRAN_PROB)
        optx = MIN_TRAN_PROB;
    if (optx < Params::getInstance().HMM_min_stran)
        optx = Params::getInstance().HMM_min_stran;
    if (optx > 1.0 - MIN_TRAN_PROB)
        optx = 1.0 - MIN_TRAN_PROB;
    return -computeFunction(optx);
}

/**
 Show parameters
 */
void ModelHmm::showParameters(ostream& out) {
    size_t i, j;
    double tranDiffCat = (1.0 - tranSameCat)/(ncat-1);
    out << "Estimated HMM transition matrix :" << endl;
    for (i = 0; i < ncat; i++) {
        for (j = 0; j < ncat; j++) {
            if (j > 0)
                out << "\t";
            if (i == j)
                out << fixed << setprecision(5) << tranSameCat;
            else
                out << fixed << setprecision(5) << tranDiffCat;
        }
        out << endl;
    }
}

int ModelHmm::getNParameters() {
    return 1;
}

// compute the log values of transition matrix
void ModelHmm::computeLogTransits() {
    double logTranSameCat = log(tranSameCat);
    double logTranDiffCat = log((1.0 - tranSameCat)/(ncat-1));
    for (size_t i = 0; i < ncat * ncat; i++)
        transitLog[i] = logTranDiffCat;
    for (size_t i = 0; i < ncat; i++)
        transitLog[i*ncat + i] = logTranSameCat;
}

// for optimization
double ModelHmm::computeFunction(double tran_same_cat) {
    tranSameCat = tran_same_cat;
    computeLogTransits();
    return -phylo_hmm->computeBackLike();
}
