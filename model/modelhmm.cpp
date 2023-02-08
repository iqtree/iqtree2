//
//  modelhmm.h
//  model
//
//  Created by Thomas Wong on 02/02/23.
//

#include "modelhmm.h"

ModelHmm::ModelHmm(int numcat) {
    ncat = numcat;
    // memory allocation of the array
    size_t transit_size = get_safe_upper_limit(ncat * ncat);
    transitLog = aligned_alloc<double>(transit_size);
    
    if (getName() == "SM")
        initialize_param();
}

/**
 destructor
 */
ModelHmm::~ModelHmm() {
    aligned_free(transitLog);
}

// initialize parameters
void ModelHmm::initialize_param() {
    // initialization of the transition probability between the same category
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
    double lower_bound = MIN_TRAN_PROB;
    double upper_bound = 1.0 - MIN_TRAN_PROB;
    
    // optimization of transition matrix
    double negative_lh;
    double ferror,optx;
    optx = minimizeOneDimen(lower_bound, tranSameCat, upper_bound, gradient_epsilon, &negative_lh, &ferror);
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
