//
//  ModelHmmTm.h
//  model
//
//  Created by Thomas Wong on 13/10/23.
//

#include "modelhmmtm.h"

/**
 constructor
 @param ncat : number of categories
 @param ntype : number of types
 */
ModelHmmTm::ModelHmmTm(int numcat, int numtype, int* siteTypes, int numsite, string* typeDescs) : ModelHmm(numcat) {
    ntype = numtype;
    ntypepair = (numtype + 1) * numtype / 2;
    nsite = numsite;
    size_t transit_size = get_safe_upper_limit(numcat * numcat * ntypepair);
    transit = aligned_alloc<double>(transit_size);
    transit_normalize = aligned_alloc<double>(transit_size);
    tranSameCats = aligned_alloc<double>(get_safe_upper_limit(ntypepair));

    transit_fwd = NULL; // site i-1 -> site i
    // create the mapping between site ID and which matrix should be used for the transition between the previous site and the current site
    initialize_transit_id(siteTypes, nsite, typeDescs);
}

/**
 destructor
 */
ModelHmmTm::~ModelHmmTm() {
    if (transit != NULL)
        aligned_free(transit);
    if (transit_normalize != NULL)
        aligned_free(transit_normalize);
    if (transit_fwd != NULL)
        aligned_free(transit_fwd);
    if (tranSameCats != NULL)
        aligned_free(tranSameCats);
}

// initialize transitLog array
void ModelHmmTm::initialize_transitLog() {
    size_t i, j, k;

    // memory allocation of the array
    size_t transit_size = get_safe_upper_limit(ncat * ncat * ntypepair);
    if (transitLog != NULL)
        aligned_free(transitLog);
    transitLog = aligned_alloc<double>(transit_size);

    for (k = 0; k < ntypepair; k++) {
        tranSameCats[k] = INITIAL_PROB_SAME_CAT;
    }
    updateTransits();
}

// create the mapping between site ID and which matrix should be used for the transition between the previous site and the current site
void ModelHmmTm::initialize_transit_id(int* siteTypes, int nsite, string* typeDescs) {
    int* typepair2id = new int[ntype * ntype];
    int i, j, k;
    
    transit_id_descs.clear();
    // compute the typepair -> ID
    k = 0;
    for (i = 0; i < ntype; i++) {
        // j = i
        typepair2id[i * ntype + i] = k;
        transit_id_descs.push_back(typeDescs[i] + " <-> " + typeDescs[i]);
        k++;
        // j > i
        for (j = i+1; j < ntype; j++) {
            typepair2id[i * ntype + j] = k;
            transit_id_descs.push_back(typeDescs[i] + " <-> " + typeDescs[j]);
            typepair2id[j * ntype + i] = k;
            k++;
        }
    }
    
    size_t array_size = get_safe_upper_limit(nsite);
    if (transit_fwd != NULL)
        aligned_free(transit_fwd);
    
    // which transit matrix should be used for each site
    transit_fwd = aligned_alloc<int>(array_size); // site i-1 -> site i
    int curr_type;
    int transit_id;
    int pre_type;
    // forward
    pre_type = siteTypes[0];
    transit_fwd[0] = -1; // not applicable
    for (i = 1; i < nsite; i++) {
        curr_type = siteTypes[i];
        transit_id = typepair2id[pre_type * ntype + curr_type];
        transit_fwd[i] = transit_id;
        pre_type = curr_type;
    }
}

/**
 Optimize the model parameters
 @param gradient_epsilon: epsilon for optimization
 @return log-likelihood value
 */
double ModelHmmTm::optimizeParameters(double gradient_epsilon) {

    // changed to EM algorithm
    return optimizeParametersByEM();

    // use BFGS algorithm
    int ndim = ntypepair;
    double *variables = new double[ndim + 1];
    double *upper_bound = new double[ndim + 1];
    double *lower_bound = new double[ndim + 1];
    bool *bound_check = new bool[ndim + 1];
    double logLike;
    
    // optimization of transition matrix
    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    logLike = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
    getVariables(variables);
    
    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    return logLike;
}

/**
 Optimize the transition matrix by EM algorithm
 @return log-likelihood value
 */
double ModelHmmTm::optimizeParametersByEM() {
    // Optimize the transition matrix by EM algorithm
    double* marginalTran;
    int dim = phylo_hmm->nsite-1;
    int i,j,k;
    
    // compute the expected value of transition probability between the same categories for all same typepair
    phylo_hmm->computeMarginalTransitProb();
    marginalTran = phylo_hmm->marginal_tran;
    
    // create tmp arrays
    double* tmpProbs = new double[ntypepair];
    int* tmpNums = new int[ntypepair];
    memset(tmpProbs, 0, sizeof(double) * ntypepair);
    memset(tmpNums, 0, sizeof(int) * ntypepair);
    for (k=1; k<nsite; k++) {
        int currTypePair = transit_fwd[k];
        tmpNums[currTypePair]++;
        /*
        // show marginalTran array
        cout << "k=" << k << endl;
        int kk=0;
        for (int ii = 0; ii < ncat; ii++) {
            for (int jj = 0; jj < ncat; jj++) {
                if (jj > 0)
                    cout << "\t";
                cout << marginalTran[kk];
                kk++;
            }
            cout << endl;
        }
        */
        // only consider the transition between the same categories
        j = 0;
        for (i=0; i<ncat; i++) {
            tmpProbs[currTypePair] += marginalTran[j];
            j += (ncat + 1);
        }
        marginalTran += sq_ncat;
    }
    for (k=0; k<ntypepair; k++) {
        if (tmpNums[k] > 0) {
            tranSameCats[k] = tmpProbs[k] / tmpNums[k];
            // cout << "tranSameCats[" << k << "]=" << tranSameCats[k] << endl;
            if (tranSameCats[k] < MIN_TRAN_PROB)
                tranSameCats[k] = MIN_TRAN_PROB;
            if (tranSameCats[k] < Params::getInstance().HMM_min_stran)
                tranSameCats[k] = Params::getInstance().HMM_min_stran;
            if (tranSameCats[k] > 1.0 - MIN_TRAN_PROB)
                tranSameCats[k] = 1.0 - MIN_TRAN_PROB;
        }
    }

    // release the memory
    delete[] tmpProbs;
    delete[] tmpNums;

    updateTransits();
    return phylo_hmm->computeBackLike();
}

/**
 Show parameters
 */
void ModelHmmTm::showParameters(ostream& out) {
    size_t i, j, k, l;
    double* curr_transit_norm = transit_normalize;
    for (l = 0; l < ntypepair; l++) {
        out << "Estimated HMM transition matrix (" << transit_id_descs[l] << "):" << endl;
        k = 0;
        for (i = 0; i < ncat; i++) {
            for (j = 0; j < ncat; j++) {
                if (j > 0)
                    out << "\t";
                out << fixed << setprecision(5) << curr_transit_norm[k];
                k++;
            }
            out << endl;
        }
        curr_transit_norm += sq_ncat;
    }
}

int ModelHmmTm::getNParameters() {
    return ntypepair;
}

void ModelHmmTm::setVariables(double *variables) {
    double* var = variables + 1;
    memcpy(var, tranSameCats, sizeof(double) * ntypepair);
}

void ModelHmmTm::getVariables(double *variables) {
    double* var = variables + 1;
    memcpy(tranSameCats, var, sizeof(double) * ntypepair);
    updateTransits();
}

void ModelHmmTm::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {

    size_t i;

    for (i = 1; i <= ntypepair; i++) {
        lower_bound[i] = MIN_VALUE;
        upper_bound[i] = MAX_VALUE;
        bound_check[i] = false;
    }
}

// update the transition matrix according to the array tranSameCats
void ModelHmmTm::updateTransits() {
    int i,j,k;
    double* curr_transit = transit;
    // the other probabilities are set as equally distributed
    for (k = 0; k < ntypepair; k++) {
        double other_tran = (1.0 - tranSameCats[k]) / ((double) ncat - 1.0);
        for (i = 0; i < ncat; i++) {
            for (j = 0; j < i; j++)
                curr_transit[i * ncat + j] = other_tran;
            curr_transit[i * ncat + i] = tranSameCats[k];
            for (j = i+1; j < ncat; j++)
                curr_transit[i * ncat + j] = other_tran;
        }
        curr_transit += sq_ncat;
    }
    computeNormalizedTransits();
    computeLogTransits();
}

// compute the normalized values of transition matrix
void ModelHmmTm::computeNormalizedTransits() {
    double sum;
    size_t i,j,k,l;
    double* curr_transit = transit;
    double* curr_transit_normalize = transit_normalize;
    for (l = 0; l < ntypepair; l++) {
        for (i = 0; i < ncat; i++) {
            sum = 0.0;
            k = i*ncat;
            for (j = 0; j < ncat; j++) {
                sum += curr_transit[k + j];
            }
            for (j = 0; j < ncat; j++) {
                curr_transit_normalize[k + j] = curr_transit[k + j] / sum;
            }
        }
        curr_transit += sq_ncat;
        curr_transit_normalize += sq_ncat;
    }
}

// compute the log values of transition matrix
void ModelHmmTm::computeLogTransits() {
    size_t i, k;
    double* curr_transit_log = transitLog;
    double* curr_transit_normalize = transit_normalize;
    for (k = 0; k < ntypepair; k++) {
        for (i = 0; i < ncat * ncat; i++)
            curr_transit_log[i] = log(curr_transit_normalize[i]);
        curr_transit_log += sq_ncat;
        curr_transit_normalize += sq_ncat;
    }
}

/**
 * @return log values of transition matrix
 */
double* ModelHmmTm::getTransitLog(int site_i) {
    if (transit_fwd == NULL || transit_fwd[site_i] == -1)
        return NULL;
    return transitLog + transit_fwd[site_i] * sq_ncat;
}
