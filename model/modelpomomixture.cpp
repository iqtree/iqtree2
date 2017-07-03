//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"
#include "rategamma.h"

ModelPoMoMixture::ModelPoMoMixture(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     string pomo_params, string pomo_rate_str)
	:
        ModelMarkov(tree),
        ModelPoMo(model_name, model_params, freq_type, freq_params, tree, pomo_params),
        ModelMixture(tree)

{
    optimizing_ratehet = false;

    // get number of categories
    int m, num_rate_cats = 4;
    if (pomo_rate_str.length() > 2 && isdigit(pomo_rate_str[2])) {
        int end_pos;
        num_rate_cats = convert_int(pomo_rate_str.substr(2).c_str(), end_pos);
        if (num_rate_cats < 1) outError("Wrong number of rate categories");
    }

    // initialize rate heterogeneity
    ratehet = new RateGamma(num_rate_cats, Params::getInstance().gamma_shape, Params::getInstance().gamma_median, tree);

    // initialize mixture
    prop = aligned_alloc<double>(num_rate_cats);

    // creating mixture components
    for (m = 0; m < num_rate_cats; m++) {
        ModelMarkov* model = new ModelMarkov(tree);
        model->init(FREQ_USER_DEFINED);
//        model->total_num_subst = ratehet->getRate(m);
        push_back(model);
        prop[m] = ratehet->getProp(m);
    }

    // allocate memory for mixture components so that they are continuous in RAM
    initMem();

    // TODO: why calling this init here
    ModelMarkov::init(FREQ_USER_DEFINED);

}




ModelPoMoMixture::~ModelPoMoMixture() {
}

void ModelPoMoMixture::saveCheckpoint() {
    ModelPoMo::saveCheckpoint();
}

void ModelPoMoMixture::restoreCheckpoint() {
    ModelPoMo::restoreCheckpoint();
}


int ModelPoMoMixture::getNDim() {
    if (optimizing_ratehet)
        return ratehet->getNDim();
    return ModelPoMo::getNDim();
}


int ModelPoMoMixture::getNDimFreq() {
    return ModelPoMo::getNDimFreq();
}


double ModelPoMoMixture::targetFunk(double x[]) {
    if (optimizing_ratehet) {
        getVariables(x);
        phylo_tree->clearAllPartialLH();
        return -phylo_tree->computeLikelihood();
    }
    return ModelPoMo::targetFunk(x);
}



void ModelPoMoMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    if (optimizing_ratehet) {
//        ratehet->setBounds(lower_bound, upper_bound, bound_check);
        lower_bound[1] = 0.2;
        upper_bound[1] = 100.0;
        bound_check[1] = false;
        return;
    }
    ModelPoMo::setBounds(lower_bound, upper_bound, bound_check);
}


void ModelPoMoMixture::writeInfo(ostream &out) {
    ModelPoMo::writeInfo(out);
}

void ModelPoMoMixture::decomposeRateMatrix() {
    // propagate eigenvalues and eigenvectors
    int m, nmix = getNMixtures(), num_states_2 = num_states*num_states;
    double saved_mutation_rate_matrix[n_alleles*n_alleles];
    memcpy(saved_mutation_rate_matrix, mutation_rate_matrix, sizeof(double)*n_alleles*n_alleles);

    // trick: reverse loop to retain eigenvalues and eigenvectors of the 0th mixture class
    for (m = nmix-1; m >= 0; m--) {
        // rescale mutation_rates
        scaleMutationRatesAndUpdateRateMatrix(ratehet->getRate(m));
        ModelPoMo::decomposeRateMatrix();
        // copy eigenvalues and eigenvectors
        if (m > 0) {
            memcpy(eigenvalues+m*num_states, eigenvalues, sizeof(double)*num_states);
            memcpy(eigenvectors+m*num_states_2, eigenvectors, sizeof(double)*num_states_2);
            memcpy(inv_eigenvectors+m*num_states_2, inv_eigenvectors, sizeof(double)*num_states_2);
        }
        // restore mutation_rate matrix
        memcpy(mutation_rate_matrix, saved_mutation_rate_matrix, sizeof(double)*n_alleles*n_alleles);
    }
    updatePoMoStatesAndRateMatrix();
}


void ModelPoMoMixture::setVariables(double *variables) {
    if (optimizing_ratehet) {
        ratehet->setVariables(variables);
        return;
    }
    ModelPoMo::setVariables(variables);
}


bool ModelPoMoMixture::getVariables(double *variables) {
    if (optimizing_ratehet) {
        bool changed = ratehet->getVariables(variables);
        if (changed) {
            decomposeRateMatrix();
        }
        return changed;
    }
    return ModelPoMo::getVariables(variables);
}


double ModelPoMoMixture::optimizeParameters(double gradient_epsilon) {

    // first optimize pomo model parameters
    double score = ModelPoMo::optimizeParameters(gradient_epsilon);

    // then optimize rate heterogeneity
    if (ratehet->getNDim() > 0) {
        optimizing_ratehet = true;
        double score_ratehet = ModelPoMo::optimizeParameters(gradient_epsilon);
        ratehet->writeInfo(cout);
        optimizing_ratehet = false;
        ASSERT(score_ratehet >= score-0.1);
        return score_ratehet;
    }
    return score;
}

void reportRate(ostream &out, PhyloTree &tree);

void ModelPoMoMixture::report(ostream &out) {
    ModelPoMo::report(out);
    RateHeterogeneity *saved_rate = phylo_tree->getRate();
    phylo_tree->setRate(ratehet);
    reportRate(out, *phylo_tree);
    phylo_tree->setRate(saved_rate);
}

bool ModelPoMoMixture::isUnstableParameters() {
    if (ModelPoMo::isUnstableParameters())
        return true;
    if (ModelMixture::isUnstableParameters())
        return true;
    return false;
}
