//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"
#include "rategamma.h"
#include "utils/tools.h"

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
    opt_mode = OPT_NONE;

    // get number of categories
    int m, num_rate_cats = 4;
    if (pomo_rate_str.length() > 2 && isdigit(pomo_rate_str[2])) {
        int end_pos;
        num_rate_cats = convert_int(pomo_rate_str.substr(2).c_str(), end_pos);
        if (num_rate_cats < 1) outError("Wrong number of rate categories");
    }

    // initialize rate heterogeneity
    ratehet = new RateGamma(num_rate_cats, Params::getInstance().gamma_shape, Params::getInstance().gamma_median, tree);

    // Adjust name.
    // this->name += pomo_rate_str;
    // this->full_name += " Gamma rate heterogeneity with " + convertIntToString(num_rate_cats) + " components;";
    this->name += ratehet->name;
    this->full_name += ratehet->full_name;

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

    // TODO: why calling this here?
    ModelMarkov::init(freq_type);

}


string ModelPoMoMixture::getName() {
    return ModelPoMo::getName();
}

ModelPoMoMixture::~ModelPoMoMixture() {
}

void ModelPoMoMixture::setCheckpoint(Checkpoint *checkpoint) {
	ModelPoMo::setCheckpoint(checkpoint);
    ratehet->setCheckpoint(checkpoint);
}

void ModelPoMoMixture::startCheckpoint() {
    checkpoint->startStruct("ModelPoMoMixture");
}

void ModelPoMoMixture::saveCheckpoint() {
    ModelPoMo::saveCheckpoint();
    startCheckpoint();
    ratehet->saveCheckpoint();
    endCheckpoint();
}

void ModelPoMoMixture::restoreCheckpoint() {
    // ratehet needs to be restored first, so that decomposeRateMatrix works properly
    startCheckpoint();
    ratehet->restoreCheckpoint();
    endCheckpoint();
    ModelPoMo::restoreCheckpoint();
}


int ModelPoMoMixture::getNDim() {
    if (opt_mode == OPT_RATEHET)
        return ratehet->getNDim();
    else if (opt_mode == OPT_POMO)
        return ModelPoMo::getNDim();
    else return ratehet->getNDim()+ModelPoMo::getNDim();
}


int ModelPoMoMixture::getNDimFreq() {
    return ModelPoMo::getNDimFreq();
}


double ModelPoMoMixture::targetFunk(double x[]) {
    if (opt_mode == OPT_RATEHET) {
        getVariables(x);
        phylo_tree->clearAllPartialLH();
        return -phylo_tree->computeLikelihood();
    }
    return ModelPoMo::targetFunk(x);
}



void ModelPoMoMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    if (opt_mode == OPT_RATEHET) {
//        ratehet->setBounds(lower_bound, upper_bound, bound_check);
        lower_bound[1] = max(POMO_GAMMA_MIN, Params::getInstance().min_gamma_shape);
        upper_bound[1] = POMO_GAMMA_MAX;
        // Boundary checking is the preferred solution to warn the user if the
        // shape parameter hits the boundary, but it seems to be too verbose.
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
        setScale(ratehet->getRate(m));
        ModelPoMo::decomposeRateMatrix();
        // TODO Check! TEST: copy state frequency
        ModelPoMo::getStateFrequency(at(m)->state_freq);
        // copy eigenvalues and eigenvectors
        if (m > 0) {
            memcpy(eigenvalues+m*num_states, eigenvalues, sizeof(double)*num_states);
            memcpy(eigenvectors+m*num_states_2, eigenvectors, sizeof(double)*num_states_2);
            memcpy(inv_eigenvectors+m*num_states_2, inv_eigenvectors, sizeof(double)*num_states_2);
            memcpy(inv_eigenvectors_transposed+m*num_states_2
                   , inv_eigenvectors_transposed, sizeof(double)*num_states_2);
        }
        // restore mutation_rate matrix
        memcpy(mutation_rate_matrix, saved_mutation_rate_matrix, sizeof(double)*n_alleles*n_alleles);
    }
    // // Reset scale.
    setScale(1.0);
    updatePoMoStatesAndRateMatrix();
    ModelPoMo::getStateFrequency(state_freq);
}


void ModelPoMoMixture::setVariables(double *variables) {
    if (opt_mode == OPT_RATEHET) {
        ratehet->setVariables(variables);
        return;
    }
    ModelPoMo::setVariables(variables);
}


bool ModelPoMoMixture::getVariables(double *variables) {
    if (opt_mode == OPT_RATEHET) {
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
    opt_mode = OPT_POMO;
    double score = ModelPoMo::optimizeParameters(gradient_epsilon);
    opt_mode = OPT_NONE;

    // then optimize rate heterogeneity
    if (ratehet->getNDim() > 0) {
        opt_mode = OPT_RATEHET;
        double score_ratehet = ModelPoMo::optimizeParameters(gradient_epsilon);
        if (verbose_mode >= VB_MIN) {
          double shape = ratehet->getGammaShape();
          if (shape <= POMO_GAMMA_MIN)
            outWarning("The shape parameter of the gamma rate heterogeneity is hitting the lower boundary.");
          ratehet->writeInfo(cout);
        }
        opt_mode = OPT_NONE;
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

int ModelPoMoMixture::get_num_states_total() {
  // We assume that all mixture model components have the same number of states.
  return num_states * getNMixtures();
}

void ModelPoMoMixture::update_eigen_pointers(double *eval, double *evec
                                             , double *inv_evec, double* inv_evec_transposed) {
    eigenvalues = eval;
    eigenvectors = evec;
    inv_eigenvectors = inv_evec;
    inv_eigenvectors_transposed = inv_evec_transposed;
    
    // We assume that all mixture model components have the same number of states.
    size_t rowOffset = 0;
    size_t matrixOffset = 0; //into matrices
    size_t num_states_squared = num_states * num_states;
    
    for (iterator it = begin(); it != end();
         it++, rowOffset+=num_states, matrixOffset+=num_states_squared) {
        (*it)->update_eigen_pointers(eval + rowOffset,
                                     evec + matrixOffset,
                                     inv_evec + matrixOffset,
                                     inv_evec_transposed + matrixOffset);
    }
    return;
}

bool ModelPoMoMixture::isUnstableParameters() {
    if (ModelPoMo::isUnstableParameters())
        return true;
    if (ModelMixture::isUnstableParameters())
        return true;
    return false;
}

// I had to write this function because of a compiler error. ModelPoMoMixture is
// inheriting functions from ModelMixture and from ModelPoMo. I defined
// computeTransMatrix for ModelPoMo because I thought that the Modelmarkov
// version did not work for non-reversible substitution models. However, this
// led to a clash because then computeTransMatrix is defined in both,
// ModelMixture and ModelPoMo and inheritance is flawed.
void ModelPoMoMixture::computeTransMatrix(double time, double *trans_matrix, int mixture, int selected_row) {
  ASSERT(mixture < getNMixtures());
  at(mixture)->computeTransMatrix(time, trans_matrix, 0, selected_row);
}
