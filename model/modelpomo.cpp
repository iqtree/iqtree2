#include "modelpomo.h"
#include "modeldna.h"
#include "modelmixture.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

ModelPoMo::ModelPoMo(PhyloTree *tree) : ModelMarkov(tree) {
}

ModelPoMo::ModelPoMo(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     string pomo_heterozygosity)
    // Set reversibility to true to allocate memory for objects (like
    // eigenvalues) necessary when the model is reversible.  In case
    // the model is not reversible memory for different object has to
    // be allocated separately.  This is done during the
    // initialization in ModelPoMo::init_mutation_model().
    : ModelMarkov(tree, true) {
    init(model_name, model_params, freq_type, freq_params, pomo_heterozygosity);
}

void ModelPoMo::init_mutation_model(const char *model_name,
                                        string model_params,
                                        StateFreqType freq_type,
                                        string freq_params)
{
    // Trick ModelDNA constructor by setting the number of states to 4 (DNA).
    phylo_tree->aln->num_states = n_alleles;
    // This would be ideal but the try and catch statement does not
    // work yet.  I guess, for the moment, the user has to find out
    // what went wrong during DNA model initialization.
    try {
        if (verbose_mode >= VB_MED)
            cout << "Initialize PoMo DNA mutation model." << endl;
        string model_str = model_name;
        if (ModelMarkov::validModelName(model_str))
            mutation_model = ModelMarkov::getModelByName(model_str, phylo_tree, model_params, freq_type, freq_params);
        else
            mutation_model = new ModelDNA(model_name, model_params, freq_type, freq_params, phylo_tree);
    }
    catch (string str) {
        cout << "Error during initialization of the underlying mutation model of PoMo." << endl;
        cout << "PoMo only works with DNA models at the moment." << endl;
        outError(str);
    }

    // Reset the number of states.
    phylo_tree->aln->num_states = num_states;

    // Set reversibility state.
    is_reversible = mutation_model->is_reversible;
    if (!is_reversible)
        setReversible(is_reversible);
}

string ModelPoMo::getName() {
  return this->name;
}

void ModelPoMo::init_sampling_method()
{
    sampling_method = phylo_tree->aln->pomo_sampling_method;
    string sampling_method_str;
    if (sampling_method == SAMPLING_SAMPLED) {
        this->name += "+S";
        sampling_method_str = "Sampled";
    }
    else if (sampling_method == SAMPLING_WEIGHTED_BINOM) {
        this->name += "+WB";
        sampling_method_str = "Weighted binomial";
    }
    else if (sampling_method == SAMPLING_WEIGHTED_HYPER) {
      this->name += "+WH";
      sampling_method_str = "Weighted hypergeometric";
    }
    else outError("Sampling type is not supported.");

    this->full_name =
        "PoMo with N=" + convertIntToString(N) + " and " +
        mutation_model->full_name + " mutation model; " +
        "Sampling method: " + sampling_method_str + "; " +
        convertIntToString(num_states) + " states in total;";
}

void ModelPoMo::init_boundary_frequencies()
{
    // The boundary state frequencies are just the state frequencies from the
    // underlying mutation model.
    freq_boundary_states = mutation_model->state_freq;
    // The empirical boundary state frequencies are computed from the data.
    freq_boundary_states_emp = new double[4];
    estimateEmpiricalBoundaryStateFreqs(freq_boundary_states_emp);
    // The frequency type is the one from the mutation model.
    freq_type = mutation_model->freq_type;

    // Handle frequency type.  This cannot be done by the underlying
    // mutation model because interpretation of polymorphic states is
    // undefined.
    switch (freq_type) {
    case FREQ_EQUAL:
        // '+FQ'.
        for (int i = 0; i < 4; i++)
            freq_boundary_states[i] = 1.0/ (double)n_alleles;
        break;
    case FREQ_ESTIMATE:
        // '+FO'.  Start estimation at empirical frequencies.
        for (int i = 0; i < 4; i++)
            freq_boundary_states[i] = freq_boundary_states_emp[i];
        break;
    case FREQ_EMPIRICAL:
        // '+F'.
        for (int i = 0; i < 4; i++)
            freq_boundary_states[i] = freq_boundary_states_emp[i];
        break;
    case FREQ_USER_DEFINED:
        // '+FU'. ModelDNA should have set them already.
        if (freq_boundary_states[0] == 0.0)
            outError("State frequencies not specified");
        break;
    case FREQ_UNKNOWN:
        outError("No frequency type given.");
        break;
    default:
        outError("Unknown frequency type.");
        break;
    }
}

void ModelPoMo::init_fixed_parameters(string model_params,
                                      string pomo_heterozygosity)
{
    fixed_model_params = false;
    fixed_heterozygosity_emp = false;
    fixed_heterozygosity_usr = false;
    fixed_heterozygosity = false;
    if (model_params.length() > 0)
        // The rest ist done by the underlying mutation model.
        fixed_model_params = true;
    if (pomo_heterozygosity.length() > 0) {
        fixed_heterozygosity = true;
        cout << setprecision(5);
        if (pomo_heterozygosity == "EMP") {
            cout << "Level of polymorphism is fixed to the estimate from the data: ";
            cout << heterozygosity << "." << endl;
            fixed_heterozygosity_emp = true;
            // No need to set the level of polymorphism here because
            // of initialization.
        }
        else {
            cout << "Heterozygosity is fixed to the value given by the user: ";
            heterozygosity = convert_double(pomo_heterozygosity.c_str());
            cout << heterozygosity << "." << endl;
            fixed_heterozygosity_usr = true;
        }
    }
}


void ModelPoMo::init(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     string pomo_heterozygosity) {
    // Initialize model constants.
    N = phylo_tree->aln->virtual_pop_size;
    n_alleles = 4;
    n_connections = n_alleles * (n_alleles-1) / 2;
    eps = 1e-8;
    // Check if number of states of PoMo match the provided data.
    ASSERT(num_states == (n_alleles + (n_alleles*(n_alleles-1)/2 * (N-1))) );

    // Main initialization of model and parameters.
    init_mutation_model(model_name,
                        model_params,
                        freq_type,
                        freq_params);
    init_sampling_method();
    init_boundary_frequencies();
    // Initialize heterozygosity and the scale factor of the mutation rates
    // which is only used for Gamma rate heterogeneity at the moment.
    heterozygosity = estimateEmpiricalWattersonTheta();
    scale = 1.0;
    init_fixed_parameters(model_params, pomo_heterozygosity);
    set_heterozygosity_boundaries();
    setInitialMutCoeff();
    //rate_matrix = new double[num_states*num_states];
    ignore_state_freq = true;
    normalize_matrix = false;
    half_matrix = false;
    delete [] rates;
    rates = new double[num_states*num_states];
    updatePoMoStatesAndRateMatrix();
    decomposeRateMatrix();


    this->name = mutation_model->getName();
//    if (model_params.length() > 0)
//        this->name += "{" + model_params + "}";
    this->name += "+P";
//    if (pomo_heterozygosity.length() > 0)
//        this->name += "{" + pomo_heterozygosity + "}";
    this->name += "+N" + convertIntToString(N);

    if (verbose_mode >= VB_MED) {
        cout << "Initialized PoMo model." << endl;
        cout << "Model name: " << this->name << "." << endl;
        cout << this->full_name << endl;
    }
    if (verbose_mode >= VB_MAX)
        writeInfo(cout);
}

ModelPoMo::~ModelPoMo() {
    // Rate matrix is deleted by ~ModelMarkov().
    // delete [] rate_matrix;
    delete mutation_model;
    delete [] freq_boundary_states_emp;
    delete [] mutation_rate_matrix;
}

double ModelPoMo::computeSumFreqBoundaryStates() {
    int i;
    double norm_boundary = 0.0;
    for (i = 0; i < n_alleles; i++)
        norm_boundary += freq_boundary_states[i];
    // Should be 1.0!
    if ((norm_boundary > 1.0 + eps) || (norm_boundary < 1.0 - eps))
      outError("Calculation of boundary state frequencies faulty (maybe a numerical problem).");
    return norm_boundary;
}

double harmonic(int n) {
    double harmonic = 0.0;
    for (int i = 1; i <= n; i++)
        harmonic += 1.0/(double)i;
    return harmonic;
}

void ModelPoMo::setInitialMutCoeff() {
    mutation_rate_matrix = new double[n_alleles*n_alleles];
    memset(mutation_rate_matrix, 0, n_alleles*n_alleles*sizeof(double));

    // Check if polymorphism data is available.
    double lambda_poly_sum_no_mu = computeSumFreqPolyStatesNoMut();
    if (lambda_poly_sum_no_mu <= 0) {
      outWarning("We discourage usage of PoMo on data without polymorphisms.");
      if (!fixed_heterozygosity_usr)
        outError("Please fix the heterozygosity when population data is unavailable.");
    }

    normalizeMutationRates();
}

double ModelPoMo::computeSumFreqPolyStatesNoMut() {
    double norm_polymorphic = 0.0;
    int i, j;
    for (i = 0; i < n_alleles; i++) {
        for (j = 0; j < i; j++)
            norm_polymorphic +=
                2 * freq_boundary_states[i] * freq_boundary_states[j];
    }
    norm_polymorphic *= harmonic(N-1);
    return norm_polymorphic;
}

double ModelPoMo::computeSumFreqPolyStates() {
    double norm_polymorphic = 0.0;
    double dpoly;
    int i, j;
    for (i = 0; i < n_alleles; i++) {
      for (j = 0; j < n_alleles; j++) {
        if (i != j) {
          dpoly = freq_boundary_states[i] * mutation_rate_matrix[i*n_alleles+j];
          norm_polymorphic += dpoly;
        }
      }
    }
    norm_polymorphic *= harmonic(N-1);
    return norm_polymorphic;
}

double ModelPoMo::computeNormConst() {
    double norm_boundary = computeSumFreqBoundaryStates();
    double norm_polymorphic = computeSumFreqPolyStates();
    return 1.0/(norm_boundary + norm_polymorphic);
}

void ModelPoMo::computeStateFreq () {
    double norm = computeNormConst();
    int state;

    double * m = mutation_rate_matrix;
    // double * r = mutation_rate_matrix_sym;
    // double * f = mutation_rate_matrix_asy;
    double * pi = freq_boundary_states;
    int n = n_alleles;

    for (state = 0; state < num_states; state++) {
        if (isBoundary(state))
            state_freq[state] = pi[state]*norm;
        else {
            // Allele count, first and second type.
            int i, a, b;
            decomposeState(state, i, a, b);
            // state_freq[state] = norm * pi[a] * pi[b];
            // state_freq[state] *= (m[a*n + b]*1.0/(N-i) + m[b*n + a]*1.0/i);
            state_freq[state] = pi[a]*m[a*n+b]/(N-i) + pi[b]*m[b*n+a]/i;
            state_freq[state] *= norm;
            // double sym =  r[a*n + b]*(1.0/i + 1.0/(N-i));
            // double asy = -f[a*n + b]*(1.0/i - 1.0/(N-i));
            // state_freq[state] *= (sym + asy);
        }
    }
}

void ModelPoMo::updatePoMoStatesAndRateMatrix () {
    computeStateFreq();

    // Compute rate matrix.
    int i, j;
    double tot_sum = 0.0;
    for (i = 0; i < num_states; i++) {
        double row_sum = 0.0;
        // Loop over columns in row state1 (transition to state2).
        for (j = 0; j < num_states; j++)
            if (i != j) {
                row_sum +=
                    (rates[i*num_states+j] =
                     computeProbBoundaryMutation(i, j));
            }
        tot_sum += state_freq[i]*row_sum;
        // diagonal will be handled later, should not assign now
        //rates[i*num_states+i] = -(row_sum);
        rates[i*num_states+i] = 0.0;
    }
    // Thu Aug 17 16:11:19 BST 2017; Dom. Normalization is preferred. Then,
    // branch lengths can be interpreted in an easy way (the length equals the
    // estimated number of events). Without normalization, branch lengths are
    // harder to interpret although from a biological point of view they make
    // more sense (see also discussion below).

    // Normalize rate matrix such that one event happens per unit time.
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            rates[i*num_states+j] /= tot_sum;
        }
    }

    // Sun Jul 16 17:43:30 BST 2017; Dom. I removed normalization, it should not
    // change output nor stability.

    // Tue Jul 18 12:32:00 BST 2017; Dom. However, it changes the branch
    // lengths, which should be normalized to one event per unit time.

    // Fri Jul 21 16:56:40 BST 2017; Dom. Going back to having no normalization,
    // because I am comparing trees to the boundary mutation model simulator and
    // they should have the same lengths.
}

void ModelPoMo::decomposeState(int state, int &i, int &nt1, int &nt2) {
    if (state < 4) {
        // Boundary A, C, G or T
        i = N;
        nt1 = state;
        nt2 = -1; // -1 for unknown nt
    } else if (state < 4+(N-1)) {
        // (iA,N-iC)
        i = state-3;
        nt1 = 0; // A
        nt2 = 1; // C
    } else if (state < 4+2*(N-1)) {
        // (iA,N-iG)
        i = state-3-(N-1);
        nt1 = 0; // A
        nt2 = 2; // G
    } else if (state < 4+3*(N-1)) {
        // (iA,N-iT)
        i = state-3-2*(N-1);
        nt1 = 0; // A
        nt2 = 3; // T
    } else if (state < 4+4*(N-1)) {
        // (iC,N-iG)
        i = state-3-3*(N-1);
        nt1 = 1; // C
        nt2 = 2; // G
    } else if (state < 4+5*(N-1)) {
        // (iC,N-iT)
        i = state-3-4*(N-1);
        nt1 = 1; // C
        nt2 = 3; // T
    } else if (state < 4+6*(N-1)) {
        // (iG,N-iT)
        i = state-3-5*(N-1);
        nt1 = 2; // G
        nt2 = 3; // T
    } else {
        outError("State exceeds limit");
    }
}

bool ModelPoMo::isBoundary(int state) {
    return (state < 4);
}

bool ModelPoMo::isPolymorphic(int state) {
    return (!isBoundary(state));
}

double ModelPoMo::computeProbBoundaryMutation(int state1, int state2) {
    // The transition rate to the same state will be calculated by
    // setting the row sum to 0.
    ASSERT(state1 != state2);

    // Both states are decomposed into the abundance of the first
    // allele as well as the nucleotide of the first and the second
    // allele.
    int i1=0, i2=0, nt1=-1, nt2=-1, nt3=-1, nt4=-1;
    decomposeState(state1, i1, nt1, nt2);
    decomposeState(state2, i2, nt3, nt4);

    // Either the first nucleotides match or the first of state 1 with
    // the second of state 2 or the first of state 2 with the second
    // of state 1.  Additionally, we have to consider boundary states as
    // special cases.
    if (nt1 == nt3 && (nt2==nt4 || nt2==-1 || nt4 == -1)) {
        ASSERT(i1 != i2); // because state1 != state2
        if (i1+1==i2)
            // e.g.: 2A8C -> 3A7C or 9A1C -> 10A
            // Changed Dom Tue Sep 29 13:31:02 CEST 2015
            // return double(i1*(N-i1)) / double(N*N);
            return double(i1*(N-i1)) / double(N);
        else if (i1-1 == i2)
            // e.g.: 3A7C -> 2A8C or 10A -> 9A1C
            if (nt2 == -1)
                // e.g. 10A -> 9A1C
                // return mutCoeff(nt1,nt4) * state_freq[nt4];
                return mutation_rate_matrix[nt1*n_alleles+nt4];
            else
                // e.g. 9A1C -> 8A2C
                // Changed Dom Tue Sep 29 13:30:43 CEST 2015
                // return double(i1*(N-i1)) / double(N*N);
                return double(i1*(N-i1)) / double(N);
        else
            return 0.0;
    } else if (nt1 == nt4 && nt2 == -1 && i2 == 1)  {
        // e.g.: 10G -> 1A9G
        //return mutCoeff(nt1,nt3) * state_freq[nt3];
        return mutation_rate_matrix[nt1*n_alleles+nt3];
    } else if (nt2 == nt3  && i1 == 1 && nt4 == -1) {
        // E.g.: 1A9G -> 10G
        // Changed Dom Tue Sep 29 13:30:25 CEST 2015
        // return double(i1*(N-i1)) / double(N*N);
        return double(i1*(N-i1)) / double(N);
    } else
        // 0 for all other transitions
        return 0.0;
}

int ModelPoMo::getNDim() {
    if (fixed_heterozygosity)
        return mutation_model->getNDim();
    else
        return mutation_model->getNDim()+1;
}

int ModelPoMo::getNDimFreq() {
    return mutation_model->getNDimFreq();
}

void ModelPoMo::setBounds(double *lower_bound,
                          double *upper_bound,
                          bool *bound_check) {
    // Set boundaries of underlying mutation model.
    mutation_model->setBounds(lower_bound, upper_bound, bound_check);

    // Level of polymorphism.
    if (!fixed_heterozygosity) {
        int ndim = getNDim();
        lower_bound[ndim] = min_heterozygosity;
        upper_bound[ndim] = max_heterozygosity;
        bound_check[ndim] = false;
    }
}

// Get rates from underlying mutation model and normalize them such
// that the heterozygosity matches.
void ModelPoMo::normalizeMutationRates() {
    double * m = mutation_rate_matrix;
    mutation_model->getQMatrix(m);

    // Normalize the mutation probability so that they resemble the given level
    // of polymorphism. The way of normalization determines heterozygosity.
    computeStateFreq();
    double poly = computeSumFreqPolyStates();
    double theta_bm = poly/harmonic(N-1);

    // The correction factor is exactly the difference between sampling with and
    // without replacement. Without replacement, the correction factor is 1.0
    // and the heterozygosity values are very far off if N is low. Even with the
    // correction, the estimated heterozygosity is still too high when N is low.
    // A correct calculation of the heterozygosity requires a rethinking of the
    // sampling step together with a correction of heterozygosity for low N.

    // No correction, sampling without replacement:
    double correction = 1.0;

    // Correction, sampling with replacement (see above), seems to
    // work better to estimate level of polymorphism but deteriorates
    // branch scrore distance, so deactivated.
    // double correction = (double) (N-1) / double (N);

    // Interestingly, a correction factor of (N-1)/(N+1) gives very
    // good results but I cannot derive it and do not know why.
    // double correction = (double) (N-1) / double (N+1);

    // Mon Apr 17 10:21:09 BST 2017. See Eq. (12.14) in Dominik's thesis.

    // heterozygosity: desired heterozygosity or level of polymorphism (Watterson's
    // estimator is mostly used to estimate this quantity).

    // theta_bm: the scaled mutation rate parameter of the PoMo model.

    // scale: a parameter that is used for the gamma rate heterogeneity model.

    double m_norm = scale * (heterozygosity / (theta_bm * (correction - harmonic(N-1) * heterozygosity)));

    if (verbose_mode >= VB_MAX) {
      cout << "Normalization constant of mutation rates: " << m_norm << endl;
      // cout << "Heterozygosity: " << heterozygosity << endl;
    }

    for (int i = 0; i < n_alleles; i++)
      for (int j = 0; j < n_alleles; j++) {
        m[i*n_alleles+j] *= m_norm;
        // DEBUG.
        // cout << setprecision(15);
        // cout << m[i*n_alleles+j] << endl;
      }

    // DEBUG.
    if (verbose_mode >= VB_MED) {
      cout << "theta_bm before normalization is " << theta_bm << endl;
      cout << "heterozygosity is " << heterozygosity << endl;
      for (int i = 0; i < n_alleles; i++) {
        for (int j = i+1; j < n_alleles; j++) {
          // The mutation rate matrix entry is the exchangeability times the
          // target allele frequency.
          double ex = m[i*n_alleles+j] / freq_boundary_states[j];
          cout << setprecision(15);
          cout << "Exchangeability " << i << " to " << j << " is " << ex << endl;
        }
      }
      computeStateFreq();
      double normc = computeNormConst();
      theta_bm = (1.0 - normc) / harmonic(N-1);
      cout << "theta_bm after normalization is " << theta_bm << endl;
    }
}

void ModelPoMo::setScale(double new_scale) {
  scale = new_scale;
  // DEBUG.
  // if (heterozygosity*new_scale < min_heterozygosity || heterozygosity*new_scale > max_heterozygosity) {
  //   cout << "Scale: " << scale << endl;
  //   outWarning("After rescaling, heterozygosity out of range. Numerical instabilities expected.");
  // }

  normalizeMutationRates();
}

double ModelPoMo::getScale() {
  return scale;
}

bool ModelPoMo::getVariables(double *variables) {
    bool changed = false;
    changed = mutation_model->getVariables(variables);

    if (!fixed_heterozygosity) {
        int ndim = getNDim();
        changed |= (heterozygosity != variables[ndim]);
        heterozygosity = variables[ndim];
    }

    normalizeMutationRates();
    updatePoMoStatesAndRateMatrix();
    return changed;
}


// Dummy function, see declaration in header file.
void ModelPoMo::setRates() {
    return;
}

void ModelPoMo::setVariables(double *variables) {
    mutation_model->setVariables(variables);

    if (!fixed_heterozygosity) {
        int ndim = getNDim();
        variables[ndim] = heterozygosity;
    }
}

void ModelPoMo::writeInfo(ostream &out) {
  report(out);
}

// TODO: s_freqs is not used.
void ModelPoMo::computeRateMatrix(double **r_matrix, double *s_freqs, int n_states) {
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            r_matrix[i][j] = rates[i*n_states+j];
        }
    }
}

double ModelPoMo::targetFunk(double x[]) {
    getVariables(x);
    // Disable test for low stationary frequency in PoMo.
    // if (state_freq[num_states-1] < 1e-4) return 1.0e+12;
    decomposeRateMatrix();
    ASSERT(phylo_tree);
    phylo_tree->clearAllPartialLH();
    return -phylo_tree->computeLikelihood();
}

bool ModelPoMo::isUnstableParameters() {
    // More checking could be done.
    for (int i = 0; i < num_states; i++) {
        if (state_freq[i] < eps) return true;
    }
    return false;
}

void ModelPoMo::normalize_boundary_freqs(double * bfs) {
    // Normalize frequencies so that they sum to 1.0.
    double sum = 0.0;
    for (int i = 0; i < n_alleles; i++)
        sum += bfs[i];
    for (int i = 0; i < n_alleles; i++)
        bfs[i] /= sum;
    if (verbose_mode >= VB_MAX) {
        std::cout << "The empirical frequencies of the boundary states are:" << std::endl;
        for (int i = 0; i < n_alleles; i++)
            std::cout << bfs[i] << " ";
        std::cout << std::endl;
    }
    check_boundary_freqs(bfs);
}

void ModelPoMo::check_boundary_freqs(double * bfs) {
    // Check if boundary frequencies are within bounds.
    bool change = false;
    for (int i = 0; i < n_alleles; i++) {
        if (bfs[i] < POMO_MIN_BOUNDARY_FREQ) {
            bfs[i] = POMO_MIN_BOUNDARY_FREQ;
            outWarning("A boundary state has very low frequency.");
            cout << "Frequency set to." << POMO_MIN_BOUNDARY_FREQ;
            change = true;
        }
        if (bfs[i] > POMO_MAX_BOUNDARY_FREQ) {
            bfs[i] = POMO_MAX_BOUNDARY_FREQ;
            outWarning("A boundary state has very high frequency.");
            cout << "Frequency set to." << POMO_MAX_BOUNDARY_FREQ;
            change = true;
        }
    }
    if (change)
        normalize_boundary_freqs(bfs);
}

void
ModelPoMo::estimateEmpiricalBoundaryStateFreqs(double * freq_boundary_states)
{
    memset(freq_boundary_states, 0, sizeof(double)*n_alleles);

    if (sampling_method == SAMPLING_SAMPLED) {
        unsigned int abs_state_freq[num_states];
        memset(abs_state_freq, 0, sizeof(unsigned int)*num_states);
        phylo_tree->aln->computeAbsoluteStateFreq(abs_state_freq);
        int n;
        int x;
        int y;

        int sum[n_alleles];
        int tot_sum = 0;
        memset (sum, 0, n_alleles * sizeof(int));

        for (int i = 0; i < num_states; i++) {
            decomposeState(i, n, x, y);
            sum[x]+= n*abs_state_freq[i];
            if (y >= 0) sum[y]+= (N-n)*abs_state_freq[i];
        }
        for (int i = 0; i < n_alleles; i++) {
            tot_sum += sum[i];
        }
        for (int i = 0; i < n_alleles; i++) {
            freq_boundary_states[i] = (double) sum[i]/tot_sum;
        }
        // Output vector if verbose mode.
        if (verbose_mode >= VB_MAX) {
            std::cout << "Absolute empirical state frequencies:" << std::endl;
            for (int i = 0; i < num_states; i++)
                std::cout << abs_state_freq[i] << " ";
            std::cout << std::endl;
        }
        // Set highest_freq_state.
        for (int i = 0; i < num_states; i++)
            if (abs_state_freq[i] > abs_state_freq[highest_freq_state])
                highest_freq_state = i;
    } else {
        for (Alignment::iterator it = phylo_tree->aln->begin();
             it != phylo_tree->aln->end(); it++) {
            for (Pattern::iterator it2 = it->begin(); it2 != it->end(); it2++) {
                int state = (int)*it2;
                if (state < num_states)
                    outError("Unknown PoMo state in pattern.");
                else if ((unsigned int)state == phylo_tree->aln->STATE_UNKNOWN)
                    continue;
                state -= num_states;
                ASSERT((unsigned int)state < phylo_tree->aln->pomo_sampled_states.size());
                // Decode the id and counts.
                int id1 = phylo_tree->aln->pomo_sampled_states[state] & 3;
                int id2 = (phylo_tree->aln->pomo_sampled_states[state] >> 16) & 3;
                int j1 = (phylo_tree->aln->pomo_sampled_states[state] >> 2) & 16383;
                int j2 = (phylo_tree->aln->pomo_sampled_states[state] >> 18);
                freq_boundary_states[id1] += j1*(it->frequency);
                freq_boundary_states[id2] += j2*(it->frequency);
            }
        }
    }
    normalize_boundary_freqs(freq_boundary_states);
    if (verbose_mode >= VB_MAX) {
        cout << "Empirical boundary state frequencies: ";
        for (int i = 0; i< n_alleles; i++)
            cout << freq_boundary_states[i] << " ";
        cout << endl;
    }
}

double
ModelPoMo::estimateEmpiricalWattersonTheta()
{
    double theta_p = 0.0;
    int sum_pol = 0;
    int sum_fix = 0;
    double sum_theta_w = 0.0;

    if (sampling_method == SAMPLING_SAMPLED) {
        unsigned int abs_state_freq[num_states];
        memset(abs_state_freq, 0, sizeof(unsigned int)*num_states);
        phylo_tree->aln->computeAbsoluteStateFreq(abs_state_freq);
        for (int i = 0; i < n_alleles; i++) sum_fix += abs_state_freq[i];
        for (int i = n_alleles; i < num_states; i++) sum_pol += abs_state_freq[i];
        theta_p = (double) sum_pol / (double) (sum_fix + sum_pol);
        // TODO DS: This is wrong because Watterson's estimator is
        // expected to decrease when sampling step is performed.  Some
        // sequences will be taken more often and necessarily,
        // polymorphism is lost.
    } else {
        for (Alignment::iterator it = phylo_tree->aln->begin();
             it != phylo_tree->aln->end(); it++) {
            for (Pattern::iterator it2 = it->begin(); it2 != it->end(); it2++) {
                int state = (int)*it2;
                if (state < num_states)
                    outError("Unknown PoMo state in pattern.");
                else if ((unsigned int)state == phylo_tree->aln->STATE_UNKNOWN)
                    continue;
                state -= num_states;
                ASSERT((unsigned int)state < phylo_tree->aln->pomo_sampled_states.size());
                // Decode counts.
                int j1 = (phylo_tree->aln->pomo_sampled_states[state] >> 2) & 16383;
                int j2 = (phylo_tree->aln->pomo_sampled_states[state] >> 18);
                if (j2 == 0) sum_fix += it->frequency;
                else {
                    // Have to use Watterson Theta because sample size may be different.
                    sum_pol += it->frequency;
                    sum_theta_w += (double) it->frequency / harmonic(j1 + j2 - 1);
                }
            }
        }
        // Calculate Watterson Theta per site.
        double theta_w_temp = sum_theta_w;
        sum_theta_w = theta_w_temp / (double) (sum_fix + sum_pol);
        theta_p = sum_theta_w;
    }
    // Output vector if verbose mode.
    if (verbose_mode >= VB_MAX) {
        cout << setprecision(8);
        cout << "Estimated relative frequency of polymorphic states:" << std::endl;
        cout << theta_p << std::endl;
        cout << setprecision(5);
    }
    // Normalize frequencies so that the last entry is 1.0.
    return theta_p;
}

void ModelPoMo::report_model_params(ostream &out, bool reset_scale) {
  if (fixed_model_params)
    out << "User-defined model parameters." << endl;
  else
    out << "Model parameters." << endl;
  // By default, reset scale before reporting.
  if (reset_scale)
    setScale(1.0);
  else
    out << "The reported rates are scaled by a factor of " << scale << "." << endl;;

  // TODO: If verbose, output rate matrix.
  if (verbose_mode >= VB_MED) {
    out << "Rate matrix: " << endl;
    for (int i = 0; i < num_states; i++) {
      for (int j = 0; j < num_states; j++) {
        out << setprecision(8);
        out << setw(8) << rates[i*num_states+j] << " ";
      }
      out << endl;
    }
  }

  // Report rates.
  // Mutation rates.
  double *rs = NULL;
  // Exchangeabilities.
  double *es = NULL;
  if (is_reversible) {
    rs = new double[n_connections];
    es = new double[n_connections];
  }
  else if (!is_reversible) {
    rs = new double[2*n_connections];
    es = new double[2*n_connections];
  }
  else outError("Model has to be either reversible or non-reversible.");
  rate_matrix_to_rates(mutation_rate_matrix, rs);
  rate_matrix_to_exchangeabilities(mutation_rate_matrix, es);
  // mutation_model->writeInfo(out);
  mutation_model->report_rates(out, "Mutation    rates", rs);
  mutation_model->report_rates(out, "Exchangeabilities", es);
  delete [] rs;
  // Report stationary frequencies.
  mutation_model->report_state_freqs(out);

  int n = n_alleles;
  if (!is_reversible) {
    out << setprecision(5);
    out << "Mutation rate matrix: " << endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        out << setw(8) << mutation_rate_matrix[i*n+j] << " ";
      out << endl;
    }

    // Exchangeability matrix.
    double * e = new double[n_alleles*n_alleles];
    memcpy(e, mutation_rate_matrix, n_alleles*n_alleles * sizeof (double));
    double * pi = freq_boundary_states;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        e[i*n+j] /= pi[j];
      }
    }
    // Reversible part of exchangeability matrix (Q^GTR or Q^REV / PI).
    double * r = new double[n_alleles*n_alleles];
    // Non-reversible part of mutation rate matrix (Q^FLUX / PI).
    double * f = new double[n_alleles*n_alleles];
    memset(r, 0, n_alleles*n_alleles*sizeof(double));
    memset(f, 0, n_alleles*n_alleles*sizeof(double));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        r[i*n+j] = e[i*n+j] + e[j*n+i];
        r[i*n+j] /= 2.0;
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i==j)
          f[i*n+j] = 0;
        else {
          f[i*n+j] = e[i*n+j] - e[j*n+i];
          f[i*n+j] /= 2.0;
        }
      }
    }

    out << "Reversible part of the exchangeability matrix:" << endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        out << setw(8) << r[i*n+j] << " ";
      out << endl;
    }
    out << "Non-reversible part of the exchangeability rate matrix:" << endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        out << setw(8) << f[i*n+j] << " ";
      out << endl;
    }

    out << setprecision(4);
    out << "Comparison of Frobenius norms." << endl;
    double frob_norm_e = frob_norm(e, n);
    double frob_norm_f = frob_norm(f, n);
    out << "Exchangeability matrix: " << setw(6) << frob_norm_e << endl;
    out << "Non-reversible part:    " << setw(6) << frob_norm_f << endl;
    out << "Ratio:                  " << setw(6) << frob_norm_f / frob_norm_e << endl;

    delete [] r;
    delete [] f;
    delete [] e;
  }

  // Report heterozygosity.
  out << setprecision(4);
  if (!fixed_heterozygosity) {
    out << "Estimated heterozygosity: " << heterozygosity << endl;
    if (sampling_method == SAMPLING_WEIGHTED_BINOM)
      out << "A slight overestimation is expected and an effect of weighted binomial sampling; see manual." << endl;
  }
  else if (fixed_heterozygosity_emp)
    out << "Empirical heterozygosity: " << heterozygosity << endl;
  else if (fixed_heterozygosity_usr)
    out << "User-defined heterozygosity: " << heterozygosity << endl;
  else
    outError("It is undefined how the heterozygosity was determined.");

  return;
}

void ModelPoMo::rate_matrix_to_exchangeabilities(double *m, double *r) {
  int c = 0;
  for (int i = 0; i < n_alleles; i++) {
    for (int j = i+1; j < n_alleles; j++) {
      r[c] = m[i*n_alleles+j] / freq_boundary_states[j];
      c++;
    }
  }
  if (!is_reversible) {
    for (int i = 0; i < n_alleles; i++) {
      for (int j = 0; j < i; j++) {
        r[c] = m[i*n_alleles+j] / freq_boundary_states[j];
        c++;
      }
    }
  }
}

void ModelPoMo::rate_matrix_to_rates(double *m, double *r) {
  int c = 0;
  for (int i = 0; i < n_alleles; i++) {
    for (int j = i+1; j < n_alleles; j++) {
      r[c] = m[i*n_alleles+j];
      c++;
    }
  }
  if (!is_reversible) {
    for (int i = 0; i < n_alleles; i++) {
      for (int j = 0; j < i; j++) {
        r[c] = m[i*n_alleles+j];
        c++;
      }
    }
  }
}

void ModelPoMo::report(ostream &out) {
  ios  state(NULL);
  out << this->full_name << endl;
  out << endl;

  // Model parameters.
  out << "--" << endl;
  report_model_params(out);

  // Empirical quantities.
  out << "--" << endl;
  out << "Empirical quantities." << endl;
  mutation_model->report_state_freqs(out, freq_boundary_states_emp);
  out << setprecision(4);
  out << "Watterson's estimator of heterozygosity: " << estimateEmpiricalWattersonTheta() << endl;
  out << "--" << endl << endl;
}

void ModelPoMo::setCheckpoint(Checkpoint *checkpoint) {
	ModelMarkov::setCheckpoint(checkpoint);
    mutation_model->setCheckpoint(checkpoint);
}

void ModelPoMo::startCheckpoint() {
    checkpoint->startStruct("ModelPoMo");
}

void ModelPoMo::saveCheckpoint() {
//    int n_rates = n_alleles * (n_alleles-1) / 2;
    startCheckpoint();
//    CKP_ARRAY_SAVE(n_rates, mutation_model->rates);
//    CKP_ARRAY_SAVE(n_alleles, mutation_model->state_freq);
    mutation_model->saveCheckpoint();
    CKP_SAVE(heterozygosity);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

// TODO DS: Check checkpointing :-).
void ModelPoMo::restoreCheckpoint() {
//    int n_rates = n_alleles * (n_alleles-1) / 2;
    // First, get variables from checkpoint.
    startCheckpoint();
//    CKP_ARRAY_RESTORE(n_rates, mutation_model->rates);
//    CKP_ARRAY_RESTORE(n_alleles, mutation_model->state_freq);
    mutation_model->restoreCheckpoint();
    CKP_RESTORE(heterozygosity);
    endCheckpoint();
    // Second, restore underlying mutation model.
    ModelMarkov::restoreCheckpoint();
    normalizeMutationRates();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

// Declaration of helper function; needed by decomposeRateMatrix().
//int computeStateFreqFromQMatrix (double Q[], double pi[], int n, double space[]);

// void ModelPoMo::decomposeRateMatrix() {
//     updatePoMoStatesAndRateMatrix();
//     // Non-reversible.
//     if (!is_reversible) {
//         if (phylo_tree->params->matrix_exp_technique == MET_EIGEN_DECOMPOSITION) {
//             eigensystem_nonrev(rate_matrix, state_freq, eigenvalues, eigenvalues_imag, eigenvectors, inv_eigenvectors, num_states);
//             return;
//         }
//         else if (phylo_tree->params->matrix_exp_technique == MET_SCALING_SQUARING) {
//             return;
//         }
//         else if (phylo_tree->params->matrix_exp_technique == MET_EIGEN3LIB_DECOMPOSITION) {
//             // Not (yet?) implemented.
//             // decomposeRateMatrixEigen3lib();
//             outError("MET_EIGEN3LIB_DECOMPOSITION does not work with PoMo.");
//         }
//         else if (phylo_tree->params->matrix_exp_technique == MET_LIE_MARKOV_DECOMPOSITION)
//             // Not possible?
//             // decomposeRateMatrixClosedForm();
//             outError("Matrix decomposition in closed form not available for PoMo.");
//         else
//             outError("Matrix decomposition method unknown.");
//     }
//     // Reversible.  Alogrithms for symmetric matrizes can be used.
//     else {
//         // TODO DS: This leaves room for speed improvements.
//         // EigenDecomposition::eigensystem_sym() expects a matrix[][]
//         // object with two indices.  However, it is not used, because
//         // ModelPoMo::computeRateMatrix() is called anyways from
//         // within eigensystem_sym().
// 		double **temp_matrix = new double*[num_states];
// 		for (int i = 0; i < num_states; i++)
// 			temp_matrix[i] = new double[num_states];
// 		eigensystem_sym(temp_matrix, state_freq, eigenvalues, eigenvectors, inv_eigenvectors, num_states);
// 		for (int i = num_states-1; i >= 0; i--)
// 			delete [] temp_matrix[i];
// 		delete [] temp_matrix;
//     return;
//     }
// }

void ModelPoMo::decomposeRateMatrix() {
    updatePoMoStatesAndRateMatrix();
    ModelMarkov::decomposeRateMatrix();
    return;
}

void ModelPoMo::set_heterozygosity_boundaries() {
  min_heterozygosity = POMO_MIN_REL_HETEROZYGOSITY * heterozygosity;
  max_heterozygosity = POMO_MAX_REL_HETEROZYGOSITY * heterozygosity;
  if (min_heterozygosity < POMO_MIN_HETEROZYGOSITY)
    outWarning("The polymorphism level in the data is very low.");
  if (max_heterozygosity > POMO_MAX_HETEROZYGOSITY)
    outWarning("The polymorphism level in the data is very high.");
}

// TODO DS or Minh: This is a copy of ModelLieMarkov::computeTransMatrix(). I do
// not see why every model needs its own function, especially when it is so
// generic like this one. However, I thought my problems may be related to the
// implementation in ModelMarkov::computeTransMatrix(). Now I am not so sure
// anymore.

// TODO DS: The parameter mixture is unused at the moment.
void ModelPoMo::computeTransMatrix(double time, double *trans_matrix, int mixture, int selected_row) {
  MatrixExpTechnique technique = phylo_tree->params->matrix_exp_technique;
  if (technique == MET_SCALING_SQUARING || !is_reversible) {
    // Do not change the object rate_matrix, but only trans_matrix.
    Eigen::Map<Eigen::MatrixXd> A(rate_matrix, num_states, num_states);
    Eigen::Map<Eigen::MatrixXd> P(trans_matrix, num_states, num_states);
    P = ((A.transpose() * time).exp()).transpose();
    int i, j;
    for (i = 0; i < num_states; i++) {
      double sum = 0.0;
      for (j = 0; j < num_states; j++)
        sum += (trans_matrix[i*num_states+j]);
      ASSERT(fabs(sum-1.0) < 1e-4);
    }
  }
  else if (technique == MET_EIGEN3LIB_DECOMPOSITION) {
    outError("TODO DS: EIGEN3LIB DECOMPOSITION not yet tested.");
    // and nondiagonalizable == false, else we used scaled squaring
    int i;
    Eigen::VectorXcd ceval_exp;
    ceval_exp.resize(num_states);
    for (i = 0; i < 4; i++)
      ceval_exp(i) = exp(ceval[i]*time);
    Eigen::Map<Eigen::MatrixXcd> cevectors(cevec, num_states, num_states);
    Eigen::Map<Eigen::MatrixXcd> cinv_evectors(cevec, num_states, num_states);
    Eigen::MatrixXcd res;
    res.resize(num_states, num_states);
    res = cevectors * ceval_exp.asDiagonal() * cinv_evectors;
    // if assertions fail, it may be due to cevec having near-zero
    // determinant, and a fix could be to relax the test for
    // nondiagonalizable in ModelLieMarkov::decomposeRateMatrixEigen3lib()
    if (technique == MET_EIGEN3LIB_DECOMPOSITION) {
			for (i = 0; i < num_states; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < num_states; j++) {
					trans_matrix[i*num_states+j] = res(j, i).real();
					ASSERT(fabs(res(j,i).imag()) < 1e-6);
					ASSERT(trans_matrix[i*num_states+j] >= -0.000001);
					ASSERT(trans_matrix[i*num_states+j] <=  1.000001);
					if (trans_matrix[i*num_states+j] < 0)
						trans_matrix[i*num_states+j] = 0.0;
					if (trans_matrix[i*num_states+j] > 1)
						trans_matrix[i*num_states+j] = 1.0;
          row_sum += trans_matrix[i*num_states+j];
				}
				ASSERT(fabs(row_sum-1.0) < 1e-4);
			}
    }
  }

  else ModelMarkov::computeTransMatrix(time, trans_matrix, 0, selected_row);
}

void ModelPoMo::computeTipLikelihood(PML::StateType state, double *lh) {
    Alignment *aln = phylo_tree->aln;
    if (state < num_states || state >= num_states+aln->pomo_sampled_states.size()) {
        ModelSubst::computeTipLikelihood(state, lh);
        return;
    }
    state = state - num_states;

  bool hypergeometric = (aln->pomo_sampling_method == SAMPLING_WEIGHTED_HYPER);
  int nstates = aln->num_states;
  int N = aln->virtual_pop_size;

  memset(lh, 0, sizeof(double)*nstates);

  // decode the id and value
  int id1 = aln->pomo_sampled_states[state] & 3;
  int id2 = (aln->pomo_sampled_states[state] >> 16) & 3;
  int j = (aln->pomo_sampled_states[state] >> 2) & 16383;
  int M = j + (aln->pomo_sampled_states[state] >> 18);

  // Number of alleles is hard coded here, change if generalization is needed.
  int nnuc = 4;

  // TODO DS: Implement down sampling or a better approach.
  if (hypergeometric && M > N)
    outError("Down sampling not yet supported.");

  // Check if observed state is a fixed one.  If so, many
  // PoMo states can lead to this data.  E.g., even (2A,8T)
  // can lead to a sampled data of 7A.
  if (j == M) {
    lh[id1] = 1.0;
    // Second: Polymorphic states.
    for (int s_id1 = 0; s_id1 < nnuc-1; s_id1++) {
      for (int s_id2 = s_id1+1; s_id2 < nnuc; s_id2++) {
        if (s_id1 == id1) {
          // States are in the order {FIXED,
          // 1A(N-1)C, ..., (N-1)A1C, ...}.
          int k;
          if (s_id1 == 0) k = s_id2 - 1;
          else k = s_id1 + s_id2;
          // Start one earlier because increment
          // happens after execution of for loop
          // body.
          int real_state = nnuc - 1 + k*(N-1) + 1;
          for (int i = 1; i < N; i++, real_state++) {
            ASSERT(real_state < nstates);
            if (!hypergeometric)
              lh[real_state] = std::pow((double)i/(double)N,j);
            else {
              lh[real_state] = 1.0;
              for (int l = 0; l<j; l++)
                lh[real_state] *= (double) (i-l) / (double) (N-l);
            }
          }
        }
        // Same but fixed allele is the second one
        // in polymorphic states.
        else if (s_id2 == id1) {
          int k;
          if (s_id1 == 0) k = s_id2 - 1;
          else k = s_id1 + s_id2;
          int real_state = nnuc - 1 + k*(N-1) + 1;
          for (int i = 1; i < N; i++, real_state++) {
            ASSERT(real_state < nstates);
            if (!hypergeometric)
              lh[real_state] = std::pow((double)(N-i)/(double)N,j);
            else {
              lh[real_state] = 1.0;
              for (int l = 0; l<j; l++)
                lh[real_state] *= (double) (N-i-l) / (double) (N-l);
            }
          }
        }
      }
    }
  }
  // Observed state is polymorphic.  We only need to set the
  // partial likelihoods for states that are also
  // polymorphic for the same alleles.  E.g., states of type
  // (ix,(N-i)y) can lead to the observed state (jx,(M-j)y).
  else {
    int k;
    if (id1 == 0) k = id2 - 1;
    else k = id1 + id2;
    int real_state = nnuc + k*(N-1);
    for (int i = 1; i < N; i++, real_state++) {
      ASSERT(real_state < nstates);
      if (hypergeometric) {
          lh[real_state] = hypergeometric_dist(j, M, i, N);
      }
      else {
          lh[real_state] = binomial_dist(j, M, (double)i / (double)N);
      }
    }
  }

  // cout << "Sample M,j,id1,id2: " << M << ", " << j << ", " << id1 << ", " << id2 << endl;
  // cout << "Real partial likelihood vector: ";
  // for (i=0; i<4; i++)
  //   cout << real_partial_lh[i] << " ";
  // cout << endl;
  // for (i=4; i<nstates; i++) {
  //   cout << real_partial_lh[i] << " ";
  //   if ((i-3)%(N-1) == 0)
  //     cout << endl;
  // }
  // double sum = 0.0;
  // for (i=0; i<nstates; i++)
  //   sum += real_partial_lh[i];
  // cout << sum << endl;
}
