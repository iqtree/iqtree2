#include "modelpomo.h"
#include "modeldna.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// TODO DS: Refactor fixed -> boundary (clashes with fixing a
// parameter).

ModelPoMo::ModelPoMo(PhyloTree *tree) : ModelMarkov(tree) {
}

ModelPoMo::ModelPoMo(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params)
    // Do not count rates; does not make sense for PoMo.
    : ModelMarkov(tree) {
    init(model_name, model_params, freq_type, freq_params, is_reversible, pomo_params);
}

void ModelPoMo::init_substitution_model(const char *model_name,
                                        string model_params,
                                        StateFreqType freq_type,
                                        string freq_params,
                                        string pomo_params)
{
    // Get DNA model info from model_name.  Use ModelDNA for this
    // purpose.  It acts as the basis of the `ModelPoMo' (the mutation
    // coefficients point to the rates of ModelDNA, the fixed state
    // frequencies to the state frequencies and so on).

    // Trick ModelDNA constructor by setting the number of states to 4 (DNA).
    phylo_tree->aln->num_states = nnuc;

    // TODO DS: Add support for general Markov models.
    // string model_str = model_name;
    // if (ModelMarkov::validModelName(model_str))
    //     dna_model = ModelMarkov::getModelByName(model_str, phylo_tree, model_params, freq_type, freq_params);
    // else
    dna_model = new ModelDNA(model_name, model_params, freq_type, freq_params, phylo_tree);
    // Reset the number of states.
    phylo_tree->aln->num_states = num_states;

    this->name = dna_model->name;
    if (model_params.length() > 0)
        this->name += "{" + model_params + "}";
    this->name += "+rP";
    if (pomo_params.length() > 0)
        this->name += "{" + pomo_params + "}";
    this->name += "+N" + convertIntToString(N);
}

void ModelPoMo::init_sampling_type()
{
    sampling_type = phylo_tree->aln->pomo_sampling_type;
    if (sampling_type == SAMPLING_SAMPLED)
        this->name += "+S";
    else if (sampling_type == SAMPLING_WEIGHTED)
        this->name += "+W";
    else outError("Sampling type is not supported.");

    string sampling_method;
    if (sampling_type == SAMPLING_SAMPLED)
        sampling_method = "Sampled";
    else if (sampling_type == SAMPLING_WEIGHTED)
        sampling_method = "Weighted";
    else outError("Sampling type is not supported.");
    this->full_name =
        "Reversible PoMo with N=" +
        convertIntToString(N) + " and " +
        dna_model->full_name + " substitution model; " +
        "Sampling method: " + sampling_method + "; " +
        convertIntToString(num_states) + " states in total.";
}

void ModelPoMo::init_freq()
{
    // Frequencies of the boundary states (fixed states, e.g., 10A).
    // These correspond to the state frequencies in the DNA
    // substitution models.
    freq_fixed_states = dna_model->state_freq;
    freq_fixed_states_emp = new double[4];
    // Get the fixed state frequencies from the data.
    estimateEmpiricalFixedStateFreqs(freq_fixed_states_emp);

    // Get frequency type from substitution model.
    freq_type = dna_model->freq_type;

    // Handle frequency type.  This cannot be done by the underlying
    // substitution model because interpretation of polymorphic states
    // is undefined.
    switch (freq_type) {
    case FREQ_EQUAL:
        // '+FQ'.
        for (int i = 0; i < 4; i++)
            freq_fixed_states[i] = 1.0/ (double)nnuc;
        break;
    case FREQ_ESTIMATE:
        // '+FO'.  Start estimation at empirical frequencies.
        for (int i = 0; i < 4; i++)
            freq_fixed_states[i] = freq_fixed_states_emp[i];
        break;
    case FREQ_EMPIRICAL:
        // '+F'.
        for (int i = 0; i < 4; i++)
            freq_fixed_states[i] = freq_fixed_states_emp[i];
        break;
    case FREQ_USER_DEFINED:
        // '+FU'. ModelDNA should have set them already.
        if (freq_fixed_states[0] == 0.0)
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
                                      string pomo_params)
{
    // TODO DS: Enable constrained maximization of likelihood and
    // fixed level of polymorphisms.  So far either all parameters are
    // fixed or none.  If separate treatment is possible, this
    // function may be split into two functions: (1) initialization of
    // parameters of underlying substitution model, (2) initizliation
    // of level of polymorphism.

    fixed_model_params = false;
    fixed_theta = false;
    if (model_params.length() > 0)
        fixed_model_params = true;
    if (pomo_params.length() > 0)
        fixed_theta = true;

    if (fixed_model_params != fixed_theta) {
        cout << "Either fix all parameters, e.g., with \"-m HKY{2.0}+rP{0.0025}\" or none." << endl;
        outError("Abort.");
    }

    // Should be set by underlying Markov model anyways.
    // if (!fixed_model_params) {
    //     // In case the fixation of parameters can be separated,
    //     // the first element should not be set here because it is
    //     // the level of polymorphism.
    //     for (vector<bool>::iterator p_fixed = dna_model->param_fixed.begin();
    //          p_fixed != dna_model->param_fixed.end(); p_fixed++)
    //         *p_fixed = false;
    // }

    if (fixed_theta) {
        cout << setprecision(5);
        if (pomo_params == "EMP") {
            cout << "Level of polymorphism will be fixed to the estimate from the data: ";
            cout << theta << "." << endl;
            // No need to set the level of polymorphism here because
            // of initialization.
        }
        else {
            cout << "Level of polymorphism will be fixed to the value given by the user: ";
            theta = convert_double(pomo_params.c_str());
            cout << theta << "." << endl;
        }
        // TODO DS: Is this needed?  It is better to save the level of
        // polymorphism not with the dna model class but in the PoMo
        // class.
        dna_model->param_fixed.front() = true;
    }
}


void ModelPoMo::init(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     bool is_reversible,
                     string pomo_params) {
    // Initialize model constants.
    N = phylo_tree->aln->virtual_pop_size;
    nnuc = 4;
    n_connections = nnuc * (nnuc-1) / 2;
    eps = 1e-6;
    // Check if number of states of PoMo match the provided data.
    assert(num_states == (nnuc + (nnuc*(nnuc-1)/2 * (N-1))) );

    // TODO DS: Remove reversibility flag.  The reversibility is
    // determined by the underlying substitution model.
    if (is_reversible != true) throw "Non-reversible PoMo not supported yet.";

    // Main initialization of model and parameters.
    init_substitution_model(model_name,
                            model_params,
                            freq_type,
                            freq_params,
                            pomo_params);
    init_sampling_type();
    init_freq();
    rate_matrix = new double[num_states*num_states];
    theta = estimateEmpiricalWattersonTheta();
    init_fixed_parameters(model_params, pomo_params);
    setInitialMutCoeff();
    updatePoMoStatesAndRateMatrix();
    decomposeRateMatrix();

    cout << "Initialized PoMo model." << endl;
    cout << "Model name: " << this->name << "." << endl;
    cout << this->full_name << endl;
    if (verbose_mode >= VB_MAX)
        writeInfo(cout);
}

ModelPoMo::~ModelPoMo() {
    delete [] rate_matrix;
//  delete [] freq_fixed_states;
    delete dna_model;
    delete [] freq_fixed_states_emp;
    // delete [] fixed_model_params_ratio;
    }

double ModelPoMo::computeSumFreqFixedStates() {
    int i;
    double norm_fixed = 0.0;
    for (i = 0; i < 4; i++)
        norm_fixed += freq_fixed_states[i];
    assert(norm_fixed == 1.0);
    return norm_fixed;
}

// Give back the harmonic number of n.
double harmonic(int n) {
    double harmonic = 0.0;
    for (int i = 1; i <= n; i++)
        harmonic += 1.0/(double)i;
    return harmonic;
}

void ModelPoMo::setInitialMutCoeff() {
    // Mutation rates point to the rates of the DNA model.
    mutation_rates = new double[n_connections];
    for (int i = 0; i < n_connections; i++)
        mutation_rates[i] = dna_model->rates[i];
    // for (int i = 0; i < 6; i++) mutation_rates[i] = POMO_INIT_RATE;
    // double m_init = 0;
    // double theta_p = theta;
    // double lambda_fixed_sum = computeSumFreqFixedStates();
    double lambda_poly_sum_no_mu = computeSumFreqPolyStatesNoMut();
    // // cout << "DEBUG" << endl;
    // // cout << theta_p << endl;
    // // cout << lambda_fixed_sum << endl;
    // // cout << lambda_poly_sum_no_mu << endl;

    if (!fixed_theta && lambda_poly_sum_no_mu <= 0) {
        // TODO DS: Obsolete because mutation rates are scaled
        // according to theta.  What needs to be done is setting an
        // predefined theta.
        outWarning("We strongly discourage to use PoMo on data without polymorphisms.");
        outWarning("Set initial rates to predefined values.");
        for (int i = 0; i < 6; i++) mutation_rates[i] = POMO_INIT_RATE;
        return;
    }

    normalizeMutationRates();
    // m_init = theta_p * lambda_fixed_sum / (lambda_poly_sum_no_mu * (1.0 - theta_p));
    // if (m_init < POMO_MIN_RATE || m_init > POMO_MAX_RATE)
    //     outError("Initial rate not within boundaries.  Please check data.");
    // // cout << "DEBUG" << endl;
    // // cout << m_init << endl;
    // // Honor fixed rate specifications.
    // double sum = 0;
    // int n_mu = 6;
    // for (int i = 0; i < n_mu; i++) sum += mutation_rates[i];
    // for (int i = 0; i < n_mu; i++) {
    //     double new_mut_prob = m_init*mutation_rates[i]*n_mu/sum;
    //     mutation_rates[i] = new_mut_prob;
    // }
}

double ModelPoMo::computeSumFreqPolyStatesNoMut() {
    double norm_polymorphic = 0.0;
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < i; j++)
            norm_polymorphic +=
                2 * freq_fixed_states[i] * freq_fixed_states[j];
    }
    norm_polymorphic *= harmonic(N-1);
    return norm_polymorphic;
}

double ModelPoMo::computeSumFreqPolyStates() {
    double norm_polymorphic = 0.0;
    int i, j;
    for (i = 0; i < nnuc; i++) {
        for (j = 0; j < i; j++)
            norm_polymorphic +=
                2 * freq_fixed_states[i] * freq_fixed_states[j] * mutCoeff(i, j);
    }
    norm_polymorphic *= harmonic(N-1);
    return norm_polymorphic;
}

double ModelPoMo::computeNormConst() {
    double norm_fixed = computeSumFreqFixedStates();
    double norm_polymorphic = computeSumFreqPolyStates();
    return 1.0/(norm_fixed + norm_polymorphic);
}

// TODO DS: Add non-reversible stuff.
void ModelPoMo::computeStateFreq () {
    double norm = computeNormConst();
    int state;

    for (state = 0; state < num_states; state++) {
        if (isFixed(state))
            state_freq[state] = freq_fixed_states[state]*norm;
        else {
            int k, X, Y;
            decomposeState(state, k, X, Y);
            state_freq[state] =
                norm * freq_fixed_states[X] * freq_fixed_states[Y] *
                mutCoeff(X, Y)*N / (k*(N-k));
        }
    }
}

void ModelPoMo::updatePoMoStatesAndRateMatrix () {
    // TODO DS: Check this. This comment is obsolete.
    // Activate this if frequencies of fixed states sum up to 1.0.
    // updateFreqFixedState();
    computeStateFreq();

    // Compute and normalzie the rate matrix such that on average one
    // event happens per delta_t = 1.0.  This seems to be stable.
    int i, j;
    double tot_sum = 0.0;
    for (i = 0; i < num_states; i++) {
        double row_sum = 0.0;
        // Loop over columns in row state1 (transition to state2).
        for (j = 0; j < num_states; j++)
            if (i != j) {
                row_sum +=
                    (rate_matrix[i*num_states+j] =
                     computeProbBoundaryMutation(i, j));
            }
        tot_sum += state_freq[i]*row_sum;
        rate_matrix[i*num_states+i] = -(row_sum);
    }
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            rate_matrix[i*num_states+j] /= tot_sum;
        }
    }
}

void ModelPoMo::decomposeState(int state, int &i, int &nt1, int &nt2) {
    if (state < 4) {
        // Fixed A, C, G or T
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

bool ModelPoMo::isFixed(int state) {
    return (state < 4);
}

bool ModelPoMo::isPolymorphic(int state) {
    return (!isFixed(state));
}

double ModelPoMo::mutCoeff(int nt1, int nt2) {
    assert(nt1!=nt2 && nt1<4 && nt2<4);
    if (nt2 < nt1) {
        int tmp=nt1;
        nt1=nt2;
        nt2=tmp;
    }
    if (nt1==0) return mutation_rates[nt2-1];
    if (nt1==1) return mutation_rates[nt2+1];
    if (nt1==2) return mutation_rates[5];
    assert(0);
}

double ModelPoMo::computeProbBoundaryMutation(int state1, int state2) {
    // The transition rate to the same state will be calculated by
    // setting the row sum to 0.
    assert(state1 != state2);

    // Both states are decomposed into the abundance of the first
    // allele as well as the nucleotide of the first and the second
    // allele.
    int i1=0, i2=0, nt1=-1, nt2=-1, nt3=-1, nt4=-1;
    decomposeState(state1, i1, nt1, nt2);
    decomposeState(state2, i2, nt3, nt4);

    // Either the first nucleotides match or the first of state 1 with
    // the second of state 2 or the first of state 2 with the second
    // of state 1.  Additionally, we have to consider fixed states as
    // special cases.
    if (nt1 == nt3 && (nt2==nt4 || nt2==-1 || nt4 == -1)) {
        assert(i1 != i2); // because state1 != state2
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
                return mutCoeff(nt1,nt4) * freq_fixed_states[nt4];
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
        return mutCoeff(nt1,nt3) * freq_fixed_states[nt3];
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
    if (fixed_theta)
        return dna_model->getNDim();
    else
        return dna_model->getNDim()+1;
}

int ModelPoMo::getNDimFreq() {
    return dna_model->getNDimFreq();
}

void ModelPoMo::setBounds(double *lower_bound,
                          double *upper_bound,
                          bool *bound_check) {
    // int i, ndim = getNDim();
    // if (verbose_mode >= VB_MAX)
    //     cout << "Set new bounds." << endl;

    // // Number of inferred substitution rates of underlying
    // // substitution model.
    // int n_rates = dna_model->num_params;    
    // for (int i = 1; i <= n_rates; i++) {
    //     lower_bound[i] = POMO_MIN_RATE;
    //     upper_bound[i] = POMO_MAX_RATE;
    //     bound_check[i] = false;
    // }

    // // Frequencies of fixed states.
    // if (freq_type == FREQ_ESTIMATE) {
    //     for (int i = n_rates+1; i <= n_rates+1+nnuc-1; i++) {
    //         lower_bound[i] = POMO_MIN_REL_FREQ * freq_fixed_states[i-n_rates-1];
    //         upper_bound[i] = POMO_MAX_REL_FREQ * freq_fixed_states[i-n_rates-1];
    //         bound_check[i] = false;
    //     }
    // }

    // Does not work because boundaries are too high.
    // Set boundaries of underlying substitution model.
    dna_model->setBounds(lower_bound, upper_bound, bound_check);
    
    // Level of polymorphism.
    if (!fixed_theta) {
        int ndim = getNDim();
        lower_bound[ndim] = POMO_MIN_THETA;
        upper_bound[ndim] = POMO_MAX_THETA;
        bound_check[ndim] = false;
    }
}

void ModelPoMo::normalizeMutationRates() {
    // Get rates from underlying substitution model.    
    double sum = 0.0;
    for (int i = 0; i < n_connections; i++) {
        mutation_rates[i] = dna_model->rates[i];
        sum += mutation_rates[i];
    }
    // TODO DS: This may not be necessary.  Remove theta bug first.
    for (int i = 0; i < n_connections; i++) {
        mutation_rates[i] /= sum;
        cout << mutation_rates[i] << endl;
    }
        
    // Normalize the mutation probability so that they resemble the
    // given level of polymorphism.
    computeStateFreq();
    double sum_pol = computeSumFreqPolyStates();
    double sum_fix = computeSumFreqFixedStates();
    double m_norm  = sum_pol * (1.0 - theta) / (sum_fix * theta);
    if (verbose_mode >= VB_MAX)
        cout << "Normalization constant of mutation rates: " << m_norm << endl;
    cout << "Normalization constant of mutation rates: " << m_norm << endl;
    for (int i = 0; i < n_connections; i++) {
        mutation_rates[i] /= m_norm;
        if (mutation_rates[i] <= POMO_MIN_RATE) {
            cout << mutation_rates[i] << endl;
            outWarning("Mutation rate below boundary after normalization.");
        }
        if (mutation_rates[i] >= POMO_MAX_RATE)
            outWarning("Mutation rate above boundary after normalization.");
    }
    // Reocmpute stationary frequency vector with updated mutation
    // rates.
    computeStateFreq();
}

void ModelPoMo::scaleMutationRatesAndUpdateRateMatrix(double scale) {
    for (int i = 0; i < n_connections; i++) {
        mutation_rates[i] = mutation_rates[i]*scale;
    }
    updatePoMoStatesAndRateMatrix();
}

bool ModelPoMo::getVariables(double *variables) {
    // A drawback: IQ-TREE looks at the values at the boundaries a
    // lot.  E.g., if HKY{3.0} is specified, the mutation_rates[0] will
    // be set to 0.01=BOUNDARY_MAX a lot, but then, e.g.,
    // mutation_rates[1] will be set to 3*0.01 which is outside the
    // boundaries.  How should this be handled?

    // int i;
    // Use a zero indexed array (variables starts at one).
    // double *vars = variables+1;

    bool changed = false;
    changed = dna_model->getVariables(variables);

    // for (i = 0; i < num_all; i++) {
    //     if (!dna_model->param_fixed[dna_model->param_spec[i]]) {
    //         changed |= (mutation_rates[i] != vars[(int)dna_model->param_spec[i]]);
    //         mutation_rates[i] = vars[(int)dna_model->param_spec[i]];
    //     }
    // }

    // int ndim = getNDim();
    // if (freq_type == FREQ_ESTIMATE) {
    //     changed |= memcmpcpy(freq_fixed_states, vars+(ndim-nnuc+1), (nnuc-1)*sizeof(double));
    //     if (verbose_mode >= VB_MAX) {
    //         for (i = 0; i < nnuc-1; i++) {
    //             cout << setprecision(8);
    //             cout << "  Estimated fixed frequencies[" << i << "] = ";
    //             cout << vars[ndim-nnuc+1+i] << endl;
    //         }
    //     }
    // }

    if (!fixed_theta) {
        int ndim = getNDim();
        changed |= (theta != variables[ndim+1]);
        theta = variables[ndim];
    }

    normalizeMutationRates();
    updatePoMoStatesAndRateMatrix();
    return changed;
}

void ModelPoMo::setVariables(double *variables) {
    dna_model->setVariables(variables);

    if (!fixed_theta) {
        int ndim = getNDim();
        variables[ndim] = theta;
    }
    // // Use a zero indexed array (variables starts at one).
    // double * vars = variables+1;
    // int num_all = dna_model->param_spec.length();
    // for (int i = 0; i < num_all; i++)
    //     if (!dna_model->param_fixed[dna_model->param_spec[i]])
    //         vars[(int)dna_model->param_spec[i]] = mutation_rates[i];

    // if (freq_type == FREQ_ESTIMATE) {
    //     int ndim = getNDim();
    //     memcpy(vars+(ndim-nnuc+1), freq_fixed_states, (nnuc-1)*sizeof(double));
    // }
}

void ModelPoMo::writeInfo(ostream &out) {
    int i;
    ios  state(NULL);
    state.copyfmt(out);

    out << setprecision(8);

    out << "Frequency of fixed states: ";
    for (i = 0; i < 4; i++)
        out << freq_fixed_states[i] << " ";
    out << endl;
    out << "Mutation rates: ";
    for (i = 0; i < 6; i++)
        out << mutation_rates[i] << " ";
    out << endl;

    out.copyfmt(state);
}

// TODO: ModelGTR::decomposeRateMatrix() calls
// EigenDecomposition::eigensystem_sym() which in turn calls
// ModelPoMo::computeRateMatrix().  I believe that there is
// unnecessary allocation of rate parameters involved in
// ModelGTR::decomposeRateMatrix() which can be removed.
void ModelPoMo::computeRateMatrix(double **r_matrix, double *s_freqs, int n_states) {
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            r_matrix[i][j] = rate_matrix[i*n_states+j];
        }
    }
}

double ModelPoMo::targetFunk(double x[]) {
    getVariables(x);
    // Disable test for low stationary frequency in PoMo.
    // if (state_freq[num_states-1] < 1e-4) return 1.0e+12;
    decomposeRateMatrix();
    assert(phylo_tree);
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

void
ModelPoMo::estimateEmpiricalFixedStateFreqs(double * freq_fixed_states)
{
    memset(freq_fixed_states, 0, sizeof(double)*nnuc);

    if (sampling_type == SAMPLING_SAMPLED) {
        unsigned int abs_state_freq[num_states];
        memset(abs_state_freq, 0, sizeof(unsigned int)*num_states);
        phylo_tree->aln->computeAbsoluteStateFreq(abs_state_freq);
        int n;
        int x;
        int y;

        int sum[nnuc];
        int tot_sum = 0;
        memset (sum, 0, nnuc * sizeof(int));

        for (int i = 0; i < num_states; i++) {
            decomposeState(i, n, x, y);
            sum[x]+= n*abs_state_freq[i];
            if (y >= 0) sum[y]+= (N-n)*abs_state_freq[i];
        }
        for (int i = 0; i < nnuc; i++) {
            tot_sum += sum[i];
        }
        for (int i = 0; i < nnuc; i++) {
            freq_fixed_states[i] = (double) sum[i]/tot_sum;
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
                assert((unsigned int)state < phylo_tree->aln->pomo_states.size());
                // Decode the id and counts.
                int id1 = phylo_tree->aln->pomo_states[state] & 3;
                int id2 = (phylo_tree->aln->pomo_states[state] >> 16) & 3;
                int j1 = (phylo_tree->aln->pomo_states[state] >> 2) & 16383;
                int j2 = (phylo_tree->aln->pomo_states[state] >> 18);
                freq_fixed_states[id1] += j1*(it->frequency);
                freq_fixed_states[id2] += j2*(it->frequency);
            }
        }
    }
    // Normalize frequencies so that they sum to 1.0.
    double sum = 0.0;
    for (int i = 0; i < nnuc; i++)
        sum += freq_fixed_states[i];
    for (int i = 0; i < nnuc; i++)
        freq_fixed_states[i] /= sum;
    if (verbose_mode >= VB_MAX) {
        std::cout << "The empirical frequencies of the fixed states are:" << std::endl;
        for (int i = 0; i < nnuc; i++)
            std::cout << freq_fixed_states[i] << " ";
        std::cout << std::endl;
    }
}

double
ModelPoMo::estimateEmpiricalWattersonTheta()
{
    double theta_p = 0.0;
    int sum_pol = 0;
    int sum_fix = 0;
    double sum_theta_w = 0.0;

    if (sampling_type == SAMPLING_SAMPLED) {
        unsigned int abs_state_freq[num_states];
        memset(abs_state_freq, 0, sizeof(unsigned int)*num_states);
        phylo_tree->aln->computeAbsoluteStateFreq(abs_state_freq);
        for (int i = 0; i < nnuc; i++) sum_fix += abs_state_freq[i];
        for (int i = nnuc; i < num_states; i++) sum_pol += abs_state_freq[i];
        theta_p = (double) sum_pol / (double) (sum_fix + sum_pol);
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
                assert((unsigned int)state < phylo_tree->aln->pomo_states.size());
                // Decode the id and counts.
                // int id1 = phylo_tree->aln->pomo_states[state] & 3;
                // int id2 = (phylo_tree->aln->pomo_states[state] >> 16) & 3;
                int j1 = (phylo_tree->aln->pomo_states[state] >> 2) & 16383;
                int j2 = (phylo_tree->aln->pomo_states[state] >> 18);
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
        // // DEBUG.
        // cout << setprecision(8);
        // cout << "DEBUG ==========: " << endl;
        // cout << "Estimated Watterson's Theta: ";
        // cout << sum_theta_w << endl;

        // Wed Jul 6 15:03:06 CEST 2016; I think Watterson's theta is
        // without * harmonic(N-1).
        // theta_p = sum_theta_w * harmonic(N-1);
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

void ModelPoMo::report_rates(ostream &out) {
    out << setprecision(8);
    out << "Mutation rates (in the order AC, AG, AT, CG, CT, GT):" << endl;
    double sum = 0.0;
    for (int i = 0; i < 6; i++)
        out << mutation_rates[i] << " ";
    out << endl;
}

void ModelPoMo::report(ostream &out) {
    out << "Reversible PoMo." << endl;
    out << "Virtual population size N: " << N << endl;
    if (sampling_type == SAMPLING_SAMPLED)
        out << "Sampling method: Sampled." << endl;
    else
        out << "Sampling method: Weighted." << endl;

    out << endl;
    out << "Estimated quantities" << endl;
    out << "--------------------" << endl;

    if (freq_type == FREQ_ESTIMATE) {
        out << "Frequencies of fixed states (in the order A, C, G T):" << endl;
        for (int i = 0; i < nnuc; i++)
            out << freq_fixed_states[i] << " ";
        out << endl;
    }
    report_rates(out);

    double poly = computeSumFreqPolyStates();
    double fixed = computeSumFreqFixedStates();
    double prop_poly = poly / (poly + fixed);
    // Thu Aug 18 22:00:49 CEST 2016; I realized that prop_poly is
    // already the Watterson's theta, so the following calculation was
    // wrong!
    // double watterson_theta = prop_poly / harmonic(N-1);
    double emp_watterson_theta = estimateEmpiricalWattersonTheta();
    out << setprecision(8);
    // out << "Estimated sum of fixed states:" << endl;
    // out << fixed << endl;
    // out << "Estimated sum of polymorphic states" << endl;
    // out << poly << endl;

    // Thu Aug 18 22:01:35 CEST 2016; see about ten lines above.
    if (!fixed_theta)
        out << "Watterson's theta: " << prop_poly << endl;

    // out << "(Estimated) proportion of polymorphic states:" << endl;
    // out << prop_poly << endl;
    // out << "Watterson Theta: " << watterson_theta << endl;

    out << endl;
    out << "Empirical quantities" << endl;
    out << "--------------------" << endl;

    out << "Frequencies of fixed states (in the order A, C, G, T):" << endl;
    for (int i = 0; i < nnuc; i++)
        out << freq_fixed_states_emp[i] << " ";
    out << endl;

    out << "Watterson's Theta: " << emp_watterson_theta << endl;
    out << endl;
    // out << "(Estimated) frequencies of fixed states:" << endl;
    // for (int i = 0; i < nnuc; i++)
    //     out << freq_fixed_states[i] << " ";
    // out << endl << endl;

    // out << "Total rate: ";
    // double tot_sum = 0.0;
    // for (int i = 0; i < num_states; i++) {
    //     double row_sum;
    //     row_sum = 0.0;
    //     for (int j = 0; j < num_states; j++) {
    //         if (i != j) row_sum += rate_matrix[i*num_states + j];
    //     }
    //     tot_sum += state_freq[i]*row_sum;
    // }
    // out << tot_sum << endl;

    // // Output expected number of substitutions without
    // // total_tree_length/t_fix.
    // out << "Total rate (mutations only): ";
    // double mu_tot = 0.0;
    // // for (int i = 0; i < nnuc; i++)
    // //     for (int j = 0; j < nnuc; j++)
    // //         if (i != j) mu_tot1 += (double) freq_fixed_states[i]*freq_fixed_states[j]*mutCoeff(i,j);
    // for (int i = 0; i < nnuc; i++)
    //     for (int j = 0; j < num_states; j++)
    //         // if (i != j) mu_tot += (double) freq_fixed_states[i]*rate_matrix[i*num_states + j];
    //         if (i != j) mu_tot += (double) state_freq[i]*rate_matrix[i*num_states + j];
    // out << mu_tot << endl;

    // double fr_tot = 0.0;
    // out << "Total rate (drift or frequency shifts only): ";
    // // fr_tot = tot_sum-mu_tot;
    // for (int i = nnuc; i < num_states; i++) {
    //     double row_sum;
    //     row_sum = 0.0;
    //     for (int j = 0; j < num_states; j++) {
    //         if (i != j) row_sum += rate_matrix[i*num_states +j];
    //     }
    //     fr_tot += state_freq[i]*row_sum;
    // }
    // out << fr_tot << endl;

    // out << "Ratio of mutation to drift (indicates if boundary mutation model is adequate): ";
    // out << mu_tot / fr_tot << endl;

    // out << "Ratio of mutation rate to total rate: ";
    // out << mu_tot / (fr_tot + mu_tot) << endl << endl;

}

void ModelPoMo::saveCheckpoint() {
    int n_rates = nnuc * (nnuc-1) / 2;
    checkpoint->startStruct("ModelPoMo");
    CKP_ARRAY_SAVE(n_rates, dna_model->rates);
    CKP_ARRAY_SAVE(nnuc, dna_model->state_freq);
    checkpoint->endStruct();
    ModelMarkov::saveCheckpoint();
}

void ModelPoMo::restoreCheckpoint() {
    int n_rates = nnuc * (nnuc-1) / 2;
    // First, get variables from checkpoint.
    checkpoint->startStruct("ModelPoMo");
    CKP_ARRAY_RESTORE(n_rates, dna_model->rates);
    CKP_ARRAY_RESTORE(nnuc, dna_model->state_freq);
    checkpoint->endStruct();
    // Second, update states and rate matrix.
    updatePoMoStatesAndRateMatrix();
    // Third, restore ModelGTR.
    ModelMarkov::restoreCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}
