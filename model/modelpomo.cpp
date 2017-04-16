#include "modelpomo.h"
#include "modeldna.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

ModelPoMo::ModelPoMo(PhyloTree *tree) : ModelMarkov(tree) {
}

ModelPoMo::ModelPoMo(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params)
    : ModelMarkov(tree) {
    init(model_name, model_params, freq_type, freq_params, is_reversible, pomo_params);
}

void ModelPoMo::init_mutation_model(const char *model_name,
                                        string model_params,
                                        StateFreqType freq_type,
                                        string freq_params,
                                        string pomo_params)
{
    // TODO DS: Add support for general Markov models.
    // Trick ModelDNA constructor by setting the number of states to 4 (DNA).
    phylo_tree->aln->num_states = nnuc;
    // string model_str = model_name;
    // if (ModelMarkov::validModelName(model_str))
    //     dna_model = ModelMarkov::getModelByName(model_str, phylo_tree, model_params, freq_type, freq_params);
    // else
    dna_model = new ModelDNA(model_name, model_params, freq_type, freq_params, phylo_tree);
    // Reset the number of states.
    phylo_tree->aln->num_states = num_states;

    // TODO DS: Remove rP flag.
    this->name = dna_model->name;
    if (model_params.length() > 0)
        this->name += "{" + model_params + "}";
    this->name += "+rP";
    if (pomo_params.length() > 0)
        this->name += "{" + pomo_params + "}";
    this->name += "+N" + convertIntToString(N);
}

void ModelPoMo::init_sampling_method()
{
    sampling_method = phylo_tree->aln->pomo_sampling_method;
    string sampling_method_str;
    if (sampling_method == SAMPLING_SAMPLED) {
        this->name += "+S";
        sampling_method_str = "Sampled";
    }
    else if (sampling_method == SAMPLING_WEIGHTED) {
        this->name += "+W";
        sampling_method_str = "Weighted";
    }
    else outError("Sampling type is not supported.");

    this->full_name =
        "Reversible PoMo with N=" +
        convertIntToString(N) + " and " +
        dna_model->full_name + " mutation model; " +
        "Sampling method: " + sampling_method_str + "; " +
        convertIntToString(num_states) + " states in total.";
}

void ModelPoMo::init_boundary_frequencies()
{
    // Get boundary state frequencies from underlying mutation model.
    freq_boundary_states = dna_model->state_freq;
    // Get the empirical boundary state frequencies from the data.
    freq_boundary_states_emp = new double[4];
    estimateEmpiricalBoundaryStateFreqs(freq_boundary_states_emp);
    // Get frequency type from mutation model.
    freq_type = dna_model->freq_type;

    // Handle frequency type.  This cannot be done by the underlying
    // mutation model because interpretation of polymorphic states is
    // undefined.
    switch (freq_type) {
    case FREQ_EQUAL:
        // '+FQ'.
        for (int i = 0; i < 4; i++)
            freq_boundary_states[i] = 1.0/ (double)nnuc;
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
                                      string pomo_params)
{
    // TODO DS: Enable constrained maximization of likelihood and
    // fixed level of polymorphisms.  So far either all parameters are
    // fixed or none.  If separate treatment is possible, this
    // function may be split into two functions: (1) initialization of
    // parameters of underlying mutation model, (2) initizliation
    // of level of polymorphism.
    fixed_model_params = false;
    fixed_theta_emp = false;
    fixed_theta_usr = false;
    fixed_theta = false;
    if (model_params.length() > 0)
        fixed_model_params = true;
    if (pomo_params.length() > 0) {
        cout << setprecision(5);
        fixed_theta = true;
        if (pomo_params == "EMP") {
            cout << "Level of polymorphism is fixed to the estimate from the data: ";
            cout << theta << "." << endl;
            fixed_theta_emp = true;
            // No need to set the level of polymorphism here because
            // of initialization.
        }
        else {
            cout << "Level of polymorphism is fixed to the value given by the user: ";
            theta = convert_double(pomo_params.c_str());
            cout << theta << "." << endl;
            fixed_theta_usr = true;
        }
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
    // determined by the underlying mutation model.
    if (is_reversible != true) throw "Non-reversible PoMo not supported yet.";

    // Main initialization of model and parameters.
    init_mutation_model(model_name,
                        model_params,
                        freq_type,
                        freq_params,
                        pomo_params);
    init_sampling_method();
    init_boundary_frequencies();
    theta = estimateEmpiricalWattersonTheta();
    init_fixed_parameters(model_params, pomo_params);
    setInitialMutCoeff();
    rate_matrix = new double[num_states*num_states];
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
    delete dna_model;
    delete [] freq_boundary_states_emp;
    delete [] mutation_rates;
}

double ModelPoMo::computeSumFreqBoundaryStates() {
    int i;
    double norm_boundary = 0.0;
    for (i = 0; i < 4; i++)
        norm_boundary += freq_boundary_states[i];
    return norm_boundary;
}

double harmonic(int n) {
    double harmonic = 0.0;
    for (int i = 1; i <= n; i++)
        harmonic += 1.0/(double)i;
    return harmonic;
}

void ModelPoMo::setInitialMutCoeff() {
    // Mutation rates point to the rates of the DNA model.
    mutation_rates = new double[n_connections];

    // TODO DS: Check if this works.
    // Check if polymorphism data is available.
    double lambda_poly_sum_no_mu = computeSumFreqPolyStatesNoMut();
    if (!fixed_theta && lambda_poly_sum_no_mu <= 0) {
        // TODO DS: Obsolete because mutation rates are scaled
        // according to theta.  What needs to be done is setting a
        // predefined theta.
        outWarning("We strongly discourage to use PoMo on data without polymorphisms.");
        outError("Setting the level of polymorphism without population data is not yet supported.");
    }

    normalizeMutationRates();
}

double ModelPoMo::computeSumFreqPolyStatesNoMut() {
    double norm_polymorphic = 0.0;
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < i; j++)
            norm_polymorphic +=
                2 * freq_boundary_states[i] * freq_boundary_states[j];
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
                2 * freq_boundary_states[i] * freq_boundary_states[j] * mutCoeff(i, j);
    }
    norm_polymorphic *= harmonic(N-1);
    return norm_polymorphic;
}

double ModelPoMo::computeNormConst() {
    double norm_boundary = computeSumFreqBoundaryStates();
    double norm_polymorphic = computeSumFreqPolyStates();
    return 1.0/(norm_boundary + norm_polymorphic);
}

// TODO DS: Add non-reversible stuff.
void ModelPoMo::computeStateFreq () {
    double norm = computeNormConst();
    int state;

    for (state = 0; state < num_states; state++) {
        if (isBoundary(state))
            state_freq[state] = freq_boundary_states[state]*norm;
        else {
            int k, X, Y;
            decomposeState(state, k, X, Y);
            state_freq[state] =
                norm * freq_boundary_states[X] * freq_boundary_states[Y] *
                mutCoeff(X, Y)*N / (k*(N-k));
        }
    }
}

void ModelPoMo::updatePoMoStatesAndRateMatrix () {
    // TODO DS: Check this. This comment is obsolete.
    // Activate this if frequencies of boundary states sum up to 1.0.
    // updateFreqBoundaryState();
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
    // of state 1.  Additionally, we have to consider boundary states as
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
                return mutCoeff(nt1,nt4) * freq_boundary_states[nt4];
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
        return mutCoeff(nt1,nt3) * freq_boundary_states[nt3];
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
    // Set boundaries of underlying mutation model.
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
    // Get rates from underlying mutation model.    
    for (int i = 0; i < n_connections; i++) {
        mutation_rates[i] = dna_model->rates[i];
    }

    // Normalize the mutation probability so that they resemble the
    // given level of polymorphism.
    computeStateFreq();
    double poly = computeSumFreqPolyStates();
    double theta_bm = poly/harmonic(N-1);

    // See Eq. (12.14) in my thesis.

    // TODO DS: Why do I need this correction?  I do know that
    // sampling with replacement gives better results (in the thesis I
    // use sampling without replacement).
    double test = (double) (N-1) / double (N+1);
    double m_norm = theta / (theta_bm * (test - harmonic(N-1) * theta));

    if (verbose_mode >= VB_MAX)
        cout << "Normalization constant of mutation rates: " << m_norm << endl;

    for (int i = 0; i < n_connections; i++)
        mutation_rates[i] *= m_norm;
    
    // // Recompute stationary frequency vector with updated mutation
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
    bool changed = false;
    changed = dna_model->getVariables(variables);

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
}

void ModelPoMo::writeInfo(ostream &out) {
    int i;
    ios  state(NULL);
    state.copyfmt(out);

    out << setprecision(8);

    out << "Frequency of boundary states: ";
    for (i = 0; i < 4; i++)
        out << freq_boundary_states[i] << " ";
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
ModelPoMo::estimateEmpiricalBoundaryStateFreqs(double * freq_boundary_states)
{
    memset(freq_boundary_states, 0, sizeof(double)*nnuc);

    if (sampling_method == SAMPLING_SAMPLED) {
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
                assert((unsigned int)state < phylo_tree->aln->pomo_states.size());
                // Decode the id and counts.
                int id1 = phylo_tree->aln->pomo_states[state] & 3;
                int id2 = (phylo_tree->aln->pomo_states[state] >> 16) & 3;
                int j1 = (phylo_tree->aln->pomo_states[state] >> 2) & 16383;
                int j2 = (phylo_tree->aln->pomo_states[state] >> 18);
                freq_boundary_states[id1] += j1*(it->frequency);
                freq_boundary_states[id2] += j2*(it->frequency);
            }
        }
    }
    // Normalize frequencies so that they sum to 1.0.
    double sum = 0.0;
    for (int i = 0; i < nnuc; i++)
        sum += freq_boundary_states[i];
    for (int i = 0; i < nnuc; i++)
        freq_boundary_states[i] /= sum;
    if (verbose_mode >= VB_MAX) {
        std::cout << "The empirical frequencies of the boundary states are:" << std::endl;
        for (int i = 0; i < nnuc; i++)
            std::cout << freq_boundary_states[i] << " ";
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

    if (sampling_method == SAMPLING_SAMPLED) {
        unsigned int abs_state_freq[num_states];
        memset(abs_state_freq, 0, sizeof(unsigned int)*num_states);
        phylo_tree->aln->computeAbsoluteStateFreq(abs_state_freq);
        for (int i = 0; i < nnuc; i++) sum_fix += abs_state_freq[i];
        for (int i = nnuc; i < num_states; i++) sum_pol += abs_state_freq[i];
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
                assert((unsigned int)state < phylo_tree->aln->pomo_states.size());
                // Decode counts.
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
    if (sampling_method == SAMPLING_SAMPLED)
        out << "Sampling method: Sampled." << endl;
    else
        out << "Sampling method: Weighted." << endl;

    out << endl;
    out << "Estimated quantities" << endl;
    out << "--------------------" << endl;

    if (freq_type == FREQ_ESTIMATE) {
        out << "Frequencies of boundary states (in the order A, C, G T):" << endl;
        for (int i = 0; i < nnuc; i++)
            out << freq_boundary_states[i] << " ";
        out << endl;
    }
    report_rates(out);

    out << setprecision(8);
    
    if (!fixed_theta)
        out << "Estimated heterozygosity: ";
    else if (fixed_theta_emp)
        out << "Empirical heterozygosity: ";
    else if (fixed_theta_usr)
        out << "User-defined heterozygosity: ";
    out << theta << endl;
    
    out << endl;
    out << "Empirical quantities" << endl;
    out << "--------------------" << endl;

    out << "Frequencies of boundary states (in the order A, C, G, T):" << endl;
    for (int i = 0; i < nnuc; i++)
        out << freq_boundary_states_emp[i] << " ";
    out << endl;

    double emp_watterson_theta = estimateEmpiricalWattersonTheta();
    out << "Watterson's Theta: " << emp_watterson_theta << endl;
    out << endl;
    // out << "(Estimated) frequencies of boundary states:" << endl;
    // for (int i = 0; i < nnuc; i++)
    //     out << freq_boundary_states[i] << " ";
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

    // // Output expected number of mutations without
    // // total_tree_length/t_fix.
    // out << "Total rate (mutations only): ";
    // double mu_tot = 0.0;
    // // for (int i = 0; i < nnuc; i++)
    // //     for (int j = 0; j < nnuc; j++)
    // //         if (i != j) mu_tot1 += (double) freq_boundary_states[i]*freq_boundary_states[j]*mutCoeff(i,j);
    // for (int i = 0; i < nnuc; i++)
    //     for (int j = 0; j < num_states; j++)
    //         // if (i != j) mu_tot += (double) freq_boundary_states[i]*rate_matrix[i*num_states + j];
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
