#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modeldna.h"

/* PoMo mutation rate boundaries.  This is not necessary anymore
   because rates are bounded by underlying Markov model*/
/* const double POMO_MIN_RATE =  5e-5; */
/* const double POMO_INIT_RATE = 1e-3; */
/* const double POMO_MAX_RATE =  1e-1; */

// HETEROZYGOSITY BOUNDARIES ARE NOW SET IN A DYNAMIC WAY TO IMPROVE NUMERICAL
// STABILITY. EXCEEDING THESE VALUES ONLY INDUCES WARNINGS.
/* Boundaries for level of polymorphism or heterozygosity. The heterozygosity
   boundaries strongly affect numerical stability. I set them so that errors are
   very seldom but there may be data that requires different boundaries. Maybe
   they should be set variable, depending on data or user input. */
const double POMO_MIN_HETEROZYGOSITY =  1e-5;
const double POMO_MAX_HETEROZYGOSITY =  1e-1;
/* Not so stringent values. However, crashes are expected, especially because of
   the maximum value (?).*/
/* const double POMO_MIN_HETEROZYGOSITY =  1e-4; */
/* const double POMO_MAX_HETEROZYGOSITY =  1e-1; */

// Relative boundaries for the heterozygosity.
const double POMO_MIN_REL_HETEROZYGOSITY = 0.5;
const double POMO_MAX_REL_HETEROZYGOSITY = 3;

/* Boundaries for boundary frequencies.  This is not necessary anymore
   because those are set by the underlying Markov model.  The actual
   boundaries will be set, e.g., to
   freq_boundary_states[i]*POMO_MIN_REL_FREQ. */
/* const double POMO_MIN_REL_FREQ = 0.5; */
/* const double POMO_MAX_REL_FREQ = 2.0; */

// Boundaries for each of the boundary states.
const double POMO_MIN_BOUNDARY_FREQ = 0.05;
const double POMO_MAX_BOUNDARY_FREQ = 0.95;

class ModelPoMo : virtual public ModelMarkov
{
 public:
    /**
     * Constructor.
     * ModelMarkov() constructor calls ModelSubst() constructor.
     * ModelSubst():
     * - allocates state_freq[tree->aln->num_states]
     * ModelMarkov():
     * - allocates rates[getNumRateEntries()] = rates[n*(n-1)/2];
     * - allocates eigenvalues and eigenvectors.
     *
     * @param model_name The name of the model (e.g., "HKY+P").
     * @param model_params The parameters of the model (user defined models.).
     * @param freq_type
     * @param freq_params
     * @param tree Associated tree for the model.
     * @param pomo_heterozygosity heterozygosity of PoMo
     *
     * @return
     */
    ModelPoMo(const char *model_name, string model_params, StateFreqType freq_type, string freq_params,
              PhyloTree *tree, string pomo_heterozygosity);

    ModelPoMo(PhyloTree *tree);

    ~ModelPoMo();


    // Tell the compiler we want both the init functions (resolve
    // warning from clang about PoMo hiding the overloaded init
    // function).
    using ModelMarkov::init;
    /**
     * Initialize the PoMo model. Run by constructor.
     *
     * @param model_name
     * @param model_params
     * @param freq_type
     * @param freq_params
     */
    virtual void init(const char *model_name,
                      string model_params,
                      StateFreqType freq_type,
                      string freq_params,
                      string pomo_heterozygosity);

    /**
     *  \brief Initialize underlying mutation model.
     *
     * PoMo is built on top of any Markov mutation model which is of class
     * ModelMarkov in IQ-TREE and is denominated "underlying" Markov model or
     * "underlying" mutation model. The idea is to use the machinery of the
     * underlying Markov model and only introduce an additional layer that adds
     * polymorphic states and one parameter, namely the heterozygosity or level
     * of polymorphism.
     *
     */
    void init_mutation_model(const char *model_name,
                             string model_params,
                             StateFreqType freq_type,
                             string freq_params);

    /**
     *  \brief Initialize sampling type.
     *
     *  Weighted (standard) or sampled.
     *
     */
    void init_sampling_method();

    /**
     *  \brief Initialize state frequencies.
     *
     *  Use the machinery of the underlying mutation model.
     *
     */
    void init_boundary_frequencies();

    /**
     *  \brief Check if parameters are optimized or not.
     *
     *  The parameters of PoMo are the ones from the underlying mutation model
     *  plus the the heterozygosity or the level of polymorphism.
     *
     */
    void init_fixed_parameters(string model_params,
                               string pomo_heterozygosity);

    /**
     * Initialize rate_matrix and state_freq for boundary mutation model.
     */
    void updatePoMoStatesAndRateMatrix();

	/**
	 * @return model name
	 */
	virtual string getName();

    /**
     *  @return Number of free parameters.
     */
    virtual int getNDim();

	/**
		@return the number of dimensions corresponding to state frequencies
	*/
	virtual int getNDimFreq();


    /**
     * Set bounds for joint optimization with BFGS.
     */
    virtual void setBounds(double *lower_bound,
                           double *upper_bound,
                           bool *bound_check);

    /**
     *  Write information to output stream (only with -vv).
     *  @param out Output stream.
     */
    virtual void writeInfo(ostream &out);

    /**
     *  The target function which needs to be optimized
     *  @param x the input vector x
     *  @return the function value at x
    */
    virtual double targetFunk(double x[]);

    /**
     *  @return TRUE if parameters are at the boundary that may cause
     *  numerical unstability
     */
    virtual bool isUnstableParameters();

    virtual bool isPolymorphismAware() { return true; };

    /**
     *  @return the number of rate entries
     */
    virtual int getNumRateEntries() { return n_alleles*(n_alleles-1)/2; };

    /**
     *  \brief Normalize boundary frequencies so that they sum to 1.0.
     *
     */
    void normalize_boundary_freqs(double * bfs);

    /**
     *  \brief Check if boundary frequencies are within bounds.
     *
     */
    void check_boundary_freqs(double * bfs);

    /**
     * Estimate the empirical (relative) boundary state frequencies.  The
     * number of A, C, G, and T in the data is summed up and the
     * relative proportion of all bases is calculated.  The empirical
     * boundary state frequencies are set to these relative proportions.
     *
     * @param freq_boundary_states (OUT) The estimated boundary frequencies,
     * size num_states.
     */
    void estimateEmpiricalBoundaryStateFreqs(double * freq_boundary_states);

    /**
     * Estimate the empirical (relative) sum of polymorhic states.
     * The number of polymorphic sites in the data is summed up and
     * the relative proportion of all sites is calculated and returned.
     *
     */
    double estimateEmpiricalWattersonTheta();

  // Extract the rate entries of a rate matrix (basically remove the diagonal
  // and align upper triangle before lower triangle). IN: rate matrix m; OUT:
  // rates r (array of size n_connections or 2*n_connections).
  void rate_matrix_to_rates(double *m, double *r);

  // Extract the exchangeability entries of a rate matrix (basically remove the
  // diagonal and align upper triangle before lower triangle; divide through
  // stationary frequencies of target alleles). IN: rate matrix m; OUT: rates r
  // (array of size n_connections or 2*n_connections).
  void rate_matrix_to_exchangeabilities(double *m, double *r);

    /**
     * Report the model rates to the output file stream 'out'.
     *
     * @param out Output file stream.
     */
    void report_model_params(ostream &out, bool reset_scale = true);

    /**
     * Report the state frequencies to the output file stream 'out'.
     *
     * @param out Output file stream.
     */
    virtual void report(ostream &out);

    /**
     * Normalize the mutation probabilities such that the given heterozygosity
     * or level of polymorphism is honored.
     *
     */
    void normalizeMutationRates();

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        Restore object from the checkpoint.
    */
    virtual void restoreCheckpoint();

	/**
       Decompose the rate matrix into eigenvalues and eigenvectors.
       This function is necessary, because handling of reversible and
       non-reversible models differs in ModelMarkov, but should not
       differ for ModelPoMo.
	*/
	virtual void decomposeRateMatrix();

  // I had serious problems with numerical instabilities because of fixed
  // boundaries for the heterozygosity. However, it is easy to get an empirical
  // value of the heterozygosity and set boundaries according to this value.
  // This improves numerical stability. This step does not make sense and,
  // hence, is not done, if the heterozygosity is set by the user or to the
  // empirical value.
  void set_heterozygosity_boundaries();

	/**
     compute the transition probability matrix.
     @param time time between two events
     @param mixture (optional) class for mixture model
     @param selected_row (optional) only compute the entries of one selected row. By default, compute all rows
     @param trans_matrix (OUT) the transition matrix between all pairs of states.
     Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, int mixture = 0, int selected_row = -1);

    /**
     *  Set the scale factor of the mutation rates to NEW_SCALE.
     *
     *  @param scale (IN).
     */
  void setScale(double new_scale);

  /**
   * get the underlying mutation model, used with PoMo model
   */
  virtual ModelSubst *getMutationModel() { return mutation_model; }

    /** compute the tip likelihood vector of a state for Felsenstein's pruning algorithm
     @param state character state
     @param[out] state_lk state likehood vector of size num_states
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);

 protected:

    ModelMarkov *mutation_model;

    /**
        Compute the rate matrix and then normalize it such that the total number of mutations is 1.
        @param rate_matrix (OUT).  It is filled with rate matrix entries
        @param state_freq state frequencies
        @param num_state number of states
    */
    virtual void computeRateMatrix(double **rate_matrix, double *state_freq, int num_state);

  // Get the current scale factor of the mutation rates.
  double getScale();
    /**
     * This function is served for the multi-dimension
     * optimization. It should pack the model parameters into a vector
     * that is index from 1 (NOTE: not from 0)
     *
     * @param variables (OUT) Vector of variables, indexed from 1.
     */
    virtual void setVariables(double *variables);

    /**
     * This function is served for the multi-dimension
     * optimization. It should assign the model parameters from a
     * vector of variables that is index from 1 (NOTE: not from 0)
     *
     * @param variables Vector of variables, indexed from 1.
     * @return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
     */
    virtual bool getVariables(double *variables);

    /**
	 * Called from getVariables() to update the rate matrix for the
	 * new model parameters.  For ModelPoMo this is only a dummy
	 * function, which has to be declared empty because otherwise,
	 * restoreCheckpoint() calls the ModelMarkov::setRates() which in
	 * turn throws an error.
	 */
	virtual void setRates();

 protected:
    /*!<  Virtual population size of the PoMo model. */
    int N;

    /**
     * Full mutation rate matrix (Q^REV + Q^NONREV; or, Q^GTR + Q^FLUX).
     */
    double *mutation_rate_matrix;

    /**
     * 4 unnormalized stationary frequencies of boundary states.
     */
    double *freq_boundary_states;

    /**
     * Normalized empirical stationary frequencies.
     */
    double *freq_boundary_states_emp;

    /* /\** */
    /*  * The rate matrix of the PoMo model. */
    /*  *\/ */
    /* double *rate_matrix; */

    /**
     * Decompose state (0..57) into abundance of two nucleotides.
     *
     * @param i (OUT) Abundance of nucleotide 1.
     * @param nt1 (OUT) Nucleotide 1 (0: A, 1: C, 2: G, 3: T).
     * @param nt2 (OUT) Nucleotide 2 (0: A, 1: C, 2: G, 3: T).
     */
    void decomposeState(int state, int &i, int &nt1, int &nt2);

    /**
     * Compute the normalized stationary frequencies that fulfill the
     * detailed balance condition.
     */
    void computeStateFreq();

    /**
     * Compute probability of change from state1 to state2 in one
     * Moran with boundary mutation model generation.
     */
    double computeProbBoundaryMutation(int state1, int state2);

    bool isBoundary(int state);

    bool isPolymorphic(int state);

    /**
     * Get symmetric mutation coefficient from nt1 to nt2 (mij = mji).
     *
     * @param nt1
     * @param nt2
     *
     * @return
     */
    double mutCoeff(int nt1, int nt2);

    /**
     * Set the initial mutation coefficients which resemble the level
     * of polymorphism in the data.
     *
     * @return
     */
    void setInitialMutCoeff();

    /**
     * Compute the sum of the frequencies of the boundary states.
     *
     * @return
     */
    double computeSumFreqBoundaryStates();

    /**
     * Compute part of normalization constant of polymorphic states.
     *
     * @return
     */
    double computeSumFreqPolyStates();

    /**
     * Compute the sum over lamda_pol without mutliplying with mutation
     * coefficients. This is useful if the mutation coefficient is constant,
     * e.g., when the initial value is set.
     *
     * @return
     */
    double computeSumFreqPolyStatesNoMut();

    /**
     * Compute the normalization constant.  This constant ensures that
     * the stationary frequencies (invariant measure) sum up to 1
     * (invariant distribution).
     *
     * @return
     */
    double computeNormConst();

    /**
     * Set the boundary frequency of T such that all boundary frequencies
     * sum up to one.  This is done, so that they can be compared to
     * the frequencies of the GTR model.
     */
    void updateFreqBoundaryState ();

    /**
     * Precision and treshold value for mathematical computations and
     * numerical stability checks.
     *
     * Feel free to set to a global value if there is one available.
     *
     */
    double eps;

    /// Number of nucleotides (alleles).  This might be useful in the
    /// future, when we do not restrict PoMo to DNA models only.
    /// Eventual todo: do not hardcode this.
    int n_alleles;

    /// True if heterozygosity has been fixed.
    bool fixed_heterozygosity;

    /// True if heterozygosity has been fixed by user.
    bool fixed_heterozygosity_usr;

    /// True if heterozygosity has been fixed to empirical estimate
    /// from data.
    bool fixed_heterozygosity_emp;

    /**
     * Heterozygosity or level of polymorphism. Will be set by init() and is
     * used to normalize the mutation coefficients. This is needed to be done
     * because PoMo had trouble to estimate this (especially, when N was large).
     * See also ModelPoMo::normalizeMutationRates().
     */
    double heterozygosity;

    /// True if the model parameters are fixed (e.g., if the
    /// transition to transversion ratio is set in the HKY model).
    bool fixed_model_params;

    /// This array contains the ratio of the fixed mutation
    /// probability to the base mutation rate that is still estimated
    /// and describes the amount of polymorphism.  I.e., for the HKY
    /// model with transition to transversion ratio 2.0, this vector
    /// would look like [1, 2, 1, 1, 2, 1]
    double * fixed_model_params_ratio;

    /// The number of connections between nucleotides.  If n_alleles=4,
    /// there are 6 connections.  Set in ModelPoMo::init().
    int n_connections;

    /// Random binomial sampling or weighted; specified when alignment
    /// is created already (Alignment::readCountsFormat()).
    SamplingType sampling_method;

  // Minimum heterozygosity, set by `set_heterozygosity_boundaries()`.
  double min_heterozygosity;

  // Maximum heterozygosity, set by `set_heterozygosity_boundaries()`.
  double max_heterozygosity;

  // The scale factor of the mutation rates. Important for Gamma rate
  // heterogeneity.
  double scale;
};
#endif /* _MODELPOMO_H_ */
