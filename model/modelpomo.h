#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modeldna.h"

const double POMO_MIN_RATE =  5e-5;
const double POMO_INIT_RATE = 1e-3;
const double POMO_MAX_RATE =  1e-2;
/* The actual boundaries will be set, e.g., to
   #freq_fixed_states[i]*POMO_MIN_REL_FREQ. */
const double POMO_MIN_REL_FREQ = 0.5;
const double POMO_MAX_REL_FREQ = 2.0;

class ModelPoMo : virtual public ModelGTR
{
 public:
    /**
     * Constructor.
     * ModelGTR() constructor calls ModelSubst() constructor.
     * ModelSubst():
     * - allocates state_freq[tree->aln->num_states]
     * ModelGtr():
     * - allocates rates[getNumRateEntries()] = rates[n*(n-1)/2];
     *   cf. modelsubst.h
     * - allocates eigenvalues and eigenvectors.
     *
     * @param model_name The name of the model (e.g., "HKY+P").
     * @param model_params The parameters of the model (user defined models.).
     * @param freq_type
     * @param freq_params
     * @param tree Associated tree for the model.
     * @param is_reversible True if the reversible model should be used.
     * @param pomo_params parameters for PoMo
     *
     * @return
     */
    ModelPoMo(const char *model_name, string model_params, StateFreqType freq_type, string freq_params, 
        PhyloTree *tree, bool is_reversible, string pomo_params);

    ModelPoMo(PhyloTree *tree);

    ~ModelPoMo();


    // Tell the compiler we want both the init functions (resolve
    // warning from clang about PoMo hiding the overloaded init
    // function).
    using ModelGTR::init;
    /**
     * Initialize the PoMo model. Run by constructor.
     *
     * @param model_name
     * @param model_params
     * @param freq_type
     * @param freq_params
     * @param is_reversible True if the reversible model should be used.
     */
    virtual void init(const char *model_name,
                      string model_params,
                      StateFreqType freq_type,
                      string freq_params,
                      bool is_reversible,
                      string pomo_params);

    /* /\** */
    /*  *  Deprecated!  Unreversible.  Initialize rate_matrix and */
    /*  *  state_freq. */
    /*  *\/ */
    /* void initMoranWithMutation(); */

    /**
     * Initialize rate_matrix and state_freq for boundary mutation model.
     */
    void updatePoMoStatesAndRates();

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
     * Write information to output stream (only with -vv).
     * @param out Output stream.
     */
    virtual void writeInfo(ostream &out);

    /**
     *  the target function which needs to be optimized
     *  @param x the input vector x
     *  @return the function value at x
    */
    virtual double targetFunk(double x[]);

    /**
     * @return TRUE if parameters are at the boundary that may cause
     * numerical unstability
     */
    virtual bool isUnstableParameters();

    virtual bool isPolymorphismAware() { return true; };

    /* /\**  */
    /*  * Set the substitution rate parameters by a specification.  From */
    /*  * ModelDNA::setRateType(). */
    /*  * */
    /*  * Sets the array #mutation_prob and the vector #param_fixed. */
    /*  *  */
    /*  * @param rate_spec a string of six letters describing how rates are related */
    /*  * @return TRUE if successful, FALSE otherwise */
    /*  *\/ */
    /* bool setRateType(const char *rate_spec); */

    /**
     *  @return the number of rate entries
     */
    virtual int getNumRateEntries() { return nnuc*(nnuc-1)/2; };

    /* /\**  */
    /*  * Read state frequencies from an input stream.  If it fails, */
    /*  * throw error message. */
    /*  *  */
    /*  * @param in input stream */
    /*  *\/ */
    /* void readFixedStateFreq(istream &in); */

    /* /\**  */
    /*  * Read state frequencies from comma-separated string.  Might */
    /*  * throw error message. */
    /*  *  */
    /*  * @param str  */
    /*  *\/ */
    /* void readFixedStateFreq(string str); */

    /* /\**  */
    /*  * Read model parameters from a file.  The file needs to contain */
    /*  * the upper-triangle rate matrix and the state frequencies. */
    /*  *  */
    /*  * @param file_name  */
    /*  *\/ */
    /* void readMutationParameters(const char *file_name); */

    /* /\**  */
    /*  * Read the upper-triangle rate matrix from an input stream. */
    /*  * Throw error message if failing. */
    /*  *  */
    /*  * @param in input stream */
    /*  *\/ */
    /* void readMutationRates(istream &in); */

    /* /\**  */
    /*  * Read rate parameters from a comma-separated string.  Throw */
    /*  * error message if failing. */
    /*  *  */
    /*  * @param str input string */
    /*  *\/ */
    /* void readMutationRates(string str); */

    /**
     * Estimate the empirical (relative) fixed state frequencies.  The
     * number of A, C, G, and T in the data is summed up and the
     * relative proportion of all bases is calculated.  The empirical
     * fixed state frequencies are set to these relative proportions.
     *
     * @param freq_fixed_states (OUT) The estimated fixed frequencies,
     * size num_states.
     */
    void estimateEmpiricalFixedStateFreqs(double * freq_fixed_states);

    /**
     * Estimate the empirical (relative) sum of polymorhic states.
     * The number of polymorphic sites in the data is summed up and
     * the relative proportion of all sites is calculated and returned.
     *
     */
    double estimateEmpiricalWattersonTheta();

    /**
     * Report the model rates to the output file stream 'out'.
     *
     * @param out Output file stream.
     */
    void report_rates(ostream &out);

    /**
     * Report the state frequencies to the output file stream 'out'.
     *
     * @param out Output file stream.
     */
    virtual void report(ostream &out);

    /** 
     * Normalize the mutation probabilities such that the given level of
     * polymorphism is honored (level_of_polymorphism).
     * 
     */
    void normalizeMutationProbs();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

 protected:

    ModelDNA *dna_model;

    /**
        compute the rate matrix and then normalize it such that the total number of substitutions is 1.
        @param rate_matrix (OUT).  It is filled with rate matrix entries
        @param state_freq state frequencies
        @param num_state number of states
    */
    virtual void computeRateMatrix(double **rate_matrix, double *state_freq, int num_state);

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


 private:

    /// Virtual population size of the PoMo model.
    int N;

    /**
     * Mutation probabilities, 6 entries for reversible model.
     */
    double *mutation_prob;

    /**
     * 4 unnormalized stationary frequencies of fixed states.
     */
    double *freq_fixed_states;

    /**
     * Normalized empirical stationary frequencies.
     */
    double *freq_fixed_states_emp;

    /**
     * The rate matrix of the PoMo model.
     */
    double *rate_matrix;

    /* /\** */
    /*  * Deprecated!  Unreversible. */
    /*  * */
    /*  * P(i,major,minor) is the probability to increase the number of */
    /*  * major alleles from i to i+1. */
    /*  * */
    /*  * @param i abundance of major allele */
    /*  * @param major major allele (0: A, 1: C, 2: G, 3: T) */
    /*  * @param minor minor allele (0: A, 1: C, 2: G, 3: T) */
    /*  *\/ */
    /* double computeP(int i, int major, int minor); */

    /* /\** */
    /*  * Deprecated!  Unreversible. */
    /*  * */
    /*  * R(i,major,minor) is the probability of no frequency change in the */
    /*  * Moran model with mutation at one locus with two alleles. */
    /*  * */
    /*  * @param i abundance of major allele */
    /*  * @param major major allele (0: A, 1: C, 2: G, 3: T) */
    /*  * @param minor minor allele (0: A, 1: C, 2: G, 3: T) */
    /*  *\/ */
    /* double computeR(int i, int major, int minor); */

    /**
     * Decompose state (0..57) into abundance of two nucleotides.
     *
     * @param i (OUT) Abundance of nucleotide 1.
     * @param nt1 (OUT) Nucleotide 1 (0: A, 1: C, 2: G, 3: T).
     * @param nt2 (OUT) Nucleotide 2 (0: A, 1: C, 2: G, 3: T).
     */
    void decomposeState(int state, int &i, int &nt1, int &nt2);

    /* /\** */
    /*  * Deprecated!  Unreversible. */
    /*  * */
    /*  * Compute probability of change from state1 to state2 in one Moran */
    /*  * model generation. */
    /*  *\/ */
    /* double computeProb(int state1, int state2); */

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

    bool isFixed(int state);

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
     * Compute the sum of the frequencies of the fixed states.
     *
     * @return
     */
    double computeSumFreqFixedStates();

    /**
     * Compute part of normalization constant of polymorphic states.
     *
     * @return
     */
    double computeSumFreqPolyStates();
    
    /**
     * Computes the sum over lamda_pol without mutliplying with
     * mutation coefficients.  This is useful if the mutation
     * coefficient is constant, e.g., when the initial value is set.
     * Compute part of normalization constant of polymorphic states.
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
     * Set the fixed frequency of T such that all fixed frequencies
     * sum up to one.  This is done, so that they can be compared to
     * the frequencies of the GTR model.
     */
    void updateFreqFixedState ();

    /**
     * Precision and treshold value for mathematical computations and
     * numerical stability checks.
     *
     * Feel free to set to a global value if there is one available.
     *
     */
    double eps;

    /* /// Rate parameter specification, a string of 6 characters.  E.g., */
    /* /// for the HKY model, it will be set to '010010' by */
    /* /// SetRateType(). */
    /* string param_spec; */

    /* /// Vector of boolean, TRUE if correspodning parameter is fixed */
    /* /// and FALSE otherwise.  Set by SetRateType(). */
    /* vector<bool> param_fixed; */

    /// Number of nucleotides (alleles).  This might be useful in the
    /// future, when we do not restrict PoMo to DNA models only.
    /// Eventual todo: do not hardcode this.
    int nnuc;

    /// True if the level of polymorphism has been fixed (either to a
    /// user given value or to the empirically observed one).
    bool fixed_level_of_polymorphism;

    /**
     * level of polymorphism.  Will be set by init() and is
     * used to normalize the mutation coefficients.  This is needed to
     * be done because PoMo had trouble to estimate this (especially,
     * when N was large).
     */
    double level_of_polymorphism;

    /// True if the model parameters are fixed (e.g., if the
    /// transition to transversion ratio is set in the HKY model).
    bool fixed_model_params;

    /// This array contains the ratio of the fixed mutation
    /// probability to the base mutation rate that is still estimated
    /// and describes the amount of polymorphism.  I.e., for the HKY
    /// model with transition to transversion ratio 2.0, this vector
    /// would look like [1, 2, 1, 1, 2, 1]
    double * fixed_model_params_ratio;

    /// The number of connections between the alleles (e.g., 6 for 4
    /// nucleotides).
    int n_connections;

    /// Random binomial sampling or weighted; specified when alignment
    /// is created already (Alignment::readCountsFormat()).
    SamplingType sampling_type;
};

#endif /* _MODELPOMO_H_ */
