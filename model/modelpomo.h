#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modeldna.h"

/* TODO: Make it work with N arbitrary, this should not be too hard */
/*       (phylotreesse.cpp has to be changed to). */
/* TODO: Implement different DNA subsitution models. */

const double POMO_MIN_RATE = 1e-6;
const double POMO_MAX_RATE = 1e-3;
const double POMO_MIN_REL_FREQ = 0.5;
const double POMO_MAX_REL_FREQ = 2.0;

class ModelPoMo : public ModelGTR
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
     * @todo Implement `model_params'.
     * @todo Implement `freq_params'.
     * 
     * @param model_name The name of the model (e.g., "HKY+P").
     * @param model_params The parameters of the model (user defined models.).
     * @param freq_type
     * @param freq_params
     * @param tree Associated tree for the model.
     * @param is_reversible True if the reversible model should be used.
     *
     * @return
     */
    ModelPoMo(const char *model_name, string model_params, StateFreqType freq_type, string freq_params, PhyloTree *tree, bool is_reversible);

    ~ModelPoMo();

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
                      bool is_reversible);

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
     *	@return the function value at x
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
     * @param abs_state_freq An array of the absolute state frequencies, size N.
     * @param freq_fixed_states (OUT) The estimated fixed frequencies,
     * size num_states.
     */
    void estimateEmpiricalFixedStateFreqs(unsigned int * abs_state_freq, double * freq_fixed_states);
        
 protected:

    ModelDNA *dna_model;

 	/**
 		compute the rate matrix and then normalize it such that the total number of substitutions is 1.
 		@param rate_matrix (IN/OUT) As input, it contains rate parameters. On output it is filled with rate matrix entries
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
     */
	virtual void getVariables(double *variables);


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
};

#endif /* _MODELPOMO_H_ */
