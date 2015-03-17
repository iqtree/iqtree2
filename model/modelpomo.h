#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modelgtr.h"

class ModelPoMo : public ModelGTR
{
 public:
    /**
     * Constructor
     *
     * @param tree Associated tree for the model.
     * @param count_rate
     */
    ModelPoMo(PhyloTree *tree, bool count_rate = true);

    ~ModelPoMo();

    /**
     * Initialize the PoMo model. Run by constructor.
     *
     * @param model_name
     * @param model_params
     * @param freq
     * @param freq_params
     */
    virtual void init(const char *model_name,
                      string model_params,
                      StateFreqType freq,
                      string freq_params);

    /* /\** */
    /*  *  Deprecated!  Unreversible.  Initialize rate_matrix and */
    /*  *  state_freq. */
    /*  *\/ */
    /* void initMoranWithMutation(); */

    /**
     * Initialize rate_matrix and state_freq for boundary mutation model.
     */
    void initMoranWithBoundaryMutation();

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

 protected:
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
    void setFreqFixedStates ();

};

#endif /* _MODELPOMO_H_ */
