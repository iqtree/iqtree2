#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modelgtr.h"

class ModelPoMo : public ModelGTR
{     
  
public:  
  /** 
   * Constructor
   * 
   * @param tree associated tree for the model
   * 
   */
  ModelPoMo(PhyloTree *tree, bool count_rate = true);
  
  ~ModelPoMo();
  /* TODO ModelSubst::rates */

  virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

  /**
   * initialize rate_matrix and state_freq
   */
  void initMoranWithMutation();

  /**
   * initialize rate_matrix and state_freq for boundary mutation model
   */
  void initMoranWithBoundaryMutation();

	/**
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables);

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
	*/
	virtual void getVariables(double *variables);


private:

  /** mutation probabilities, either 4x4 or 6 entries for non-reversible or reversible model, respectively */
  double *mutation_prob;

  /** 4 unnormalized stationary frequencies of fixed states (reversible model only) */
  double *freq_fixed_states;

  double *rate_matrix;

  /**
   * P(i,major,minor) is the probability to increase the number of
   * major alleles from i to i+1.
   * 
   * @param i abundance of major allele
   * @param major major allele (0: A, 1: C, 2: G, 3: T)
   * @param minor minor allele (0: A, 1: C, 2: G, 3: T)
   */
  double computeP(int i, int major, int minor);

  /**
   * R(i,major,minor) is the probability of no frequency change in the
   * Moran model with mutation at one locus with two alleles.
   * 
   * @param i abundance of major allele
   * @param major major allele (0: A, 1: C, 2: G, 3: T)
   * @param minor minor allele (0: A, 1: C, 2: G, 3: T)
   */
  double computeR(int i, int major, int minor);

  /**
   * Decompose state (0..57) into abundance of two nucleotides.
   * 
   * @param i (OUT) abundance of nucleotide 1
   * @param nt1 (OUT) nucleotide 1 (0: A, 1: C, 2: G, 3: T)
   * @param nt2 (OUT) nucleotide 2 (0: A, 1: C, 2: G, 3: T)
   */
  void decomposeState(int state, int &i, int &nt1, int &nt2);

  /**
   * Compute probability of change from state1 to state2 in one Moran
   * model generation.
   */
  double computeProb(int state1, int state2);

  void computeStateFreq();

  /**
   * Compute probability of change from state1 to state2 in one Moran with boundary mutation
   * model generation.
   */
  double computeProbBoundaryMutation(int state1, int state2);

  bool isFixed(int state);

  bool isPolymorphic(int state);

  double mutCoeff(int nt1, int nt2);

  double computeNormConst();

};

#endif /* _MODELPOMO_H_ */
