#ifndef _MODELPOMO_H_
#define _MODELPOMO_H_

#include "modelnonrev.h"

class ModelPoMo : public ModelNonRev
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

  void initMoranWithMutation();

  virtual int getNDim();

	/**
		optimize model parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double epsilon);


private:

  /** 4*4 mutation probabilities */
  double *mutation_prob;

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

};

#endif /* _MODELPOMO_H_ */
