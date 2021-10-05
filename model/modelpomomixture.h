//
//  modelpomomixture.h
//  iqtree
//  Mixture PoMo models to include e.g. Gamma-rate heterogeneity
//
//  Created by Minh Bui on 7/22/16.
//
//

#ifndef modelpomomixture_h
#define modelpomomixture_h

#include <stdio.h>
#include "modelpomo.h"
#include "modelmixture.h"
#include "rateheterogeneity.h"

// This would be the preferred lower bound but eigendecomposition leads to
// numerical errors.
// const double POMO_GAMMA_MIN = 0.02;
const double POMO_GAMMA_MIN = 0.05;
const double POMO_GAMMA_MAX = 100;

enum PomoMixtureOptMode {OPT_NONE, OPT_RATEHET, OPT_POMO};

/**
    Mixture PoMo models
*/
class ModelPoMoMixture : public ModelPoMo {

public:
	typedef ModelPoMo super;
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelPoMoMixture(const char *model_name,
                     const string& model_params,
                     StateFreqType freq_type,
                     const string& freq_params,
                     PhyloTree *tree,
                     string pomo_params,
                     string pomo_rate_str,
                     PhyloTree* report_to_tree);

    virtual ~ModelPoMoMixture();

	/**
	 * @return model name
	 */
	virtual string getName() const override;

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint) override;

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint() override;

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint() override;

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint() override;

	/**
		@return the number of dimensions
	*/
	virtual int getNDim() const override;

	/**
		@return the number of dimensions corresponding to state frequencies
	*/
	virtual int getNDimFreq() const override;


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) override;

    /**
     * @return TRUE if parameters are at the boundary that may cause
     * numerical unstability
     */
    virtual bool isUnstableParameters() const override;

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, 
	                       bool *bound_check) override;

	/**
		optimize model parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon,
                                      PhyloTree* report_to_tree) override;

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) override;

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix() override;

    /**
     * Report the state frequencies to the output file stream 'out'.
     *
     * @param out Output file stream.
     */
    virtual void report(ostream &out) override;


    /**
        rate heterogeneity among sites, TODO: redesign this as separate
    */
    RateHeterogeneity *ratehet;

  // Mon Jul 3 14:47:08 BST 2017; added by Dominik. I had problems with mixture
  // models together with PoMo and rate heterogeneity. E.g., a model
  // "MIX{HKY+P+N9+G2,GTR+P+N9+G2}" leads to segmentation faults because the
  // `ModelPoMoMixture` reports a /wrong/ number of states (i.e., it reports 52
  // instead of 104). Consequently, the `initMem()` function of ModelMixture,
  // messes up the `eigenvalues`, etc., variables of the `ModelPoMoMixture`s. I
  // circumvent this, by adding this virtual function; for normal models, it
  // just returns `num_states`, however, for mixture models, it returns
  // `num_states*nmixtures`.
  virtual int get_num_states_total() override;

  // Mon Jul 3 15:53:00 BST 2017; added by Dominik. Same problem as with
  // `get_num_states_total()`. The pointers to the eigenvalues and eigenvectors
  // need to be updated recursively, if the model is a mixture model. For a
  // normal Markov model, only the standard pointers are set. This was done in
  // `ModelMixture::initMem()` before.
  virtual void update_eigen_pointers(double *eval, double *evec, double *inv_evec, 
									 double *inv_evec_transposed) override;


	/**
     compute the transition probability matrix.
     @param time time between two events
     @param mixture (optional) class for mixture model
     @param trans_matrix (OUT) the transition matrix between all pairs of states.
     Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, 
	                                int mixture = 0) const override;


protected:

	ModelMixture mixture;

    /** normally false, set to true while optimizing rate heterogeneity */
    PomoMixtureOptMode opt_mode;

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables) override;

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(const double *variables) override;

#ifdef _MSC_VER
	//MSVC generates warning messages about these member functions
	//being inherited "via dominance". Explictly declaring them
	//instead shuts those warnings up.
public:
	virtual string getNameParams() const override { return mixture.getNameParams(); }
	virtual bool isMixture()       const override { return mixture.isMixture(); }
	virtual int getNMixtures()     const override { return mixture.getNMixtures(); }
	virtual double getMixtureWeight(int cat)        const override { return mixture.getMixtureWeight(cat); }
	virtual void setMixtureWeight(int cat, double weight) override { mixture.setMixtureWeight(cat, weight); }
	virtual void setFixMixtureWeight(bool fix_weight)     override { mixture.setFixMixtureWeight(fix_weight); }
	virtual ModelSubst* getMixtureClass(int cat)    const override { return mixture.getMixtureClass(cat); }
	virtual void setMixtureClass(int cat, ModelSubst* m)  override { mixture.setMixtureClass(cat, m); }
	virtual void getStateFrequency(double* state_freq,
	                               int mixture_num = 0) const override 
	{ 
		mixture.getStateFrequency(state_freq, mixture_num);
	}
	virtual void computeTransDerv(double time, double* trans_matrix,
		double* trans_derv1, double* trans_derv2, int mixture_num = 0) const {
		mixture.computeTransDerv(time, trans_matrix,
			trans_derv1, trans_derv2, mixture_num) override;
	}
	virtual void setOptimizeSteps(int steps)   override { mixture.setOptimizeSteps(steps); }
	virtual uint64_t getMemoryRequired() const override { return mixture.getMemoryRequired(); }
	virtual void writeParameters(ostream& out) override { mixture.writeParameters(out); }

    virtual bool isPolymorphismAware()   const override { return super::isPolymorphismAware(); };
    virtual int getNumRateEntries()      const override { return super::getNumRateEntries();  }
    virtual void computeTipLikelihood
					(PML::StateType state, double* state_lk) 
					const override {
		super::computeTipLikelihood(state, state_lk);
    }
    virtual ModelSubst* getMutationModel() override { return super::getMutationModel(); }
    virtual void setRates()                override { ModelPoMo::setRates(); }
    virtual void computeRateMatrix(double** rate_matrix, 
	                               double* state_freq, 
								   int num_state) const override{
		super::computeRateMatrix(rate_matrix, state_freq, num_state);
    }
#endif

};

#endif /* modelpomomixture_h */
