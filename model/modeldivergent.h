
#pragma once
#ifndef model_divergent_h
#define model_divergent_h

#include "modelsubst.h"
#include "modelinfofromyamlfile.h" //for ModelVariable class

class ModelDivergent: public ModelSubst {

    //Consistency restrictions: Ascertainment biases have
    //to agree, all must be reversible (or all not).

protected:
    std::vector<ModelSubst*>   subtree_models;
    std::vector<ModelVariable> own_parameters;
    PhyloTree* phylo_tree; //Set by setTree(), and used in
                           //decomposeRateMatrix() member function.
    int optimize_steps;

public:
    typedef ModelSubst super;
    explicit ModelDivergent(int nstates);
    virtual ~ModelDivergent();
	virtual void setTree(PhyloTree *tree) override;      //for each subtree model
	virtual void setNumberOfStates(int states) override; //for each subtree model
	virtual bool getSpecifiedAscertainmentBiasCorrection(ASCType& asc_type) override;
	virtual int  getNDim()           const override;
	virtual int  getNDimFreq()       const override;
    virtual bool isReversible()      const override;
    virtual int  getNumRateEntries() const override;
	virtual double computeTrans(double time, int model_id, 
                                int state1, int state2) const override;
	virtual double computeTrans(double time, int model_id, 
                                int state1, int state2, 
                                double &derv1, double &derv2) const override;
	virtual void   getStateFrequency(double *state_freq, 
	                                 int mixture = 0) const override;
    virtual void   setStateFrequency(double *state_freq) override;
    virtual void   decomposeRateMatrix() override;
    virtual void   setOptimizeSteps(int optimize_steps) override;
	virtual double optimizeParameters(double gradient_epsilon,
                                      PhyloTree* report_to_tree) override;
	virtual void setBounds(double *lower_bound, double *upper_bound, 
                           bool *bound_check) override;
	virtual void setVariables(double *variables) override;
	virtual bool getVariables(const double *variables) override;
	virtual void afterVariablesChanged() override;
	virtual bool isUnstableParameters()  const override;
	virtual void writeInfo(ostream &out) override;
	virtual void report(ostream &out)    override;
    virtual uint64_t getMemoryRequired() const override;
    virtual void startCheckpoint()       override;
    virtual void saveCheckpoint()        override;
    virtual void restoreCheckpoint()     override;

    //Member functions that are not (as yet?) overridden:
    #if (0)
	virtual RateHeterogeneity* getSpecifiedRateModel(PhyloTree* tree) override;
    virtual string getName() const override;
    virtual std::string getNameParams() const override;
    virtual bool fixParameters(bool fix) override;
    virtual bool isSiteSpecificModel() const override; //Hell no!
    virtual bool isMixture() const override; //Hell no!
    virtual bool isPolymorphismAware() override; //Hell no!
	virtual int  getNMixtures() const; //Always 1
	virtual double getMixtureWeight(int cat) const override; //Always 1.0
    virtual ModelSubst* getMixtureClass(int cat) const override;
    virtual int getNumRateEntries() const override; //num_states squared
	virtual int getPtnModelID(int ptn) const; //might as well return 0
	virtual StateFreqType getFreqType() const;
        //Currently always returns StateFreqType::FREQ_EQUAL; }
        //Should perhaps return the frequency type of the first 
        //divergent model.
    virtual ModelSubst *getMutationModel();

    #endif

    //Member functions that *can't* legitimately be called 
    //(and will result in assertion failures if they are)

    virtual void    setMixtureWeight(int cat, double weight) override;
    virtual void    setFixMixtureWeight(bool fix_prop) override;
    virtual void    setMixtureClass(int cat, ModelSubst* m) override;
	virtual int     getTransMatrixSize() const override;
	virtual void    computeTransMatrix(double time, double *trans_matrix, 
                                       int mixture = 0) const override;
	virtual double  computeTrans    (double time, 
                                     int state1, int state2) 
                                     const override;
	virtual double  computeTrans    (double time, 
                                     int state1, int state2, 
                                     double &derv1, double &derv2) 
                                     const override;
    virtual void    getRateMatrix   (double *rate_mat) const override;
	virtual void    getQMatrix      (double *q_mat) const override;
	virtual void    computeTransDerv(double time, double *trans_matrix, 
		                             double *trans_derv1, double *trans_derv2, 
					  			     int mixture = 0) const override;
	virtual double* getEigenvalues()                   const override;
	virtual double* getEigenvectors()                  const override;
	virtual double* getInverseEigenvectors()           const override;
    virtual double* getInverseEigenvectorsTransposed() const override;


};

#endif