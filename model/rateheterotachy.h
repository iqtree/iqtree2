//
//  rateheterotachy.hpp
//  iqtree
//
//  Created by Minh Bui on 11/8/16.
//
//

#ifndef rateheterotachy_hpp
#define rateheterotachy_hpp


#include "rateheterogeneity.h"

class PhyloTree;

/**
    rate-heterotachy model, allowing for mixed branch lengths
*/
class RateHeterotachy: public RateHeterogeneity {

    friend class ModelFactoryMixlen;

public:
	typedef RateHeterogeneity super;
	/**
		constructor
		@param ncat number of rate categories
        @param sorted_rates true to sort the rate in increasing order
		@param tree associated phylogenetic tree
	*/
    RateHeterotachy(int ncat, string params, PhyloTree* tree);

	/**
		constructor (used by, for example, YAMLRateModelWrapper)
		@param ncat number of rate categories
		@param tree associated phylogenetic tree
		@param report_to_tree send any log messages to this tree.
	*/
	RateHeterotachy(int ncat, PhyloTree* tree, PhyloTree* report_to_tree);

	/**
		destructor
	*/
    virtual ~RateHeterotachy();

    /**
        @return TRUE if this is a heterotachy model, default: FALSE
    */
    virtual bool isHeterotachy() const override { return true; }

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
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual std::string getNameParams() const override;


    /**
        fix parameters, so that no optimization done
        @param mode some input mode
    */
    virtual int getFixParams() const override { return fix_params; }

    /**
        fix parameters, so that no optimization done
        @param mode some input mode
    */
    virtual void setFixParams(int mode) override { fix_params = mode; }

	/**
		return the number of dimensions
	*/
	virtual int getNDim() const override;

	/**
		@return the number of rate categories
	*/
	virtual int getNRate() const override { return ncategory; }

	/**
		get the number of rate categories for site-specific category model
		@return the number of rate categories
	*/
	virtual int getNDiscreteRate() const override { return ncategory; }

	/**
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) const override {
        return 1.0;
    }

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) const override { return prop[category]; }

	/**
		set the proportion of a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual void setProp(int category, double value) override { prop[category] = value; }


    /** 
        set number of optimization steps
        @param opt_steps number of optimization steps
    */
    virtual void setOptimizeSteps(int steps) override { this->optimize_steps = steps; }

	/**
		optimize parameters. Default is to optimize gamma shape
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon,
                                      PhyloTree* report_to_tree) override;

    /** 
        optimize rate parameters using EM algorithm 
        @return log-likelihood of optimized parameters
    */
    double optimizeWithEM(PhyloTree* report_to_tree);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) override;

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out) override;

    /**
        set number of rate categories
        @param ncat #categories
    */
	virtual void setNCategory           (int ncat) override;

	virtual bool isOptimizingProportions() const;
	virtual bool isOptimizingRates      () const;
	virtual bool isOptimizingShapes     () const;
	virtual void sortUpdatedRates       ();

	virtual void setFixProportions      (bool   fix);
	virtual void setFixRates            (bool   fix);
	virtual void setProportionTolerance (double tol) override;
	virtual void setRateTolerance       (double tol) override;

protected:

	/**
		number of rate categories
	*/
	int ncategory;

	/**
	 * proportion of sites for each rate categories
	 */
    double *prop;

	/** TRUE to fix parameters */
	int fix_params;

    /** number of optimization steps, default: ncategory*2 */
    int optimize_steps;

    /** tolerance (for deciding when proportions have converged) */
	double prop_tolerance;

    /** minimum value for a proportion */
	double prop_minimum;

    /** tolerance (for deciding when rates have converged) */
	double rate_tolerance;

};




#endif /* rateheterotachy_hpp */
