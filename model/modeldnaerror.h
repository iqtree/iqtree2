//
//  modeldnaerror.hpp
//  model
//
//  Created by Minh Bui on 26/6/20.
//

#ifndef modeldnaerror_hpp
#define modeldnaerror_hpp

#include "modeldna.h"

/**
DNA models with sequecing error

    @author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelDNAError : public ModelDNA
{
public:

    /**
        constructor
        @param tree associated tree for the model
    */
    ModelDNAError(PhyloTree *tree);

    /**
        constructor
        @param model_name model name, e.g., JC, HKY.
        @param freq state frequency type
        @param tree associated phylogenetic tree
    */
    ModelDNAError(const char *model_name, string model_params, StateFreqType freq, string freq_params, string seqerr, PhyloTree *tree);

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    /**
     * @return model name
     */
    virtual string getName();

    /**
     * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
     */
    virtual string getNameParams(bool show_fixed_params = false);

    /**
        write information
        @param out output stream
    */
    virtual void writeInfo(ostream &out);


    /**
        return the number of dimensions
    */
    virtual int getNDim();

    /**
     * setup the bounds for joint optimization with BFGS
     */
    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);


    /** compute the tip likelihood vector of a state for Felsenstein's pruning algorithm
     @param state character state
     @param[out] state_lk state likehood vector of size num_states
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);
    
    /**
     * @return TRUE if this is a DNA error model, FALSE otherwise
     */
    virtual bool containDNAerror() { return true; }
    
    /**
     * get the dna error probability, by default error probability = 0
     */
    virtual double getDNAErrProb(int mixture_index = 0) { return epsilon; }

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
        @return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
    */
    virtual bool getVariables(double *variables);

private:
    
    /** sequencing error */
    double epsilon;
    
    /** true to fix epsilon */
    bool fix_epsilon;
    
    /** error model name */
    string seqerr_name;
    
};

#endif /* modeldnaerror_hpp */
