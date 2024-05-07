//
// C++ Interface: substmodel
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SUBSTMODEL_H
#define SUBSTMODEL_H

#include <string>
#include "utils/tools.h"
#include "utils/optimization.h"
#include "utils/checkpoint.h"
#include "phylo-yaml/statespace.h"

using namespace std;

const char OPEN_BRACKET = '{';
const char CLOSE_BRACKET = '}';

class PhyloTree;

/**
Substitution model abstract class

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class ModelSubst: public Optimization, public CheckpointFactory
{
	friend class ModelFactory;
    friend class PartitionModel;
    friend class IQTreeMix;

public:
	/**
		constructor
		@param nstates number of states, e.g. 4 for DNA, 20 for proteins.
	*/
    ModelSubst(int nstates);


	/**
		@return the number of dimensions
	*/
	virtual int getNDim() { return 0; }

	/**
		@return the number of dimensions corresponding to state frequencies
	*/
	virtual int getNDimFreq() { return 0; }
	
	/**
	 * @return model name
	 */
	virtual string getName() { return name; }

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams(bool show_fixed_params = false) { return name; }

	/**
		@return TRUE if model is time-reversible, FALSE otherwise
	*/
	virtual bool isReversible() { return true; };

    /** return true if using reversible likelihood kernel, false for using non-reversible kernel */
    bool useRevKernel() {
        return isReversible() && !Params::getInstance().kernel_nonrev;
    };

    /**
        fix parameters of the model
        @param fix true to fix, false to not fix
        @return the current state of fixing parameters
     */
    virtual bool fixParameters(bool fix) {
        bool current = fixed_parameters;
        fixed_parameters = fix;
        return current;
    }
    
	/**
	 * @return TRUE if this is a site-specific model, FALSE otherwise
	 */
	virtual bool isSiteSpecificModel() { return false; }

	/**
	 * @return TRUE if this is a mixture model, FALSE otherwise
	 */
	virtual bool isMixture() { return false; }
    
    /**
     * @return TRUE if this is a liemarkov model, FALSE otherwise
     */
    virtual bool isLieMarkov() { return false; }
    
    /**
     * @return TRUE if this is a mixture model and all model components share the same rate matrix, FALSE otherwise
     */
    virtual bool isMixtureSameQ() { return false; }
    
    /**
     * @return TRUE if this is a DNA error model, FALSE otherwise
     */
    virtual bool containDNAerror() { return false; }
    
    /**
     * get the dna error probability, by default error probability = 0
     */
    virtual double getDNAErrProb(int mixture_index = 0) { return 0; }

    /** 
     * Confer to modelpomo.h.
     * 
     * @return TRUE if PoMo is being used, FALSE otherise.
     */
    virtual bool isPolymorphismAware() { return false; }

	/**
	 * @return the number of mixture model components
	 */
	virtual int getNMixtures() { return 1; }

 	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual double getMixtureWeight(int cat) { return 1.0; }

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual void setMixtureWeight(int cat, double weight) {}

	/**
	 * @param cat mixture class
	 * @return weight of a mixture model component
	 */
	virtual void setFixMixtureWeight(bool fix_prop) {}

	/**
	 * @param cat mixture class ID
	 * @return corresponding mixture model component
	 */
    virtual ModelSubst* getMixtureClass(int cat) { return NULL; }

	/**
	 * @param cat mixture class ID
	 * @param m mixture model class to set
	 */
    virtual void setMixtureClass(int cat, ModelSubst* m) { }

	/**
		@return the number of rate entries, equal to the number of elements
			in the upper-diagonal of the rate matrix (since model is reversible)
	*/
	virtual int getNumRateEntries() { return num_states*(num_states-1)/2; }
    
    /**
     set num_params variable
     */
    virtual void setNParams(int num_params) {}
    
    /**
     get num_params variable
     */
    virtual int getNParams() {
        return 0;
    }

	/**
	 * get the size of transition matrix, default is num_states*num_states.
	 * can be changed for e.g. site-specific model
	 */
	virtual int getTransMatrixSize() { return num_states * num_states; }

	/**
		compute the transition probability matrix. One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
        @param mixture (optional) class for mixture model
        @param selected_row (optional) only compute the entries of one selected row. By default, compute all rows
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, int mixture = 0, int selected_row = -1);

	/**
		compute the transition probability between two states. 
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2);

	/**
		compute the transition probability between two states at a specific model ID, useful for partition model 
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param model_id model ID
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int model_id, int state1, int state2);

	/**
		compute the transition probability and its 1st and 2nd derivatives between two states. 
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int state1, int state2, double &derv1, double &derv2);

	/**
		compute the transition probability and its 1st and 2nd derivatives between two states at a specific model ID
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param model_id model ID
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int model_id, int state1, int state2, double &derv1, double &derv2);

	/**
	 * @return pattern ID to model ID map, useful for e.g., partition model
	 * @param ptn pattern ID of the alignment
	 */
	virtual int getPtnModelID(int ptn) { return 0; }
	

	/**
	 * Get the rate parameters like a,b,c,d,e,f for DNA model!!!
		Get the above-diagonal entries of the rate matrix, assuming that the last element is 1.
		ONE SHOULD OVERRIDE THIS FUNCTION WHEN DEFINING NEW MODEL!!!
		The default is equal rate of 1 (JC Model), valid for all kind of data.
		@param rate_mat (OUT) upper-triangle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	
	virtual void getRateMatrix(double *rate_mat);

	/**
		Get the rate matrix Q. One should override this function when defining new model.
		The default is equal rate of 1 (JC Model), valid for all kind of data.
		@param rate_mat (OUT) upper-triagle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void getQMatrix(double *q_mat, int mixture = 0);

	/**
		compute the state frequency vector. One should override this function when defining new model.
		The default is equal state sequency, valid for all kind of data.
        @param mixture (optional) class for mixture model
		@param[out] state_freq state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq, int mixture = 0);

    /**
     set the state frequency vector.
     @param state_freq state frequency vector. Assume state_freq has size of num_states
     */
    virtual void setStateFrequency(double *state_freq);

	/**
		get frequency type
		@return frequency type
	*/
	virtual StateFreqType getFreqType() { return FREQ_EQUAL; }

    /**
        set the associated tree
        @param tree the associated tree
    */
    virtual void setTree(PhyloTree *tree) {}
    
    /** for reversible models, multiply likelihood with inverse eigenvectors for fast pruning algorithm
            @param[in/out] state_lk state likelihood multiplied with inverse eigenvectors
     */
    void multiplyWithInvEigenvector(double *state_lk);

    /** compute the tip likelihood vector of a state for Felsenstein's pruning algorithm
     @param state character state
     @param[out] state_lk state likehood vector of size num_states
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);
    
	/**
		allocate memory for a transition matrix. One should override this function when defining new model
		such as Gamma model. The default is to allocate a double vector of size num_states * num_states. This
		is equivalent to the memory needed by a square matrix.
		@return the pointer to the newly allocated transition matrix
	*/
	virtual double *newTransMatrix();


	/**
		compute the transition probability matrix.and the derivative 1 and 2
		@param time time between two events
        @param mixture (optional) class for mixture model
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	virtual void computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2, int mixture = 0);

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix() {}


    /** 
        set number of optimization steps
        @param opt_steps number of optimization steps
    */
    virtual void setOptimizeSteps(int optimize_steps) { }

	/**
		optimize model parameters. One should override this function when defining new model.
		The default does nothing since it is a Juke-Cantor type model, hence no parameters involved.
		@param epsilon accuracy of the parameters during optimization
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double gradient_epsilon) { return 0.0; }

	/**
	 * @return TRUE if parameters are at the boundary that may cause numerical unstability
	 */
	virtual bool isUnstableParameters() { return false; }

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) {}

	/**
		report model
		@param out output stream
	*/
    virtual void report(ostream &out) {}

	virtual double *getEigenvalues() const {
		return NULL;
	}

	virtual double *getEigenvectors() const {
		return NULL;
	}

	virtual double *getInverseEigenvectors() const {
		return nullptr;
	}

    virtual double *getInverseEigenvectorsTransposed() const {
        return nullptr;
    }

    
    /**
     * compute the memory size for the model, can be large for site-specific models
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired() {
    	return num_states*sizeof(double);
    }
    
    /** @return true if model is a mixture model and it's fused with site_rate */
    virtual bool isFused(){
        return false;
    };

    /**
    * get the underlying mutation model, used with PoMo model
    */
    virtual ModelSubst *getMutationModel() { return this; }

	/*****************************************************
		Checkpointing facility
	*****************************************************/

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
		number of states
	*/
	int num_states;

	/**
		name of the model
	*/
	string name;


	/**
		full name of the model
	*/
	string full_name;
    
    /** true to fix parameters, otherwise false */
    bool fixed_parameters;

	/**
	 state frequencies
	 */
	double *state_freq;
	

	/**
		state frequency type
	*/
	StateFreqType freq_type;

    /** state set for each sequence in the alignment */
    //vector<vector<int> > seq_states;

	/**
		destructor
	*/
    virtual ~ModelSubst();

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables) {}

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(double *variables) { return false; }

    
};

#endif
