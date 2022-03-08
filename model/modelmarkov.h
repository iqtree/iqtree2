/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef MODELMARKOV_H
#define MODELMARKOV_H

#include "tree/phylotree.h"
#include "modelsubst.h"
#include "utils/optimization.h"
#include "alignment/alignment.h"
#include "utils/eigendecomposition.h"
#include <complex>

const double MIN_RATE = 1e-4;
const double TOL_RATE = 1e-4;
const double MAX_RATE = 100;

string freqTypeString(StateFreqType freq_type, SeqType seq_type, bool full_str);

/**
General Markov model of substitution (reversible or non-reversible)
This works for all kind of data

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelMarkov : public ModelSubst, public EigenDecomposition
{
	
	friend class ModelSet;
	friend class ModelMixture;
    friend class ModelPoMo;
    friend class PartitionModel;
    friend class PartitionModelPlen;
	
public:
	/**
		constructor
		@param tree associated tree for the model
        @param reversible TRUE (default) for reversible model, FALSE for non-reversible
        @param adapt_tree TRUE (default) to convert rooted<->unrooted tree
	*/
    ModelMarkov(PhyloTree *tree, bool reversible = true, bool adapt_tree = true);

	/**
		@return TRUE if model is time-reversible, FALSE otherwise
	*/
	virtual bool isReversible() { return is_reversible; };

    /**
        set the reversibility of the model
        @param reversible TRUE to make model reversible, FALSE otherwise
        @param adapt_tree TRUE (default) to convert between rooted and unrooted tree
    */
    virtual void setReversible(bool reversible, bool adapt_tree = true);


	/**
		init the model and decompose the rate matrix. This function should always be called
		after creating the class. Otherwise it will not work properly.
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	void init(StateFreqType freq_type);

	/**
	   initializes ModelSubst::freq_type array according to freq_type
           (can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE)
	 */
	void init_state_freq(StateFreqType freq_type);

	/**
		this function is served for ModelDNA and ModelProtein
		@param model_name name of the model
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params) {}

	/**
		destructor
	*/
    virtual ~ModelMarkov();

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
        internal function: return string for frequency
        @param retname output stream
    */
    void getNameParamsFreq(ostream &retname);

	/**
		@return the number of rate entries, equal to the number of non-diagonal elements
			of the rate matrix (since model is NOT reversible)
	*/
	virtual int getNumRateEntries();

	/**
		set the associated tree
		@param tree the associated tree
	*/
    void setTree(PhyloTree *tree);

	/**
		Read the upper-triangle rate matrix from an input stream. 
		It will throw error messages if failed
		@param in input stream
	*/
	virtual void readRates(istream &in) throw(const char*, string);

	/**
		Read the rate parameters from a comma-separated string
		It will throw error messages if failed
		@param in input stream
	*/
	virtual void readRates(string str) throw(const char*);

	/**
		Read state frequencies from an input stream. 
		It will throw error messages if failed
		@param in input stream
	*/
	virtual void readStateFreq(istream &in) throw(const char*);

	/**
		Read state frequencies from comma-separated string
		It will throw error messages if failed
		@param str input string
	*/
	virtual void readStateFreq(string str) throw(const char*);

	/**
		read model parameters from a file
		@param file_name file containing rate matrix and state frequencies
	*/
	void readParameters(const char *file_name, bool adapt_tree = true);

	/**
		read model parameters from a string
		@param model_str string containing rate matrix and state frequencies
	*/
	void readParametersString(string &model_str, bool adapt_tree = true);

	/**
		compute the transition probability matrix.
		@param time time between two events
        @param mixture (optional) class for mixture model
        @param selected_row (optional) only compute the entries of one selected row. By default, compute all rows
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, int mixture = 0, int selected_row = -1);

    /**
     compute the transition probability matrix for non-reversible model
     @param time time between two events
     @param mixture (optional) class for mixture model
     @param trans_matrix (OUT) the transition matrix between all pairs of states.
     Assume trans_matrix has size of num_states * num_states.
     */
    virtual void computeTransMatrixNonrev(double time, double *trans_matrix, int mixture = 0);

	/**
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2);

	/**
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int state1, int state2, double &derv1, double &derv2);

	/**
		Get the rate matrix.
		@param rate_mat (OUT) upper-triagle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void getRateMatrix(double *rate_mat);

	/**
		Set the rate matrix.
		@param rate_mat upper-triagle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void setRateMatrix(double *rate_mat);

    /**
     Set the full rate matrix of size num_states*num_states
     @param rate_mat full rate matrix
     @param freq state frequency
     */
    virtual void setFullRateMatrix(double *rate_mat, double *freq);

	/**
		compute the state frequency vector
        @param mixture (optional) class for mixture model
		@param state_freq (OUT) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq, int mixture = 0);

	/**
		set the state frequency vector
		@param state_freq (IN) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void setStateFrequency(double *state_freq);

    /**
     set the state frequency vector
     @param state_freq (IN) state frequency vector. Assume state_freq has size of num_states
     */
    virtual void adaptStateFrequency(double *state_freq);

	/**
	 * compute Q matrix 
	 * @param q_mat (OUT) Q matrix, assuming of size num_states * num_states
	 */
	virtual void getQMatrix(double *q_mat, int mixture = 0);

	/**
		rescale the state frequencies
		@param sum_one TRUE to make frequencies sum to 1, FALSE to make last entry equal to 1
	*/
	void scaleStateFreq(bool sum_one);

	/**
		get frequency type
		@return frequency type
	*/
	virtual StateFreqType getFreqType() { return freq_type; }


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
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
		@return the number of dimensions corresponding to state frequencies, which is 
            not counted in getNDim(). This serves e.g. for computing AIC, BIC score
	*/
	virtual int getNDimFreq();
	

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		optimize model parameters
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double gradient_epsilon);

	/**
	 * @return TRUE if parameters are at the boundary that may cause numerical unstability
	 */
	virtual bool isUnstableParameters();

  // A simple helper function that prints the rates in a nice way and can be
  // reused by children. The title is necessary, because, e.g., for PoMo, the
  // rates are mutation rates and not substitution rates, and also
  // exchangeabilities may be reported.
  void report_rates(ostream &out, string title, double *r);

  // Report the stationary frequencies state_freq or custom_state_freq (if
  // given) to output stream out.
  void report_state_freqs(ostream &out, double *custom_state_freq=NULL);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out){}

    /** decompose rate matrix for non-reversible models */
    virtual void decomposeRateMatrixNonrev();

    /** old version of decompose rate matrix for reversible models */
    void decomposeRateMatrixRev();

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

//	double *getEigenCoeff() const;

	virtual double *getEigenvalues() const;

	virtual double *getEigenvectors() const;
	virtual double *getInverseEigenvectors() const;
    virtual double *getInverseEigenvectorsTransposed() const;

//	void setEigenCoeff(double *eigenCoeff);

	void setEigenvalues(double *eigenvalues);

	void setEigenvectors(double *eigenvectors);

	void setInverseEigenvectors(double *inv_eigenvectors);

    void setInverseEigenvectorsTransposed(double *inv_eigenvectors);

    static void calculateExponentOfScalarMultiply(const double* source, int size
                                                  , double scalar, double* dest);
    
    static void calculateHadamardProduct(const double* first, const double* second
                                         , int size, double *dest);
    
    static double dotProduct(const double* first, const double* second, int size);
    
    static void calculateSquareMatrixTranspose(const double* original, int rank
                                               , double* transpose);
    
    static void aTimesDiagonalBTimesTransposeOfC
        ( const double* matrixA, const double* rowB
         , const double* matrixCTranspose, int rank
         , double* dest, int selected_row = -1);
    /**
     * compute the memory size for the model, can be large for site-specific models
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired() {
    	return ModelSubst::getMemoryRequired() + sizeof(double)*num_states*num_states*3;
    }

    /** default TRUE: store only upper half of the rate matrix */
    bool half_matrix;

    /****************************************************/
    /*      NON-REVERSIBLE STUFFS                       */
    /****************************************************/

    /**
     * Return a model of type given by model_name. (Will be some subclass of ModelMarkov.)
     */
    static ModelMarkov* getModelByName(string model_name, PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params);

    /**
     * true if model_name is the name of some known non-reversible model
     */
    static bool validModelName(string model_name);

  // Mon Jul 3 14:47:08 BST 2017; added by Dominik. I had problems with mixture
  // models together with PoMo and rate heterogeneity. E.g., a model
  // "MIX{HKY+P+N9+G2,GTR+P+N9+G2}" leads to segmentation faults because the
  // `ModelPoMoMixture` reports a /wrong/ number of states (i.e., it reports 52
  // instead of 104). Consequently, the `initMem()` function of ModelMixture,
  // messes up the `eigenvalues`, etc., variables of the `ModelPoMoMixture`s. I
  // circumvent this, by adding this virtual function; for normal models, it
  // just returns `num_states`, however, for mixture models, it returns
  // `num_states*nmixtures`.
  virtual int get_num_states_total();

  // Mon Jul 3 15:53:00 BST 2017; added by Dominik. Same problem as with
  // `get_num_states_total()`. The pointers to the eigenvalues and eigenvectors
  // need to be updated recursively, if the model is a mixture model. For a
  // normal Markov model, only the standard pointers are set. This was done in
  // `ModelMixture::initMem()` before.
  virtual void update_eigen_pointers(double *eval, double *evec
                                     , double *inv_evec, double *inv_evec_transposed);


    /**
        set num_params variable
     */
    virtual void setNParams(int num_params) {
        this->num_params = num_params;
    }
    
    /**
        get num_params variable
     */
    virtual int getNParams() {
        return num_params;
    }
    
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


	/**
	 * Called from getVariables to update the rate matrix for the new
	 * model parameters.
	 */
	virtual void setRates();

    /**
        free all allocated memory 
		(internal version not affected by virtual dispatch; safe to call from destructor)
    */
	void internalFreeMem();

	/**
		free all allocated memory
	*/
	virtual void freeMem();

    /** TRUE if model is reversible */
    bool is_reversible;

	/**
		phylogenetic tree associated
	*/
	PhyloTree *phylo_tree;

    /**
		rates between pairs of states of the unit rate matrix Q.
		In order A-C, A-G, A-T, C-G, C-T (rate G-T = 1 always)
	*/
	double *rates;

	/**
		the number of free rate parameters
	*/
	int num_params;

	/**
		eigenvalues of the rate matrix Q
	*/
	double *eigenvalues;

	/**
		eigenvectors of the rate matrix Q 
	*/
	double *eigenvectors;

	/**
		inverse eigenvectors of the rate matrix Q 
	*/
	double *inv_eigenvectors;

    /**
        transpose of the matrix of the inverse eigenvectors of the rate matrix Q
    */
    double *inv_eigenvectors_transposed;
    
	/** state with highest frequency, used when optimizing state frequencies +FO */
	int highest_freq_state;

    /****************************************************/
    /*      NON-REVERSIBLE STUFFS                       */
    /****************************************************/

	/**
		compute the transition probability matrix using (complex) eigenvalues
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	void computeTransMatrixEigen(double time, double *trans_matrix);

	/**
		unrestricted Q matrix. Note that Q is normalized to 1 and has row sums of 0.
		no state frequencies are involved here since Q is a general matrix.
	*/
	double *rate_matrix;

	/** imaginary part of eigenvalues */
	double *eigenvalues_imag;
	
    /**
        complex eigenvalues and eigenvectors, pointing to the same pointer 
        to the previous double *eigenvalues and double *eigenvectors
    */
    std::complex<double> *ceval, *cevec, *cinv_evec;

    /** will be set true for nondiagonalizable rate matrices,
     then will use scaled squaring method for matrix exponentiation.
    */
    bool nondiagonalizable;

};

#endif
