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
#ifndef GTRMODEL_H
#define GTRMODEL_H

#define EIGEN
#include "phylotree.h"
#include "modelsubst.h"
#include "optimization.h"
#include "alignment.h"
#include "eigendecomposition.h"

const double MIN_RATE = 1e-4;
const double TOL_RATE = 1e-4;
const double MAX_RATE = 100;


/**
General Time Reversible (GTR) model of substitution.
This works for all kind of data, not only DNA

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelGTR : public ModelSubst, public EigenDecomposition
{
	
	friend class ModelSet;
	friend class ModelMixture;
    friend class ModelPoMo;
	
public:
	/**
		constructor
		@param tree associated tree for the model
	*/
    ModelGTR(PhyloTree *tree, bool count_rates = true);

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
    virtual ~ModelGTR();

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
	virtual string getNameParams();

    /**
        internal function: return string for frequency
        @param retname output stream
    */
    void getNameParamsFreq(ostream &retname);

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
		@param file_name file containing upper-triangle rate matrix and state frequencies
	*/
	void readParameters(const char *file_name);


	/**
		compute the transition probability matrix.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix);

	
	/**
	 * wrapper for computing transition matrix times state frequency vector
	 * @param time time between two events
	 * @param trans_matrix (OUT) the transition matrix between all pairs of states.
	 * 	Assume trans_matrix has size of num_states * num_states.
	 */
	virtual void computeTransMatrixFreq(double time, double *trans_matrix);

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
		compute the state frequency vector
		@param state_freq (OUT) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq);

	/**
		set the state frequency vector
		@param state_freq (IN) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void setStateFrequency(double *state_freq);

	/**
	 * compute Q matrix 
	 * @param q_mat (OUT) Q matrix, assuming of size num_states * num_states
	 */
	virtual void getQMatrix(double *q_mat);

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
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	virtual void computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2);

	/**
		compute the transition probability matrix.and the derivative 1 and 2 times state frequency vector
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	virtual void computeTransDervFreq(double time, double rate_val, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2);

	/**
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
		@return the number of dimensions corresponding to state frequencies
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


	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

//	double *getEigenCoeff() const;

	virtual double *getEigenvalues() const;

	virtual double *getEigenvectors() const;
	virtual double *getInverseEigenvectors() const;

//	void setEigenCoeff(double *eigenCoeff);

	void setEigenvalues(double *eigenvalues);

	void setEigenvectors(double *eigenvectors);

    /**
     * compute the memory size for the model, can be large for site-specific models
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired() {
    	return ModelSubst::getMemoryRequired() + sizeof(double)*num_states*num_states*3;
    }

    /** default TRUE: store only upper half of the rate matrix */
    bool half_matrix;

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

	virtual void freeMem();

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
		coefficient cache, served for fast computation of the P(t) matrix
	*/
//	double *eigen_coeff;

	/** state with highest frequency, used when optimizing state frequencies +FO */
	int highest_freq_state;

};

#endif
