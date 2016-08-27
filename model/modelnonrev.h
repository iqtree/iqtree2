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
#ifndef MODELNONREV_H
#define MODELNONREV_H

#include "phylotree.h"
#include "modelgtr.h"
#include <complex>

/**
The general non-reversible model

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelNonRev : public ModelGTR
{
public:
    ModelNonRev(PhyloTree *tree);

    virtual ~ModelNonRev();
    /**
     * Return a model of type given by model_name. (Will be some subclass of ModelNonRev.)
     */
    static ModelNonRev* getModelByName(string model_name, PhyloTree *tree, string model_params, bool count_rates);

    /**
     * true if model_name is the name of some known non-reversible model
     */
    static bool validModelName(string model_name);

    /** 
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

	/**
		@return TRUE if model is time-reversible, FALSE otherwise
	*/
	virtual bool isReversible() { return false; };

	/**
		@return the number of dimensions
	*/
	virtual int getNDim();

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

	/**
		optimize model parameters
		@return the best likelihood 
	*/
	virtual double optimizeParameters(double gradient_epsilon);

	/**
		@return the number of rate entries, equal to the number of non-diagonal elements
			of the rate matrix (since model is NOT reversible)
	*/
	virtual int getNumRateEntries() { return num_states*(num_states-1); }

	virtual void getQMatrix(double *rate_mat);
	
	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	/**
		compute the transition probability matrix.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix);

	/**
		compute the transition probability matrix using (complex) eigenvalues
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	void computeTransMatrixEigen(double time, double *trans_matrix);

	/**
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

protected:

	/*
	 * MDW:
	 * I strongly object to the function naming here - it completely threw me.
	 * a 'get' method should copy data from the object, and a 'set' method should
	 * write data into the object. These are named completely the wrong way around!
	 */

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

	virtual void freeMem();

	/**
	    Model parameters - cached so we know when they change, and thus when
	    recalculations are needed.

	 */
	double *model_parameters;

	/** true to fix parameters, otherwise false */
	bool fixed_parameters;

	/**
		unrestricted Q matrix. Note that Q is normalized to 1 and has row sums of 0.
		no state frequencies are involved here since Q is a general matrix.
	*/
	double *rate_matrix;

	/** imaginary part of eigenvalues */
	double *eigenvalues_imag;
	
	/**
		temporary working space
	*/
	double *temp_space;
    
    /**
        complex eigenvalues and eigenvectors, pointing to the same pointer 
        to the previous double *eigenvalues and double *eigenvectors
    */
    std::complex<double> *ceval, *cevec, *cinv_evec;
    
};

#endif
