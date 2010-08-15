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

#include "substmodel.h"
#include "optimization.h"
#include "phylotree.h"
#include "alignment.h"
#include "eigendecomposition.h"


/**
General Time Reversible (GTR) model of substitution.
This works for all kind of data, not only DNA

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class GTRModel : public SubstModel, public Optimization, public EigenDecomposition
{
public:
	/**
		constructor
		@param tree associated tree for the model
	*/
    GTRModel(PhyloTree *tree);
	
	/**
		init the model and decompose the rate matrix. This function should always be called
		after creating the class. Otherwise it will not work properly.
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	void init(StateFreqType freq_type);

	/**
		this function is served for ModelDNA and ModelProtein
		@param model_name name of the model
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	virtual void init(const char *model_name, StateFreqType freq_type) {}

	/**
		destructor
	*/
    virtual ~GTRModel();

	/**
		set the associated tree
		@param tree the associated tree
	*/
    void setTree(PhyloTree *tree);


	/**
		compute the transition probability matrix.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix);

	/**
		Get the rate matrix.
		@param rate_mat (OUT) upper-triagle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void getRateMatrix(double *rate_mat);

	/**
		compute the state frequency vector
		@param state_freq (OUT) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq);


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
		get frequency type
		@return frequency type
	*/
	StateFreqType getFreqType() { return freq_type; }

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
	virtual double optimizeParameters();


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
	void decomposeRateMatrix();

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
		state frequencies
	*/
	double *state_freq;
	
	/**
		state frequency type
	*/
	StateFreqType freq_type;

	/**
		eigenvalues of the rate matrix Q
	*/
	double *eigenvalues;

	/**
		eigenvectors of the rate matrix Q 
	*/
	double **eigenvectors;

	/**
		inverse eigenvectors of the rate matrix Q 
	*/
	double **inv_eigenvectors;

	/**
		coefficient cache, served for fast computation of the P(t) matrix
	*/
	double *eigen_coeff;


};

#endif
