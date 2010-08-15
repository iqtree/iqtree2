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

using namespace std;

/**
Substitution model abstract class

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class SubstModel
{
public:
	/**
		constructor
		@param nstates number of states, e.g. 4 for DNA, 20 for proteins.
	*/
    SubstModel(int nstates);

	/**
		compute the transition probability matrix. One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix);


	/**
		Get the rate matrix. One should override this function when defining new model.
		The default is equal rate of 1 (JC Model), valid for all kind of data.
		@param rate_mat (OUT) upper-triagle rate matrix. Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void getRateMatrix(double *rate_mat);

	/**
		compute the state frequency vector. One should override this function when defining new model.
		The default is equal state sequency, valid for all kind of data.
		@param state_freq (OUT) state frequency vector. Assume state_freq has size of num_states
	*/
	virtual void getStateFrequency(double *state_freq);

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
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	virtual void computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2);


	/**
		optimize model parameters. One should override this function when defining new model.
		The default does nothing since it is a Juke-Cantor type model, hence no parameters involved.
		@return the best likelihood 
	*/
	virtual double optimizeParameters() { return 0.0; }

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) {}

	/**
		number of states
	*/
	int num_states;

	/**
		name of the model
	*/
	string name;

	/**
		destructor
	*/
    virtual ~SubstModel();

};

#endif
