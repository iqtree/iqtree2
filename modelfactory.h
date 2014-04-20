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
#ifndef MODELFACTORY_H
#define MODELFACTORY_H

#include "tools.h"
#include "modelsubst.h"
#include "rateheterogeneity.h"

/**
Store the transition matrix corresponding to evolutionary time so that one must not compute again. 
For efficiency purpose esp. for protein (20x20) or codon (61x61).
The values of the map contain 3 matricies consecutively: transition matrix, 1st, and 2nd derivative

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelFactory : public unordered_map<int, double*>, public Optimization
{
public:

	/**
		constructor
		create substitution model with possible rate heterogeneity. Create proper class objects
		for two variables: model and site_rate. It takes the following field of params into account:
			model_name, num_rate_cats, freq_type, store_trans_matrix
		@param params program parameters
		@param tree associated phylogenetic tree
	*/
	ModelFactory(Params &params, PhyloTree *tree);

	/**
		blank constructor
	*/
	
	ModelFactory();

	/**
	 * create a substitution model
	 * @param model_str model nme
	 * @param freq_type state frequency type
	 * @param freq_params frequency parameters
	 * @param tree associated phylo tree
	 * @param count_rates TRUE to assign rates counted from alignment, FALSE to not initialize rates
	 * @return substitution model created
	 */
	ModelSubst *createModel(string model_str, StateFreqType freq_type, string freq_params,
			PhyloTree *tree, bool count_rates = true);

	/**
	 * read site specific state frequency vectors from a file to create corresponding model (Ingo's idea)
	 * @param aln input alignment
	 * @param site_freq_file file name
	 * @param site_model (OUT) site to model ID map
	 * @param freq_vec (OUT) vector of frequency vectors
	 */
	void readSiteFreq(Alignment *aln, char* site_freq_file, IntVector &site_model, vector<double*> &freq_vec);

	/**
		get the name of the model
	*/
	//string getModelName();

	void writeInfo(ostream &out);

	/**
		Start to store transition matrix for efficiency
	*/
	void startStoringTransMatrix();

	/**
		Stop storing transition matrix, e.g., when optimizing model parameters
	*/
	void stopStoringTransMatrix();

	/**
		Wrapper for computing the transition probability matrix from the model. It use ModelFactory
		that stores matrix computed before for effiency purpose.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
	*/
	void computeTransMatrix(double time, double *trans_matrix);

	/**
	 * wrapper for computing transition matrix times state frequency vector
	 * @param time time between two events
	 * @param state_freq state frequency vector
	 * @param trans_matrix (OUT) the transition matrix between all pairs of states.
	 * 	Assume trans_matrix has size of num_states * num_states.
	 */
	void computeTransMatrixFreq(double time, double *state_freq, double *trans_matrix);

	/**
		Wrapper for computing the transition probability between two states.
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	double computeTrans(double time, int state1, int state2);

	/**
		Wrapper for computing the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int state1, int state2, double &derv1, double &derv2);

	/**
		Wrapper for computing the transition probability matrix and the derivative 1 and 2 from the model.
		It use ModelFactory that stores matrix computed before for effiency purpose.
		@param time time between two events
		@param trans_matrix (OUT) the transition matrix between all pairs of states. 
			Assume trans_matrix has size of num_states * num_states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states. 
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states. 
	*/
	void computeTransDerv(double time, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2);

	void computeTransDervFreq(double time, double rate_val, double *state_freq, double *trans_matrix, 
		double *trans_derv1, double *trans_derv2);

	/**
		 destructor
	*/
    virtual ~ModelFactory();

    /**
     * @return #parameters of the model + # branches
     */
    virtual int getNParameters();

	/**
		optimize model parameters and tree branch lengths
		@param fixed_len TRUE to fix branch lengths, default is false
		@return the best likelihood 
	*/
	virtual double optimizeParameters(bool fixed_len = false, bool write_info = true, double epsilon = 0.001);

	/**
		pointer to the model, will not be deleted when deleting ModelFactory object
	*/
	ModelSubst *model;


	/**
		pointer to the site-rate heterogeneity, will not be deleted when deleting ModelFactory object
	*/
	RateHeterogeneity *site_rate;

	/**
		TRUE to store transition matrix into this hash table for computation efficiency
	*/
	bool store_trans_matrix;

	/**
		TRUE for storing process
	*/
	bool is_storing;

	/**
	 * encoded constant sites that are unobservable and added in the alignment
	 * this involves likelihood function for ascertainment bias correction for morphological or SNP data (Lewis 2001)
	 */
	string unobserved_ptns;

	/**
	 * optimize model and site_rate parameters
	 */
	double optimizeParametersOnly(double epsilon);

	/************* FOLLOWING FUNCTIONS SERVE FOR JOINT OPTIMIZATION OF MODEL AND RATE PARAMETERS *******/

	/**
	 * TRUE to optimize all parameters simultaneously, default: FALSE
	 */
	bool joint_optimize;
	/**
		return the number of dimensions
	*/
	virtual int getNDim();

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

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

};

#endif
