/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  BUI Quang Minh <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef MODELSET_H
#define MODELSET_H

#include "modelmarkov.h"

/**
 * A set of substitution models to implement site-specific models
 */
class ModelSet : public ModelMarkov, public vector<ModelMarkov*>
{

public:
	ModelSet(const string model_name, ModelsBlock *models_block,
		StateFreqType freq, string freq_params, PhyloTree *tree);

	virtual ~ModelSet();

	/**
		set checkpoint object
		@param checkpoint
	*/
	virtual void setCheckpoint(Checkpoint *checkpoint);

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
	 * @return TRUE if this is a site-specific model, FALSE otherwise
	 */
	virtual bool isSiteSpecificModel() { return true; }

	/**
	 * @return TRUE if the model has site-specific frequencies, FALSE otherwise
	 */
	virtual bool isSSF() { return phylo_tree->aln->isSSF(); }

	/**
	 * @return TRUE if the model has site-specific rates, FALSE otherwise
	 */
	virtual bool isSSR() { return phylo_tree->aln->isSSR(); }

	/**
	 * get the size of transition matrix, default is num_states*num_states.
	 * can be changed for e.g. site-specific model
	 */
	virtual int getTransMatrixSize() { return num_states * num_states * size(); }

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
		get the rate matrix.
		@param rate_mat (OUT) upper-triagle rate matrix.
		Assume rate_mat has size of num_states*(num_states-1)/2
	*/
	virtual void getRateMatrix(double *rate_mat);

	/**
		compute the state frequency vector. One should override this function when defining new model.
		The default is equal state sequency, valid for all kind of data.
		@param mixture (optional) class for mixture model
		@param[out] state_freqs state frequency vector. Assume state_freqs has size of num_states
	*/
	virtual void getStateFrequency(double *state_freqs, int mixture = 0);

	/**
		compute Q matrix
		@param q_mat (OUT) Q matrix, assuming of size num_states * num_states
	*/
	virtual void getQMatrix(double *q_mat, int mixture = 0);

	/**
		get frequency type
		@return frequency type
	*/
	virtual StateFreqType getFreqType();

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
		@return TRUE if parameters are at the boundary that may cause numerical unstability
	*/
	virtual bool isUnstableParameters();

	/**
		setup the bounds for joint optimization with BFGS
	*/
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		rescale the state frequencies
		@param sum_one TRUE to make frequencies sum to 1, FALSE to make last entry equal to 1
	*/
	virtual void scaleStateFreq(bool sum_one);

	/**
		optimize model parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon);

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

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
		compute the transition probability matrix.and the derivative 1 and 2
		@param time time between two events
		@param mixture (optional) class for mixture model
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
		@param trans_derv1 (OUT) the 1st derivative matrix between all pairs of states.
		@param trans_derv2 (OUT) the 2nd derivative matrix between all pairs of states.
		Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransDerv(double time, double *trans_matrix,
		double *trans_derv1, double *trans_derv2, int mixture = 0);

	/**
		To AVOID 'hides overloaded virtual functions
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2) { return 0; }

	/**
		To AVOID 'hides overloaded virtual functions
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int state1, int state2, double &derv1, double &derv2) { return 0; }

	/**
		compute the transition probability between two states at a specific site
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param model_id model ID
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int model_id, int state1, int state2);

	/**
		compute the transition probability and its 1st and 2nd derivatives between two states at a specific site
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
	virtual int getPtnModelID(int ptn);

	/**
		compute the memory size for the model, can be large for site-specific models
		@return memory size required in bytes
	*/
	virtual uint64_t getMemoryRequired();

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

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
		join memory for eigen into one chunk
	*/
	void joinEigenMemory();

};

#endif // MODELSET_H
