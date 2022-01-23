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
 * a set of substitution models, used eg for site-specific state frequency model or 
 * partition model with joint branch lengths
 */
class ModelSet : public ModelMarkov
{
public:
	typedef ModelMarkov super;
    ModelSet(const char* model_name, PhyloTree* tree);
	/**
	 * @return TRUE if this is a site-specific model, FALSE otherwise
	 */
	virtual bool isSiteSpecificModel() const override { return true; }

	/**
	 * get the size of transition matrix, default is num_states*num_states.
	 * can be changed for e.g. site-specific model
	 */
	virtual int getTransMatrixSize() const override { 
		return num_states * num_states 
	           * static_cast<int>(models.size()); 
	}

	/**
		compute the transition probability matrix.
		@param time time between two events
        @param mixture (optional) class for mixture model
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, 
	                                int mixture = 0) const override;

	
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
		                          double *trans_derv1, double *trans_derv2,
								  int mixture = 0) const override;

	/**
		To AVOID 'hides overloaded virtual functions
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2) 
	               const override { return 0; }

	/**
		To AVOID 'hides overloaded virtual functions
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
		@param derv1 (OUT) 1st derivative
		@param derv2 (OUT) 2nd derivative
	*/
	virtual double computeTrans(double time, int state1, int state2, 
	                            double &derv1, double &derv2) 
								const override { return 0; }



	/**
		compute the transition probability between two states at a specific site 
		One should override this function when defining new model.
		The default is the Juke-Cantor model, valid for all kind of data (DNA, AA, Codon, etc)
		@param time time between two events
		@param model_id model ID
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int model_id, 
	                            int state1, int state2) const override;

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
	virtual double computeTrans(double time, int model_id, 
	                            int state1, int state2, 
	                            double &derv1, double &derv2) const override;

	/**
	 * @return pattern ID to model ID map, useful for e.g., partition model
	 * @param ptn pattern ID of the alignment
	 */
	virtual int getPtnModelID(int ptn) const override;
	
	/**
		return the number of dimensions
	*/
	virtual int getNDim() const override;
	

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out) override;

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix() override;

    virtual ~ModelSet();

    /**
     * compute the memory size for the model, can be large for site-specific models
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired() const override {
    	uint64_t mem = ModelMarkov::getMemoryRequired();
    	for (auto it = models.begin(); it != models.end(); ++it) {
    		mem += (*it)->getMemoryRequired();
		}
    	return mem;
    }

	/** map from pattern ID to model ID */
	IntVector pattern_model_map;

    /**
        join memory for eigen into one chunk
    */
    void joinEigenMemory();

	ModelMarkov* firstModel() const;
	void addModel(ModelMarkov* model_to_add);


protected:
	vector<ModelMarkov*> models;
	
	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters 
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables) override;

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters 
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(const double *variables) override;

	virtual void getDivergentModels
             (DivergentModels& div_models) override;

};

#endif // MODELSET_H
