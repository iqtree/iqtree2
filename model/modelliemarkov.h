/*
 * modelliemarkov.h
 *
 *  Created on: 24/05/2016
 *      Author: Michael Woodhams
 */

#ifndef MODELLIEMARKOV_H_
#define MODELLIEMARKOV_H_

#include "modelmarkov.h"

class ModelLieMarkov: public ModelMarkov {
public:
        ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params);
        virtual ~ModelLieMarkov();

	/**
		this function is served for model testing
		@param model_name name of the model
		@param freq_type state frequency type, can be FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, or FREQ_ESTIMATE
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

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
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	static void getLieMarkovModelInfo(string model_name, string &name, string &full_name, int &model_num, int &symmetry, StateFreqType &def_freq);

	static string getModelInfo(string model_name, string &full_name, StateFreqType &def_freq);

	// DO NOT override this function, because
    // BQM, 2017-05-02: getNDimFreq should return degree of freedom, which is not included in getNDim()
    // That's why 0 is returned for FREQ_ESTIMATE, num_states-1 for FREQ_EMPIRICAL
//	virtual int getNDimFreq();

    // this is redundant, there is already the same function below
//	bool isTimeReversible();

	/**
		@return TRUE if model is time-reversible, FALSE otherwise
	*/
	virtual bool isReversible();
    
    /**
     * @return TRUE if this is a liemarkov model, FALSE otherwise
     */
    virtual bool isLieMarkov() { return true; }

    /**
         initialize random state frequencies when running AliSim without inference mode
    */
    void initStateFreqsAliSim(StateFreqType expected_freq_type);
    
    /**
         read user-specified state frequencies
    */
    void readFreqs(StateFreqType expected_freq_type, string freq_params);
    
    /**
         mapping state frequencies from user-specified/random frequencies
    */
    void mappingFreqs(StateFreqType expected_freq_type, double *freqs);

	static bool validModelName(string model_name);
	void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);
	
	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

    /** decompose rate matrix using closed formula derived by Michael Woodhams */
    void decomposeRateMatrixClosedForm();

    /** decompose rate matrix using Eigen library */
    virtual void decomposeRateMatrixEigen3lib();

	/**
		compute the transition probability matrix.
		@param time time between two events
        @param mixture (optional) class for mixture model
        @param selected_row (optional) only compute the entries of one selected row. By default, compute all rows
		@param trans_matrix (OUT) the transition matrix between all pairs of states.
			Assume trans_matrix has size of num_states * num_states.
	*/
	virtual void computeTransMatrix(double time, double *trans_matrix, int mixture = 0, int selected_row = -1);
	// overrides Optimization::restartParameters
	bool restartParameters(double guess[], int ndim, double lower[], double upper[], bool bound_check[], int iteration);

protected:
	/**
	    Model parameters - cached so we know when they change, and thus when
	    recalculations are needed.

	 */
	double *model_parameters;


	double **basis;
	int symmetry; // RY->0, WS->1, MK->2
	int model_num; // 0->1.1, etc to 36->12.12
	void setBasis();
	virtual void setRates();

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

	static void parseModelName(string model_name, int* model_num, int* symmetry);
	/*
         * Overrides ModelMarkov::getName().
	 * Avoids appending +FO to name, as this is implied by how LM models 
	 * work.
         * Minh: you might chose to remove this override, if you like "+FO"
	 * to be on LM model names.
         */
	string getName();

	/*
	const static double ***BASES;
	const static int *MODEL_PARAMS;
	const static string *SYMMETRY;
	const static string *MODEL_NAMES;
	const static int NUM_RATES;
	const static int NUM_LM_MODELS;
	*/
        bool validFreqType();

};
#endif /* MODELLIEMARKOV_H_ */
