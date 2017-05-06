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
#ifndef STOPRULE_H
#define STOPRULE_H

#include "tools.h"
#include "checkpoint.h"

/**
Stopping rule
	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class StopRule : public CheckpointFactory
{
public:

	/**
		constructor
	*/
    StopRule();

    void initialize(Params &params);
	/**
		destructor
	*/
    ~StopRule();

    void getUFBootCountCheck(int &ufboot_count, int &ufboot_count_check);

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

	/**
		read improved iteration number from a file
		@param fileName file name
	*/
	void readFile (const char *fileName);

	/**
		Add the iteration number that improve trees
		@param iteration improved iteration number
	*/
	void addImprovedIteration(int iteration);

	/**
		Get the last iteration number that improved trees
		@return the last iteration number that improved trees
	*/
	int getLastImprovedIteration();

	/**
		main function to check the stop condition
		@param current_iteration current iteration number
		@param cur_correlation current correlation coefficient for bootstrap convergence
		@return TRUE if stop condition is met, FALSE otherwise
	*/
	bool meetStopCondition(int cur_iteration, double cur_correlation);

    /**
        return TRUE if cur_correlation is high enough
        @param cur_correlation correlation coefficient
    */
    bool meetCorrelation(double cur_correlation) {
        return cur_correlation >= min_correlation;
    }
	
	/** get the remaining time to converge, in seconds */
	double getRemainingTime(int cur_iteration);

	/**
		@return the number of iterations required to stop the search
	*/
//	int getNumIterations();

	/**
		@return predicted iteration, 0 if no prediction has been made
	*/
//	int getPredictedIteration(int cur_iteration);


    int getCurIt() const {
        return curIteration;
    }

    void setCurIt(int curIteration) {
        StopRule::curIteration = curIteration;
    }

    void shouldStop() {
        should_stop = true;
    }

private:

    /**
	 *  Current iteration number
	 */
	int curIteration;

	double predict (double &upperTime);

	/**
		stop condition 
	*/
	STOP_CONDITION stop_condition;	
	
	/**
		confidence value of prediction
	*/
	double confidence_value;

	/**
		minimum number of iterations
	*/
	int min_iteration;

	/**
		maximum number of iterations
	*/
	int max_iteration;

	/**
		predicted number of iterations
	*/
	int predicted_iteration;

	/** number of unsuccessful iterations to stop the search */
	int unsuccess_iteration;

	/** bootstrap correlation threshold to stop */
	double min_correlation;

	/** step size for checking bootstrap convergence */
	int step_iteration;

	/** max wall-clock running time to stop */
	double max_run_time;

    /** starting real time of the program */
    double start_real_time;

    /** TRUE to override stop condition */
    bool should_stop;

	/* FOLLOWING CODES ARE FROM IQPNNI version 3 */	

//	int nTime_;
	DoubleVector time_vec;

	void cmpInvMat (DoubleMatrix &oriMat, DoubleMatrix &invMat, int size);

	void readMat (char *fileName, DoubleMatrix &oriMat, int &size);

	void multiple (DoubleMatrix &mat1, DoubleMatrix &mat2, DoubleMatrix &proMat);


	void multiple (DoubleMatrix &mat1, DoubleVector &vec2, DoubleVector &proVec);

	void multiple (DoubleVector &vec1, DoubleMatrix &mat2, DoubleVector &proVec);
	void multiple (DoubleVector &vec1, DoubleVector &vec2, DoubleMatrix &proMat);
	double multiple (DoubleVector &vec1, DoubleVector &vec2);

	void readVector(DoubleVector &tmpTimeVec_);

	/* THE FOLLOWING CODE COMES FROM tools.c in Yang's PAML package */
	//----------------------------------------------------------------------------------------
	double cmpLnGamma (double alpha);

	double cmpMuy (int k);

	void cmpLamdaMat (int k, DoubleMatrix &lamdaMat);

	void cmpVecA (int k, DoubleVector &aVec);

	double cmpExtinctTime (int k);
	double cmpUpperTime (int k, double alpha);
};

#endif
