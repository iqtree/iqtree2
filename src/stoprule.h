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


/**
Stopping rule
	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class StopRule
{
public:

	/**
		constructor
	*/
    StopRule();

	/**
		destructor
	*/
    ~StopRule();

	/**
		@param sc stop condition
	*/
	void setStopCondition(STOP_CONDITION sc);

	/**
		@param confidence_val confidence value for the prediction
	*/
	void setConfidenceValue(double confidence_val);

	/**
		set the number of iterations
		@param min_it minimum iteration
		@param max_it maximum iteration
	*/
	void setIterationNum(int min_it, int max_it);

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
		main function to check the stop condition
		@param current_iteration current iteration number
		@return TRUE if stop condition is met, FALSE otherwise
	*/
	bool meetStopCondition(int current_iteration);

private:

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

	/* FOLLOWING CODES ARE FROM IQPNNI version 3 */	

	int nTime_;
	DoubleVector timeVec_;

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
