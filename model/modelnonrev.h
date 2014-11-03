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

/**
The general non-reversible model

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelNonRev : public ModelGTR
{
public:
    ModelNonRev(PhyloTree *tree, bool count_rates = true);

	/**
		@return TRUE if model is time-reversible, FALSE otherwise
	*/
	virtual bool isReversible() { return false; };

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
		compute the transition probability between two states
		@param time time between two events
		@param state1 first state
		@param state2 second state
	*/
	virtual double computeTrans(double time, int state1, int state2);

protected:

	virtual void freeMem();

	/**
		unrestricted Q matrix. Note that Q is normalized to 1 and has row sums of 0.
		no state frequencies are involved here since Q is a general matrix.
	*/
	double *rate_matrix;
	
	/**
		temporary working space
	*/
	double *temp_space;
};

#endif
