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
#ifndef SUPERALIGNMENTPAIRWISE_H
#define SUPERALIGNMENTPAIRWISE_H

#include "alignmentpairwise.h"
#include "superalignment.h"
#include "tree/phylosupertree.h"

/**
	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class SuperAlignmentPairwise : public AlignmentPairwise
{
public:
    SuperAlignmentPairwise();
	/**
		construct the pairwise alignment from two sequences of a multiple alignment
		@param aln input multiple alignment
		@param seq_id1 ID of the first sequence
		@param seq_id2 ID of the second sequence
	*/
	SuperAlignmentPairwise(PhyloSuperTree *atree, int seq1, int seq2);

    ~SuperAlignmentPairwise();

	/**
		compute the likelihood for a distance between two sequences. Used for the ML optimization of the distance.
		@param value x-value of the function
		@return log-likelihood 
	*/
	virtual double computeFunction(double value);


	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		used by Newton raphson method to minimize the function.
		@param value x-value of the function
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
		@return f(value) of function f you want to minimize
	*/
	virtual void computeFuncDerv(double value, double &df, double &ddf);

	vector<AlignmentPairwise> partitions;

};

#endif
