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
#ifndef ALIGNMENTPAIRWISE_H
#define ALIGNMENTPAIRWISE_H

#include "alignment.h"
#include "optimization.h"
#include "modelfactory.h"
#include "rateheterogeneity.h"

/**
Pairwise alignment

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class AlignmentPairwise : public Alignment, public Optimization
{
public:
    AlignmentPairwise();

	/**
		construct the pairwise alignment from two sequences of a multiple alignment
		@param aln input multiple alignment
		@param seq_id1 ID of the first sequence
		@param seq_id2 ID of the second sequence
	*/
    AlignmentPairwise(Alignment *aln, int seq_id1, int seq_id2);

	/**
		compute the likelihood for a distance between two sequences. Used for the ML optimization of the distance.
		@param value x-value of the function
		@return log-likelihood 
	*/
	virtual double computeFunction(double value);

	/**
		compute the ML distance between two sequences
		@return the ML distance
	*/
	double optimizeDist();


	/**
		destructor
	*/
    ~AlignmentPairwise();

	/**
		pairwise state frequencies
	*/
	int *pair_freq;

	/**
		stores transition matrices computed before for efficiency purpose, eps. AA or CODON model.
	*/
	ModelFactory *model_factory;

	/**
		among-site rates 
	*/
	RateHeterogeneity *site_rate;

};

#endif
