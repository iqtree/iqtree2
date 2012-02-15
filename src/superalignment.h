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
#ifndef SUPERALIGNMENT_H
#define SUPERALIGNMENT_H

#include "alignment.h"

class PhyloSuperTree;

/**
Super alignment representing m partitions for a total of n sequences. It has the form:
		Site_1 Site_2 ... Site_m
Seq_1     1      0    ...   1
Seq_2     0      1    ...   0
...      ...
Seq_n     1      1    ...   0

Where (i,j)=1 means Seq_i is present in partition j, 0 otherwise

So data is binary.

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/

class SuperAlignment : public Alignment
{
public:
    SuperAlignment(PhyloSuperTree *super_tree);
    SuperAlignment();

    ~SuperAlignment();

	virtual bool isSuperAlignment() { return true; }

	/**
		Quit if some sequences contain only gaps or missing data
	*/
	virtual void checkGappySeq();

	/**
		create a non-parametric bootstrap alignment from an input alignment
		@param aln input alignment
	*/
	virtual void createBootstrapAlignment(Alignment *aln);

	/**
		compute the observed distance (number of different pairs of positions per site) 
			between two sequences
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2
		@return the observed distance between seq1 and seq2 (between 0.0 and 1.0)
	*/
	virtual double computeObsDist(int seq1, int seq2);

	/**
		compute the Juke-Cantor corrected distance between 2 sequences over all partitions
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2		
		@return any distance between seq1 and seq2
	*/
	virtual double computeDist(int seq1, int seq2);

	void printCombinedAlignment(const char *filename, bool append = false);

	/**
		@return unconstrained log-likelihood (without a tree)
	*/
	virtual double computeUnconstrainedLogL();

	/**
		actual partition alignments
	*/
	vector<Alignment*> partitions;

	/**
		matrix represents the index of taxon i in partition j
	*/
	vector<IntVector> taxa_index;

};

#endif
