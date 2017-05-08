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
#include "superalignmentpairwise.h"

SuperAlignmentPairwise::SuperAlignmentPairwise()
 : AlignmentPairwise()
{
}

SuperAlignmentPairwise::SuperAlignmentPairwise(PhyloSuperTree *atree, int seq1, int seq2) 
 : AlignmentPairwise() 
{
	tree = atree;
	seq_id1 = seq1;
	seq_id2 = seq2;
	SuperAlignment *aln = (SuperAlignment*) atree->aln;
	int part = 0;
	for (PhyloSuperTree::iterator it = atree->begin(); it != atree->end(); it++, part++) {
		int id1 = aln->taxa_index[seq1][part];
		int id2 = aln->taxa_index[seq2][part];
		if (id1 >= 0 && id2 >= 0)
		partitions.push_back(new AlignmentPairwise((*it), id1, id2));
	}
}

double SuperAlignmentPairwise::computeFunction(double value) {
	double lh = 0.0;
	for (vector<AlignmentPairwise*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		lh += (*it)->computeFunction(value);
	}
	return lh;
}


void SuperAlignmentPairwise::computeFuncDerv(double value, double &df, double &ddf) {
//	double lh = 0.0;
	df = 0.0;
	ddf = 0.0;
	for (vector<AlignmentPairwise*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		double d1, d2;
//		lh += (*it)->computeFuncDerv(value, d1, d2);
		(*it)->computeFuncDerv(value, d1, d2);
		df += d1;
		ddf += d2;
	}
//	return lh;
}


SuperAlignmentPairwise::~SuperAlignmentPairwise()
{
	for (vector<AlignmentPairwise*>::reverse_iterator it = partitions.rbegin(); it != partitions.rend(); it++)
		delete (*it);
	partitions.clear();
}


