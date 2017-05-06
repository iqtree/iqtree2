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

/*
	Geneset selection (GSS) for Roland 
*/

#ifndef GSS_H
#define GSS_H


#include "tools.h"
#include "pda/pdnetwork.h"

class GSSNetwork : public PDNetwork {
public:
	/**
		construct PD network from a NEXUS or NEWICK file, e.g. produced by SplitsTree
		@param params program parameters
	*/
    GSSNetwork(Params &params);

	/**
		transform the problem into an Integer Linear Programming and write to .lp file
		@param params program parameters
		@param outfile name of output file in LP format
		@param total_size k for PD_k or total budget
		@param make_bin TRUE if creating binary programming
	*/
	void transformLP_GSS(Params &params, const char *outfile, int total_size, bool make_bin);

	/**
		main function to search for maximal phylogenetic diversity
		@param params program parameters
		@param taxa_set (OUT) the vector of set of taxa in the maximal PD set
		@param taxa_order (OUT) order of inserted taxa
	*/
	virtual void findPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order);

	/**
		@return TRUE if we are doing PD area optimization
	*/
	virtual bool isPDArea();

	void readGenePValues(Params &params);

protected:

	/**
		names of the genes
	*/
	StrVector genes;

	map<string, int> gene_index;
	
	/**
		p-values of the genes
	*/
	DoubleVector gene_pvalues;


	/**
		z variables for genes in the LP formulation, check if it can be dropped or equals some x variable.
		@param total_size k for PD_k or total budget
		@param z_value (OUT): vector of: -1 if cannot reduce, 1 if equals 1, or id+2 where id is the trivial split id 
	*/
	void checkZValue(int total_size, vector<int> &z_value);

	void lpObjectiveGSS(ostream &out, Params &params, IntVector &y_value, IntVector &z_value, int total_size);

	void lpVariableBound(ostream &out, Params &params, Split &included_vars, IntVector &y_value, IntVector &z_value);

	void lpGeneConstraint(ostream &out, Params &params, IntVector &z_value);

};

void runGSSAnalysis(Params &params);

#endif
