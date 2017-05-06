/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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
#ifndef CIRCULARNETWORK_H
#define CIRCULARNETWORK_H

#include "pdnetwork.h"

/**
Circular Network for PDA algorithm

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class CircularNetwork : public PDNetwork
{
public:

	/**
		empty constructor
	*/
    CircularNetwork();

	/**
		construct network from a NEXUS file, e.g. produced by SplitsTree
		@param params program parameters
	*/
    CircularNetwork(Params &params);


	/**
		MAIN FUNCTION which will call other other findPD depending on the input splits.
		Search for maximal phylogenetic diversity of a given size 
		@param params program parameters
		@param taxa_set (OUT) the set of taxa in the maximal PD set
		@param taxa_order (OUT) order of inserted taxa
	*/
	virtual void findPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order);

/********************************************************
	Dynamic programming strategy
********************************************************/

	/**
		dynamic programming algorithm in UNROOTED circular splits graph 
			for phylogenetic diversity of a given size 
		@param params parameters
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (IN) order of inserted taxa
	*/
	void findCircularPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order);

	/**
		dynamic programming algorithm in ROOTED circular splits graph 
			for phylogenetic diversity of a given size 
		@param params parameters
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (IN) order of inserted taxa
		@param root index of the root taxon
	*/
	void findCircularRootedPD(Params &params, vector<SplitSet> &taxa_set, 
		vector<int> &taxa_order, int root);

	/**
		dynamic programming algorithm with cost-constrained in UNROOTED circular splits graph 
			for phylogenetic diversity under budget constraint
		@param params program parameters
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (IN) order of inserted taxa
		@return the PD score of the maximal set, also returned in taxa_set.weight
	*/
	void findCircularPDBudget(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order);
	
	/**
		dynamic programming algorithm with cost-constrained in ROOTED circular splits graph 
			for phylogenetic diversity under budget constraint
		@param params program parameters
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (IN) order of inserted taxa
		@return the PD score of the maximal set, also returned in taxa_set.weight
		@param root index of the root taxon
	*/
	void findCircularRootedPDBudget(Params &params, vector<SplitSet> &taxa_set, 
		vector<int> &taxa_order, int root);


protected:

/********************************************************
	CIRCULAR NETWORKS
********************************************************/

	/**
		compute the PD information table
		@param params program parameters
		@param table (OUT) computed information
		@param dist distance matrix
		@param root index of the root taxon
	*/
	void computePDInfo(Params &params, DoubleMatrix &table, DoubleMatrix  &dist, int root);

	/**
		compute the PD score
		@param sub_size the subset size
		@param table computed information
		@param root index of the root taxon
	*/
	double computePDScore(int sub_size, DoubleMatrix &table, int root);


	/**
		construct optimal PD set from computed information for ROOTED case circular network
		@param sub_size the subset size
		@param find_all TRUE of want to find all PD sets
		@param pd_limit maximum number of returned PD sets
		@param table computed information
		@param dist distance matrix
		@param taxa_set (OUT) sets of taxa with optimal PD
		@param taxa_order circular order
		@param root the root
	*/
	void constructPD(int sub_size, bool find_all, int pd_limit, DoubleMatrix &table, DoubleMatrix  &dist, 
		SplitSet &taxa_set, vector<int> &taxa_order, int root);

	/**
		construct optimal PD set from computed information for ROOTED case circular network
		@param sub_size the subset size
		@param max_v end taxon
		@param pd_limit maximum number of returned PD sets
		@param pd_set the current constructed PD set
		@param table computed information
		@param dist distance matrix
		@param taxa_set (OUT) sets of taxa with optimal PD
		@param taxa_order circular order
		@param root the root
	*/
	void constructPD(int sub_size, int max_v, int pd_limit, Split *pd_set, DoubleMatrix &table, 
		DoubleMatrix  &dist, SplitSet &taxa_set, vector<int> &taxa_order, int root);


/********************************************************
	CIRCULAR NETWORKS WITH BUDGET CONSTRAINT
********************************************************/


	/**
		calculate the maximum budget required from u to v (excluding u and v)
		@param budget total budget
		@param max_b (OUT) max budget matrix between taxa
		@param taxa_order circular order		
	*/
	void calcMaxBudget(int budget, mmatrix(int) &max_b, vector<int> &taxa_order);

	/**
		construct optimal PD set from computed information for budget constraint (ROOTED case)
		@param budget total budget
		@param find_all TRUE of want to find all PD sets
		@param table computed information
		@param dist distance matrix
		@param taxa_set (OUT) sets of taxa with optimal PD
		@param taxa_order circular order
		@param max_b max budget matrix between taxa
		@param root the root
	*/
	void constructPDBudget(int budget, bool find_all, mmatrix(double) &table, 
		mmatrix(double) &dist,SplitSet &taxa_set, 
		vector<int> &taxa_order, mmatrix(int) &max_b, int root);


	/**
		construct optimal PD set from computed information for budget constraint (ROOTED case)
		@param budget total budget
		@param max_v end taxon
		@param pd_set the current constructed PD set
		@param table computed information
		@param dist distance matrix
		@param taxa_set (OUT) sets of taxa with optimal PD
		@param taxa_order circular order
		@param max_b max budget matrix between taxa
		@param root the root
	*/
	void constructPDBudget(int budget, int max_v, Split *pd_set, 
		mmatrix(double) &table, mmatrix(double) &dist, SplitSet &taxa_set, 
		vector<int> &taxa_order, mmatrix(int) &max_b, int root);

	/**
		compute the PD information table with budget
		@param params program parameters
		@param table (OUT) computed information
		@param id (OUT) computed information
		@param dist distance matrix
		@param taxa_order circular order
		@param max_b max budget matrix between taxa
		@param root index of the root taxon
	*/
	void computePDBudgetInfo(Params &params, mmatrix(double) &table, mmatrix(int) &id, 
		mmatrix(double) &dist, vector<int> &taxa_order, mmatrix(int) &max_b, int root);

	/**
		compute the PD score with budget
		@param budget total budget
		@param table (OUT) computed information
		@param dist distance matrix
		@param taxa_order circular order
		@param max_b max budget matrix between taxa
		@param root index of the root taxon
	*/
	double computePDBudgetScore(int budget, mmatrix(double) &table,
		mmatrix(double) &dist, vector<int> &taxa_order, mmatrix(int) &max_b, int root);


};

/**
	display the matrix into out
*/
template <class T>
ostream &operator<<(ostream &out, mmatrix(T) &mat) {
	unsigned int i, j;
	for (i = 0; i < mat.size(); i++) {
		for (j = 0; j < mat[i].size(); j++) {
			if (j < i) 
				out << " &  "; 
			else if (j < mat[i].size()-1) 
				out << mat[i][j] << " & ";
			else
				out << mat[i][j];
		}
		if (i < mat.size()-1)
			out << " \\\\";
		out << endl;
	}
	return out;
} 



#endif
