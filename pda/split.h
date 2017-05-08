/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#ifndef SPLIT_H
#define SPLIT_H

#include <vector>
#include <string>
#include "utils/tools.h"
#include "nclextra/msplitsblock.h"

using namespace std;


const int UINT_BITS = sizeof(UINT) * 8;
const int BITS_DIV = (sizeof(int) == 2) ? 4 : ((sizeof(int)==4) ? 5 : 6);
const int BITS_MODULO = UINT_BITS-1;

/**
Defining a split, also a set of taxa.

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class Split : public vector<UINT>
{
public:

	friend class MSplitsBlock;
	friend class SplitGraph;
	friend class PDNetwork;
	friend class CircularNetwork;
	friend class PDTree;
    friend class MTreeSet;

	/**
		empty constructor
	*/
    Split();

	/**
		constructor
		@param antaxa number of taxa
		@param aweight weight of split
	*/
	Split(int antaxa, double aweight = 0.0);


	/**
		constructor copy from another split
		@param sp split to be copied from
	*/
	Split(const Split &sp);

	/**
		construct the split from a taxa list
		@param antaxa number of taxa
		@param aweight weight of split
		@param taxa_list list of taxa in one side of split
	*/
    Split(int antaxa, double aweight, vector<int> taxa_list);

	/**
		print infos of split graph
		@param out the output stream
	*/	
	void report(ostream &out);

	/**
		destructor
	*/
    ~Split();

	/**
		get number of taxa
		@return number of taxa
	*/
	inline int getNTaxa() const {
		return ntaxa;
	}

	/**
		get number of taxa being in the split
		@return number of taxa
	*/
	int countTaxa() const;

	/**
		copy from another split
	*/
	//void copy(Split &sp); 

	/**
		set number of taxa
		@param antaxa number of taxa
	*/
	void setNTaxa(int antaxa);

	/**
		get the first taxon in the set
		@return the first taxon or -1 if empty
	*/
	int firstTaxon();

	/**
		@return TRUE if the set is empty
	*/
	bool isEmpty();

	/**
		get weight
		@return weight
	*/
	inline double getWeight() const {
		return weight;
	}

	/**
		set weight
		@param aweight the new weight
	*/
	inline void setWeight(double aweight) {
		weight = aweight;
	}

	/**
		check whether the split should be inverted (number of taxa > ntaxa / 2)
		@return TRUE yes, should be inverted
	*/
	bool shouldInvert();

	/**
		invert the split (0->1, 1->0)
	*/
	void invert();

	/**
		@param sp the other split
		@return true if this split is compatible with sp
	*/
	bool compatible(Split &sp);

	/**	
		@param taxa_set set of taxa
		@return true if this split is preserved in the set taxa_set
	*/
	bool preserved(Split &taxa_set);

	/**
		if the split is trivial (contains 1 taxon in 1 side), return the taxon id, 
		otherwise return -1
		@return the taxon id if the split is trivial 
	*/
	int trivial();

	/**
		add a taxon into the split
		@param tax_id id of taxon from 0..ntaxa-1
	*/
	void addTaxon(int tax_id); 

	/**
		remove a taxon from the split
		@param tax_id id of taxon from 0..ntaxa-1
	*/
	void removeTaxon(int tax_id); 

	/**
		@param tax_id id of taxon from 0..ntaxa-1
		@return true if tax_id is in the set
	*/
	bool containTaxon(int tax_id); 

	/**
		@param tax_id vector of id of taxa from 0..ntaxa-1
		@return true if SOME taxon in tax_id is in the set
	*/
	bool containAny(IntVector &tax_id); 

	/**
		@param tax_id vector of id of taxa from 0..ntaxa-1
		@return true if ALL taxa in tax_id is in the set
	*/
	bool containAll(IntVector &tax_id); 

	/**
		get the list of taxa contained in split
		@param invec (OUT) taxa in this side of split
	*/
	void getTaxaList(vector<int> &invec);

	/**
		get the list of taxa contained in split and not contained in split
		@param invec (OUT) taxa in this side of split
		@param outvec (OUT) taxa on the other side
	*/
	void getTaxaList(vector<int> &invec, vector<int> &outvec);

	/**
	 *  Test whether the current split is smaller than \a sp
	 *  @param sp the other split to compare
	 *  @return true if the current split contains less taxa than \a sp
	 */
    bool operator<(const Split &sp) const;

	/**
		compare two split, do not compare the weight
		@param sp the target split
		@return TRUE if equal, FALSE otherwise
	*/
	bool operator==(const Split &sp) const;

	/**
		add all taxa from another split into this split (union)
		@param sp a split
	*/
	Split &operator+=(Split &sp);

	/**
		get the intersection with another split
		@param sp a split
	*/
	Split &operator*=(Split &sp);

	/**
		get the set difference with another split
		@param sp a split
	*/
	Split &operator-=(Split &sp);

	/**
		@return TRUE if there is overlapped taxon with sp, FALSE otherwise
		@param sp a split
	*/
	bool overlap(Split &sp);

	/**
		assignment
		@param sp a split
	*/
	Split &operator= (const Split &sp);

	/**
		subset operator
		@param sp a split
		@return TRUE of this set is a subset of sp
	*/
	bool subsetOf (Split &sp);


	/**
		randomize the set of a specific size
		@param size number of taxa in the resulting set
	*/
	void randomize(int size);

	Split *extractSubSplit(Split &taxa_mask);

	string &getName() { return name; }

protected:
	/**
		number of taxa
	*/
	int ntaxa; 

	/**
		weight of split
	*/
	double weight;
    
    /** 2018-08-23: split name */
    string name;

};

inline int splitweightcmp(const Split* a, const Split* b)
{
	return (a->getWeight() > b->getWeight());
}

typedef Split TaxaSet;

#endif
