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
#include "split.h"

Split::Split()
		: vector<UINT>()
{
	ntaxa = 0;
	weight = 0.0;
}

Split::Split(int antaxa, double aweight)
		: vector<UINT>()
{
	weight = aweight;
	setNTaxa(antaxa);
}

Split::Split(const Split &sp)
		: vector<UINT>(sp) 
{
	weight = sp.weight;
	ntaxa = sp.ntaxa;
/*
	setNTaxa(sp.ntaxa);
	int i = 0;
	for (iterator it = begin(); it != end(); it++, i++)
		(*it) = sp[i];
*/
}



Split::Split(int antaxa, double aweight, vector<int> taxa_list)
		: vector<UINT>()
{
	ntaxa = antaxa;
	weight = aweight;
	vector<int>::iterator it;

	// inverted mode: include the remaining part into the split

	// if taxa_list contains more than half of taxa, turn on the inverted mode
	bool inverted = (taxa_list.size() > ntaxa / 2);
	/* if taxa_list contains exactly one haft of taxa, only turn on the inverted mode
	 if taxon 0 is not in the list */
	if (taxa_list.size() == ntaxa / 2) {
		for (it = taxa_list.begin(); it != taxa_list.end(); it++)
			if ((*it) == 0) {
				break;
			}
		inverted = true;
	}

	// resize the split size
	resize((ntaxa + UINT_BITS -1) / UINT_BITS, 0);

	for (it = taxa_list.begin(); it != taxa_list.end(); it++)
	{
		int value = *it;
		int bit_pos = value / UINT_BITS;
		int bit_off = value % UINT_BITS;
		(*this)[bit_pos] |= (UINT) (1 << bit_off);
	}


	if (inverted) invert();

}

void Split::invert() {
	for (iterator uit = begin(); uit != end(); uit++)
	{
		int num_bits = (uit+1 == end()) ? ntaxa % UINT_BITS : UINT_BITS;

		*uit = (1 << (num_bits-1)) - 1 + (1 << (num_bits-1)) - (*uit);
	}
}


bool Split::shouldInvert() {
	int count = countTaxa();
	if (count * 2 < ntaxa) 
		return false;
	if (count * 2 > ntaxa)
		return true;
	return !containTaxon(0);
}


/**
	set number of taxa
	@param antaxa number of taxa
*/
void Split::setNTaxa(int antaxa)
{
	ntaxa = antaxa;
	resize((ntaxa + UINT_BITS - 1) / UINT_BITS, 0);
	for (iterator it = begin(); it != end(); it++)
		(*it) = 0;
}

int Split::countTaxa() const {
	int count=0;
	for (int i = 0; i < size(); i++)
		for (UINT j = 0; j < UINT_BITS && (i*UINT_BITS+j < getNTaxa()); j++)
			if ((*this)[i] & (1 << j))
			{
				count++;
			}
	return count;
}

void Split::report(ostream &out)
{

	out << getWeight() << '\t';
	for (int i = 0; i < size(); i++)
		for (UINT j = 0; j < UINT_BITS && (i*UINT_BITS+j < getNTaxa()); j++)
			if ((*this)[i] & (1 << j))
			{
				out << i * UINT_BITS + j + 1 << " ";
			}
	out << endl;
}


int Split::firstTaxon() {
	for (int i = 0; i < size(); i++)
		if ((*this)[i] != 0) {
			for (UINT j = 0; j < UINT_BITS && (i*UINT_BITS+j < getNTaxa()); j++)
				if ((*this)[i] & (1 << j)) {
					return (i * UINT_BITS + j);
				}
		}
	return -1;
}

bool Split::isEmpty() {
	for (iterator it = begin(); it != end(); it++)
		if (*it != 0) return false;
	return true;
}


/**
	@param sp the other split
	@return true if this split is compatible with sp
*/
bool Split::compatible(Split &sp)
{
	// be sure that	the two split has the same size
	assert(sp.size() == size() && sp.ntaxa == ntaxa);

	UINT res = 0, res2 = 0, res3 = 0, res4 = 0;
	for (iterator it = begin(), sit = sp.begin(); it != end(); it++, sit++)
	{
		int num_bits = (it+1 == end()) ? ntaxa % UINT_BITS : UINT_BITS;
		UINT it2 = (1 << (num_bits-1)) - 1 + (1 << (num_bits-1)) - (*it);
		UINT sit2 = (1 << (num_bits-1)) - 1 + (1 << (num_bits-1)) - (*sit);

		res |= (*it) & (*sit);
		res2 |= (it2) & (sit2);
		res3 |= (*it) & (sit2);
		res4 |= (it2) & (*sit);
		if (res != 0 && res2 != 0 && res3 != 0 && res4 != 0)
			return false;
		//if (res != 0 && res != (*it) && res != (*sit) && res2 != 0)
			//return false;
	}
	return true;
	//return (res == 0) || (res2 == 0) || (res3 == 0) || (res4 == 0);
}

/**
	@param taxa_set set of taxa
	@return true if this split is preserved in the set taxa_set
*/
bool Split::preserved(Split &taxa_set)
{
	// be sure that	the two split has the same size
	assert(taxa_set.size() == size() && taxa_set.ntaxa == ntaxa);

	int time_zero = 0, time_notzero = 0;
	
	for (iterator it = begin(), sit = taxa_set.begin(); it != end(); it++, sit++)
	{
		UINT res = (*it) & (*sit);
		if (res != 0 && res != (*sit))
			return true;
		if (*sit != 0) {
			if (res == 0) time_zero++; else time_notzero++;
			if (res == 0 && time_notzero > 0) return true;
			if (res != 0 && time_zero > 0) return true;
		}
	}
	return false;
}

int Split::trivial() {
/*
	int num = countTaxa();
	if (num == 1) {
		// trivial split, fetch the first bit-1
		int tax = 0;
		for (iterator it = begin(); it != end(); it++) {
			for (int i = 0; i < UINT_BITS && tax < ntaxa; i++, tax++)
				if (((*it) & (1 << i)) != 0)
					return tax;
		}
	} else if (num == ntaxa - 1) {
		// trivial split, fetch the first bit-0
		int tax = 0;
		for (iterator it = begin(); it != end(); it++) {
			for (int i = 0; i < UINT_BITS && tax < ntaxa; i++, tax++)
				if (((*it) & (1 << i)) == 0)
					return tax;
		}
	} else 
		// not a trivial split
		return -1;
*/
	int id0 = 0, id1 = 0, pos = 0;
	int bit0s = 0, bit1s = 0;
	for (iterator it = begin(); it != end(); it++, pos++) {
		UINT content = *it;
		int max_step;
		if ((it + 1) == end()) 
			max_step = ntaxa % UINT_BITS;
		else
			max_step = UINT_BITS;

		for (int i = 0; i < max_step; i++) {
			if ((content & ( 1 << i)) != 0) {
				bit1s ++;
				if (bit1s == 1) 
					id1 = pos * UINT_BITS + i;
			}
			else {
				bit0s ++;
				if (bit0s == 1) 
					id0 = pos * UINT_BITS + i;
			}
			// if both number of bit 0 and 1 greater than 1, return -1 (not trivial)
			if (bit1s > 1 && bit0s > 1) 
				return -1;
		}
	}
	if (bit1s == 1)
		return id1;
	else if (bit0s == 1)
		return id0;
	else
		return -1;
}


/**
	add a taxon into the split
	@param tax_id id of taxon from 0..ntaxa-1
*/
void Split::addTaxon(int tax_id)
{
	assert(tax_id >= 0 && tax_id < ntaxa);
	int pos = tax_id / UINT_BITS, off = tax_id % UINT_BITS;
	(*this)[pos] |= 1 << off;
}

/**
	remove a taxon from the split
	@param tax_id id of taxon from 0..ntaxa-1
*/
void Split::removeTaxon(int tax_id)
{
	assert(tax_id >= 0 && tax_id < ntaxa);
	int pos = tax_id / UINT_BITS, off = tax_id % UINT_BITS;

	(*this)[pos] &= -1 - (1 << off);
}

/**
	@param tax_id id of taxon from 0..ntaxa-1
	@return true if tax_id is in the set
*/
bool Split::containTaxon(int tax_id)
{
	assert(tax_id >= 0 && tax_id < ntaxa);
	int pos = tax_id / UINT_BITS, off = tax_id % UINT_BITS;
	return ((*this)[pos] & ( 1 << off)) != 0;

}

void Split::getTaxaList(vector<int> &invec) {
	int tax = 0;
	invec.clear();
	for (iterator it = begin(); it != end(); it++) {
		for (int i = 0; i < UINT_BITS && tax < ntaxa; i++, tax++)
			if (((*it) & (1 << i)) != 0) // inside the split
				invec.push_back(tax);
	}
}


void Split::getTaxaList(vector<int> &invec, vector<int> &outvec) {
	int tax = 0;
	invec.clear();
	outvec.clear();
	for (iterator it = begin(); it != end(); it++) {
		for (int i = 0; i < UINT_BITS && tax < ntaxa; i++, tax++)
			if (((*it) & (1 << i)) != 0) // inside the split
				invec.push_back(tax);
			else
				outvec.push_back(tax);
	}
}

Split &Split::operator+=(Split &sp) {
	assert(sp.ntaxa == ntaxa);
	iterator it1, it2;
	for (it1 = begin(), it2 = sp.begin(); it1 != end(); it1++, it2++) {
		(*it1) |= (*it2);
	}
	return *this;
}

Split &Split::operator*=(Split &sp) {
	assert(sp.ntaxa == ntaxa);
	iterator it1, it2;
	for (it1 = begin(), it2 = sp.begin(); it1 != end(); it1++, it2++) {
		(*it1) &= (*it2);
	}
	return *this;	
}

Split &Split::operator-=(Split &sp) {
	assert(sp.ntaxa == ntaxa);
	iterator it1, it2;
	for (it1 = begin(), it2 = sp.begin(); it1 != end(); it1++, it2++) {
		(*it1) -= (*it1) & (*it2);
	}
	return *this;	
}

bool Split::operator==(const Split &sp) const{
	if (ntaxa != sp.ntaxa) return false;
	for (const_iterator it = begin(), it2 = sp.begin(); it != end(); it++, it2++)
		if ((*it) != (*it2))
			return false;
	return true;
}

bool Split::subsetOf (Split &sp) {
	assert(ntaxa == sp.ntaxa);
	for (iterator it = begin(), it2 = sp.begin(); it != end(); it++, it2++)
		if ( ((*it) & (*it2)) != (*it) )
			return false;
	return true;
}

Split &Split::operator= (const Split &sp) {
	assert(ntaxa == sp.ntaxa);
	vector<UINT>::operator= (sp);
	weight = sp.weight;
	return *this;
}

/*
void Split::copy(const Split &sp) {
	assert(ntaxa == sp.ntaxa);
	for (iterator it = begin(), it2 = sp.begin(); it != end(); it++, it2++)
		(*it) = (*it2);
	weight = sp.weight;
}
*/

void Split::randomize(int size) {
	assert(size < ntaxa);
	int num = countTaxa();
	int cnt;
	// repeat at most 10 times
	const int MAX_STEP = 20;
	const int PROB_STEP = 5;
	for (int step = 0; step < MAX_STEP && num < size; step++) {
		// probability of including a taxon
		double prob = (double)(size - num) / ntaxa;
		// increase the probability if passing too many iterations
		if (step >= PROB_STEP) prob *= 2.0;
		if (step >= PROB_STEP*2) prob *= 2.0;
		if (step == MAX_STEP - 1) prob = 1.0;
		// now scan through all elements, pick up at random
		for (cnt = 0; cnt < ntaxa && num < size; cnt++)
			if (!containTaxon(cnt) && ( ((double)(rand()) / RAND_MAX) <= prob )) {
				addTaxon(cnt);
				num++;
			}
	}
	//report(cout);
	if (num >= size) return;
	cerr << "Warning: random set has less than " << size << "taxa." << endl;
}


bool Split::overlap(Split &sp) {
	assert(ntaxa == sp.ntaxa);
	iterator it, it2;
	for (it = begin(), it2 = sp.begin(); it != end(); it++, it2++)
		if ((*it) & (*it2)) return true;
	return false;
	
}


Split::~Split()
{}

bool Split::containAny(IntVector &tax_id) {
	for (IntVector::iterator it = tax_id.begin(); it != tax_id.end(); it++)
		if (containTaxon(*it)) return true;
	return false;
}
