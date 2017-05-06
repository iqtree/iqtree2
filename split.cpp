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
	/* if taxa_list contains exactly one haft of taxa, only turn on the inverted mode
	 if taxon 0 is not in the list */
/*	bool inverted = (taxa_list.size()*2 > ntaxa);
	if (taxa_list.size()*2 == ntaxa) {
		inverted = true;
		for (it = taxa_list.begin(); it != taxa_list.end(); it++)
			if ((*it) == 0) {
				inverted = false;
				break;
			}
	}*/

	// resize the split size
	resize((ntaxa + UINT_BITS -1) / UINT_BITS, 0);

	for (it = taxa_list.begin(); it != taxa_list.end(); it++)
	{
		int value = *it;
		int bit_pos = value / UINT_BITS;
		int bit_off = value % UINT_BITS;
		(*this)[bit_pos] |= (UINT) (1 << bit_off);
	}


	//if (inverted) invert();
	if (shouldInvert()) invert();

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
				//out << i * UINT_BITS + j + 1 << " ";
				out << i * UINT_BITS + j << " ";
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
	ASSERT(sp.size() == size() && sp.ntaxa == ntaxa);

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
	ASSERT(taxa_set.size() == size() && taxa_set.ntaxa == ntaxa);

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
		if ((it + 1) == end()) {
			max_step = ntaxa % UINT_BITS;
			if (!max_step) max_step = UINT_BITS;
		}
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
	ASSERT(tax_id >= 0 && tax_id < ntaxa);
	int pos = tax_id / UINT_BITS, off = tax_id % UINT_BITS;
	(*this)[pos] |= 1 << off;
}

/**
	remove a taxon from the split
	@param tax_id id of taxon from 0..ntaxa-1
*/
void Split::removeTaxon(int tax_id)
{
	ASSERT(tax_id >= 0 && tax_id < ntaxa);
	int pos = tax_id / UINT_BITS, off = tax_id % UINT_BITS;

	(*this)[pos] &= -1 - (1 << off);
}

/**
	@param tax_id id of taxon from 0..ntaxa-1
	@return true if tax_id is in the set
*/
bool Split::containTaxon(int tax_id)
{
	ASSERT(tax_id >= 0 && tax_id < ntaxa);
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

bool Split::operator<(const Split &sp) const {
	return countTaxa() < sp.countTaxa();
}

Split &Split::operator+=(Split &sp) {
	ASSERT(sp.ntaxa == ntaxa);
	iterator it1, it2;
	for (it1 = begin(), it2 = sp.begin(); it1 != end(); it1++, it2++) {
		(*it1) |= (*it2);
	}
	return *this;
}

Split &Split::operator*=(Split &sp) {
	ASSERT(sp.ntaxa == ntaxa);
	iterator it1, it2;
	for (it1 = begin(), it2 = sp.begin(); it1 != end(); it1++, it2++) {
		(*it1) &= (*it2);
	}
	return *this;	
}

Split &Split::operator-=(Split &sp) {
	ASSERT(sp.ntaxa == ntaxa);
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
	ASSERT(ntaxa == sp.ntaxa);
	for (iterator it = begin(), it2 = sp.begin(); it != end(); it++, it2++)
		if ( ((*it) & (*it2)) != (*it) )
			return false;
	return true;
}

Split &Split::operator= (const Split &sp) {
	ASSERT(ntaxa == sp.ntaxa);
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
	ASSERT(size < ntaxa);
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
			if (!containTaxon(cnt) && ( random_double() <= prob )) {
				addTaxon(cnt);
				num++;
			}
	}
	//report(cout);
	if (num >= size) return;
	cerr << "WARNING: random set has less than " << size << "taxa." << endl;
}


bool Split::overlap(Split &sp) {
	ASSERT(ntaxa == sp.ntaxa);
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

Split *Split::extractSubSplit(Split &taxa_mask) {
	ASSERT(taxa_mask.getNTaxa() == getNTaxa());
	Split *sp = new Split(taxa_mask.countTaxa());
	int id = 0;
	for (int tax = 0; tax < ntaxa; tax++)
	if (taxa_mask.containTaxon(tax)) {
		if (containTaxon(tax))
			sp->addTaxon(id);
		id++;
	}
	ASSERT(id == sp->getNTaxa());
	return sp;
}


/**
	Solve k-means problem for one-dimensional data with dynamic programming
	@param n number of data points
	@param ncat number of clusters
	@param data data point of size n: x[0..n-1]
	@param center (OUT) output k centers of k clusters: center[0...k-1] will be filled
	@param cluster (OUT) cluster assignments for each data point: cluster[0...n-1] will be filled
	@return the minimum sum of squares over all k clusters
*/
double kMeansOneDim(int n, int ncat, double *data, double *center, int *cluster) {
	int i, j, m, k = ncat;
	if (ncat == 0) k = n;
	/**
		dynamic programming cost matrix, c[i][j] = cost of i clusters for {x1...xj}
	*/
	double **c = (double**) new double[k]; 
	/**
		id is used to trace back the minimal solution
	*/
	double **id = (double**) new double[k]; 
	/** 
		c1[i][j] = cost of 1 cluster for {xi...xj}
	*/
	double **c1 = (double**) new double[n];
	/** 
		m1[i][j] = mean of {xi...xj}
	*/
	double **m1 = (double**) new double[n];
	
	double *x = new double[n]; // sorted data points

	double *h = new double[n]; // Hartigan index

	// allocate memory 
	for (i = 0; i < k; i++) c[i] = new double[n];
	for (i = 0; i < k; i++) id[i] = new double[n];
	for (i = 0; i < n; i++) c1[i] = new double[n];
	for (i = 0; i < n; i++) m1[i] = new double[n];

	// first sort data into x
	memmove(x, data, sizeof(double)*n);
	std::sort(x, x+n);
	// first compute c1 matrix
	for (i = 0; i < n; i++) {
		double sum = 0.0;
		for (j = i; j < n; j++) {
			sum += x[j];
			double mean = sum / (j-i+1);
			m1[i][j] = mean;
			double ss = 0; 
			for (m = i; m <= j; m++) 
				ss += (x[m]-mean)*(x[m]-mean); // sum of squared difference
				//ss += fabs(x[m]-mean); // sum of absolute difference
			c1[i][j] = ss;
		}
	}

	/* now compute dynamic programming matrix */
	// initialize the 1st row
	for (j = 0; j < n; j++) {
		c[0][j] = c1[0][j];
		id[0][j] = -1;
	}
	for (i = 1; i < k; i++) {
		// no i clusters exist for less than i data points
		for (j = 0; j < i; j++) { c[i][j] = INFINITY; id[i][j] = -1; }
		for (j = i; j < n; j++) {
			c[i][j] = INFINITY;
			for (m = i-1; m < j; m++)
				if (c[i][j] > c[i-1][m] + c1[m+1][j]) {
					c[i][j] = c[i-1][m] + c1[m+1][j];
					id[i][j] = m;
				}
		}
		// compute Hartigan index
		h[i-1] = (n-i-1)*(c[i-1][n-1]-c[i][n-1]) / c[i][n-1];
		//cout << i << " clusters " << h[i-1] << endl;
	}

	double min_cost = c[k-1][n-1];
	int *bound = new int[k+1];
	// now trace back
	bound[k] = n-1;
	for (i = k-1; i >= 0; i--) {
		bound[i] = id[i][bound[i+1]];
	}

	for (i = 0; i < k; i++) {
		center[i] = m1[bound[i]+1][bound[i+1]];
		for (j = 0; j < n; j++)
			if (data[j] <= x[bound[i+1]] && data[j] >= x[bound[i]+1])
				cluster[j] = i;
	}

	// free memory
	delete [] bound;
	for (i = n-1; i >= 0; i--) delete [] m1[i];
	for (i = n-1; i >= 0; i--) delete [] c1[i];
	for (i = k-1; i >= 0; i--) delete [] id[i];
	for (i = k-1; i >= 0; i--) delete [] c[i];

	delete [] h;
	delete [] x;
	delete [] m1;
	delete [] c1;
	delete [] id;
	delete [] c;
		
	return min_cost;
}
