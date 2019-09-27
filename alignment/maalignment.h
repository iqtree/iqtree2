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
#ifndef MAALIGNMENT_H
#define MAALIGNMENT_H

#include "alignment.h"

typedef vector< Alignment > AlignmentVector;

/**
Extended Alignment class to serve some analysis, created by MA

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class MaAlignment : public Alignment
{
public:
    MaAlignment() : Alignment() {};

    MaAlignment(char *filename,  char *sequence_type, InputType &intype, string model) : Alignment(filename, sequence_type, intype, model){};
	
	MaAlignment(Alignment &align) : Alignment(align){};

	/**
		To generate a new alignment from a given alignment (with the Expected Normalized Frequency)
		@param inputAlign the input alignment for which we can derive the expected normalized requency (CONSTANT)
		@param prop (OUT) the probability of the new alignment given the observed frequency of patterns in the input alignment (inputAlign).
		THEN THIS ALIGNMENT IS UPDATED!
		prop is computed as follows:
		- We have pattern 1 ... k in the inputAlign with observed freq. d_1 ... d_k (d_1+..+d_k = ell)
		==> The observed (relative) frequencies are p_1 ... p_k, p_i = d_i/ell
		- From some tree T we know the likelihood of each pattern given the tree and we derive the expected frequency (the expected alignment) d1(T) ... dk(T) where sum d_i(T) = ell
		===> prop = [ell!/product(d_i(T)!)] * product(p_i^d_i(T)).
		
	*/
	void generateExpectedAlignment(MaAlignment *inputAlign, double &prop);

	/**
		To generate a new alignment with the Expected Normalized Frequency
	*/
	//void generateExpectedAlignment(Alignment &returnAlign, const IntVector expectedNorFre);

	/**
		To print a list containing: patterns and the corresponding observed, expected frequencies
		@param fileName a file to store the information		
	*/
	void printPatObsExpFre(const char *fileName);

	/**
		To print a list containing: patterns and the corresponding observed, expected frequencies
		@param fileName a file to store the information
		@param expectedNorFre a vector containing the expected frequencies
	*/
	void printPatObsExpFre(const char *fileName, const IntVector expectedNorFre);
	/**		
		To compute the Expected Normalized Frequencies of the patterns in the alignment.
		The values in  this vector should be in the same order as the patterns in the pattern vector.
		These values are computed based on
			+ the length of the alignment (ell)
			+ the logLL vector
		How?
			(do not need but to be clear: observed frequencies d1, d2, ..., dk)
			logLL --> likelihood: p1, p2 ... pk
			Because we have in total 4^n patterns but may observe k < 4^n patterns ==> p1 + p2 + ... + pk <= 1.
			This also means (p1 + p2 + ... + pk)*ell <= ell
			Now we want to derive expected frequencies ^d1, ^d2, ..., ^dk such that ^d1 + ^d2 + ... + ^dk = ell based on p1, p2, ..., pk.
			We do the followings:
			+ Compute li = ell*pi / sum_i(pi)  ==> sum_i (li) = ell
			+ Because li is usually not an integer, we now have to round li as below:
				* r1 = l1
				* r_{i+1} = l_{i+1} + r_{i} - [r_i]
				* Finally set: ^d_i = [r_i]
				* where [.] denotes ordinary rounding
	*/
	IntVector computeExpectedNorFre();
	/**
		To read the log likelihood of the patterns from a file
		@param filename file contains the site and log likelihood
		@result: the vector logLL will be changed by this function
	*/
	void readLogLL(char *filename);

	/**
		Compute the multinomial probabilities for a vector of object alignments according to the parameters determined by THIS alignment
		@param objectAligns vector containing the object alignments
		@param probs (OUT) returned vector containing the probabilities (double)
	*/
	//void multinomialProb (AlignmentVector objectAligns, DoubleVector &probs);
private:
	/*
		Log likelihood of the patterns.
	 	The values in this vector should be in the same order as the patterns in the pattern vector. 
		This must be made sure while reading the file containing these log likelihood values.
	*/
	DoubleVector logLL;	
	//IntVector expectedNorFre;
};

#endif
