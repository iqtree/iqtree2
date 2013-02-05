/* 
 * DPPDiv version 1.0b source code (git: 9c0ac3d2258f89827cfe9ba2b5038f0f656b82c1)
 * Copyright 2009-2011
 * Tracy Heath(1,2,3) (NSF postdoctoral fellowship in biological informatics DBI-0805631)
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck
 */
 
#include "Alignment.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <vector>
#include <cstdlib>

using namespace std;

Alignment::Alignment(string fn) {
	
	numTaxa = numChar = numPatterns = 0;
	matrix = compressedMatrix = NULL;
	patternCount = NULL;
	isCompressed = false;

	ifstream seqStream(fn.c_str());
	if (!seqStream) 
		{
		cerr << "Cannot open file \"" + fn + "\"\n";
		exit(1);
		}

	string lineString = "";
	int lineNum = 0;
	int taxonNum = 0;
	while ( getline(seqStream, lineString).good() )
		{
		istringstream linestream(lineString);
		int ch;
		string word = "";
		int wordNum = 0;
		int siteNum = 0;
		string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			if (lineNum == 0)
				{
				int x;
				istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numChar = x;
				if (numTaxa > 0 && numChar > 0 && matrix == NULL)
					{	
					cout << "   Number of taxa  = " << numTaxa << '\n';
					cout << "   Number of sites = " << numChar << '\n';
					cout.flush();
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa*numChar];
					for (int i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numChar;
					}
				}
			else
				{
				if (wordNum == 1)
					{
					taxonNames.push_back(word);
					taxonNum++;
					}
				else
					{
					for (unsigned i=0; i<word.length(); i++)
						{
						char site = word.at(i);
						matrix[taxonNum-1][siteNum++] = nucID(site);
						}
					}
				}
			} while ( (ch=linestream.get()) != EOF );		
		lineNum++;
		}
	seqStream.close();
}

Alignment::~Alignment(void) {

	if (matrix != NULL)
		{
		delete [] matrix[0];
		delete [] matrix;
		}
	if (compressedMatrix != NULL)
		{
		delete [] compressedMatrix[0];
		delete [] compressedMatrix;
		}
	if (patternCount != NULL)
		delete [] patternCount;
}

void Alignment::compress(void) {

	if (isCompressed == false)
		{
		int *tempCnt = new int[numChar];
		for (int i=0; i<numChar; i++)
			tempCnt[i] = 1;
		for (int i=0; i<numChar; i++)
			{
			if (tempCnt[i] > 0)
				{
				numPatterns++;
				for (int j=i+1; j<numChar; j++)
					{
					if (tempCnt[j] > 0)
						{
						bool isSame = true;
						for (int k=0; k<numTaxa; k++)
							{
							if (matrix[k][i] != matrix[k][j])
								{
								isSame = false;
								break;
								}
							}
						if (isSame == true)
							{
							tempCnt[i]++;
							tempCnt[j] = 0;
							}
						}
					}
				}
			}
		
		compressedMatrix = new int*[numTaxa];
		compressedMatrix[0] = new int[numTaxa * numPatterns];
		for (int i=1; i<numTaxa; i++)
			compressedMatrix[i] = compressedMatrix[i-1] + numPatterns;
		patternCount = new int[numPatterns];
		
		for (int i=0, j=0; i<numChar; i++)
			{
			if (tempCnt[i] > 0)
				{
				for (int k=0; k<numTaxa; k++)
					compressedMatrix[k][j] = matrix[k][i];
				patternCount[j] = tempCnt[i];
				j++;
				}
			}
		
		isCompressed = true;
		
		delete[] tempCnt;
		}
}

int Alignment::getIndexForTaxonNamed(string nm) {

	int idx = -1;
	for (unsigned i=0; i<taxonNames.size(); i++)
		{
		if ( nm == taxonNames[i] )
			{
			idx = i;
			break;
			}
		}
	return idx;
}

int Alignment::getNumChar(void) {

	if (isCompressed == false)
		return numChar;
	return numPatterns;
}

int Alignment::getNucleotide(int i, int j) {

	if (isCompressed == false)
		return matrix[i][j];
	return compressedMatrix[i][j];
}

int Alignment::getNumSitesOfPattern(int i) {

	if (isCompressed == false)
		return 1;
	return patternCount[i];
}

/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs (int nucCode, int *nuc) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
}

bool Alignment::isTaxonPresent(string nm) {

	bool isPresent = false;
	for (unsigned i=0; i<taxonNames.size(); i++)
		{
		if ( nm == taxonNames[i] )
			{
			isPresent = true;
			break;
			}
		}
	return isPresent;
}

/*-------------------------------------------------------------------
|
|   NucID: 
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1 
|        C            2     
|        G            4      
|        T U          8     
|        R            5      
|        Y           10       
|        M            3      
|        K           12   
|        S            6     
|        W            9      
|        H           11      
|        B           14     
|        V            7      
|        D           13  
|        N - ?       15       
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else
		return -1;
}

void Alignment::print(std::ostream & o) const {

	if (isCompressed == false)
		{
		for (int j=0; j<numChar; j++)
			{
			o << setw(10) << j+1 << " -- ";
			for (int i=0; i<numTaxa; i++)
				{
				int state = matrix[i][j];
				o << setw(5) << state << " ";
				}
			o << '\n';
			}
		}
	else
		{
		for (int j=0; j<numPatterns; j++)
			{
			o << setw(10) << j+1 << " -- " << setw(4) << patternCount[j] << " -- ";
			for (int i=0; i<numTaxa; i++)
				{
				int state = compressedMatrix[i][j];
				o << setw(5) << state << " ";
				}
			o << '\n';
			}
		}
	o.flush();
}
