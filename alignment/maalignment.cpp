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
#include "maalignment.h"

void MaAlignment::readLogLL(char *fileName)
{
	//First read the values from inFile to a DoubleVector
	DoubleVector _logllVec;
	int siteNum = -1;
	string currentString;
	cout << "\nReading file containing site's loglikelihood: " << fileName << "...." << endl;
    ifstream inFile;
	try{
		inFile.exceptions (ios::failbit | ios::badbit);
		inFile.open(fileName);
		/**really start reading*/
		//read number of sites
		inFile >> currentString;
		siteNum = convert_int(currentString.c_str());
		//ignore "Site_Lh"		
		inFile >> currentString;		
		while (!inFile.eof())
		{
			//reading each line of the file
			//remove the badbit
			inFile.exceptions (ios::badbit);
			if ( !(inFile >> currentString) ) break;
			//set the failbit again
			inFile.exceptions (ios::failbit | ios::badbit);
			_logllVec.push_back(convert_double(currentString.c_str()));
		}/**finish reading*/
		inFile.clear();
		inFile.exceptions (ios::failbit | ios::badbit);
		inFile.close();
	} catch(bad_alloc){
			outError(ERR_NO_MEMORY);
	} catch (const char *str){
			outError(str);
	} catch (string str){
			outError(str);
	} catch (ios::failure){
			outError(ERR_READ_INPUT);
	} catch (...){
			outError(ERR_READ_ANY);
	}
	if (siteNum != _logllVec.size())
		outError("Actual number of site's likelihoods is not consistent with the announced number in the first line.");
	cout << "Finish reading, now assign the logLL to the pattern:" << endl;

	logLL.resize(getNPattern(),0.0);
	for (int i = 0; i < siteNum; i++)
	{
		int patIndex = getPatternID(i);
		if ( logLL[patIndex] == 0 )
			logLL[patIndex] = _logllVec[i];
		else
			if ( logLL[patIndex] != _logllVec[i] )
//				outError("Conflicting between the likelihoods reported for pattern", (*this)[i]);
                outError("Conflicting between the likelihoods reported for pattern");
	}
//	int npat = getNPattern();
//	cout << "Number of patterns: " << npat << endl;
//	for ( int j = 0; j < npat; j++ )
//		cout << j << "\t" << at(j) << "\t" << logLL[j] << endl;
	cout << "Finish assigning logLL to the patterns!" << endl;	 
}

IntVector MaAlignment::computeExpectedNorFre()
{
	IntVector expectedNorFre;
	if ( logLL.empty()) 
		outError("Error: log likelihood of patterns are not given!");

	size_t patNum = getNPattern();
	size_t alignLen = getNSite();		
	//resize the expectedNorFre vector
	expectedNorFre.resize(patNum,-1);

	//Vector containing the likelihood of the pattern p_i
	DoubleVector LL(patNum,-1.0);
	double sumLL = 0; //sum of the likelihood of the patterns in the alignment

	//Compute the likelihood from the logLL
	for ( int i = 0; i < patNum; i++ )
	{
		LL[i] = exp(logLL[i]);
		sumLL += LL[i];
	}

	//Vector containing l_i = p_i*ell/sum_i(p_i)
	DoubleVector ell(patNum, -1.0);
	//Compute l_i
	for ( int i = 0; i < patNum; i++ )
	{
		ell[i] = (double)alignLen * LL[i] / sumLL;
	}


	//Vector containing r_i where r_0 = ell_0; r_{i+1} = ell_{i+1} + r_i - ordinaryRounding(r_i)
	DoubleVector r(patNum, -1.0);
	//Compute r_i and the expected normalized frequencies
	r[0] = ell[0];
	expectedNorFre[0] = (int)floor(ell[0]+0.5); //note that floor(_number+0.5) returns the ordinary rounding of _number
	int sum = expectedNorFre[0];
	for (int j = 1; j < patNum; j++ )
	{
		r[j] = ell[j] + r[j-1] - floor(r[j-1]+0.5);
		expectedNorFre[j] = (int)floor(r[j]+0.5);
		sum += expectedNorFre[j];
	}
	
	//cout << "Number of patterns: " << patNum << ", sum of expected sites: " << sum << endl;
	return expectedNorFre;
}

void MaAlignment::printPatObsExpFre(const char *fileName)
{
	IntVector expectedNorFre = computeExpectedNorFre();
	printPatObsExpFre(fileName, expectedNorFre);
}

void MaAlignment::printPatObsExpFre(const char *fileName, const IntVector expectedNorFre)
{	
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(fileName);
		out << "Pattern\tLogLL\tObservedFre\tExpectedFre" << endl;

		size_t patNum = getNPattern();
		size_t seqNum = getNSeq();

		for ( size_t i = 0; i < patNum; ++i )
		{
			for ( size_t seqID = 0; seqID < seqNum; ++seqID ){
				out << convertStateBackStr(at(i)[seqID]);
			}
			out << "\t" << logLL[i] << "\t" << (*this)[i].frequency << "\t" << expectedNorFre[i] << endl;
		}
		out.close();
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, fileName);
	}
}

void MaAlignment::generateExpectedAlignment(MaAlignment *aln, double &prob)
{
	//cout << "In function: generating expected alignment!" << endl;
	IntVector expectedNorFre = aln->computeExpectedNorFre();
	
	int nsite = aln->getNSite();
	seq_names.insert(seq_names.begin(), aln->seq_names.begin(), aln->seq_names.end());
	num_states = aln->num_states;
	site_pattern.resize(nsite, -1);
	clear();
	pattern_index.clear();
	VerboseMode save_mode = verbose_mode; 
	verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern

	int site = 0;
	int npat = aln->getNPattern();

	double sumFac = 0;
	double sumProb = 0;
	double fac = logFac(nsite);

	double sumFacMax = 0;
	double sumProbMax = 0;

	for (int patID = 0; patID < npat; ++patID) {
		int patFre = expectedNorFre[patID];
		for ( int patSite = 0; patSite < patFre; patSite++)
		{			
			Pattern pat = aln->at(patID);
			addPattern(pat,site);
			site++;	
		}

		//to compute the probability of the new alignment given the multinomial distribution
		sumFac += logFac(patFre);
		sumProb += (double)patFre*log((double)aln->at(patID).frequency/(double)nsite);

		//for the unconstraint maximum log likelihood
		sumFacMax += logFac(aln->at(patID).frequency);
		sumProbMax += (double)aln->at(patID).frequency*log((double)aln->at(patID).frequency/(double)nsite);
	}
	prob = fac - sumFac + sumProb;

	double probMax = fac - sumFacMax + sumProbMax;
//	cout << "total number of sites: " << site << endl;
	verbose_mode = save_mode;
	countConstSite();
	//cout << "Finish generating expected alignment!" << endl;
	cout << "Logarithm of the probability of the new alignment given the multinomial distribution of the input alignment is: " << prob << endl;
	cout << "Maximum unconstraint (log) likelihood of the input alignment: " << probMax << endl;
// 	cout << "Maximum unconstraint likelihood: " << exp(probMax) << endl;
}

/*void MaAlignment::multinomialProb(Alignment objectAlign, double &prob)
{
	cout << "Computing the multinomial probability of an object alignment given a reference alignment ..." << endl;
	//should we check for compatibility of sequence's names and sequence's order in THIS alignment and in the objectAlign??
	//check alignment length
	int nsite = getNSite();
	assert(nsite == objectAlign.getNSite());
	double sumFac = 0;
	double sumProb = 0;
	double fac = logFac(nsite);
	int index;
	for ( Alignment::iterator objectIt = objectAlign.begin(); objectIt != objectAlign.end() ; objectIt++)
	{
		PatternIntMap::iterator pat_it = pattern_index.find((*objectIt));
		if ( pat_it == pattern_index.end() ) //not found ==> error
			outError("Pattern in the object alignment is not found in the reference alignment!");
		sumFac += logFac((*objectIt).frequency);
		index = pat_it->second;
		sumProb += (double)(*objectIt).frequency*log((double)at(index).frequency/(double)nsite);
	}
	prob = fac - sumFac + sumProb;
}*/

/*void MaAlignment::multinomialProb(AlignmentVector objectAligns, DoubleVector &probs)
{
	int num = objectAligns.size();
	double curProb;
	if (num > 0)
	{
		probs.resize(num,0);
		for ( int i = 0; i < num; i++ )
		{
			(*this).multinomialProb(objectAligns[i], curProb);
			probs[i] = curProb;
		}
	}
}*/
