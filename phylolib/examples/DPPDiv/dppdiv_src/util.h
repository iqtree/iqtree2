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

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

inline std::string getLineFromFile(std::string fileName, int lineNum) {

	std::ifstream fileStream(fileName.c_str());
	if (!fileStream){
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
	}
	std::string linestring = "";
	int line = 0;
	while( getline(fileStream, linestring).good() ){
		line++;
		if (line == lineNum)
			break;
	}
	fileStream.close();
	if (line != lineNum){
		std::cerr << "The file \"" + fileName + "\" has " << line << " lines. Could not find line " << lineNum << std::endl;
		exit(1);
	}
	return linestring;
}

inline double expNumTables(double a, int n) {
	
	double expectedNum = 0.0;
	for (int i=1; i<=n; i++)
		expectedNum += ( 1.0 / (i - 1.0 + a) );
	expectedNum *= a;
	return expectedNum;
}

inline double calculateFromPriorMean(double target, int n) {
	
	double a = 0.000001;
	double ea = expNumTables(a, n);
	bool goUp;
	if (target <= 1.0) 
		target = 1.01;
	if (ea < target)
		goUp = true;
	else
		goUp = false;
	double increment = 0.1;
	while ( fabs(ea - target) > 0.000001 )
	{
		if (ea < target && goUp == true)
		{
			a += increment;
		}
		else if (ea > target && goUp == false)
		{
			a -= increment;
		}
		else if (ea < target && goUp == false)
		{
			increment /= 2.0;
			goUp = true;
			a += increment;
		}
		else
		{
			increment /= 2.0;
			goUp = false;
			a -= increment;
		}
		ea = expNumTables(a, n);
	}
	return a;
}

inline void normalizeVector(std::vector<double> &v) {
	
	int n = (int)v.size();
	double lnC = v[0];
	for (int i=1; i<n; i++)
	{
		if (v[i] > lnC)
			lnC = v[i];
	}
	
	for (int i=0; i<n; i++)
		v[i] -= lnC;
	
	double sum = 0.0;
	for (int i=0; i<n; i++)
	{
		if ( v[i] < -300.0 )
			v[i] = 0.0;
		else
			v[i] = exp( v[i] );
		sum += v[i];
	}
	
	for (int i=0; i<n; i++)
		v[i] /= sum;
}

#endif
