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

#include "Calibration.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

Calibration::Calibration(string calstr){
	
	isRootCal = false;
	nodeIDX = -1;
	prDistType = 1;
	exponRate = -1.0;
	exponMean = -1.0;
	stringstream ss;
	string tmp = "";
	ss << calstr;
	ss >> tmp;
	if(tmp[0] == '-'){
		if(tmp[1] == 'U' || tmp[1] == 'u')
			prDistType = 1;
		else if(tmp[1] == 'E' || tmp[1] == 'e')
			prDistType = 2;
		else{
			cerr << "ERROR: There's a problem with the calibration file " << endl;
			exit(1);
		}
		
		ss >> txn1;
		if(txn1 != "root")
			ss >> txn2;
		else { 
			isRootCal = true;
			txn2 = "root";
		}
		ss >> tmp;
		if(prDistType == 1){
			youngtime = atof(tmp.c_str());
			ss >> tmp;
			oldtime = atof(tmp.c_str());
			cout << "   Uniform calibration on MRCA[" << txn1 << ", " << txn2 << "] --> ("
				<< youngtime << ", " << oldtime << ")" << endl;
		}
		else if(prDistType == 2){
			youngtime = atof(tmp.c_str());
			ss >> tmp;
			if(tmp[0] == '-'){
				if(tmp == "-r"){
					ss >> tmp;
					exponRate = atof(tmp.c_str());
					exponMean = 1.0 / exponRate;
				}
				else if(tmp == "-m"){
					ss >> tmp;
					exponMean = atof(tmp.c_str()) - youngtime;
					exponRate = 1.0 / exponMean;
				}
			}
			else{
				if(isRootCal){
					exponRate = 1.0 / (youngtime * 0.4);
					exponMean = 1.0 / exponRate;
				}
				else{
					exponRate = 1.0 / (youngtime * 0.25);
					exponMean = 1.0 / exponRate;
				}
			}
			cout << "   Offset-exponential calibration on MRCA[" << txn1 << ", " << txn2 << "] --> (offset="
			<< youngtime << ", lambda=" << exponRate << ", mean=" << exponMean + youngtime << ")" << endl;
		}
	}
	else{
		txn1 = tmp;
		isRootCal = false;
		if(txn1 != "root")
			ss >> txn2;
		else { 
			isRootCal = true;
			txn2 = "root";
		}
		ss >> tmp;
		youngtime = atof(tmp.c_str());
		ss >> tmp;
		oldtime = atof(tmp.c_str());
		cout << "   Uniform calibration on MRCA[" << txn1 << ", " << txn2 << "] --> ("
		<< youngtime << ", " << oldtime << ")" << endl;
	}
}

