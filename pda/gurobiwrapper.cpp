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


#include <string.h>
#include <sstream>
#include "utils/tools.h"
#include "gurobiwrapper.h"

#define tolerance 0.000001

/**
	interface to call GUROBI LP solver
	@param filename name of input lp file
	@param ntaxa number of taxa
	@param score (OUT) returned optimal score
	@param variables (OUT) array of returned solution
	@param verbose_mode verbose mode
	@return 
		-1 if gurobi was not installed properly or does not exist at all
		0 if everything works file, 
		5 if solution is not optimal, 
		6 if some variable has wrong name, 
		7 if returned solution is not binary. In this case, one should run the solver 
		again with strict binary variable constraint.
*/
int gurobi_solve(char *filename, int ntaxa, double *score, double *variables, int verbose_mode, int num_threads) {
	int ret = 0;
	*score = -1;
	string command;
	ostringstream ss;

	ss << "gurobi_cl Threads=" << num_threads << " ResultFile=" << filename
		<< ".sol MIPGap=0 "<< filename  << " >" << filename << ".log ";
	command = ss.str();
	if (verbose_mode >= VB_MED)
		cout << command << endl;
	int sys_ret = system(command.c_str());
	if (sys_ret != 0) {
		cout << "gurobi_cl could not be executed. Make sure it was installed with proper license." << endl;
		cout << command << endl;
		return -1;
	}

	command = filename;
	command += ".sol";

	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(command.c_str());
		string str;

		while (!in.eof()) {
			// remove the failbit
			in.exceptions(ios::badbit);
			if(!(in >> str)) break;
			// set the failbit again
			in.exceptions(ios::failbit | ios::badbit);
			if (str[0] != 'x') continue;
			int index = convert_int(str.substr(1).c_str());
			if (index < 0 || index >= ntaxa) {
				cout << "Index x_" << index << " is not in the range!" << endl;
				ret = 6;
				break;
			}
			double value;
			in >> value;
			if (value > tolerance && (1.0 - value) > tolerance) {
				if (verbose_mode >= VB_MED) cout << endl << str << " = " << value;
				ret = 7;
				if (!verbose_mode) break;
			}
			variables[index] = value;
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (string str) {
		outError(str);
	}

	command = filename;
	command += ".log";
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(command.c_str());
		string str;

		while (!in.eof()) {
			in.exceptions(ios::badbit);
			if(!(in >> str)) break;
			// set the failbit again
			in.exceptions(ios::failbit | ios::badbit);
			if (str != "Best" && str != "Optimal") continue;
			in >> str;
			if (str != "objective") continue;
			in >> str;
			// remove the ending comma ,
			if (*str.rbegin() == ',') str.erase(str.length()-1);
			*score = convert_double(str.c_str());
			break;
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (string str) {
		outError(str);
	}
	return ret;
}
