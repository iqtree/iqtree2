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

#ifndef _LP_WRAPPER
#define _LP_WRAPPER

#define tolerance 0.000001


#ifdef __cplusplus
extern "C" {
#endif

/**
	interface to call LP_SOLVE
	@param filename name of input lp file
	@param ntaxa number of taxa
	@param score (OUT) returned optimal score
	@param variables (OUT) array of returned solution
	@param verbose_mode verbose mode
	@return 
		0 if everything works file, 
		5 if solution is not optimal, 
		6 if some variable has wrong name, 
		7 if returned solution is not binary. In this case, one should run the solver 
		again with strict binary variable constraint.
*/
int lp_solve(char *filename, int ntaxa, double *score, double *variables, int verbose_mode);

/*int lp_demo();*/

void lp_solve_version_info(int *majorversion, int *minorversion, int *release, int *build);

#ifdef __cplusplus
}
#endif


#endif
