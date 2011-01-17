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
#include "rategammainvar.h"

RateGammaInvar::RateGammaInvar(int ncat, double shape, double p_invar_sites, PhyloTree *tree)
: RateInvar(p_invar_sites, tree), RateGamma(ncat, shape, tree)
{
	name = "+I" + name;
	full_name = "Invar+" + full_name;
}


double RateGammaInvar::computeFunction(double value) {
	if (cur_optimize == 0)
		return RateGamma::computeFunction(value);
	else 
		return RateInvar::computeFunction(value);
	
}

double RateGammaInvar::optimizeParameters() {
	double tree_lh;
	cur_optimize = 0;
	tree_lh = RateGamma::optimizeParameters();
	cur_optimize = 1;
	tree_lh = RateInvar::optimizeParameters();
	return tree_lh;
}

void RateGammaInvar::writeInfo(ostream &out) {
	RateInvar::writeInfo(out);
	RateGamma::writeInfo(out);
}

void RateGammaInvar::writeParameters(ostream &out) {
	RateInvar::writeParameters(out);
	RateGamma::writeParameters(out);
}


