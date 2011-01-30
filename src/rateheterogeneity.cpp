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

#include "phylotree.h"
#include "rateheterogeneity.h"


RateHeterogeneity::RateHeterogeneity()
 : Optimization()
{
	name = "";
	full_name = "Uniform";
	phylo_tree = NULL;
}

void RateHeterogeneity::setTree(PhyloTree *tree) {
	phylo_tree = tree;
}

RateHeterogeneity::~RateHeterogeneity()
{
}


void RateHeterogeneity::writeSiteRates(DoubleVector &pattern_rates, const char *file_name) {
	int nsite = phylo_tree->aln->getNSite();
	int i;
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name);
		for (i = 0; i < nsite; i++) 
			out << i+1 << "\t" << pattern_rates[phylo_tree->aln->getPatternID(i)] << endl;
		out.close();
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}
}
