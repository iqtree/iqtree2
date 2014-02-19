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

void RateHeterogeneity::writeSiteRates(ostream &out, DoubleVector &pattern_rates, IntVector &pattern_cat) {
	int nsite = phylo_tree->aln->getNSite();
	int i;
	
	out.setf(ios::fixed,ios::floatfield);
	out.precision(5);
	out << "Site\tRate";
	if (!pattern_cat.empty()) out << "\tCategory\tCategorized_rate";
	out << endl;
	//cout << __func__ << endl;
	for (i = 0; i < nsite; i++) {
		int ptn = phylo_tree->aln->getPatternID(i);
		out << i+1 << "\t";
		if (pattern_rates[ptn] >= MAX_SITE_RATE) out << "100.0"; else out << pattern_rates[ptn];
		//cout << i << " "<< ptn << " " << pattern_cat[ptn] << endl;
		if (!pattern_cat.empty()) out << "\t" << pattern_cat[ptn]+1 << "\t" << getRate(pattern_cat[ptn]);
		out << endl;
	}
}

void RateHeterogeneity::writeSiteRates(ostream &out) {
	DoubleVector pattern_rates;
	IntVector pattern_cat;
	computePatternRates(pattern_rates, pattern_cat);
	if (pattern_rates.empty()) return;
	writeSiteRates(out, pattern_rates, pattern_cat);
}

void RateHeterogeneity::writeSiteRates(const char *file_name) {
	DoubleVector pattern_rates;
	IntVector pattern_cat;
	computePatternRates(pattern_rates, pattern_cat);
	if (pattern_rates.empty()) return;

	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name);
		writeSiteRates(out, pattern_rates, pattern_cat);
		out.close();
		cout << "Site rates printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}
}

double RateHeterogeneity::targetFunk(double x[]) {
	return -phylo_tree->computeLikelihood();
}
