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
#include "modeluser.h"

ModelUser::ModelUser(PhyloTree *tree)
 : GTRModel(tree)
{
}

void ModelUser::readRates(istream &in) throw(const char*) {
	int nrates = num_states*(num_states-1)/2;
	for (int i = 0; i < nrates; i++) {
		if (!(in >> rates[i]))
			throw "Rate entries could not be read";
		if (rates[i] < 0.0)
			throw "Negative rates found";
	}
}

void ModelUser::readStateFreq(istream &in) throw(const char*) {
	int i;
	for (i = 0; i < num_states; i++) {
		if (!(in >> state_freq[i])) 
			throw "State frequencies could not be read";
		if (state_freq[i] < 0.0)
			throw "Negative state frequencies found";
	}
	double sum = 0.0;
	for (i = 0; i < num_states; i++) sum += state_freq[i];
	if (fabs(sum-1.0) > 1e-5)
		throw "State frequencies do not sum up to 1.0";
}

int ModelUser::getNDim() { 
	return 0;
	assert(freq_type != FREQ_UNKNOWN);
	int ndim = 0; 
	if (freq_type == FREQ_ESTIMATE) 
		ndim += num_states-1;
	return ndim;
}
