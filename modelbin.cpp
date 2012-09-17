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
#include "modelbin.h"

ModelBIN::ModelBIN(const char *model_name, StateFreqType freq, PhyloTree *tree, bool count_rates)
: GTRModel(tree, count_rates) {
	assert(num_states == 2); // make sure that you create model for DNA
	StateFreqType def_freq = FREQ_UNKNOWN;
	name = model_name;
	full_name = model_name;

	//cout << "User-specified model "<< model_name << endl;
	readParameters(model_name);
	//name += " (user-defined)";
	if (freq == FREQ_UNKNOWN || def_freq == FREQ_EQUAL) freq = def_freq;
	GTRModel::init(freq);

}



