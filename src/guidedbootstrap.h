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

#ifndef GUIDED_BOOTSTRAP_H
#define GUIDED_BOOTSTRAP_H

#include "tools.h"
#include "iqptree.h"
#include "alignment.h"

#if defined(WIN32)
namespace stdext {
#else
namespace __gnu_cxx {
#endif

	template<>
	struct hash<IntVector*> {
		size_t operator()(const IntVector* sp) const {
			size_t sum = 0;
			for (IntVector::const_iterator it = sp->begin(); it != sp->end(); it++)
				sum = (*it) + (sum << 6) + (sum << 16) - sum;
			return sum;
		}
	};
} // namespace __gnu_cxx

typedef hash_map<IntVector*, int> IntVectorMap;

typedef vector<IntVector*> IntVectorCollection;

/**
	run guided bootstrap
*/
void runGuidedBootstrap(Params &params, Alignment *alignment, IQPTree &tree);

void runAvHTest(Params &params, Alignment *alignment, IQPTree &tree);

#endif