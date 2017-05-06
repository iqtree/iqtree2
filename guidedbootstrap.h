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
#include "iqtree.h"
#include "alignment.h"

struct hashfunc_IntVector {
	size_t operator()(const IntVector* sp) const {
		size_t sum = 0;
		for (IntVector::const_iterator it = sp->begin(); it != sp->end(); it++)
			sum = (*it) + (sum << 6) + (sum << 16) - sum;
		return sum;
	}
};

namespace std {
	/**
		Define equal_to of two IntVector, used for hash_set (or hash_map) template
	*/
	template<>
	struct equal_to<IntVector*> {
		/**
			@return true if *s1 == *s2
			@param s1 first IntVector
			@param s2 second IntVector
		*/
		bool operator()(const IntVector* s1, const IntVector* s2) const{
			return *s1 == *s2;
		}
	};
	/**
		Define less than relationship of two IntVector, used for set (or map) template
	*/
	template<>
	struct less<IntVector*> {
		/**
			@return true if *s1 < *s2 alphabetically
			@param s1 first IntVector
			@param s2 second IntVector
		*/
		bool operator()(const IntVector *s1, const IntVector *s2) const {
			ASSERT(s1->size() == s2->size());
			for (int i = 0; i < s1->size(); i++)
				if ((*s1)[i] < (*s2)[i]) 
					return true;
				else if ((*s1)[i] > (*s2)[i]) return false;
			return false;
		}
	};
} // namespace std


#ifdef USE_HASH_MAP
typedef unordered_map<IntVector*, int, hashfunc_IntVector> IntVectorMap;
#else
typedef map<IntVector*, int> IntVectorMap;
#endif

typedef vector<IntVector*> IntVectorCollection;

/**
	OBSOLETE: run guided bootstrap (this function was only used at the beginning of the UFBoot project
*/
void runGuidedBootstrap(Params &params, Alignment *alignment, IQTree &tree);

void runAvHTest(Params &params, Alignment *alignment, IQTree &tree);

/**
 * make the plot with x-axis being the alignments and y-axis being the likelihood of all trees
 * Right now we use Kullback-Leibler distance to arrange alignments on the x-axis
 */
void runBootLhTest(Params &params, Alignment *alignment, IQTree &tree);

#endif
