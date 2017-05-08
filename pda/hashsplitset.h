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

#ifndef HASHSPLITSET_H
#define HASHSPLITSET_H


#include "utils/tools.h"
#include "split.h"

using namespace std;

class SplitGraph;

#ifdef USE_HASH_MAP
/*
	Define the hash function of Split
*/
struct hashfunc_Split {
	size_t operator()(const Split* sp) const {
		size_t sum = 0;
		for (Split::const_iterator it = sp->begin(); it != sp->end(); it++)
			sum = (*it) + (sum << 6) + (sum << 16) - sum;
		return sum;
	}
};
#endif // USE_HASH_MAP

namespace std {
	/**
		Define equal_to of two splits, used for hash_set (or hash_map) template
	*/
	template<>
	struct equal_to<Split*> {
		/**
			@return true if *s1 == *s2
			@param s1 first split
			@param s2 second split
		*/
		bool operator()(const Split* s1, const Split* s2) const{
			return *s1 == *s2;
		}
	};
	/**
		Define less than relationship of two splits, used for set (or map) template
	*/
	template<>
	struct less<Split*> {
		/**
			@return true if *s1 < *s2 alphabetically
			@param s1 first split
			@param s2 second split
		*/
		bool operator()(const Split *s1, const Split *s2) const {
			ASSERT(s1->size() == s2->size());
			for (int i = 0; i < s1->size(); i++)
				if ((*s1)[i] < (*s2)[i]) 
					return true;
				else if ((*s1)[i] > (*s2)[i]) return false;
			return false;
		}
	};
} // namespace std


/**
SplitSet for quick search purpose

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
#ifdef USE_HASH_MAP
class SplitIntMap : public unordered_map<Split*, int, hashfunc_Split>
#else
class SplitIntMap : map<Split*, int, std::less<Split*> > 
#endif
{
public:

	/**
		find a split
		@param sp target split
		@return the split containing the same set of taxa with sp, NULL if not found
	*/
	Split *findSplit(Split *sp);

	/**
		find a split
		@param sp target split
		@param value (OUT) associated value
		@return the split containing the same set of taxa with sp, NULL if not found
	*/
	Split *findSplit(Split *sp, int &value);

	int getValue(Split *sp);

	void setValue(Split *sp, int value);

	void eraseSplit(Split *sp);

	void insertSplit(Split *sp, int value);

	/**
	 * build a map from the input split graph
	 * @param sg input split graph
	 * @param use_index TRUE to map to index of splits in sg, FALSE to map to split weights
	*/
	void buildMap(SplitGraph &sg, bool use_index = true);

	int getNumTree() {
		return numTree;
	}

	void setNumTree(int maxValue) {
		this->numTree = maxValue;
	}

private:

	/**
	 *  The maximum weight value. If the splits are generated from n trees and splits of every tree
	 *  all have weight = 1, then maxValue = n
	 *  This variable is used to determine whether a split appear on all input trees.
	 */
	int numTree;

};

#endif
