/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include "splitset.h"

SplitSet::SplitSet()
 : vector<Split*>()
{
}

double SplitSet::getWeight() {
	if (empty())
		return 0;
	//assert(!empty());
	return (*begin())->getWeight();
}


/**
	release the memory of all element splits, then resize to 0
*/
void SplitSet::removeAll() {
	for (reverse_iterator it = rbegin(); it != rend(); it++)
		if (*it) delete *it;
	clear();
}


bool SplitSet::compatible(Split *sp) {
	for (iterator it = begin(); it != end(); it++)
		if (!(*it)->compatible(*sp))
			return false;
	return true;
}


SplitSet::~SplitSet()
{
	removeAll();
	//cout << "deleted" << endl;
}


