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
#ifndef SPLITSET_H
#define SPLITSET_H

#include <vector>
#include "split.h"

/**
Vector of Splits

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class SplitSet : public vector<Split*>
{
public:
    SplitSet();
	
	/**
		release the memory of all element splits, then resize to 0
	*/
	void removeAll();
	
	/**
		get the weight of the first split in the vector
	*/
	double getWeight();

	/**
 		check the compatibility of sp against all splits in this set
		@param sp the target split
		@return TRUE if sp is compatible with all splits here, otherwise FALSE
	*/
	bool compatible(Split *sp);

    virtual ~SplitSet();

};

#endif
