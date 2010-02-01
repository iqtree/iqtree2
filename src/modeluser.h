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
#ifndef MODELUSER_H
#define MODELUSER_H

#include "gtrmodel.h"
#include <iostream>

/**
User-defined models: The rates will be fixed.

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelUser : public GTRModel
{
public:
	/**
		construction
		@param tree associated phylogenetic tree
	*/
    ModelUser(PhyloTree *tree);

	/**
		read the rates from an input stream. it will throw error messages if failed
		@param in input stream
	*/
	virtual void readRates(istream &in) throw(const char*);
	
	/**
		read state frequencies from an input stream. it will throw error messages if failed
		@param in input stream
	*/
	virtual void readStateFreq(istream &in) throw(const char*);

	/**
		@return the number of dimensions
	*/
	virtual int getNDim();
	
};

#endif
