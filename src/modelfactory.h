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
#ifndef MODELFACTORY_H
#define MODELFACTORY_H

#include "tools.h"

/**
Store the transition matrix corresponding to evolutionary time so that one must not compute again. 
For efficiency purpose esp. for protein (20x20) or codon (61x61)

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelFactory : public hash_map<int, double*>
{
public:


	/**
		get the transition matrix
		@param time branch length
		@return transition matrix for time
	*/
	double *getTransMatrix(double time);

	double *getTransMatrixDerv(double time, double *trans_derv1, double trans_derv2);

	void addTransMatrix(double time, double *trans_mat);

    ~ModelFactory();

};

#endif
