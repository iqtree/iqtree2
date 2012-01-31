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
#ifndef PHYLOSUPERTREE_H
#define PHYLOSUPERTREE_H

#include <phylotree.h>

/**
Phylogenetic tree for partition model (multi-gene alignment)

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class PhyloSuperTree : public PhyloTree, public vector<PhyloTree* >
{
public:
	/** 
		constructor
	*/
    PhyloSuperTree();

	/** 
		constructor
	*/
    PhyloSuperTree(Params &params);


    virtual ~PhyloSuperTree();

protected:
	/**
		partition names
	*/
	StrVector part_names; 

};

#endif
