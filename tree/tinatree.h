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
#ifndef TINATREE_H
#define TINATREE_H

#include "phylotree.h"

/**
	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class TinaTree : public PhyloTree
{
public:
    TinaTree();
    /**
     * Constructor with given alignment
     * @param alignment
     */
    TinaTree(Alignment *alignment);

    ~TinaTree();
    /**
            SLOW VERSION: compute the parsimony score of the tree, given the alignment
            @return the parsimony score
     */
    int computeParsimonyScore();

    /**
            SLOW VERSION: compute the parsimony score of the tree, given the alignment
            @return the parsimony score
            @param node the current node
            @param dad dad of the node, used to direct the search
            @param ptn pattern ID
            @param states set of admissible states at the current node (in binary code)
     */
    int computeParsimonyScore(int ptn, int &states, PhyloNode *node = NULL, PhyloNode *dad = NULL);

	virtual void initializeAllPartialLh();

	virtual void initializeAllPartialLh(int &index, int &indexlh,
                                        PhyloNode *node = NULL, PhyloNode *dad = NULL);

};

#endif
