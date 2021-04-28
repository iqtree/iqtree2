/***************************************************************************
 *   Copyright (C) 2018 by Lukasz Reszczynski                              *
 *   lukasz.reszczynski@univie.ac.at                                       *
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
#ifndef TERRACETPHAST_H
#define TERRACETPHAST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "tree/phylotree.h"
#include "alignment/superalignment.h"
#include "utils/tools.h"

#include "terraphast/include/terraces/errors.hpp"
#include "terraphast/include/terraces/parser.hpp"
#include "terraphast/include/terraces/simple.hpp"
#include "terraphast/include/terraces/advanced.hpp"
#include "terraphast/include/terraces/bitmatrix.hpp"

/**
    A phylogenetic terrace
    @author Lukasz Reszczynski <lukasz.reszczynski@univie.ac.at>
*/
class TerraceTP {

public:
	/**
	 *  Constructor
	 *  @param tree tree
	 *  @param saln superalignment
	 */
    TerraceTP(PhyloTree &tree, SuperAlignment* saln);

    /**
     *  @return The terrace size
     */
    uint64_t getSize();

    /**
     *  Print trees from the terrace in the newick format.
     *  @param out the output stream
     */
    void printTrees(ostream &out);    

    /**
     *  Print trees from the terrace in the compressed newick format.
     *  @param out the output stream 
     */
    void printTreesCompressed(ostream &out);

	void init();

    ~TerraceTP();

private:
  	terraces::bitmatrix coverage;
  	terraces::name_map names;
  	terraces::index_map indices;
  	terraces::supertree_data supertree;
};


#endif
