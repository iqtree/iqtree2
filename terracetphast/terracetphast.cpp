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
#include "terracetphast.h"

TerraceTP::TerraceTP(PhyloTree &tree, SuperAlignment* saln) :
    coverage(tree.aln->getSeqNames().size(), saln->taxa_index[0].size())
{
    stringstream nwk;
    tree.printTree(nwk, 0);
    
    terraces::index cols{};
    terraces::index rows{};

    terraces::bitmatrix coverage_matrix{rows, cols};

    vector<string> labels = tree.aln->getSeqNames();

    names.resize(labels.size());

    for (int i=0; i<labels.size(); i++) {
        string label = labels[i];

        names[i] = label;
        indices[label] = i;
    }

    auto terraphast_nwk = terraces::parse_nwk(nwk.str(), indices);

    int n_partitions = saln->taxa_index[0].size();

    for (int i=0; i<labels.size(); i++) {
        for (int j=0; j<n_partitions; j++) {
            bool value = (saln->taxa_index[i][j] != -1);
            coverage.set(i, j, value);
        }
    }

    supertree = terraces::create_supertree_data(terraphast_nwk, coverage);

    init();  
}

void TerraceTP::init()
{

}

uint64_t TerraceTP::getSize()
{
    return terraces::count_terrace(supertree);
}

void TerraceTP::printTrees(ostream &out)
{
    terraces::print_terrace(supertree, names, out);
}

void TerraceTP::printTreesCompressed(ostream &out)
{
    terraces::print_terrace_compressed(supertree, names, out);
}

TerraceTP::~TerraceTP()
{
}


