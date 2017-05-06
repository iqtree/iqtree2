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
#include "pdtreeset.h"

PDTreeSet::PDTreeSet()
 : MTreeSet()
{
}


PDTreeSet::PDTreeSet(Params &params) {
	init(params);
}

void PDTreeSet::init(Params &params) {
	MTreeSet::init(params.user_file, params.is_rooted, params.tree_burnin, params.tree_max_count);

	if (isRootedTrees()) {
		params.sub_size++;
		params.min_size++;
	}
	if (isRootedTrees() && params.root != NULL) {
		outError(ERR_CONFLICT_ROOT);
	}

	if (params.sub_size > getNTaxa()) {
		ostringstream err;
		err << "Subset size k = " << params.sub_size - params.is_rooted <<
			" is greater than the number of taxa = " << getNTaxa() - params.is_rooted;
		outError(err.str());
	}

	if (isRootedTrees()) {
		char *rname = (char*)ROOT_NAME;
		readRootNode(rname);
	}
	// read the parameter file
	if (params.param_file != NULL) {
		readParams(params);
	}
	// identify the root
	if (params.root != NULL) 
		readRootNode(params.root);

	// read the initial set of taxa, incoporate info into the split system
	if (params.initial_file != NULL) {
		readInitialSet(params);
	}
}

bool PDTreeSet::isRootedTrees() {
	ASSERT(size() > 0);
	return front()->rooted;
}

int PDTreeSet::getNTaxa() {
	ASSERT(size() > 0);
	return front()->leafNum;
}

void PDTreeSet::readRootNode(const char *root_name) {
	string name = root_name;
	init_taxa.push_back(name);
	for (iterator it = begin(); it != end(); it++)
		((PDTree*)(*it))->readRootNode(root_name);
}

void PDTreeSet::readParams(Params &params) {

	int ntaxa = getNTaxa() - params.is_rooted;

	// read parameters from file
	double scale;
	StrVector tax_name;
	DoubleVector ori_weight;
	readWeightFile(params, ntaxa, scale, tax_name, ori_weight);

	for (iterator it = begin(); it != end(); it++) {
		// now convert the weights
		PDTree *mytree = (PDTree*)(*it);
		LeafMapName lsn;
		mytree->buildLeafMapName(lsn);
		DoubleVector tax_weight;
		tax_weight.resize(ntaxa, 0);
		for (int i = 0; i < tax_name.size(); i++) {
			LeafMapName::iterator nameit = lsn.find(tax_name[i]);
			if (nameit == lsn.end())
				outError(ERR_NO_TAXON, tax_name[i]);
			tax_weight[(*nameit).second->id] = ori_weight[i];
		}
	
		// incoporate them into the tree
		mytree->incoporateParams(scale, tax_weight);
	}
}

/**
	read the initial set of taxa to be included into PD-tree
*/
void PDTreeSet::readInitialSet(Params &params) {
	int ntaxa = getNTaxa() - params.is_rooted;
	StrVector tax_name;
	readInitTaxaFile(params, ntaxa, tax_name);
	init_taxa.insert(init_taxa.end(), tax_name.begin(), tax_name.end());

	for (iterator itree = begin(); itree != end(); itree++) {
		PDTree *mytree = (PDTree*)(*itree);
		LeafMapName lsn;
		mytree->buildLeafMapName(lsn);
		for (StrVector::iterator it = tax_name.begin(); it != tax_name.end(); it++) {
			LeafMapName::iterator nameit = lsn.find((*it));
			if (nameit == lsn.end()) {
				outError(ERR_NO_TAXON, *it);
			}
			mytree->initialset.push_back((*nameit).second);
		}
	}
}
