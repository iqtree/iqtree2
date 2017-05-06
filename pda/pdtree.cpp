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
#include "ncl/ncl.h"
#include "utils/tools.h"
#include "pdtree.h"
#include "nclextra/myreader.h"

/*********************************************
	class PDTree
*********************************************/
PDTree::PDTree(Params &params)
{
	init(params);
}

void PDTree::init(Params &params) {
	MTree::init(params.user_file, params.is_rooted);
	if (params.is_rooted) {
		params.sub_size++;
		params.min_size++;
	}
	if (params.is_rooted && params.root != NULL) {
		outError(ERR_CONFLICT_ROOT);
	}

	if (params.sub_size > leafNum) {
		ostringstream err;
		err << "Subset size k = " << params.sub_size-params.is_rooted << 
			" is greater than the number of taxa = " << leafNum-params.is_rooted;
		outError(err.str());
	}

	if (params.is_rooted) {
		initialset.push_back(root);
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


/**
	constructor
*/
PDTree::PDTree(PDTree &tree)
{
	init(tree);
}

void PDTree::init(PDTree &tree) {
	MTree::init(tree);
	//subsize = tree.subsize;
	initialset = tree.initialset;
}


void PDTree::buildLeafMapName(LeafMapName &lsn, Node *node, Node* dad) {
	if (!node) node = root;
	if (node->isLeaf()) {
		if (lsn.find(node->name) != lsn.end()) 
			outError(ERR_DUPLICATED_TAXA);
		lsn[node->name] = node;
	}
	//for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
		//if ((*it)->node != dad)
	FOR_NEIGHBOR_IT(node, dad, it)
		buildLeafMapName(lsn, (*it)->node, node);
}

/*
Node *PDTree::findNode(char *name, Node *node, Node *dad) {
	if (!node) node = root;
	// check the name if a leaf
	if (node->isLeaf()) {
		if (node->name == name)
			return node;
	}
	// recursive search
	//for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++) 
		//if ((*it)->node != dad) {
	FOR_NEIGHBOR_IT(node, dad, it) {
		Node *res = findNode(name, (*it)->node, node);
		if (res != NULL)
			return res;
	}
	return NULL;
}
*/

void PDTree::readRootNode(const char *root_name) {
	string name = root_name;
	Node *node = findNodeName(name);
	if (node == NULL)
		outError(ERR_NO_ROOT, root_name);
	initialset.push_back(node);
}


/**
	read the initial set of taxa to be included into PD-tree
*/
void PDTree::readInitialSet(Params &params) {
	LeafMapName lsn;
	buildLeafMapName(lsn);
	int ntaxa = leafNum - params.is_rooted;
	StrVector tax_name;
	readInitTaxaFile(params, ntaxa, tax_name);
	for (StrVector::iterator it = tax_name.begin(); it != tax_name.end(); it++) {
		LeafMapName::iterator nameit = lsn.find((*it));
		if (nameit == lsn.end()) {
			Node *node = findNodeName(*it);
			if (!node)
				cout << "Find no taxon with name " << *it << endl;
			else {
				Node *taxon;
				int distance = findNearestTaxon(taxon, node);
				cout << "Replace internal node " << node->name << " by taxon " 
					 << taxon->name << " (" << distance << " branches away)" << endl;
				initialset.push_back(taxon);
			}
		} else
		initialset.push_back((*nameit).second);
	}
	cout << initialset.size() - rooted << " initial taxa" << endl;
}


void PDTree::readParams(Params &params) {
	int ntaxa = leafNum - params.is_rooted;

	// read parameters from file
	double scale;
	StrVector tax_name;
	DoubleVector ori_weight, tax_weight;
	readWeightFile(params, ntaxa, scale, tax_name, ori_weight);

	// now convert the weights
	LeafMapName lsn;
	buildLeafMapName(lsn);
	tax_weight.resize((unsigned long) ntaxa, 0);
	for (int i = 0; i < tax_name.size(); i++) {
		LeafMapName::iterator nameit = lsn.find(tax_name[i]);
		if (nameit == lsn.end())
			outError(ERR_NO_TAXON, tax_name[i]);
		tax_weight[(*nameit).second->id] = ori_weight[i];
	}

	if (params.scaling_factor >= 0) {
		if (params.scaling_factor > 1) outError("Scaling factor must be between 0 and 1");
		cout << "Rescaling branch lengths with " << params.scaling_factor << 
			" and taxa weights with " << 1 - params.scaling_factor << endl;
		scale = params.scaling_factor;
		for (DoubleVector::iterator it = tax_weight.begin(); it != tax_weight.end(); it++)
			(*it) *= (1 - scale);
	}

	// incoporate them into the tree
	incoporateParams(scale, tax_weight);
}

void PDTree::incoporateParams(double &scale, DoubleVector &tax_weight, Node* node, Node* dad) {
	if (!node) node = root;
	FOR_NEIGHBOR_DECLARE(node, NULL, it) {
		double newlen;
		newlen = (*it)->length * scale;
		if (node->isLeaf())
			newlen += tax_weight[node->id];
		else if ((*it)->node->isLeaf())
			newlen += tax_weight[(*it)->node->id];
		(*it)->length = newlen;
	}
	FOR_NEIGHBOR(node, dad, it)
		incoporateParams(scale, tax_weight, (*it)->node, node);
			
}

void PDTree::computePD(Params &params, vector<PDTaxaSet> &taxa_set, PDRelatedMeasures &pd_more) {
	LeafMapName lsn;
	buildLeafMapName(lsn);

	MSetsBlock *sets;
	TaxaSetNameVector *allsets;
	sets = new MSetsBlock();

 	cout << "Reading taxa sets in file " << params.pdtaxa_file << "..." << endl;

	bool nexus_formated = (detectInputFile(params.pdtaxa_file) == IN_NEXUS);
	if (nexus_formated) {
		MyReader nexus(params.pdtaxa_file);
		nexus.Add(sets);
		MyToken token(nexus.inf);
		nexus.Execute(token);
	} else {
		readTaxaSets(params.pdtaxa_file, sets);
	}

	allsets = sets->getSets();

	//sets->Report(cout);

	taxa_set.resize((unsigned long) sets->getNSets());

	vector<PDTaxaSet>::iterator it_ts;
	TaxaSetNameVector::iterator i;

	for (i = allsets->begin(), it_ts = taxa_set.begin(); i != allsets->end(); i++, it_ts++) {
		set<string> taxa_name;
		for (NodeVector::iterator it = initialset.begin(); it != initialset.end(); it++)
			taxa_name.insert((*it)->name);
		for (vector<string>::iterator it2 = (*i)->taxlist.begin(); it2 != (*i)->taxlist.end(); it2++) {
			LeafMapName::iterator nameit = lsn.find(*it2);
			if (nameit == lsn.end())
				outError(ERR_NO_TAXON, *it2);
			taxa_name.insert(*it2);
		}

		Split id_set;
		makeTaxaSet(taxa_name, *it_ts);
		(*it_ts).makeIDSet(leafNum, id_set);
		if (params.exclusive_pd) {
			calcExclusivePD(id_set);
			pd_more.exclusivePD.push_back(id_set.getWeight());
		}
		calcPD(id_set);
		(*it_ts).score = id_set.getWeight();
		(*it_ts).name = (*i)->name;
		pd_more.PDScore.push_back(id_set.getWeight());
		pd_more.setName.push_back((*i)->name);
	}

	delete sets;
}



void PDTree::makeTaxaSet(set<string> &taxa_name, PDTaxaSet &taxa_set, Node *node, Node *dad) {
	if (!node) node = root;
	if (node->isLeaf() && taxa_name.find(node->name) != taxa_name.end()) {
		taxa_set.push_back(node);
	}
	FOR_NEIGHBOR_IT(node, dad, it) {
		makeTaxaSet(taxa_name, taxa_set, (*it)->node, node);
	}
}

bool PDTree::calcPD(Split &id_set, double cur_len, Node *node, Node *dad) {
	if (!node) { 
		node = root; 
		id_set.weight = 0.0;
		if (!rooted && !id_set.containTaxon(node->id)) {
			int id = id_set.firstTaxon();
			if (id < 0) return false;
			node = findNodeID(id);
		}
	}

	bool resval = false;

	if (node->isLeaf() && id_set.containTaxon(node->id)) {
		id_set.weight += cur_len;
		resval = true;
	}
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (calcPD(id_set, cur_len + (*it)->length, (*it)->node, node)) {
			cur_len = 0.0;
			resval = true;
		}
	}
	return resval;
}

void PDTree::calcExclusivePD(Split &id_set) {
	id_set.invert();
	calcPD(id_set);
	id_set.invert();
	id_set.weight = treeLength() - id_set.weight;
}


void PDTree::calcPDEndemism(vector<PDTaxaSet> &area_set, DoubleVector &pd_endem) {
	vector<Split> id_sets;
	vector<Split>::iterator it_s;
	vector<PDTaxaSet>::iterator it_a;

	// convert taxa set to id set
	id_sets.resize(area_set.size());
	for (it_a = area_set.begin(), it_s = id_sets.begin(); it_a != area_set.end(); it_a++, it_s++) 
		(*it_a).makeIDSet(leafNum, *it_s);

	// make union of all id_sets
	Split id_union(leafNum);
	for (it_s = id_sets.begin(); it_s != id_sets.end(); it_s++) 
		id_union += *it_s;
	
	// calculate PD of union 
	calcPD(id_union);

	// now calculate PD endemism
	pd_endem.clear();
	for (it_s = id_sets.begin(); it_s != id_sets.end(); it_s++) {
		// make union of all other set
		Split id_other(leafNum);
		for (vector<Split>::iterator it_s2 = id_sets.begin(); it_s2 != id_sets.end(); it_s2++)
			if (it_s2 != it_s) id_other += *it_s2;
		// calculate PD of all other sets
		calcPD(id_other);

		// calc PD endemism
		pd_endem.push_back(id_union.weight - id_other.weight);
	}
}


void PDTree::calcPDComplementarity(vector<PDTaxaSet> &area_set, char *area_names, DoubleVector &pd_comp) {

	set<string> given_areas;

	parseAreaName(area_names, given_areas);

/*
	for (set<string>::iterator it = given_areas.begin(); it != given_areas.end(); it++)
		cout << (*it) << "!";
	cout << endl;
*/
	vector<Split> id_sets;
	vector<Split>::iterator it_s;
	vector<PDTaxaSet>::iterator it_a;

	Split given_id(leafNum);

	// convert taxa set to id set
	id_sets.resize(area_set.size());
	for (it_a = area_set.begin(), it_s = id_sets.begin(); it_a != area_set.end(); it_a++, it_s++) {
		(*it_a).makeIDSet(leafNum, *it_s);
		if (given_areas.find((*it_a).name) != given_areas.end())
			given_id += *it_s;
	}
	
	if (given_id.countTaxa() == 0)
		outError("Complementary area name(s) not correct");
	calcPD(given_id);

	

	// now calculate PD complementarity
	pd_comp.clear();
	for (it_s = id_sets.begin(); it_s != id_sets.end(); it_s++) {
		// make union the two sets
		Split id_both(*it_s);
		id_both += given_id;
		// calculate PD of both sets
		calcPD(id_both);
		// calc PD complementarity
		pd_comp.push_back(id_both.weight - given_id.weight);
	}

}
int PDTree::findNearestTaxon(Node* &taxon, Node *node, Node *dad) {
	if (node->isLeaf()) {
		taxon = node;
		return 0;
	}
	int distance = 10000000;
	taxon = NULL;
	FOR_NEIGHBOR_IT(node, dad, it) {
		Node *mytaxon;
		int mydistance = findNearestTaxon(mytaxon, (*it)->node, node);
		if (mydistance < distance) {
			distance = mydistance;
			taxon = mytaxon;
		}
	}
	return distance+1;
}
