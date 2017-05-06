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
#include "pruning.h"
#include <algorithm>

/*********************************************
	class Pruning
*********************************************/
/**
	Steffen Klaere's pruning algorithm
*/
void Pruning::run(Params &params, vector<PDTaxaSet> &taxa_set)
{
	//if (params.min_size < 2) 
		params.min_size = params.sub_size;

	list_size = 2*(leafNum-params.sub_size)-1;
	if (!initialset.empty()) {
		doInitialSet();
	}
	buildLeaves();
	for (int step = leafNum; step > params.sub_size; step--)
	{
		deleteExNode(nearestLeaf());
		list_size -= 2;
	}
	taxa_set.resize(1);
	taxa_set[0].setTree(*this);

}

void Pruning::doInitialSet() {
	for (NodeVector::iterator it = initialset.begin(); it != initialset.end(); it++) {
		(*it)->height = 1;
	}
}

void Pruning::printLeaves()
{
	// print info
	for (LeafSet::iterator it = leaves.begin(); it != leaves.end(); it++)
	{
		Node *node = *it;
		cout << node->id << " " << node->neighbors[0]->length << endl;
	}
}


void Pruning::buildLeaves(Node *node, Node *dad)
{
	if (!node) node = root;
	if (node->isLeaf())
		addLeaf(node);
	FOR_NEIGHBOR_IT(node, dad, it)
		buildLeaves((*it)->node, node);
}

LeafSet::iterator Pruning::findNode(Node *node)
{
	pair<LeafSet::iterator, LeafSet::iterator>
	range = leaves.equal_range(node);
	for (LeafSet::iterator it = range.first; it != range.second; it++)
		if (*it == node)
			return it;
	return leaves.end();
}

void Pruning::deleteExNode(LeafSet::iterator pos)
{
	// delete from the tree
	Node *node = *pos;
	Node *innode = node->neighbors[0]->node;
	Node *othernodes[2] = { NULL, NULL };
	int i;
	NeighborVec::iterator it;
	double length = 0;

	bool should_merge = true;
	//for (it = innode->neighbors.begin(); it != innode->neighbors.end(); it++)
		//if ((*it)->node != node)
	FOR_NEIGHBOR(innode, node, it)	{
		length += (*it)->length;
		if (othernodes[0] == NULL)
			othernodes[0] = (*it)->node;
		else if (othernodes[1] == NULL)
			othernodes[1] = (*it)->node;
		else
			should_merge = false;
	}

	if (should_merge)
	{	
		// merge two branches
		for (i = 0; i < 2; i++)
			if (othernodes[i]->isLeaf()) {
				LeafSet::iterator temp = findNode(othernodes[i]);
				if (temp != leaves.end())
					leaves.erase(temp);
			}

		for (i = 0; i < 2; i++)
			for (it = othernodes[i]->neighbors.begin(); it != othernodes[i]->neighbors.end(); it++)
				if ((*it)->node == innode)
				{
					(*it)->node = othernodes[1-i];
					(*it)->length = length;
				}
	} else {
		// simple delete the neighbor of innode
		for (it = innode->neighbors.begin(); it != innode->neighbors.end(); it++)
			if ((*it)->node == node) {
				innode->neighbors.erase(it);
				break;
			}
	}
	// delete the first element
	leaves.erase(pos);
	if (should_merge)
	for (i = 0; i < 2; i++)
		if (othernodes[i]->isLeaf())
			addLeaf(othernodes[i]);

	// also delete the last element if necessary
	if (leaves.size() > list_size && leaves.size() > 1) {
		LeafSet::iterator last = leaves.end();
		last--;
		leaves.erase(last);
	}	

	if (node == root)
	{	// find another root
		root = *(leaves.begin());
	}
}

LeafSet::iterator Pruning::nearestLeaf()
{
	return leaves.begin();
}


/**
	insert a leaf into the LeafSet
*/
void Pruning::addLeaf(Node* leaf) {

	// if leaf in the initial set
	if (leaf->height == 1) {
		return;
	}

	if (list_size <= 0)
		return;
	if (leaves.size() < list_size) 
		leaves.insert(leaf);
	else {
		LeafSet::iterator last = leaves.end();
		last--;
		Node* endp = *last;
		if (leaf->neighbors[0]->length < endp->neighbors[0]->length) {
			leaves.erase(last);
			leaves.insert(leaf);
		}
	}
}

