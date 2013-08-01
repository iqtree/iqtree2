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
#include "greedy.h"

/*********************************************
	class Greedy
*********************************************/
/**
	run the algorithm
*/
void Greedy::run(Params &params, vector<PDTaxaSet> &taxa_set)
{
	Node *node1, *node2;
	NodeVector subtree;
	subtree.resize(nodeNum, NULL);

	//if (params.is_rooted) subsize++;

	if (params.min_size < 2) 
		params.min_size = params.sub_size;

	taxa_set.resize(params.sub_size - params.min_size + 1);

	if (initialset.empty()) {
		taxa_set[0].score = root->longestPath2(node1, node2);
		root = node1;
		// initialize the subtree length
		subtree[root->id] = root;
		// update the current PD set
		taxa_set[0].push_back(root);

		list_size = params.sub_size-2;
		updateOnLongestPath(root->highestNei->node, subtree, taxa_set[0]);
	} else if (initialset.size() == 1) {
		root = initialset[0];
		subtree[root->id] = root;
		root->calcHeight();
		// initialize the subtree length
		taxa_set[0].score = root->height;
		taxa_set[0].push_back(root);
		list_size = params.sub_size-2;
		updateOnLongestPath(root->highestNei->node, subtree, taxa_set[0]);
	} else {
		root = initialset[0];
		int included = initialset.size();
		// put all taxa on the initial set into subtree
		for (NodeVector::iterator it = initialset.begin(); it != initialset.end(); it++) {
			if (subtree[(*it)->id]) {
				cout << "Duplicated " << (*it)->name << endl;
				included--;
				continue;
			}
			subtree[(*it)->id] = (*it);
			taxa_set[0].push_back(*it);
		}
		list_size = params.sub_size - included;
		cout << included - rooted << " distinct taxa included, adding " << list_size << " more taxa" << endl;
		// initialize maximal distance set
		root->calcHeight();

		NodeVector nodestack;
		buildOnInitialSet(subtree, nodestack);
		taxa_set[0].score = updateOnInitialSet(subtree);
		//taxa_set[0].insert(taxa_set[0].end(), initialset.begin(), initialset.end());
	}

	// greedy step

	if (list_size < 0) outError("Too small k");

	int ts;
	for (ts = 0; list_size > 0; list_size--)
	{
		if (params.sub_size - list_size >= params.min_size) {
			taxa_set[ts].setSubTree(*this, subtree);
			ts++;
			taxa_set[ts] = taxa_set[ts-1];
		}
		NeighborSet::iterator itneigh = highestNeighbor();
		Neighbor* neigh = *itneigh;
		neighset.erase(itneigh);
		// update the subtree length
		taxa_set[ts].score += neigh->length + neigh->node->height;
		// update the subtree
		updateOnLongestPath(neigh->node, subtree, taxa_set[ts]);
	}

	taxa_set[ts].setSubTree(*this, subtree);

}

/**
	initialize the ordered list based on the initial tree structure
*/
void Greedy::buildOnInitialSet(NodeVector &subtree, NodeVector &nodestack, Node *node, Node *dad) {
	if (!node) node = root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		Node *next = (*it)->node;
		nodestack.push_back(next);

		if (next->isLeaf() && subtree[next->id] != NULL) {
			// the next node is a leaf and is in the initial set
			// put all node on the stack into the subtree
			for (NodeVector::iterator itnode = nodestack.begin(); itnode != nodestack.end(); itnode++) {
				subtree[(*itnode)->id] = (*itnode);
			}
		}
		buildOnInitialSet(subtree, nodestack, next, node);
		nodestack.pop_back();
	}
}

/**
	initialize the ordered list based on the initial subtree structure
	@param subtree vector containing nodes in the subtree
	@return the subtree length
*/
double Greedy::updateOnInitialSet(NodeVector &subtree) {
	int i;
	// scan through interior nodes
	for (i = leafNum; i < nodeNum; i++) 
		if (subtree[i] != NULL) {
			Node *node = subtree[i];
			for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++) 
				if (subtree[(*it)->node->id] == NULL){
					addNeighbor((*it));
				}
		}
	double len = 0.0;
	for (i = 0; i < nodeNum; i++) 
		if (subtree[i] != NULL) {
			Node *node = subtree[i];
			for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++) 
				if (subtree[(*it)->node->id] != NULL){
					len += (*it)->length;
				}
		}
	return len / 2.0;
}

void Greedy::updateOnLongestPath(Node *node, NodeVector &subtree, PDTaxaSet &cur_set)
{

	Node* next;
	Node* current;

	for (current = node; !current->isLeaf(); current = next)
	{
		subtree[current->id] = current;
		next = current->highestNei->node;
		// redirect the highest neighbor of the current
		//for (int i = 0; i < current->neighbors.size(); i++)
			//if (subtree[current->neighbors[i]->node->id] == NULL && current->neighbors[i]->node != next)
		FOR_NEIGHBOR_IT(current, next, it)
			if (subtree[(*it)->node->id] == NULL) {
				addNeighbor((*it));
			}
	}
	subtree[current->id] = current;
	cur_set.push_back(current);
}

NeighborSet::iterator Greedy::highestNeighbor()
{
	return neighset.begin();
}

/**
	add an edge into the NeighborSet
*/
void Greedy::addNeighbor(Neighbor* neigh) {
	if (list_size <= 0)
		return;
	if (neighset.size() < list_size)
		neighset.insert(neigh);
	else {
		NeighborSet::iterator last = neighset.end();
		last--;
		Neighbor* endn = *last;
    	if ((neigh->length + neigh->node->height) > (endn->length + endn->node->height)) {
			neighset.erase(last);
			neighset.insert(neigh);
		}
 	}
}
