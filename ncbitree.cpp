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
#include "ncbitree.h"

NCBITree::NCBITree()
        : MTree()
{
}


NCBITree::~NCBITree()
{
}

void NCBITree::readNCBINames(const char* infile, const char *name_type) {
    ifstream in;
    cout << "Reading NCBI names file " << infile << endl;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        readNCBINames(in, name_type);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void NCBITree::readNCBINames(ifstream &in, const char *name_type) {
    ASSERT(!nodes.empty());
    char ch;
    int node_id;
    string node_name, unique_name;

    in_line = in_column = 0;

    while (!in.eof()) {
        node_id = 0;
        if (!(in >> node_id)) break;
        in_line++;
        if (node_id <= 0) throw "Wrong node ID";
        if (node_id > nodes.size()) throw "Too large node ID";
        if (nodes[node_id]) {
            in >> ch;
            if (ch != '|') throw "No | between node ID and name";
            in.get(ch);
            getline(in, node_name, '\t');
            if (node_name == "") throw "Empty node name";
            in >> ch;
            if (ch != '|') throw "No | between name and unique name";
            in.get(ch);
            getline(in, unique_name, '\t');
            if (unique_name != "") node_name = unique_name;

            for (string::iterator i = node_name.begin(); i != node_name.end(); i++) {
                if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.') {
                    (*i) = '_';
                }
            }
            nodes[node_id]->name = node_name;
        }
        // get the rest of the line
        string str;
        getline(in, str);
    }

}

Node *NCBITree::readNCBITree(const char *infile, int root_id, const char* taxon_level, const char *ignore_level) {
    ifstream in;
    cout << "Reading NCBI nodes file " << infile << endl;
    Node *parent = NULL;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        parent = readNCBITree(in, root_id, taxon_level, ignore_level);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }

    return parent;
}

Node* NCBITree::readNCBITree(istream &in, int root_id, const char* taxon_level, const char *ignore_level) {
    //IntVector parents_id;
    nodes.resize(MAX_TAXONOMY_ID, NULL);
    node_levels.resize(MAX_TAXONOMY_ID);
    string node_level;
    int node_id, parent_id, max_node_id = 0, num_nodes = 0;
    char ch;
    in_line = in_column = 0;

    while (!in.eof()) {
        node_id = parent_id = 0;
        if (!(in >> node_id)) break;
        in_line++;
        num_nodes ++;
        if (node_id <= 0) throw "Wrong node ID";
        if (node_id >= nodes.size()) throw "Too large node ID";
        in >> ch;
        if (ch != '|') throw "No | between node ID and parent ID";
        in >> parent_id;
        if (parent_id <= 0) throw "Wrong parent ID";
        if (parent_id >= nodes.size()) throw "Too large parent ID";
        in >> ch;
        if (ch != '|') throw "No | between parent ID and node rank";
        in.get(ch);
        getline(in,node_level,'\t');

        string str;
        getline(in, str);
        if (node_id > max_node_id) max_node_id = node_id;
        if (nodes[node_id]) throw "Duplicated node ID";
        nodes[node_id] = newNode(node_id, node_id);
        nodes[node_id]->height = parent_id; // use height temporarily for parent_id
        node_levels[node_id] = node_level;
    }

    nodes.resize(max_node_id+1);
    node_levels.resize(max_node_id+1);
    int ignored = 0;

    for (node_id = 0; node_id <= max_node_id; node_id++)
        if (nodes[node_id]) {
            parent_id = nodes[node_id]->height;
            if (!nodes[parent_id]) throw "Parent ID not found";
            if (parent_id == node_id) {
                cout << "Ignore " << node_id << " | " << parent_id << endl;
                continue;
            }
            double len = 1.0;
            if (ignore_level && node_levels[node_id] == ignore_level) {
                len = 0.0;
                ignored++;
            }
            nodes[node_id]->addNeighbor(nodes[parent_id], len);
            nodes[parent_id]->addNeighbor(nodes[node_id], len);
        }

    if (ignore_level)
        cout << ignored << " branches are set to zero because of " << ignore_level << endl;

    rooted = true;
    if (!nodes[root_id]) throw "Root node not available";
    root = nodes[root_id];

    if (taxon_level) {
        int pruned = pruneTaxa(node_levels, taxon_level, root, nodes[root->height]);
        cout << pruned << " nodes below " << taxon_level << " are pruned" << endl;
    }

//	int pruned = pruneBridgeNodes(root, nodes[root->height]);
//	cout << pruned << " nodes of degree 2 are pruned" << endl;

    leafNum = nodeNum = branchNum = 0;
    countNodeNum(root, nodes[root->height]);

    /*	for (node_id = 0; node_id <= max_node_id; node_id++)
    	if (nodes[node_id] && nodes[node_id]->isLeaf()) {
    		Node *taxon = nodes[node_id];
    		taxon->id = leafNum;
    		leafNum++;
    	}
    	initializeTree();*/


    cout << num_nodes << " NCBI nodes, " << nodeNum << " tree nodes, " << leafNum << " leaves, " << branchNum << " branches" << endl;
    return nodes[nodes[root_id]->height];
}

int NCBITree::pruneTaxa(StrVector &node_levels, const char* taxon_level, Node *node, Node *dad) {
    int num_nodes = 0;
    //if (node_levels[node->id].find(taxon_level) != string::npos) {
    if (node_levels[node->id] == taxon_level) {
        // prune subtree below node
        Neighbor *node_nei = node->findNeighbor(dad);
        FOR_NEIGHBOR_IT(node, dad, it) {
            num_nodes += freeNode((*it)->node, node);
            delete (*it);
        }
        node->neighbors.resize(1);
        node->neighbors[0] = node_nei;
        return num_nodes;
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    num_nodes += pruneTaxa(node_levels, taxon_level, (*it)->node, node);
    return num_nodes;
}


void NCBITree::countNodeNum(Node *node, Node *dad) {
    nodeNum++;
    if (node->isLeaf()) leafNum++;
    FOR_NEIGHBOR_IT(node, dad, it) {
        branchNum++;
        countNodeNum((*it)->node, node);
    }
}

int NCBITree::pruneBridgeNodes(Node *node, Node *dad) {
    int num_nodes = 0;
    FOR_NEIGHBOR_IT(node, dad, it)
    num_nodes += pruneBridgeNodes((*it)->node, node);
    if (node->neighbors.size() == 2) {
        Node *child;
        if (node->neighbors[0]->node == dad)
            child = node->neighbors[1]->node;
        else
            child = node->neighbors[0]->node;
        double len = node->neighbors[0]->length + node->neighbors[1]->length;
        dad->updateNeighbor(node, child, len);
        child->updateNeighbor(node, dad, len);
        nodes[node->id] = NULL;
        delete node;
        num_nodes++;
    }
    return num_nodes;
}

int NCBITree::freeNode(Node *node, Node *dad)
{
    if (!node) node = root;
    NeighborVec::reverse_iterator it;
    int num_nodes = 1;
    for (it = node->neighbors.rbegin(); it != node->neighbors.rend(); it++)
        if ((*it)->node != dad) {
            num_nodes += freeNode((*it)->node, node);
        }
    nodes[node->id] = NULL;
    delete node;
    return num_nodes;
}
