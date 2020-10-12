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
#include "mtree.h"
#include <iostream>
//#include <fstream>
#include <iterator>
//#include <mtree.h>
#include "pda/splitgraph.h"
#include "utils/tools.h"
#include "mtreeset.h"
using namespace std;

/*********************************************
	class MTree
*********************************************/

MTree::MTree() {
    root = NULL;
    leafNum = 0;
    nodeNum = 0;
    rooted = false;
    if (Params::getInstance().min_branch_length <= 0)
        num_precision = 6;
    else
        num_precision = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
    len_scale = 1.0;
	fig_char = "|-+++";
}

MTree::MTree(const char *userTreeFile, bool &is_rooted)
{
    init(userTreeFile, is_rooted);
}

void MTree::init(const char *userTreeFile, bool &is_rooted) {
    if (Params::getInstance().min_branch_length <= 0)
        num_precision = 6;
    else
        num_precision = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
    len_scale = 1.0;
    readTree(userTreeFile, is_rooted);
    //printInfo();
	fig_char = "|-+++";
}


/**
	constructor, get from another tree
*/
MTree::MTree(MTree &tree) {
    init(tree);
}

MTree::MTree(string& treeString, vector<string>& taxaNames, bool isRooted) {
    stringstream str;
    str << treeString;
    str.seekg(0, ios::beg);
    readTree(str, isRooted);
    assignIDs(taxaNames);
    assignLeafID();
}

MTree::MTree(string& treeString, bool isRooted) {
    stringstream str;
    str << treeString;
    str.seekg(0, ios::beg);
    readTree(str, isRooted);
    assignLeafID();
}

void MTree::init(MTree &tree) {
    root = tree.root;
    leafNum = tree.leafNum;
    nodeNum = tree.nodeNum;
    rooted = tree.rooted;
    //userFile = tree.userFile;
    // have to delete the root when exchange to another object
    tree.root = NULL;
    num_precision = tree.num_precision;
    len_scale = tree.len_scale;
    fig_char = tree.fig_char;
}

void MTree::assignIDs(vector<string>& taxaNames) {
    bool err = false;
    int nseq = taxaNames.size();
    for (int seq = 0; seq < nseq; seq++) {
        string seq_name = taxaNames[seq];
        Node *node = findLeafName(seq_name);
        if (!node) {
            string str = "Sequence ";
            str += seq_name;
            str += " does not appear in the tree";
            err = true;
            outError(str, false);
        } else {
            ASSERT(node->isLeaf());
            node->id = seq;
        }
    }
    StrVector taxname;
    getTaxaName(taxname);
    for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
        bool foundTaxa = false;
        for (vector<string>::iterator it2 = taxaNames.begin(); it2 != taxaNames.end(); it2++) {
            if ( *it == *it2 ) {
                foundTaxa = true;
                break;
            }
        }
        if (!foundTaxa) {
            outError((string) "Tree taxon " + (*it) + " does not appear in the input taxa names", false);
            err = true;
        }
    }
    if (err) outError("Tree taxa and input taxa names do not match (see above)");
}

void MTree::copyTree(MTree *tree) {
    if (root) freeNode();
    stringstream ss;
    tree->printTree(ss);
    readTree(ss, tree->rooted);
    rooted = tree->rooted;
}

void MTree::copyTree(MTree *tree, string &taxa_set) {
    rooted = tree->rooted;
    if (rooted) {
        ASSERT(tree->rooted);
        taxa_set.push_back(1);
    }
//    if (tree->leafNum != taxa_set.length())
//        outError("#leaves and taxa_set do not match!");
    leafNum = nodeNum = branchNum = 0;
    for (string::iterator it = taxa_set.begin(); it != taxa_set.end(); it++)
        nodeNum += (*it);
    double new_len;
    if (root) freeNode();
    root = NULL;
    root = copyTree(tree, taxa_set, new_len);
    if (rooted)
        ASSERT(root->name == ROOT_NAME);
}

Node* MTree::copyTree(MTree *tree, string &taxa_set, double &len, Node *node, Node *dad) {
    if (!node) {
        if (taxa_set[tree->root->id]) {
            node = tree->root;
        } else {
            for (int i = 0; i < tree->leafNum; i++)
                if (taxa_set[i]) {
                    node = tree->findNodeID(i);
                    break;
                }
        }
    }
    Node *new_node = NULL;
    NodeVector new_nodes;
    DoubleVector new_lens;
    if (node->isLeaf()) {
        len = 0.0;
        if (taxa_set[node->id]) {
            new_node = newNode(leafNum++, node->name.c_str());
        }
        if (dad) return new_node;
    }
    if (new_node) {
        new_nodes.push_back(new_node);
        new_lens.push_back(len);
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        double new_len;
        new_node = copyTree(tree, taxa_set, new_len, (*it)->node, node);
        if (new_node) {
            new_nodes.push_back(new_node);
            new_lens.push_back((*it)->length + new_len);
        }
    }
    if (new_nodes.empty()) return NULL;
    if (new_nodes.size() == 1) {
        len = new_lens[0];
        return new_nodes[0];
    }
    if (!dad && new_nodes.size() == 2) {
        double sum_len = new_lens[0] + new_lens[1];
        new_nodes[0]->addNeighbor(new_nodes[1], sum_len, branchNum);
        new_nodes[1]->addNeighbor(new_nodes[0], sum_len, branchNum);
        branchNum++;
        return new_nodes[0];
    }
    Node* int_node = newNode(nodeNum++, node->name.c_str());
    len = 0.0;
    for (int i = 0; i < new_nodes.size(); i++) {
        int_node->addNeighbor(new_nodes[i], new_lens[i], branchNum);
        new_nodes[i]->addNeighbor(int_node, new_lens[i], branchNum);
        branchNum++;
    }
    return int_node;
}

void MTree::extractBifurcatingSubTree(Node *node, Node *dad) {
    if (!node) node = root;
    if (node->degree() > 3) {
        int id1, id2, id3;
        id1 = node->findNeighborIt(dad) - node->neighbors.begin();
        do {
            id2 = random_int(node->degree());
        } while (id2 == id1);
        
        // make sure that id1 < id2
        if (id1 > id2) {
            int tmp = id1;
            id1 = id2;
            id2 = tmp;
        }
        do {
            id3 = random_int(node->degree());
        } while (id3 == id1 || id3 == id2);
        //make sure that id1 < id2 < id3
        if (id3 < id2) {
            if (id3 < id1) {
                // id3 < id1 < id2
                int tmp = id1;
                id1 = id3;
                id3 = id2;
                id2 = tmp;
            } else {
                // id1 < id3 < id2
                int tmp = id2;
                id2 = id3;
                id3 = tmp;
            }
        }
        // remove all neighbors except id1, id2, id3
        for (int i = 0; i != node->neighbors.size(); i++)
            if (i != id1 && i != id2 && i != id3) {
                freeNode(node->neighbors[i]->node, node);
                delete node->neighbors[i];
            }
        node->neighbors[0] = node->neighbors[id1];
        node->neighbors[1] = node->neighbors[id2];
        node->neighbors[2] = node->neighbors[id3];
        node->neighbors.erase(node->neighbors.begin()+3, node->neighbors.end());
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->node->isLeaf())
            extractBifurcatingSubTree((*it)->node, node);
    }
}

void MTree::resolveMultifurcation() {
    // randomly resolve multifurcating node

    NodeVector nodes;
    getInternalNodes(nodes);
    for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++)
        while ((*it)->degree() > 3) {
            Node *new_node = newNode();
            int id1 = random_int((*it)->degree());
            int id2;
            do {
                id2 = random_int((*it)->degree());
            } while (id2 == id1);
            
            // make sure that id1 < id2
            if (id1 > id2) {
                int tmp = id1;
                id1 = id2;
                id2 = tmp;
            }
            Neighbor *nei1 = (*it)->neighbors[id1];
            Neighbor *nei2 = (*it)->neighbors[id2];
            
            // connect id1 with new_node
            nei1->node->updateNeighbor((*it), new_node);
            new_node->neighbors.push_back(nei1);
            
            // connect id2 with new_node
            nei2->node->updateNeighbor((*it), new_node);
            new_node->neighbors.push_back(nei2);
            
            // connect new_node with old node
            new_node->addNeighbor((*it), -1.0);
            (*it)->neighbors.erase((*it)->neighbors.begin() + id2);
            (*it)->neighbors.erase((*it)->neighbors.begin() + id1);
            (*it)->addNeighbor(new_node, -1.0);
        }
}

Node* MTree::newNode(int node_id, const char* node_name) {
    return new Node(node_id, node_name);
}

Node* MTree::newNode(int node_id, int node_name) {
    return new Node(node_id, node_name);
}

bool MTree::isBifurcating(Node *node, Node *dad) {
	if (!node) node = root;
	if (!node->isLeaf() && node->degree() != 3) return false;
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!(*it)->node->isLeaf() && (*it)->node->degree() != 3) return false;
		if (!isBifurcating((*it)->node, node)) return false;
	}
	return true;
}

void MTree::printBranchLengths(ostream &out, Node *node, Node *dad)
{
    if (node == NULL) {
    	node = root;
    	sortTaxa();
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (node->name != "") out << node->name; else out << node->id;
        out << "\t";
        if ((*it)->node->name != "") out << (*it)->node->name; else out << (*it)->node->id;
        out << "\t" << (*it)->length << endl;
        printBranchLengths(out, (*it)->node, node);
    }
}

int MTree::countZeroBranches(Node *node, Node *dad, double epsilon) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length <= epsilon) count++;
        count += countZeroBranches((*it)->node, node, epsilon);
    }
    return count;
}

int MTree::countZeroInternalBranches(Node *node, Node *dad, double epsilon) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length <= epsilon && !(*it)->node->isLeaf() && !node->isLeaf()) count++;
        count += countZeroInternalBranches((*it)->node, node, epsilon);
    }
    return count;

}

int MTree::countLongBranches(Node *node, Node *dad, double upper_limit) {
    int count = 0;
    if (node == NULL) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length >= upper_limit) count++;
        count += countLongBranches((*it)->node, node, upper_limit);
    }
    return count;
}


void MTree::printTree(const char *ofile, int brtype)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (brtype & WT_APPEND)
            out.open(ofile, ios_base::out | ios_base::app);
        else
            out.open(ofile);
        printTree(out, brtype);
        out.close();
        if (verbose_mode >= VB_DEBUG)
            cout << "Tree was printed to " << ofile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, ofile);
    }
}

void MTree::printNexus(string ofile, int brtype, string nexus_comment)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (brtype & WT_APPEND)
            out.open(ofile, ios_base::out | ios_base::app);
        else
            out.open(ofile);
        out << "#NEXUS" << endl;
        if (!nexus_comment.empty())
            out << "[ " << nexus_comment << " ]" << endl;
        out << "begin trees;" << endl;
        out << "  tree tree_1 = ";
        printTree(out, brtype | WT_BR_ATTR);
        out << endl;
        out << "end;" << endl;
        out.close();
        if (verbose_mode >= VB_DEBUG)
            cout << "Tree was printed to " << ofile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, ofile);
    }
}

void MTree::printTree(ostream &out, int brtype) {
    if (root->isLeaf()) {
        if (root->neighbors[0]->node->isLeaf()) {
            // tree has only 2 taxa!
            out << "(";
            printTree(out, brtype, root);
            out << ",";
            if (brtype & WT_TAXON_ID)
                out << root->neighbors[0]->node->id;
            else
                out << root->neighbors[0]->node->name;

            if (brtype & WT_BR_LEN)
                out << ":0";
            out << ")";
        } else
            // tree has more than 2 taxa
            printTree(out, brtype, root->neighbors[0]->node);
    } else
        printTree(out, brtype, root);

    out << ";";
    if (brtype & WT_NEWLINE) out << endl;
}

struct IntString {
    int id;
    string str;
};

/**
	nodecmp, for pruning algorithm
*/
struct IntStringCmp
{
    /**
    	nodecmp, for pruning algorithm
    */
    bool operator()(const IntString* s1, const IntString* s2) const
    {
        return (s1->id) < (s2->id);
    }
};

typedef set<IntString*, IntStringCmp> IntStringSet;

void MTree::printBranchLength(ostream &out, int brtype, bool print_slash, Neighbor *length_nei) {
    if (length_nei->length == -1.0)
        return; // NA branch length
    int prec = 10;
	double length = length_nei->length;
    if (brtype & WT_BR_SCALE) length *= len_scale;
    if (brtype & WT_BR_LEN_SHORT) prec = 6;
    if (brtype & WT_BR_LEN_ROUNDING) length = round(length);
    out.precision(prec);
    if ((brtype & WT_BR_ATTR) && !length_nei->attributes.empty()) {
        // print branch attributes
        out << "[&";
        bool first = true;
        for (auto attr : length_nei->attributes) {
            if (!first)
                out << ",";
            out << attr.first << "=\"" << attr.second << '"';
            first = false;
        }
        out << "]";
    }
    
    if (brtype & WT_BR_LEN) {
        if (brtype & WT_BR_LEN_FIXED_WIDTH)
            out << ":" << fixed << length;
        else
            out << ":" << length;
    } else if (brtype & WT_BR_CLADE && length_nei->node->name != ROOT_NAME) {
    	if (print_slash)
    		out << "/";
        out << length;
    }
}

int MTree::printTree(ostream &out, int brtype, Node *node, Node *dad)
{
    int smallest_taxid = leafNum;
    out.precision(num_precision);
    if (!node) node = root;
    if (node->isLeaf()) {
        smallest_taxid = node->id;
        if (brtype & WT_TAXON_ID)
            out << node->id;
        else
            out << node->name;

        if (brtype & WT_BR_LEN) {
        	out.setf( std::ios::fixed, std:: ios::floatfield ); // some sofware does handle number format like '1.234e-6'
//            out.precision(10); // increase precision to avoid zero branch (like in RAxML)
            printBranchLength(out, brtype, false, node->neighbors[0]);
//        	double len = node->neighbors[0]->length;
//            if (brtype & WT_BR_SCALE) len *= len_scale;
//            if (brtype & WT_BR_LEN_ROUNDING) len = round(len);
//            if (brtype & WT_BR_LEN_FIXED_WIDTH)
//                out << ":" << fixed << len;
//            else
//                out << ":" << len;
        }
    } else {
        // internal node
        out << "(";
        bool first = true;
        Neighbor *length_nei = NULL;
        //for (int i = 0; i < node->neighbors.size(); i++)
        //if (node->neighbors[i]->node != dad)
        if (! (brtype & WT_SORT_TAXA)) {
            FOR_NEIGHBOR_IT(node, dad, it) {
                if ((*it)->node->name != ROOT_NAME) {
                    if (!first)
                        out << ",";
                    int taxid = printTree(out, brtype, (*it)->node, node);
                    if (taxid < smallest_taxid) smallest_taxid = taxid;
                    first = false;
                } else
                    length_nei = (*it);
            } else {
                length_nei = (*it);
            }
        } else {
            IntStringSet strout;
            FOR_NEIGHBOR_IT(node, dad, it) {
                if ((*it)->node->name != ROOT_NAME) {
                    ostringstream ss;
                    IntString *str = new IntString;
                    str->id = printTree(ss, brtype, (*it)->node, node);
                    //ss.flush();
                    str->str = ss.str();
                    strout.insert(str);
                } else
                	length_nei = (*it);
            } else {
            	length_nei = (*it);
            }
            smallest_taxid = (*strout.begin())->id;
            IntStringSet::iterator iss;
            for (iss = strout.begin(); iss != strout.end(); iss++) {
                if (!first) out << ",";
                out << (*iss)->str;
                first = false;
            }
            for (iss = strout.begin(); iss != strout.end(); iss++)
                delete (*iss);
        }
        out << ")";
        if (brtype & WT_INT_NODE)
            out << node->id;
        else if (!node->name.empty() && (brtype & WT_BR_ATTR) == 0)
            out << node->name;
        if (dad != NULL || length_nei) {
        	printBranchLength(out, brtype, !node->name.empty(), length_nei);
        }
    }
    return smallest_taxid;
}


void MTree::printSubTree(ostream &out, NodeVector &subtree) {
    if (root->isLeaf())
        printSubTree(out, subtree, root->neighbors[0]->node);
    else
        printSubTree(out, subtree, root);
    out << ";";
}

void MTree::printSubTree(ostream &out, NodeVector &subtree, Node *node, Node *dad) {
    if (!node) node = root;

    NeighborVec::iterator it;
    double length = 0.0, dad_length = 0.0;
    // go down if only 1 child available
    Node *child = NULL;
    int degree;
    do {
        degree = 0;
        FOR_NEIGHBOR(node, dad, it) {
            if (subtree[(*it)->node->id] != NULL) {
                degree++;
                child = (*it)->node;
            }
        } else dad_length = (*it)->length;

        if (degree == 1) {
            dad = node;
            node = child;
            length += dad_length;
        }
    } while (degree == 1 && !node->isLeaf());

    if (node->isLeaf())
        out << node->name << ":" << node->neighbors[0]->length + length;
    else
    {
        // internal node
        out << "(";
        bool first = true;

        FOR_NEIGHBOR(node, dad, it)	{
            if (subtree[(*it)->node->id] != NULL) {
                if ((*it)->node->name != ROOT_NAME) {
                    if (!first)
                        out << ",";
                    printSubTree(out, subtree, (*it)->node, node);
                    first = false;
                } else
                    length += (*it)->length;
            }
        } else {
            length += (*it)->length;
        }
        out << ")";
        if (!node->name.empty())
            out << node->name;
        if (dad != NULL || length > 1e-20)
            out << ":" << length;
    }
}


void MTree::printTaxa(const char *ofile)
{
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(ofile);
        if (root->isLeaf())
            printTaxa(out, root->neighbors[0]->node);
        else
            printTaxa(out);
        out.close();
        cout << "Taxa list was printed to " << ofile << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, ofile);
    }
}

void MTree::printTaxa(ostream &out, Node *node, Node *dad)
{
    if (!node) node = root;
    if (node->isLeaf())
        out << node->name << endl;
    else
    {
        // internal node
        //for (int i = 0; i < node->neighbors.size(); i++)
        //if (node->neighbors[i]->node != dad)
        FOR_NEIGHBOR_IT(node, dad, it)	{
            printTaxa(out, (*it)->node, node);
        }
    }
}

void MTree::printTaxa(ostream &out, NodeVector &subtree) {
    for (int i = 0; i < leafNum; i++)
        if (subtree[i] != NULL) {
            out << subtree[i]->name << endl;
        }
}

void MTree::readTree(const char *infile, bool &is_rooted) {
    ifstream in;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        readTree(in, is_rooted);
        in.close();
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
    rooted = is_rooted;
    if (verbose_mode >= VB_MED) {
        cout << "Tree contains " << leafNum - is_rooted <<
             " taxa and " << nodeNum-1-is_rooted << " branches" << (is_rooted ? " (rooted)" : "") << endl;
    }
}

/*
void MTree::readTreeString(string tree_string, bool is_rooted) {
	stringstream str;
	str << tree_string;
	str.seekg(0, ios::beg);
	freeNode();
	readTree(str, is_rooted);
}
*/


void MTree::readTree(istream &in, bool &is_rooted)
{
    in_line = 1;
    in_column = 1;
    in_comment = "";
    try {
        char ch;
        ch = readNextChar(in);
        if (ch != '(') {
        	cout << in.rdbuf() << endl;
            throw "Tree file does not start with an opening-bracket '('";
        }

        leafNum = 0;

        DoubleVector branch_len;
        Node *node;
        parseFile(in, ch, node, branch_len);
        // 2018-01-05: assuming rooted tree if root node has two children
        if (is_rooted || (!branch_len.empty() && branch_len[0] != 0.0) || node->degree() == 2) {
            if (branch_len.empty())
                branch_len.push_back(-1.0);
            if (branch_len[0] == -1.0) branch_len[0] = 0.0;
            if (branch_len[0] < 0.0)
                throw ERR_NEG_BRANCH;
            rooted = is_rooted = true;
            root = newNode(leafNum, ROOT_NAME);
            root->addNeighbor(node, branch_len);
            node->addNeighbor(root, branch_len);
            leafNum++;
            rooted = true;
        } else { // assign root to one of the neighbor of node, if any
            FOR_NEIGHBOR_IT(node, NULL, it)
            if ((*it)->node->isLeaf()) {
                root = (*it)->node;
                break;
            }
        }
        // make sure that root is a leaf
        ASSERT(root->isLeaf());

        if (in.eof() || ch != ';')
            throw "Tree file must be ended with a semi-colon ';'";
    } catch (bad_alloc) {
        outError(ERR_NO_MEMORY);
    } catch (const char *str) {
        outError(str, reportInputInfo());
    } catch (string str) {
        outError(str.c_str(), reportInputInfo());
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, reportInputInfo());
    } catch (...) {
        // anything else
        outError(ERR_READ_ANY, reportInputInfo());
    }

    nodeNum = leafNum;
    initializeTree();

    //bool stop = false;
    //checkValidTree(stop);
}

void MTree::initializeTree(Node *node, Node* dad)
{
    if (!node) {
        node = root;
        nodeNum = leafNum;
        branchNum = 0;
    }
    if (!node->isLeaf())
    {
        node->id = nodeNum;
        nodeNum++;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        (*it)->id = branchNum;
        (*it)->node->findNeighbor(node)->id = branchNum;
        branchNum++;
        initializeTree((*it)->node, node);
    }
}

void MTree::parseBranchLength(string &lenstr, DoubleVector &branch_len) {
//    branch_len.push_back(convert_double(lenstr.c_str()));
    double len = convert_double(lenstr.c_str());
    if (in_comment.empty()) {
        branch_len.push_back(len);
        return;
    }
    convert_double_vec(in_comment.c_str(), branch_len, BRANCH_LENGTH_SEPARATOR);
//    char* str = (char*)in_comment.c_str() + 1;
//    int pos;
//    for (int i = 1; str[0] == 'L'; i++) {
//        str++;
//        int id = convert_int(str, pos);
//        if (id != i)
//            throw "Wrong ID in " + string(str);
//        if (str[pos] != '=')
//            throw "= is expected in " + string(str);
//        str += pos+1;
//        double val = convert_double(str, pos);
//        branch_len.push_back(val);
//        if (str[pos] == ',') {
//            str += pos+1;
//            continue;
//        } else
//            break;
//    }
}


void MTree::parseFile(istream &infile, char &ch, Node* &root, DoubleVector &branch_len)
{
    Node *node;
    int maxlen = 1000;
    string seqname;
    int seqlen;
    DoubleVector brlen;
    branch_len.clear();

    root = newNode();

    if (ch == '(') {
        // internal node
        ch = readNextChar(infile);
        while (ch != ')' && !infile.eof())
        {
            node = NULL;
            parseFile(infile, ch, node, brlen);
            //if (brlen == -1.0)
            //throw "Found branch with no length.";
            //if (brlen < 0.0)
            //throw ERR_NEG_BRANCH;
            root->addNeighbor(node, brlen);
            node->addNeighbor(root, brlen);
            if (infile.eof())
                throw "Expecting ')', but end of file instead";
            if (ch == ',')
                ch = readNextChar(infile);
            else if (ch != ')') {
                string err = "Expecting ')', but found '";
                err += ch;
                err += "' instead";
                throw err;
            }
        }
        if (!infile.eof()) ch = readNextChar(infile);
    }
    // now read the node name
    seqlen = 0;
    char end_ch = 0;
    if (ch == '\'' || ch == '"') end_ch = ch;
    seqname = "";

    while (!infile.eof() && seqlen < maxlen)
    {
        if (end_ch == 0) {
            if (is_newick_token(ch) || controlchar(ch)) break;
        }
        seqname += ch;
        seqlen++;
//        seqname[seqlen++] = ch;
        ch = infile.get();
        in_column++;
        if (end_ch != 0 && ch == end_ch) {
            seqname += ch;
            seqlen++;
//            seqname[seqlen++] = ch;
            break;
        }
    }
    if ((controlchar(ch) || ch == '[' || ch == end_ch) && !infile.eof())
        ch = readNextChar(infile, ch);
    if (seqlen == maxlen)
        throw "Too long name ( > 1000)";
    if (root->isLeaf())
        renameString(seqname);
//    seqname[seqlen] = 0;
    if (seqlen == 0 && root->isLeaf())
        throw "Redundant double-bracket ‘((…))’ with closing bracket ending at";
    if (seqlen > 0)
        root->name.append(seqname);
    if (root->isLeaf()) {
        // is a leaf, assign its ID
        root->id = leafNum;
        if (leafNum == 0)
            MTree::root = root;
        leafNum++;
    }

    if (ch == ';' || infile.eof())
        return;
    if (ch == ':')
    {
        string saved_comment = in_comment;
        ch = readNextChar(infile);
        if (in_comment.empty())
            in_comment = saved_comment;
        seqlen = 0;
        seqname = "";
        while (!is_newick_token(ch) && !controlchar(ch) && !infile.eof() && seqlen < maxlen)
        {
//            seqname[seqlen] = ch;
            seqname += ch;
            seqlen++;
            ch = infile.get();
            in_column++;
        }
        if ((controlchar(ch) || ch == '[') && !infile.eof())
            ch = readNextChar(infile, ch);
        if (seqlen == maxlen || infile.eof())
            throw "branch length format error.";
//        seqname[seqlen] = 0;
        parseBranchLength(seqname, branch_len);
//        convert_double_vec(seqname.c_str(), branch_len, BRANCH_LENGTH_SEPARATOR);
    }
}



/**
	check tree is bifurcating tree (every leaf with level 1 or 3)
*/
void MTree::checkValidTree(bool &stop, Node *node, Node *dad)
{
    if (!node) node = root;
    if (node->degree() != 1 && node->degree() != 3) {
        cout << "Tree is not bifurcating." << endl;
        stop = true;
        return;
    }
    //for (int i = 0; i < node->neighbors.size(); i++)
    //if (node->neighbors[i]->node != dad) {
    FOR_NEIGHBOR_IT(node, dad, it) {
        checkValidTree(stop, (*it)->node, node);
        if (stop)
            return;
    }
}

double MTree::treeLength(Node *node, Node *dad)
{
    if (!node) node = root;
    double sum = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        sum += (*it)->length + treeLength((*it)->node, node);
    }
    return sum;
}

double MTree::treeLengthInternal( double epsilon, Node *node, Node *dad)
{
    if (!node) node = root;
    double sum = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
    	if (!(*it)->node->isLeaf() && !node->isLeaf())
    	{
    		if (treeLength((*it)->node, node) > epsilon) {
    			sum += (*it)->length + treeLengthInternal(epsilon, (*it)->node, node);
    		}
    	}
    	else {
    		if (treeLength((*it)->node, node) > epsilon) {
    			sum += treeLengthInternal(epsilon, (*it)->node, node);
    		}
    	}
    }
    return sum;
}

double MTree::treeDepth(Node *node, Node *dad)
{
    if (!node) node = root;
    double maxsum = 0.0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        double len = (*it)->length;
        if (len < 0.0) len = 0.0;
        double sum = len + treeDepth((*it)->node, node);
        if (sum > maxsum) maxsum = sum;
    }
    return maxsum;
}

void MTree::getNonCherryLeaves(NodeVector &noncherry, NodeVector &cherry, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
    	if (node->isInCherry()) {
    		cherry.push_back(node);
    	} else {
            noncherry.push_back(node);
    	}
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
    	getNonCherryLeaves(noncherry, cherry, (*it)->node, node);
    }
}

void MTree::getTaxa(NodeVector &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        taxa.push_back(node);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxa(taxa, (*it)->node, node);
    }
}

void MTree::getAllNodesInSubtree(Node *node, Node *dad, NodeVector &nodeList) {
    ASSERT(node);
    nodeList.push_back(node);
    if (node->isLeaf()) {
        return;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        getAllNodesInSubtree((*it)->node, node, nodeList);
    }
}

int MTree::getNumTaxa(Node *node, Node *dad) {
    int numLeaf = 0;
    if (!node) {
    	node = root;
    	numLeaf = 1;
    } else {
        if (node->isLeaf()) {
            return 1;
        }
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        numLeaf += getNumTaxa((*it)->node, node);
    }
    return numLeaf;
}

void MTree::getInternalNodes(NodeVector &nodes, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!(*it)->node->isLeaf()) {
        getInternalNodes(nodes, (*it)->node, node);
        nodes.push_back((*it)->node);
    }
}

void MTree::getMultifurcatingNodes(NodeVector &nodes, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)    {
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!(*it)->node->isLeaf()) {
        if ((*it)->node->degree() > 3)
            nodes.push_back((*it)->node);
        getMultifurcatingNodes(nodes, (*it)->node, node);
    }
}

void MTree::generateNNIBraches(NodeVector &nodes1, NodeVector &nodes2, SplitGraph* excludeSplits, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!(*it)->node->isLeaf()) {
        generateNNIBraches(nodes1, nodes2, excludeSplits, (*it)->node, node);
        if (!node->isLeaf()) {
        	if (excludeSplits != NULL && excludeSplits->size() != 0) {
        		Split* sp = getSplit(node, (*it)->node);
        		if (excludeSplits->containSplit(*sp)) {
        			delete sp;
        			continue;
        		}
        		delete sp;
        	}
			if (node->id < (*it)->node->id) {
				nodes1.push_back(node);
				nodes2.push_back((*it)->node);
			} else {
				nodes1.push_back((*it)->node);
				nodes2.push_back(node);
			}
        }
    }
}

//bool MTree::branchExist(Node* node1, Node* node2, NodeVector& nodes1, NodeVector& nodes2) {
//	assert(nodes1.size() == nodes2.size());
//	bool existed = false;
//	for (int i = 0; i < nodes1.size(); i++) {
//		if (nodes1[i] == node1) {
//			if (nodes2[i] == node2) {
//				existed = true;
//				break;
//			}
//		}
//		if (nodes1[i] == node2) {
//			if (nodes2[i] == node1) {
//				existed = true;
//				break;
//			}
//		}
//	}
//	return existed;
//}

void MTree::getSurroundingInnerBranches(Node *node, Node *dad, int depth, Branches &surrBranches) {
    if (depth == 0)
      return;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->node->isLeaf()) {
            Branch curBranch;
            curBranch.first = node;
            curBranch.second = (*it)->node;
            int branchID = pairInteger(node->id, (*it)->node->id);
            if (surrBranches.find(branchID) == surrBranches.end())
                surrBranches.insert(pair<int,Branch>(branchID, curBranch));
            getSurroundingInnerBranches((*it)->node, node, depth-1, surrBranches);
        }
    }
}

bool MTree::isInnerBranch(Node* node1, Node* node2) {
    return(node1->degree() >= 3 && node2->degree() >= 3 && isABranch(node1, node2));
}

bool MTree::isABranch(Node* node1, Node* node2) {
    return (node1->findNeighbor(node2) != NULL && node2->findNeighbor(node1) != NULL);
}

void MTree::getBranches(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad, bool post_traversal) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)   {
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!post_traversal) {
            if (node->id < (*it)->node->id) {
                nodes.push_back(node);
                nodes2.push_back((*it)->node);
            } else {
                nodes.push_back((*it)->node);
                nodes2.push_back(node);
            }
        }
        getBranches(nodes, nodes2, (*it)->node, node, post_traversal);
        if (post_traversal) {
            if (node->id < (*it)->node->id) {
                nodes.push_back(node);
                nodes2.push_back((*it)->node);
            } else {
                nodes.push_back((*it)->node);
                nodes2.push_back(node);
            }
        }
    }
}

void MTree::getBranches(int max_dist, NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad) {
    if (!node) node = root;
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)   {
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (node->id < (*it)->node->id) {
            nodes.push_back(node);
            nodes2.push_back((*it)->node);
        } else {
            nodes.push_back((*it)->node);
            nodes2.push_back(node);
        }
        if (max_dist > 1)
            getBranches(max_dist-1, nodes, nodes2, (*it)->node, node);
    }
}

void MTree::getBranches(BranchVector& branches, Node *node, Node *dad, bool post_traversal) {
    if (!node) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!post_traversal) {
            Branch branch;
            branch.first = node;
            branch.second = (*it)->node;
            branches.push_back(branch);
        }
        getBranches(branches, (*it)->node, node, post_traversal);
        if (post_traversal) {
            Branch branch;
            branch.first = node;
            branch.second = (*it)->node;
            branches.push_back(branch);
        }
    }
}

void MTree::getInnerBranches(Branches& branches, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
    	if (isInnerBranch((*it)->node, node)) {
            Branch branch;
            branch.first = node;
            branch.second = (*it)->node;
            branches.insert(pair<int, Branch>(pairInteger(branch.first->id, branch.second->id), branch));
    	}
    	getInnerBranches(branches, (*it)->node, node);
    }
}

void MTree::getInnerBranches(BranchVector& branches, Node *node, Node *dad, bool post_traversal) {
    if (!node) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!node->isLeaf() && !(*it)->node->isLeaf() && !post_traversal) {
            Branch branch;
            branch.first = node;
            branch.second = (*it)->node;
            branches.push_back(branch);
        }
        getInnerBranches(branches, (*it)->node, node, post_traversal);
        if (!node->isLeaf() && !(*it)->node->isLeaf() && post_traversal) {
            Branch branch;
            branch.first = node;
            branch.second = (*it)->node;
            branches.push_back(branch);
        }
    }
}

void MTree::getBranchLengths(vector<DoubleVector> &len, Node *node, Node *dad) {
    if (!node) {
        node = root;
        ASSERT(len.size() == branchNum);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        (*it)->getLength(len[(*it)->id]);
        getBranchLengths(len, (*it)->node, node);
    }
}

void MTree::setBranchLengths(vector<DoubleVector> &len, Node *node, Node *dad) {
    if (!node) {
        node = root;
        ASSERT(len.size() == branchNum);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        (*it)->setLength(len[(*it)->id]);
        (*it)->node->findNeighbor(node)->setLength(len[(*it)->id]);
        setBranchLengths(len, (*it)->node, node);
    }
}

void MTree::getOrderedTaxa(NodeVector &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        if (taxa.empty()) taxa.resize(leafNum);
        taxa[node->id] = node;
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getOrderedTaxa(taxa, (*it)->node, node);
    }
}

void MTree::getTaxaName(vector<string> &taxname, Node *node, Node *dad) {
    if (!node) {
        node = root;
    }
    if (node->isLeaf()) {
        if (taxname.empty()) {
            taxname.resize(leafNum);
        }
        taxname[node->id] = node->name;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxaName(taxname, (*it)->node, node);
    }
}

void MTree::getMapOfTaxonNameToNode(Node* node, Node* dad
                                    , map<string, Node*> &map) {
    if (!node) {
        node = root;
    }
    if (node->isLeaf()) {
        map[node->name] = node;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        getMapOfTaxonNameToNode((*it)->node, node, map);
    }
}

void MTree::getArrayOfTaxaNodesById(Node* node, Node* dad,
                             NodeVector& array) {
    if (!node) {
        node = root;
    }
    if (node->isLeaf() && 0<node->id) {
        if (array.size() <= node->id) {
            array.resize(node->id+1, nullptr);
        }
        array[node->id] = node;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        getArrayOfTaxaNodesById((*it)->node, node, array);
    }
}


void MTree::getNodeName(set<string> &nodename, Node *node, Node *dad) {
    if (!node) node = root;
    if (!node->name.empty())
        nodename.insert(node->name);
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)    {
    FOR_NEIGHBOR_IT(node, dad, it) {
        getNodeName(nodename, (*it)->node, node);
    }
}

void MTree::getUnorderedTaxaName(vector<string> &taxname, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
    	taxname.push_back(node->name);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getUnorderedTaxaName(taxname, (*it)->node, node);
    }
}

void MTree::getTaxaID(vector<int> &taxa, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        taxa.push_back(node->id);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        getTaxaID(taxa, (*it)->node, node);
    }
}

bool MTree::containsSplits(SplitGraph& splits) {
	SplitGraph treeSplits;
	convertSplits(treeSplits);
	//check if treeSplits contains all splits in splits
	for (SplitGraph::iterator it = splits.begin(); it != splits.end(); it++) {
		if (!treeSplits.containSplit(**it))
			return false;
	}
	//treeSplits.report(cout);
	//splits.report(cout);
	return true;
}

Split* MTree::getSplit(Node* node1, Node* node2) {
    Neighbor* node12 = node1->findNeighbor(node2);
    return node12->split;
}

Split* MTree::_getSplit(Node* node1, Node* node2) {
    Split* sp = new Split(leafNum);
    getTaxa(*sp, node1, node2);
    if (sp->shouldInvert())
        sp->invert();
    return sp;
}

void MTree::convertSplits(SplitGraph &sg, Split *resp, NodeVector *nodes, Node *node, Node *dad) {
    if (!node) node = root;
    ASSERT(resp->getNTaxa() == leafNum);
    bool has_child = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
        //vector<int> taxa;
        //getTaxaID((*it)->node, node, taxa);

        Split *sp = new Split(leafNum, (*it)->length);
        convertSplits(sg, sp, nodes, (*it)->node, node);
        *resp += *sp;
        if (sp->shouldInvert())
            sp->invert();
		 /* ignore nodes with degree of 2 because such split will be added before */
        if (node->degree() != 2) {
        	sg.push_back(sp);
        	if (nodes) nodes->push_back((*it)->node);
        }
        has_child = true;
    }
    if (!has_child)
        resp->addTaxon(node->id);
}

void MTree::convertSplits(vector<string> &taxname, SplitGraph &sg, NodeVector *nodes, Node *node, Node *dad) {
    if (!sg.taxa) {
        sg.taxa = new NxsTaxaBlock();
        for (vector<string>::iterator it = taxname.begin(); it != taxname.end(); it++)
            sg.taxa->AddTaxonLabel(NxsString(it->c_str()));
    }
    if (!sg.splits)
        sg.splits = new MSplitsBlock(&sg);
    if (!sg.pda)
        sg.pda = new MPdaBlock(&sg);

    // make the cycle
    getTaxaID(sg.splits->cycle);
    // make the splits
    Split sp(leafNum);
    convertSplits(sg, &sp, nodes, node, dad);
}

void MTree::convertSplits(SplitGraph &sg, NodeVector *nodes, Node *node, Node *dad) {

    // make the taxa name
    vector<string> taxname;
    taxname.resize(leafNum);
    getTaxaName(taxname);

    convertSplits(taxname, sg, nodes, node, dad);
}

inline int splitnumtaxacmp(const Split* a, const Split* b)
{
    return (a->countTaxa() < b->countTaxa());
}

void MTree::convertToTree(SplitGraph &sg) {
    SplitGraph::iterator it;
    int taxid;
    int count;
	BoolVector has_tax;
	has_tax.resize(sg.getNTaxa(), false);
	// first add trivial splits if not existed
	for (it = sg.begin(); it != sg.end(); it++) {
		taxid = (*it)->trivial();
		if (taxid >= 0) has_tax[taxid] = true;
	}
	for (count = 0; count < has_tax.size(); count++)
		if (!has_tax[count]) {
			Split *sp = new Split(sg.getNTaxa());
			sp->addTaxon(count);
			sg.push_back(sp);
		}
    // sort splits by the number of taxa they contain
    sort(sg.begin(), sg.end(), splitnumtaxacmp);

    // initialize the tree
    rooted = false;
    leafNum = sg.getNTaxa();
    nodeNum = leafNum;

    // create the ground nodes, first as the leaves
    NodeVector leaves;
    vector<Split*> cladetaxa;
    leaves.resize(leafNum, NULL);
    cladetaxa.resize(leafNum, NULL);
    // first add all trivial splits into tree
    for (it = sg.begin(), count = 0; it != sg.end(); it++, count++) {
        //(*it)->report(cout);
        taxid = (*it)->trivial();
        if (taxid < 0) break;
        ASSERT(leaves[taxid] == NULL);
        leaves[taxid] = newNode(taxid, sg.getTaxa()->GetTaxonLabel(taxid).c_str());
        leaves[taxid]->addNeighbor(NULL, (*it)->getWeight());
        cladetaxa[taxid] = (*it);
    }
    // now fill in all missing taxa with zero terminal branch
    for (taxid = 0; taxid < leafNum; taxid++)
        ASSERT(leaves[taxid]);

    // now add non-trivial splits, cotinue with the interrupted iterator
    for (/*it = sg.begin()*/; it != sg.end(); it++) {
        //(*it)->report(cout);
        Split *mysp = *it;
        Node *newnode = newNode(nodeNum);
        int count = 0;

        for (taxid = 0; taxid < leaves.size(); )
            if (cladetaxa[taxid]->subsetOf(*mysp)) // clade is a subset of current split
            {
                count += cladetaxa[taxid]->countTaxa();
                double len = leaves[taxid]->updateNeighbor(NULL, newnode);
                newnode->addNeighbor(leaves[taxid], len);
                leaves[taxid] = leaves.back();
                leaves.pop_back();
                cladetaxa[taxid] = cladetaxa.back();
                cladetaxa.pop_back();
            } else taxid++;
        ASSERT(count == mysp->countTaxa());
        cladetaxa.push_back(mysp);
        leaves.push_back(newnode);

        newnode->addNeighbor(NULL, mysp->getWeight());
        nodeNum++;
    }
    ASSERT(leaves.size() >= 3);
    Node *newnode = newNode(nodeNum);
    for (taxid = 0; taxid < leaves.size(); taxid++) {
        double len = leaves[taxid]->updateNeighbor(NULL, newnode);
        newnode->addNeighbor(leaves[taxid], len);
    }
    root = newnode;
    nodeNum++;
    cladetaxa.clear();
    string root_name = ROOT_NAME;
    newnode = findLeafName(root_name);
    if (newnode) {
        rooted = true;
        root = newnode;
    }
    
}

Node *MTree::findNodeName(string &name, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->name == name) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findNodeName(name, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}

bool MTree::findNodeNames(unordered_set<string> &taxa_set, pair<Node*,Neighbor*> &res, Node *node, Node *dad) {
    int presence = 0;
    Neighbor *target = NULL;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->node->isLeaf()) {
            if (taxa_set.find((*it)->node->name) != taxa_set.end()) {
                presence++;
            } else target = (*it);
        } else {
            if (findNodeNames(taxa_set, res, (*it)->node, node))
                presence++;
            else
                target = *it;
            if (res.first)
                return false;
        }
    }
    // all presence or absence
    if (presence == node->neighbors.size()-1)
        return true;
    if (presence == 0)
        return false;
    // inbetween: detect it!
    res.first = node;
    res.second = target;
    if (target != node->neighbors[0]) {
        // move target into the first neighbor
        FOR_NEIGHBOR_IT(node, NULL, it)
        if ((*it) == target) {
            (*it) = node->neighbors[0];
            node->neighbors[0] = target;
            break;
        }
    }
    return false;
}

Node *MTree::findLeafName(const string &name, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf() && node->name == name) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findLeafName(name, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}

Node *MTree::findNodeID(int id, Node *node, Node* dad) {
    if (!node) node = root;
    if (node->id == id) return node;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *res = findNodeID(id, (*it)->node, node);
        if (res) return res;
    }
    return NULL;
}


void MTree::scaleLength(double norm, bool make_int, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_DECLARE(node, NULL, it) {
        (*it)->length *= norm;
        if (make_int)
            (*it)->length = round((*it)->length);
    }

    FOR_NEIGHBOR(node, dad, it) {
        scaleLength(norm, make_int, (*it)->node, node);
    }
}

void MTree::transformBranchLenRAX(double factor, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_DECLARE(node, NULL, it) {
        (*it)->length /= factor;
        (*it)->length = exp(-(*it)->length);
    }

    FOR_NEIGHBOR(node, dad, it) {
    	transformBranchLenRAX(factor, (*it)->node, node);
    }
}

void MTree::scaleCladeSupport(double norm, bool make_int, Node *node, Node *dad) {
    if (!node) node = root;
    if (!node->isLeaf() && !node->name.empty()) {
        double supp = 0.0;
        try {
            supp = convert_double(node->name.c_str());
        } catch (string str) {
            outError(str);
        }
        supp *= norm;
        if (make_int)
            supp = round(supp);
        node->name = "";
        node->name += supp;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        scaleCladeSupport(norm, make_int, (*it)->node, node);
    }
}


MTree::~MTree()
{
    if (root != NULL)
        freeNode();
    root = NULL;
}

int MTree::freeNode(Node *node, Node *dad)
{
	if ( root == NULL )
		return 0;
    if (!node) node = root;
    NeighborVec::reverse_iterator it;
    int num_nodes = 1;
    for (it = node->neighbors.rbegin(); it != node->neighbors.rend(); it++)
        if ((*it)->node != dad) {
            num_nodes += freeNode((*it)->node, node);
        }
    delete node;
    return num_nodes;
}

char MTree::readNextChar(istream &in, char current_ch) {
    char ch;
    if (current_ch == '[')
        ch = current_ch;
    else {
        in.get(ch);
        in_column++;
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
    }
    while (controlchar(ch) && !in.eof()) {
        in.get(ch);
        in_column++;
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
    }
    in_comment = "";
    // ignore comment
    while (ch=='[' && !in.eof()) {
        while (ch!=']' && !in.eof()) {
            in.get(ch);
            if (ch != ']')
                in_comment += ch;
            in_column++;
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
        if (ch != ']') throw "Comments not ended with ]";
        in_column++;
        in.get(ch);
        if (ch == 10) {
            in_line++;
            in_column = 1;
        }
        while (controlchar(ch) && !in.eof()) {
            in_column++;
            in.get(ch);
            if (ch == 10) {
                in_line++;
                in_column = 1;
            }
        }
    }
    return ch;
}

string MTree::reportInputInfo() {
    string str = " (line ";
    str += convertIntToString(in_line) + " column " + convertIntToString(in_column-1) + ")";
    return str;
}

void MTree::convertToUnrooted() {
    ASSERT(rooted && root);
    ASSERT(root->isLeaf() && root->id == leafNum-1);
    Node *node = root->neighbors[0]->node;
    Node *taxon = findFirstTaxon();
    
    rooted = false;
    leafNum--;
    
    // delete root node
    if (node->degree() == 3) {
        // delete and join adjacent branches
        Node *node1 = NULL, *node2 = NULL;
        double len = 0.0;
        FOR_NEIGHBOR_IT(node, root, it) {
            if (!node1) node1 = (*it)->node; else node2 = (*it)->node;
            len += (*it)->length;
        }
        node1->updateNeighbor(node, node2, len);
        node2->updateNeighbor(node, node1, len);
        delete node;
    } else {
        // only delete root node
        auto it = node->findNeighborIt(root);
        delete *it;
        node->neighbors.erase(it);
        
    }
    
    delete root;
    // set a temporary taxon so that tree traversal works
    root = taxon;
    
    initializeTree();
    //    computeBranchDirection();
}

typedef map<int, Neighbor*> IntNeighborMap;

int MTree::sortTaxa(Node *node, Node *dad) {
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    }
    if (node->isLeaf())
        return node->id;
    IntNeighborMap taxid_nei_map;
    FOR_NEIGHBOR_IT(node, dad, it) {
        int taxid = sortTaxa((*it)->node, node);
        taxid_nei_map.insert(IntNeighborMap::value_type(taxid, (*it)));
    }
    ;
    int i = 0;
    for (IntNeighborMap::iterator it = taxid_nei_map.begin(); it != taxid_nei_map.end(); it++, i++) {
        if (node->neighbors[i]->node == dad) i++;
        node->neighbors[i] = it->second;
    }

    return taxid_nei_map.begin()->first;
}

void MTree::setExtendedFigChar() {
	//fig_char[0] = 179;
	//fig_char[1] = 196;
	fig_char[2] = '/';
	//fig_char[3] = 195;
	fig_char[4] = '\\';
}

void MTree::drawTree(ostream &out, int brtype, double zero_epsilon) {
    IntVector sub_tree_br;
    if (verbose_mode >= VB_DEBUG) {
        printTree(cout);
        cout << endl;
    }
    Node *node = root;
    if (node->isLeaf()) node = node->neighbors[0]->node;
    double scale = 60.0/treeDepth(node);
    //if (verbose_mode >= VB_DEBUG)
    //cout << "Tree depth: " << scale<< endl;
    drawTree2(out, brtype, scale, sub_tree_br, zero_epsilon);
    /*
    if (brtype & WT_INT_NODE)
        drawTree2(out, brtype, scale, sub_tree_br, zero_epsilon);
    else
        drawTree(out, brtype, scale, sub_tree_br, zero_epsilon);
    */
    out << endl;
}

/*
void MTree::drawTree(ostream &out, int brtype, double brscale, IntVector &subtree_br, double zero_epsilon, Node *node, Node *dad) {
    int i, br_len = 3;
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    } else {

        if (brtype & WT_BR_SCALE) {
            br_len = floor(node->findNeighbor(dad)->length * brscale)-1;
            if (br_len < 3) br_len = 3;
            //if (!node->isLeaf() && br_len < 4) br_len = 4;
        }
        out << '+';
        if ((brtype & WT_INT_NODE) && !node->isLeaf()) {
            string str = convertIntToString(node->id);
            for (i = 0; i < br_len-str.length(); i++) out << '-';
            out << node->id;
        } else
            for (i = 0; i < br_len; i++) out << '-';
    }
    if (node->isLeaf()) {
        out << node->name;
        if (brtype & WT_TAXON_ID)
            out << " (" << node->id << ")";
        out << endl;
        return;
    }
    int descendant_cnt = node->degree();
    if (dad) descendant_cnt--;
    int cnt = 0;
    subtree_br.push_back(br_len);
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (cnt == descendant_cnt-1)
            subtree_br.back() = -subtree_br.back();

        drawTree(out, brtype, brscale, subtree_br, zero_epsilon, (*it)->node, node);
        cnt++;
        if (cnt == descendant_cnt) break;
        for (IntVector::iterator it = subtree_br.begin()+1; it != subtree_br.end(); it++)
        {
            if ((*(it-1)) > 0) out << '|';
            else out << ' ';
            for (i = 0; i < abs(*it); i++) out << ' ';
        }
    }
    subtree_br.pop_back();
}
*/

void MTree::drawTree2(ostream &out, int brtype, double brscale, IntVector &subtree_br, double zero_epsilon, Node *node, Node *dad) {
    int i, br_len = 3;
    IntVector::iterator ii;
    bool zero_length = false;

    //cout << "DrawTree2!" << endl;
    if (!node) {
        node = root;
        if (node->isLeaf()) node = node->neighbors[0]->node;
    } else {
        if (brtype & WT_BR_SCALE) {
            br_len = floor(node->findNeighbor(dad)->length * brscale)-1;
            if (br_len < 2) br_len = 2;
        }
        if (node->findNeighbor(dad)->length <= zero_epsilon) zero_length = true;
    }
    if (node->isLeaf()) {
        for (ii = subtree_br.begin()+1; ii != subtree_br.end(); ii++) {
            if (abs(*(ii-1)) > 1000) out << ' ';
            else out << fig_char[0];
            int num = abs(*ii);
            if (num > 1000) num -= 1000;
            for (i = 0; i < num; i++) out << ' ';
        }
        out << ((node==dad->neighbors.front()->node) ? fig_char[2] : ((node==dad->neighbors.back()->node) ? fig_char[4] : fig_char[3]));
        for (i = 0; i < br_len; i++)
            out << ((zero_length) ? '*' : fig_char[1]);
        out << node->name;
        if (brtype & WT_TAXON_ID)
            out << " (" << node->id << ")";
        if (brtype & WT_BR_ID)
            out << " [" << node->neighbors[0]->id << "]";
        if (brtype & WT_BR_LEN)
            out << " " << node->neighbors[0]->length;
        //out << " ";
        //copy (subtree_br.begin(), subtree_br.end(), ostream_iterator<int> (out, " "));
        out << endl;
        return;
    }
    int descendant_cnt = node->degree();
    if (dad) descendant_cnt--;
    int cnt = 0;
    bool first = true;

    br_len = br_len+1000;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (cnt == descendant_cnt-1)
            br_len = -br_len;
        subtree_br.push_back(br_len);

        drawTree2(out, brtype, brscale, subtree_br, zero_epsilon, (*it)->node, node);
        subtree_br.pop_back();
        if (br_len > 1000) br_len -= 1000;
        cnt++;
        if (cnt == descendant_cnt) break;
        if (subtree_br.size() > 1)
            for (ii = subtree_br.begin()+1; ii != subtree_br.end(); ii++) {
                if (abs(*(ii-1)) > 1000) out << ' ';
                else out << fig_char[0];
                if (ii == subtree_br.begin()) continue;
                int num = abs(*ii);
                if (num > 1000) num -= 1000;
                for (i = 0; i < num; i++) out << ' ';
            }
        if (first) {
            if (dad) {
				out << ((node==dad->neighbors.front()->node) ? fig_char[2] : ((node==dad->neighbors.back()->node) ? fig_char[4] : fig_char[3]));
                for (i = 0; i < abs(br_len); i++)
                    out << ((zero_length) ? '*' : fig_char[1]);
            }
            if (brtype & WT_INT_NODE)
            	out << node->id;
            else
            	out << fig_char[0];
            if (!node->name.empty())
                out << " (" << node->name << ")";
            if (brtype & WT_BR_LEN && dad)
                out << " " << node->findNeighbor(dad)->length;
            if (brtype & WT_BR_ID && dad)
                out << " [" << node->findNeighbor(dad)->id << "]";
            if (!subtree_br.empty()) {
                if (subtree_br.back() >1000)
                    subtree_br.back() -= 1000;
                else if (subtree_br.back() < 0)
                    subtree_br.back() -= 1000;
            }
        } else {
            if (dad) {
                if (abs(subtree_br.back()) > 1000) out << ' ';
                else out << fig_char[0];
                for (i = 0; i < abs(br_len); i++)
                    out << ' ';
            }
            out << fig_char[0];
        }
        //out << " ";
        //copy (subtree_br.begin(), subtree_br.end(), ostream_iterator<int> (out, " "));
        out << endl;
        first = false;
    }
}

bool MTree::equalTopology(MTree *tree) {
	ASSERT(root->isLeaf());
	Node *root2 = tree->findLeafName(root->name);
	if (!root2) return false;
	ostringstream ostr, ostr2;
	printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
	tree->printTree(ostr2, WT_TAXON_ID | WT_SORT_TAXA, root2);
	return ostr.str() == ostr2.str();
}

void MTree::calcDist(char *filename) {
    vector<string> taxname;
    int i, j;

    // allocate memory
    taxname.resize(leafNum);
    double *dist = new double [leafNum * leafNum];
    // calculate the distances
    calcDist(dist);
    // get the taxa name
    getTaxaName(taxname);

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);

        // now write the distances in phylip .dist format
        out << leafNum << endl;

        for (i = 0; i < leafNum; i++) {
            out << taxname[i] << "   ";
            for (j = 0; j < leafNum; j++) {
                out << dist[i*leafNum + j] << "  ";
            }
            out << endl;
        }
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    delete [] dist;
}

void MTree::calcDist(double* &dist, Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        calcDist(node, 0.0, dist, node, NULL);
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
        calcDist(dist, (*it)->node, node);
    }
}

void MTree::calcDist(Node *aroot, double cur_len, double* &dist, Node *node, Node *dad) {
    double branch_length;
	if (!node) node = root;
    if (node->isLeaf()) {
        dist[aroot->id * leafNum + node->id] = cur_len;
        dist[node->id * leafNum + aroot->id] = cur_len;
    }
    //for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++)
    //if ((*it)->node != dad)	{
    FOR_NEIGHBOR_IT(node, dad, it) {
    	branch_length = (*it)->length;
        calcDist(aroot, cur_len + branch_length, dist, (*it)->node, node);
    }

}


/*********************************************
	class PDTaxaSet
*********************************************/

void PDTaxaSet::setSubTree(MTree &tree, NodeVector &subtree) {
    stringstream ostr;
    tree.printSubTree(ostr, subtree);
    tree_str = ostr.str();
}

void PDTaxaSet::setTree(MTree &tree) {
    // assign the taxa set
    tree.getTaxa(*this);
    // assign the score
    score = tree.treeLength();

    // assign tree_str
    stringstream ostr;
    tree.printTree(ostr);
    tree_str = ostr.str();
}


void PDTaxaSet::printTaxa(ostream &out) {
    for (iterator it = begin(); it != end(); it++)
        if ((*it)->name != ROOT_NAME)
            out << (*it)->name << endl;
}

void PDTaxaSet::printTaxa(char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        printTaxa(out);
        out.close();
        cout << "Taxa list was printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void PDTaxaSet::printTree(ostream &out) {
    if (!tree_str.empty())
        out << tree_str << endl;
}

void PDTaxaSet::printTree(char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        printTree(out);
        out.close();
        cout << "Tree was printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}


void PDTaxaSet::makeIDSet(int ntaxa, Split &id_set) {
    id_set.setNTaxa(ntaxa);
    id_set.setWeight(score);
    for (iterator it = begin(); it != end(); it++)
        id_set.addTaxon((*it)->id);
}

void MTree::writeInternalNodeNames(string &out_file) {
    try {
        ofstream out(out_file.c_str());
        NodeVector nodes;
        getInternalNodes(nodes);
        for (NodeVector::iterator nit = nodes.begin(); nit != nodes.end(); nit++) {
            out  << " " << (*nit)->name;
        }
        out << endl;
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, out_file);
    }
}

void MTree::assignLeafID(Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        node->id = atoi(node->name.c_str());
        ASSERT(node->id >= 0 && node->id < leafNum);
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    assignLeafID((*it)->node, node);
}

void MTree::assignLeafNameByID(Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
//        node->id = atoi(node->name.c_str());
//        assert(node->id >= 0 && node->id < leafNum);
        node->name = convertIntToString(node->id);
    }
    FOR_NEIGHBOR_IT(node, dad, it)
        assignLeafNameByID((*it)->node, node);
}

void MTree::getTaxa(Split &taxa, Node *node, Node *dad) {
	if (!node) node = root;
	if (node->isLeaf()) {
		taxa.addTaxon(node->id);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
		getTaxa(taxa, (*it)->node, node);
}


void MTree::extractQuadSubtrees(vector<Split*> &subtrees, BranchVector &branches, Node *node, Node *dad) {
	if (!node) node = root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		extractQuadSubtrees(subtrees, branches, (*it)->node, node);
		if ((*it)->node->isLeaf()) continue;
		// internal branch
		ASSERT(node->degree() == 3 && (*it)->node->degree() == 3);
		int cnt = 0;
		Node *child = (*it)->node;
        string treestrings[4];
        int nodeid = 0;
		FOR_NEIGHBOR_DECLARE(child, node, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, child);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
            if (Params::getInstance().print_df1_trees) {
                ostringstream ostr;
                printTree(ostr, WT_BR_LEN, (*it2)->node, child);
                treestrings[nodeid] = ostr.str();
                treestrings[nodeid] = treestrings[nodeid].substr(0, treestrings[nodeid].length()-1);
                nodeid++;
            }
		}
		FOR_NEIGHBOR(node, child, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, node);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
            if (Params::getInstance().print_df1_trees) {
                ostringstream ostr;
                printTree(ostr, WT_BR_LEN, (*it2)->node, node);
                treestrings[nodeid] = ostr.str();
                treestrings[nodeid] = treestrings[nodeid].substr(0, treestrings[nodeid].length()-1);
                nodeid++;
            }
		}
		ASSERT(cnt == leafNum);
        branches.push_back({node, child});
        if (Params::getInstance().print_df1_trees) {
            ASSERT(nodeid == 4);
            // output NNI-1 tree
            string treeDF1 = "(" + treestrings[0] + "," + treestrings[2] + ",(" + treestrings[1] + "," + treestrings[3] + "));";
            string treeDF2 = "(" + treestrings[0] + "," + treestrings[3] + ",(" + treestrings[1] + "," + treestrings[2] + "));";
            PUT_ATTR(child->findNeighbor(node), treeDF1);
            PUT_ATTR(child->findNeighbor(node), treeDF2);
        }
	}
}

/*
void MTree::assignBranchSupport(const char *trees_file, map<int,BranchSupportInfo> &branch_supports) {
	cout << "Reading input trees file " << trees_file << endl;
	try {
		ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees_file);
        assignBranchSupport(in, branch_supports);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, trees_file);
	}
}

void MTree::assignBranchSupport(istream &in, map<int,BranchSupportInfo> &branch_supports) {
	SplitGraph mysg;
	NodeVector mynodes;
	convertSplits(mysg, &mynodes, root->neighbors[0]->node);
	vector<Split*> subtrees;
	extractQuadSubtrees(subtrees, root->neighbors[0]->node);
	IntVector decisive_counts;
	decisive_counts.resize(mynodes.size(), 0);
	StrVector occurence_trees; // list of tree IDs where each split occurs
	if (verbose_mode >= VB_MED)
		occurence_trees.resize(mynodes.size());
	SplitGraph::iterator sit;
	for (sit = mysg.begin(); sit != mysg.end(); sit++)
		(*sit)->setWeight(0.0);
	int ntrees, taxid;
	for (ntrees = 1; !in.eof(); ntrees++) {
		MTree tree;
		bool is_rooted = false;

		// read in the tree and convert into split system for indexing
		tree.readTree(in, is_rooted);
		if (verbose_mode >= VB_DEBUG)
			cout << ntrees << " " << endl;
		StrVector taxname;
		tree.getTaxaName(taxname);
		// create the map from taxa between 2 trees
		Split taxa_mask(leafNum);
		for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
			taxid = mysg.findLeafName(*it);
			if (taxid < 0)
				outError("Taxon not found in full tree: ", *it);
			taxa_mask.addTaxon(taxid);
		}
		// make the taxa ordering right before converting to split system
		taxname.clear();
		int smallid;
		for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
			if (taxa_mask.containTaxon(taxid)) {
				taxname.push_back(mysg.getTaxa()->GetTaxonLabel(taxid));
				string name = (string)mysg.getTaxa()->GetTaxonLabel(taxid);
				tree.findLeafName(name)->id = smallid++;
			}
		ASSERT(taxname.size() == tree.leafNum);

		SplitGraph sg;
		//NodeVector nodes;
		tree.convertSplits(sg);
		SplitIntMap hash_ss;
		for (sit = sg.begin(); sit != sg.end(); sit++)
			hash_ss.insertSplit((*sit), 1);

		// now scan through all splits in current tree
		int id, qid;
		for (sit = mysg.begin(), id = 0, qid = 0; sit != mysg.end(); sit++, id++)
		if ((*sit)->trivial() < 0) // it is an internal split
		{

			bool decisive = true;
			for (int i = 0; i < 4; i++) {
				if (!taxa_mask.overlap(*subtrees[qid+i])) {
					decisive = false;
					break;
				}
			}
			qid += 4;
			if (!decisive) continue;

			decisive_counts[id]++;
			Split *subsp = (*sit)->extractSubSplit(taxa_mask);
			if (subsp->shouldInvert())
				subsp->invert();
			Split *sp = hash_ss.findSplit(subsp);
			if (sp && sp->trivial() < 0) {
				(*sit)->setWeight((*sit)->getWeight()+1.0);
				if (verbose_mode >= VB_MED)
					occurence_trees[id] += convertIntToString(ntrees) + " ";
				if (verbose_mode >= VB_MAX) {
					for (taxid = 0; taxid < (*sit)->getNTaxa(); taxid++)
						if ((*sit)->containTaxon(taxid))
							cout << " " << mysg.getTaxa()->GetTaxonLabel(taxid);
					cout << " --> ";
					for (taxid = 0; taxid < sp->getNTaxa(); taxid++)
						if (sp->containTaxon(taxid))
							cout << " " << taxname[taxid];
					cout << endl;
				}
			}
			delete subsp;
		}

		char ch;
		in.exceptions(ios::goodbit);
		(in) >> ch;
		if (in.eof()) break;
		in.unget();
		in.exceptions(ios::failbit | ios::badbit);

	}

	cout << ntrees << " trees read" << endl;

	for (int i = 0; i < mysg.size(); i++)
	if (!mynodes[i]->isLeaf())
	{
        BranchSupportInfo brsup;
        brsup.id = mynodes[i]->id;
        brsup.name = mynodes[i]->name;
        brsup.geneCF = mysg[i]->getWeight()/decisive_counts[i];
        brsup.geneN = decisive_counts[i];
        branch_supports[brsup.id] = brsup;
        
		stringstream tmp;
		if (mysg[i]->getWeight() == 0.0)
			tmp << "0";
		else
			tmp << round((mysg[i]->getWeight()/decisive_counts[i])*1000)/10;
		if (verbose_mode >= VB_MED)
			tmp << "%" << decisive_counts[i];

        if (Params::getInstance().newick_extended_format) {
            if (mynodes[i]->name.empty() || mynodes[i]->name.back() != ']')
                mynodes[i]->name += "[&CF=" + tmp.str() + "]";
            else
                mynodes[i]->name = mynodes[i]->name.substr(0, mynodes[i]->name.length()-1) + ",!CF=" + tmp.str() + "]";
        } else {
            if (!mynodes[i]->name.empty())
                mynodes[i]->name.append("/");
            mynodes[i]->name.append(tmp.str());
        }
		if (verbose_mode >= VB_MED) {
			cout << mynodes[i]->name << " " << occurence_trees[i] << endl;
		}
	}
	for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
		delete (*it);
}
*/

void MTree::computeRFDist(const char *trees_file, DoubleVector &dist, int assign_sup) {
	cout << "Reading input trees file " << trees_file << endl;
	try {
		ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(trees_file);
        computeRFDist(in, dist, assign_sup);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, trees_file);
	}
}

void MTree::computeRFDist(istream &in, DoubleVector &dist, int assign_sup, bool one_tree) {
	SplitGraph mysg;
    NodeVector nodes;
	convertSplits(mysg, &nodes, root->neighbors[0]->node);
    StringIntMap name_index;
    int ntrees, taxid;
    for (taxid = 0; taxid < mysg.getNTaxa(); taxid++) {
        name_index[mysg.getTaxa()->GetTaxonLabel(taxid)] = taxid;
    }
    NodeVector::iterator nit;
    if (assign_sup) {
        for (nit = nodes.begin(); nit != nodes.end(); nit++) {
            (*nit)->height = 0.0;
        }
    }    
	SplitGraph::iterator sit;
    for (sit = mysg.begin(); sit != mysg.end(); sit++) {
        (*sit)->setWeight(0.0);
    }
	for (ntrees = 1; !in.eof(); ntrees++) {
		MTree tree;
		bool is_rooted = false;

		// read in the tree and convert into split system for indexing
		tree.readTree(in, is_rooted);
		if (verbose_mode >= VB_DEBUG)
			cout << ntrees << " " << endl;
		StrVector taxname;
		tree.getTaxaName(taxname);
		// create the map from taxa between 2 trees
		Split taxa_mask(leafNum);
		for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
            if (name_index.find(*it) == name_index.end())
                outError("Taxon not found in full tree: ", *it);
			taxid = name_index[*it];
			taxa_mask.addTaxon(taxid);
		}
		// make the taxa ordering right before converting to split system
		taxname.clear();
		int smallid;
		for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
			if (taxa_mask.containTaxon(taxid)) {
				taxname.push_back(mysg.getTaxa()->GetTaxonLabel(taxid));
				string name = (string)mysg.getTaxa()->GetTaxonLabel(taxid);
				tree.findLeafName(name)->id = smallid++;
			}
		ASSERT(taxname.size() == tree.leafNum);

		SplitGraph sg;
		//NodeVector nodes;
		tree.convertSplits(sg);
		SplitIntMap hash_ss;
		for (sit = sg.begin(); sit != sg.end(); sit++)
			hash_ss.insertSplit((*sit), 1);

		// now scan through all splits in current tree
		int common_splits = 0;
		for (sit = mysg.begin(); sit != mysg.end(); sit++)
		if ((*sit)->trivial() < 0) // it is an internal split
		{

			Split *subsp = (*sit)->extractSubSplit(taxa_mask);
			if (subsp->shouldInvert())
				subsp->invert();
			Split *sp = hash_ss.findSplit(subsp);
			if (sp) {
				common_splits++;
				//(*sit)->setWeight((*sit)->getWeight()+1.0);
				if (verbose_mode >= VB_MAX) {
					for (taxid = 0; taxid < (*sit)->getNTaxa(); taxid++)
						if ((*sit)->containTaxon(taxid))
							cout << " " << mysg.getTaxa()->GetTaxonLabel(taxid);
					cout << " --> ";
					for (taxid = 0; taxid < sp->getNTaxa(); taxid++)
						if (sp->containTaxon(taxid))
							cout << " " << taxname[taxid];
					cout << endl;
				}
                if (assign_sup && subsp->trivial() < 0)
                    nodes[sit-mysg.begin()]->height++;
			}
			delete subsp;
		}

		//cout << "common_splits = " << common_splits << endl;
        double max_dist = branchNum-leafNum + tree.branchNum-tree.leafNum;
        double rf_val = max_dist - 2*common_splits;
        if (Params::getInstance().normalize_tree_dist) {
            rf_val = rf_val / max_dist;
        }
		dist.push_back(rf_val);
		char ch;
		in.exceptions(ios::goodbit);
		(in) >> ch;
		if (in.eof()) break;
		in.unget();
		in.exceptions(ios::failbit | ios::badbit);
        
        if (one_tree)
            break;

	}
    if (assign_sup)
        for (nit = nodes.begin(); nit != nodes.end(); nit++)
            if (!(*nit)->isLeaf())
                (*nit)->name = convertIntToString((*nit)->height);

//	cout << ntrees << " trees read" << endl;


}

void MTree::reportDisagreedTrees(vector<string> &taxname, MTreeSet &trees, Split &mysplit) {
	for (MTreeSet::iterator it = trees.begin(); it != trees.end(); it++) {
		MTree *tree = (*it);
		SplitGraph sg;
		tree->convertSplits(taxname, sg);
		if (!sg.containSplit(mysplit)) {
			tree->printTree(cout, 0); // don't print branch lengths
			cout << endl;
		}
	}
}


void MTree::createBootstrapSupport(vector<string> &taxname, MTreeSet &trees, SplitIntMap &hash_ss,
    char *tag, Node *node, Node *dad) {
	if (!node) node = root;	
	FOR_NEIGHBOR_IT(node, dad, it) {
		if (!node->isLeaf() && !(*it)->node->isLeaf()) {
			vector<int> taxa;
			getTaxaID(taxa, (*it)->node, node);
			Split mysplit(leafNum, 0.0, taxa);
			if (mysplit.shouldInvert())
				mysplit.invert();
			//mysplit.report(cout);
			//SplitIntMap::iterator ass_it = hash_ss.find(&mysplit);
			Split *sp = hash_ss.findSplit(&mysplit);
			// if found smt
			if (sp != NULL) {
				//Split *sp = ass_it->first;
				/*char tmp[100];
				if ((*it)->node->name.empty()) {
					sprintf(tmp, "%d", round(sp->getWeight()));
				} else
					sprintf(tmp, "/%d", round(sp->getWeight()));*/
				stringstream tmp;
				if ((*it)->node->name.empty())
				  tmp << sp->getWeight();
				else
				  tmp << "/" << sp->getWeight();
                  
                // assign tag
                if (tag && (iEquals(tag, "ALL") || (*it)->node->name == tag))
                    tmp << sp->getName();                
				(*it)->node->name.append(tmp.str());
			} else {
				if (!(*it)->node->name.empty()) (*it)->node->name.append("/");
				(*it)->node->name.append("0");
				if (verbose_mode >= VB_MED) {
					cout << "split not found:" << endl;
					mysplit.report(cout);
				}
			} 
			/* new stuff: report trees that do not contain the split */
			if (strncmp((*it)->node->name.c_str(), "INFO", 4) == 0) {
				cout << "Reporting trees not containing the split " << (*it)->node->name << endl;
				reportDisagreedTrees(taxname, trees, mysplit);
			}
		}
		createBootstrapSupport(taxname, trees, hash_ss, tag, (*it)->node, node);
	}	
}

void MTree::removeNode(Node *dad, Node *node) {
//    Node *child = (*it)->node;
    bool first = true;
    FOR_NEIGHBOR_IT(node, dad, it2) {
        if (first)
            dad->updateNeighbor(node, (*it2)->node, (*it2)->length);
        else
            dad->addNeighbor((*it2)->node, (*it2)->length);
        (*it2)->node->updateNeighbor(node, dad);
        first = false;
    }
    delete node;

//    Node *child = (*it)->node;
//    bool first = true;
//    FOR_NEIGHBOR_IT(child, node, it2) {
//        if (first)
//            node->updateNeighbor(child, (*it2)->node, (*it2)->length);
//        else
//            node->addNeighbor((*it2)->node, (*it2)->length);
//        (*it2)->node->updateNeighbor(child, node);
//        first = false;
//    }
//    delete child;
}


int MTree::collapseZeroBranches(Node *node, Node *dad, double threshold) {
	if (!node) node = root;
    int count = 0;
	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		count += collapseZeroBranches((*it)->node, node, threshold);
	}
	NeighborVec nei_vec;
	nei_vec.insert(nei_vec.begin(), node->neighbors.begin(), node->neighbors.end());
	for (it = nei_vec.begin(); it != nei_vec.end(); it++) 
	if ((*it)->node != dad && (*it)->length <= threshold) {
		// delete the child node
        removeNode(node, (*it)->node);
        count++;
	}
    return count;
}

int MTree::collapseInternalBranches(Node *node, Node *dad, double threshold) {
	if (!node) node = root;
    int count = 0;
	FOR_NEIGHBOR_DECLARE(node, dad, it) {
		count += collapseInternalBranches((*it)->node, node, threshold);
	}
    if (node->isLeaf()) {
        return count;
    }
	NeighborVec nei_vec;
	nei_vec.insert(nei_vec.begin(), node->neighbors.begin(), node->neighbors.end());
	for (it = nei_vec.begin(); it != nei_vec.end(); it++) 
	if ((*it)->node != dad && !(*it)->node->isLeaf() && (*it)->length <= threshold) {
		// delete the child node
        removeNode(node, (*it)->node);
        count++;
	}
    return count;
}

void MTree::insertTaxa(StrVector &new_taxa, StrVector &existing_taxa) {
	if (new_taxa.empty()) return;
	IntVector id;
	int i;
	id.resize(new_taxa.size());
	for (i = 0; i < id.size(); i++)
		id[i] = i;
	// randomize order before reinsert back into tree
	my_random_shuffle(id.begin(), id.end());

	for (int i = 0; i < new_taxa.size(); i++) {
		Node *old_taxon = findLeafName(existing_taxa[id[i]]);
		ASSERT(old_taxon);
		double len = old_taxon->neighbors[0]->length;
		Node *old_node = old_taxon->neighbors[0]->node;
		Node *new_taxon = newNode(leafNum+i, new_taxa[id[i]].c_str());
		Node *new_node = newNode();
		// link new_taxon - new_node
		new_taxon->addNeighbor(new_node, 0.0);
		new_node->addNeighbor(new_taxon, 0.0);
		// link old_taxon - new_node
		new_node->addNeighbor(old_taxon, 0.0);
		old_taxon->updateNeighbor(old_node, new_node, 0.0);
		// link old_node - new_node
		new_node->addNeighbor(old_node, len);
		old_node->updateNeighbor(old_taxon, new_node, len);
	}

    leafNum = leafNum + new_taxa.size();
    initializeTree();
}

Node *MTree::findFirstTaxon(Node *node, Node *dad) {
	if (!node) node = root;
//	Node *next;
	for (int i = 0; i < nodeNum; i++)
		FOR_NEIGHBOR_IT(node, dad, it) {
			if ((*it)->node->isLeaf()) return (*it)->node;
			dad = node;
			node = (*it)->node;
            break;
		}
	return NULL;
}

int MTree::removeTaxa(const StrVector &taxa_names,
                      bool reassignNodeIDs, const char* context) {
    if (taxa_names.empty()) {
        return 0;
    }
    int count = 0;
    bool showingProgress = (context!=nullptr && *context!='\0');
    if (showingProgress) {
        initProgress(taxa_names.size(), context, "removed", "taxon");
    }
    for (auto sit = taxa_names.begin(); sit != taxa_names.end(); sit++) {
        Node *node = findLeafName(*sit);
        if (!node) {
            if (context!=nullptr) {
                trackProgress(1.0);
            }
            continue;
        }
        count++;
        if (node == root) {
            // find another root
            root = findFirstTaxon(root);
        }
		Node *innode = node->neighbors[0]->node;
		Node *othernodes[2] = { NULL, NULL };
		int i;
		double length = 0;

		bool should_merge = true;

		FOR_NEIGHBOR_DECLARE(innode, node, it)	{
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
				for (it = othernodes[i]->neighbors.begin(); it != othernodes[i]->neighbors.end(); it++)
					if ((*it)->node == innode)
					{
						(*it)->node = othernodes[1-i];
						(*it)->length = length;
					}
            delete innode;
		} else {
			// simple delete the neighbor of innode
			for (it = innode->neighbors.begin(); it != innode->neighbors.end(); it++)
				if ((*it)->node == node) {
					innode->neighbors.erase(it);
					break;
				}
		}
		delete node;
        if (showingProgress) {
            trackProgress(1.0);
        }
    }
    if (showingProgress) {
        doneProgress();
    }
    if (count==0 || !reassignNodeIDs) {
        return count;
    }
    NodeVector taxa;
    getTaxa(taxa);
    ASSERT(taxa.size() > 0);
    // reassign taxon IDs
    int id = 0;
    for (NodeVector::iterator nit = taxa.begin(); nit != taxa.end(); nit++) {
        if (*nit == root && rooted) {
            (*nit)->id = taxa.size()-1;
        } else {
            (*nit)->id = id++;
        }
    }
    leafNum = taxa.size();
    initializeTree();
    return count;
}

void MTree::getSplits(SplitGraph &splits, Node* node, Node* dad) {
   if (!node) {
       node = root;
   }
   FOR_NEIGHBOR_IT(node, dad, it) {
           getSplits(splits, (*it)->node, node);
           Split* mySplit = new Split(*((*it)->split));
           if (mySplit->shouldInvert())
               mySplit->invert();
           splits.push_back(mySplit);
       }
}

void MTree::buildNodeSplit(Split *resp, Node *node, Node *dad) {
    if (!node) {
        node = root;
        // The neighbor that represents root
        Neighbor* rootNei = root->neighbors[0]->node->findNeighbor(root);
        if (rootNei->split == NULL) {
            rootNei->split = new Split(leafNum);
        } else {
            delete rootNei->split;
            rootNei->split = new Split(leafNum);
        }
        resp = rootNei->split;
    }
    bool has_child = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
            if ((*it)->split == NULL) {
                (*it)->split = new Split(leafNum);
            } else {
                delete (*it)->split;
                (*it)->split = new Split(leafNum);
            }
            buildNodeSplit((*it)->split, (*it)->node, node);
            //(*it)->split->report(cout);
            *resp += *((*it)->split);
            has_child = true;
        }

    if (dad != NULL) {
        Neighbor* dadNei = node->findNeighbor(dad);
        dadNei->split = new Split(*resp);
        dadNei->split->invert();
    }

    if (!has_child) {
        resp->addTaxon(node->id);
    }
}

void MTree::initializeSplitMap(Split *resp, Node *node, Node *dad) {
    if (!node) node = root;
    if (!resp) {
        resp = new Split(leafNum);
    }
    bool has_child = false;
    FOR_NEIGHBOR_IT(node, dad, it) {
            Split *sp = new Split(leafNum);
            initializeSplitMap(sp, (*it)->node, node);
            *resp += *sp;
            if (sp->shouldInvert())
                sp->invert();
            /* ignore nodes with degree of 2 because such split will be added before */
            if (node->degree() != 2) {
                Branch curBranch((*it)->node, node);
                splitBranchMap.insert(make_pair(sp, curBranch));
            }
            has_child = true;
        }
    if (!has_child) {
        resp->addTaxon(node->id);
    }
}

Node *MTree::findFarthestLeaf(Node *node, Node *dad) {
    if (rooted) // special treatment for rooted tree
        return root;
    if (!node) 
        node = root;
    
    if (dad && node->isLeaf()) {
        node->height = 0.0;
        return node;
    }
    Node *res = NULL;
    node->height = 0.0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *leaf = findFarthestLeaf((*it)->node, node);
        if (node->height < (*it)->node->height+1) {
            node->height = (*it)->node->height+1;
            res = leaf;
        }
    }
    return res;
}

void MTree::getPreOrderBranches(NodeVector &nodes, NodeVector &nodes2, Node *node, Node *dad) {
    if (dad) {
        nodes.push_back(node);
        nodes2.push_back(dad);
    }

    NeighborVec neivec = node->neighbors;
    NeighborVec::iterator i1, i2;
    for (i1 = neivec.begin(); i1 != neivec.end(); i1++)
        for (i2 = i1+1; i2 != neivec.end(); i2++)
            if ((*i1)->node->height > (*i2)->node->height) {
                Neighbor *nei = *i1;
                *i1 = *i2;
                *i2 = nei;
            }
    for (i1 = neivec.begin(); i1 != neivec.end(); i1++)
        if ((*i1)->node != dad)
            getPreOrderBranches(nodes, nodes2, (*i1)->node, node);
//    FOR_NEIGHBOR_IT(node, dad, it) 
//        getPreOrderBranches(nodes, nodes2, (*it)->node, node);
}
