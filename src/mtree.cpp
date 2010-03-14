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
#include "splitgraph.h"

/*********************************************
	class MTree
*********************************************/

MTree::MTree() { 
	root = NULL; 
	leafNum = 0;
	nodeNum = 0;
	rooted = false;
}

MTree::MTree(const char *userTreeFile, bool &is_rooted)
{
	init(userTreeFile, is_rooted);
}

void MTree::init(const char *userTreeFile, bool &is_rooted) {
	readTree(userTreeFile, is_rooted);
	//printInfo();
}


/**
	constructor, get from another tree
*/
MTree::MTree(MTree &tree) {
	init(tree);
}

void MTree::init(MTree &tree) {
	root = tree.root;
	leafNum = tree.leafNum;
	nodeNum = tree.nodeNum;
	rooted = tree.rooted;
	//userFile = tree.userFile;
	// have to delete the root when exchange to another object
	tree.root = NULL;
}

void MTree::copyTree(MTree *tree) {
	if (root) freeNode();
	stringstream ss;
	tree->printTree(ss);
	readTree(ss, tree->rooted);
}


Node* MTree::newNode(int node_id, const char* node_name) {
	return new Node(node_id, node_name);
}

Node* MTree::newNode(int node_id, int node_name) {
	return new Node(node_id, node_name);
}


void MTree::printInfo(Node *node, Node *dad)
{
	if (node == NULL) node = root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		cout << node->name << " " << (*it)->node->name << " " << (*it)->length << endl;
		printInfo((*it)->node, node);
	}
}

void MTree::printTree(const char *ofile, int brtype)
{
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(ofile);
		printTree(out, brtype);
		out.close();
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
}

void MTree::printTree(ostream &out, int brtype, Node *node, Node *dad)
{
	out.precision(6);
	if (!node) node = root;
	if (node->isLeaf()) {
		if (brtype & WT_TAXON_ID) 
			out << node->id;
		else
			out << node->name;

		if (brtype & WT_BR_LEN)
			out << ":" << node->neighbors[0]->length; 
	} else {
		// internal node
		out << "(";
		bool first = true;
		double length = 0.0;
		//for (int i = 0; i < node->neighbors.size(); i++)
			//if (node->neighbors[i]->node != dad)
		FOR_NEIGHBOR_IT(node, dad, it) {
			if ((*it)->node->name != ROOT_NAME) {
				if (!first)
					out << ",";
				printTree(out, brtype, (*it)->node, node);
				first = false;
			} else
				length = (*it)->length;
		} else {
			length = (*it)->length;
		}
		out << ")";
		if (!node->name.empty())
			out << node->name;
		if (dad != NULL || length > 0.0)
			if (brtype & WT_BR_LEN)
				out << ":" << length; 
			else if (brtype & WT_BR_CLADE) {
				if (! node->name.empty()) out << "/";
				out << length; 
			}
	}
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
			if (subtree[(*it)->node->id] != NULL)
				if((*it)->node->name != ROOT_NAME) {
					if (!first)
						out << ",";
					printSubTree(out, subtree, (*it)->node, node);
					first = false;
				} else
					length += (*it)->length;
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
	cout << "Reading tree file " << infile << " ..." << endl;
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

	cout << "Tree contains " << leafNum - is_rooted << 
		" taxa and " << nodeNum-1-is_rooted << " branches" << endl;
}


void MTree::readTree(istream &in, bool &is_rooted)
{
	in_line = 1;
	in_column = 1;
	try {
		char ch;
		ch = readNextChar(in);
		if (ch != '(')
			throw "Tree file not started with an opening-bracket '('";
	
		leafNum = 0;
	
		double branch_len;
		Node *node;
		parseFile(in, ch, node, branch_len);
		if (is_rooted || branch_len > 0.0) {
			if (branch_len == -1.0) branch_len = 0.0;
			if (branch_len < 0.0) 
				throw ERR_NEG_BRANCH; 
			is_rooted = true;
			root = newNode(leafNum, ROOT_NAME);
			root->addNeighbor(node, branch_len);
			node->addNeighbor(root, branch_len);
			leafNum++;
		}
		// make sure that root is a leaf
		assert(root->isLeaf());
	
		if (in.eof() || ch != ';')
			throw "Tree file must be ended with a semi-colon ';'";
	} catch (bad_alloc) {
		outError(ERR_NO_MEMORY);
	} catch (const char *str) {
		outError(str, reportInputInfo());
	} catch (char *str) {
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
	if (!node) node = root;
	if (!node->isLeaf())
	{
		node->id = nodeNum;
		nodeNum++;
		//node->name = node->id;

	}
	//for (int i = 0; i < node->neighbors.size(); i++)
		//if (node->neighbors[i]->node != dad)
	FOR_NEIGHBOR_IT(node, dad, it)
		initializeTree((*it)->node, node);
}


void MTree::parseFile(istream &infile, char &ch, Node* &root, double &branch_len)
{
	Node *node;
	int maxlen = 100;
	char seqname[100];
	int seqlen;
	double brlen;
	branch_len = -1.0;

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
	while (!is_newick_token(ch) && !controlchar(ch) && !infile.eof() && seqlen < maxlen)
	{
		seqname[seqlen] = ch;
		seqlen++;
		ch = infile.get();
		in_column++;
	}
	if ((controlchar(ch) || ch == '[') && !infile.eof()) 
		ch = readNextChar(infile, ch);
	if (seqlen == maxlen)
		throw "Too long name ( > 100)";
	seqname[seqlen] = 0;
	if (seqlen == 0 && root->isLeaf())
		throw "A taxon has no name.";
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
		ch = readNextChar(infile);
		seqlen = 0;
		while (!is_newick_token(ch) && !controlchar(ch) && !infile.eof() && seqlen < maxlen)
		{
			seqname[seqlen] = ch;
			seqlen++;
			ch = infile.get();
			in_column++;
		}
		if ((controlchar(ch) || ch == '[') && !infile.eof()) 
			ch = readNextChar(infile, ch);
		if (seqlen == maxlen || infile.eof())
			throw "branch length format error.";
		seqname[seqlen] = 0;
		branch_len = convert_double(seqname);
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

void MTree::getInternalNodes(NodeVector &nodes, Node *node, Node *dad) {
	if (!node) node = root;
	//for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++) 
		//if ((*it)->node != dad)	{
	FOR_NEIGHBOR_IT(node, dad, it)
	if (!(*it)->node->isLeaf()) {
		nodes.push_back((*it)->node);
		getInternalNodes(nodes, (*it)->node, node);
	}
}


void MTree::getTaxaName(vector<NxsString> &taxname, Node *node, Node *dad) {
	if (!node) node = root;
	if (node->isLeaf()) {
		taxname[node->id] = node->name;
	}
	//for (NeighborVec::iterator it = node->neighbors.begin(); it != node->neighbors.end(); it++) 
		//if ((*it)->node != dad)	{
	FOR_NEIGHBOR_IT(node, dad, it) {
		getTaxaName(taxname, (*it)->node, node);
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


void MTree::convertSplits(SplitGraph &sg, Split *resp, Node *node, Node *dad) {
	if (!node) node = root;
	assert(resp->getNTaxa() == leafNum);
	bool has_child = false;
	FOR_NEIGHBOR_IT(node, dad, it) {
		//vector<int> taxa;
		//getTaxaID((*it)->node, node, taxa);
		
		Split *sp = new Split(leafNum, (*it)->length);
		convertSplits(sg, sp, (*it)->node, node);
		*resp += *sp;
		if (sp->shouldInvert())
			sp->invert();
		sg.push_back(sp);
		has_child = true;
	}
	if (!has_child)
		resp->addTaxon(node->id);
}

void MTree::convertSplits(vector<NxsString> &taxname, SplitGraph &sg) {
	if (!sg.taxa)
		sg.taxa = new NxsTaxaBlock();
	if (!sg.splits)
		sg.splits = new MSplitsBlock(&sg);
	if (!sg.pda)
		sg.pda = new MPdaBlock(&sg);

	for (vector<NxsString>::iterator it = taxname.begin(); it != taxname.end(); it++)
		sg.taxa->AddTaxonLabel(*it);
	
	// make the cycle
	getTaxaID(sg.splits->cycle);
	// make the splits
	Split sp(leafNum);
	convertSplits(sg, &sp);
}

void MTree::convertSplits(SplitGraph &sg) {

	// make the taxa name
	vector<NxsString> taxname;
	taxname.resize(leafNum);
	getTaxaName(taxname);

	convertSplits(taxname, sg);
}

inline int splitnumtaxacmp(const Split* a, const Split* b)
{
	return (a->countTaxa() < b->countTaxa());
}

void MTree::convertToTree(SplitGraph &sg) {
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
	SplitGraph::iterator it;
	int taxid;
	int count = 0;
	// first add all trivial splits into tree
	for (it = sg.begin(); it != sg.end(); it++, count++) {
		//(*it)->report(cout);
		taxid = (*it)->trivial();
		if (taxid < 0) break;
		assert(leaves[taxid] == NULL);
		leaves[taxid] = newNode(taxid, sg.getTaxa()->GetTaxonLabel(taxid).c_str());
		leaves[taxid]->addNeighbor(NULL, (*it)->getWeight());
		cladetaxa[taxid] = (*it);
	}
	// now fill in all missing taxa with zero terminal branch
	for (taxid = 0; taxid < leafNum; taxid++) 
		assert(leaves[taxid]);

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
		assert(count == mysp->countTaxa());
		cladetaxa.push_back(mysp);
		leaves.push_back(newnode);

		newnode->addNeighbor(NULL, mysp->getWeight());
		nodeNum++;
	}
	assert(leaves.size() >= 3);
	Node *newnode = newNode(nodeNum);
	for (taxid = 0; taxid < leaves.size(); taxid++) {
		double len = leaves[taxid]->updateNeighbor(NULL, newnode);
		newnode->addNeighbor(leaves[taxid], len);
	}
	root = newnode;
	nodeNum++;
	cladetaxa.clear();
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

void MTree::freeNode(Node *node, Node *dad)
{
	if (!node) node = root;
	NeighborVec::reverse_iterator it;
	for (it = node->neighbors.rbegin(); it != node->neighbors.rend(); it++)
		if ((*it)->node != dad) {
			freeNode((*it)->node, node);
		}
	delete node;
}

char MTree::readNextChar(istream &in, char current_ch) {
	char ch;
	if (current_ch == '[') 
		ch = current_ch; 
	else { 
		in.get(ch); 
		in_column++;
		if (ch == 10) { in_line++; in_column = 1; }
	}
	while (controlchar(ch) && !in.eof()) { 
		in.get(ch); 
		in_column++;
		if (ch == 10) { in_line++; in_column = 1; }
	}
	// ignore comment
	while (ch=='[' && !in.eof()) {
		while (ch!=']' && !in.eof()) {
			in.get(ch); 
			in_column++;
			if (ch == 10) { in_line++; in_column = 1; }
		}
		if (ch != ']') throw "Comments not ended with ]";
		in_column++;
		in.get(ch); 
		if (ch == 10) { in_line++; in_column = 1; }
		while (controlchar(ch) && !in.eof()) {
			in_column++;
			in.get(ch); 
			if (ch == 10) { in_line++; in_column = 1; }
		}
	}
	return ch;
}

string MTree::reportInputInfo() {
	string str = " (line ";
	str += convertIntToString(in_line) + " column " + convertIntToString(in_column-1) + ")";
	return str;
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
