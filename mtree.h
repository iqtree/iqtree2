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
#ifndef MTREE_H
#define MTREE_H

#include "node.h"
//#include "splitgraph.h"
#include "split.h"
#include <iostream>
#include <sstream>
#include "hashsplitset.h"
#include "splitset.h"

const char ROOT_NAME[] = "_root";

class SplitGraph;

/**
General-purposed tree
@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
 */
class MTree {
public:

    /********************************************************
            CONSTRUCTORs, INITIALIZATION AND DESTRUCTORs
     ********************************************************/

    /**
            constructor, read tree from user file
            @param userTreeFile the name of the user tree
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    MTree(const char *userTreeFile, bool &is_rooted);

    /**
            constructor, get from another tree
            @param tree another MTree
     */
    MTree(MTree &tree);

    /**
            constructor
     */
    MTree();

    /**
            copy the tree structure into this tree
            @param tree the tree to copy
     */
    virtual void copyTree(MTree *tree);

    /**
            copy the sub-tree structure into this tree
            @param tree the tree to copy
            @param taxa_set 0-1 string of length leafNum (1 to keep the leaf)
     */
    virtual void copyTree(MTree *tree, string &taxa_set);

    Node* copyTree(MTree *tree, string &taxa_set, double &len, Node *node = NULL, Node *dad = NULL);

    /**
            initialize the tree from a NEWICK tree file
            @param userTreeFile the name of the user tree
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    void init(const char *userTreeFile, bool &is_rooted);

    /**
            initialize tree, get from another tree
            @param tree another MTree
     */
    void init(MTree &tree);


    /**
            destructor
     */
    virtual ~MTree();

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name
            @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = NULL);

    virtual Node* newNode(int node_id, int node_name);


    /**
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return the number of branches with zero length ( <= epsilon)
     */
    int countZeroBranches(Node *node = NULL, Node *dad = NULL, double epsilon = 0.000001);

    /**
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return the number of internal branches with zero length ( <= epsilon)
     */
    int countZeroInternalBranches(Node *node = NULL, Node *dad = NULL, double epsilon = 0.000001);

	/**
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@return the number of long branches
	*/
	int countLongBranches(Node *node = NULL, Node *dad = NULL, double epsilon = 8.8);
    /********************************************************
            PRINT INFORMATION
     ********************************************************/

	/** @return true if tree is bifurcating, false otherwise */
	bool isBifurcating(Node *node = NULL, Node *dad = NULL);
    /**
            print information
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void printBranchLengths(ostream &out, Node *node = NULL, Node *dad = NULL);

    /**
            print the tree to the output file in newick format
            @param outfile the output file.
            @param brtype type of branch to print
     */
    void printTree(const char *outfile, int brtype = WT_BR_LEN);

    /**
            print the tree to the output file in newick format
            @param out the output stream.
            @param brtype type of branch to print
     */
    void printTree(ostream & out, int brtype = WT_BR_LEN);


    string getTreeString();

    /**
            print the tree to the output file in newick format
            @param out the output file.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param brtype type of branch to print
            @return ID of the taxon with smallest ID
     */
    int printTree(ostream &out, int brtype, Node *node, Node *dad = NULL);


    /**
            print the sub-tree to the output file in newick format
            @param out the output file.
            @param subtree list of nodes (internal & external) contained in the new tree
     */
    void printSubTree(ostream &out, NodeVector &subtree);

    /**
            print the sub-tree to the output file in newick format
            @param out the output file.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param subtree list of nodes (internal & external) contained in the new tree
     */
    void printSubTree(ostream &out, NodeVector &subtree, Node *node, Node *dad = NULL);


    /**
            print the taxa set to the output file
            @param outfile the output file.
     */
    void printTaxa(const char *outfile);

    /**
            print the taxa set to the output file
            @param out the output file.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void printTaxa(ostream &out, Node *node = NULL, Node *dad = NULL);

    /**
            print the taxa set of a given subtree
            @param out the output file.
            @param subtree the subtree vector
     */
    void printTaxa(ostream &out, NodeVector &subtree);

    void writeInternalNodeNames(string &out_file);

    /********************************************************
            DRAW TREE
     ********************************************************/

    /**
            Sort the taxa by their IDs
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return smallest taxon ID of the subtree
     */
    int sortTaxa(Node *node = NULL, Node *dad = NULL);

	void drawTree(ostream &out, int brtype = WT_BR_SCALE + WT_INT_NODE, double zero_epsilon = 2e-6);

	/** OBSOLETE:
	void drawTree(ostream &out, int brtype, double brscale, IntVector &sub_tree_br, double zero_epsilon,
            Node *node = NULL, Node *dad = NULL);
    */

	void drawTree2(ostream &out, int brtype, double brscale, IntVector &sub_tree_br, double zero_epsilon,
            Node *node = NULL, Node *dad = NULL);

    /**
     * @param tree the other tree to compare with
     * @return TRUE if this tree is topologically equal to tree
     */
    bool equalTopology(MTree *tree);

    /********************************************************
            READ TREE FROM FILE
     ********************************************************/

    /**
            read the tree from the input file in newick format
            @param infile the input file file.
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    void readTree(const char *infile, bool &is_rooted);

    /**
            read the tree from the ifstream in newick format
            @param in the input stream.
            @param is_rooted (IN/OUT) true if tree is rooted
     */
    void readTree(istream &in, bool &is_rooted);

    /**
            parse the tree from the input file in newick format
            @param infile the input file
            @param ch (IN/OUT) current char
            @param root (IN/OUT) the root of the (sub)tree
            @param branch_len (OUT) branch length associated to the current root
		
     */
    void parseFile(istream &infile, char &ch, Node* &root, double &branch_len);


    /**
            initialize tree, set node structure
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void initializeTree(Node *node = NULL, Node* dad = NULL);


    /********************************************************
            GET INFORMATION
     ********************************************************/

    /**
            @return sum of all branch lengths
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    double treeLength(Node *node = NULL, Node *dad = NULL);

    /**
            @return sum length of all internal branches
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    double treeLengthInternal(double epsilon, Node *node = NULL, Node *dad = NULL);

    /**
            @return maximum path length from root node to taxa
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    double treeDepth(Node *node = NULL, Node *dad = NULL);
    /**
        get the descending taxa ID list below the node
        @param node the starting node, NULL to start from the root
        @param dad dad of the node, used to direct the search
        @param taxa (OUT) taxa ID
     */
    void getTaxaID(vector<int> &taxa, Node *node = NULL, Node *dad = NULL);

    /**
     * get all node within a subtree
     * TODO: This is probably identical with getTaxa
     * @param node root of the subtree
     * @param dad node to define the subtree
     * @param nodeList (OUT) vector containing all nodes of the subtree
     */
    void getAllNodesInSubtree(Node *node, Node *dad, NodeVector &nodeList);

    /**
     * get number of taxa below the node
     * @param node the starting node, NULL to start from the root
     * @param dad dad of the node, used to direct the search
     * @return number of taxa
     */
    int getNumTaxa(Node *node = NULL, Node *dad = NULL);

    /**
            get the descending taxa below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param taxa (OUT) vector of taxa
     */
    void getTaxa(NodeVector &taxa, Node *node = NULL, Node *dad = NULL);

    /**
     	get all descending taxa which are in non-cherry position
  		@param node the starting node, NULL to start from the root
        @param dad dad of the node, used to direct the search
        @param noncherry (OUT) vector of non-cherry taxa
        @param cherry (OUT) vector of cherry taxa
     */
    void getNonCherryLeaves(NodeVector &noncherry, NodeVector &cherry, Node *node = NULL, Node *dad = NULL);

	/**
		get the descending taxa below the node
		@param node the starting node, NULL to start from the root
		@param dad dad of the node, used to direct the search
		@param taxa (OUT) vector of taxa
	*/
	void getTaxa(Split &taxa, Node *node = NULL, Node *dad = NULL);

    /**
            get the descending taxa below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param taxa (OUT) vector of taxa
     */
    void getOrderedTaxa(NodeVector &taxa, Node *node = NULL, Node *dad = NULL);

    /**
            get the descending taxa names below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param taxname (OUT) taxa name
     */
    void getTaxaName(vector<string> &taxname, Node *node = NULL, Node *dad = NULL);

    /**
            get the descending internal nodes below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param nodes (OUT) vector of internal nodes
     */
    void getInternalNodes(NodeVector &nodes, Node *node = NULL, Node *dad = NULL);

    /**
            get the descending internal branches below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param nodes (OUT) vector of one end node of branch
            @param nodes2 (OUT) vector of the other end node of branch
     */
    void getInternalBranches(NodeVector &nodes, NodeVector &nodes2, Node *node = NULL, Node *dad = NULL);

    /**
            get all descending branches below the node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param nodes (OUT) vector of one end node of branch
            @param nodes2 (OUT) vector of the other end node of branch
     */
    void getBranches(NodeVector &nodes, NodeVector &nodes2, Node *node = NULL, Node *dad = NULL);

    /**
     *      get all descending internal branches below \a node and \a dad up to depth \a depth
     *      @param[out] brans the branches in question
     */
    void getInBranches(map<string, Branch> &brans, int depth, Node *node, Node *dad);

    /**
     * @brief: check if the branch is internal
     * @param[in] node1 one end of the branch
     * @param[in] node2 the other end of the branch
     */
    bool isInBran(Node* node1, Node* node2);

    void getBranchLengths(DoubleVector &len, Node *node = NULL, Node *dad = NULL);

    void setBranchLengths(DoubleVector &len, Node *node = NULL, Node *dad = NULL);

    /**
            find a node with corresponding name
            @param name node name
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return node if found, otherwise NULL
     */
    Node *findNodeName(string &name, Node *node = NULL, Node* dad = NULL);

    /**
            find a leaf with corresponding name
            @param name leaf name
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return node if found, otherwise NULL
     */
    Node *findLeafName(string &name, Node *node = NULL, Node* dad = NULL);

    /**
            find a node with corresponding ID
            @param id node ID
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @return node if found, otherwise NULL
     */
    Node *findNodeID(int id, Node *node = NULL, Node* dad = NULL);


    /**
            scale the length of all branches to a norm factor
            @param norm normalized factor
            @param make_int TRUE to round lengths to int, FALSE otherwise
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void scaleLength(double norm, bool make_int = false, Node *node = NULL, Node *dad = NULL);

    /**
            scale the length of all branches for RAxML internal presentation
            @param norm normalized factor
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void transformBranchLenRAX(double factor, Node *node = NULL, Node *dad = NULL);

    /**
            scale the clade supports of all internal nodes to a norm factor
            @param norm normalized factor
            @param make_int TRUE to round lengths to int, FALSE otherwise
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void scaleCladeSupport(double norm, bool make_int = false, Node *node = NULL, Node *dad = NULL);

    /**
            assign the leaf IDs with their names
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search

     */
    void assignLeafID(Node *node = NULL, Node *dad = NULL);


    /********************************************************
            CONVERT TREE INTO SPLIT SYSTEM
     ********************************************************/

    /**
            convert the tree into the split system
            @param sg (OUT) resulting split graph
     */
	void convertSplits(SplitGraph &sg, NodeVector *nodes = NULL, Node *node = NULL, Node *dad = NULL);

    /**
            convert the tree into the split system
            @param taxname certain taxa name
            @param sg (OUT) resulting split graph
     */
	void convertSplits(vector<string> &taxname, SplitGraph &sg, NodeVector *nodes = NULL, Node *node = NULL, Node *dad = NULL);

    /**
            convert the tree into the split system, iterative procedure
            @param sg (OUT) resulting split graph
            @param resp (internal) set of taxa below node
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    void convertSplits(SplitGraph &sg, Split *resp, NodeVector *nodes = NULL, Node *node = NULL, Node *dad = NULL);

    /********************************************************
            CONVERT SPLIT SYSTEM INTO TREE
     ********************************************************/
    /**
            convert compatible split set into tree
            @param sg source split graph
     */
    void convertToTree(SplitGraph &sg);


    /********************************************************
            calculate distance matrix
     ********************************************************/


    /**
            calculate the pairwise distances on the tree, print the matrix to file (in phylip format)
            @param filename file name
     */
    void calcDist(char *filename);

    /**
            calculate the pairwise distances on the tree
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param dist (OUT) distance matrix
     */
    void calcDist(double* &dist, Node *node = NULL, Node *dad = NULL);

    /**
            calculate the pairwise distances on the tree
            @param aroot the starting root
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param cur_len current length from aroot to node
            @param dist (OUT) distance matrix
     */
    void calcDist(Node *aroot, double cur_len, double* &dist, Node *node, Node *dad);
/********************************************************
	STATISTICS
********************************************************/

	void extractQuadSubtrees(vector<Split*> &subtrees, Node *node = NULL, Node *dad = NULL);

	/**
	 * for each branch, assign how many times this branch appears in the input set of trees.
	 * Work fine also when the trees do not have the same taxon set.
	 * @param trees_file set of trees in NEWICK
	 */
	void assignBranchSupport(const char *trees_file);

	void assignBranchSupport(istream &in);

	/**
	 * compute robinson foulds distance between this tree and a set of trees.
	 * Work fine also when the trees do not have the same taxon set.
	 * @param trees_file set of trees in NEWICK
	 * @param dist (OUT) distance vector
	 */
	void computeRFDist(const char *trees_file, IntVector &dist);

	void computeRFDist(istream &in, IntVector &dist);

	/**
	 * insert new taxa next to the existing taxa in the tree
	 * @param new_taxa name of new taxa to be inserted
	 * @param existing_taxa names of existing taxa in the tree
	 */
	void insertTaxa(StrVector &new_taxa, StrVector &existing_taxa);


    /********************************************************
            PROPERTIES OF TREE
     ********************************************************/
    /**
            root node.
     */
    Node *root;

    /**
            number of leaves
     */
    int leafNum;

    /**
            total number of nodes in the tree
     */
    int nodeNum;

    /**
            total number of branches in the tree
     */
    int branchNum;

    /**
            user tree file name
     */
    //char *userFile;

    /**
            TRUE if the tree is rooted
     */
    bool rooted;

    /**
            precision to print branch lengths, default: 6
     */
    int num_precision;

    /** if WT_BR_SCALE turned on, printTree will scale branch length with this factor */
    double len_scale;

    /**
            release the nemory.
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
     */
    int freeNode(Node *node = NULL, Node *dad = NULL);

    void setExtendedFigChar();

protected:

    /**
            line number of the input file, used to output errors in input file
     */
    int in_line;
    /**
            column number of the input file, used to output errors in input file
     */
    int in_column;


    /**
     * special character for drawing tree figure
     * 0: vertical line
     *  1: horizontal line
     *  2: top corner
     *  3: middle corner
     *  4: bottom corner
     */
    string fig_char;

    /**
            check tree is bifurcating tree (every leaf with level 1 or 3)
            @param node the starting node, NULL to start from the root
            @param dad dad of the node, used to direct the search
            @param stop (IN/OUT) set = true to stop the search
     */
    void checkValidTree(bool& stop, Node *node = NULL, Node *dad = NULL);

    /**
            read the next character from a NEWICK file. Ignore comments [...]
            @param in input stream
            @param current_ch current character in the stream
            @return next character read from input stream
     */
    char readNextChar(istream &in, char current_ch = 0);

    string reportInputInfo();

    /**
     * Convert node IDs of a pair of nodes to a string in form "id1-id2"
     * where id1 is smaller than id2. This is done to create a key for the map data structure
     * @param node1
     * @param node2
     * @return the string key for the node pair
     */
    inline string nodePair2String(Node* node1, Node* node2) {
        string key("");
        if (node1->id < node2->id) {
            key += convertIntToString(node1->id) + "-"
                    + convertIntToString(node2->id);
        } else {
            key += convertIntToString(node2->id) + "-"
                    + convertIntToString(node1->id);
        }
        return key;
    }
};

/**
PD set
@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
 */
class PDTaxaSet : public vector<Node*> {
public:
    /**
            PD score
     */
    double score;

    /**
            string representing subtree connecting taxa in the PD set
     */
    string tree_str;

    /**
            name of this taxa set
     */
    string name;

    /**
            assign subtree string
            @param tree a MTree class
            @param subtree list of nodes (internal & external) contained in the tree
            @return score and tree_str variables of this class
     */
    void setSubTree(MTree &tree, NodeVector &subtree);


    /**
            assign the taxa, score and subtree string
            @param tree a MTree class
     */
    void setTree(MTree &tree);

    /**
            print taxa to stream
            @param out output stream
     */
    void printTaxa(ostream &out);

    /**
            print taxa to file
            @param filename output file name
     */
    void printTaxa(char *filename);

    /**
            print tree to stream
            @param out output stream
     */
    void printTree(ostream &out);

    /**
            print tree to file
            @param filename output file name
     */
    void printTree(char *filename);

    /**
            convert from the taxa node vector to set of their IDs
            @param ntaxa total number of taxa
            @param id_set (OUT) set of their IDs
     */
    void makeIDSet(int ntaxa, Split &id_set);

};


#endif
