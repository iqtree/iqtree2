//
// C++ Interface: phylonode
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PHYLONODE_H
#define PHYLONODE_H

#include <node.h>

/**
A neighbor in a phylogenetic tree

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class PhyloNeighbor : public Neighbor {

	friend class PhyloNode;
	friend class PhyloTree;
	friend class IQPTree;

public:
	/**
		construct class with a node and length
		@param anode the other end of the branch
		@param alength length of branch
	*/
	PhyloNeighbor(Node *anode, double alength) : Neighbor(anode, alength) {	
		partial_lh = NULL; 
		partial_lh_computed = false;
		lh_scale_factor = 0.0;
	}

	/**
		destructor
	*/
	virtual ~PhyloNeighbor() {
		if (!partial_lh) delete partial_lh;
		partial_lh = NULL;
	}

	/**
		tell that the partial likelihood vector is not computed
	*/
	void clearPartialLh() { partial_lh_computed = false; }

private:

	/**
		true if the partial likelihood was computed
	*/
	bool partial_lh_computed;

	/**
		vector containing the partial likelihoods
	*/
	double *partial_lh;

	/**
		likelihood scaling factor
	*/
	double lh_scale_factor;
};

/**
A node in a phylogenetic tree

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class PhyloNode : public Node
{
	friend class PhyloTree;

public:
	/**
		constructor 
	*/
    PhyloNode();

	/**
		constructor 
		@param aid id of this node
	*/
	PhyloNode(int aid);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	PhyloNode(int aid, int aname);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	PhyloNode(int aid, const char *aname);

	/**
		initialization
	*/
	void init();

	/**
		add a neighbor
		@param node the neighbor node
		@param length branch length
	*/
	virtual void addNeighbor(Node *node, double length);



	/**
		tell that all partial likelihood vectors below this node are not computed
	*/
	void clearAllPartialLh(PhyloNode *dad);

	/**
		tell that all partial likelihood vectors (in reverse direction) below this node are not computed
	*/
	void clearReversePartialLh(PhyloNode *dad);

};


/**
	Node vector
*/
typedef vector<PhyloNode*> PhyloNodeVector;


#endif
