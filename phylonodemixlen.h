//
//  phylonodemixlen.h
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#ifndef __iqtree__phylonodemixlen__
#define __iqtree__phylonodemixlen__

#include <stdio.h>
#include "phylonode.h"

/**
A neighbor in a phylogenetic tree with mixture branch lengths

    @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloNeighborMixlen : public PhyloNeighbor {
public:

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param aid branch ID
     */
    PhyloNeighborMixlen(Node *anode, double alength, int aid = -1) : PhyloNeighbor(anode, alength, aid) {
        lengths.clear();
    }

    PhyloNeighborMixlen(Node *anode, DoubleVector &alength, int aid = -1) : PhyloNeighbor(anode, -1.0, aid) {
        lengths = alength;
        if (!lengths.empty()) {
            length = 0.0;
            for (int i = 0; i < lengths.size(); i++)
                length += lengths[i];
            length /= lengths.size();
        }
    }

    /** branch lengths for mixture */
    DoubleVector lengths;

    /**
        get branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual double getLength(int c) { 
        if (lengths.empty())
            return length; 
        assert(c < lengths.size());
        return lengths[c];
    }

protected:

};

class PhyloNodeMixlen : public PhyloNode {
public:
    /**
        constructor
     */
    PhyloNodeMixlen() : PhyloNode() {}

    /**
        constructor
        @param aid id of this node
     */
    PhyloNodeMixlen(int aid) : PhyloNode(aid) {}

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNodeMixlen(int aid, int aname) : PhyloNode(aid, aname) {}

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNodeMixlen(int aid, const char *aname) : PhyloNode(aid, aname) {}

    /**
        add a neighbor
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);

    /**
        add a neighbor for heterotachy model
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, DoubleVector &length, int id = -1);

protected:
};

#endif /* defined(__iqtree__phylonodemixlen__) */
