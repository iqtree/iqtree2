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

class PhyloNodeMixlen;
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
    PhyloNeighborMixlen(Node *anode, double alength, int aid = -1);

    PhyloNeighborMixlen(Node *anode, DoubleVector &alength, int aid = -1);

    PhyloNeighborMixlen(const PhyloNeighborMixlen &nei);
    
    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual PhyloNeighborMixlen* newNeighbor() const;

    /** branch lengths for mixture */
    DoubleVector lengths;

    /**
        get branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual double getLength(int c);

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void getLength(DoubleVector &vec);

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec (OUT) destination branch length vector
        @param start_pos starting position in vec to copy to
    */
    virtual void getLength(DoubleVector &vec, int start_pos);


    /**
        set branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual void setLength(int c, double len);

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void setLength(const DoubleVector &vec);

    /**
        set branch length by length of a Neighbor, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param nei source neigbor to copy branch lengths
        @return branch length for class c
    */
    virtual void setLength(Neighbor *nei);
    
    /**
        set branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec source branch length vector
        @param start_pos starting position in vec to copy from
    */
    virtual void setLength(const DoubleVector &vec, int start_pos, int num_elem);
    
    virtual PhyloNodeMixlen* getNode() const;
};

class PhyloNodeMixlen : public PhyloNode {
public:
    /**
        constructor
     */
    PhyloNodeMixlen();

    /**
        constructor
        @param aid id of this node
     */
    PhyloNodeMixlen(int aid);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNodeMixlen(int aid, int aname);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNodeMixlen(int aid, const char *aname);

    /**
        add a neighbor
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    void addNeighbor(Node *node, double length, int id = -1);

    /**
        add a neighbor for heterotachy model
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    void addNeighbor(Node *node, DoubleVector &length, int id = -1);

    virtual PhyloNeighborMixlen* firstNeighbor() const;

    virtual PhyloNeighborMixlen* findNeighbor(PhyloNodeMixlen* node) const;

};

typedef SubclassPointerVector<PhyloNodeMixlen, NodeVector> PhyloNodeMixlenVector;


#endif /* defined(__iqtree__phylonodemixlen__) */
