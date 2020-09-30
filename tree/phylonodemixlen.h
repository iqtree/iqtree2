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

    PhyloNeighborMixlen(PhyloNeighborMixlen *nei) : PhyloNeighbor(nei) {
        lengths = nei->lengths;
    }
    
    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual PhyloNeighborMixlen* newNeighbor() {
        return (new PhyloNeighborMixlen(this));
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
        ASSERT(c < lengths.size());
        return lengths[c];
    }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void getLength(DoubleVector &vec) { 
        vec = lengths;
    }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec (OUT) destination branch length vector
        @param start_pos starting position in vec to copy to
    */
    virtual void getLength(DoubleVector &vec, int start_pos) { 
        ASSERT(start_pos+lengths.size() <= vec.size());
        for (int i = 0; i < lengths.size(); i++)
            vec[start_pos+i] = lengths[i];
    }


    /**
        set branch length for a mixture class c, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param c class index
        @return branch length for class c
    */
    virtual void setLength(int c, double len) {
        if (lengths.empty()) {
            length = len;
            return;
        }
        ASSERT(c < lengths.size());
        lengths[c] = len;
    }

    /**
        get branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @return branch length for class c
    */
    virtual void setLength(DoubleVector &vec) { 
        lengths = vec;
    }

    /**
        set branch length by length of a Neighbor, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param nei source neigbor to copy branch lengths
        @return branch length for class c
    */
    virtual void setLength(Neighbor *nei) { 
        length = nei->length; 
        lengths = ((PhyloNeighborMixlen*)nei)->lengths;
    }
    
    /**
        set branch lengths, used by heterotachy model (PhyloNeighborMixlen)
        the default is just to return a single branch length
        @param vec source branch length vector
        @param start_pos starting position in vec to copy from
    */
    virtual void setLength(DoubleVector &vec, int start_pos, int num_elem) { 
        ASSERT(start_pos+num_elem <= vec.size());
        lengths.clear();
        lengths.insert(lengths.begin(), vec.begin()+start_pos, vec.begin()+start_pos+num_elem);
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
