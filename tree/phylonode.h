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

#include "node.h"

#define FOR_EACH_ADJACENT_PHYLO_NODE(mynode, mydad, it, mychild) \
for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
    for (PhyloNode* mychild = (PhyloNode*)(*it)->node ; mychild != mydad; mychild = mydad )

#define FOR_EACH_PHYLO_NEIGHBOR(mynode, mydad, it, nei) \
for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
    for (PhyloNeighbor* nei = (PhyloNeighbor*)(*it); nei!=nullptr && nei->getNode() != mydad; nei=nullptr )


std::string pointer_to_hex(void *ptr);

typedef unsigned short UBYTE;

/**
 * direction of a Neighbor from the root, for rooted tree only
 */
enum RootDirection {UNDEFINED_DIRECTION, TOWARD_ROOT, AWAYFROM_ROOT};

class PhyloNode;
class BlockAllocator;
/**
A neighbor in a phylogenetic tree

    @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloNeighbor : public Neighbor {
    friend class PhyloNode;
    friend class PhyloTree;
    friend class IQTree;
    friend class PhyloSuperTree;
    friend class PhyloTreeMixlen;
    friend class MemSlotVector;
    friend class ParsTree;
    friend class BlockAllocator;
    friend class TaxonToPlace;
    friend class TargetBranch;

public:
    friend class TinaTree;
    friend class PhyloSuperTreePlen;

    /**
        construct class with a node and length		
        @param anode the other end of the branch
        @param alength length of branch
     */
    PhyloNeighbor(Node *anode, double alength) : Neighbor(anode, alength) {
        partial_lh = NULL;
        scale_num = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
        direction = UNDEFINED_DIRECTION;
        size = 0;
    }

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param aid branch ID
     */
    PhyloNeighbor(Node *anode, double alength, int aid) : Neighbor(anode, alength, aid) {
        partial_lh = NULL;
        scale_num = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
        direction = UNDEFINED_DIRECTION;
        size = 0;
    }

    /**
     construct class with another Neighbor
     @param nei another Neighbor
     */
    PhyloNeighbor(PhyloNeighbor *nei) : Neighbor(nei) {
        partial_lh = NULL;
        scale_num = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
        direction = nei->direction;
        size = nei->size;
    }

    
    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual Neighbor* newNeighbor() {
        return (new PhyloNeighbor(this));
    }

    /**
        tell that the partial likelihood vector is not computed
     */
    inline void clearPartialLh() {
        partial_lh_computed = 0;
    }

    /**
     *  tell that the partial likelihood vector is computed
     */
    inline void unclearPartialLh() {
        partial_lh_computed = 1;
    }

    /**
        clear all partial likelihood recursively in forward direction
        @param dad dad of this neighbor
     */
    void clearForwardPartialLh(PhyloNode *dad);

    /**
        DEPRECATED, moved to PhyloTree
        if partial_lh is NULL, reorient partial_lh (LM_PER_NODE technique)
        @param dad dad of this neighbor
    */
//    void reorientPartialLh(PhyloNode *dad);

	/**
	* For Upper Bounds analysis: get partial likelihood and lh scale factor
	*/
	double* get_partial_lh(){
	return partial_lh;
	}

	double get_lh_scale_factor(){
	return lh_scale_factor;
	}

	int get_partial_lh_computed(){
	return partial_lh_computed;
	}

	/**
	 * true if this Neighbor is directed towards the root
	 */
	bool isTowardsRoot() {
		ASSERT(direction != UNDEFINED_DIRECTION);
		return (direction == TOWARD_ROOT);
	}

    int getSize() {
        return size;
    }
    
    PhyloNode* getNode() {
        return (PhyloNode*)node;
    }
    
    bool isLikelihoodComputed() const {
        return ( partial_lh_computed & LIKELIHOOD_IS_COMPUTED ) != 0;
    }
    
    void setLikelihoodComputed(bool set) {
        if (set) {
            partial_lh_computed |= LIKELIHOOD_IS_COMPUTED;
        } else {
            partial_lh_computed &= ~LIKELIHOOD_IS_COMPUTED;
        }
    }
    
    bool isParsimonyComputed() const {
       return ( partial_lh_computed & PARSIMONY_IS_COMPUTED ) != 0;
    }
           
    void setParsimonyComputed(bool set) {
        if (set) {
            partial_lh_computed |= PARSIMONY_IS_COMPUTED;
        } else {
            partial_lh_computed &= ~PARSIMONY_IS_COMPUTED;
        }
    }
    
    void clearComputedFlags() {
        partial_lh_computed = 0;
    }


private:

    /**
        indicates whether the partial likelihood (and/or) the parsimony,
        is computed (and up to date).
     */
    int partial_lh_computed;
    
    static const int LIKELIHOOD_IS_COMPUTED = 1;
    static const int PARSIMONY_IS_COMPUTED  = 2;

private:
    /**
        vector containing the partial likelihoods
     */
    double *partial_lh;

    /**
        likelihood scaling factor
     */
    double lh_scale_factor;

    /**
        vector containing number of scaling events per pattern // NEW!
     */
    UBYTE *scale_num;

    /**
        vector containing the partial parsimony scores
     */
    UINT *partial_pars;

    /**
     * direction of the Neighbor in a rooted tree
     */
    RootDirection direction;

    /** size of subtree below this neighbor in terms of number of taxa */
    int size;

};

/**
A node in a phylogenetic tree

    @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloNode : public Node {
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
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);

    /**
        tell that all partial likelihood vectors below this node are not computed
     */
    void clearAllPartialLh(bool set_to_null, PhyloNode *dad);
    
    /**
        forget all scale_num vectors below this node
     */
    void clearAllScaleNum(bool set_to_null, PhyloNode* dad);

    /**
        tell that all partial parsimony vectors below this node are not computed
        @param dad the node below which, all partial parsimony vectors are to be marked as uncomputed
     */
    void clearAllPartialParsimony(bool set_to_null, PhyloNode *dad);
    
    /**
        tell that all partial parsimony vectors (in reverse direction), reach from this node
        are not computed (e.g. dad might be a newly inserted node).
     */
    void clearReversePartialParsimony(PhyloNode* dad);
    
    /**
        tell that all partial likelihood vectors (in reverse direction) below this node are not computed
     */
    void clearReversePartialLh(PhyloNode *dad);

    void computeReversePartialLh(PhyloNode *dad);

    /** 
        compute the size (#taxa) of the subtree rooted at this node
        using buffered 'size' attribute if computed beforehand
        @param dad dad of this node
    */
    int computeSize(PhyloNode *dad);

    PhyloNeighbor* findNeighbor(Node* node);
    
    PhyloNeighbor* firstNeighbor();

    PhyloNeighbor* getNeighborByIndex(size_t index);
    
};

template <class T, class S> class SubclassPointerVector: public S {
    //
    //A subclass of vector class S, of pointers,
    //where the entries in S are pointer types of T
    //(where T is a subclass of S::value_type).
    //
    //eg. SubclassPointerVector<PhyloNode, NodeVector>
    //    is a subclass of NodeVector, such that operator[] and
    //    begin() and end() work in terms of PhyloNode* rather
    //    than Node*.  So client code won't need to cast
    //    iterator dereferences or [] entry-lookups to PhyloNode*.
    //
    //
public:
    typedef S super;
    using typename S::size_type;
    
    class reference {
        private:
            S& to_vector;
            size_type at_index;
        public:
            reference(S& vector, size_type index):
                to_vector(vector), at_index(index) {}
            operator T*()   { return dynamic_cast<T*> ( to_vector[at_index] ) ; }
            T* operator->() { return dynamic_cast<T*> ( to_vector[at_index] ) ; }
            reference& operator= (T* new_value) {
                to_vector[at_index] = new_value;
                return *this;
            }
            reference& operator= (const reference& elsewhere) {
                to_vector[at_index] = elsewhere.to_vector[elsewhere.at_index];
                return *this;
            }
    };

    class iterator : public super::iterator {
        public:
            iterator( typename super::iterator i) : super::iterator(i) {};
            T* operator*() { return dynamic_cast<T*>( super::iterator::operator*() ); }
    };
    class const_iterator: public super::const_iterator {
        public:
            const_iterator( typename super::const_iterator i) : super::const_iterator(i) {};
            T* operator*() { return dynamic_cast<T*>( super::const_iterator::operator*() ); }
    };
    iterator begin() { return iterator( super::begin() ); }
    iterator end()   { return iterator( super::end() ); }
    const_iterator begin() const { return const_iterator ( super::begin() ); }
    const_iterator end() const { return const_iterator ( super::end() ); }
    T* operator[] (typename super::size_type i) const {
        return dynamic_cast<T*>( super::operator[] (i) );
    }
    reference operator[] (size_type i) {
        return reference(*(dynamic_cast<S*>(this)), i );
    }
};

/**
    Node vector
 */
typedef SubclassPointerVector<PhyloNode, NodeVector> PhyloNodeVector;

typedef pair<PhyloNode*, PhyloNode*> PhyloBranch;
typedef CastingVector<PhyloBranch, BranchVector> PhyloBranchVector;

#endif
