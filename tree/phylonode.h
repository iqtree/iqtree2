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
#include <utils/subclasspointervector.h>

#define FOR_EACH_ADJACENT_PHYLO_NODE(mynode, mydad, it, mychild) \
    for (PhyloNode* mychild=nullptr, *child2x=(mynode); child2x!=nullptr; child2x=nullptr) \
        for (auto it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
            if ((mychild = dynamic_cast<PhyloNode*>((*it)->node)) && mychild != (mydad) )

#define FOR_EACH_PHYLO_NEIGHBOR(mynode, mydad, it, nei) \
    for (PhyloNeighbor* nei=nullptr, *nei2x=(dynamic_cast<PhyloNode*>(mynode))->firstNeighbor(); nei2x!=nullptr ; nei2x=nullptr) \
        for (auto it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
                if ((nei = dynamic_cast<PhyloNeighbor*>(*it)) && nei->getNode() != (mydad) )

#define GET_OTHER_ADJACENT_PHYLO_NODES(node, dad, child_one, child_two) \
if (1) { \
    child_one = child_two = nullptr; \
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) { \
        if (child_one==nullptr) child_one = child; else child_two = child; \
    } \
} else 0

std::string pointer_to_hex(const void *ptr);

template <class T> std::string array_to_string(T* base, size_t N) {
    std::stringstream s;
    s << "{ ";
    for (size_t r = 0; r<N; ++r) {
        if (0<r) s << ", ";
        s << base[r];
    }
    s << " }";
    return s.str();
}

typedef unsigned short UBYTE;

/**
 * direction of a Neighbor from the root, for rooted tree only
 */
enum class RootDirection {
    UNDEFINED_DIRECTION, 
    TOWARD_ROOT, 
    AWAY_FROM_ROOT
};

class PhyloNeighbor;
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
    explicit PhyloNode(int aid);

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
    virtual void addNeighbor(Node *node, double length, int id = -1) override;

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

    bool hasNeighbor(PhyloNode* node);
    
    virtual PhyloNeighbor* firstNeighbor() const;
    
    virtual PhyloNeighbor* lastNeighbor() const;

    PhyloNeighbor* getNeighborByIndex(size_t index);
    
    int  getSubsetNumber() const;
    void setSubsetNumber(int subset_number);

protected:
    int subset;
};

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
    friend class LikelihoodBlockAllocator;
    friend class LikelihoodBlockPair;
    friend class LikelihoodCostCalculator;
    friend class PlacementTraversalInfo;
    friend class ParsimonyMatrix;
    friend class NNIContext;
    friend class ParsimonyNNIMove;
    friend class ParsimonyRouter;
    
public:
    friend class TinaTree;
    friend class PhyloSuperTreePlen;

    /**
        construct class with a node and length		
        @param anode the other end of the branch
        @param alength length of branch
     */
    PhyloNeighbor(Node *anode, double alength);

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param aid branch ID
     */
    PhyloNeighbor(Node *anode, double alength, int aid);

    /**
     construct class with another Neighbor
     @param nei another Neighbor
     */
    PhyloNeighbor(const PhyloNeighbor& nei);

    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual PhyloNeighbor* newNeighbor() const override;

    /**
        tell that the partial likelihood vector is not computed
     */
    inline void clearPartialLh() {
        partial_lh_computed &= ~LIKELIHOOD_IS_COMPUTED;
    }

    /**
     *  tell that the partial likelihood vector is computed
     */
    inline void unclearPartialLh() {
        partial_lh_computed |= LIKELIHOOD_IS_COMPUTED;
    }

    /**
        clear all partial likelihood recursively in forward direction
        @param dad dad of this neighbor
     */
    void clearForwardPartialLh(PhyloNode *dad);
    
	/**
	* For Upper Bounds analysis: get partial likelihood and lh scale factor
	*/
	inline double* get_partial_lh(){
        return partial_lh;
	}

	inline double get_lh_scale_factor(){
	    return lh_scale_factor;
	}

	inline int get_partial_lh_computed(){
        return partial_lh_computed;
	}

    /**
        tell that the partial parsimony vector is not computed
     */
    inline void clearPartialParsimony() {
        partial_lh_computed &= ~PARSIMONY_IS_COMPUTED;
    }

    /**
     *  tell that the partial likelihood vector *is* computed
     */
    inline void unclearPartialParsimony() {
        partial_lh_computed |= PARSIMONY_IS_COMPUTED;
    }

	/**
	 * true if this Neighbor is directed towards the root
	 */
	inline bool isTowardsRoot() {
		ASSERT(direction != RootDirection::UNDEFINED_DIRECTION);
		return (direction == RootDirection::TOWARD_ROOT);
	}

    inline int getSize() {
        return size;
    }
    
    PhyloNode* getNode() const;
    
    bool isLikelihoodComputed() const;
    
    void setLikelihoodComputed(bool set);
    
    bool isParsimonyComputed() const;
           
    void setParsimonyComputed(bool set);
    
    UINT* get_partial_pars();
    
    void clearComputedFlags();

    /**
     copy partial likelihood and partial parsimony information
     from another PhyloNeighbor.
     */
    void copyComputedState(const PhyloNeighbor* donor);

private:

    /**
        indicates whether the partial likelihood (and/or) the parsimony,
        is computed (and up to date).
     */
    int partial_lh_computed;
    
    static const int NOTHING_IS_COMPUTED    = 0;
    static const int LIKELIHOOD_IS_COMPUTED = 1;
    static const int PARSIMONY_IS_COMPUTED  = 2;

private:
    /**
        vector containing the partial likelihoods
     */

    double* partial_lh;

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
    PhyloNeighbor and PhyloNode vectors
 */
typedef SubclassPointerVector<PhyloNeighbor, NeighborVec> PhyloNeighborVec;
typedef SubclassPointerVector<PhyloNode, NodeVector> PhyloNodeVector;

struct PhyloBranch: public pair<PhyloNode*, PhyloNode*> {
    typedef pair<PhyloNode*, PhyloNode*> super;
    PhyloBranch();
    PhyloBranch(Node* left, Node* right);
    PhyloBranch(PhyloNode* left, PhyloNode* right);
    PhyloBranch& operator=(const PhyloBranch& rhs) = default;
    explicit PhyloBranch(const Branch& copyMe );
    PhyloBranch& operator=(const Branch& copyMe );
    operator Branch() const;

    int            getBranchID()       const;
    PhyloNeighbor* getLeftNeighbor()   const;
    PhyloNeighbor* getRightNeighbor()  const;
    bool           stillExists()       const;
    bool           isInnerBranch()     const;
    bool           isABranch()         const;
    bool           isDivergentBranch() const;
};

typedef CastingVector<PhyloBranch, BranchVector> PhyloBranchVector;

#endif
