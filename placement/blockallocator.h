//
// blockallocator.h
// Does allocation (and tracking of) parsimony, likelihood, and
// scale vectors.
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef blockallocator_hpp
#define blockallocator_hpp
#include <tree/phylotree.h>
#include <vector>

class BlockAllocator {
protected:
    PhyloTree& phylo_tree;           //the tree that owns the memory (BlockAllocator
                                     //doesn't own any memory)
    int        index_parsimony;      //the number of parsimony blocks allocated thus far
    size_t     nptn;                 //
    std::vector<UINT*> uint_blocks;  //blocks of locally-allocated unsigned integer to release

public:
    BlockAllocator(PhyloTree& tree,
                   int parsimonyIndex);
    virtual ~BlockAllocator();
    
    /** request a partial parsimony block
     @param[out] partial_pars - pointer to the block (to be set to the newly allocated block)*/
    void         allocateParsimonyBlock(UINT*& partial_pars);

    /** allocate partial parsimony (and perhaps likelihood and scalenum) blocks
     to a PhyloNeighbor instance (in the tree).
     @param[out] nei - pointer to the PhyloNeighbor instance*/
    virtual void allocateMemoryFor(PhyloNeighbor* nei);

    /** returns a reference to the PhyloTree,
     for which the BlockAllocator instance is tracking
     parsimony and/or likelihood and scalenum vectors
     @return a reference to the tree*/
    PhyloTree&   getTree();
    
    /** returns the number of partial parsimony blocks that
     have been assigned to PhyloNeighbor instances;
     @return the number of assigned parsimony blocks*/
    int          getParsimonyBlockCount() const;
    
    /** hands over computed partial parsimony and/or likelihood
        vector block from one PhyloNeighbor instance to another.
        @param from_nei - the donor
        @param to_nei - the recipient */
    virtual void handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei);
    
    /** returns the number of partial parsimony blocks that
     have been assigned to PhyloNeighbor instances
     @return the number of assigned parsimony blocks*/
    virtual int  getLikelihoodBlockCount() const;
    
    /** may initialize pointers to partial likelihood and scale vectors
        (if this block allocator handles likelihood resources).
        @param partial_lh address of a partial likelihood vector pointer
                (only to be initialized if the pointer is null)
        @param scale_num address of a scale vector pointer
                (only to be initialized if *partial_lh is null)*/
    virtual void allocateLikelihoodBlocks(double*& partial_lh, UBYTE*& scale_num);

    /** indicates if this block allocator is keeping track of
        likelihood and scalenum vectors
         @return true if it is, false if not*/
    virtual bool usesLikelihood();

    /** readies the entire tree for the calculation of partial parsimony
        or partial likelihood, from first oriented toward second, and vice
        versa, for determining branch parsimony, branch length, or
        branch likelihood.
         @param first     the first node
         @param second   the second node*/
    virtual void makeTreeReady(PhyloNode* first, PhyloNode* second);
};

class LikelihoodBlockPair {
    public:
        double* partial_lh;
        UBYTE*  scale_num;
        bool    blocks_are_owned;
        LikelihoodBlockPair();
        LikelihoodBlockPair(const LikelihoodBlockPair& rhs);
        LikelihoodBlockPair& operator = (const LikelihoodBlockPair& rhs);
        LikelihoodBlockPair(double* lh, UBYTE* scale, bool owned=false);
        void    clear();
        virtual ~LikelihoodBlockPair();
        void    copyFrom(PhyloTree& tree, PhyloNeighbor* source);
        void    allocate(PhyloTree& tree);
        void    lendTo(PhyloNeighbor* borrower);
};

typedef std::vector<LikelihoodBlockPair> LikelihoodBlockPairs;

class LikelihoodBlockAllocator: public BlockAllocator {
protected:
    int                  index_lh;        //the number of likelihood blocks allocated so far
    std::vector<double*> double_blocks;
    std::vector<UBYTE*>  ubyte_blocks;
    std::vector<LikelihoodBlockPair> spare_block_pairs;

public:
    typedef BlockAllocator super;
    LikelihoodBlockAllocator(PhyloTree& tree, int parsimonyIndex, int likelihoodIndex);
    virtual      ~LikelihoodBlockAllocator();
    int          getLikelihoodBlockCount() const;
    virtual void allocateLikelihoodBlocks(double*& partial_lh, UBYTE*& scale_num);
    virtual void allocateMemoryFor(PhyloNeighbor* nei);
    virtual void handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei);

    /** indicates if this block allocator is keeping track of
        likelihood and scalenum vectors
         @return true if it is (for LikelihoodBlocAllocator, it always is), 
         false if not*/
    virtual bool usesLikelihood();
    
    /** Searches a subtree (on the side of second opposite first)
        for PhyloNeighbor instances that do not have
        likelihood vectors allocated, that will need them,
        if the partial likelihood of first->findNeighbor(second)
        is to be calculated.
     
     @param first    a node outside the subtree
     @param second  the node that marks the top of the subtree
     @param[out] list - the list of PhyloNeighbor
     */
    void findMissingInputs(PhyloNode* first, PhyloNode* second,
                           PhyloNeighborVec& list);

    /** Searches a subtree (on the side of second opposite first)
        for PhyloNeighbor instances that do have
        likelihood vectors allocated, that will *not* need them,
        during the calculation of the partial likelihood of
        first->findNeighbor(second)
     
     @param first    a node outside the subtree
     @param second  the node that marks the top of the subtree
     @param[out] list - the list of PhyloNeighbor
     */
    void findUselessOutputs(PhyloNode* first, PhyloNode* second,
                            PhyloNeighborVec& list);
    
    /**
     Ensure that the partial_lh (likelihood) and scale_num vectors
     that are needed (in a tree that uses per-node allocation) are
     found on the PhyloNeighbor instances that will want them
     (and nowhere else).
     (Sort of like what reorientPartialLh does, only tree-wide).

     perform tree search
     @return best likelihood found
     */

    virtual void makeTreeReady(PhyloNode* first, PhyloNode* second);
};

#endif /* blockallocator_h */
