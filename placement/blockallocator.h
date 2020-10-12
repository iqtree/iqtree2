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
    uint64_t   parsimony_block_size; //the size of an (interior) parsimony vector (in UINTs)
    uint64_t   lh_block_size;        //the size of an (interior) likelihood vector (in doubles)
    uint64_t   scale_block_size;     //the size of an (interior) scale num vector (in UBYTEs)
    int        index_parsimony;      //the number of parsimony blocks allocated thus far
    size_t     nptn;                 //
    std::vector<UINT*> uint_blocks;  //blocks of locally-allocated unsigned integer to release

public:
    BlockAllocator(PhyloTree& tree,
                   int parsimonyIndex);
    virtual ~BlockAllocator();
    
    /** request a partial parsimony block
     @params[out] pointer to the block (to be set to the newly allocated block)*/
    void         allocateParsimonyBlock(UINT*& partial_pars);

    /** allocate partial parsimony (and perhaps likelihood and scalenum) blocks
     to a PhyloNeighbor instance (in the tree).
     @params[out] pointer to the PhyloNeighbor instance*/
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
    
    virtual void handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei);
    
    /** returns the number of partial parsimony blocks that
     have been assigned to PhyloNeighbor instances
     @return the number of assigned parsimony blocks*/
    virtual int  getLikelihoodBlockCount() const;
    
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

class LikelihoodBlockAllocator: public BlockAllocator {
protected:
    int                  index_lh;        //the number of likelihood blocks allocated so far
    std::vector<double*> double_blocks;
    std::vector<UBYTE*>  ubyte_blocks;
    struct spare_block_pair {
        public:
            double* partial_lh;
            UBYTE*  scale_num;
            spare_block_pair(double* lh, UBYTE* scale);
    };
    std::vector<spare_block_pair> spare_block_pairs;

public:
    typedef BlockAllocator super;
    LikelihoodBlockAllocator(PhyloTree& tree, int parsimonyIndex, int likelihoodIndex);
    virtual      ~LikelihoodBlockAllocator();
    int          getLikelihoodBlockCount() const;
    virtual void allocateLikelihoodBlocks(double*& partial_lh, UBYTE*& scale_num);
    virtual void allocateMemoryFor(PhyloNeighbor* nei);
    virtual void handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei);
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
