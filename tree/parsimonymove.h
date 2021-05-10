//
//  parsimonymove.h
//  Describes a possible parsimony move (in practice, one
//  found during parsimony SPR or parsimony TBR).
//
//  Created by James Barbetti on 18-Jan-2021.
//

#ifndef parsimonymove_h
#define parsimonymove_h

#include <stdio.h>
#include "phylotree.h"
#include "iqtree.h"
#include <placement/targetbranch.h> //for TargetBranchRange

class ParsimonyPathVector: public std::vector< std::vector<UINT*> > {
    intptr_t pv_per_thread;
    intptr_t min_thread_count;
    intptr_t thread_count;
public:
    ParsimonyPathVector() = delete;
    ParsimonyPathVector(intptr_t blocks, intptr_t min_threads, 
                        intptr_t threads_to_use);
    ~ParsimonyPathVector() = default;
    intptr_t getBlocksPerThread() const;
    intptr_t getNumberOfPathsRequired() const;
    intptr_t getTotalNumberOfBlocksRequired() const;
};

class ParsimonyMove {
public:
    bool     lazy;
    double   benefit;
    intptr_t source_branch_id;
    int64_t  positions_considered;

    ParsimonyMove();
    ParsimonyMove(const ParsimonyMove& rhs) = default;
    ParsimonyMove& operator=(const ParsimonyMove& rhs) = default;
    virtual ~ParsimonyMove() = default;
    
    bool   operator <  (const ParsimonyMove& rhs) const; //order by ascending benefit
    bool   operator <= (const ParsimonyMove& rhs) const;
    double getBenefit() const;
    
    /** indicates how many additional parsimony vectors need to be
     allocated, and made available to each thread, and passed in via
     path_parsimony, when findMove is called
     @param radius the search radius
     @return the number of parsimony vectors to allocate for each thread
     
        This needs to be replaced in your implementation.  The version
        in ParsimonyMove will generate an assertion failure via ASSERT(0).
     */
    static intptr_t getParsimonyVectorSize(intptr_t radius);

    static intptr_t getMinimumPathVectorCount();
    
    /** set up to search for a partial tree rearrangement (or move),
        involving the specified brach
     @param source_branch the id number of the branch (which will be
     used as one of the starting points for the search for a beneficial move)
     @param be_lazy true, if the search is to be lazy (is to skimp on
     parsimony scoring, and use heuristic cost/benefit functions instead).
     */
    virtual void   initialize(intptr_t source_branch_id,
                              bool be_lazy) = 0;
    /** search for a partial tree rearrangement (or move).
     @param tree a phylo tree
     @param branches an up-to-date (synchronized) vector of TargetBranch
            instances that supports random access (id #s of PhyloNeighbor
            instances agree with indexes of corresponding entries in the vector).
            The implementation may assume that computeState has been called
            on each branch *since* the last time the tree structure was changed.
     @param radius the search radius
     @param path_parsimony a vector (private to the currently executing thread
            that is large enough to store the partial parsimony information (if any)
            that is being calculated during local tree traversal.  how large, "large
            enough" is is indicated by getParsimonyVectorSize().
     @param parsimony_score the current parsimony score for the tree
     */
    virtual void   findMove(const PhyloTree& tree,
                            const TargetBranchRange& branches,
                            int radius,
                            std::vector<UINT*> &path_parsimony,
                            double parsimony_score) = 0;
    /** return a description of the most beneficial move found
     (implementations can assume that, this will only be called if
     a >0 benefit move *has* already been found; getBenefit() > 0). */
    virtual std::string getDescription() const = 0;
    
    /** called when searching (for partial tree arrangements) is done
     * (this is an opportunity for the ParsimonyMove implementation to record
     * information, about the current structure of the phylo tree, that it can
     * make use of in its isStillPossible() implementation).
     * @param tree a phylo tree
     * @param branches an up-to-date vector of TargetBranch instances
     *        that correspond 1:1 with the branches in the tree.
     */
    virtual void   finalize(PhyloTree& tree,
                            const TargetBranchRange& branches) = 0;
    
    /** determines whether the most beneficial move (if there was one) is still
     *  possible.  This should always return false if no beneficial move was found
     *  @param branches an up-to-date vector of TargetBranch instances
     *         that correspond 1:1 with the branches in the tree.
     *  @param path a PhyloBranchVector into which the "certifying path"
     *         (that proves the move is still possible) is to be written.
     *         (cleared if the return the move is no longer possible)
     *  @return true if the move is still possible
     *  @note I'm not sure that the path parameter is really needed by any caller.
     *        it might be worth removing it. -James B. 21-Jan-2021
     */
    virtual bool   isStillPossible(const TargetBranchRange& branches,
                                   PhyloBranchVector& path) const = 0;

    /** Recalculate the net benefit of a move
     * @param tree the phylo tree
     * @param tree_parsimony_score the current parsimony score of the
     *        tree (or -1 if that is not known).
     * @param branches an up-to-date vector of TargetBranch instances
     *        that correspond 1:1 with the branches in the tree.
     * @param blocks LikelihoodBlockPairs the might be needed to
     *        if likelihood is being taken into account (as yet, it never is but
     *        TargetBranch::computeState expects a LikelihoodBlockPairs
     *        passed to it, in case likelihood calculations are being taken into
     *        account) (a LikelihoodSPRMove class, if there were one, would
     *        genuinely need to be supplied some likelihood and scale vectors,
     *        via a "real" LikelihoodBlockPairs).
     *
     * @return the net benefit of the move, if it were applied to the tree as it now is
     */
    virtual double recalculateBenefit
                   ( PhyloTree& tree,
                     double tree_parsimony_score,
                     TargetBranchRange& branches,
                     LikelihoodBlockPairs &blocks,
                     ParsimonyPathVector& parsimony_path_vectors) const = 0;

    /** Apply a move to the tree (updating branches to match), in such a way
     *  that calling apply() a second time, with the same parameters, immediately
     *  afterwards, will exactly reverse the effect of the move.
     *  Callers may expect that parsimony scores will either be updated,
     *  or marked as being out of date, and should be able to trust the parsimony
     *  score that this member function returns.
     * @param tree the phylo tree
     * @param branches an up-to-date vector of TargetBranch instances
     *        that correspond 1:1 with the branches in the tree with
     *        PhyloNeighbor::id numbers matching the indexes into branches.
     * @param blocks (as for recalculateBenefit) LikelihoodBlockPairs
     * @return a parsimony score (it's declared double, because someday it might
     *         be minus one times the log of the tree likelihood, in subclasses that
     *         take account of and update likelihood)
     */
    virtual double apply(PhyloTree& tree,
                         double parsimony_score,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks,
                         ParsimonyPathVector& parsimony_path_vectors) = 0;

    /** Indicates whether a move can be applied (returns the true/false opposite
     *  of isStillPossible().
     * @param branches an up-to-date vector of TargetBranch instances
     *        that correspond 1:1 with the branches in the tree.
     * @param path a PhyloBranchVector into which the "certifying path"
     *        (that proves the move is still possible) is to be written
     *        (cleared if the return code is true, set if the return code is false)
     * @return false if the move can be applied, true if it can't
     */
    bool isNoLongerPossible(const TargetBranchRange& branches,
                            PhyloBranchVector& path) const;
    
    /** Indicates whether there is a path from node a to node c, through node b,
     *  and if there is one, returns it in path.  If there isn't path is cleared.
     *  @param a the starting node
     *  @param b a node that must be on the path adjacent to a
     *  @param c the destination node
     *  @param path the first path found (and, we can hope, the only path)
     *  @return true if there is a path
     *  @note the implementation should use a breadth-first, rather than a
     *  depth-first search, because movement radii tend to be limited
     *  (and, so, a *successful* breadth-first search for a short path should
     *  on average, return much more quickly than a *successful* depth-first
     *  search would) (performance in the "not connected case" is less important).
     */
    bool isAConnectedThroughBToC(PhyloNode* a,
                                 PhyloNode* b, PhyloNode* c,
                                 PhyloBranchVector& path) const;
};

#endif /* parsimonymove_h */
