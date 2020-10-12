//
// targetbranch.h
// Defines the TargetBranchRange class, which describes
// a branch (in a PhyloTree) to which new taxa might be
// linked, during taxon placement.
//
// Defines the TargetBranchRange and TargetBranchRef classes.
// TargetBranchRange is a vector of target branches, and
// TargetBranchRef is a "by-index" reference to a target
// branch that is contained *in* a TargetBranchRange instance.
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef targetbranch_h
#define targetbranch_h

#include <vector>
#include <tree/phylonode.h>
#include <tree/phylotree.h>

class TargetBranch;
class TargetBranchRange;
class TargetBranchRef {
private:
    TargetBranchRange* target_range;
    size_t             target_index;
public:
    TargetBranchRef ( );
    TargetBranchRef ( const TargetBranchRef& r );
    TargetBranchRef ( TargetBranchRange* range, size_t index );
    TargetBranchRef& operator= ( const TargetBranchRef& r );
    bool          isUsedUp()       const;
    PhyloNode*    getFirst()       const;
    PhyloNode*    getSecond()      const;
    TargetBranch* getTarget()      const;
    size_t        getTargetIndex() const;
};

typedef std::vector<TargetBranchRef>
    ReplacementBranchList;

class SearchHeuristic;
class PlacementCostCalculator;
class TaxonToPlace;

class TargetBranch : public std::pair<PhyloNode*, PhyloNode*> {
    //A place where a node could be inserted, with likelihood and
    //partial parsimony determined, looking into the tree from the
    //insertion point.
    BlockAllocator* blocker;
    UINT*   partial_pars;
    double* partial_lh;
    UBYTE*  scale_num;
    bool    used;
    ReplacementBranchList* replacements;
    friend class TargetBranchRef;
public:
    typedef std::pair<PhyloNode*, PhyloNode*> super;
    TargetBranch();
    ~TargetBranch();
    TargetBranch(const TargetBranch& rhs);
    TargetBranch& operator= (const TargetBranch& rhs);
    TargetBranch(BlockAllocator* allocator,
                 PhyloNode* node1, PhyloNode* node2,
                 bool likelihood_wanted);
    void computeState(PhyloTree& phylo_tree) const;
    void forgetState() const;
    bool isUsedUp() const;
    void handOverComputedStateTo(PhyloNeighbor* nei) ;
    const UINT* getParsimonyBlock() const;
    const double* getLikelihoodBlock() const;
    const UBYTE* getScaleNumBlock() const;
    
    template <class T>
    void costPlacementOfTaxa(PhyloTree& tree,
                             TargetBranchRange* targets,
                             size_t targetNumber,
                             T* candidateStart,
                             T* candidateStop,
                             SearchHeuristic*   heuristic,
                             PlacementCostCalculator* calculator,
                             bool isFirstTargetBranch) const;
    
    void takeOwnershipOfReplacementVector(ReplacementBranchList* branches);
    ReplacementBranchList* getReplacements();
};

class TargetBranchRef;

class TargetBranchRange : public vector<TargetBranch> {
public:
    typedef  vector<TargetBranch> super;
    TargetBranchRange(PhyloTree& phylo_tree, BlockAllocator* b,
                      PlacementCostCalculator* calculator);
    void removeUsed();
    TargetBranchRef addNewRef(BlockAllocator* allocator,
                              PhyloNode* node1, PhyloNode* node2,
                              bool likelihood_wanted);
};






#endif /* targetbranch_h */
