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
#include "blockallocator.h"

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
class TaxaToPlace;

struct BenefitPair {
public:
    bool   hasForwardBenefit;
    double forwardBenefit;
    bool   hasBackwardBenefit;
    double backwardBenefit;
    BenefitPair(): hasForwardBenefit(false),  forwardBenefit(0)
              , hasBackwardBenefit(false), backwardBenefit(0) {}
};

class TargetBranch : public std::pair<PhyloNode*, PhyloNode*> {
    //A place where a node could be inserted, with likelihood and
    //partial parsimony determined, looking into the tree from the
    //insertion point.
    BlockAllocator* blocker;
    UINT*           partial_pars;
    double          connection_cost;     //The cost of adding a node in the middle of the branch
    double          branch_cost;         //The branch cost (typically, the parsimony score) of
                                         //the first<-->second branch).
    double*         partial_lh;
    UBYTE*          scale_num;
    mutable double  branch_lh_scale_factor;
    bool            used;
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
    double computeState (PhyloTree& phylo_tree,
                         intptr_t target_branch_index,
                         LikelihoodBlockPairs &blocks);
    void   dumpNeighbor (VerboseMode level, const char* prefix,
                         PhyloTree& phylo_tree, PhyloNeighbor* nei) const;
    void   updateMapping(intptr_t branch_id,
                         PhyloNode* updated_first,
                         PhyloNode* updated_second);
    void forgetState()            const;
    bool isUsedUp()               const;
    void handOverComputedStateTo(PhyloNeighbor* nei) ;
    UINT*   getParsimonyBlock()   const;
    double* getLikelihoodBlock()  const;
    UBYTE*  getScaleNumBlock()    const;
    double  getLhScaleFactor()    const;
    void    setLhScaleFactor(double v);

    void costPlacementOfTaxa(PhyloTree& tree,
                             TargetBranchRange& targets,
                             size_t targetNumber,
                             TaxaToPlace& candidates,
                             intptr_t candidateStartIndex,
                             intptr_t candidateStopIndex,
                             SearchHeuristic*   heuristic,
                             PlacementCostCalculator* calculator,
                             bool isFirstTargetBranch) const;
    
    void takeOwnershipOfReplacementVector(ReplacementBranchList* branches);
    ReplacementBranchList* getReplacements();
    
    //The following functions are used for parsimony rearrangement
    
    double getBranchCost()        const;
    double getFirstSubTreeCost()  const;
    double getSecondSubTreeCost() const;
    double getConnectionCost()    const;

    /** Calculate the benefit (in terms of reduced state changes required
     according to maximum parsimony), if this branch is disconnected
     @param phylo_tree
     @param other_branch
     @return  If this branch is AB, and the other is EF, returns
            cost(AC)+cost(BC)+cost(CD)+cost(CE)+cost(CF) - cost(AB) - cost(EF)
     
     */
    double      getFullDisconnectionBenefit    (const PhyloTree& phylo_tree,
                                                const TargetBranchRange& branches) const;
    BenefitPair getPartialDisconnectionBenefit (const PhyloTree& phylo_tree,
                                                const TargetBranchRange& branches) const;
    
    double      getFullConnectionCost    (const PhyloTree& phylo_tree,
                                          const TargetBranch& other_branch) const;
    double      getForwardConnectionCost (const PhyloTree& phylo_tree,
                                          const TargetBranch& other_branch) const;
    double      getBackwardConnectionCost(const PhyloTree& phylo_tree,
                                          const TargetBranch& other_branch) const;
    bool        isExternalBranch() const;
};

class TargetBranchRef;

class TargetBranchRange : public vector<TargetBranch> {
public:
    typedef  vector<TargetBranch> super;
    TargetBranchRange(PhyloTree& phylo_tree, BlockAllocator* b,
                      PlacementCostCalculator* calculator,
                      bool match_branch_numbers);
    TargetBranch* getTargetBranch(size_t i) {
        return &at(i);
    }
    void removeUsed();
    TargetBranchRef addNewRef(BlockAllocator& allocator,
                              LikelihoodBlockPairs& blocks,
                              PhyloNode* node1, PhyloNode* node2,
                              bool likelihood_wanted);
};

#endif /* targetbranch_h */
