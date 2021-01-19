//
//  parsimonymove.h
//  Created by James Barbetti on 18/1/21.
//

#ifndef parsimonymove_h
#define parsimonymove_h

#include <stdio.h>
#include "phylotree.h"
#include "iqtree.h"
#include <placement/targetbranch.h> //for TargetBranchRange

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
    
    bool   operator < (const ParsimonyMove& rhs)  const;
    bool   operator <= (const ParsimonyMove& rhs) const;
    double getBenefit() const;

    virtual void   initialize(intptr_t source_branch, bool beLazy) = 0;
    virtual bool   isStillPossible(const TargetBranchRange& branches,
                                 PhyloBranchVector& path) const = 0;
    virtual double recalculateBenefit(PhyloTree& tree,
                                      TargetBranchRange& branches,
                                      LikelihoodBlockPairs &blocks) = 0;
    virtual double apply(PhyloTree& tree, LikelihoodBlockPairs blocks,
                         TargetBranchRange& branches) = 0;
    virtual void   finalize(PhyloTree& tree,
                            const TargetBranchRange& branches) = 0;

    bool isNoLongerPossible(const TargetBranchRange& branches,
                            PhyloBranchVector& path) const;
    bool isAConnectedThroughBToC(PhyloNode* a, PhyloNode* b,
                                 PhyloNode* c, PhyloBranchVector& path) const;
};

#endif /* parsimonymove_hpp */
