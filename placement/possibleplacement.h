//
// possibleplacement.h
// Defines the PossiblePlacement class, which records
// how a TaxonToPlace would be attached to a tree, with
// its interior node linked into the middle of a specific
// target branch (target_branch member).
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef possibleplacement_hpp
#define possibleplacement_hpp

#include "targetbranch.h"

class PossiblePlacement {
public:
    TargetBranchRef  target_branch;  //
    int              parsimony_score;//the parsimony score (+ve)
    double           score;          //either +ve parsimony score, or minus
                                     //the log-likelihood (lower is better).
    double           lenToNewTaxon;  //(best scoring) length of the edge between
                                     //new_taxon and added_node
    double           lenToNode1;     //(best scoring) length of edge between
                                     //target_dad and added_node
    double           lenToNode2;     //(best scoring) length of edge between
                                     //target_child and added_node
    
    PossiblePlacement();
    PossiblePlacement& operator =  ( const PossiblePlacement& rhs );
    bool               operator <  ( const PossiblePlacement& rhs ) const;
    bool               operator <= ( const PossiblePlacement& rhs ) const;
    void               setTargetBranch(TargetBranchRange* targetRange, size_t index);
    void               setTargetBranch(TargetBranchRef& branch_ref);
    bool               canStillUse()    const;
    TargetBranch*      getTarget()      const;
    size_t             getTargetIndex() const;
    void               forget();
};

#include <stdio.h>

#endif /* possibleplacement_h */
