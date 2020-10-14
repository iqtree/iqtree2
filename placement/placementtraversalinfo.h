//
// placementtraversalinfo.h
// The PlacementTraversalInfo class is a subclass of
// TraversalInfo (which is used for keeping track of
// branches that need to be processed during the
// calculation of a partial likelihood).
//
// During Likelihood placement, partial likelihoods
// have to be calculated "looking outward" from many
// different target branches (rather than, predominantly,
// "down" from a root), with the "point of perspective"
// potentially jumping around in the tree (not necessarily
// moving from one target branch to an adjacent branch).
// PlacementTraversalInfo helps with those "bigger jumps".
//
// Additionally, during likelihood placement, likelihoods
// need to be calculated for "grafted" trees, where a few
// (typically two) "new" nodes are "sort of" attached to the
// tree, and likelihood is calculated for the composite
// tree (from the point of the view of the "new" nodes).
// PlacementTraversalInfo provides member functions that
// help with those calculations.
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef placementtraversalinfo_hpp
#define placementtraversalinfo_hpp

#include <tree/phylotree.h>
#include <tree/phylonode.h>

class PlacementTraversalInfo: public TraversalInfo {
protected:
    PhyloTree& phylo_tree;
    LikelihoodBufferSet& buffers;
public:
    typedef TraversalInfo super;
    PlacementTraversalInfo(PhyloTree& tree,
                           LikelihoodBufferSet& buffersToUse,
                           PhyloNeighbor* dad_branch,
                           PhyloNode* dad);
    void computePartialLikelihood(PhyloNeighbor* nei=nullptr,
                                  PhyloNode* node=nullptr) ;
    double getBranchScore();
    ~PlacementTraversalInfo();
};

#endif /* placementtraversalinfo_h */
