//
// searchheuristic.h
// Defines the SearchHeuristic class.  Subclasses will
// decide which taxa (TaxonToPlace instances) should be
// scored for placement against which target branches
// (TargetBranchRef instances).
//
//  Created by James Barbetti on 8-Oct-2020.
//

#ifndef searchheuristic_h
#define searchheuristic_h

#include <stdlib.h>
#include <utils/distancematrix.h>

class TaxonToPlace;
class TargetBranchRef;
class TargetBranchRange;
class TaxaToPlace;
class PhyloTree;

class SearchHeuristic {
public:
    virtual ~SearchHeuristic();
    virtual bool isGlobalSearch();
    
    /** called to give the SearchHeuristic a chance to set up (or allocate)
        internal structures, before it is asked whether placements should be
        considered.
     @param tree the phylo tree
     @param targets a range of target branches to be searched
     @param startTarget the index (in targets) of the first target branch to be consiered
     @param stopTarget one more than the index of the last target branch to be considered
     @param taxa a container (a TaxaToPlace) instance of taxa to insert
     @param startTaxon the index (in taxa) of the first taxon to be considered
     @param stopTaxon one more than the index of the last taxon to be considered*/
    virtual void prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                 size_t startTarget, size_t stopTarget,
                                 TaxaToPlace& taxa, size_t startTaxon, size_t stopTaxon);
    
    
    virtual bool isPlacementWorthTrying(const TaxonToPlace& taxon,
                                        size_t taxon_index,
                                        const TargetBranchRef& target );

    /** called to let the search heuristic know isPlacementWorthTrying won't be called
        until after prepareToFilter is called again (an opportunity for the SearchHeuristic
        to deallocate any resources it allocated when prepareToFilter was called.*/
    virtual void doneFiltering();
    static SearchHeuristic* getSearchHeuristic();
};

class PlacementCostCalculator;

class BaseballSearchHeuristic : public SearchHeuristic {
    
private:
    PlacementCostCalculator* calculator; //Owned. Deleted in ~SearchHeuristic destructor.
    size_t       taxon_base;
    size_t       target_base;
    Matrix<bool> is_worth_trying; //Rows are target branches (by index into a TargetBranchRange);
                                  //Columns are taxa (by index into TaxaToPlace).
    PhyloTree*   tree_in_use;

    //
    //1. prepareToFilter ... runs the calculator (ideally a cheap one! Parsimony?)
    //   for each combination of target branch & taxon to place,
    //   identifies the best (lowest) scoring target branches for each taxon, and
    //   marks is_worth_trying[target_index-target-base][taxon_index-taxon_base]
    //   with true for those (and false for the others).
    //2. isPlacementWorthTrying ... looks up is_worth_trying.
    //3. doneFiltering discards is_worth_trying.
    //

public:
    /**
     @param calculator a new PlacementCostCalulator instance
            (BaseballSearchHeuristic's destructor will delete it) */
    BaseballSearchHeuristic(PlacementCostCalculator* calculatorToUse);
    virtual ~BaseballSearchHeuristic();
    virtual bool isGlobalSearch();
    virtual void prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                 size_t startTarget, size_t stopTarget,
                                 TaxaToPlace& taxa, size_t startTaxon, size_t stopTaxon);
    virtual bool isPlacementWorthTrying(const  TaxonToPlace& taxon,
                                        size_t taxonIndex,
                                        const  TargetBranchRef& target );
    virtual void doneFiltering();
};

#endif /* searchheuristic_h */
