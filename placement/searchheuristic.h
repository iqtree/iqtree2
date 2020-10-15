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
    virtual void prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                 size_t startTarget, size_t stopTarget,
                                 TaxaToPlace& taxa, size_t startTaxon, size_t stopTaxon);
    virtual bool isPlacementWorthTrying(const TaxonToPlace* taxon,
                                        size_t taxon_index,
                                        const TargetBranchRef& target );
    virtual void doneFiltering();
    static SearchHeuristic* getSearchHeuristic();
};

class PlacementCostCalculator;

class BaseballSearchHeuristic : public SearchHeuristic {
private:
    PlacementCostCalculator* calculator;
    size_t       taxon_base;
    size_t       target_base;
    Matrix<bool> is_worth_trying;
public:
    BaseballSearchHeuristic(PlacementCostCalculator* calculatorToUse);
    virtual ~BaseballSearchHeuristic();
    virtual bool isGlobalSearch();
    virtual void prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                 size_t startTarget, size_t stopTarget,
                                 TaxaToPlace& taxa, size_t startTaxon, size_t stopTaxon);
    virtual bool isPlacementWorthTrying(const  TaxonToPlace* taxon,
                                        size_t taxonIndex,
                                        const  TargetBranchRef& target );
    virtual void doneFiltering();
};

#endif /* searchheuristic_h */
