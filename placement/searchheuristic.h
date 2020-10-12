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

class TaxonToPlace;
class TargetBranchRef;
class SearchHeuristic {
public:
    virtual ~SearchHeuristic();
    virtual bool isPlacementWorthTrying(const TaxonToPlace* taxon,
                                        const TargetBranchRef& target );
    virtual bool isGlobalSearch();

    static SearchHeuristic* getSearchHeuristic();
};

#endif /* searchheuristic_h */
