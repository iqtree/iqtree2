//
// searchheuristic.cpp
// Implementation of the SearchHeuristic class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "searchheuristic.h"

SearchHeuristic::~SearchHeuristic() = default;
bool SearchHeuristic::isPlacementWorthTrying
    ( const TaxonToPlace* taxon,
      const TargetBranchRef& target ) {
        return true;
}
bool SearchHeuristic::isGlobalSearch() {
    return true;
}
SearchHeuristic* SearchHeuristic::getSearchHeuristic() {
    return new SearchHeuristic;
}

