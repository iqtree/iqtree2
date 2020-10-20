//
// searchheuristic.cpp
// Implementation of the SearchHeuristic class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "searchheuristic.h"
#include "targetbranch.h"
#include "placementcostcalculator.h"
#include "taxontoplace.h"             //for TaxaToPlace
#include <utils/heapsort.h>

SearchHeuristic::~SearchHeuristic() = default;

void SearchHeuristic::prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                      size_t startTarget, size_t stopTarget,
                                      TaxaToPlace& taxa,
                                      size_t startTaxon, size_t stopTaxon) {
}

bool SearchHeuristic::isPlacementWorthTrying(const TaxonToPlace& taxon,
                                             size_t taxonIndex, /* not id, index into TaxaToPlace */
                                             const TargetBranchRef& target ) {
        return true;
}

void SearchHeuristic::doneFiltering() {
}

bool SearchHeuristic::isGlobalSearch() {
    return true;
}
SearchHeuristic* SearchHeuristic::getSearchHeuristic() {
    auto heuristic = Placement::getIncrementalParameter('H', "");
    if (heuristic=="") {
        return new SearchHeuristic;
    }
    else if (heuristic=="MP") {
        return new BaseballSearchHeuristic(new ParsimonyCostCalculator(false));
    } else {
        std::stringstream s;
        s << "Did not recognize heuristic " << heuristic;
        outError(s.str());
        return nullptr;
    }
}

BaseballSearchHeuristic::BaseballSearchHeuristic(PlacementCostCalculator* calculatorToUse)
    : calculator(calculatorToUse), tree_in_use(nullptr) {
}

BaseballSearchHeuristic::~BaseballSearchHeuristic() {
    delete calculator;
}

bool BaseballSearchHeuristic::isGlobalSearch() {
    return false;
}

void BaseballSearchHeuristic::prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                              size_t startTarget, size_t stopTarget,
                                              TaxaToPlace& taxa,
                                              size_t startTaxon, size_t stopTaxon) {
    target_base = startTarget;
    taxon_base  = startTaxon;
    is_worth_trying.setDimensions( stopTarget-startTarget, stopTaxon-startTaxon) ;
    Matrix<double> scores;
    scores.setDimensions( stopTarget-startTarget, stopTaxon-startTaxon);
    for (size_t b = startTarget; b<stopTarget; ++b ) { //branch
        targets.getTargetBranch(b)->computeState(tree);
        double* scoreRow = scores.getRow(b-startTarget);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t c = startTaxon; c<stopTaxon; ++c ) { //candidate taxon
            PossiblePlacement p;
            p.setTargetBranch(&targets, b);
            calculator->assessPlacementCost(tree, taxa.getTaxonByIndex(c), p);
            scoreRow[c] = p.score;
        }
    }
    size_t taxon_count  = stopTaxon  -  startTaxon;
    size_t target_count = stopTarget -  startTarget;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (size_t c = 0; c<taxon_count; ++c ) { //candidate taxon (0-based)
        std::vector<double> scoresForTaxon;
        scores.appendColumnToVector(c, scoresForTaxon);
        std::vector<size_t> targetIndices;
        for (size_t b = 0; b < target_count ; ++b ) {
            targetIndices.emplace_back(b);
        }
        mirroredHeapsort( scoresForTaxon, targetIndices );
        size_t take = static_cast<size_t> ( floor(sqrt(target_count)) );
        //Always at least 1, never more than target_count.
        
        if ( scoresForTaxon[0] == scoresForTaxon[take-1] ) {
            //Might need more (if the best (take) are all tied)
            double same = scoresForTaxon[0];
            for ( ; take < target_count && scoresForTaxon[take] == same; ++take) {}
        } else {
            //Can take fewer if there are lots of high-scoring ties.
            double same = scoresForTaxon[take-1];
            for (; 1<take && scoresForTaxon[take-2] == same; --take) {}
        }
        for (size_t t = 0; t < take ; ++t) {
            size_t b = targetIndices[t];
            is_worth_trying.cell(b, c ) = true;
        }
        if (VB_DEBUG <= verbose_mode) {
            std::stringstream s;
            s << taxa.getTaxonByIndex(c).taxonName << " took top " << take << "target branches";
            tree.logLine(s.str());
            s.clear();
            if (take>3) {
                take=3;
                s<< "The top " << take << " were: ";
            } else {
                s<< "They were: ";
            }
            for (size_t t=0; t<take; ++t) {
                s << " " << targetIndices[t] << "(score " << scoresForTaxon[t] << ")";
            }
            tree.logLine(s.str());
        }
    }
    //Todo: need to tell the tree that we're doing stuff (report progress).
    //      The issue is, how do we count progress here, versus progress in the
    //      more expensive cost-calculation to which this is feeding "combinations
    //      worth trying".
    //
    tree_in_use = &tree;
}

bool BaseballSearchHeuristic::isPlacementWorthTrying(const TaxonToPlace& taxon, size_t taxonIndex,
                                    const TargetBranchRef& target ) {
    bool tryIt = is_worth_trying.cell(target.getTargetIndex() - target_base, taxonIndex-taxon_base);
    if (tryIt && tree_in_use!=nullptr) {
        TREE_LOG_LINE(*tree_in_use, VB_DEBUG, "Will try " << taxon.taxonName
                      << " against target branch " << target.getTargetIndex());
    }
    return tryIt;
}

void BaseballSearchHeuristic::doneFiltering() {
    is_worth_trying.clear();
    tree_in_use = nullptr;
}
