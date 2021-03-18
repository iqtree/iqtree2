//
// searchheuristic.cpp
// Implementation of the SearchHeuristic class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include <utils/heapsort.h>
#include "searchheuristic.h"
#include "targetbranch.h"
#include "placementcostcalculator.h"
#include "taxontoplace.h"             //for TaxaToPlace

SearchHeuristic::~SearchHeuristic() = default;

bool SearchHeuristic::isGlobalSearch() const {
    return true;
}

bool SearchHeuristic::usesLikelihood() const {
    return false;
}

void SearchHeuristic::prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                      intptr_t startTarget, intptr_t stopTarget,
                                      TaxaToPlace& taxa,
                                      intptr_t startTaxon, intptr_t stopTaxon) {
}

bool SearchHeuristic::isPlacementWorthTrying(const TaxonToPlace& taxon,
                                             size_t taxonIndex, /* not id, index into TaxaToPlace */
                                             const TargetBranchRef& target ) {
        return true;
}

void SearchHeuristic::doneFiltering() {
}

SearchHeuristic* SearchHeuristic::getSearchHeuristic
                 (const PlacementParameters &placement_params) {
    auto heuristic = placement_params.getIncrementalParameter('H', "");
    if (heuristic=="" || heuristic=="0") {
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

bool BaseballSearchHeuristic::isGlobalSearch() const {
    return false;
}

bool BaseballSearchHeuristic::usesLikelihood() const {
    return true;
}

BaseballSearchHeuristic::~BaseballSearchHeuristic() {
    delete calculator;
}

void BaseballSearchHeuristic::prepareToFilter(PhyloTree& tree, TargetBranchRange& targets,
                                              intptr_t startTarget, intptr_t stopTarget,
                                              TaxaToPlace& taxa,
                                              intptr_t startTaxon, intptr_t stopTaxon) {
    target_base = startTarget;
    taxon_base  = startTaxon;
    is_worth_trying.setDimensions( stopTarget-startTarget, stopTaxon-startTaxon) ;
    Matrix<double> scores;
    scores.setDimensions( stopTarget-startTarget, stopTaxon-startTaxon);
    LikelihoodBlockPairs blocks(2);
    double parsimony_score = -1.0;
    for (intptr_t t = startTarget; t<stopTarget; ++t ) { //branch
        targets.getTargetBranch(t)->computeState(tree, parsimony_score, t, blocks);
        double* scoreRow = scores.getRow(t-startTarget);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t c = startTaxon; c<stopTaxon; ++c ) { //candidate taxon
            PossiblePlacement p;
            p.setTargetBranch(&targets, t);
            calculator->assessPlacementCost(tree, taxa.getTaxonByIndex(c), p);
            scoreRow[c] = p.score;
        }
    }
    intptr_t taxon_count  = stopTaxon  -  startTaxon;
    intptr_t target_count = stopTarget -  startTarget;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t c = 0; c<taxon_count; ++c ) { //candidate taxon (0-based)
        DoubleVector scoresForTaxon;
        scores.appendColumnToVector(c, scoresForTaxon);
        std::vector<intptr_t> targetIndices;
        for (intptr_t b = 0; b < target_count ; ++b ) {
            targetIndices.emplace_back(b);
        }
        mirroredHeapsort( scoresForTaxon, targetIndices );
        intptr_t take = static_cast<intptr_t> ( floor(sqrt(target_count)) );
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
        for (intptr_t t = 0; t < take ; ++t) {
            intptr_t b = targetIndices[t];
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
            for (intptr_t  t=0; t<take; ++t) {
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
