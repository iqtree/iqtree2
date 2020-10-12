//
// placement.h
// Declares enumerated types and functions that are used, in
// multiple classes, for determining parsimony placement, and
// likelihood placement.
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef placement_h
#define placement_h

#include <string>

namespace Placement {

    enum CostFunction {
        MAXIMUM_PARSIMONY,           //maximum parsimony  (each taxon, each possible insertion place)
        SANKOFF_PARSIMONY,           //ditto, but using Sankoff parsimony
        MAXIMUM_LIKELIHOOD_MIDPOINT, //maximum likelihood (at midpoint of existing branch)
        MAXIMUM_LIKELIHOOD_ANYWHERE  //maximum likelihood (anyhwere in existing branch)
    };
    enum LocalCleanup {
        NO_LOCAL_CLEANUP
    };
    enum BatchCleanup {
        NO_BATCH_CLEANUP
    };
    enum GlobalCleanup {
        NO_GLOBAL_CLEANUP
    };

    std::string getIncrementalParameter(const char letter, const char* defaultValue);
    size_t getIncrementalParameter(const char letter, size_t defaultValue);
    size_t getNumberOfTaxaToRemove(size_t countOfTaxa);
    CostFunction getCostFunction();
    LocalCleanup getLocalCleanupAlgorithm();
    size_t getTaxaPerBatch(size_t totalTaxa);
    size_t getInsertsPerBatch(size_t totalTaxa, size_t taxaPerBatch) ;
    BatchCleanup getBatchCleanupAlgorithm();
    GlobalCleanup getGlobalCleanupAlgorithm();

}
#endif /* placement_h */
