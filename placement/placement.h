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
    enum LocalOptimization {
        NO_LOCAL_OPTIMIZATION
    };
    enum BatchOptimization {
        NO_BATCH_OPTIMIZATION
    };
    enum GlobalOptimization {
        NO_GLOBAL_OPTIMIZATION
    };

    /** determine the value of a incremental placement sub-parameter, that is part 
        of the -incremental parameter passed to the command line.  Each sub parameter 
        is a single letter, followed by a value.  For example, 'B' for batch size, 
        followed by a number (of taxa to add per batch) or a percentage (of remaining
        taxa to add in the next batch.  Sub-parameters are separated by commas.
        If a subparameter contains a comma, curly braces { } can be used to ensure
        that any enclosed commas aren't treated as a sub-parameter separator
        (the outermost pair of curly braces is removed).
        @param letter       - the letter that indicates the parameter
        @param defaultValue - the default value of the parameter
        @return the value of the parameter*/
    std::string   getIncrementalParameter(const char letter, const char* defaultValue);

    /** determine the value of a numeric incremental placement sub-parameter       
        @param letter       - the letter that indicates the parameter
        @param defaultValue - the default value of the parameter
        @return the value of the parameter*/
    size_t        getIncrementalParameter(const char letter, size_t defaultValue);

    /** determine how many (if any) taxa are to be removed, and reinserted
        @param countOfTaxa - the total number of taxa in the tree
        @return the number of taxa to be removed, based on the R sub-parameter,
                which can be specified as a number or a percentage.*/
    size_t        getNumberOfTaxaToRemoveAndReinsert(size_t countOfTaxa);

    /** determine what cost function is to be used for taxon placement
        @return the cost function*/
    CostFunction  getCostFunction();

    /** determine how many taxa should be processed in the next batch.
        (this *cannot* be specified as a percentage of totalTaxa)
        @param  totalTaxa - the number of taxa still to be processed
        @return the number of taxa to be processed in the next batch
                (always at least 1, even if totalTaxa was zero)*/
    size_t        getTaxaPerBatch(size_t totalTaxa);

    /** determine how many taxa should be inserted in the next batch
        (this can be specified as either a number or a percentage)
        @param  totalTaxa - the number of taxa still to be processed
                (including those in the current batch)
        @param  taxaPerBatch - the number of taxa to be processed 
                (in the current batch) (it is assumed this as least 1)
        @return the number of taxa, from the current batch, that are to
                be added to the phylogenetic tree.  Always at least 1*/
    size_t        getInsertsPerBatch(size_t totalTaxa, size_t taxaPerBatch) ;

    LocalOptimization  getLocalOptimizationAlgorithm();
    BatchOptimization  getBatchOptimizationAlgorithm();
    GlobalOptimization getGlobalOptimizationAlgorithm();
}
#endif /* placement_h */
