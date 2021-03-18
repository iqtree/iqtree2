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

class PlacementParameters {
protected:
    std::string incremental_method;
        //Letter Meaning             Examples (comments)
        //A      batch optimisation  "SPR"
        //B      batch size          "100%"   (or a number)
        //C      cost function       "MP"     (maximum parsimony) (SMP == force sankoff)
        //H      heuristic method    "MP"     (cost function to use as heuristic)
        //                                    (only worth setting this if 'C' parameter is "ML")
        //I      inserts per batch   "100%"   (or a number)
        //L      local optimization           (not honoured)
        //R      taxa to remove      "10%"    (or a number) (only used for testing)
        //T      global optimisation "SPR"    (only blank and SPR (MP only) supported as yet).
        //
    
public:
    enum LocalOptimization {
        NO_LOCAL_OPTIMIZATION
    };
    enum BatchOptimization {
        NO_BATCH_OPTIMIZATION
    };
    enum GlobalOptimization {
        NO_GLOBAL_OPTIMIZATION
    };
    
    
    PlacementParameters();
    explicit PlacementParameters(const char* incremental_method_text);
    
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
    std::string   getIncrementalParameter(const char letter, const char* defaultValue) const;

    /** determine the value of a numeric incremental placement sub-parameter       
        @param letter       - the letter that indicates the parameter
        @param defaultValue - the default value of the parameter
        @return the value of the parameter*/
    size_t        getIncrementalSizeParameter(const char letter, size_t defaultValue) const;

    /** determine how many (if any) taxa are to be removed, and reinserted
        @param countOfTaxa - the total number of taxa in the tree
        @return the number of taxa to be removed, based on the R sub-parameter,
                which can be specified as a number or a percentage.*/
    size_t        getNumberOfTaxaToRemoveAndReinsert(size_t countOfTaxa) const;

    /** indicates whether placement will use parsimony (more to the point:
        will it need to use partial parsimony vectors
     @return true if it will*/
    bool doesPlacementUseParsimony() const;

    /** indicates whether placement will use Sankoff parsimony
     @return true if it will*/
    bool doesPlacementUseSankoffParsimony() const;

    /** indicates whether placement will use likelihood (more to the point:
     will it need to use partial likehood and scale num vectors).
     @return true if it will*/
    bool doesPlacementUseLikelihood() const;

    /** determine how many taxa should be processed in the next batch.
        (this *cannot* be specified as a percentage of totalTaxa)
        @param  totalTaxa - the number of taxa still to be processed
        @return the number of taxa to be processed in the next batch
                (always at least 1, even if totalTaxa was zero)*/
    size_t        getTaxaPerBatch(size_t totalTaxa) const;

    /** determine how many taxa should be inserted in the next batch
        (this can be specified as either a number or a percentage)
        @param  totalTaxa - the number of taxa still to be processed
                (including those in the current batch)
        @param  taxaPerBatch - the number of taxa to be processed 
                (in the current batch) (it is assumed this as least 1)
        @return the number of taxa, from the current batch, that are to
                be added to the phylogenetic tree.  Always at least 1*/
    size_t        getInsertsPerBatch(size_t totalTaxa, size_t taxaPerBatch) const;

    LocalOptimization  getLocalOptimizationAlgorithm()  const;
    BatchOptimization  getBatchOptimizationAlgorithm()  const;
    GlobalOptimization getGlobalOptimizationAlgorithm() const;
};
#endif /* placement_h */
