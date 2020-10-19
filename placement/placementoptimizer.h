//
// placementoptimizer.h
// Defines the TaxonPlacementOptimizer class (which is
// called after each taxon has been placed), the
// BatchPlacementOptimizer class (called after each batch
// of new taxa has been placed), and the
// GlobalPlacementOptimizer class (called after all new
// taxa have been placed).
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef placementoptimizer_h
#define placementoptimizer_h

#include <tree/phylotree.h>
#include "taxontoplace.h"

class TaxonPlacementOptimizer {
public:
    TaxonPlacementOptimizer();
    virtual ~TaxonPlacementOptimizer();
    void cleanUpAfterTaxonPlacement(TaxaToPlace& taxa,
                                    size_t taxon_index,
                                    TargetBranchRange& targets,
                                    PhyloTree& tree);
    
    /** allocate a new TaxonPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @return a new TaxonPlacementOptimizer instance (it is up to the caller to delete it) */
    static TaxonPlacementOptimizer* getNewTaxonPlacementOptimizer();
};

class BatchPlacementOptimizer {
public:
    BatchPlacementOptimizer();
    virtual ~BatchPlacementOptimizer();
    virtual void cleanUpAfterBatch(TaxaToPlace& taxa,
                                   int start_taxon_index,
                                   int stop_taxon_index,
                                   TargetBranchRange& targets,
                                   PhyloTree& tree);
    
    /** allocate a new BatchPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @return a new BatchPlacementOptimizer instance (it is up to the caller to delete it) */
    static BatchPlacementOptimizer* getNewBatchPlacementOptimizer();
};

class GlobalPlacementOptimizer {
public:
    GlobalPlacementOptimizer();
    virtual ~GlobalPlacementOptimizer();
    void cleanUpAfterPlacement(PhyloTree& tree);
    
    /** allocate a new GlobalPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @return a new GlobalPlacementOptimizer instance (it is up to the caller to delete it) */
    static GlobalPlacementOptimizer* getNewGlobalPlacementOptimizer();
};

#endif /* placementoptimizer_h */
