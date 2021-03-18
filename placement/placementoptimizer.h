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
    void optimizeAfterTaxonPlacement(TaxaToPlace& taxa,
                                    size_t taxon_index,
                                    TargetBranchRange& targets,
                                    PhyloTree& tree);
    
    /** allocate a new TaxonPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @return a new TaxonPlacementOptimizer instance (it is up to the caller to delete it) */
    static TaxonPlacementOptimizer* getNewTaxonPlacementOptimizer();
};

class BatchPlacementOptimizer {
protected:
    bool be_quiet;
public:
    BatchPlacementOptimizer(bool be_quiet);
    virtual ~BatchPlacementOptimizer();
    virtual void optimizeAfterBatch(TaxaToPlace& taxa,
                                   intptr_t start_taxon_index,
                                   intptr_t stop_taxon_index,
                                   TargetBranchRange& targets,
                                   PhyloTree& tree);
    
    /** allocate a new BatchPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @params be_silent if true, don't log anything
        @return a new BatchPlacementOptimizer instance (it is up to the caller to delete it) */
    static BatchPlacementOptimizer* getNewBatchPlacementOptimizer(bool be_silent);
};

class GlobalPlacementOptimizer {
protected:
    bool be_quiet;
public:
    GlobalPlacementOptimizer(bool be_silent);
    virtual ~GlobalPlacementOptimizer();
    virtual void optimizeAfterPlacement(PhyloTree& tree);
    
    /** allocate a new GlobalPlacementOptimizer that matches what
        has been asked for in the -incremental parameter
        @return a new GlobalPlacementOptimizer instance
                (it is up to the caller to delete it) */
    static GlobalPlacementOptimizer* getNewGlobalPlacementOptimizer(bool use_likelihood, bool be_silent);
};

class GlobalLikelihoodPlacementOptimizer: public GlobalPlacementOptimizer {
public:
    typedef GlobalPlacementOptimizer super;
    GlobalLikelihoodPlacementOptimizer(bool be_silent);
    virtual ~GlobalLikelihoodPlacementOptimizer();
    virtual void optimizeAfterPlacement(PhyloTree& tree);
};

#endif /* placementoptimizer_h */
