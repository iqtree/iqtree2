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
    void cleanUpAfterTaxonPlacement(const TaxonToPlace& taxon,
                                    PhyloTree* tree);
    static TaxonPlacementOptimizer* getNewTaxonPlacementOptimizer();
};

class BatchPlacementOptimizer {
public:
    BatchPlacementOptimizer();
    virtual ~BatchPlacementOptimizer();
    virtual void cleanUpAfterBatch(TaxaToPlace& taxa,
                                   int firstTaxon, int lastTaxon,
                                   PhyloTree* tree);
    static BatchPlacementOptimizer* getNewBatchPlacementOptimizer();
};

class GlobalPlacementOptimizer {
public:
    GlobalPlacementOptimizer();
    virtual ~GlobalPlacementOptimizer();
    void cleanUpAfterPlacement(PhyloTree* tree);
    static GlobalPlacementOptimizer* getNewGlobalPlacementOptimizer();
};


#endif /* placementoptimizer_h */
