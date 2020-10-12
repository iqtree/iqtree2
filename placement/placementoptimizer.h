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
    static TaxonPlacementOptimizer* getTaxonPlacementOptimizer();
};

class BatchPlacementOptimizer {
public:
    BatchPlacementOptimizer();
    virtual ~BatchPlacementOptimizer();
    template <class T> void cleanUpAfterBatch(TaxaToPlace<T>& taxa,
                           int firstTaxon, int lastTaxon,
                           PhyloTree* tree) {
        if (VB_MIN <= verbose_mode) {
            std::stringstream s;
            s << "Processed batch of "
              << (lastTaxon - firstTaxon) << " taxa";
            tree->logLine(s.str() );
        }
    }
    static BatchPlacementOptimizer* getBatchPlacementOptimizer();
};

class GlobalPlacementOptimizer {
public:
    GlobalPlacementOptimizer();
    virtual ~GlobalPlacementOptimizer();
    void cleanUpAfterPlacement(PhyloTree* tree);
    static GlobalPlacementOptimizer* getGlobalPlacementOptimizer();
};


#endif /* placementoptimizer_h */
