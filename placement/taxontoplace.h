//
// taxontoplace.h
// Defines the TaxonToPlace and LessFussyTaxon classes.
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef taxontoplace_hpp
#define taxontoplace_hpp

#include <string>
#include "possibleplacement.h"
#include <tree/phylotree.h>

class SearchHeuristic;
class PlacementCostCalculator;
class TaxonToPlace {
    //A taxon that could be added to a tree.
protected:
    PossiblePlacement bestPlacement;
public:
    int               taxonId;       //The id of the taxon   (new_leaf->id)
    std::string       taxonName;     //The name of the taxon (new_leaf->name)
    bool              inserted;      //true if this taxon has been inserted
    PhyloNode*        new_leaf;      //The new leaf node for this taxon
    PhyloNode*        new_interior;  //The new interior node that the leaf node
                                     //will be linked with.
    const UINT*       partial_pars;  //partial parsimony for new leaf, seen
                                     //from the new interior
    double*           partial_lh;    //partial likelihood and scaling info, also
    UBYTE*            scale_num;     //as seen from the new interior
    
    TaxonToPlace();
    TaxonToPlace(const TaxonToPlace& rhs); //copy everything!
    TaxonToPlace(BlockAllocator* ba, int id, std::string name);
    virtual      ~TaxonToPlace();
        
    const UINT*   getParsimonyBlock()  const;
    const double* getLikelihoodBlock() const;
    const UBYTE*  getScaleNumBlock()   const;
    int           getTaxonId()         const;
    
    virtual size_t considerPlacements(const PossiblePlacement* placements,
                                      size_t count);
    virtual bool   considerAdditionalPlacement(const PossiblePlacement& placement) ;
    virtual const  PossiblePlacement& getBestPlacement() const;
    virtual bool   canInsert() const;
    void findPlacement ( PhyloTree& phylo_tree,
                         size_t taxonIndex,
                         TargetBranchRange& range,
                         SearchHeuristic* heuristic,
                         const PlacementCostCalculator* calculator );
    
    bool operator <  ( const TaxonToPlace& rhs ) const;
    bool operator <= ( const TaxonToPlace& rhs ) const;
    void insertIntoTree ( PhyloTree& phylo_tree, BlockAllocator* b,
                          TargetBranchRange& dest,
                          PlacementCostCalculator& calculator);
    virtual void forgetGazumpedPlacements();
    bool insertNearby ( PhyloTree& phylo_tree, BlockAllocator* b,
                        TargetBranchRange& dest,
                        PlacementCostCalculator& calculator );
    void assessNewTargetBranches ( PhyloTree& phylo_tree,
                                   PlacementCostCalculator& calculator,
                                   TargetBranch* tb,
                                   std::vector<PossiblePlacement>& scores);
};

class TaxaToPlace {
public:
    virtual ~TaxaToPlace() = default;
    virtual TaxonToPlace* getTaxonByIndex(size_t index) = 0;
    virtual bool isEmpty() const = 0;
    virtual void sortBatch(size_t batchStart, size_t batchStop) = 0;
};

template <class T=TaxonToPlace> class TypedTaxaToPlace: public TaxaToPlace, public std::vector<T> {
public:
    typedef std::vector<T> super;
    explicit TypedTaxaToPlace(size_t reservation) {
        super::reserve(reservation);
    }
    virtual bool isEmpty() const                        { return std::vector<T>::empty(); }
    virtual TaxonToPlace* getTaxonByIndex(size_t index) { return & std::vector<T>::at(index); }
    virtual void sortBatch(size_t batchStart, size_t batchStop) {
        std::sort( std::vector<T>::begin() + batchStart, std::vector<T>::begin() + batchStop);

    }
    virtual ~TypedTaxaToPlace() = default;
};

class LessFussyTaxon: public TaxonToPlace {
protected:
    std::vector<PossiblePlacement> placementStore;
    static const size_t max_placements_to_keep = 5;
    //Note: this had better not be more than about 10, because
    //      if it were large, you'd want to maintain a heap.
public:
    typedef TaxonToPlace super;
    LessFussyTaxon ( );
    LessFussyTaxon ( const LessFussyTaxon& rhs );
    LessFussyTaxon ( BlockAllocator* ba, int id, std::string name);
    virtual size_t considerPlacements(const PossiblePlacement* placements,
                                      size_t count);
    virtual bool   considerAdditionalPlacement(const PossiblePlacement& placement);
    virtual void   forgetGazumpedPlacements();
};
#endif /* taxontoplace_h */
