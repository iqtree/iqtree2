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

    //
    //Note 1: new_leaf, new_interior are "owned", but they are "given away" to
    //        the PhyloTree when the taxon is inserted into  the tree, so they
    //        aren't destroyed in the ~TaxonToPlace destructor.
    //        partial_pars, partial_lh, and scale_num are just copies of the
    //        same-named pointers on to/on new_leaf->findNeighbor(new_interior).
    //Note 2: taxon ID is a sequence number in an alignment.
    //        taxon index (often mentioned when talking about taxa in a
    //        TaxaToPlace container) is an index into *that* container.
    //        Not the same thing.
    //
    
    TaxonToPlace();
    TaxonToPlace(const TaxonToPlace& rhs); //copy everything!
    TaxonToPlace(BlockAllocator* ba, int id, std::string name);
    virtual      ~TaxonToPlace();
        
    const UINT*   getParsimonyBlock()  const; //returns partial_pars
    const double* getLikelihoodBlock() const; //returns partial_lh
    const UBYTE*  getScaleNumBlock()   const; //returns scale_num
    int           getTaxonId()         const; //returns taxon id
    
    /** Given a contiguous range of placements, identify the one with the best (read: lowest) score.
     @param placements - pointer to the first possible placement
     @param count - the number of possible placements
     @return the index of the placement that had the lowest score*/
    virtual size_t considerPlacements(const PossiblePlacement* placements,
                                      size_t count);
    
    /** After a placement has already been chosen (as the best - lowest score - placement),
        consider another placement (which might have a lower score).
     @param placement - the placement to consider
     @return true if the placement is the best one seen so far (note: if its score is equal
     to that of the previous placement, false will be returned).*/
    virtual bool   considerAdditionalPlacement(const PossiblePlacement& placement) ;
    
    /** return a reference to the best placement found, so far, for this taxon.
     @return best placement*/
    virtual const  PossiblePlacement& getBestPlacement() const;
    
    /** indicates whether this taxon can be inserted in its chosen (or preferred) place.
     @return true if it can*/
    virtual bool   canInsert() const;
    
    /** find the best placement, given a range of possible target branches, a search heuristic,
        and a placement cost calculator.
     @param phylo_tree the tree
     @param taxonIndex the index, of the taxon, in a TaxaToPlace container
            (this is supplied as a parameter as a SearchHeuristic might need need to know it)
     @param range the possible target branches (all will be tried)
     @param heuristic a SearchHeuristic (which has been set up, and is ready to have its
            isPlacementWorthTrying member function called)
     @param calculator the PlacementCostCalculator */
    void findPlacement ( PhyloTree& phylo_tree,
                         size_t taxonIndex,
                         TargetBranchRange& range,
                         SearchHeuristic* heuristic,
                         const PlacementCostCalculator* calculator );

    /** comparison operator (compares the scores of the best placement, for this taxon
     with the best placement, for another
     @param rhs the other taxon
     @return true if this taxon's best placement's score is less than that of the other*/
    bool operator <  ( const TaxonToPlace& rhs ) const;
    bool operator <= ( const TaxonToPlace& rhs ) const;
    
    /** insert this taxon into the tree, at its preferred location, mark the
        target branch of its preferred location as used, and add additional
        target branches (corresponding to each of the branches that are being
        added to the phylo tree by the insertion).
     @param phylo_tree the phylo tree
     @param b the BlockAllocator (needed, when adding new branches
     @param dest the target branch range (new target branches)
     @param calculator the PlacementCostCalculator (at present, only asked if uses likelihood!)*/
    void insertIntoTree ( PhyloTree& phylo_tree, BlockAllocator* b,
                          TargetBranchRange& dest,
                          PlacementCostCalculator& calculator);
    
    virtual void forgetGazumpedPlacements();
    
    /** insert this taxon into the tree, NEAR its preferred location, marking the
        target branch of the chosen location as used, and add additional
        target branches (corresponding to each of the branches that are being
        added to the phylo tree by the insertion).
     @param phylo_tree the phylo tree
     @param b the BlockAllocator (needed, when adding new branches
     @param dest the target branch range (new target branches)
     @param calculator the PlacementCostCalculator (at present, only asked if uses likelihood!)*/
    bool insertNearby ( PhyloTree& phylo_tree, BlockAllocator* b,
                        TargetBranchRange& dest,
                        PlacementCostCalculator& calculator );
    
    /** (this is a supporting member function that helps out for insertNearby)
        Search, in the target branches that have replaced a "used" target branch
        (and possibly recursively, in the replacements for those, if they in turn have
     been used), for possible placements (and score them)
     @param phylo_tree the phylo tree
     @param calculator the PlacementCostCalculator to score each possible placement
     @param tb                   the *used* target branch, whose replacement target branches
                        are to be considered
     @param placements the vector of PossiblePlacement instances,
     to which to append a PossiblePlacement for each
     unused replacement target branch.*/
    void assessNewTargetBranches ( PhyloTree& phylo_tree,
                                   PlacementCostCalculator& calculator,
                                   TargetBranch* tb,
                                   std::vector<PossiblePlacement>& placements);
};

class TaxaToPlace {
public:
    virtual ~TaxaToPlace() = default;
    
    /** fetch a taxon by index */
    virtual TaxonToPlace& getTaxonByIndex(size_t index) = 0;
    
    /** returns true if and only if there are no taxa left*/
    virtual bool isEmpty() const = 0;
    
    /** sort a batch by score
     @param batchStart the index of the first taxon in the batch
     @param batchStop one more than the index of the last taxon in the batch*/
    virtual void sortBatch(size_t batchStart, size_t batchStop) = 0;
};

template <class T=TaxonToPlace> class TypedTaxaToPlace: public TaxaToPlace, public std::vector<T> {
public:
    typedef std::vector<T> super;
    explicit TypedTaxaToPlace(size_t reservation) {
        super::reserve(reservation);
    }
    virtual bool isEmpty() const                        { return std::vector<T>::empty(); }
    virtual TaxonToPlace& getTaxonByIndex(size_t index) { return std::vector<T>::at(index); }
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
