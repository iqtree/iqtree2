//
// placementcalculator.h
// Declarations of classes for calculating placement costs:
// 1. PlacementCostCalculator  (the superclass)
// 2. ParsimonyCostCalculator  (for costing placements by parsimony)
// 3. LikelihoodCostCalculator (for costing placements by (un)likelihood)
//
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef placementcostcalculator_hpp
#define placementcostcalculator_hpp

#include <tree/phylotree.h>
#include "possibleplacement.h"
#include "placement.h" //for Placement::doesPlacementUseParsimony, doesPlacementUseLikelihood

class TaxonToPlace;

class PlacementCostCalculator {
public:
    virtual ~PlacementCostCalculator();

    /** Determine what the cost is, of placing a taxon on a particular target branch.
        (lower cot is better).
     @param tree the PhyloTree
     @param taxon The taxon to be placed
     @param p A place it might go*/
    virtual void assessPlacementCost(PhyloTree& tree, const TaxonToPlace& taxon,
                                     PossiblePlacement& p) const = 0;

    /** Indicate whether this placement cost calculator does likelihood calculations
     (and needs PhyloNeighbor instances, for example, to have partial likelihood vectors allocated)
     @returns true if it does likelihood calculations, false if not
     */
    virtual bool usesLikelihood();

    /** Indicate whether this placement cost calculator uses sankoff parsimony
     (and needs PhyloNeighbor instances, for example, to have partial parsimony vectors allocated)
     @returns true if it does SANKOFF parsimony calculations, false if not
     (so true if uses regular parsimony)
     */
    virtual bool usesSankoffParsimony();

    /**
     @param  the placement const function 
     @return the placement cost calculator to use, allocated via new
             (it is up to the caller to delete it)
     */
    static PlacementCostCalculator* getNewCostCalculator();
};

class ParsimonyCostCalculator : public PlacementCostCalculator {
private:
    bool sankoff;
public:
    ParsimonyCostCalculator(bool usesSankoff);
    virtual bool usesSankoffParsimony();
    virtual void assessPlacementCost(PhyloTree& phylo_tree,
                                     const TaxonToPlace& taxon,
                                     PossiblePlacement& placement) const;
};

class LikelihoodCostCalculator : public ParsimonyCostCalculator {
private:
    bool midpoint;
public:
    typedef ParsimonyCostCalculator super;
    LikelihoodCostCalculator(bool useMidpoint);
    virtual bool usesLikelihood();
    virtual void assessPlacementCost(PhyloTree& tree, const TaxonToPlace& taxon,
                                     PossiblePlacement& placement) const;
};


#endif /* placementcalculator_h */
