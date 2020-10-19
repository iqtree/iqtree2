//
// placementcostcalculator.cpp
// Implementations of the PlacementCostCalculator,
// ParsimonyCostCalculator, and LikelihoodCostCalculator classes.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "placementcostcalculator.h"
#include "taxontoplace.h"
#include "placementtraversalinfo.h"

bool PlacementCostCalculator::usesLikelihood() { return false; }
PlacementCostCalculator::~PlacementCostCalculator() = default;

void ParsimonyCostCalculator::assessPlacementCost(PhyloTree& phylo_tree,
                                 const TaxonToPlace& taxon,
                                 PossiblePlacement& placement)  const {
    auto target = placement.getTarget();
    phylo_tree.computeParsimonyOutOfTree( target->getParsimonyBlock(),
                                          taxon.getParsimonyBlock(),
                                         &placement.parsimony_score);
    placement.score      = placement.parsimony_score;
    placement.lenToNode1 = target->first->findNeighbor(target->second)->length * 0.5;
    placement.lenToNode2 = placement.lenToNode1;
    if ( verbose_mode >= VB_MAX ) {
        std::stringstream s2;
        s2  << "Parsimony score for taxon " << taxon.getTaxonId()
        << " at target branch " << placement.getTargetIndex()
        << " would be " << placement.parsimony_score;
        phylo_tree.logLine(s2.str());
    }
    
    //Todo: Couldn't we get a better handle on lenToNode1 and
    //      lenToNode2. If not here, then later.
    //      After all, computeParsimonyOutOfTree could be used:
    //      pass it ( firstNei->partial_pars,  taxon->getParsimonyBlock() )
    //      and     ( secondNei->partial_pars, taxon->getParsimonyBlock() )
    //      where firstNei == first->findNeighbor(second),
    //      and use the ratio of those parsimony scores to get
    //      better estimates of lenToNode1 to lenToNode2?
    //
}

bool LikelihoodCostCalculator::usesLikelihood() {
    return true;
}
void LikelihoodCostCalculator::assessPlacementCost(PhyloTree& tree, const TaxonToPlace& taxon,
                                 PossiblePlacement& placement) const {
    
    super::assessPlacementCost(tree, taxon, placement);
    
    auto   target            = placement.getTarget();
    double score             = (placement.parsimony_score > 1) ? placement.parsimony_score : 1;
    double normalized        = score / tree.getAlnNSite();
    double alpha             = tree.getGammaShape();
    double parsimony_length  = tree.correctBranchLengthF81( normalized, alpha );
    placement.lenToNewTaxon = parsimony_length;
            
    PhyloNode fakeInterior;
    fakeInterior.is_floating_interior = true;
    
    PhyloNode fakeLeaf(taxon.new_leaf->id, taxon.new_leaf->name.c_str());
    fakeLeaf.addNeighbor(&fakeInterior, placement.lenToNewTaxon );
    
    PhyloNeighbor* neighUp = fakeLeaf.findNeighbor(&fakeInterior);
    neighUp->partial_lh    = const_cast<double*> ( target->getLikelihoodBlock());
    neighUp->scale_num     = const_cast<UBYTE*>  ( target->getScaleNumBlock() );
    neighUp->setLikelihoodComputed(true);

    double halfLen = 0.5 * target->first->findNeighbor(target->second)->length;
    fakeInterior.addNeighbor(&fakeLeaf, placement.lenToNewTaxon );
    fakeInterior.addNeighbor(target->first,  halfLen );
    fakeInterior.addNeighbor(target->second, halfLen );

    PhyloNeighbor* neighLeft    = fakeInterior.findNeighbor(target->first);
    PhyloNeighbor* realNeiLeft  = target->second->findNeighbor(target->first);
    neighLeft->partial_lh       = realNeiLeft->partial_lh;
    neighLeft->scale_num        = realNeiLeft->scale_num;
    neighLeft->setLikelihoodComputed(true);

    PhyloNeighbor* neighRight   = fakeInterior.findNeighbor(target->second);
    PhyloNeighbor* realNeiRight = target->first->findNeighbor(target->second);
    neighRight->partial_lh      = realNeiRight->partial_lh;
    neighRight->scale_num       = realNeiRight->scale_num;
    neighRight->setLikelihoodComputed(true);
    
    PhyloNeighbor* neighDown    = fakeInterior.findNeighbor(&fakeLeaf);
    
    //TREE_LOG_LINE(tree, VB_MIN, "neighLeft  is " << neighLeft);
    //TREE_LOG_LINE(tree, VB_MIN, "neighRight is " << neighRight);
    //TREE_LOG_LINE(tree, VB_MIN, "neighDown  is " << neighDown);
    
    struct LikelihoodOptimizer : public Optimization
    {
        PhyloTree&          tree;
        PhyloNode*          leaf_node;
        PhyloNeighbor*      neigh_up;
        PhyloNeighbor*      neigh_down;
        LikelihoodBufferSet buffers;
        double              initial_length;
        double              initial_score;
        
        LikelihoodOptimizer(PhyloTree&     phyloTree, PhyloNode*     leafNode,
                            PhyloNeighbor* neighUp,   PhyloNeighbor* neighDown)
            : Optimization(), tree(phyloTree), leaf_node(leafNode)
            , neigh_up(neighUp), neigh_down(neighDown)
            , buffers(phyloTree.tree_buffers) //new buffers, same sizes
            , initial_length(-1), initial_score(0) {
        }
        virtual ~LikelihoodOptimizer() = default;
        virtual double computeFunction(double value) {
            neigh_up->length   = value;
            neigh_down->length = value;
            double score = -tree.computeLikelihoodBranch(neigh_up, leaf_node, buffers);
            TREE_LOG_LINE(tree, VB_MAX, "  length " << value
                          << ", f " << score
                          << " for taxon " << leaf_node->id);
            return score;
        }
        virtual void computeFuncDerv(double value, double &df, double &ddf) {
            neigh_up->length   = value;
            neigh_down->length = value;
            tree.computeLikelihoodDerv(neigh_up, leaf_node, &df, &ddf, buffers);
            df  = -df;
            ddf = -ddf;
            TREE_LOG_LINE(tree, VB_MAX, "  length " << value
                          << ", df " << df << ", ddf " << ddf
                          << " for taxon " << leaf_node->id);
        }
        void optimize(PossiblePlacement& p) {
            PlacementTraversalInfo info(tree, buffers, neigh_up, leaf_node );
            info.computePartialLikelihood();
            
            initial_length  = neigh_up->length;
            initial_score   = computeFunction(initial_length);
                            
            if (tree.optimize_by_newton) {
                double optx, derivative_of_likelihood_wrt_length=0;
                // Newton-Raphson method (note: 10 = max # steps hardcoded for now
                optx = minimizeNewton ( tree.params->min_branch_length, initial_length,
                                        tree.params->max_branch_length,
                                        tree.params->min_branch_length,
                                        derivative_of_likelihood_wrt_length, 10);
                p.lenToNewTaxon = optx;
                p.score         = computeFunction(optx);
                if (optx > tree.params->max_branch_length*0.95 && !tree.isSuperTree()) {
                    // newton raphson diverged, reset
                    if (initial_score < p.score) {
                        p.lenToNewTaxon = initial_length;
                        p.score         = initial_score;
                    }
                }
            } else {
                // Brent method
                double optx, negative_lh = 0, ferror;
                optx = minimizeOneDimen(tree.params->min_branch_length,
                                        initial_length,
                                        tree.params->max_branch_length,
                                        tree.params->min_branch_length,
                                        &negative_lh, &ferror);
                p.lenToNewTaxon = optx;
                p.score         = negative_lh;
            }
        }
        double getInitialScore() {
            return initial_score;
        }
    };
    LikelihoodOptimizer opt(tree, &fakeLeaf, neighUp, neighDown);
    opt.optimize(placement);
    TREE_LOG_LINE(tree, VB_DEBUG, taxon.taxonName
                  << " at target " << placement.getTargetIndex()
                  << " had parsimony score " << placement.parsimony_score
            << ", len "    << parsimony_length
            << ", old lh " << opt.getInitialScore()
                  << ", new lh " << placement.score
                  << ", len "    << placement.lenToNewTaxon);
}

PlacementCostCalculator* PlacementCostCalculator::getNewCostCalculator(Placement::CostFunction fun) {
    switch (fun) {
        case Placement::MAXIMUM_PARSIMONY: /* fall-through */
        case Placement::SANKOFF_PARSIMONY:
            return new ParsimonyCostCalculator();
        case Placement::MAXIMUM_LIKELIHOOD_ANYWHERE:
        case Placement::MAXIMUM_LIKELIHOOD_MIDPOINT:
            return new LikelihoodCostCalculator();
        default:
            outError("No cost calculator class available for specified cost function");
            return nullptr;
    }
}
