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
#include <utils/hammingdistance.h>  //for hammingDistance() template function

PlacementCostCalculator::~PlacementCostCalculator() = default;

bool PlacementCostCalculator::usesHammingDistance()  const { return false; }

bool PlacementCostCalculator::usesLikelihood()       const { return false; }

bool PlacementCostCalculator::usesParsimony()        const { return false; }

bool PlacementCostCalculator::usesSankoffParsimony() const { return false; }

void PlacementCostCalculator::assessPlacementCost(PhyloTree& phylo_tree,
                                                  const TaxonToPlace& taxon,
                                                  PossiblePlacement& placement)  const {
    auto target = placement.getTarget();
    placement.parsimony_score = -1;
    placement.score = -1;
    placement.lenToNode1 = target->first->findNeighbor(target->second)->length * 0.5;
    placement.lenToNode2 = placement.lenToNode1;
}

HammingDistanceCostCalculator::HammingDistanceCostCalculator() = default;

bool HammingDistanceCostCalculator::usesHammingDistance() const {
    return true;
}

void HammingDistanceCostCalculator::assessPlacementCost
        (PhyloTree &phylo_tree, const TaxonToPlace &taxon, 
         PossiblePlacement &placement) const {
    auto target = placement.getTarget();

    char     unk         = static_cast<char>(phylo_tree.aln->STATE_UNKNOWN);
    auto     freq_vector = phylo_tree.getConvertedSequenceFrequencies();
    intptr_t seq_len     = phylo_tree.getConvertedSequenceLength();
    double   denominator = (double)phylo_tree.aln->getNSite();
    const char* that_seq = nullptr;

    if (target->first->isLeaf()) {
        that_seq = phylo_tree.getConvertedSequenceByNumber(target->first->id);
    } else if (target->second->isLeaf()) {
        that_seq = phylo_tree.getConvertedSequenceByNumber(target->first->id);
    }

    if (that_seq==nullptr) {
        placement.score = 10.0; //Don't want it here!
    }
    else {
        const char* this_seq = phylo_tree.getConvertedSequenceByNumber(taxon.getTaxonId());
        double unknown_freq  = 0;
        double hamming_dist  = hammingDistance(unk, this_seq, that_seq,
                                            seq_len, freq_vector, unknown_freq);
        if (unknown_freq < denominator) {
            placement.score = hamming_dist / (denominator - unknown_freq);
        } else {
            placement.score = 0;
        }
    }
    placement.lenToNewTaxon = -1; //Let Parsimony &/or Likelihood decide
    placement.lenToNode1    = -1; //Ditto
    placement.lenToNode2    = -1; //Double Ditto
}

ParsimonyCostCalculator::ParsimonyCostCalculator(bool useSankoff)
: sankoff(useSankoff) {
}

bool ParsimonyCostCalculator::usesParsimony()        const { return true; }

bool ParsimonyCostCalculator::usesSankoffParsimony() const { return sankoff; }

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
    if ( verbose_mode >= VerboseMode::VB_MAX ) {
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


LikelihoodCostCalculator::LikelihoodCostCalculator(bool useMidpoint)
    : ParsimonyCostCalculator(false), midpoint(useMidpoint) {
}

LikelihoodCostCalculator::~LikelihoodCostCalculator()
{
}

bool LikelihoodCostCalculator::usesLikelihood() const {
    return true;
}

void LikelihoodCostCalculator::assessPlacementCost(PhyloTree& tree, const TaxonToPlace& taxon,
                                                   PossiblePlacement& placement) const {
    super::assessPlacementCost(tree, taxon, placement);
    
    const TargetBranch* target = placement.getTarget();
    double score               = placement.parsimony_score;
    double adjustment          = 0;
    if (tree.params->add_uninformative_sites_to_parsimony_length) {
        if (taxon.taxonId < tree.aln->singleton_parsimony_states.size()) {
            adjustment = static_cast<double>(tree.aln->singleton_parsimony_states[taxon.taxonId]);
        }
    }

    double parsimony_length    = (score+adjustment) / tree.getAlnNSite();
    tree.correctBranchLengthIfNeedBe(parsimony_length);
    
    placement.lenToNewTaxon    = parsimony_length;
            
    PhyloNode fakeInterior;
    fakeInterior.is_floating_interior = true;
    
    PhyloNode fakeLeaf(taxon.new_leaf->id, taxon.new_leaf->name.c_str());
    fakeLeaf.addNeighbor(&fakeInterior, placement.lenToNewTaxon );
    
    PhyloNeighbor* neighUp   = fakeLeaf.findNeighbor(&fakeInterior);
    neighUp->partial_lh      = target->getLikelihoodBlock();
    neighUp->scale_num       = target->getScaleNumBlock();
    neighUp->lh_scale_factor = target->getLhScaleFactor();
    neighUp->setLikelihoodComputed(true);

    fakeInterior.addNeighbor(&fakeLeaf, placement.lenToNewTaxon );
    PhyloNeighbor* neighDown = fakeInterior.findNeighbor(&fakeLeaf);
    
    if (!this->midpoint) {
        double halfLen = 0.5 * target->first->findNeighbor(target->second)->length;
        fakeInterior.addNeighbor(target->first,  halfLen );
        fakeInterior.addNeighbor(target->second, halfLen );

        PhyloNeighbor* neighLeft    = fakeInterior.findNeighbor(target->first);
        PhyloNeighbor* realNeiLeft  = target->second->findNeighbor(target->first);
        neighLeft->partial_lh       = realNeiLeft->partial_lh;
        neighLeft->scale_num        = realNeiLeft->scale_num;
        neighLeft->lh_scale_factor  = realNeiLeft->lh_scale_factor;
        neighLeft->setLikelihoodComputed(true);

        PhyloNeighbor* neighRight   = fakeInterior.findNeighbor(target->second);
        PhyloNeighbor* realNeiRight = target->first->findNeighbor(target->second);
        neighRight->partial_lh      = realNeiRight->partial_lh;
        neighRight->scale_num       = realNeiRight->scale_num;
        neighRight->lh_scale_factor = realNeiRight->lh_scale_factor;
        neighRight->setLikelihoodComputed(true);
    }    
    
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
#if (0)
            TREE_LOG_LINE(tree, VerboseMode::VB_MIN, "neigh_up   is " << pointer_to_hex(neighUp));
            TREE_LOG_LINE(tree, VerboseMode::VB_MIN, "neigh_down is " << pointer_to_hex(neighDown));
            std::stringstream s;
            const char* sep = "nu ";
            auto lh_block = tree.getLhBlockSize();
            for (int i=0; i<lh_block; ++i) {
                s << (((i&7)==0) ? sep : " ");
                s << neighUp->partial_lh[i];
                sep = "\nnu ";
            }
            std::cout << s.str() << std::endl;
#endif
        }
        virtual ~LikelihoodOptimizer() = default;
        virtual double computeFunction(double value) {
            neigh_up->length   = value;
            neigh_down->length = value;
            double score = -tree.computeLikelihoodBranch(neigh_up, leaf_node, buffers);
            TREE_LOG_LINE(tree, VerboseMode::VB_MAX, "  length " << value
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
            TREE_LOG_LINE(tree, VerboseMode::VB_MAX, "  length " << value
                          << ", df " << df << ", ddf " << ddf
                          << " for taxon " << leaf_node->id);
        }
        void optimize(PossiblePlacement& p) {
            #if (0)
                PlacementTraversalInfo info(tree, buffers, neigh_up, leaf_node );
                info.computePartialLikelihood();
            #endif
            
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
    const_cast<TargetBranch*>(target)->setLhScaleFactor(neighUp->lh_scale_factor);
    
    TREE_LOG_LINE(tree, VerboseMode::VB_MED, taxon.taxonName
                  << " (taxon " << taxon.taxonId << ")"
                  << " at target branch " << placement.getTargetIndex()
                  << " had parsimony score " << placement.parsimony_score
                  << ", len "    << parsimony_length
                  << ", old lh " << opt.getInitialScore()
                  << ", new lh " << placement.score
                  << ", len "    << placement.lenToNewTaxon);
}

PlacementCostCalculator* PlacementCostCalculator::getNewCostCalculator
                         (const PlacementParameters& params) {
    
    auto cf = params.getIncrementalParameter('C', DEFAULT_COST_FUNCTION);
    if (cf=="ML") {
        return new LikelihoodCostCalculator(true);
    } else if (cf=="FML") {
        return new LikelihoodCostCalculator(false);
    } else if (cf=="MP") {
        return new ParsimonyCostCalculator(false);
    } else if (cf=="SMP") {
        return new ParsimonyCostCalculator(true);
    } else if (cf=="HD") {
        return new HammingDistanceCostCalculator();
    }
    outError("No cost calculator class available for specified cost function");
    return nullptr;
}
