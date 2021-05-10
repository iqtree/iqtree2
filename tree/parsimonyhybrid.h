//
//parsimonyhybrid.h
//Describes a possible move, that is the union of two different 
//types of possible move (e.g. an NNI *or* an SPR move).
//
//Created by James Barbetti on 05-Mar-2021.
//

#ifndef parsimonyhybrid_h
#define parsimonyhybrid_h

#include "parsimonymove.h"
#include "parsimonynni.h"
#include "parsimonyspr.h"
#include "parsimonytbr.h"

template <class A=ParsimonySPRMove,
          class B=ParsimonyNNIMove,
          bool halfRadiusForB=false>
class ParsimonyHybridMove: public ParsimonyMove {
protected:
    A    alpha_move;
    B    beta_move;
    mutable bool alpha_move_was_better;
    void copyFrom(const ParsimonyMove& original) {
        benefit = original.benefit;
        source_branch_id = original.source_branch_id;
    }
    
public:    
    static intptr_t getParsimonyVectorSize(intptr_t radius) {
        intptr_t alpha_size = A::getParsimonyVectorSize(radius);
        intptr_t beta_size  = B::getParsimonyVectorSize(radius);
        return (alpha_size < beta_size) ? beta_size : alpha_size;
    }

    static intptr_t getMinimumPathVectorCount() {
        intptr_t alpha_count = A::getMinimumPathVectorCount();
        intptr_t beta_count  = B::getMinimumPathVectorCount();
        return (alpha_count < beta_count) ? beta_count : alpha_count;
    }

    virtual void   initialize(intptr_t source_branch_id_to_use,
                              bool be_lazy) {
        benefit          = 0;
        source_branch_id = source_branch_id_to_use;
        lazy             = be_lazy;
        alpha_move.initialize(source_branch_id, be_lazy);
        beta_move.initialize(source_branch_id, be_lazy);
    }
    
    virtual void   findMove(const PhyloTree& tree,
                            const TargetBranchRange& branches,
                            int radius,
                            std::vector<UINT*> &path_parsimony,
                            double parsimony_score) {
        alpha_move.findMove(tree, branches, radius,
                            path_parsimony, parsimony_score);
        beta_move.findMove(tree, branches,
                           halfRadiusForB ? (radius/2) : radius,
                           path_parsimony, parsimony_score);
    }
    
    virtual void   finalize(PhyloTree& tree,
                            const TargetBranchRange& branches) {
        alpha_move.finalize(tree, branches);
        beta_move.finalize(tree, branches);
        alpha_move_was_better = beta_move <= alpha_move;
        if (alpha_move_was_better) {
            copyFrom(alpha_move);
        }
        else {
            copyFrom(beta_move);
        }
        positions_considered = alpha_move.positions_considered
                             + beta_move.positions_considered;
    }
    
    virtual std::string getDescription() const {
        return alpha_move_was_better
            ? alpha_move.getDescription()
            : beta_move.getDescription();
    }
    
    virtual bool   isStillPossible(const TargetBranchRange& branches,
                                   PhyloBranchVector& path) const {
        if (alpha_move_was_better) {
            if (alpha_move.isStillPossible(branches, path)) {
                return true;
            }
            alpha_move_was_better = false;
            const_cast<ParsimonyHybridMove*>(this)->copyFrom(beta_move);
            return 0 < beta_move.benefit && beta_move.isStillPossible(branches, path);
        } else {
            if (beta_move.isStillPossible(branches, path)) {
                return true;
            }
            alpha_move_was_better = true;
            const_cast<ParsimonyHybridMove*>(this)->copyFrom(alpha_move);
            return 0 < alpha_move.benefit && alpha_move.isStillPossible(branches, path);
        }
    }
    
    virtual double recalculateBenefit(PhyloTree& tree, double tree_parsimony_score,
                                      TargetBranchRange& branches,
                                      LikelihoodBlockPairs &blocks,
                                      ParsimonyPathVector& parsimony_path_vectors) const {
        return alpha_move_was_better
            ? alpha_move.recalculateBenefit(tree, tree_parsimony_score,
                                            branches, blocks,
                                            parsimony_path_vectors)
            : beta_move.recalculateBenefit(tree, tree_parsimony_score,
                                           branches, blocks,
                                           parsimony_path_vectors);
    }
    
    virtual double apply(PhyloTree& tree,
                         double parsimony_score,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks,
                         ParsimonyPathVector& parsimony_path_vectors) {
        return (alpha_move_was_better)
            ? alpha_move.apply(tree, parsimony_score, branches, blocks, parsimony_path_vectors)
            : beta_move.apply (tree, parsimony_score, branches, blocks, parsimony_path_vectors);
    }
};

#endif /* parsimonyhybrid_h */
