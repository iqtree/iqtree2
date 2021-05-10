//
//parsimonynni.h
//Declarations of ParsimonyNNIMove (which describes a possible NNI move)
//Created by James Barbetti on 05-Mar-2021 (declarations moved here).
//

#ifndef parsimonynni_h
#define parsimonynni_h

#include "parsimonymove.h"

class ParsimonyNNIMove : public ParsimonyMove {
public:
    typedef ParsimonyMove super;
    PhyloNode*  left;
    PhyloBranch middle;
    PhyloNode*  right;
    ParsimonyNNIMove();
    ParsimonyNNIMove(const ParsimonyNNIMove& rhs);
    ParsimonyNNIMove& operator=(const ParsimonyNNIMove& rhs);
    virtual ~ParsimonyNNIMove() = default;
    virtual void initialize(intptr_t source_branch, bool beLazy);
    virtual std::string getDescription() const;
    static intptr_t getParsimonyVectorSize(intptr_t radius);
    static intptr_t getMinimumPathVectorCount();
    virtual void   findMove(const PhyloTree& tree, 
                            const TargetBranchRange& branches,
                            int /*radius*/ /* ignored; just part of signature */,
                            std::vector<UINT*> &path_parsimony,
                            double parsimony_score);
    virtual void finalize(PhyloTree& tree,
                          const TargetBranchRange& branches);
    virtual bool isStillPossible(const TargetBranchRange& branches,
                                 PhyloBranchVector& path) const;
    virtual double recalculateBenefit
                   ( PhyloTree& tree, double tree_parsimony_score,
                     TargetBranchRange& branches,
                     LikelihoodBlockPairs &blocks,
                     ParsimonyPathVector& parsimony_path_vectors) const;
    virtual double apply(PhyloTree& tree,
                         double parsimony_score,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks,
                         ParsimonyPathVector& parsimony_path_vectors);
protected:
    void consider(intptr_t branch_id, PhyloNode* leftNode,
                  const PhyloBranch& middleBranch,
                  PhyloNode* rightNode, double move_benefit);
    PhyloNode* getOtherLeftNode() const;
    PhyloNode* getOtherRightNode() const;
};


#endif /* parsimonynni_h */
