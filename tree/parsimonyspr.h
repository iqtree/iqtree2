//
//parsimonyspr.h
//Declarations of ParsimonyLazySPRMove and ParsimonySPRMove.
//Created by James Barbetti on 05-Mar-2021 (declarations moved here).
//

#ifndef parsimonyspr_h
#define parsimonyspr_h

#include "parsimonymove.h"

struct ParsimonyLazySPRMove : public ParsimonyMove {
    public:
        typedef  ParsimonyLazySPRMove this_type;
        typedef  ParsimonyMove        super;

        intptr_t target_branch_id;
        bool     isForward;        //branch.first moves, branch.second does not
        
        //Members used for checking whether an SPR is still valid
        //after other changes might have been made to the tree.
    private:
        PhyloNode* source_first;
        PhyloNode* source_second;
        PhyloNode* target_first;
        PhyloNode* target_second;
        
    public:
        ParsimonyLazySPRMove(const this_type& rhs) = default;
        ParsimonyLazySPRMove& operator=(const this_type& rhs) = default;
        
        struct LazySPRSearch {
            public:
            const PhyloTree&         tree;
            const TargetBranchRange& branches;
            const TargetBranch&      source;
            double                   discon;
            ParsimonyLazySPRMove&    put_answer_here;
            
            LazySPRSearch(const PhyloTree& phylo_tree,
                          const TargetBranchRange& target_branches,
                          double disconnection_benefit,
                          ParsimonyLazySPRMove& output);
            
            void searchForForwardsSPR(PhyloNode* current, PhyloNode* prev,
                                      int radius);
            void searchForBackwardsSPR(PhyloNode* current, PhyloNode* prev,
                                       int radius);
        };
public:
    ParsimonyLazySPRMove();
    static intptr_t getParsimonyVectorSize(intptr_t radius);
    virtual void initialize(intptr_t source_branch, bool beLazy);
    virtual std::string getDescription() const;
    virtual void finalize(PhyloTree& tree,
                          const TargetBranchRange& branches);
    void findForwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                            int radius, double disconnection_benefit);
    void findBackwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                             int radius, double disconnection_benefit);
    virtual bool isStillPossible(const TargetBranchRange& branches,
                                 PhyloBranchVector& path) const;
    virtual double recalculateBenefit
                   ( PhyloTree& tree, double parsimony_score,
                     TargetBranchRange& branches,
                     LikelihoodBlockPairs &blocks,
                     ParsimonyPathVector& parsimony_path_vectors) const;
    virtual double apply(PhyloTree& tree,
                         double parsimony_score,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks,
                         ParsimonyPathVector& parsimony_path_vectors);
}; //ParsimonyLazySPRMove

class ParsimonySPRMove: public ParsimonyLazySPRMove {
public:
    typedef  ParsimonySPRMove this_type;
    typedef  ParsimonyLazySPRMove super;
    ParsimonySPRMove(const this_type& rhs) = default;
    ParsimonySPRMove& operator=(const this_type& rhs) = default;
    ParsimonySPRMove() = default;

    static intptr_t getParsimonyVectorSize(intptr_t radius);

    struct ProperSPRSearch: public LazySPRSearch {
        public:
        typedef LazySPRSearch super;
        std::vector<UINT*>& path_parsimony;
            //This is a vector of size at least radius + 2, for any
            //radius passed to searchForForwardsSPR or to
            //searchForBackwardsSPR.  The [tree->params->spr_radius+1]
            //entry in the vector can be a copy (when disconnecting
            //x from left and right), it is either the parsimony
            //vector of the view, from x, towards the subtree
            //containing left, OR the parsimony vector of the view
            //from x, towards the subtree containing right.
            //The other radius + 1 entries are all "private" to the
            //search and will be calculated (and repeatedly overwritten)
            //during the search.
        ProperSPRSearch(const PhyloTree& phylo_tree,
                        const TargetBranchRange& target_branches,
                        double disconnection_benefit,
                        std::vector<UINT*>& path_parsimony_to_use,
                        ParsimonyLazySPRMove& output);
        PhyloNode* other_adjacent_node(PhyloNode*a, PhyloNode*b, PhyloNode*c);
        void prepareToSearch(PhyloNode* left, PhyloNode* right,
                             PhyloNode* snipped, int radius) ;
        void searchForForwardsSPR(PhyloNode* current, PhyloNode* prev,
                                  int radius, double parsimony);
        void searchForBackwardsSPR(PhyloNode* current, PhyloNode* prev,
                                   int radius, double parsimony);
    };
    void findForwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                        int radius, double disconnection_benefit,
                        std::vector<UINT*> &path_parsimony,
                        double parsimony_score);
    void findBackwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                         int radius, double disconnection_benefit,
                         std::vector<UINT*> &path_parsimony,
                         double parsimony_score);
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius,
                          std::vector<UINT*> &path_parsimony,
                          double parsimony_score);
}; //ParsimonySPRMove

#endif /* parsimonyspr_h */
