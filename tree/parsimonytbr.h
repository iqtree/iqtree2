//
//parsimonytbr.h
//Declarations of ParsimonyLazyTBRMove and ProperParsimonyTBRMove.
//Created by James Barbetti on 05-Mar-2021 (declarations moved here).
//

#ifndef parsimonytbr_h
#define parsimonytbr_h

#include "parsimonymove.h"

struct ParsimonyLazyTBRMove : public ParsimonyMove {
public:
    typedef  ParsimonyLazyTBRMove this_type;
    typedef  ParsimonyMove        super;

    double   disconnection_benefit; //benefit of snipping the source branch
    intptr_t first_target_branch_id;
    intptr_t second_target_branch_id;
    intptr_t better_positions;
    
    PhyloBranch copy_of_source;
    PhyloBranch copy_of_first_target;
    PhyloBranch copy_of_second_target;
    
    int depth; //search depth is used for tie breaks
               //when two TBR moves have equal benefit. Short range
               //tbr moves mess up less of the tree, and are to
               //be preferred over long-range moves, for that reason.

    ParsimonyLazyTBRMove(const this_type& rhs)            = default;
    ParsimonyLazyTBRMove& operator=(const this_type& rhs) = default;
    ParsimonyLazyTBRMove();
    
    static intptr_t getParsimonyVectorSize(intptr_t radius);
    
    virtual void initialize(intptr_t id_of_source_branch, bool beLazy);
    
    virtual std::string getDescription() const;
    
    virtual void finalize(PhyloTree& tree,
                          const TargetBranchRange& branches) ;
    
    bool doBranchesTouch(const TargetBranchRange& branches,
                         intptr_t id_1, intptr_t id_2) const;
    
    virtual bool isStillPossible(const TargetBranchRange& branches,
                                 PhyloBranchVector& path) const;
    virtual double recalculateBenefit
                   ( PhyloTree& tree, double parsimony_score,
                     TargetBranchRange& branches,
                     LikelihoodBlockPairs &blocks,
                     ParsimonyPathVector& parsimony_path_vectors) const;
    
    class LazyTBRSearch {
    public:
        const PhyloTree&         tree;
        const TargetBranchRange& branches;
        int                      max_radius;
        std::vector<UINT*>&      path_parsimony;
        double                   parsimony_score; //of the tree as it is *before* any TBR move
        ParsimonyLazyTBRMove&    move;
        PhyloNode*               front;
        PhyloNode*               back;
        intptr_t                 first_id;
        int                      first_depth;
        
        LazyTBRSearch(const PhyloTree& t, const TargetBranchRange& b, int r,
                      std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m );
         /**
          * @param from    where we are
          * @param prev    where we were
          * @param depth  how far we have gone (0 if we are on a branch adjacent
          *              to the source branch).
          * @note  this should not be declared virtual, as it has the same name
          *        as a member function of one of its subclasses and we want
          *        LazyTBRSearch::findMove to see *this* function, not the version
          *        in the subclass.
          */
        void searchPart1(PhyloNode* from, PhyloNode* prev, int depth);
        /**
         * @param from    where we are
         * @param prev    where we were
         * @param depth  how far we have searched on both sides of the
         *              source branch (==first_depth+1, if we are looking at
         *              a branch adjacent to the "back" end of the source branch).
         * @note  this should not be declared virtual, as it has the same name
         *        as a member function of one of its subclasses and we want
         *        LazyTBRSearch::findMove to see *this* function, not the version
         *        in the subclass.
         */
        void searchPart2(PhyloNode* from, PhyloNode* prev,
                         int depth);
        void considerMove(intptr_t first_id, intptr_t second_id,
                          double gain, int depth);
    };
    
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius,
                          std::vector<UINT*>& path_parsimony,
                          double parsimony_score);
    
    void getOtherNeighbors(PhyloNode* of, PhyloNode* but_not,
                           PhyloNode** put_here, intptr_t* branch_ids);
    
    void getBranchNodes(const TargetBranch& b, PhyloNode** put_here);
    
    void disconnect(PhyloNode* first, PhyloNode* second, PhyloNode* third);
    
    void reconnect(PhyloNode* first, PhyloNode* second,
                   PhyloNode* third, PhyloNode* fourth);
    
    void updateBranch(TargetBranchRange& branches, intptr_t id,
                      PhyloNode* left, PhyloNode* right);
    
    virtual double apply(PhyloTree& tree,
                         double parsimony_score,
                         TargetBranchRange& branches,
                         LikelihoodBlockPairs blocks,
                         ParsimonyPathVector& parsimony_path_vectors);
};

struct ProperParsimonyTBRMove : public ParsimonyLazyTBRMove {
public:
    typedef ParsimonyLazyTBRMove super;
    static intptr_t getParsimonyVectorSize(intptr_t radius) ;
    class ProperTBRSearch: public LazyTBRSearch {
    public:
        typedef LazyTBRSearch super;
        std::vector<PhyloNode*> path;
        UINT*  front_parsimony;
        double front_score;
        UINT*  back_parsimony;
        double back_score;
        ProperTBRSearch(const PhyloTree& t, const TargetBranchRange& b, int r,
                        std::vector<UINT*>& p, double s, ParsimonyLazyTBRMove& m );
        inline UINT* offPathParsimony(PhyloNode* a, PhyloNode* b, PhyloNode* c ) ;
        
        /**
         * @param from    where we are
         * @param prev    where we were
         * @param depth  how far we have gone (0 if we are on a branch adjacent
         *              to the source branch).
         * @note  path_parsimony[max_radius-1] is the calculated vector
         *        for the view of the subtree from the first target branch.
         *        path[0..depth-1] are the nodes visited in the path to "prev".
         */
        void searchPart1(PhyloNode* from, PhyloNode* prev, int depth);
        /**
         * @param from    where we are (on the back side of the source branch)
         * @param prev    where we were
         * @param depth  how far we have searched on both sides of the
         *              source branch (==first_depth+1, if we are looking at
         *              a branch adjacent to the "back" end of the source branch).
         * @note  path_parsimony[max_radius] is the calculated vector
         *        for the view of the subtree from the second target branch.
         *        path[0..first_depth]] are the nodes visited in the path from
         *        front to the first target branch.
         *        path[first_depth+1..depth-1] are the nodes visited in the path
         *        from back to where we are search now.
         */
        void searchPart2(PhyloNode* from, PhyloNode* prev, int depth) ;
    };
    
    virtual void findMove(const PhyloTree& tree,
                          const TargetBranchRange& branches,
                          int radius,
                          std::vector<UINT*>& path_parsimony,
                          double parsimony_score);
};

#endif /* parsimonytbr_h */
