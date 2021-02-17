//
// parallelparsimonycalculator.h
// Created by James Barbetti on 08-Oct-2020.
//

#ifndef parallelparsimonycalculator_h
#define parallelparsimonycalculator_h

#include <tree/phylotree.h>
#include <tree/phylonode.h>

class ParallelParsimonyCalculator {
private:
    PhyloTree& tree; //The PhyloTree
    typedef std::pair<PhyloNeighbor*, PhyloNode*> WorkItem;
    std::vector <WorkItem> workToDo;  //A *stack* (rather than a queue)
                                      //of work to do ( partial parsimonies to calculate).
    std::vector <int> workDone;       //A *stack* of boundaries
    
    const char* task_to_start;
    const char* task_in_progress;
    bool        report_progress_to_tree;
public:
    explicit ParallelParsimonyCalculator(PhyloTree& phylo_tree, bool report_back=false);

    /**
     Indicate a PhyloNeighbor whose partial parsimony is to be calculated
     (but don't request its calculation yet - see the calculate() method).
     @param dad_branch the branch for which parsimony is to be calculated
     @param dad the PhyloNode for which dad_branch is a neighor
     @return 1 if it needs to be recalculated, 0 if not
     */
    int schedulePartialParsimony(PhyloNeighbor* dad_branch, PhyloNode* dad) ;

    /**
     Indicate a PhyloNeighbor whose partial parsimony is to be calculated
     (but don't request its calculation yet - see the calculate() method).
     @param branch the branch (parsimony is to be calculated for
         the view from branch.first to branch.second, but not the converse).
     @return 1 if it needs to be recalculated, 0 if not
     */
    int schedulePartialParsimony(const PhyloBranch& branch);

    /**
     Calculate the partial parsimony for a branch (by calling computePartialParsimony for each
     Neighbor (the branch seen "in both orientations", and then calling calculate()).
     @param dad_branch the branch for which parsimony is to be calculated
     @param dad the PhyloNode for which dad_branch is a neighor
     @param taskDescription - description to display in a progress bar
        (if progress bars are being displayed)
     */
    int  computeParsimonyBranch(PhyloNeighbor* dad_branch, PhyloNode* dad,
                                const char* taskDescription="");

    /**
     Calculate reverse partial parsimony, back towards a specify branch, for every
     other node in the tree.
     @param first    one of the nodes in the branch
     @param second  the other nodes in the branch
     */
    void computeReverseParsimony(PhyloNode* first, PhyloNode* second);

    /**
     Calculate partial parsimony in both directions, along every branch in the tree
     @param dad_branch    the branch for which parsimony is to be calculated
     @param dad the PhyloNode for which dad_branch is a neighor
     @param taskDescription - description to display in a progress bar
        (if progress bars are being displayed)
     */
    int  computeAllParsimony(PhyloNeighbor* dad_branch, PhyloNode* dad,
                             const char* taskDescription="");

    /**
     Calculate the partial parsimonies of the PhyloNeighbor instances that have
     been passed to computePartialParsimony(), since the last time calculate() was called
     (if it ever was).
     @param start_index the "top" of the stack (the point at which to begin calculation)
     @param taskDescription a description of the task; may point to an empty string
     */
    void calculate(intptr_t start_index = 0, const char* taskDescription="");
    
    /**
     * Calculate the cost of linking the nodes a, through f, in the configuration shown at
     *           left.  On entry, it is assumed that a, b, e, and f are connected
     *  a      b          somehow to c and d. Any of the following configurations is valid
     *   \   /            (a,e) connected to c, (a,f) connected to c, (a,b) connected to c.
     *    c             (b,e) connected to c, (b,f) connected to c, (e,f) connected to c.
     *    |
     *    d            The subtree costs (viewing out from c and d) are calculated,
     *   /    \          from the existing tree structure, then the cost of connection the a, and b
     *  e       f        subtrees (to c), the cost of connecting the e and f (to d), and lastly the
     *          cost of connecting c and d.
     *
     * @param a a node
     * @param b a node (its subtree is to be connected, first, to the subtree rooted on a)
     * @param c an interior node (to which any two of a,b, e, f are connected)
     * @param d an interior node (to which the other two of a,b, e, f are connected)
     * @param e a node
     * @param f a node (its subtree is to be connected, first, to the subtree rooted on e
     * @param buffer1 sufficient storage for a parsimony vector for the "ab" subtree,
     *                viewed from d toward c.
     * @param buffer2 sufficient storage for a parsimony vector for the "ef" subtree
     *                viewed from c toward a.
     */
    double parsimonyLink4Cost(PhyloNode* a, PhyloNode* b, PhyloNode* c,
                              PhyloNode* d, PhyloNode* e, PhyloNode* f,
                              UINT* buffer1, UINT* buffer2);
    
    static double parsimonyLink4CostOutOfTree
                  ( const PhyloTree& tree, PhyloNode* a, PhyloNode* b, PhyloNode* c,
                    PhyloNode* d, PhyloNode* e, PhyloNode* f,
                    UINT* buffer1, UINT* buffer2);
    
};


#endif /* parallelparsimonycalculator_h */
