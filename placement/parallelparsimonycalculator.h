//
// parallelparsimonycalculator.h
//
//
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
    const char* task_to_start;
    const char* task_in_progress;
public:
    ParallelParsimonyCalculator(PhyloTree& phylo_tree);

    /**
     Indicate a PhyloNeighbor whose partial parsimony is to be calculated
     (but don't requiest its calculation yet - see the calculate() method).
     @param dad_branch the branch for which parsimony is to be calculated
     @param dad the PhyloNode for which dad_branch is a neighor
     */
    void computePartialParsimony(PhyloNeighbor* dad_branch, PhyloNode* dad) ;

    /**
     Calculate the partial parsimony for a branch (by calling computePartialParsimony for each
     Neighbor (the branch seen "in both orientations", and then calling calculate()).
     @param dad_branch the branch for which parsimony is to be calculated
     @param dad the PhyloNode for which dad_branch is a neighor
     @param taskDescription - description to display in a progress bar
        (if progress bars are being displayed)
     */
    int  computeParsimonyBranch(PhyloNeighbor* dad_branch, PhyloNode* dad, const char* taskDescription="");

    /**
     Calculate the partial parsimonies of the PhyloNeighbor instances that have
     been passed to computePartialParsimony(), since the last time calculate() was called
     (if it ever was).
     @param start_index the "top" of the stack (the point at which to begin calculation)
     @param taskDescription a description of the task; may point to an empty string
     */
    void calculate(int start_index = 0, const char* taskDescription="");
};


#endif /* parallelparsimonycalculator_h */
