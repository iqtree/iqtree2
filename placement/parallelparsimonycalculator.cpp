//
// parallelparsimonycalculator.cpp
//
// Parallel calculation of parsimony (parallelizes across nodes a given
// "distance" from the branch for which parsimny is being calculated)
// (so, in large trees, close-to-linear scaling should be possible,
// regardless of the number of threads of execution).
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "parallelparsimonycalculator.h"

ParallelParsimonyCalculator::ParallelParsimonyCalculator(PhyloTree& phylo_tree)
    : tree(phylo_tree), task_to_start(nullptr), task_in_progress(nullptr) {}

void ParallelParsimonyCalculator::computePartialParsimony
    ( PhyloNeighbor* dad_branch, PhyloNode* dad ) {
    if (!dad_branch->isParsimonyComputed()) {
        workToDo.emplace_back(dad_branch, dad);
    }
}

int  ParallelParsimonyCalculator::computeParsimonyBranch
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      const char* task_description ) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);

    int startIndex = workToDo.size();
    computePartialParsimony(dad_branch,  dad);
    computePartialParsimony(node_branch, node);
    calculate(startIndex, task_description);
    return tree.computeParsimonyBranch(dad_branch, dad);
}

void ParallelParsimonyCalculator::calculate(int start_index,
                                            const char* task_description) {
    int stop_index = workToDo.size();
    bool tasked = (task_description!=nullptr && task_description[0]!='\0');
    if (stop_index <= start_index) {
        //Bail, if nothing to do
        return;
    }
    
    if (tasked && task_to_start==nullptr) {
        task_to_start   = task_description;
        task_in_progress = task_description;
    }
    
    //1. Find work to do at a lower level
    for (int i=stop_index-1; i>=start_index; --i) {
        auto item = workToDo[i];
        PhyloNeighbor* dad_branch = item.first;
        PhyloNode*     dad        = item.second;
        PhyloNode*     node       = dad_branch->getNode();
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
            computePartialParsimony(nei, node);
        }
    }
    
    //2. Do it, and then forget about it
    calculate(stop_index, "");
    workToDo.resize(stop_index);
    
    //3. Do the actual parsimony calculations at the
    //   current level (this doesn't change the content
    //   of workToDo so workToDo's contents can be
    //   processed with a parallel pointer loop).
    WorkItem* startItem = workToDo.data() + start_index;
    WorkItem* stopItem  = workToDo.data() + stop_index;
    if (task_to_start != nullptr) {
        tree.initProgress( workToDo.size(), task_to_start, "", "" );
        task_to_start = nullptr;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (WorkItem* item = startItem; item<stopItem; ++item) {
        PhyloNeighbor* dad_branch = item->first;
        PhyloNode*     dad        = item->second;
        tree.computePartialParsimony(dad_branch, dad);
    }
    if (task_in_progress != nullptr) {
        tree.trackProgress(stopItem-startItem);
    }
    workToDo.resize(start_index);
    if (tasked) {
        tree.doneProgress();
        task_in_progress = nullptr;
    }
}
