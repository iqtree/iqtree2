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

ParallelParsimonyCalculator::ParallelParsimonyCalculator(PhyloTree& phylo_tree,
                                                         bool report_progress)
    : tree(phylo_tree), task_to_start(nullptr)
    , task_in_progress(nullptr), report_progress_to_tree(report_progress)
    {}

int ParallelParsimonyCalculator::schedulePartialParsimony
    ( PhyloNeighbor* dad_branch, PhyloNode* dad ) {
    if (dad_branch->isParsimonyComputed()) {
        return 0;
    }
    workToDo.emplace_back(dad_branch, dad);
    return 1;
}

int ParallelParsimonyCalculator::schedulePartialParsimony(const PhyloBranch& branch) {
    PhyloNeighbor* dad_branch = branch.first->findNeighbor(branch.second);
    if (dad_branch->isParsimonyComputed()) {
        return 0;
    }
    workToDo.emplace_back(dad_branch, branch.first);
    return 1;
}

int ParallelParsimonyCalculator::parsimonyLink4Cost(PhyloNode* a, PhyloNode* b, PhyloNode* c,
                                                    PhyloNode* d, PhyloNode* e, PhyloNode* f,
                                                    UINT* buffer1, UINT* buffer2) {
    PhyloBranchVector branches;
    branches.resize(4);
    FOR_EACH_ADJACENT_PHYLO_NODE(c, d, it, x) {
        if      (x==a) branches[0] = PhyloBranch(c, a);
        else if (x==b) branches[1] = PhyloBranch(c, b);
        else if (x==e) branches[2] = PhyloBranch(c, e);
        else if (x==f) branches[3] = PhyloBranch(c, f);
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(d, c, it, x) {
        if      (x==a) branches[0] = PhyloBranch(d, a);
        else if (x==b) branches[1] = PhyloBranch(d, b);
        else if (x==e) branches[2] = PhyloBranch(d, e);
        else if (x==f) branches[3] = PhyloBranch(d, f);
    }
    int computes=0;
    std::vector<UINT*> inputs;
    PhyloNeighborVec neighbors;
    for (int i=0; i<4; ++i) {
        PhyloBranch b (branches[i]);
        ASSERT(b.first != nullptr && b.second != nullptr);
        PhyloNeighbor* nei = b.first->findNeighbor(b.second);
        computes += schedulePartialParsimony(nei, b.first);
        inputs.push_back(nei->get_partial_pars());
    }
    if (0<computes) {
        calculate();
    }
    double score = 0;
    score += tree.computePartialParsimonyOutOfTree(inputs[0], inputs[1], buffer1);
    score += tree.computePartialParsimonyOutOfTree(inputs[2], inputs[3], buffer2);
    int branchCost;
    tree.computeParsimonyOutOfTree(buffer1, buffer2, &branchCost);
    score += branchCost;
    return score;
}

int  ParallelParsimonyCalculator::computeParsimonyBranch
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      const char* task_description ) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);

    int startIndex = workToDo.size();
    schedulePartialParsimony(dad_branch,  dad);
    schedulePartialParsimony(node_branch, node);
    calculate(startIndex, task_description);
    return tree.computeParsimonyBranch(dad_branch, dad);
}

int  ParallelParsimonyCalculator::computeAllParsimony
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      const char* task_description ) {
    int result = computeParsimonyBranch(dad_branch, dad, task_description);
    computeReverseParsimony(dad_branch, dad);
    PhyloNode*     node     = dad_branch->getNode();
    PhyloNeighbor* back_nei = node->findNeighbor(dad);
    computeReverseParsimony(back_nei, node);
    return result;
}

void ParallelParsimonyCalculator::computeReverseParsimony(PhyloNeighbor* dad_branch,
                                                          PhyloNode* dad) {
    std::vector <WorkItem> stuffToDo;
    stuffToDo.emplace_back(dad_branch, dad);
    while (!stuffToDo.empty()) {
        std::vector <WorkItem> stuffToDoNext;
        size_t r = 0;
        size_t w = 0;
        for (;r<stuffToDo.size(); ++r) {
            PhyloNeighbor* nei   = stuffToDo[r].first;
            PhyloNode*     node  = stuffToDo[r].second;
            PhyloNode*     child = nei->getNode();
            FOR_EACH_PHYLO_NEIGHBOR(child, node, it, nei2) {
                stuffToDoNext.emplace_back(nei2, child);
            }
            PhyloNeighbor* back_nei = child->findNeighbor(node);
            //Only keep the (PhyloNeighbor*, PhyloNode*)
            //pairs where we need to compute reverse parsimony now.
            stuffToDo[w].first  = back_nei; //reverse direction!
            stuffToDo[w].second = child;
            w                  += ( back_nei->isParsimonyComputed() ? 0 : 1 );
        }
        stuffToDo.resize(w);
        WorkItem* startItem = stuffToDo.data();
        WorkItem* stopItem  = stuffToDo.data() + w;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (WorkItem* item = startItem; item<stopItem; ++item) {
            PhyloNeighbor* dad_branch = item->first;
            PhyloNode*     dad        = item->second;
            tree.computePartialParsimony(dad_branch, dad);
        }
        if (report_progress_to_tree) {
            tree.trackProgress(r);
        }
        std::swap(stuffToDo, stuffToDoNext);
    }
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
        task_to_start    = task_description;
        task_in_progress = task_description;
    }
    
    //1. Find work to do at a lower level
    for (int i=stop_index-1; i>=start_index; --i) {
        auto item = workToDo[i];
        PhyloNeighbor* dad_branch = item.first;
        PhyloNode*     dad        = item.second;
        PhyloNode*     node       = dad_branch->getNode();
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
            schedulePartialParsimony(nei, node);
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
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (WorkItem* item = startItem; item<stopItem; ++item) {
        PhyloNeighbor* dad_branch = item->first;
        PhyloNode*     dad        = item->second;
        tree.computePartialParsimony(dad_branch, dad);
    }
    if (task_in_progress != nullptr || report_progress_to_tree) {
        tree.trackProgress(stopItem-startItem);
    }
    workToDo.resize(start_index);
    if (tasked) {
        tree.doneProgress();
        task_in_progress = nullptr;
    }
}
