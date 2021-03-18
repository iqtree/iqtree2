//
// parallelparsimonycalculator.cpp
// ===============================
// Parallel calculation of parsimony (parallelizes across nodes a given
// "distance" from the branch for which parsimony is being calculated)
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

double ParallelParsimonyCalculator::parsimonyLink4Cost(PhyloNode* a, PhyloNode* b, PhyloNode* c,
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
        PhyloBranch branch (branches[i]);
        ASSERT(branch.first != nullptr && branch.second != nullptr);
        PhyloNeighbor* nei = branch.first->findNeighbor(branch.second);
        computes += schedulePartialParsimony(nei, branch.first);
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

double ParallelParsimonyCalculator::parsimonyLink4CostOutOfTree
       ( const PhyloTree& tree, PhyloNode* a, PhyloNode* b, PhyloNode* c,
         PhyloNode* d, PhyloNode* e, PhyloNode* f,
         UINT* buffer1, UINT* buffer2) {
    std::vector<UINT*> inputs;
    inputs.resize(4, nullptr);

    FOR_EACH_ADJACENT_PHYLO_NODE(c, d, it, x) {
        if      (x==a) inputs[0] = c->findNeighbor(a)->get_partial_pars();
        else if (x==b) inputs[1] = c->findNeighbor(b)->get_partial_pars();
        else if (x==e) inputs[2] = c->findNeighbor(e)->get_partial_pars();
        else if (x==f) inputs[3] = c->findNeighbor(f)->get_partial_pars();
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(d, c, it, x) {
       if      (x==a) inputs[0] = d->findNeighbor(a)->get_partial_pars();
       else if (x==b) inputs[1] = d->findNeighbor(b)->get_partial_pars();
       else if (x==e) inputs[2] = d->findNeighbor(e)->get_partial_pars();
       else if (x==f) inputs[3] = d->findNeighbor(f)->get_partial_pars();
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

    intptr_t startIndex = workToDo.size();
    schedulePartialParsimony(dad_branch,  dad);
    schedulePartialParsimony(node_branch, node);
    calculate(startIndex, task_description);
    return tree.computeParsimonyBranch(dad_branch, dad);
}

int  ParallelParsimonyCalculator::computeAllParsimony
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      const char* task_description ) {
    int        result = computeParsimonyBranch(dad_branch, dad, task_description);
    PhyloNode* node   = dad_branch->getNode();
    computeReverseParsimony(dad, node);
    return result;
}

void ParallelParsimonyCalculator::computeReverseParsimony(PhyloNode* first,
                                                          PhyloNode* second) {
    std::vector <PhyloBranch> stuffToDo;
    stuffToDo.emplace_back(first, second);
    stuffToDo.emplace_back(second, first);
    while (!stuffToDo.empty()) {
        std::vector <PhyloBranch> stuffToDoNext;
        intptr_t r = 0;
        intptr_t w = 0;
        for (;r<static_cast<intptr_t>(stuffToDo.size()); ++r) {
            PhyloNode* first  = stuffToDo[r].first;
            PhyloNode* second = stuffToDo[r].second;
            FOR_EACH_ADJACENT_PHYLO_NODE(first, second, it, back) {
                stuffToDoNext.emplace_back(back, first);
            }
            //Only keep the (PhyloNeighbor*, PhyloNode*)
            //pairs where we need to compute reverse parsimony now.
            stuffToDo[w]       = stuffToDo[r];
            PhyloNeighbor* nei = first->findNeighbor(second);
            w                 += ( nei->isParsimonyComputed() ? 0 : 1 );
        }
        stuffToDo.resize(w);
        PhyloBranch* firstItem = stuffToDo.data();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t i = 0; i < w; ++i) {
            PhyloBranch*   item = firstItem + i;
            PhyloNeighbor* nei  = item->first->findNeighbor(item->second);
            tree.computePartialParsimony(nei, item->first);
            if (report_progress_to_tree && (i%1000) == 999) {
                tree.trackProgress(1000.0);
            }
        }
        if (report_progress_to_tree) {
            tree.trackProgress(static_cast<double>(w%1000));
        }
        std::swap(stuffToDo, stuffToDoNext);
    }
}

void ParallelParsimonyCalculator::calculate(intptr_t start_index,
                                            const char* task_description) {
    intptr_t stop_index = workToDo.size();
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
    for (intptr_t  i=stop_index-1; i>=start_index; --i) {
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
    WorkItem* item_data = workToDo.data();
    if (task_to_start != nullptr) {
        double estimate = static_cast<double>(workToDo.size());
        tree.initProgress( estimate, task_to_start, "", "" );
        task_to_start = nullptr;
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (intptr_t i = start_index; i < stop_index; ++i) {
        WorkItem*      item       = item_data + i;
        PhyloNeighbor* dad_branch = item->first;
        PhyloNode*     dad        = item->second;
        
        tree.computePartialParsimony(dad_branch, dad);
        
        if (task_in_progress != nullptr || report_progress_to_tree) {
            intptr_t j = i - start_index;
            if ((j%1000)==999) {
                tree.trackProgress(1000.0);
            }
        }
    }
    if (task_in_progress != nullptr || report_progress_to_tree) {
        double work_done = (double)((stop_index - start_index)%1000);
        tree.trackProgress(work_done);
    }
    workToDo.resize(start_index);
    if (tasked) {
        tree.doneProgress();
        task_in_progress = nullptr;
    }
}
