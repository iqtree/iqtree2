//
//  phylotreethreadingcontext.cpp
//  alignment
//
//  Created by James Barbetti on 21/1/21.
//

#include "phylotreethreadingcontext.h"
#include <omp.h> //for omp_get_num_threads and omp_get_num_procs

PhyloTreeThreadingContext::PhyloTreeThreadingContext(PhyloTree& phylo_tree,
                                                     bool force_use_of_all_threads)
: tree(phylo_tree), old_num_threads(phylo_tree.num_threads)
, old_thread_count(0), was_omp_thread_count_set(false) {
#ifdef _OPENMP
    old_thread_count = omp_get_num_threads();
    auto max_cores   = omp_get_num_procs();
#endif
    if (old_num_threads==0) {
        if (0<tree.params->num_threads) {
            tree.num_threads = tree.params->num_threads;
            if (force_use_of_all_threads) {
                #ifdef _OPENMP
                    if (old_thread_count < max_cores) {
                        TREE_LOG_LINE(tree, VB_MIN, "Temporarily elevating number of threads to " << max_cores);
                        omp_set_num_threads(max_cores);
                        was_omp_thread_count_set = true;
                    }
                    tree.num_threads = omp_get_max_threads();
                #endif
            }
        } else {
            #ifdef _OPENMP
            tree.num_threads = omp_get_max_threads();
            #else
            tree.num_threads = 0;
            #endif
        }
        return;
    }
#ifdef _OPENMP
    if (!force_use_of_all_threads) {
        return;
    }
    TREE_LOG_LINE(tree, VB_MIN, "Maximum number of cores is " << max_cores);
    if (max_cores <= tree.num_threads) {
        return;
    }
    old_thread_count = omp_get_num_threads();
    if (old_thread_count < max_cores) {
        TREE_LOG_LINE(tree, VB_MIN, "Temporarily elevating number of threads to " << max_cores);
        omp_set_num_threads(max_cores);
        was_omp_thread_count_set = true;
    }
    tree.num_threads = omp_get_max_threads();
#endif
}

int PhyloTreeThreadingContext::getThreadNumber() const {
    #ifdef _OPENMP
        return omp_get_thread_num();
    #else
        return 0;
    #endif
}

int PhyloTreeThreadingContext::getThreadCount() const {
    return tree.num_threads;
}

/*static*/ int PhyloTreeThreadingContext::getMaximumThreadCount() {
    #ifdef _OPENMP
        return omp_get_max_threads();
    #else
        return 0;
    #endif
}

PhyloTreeThreadingContext::~PhyloTreeThreadingContext() {
    if (old_num_threads != tree.num_threads) {
        tree.num_threads = old_num_threads;
    }
    if (was_omp_thread_count_set) {
        omp_set_num_threads(old_thread_count);
    }
}
