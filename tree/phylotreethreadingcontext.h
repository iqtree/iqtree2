//
//  phylotreethreadingcontext.h
//  alignment
//
//  Created by James Barbetti on 21/1/21.
//

#ifndef phylotreethreadingcontext_h
#define phylotreethreadingcontext_h

#include "phylotree.h"

class PhyloTreeThreadingContext {
public:
    PhyloTree& tree;
    int  old_num_threads;
    int  old_thread_count;         //OMP's thread count
    bool was_omp_thread_count_set; //True if it was adjusted upward
    
    PhyloTreeThreadingContext(PhyloTree& tree, bool force_use_of_all_threads);
    int getThreadNumber() const;
    int getThreadCount()  const;
    static int getMaximumThreadCount();

    ~PhyloTreeThreadingContext();
};

#endif /* phylotreethreadingcontext_h */
