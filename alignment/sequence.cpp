#include "sequence.h"

Sequence::Sequence() {
    nums_children_done_simulation.resize(1);
    sequence_chunks.resize(1);
    num_threads_done_simulation = 0;
    num_threads_reach_barrier = 0;
    num_gaps = 0;
    depth = 0;
    insertion_pos = NULL;
    parent = NULL;
}

