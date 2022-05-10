#include "sequence.h"

Sequence::Sequence() {
    num_children_done_simulation = 0;
    num_threads_done_simulation = 0;
    num_gaps = 0;
    insertion_pos = NULL;
    parent = NULL;
}

