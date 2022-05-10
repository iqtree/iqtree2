#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "utils/tools.h"
#include "tree/genometree.h"

using namespace std;

class Sequence {
public:
    
    /**
     chunks of sequence
     */
    vector<short int> sequence_chunks;
    
    /**
        chunks of sequence (in string)
     */
    string sequence_str_chunks;
    
    /**
        number of children which have completed simulating the sequence (for AliSim)
     */
    short int num_children_done_simulation;
    
    /**
        number of OPENMP threads that completed simulating the sequence (for AliSim)
     */
    short int num_threads_done_simulation;
    
    /**
        pointer to the position of the insertion event that occurs after simulating sequence at this node
     */
    Insertion* insertion_pos;
    
    /**
        parent node of the current node (only use when simulating Indels with AliSim)
     */
    Node* parent;
    
    /**
        number of gaps in the sequence
     */
    int num_gaps;

    /**
        constructor
     */
    Sequence();
    
    /**
        deconstructor
     */
    ~Sequence() {
        // do nothing, pointers should be deallocated by other tasks
    }
};

#endif
