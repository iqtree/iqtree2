#ifndef SEQUENCECHUNK_H
#define SEQUENCECHUNK_H

#include "utils/tools.h"

/** class storing a chunk of readable sequence waiting in writing queue */
class SequenceChunkStr {
public:
    SEQ_CHUNK_STATUS chunk_status;
    int64_t pos;
    string chunk_str;
    
    /**
        constructor
     */
    SequenceChunkStr();
    
    /**
        deconstructor
     */
    ~SequenceChunkStr();
    
};

#endif
