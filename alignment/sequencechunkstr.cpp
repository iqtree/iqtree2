#include "sequencechunkstr.h"

SequenceChunkStr::SequenceChunkStr()
{
    chunk_status = EMPTY;
    pos = 0;
    chunk_str = "";
}

SequenceChunkStr::~SequenceChunkStr()
{
    string().swap(chunk_str);
}
