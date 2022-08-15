#include "sequencechunk.h"

SequenceChunk::SequenceChunk()
{
    chunk_status = EMPTY;
    pos = 0;
    chunk_str = "";
}

SequenceChunk::~SequenceChunk()
{
    string().swap(chunk_str);
}
