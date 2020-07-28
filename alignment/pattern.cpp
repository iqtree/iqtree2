//
// C++ Implementation: pattern
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pattern.h"
#include "alignment/alignment.h"
#include <vectorclass/vectorclass.h>

Pattern::Pattern()
        : vector<StateType>()
{
    frequency = 0;
//    is_const = false;
//    is_informative = false;
    flag = 0;
    const_char = -1;
    num_chars = 0;
}

Pattern::Pattern(int nseq, int freq)
: vector<StateType>(nseq)
{
    frequency = freq;
    //    is_const = false;
    //    is_informative = false;
    flag = 0;
    const_char = -1;
    num_chars = 0;
}

Pattern::Pattern(const Pattern &pat)
        : vector<StateType>(pat)
{
    frequency = pat.frequency;
//    is_const = pat.is_const;
//    is_informative = pat.is_informative;
    flag = pat.flag;
    const_char = pat.const_char;
    num_chars = pat.num_chars;
}

Pattern::~Pattern()
{
}

int Pattern::computeAmbiguousChar(int num_states) {
    int num = 0;
    for (iterator i = begin(); i != end(); i++)
        if (*i >= num_states) num++;
    return num;
}

#define VECTORIZE_GAPCHAR_COUNT 1
int Pattern::computeGapChar(int num_states, int STATE_UNKNOWN) const {
    int num = 0;
#if VECTORIZE_GAPCHAR_COUNT
    //This won't compile unless value_type is based on uint32_t
    //(nor should it! You'd need to use different vector types!)
    const uint32_t* dataStart = data();
    size_type  count   = size();
    size_type  vecSize = Vec8ui::size();
    Vec8ui     unknown = STATE_UNKNOWN;
    const uint32_t* dataStop   = dataStart + count;
    const uint32_t* blockStop  = dataStop - (count & (vecSize-1));
    for (const uint32_t* block=dataStart; block<blockStop; block+=vecSize) {
        Vec8ui a;
        a.load(block);
        num -= horizontal_add( Vec8ui(a == unknown) );
    }
    for (const uint32_t* single=blockStop; single<dataStop; ++single) {
        if (*single == STATE_UNKNOWN) {
            ++num;
        }
    }
#else
    for (iterator i = begin(); i != end(); i++)
        if (*i == STATE_UNKNOWN) num++;
#endif
    return num;
}

//Pattern &Pattern::operator= (Pattern pat) {
//    assign(pat);
//    frequency = pat.frequency;
//    is_const = pat.is_const;
//    const_char = pat.const_char;
//    return *this;
//}
