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

Pattern::Pattern()
        : vector<StateType>()
{
    frequency = 0;
//    is_const = false;
//    is_informative = false;
    flag = 0;
    const_char = 255;
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

int Pattern::computeGapChar(int num_states, int STATE_UNKNOWN) {
    int num = 0;
    for (iterator i = begin(); i != end(); i++)
        if (*i == STATE_UNKNOWN) num++;
    return num;
}

//Pattern &Pattern::operator= (Pattern pat) {
//    assign(pat);
//    frequency = pat.frequency;
//    is_const = pat.is_const;
//    const_char = pat.const_char;
//    return *this;
//}
