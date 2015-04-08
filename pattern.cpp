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
#include "alignment.h"

Pattern::Pattern()
        : string()
{
    frequency = 0;
    is_const = false;
    const_char = 255;
}

Pattern::Pattern(const Pattern &pat)
        : string(pat)
{
    frequency = pat.frequency;
    is_const = pat.is_const;
    const_char = pat.const_char;
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
