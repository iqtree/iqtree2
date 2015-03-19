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


