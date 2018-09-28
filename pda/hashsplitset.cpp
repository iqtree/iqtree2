/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "hashsplitset.h"
#include "splitgraph.h"

Split *SplitIntMap::findSplit(Split *sp) {
    iterator ass_it = find(sp);
    if (ass_it != end()) {
        return ass_it->first;
    }
    return NULL;
}


Split *SplitIntMap::findSplit(Split *sp, int &value) {
    iterator ass_it = find(sp);
    if (ass_it != end()) {
        value = ass_it->second;
        return ass_it->first;
    }
    value = 0;
    return NULL;
}

int SplitIntMap::getValue(Split *sp) {
    int value;
    Split* findsp = findSplit(sp, value);
    ASSERT(findsp);
    return value;
}

void SplitIntMap::setValue(Split *sp, int value) {
    ASSERT(findSplit(sp));
    (*this)[sp] = value;
}

void SplitIntMap::eraseSplit(Split *sp) {
    ASSERT(findSplit(sp));
    erase(sp);
}

void SplitIntMap::insertSplit(Split *sp, int value) {
    ASSERT(!findSplit(sp));
    if (verbose_mode >= VB_MAX) sp->report(cout);
    (*this)[sp] = value;
}

void SplitIntMap::buildMap(SplitGraph &sg, bool use_index) {
    clear();
    for (int i = 0; i < sg.size(); i++) {
        if (use_index) 
            insertSplit(sg[i], i);
        else
            insertSplit(sg[i], sg[i]->getWeight());
    }
}
