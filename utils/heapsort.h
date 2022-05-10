//
//  heapsort.h
//  iqtree
//
//  Copyright (C) 2020, James Barbetti.
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//  Created by James Barbetti on 22/6/20.
//
// The heapsort algorithm here differs from Floyd's TreeSort3 [RF1964]
// (the algorithm that people mean, when they claim they've
//  based their code on J.W.J.William's Heapsort [JW1964] - when they
//  haven't, since [JW1964] was actually an out-of-place heapsort(!)):
// There is no "top" of the heap at vaueArray[start]; each other
// item in the array is, in turn, treated as though it were the top
// of the heap.  This trick makes the algorithm slightly shorter, and
// saves O(n) moves and O(log(n)) comparisons.
//

#ifndef heapsort_h
#define heapsort_h

#include <cstddef>

template <class V, class S>
void constructMirroredHeap ( V* valueArray
                            , ptrdiff_t start, ptrdiff_t stop
                            , S* satteliteArray) {
    //
    //Constructs a "top-less" radix 2 max heap
    //of valueArray[start..stop] and "mirrors"
    //the permutation of the values
    //on satteliteArray[start..stop].
    //
    ptrdiff_t fudge = 2 - start;
    for (ptrdiff_t h = start+(stop-start)/2; h >= start; --h) {
        ptrdiff_t i = h;
        V         v = valueArray[i];
        S         s = satteliteArray[i];
        ptrdiff_t j = i + i + fudge;
        while ( j < stop ) {
            if ( j + 1 < stop ) {
                j += ( valueArray[j] < valueArray[j+1] ) ? 1 : 0;
            }
            if ( valueArray[j] <= v ) {
                break;
            }
            valueArray[i] = valueArray[j];
            satteliteArray[i] = satteliteArray[j];
            i = j;
            j = i + i + fudge;
        }
        valueArray[i]     = v;
        satteliteArray[i] = s;
    }
}
template <class V, class S>
void extractFromMirroredHeap ( V* valueArray
                              , ptrdiff_t start, ptrdiff_t stop
                              , S* satteliteArray) {
    //
    //Extracts a value (and some sattelite information)
    //from a radix 2 max heap.
    //
    ptrdiff_t fudge = 2 - start;
    for (--stop; start<=stop; --stop) {
        ptrdiff_t i = stop;
        V         v = valueArray[i];
        S         s = satteliteArray[i];
        ptrdiff_t j = start;
        while ( j < stop ) {
            if ( j + 1 < stop ) {
                j += ( valueArray[j] < valueArray[j+1] ) ? 1 : 0;
            }
            if ( valueArray[j] <= v ) {
                break;
            }
            valueArray[i] = valueArray[j];
            satteliteArray[i] = satteliteArray[j];
            i = j;
            j = i + i + fudge;
        }
        valueArray[i]     = v;
        satteliteArray[i] = s;
    }
}

template <class V, class S>
void mirroredHeapsort ( V* valueArray
                       , size_t start, size_t stop
                       , S* satteliteArray) {
    //
    //Sorts valueArray[start..stop] and "mirrors"
    //the permutation of the values
    //on satteliteArray[start..stop].
    //
    if ( start + 1 < stop ) {
        constructMirroredHeap(valueArray, start, stop, satteliteArray);
        extractFromMirroredHeap(valueArray, start, stop, satteliteArray);
    }
}

#endif /* heapsort_h */
