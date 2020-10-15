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
// There is no "top" of the heap at valueArray[start]; each other
// item in the array is, in turn, treated as though it were the top
// of the heap.  This trick makes the algorithm slightly shorter, and
// saves O(n) moves and O(log(n)) comparisons.
//

#ifndef heapsort_h
#define heapsort_h

#include "progress.h"

template <class V>
void constructMinHeapLayer ( V* valueArray, ptrdiff_t start,
                        ptrdiff_t layerStart, ptrdiff_t layerStop,
                        ptrdiff_t stop, progress_display* progress ) {
    //
    //Note: For lower levels in the heap, it would be better (due to the
    //      much better temporal locality of reference), to construct the
    //      heap "top down" (as per J.W.J. Williams' original Heapsort)
    //      rather than "bottom up" (as per Floyd's Treesort3) (the
    //      *later* algorithm that "stole" the original Heapsort's name).
    //      But the difference doesn't become important unless there are
    //      many millions of items in the heap.  So this implementation
    //      does things the Treesort3 way.
    //
    ptrdiff_t fudge         = 2 - start;
    ptrdiff_t nextLayerStop = layerStop + layerStop + fudge;
    if (nextLayerStop < stop) {
        ptrdiff_t nextLayerStart = layerStart + layerStart + fudge;
        if ( start+(stop-start)/2 + 1 < nextLayerStop ) {
            nextLayerStop = start+(stop-start)/2 + 1 ;
        }
        constructMinHeapLayer( valueArray, start, nextLayerStart, nextLayerStop, stop, progress);
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024)
    #endif
    for (ptrdiff_t h = layerStop-1; h >= layerStart; --h) {
        ptrdiff_t i = h;
        V         v = valueArray[i];
        ptrdiff_t j = i + i + fudge;
        while ( j < stop ) {
            if ( j + 1 < stop ) {
                j += ( valueArray[j+1] < valueArray[j] ) ? 1 : 0;
            }
            if ( v <= valueArray[j] ) {
                break;
            }
            valueArray[i] = valueArray[j];
            i = j;
            j = i + i + fudge;
        }
        valueArray[i] = v;
        if ( progress!=nullptr && (h&1023)==0 ) {
            progress->incrementBy(1024);
        }
    }
}

template <class V>
void constructMinHeap ( V* valueArray, ptrdiff_t start, ptrdiff_t stop
                      , const char* task_name ) {
    //
    //Constructs a "top-less" radix 2 min heap
    //of valueArray[start..stop].
    //
    ptrdiff_t count      = (stop-start);
    if (task_name==nullptr || *task_name=='\0') {
        constructMinHeapLayer(valueArray, start, start, start+2, stop, nullptr);
    } else {
        progress_display progress(count/2, task_name);
        if (2<count) {
            constructMinHeapLayer(valueArray, start, start, start+2, stop, &progress);
        }
        progress.done();
    }
}

template <class V>
void extractTopFromMinHeap ( V* valueArray
                        , ptrdiff_t start, ptrdiff_t stop) {
    //
    //Extracts a *single* value from a radix 2 min heap.
    //
    ptrdiff_t fudge = 2 - start;
    ptrdiff_t i = stop;
    V         v = valueArray[i];
    ptrdiff_t j = start;
    while ( j < stop ) {
        if ( j + 1 < stop ) {
            j += ( valueArray[j+1] < valueArray[j] ) ? 1 : 0;
        }
        if ( v <= valueArray[j] ) {
            break;
        }
        valueArray[i] = valueArray[j];
        i = j;
        j = i + i + fudge;
    }
    valueArray[i] = v;
}

template <class V> class MinHeapOnArray {
protected:
    V*     data;
    size_t count;
public:
    MinHeapOnArray(V* arrayStart, size_t elementCount, const char* taskName):
        data(arrayStart), count(elementCount) {
        constructMinHeap(data, 0, count, taskName);
    }
    V pop_min() {
        //ASSERT ( count!=0 );
        --count;
        extractTopFromMinHeap( data, 0, count );
        return data[count];
    }
};

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
    //Note: This could be parallelized
    //(as per constructMinHeapLayer), but as yet,
    //it hasn't.  It'd be better to have a separate
    //parallelConstructMirroredHeap function for that.
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
    //Extracts all the values (and associated sattelite information)
    //from a radix 2 max heap (and an array of sattelite information,
    //whose entries correspond one-to-one with the entries in the heap).
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
        constructMirroredHeap   ( valueArray, start, stop, satteliteArray );
        extractFromMirroredHeap ( valueArray, start, stop, satteliteArray );
    }
}

template <class T, class U>
void mirroredHeapsort( T& values, U& sattelite) {
    mirroredHeapsort ( values.data(), 0, values.size(), sattelite.data() );
}


#endif /* heapsort_h */
