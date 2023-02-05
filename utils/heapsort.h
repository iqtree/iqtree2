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
//
// There is no "top" of the heap at valueArray[start]; each other
// item in the array is, in turn, treated as though it were the top
// of the heap.  This trick makes the algorithm slightly shorter, and
// saves O(n) moves and O(log(n)) comparisons.
//

#ifndef heapsort_h
#define heapsort_h

#include <stddef.h>          //for ptrdiff_t
#include "progress.h"
#include <utils/my_assert.h> //for ASSERT macro

/**
 * @brief  Given a partially constructed heap (where layers below the 
 *         layerStart..(layerStop-1) layer of the heap have been 
 *         constructed, and each item in the index range 
 *         layerStart..(layerStop-1) will be the "top" of its own
 *         subheap, once it is insert into that subheap)... insert the
 *         items with indices layerStart..(layerStop-1) into their
 *         subheaps.
 * @tparam V          The element type (must implement operator< and operator<=)
 * @param  valueArray The array in which the heap is being constructed
 * @param  start      The index of the first item in the heap
 * @param  layerStart The index of the first item in this layer of the heap
 * @param  layerStop  One more than the index of the last item in this layer
 * @param  stop       One more than the index of the last item in the heap
 * @param  progress   A pointer to something which keeps track of progress
 * @note   If _OPENMP is defined and nn-zero, the construction runs in
 *         parallel. Each element has its own "sub-heap" and the operations
 *         in each subheap are independent of those in all the others.
 * @note   For lower levels in the heap, it would be better (due to the
 *         much better temporal locality of reference), to construct the
 *         heap "top down" (as per J.W.J. Williams' original Heapsort)
 *         rather than "bottom up" (as per Floyd's Treesort3) (the
 *         *later* algorithm that "stole" the original Heapsort's name).
 *         But the difference doesn't become important unless there are
 *         many millions of items in the heap.  So this implementation
 *         does things the Treesort3 way.  In any case, the "Williams way"
 *         would be much more difficult to parallelize.
 * @note   The heap is a MIN heap, with each item in a subheap the *minimum*
 *         NOT the maximum (you might expected) of the items in "lower" levels
 *         of the heap.  
 */
template <class V>
void constructMinHeapLayer ( V* valueArray, ptrdiff_t start,
                            ptrdiff_t layerStart, ptrdiff_t layerStop,
                            ptrdiff_t stop,
                            progress_display_ptr progress = nullptr ) {
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
            *progress += 1024.0;
        }
    }
}

/**
 * @brief  Construct "top-less" radix 2 min heap of valuArray[start..stop]
 * @tparam V          the type of element, in the array.   Must support 
 *                    operator< and operator<=.
 * @param  valueArray the array, in which the heap is to be constructed
 * @param  start      the index of the first item in what is to be the heap.
 * @param  stop       one more than the index of the last item
 * @param  task_name 
 * @note   If _OPENMP is defined and nn-zero, the construction of each 
 *         "level" of the heap runs in parallel. Each element at the
 *         current level has its own "sub-heap" and the operations in 
 *         each subheap are independent of those in all the others.
 */
template <class V>
void constructMinHeap ( V* valueArray, ptrdiff_t start, ptrdiff_t stop
                      , const char* task_name ) {
    size_t count      = (stop-start);
    if (task_name==nullptr || *task_name=='\0') {
        constructMinHeapLayer(valueArray, start, start, start+2, stop, nullptr);
    } else {
        #if USE_PROGRESS_DISPLAY
        progress_display progress(count/2, task_name);
        #else
        double progress = 0.0;
        #endif
        if (2<count) {
            constructMinHeapLayer(valueArray, start, start, start+2, stop, &progress);
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
    }
}

/**
 * @brief  Extracts a *single* value from a radix 2 min heap,
 *         in valueArray[start..(stop-1)], and - at the same time - 
 *         inserts the value, v, that was at valueArray[start] in
 *         that min heap, it is place.
 * @tparam V          - the element type (must support operator< and operator<=)
 * @param  valueArray - the array in which the heap has been built/maintained
 * @param  start      - the index of the first element in the heap
 * @param  stop       - one more than the last
 */
template <class V>
void extractTopFromMinHeap ( V* valueArray
                        , ptrdiff_t start, ptrdiff_t stop) {
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

/**
 * @brief  A min-heap, on a naked array, of V
 * @tparam V the type, of the elements, in the array
 */
template <class V> class MinHeapOnArray {
protected:
    V*     data;  //pointer to the start of the heap
    size_t count; //the number of items currently in the heap
                  //(initially, the number of items in the naked array range 
                  //to be assembled into a MIN heap; later, the number of
                  //items left in the MIN heap) (items are, as the extracted,
                  //written to the "right" of what remains of the heap).
public:
    MinHeapOnArray(V* arrayStart, size_t elementCount, const char* taskName):
        data(arrayStart), count(elementCount) {
        constructMinHeap(data, 0, count, taskName);
    }
    /**
     * @brief  Extract the minimum item, currently in the heap, from
     *         the heap (preserving the heap properties, but reducing
     *         the size of the heap by 1).
     * @return V - the (minimum) item, popped from the heap.
     * @note   No index-range checking is performed.  Nor is there any 
     *         check that 0<count. But it is assumed that 0<count.
     */
    V pop_min() {
        //ASSERT ( count!=0 );
        --count;
        extractTopFromMinHeap( data, 0, count );
        return data[count];
    }

    /**
     * @brief  As for pop_min(), except that there's a check that the
     *         heap is not empty.
     * @return V - the (minimum) item, popped from the heap.
     */
    V safe_pop_min() {
        ASSERT( 0 < count );
        --count;
        extractTopFromMinHeap( data, 0, count );
        return data[count];
    }

    /**
     * @brief  Either pop the minimum item from the heap, or return false
     *         to indicate that the heap is empty
     * @param  v     - a reference, to a variable, which is to hold
     *                 what was the minimum item in the heap, if there
     *                 were any items in the heap.
     * @return true  - if the heap wasn't empty
     * @return false - if the heap was empty
     */
    bool pop_min_into(V& v) {
        if (0==count) {
            return false;
        }
        v = pop_min();
        return true;
    }

    bool empty() const {
        return count==0;
    }

    size_t size() const {
        return count;
    }
};


/**
 * @brief  Construct a "mirrored" heap (suitable for sorting one array, and
 *         "mirroring" what happens on that array, in *another* array, to
 *         maintain the existing one to one correspondence between the items
 *         in the two arrays, by reordering both in exactly the same way).
 * @tparam V - the type of the items in the value array (the array of sort
 *             keys).  Must support operator< and operator<=.
 * @tparam S - the type of the items in the sattelite array (the array of
 *             values to be ordered by the keys).
 * @param  valueArray - the array of comparable value items
 * @param  start      - the index of the first item (in valueArray) to
 *                      be included in the mirrored heap
 * @param  stop       - one more than the index of the last item 
 *                      (in valueArray) to be included in the mirrored heap
 * @param  satteliteArray - the array of sattelite items
 * @note   This function constructs a "top-less" radix 2 max heap of 
 *         valueArray[start..stop] and "mirrors" the permutation of the 
 *         values on satteliteArray[start..stop].  Note: a MAX heap.
 * @note   This could be parallelized (as per constructMinHeapLayer), 
 *         but as yet, it hasn't.  It'd be better to have a separate
 *         parallelConstructMirroredHeap function for that. -James B.
 */
template <class V, class S>
void constructMirroredHeap ( V* valueArray
                            , ptrdiff_t start, ptrdiff_t stop
                            , S* satteliteArray) {
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

/**
 * @brief  Extracts all the values (and associated sattelite information)
 *         from a radix 2 max heap (and an array of sattelite information,
 *         whose entries correspond one-to-one with the entries in the heap).
 * @tparam V - the type of the items in the value array (the array of sort
 *             keys).  Must support operator< and operator<=.
 * @tparam S - the type of the items in the sattelite array (the array of
 *             values to be ordered by the keys).
 * @param valueArray - the array of comparable value items
 * @param start      - the index of the first item (in valueArray) to
 *                     be included in the mirrored heap
 * @param stop       - one more than the index of the last item 
 *                     (in valueArray) to be included in the mirrored heap
 * @param satteliteArray - the array of sattelite items
 */
template <class V, class S>
void extractFromMirroredHeap ( V* valueArray
                              , ptrdiff_t start, ptrdiff_t stop
                              , S* satteliteArray) {
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
            valueArray[i]     = valueArray[j];
            satteliteArray[i] = satteliteArray[j];
            i = j;
            j = i + i + fudge;
        }
        valueArray[i]     = v;
        satteliteArray[i] = s;
    }
}

/**
 * @brief  Sorts valueArray[start..stop] and "mirrors"
 *         the permutation of the values on satteliteArray[start..stop].
 * @tparam V - the type of the items in the value array (the array of sort
 *             keys).  Must support operator< and operator<=.
 * @tparam S - the type of the items in the sattelite array (the array of
 *             values to be ordered by the keys).
 * @param  valueArray - the array of comparable value items
 * @param  start      - the index of the first item (in valueArray) to
 *                      be included in the mirrored heap
 * @param  stop       - one more than the index of the last item 
 *                      (in valueArray) to be included in the mirrored heap
 * @param  satteliteArray 
 * @note   In practice, V and S *could* be the same, and the arrays could
 *         be the same, but what would the point be, of moving eerything
 *         twice.
 * @note   It is assumed that the two array ranges do not overlap (!).
 *         But this is not checked.
 */
template <class V, class S>
void mirroredHeapsort ( V* valueArray
                       , size_t start, size_t stop
                       , S* satteliteArray) {
    //
    //
    //
    if ( start + 1 < stop ) {
        constructMirroredHeap   ( valueArray, start, stop, satteliteArray );
        extractFromMirroredHeap ( valueArray, start, stop, satteliteArray );
    }
}

/**
 * @brief  Duck-type sorting of vectors (or other types with size() and
 *         data() member functions)
 * @tparam T - the vector type, of the first vector, which contains values;
 *             to be compared and sorted. The value type is whatever is 
 *             pointed to by the return value of T::size(). The value type
 *             must implement operator< and operator<=.
 * @tparam U - the vector type, of the second, "mirrored" vector,
 *             which contains sattelite values, which are to be rearranged
 *             in a fashion that mirrors what happens to the first vector
 *             (or: is to permuted according to the same permutation).
 *             The sattelite type is whatever is pointed to by
 */
template <class T, class U>
void mirroredHeapsort( T& values, U& sattelite) {
    mirroredHeapsort ( values.data(), 0, values.size(), sattelite.data() );
}

#endif /* heapsort_h */
