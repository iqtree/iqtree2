//
// parallel_mergesort.h
// Created by James Barbetti on 17-Mar-2021.
//

#ifndef parallel_mergesort_h
#define parallel_mergesort_h

#include <stdint.h>        //for intptr_t
#include <algorithm>       //for std::sort
#include "parallel_sort.h"

template <class T> void merge_to(const T *aStart, const T* aStop,
                                 const T* bStart, const T* bStop,
                                 T* dest) {
    //
    //Stably merges two non-overlapping, sorted ranges, 
    //writing sequentially to a (non-overlapping) output range.
    //
    //Assumes: 1. aStart<aStop, bStart<bStop.
    //         2. dest[0..(aStop-aStart + bStop-bStart)-1] writable.
    //         3. aStart[0..aStop-aStart-1] in sorted order.
    //         4. bStart[0..bStop-bStart-1] in sorted order.
    //         5. Neither of those ranges overlap with each other.
    //         6. Neither overlap with the range,
    //            dest[0..aStop-aStart+bStop-bStart-1].
    //
    //Notes:   1. When there is a tie the next value in the first range
    //            is taken (rather than the next in the second range)
    //         2. Element moves are unconditional! Pointer increments
    //            are determined via conditional moves (to reduce the
    //            number of pipeline stalls, if T is cheap to copy).
    //
    //            This makes sense when T is a primitive (or at least 
    //            a small) type. But it would be a very bad idea if T is 
    //            large and/or expensive to copy.
    //

    do {
        int aStep = ( *bStart < *aStart ) ? 0 : 1 ;
        *dest   = *aStart;
        aStart += aStep;
        dest   += aStep;
        if (aStop<=aStart) {
            //first input range exhausted.  
            //copy the remainder of the second.
            for (;bStart<bStop;++bStart,++dest) {
                *dest = *bStart;
            }
            return;
        }
        int bStep = 1 - aStep;
        *dest     = *bStart;
        bStart   += bStep;
        dest     += bStep;
    } while (bStart<bStop);
    //second input range exhausted. 
    //copy the remainder of the first.
    for (;aStart<aStop;++aStart,++dest) {
        *dest = *aStart;
    }
}

template <class T> void parallel_mergesort(T* data, T* buffer,
                                           intptr_t count) {
    if (count<256) {
        return std::sort(data, data+count);
    }

    //This is a straight mergesort, but with evenly-sized partitions
    //(as that is a more efficent way to divide the array into
    //sub-arrays, especially when count is just above a power of 2).
    //
    //Rational arithmetic is used to keep track of block sizes.
    //denominator is a power of 2.  The mean size of a range,
    //on each pass is (count/denominator).  denominator starts
    //at approximately count/16, and is halved on each pass.
    //
    //Data moves backwards and forwards between data and buffer
    //on successive passes.  There is always an even number of 
    //passes. power is the base-2 logarithm of denominator, and
    //is also the number of merging passes that remain.
    //
    //Each merge of a pair of blocks can be performed by a different
    //thread (as they are independent).  There is no parallelization
    //within individual merging operations which means that there are
    //powerof2near(count/16), powerof2near(count/8), .... , 8, 4, 2, 1
    //threads of execution (and the last merging pass is single-
    //threaded).  It *would* be possible to further parallelize the
    //last few passes, but that hasn't been done.

    intptr_t denominator  = 1;
    int      power = 0;
    for (; denominator<count; denominator+=denominator) {
        ++power;
    }
    denominator >>= 4;
    power-=4;
    if (power & 1) {
        //Want, an even number of merging passes, so the
        //sorted output ends up in data, not buffer.
        denominator >>= 1;
        --power;
    }
    intptr_t numerator_stop = count*denominator;

    //Sort the small blocks before the first merge pass
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t numerator=0; numerator<numerator_stop;
         numerator+=count) {
        intptr_t start = numerator>>power;
        intptr_t stop  = (numerator+count)>>power;
        std::sort(data+start, data+stop);
    }
                  
    T* source = data;
    T* dest   = buffer;
    do {
        intptr_t double_count = count+count;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t numerator=0; numerator<numerator_stop;
             numerator+=double_count) {
            intptr_t start1 = numerator/denominator;
            intptr_t start2 = (numerator+count)/denominator;
            intptr_t stop   = (numerator+double_count)/denominator;
            if (count<stop) {
                stop=count;
            }
            if (start2<count) {
                merge_to(source+start1, source+start2,
                         source+start2, source+stop,
                         dest+start1 );
            }
        }
        std::swap(source,dest);
        denominator >>= 1;
        numerator_stop = count*denominator;
        --power;
    } while (denominator>1); //or 0<power if you prefer. They're equivalent.
}

template <class T> void single_thread_mergesort(T* data, T* buffer,
                                                intptr_t count) {
    //The same as the parallel version minus the omp directives
    //(if and when the parallel one does parallel merging -
    //parallelizing individual block merges - this won't
    //look quite so similar).
    if (count<256) {
        return std::sort(data, data+count);
    }
    intptr_t denominator  = 1;
    int      power = 0;
    for (; denominator<count; denominator+=denominator) {
        ++power;
    }
    denominator >>= 4;
    power-=4;
    if (power & 1) {
        //Want, an even number of merging passes, so the
        //sorted output ends up in data, not buffer.
        denominator >>= 1;
        --power;
    }
    intptr_t numerator_stop = count*denominator;

    for (intptr_t numerator=0; numerator<numerator_stop;
         numerator+=count) {
        intptr_t start = numerator>>power;
        intptr_t stop  = (numerator+count)>>power;
        std::sort(data+start, data+stop);
    }
                  
    T* source = data;
    T* dest   = buffer;
    do {
        intptr_t double_count = count+count;
        for (intptr_t numerator=0; numerator<numerator_stop;
             numerator+=double_count) {
            intptr_t start1 = numerator/denominator;
            intptr_t start2 = (numerator+count)/denominator;
            intptr_t stop   = (numerator+double_count)/denominator;
            if (count<stop) {
                stop=count;
            }
            if (start2<count) {
                merge_to(source+start1, source+start2,
                         source+start2, source+stop,
                         dest+start1 );
            }
        }
        std::swap(source,dest);
        denominator >>= 1;
        numerator_stop = count*denominator;
        --power;
    } while (denominator>1);
}

template <class T> class MergeSorter: public ParallelSorter<T> {
protected:
    std::vector<T> buffer;
    void allocateBuffer(intptr_t count) {
        buffer.resize(count);
    }
public:
    void parallel_sort(T* data, intptr_t count) override {
        allocateBuffer(count);
        parallel_mergesort(data, buffer.data(), count);
    }
    void single_thread_sort(T* data, intptr_t count) override {
        allocateBuffer(count);
        single_thread_mergesort(data, buffer.data(), count);
    }
};

template <class T, class U, class P=ValueAndSattelite<T,U>,
          class S=MergeSorter<P> >
class MirrorMergeSorter: public ParallelMirrorSorterBase<T,U,P,S>  {
};

#endif /* parallel_mergesort_h */
