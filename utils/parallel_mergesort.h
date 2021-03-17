//
// parallel_mergesort.h
// Created by James Barbetti on 17-Mar-2021.
//

#ifndef parallel_mergesort_h
#define parallel_mergesort_h

#include "parallel_sort.h"

template <class T> void merge_to(T *aStart, T* aStop,
                                 T* bStart, T* bStop,
                                 T* dest) {
    //
    //Assumes: aStart<aStop, bStart<bStop.
    //         dest[0..(aStop-aStart + bStop-bStart)-1] writable.
    //         aStart[0..aStop-aStart-1] in sorted order.
    //         bStart[0..bStop-bStart-1] in sorted order.
    do {
        int aStep = ( *aStart < *bStart ) ? 1 : 0 ;
        *dest   = *aStart;
        aStart += aStep;
        dest   += aStep;
        if (aStop<=aStart) {
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
    for (;aStart<aStop;++aStart,++dest) {
        *dest = *aStart;
    }
}

template <class T> void parallel_mergesort(T* data, T* buffer,
                                           intptr_t count) {
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
    } while (denominator>1);
}

template <class T> void single_thread_mergesort(T* data, T* buffer,
                                                intptr_t count) {
    //The same as the parallel version minus the omp directives
    //(if when the parallel one does parallel merging this won't
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
    void parallel_sort(T* data, intptr_t count) {
        allocateBuffer(count);
        parallel_mergesort(data, buffer.data(), count);
    }
    void single_thread_sort(T* data, intptr_t count) {
        allocateBuffer(count);
        single_thread_mergesort(data, buffer.data(), count);
    }
};

template <class T, class U, class P=ValueAndSattelite<T,U>,
          class S=MergeSorter<P> >
class MirrorMergeSorter: public ParallelMirrorSorterBase<T,U,P,S>  {
};

#endif /* parallel_mergesort_h */
