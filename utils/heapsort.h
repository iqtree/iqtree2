//
//  heapsort.h
//  iqtree
//
//  Created by James Barbetti on 22/6/20.
//

#ifndef heapsort_h
#define heapsort_h

template <class V, class S>
void constructMirroredHeap ( V* valueArray
                            , int start, int stop
                            , S* satteliteArray) {
    //
    //Constructs a "top-less" radix 2 max heap
    //of valueArray[start..stop] and "mirrors"
    //the permutation of the values
    //on satteliteArray[start..stop].
    //
    int fudge = 2 - start;
    for (int h = start+(stop-start)/2; h >= start; --h) {
        int i = h;
        V   v = valueArray[i];
        S   s = satteliteArray[i];
        int j = i + i + fudge;
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
                              , int start, int stop
                              , S* satteliteArray) {
    //
    //Extracts a value (and some sattelite information)
    //from a radix 2 max heap.
    //
    int fudge = 2 - start;
    for (--stop; start<=stop; --stop) {
        int i = stop;
        V   v = valueArray[i];
        S   s = satteliteArray[i];
        int j = start;
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
                       , int start, int stop
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
