//
//  parallel_sort.h
//  iqtree
//
//  Created by James Barbetti on 17/3/21.
//

#ifndef parallel_sort_h
#define parallel_sort_h

#include <stdint.h> //for definition of intptr_t
#include <vector>   //for definition of std::vector template class

template <class T> class ParallelSorter {
public:
    virtual void parallel_sort     (T* data, intptr_t count) = 0;
    virtual void single_thread_sort(T* data, intptr_t count) = 0;
    virtual ~ParallelSorter() = default;
};

template <class T, class U> class ParallelMirrorSorter {
public:
    virtual void parallel_mirror_sort     (T* data, intptr_t count, U* auxiliary) = 0;
    virtual void single_thread_mirror_sort(T* data, intptr_t count, U* auxiliary) = 0;
    virtual ~ParallelMirrorSorter() = default;
};

template <class T, class U> class ValueAndSattelite {
public:
    T t;
    U u;
    ValueAndSattelite() = default;
    ValueAndSattelite(const ValueAndSattelite& rhs) = default;
    ValueAndSattelite(const T& first, const U& second) : t(first), u(second) {}
    ValueAndSattelite& operator=(const ValueAndSattelite& rhs) = default;
    bool operator <  ( const ValueAndSattelite& rhs ) const { return t <  rhs.t; }
    bool operator <= ( const ValueAndSattelite& rhs ) const { return t <= rhs.t; }
};

template <class T,class U, class P=ValueAndSattelite<T,U>,
          class S=ParallelSorter<P> >
class ParallelMirrorSorterBase: public S, public ParallelMirrorSorter<T,U> {
protected:
    std::vector<P> copy;
    void loadCopy(T* data, intptr_t count, U* auxiliary) {
        copy.resize(count);
        for (int i=0;i<count;++i) {
            copy[i]=P(data[i], auxiliary[i]);
        }
    }
    void unloadCopy(T* data, intptr_t count, U* auxiliary) {
        P* copy_ptr = copy.data();
        for (int i=0;i<count;++i) {
            data[i]      = copy_ptr[i].t;
            auxiliary[i] = copy_ptr[i].u;
        }
    }
public:
    typedef S super;
    using super::parallel_sort;
    using super::single_thread_sort;
    virtual void parallel_mirror_sort(T* data, intptr_t count, U* auxiliary) override {
        loadCopy(data, count, auxiliary);
        super::parallel_sort(copy.data(), count);
        unloadCopy(data, count, auxiliary);
    }
    virtual void single_thread_mirror_sort(T* data, intptr_t count, U* auxiliary) override {
        loadCopy(data, count, auxiliary);
        super::single_thread_sort(copy.data(), count);
        unloadCopy(data, count, auxiliary);
    }
};

#endif /* parallel_sort_h */
