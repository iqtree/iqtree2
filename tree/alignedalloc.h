//
//  alignedalloc.h
//  iqtree
//
//  Moved to this file (from phylotree.h) by James Barbetti on 9/10/20.
//  phylotree.h due to: BUI Quang Minh <minh.bui@univie.ac.at>
//                      Lam-Tung Nguyen <nltung@gmail.com>
//
//  (The functions in this file were moved here, so they are visible
//   to member functions implemented in likelihoodbufferset.cpp).
//

#ifndef alignedalloc_h
#define alignedalloc_h

#include <utils/tools.h>    //for Params::getInstance()

#ifndef ROUND_UP_TO_MULTIPLE
#define ROUND_UP_TO_MULTIPLE
template<class C, class G> C roundUpToMultiple(C count, G grain) {
    //Assumed: division rounds down, for type C.
    return ((count+grain-1)/grain)*grain;
}
#endif

template< class T>
inline T *aligned_alloc(size_t size) {
    size_t MEM_ALIGNMENT = (Params::getInstance().SSE >= LK_AVX512) ? 64 : ((Params::getInstance().SSE >= LK_AVX) ? 32 : 16);
    void *mem;

#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined(WIN64)
    #if (defined(__MINGW32__) || defined(__clang__)) && defined(BINARY32)
        mem = __mingw_aligned_malloc(size*sizeof(T), MEM_ALIGNMENT);
    #else
        mem = _aligned_malloc(size*sizeof(T), MEM_ALIGNMENT);
    #endif
#else
    int res = posix_memalign(&mem, MEM_ALIGNMENT, size*sizeof(T));
    if (res == ENOMEM) {
#if (defined(__GNUC__) || defined(__clang__)) && !defined(WIN32) &&!defined(WIN64) && !defined(__CYGWIN__)
        print_stacktrace(cerr);
#endif
        outError("Not enough memory, allocation of " + convertInt64ToString(size*sizeof(T)) + " bytes failed (bad_alloc)");
    }
#endif
    if (mem == NULL) {
#if (defined(__GNUC__) || defined(__clang__)) && !defined(WIN32) &&!defined(WIN64) && !defined(__CYGWIN__)
        print_stacktrace(cerr);
#endif
        outError("Not enough memory, allocation of " + convertInt64ToString(size*sizeof(T)) + " bytes failed (bad_alloc)");
    }
    return (T*)mem;
}

template <class T> T* ensure_aligned_allocated(T* & ptr, size_t size) {
    if (ptr == nullptr) {
        ptr = aligned_alloc<T>(size);
    }
    return ptr;
}

 template <class T> void aligned_free(T* & mem) {
     if (mem == nullptr) {
         return;
     }
#if defined WIN32 || defined _WIN32 || defined __WIN32__
    #if (defined(__MINGW32__) || defined(__clang__)) && defined(BINARY32)
        __mingw_aligned_free(mem);
    #else
        _aligned_free(mem);
    #endif
#else
    free(mem);
#endif
    mem = nullptr;
}

#endif /* alignedalloc_h */
