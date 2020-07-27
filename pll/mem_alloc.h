#ifndef __mem_alloc_h
#define __mem_alloc_h

#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
#include <stdlib.h>
#include <malloc.h>
#endif

#include <stddef.h>
#include <stdlib.h>
#ifdef __linux__
#include <malloc.h>
#endif
#include "pll.h"
#include <string.h>
#include <iqtree_config.h>

//#define rax_memalign memalign
//#define rax_malloc malloc
//#define rax_calloc calloc
//#define rax_realloc realloc


#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined(WIN64)
    #if (defined(__MINGW32__) || defined(__clang__)) && defined(BINARY32)
        #define rax_posix_memalign(ptr,alignment,size) *(ptr) = __mingw_aligned_malloc((size),(alignment))
        #define rax_malloc(size) __mingw_aligned_malloc((size), PLL_BYTE_ALIGNMENT)
        void *rax_calloc(size_t count, size_t size);
        #define rax_free __mingw_aligned_free
    #else
        #define rax_posix_memalign(ptr,alignment,size) *(ptr) = _aligned_malloc((size),(alignment))
        #define rax_malloc(size) _aligned_malloc((size), PLL_BYTE_ALIGNMENT)
        void *rax_calloc(size_t count, size_t size);
        #define rax_free _aligned_free
    #endif
#else
    #define rax_posix_memalign posix_memalign
    #define rax_malloc malloc
    #define rax_calloc calloc
    #define rax_free free
#endif

//#define rax_malloc_aligned(x) memalign(PLL_BYTE_ALIGNMENT,x)

//void *rax_memalign(size_t align, size_t size);
//void *rax_malloc(size_t size);
//void *rax_realloc(void *p, size_t size);
//void rax_free(void *p);
//int rax_posix_memalign(void **p, size_t align, size_t size);
//void *rax_calloc(size_t n, size_t size);
//
//void *rax_malloc_aligned(size_t size);


/* for strndup stuff */
char *my_strndup(const char *s, size_t n);

#ifndef HAVE_STRTOK_R
char *strtok_r (char * s, const char * delim, char **save_ptr);
#endif

#if 0
// using the following contraption to trigger a compile-time error does not work on some gcc versions. It will trigger a confising linker error in the best case, so it is deativated.

#if defined(RAXML_USE_LLALLOC) && !defined(MEM_ALLOC_NO_GUARDS)
#define malloc(x) XXX_DONT_USE_MALLOC_WITHOUT_RAX_PREFIX_XXX
#define free(x) XXX_DONT_USE_FREE_WITHOUT_RAX_PREFIX_XXX
#define calloc(x,y) XXX_DONT_USE_CALLOC_WITHOUT_RAX_PREFIX_XXX
#define realloc(x,y) XXX_DONT_USE_REALLOC_WITHOUT_RAX_PREFIX_XXX
#define malloc_aligned(x) XXX_DONT_USE_MALLOC_ALIGNED_WITHOUT_RAX_PREFIX_XXX
#define posix_memalign(x,y,z) XXX_DONT_USE_POSIX_MEMALIGN_ALIGNED_WITHOUT_RAX_PREFIX_XXX
#endif
#endif

#endif
