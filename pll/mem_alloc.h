#ifndef __mem_alloc_h
#define __mem_alloc_h
#include <stddef.h>
#include <stdlib.h>
#ifdef __linux__
#include <malloc.h>
#endif
#include "pll.h"

#define rax_memalign memalign
#define rax_malloc malloc
#define rax_realloc realloc
#define rax_free free
#define rax_posix_memalign posix_memalign
#define rax_calloc calloc
//#define rax_malloc_aligned(x) memalign(PLL_BYTE_ALIGNMENT,x)

//void *rax_memalign(size_t align, size_t size);
//void *rax_malloc(size_t size);
//void *rax_realloc(void *p, size_t size);
//void rax_free(void *p);
//int rax_posix_memalign(void **p, size_t align, size_t size);
//void *rax_calloc(size_t n, size_t size);
//
//void *rax_malloc_aligned(size_t size);


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
