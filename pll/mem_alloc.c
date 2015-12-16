
#define MEM_ALLOC_NO_GUARDS 1

#include "mem_alloc.h"

#include <stdio.h>
#include <stdlib.h>
#ifndef __APPLE__
#include <malloc.h>             // this is probably not necessary
#endif

#ifdef RAXML_USE_LLALLOC

// the llalloc library implementation in lockless_alloc/ll_alloc.c exports the alloction functions prefixed
// with 'llalloc'. The following are the forward declarations of the llalloc* functions 

#define PREFIX(X)   llalloc##X

void *PREFIX(memalign)(size_t align, size_t size);
void *PREFIX(malloc)(size_t size);
void *PREFIX(realloc)(void *p, size_t size);
int PREFIX(posix_memalign)(void **p, size_t align, size_t size);
void *PREFIX(calloc)(size_t n, size_t size);
void PREFIX(free)(void *p);


// wrappers that forward the rax_* functions to the corresponding llalloc* functions


void *rax_memalign(size_t align, size_t size) {
  return PREFIX(memalign)(align, size);
}

void *rax_malloc( size_t size ) {
  return PREFIX(malloc)(size);
}
void *rax_realloc( void *p, size_t size ) {
  return PREFIX(realloc)(p, size);
}


void rax_free(void *p) {
  PREFIX(free)(p);
}

int rax_posix_memalign(void **p, size_t align, size_t size) {
  return PREFIX(posix_memalign)(p, align, size);
}
void *rax_calloc(size_t n, size_t size) {
  return PREFIX(calloc)(n,size);
}

void *rax_malloc_aligned(size_t size) 
{
  const size_t PLL_BYTE_ALIGNMENT = 32;
  return rax_memalign(PLL_BYTE_ALIGNMENT, size);
  
}

#else // RAXML_USE_LLALLOC
// if llalloc should not be used, forward the rax_* functions to the corresponding standard function

void *rax_memalign(size_t align, size_t size) {
#if defined (__APPLE__)
    void * mem;
    if (posix_memalign (&mem, align, size))
      return (NULL);
    else
      return (mem);
#else
    return memalign(align, size);
#endif
    
}

void *rax_malloc( size_t size ) {
  return malloc(size);
}
void *rax_realloc( void *p, size_t size ) {
  return realloc(p, size);
}


void rax_free(void *p) {
  free(p);
}

int rax_posix_memalign(void **p, size_t align, size_t size) {
  return posix_memalign(p, align, size);
}
void *rax_calloc(size_t n, size_t size) {
  return calloc(n,size);
}

void *rax_malloc_aligned(size_t size) 
{
  const size_t PLL_BYTE_ALIGNMENT = 32;
  return rax_memalign(PLL_BYTE_ALIGNMENT, size);
  
}

#endif



#if 0
//
// two test cases to check if the default malloc plays along well with lockless malloc. Normally there shoudl not be a
// problem as long as everyone handles 'foreign' sbrk calls gracefully (as lockless and glibc seem to do). 
// WARNING: there is a slightly worrying comment in glibc malloc, which seems to assume that magically no foreign sbrks
// happen between two consecutive sbrk calls while re-establishing page alignment in some obscure special case. IMHO, this
// is clearly an error (race) in multithreaded programs, as there is no way how a foreign sbrk user can properly lock anything.
// see: http://sourceware.org/git/?p=glibc.git;a=blob;f=malloc/malloc.c;h=0f1796c9134ffef289ec31fb1cd538f3a9490ae1;hb=HEAD#l2581
//
// If all threads consistently only use the rax_* wrappers this is not a problem, but as this is a library, we can not be sure 
// that no other thread uses default malloc... note that lockless malloc only uses sbrk for the slab (=small block) area, while 
// raxml heavy uses malloc/free only on much larger blocks...
// If anything ever goes wrong while using mixed glibc/lockless malloc, this should be investigated.
//
// TODO: the potential race seems to be related to handling the case where a 'foreign sbrk' adjusted the break to a non page-boundary.
// check if lockless malloc actually ever adjusts to non page-boundaries.


void check_block( void *p, size_t size ) {
    size_t i;
    char *cp = (char*)p;
    
    for( i = 0; i < size; ++i ) {
        
        if( cp[i] != (char)i ) {
            printf( "MEEEEEEEEEEEEEEEEEEEEP\n" );
            abort();
        }
    }
    
}


void fill_block( void *p, size_t size ) {
    size_t i;
    char *cp = (char*)p;
    
    for( i = 0; i < size; ++i ) {
        cp[i] = (char)i;
    }
}


void malloc_stress() {
    const int n_slots = 100000;
    
    void *blocks1[n_slots];
    size_t sizes1[n_slots];
    void *blocks2[n_slots];
    size_t sizes2[n_slots];
    
    memset( blocks1, 0, sizeof( void * ) * n_slots ); 
    memset( blocks2, 0, sizeof( void * ) * n_slots ); 
    
    memset( sizes1, 0, sizeof( size_t ) * n_slots );
    memset( sizes2, 0, sizeof( size_t ) * n_slots );
    
    
    
    while( 1 ) {
        int r = rand() % n_slots;
        
        void *bs;
        
        
        int size;
        if( rand() % 2 == 0 ) {
            size = rand() % (32 * 16); // hit slab
        } else {
            size = (rand() % 128) * 128; // not slab
        }
            
            
        if( 1 || rand() % 2 == 0 ) {
            if( blocks1[r] == 0 ) {
                blocks1[r] = malloc( size );
                sizes1[r] = size;
                fill_block( blocks1[r], sizes1[r] );
            } else {
                check_block( blocks1[r], sizes1[r] );
                free( blocks1[r] );
                blocks1[r] = 0;
            }
        } else {
            if( blocks2[r] == 0 ) {
                blocks2[r] = rax_malloc( size );
                sizes2[r] = size;
                fill_block( blocks2[r], sizes2[r] );
            } else {
                check_block( blocks2[r], sizes2[r] );
                
                rax_free( blocks2[r] );
                blocks2[r] = 0;
            }
        }
            
       
        
    }
    
}


void malloc_stress2() {
    const size_t n_slots = 1000;
    
    void *blocks[n_slots];
    size_t i;
    for( i = 0; i < n_slots; ++i ) {
        blocks[i] = malloc( (rand() % 32) * 1024 ); 
        
    }
    sbrk( 10 );
    for( i = 0; i < n_slots; ++i ) {
        free(blocks[i]);
        
    }
    
    
    
}
#endif

