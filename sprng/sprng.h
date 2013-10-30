#ifndef _sprng_h_
#define _sprng_h_

#include "interface.h"

#define SPRNG_DEFAULT 0
#define CRAYLCG 0
#define DRAND48 1
#define FISH1   2
#define FISH2   3
#define FISH3   4
#define FISH4   5
#define FISH5   6
#define LECU1   0
#define LECU2   1
#define LECU3   2
#define LAG1279  0
#define LAG17    1
#define LAG31    2
#define LAG55    3
#define LAG63    4
#define LAG127   5
#define LAG521   6
#define LAG521B  7
#define LAG607   8
#define LAG607B  9
#define LAG1279B 10

#define CHECK 1

#define MAX_PACKED_LENGTH 24000

#ifdef SPRNG_USE_MPI
#define MPINAME(A) A ## _mpi
#else
#define MPINAME(A) A
#endif

#define make_sprng_seed MPINAME(make_new_seed)

#if defined(SIMPLE_SPRNG)

/* HAS ;-) */
#if 0
#define pack_sprng pack_rng_simple
#define unpack_sprng unpack_rng_simple
#define isprng  MPINAME(get_rn_int_simple)
#define init_sprng MPINAME(init_rng_simple)
#define print_sprng print_rng_simple

#ifdef FLOAT_GEN
#define sprng  MPINAME(get_rn_flt_simple)
#else
#define sprng  MPINAME(get_rn_dbl_simple)
#endif
#endif

#elif !defined(CHECK_POINTERS)

#define free_sprng free_rng
#define pack_sprng pack_rng
#define unpack_sprng unpack_rng
#define isprng  get_rn_int
#define spawn_sprng(A,B,C) spawn_rng(A,B,C,!CHECK)
#define init_sprng init_rng
#define print_sprng print_rng

#ifdef FLOAT_GEN
#define sprng  get_rn_flt
#else
#define sprng  get_rn_dbl
#endif

#else

/* HAS ;-) */
#if 0
#define free_sprng(A) ((deleteID(A)==NULL) ? -1 : free_rng(A))
#define pack_sprng(A,B) ((checkID(A)==NULL) ? 0 : pack_rng(A,B))
#define unpack_sprng(A) addID(unpack_rng(A))
#define isprng(A)  ((checkID(A)==NULL) ? -1 : get_rn_int(A))
#define spawn_sprng(A,B,C) ((checkID(A)==NULL) ? 0 : spawn_rng(A,B,C,CHECK))
#define init_sprng(A,B,C,D) addID(init_rng(A,B,C,D))
#define print_sprng(A) ((checkID(A)==NULL) ? 0 : print_rng(A))

#ifdef FLOAT_GEN
#define sprng(A)  ((checkID(A)==NULL) ? -1.0 : get_rn_flt(A))
#else
#define sprng(A)  ((checkID(A)==NULL) ? -1.0 : get_rn_dbl(A))
#endif
#endif

#endif

#endif
