
#ifndef _interface_h_
#define _interface_h_

#ifndef ANSI_ARGS
#ifdef __STDC__
#define ANSI_ARGS(args) args
#else
#define ANSI_ARGS(args) ()
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

int get_rn_int ANSI_ARGS((int *igenptr));
float get_rn_flt ANSI_ARGS((int *igenptr));
double get_rn_dbl ANSI_ARGS((int *igenptr));
int *init_rng ANSI_ARGS(( int gennum, int total_gen,  int seed,
			  int mult));
int spawn_rng ANSI_ARGS((int *igenptr, int nspawned, int ***newgens, int checkid) );
int make_new_seed ANSI_ARGS((void));
int make_new_seed_mpi ANSI_ARGS((void));
int get_seed__rng ANSI_ARGS((int *genptr));
int free_rng ANSI_ARGS((int *genptr));
int pack_rng ANSI_ARGS(( int *genptr, char **buffer));
int *unpack_rng ANSI_ARGS(( char *packed));
int print_rng ANSI_ARGS(( int *igen));
int *checkID ANSI_ARGS(( int *igen));
int *addID ANSI_ARGS(( int *igen));
int *deleteID ANSI_ARGS(( int *igen));


/* HAS ;-) */
#if 0
int *init_rng_simple ANSI_ARGS(( int seed,  int mult));
int *init_rng_simple_mpi ANSI_ARGS(( int seed,  int mult));
int get_rn_int_simple ANSI_ARGS((void));
int get_rn_int_simple_mpi ANSI_ARGS((void));
float get_rn_flt_simple ANSI_ARGS((void));
float get_rn_flt_simple_mpi ANSI_ARGS((void));
double get_rn_dbl_simple ANSI_ARGS((void));
double get_rn_dbl_simple_mpi ANSI_ARGS((void));
int pack_rng_simple ANSI_ARGS((char **buffer));
int *unpack_rng_simple ANSI_ARGS(( char *packed));
int print_rng_simple ANSI_ARGS((void));
#endif


#ifdef __cplusplus
}
#endif


#endif
