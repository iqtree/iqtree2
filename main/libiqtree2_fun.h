#ifndef LIBIQTREE2_FUN
#define LIBIQTREE2_FUN

#include "tree/mtree.h"
#include "tree/phylotree.h"
#include "utils/tools.h"
#include "main/alisim.h"
#include "utils/starttree.h"
#include "suppFunc.h"
#include <vector>
#include <string>
#include <cstring>

using namespace std;

#ifdef _MSC_VER
#pragma pack(push, 1)
#else
#pragma pack(1)
#endif

typedef struct {
  const char** strings;
  size_t length;
} StringArray;

typedef struct {
  double* doubles;
  size_t length;
} DoubleArray;

typedef struct {
  int value;
  char* errorStr;
} IntegerResult;

typedef struct {
  char* value;
  char* errorStr;
} StringResult;

typedef struct {
  double* value;
  size_t length;
  char* errorStr;
} DoubleArrayResult;

#ifdef _MSC_VER
#pragma pack(pop)
#else
#pragma pack()
#endif

/*
 * Calculates the robinson fould distance between two trees
 */
extern "C" IntegerResult robinson_fould(const char* ctree1, const char* ctree2);

/*
 * Generates a set of random phylogenetic trees
 * tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
 * output: a newick tree (in string format)
 */
extern "C" StringResult random_tree(int num_taxa, const char* tree_gen_mode, int num_trees, int rand_seed = 0);

/*
 * Perform phylogenetic analysis on the input alignment
 * With estimation of the best topology
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the tree and the details of parameters
 */
extern "C" StringResult build_tree(StringArray& names, StringArray& seqs, const char* model, int rand_seed = 0, int bootstrap_rep = 0, int num_thres = 1);

/*
 * Perform phylogenetic analysis on the input alignment
 * With restriction to the input toplogy
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the details of parameters
 */
extern "C" StringResult fit_tree(StringArray& names, StringArray& seqs, const char* model, const char* intree, int rand_seed = 0, int num_thres = 1);

/*
 * Perform phylogenetic analysis with ModelFinder
 * on the input alignment (in string format)
 * model_set -- a set of models to consider
 * freq_set -- a set of frequency types
 * rate_set -- a set of RHAS models
 * rand_seed -- random seed, if 0, then will generate a new random seed
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: modelfinder results in YAML format
 */
extern "C" StringResult modelfinder(StringArray& names, StringArray& seqs, int rand_seed = 0,
                   const char* model_set = "", const char* freq_set = "", const char* rate_set = "", int num_thres = 1);

/*
 * Build pairwise JC distance matrix
 * output: set of distances
 * (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
 * where n is the number of sequences
 * num_thres -- number of cpu threads to be used, default: 1; 0 - use all available cpu threads on the machine
 */
extern "C" DoubleArrayResult build_distmatrix(StringArray& names, StringArray& seqs, int num_thres = 1);

/*
 * Using Rapid-NJ to build tree from a distance matrix
 * output: a newick tree (in string format)
 */
extern "C" StringResult build_njtree(StringArray& names, DoubleArray& distances);

/*
 * verion number
 */
extern "C" StringResult version();

#endif /* LIBIQTREE2_FUN */
