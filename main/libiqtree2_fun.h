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

using namespace std;

/*
 * Calculates the robinson fould distance between two trees
 */
int robinson_fould(const string& tree1, const string& tree2);

/*
 * Generates a set of random phylogenetic trees
 * tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
 * output: a newick tree (in string format)
 */
string random_tree(int num_taxa, string tree_gen_mode, int num_trees, int rand_seed = 0);

/*
 * Perform phylogenetic analysis on the input alignment
 * With estimation of the best topology
 * num_thres -- number of cpu threads to be used, default: 1
 * output: results in YAML format with the tree and the details of parameters
 */
string build_tree(vector<string>& names, vector<string>& seqs, string model, int rand_seed = 0, int bootstrap_rep = 0, int num_thres = 1);

/*
 * Perform phylogenetic analysis on the input alignment
 * With restriction to the input toplogy
 * num_thres -- number of cpu threads to be used, default: 1
 * output: results in YAML format with the details of parameters
 */
string fit_tree(vector<string>& names, vector<string>& seqs, string model, string intree, int rand_seed = 0, int num_thres = 1);

/*
 * Perform phylogenetic analysis with ModelFinder
 * on the input alignment (in string format)
 * model_set -- a set of models to consider
 * freq_set -- a set of frequency types
 * rate_set -- a set of RHAS models
 * rand_seed -- random seed, if 0, then will generate a new random seed
 * num_thres -- number of cpu threads to be used, default: 1
 * output: modelfinder results in YAML format
 */
string modelfinder(vector<string>& names, vector<string>& seqs, int rand_seed = 0,
                   string model_set = "", string freq_set = "", string rate_set = "", int num_thres = 1);

/*
 * Build pairwise JC distance matrix
 * output: set of distances
 * (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
 * where n is the number of sequences
 * num_thres -- number of cpu threads to be used, default: 1
 */
vector<double> build_distmatrix(vector<string>& names, vector<string>& seqs, int num_thres = 1);

/*
 * Using Rapid-NJ to build tree from a distance matrix
 * output: a newick tree (in string format)
 */
string build_njtree(vector<string>& names, vector<double>& distances);

#endif /* LIBIQTREE2_FUN */
