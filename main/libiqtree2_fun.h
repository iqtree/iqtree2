#ifndef LIBIQTREE2_FUN
#define LIBIQTREE2_FUN

#include "tree/mtree.h"
#include "tree/phylotree.h"
#include "utils/tools.h"
#include "main/alisim.h"
#include "suppFunc.h"
#include <vector>
#include <string>

using namespace std;

// Calculates the robinson fould distance between two trees
int robinson_fould(const string& tree1, const string& tree2);

// Generates a set of random phylogenetic trees
// tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
string random_tree(int num_taxa, string tree_gen_mode, int num_trees, int rand_seed = 0);

// Perform phylogenetic analysis on the input alignment (in string format)
// With estimation of the best topology
string build_tree(vector<string>& names, vector<string>& seqs, string model, int rand_seed = 0);

// Perform phylogenetic analysis on the input alignment (in string format)
// With restriction to the input toplogy
string fit_tree(vector<string>& names, vector<string>& seqs, string model, string intree, int rand_seed = 0);

#endif /* LIBIQTREE2_FUN */
