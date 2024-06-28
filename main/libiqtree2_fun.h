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

// Calculates the RF distance between two trees
int calculate_RF_distance(const string& tree1, const string& tree2);

// Generates a random phylogenetic tree
void generate_random_tree_file(int numtaxa, int seed, string tree_gen_mode, string outfile);

// perform phylogenetic analysis on the input alignment file
void phylogenetic_analysis(string& align_file, int ncpus = 1);

// -----------------------------
// Update functions
// -----------------------------

// Calculates the RF distances betweeen a set of trees
// Input: a string containing a set of trees, separated by ";"
// string robinson_foulds(string& trees);

// Generates a random phylogenetic tree
string generate_random_tree(int num_taxa, string tree_gen_mode, int num_trees, int rand_seed = 0);

// Perform phylogenetic analysis on the input alignment (in string format)
string build_tree(vector<string>& names, vector<string>& seqs, string model, int rand_seed = 0);

#endif /* LIBIQTREE2_FUN */
