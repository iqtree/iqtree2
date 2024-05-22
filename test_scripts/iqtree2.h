#ifndef LIBIQTREE2_FUN
#define LIBIQTREE2_FUN

#include <string>

using namespace std;

// Calculates the RF distance between two trees
int calculate_RF_distance(const string& tree1, const string& tree2);

// Generates a random phylogenetic tree
void generate_random_tree_file(int numtaxa, int seed, string tree_gen_mode, string outfile);

// perform phylogenetic analysis on the input alignment file
void phylogenetic_analysis(string& align_file, int ncpus = 1);

#endif /* LIBIQTREE2_FUN */
