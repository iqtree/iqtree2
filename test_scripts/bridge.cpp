#include<vector>
#include<string>
#include<iostream>

/*
 * Calculates the robinson fould distance between two trees
 */
extern "C" int robinson_fould(const char* tree1, const char* tree2);

/*
 * Generates a set of random phylogenetic trees
 * tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
 * output: a newick tree (in string format)
 */
std::string random_tree(int num_taxa, std::string tree_gen_mode, int num_trees, int rand_seed = 0);

/*
 * Perform phylogenetic analysis on the input alignment
 * With estimation of the best topology
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the tree and the details of parameters
 */
std::string build_tree(std::vector<std::string>& names, std::vector<std::string>& seqs, std::string model, int rand_seed = 0, int bootstrap_rep = 0, int num_thres = 1);

/*
 * Perform phylogenetic analysis on the input alignment
 * With restriction to the input toplogy
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the details of parameters
 */
std::string fit_tree(std::vector<std::string>& names, std::vector<std::string>& seqs, std::string model, std::string intree, int rand_seed = 0, int num_thres = 1);

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
std::string modelfinder(std::vector<std::string>& names, std::vector<std::string>& seqs, int rand_seed = 0,
                   std::string model_set = "", std::string freq_set = "", std::string rate_set = "", int num_thres = 1);

/*
 * Build pairwise JC distance matrix
 * output: set of distances
 * (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
 * where n is the number of sequences
 * num_thres -- number of cpu threads to be used, default: 1; 0 - use all available cpu threads on the machine
 */
std::vector<double> build_distmatrix(std::vector<std::string>& names, std::vector<std::string>& seqs, int num_thres = 1);

/*
 * Using Rapid-NJ to build tree from a distance matrix
 * output: a newick tree (in string format)
 */
std::string build_njtree(std::vector<std::string>& names, std::vector<double>& distances);

/*
 * verion number
 */
extern "C" const char* version();

int main() {
  const char* tree1 = "(a,b,(c,(d,e)));";
  const char* tree2 = "(e,b,(c,(d,a)));";

  std::cout << robinson_fould(tree1, tree2) << std::endl;
  return 0;
}
