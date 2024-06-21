#include "iqtree2.h"
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

    // ====================================================
    // for testing the available functions for PiQTree
    // ====================================================
        
        if (argc < 2) {
            // show the available options for testing piqtree
            std::cout << "available options:" << endl;
            std::cout << "calculate RF distance: " << argv[0] << " 1 [tree1 newick string] [tree2 newick string]" << endl;
            std::cout << "generate random tree: " << argv[0] << " 2 [num taxa] [tree mode (e.g. 'YULE_HARDING')] [num trees] [optional: rand_seed]" << endl;
            std::cout << "phylogenetic analysis: " << argv[0] << " 3 [alignment file]" << endl;
            return 0;
        }

        if (atoi(argv[1]) == 1) {
            // calculate the RF distance
            if (argc < 4) {
                std::cout << "Syntax: " << argv[0] << " 1 [tree1 newick string] [tree2 newick string]" << endl;
                return 0;
            }
            std::string tree1 = argv[2];
            std::string tree2 = argv[3];
            std::cout << "RF dist: " << calculate_RF_distance(tree1, tree2) << std::endl;
            return 0;
        }
        
        if (atoi(argv[1]) == 2) {
            // generate a random tree
            if (argc < 5) {
                std::cout << "Syntax: " << argv[0] << " 2 [num taxa] [mode (e.g. 'YULE_HARDING')] [num trees] [optional: rand_seed]" << endl;
                return 0;
            }
            int num_taxa = atoi(argv[2]);
            string tree_mode = argv[3];
            int num_tree = atoi(argv[4]);
            int seed = 0;
            if (argc > 5)
                seed = atoi(argv[5]);
            std::cout << generate_random_tree(num_taxa, tree_mode, num_tree, seed) << endl;
            return 0;
        }

        if (atoi(argv[1]) == 3) {
            // perform phylogenetic analysis
            if (argc < 3) {
                std::cout << "Syntax: " << argv[0] << " 3 [alignment file] [# cpu threads]" << endl;
                return 0;
            }
            string align_file = argv[2];
            if (argc < 4)
                phylogenetic_analysis(align_file);
            else {
                int ncpus = atoi(argv[3]);
                phylogenetic_analysis(align_file, ncpus);
            }
            return 0;
        }
        
        return 0;

}
