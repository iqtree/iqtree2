#include "iqtree2.h"
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

    // =============================================================
    // for testing the available functions declared inside iqtree2.h
    // =============================================================
        
        if (argc < 2) {
            // show the available options for testing piqtree
            std::cout << "available options:" << endl;
            std::cout << "calculate RF distance: " << argv[0] << " 1 [tree1 newick string] [tree2 newick string]" << endl;
            std::cout << "generate random tree: " << argv[0] << " 2 [num taxa] [seed] [mode (e.g. 'YULE_HARDING')] [out filename]" << endl;
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
            if (argc < 6) {
                std::cout << "Syntax: " << argv[0] << " 2 [num taxa] [seed] [mode (e.g. 'YULE_HARDING')] [out filename]" << endl;
                return 0;
            }
            int numtaxa = atoi(argv[2]);
            int seed = atoi(argv[3]);
            string gen_mode = argv[4];
            string outfilename = argv[5];
            generate_random_tree_file(numtaxa, seed, gen_mode, outfilename);
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
