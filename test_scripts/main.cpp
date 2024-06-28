#include "iqtree2.h"
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {

    // ====================================================
    // for testing the available functions for PiQTree
    // ====================================================
        
        if (argc < 2) {
            // show the available options for testing piqtree
            cout << "available options:" << endl;
            cout << "calculate RF distance: " << argv[0] << " 1 [tree1 newick string] [tree2 newick string]" << endl;
            cout << "generate random tree: " << argv[0] << " 2 [num taxa] [tree mode (e.g. 'YULE_HARDING')] [num trees] [optional: rand_seed]" << endl;
            cout << "phylogenetic analysis: " << argv[0] << " 3 [alignment file] [model] [optional: rand_seed]" << endl;
            cout << "fitting-tree analysis: " << argv[0] << " 4 [alignment file] [model] [tree newick string] [optional: rand_seed]" << endl;
            return 0;
        }

        if (atoi(argv[1]) == 1) {
            // calculate the RF distance
            if (argc < 4) {
                cout << "Syntax: " << argv[0] << " 1 [tree1 newick string] [tree2 newick string]" << endl;
                return 0;
            }
            string tree1 = argv[2];
            string tree2 = argv[3];
            cout << "RF dist: " << robinson_fould(tree1, tree2) << endl;
        }
        
        if (atoi(argv[1]) == 2) {
            // generate a random tree
            if (argc < 5) {
                cout << "Syntax: " << argv[0] << " 2 [num taxa] [mode (e.g. 'YULE_HARDING')] [num trees] [optional: rand_seed]" << endl;
                return 0;
            }
            int num_taxa = atoi(argv[2]);
            string tree_mode = argv[3];
            int num_tree = atoi(argv[4]);
            int seed = 0;
            if (argc > 5)
                seed = atoi(argv[5]);
            cout << random_tree(num_taxa, tree_mode, num_tree, seed) << endl;
        }

        if (atoi(argv[1]) == 3) {
            // perform phylogenetic analysis
            if (argc < 4) {
                cout << "Syntax: " << argv[0] << " 3 [alignment file] [model] [optional: rand_seed]" << endl;
                return 0;
            }
            string align_file = argv[2];
            string model = argv[3];
            int seed = 0;
            if (argc > 4)
                seed = atoi(argv[4]);

            // create the arrays for names and sequences
            ifstream fin;
            vector<string> names;
            vector<string> seqs;
            fin.open(align_file);
            string aline;
            while (getline(fin,aline)) {
               if (aline.length() == 0)
                   continue;
               if (aline[0] == '>') {
                   // sequence name
                   names.push_back(aline.substr(1));
                   seqs.push_back("");
               } else if (seqs.size() > 0) {
                   seqs[seqs.size() - 1].append(aline);
               }
            }
            fin.close();

            if (names.size() != seqs.size()) {
               cerr << "Error! Number of items inside names does not equal the number of item in seqs" << endl;
               exit(1);
            } 
            string output_str = build_tree(names, seqs, model, seed);
            cout << output_str << endl;
        }

        if (atoi(argv[1]) == 4) {
            // fitting a tree for phylogenetic analysis
            if (argc < 5) {
                cout << "Syntax: " << argv[0] << " 3 [alignment file] [model] [tree newick] [optional: rand_seed]" << endl;
                return 0;
            }
            string align_file = argv[2];
            string model = argv[3];
            string tree = argv[4];
            int seed = 0;
            if (argc > 5)
                seed = atoi(argv[5]);

            // create the arrays for names and sequences
            vector<string> names;
            vector<string> seqs;
            ifstream fin;
            fin.open(align_file);
            string aline;
            while (getline(fin,aline)) {
               if (aline.length() == 0)
                   continue;
               if (aline[0] == '>') {
                   // sequence name
                   names.push_back(aline.substr(1));
                   seqs.push_back("");
               } else if (seqs.size() > 0) {
                   seqs[seqs.size() - 1].append(aline);
               }
            }
            fin.close();

            if (names.size() != seqs.size()) {
               cerr << "Error! Number of items inside names does not equal the number of item in seqs" << endl;
               exit(1);
            } 
            string output_str = fit_tree(names, seqs, model, tree, seed);
            cout << output_str << endl;
        }
        
        return 0;

}
