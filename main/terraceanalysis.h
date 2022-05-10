/*
 *  terraceanalysis.ch
 *  Interface to work with terraces
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#ifndef terraceanalysis_hpp
#define terraceanalysis_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "terrace/terrace.hpp"

/*
 *  Main function for terrace analysis
 */
void runterraceanalysis(Params &params);

/*
 * Run terrace check:
 * - check if query trees lie on the same terrace with a representative tree. Naive pairwise comparison of induced partition trees.
 */
void run_terrace_check(Terrace *terrace, Params &params);


/**
        This is the function that calls generation of terrace trees, when the input data (representative tree and presence-absence matrix are processed already)
 */
void run_generate_trees(Terrace *terrace, Params &params,const int m);

/**
        The function is used to read a set of subtrees to be considered as partition trees for terrace analysis
 */
void read_tree_set(const char *infile, bool &is_rooted, vector<TerraceTree*> &subtrees);

#endif /* terraceanalysis_hpp */
