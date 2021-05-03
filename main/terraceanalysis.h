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

#endif /* terraceanalysis_hpp */
