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


typedef struct Task{
    
    void (*taskFunction)(Terrace*, 
                        vector<Terrace*>&, 
                        vector<string>& ,
                        int&, 
                        bool&, 
                        int&, 
                        string&,
                        bool&, 
                        vector<int>*,
                        vector<int>*);
    
    
    vector<string> *list_taxa_to_insert;
    string _taxon_name;
    int path_size;
    int _taxon_to_insert;
    vector<int>* ids1;
    vector<int>* ids2;
    bool use_dynamic_order;
    int** path_up_to_now;

} Task;

void terrace_main(Params &params);

/*

 * Function for terrace analysis
 */
void runterraceanalysis(Params &params);
void runterraceanalysis_parallel(Params &params, int thread_num);

/*
 * Run terrace check:
 * - check if query trees lie on the same terrace with a representative tree. Naive pairwise comparison of induced partition trees.
 */
void run_terrace_check(Terrace *terrace, Params &params);


/**
        This is the function that calls generation of terrace trees, when the input data (representative tree and presence-absence matrix are processed already)
 */
void run_generate_trees(Terrace* terrace, Params &params,const int m);
void run_generate_trees_parallel(Terrace* terrace, Params &params, const int m, int thread_num = 0);

/**
        The function is used to read a set of subtrees to be considered as partition trees for terrace analysis
 */
void read_tree_set(const char *infile, bool &is_rooted, vector<TerraceTree*> &subtrees);


void split_threads(Terrace* init_terrace,
                vector<string> &list_taxa_to_insert, 
                int taxon_to_insert, 
                bool use_dynamic_order,
                string &taxon_name,
                vector<int> &ids1,
                vector<int> &ids2);


void generateTerrace_wrapper(Terrace* init_terrace,
                        vector<Terrace*> &part_tree_pairs, 
                        vector<string> &list_taxa_to_insert,
                        int &taxon_to_insert, 
                        bool &thread_call, 
                        int &new_thread_num, 
                        string &taxon_name,
                        bool &use_dynamic_taxon_order, 
                        vector<int>* ids1 = nullptr,
                        vector<int>* ids2 = nullptr);


// parallelization
void submitTask(Task* task);
void executeTask(Terrace* init_terrace, vector<Terrace*> *part_tree_pairs, Task* task);
void* startThread(Terrace* init_terrace, vector<Terrace*> *part_tree_pairs, int &thread_num);

#endif /* terraceanalysis_hpp */
