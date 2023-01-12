/*
 *  terrace.hpp
 *  Terrace class is used to generate trees that have idetical set of induced partition trees.
 *  Created on: Sep 14, 2020
 *      Author: Olga
 */

#ifndef terrace_hpp
#define terrace_hpp

#include <stdio.h>
#include "terracetree.hpp"
#include "terracenode.hpp"
#include "presenceabsencematrix.hpp"
#include "ParallelContext.hpp"


class Terrace: public TerraceTree
{
public:
    /*
     * constructor
     */
    Terrace();

    /*
     * constructor
     */
    Terrace(TerraceTree &tree, std::shared_ptr<PresenceAbsenceMatrix> &m);

    /*
     * constructor
     */
    Terrace(TerraceTree &tree, vector<TerraceTree*> input_induced_trees);
    
    /*
     *  constructor
     */
    Terrace(const char *infile_tree, bool is_rooted,const char *infile_matrix);
    
    /*
     * constructor
     */
    Terrace(vector<TerraceTree*> input_induced_trees);
    
    /*
     * destructor
     */
    ~Terrace();
    
    /*
     *  Initialisator
     */
    void init();
    
    vector<TerraceTree*> induced_trees;
    std::shared_ptr<PresenceAbsenceMatrix> matrix;
    StrVector terrace_trees;
    
    int taxa_num;
    int part_num;

    /**
     to track, if the initial tree was generate using backward approach: remove m leaves instead of take induced partition tree
     */
    int rm_leaves{0};
    
    /**
        Master terrace. The pointer to a terrace, which is beeing explored using current agile terrace. (this) is derrived from master.
     */
    Terrace* master_terrace{nullptr};
    
    /*
     * Original taxon names
     */
    StrVector taxa_names_orgn;
    
    /*
     *  Number of trees on terrace
     */
    unsigned int terrace_trees_num;
    unsigned int cur_terrace_trees;
    
    /*
     *  Number of intermediate trees visited
     */
    int intermediated_trees_num;
    int cur_intermediate_trees;
    
    /*
     *  Number of dead ends encountered
     */
    unsigned int dead_ends_num;

    double start_real_time;
    double terrace_trees_time;
    double intermediate_trees_time;
    double dead_ends_time;
    
    // file to output all generated terrace trees
    string out_file;
    ofstream out;
    bool terrace_out;
    int trees_out_lim;
    
    // Stopping rules
    unsigned int terrace_max_trees;
    unsigned int intermediate_max_trees;
    int seconds_max;

    // for parallelization
    shared_ptr<ParallelContext> pContext;  // have to modify that
    
    int remaining_threads_to_assign;
    int parallel_case;
    
    int artificial_thread_num;
    int real_thread_num;

    bool initial_split_done = true;
    bool parallel_exec = false;
    bool stop = false;
    double threshold;
    int taskThreshold;
    int initial_split_done_index;

    int **path_up_to_now;
    int path_size;
    
    char** inserted_taxa;
    int max_string_size;

    /*
    * Adds/substracts one to global variables - locks and unlocks the mutex
    */
    
    /*
     *  Print terrace info: a representative tree, induced trees and presence-absence matrix
     */
    void printInfo(ostream &out=cout);
    
    /*
     *  get induced partition trees
     */
    void get_part_trees();
    
    /*
     *  set induced partition trees, if they already exist
     */
    void set_part_trees(vector<TerraceTree*> input_induced_trees);
    void unset_part_trees();
    
    /*
     * link parent tree and induced partition trees
     */
    
    void linkTrees(bool back_branch_map, bool back_taxon_map);
    
    /*
     * link parent tree and a single induced partition tree
     */
    void linkTree(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node = nullptr, TerraceNode *dad = nullptr);

    /*
     *  link one branch of parent tree on partition tree part (for internal branches only, oder?)
     */
    void linkBranch(int part, TerraceNeighbor *nei, TerraceNeighbor *dad_nei, bool back_branch_map, bool back_taxon_map);
    
    /*
     *  Local map update after insertion/deletion of the taxon
     */
    
    void update_map(int part, NodeVector &part_taxa, bool back_branch_map, bool back_taxon_map, TerraceNode *node, TerraceNode *dad = nullptr);
    
    /*
     *  Print Info about the branch between parent tree and induced partition trees
     */
    void printMapInfo(int partition = -1);
    
    /*
     *  Print Info about inverse branch and taxon images from induced partition trees to master and upper level induced trees, respectivelly
     */
    void printBackMapInfo();
    
    /*
     *  For all nodes of the tree clear information about empty branches and empty nodes, which might have remained after the previous partition
     */
    
    void clearEmptyBranchAndTaxaINFO(TerraceNode *node = nullptr, TerraceNode *dad = nullptr);
    
    /*
     *  For a given taxon name get allowed branches. aux_terrace contains pointers to top level induced partition trees (which are stored as terraces)
     */
    //void getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec, NodeVector *branch_end_1 = nullptr, NodeVector *branch_end_2 = nullptr);
    void getAllowedBranches(string taxon_name, vector<Terrace*> aux_terrace, NodeVector *node1_vec, NodeVector *node2_vec);
    
    /*
     *  Insert a new taxon to the parent tree, update induced partition trees, update mapping (locally)
     */
    void extendNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> aux_terrace, bool constructing_thread = false);
    
    /*
     *  Clean all link neighbours and taxa on parent tree and on induced partition trees
     */
    
    void cleanAllLinkNeighboursAndTaxa(bool clean_induced_part_maps = false);
    
    /*
     *  Prepare top-low induced partition tree pairs: induced tree from the terrace and a common subtree with the initial tree (to be modified by inserting new taxa). Top level provided by the passed terrace, low level by the current terrace (which is initial terrace).
     */
    
    void create_Top_Low_Part_Tree_Pairs(vector<Terrace*> &part_tree_pairs, Terrace *terrace);

    /*
     *  The main function to generate trees by recursive taxon insertion
     */


    void generateTerraceTrees(Terrace *terrace, 
                            vector<Terrace*> &part_tree_pairs, 
                            vector<string> &list_taxa_to_insert, 
                            int taxon_to_insert = -1, 
                            bool use_dynamic_order = true, 
                            bool thread_call = false,
                            string _taxon_name = "",
                            vector<int>* ids1 = nullptr,
                            vector<int>* ids2 = nullptr);
    
    /*
     *  Parallel version
     */
    
    void print_debug(NodeVector &node1_vec_branch, NodeVector &node2_vec_branch, int taxon_to_insert, bool give_birth = false, int thread_num = 0);

    void sorting_vector(NodeVector &node1_vec_branch, NodeVector &node2_vec_branch);
    void keepIndicesAssignedToThread(NodeVector &node1_vec_branch, NodeVector &node2_vec_branch); 
    void keepIndicesAssignedToThread_initial(NodeVector &node1_vec_branch, NodeVector &node2_vec_branch, int parallel_case, int thread_num, int nun_threads_remaining);
    /*
     *  Get next taxon to be inserted - a taxon with the least number of allowed branches
     */
    string getNextTaxon(vector<Terrace*> &part_tree_pairs, vector<string> *list_taxa_to_insert, int taxon_to_insert, NodeVector &node1_vec_main, NodeVector &node2_vec_main);
    

    /*
     *  Remove one taxon from the terrace tree, from induced partition trees, update the mapping
     *  Here the update is local update: only involved branches and maps of trees are beeing updated.
     */
    
    void remove_one_taxon(string taxon_name, vector<Terrace*> part_tree_pairs);

    /*
     *  Re-link
     */
    void relinkALL(vector<Terrace*> part_tree_pairs);
    
    /*
     * Print initial tree, print TOP and LOW induced trees
     */
    void print_ALL_DATA(vector<Terrace*> part_tree_pairs);
    
    /*
     *  Rename taxa on a tree and in presence-absence matrix
     */
    void renameTaxa();
    
    /*
     *  Check if the query tree has the same set of induced partition trees
     */
    bool check_two_trees(MTree* query_tree);
    
    /*
     *  Write all generated trees to file
     */
    void write_terrace_trees_to_file();
    
    /*
     *  Write summary of generating trees from a terrace
     */
    void write_summary_generation();
    
    /*
     *  Write warning depending on the activated stopping rule
     */
    void write_warning_stop(int type);
};

#endif /* terrace_hpp */
