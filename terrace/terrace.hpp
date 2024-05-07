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
    Terrace(TerraceTree tree, PresenceAbsenceMatrix *m);
    
    /*
     * constructor
     */
    Terrace(TerraceTree tree, PresenceAbsenceMatrix *m, vector<TerraceTree*> input_induced_trees);
    
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
    PresenceAbsenceMatrix *matrix;
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
    Terrace *master_terrace{nullptr};
    
    /*
     * Original taxon names
     */
    StrVector taxa_names_orgn;
    
    /*
     *  Number of trees on terrace
     */
    unsigned int terrace_trees_num;
    
    /*
     *  Number of intermediate trees visited
     */
    unsigned int intermediated_trees_num;
    
    /*
     *  Number of dead ends encountered
     */
    unsigned int dead_ends_num;
    
    // file to output all generated terrace trees
    string out_file;
    ofstream out;
    bool terrace_out;
    int trees_out_lim;
    
    // Stopping rules
    unsigned int terrace_max_trees;
    unsigned int intermediate_max_trees;
    int seconds_max;
    
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
    void extendNewTaxon(string node_name, TerraceNode *node_1_branch, TerraceNode *node_2_branch, vector<Terrace*> aux_terrace);
    
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
    
    void generateTerraceTrees(Terrace *terrace, vector<Terrace*> &part_tree_pairs, vector<string> &list_taxa_to_insert, int taxon_to_insert = -1,vector<string> *ordered_taxa_to_insert = nullptr);
    
    /*
     *  Get next taxon to be inserted - a taxon with the least number of allowed branches
     */
    string getNextTaxon(vector<Terrace*> &part_tree_pairs, vector<string> *ordered_taxa_to_insert,NodeVector &node1_vec_main, NodeVector &node2_vec_main);
    
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
