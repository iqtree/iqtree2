/*
 *  presenceabsencematrix.hpp
 *  Class of presence-absence matrices (gene pr-ab)
 *  Created on: Sep 10, 2020
 *      Author: Olga
 */

#ifndef presenceabsencematrix_hpp
#define presenceabsencematrix_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "tree/mtree.h"
#include "terracetree.hpp"

class PresenceAbsenceMatrix {
    
public:
    
    PresenceAbsenceMatrix();
    ~PresenceAbsenceMatrix();
    
    vector<IntVector> pr_ab_matrix;
    vector<string> taxa_names;
    
    vector<IntVector> overlap_matrix;
    
    int taxa_num{0};
    int part_num{0};
    int uniq_taxa_num{0};
    int uniq_taxa_to_insert_num{0};
    double missing_percent{0.0};
    
    /*
     *  initialization
     */
    void init();
    
    /*
     *  Fill in presence-absence matrix from an input alignment
     */
    void get_from_alignment(Params &params);
    
    /*
     *  Build presence-absence matrix from input set of partition trees
     */
    void get_from_subtrees(vector<TerraceTree*> subtrees);

    /*
     *  Reading a presence-absence matrix for supermatrix
     */
    void read_pr_ab_matrix(const char *infile);
    void read_pr_ab_matrix(istream &in);
    
    void print_pr_ab_matrix(ostream &out = cout);
    
    /*
     * Find taxon id by taxon name
     */
    
    int findTaxonID(string taxon_name);
    
    /*
     *  a NodeVector part_taxa of (pointers to) leaf nodes from the tree is filled according to the presence-absence info for part
     */
    void getPartTaxa(int part, MTree *tree, MTree *part_tree, NodeVector &part_taxa);
    
    /*
     *  reorder vector of taxa_names and pr_ab_matrix according to the the order of taxa on a tree (based on node->ids)
     */
    void reorderAccordingToTree(NodeVector taxa_nodes);
    
    /*
     *   A variable to keep track if the reordering of taxa according to the tree was already performed
     */
    bool flag_reorderAccordingToTree;
    
    /*
     *  Get a submatrix of presence-absence matrix corresponding to taxa passed through taxa_names.
     */
    
    void getSubPrAbMatrix(vector<string> taxa_names_subset, PresenceAbsenceMatrix *submatrix, IntVector *parts = nullptr);
    void getSubPrAbMatrix(NodeVector taxon_nodes, PresenceAbsenceMatrix *submatrix, IntVector *parts = nullptr);
    
    /*
     *  Extend presence-absence matrix by the new taxon
     */
    
    void extend_by_new_taxa(string taxon_name, IntVector pr_ab_pattern);
    
    /*
     *  Remove a taxon from presence_absence matrix
     */
    
    void remove_taxon(string taxon_name);
    
    /*
     *  Function to get informaton about initial tree and a taxon order for stepwise insertion.
     *  taxa_names_sub - list of taxa to be on the intial tree (TODO: maybe get the tree right away?).
     *  list_taxa_to_insert - taxa in the order to be inserted
     */
    
    int getINFO_init_tree_taxon_order(vector<string> &taxa_names_sub, vector<string> &list_taxa_to_insert,const int m);
    
    /*
     *  Order partitions by their overlap
     */
    void orderPartByOverlap(IntVector &ordered_partitions,IntVector &part_cov);
    
    /*
     *  Order partitions by their overlap within the group
     */
    void orderPartByOverlap_within(int upper_lim, IntVector &ordered_partitions, IntVector &group, IntVector &part_cov);
    
    /*
     *  Order partitions by their overlap with respect to preceding partitions
     */
    void orderPartByOverlap_preceding(int upper_lim, IntVector &ordered_partitions, IntVector &group, IntVector &part_cov);
    
    /*
     *  Get pairwise taxon overlap between partition part and all other partitions in a group
     */
    void getPartOverlap(int part, IntVector &group, IntVector &part_cov, IntVector &overlap);
    int get2PartOverlap(int part_1, int part_2, int max_overlap);
    
    void reordering(vector<IntVector> &new_order, vector<IntVector>::iterator it_b, vector<IntVector>::iterator it_e, IntVector &part_cov, int i=0);
    /*
     *  Order partitions based on the overlap with one preceding partition. (If a subgroup has the same overlap with the considered partition, analyse overlap with the next partition, till all preceding partitions)
     */
    
    /*
     *  Order taxa by their coverage
     */
    
    void orderTaxaByCoverage(vector<int> &taxon_ids, vector<IntVector> &coverage_info, IntVector &ordered_taxa);
    
    /*
     *  Compute amount of missing sequences (in %)
     */
    
    void percent_missing();
    
    /*
     * Compute overlap between partitions
     */
    void getPartOverlapComplete();
    void print_overlap_matrix(ostream &out);
};

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree);

#endif /* presenceabsencematrix_hpp */
