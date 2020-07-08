/*

BOOSTER: BOOtstrap Support by TransfER: 
BOOSTER is an alternative method to compute bootstrap branch supports 
in large trees. It uses transfer distance between bipartitions, instead
of perfect match.

Copyright (C) 2017 Frederic Lemoine, Jean-Baka Domelevo Entfellner, Olivier Gascuel

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef _TREE_H_
#define _TREE_H_

#include "hashtables_bfields.h"	/* for the hashtables to store taxa names on the branches */
#include "hashmap.h"
#include "io.h"
#include <ctype.h>

#define	TRUE	1
#define	FALSE	0

#define	MIN_BRLEN	1e-8
#define MAX_TREELENGTH	10000000 /* more or less 10MB for a tree file in NH format */
#define MAX_NODE_DEPTH	100000 /* max depth for nodes in the tree */

#define MAX_NAMELENGTH		255	/* max length of a taxon name */
#define MAX_COMMENTLENGTH	255	/* max length of a comment string in NHX format */

/* TYPES */
/* Every node in our binary trees has several neighbours with indices  0, 1, 2.... We allow polytomies of any degree.
   An internal node with no multifurcation has 3 outgoing directions/neighbours.

   In rooted trees, the interpretation is the following:
   - for internal nodes, direction 0 is to the father, other directions are to the sons
   - for tips (leaves), direction 0 is to the father, no other neighbours
   - the root has only two neighbours.
   So it's easy to check whether a node is a leaf (one neighbour) or the root of a rooted tree (two neighbours).

   For unrooted trees, the interpretation is the same, except that no node has two neighbours.
   The pseudo-root (from the NH import) is three-furcated, so behaves exactly like any other internal node.

   It is not advisable to have several nodes of degree two in the same tree. 

   The root or pseudo-root is ALWAYS assigned id 0 at beginning. May change later, upon rerooting.
   In any case, it is always pointed to by tree->node0. */

typedef struct __Node {
	char* name;
	char* comment;		/* for further use: store any comment (e.g. from NHX format) */
	int id;			/* unique id attributed to the node */
	short int nneigh;	/* number of neighbours */
	struct __Node** neigh;	/* neighbour nodes */
	struct __Edge** br;	/* corresponding branches going from this node */
	double depth;		/* the depth of a node is its min distance to a leaf */
} Node;


/* Every edge connects two nodes.
   By convention, all terminal branches will have the tip on their RIGHT end */


typedef struct __Edge {
	int id;
	Node *left, *right; /* in rooted trees the right end will always be the descendant.
			       In any case, a leaf is always on the right side of its branch. */
	double brlen;
	double branch_support;
	int* subtype_counts[2];			/* first index is 0 for the left of the branch, 1 for its right side */
	id_hash_table_t *hashtbl[2];		/* hashtables containing the ids of the taxa in each subtree */
						/* index 0 corresponds to the left of the branch, index 1 to its right.
						   following our implementation, we only keep hashtbl[1] populated. */
	short int had_zero_length; 		/* set at the moment when we read the tree, even though
				      		   we then immediately set the branch length to MIN_BRLEN */
	short int has_branch_support; 	
	int topo_depth;				/* the topological depth is the number of taxa on the lightest side of the bipar */
} Edge;


typedef struct __Tree {
	Node** a_nodes;	/* array of node pointers */
	Edge** a_edges;	/* array of edge pointers */
	Node* node0;	/* the root or pseudo-root node */
	int nb_nodes;
	int nb_edges;
	int nb_taxa;
	char** taxa_names; /* store only once the taxa names */
	int length_hashtables; /* the number of chained lists in the hashtables on the edges */
	int next_avail_node_id;
	int next_avail_edge_id;
	int next_avail_taxon_id;
	char** taxname_lookup_table;
} Tree;
	

/* FUNCTIONS */

/* UTILS/DEBUG: counting specific branches or nodes in the tree */
int count_zero_length_branches(Tree* tree);
int count_leaves(Tree* tree);
int count_roots(Tree* tree);
int count_multifurcations(Tree* tree);
int dir_a_to_b(Node* a, Node* b);

/* various statistics on branch support */
double mean_bootstrap_support(Tree* tree);
double median_bootstrap_support(Tree* tree);
int summary_bootstrap_support(Tree* tree, double* result);

/* parsing utils: discovering and dealing with tokens */
int index_next_toplevel_comma(char* in_str, int begin, int end);
int count_outer_commas(char* in_str, int begin, int end);
void strip_toplevel_parentheses(char* in_str, int begin, int end, int* pair);
int index_toplevel_colon(char* in_str, int begin, int end);
void parse_double(char* in_str, int begin, int end, double* location);

/* creating a node, a branch, a tree: to create a tree from scratch, not from parsing */
Node* new_node(const char* name, Tree* t, int degree);
Edge* new_edge(Tree* t);
Tree* new_tree(int nb_taxa, const char* name);
Node* graft_new_node_on_branch(Edge* target_edge, Tree* tree, double ratio_from_left, double new_edge_length, char* node_name);


/* collapsing a branch */
void collapse_branch(Edge* branch, Tree* tree);

/**
   This function removes a taxon from the tree (identified by its taxon_id)
   And recomputed the branch length of the branch it was branched on.

   Also, the bootstrap support (if any) is recomputed, taking the maximum of the
   supports of the joined branches

   Be careful: The taxnames_lookup_table is modified after this function!
   Do not use this function if you share the same taxnames_lookup_table in
   several trees.
*/
void remove_taxon(int taxon_id, Tree* tree);
/**
   If there remains 2 neighbors to connect_node
   We connect them directly and delete connect_node
   We keep l_edge and delete r_edge
   -> If nneigh de connect node != 2 : Do nothing
   This function deletes a node like that:
              connect_node
             l_edge   r_edge
     l_node *-------*--------* r_node
   => Careful: After this function, you may want to call 
   => recompute_identifiers()

*/
void remove_single_node(Tree *tree, Node *connect_node);
/**
   This method recomputes all the identifiers 
   of the nodes and of the edges
   for which the tree->a_nodes is not null
   or tree->a_edges is not null
   It also recomputes the total number of edges 
   and nodes in the tree
 */
void recompute_identifiers(Tree *tree);

/**
   This function shuffles the taxa of an input tree 
*/
void shuffle_taxa(Tree *tree);

/* (re)rooting a tree */
void reroot_acceptable(Tree* t);
void unrooted_to_rooted(Tree* t);

/* To be called after a reroot*/
void reorient_edges(Tree *t);
void reorient_edges_recur(Node *n, Node *prev, Edge *e);

/* utility functions to deal with NH files */
unsigned int tell_size_of_one_tree(const char* filename);
int copy_nh_stream_into_str(FILE* nh_stream, char* big_string);

/* actually parsing a tree */
void process_name_and_brlen(Node* son_node, Edge* edge, Tree* current_tree, char* in_str, int begin, int end);
Node* create_son_and_connect_to_father(Node* current_node, Tree* current_tree, int direction, char* in_str, int begin, int end);
void parse_substring_into_node(char* in_str, int begin, int end, Node* current_node, int has_father, Tree* current_tree);
Tree* parse_nh_string(char* in_str);

/* complete parse tree: parse NH string, update hashtables and subtype counts */
Tree *complete_parse_nh(char* big_string, char*** taxname_lookup_table);


/* taxname lookup table functions */
char** build_taxname_lookup_table(Tree* tree);
map_t build_taxid_hashmap(char** taxname_lookup_table, int nb_taxa);
void free_taxid_hashmap(map_t taxmap);
int free_hashmap_data(any_t arg,any_t key, any_t elemt);

char** get_taxname_lookup_table(Tree* tree);
Taxon_id get_tax_id_from_tax_name(char* str, char** lookup_table, int length);

/* (unnecessary/deprecated) multifurcation treatment */
void regraft_branch_on_node(Edge* branch, Node* target_node, int dir);

/***************************************************************
  ******************* neatly implementing tree traversals ******
***************************************************************/

void post_order_traversal_recur(Node* current, Node* origin, Tree* tree, void (*func)(Node*, Node*, Tree*));
void post_order_traversal(Tree* t, void (*func)(Node*, Node*, Tree*));

/* post order traversal with any data passed to the function call */
void post_order_traversal_data_recur(Node* current, Node* origin, Tree* tree, void*, void (*func)(Node*, Node*, Tree*, void*));
void post_order_traversal_data(Tree* t, void*, void (*func)(Node*, Node*, Tree*, void*));


void pre_order_traversal_recur(Node* current, Node* origin, Tree* tree, void (*func)(Node*, Node*, Tree*));
void pre_order_traversal(Tree* t, void (*func)(Node*, Node*, Tree*));

/* pre order traversal with any data passed to the function call */
void pre_order_traversal_data_recur(Node* current, Node* origin, Tree* tree, void* data, void (*func)(Node*, Node*, Tree*, void*));
void pre_order_traversal_data(Tree* t, void* data, void (*func)(Node*, Node*, Tree*, void*));

/* bootstrap values */
void update_bootstrap_supports_from_node_names(Tree* tree);
void update_bootstrap_supports_doer(Node* current, Node* origin, Tree* tree);


/* node depths */

void update_node_depths_post_doer(Node* target, Node* orig, Tree* t);
void update_node_depths_pre_doer(Node* target, Node* orig, Tree* t);
void update_node_depths_post_alltree(Tree* tree);
void update_node_depths_pre_alltree(Tree* tree);

/* topological depths */
void update_all_topo_depths_from_hashtables(Tree* tree);
int greatest_topo_depth(Tree* tree);

/* WORKING WITH HASHTABLES */

void update_hashtables_post_doer(Node* current, Node* orig, Tree* t);
void update_hashtables_pre_doer(Node* current, Node* orig, Tree* t);

void update_hashtables_post_alltree(Tree* tree);
void update_hashtables_pre_alltree(Tree* tree);


/* UNION AND INTERSECT CALCULATIONS FOR THE TRANSFER METHOD (from Bréhélin/Gascuel/Martin 2008) */
void update_i_c_post_order_ref_tree(Tree* ref_tree, Node* orig, Node* target, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix);
void update_all_i_c_post_order_ref_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix);

void update_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, Node* orig, Node* target, short unsigned** i_matrix, short unsigned** c_matrix, short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge);
void update_all_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix, short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge);


/*Generate Random Tree*/
/**
   - nbr_taxa: Number of leaves in the output tree
   - taxa_names: array of leaf names: must be NULL if you do not want to use it 
   (names will be numbered in this case)

   - The output tree has branch lengths attributed using a normal distribution  N(0.1,0.05), and any br len < 0 is set to 0
*/
Tree * gen_rand_tree(int nbr_taxa, char **taxa_names);

/* writing a tree */

void write_nh_tree(Tree* tree, FILE* stream);
void write_subtree_to_stream(Node* node, Node* node_from, FILE* stream);

/* freeing stuff */

void free_edge(Edge* edge);
void free_node(Node* node);
void free_tree(Tree* tree);
#endif /* _TREE_H_ */
