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

#include "tree.h"
#include "externs.h"

int ntax;		/* this global var is set here, in parse_nh */

/* UTILS/DEBUG: counting specific branches or nodes in the tree */

int count_zero_length_branches(Tree* tree) {
	int count = 0;
	int i, n = tree->nb_edges;
	for (i = 0; i < n; i++) if(tree->a_edges[i]->had_zero_length) count++;
	return count;
}

int count_leaves(Tree* tree) {
	int count = 0;
	int i, n = tree->nb_nodes;
	for (i = 0; i < n; i++) if(tree->a_nodes[i]->nneigh == 1) count++;
	return count;
}

int count_roots(Tree* tree) { /* to ensure there is exactly zero or one root */
	int count = 0;
	int i, n = tree->nb_nodes;
	for (i = 0; i < n; i++) if(tree->a_nodes[i]->nneigh == 2) count++;
	return count;
}

int count_multifurcations(Tree* tree) { /* to ensure there is exactly zero or one root */
	int count = 0;
	int i, n = tree->nb_nodes;
	for (i = 0; i < n; i++) if(tree->a_nodes[i]->nneigh > 3) count++;
	return count;
}

int dir_a_to_b(Node* a, Node* b) {
	/* this returns the direction from a to b when a and b are two neighbours, otherwise yields an error */
	int i, n = a->nneigh;
	for(i=0; i<n; i++) if (a->neigh[i] == b) break;
	if (i < n) return i; else {
	  fprintf(stderr,"Fatal error : nodes are not neighbours.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
	return -1;
} /* end dir_a_to_b */


/* various statistics on tree branch support */

double mean_bootstrap_support(Tree* tree) {
	/* this function returns the mean bootstrap support calculated on those branches that have a bootstrap support value */
	int i, total_num = 0;
	double accu = 0.0;
	int n_br = tree->nb_edges;
	for(i = 0; i < n_br; i++) {
		if (tree->a_edges[i]->has_branch_support) {
			accu += tree->a_edges[i]->branch_support;
			total_num++;
		}
	} /* end for */

	return accu / total_num;
} /* end mean_bootstrap_support */



double median_bootstrap_support(Tree* tree) {
	/* this function returns the median bootstrap support calculated on those branches that have a bootstrap support value */
	/* we first create an array with all bootstrap supports */
	int i, j,  total_num = 0, n_br = tree->nb_edges;
	for(i = 0; i < n_br; i++) if(tree->a_edges[i]->has_branch_support) total_num++;
        double* branch_supports = (double*) malloc (total_num * sizeof(double));

	j=0;
	for(i = 0; i < n_br; i++) if(tree->a_edges[i]->has_branch_support) branch_supports[j++] = tree->a_edges[i]->branch_support;
	double result = median_double_vec(branch_supports, total_num);
	free(branch_supports);
	return result;
} /* end median_bootstrap_support */
		

int summary_bootstrap_support(Tree* tree, double* result) {
	/* this function stores all the bootstrap values in a vector and outputs the statistical summary
	   of that vector into the result array. Same order as in R. */
	/* RESULT MUST HAVE ALLOCATED SIZE >= 6 */
	/* retcode is -1 in case no support values found */
	int i, j, num_bootstrap_values = 0, n_br = tree->nb_edges;
	for(i = 0; i < n_br; i++) if(tree->a_edges[i]->has_branch_support) num_bootstrap_values++;

	if (num_bootstrap_values == 0) return -1;

	/* allocating vector */
	double* bootstrap_vals = (double*) malloc(num_bootstrap_values * sizeof(double));
	/* filling in vector */
	for(i = j = 0; i < n_br; i++) if(tree->a_edges[i]->has_branch_support) bootstrap_vals[j++] = tree->a_edges[i]->branch_support;
	/* summary */
	summary_double_vec_nocopy(bootstrap_vals, num_bootstrap_values, result);
	free(bootstrap_vals);
	return 0;
} /* end summary_bootstrap_support */




/* parsing utils: discovering tokens */

int index_next_toplevel_comma(char* in_str, int begin, int end) {
	/* returns the index of the next toplevel comma, from position begin included, up to position end.
	   the result is -1 if none is found. */
	int level = 0, i;
	for (i = begin; i <= end; i++) {
		switch(in_str[i]) {
			case '(':
				level++;
				break;
			case ')':
				level--;
				break;
			case ',':
				if (level == 0) return i;
		} /* endswitch */
	} /* endfor */
	return -1; /* reached if no outer comma found */
} /* end index_next_toplevel_comma */



int count_outer_commas(char* in_str, int begin, int end) {
	/* returns the number of toplevel commas found, from position begin included, up to position end. */
	int count = 0, level = 0, i;
	for (i = begin; i <= end; i++) {
		switch(in_str[i]) {
			case '(':
				level++;
				break;
			case ')':
				level--;
				break;
			case ',':
				if (level == 0) count++;
		} /* endswitch */
	} /* endfor */
	return count;
} /* end count_outer_commas */



void strip_toplevel_parentheses(char* in_str, int begin, int end, int* pair) {
	/* returns the new (begin,end) pair comprising all chars found strictly inside the toplevel parentheses.
	   The input "pair" is an array of two integers, we are passing the output values through it.
	   It is intended that here, in_str[pair[0]-1] == '(' and in_str[pair[1]+1] == ')'.
	   In case no matching parentheses are simply return begin and end in pair[0] and pair[1]. It is NOT an error. */
	/* This function also tests the correctness of the NH syntax: if no balanced pars, then return an error and abort. */
	int i, found_par = 0;
	
	pair[0] = end+1; pair[1] = -1; /* to ensure termination if no parentheses are found */

	/* first seach opening par from the beginning of the string */
	for (i = begin; i <= end; i++) if (in_str[i] == '(') { pair[0] = i+1; found_par += 1; break; } 

	/* and then search the closing par from the end of the string */
	for (i = end; i >= begin; i--) if (in_str[i] == ')') { pair[1] = i-1; found_par += 1; break; } 

	switch (found_par) {
		case 0:
			pair[0] = begin;
			pair[1] = end;
			break;
		case 1:
		  fprintf(stderr,"Syntax error in NH tree: unbalanced parentheses between string indices %d and %d. Aborting.\n", begin, end);
		  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	} /* end of switch: nothing to do in case 2 (as pair[0] and pair[1] correctly set), and found_par can never be > 2 */
}


int index_toplevel_colon(char* in_str, int begin, int end) {
	/* returns the index of the (first) toplevel colon only, -1 if not found */
	int level = 0, i;
	for (i = end; i >= begin; i--) {/* more efficient to proceed from the end in this case */
		switch(in_str[i]) {
			case ')':
				level++;
				break;
			case '(':
				level--;
				break;
			case ':':
				if (level == 0) return i;
		} /* endswitch */
	} /* endfor */
	return -1;
} /* end index_toplevel_colon */


void parse_double(char* in_str, int begin, int end, double* location) {
	/* this function parses a numerical value and puts it into location. Meant to be used for branch lengths. */
	if (end < begin) {
	  fprintf(stderr,"Missing branch length at offset %d in the New Hampshire string. Branch length set to 0.\n", begin);
	  sscanf("0.0", "%lg", location);
	  return;
	}
	char numerical_string[52] = { '\0' };
	strncpy(numerical_string, in_str+begin, end-begin+1);
	int n_matches = sscanf(numerical_string, "%lg", location);
	if (n_matches != 1) {
	  fprintf(stderr,"Fatal error in parse_double: unable to parse a number out of \"%s\". Aborting.\n", numerical_string);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
} /* end parse_double */




/* CREATION OF A NEW TREE FROM SCRATCH, ADDING TAXA ONE AT A TIME */

Node* new_node(const char* name, Tree* t, int degree) {
	int i;
	Node* nn = (Node*) malloc(sizeof(Node));
	nn->nneigh = degree;
	nn->neigh = malloc(degree * sizeof(Node*));
	nn->br = malloc(degree * sizeof(Edge*));
	nn->id = t->next_avail_node_id++;
	if(degree==1 && !name) { fprintf(stderr,"Fatal error : won't create a leaf with no name. Aborting.\n"); Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);}
	if(name) { nn->name = strdup(name); } else nn->name = NULL;
	if(degree==1) { t->taxa_names[t->next_avail_taxon_id++] = strdup(name); }
	nn->comment = NULL;
	for(i=0; i < nn->nneigh; i++) { nn->neigh[i] = NULL; nn->br[i] = NULL; }
	nn->depth = MAX_NODE_DEPTH;
	t->a_nodes[nn->id] = nn; /* warning: not checking anything here! This array haas to be big enough from start */
	t->nb_nodes++;
	return nn;
}

Edge* new_edge(Tree* t) {
	Edge* ne = (Edge*) malloc(sizeof(Edge));
	ne->id = t->next_avail_edge_id++;
	ne->has_branch_support = 0;
	ne->hashtbl[0] = ne->hashtbl[1] = NULL;
	ne->subtype_counts[0] = ne->subtype_counts[1] = NULL;
	t->a_edges[ne->id] = ne;
	t->nb_edges++;
	return ne;
}


Tree* new_tree(int nb_taxa, const char* name) {
	/* allocates the space for a new tree and gives it as an output (pointer to the new tree) */
	/* optional is the name of the first taxa. If we don't provide it, there exists a risk that we will build a
	   tree with finally one leaf with no name */
	if (nb_taxa <= 0) return NULL; /* at least one node, that is node0 */
	Tree* t = (Tree*) malloc(sizeof(Tree));
	t->taxa_names = (char**) calloc(nb_taxa, sizeof(char*)); /* store only once the taxa names */
	t->next_avail_node_id = t->next_avail_edge_id = t->next_avail_taxon_id = t->nb_nodes = t->nb_edges = 0;
	t->nb_taxa = nb_taxa; /* here we don't put the actual number of taxa, but the value to be reached by growing the tree */

	t->a_nodes = (Node**) calloc(2*nb_taxa-1, sizeof(Node*));	/* array of node pointers, enough for a rooted tree */
	t->a_edges = (Edge**) calloc(2*nb_taxa-2, sizeof(Edge*));	/* array of edge pointers, enough for a rooted tree */
	
	t->node0 = new_node(name, t, 1);	/* this first node _is_ a leaf */

	t->taxname_lookup_table = NULL;
	return t;
}


/* for the moment this function is used to create binary trees (where all internal nodes have three neighbours) */
Node* graft_new_node_on_branch(Edge* target_edge, Tree* tree, double ratio_from_left, double new_edge_length, char* node_name) {
	/* this grafts a new node on an existing branch. the ratio has to be between 0 and 1, and is relative to the "left" tip of the branch */
	int orig_dir_from_node_l, orig_dir_from_node_r;

	if(tree == NULL) {
	  fprintf(stderr,"Error : got a NULL tree pointer. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	if(ratio_from_left <= 0 && ratio_from_left >= 1) {
	  fprintf(stderr,"Error : invalid ratio %.2f for branch grafting. Aborting.\n", ratio_from_left);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	if(new_edge_length <= 0) {
	  fprintf(stderr,"Error : nonpositive new branch length %.2f. Aborting.\n", new_edge_length);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}


	if(node_name == NULL) {
	  fprintf(stderr,"Error : won't create a leaf with no name. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	if(target_edge == NULL) {
		/* here we treat the special case of the insertion of the second node (creation of the very first branch) */
		if (tree->nb_edges!= 0 || tree->next_avail_node_id != 1 || tree->next_avail_edge_id != 0) {
		  fprintf(stderr,"Error : I get a NULL branch pointer while there is at least one existing branch in the tree. Aborting.\n");
		  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
		}
		Node* second_node = new_node(node_name, tree, 1); /* will be the right node, also a leaf */
		Edge* only_edge = new_edge(tree);
		only_edge->left = tree->node0;
		only_edge->right = second_node;
		only_edge->brlen = new_edge_length;
		only_edge->had_zero_length = 0;

		second_node->neigh[0] = tree->node0; tree->node0->neigh[0] = second_node;
		second_node->br[0] = tree->node0->br[0] = only_edge;

		return second_node;

	} /* end of the treatment of the insertion in the case of the second node */

	if(tree->a_edges[target_edge->id] != target_edge) {
	  fprintf(stderr,"Error : wrong edge id rel. to the tree. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
	
	/* create two new nodes in the tree: the father and the son. The father breaks the existing edge into two. */
	/* Steps:
	   (1) create a new node, the breaking point
	   (2) create a new edge (aka. right edge, because the target edge remains left of the breakpoint)
	   (3) shorten the initial edge and give length value to the new edge
	   (4) rearrange the tips and update the node
	   (4bis) VERY IMPORTANT: FOR EACH END OF THE INITIAL EDGE THAT IS A LEAF, MAKE SURE THAT THE BREAKPOINT IS IN DIR 0 FROM IT
	   (5) create the son node
	   (6) create the edge leading to this son and update it and the node
	   */

	/* record the original situation */
	Node* node_l = target_edge->left;
	Node* node_r = target_edge->right;
	orig_dir_from_node_l = dir_a_to_b(node_l,node_r);
	orig_dir_from_node_r = dir_a_to_b(node_r,node_l);

	/* (1) */
	Node* breakpoint = new_node(NULL, tree, 3); /* not a leaf, so has three neighbours */

	/* (2) */
	Edge* split_edge = new_edge(tree); /* the breakpoint sits between the target_edge and the split_edge */

	/* (3) */
	split_edge->brlen = 2.0 * (1.0 - ratio_from_left) * target_edge->brlen; /* double the length so that we never get tiny edges after multiple insertions on the same branch */
	split_edge->had_zero_length = 0;
	target_edge->brlen *= 2.0 * ratio_from_left;

	/* (4) */
	/* edge tips */

	split_edge->left = breakpoint;
	split_edge->right = node_r;
	target_edge->right = breakpoint;

	if(node_l->nneigh ==1){
	  /* Case of the first edge that connects TWO leaves
	     We need to connect the left leaf to the right side of the branch
	     to be consistent with the tree definition
	  */
	  target_edge->right = target_edge->left;
	  target_edge->left = breakpoint;
	}

	/* update the 3 nodes */
	breakpoint->neigh[0] = node_l;
	breakpoint->br[0] = target_edge;

	breakpoint->neigh[1] = node_r;
	breakpoint->br[1] = split_edge;

	/* (4bis) */
	if (node_l->nneigh == 1 && orig_dir_from_node_l != 0) { /* change direction to 0 */
		node_l->neigh[0] = breakpoint;
		node_l->br[0] = target_edge;
		node_l->neigh[orig_dir_from_node_l] = NULL;
		node_l->br[orig_dir_from_node_l] = NULL;
	} else {
		node_l->neigh[orig_dir_from_node_l] = breakpoint;
		/* target_edge was already registered as the branch in this direction */
	}

	if (node_r->nneigh == 1 && orig_dir_from_node_r != 0) { /* change direction to 0 */
		node_r->neigh[0] = breakpoint;
		node_r->br[0] = split_edge;
		node_r->neigh[orig_dir_from_node_r] = NULL;
		node_r->br[orig_dir_from_node_r] = NULL;
	} else {
		node_r->neigh[orig_dir_from_node_r] = breakpoint;
		node_r->br[orig_dir_from_node_r] = split_edge;
	}

	/* (5) */
	Node* son = new_node(node_name, tree, 1); /* a leaf */

	/* (6) */
	Edge* outer_edge = new_edge(tree);
	outer_edge->left = breakpoint;
	outer_edge->right = son; /* the leaf is right of the branch */
	outer_edge->brlen = new_edge_length;
	outer_edge->had_zero_length = (new_edge_length == 0); /* but was already ruled out, see beginning of func. */

	son->neigh[0] = breakpoint; breakpoint->neigh[2] = son; /* necessarily the father is in direction 0 from the leaf */
	son->br[0] = breakpoint->br[2] = outer_edge;

	return son;
}

/* collapsing a branch */
void collapse_branch(Edge* branch, Tree* tree) {
	/* this function collapses the said branch and creates a higher-order multifurcation (n1 + n2 - 2 neighbours for the resulting node).
	   We also have to remove the extra node from tree->a_nodes and the extra edge from t->a_edges.
	   to be done:
	   (1) create a new node with n1+n2-2 neighbours. Ultimately we will destroy the original node.
	   (2) populate its list of neighbours from the lists of neighbours corresponding to the two original nodes
	   (3) populate its list of neighbouring edges form the br lists of the two original nodes
	   (4) for each of the neighbours, set the info regarding their new neighbour (that is, our new node)
	   (5) for each of the neighbouring branches, set the info regarding their new side (that is, our new node)
	   (6) destroy the two original nodes and commit this info to a_nodes. Modify tree->nb_nodes
	   (7) destroy the original edge and commit this info to a_edges. Modify tree->nb_edges */

	/* WARNING: this function won't accept to collapse terminal edges */
	Node *node1 = branch->left, *node2 = branch->right;
	int i, j, n1 = node1->nneigh, n2 = node2->nneigh;
	if (n1 == 1 || n2 == 1) { fprintf(stderr,"Warning: %s() won't collapse terminal edges.\n",__FUNCTION__); return; }
	int degree = n1+n2-2;
	/* (1) */
	/* Node* new = new_node("collapsed", tree, n1 + n2 - 2); */ /* we cannot use that because we want to reuse n1's spot in tree->a_nodes */
	Node* new = (Node*) malloc(sizeof(Node));
	new->nneigh = degree;
	new->neigh = malloc(degree * sizeof(Node*));
	new->br = malloc(degree * sizeof(Edge*));
	new->id = node1->id; /* because we are going to store the node at this index in tree->a_nodes */
	new->name = strdup("collapsed");
	new->comment = NULL;
	new->depth = min_int(node1->depth, node2->depth);

	/* very important: set tree->node0 to new in case it was either node1 or node2 */
	if (tree->node0 == node1 || tree->node0 == node2) tree->node0 = new;


	int ind = 0; /* index in the data structures in new */
	/* (2) and (3) and (4) and (5) */
	for (i=0; i < n1; i++) {
		if (node1->neigh[i] == node2) continue;
		new->neigh[ind] = node1->neigh[i];
		/*  then change one of the neighbours of that neighbour to be the new node... */
		for (j=0; j < new->neigh[ind]->nneigh; j++) {
			if(new->neigh[ind]->neigh[j] == node1) {
				new->neigh[ind]->neigh[j] = new;
				break;
			}
		} /* end for j */

		new->br[ind] = node1->br[i];
		/* then change one of the two ends of that branch to be the new node... */
		if (new->neigh[ind] == new->br[ind]->right) new->br[ind]->left = new; else new->br[ind]->right = new; 
		ind++;
	}

	for (i=0; i < n2; i++) {
		if (node2->neigh[i] == node1) continue;
		new->neigh[ind] = node2->neigh[i];
		/*  then change one of the neighbours of that neighbour to be the new node... */
		for (j=0; j < new->neigh[ind]->nneigh; j++) {
			if(new->neigh[ind]->neigh[j] == node2) {
				new->neigh[ind]->neigh[j] = new;
				break;
			}
		} /* end for j */

		new->br[ind] = node2->br[i];
		/* then change one of the two ends of that branch to be the new node... */
		if (new->neigh[ind] == new->br[ind]->right) new->br[ind]->left = new; else new->br[ind]->right = new; 
		ind++;
	}

	/* (6) tidy up tree->a_nodes and destroy old nodes */
	assert(tree->a_nodes[new->id] == node1);
	tree->a_nodes[new->id] = new;
	/* current last node in tree->a_edges changes id and is now placed at the position were node2 was */
	int id2 = node2->id;
	assert(tree->a_nodes[id2] == node2);
	tree->a_nodes[id2] = tree->a_nodes[-- tree->next_avail_node_id]; /* moving the last node into the spot occupied by node2... */
	tree->a_nodes[id2]->id = id2;					/* and changing its id accordingly */
	tree->a_nodes[tree->next_avail_node_id] = NULL; /* not strictly necessary, but... */
	tree->nb_nodes--;
	free_node(node1);
	free_node(node2);

	/* (7) tidy up tree->a_edges and destroy the old branch */
	assert(tree->a_edges[branch->id] == branch);
	tree->a_edges[branch->id] = tree->a_edges[-- tree->next_avail_edge_id]; /* moving the last branch into the spot occupied by 'branch' */
	tree->a_edges[branch->id]->id = branch->id; 				/* ... and changing its id accordingly */
	tree->a_edges[tree->next_avail_edge_id] = NULL; /* not strictly necessary, but... */
	tree->nb_edges--;
	free_edge(branch);

} /* end collapse_branch */


/**
   This function removes a taxon from the tree (identified by its taxon_id)
   And recomputed the branch length of the branch it was branched on.

   Be careful: The taxnames_lookup_table is modified after this function!
   Do not use this function if you share the same taxnames_lookup_table in
   several trees.

               connect_node
             l_edge   r_edge
     l_node *-------*--------* r_node
                    |e_to_remove_index
                    | e_to_remove
                    |
                    *
                n_to_remove
*/
void remove_taxon(int taxon_id, Tree* tree){
  Node *n_to_remove = NULL;
  Edge *e_to_remove,  *r_edge;
  Node *connect_node, *r_node;

  int i,j;
  int e_to_remove_local_index = 0;
  int e_to_remove_global_index = 0;
  int n_to_remove_global_index = 0;
  int connect_node_global_index = -1;
  int r_edge_global_index = -1;

  char **new_taxa_names;

  /**
     initialization of nodes and edge to delete
   */
  if(taxon_id>tree->nb_taxa){
    fprintf(stderr,"Warning: %s - the given taxon_id is > the number of taxa: %d\n",__FUNCTION__,taxon_id); 
    return;
  }

  for(i=0;i<tree->nb_nodes;i++){
    if(tree->a_nodes[i]->nneigh==1 && strcmp(tree->a_nodes[i]->name,tree->taxname_lookup_table[taxon_id])==0){
      n_to_remove = tree->a_nodes[i];
    }
  }

  if(n_to_remove==NULL || n_to_remove->nneigh != 1){
    fprintf(stderr,"Warning: %s() won't remove non terminal node.\n",__FUNCTION__); 
    return;
  }

  e_to_remove = n_to_remove->br[0];
  connect_node = n_to_remove->neigh[0];

  e_to_remove_global_index = e_to_remove->id;
  n_to_remove_global_index = n_to_remove->id;
  connect_node_global_index = connect_node->id;
  
  /* We get the index of the node/edge to remove*/
  for(i=0;i<connect_node->nneigh;i++){
    if(connect_node->neigh[i] == n_to_remove){
      e_to_remove_local_index = i;
    }
  }

  /**
     We remove the branch e_to_remove from the connect_node
     And the node n_to_remove from its neighbors
  */
  for(i=e_to_remove_local_index; i < connect_node->nneigh-1;i++){
    connect_node->br[i] = connect_node->br[i+1];
    connect_node->neigh[i] = connect_node->neigh[i+1];
  }
  connect_node->nneigh--;

  new_taxa_names = malloc((tree->nb_taxa-1)*sizeof(char*));

  /**
     We remove the name of the taxon from the taxa_names array 
  */
  j=0;
  for(i=0;i<tree->nb_taxa;i++){
    if(strcmp(n_to_remove->name,tree->taxa_names[i]) != 0){
      new_taxa_names[j] = strdup(tree->taxa_names[i]);
      j++;
    }
    free(tree->taxa_names[i]);
  }
  free(tree->taxa_names);
  tree->taxa_names=new_taxa_names;
  free_node(n_to_remove);
  free_edge(e_to_remove);

  tree->a_nodes[n_to_remove_global_index] = NULL;
  tree->a_edges[e_to_remove_global_index] = NULL;

  /**
     If there remains 1 neighbor, it means that connect node is the root of
     a rooted tree
     -----*r_node
     |r_edge
     *connect_node
     |e_to_remove
     -----*n_to_remove
   */
  if(connect_node->nneigh == 1){
    r_edge = connect_node->br[0];
    r_node = connect_node->neigh[0];
    r_edge_global_index = r_edge->id;
    int index = -1;
    /**
       We remove the branch r_edge from the r_node
       And the node connect_node from its neighbors
    */
    for(i=0;i<r_node->nneigh-1;i++){
      if(r_node->neigh[i] == connect_node){
	index = i;
      }
      if(index != -1){
	r_node->br[i] = r_node->br[i+1];
	r_node->neigh[i] = r_node->neigh[i+1];
      }
    }
    r_node->nneigh--;

    /* The new root is r_node*/
    if(tree->node0 == connect_node){
      tree->node0 = r_node;
    }
    free_edge(r_edge);
    free_node(connect_node);

    tree->a_nodes[connect_node_global_index] = NULL;
    tree->a_edges[r_edge_global_index] = NULL;

  } else if(connect_node->nneigh == 2){ 
    /**
       If there remains 2 neighbors to connect_node
       We connect them directly and delete connect_node
       We keep l_edge and delete r_edge
    */
    remove_single_node(tree, connect_node);
  }
  recompute_identifiers(tree);

  /**
     We update the taxname_lookup_table
  */
  for(i=0; i < tree->nb_taxa; i++){
    free(tree->taxname_lookup_table[i]);
    if(i<(tree->nb_taxa-1))
      tree->taxname_lookup_table[i] = strdup(tree->taxa_names[i]);
  }

  /**
     We update the hashtables
   */
  for(i=0;i<tree->nb_edges;i++){
    free_id_hashtable(tree->a_edges[i]->hashtbl[1]);
  }
  tree->length_hashtables = (int)((tree->nb_taxa-1) / ceil(log10((double)(tree->nb_taxa-1))));
  for(i=0;i<tree->nb_edges;i++){
    tree->a_edges[i]->hashtbl[0] = create_id_hash_table(tree->length_hashtables);
    tree->a_edges[i]->hashtbl[1] = create_id_hash_table(tree->length_hashtables);
  }
  tree->nb_taxa--;
  ntax--;
  update_hashtables_post_alltree(tree);
  update_hashtables_pre_alltree(tree);
  update_node_depths_post_alltree(tree);
  update_node_depths_pre_alltree(tree);

  /**
     now for all the branches we can delete the **left** hashtables, because the information is redundant and
     we have the equal_or_complement function to compare hashtables
  */
  for (i = 0; i < tree->nb_edges; i++) {
    free_id_hashtable(tree->a_edges[i]->hashtbl[0]); 
    tree->a_edges[i]->hashtbl[0] = NULL;
  }

  /**
     topological depths of branches
  */
  update_all_topo_depths_from_hashtables(tree);
}

/**
   This method recomputes all the identifiers 
   of the nodes and of the edges
   for which the tree->a_nodes is not null
   or tree->a_edges is not null
   It also recomputes the total number of edges 
   and nodes in the tree
 */
void recompute_identifiers(Tree *tree){
  int new_nb_edges = 0;
  int new_nb_nodes = 0;

  Node **new_nodes;
  Edge **new_edges;

  int i, j;

  for(i=0;i<tree->nb_edges;i++){
    if(tree->a_edges[i]!=NULL){
      new_nb_edges++;
    }
  }

  for(i=0;i<tree->nb_nodes;i++){
    if(tree->a_nodes[i]!=NULL){
      new_nb_nodes++;
    }
  }

  /**
     We recompute all node identifiers 
  */
  new_nodes = malloc(new_nb_nodes*sizeof(Node*));
  new_edges = malloc(new_nb_edges*sizeof(Edge*));

  j=0;
  for(i=0;i<tree->nb_nodes;i++){
    if(tree->a_nodes[i]!=NULL){
      tree->a_nodes[i]->id=j;
      new_nodes[j] = tree->a_nodes[i];
      j++;
    }
  }

  /**
     We recompute all edge identifiers 
  */
  j=0;
  for(i=0;i<tree->nb_edges;i++){
    if(tree->a_edges[i] != NULL){
      tree->a_edges[i]->id=j;
      new_edges[j] = tree->a_edges[i];
      j++;
    }
  }
  free(tree->a_nodes);
  tree->a_nodes = new_nodes;
  tree->nb_nodes=new_nb_nodes;
  free(tree->a_edges);
  tree->a_edges = new_edges;
  tree->nb_edges=new_nb_edges;
}

/**
   If there remains 2 neighbors to connect_node
   We connect them directly and delete connect_node
   We keep l_edge and delete r_edge
   -> If nneigh de connect node != 2 : Do nothing
              connect_node
             l_edge   r_edge
     l_node *-------*--------* r_node
   => Careful: After this function, you may want to call 
   => recompute_identifiers()
*/
void remove_single_node(Tree *tree, Node *connect_node){

  Edge *l_edge = connect_node->br[0];
  Edge *r_edge = connect_node->br[1];
  int r_edge_global_index = r_edge->id;
  int connect_node_global_index = connect_node->id;

  Node *l_node = (l_edge->left == connect_node) ? l_edge->right : l_edge->left;
  Node *r_node = (r_edge->left == connect_node) ? r_edge->right : r_edge->left;

  Node *tmp;
  double sum_brlengths = 0;
  char * new_right_name = NULL;
  double new_branch_support = -1000;
  int i;

  if(connect_node->nneigh!=2){
    return;
  }

  new_right_name = NULL;
  for(i=0;i<connect_node->nneigh;i++){
    sum_brlengths+=connect_node->br[i]->brlen;
    if(connect_node->br[i]->has_branch_support
       && connect_node->br[i]->branch_support > new_branch_support){
      new_branch_support = connect_node->br[i]->branch_support;
      new_right_name = connect_node->br[i]->right->name;
    }
  }

  /**
     We replace connect_node by r_node from l_node neighbors 
  */
  for(i=0;i<l_node->nneigh;i++){
    if(l_node->neigh[i] == connect_node){
      l_node->neigh[i] = r_node;
    }
  }

  /**
     We replace connect_node by l_node from r_node neighbors 
  */
  for(i=0;i<r_node->nneigh;i++){
    if(r_node->neigh[i] == connect_node){
      r_node->neigh[i] = l_node;
      r_node->br[i] = l_edge;
    }
  }
    
  /**
     We replace the left or right of l_edge by r_edge 
  */
  if(l_edge->left == connect_node){
    l_edge->left = r_node;
  }else{
    l_edge->right = r_node;
  }
  
  /**
     We check that the left is not a tax node, otherwise, we swap them 
  */
  if(l_edge->left->nneigh==1){
    tmp = l_edge->left;
    l_edge->left = l_edge->right;
    l_edge->right = tmp;
  }

  l_edge->brlen = sum_brlengths;
  
  /**
     If right is a tax node, then no branch support anymore
  */
  if(l_edge->right->nneigh==1){
    l_edge->has_branch_support = 0;
    l_edge->branch_support = 0;
  }else{
    /**
       Otherwise we take the max branch_support computed earlier
    */
    l_edge->branch_support = new_branch_support;
    if(l_edge->right->name != new_right_name)
      strcpy(l_edge->right->name,new_right_name);
  }
  
  /**
     if the root was the deleted node, we take a new root
  */
  if(tree->node0 == connect_node){
    tree->node0 = l_edge->left;
    free(tree->node0->name);
    tree->node0->name = NULL;
  }

  tree->a_edges[r_edge_global_index] = NULL;
  tree->a_nodes[connect_node_global_index] = NULL;

  free_edge(r_edge);
  free_node(connect_node);
}

/**
   This function shuffles the taxa of an input tree 
   It takes also in argument an array of indices that
   will be shuffled, and will be used to shuffle taxa 
   names.
   - If the array is NULL: then it will init it with [0..nb_taxa]
   and then shuffle it. It is freed at the end
   - If the array is not NULL: it must contain all indices from 0 
   to nb_taxa (in any order), and it will be shuffled. It will not be freed

   if duplicate taxnames : then it will copy string from tax_name array to nodes->name
   else : it will just assign pointer from tax_name array to nodes->name

   --> it the last case : be careful of assign NULL to node->name after the function
   otherwise the memory from tax_name and node->name will be freed twice when free_tree will
   be applied
*/
void shuffle_taxa(Tree *tree){
  int * shuffled_indices = NULL;
  int i = 0;
  int node = 0;

  shuffled_indices = (int*) malloc(tree->nb_taxa * sizeof(int));
  for(i=0; i < tree->nb_taxa ; i++){
    shuffled_indices[i]=i;
  }

  for (i=0; i < tree->nb_nodes; i++) {
    if (tree->a_nodes[i]->nneigh > 1) continue;
    if(tree->a_nodes[i]->name) {
      free(tree->a_nodes[i]->name); 
      tree->a_nodes[i]->name = NULL; 
    }
  } /* end freeing all leaf names */
  
  shuffle(shuffled_indices,tree->nb_taxa, sizeof(int));
  /* and then we change accordingly all the pointers node->name for the leaves of the tree */
  node = 0;
  for (i=0; i < tree->nb_nodes; i++){
    if (tree->a_nodes[i]->nneigh == 1){
      /* if(input_tree->a_nodes[i]->name) { free(input_tree->a_nodes[i]->name); input_tree->a_nodes[i]->name = NULL; } */
      tree->a_nodes[i]->name = strdup(tree->taxa_names[shuffled_indices[node]]);
      node++;
    }
  }

  /**
     We update the hashtables
   */
  for(i=0;i<tree->nb_edges;i++){
    free_id_hashtable(tree->a_edges[i]->hashtbl[1]);
  }
  for(i=0;i<tree->nb_edges;i++){
    tree->a_edges[i]->hashtbl[0] = create_id_hash_table(tree->length_hashtables);
    tree->a_edges[i]->hashtbl[1] = create_id_hash_table(tree->length_hashtables);
  }

  update_hashtables_post_alltree(tree);
  update_hashtables_pre_alltree(tree);
  update_node_depths_post_alltree(tree);
  update_node_depths_pre_alltree(tree);
  
  /**
     now for all the branches we can delete the **left** hashtables, because the information is redundant and
     we have the equal_or_complement function to compare hashtables
  */
  for (i = 0; i < tree->nb_edges; i++) {
    free_id_hashtable(tree->a_edges[i]->hashtbl[0]); 
    tree->a_edges[i]->hashtbl[0] = NULL;
  }
  /**
     topological depths of branches
  */
  update_all_topo_depths_from_hashtables(tree);
  
  free(shuffled_indices);
}


void reroot_acceptable(Tree* t) {
	/* this function replaces t->node0 on a trifurcated node (or bigger polytomy) selected at random */
	int i, myrandom, chosen_index_in_a_nodes, nb_trifurcated = 0;
	Node *candidate, *chosen;
	/* we first create a table of all indices of the trifurcated nodes */
	int* mytable = calloc(t->nb_nodes, sizeof(int));
	for (i = 0; i < t->nb_nodes; i++) {
		candidate = t->a_nodes[i];
		if(candidate->nneigh >= 3) mytable[nb_trifurcated++] = i;
	}
	if(nb_trifurcated == 0) {
	  fprintf(stderr,"Warning: %s was not able to find a trifurcated node! No rerooting.\n", __FUNCTION__);
		return; }
	else {
		myrandom = rand_to(nb_trifurcated); /* between 0 and nb_trifurcated excluded */
		chosen_index_in_a_nodes = mytable[myrandom];
		chosen = t->a_nodes[chosen_index_in_a_nodes];
		t->node0 = chosen;
	}

	reorient_edges(t);
	
	free(mytable);
} /* end reroot_acceptable */


void reorient_edges(Tree *t){
  int i=0;
  for(i=0; i < t->node0->nneigh; i++)
    reorient_edges_recur(t->node0->neigh[i], t->node0, t->node0->br[i]);
}

void reorient_edges_recur(Node *n, Node *prev, Edge *e){
  int i;
  /* We reorient the edge */
  if(e->left == n && e->right == prev){
    e->left = prev;
    e->right= n;
  }else{
    assert(e->left == prev && e->right == n); /* descendant */
  }
    
  for(i = 0; i < n->nneigh ; i++){
    if(n->neigh[i] != prev){
      reorient_edges_recur(n->neigh[i], n, n->br[i]);
    }
  }
}


void unrooted_to_rooted(Tree* t) {
	/* this function takes an unrooted tree and simply roots it on node0:
	   at the end of the process, t->node0 has exactly two neighbours */
	/* it assumes there is enough space in the tree's node pointer and edge pointer arrays. */
	if (t->node0->nneigh == 2) {
	  fprintf(stderr,"Warning: %s was called on a tree that was already rooted! Nothing to do.\n", __FUNCTION__);
	  return;
	}
	Node* old_root = t->node0;
	Node* son0 = old_root->neigh[0];
	Edge* br0 = old_root->br[0];
	/* we create a new root node whose left son will be what was in dir0 from the old root, and right son will be the old root. */
	Node* new_root = new_node("root", t, 2); /* will have only two neighbours */
	t->node0 = new_root;
	

	Edge* new_br = new_edge(t); /* this branch will have length MIN_BRLEN and links the new root to the old root as its right son */
	new_br->left = new_root;
	new_br->right = old_root;
	new_br->brlen = MIN_BRLEN;
	new_br->had_zero_length = 1;
	new_br->has_branch_support = 0;
	/* copying hashtables */
	assert(br0->right == son0); /* descendant */
	/* the hashtable for br0 is not modified: subtree rooted on son0 remains same */
	new_br->hashtbl[1] = complement_id_hashtbl(br0->hashtbl[1], t->nb_taxa);
	/* WARNING: not dealing with subtype counts nor topological depth */

	new_root->neigh[0] = son0;
	new_root->br[0] = br0;

	new_root->neigh[1] = old_root;
	new_root->br[1] = new_br;

	assert(son0->br[0] == br0 && br0->right == son0); /* must be the case because son0 was the neighbour of the old root in direction 0 */
	son0->neigh[0] = new_root;

	br0->left = new_root;

	old_root->neigh[0] = new_root;
	old_root->br[0] = new_br;
	/* done rerooting */
}



/* THE FOLLOWING FUNCTIONS ARE USED TO BUILD A TREE FROM A STRING (PARSING) */

/* utility functions to deal with NH files */

unsigned int tell_size_of_one_tree(const char* filename) {
	/* the only purpose of this is to know about the size of a treefile (NH format) in order to save memspace in allocating the string later on */
	/* wew open and close this file independently of any other fopen */
	unsigned int mysize = 0;
	char u;
	FILE* myfile = fopen(filename, "r");
	if (myfile) {
		while ( (u = fgetc(myfile))!= ';' ) { /* termination character of the tree */
			if (u == EOF) break; /* shouldn't happen anyway */
			if (isspace(u)) continue; else mysize++;
		}
		fclose(myfile);
	} /* end if(myfile) */
	return (mysize+1);
}	


int copy_nh_stream_into_str(FILE* nh_stream, char* big_string) {
	int index_in_string = 0;
	char u;
	/* rewind(nh_stream); DO NOT go to the beginning of the stream if we want to make this flexible enough to read several trees per file */
	while ( (u = fgetc(nh_stream))!= ';' ) { /* termination character of the tree */
		if (u == EOF) { big_string[index_in_string] = '\0'; return 0; } /* error code telling that no tree has been read properly */
		if (index_in_string == MAX_TREELENGTH - 1) {
		  fprintf(stderr,"Fatal error: tree file seems too big, are you sure it is an NH tree file? Aborting.\n");
		  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
		}
		if (isspace(u)) continue;
		big_string[index_in_string++] = u; 
	}
	big_string[index_in_string++] = ';';
	big_string[index_in_string] = '\0';
	return 1; /* leaves the stream right after the terminal ';' */
} /*end copy_nh_stream_into_str */




/* actually parsing a tree */


void process_name_and_brlen(Node* son_node, Edge* edge, Tree* current_tree, char* in_str, int begin, int end) {
	/* looks into in_str[begin..end] for the branch length of the "father" edge
	   and updates the edge and node structures accordingly */
	int colon = index_toplevel_colon(in_str,begin,end);
	int closing_par = -1, opening_bracket = -1;
	int i, ignore_mode, name_begin, name_end, name_length, effective_length;
	double brlen = .0;

	/* processing the optional BRANCH LENGTH... */
	if (colon == -1) {
		edge->had_zero_length = TRUE;
		edge->brlen = MIN_BRLEN;
	} else {
		parse_double(in_str,colon+1,end,&brlen);
		edge->had_zero_length = (brlen == 0.0);
		edge->brlen = (brlen < MIN_BRLEN ? MIN_BRLEN : brlen);
	}

			
	/* then scan backwards from the colon (or from the end if no branch length) to get the NODE NAME,
	   not going further than the first closing par */
	/* we ignore the NHX-style comments for the moment, hence the detection of the brackets, which can contain anything but nested brackets */
	ignore_mode = 0;
	for (i = (colon == -1 ? end : colon - 1); i >= begin; i--) {
		if (in_str[i] == ']' && ignore_mode == 0) { ignore_mode = 1; } 
		else if (in_str[i] == ')' && ignore_mode == 0) { closing_par = i; break; }
		else if (in_str[i] == '[' && ignore_mode) { ignore_mode = 0; opening_bracket = i; }
	} /* endfor */

	name_begin = (closing_par == -1 ? begin : closing_par + 1);
	if (opening_bracket != -1) name_end = opening_bracket - 1; else name_end = (colon == -1 ? end : colon - 1);
	/* but now if the name starts and ends with single or double quotes, remove them */
	if (in_str[name_begin] == in_str[name_end] && ( in_str[name_begin] == '"' || in_str[name_begin] == '\'' )) { name_begin++; name_end--; }
	name_length = name_end - name_begin + 1;
	effective_length = (name_length > MAX_NAMELENGTH ? MAX_NAMELENGTH : name_length);
	if (name_length >= 1) {
		son_node->name = (char*) malloc((effective_length+1) * sizeof(char));
		strncpy(son_node->name, in_str+name_begin, effective_length);
		son_node->name[effective_length] = '\0'; /* terminating the string */
	}


} /* end of process_name_and_brlen */




Node* create_son_and_connect_to_father(Node* current_node, Tree* current_tree, int direction, char* in_str, int begin, int end) {
	/* This function creates (allocates) the son node in the given direction from the current node.
	   It also creates a new branch to connect the son to the father.
	   The array structures in the tree (a_nodes and a_edges) are updated accordingly.
	   Branch length and node name are processed.
	   The input string given between the begin and end indices (included) is of the type:
	   (...)node_name:length
	   OR
	   leaf_name:length
	   OR
	   a:1,b:0.31,c:1.03
	   In both cases the length is optional, and replaced by MIN_BR_LENGTH if absent. */

	if (direction < 0) {
	  fprintf(stderr,"Error in the direction given to create a son! Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	int i;
	Node* son = (Node*) malloc(sizeof(Node));
	son->id = current_tree->next_avail_node_id++;
	current_tree->a_nodes[son->id] = son;
	current_tree->nb_nodes++;

	son->name = son->comment = NULL;
	son->depth = MAX_NODE_DEPTH;

	Edge* edge = (Edge*) malloc(sizeof(Edge));
	edge->id = current_tree->next_avail_edge_id++;
	current_tree->a_edges[edge->id] = edge;
	current_tree->nb_edges++;

	edge->hashtbl[0] = create_id_hash_table(current_tree->length_hashtables);
	edge->hashtbl[1] = create_id_hash_table(current_tree->length_hashtables);

	// for (i=0; i<2; i++) edge->subtype_counts[i] = (int*) calloc(NUM_SUBTYPES, sizeof(int));
	for (i=0; i<2; i++) edge->subtype_counts[i] = NULL; /* subtypes.c will have to create that space */

	edge->right = son;
	edge->left = current_node;

	edge->has_branch_support = 0;

	current_node->neigh[direction] = son;
	current_node->br[direction] = edge; 

	/* process node name (of the son) and branch length (of the edge we just created)... */
	process_name_and_brlen(son, edge, current_tree, in_str, begin, end);

	return son;
} /* end of create_son_and_connect_to_father */



void parse_substring_into_node(char* in_str, int begin, int end, Node* current_node, int has_father, Tree* current_tree) {
	/* this function supposes that current_node is already allocated, but not the data structures in there.
	   It reads starting from character of in_str at index begin and stops at character at index end.
	   It is supposed that the input to this function is what has been seen immediately within a set of parentheses.
	   The outer parentheses themselves are not included in the range [begin, end].
	   So we expect in_str[begin, end] to contain something like:
	   MyTaxa:1.2e-3
	   OR
	   (A:4,B:6)Archae:0.45,Ctax:0.004
	   OR
	   MyTaxa
	   OR
	   (A:4,B:6),Ctax
	   OR
	   A,B,C,D,E,etc (we allow large multifurcations, with no limit on the number of sons)
	 */

	/* When called, the current node has just been created but doesn't know yet its number of neighbours. We are going to discover
	   this when counting the number of outer commas in the substring. This function:
	   (1) checks how many outer commas are here: this is the number of "sons" of this node. Add one to it if the node has a father.
	   (2) creates the stuctures (array of node pointers and array of edge pointers) accordingly (+1 for the father)
	   (3) fills them. index 0 corresponds to the "father", the other to the "sons". */

	if (begin>end) {
	  fprintf(stderr,"Error in parse_substring_into_node: begin > end. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	int i;
	int pair[2]; /* to be the beginning and end points of the substrings describing the various nodes */
	int inner_pair[2]; /* to be the beginning and end points of the substrings after removing the name and branch length */
	int nb_commas = count_outer_commas(in_str, begin, end);
	int comma_index = begin - 1;
	int direction;
	Node* son;

	/* allocating the data structures for the current node */
	current_node->nneigh = (nb_commas==0 ? 1 : nb_commas + 1 + has_father);
	current_node->neigh = malloc(current_node->nneigh * sizeof(Node*));
	current_node->br = malloc(current_node->nneigh * sizeof(Edge*));
	
	if (nb_commas == 0) { /* leaf: no recursive call */
		/* this means there is no split here, terminal node: we know that the current node is a leaf.
		   Its name is already there in node->name, we just have to update the taxname table and all info related
		   to the fact that we have a taxon here. */
		/* that's also the moment when we check that there are no two identical taxa on different leaves of the tree */
		for(i=0;i < current_tree->next_avail_taxon_id; i++) {
			if (!strcmp(current_node->name, current_tree->taxa_names[i])) {
			  fprintf(stderr,"Fatal error: duplicate taxon %s.\n", current_node->name);
			  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
				} /* end if */
		} /* end for */

		current_tree->taxa_names[current_tree->next_avail_taxon_id++] = strdup(current_node->name);

	} else { /* at least one comma, so at least two sons: */
		for (i=0; i <= nb_commas; i++) { /* e.g. three iterations for two commas */
			direction = i + has_father;
			pair[0] = comma_index + 1; /* == begin at first iteration */
			comma_index = (i == nb_commas ? end + 1 : index_next_toplevel_comma(in_str, pair[0], end));
			pair[1] = comma_index - 1;

			son = create_son_and_connect_to_father(current_node, current_tree, direction /* dir from current */,
								in_str, pair[0], pair[1]);
			/* RECURSIVE TREATMENT OF THE SON */
			strip_toplevel_parentheses(in_str,pair[0],pair[1],inner_pair); /* because name and brlen already processed by create_son */
			parse_substring_into_node(in_str,inner_pair[0],inner_pair[1], son, 1, current_tree); /* recursive treatment */
			/* after the recursive treatment of the son, the data structures of the son have been created, so now we can write
			   in it the data corresponding to its direction0 (father) */
			son->neigh[0] = current_node;
			son->br[0] = current_node->br[direction];
		} /* end for i (treatment of the various sons) */


	} /* end if/else on the number of commas */


} /* end parse_substring_into_node */



Tree* parse_nh_string(char* in_str) {
	/* this function allocates, populates and returns a new tree. */
	/* returns NULL if the file doesn't correspond to NH format */
	int in_length = (int) strlen(in_str);
	int i; /* loop counter */
	int begin, end; /* to delimitate the string to further process */
	int n_otu = 0;

	/* SYNTACTIC CHECKS on the input string */ 	
	i = 0; while (isspace(in_str[i])) i++;
	if (in_str[i] != '(') { fprintf(stderr,"Error: tree doesn't start with an opening parenthesis.\n"); return NULL; }
	else begin = i+1;
	/* begin: AFTER the very first parenthesis */

	i = in_length-1;
	while (isspace(in_str[i])) i--;
	if (in_str[i] != ';') { fprintf(stderr,"Error: tree doesn't end with a semicolon.\n"); return NULL; }
	while (in_str[--i] != ')') ;
	end = i-1;
	/* end: BEFORE the very last parenthesis, discarding optional name for the root and uncanny branch length for its "father" branch */

	/* we make a first pass on the string to discover the number of taxa. */
	/* there are as many OTUs as commas plus 1 in the nh string */
	for (i = 0; i < in_length; i++) if (in_str[i] == ',') n_otu++;
	n_otu++;

	/* immediately, we set the global variable ntax. TODO: see if we can simply get rid of this global var. */
	ntax = n_otu;



	/************************************
	initialisation of the tree structure 
	*************************************/
	Tree *t = (Tree *) malloc(sizeof(Tree));
	/* in a rooted binary tree with n taxa, (2n-2) branches and (2n-1) nodes in total.
	  this is the maximum we can have. multifurcations will reduce the number of nodes and branches, so set the data structures to the max size */
	t->nb_taxa = n_otu;

	t->a_nodes = (Node**) calloc(2*n_otu-1, sizeof(Node*));
	t->nb_nodes = 1; /* for the moment we only have the node0 node. */

	t->a_edges = (Edge**) calloc(2*n_otu-2, sizeof(Edge*));
	t->nb_edges = 0; /* none at the moment */
	
	t->node0 = (Node*) malloc(sizeof(Node));
	t->a_nodes[0] = t->node0;

	t->node0->id = 0;
	t->node0->name = NULL;
	t->node0->comment = NULL;

	t->node0->depth = MAX_NODE_DEPTH;
	t->taxa_names = (char**) malloc(n_otu * sizeof(char*));
	t->length_hashtables = (int) (n_otu / ceil(log10((double)n_otu)));

	t->taxname_lookup_table = NULL;

	t->next_avail_node_id = 1; /* root node has id 0 */
	t->next_avail_edge_id = 0; /* no branch added so far */
	t->next_avail_taxon_id = 0; /* no taxon added so far */

	/* ACTUALLY READING THE TREE... */

	parse_substring_into_node(in_str, begin, end, t->node0, 0 /* no father node */, t);

	/* SANITY CHECKS AFTER READING THE TREE */

	//printf("\n*** BASIC STATISTICS ***\n\n", in_str);
	//printf("Number of taxa in the tree read: %d\n", t->nb_taxa);
	//printf("Number of nodes in the tree read: %d\n", t->nb_nodes);
	//printf("Next available node id in the new tree: %d\n", t->next_avail_node_id);
	//printf("Number of edges in the tree read: %d\n", t->nb_edges);
	//printf("Next available edge id in the new tree: %d\n\n", t->next_avail_edge_id);
	//printf("Number of leaves according to the tree structure: %d\n", count_leaves(t));
	//printf("Number of roots in the whole tree (must be 1): %d\n", count_roots(t));
	//printf("Number of edges with zero length: %d\n", count_zero_length_branches(t));

	/* DEBUG printf("Array of node pointers:\n");
	for(i=0; i<t->nb_nodes; i++) printf("%p\t",t->a_nodes[i]); printf("\n");
	printf("Node names:\n");
	for(i=0; i<t->nb_nodes; i++) printf("%s\n",t->a_nodes[i]->name);
	*/

	return t;

} /* end parse_nh_string */


Tree *complete_parse_nh(char* big_string, char*** taxname_lookup_table) {
	/* trick: iff taxname_lookup_table is NULL, we set it according to the tree read, otherwise we use it as the reference taxname lookup table */
	int i;
 	Tree* mytree = parse_nh_string(big_string); 
	if(mytree == NULL) { fprintf(stderr,"Not a syntactically correct NH tree.\n"); return NULL; }

	if(*taxname_lookup_table == NULL)  *taxname_lookup_table = build_taxname_lookup_table(mytree);
	mytree->taxname_lookup_table = *taxname_lookup_table;

	update_bootstrap_supports_from_node_names(mytree);
	/* update_subtype_counts_post_alltree(mytree);
	update_subtype_counts_pre_alltree(mytree);
	update_branch_subtype_counts_from_nodes(mytree); */

	update_hashtables_post_alltree(mytree);
	update_hashtables_pre_alltree(mytree);

	update_node_depths_post_alltree(mytree);
	update_node_depths_pre_alltree(mytree);

	/* for all branches in the tree, we should assert that the sum of the number of taxa on the left
	   and on the right of the branch is equal to tree->nb_taxa */
	for (i = 0; i < mytree->nb_edges; i++)
		if(!mytree->a_edges[i]->had_zero_length)
			assert(mytree->a_edges[i]->hashtbl[0]->num_items
			     + mytree->a_edges[i]->hashtbl[1]->num_items
			    == mytree->nb_taxa);


	/* now for all the branches we can delete the **left** hashtables, because the information is redundant and
	   we have the equal_or_complement function to compare hashtables */

	for (i = 0; i < mytree->nb_edges; i++) {
		free_id_hashtable(mytree->a_edges[i]->hashtbl[0]); 
		mytree->a_edges[i]->hashtbl[0] = NULL;
	}

	/* topological depths of branches */
	update_all_topo_depths_from_hashtables(mytree);

	return mytree;
}




/* taxname lookup table functions */

char** build_taxname_lookup_table(Tree* tree) {
	/* this function ALLOCATES a lookup table, a mere array of strings */
	/* lookup tables are shared between trees, be able to compare hashtables (one taxon == one index in the lookup table) */
	int i;
	char** output = (char**) malloc(tree->nb_taxa * sizeof(char*));
	for(i=0; i < tree->nb_taxa; i++) output[i] = strdup(tree->taxa_names[i]);
	return output;
}

/**
   The tax_id_lookup table is useful to make the correspondance between a
   node id and a taxon id: Is avoids to look for taxon name in the lookup_table 
   which is very time consuming: traverse the whole array and compare strings
   This structure is not stored in the tree, but may be computed when needed
*/
map_t build_taxid_hashmap(char** taxname_lookup_table, int nb_taxa){
  map_t h = hashmap_new();

  int i;
  for(i=0;i<nb_taxa;i++){
    int *val = malloc(sizeof(int));
    *val=i;
    hashmap_put(h, taxname_lookup_table[i], val);
  }
  return h;
}
void free_taxid_hashmap(map_t taxmap){
  hashmap_iterate(taxmap,&free_hashmap_data,NULL);
  hashmap_free(taxmap);
}

int free_hashmap_data(any_t arg,any_t key, any_t elemt){
  free(elemt);
  return MAP_OK;
}

char** get_taxname_lookup_table(Tree* tree) {
	return tree->taxname_lookup_table;
}


Taxon_id get_tax_id_from_tax_name(char* str, char** lookup_table, int length) {
	/* just exits on an error if the taxon is not to be found by this linear search */
	int i;
	for(i=0; i < length; i++) if (!strcmp(str,lookup_table[i])) return i;
	fprintf(stderr,"Fatal error : taxon %s not found! Aborting.\n", str);
	Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	return MAX_TAXON_ID; /* just in case the compiler would complain */
} /* end get_tax_id_from_tax_name */


/* (unnecessary/deprecated) multifurcation treatment */
void regraft_branch_on_node(Edge* edge, Node* target_node, int dir) {
	/* this function modifies the given edge and target node, but nothing concerning the hashtables, subtype_counts, etc.
	   This function is meant to be called during the tree construction or right afterwards, at a moment when the complex
	   recursive structures are not yet populated. */

	/* modifying the info into the branch */
	edge->left = target_node; /* modify the ancestor, not the descendant */
	
	/* modifying the info into the node on which we graft */
	target_node->br[dir] = edge;
	target_node->neigh[dir] = edge->right;

	/* modifying the info on the node at the right end of the branch we just grafted */
	Node* son = edge->right;
	son->neigh[0] = target_node; /* the father is always in direction 0 */

} /* end regraft_branch_on_node */



/***************************************************************
  ******************* neatly implementing tree traversals ******
***************************************************************/

/* in all cases below we accept that origin can be NULL:
   this describes the situation where we are on the pseudoroot node. */

void post_order_traversal_recur(Node* current, Node* origin, Tree* tree, void (*func)(Node*, Node*, Tree*)) {
	/* does the post order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
	int i, n = current->nneigh;
	int cur_to_orig = (origin ? dir_a_to_b(current, origin) : -1); /* direction from the current node to the origin of the traversal */

	/* process children first */
	if (cur_to_orig == -1) { /* current is the pseudoroot node */
		for(i=0; i < n; i++) post_order_traversal_recur(current->neigh[i], current, tree, func);
	} else {
		for(i=1; i < n; i++) post_order_traversal_recur(current->neigh[(cur_to_orig+i)%n], current, tree, func); /* no iter when n==1 (leaf) */
	}

	/* and then in any case, call the function on the current node */
	func(current, origin /* may be NULL, it's up to func to deal with that properly */, tree);
}

void post_order_traversal(Tree* t, void (*func)(Node*, Node*, Tree*)) {
	post_order_traversal_recur(t->node0, NULL, t, func);
}

/* Post order traversal with any data that can be passed to the recur function */
void post_order_traversal_data_recur(Node* current, Node* origin, Tree* tree, void* data, void (*func)(Node*, Node*, Tree*, void*)) {
  /* does the post order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
  int i, n = current->nneigh;
  int cur_to_orig = (origin ? dir_a_to_b(current, origin) : -1); /* direction from the current node to the origin of the traversal */

  /* process children first */
  if (cur_to_orig == -1) { /* current is the pseudoroot node */
    for(i=0; i < n; i++) post_order_traversal_data_recur(current->neigh[i], current, tree, data, func);
  } else {
    for(i=1; i < n; i++) post_order_traversal_data_recur(current->neigh[(cur_to_orig+i)%n], current, tree, data, func); /* no iter when n==1 (leaf) */
  }

  /* and then in any case, call the function on the current node */
  func(current, origin /* may be NULL, it's up to func to deal with that properly */, tree, data);
}

void post_order_traversal_data(Tree* t, void* data, void (*func)(Node*, Node*, Tree*,void*)) {
  post_order_traversal_data_recur(t->node0, NULL, t, data, func);
}

void pre_order_traversal_recur(Node* current, Node* origin, Tree* tree, void (*func)(Node*, Node*, Tree*)) {
	/* does the pre order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
	int i, n = current->nneigh;
	int cur_to_orig = (origin ? dir_a_to_b(current, origin) : -1); /* direction from the current node to the origin of the traversal */

	/* in any case, call the function on the current node first */
	func(current, origin /* may be NULL, it's up to func to deal with that properly */, tree);

	/* if current is not a leaf, process its children */
	if (cur_to_orig == -1) { /* current is the pseudoroot node */
		for(i=0; i < n; i++) pre_order_traversal_recur(current->neigh[i], current, tree, func);
	} else {
		for(i=1; i < n; i++) pre_order_traversal_recur(current->neigh[(cur_to_orig+i)%n], current, tree, func); /* no iter when n==1 (leaf) */
	}
}


void pre_order_traversal(Tree* t, void (*func)(Node*, Node*, Tree*)) {
	pre_order_traversal_recur(t->node0, NULL, t, func);
}

/* Pre order traversal with any data that can be passed to the recur function */
void pre_order_traversal_data_recur(Node* current, Node* origin, Tree* tree, void* data, void (*func)(Node*, Node*, Tree*, void*)) {
	/* does the pre order traversal on current Node and its "descendants" (i.e. not including origin, who is a neighbour of current */
	int i, n = current->nneigh;
	int cur_to_orig = (origin ? dir_a_to_b(current, origin) : -1); /* direction from the current node to the origin of the traversal */

	/* in any case, call the function on the current node first */
	func(current, origin /* may be NULL, it's up to func to deal with that properly */, tree, data);

	/* if current is not a leaf, process its children */
	if (cur_to_orig == -1) { /* current is the pseudoroot node */
	  for(i=0; i < n; i++) pre_order_traversal_data_recur(current->neigh[i], current, tree, data, func);
	} else {
	  for(i=1; i < n; i++) pre_order_traversal_data_recur(current->neigh[(cur_to_orig+i)%n], current, tree, data, func); /* no iter when n==1 (leaf) */
	}
}


void pre_order_traversal_data(Tree* t, void* data, void (*func)(Node*, Node*, Tree*, void*)) {
  pre_order_traversal_data_recur(t->node0, NULL, t, data, func);
}


/* BOOTSTRAP SUPPORT UTILITIES */

void update_bootstrap_supports_from_node_names(Tree* tree) {
	/* this calls the recursive function to update all branch bootstrap supports, originally imported as internal node names from the NH file */
	pre_order_traversal(tree,&update_bootstrap_supports_doer);
}

void update_bootstrap_supports_doer(Node* current, Node* origin, Tree* tree) {
	/* a branch takes its support value from its descendant node (son).
	   The current node under examination will give its value (node name) to its father branch, if that one exists.
	   We modify here the bootstrap support on the edge between current and origin. It is assumed that the node "origin" is on
	   the path from "current" to the (pseudo-)root */
	if(!origin || current->nneigh == 1) return; /* nothing to do for a leaf or for the root */

	double value;
	Edge* edge = current->br[dir_a_to_b(current, origin)];

	if (current->name && strlen(current->name) > 0 && sscanf(current->name,"%lf", &value) == 1) { /* if succesfully parsing a number */
		edge->has_branch_support = 1;
		edge->branch_support = value; 
	} else {
		edge->has_branch_support = 0;
	}
} /* end of update_bootstrap_supports_doer */




/* CALCULATING NODE DEPTHS */

void update_node_depths_post_doer(Node* target, Node* orig, Tree* t) {
	/* here we update the depth of the target node */
	int i;
	double depth = MAX_NODE_DEPTH;
	if (target->nneigh == 1)
	       	target->depth = 0.0;
	else {
		/* the following loop also takes care of the case where origin == NULL (target is root) */
		for (i=0; i < target->nneigh; i++) {
			if (target->neigh[i] == orig) continue;
			depth = min_double(depth, target->neigh[i]->depth + (target->br[i]->had_zero_length ? 0.0 : target->br[i]->brlen));
		}
		target->depth = depth;
	}
} /* end of update_node_depths_post_doer */


void update_node_depths_pre_doer(Node* target, Node* orig, Tree* t) {
	/* when we enter this function, orig already has its depth set to its final value. Update the target if its current depth is larger
	   than the one we get taking into account the min path to a leave from target via origin */
	if (!orig) return; /* nothing to do on the root for this preorder: value is already correctly set by the postorder */
	int dir_target_to_orig = dir_a_to_b(target, orig);
	double alt_depth = orig->depth + (target->br[dir_target_to_orig]->had_zero_length ? 0.0 : target->br[dir_target_to_orig]->brlen);
	if (alt_depth < target->depth) target->depth = alt_depth;
} /* end of update_node_depths_pre_doer */


void update_node_depths_post_alltree(Tree* tree) {
	post_order_traversal(tree, &update_node_depths_post_doer);
} /* end of update_node_depths_post_alltree */


void update_node_depths_pre_alltree(Tree* tree) {
	pre_order_traversal(tree, &update_node_depths_pre_doer);
} /* end of update_node_depths_pre_alltree */



/* working with topological depths: number of taxa on the lightest side of the branch */

void update_all_topo_depths_from_hashtables(Tree* tree) {
	int i, m, n = tree->nb_taxa;
	for (i = 0; i < tree->nb_edges; i++) {
		m = tree->a_edges[i]->hashtbl[1]->num_items;
		tree->a_edges[i]->topo_depth = min_int(m, n-m);
	}

} /* end update_all_topo_depths_from_hashtables */




int greatest_topo_depth(Tree* tree) {
	/* returns the greatest branch depth in the tree */
	int i, greatest = 0;
	for (i = 0; i < tree->nb_edges; i++) {
		if (tree->a_edges[i]->topo_depth > greatest) greatest = tree->a_edges[i]->topo_depth;
	}
	return greatest;
} /* end greatest_topo_depth */



/* WORKING WITH HASHTABLES */

void update_hashtables_post_doer(Node* current, Node* orig, Tree* t) {
	/* we are going to update one of the two hashtables sitting on the branch between current and orig. */
	if (orig==NULL) return;
	int i, n = current->nneigh;
	int curr_to_orig = dir_a_to_b(current, orig); 
	Edge* br = current->br[curr_to_orig], *br2; /* br: current to orig; br2: any _other_ branch from current */

	for(i=1 ; i < n ; i++) {
		br2 = current->br[(curr_to_orig + i)%n];
		/* we are going to update the info on br with the info from br2 */
		update_id_hashtable(br2->hashtbl[current==br2->left], /* source */
					    br->hashtbl[current==br->right]); /* dest */
	}

	/* but if n = 1 we haven't done anything (leaf): we must put the info corresponding to the taxon into the branch */
	if (n == 1) {
		assert(br->right == current);
		/* add the id of the taxon to the right hashtable of the branch */
		add_id(br->hashtbl[1],get_tax_id_from_tax_name(current->name, t->taxname_lookup_table, t->nb_taxa));
	}
} /* end update_hashtables_post_doer */


void update_hashtables_pre_doer(Node* current, Node* orig, Tree* t) {
	/* we are going to update one of the two hashtables sitting on the branch between current and orig. */
	if (orig==NULL) return;
	int i, n = orig->nneigh;
	int orig_to_curr = dir_a_to_b(orig, current); 
	Edge* br = orig->br[orig_to_curr], *br2; /* br: current to orig; br2: any _other_ branch from orig */
	id_hash_table_t* hash_to_update = br->hashtbl[current==br->left];

	/* if current is a leaf we just put in the left hashtable the full hashtable minus the taxon on the leaf */
	if (current->nneigh == 1) {
		assert(current == br->right); /* leaf should be on the right of the branch */
		//fill_id_hashtable(hash_to_update, t->nb_taxa);	
		//delete_id(hash_to_update, get_tax_id_from_tax_name(current->name, t->taxname_lookup_table, t->nb_taxa));
		complement_id_hashtable(hash_to_update /*dest*/, br->hashtbl[1] /*source*/, t->nb_taxa);
		return;
	}

	/* else we are going to update that hashtable with the info from the _other_ neighbours of the origin node. Origin can never be a leaf. */
	for(i=1 ; i < n ; i++) {
		br2 = orig->br[(orig_to_curr + i)%n];
		/* we are going to update the info on br with the info from br2 */
		update_id_hashtable(br2->hashtbl[orig==br2->left], /* source */
					    hash_to_update); /* dest */
	}
} /* end update_hashtables_pre_doer */


void update_hashtables_post_alltree(Tree* tree) {
	post_order_traversal(tree, &update_hashtables_post_doer);
} /* end of update_hashtables_post_alltree */

void update_hashtables_pre_alltree(Tree* tree) {
	pre_order_traversal(tree, &update_hashtables_pre_doer);
} /* end of update_hashtables_pre_alltree */



/* UNION AND INTERSECT CALCULATIONS (FOR THE TRANSFER METHOD) */

void update_i_c_post_order_ref_tree(Tree* ref_tree, Node* orig, Node* target, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix) {
	/* this function does the post-order traversal (recursive from the pseudoroot to the leaves, updating knowledge for the subtrees)
	   of the reference tree, examining only leaves (terminal edges) of the bootstrap tree.
	   It sends a probe from the orig node to the target node (nodes in ref_tree), calculating I_ij and C_ij
	   (see Brehelin, Gascuel, Martin 2008). */
	int j, k, dir, orig_to_target, target_to_orig;
	Edge* my_br; /* branch of the ref tree connecting orig to target */
	int edge_id; /* its id */
	int edge_id2;

	/* we first have to determine which is the direction of the edge (orig -> target and target -> orig) */
	orig_to_target = dir_a_to_b(orig,target);
	target_to_orig = dir_a_to_b(target,orig);
	my_br = orig->br[orig_to_target];
	edge_id = my_br->id; /* all this is in ref_tree */
	assert(target==my_br->right); /* the descendant should always be the right side of the edge */

	if(target->nneigh == 1) {
		for (j=0; j < boot_tree->nb_edges; j++) { /* for all the terminal edges of boot_tree */ 
			if(boot_tree->a_edges[j]->right->nneigh != 1) continue;
			/* we only want to scan terminal edges of boot_tree, where the right son is a leaf */
			/* else we update all the I_ij and C_ij with i = edge_id */
			if (strcmp(target->name,boot_tree->a_edges[j]->right->name)) {
				/* here the taxa are different */
				i_matrix[edge_id][j] = 0;
				c_matrix[edge_id][j] = 1;
			} else {
				/* same taxa here in T_ref and T_boot */
				i_matrix[edge_id][j] = 1;
				c_matrix[edge_id][j] = 0;
			}
		} /* end for on all edges of T_boot, for my_br being terminal */
	} else {
		/* now the case where my_br is not a terminal edge */
		/* first initialise (zero) the cells we are going to update */
		for (j=0; j < boot_tree->nb_edges; j++)
		  /**
		     We initialize the i and c matrices for the edge edge_id with :
		     * 0 for i : because afterwards we do i[edge_id] = i[edge_id] || i[edge_id2]
		     * 1 for c : because afterwards we do c[edge_id] = c[edge_id] && c[edge_id2]
		  */
		  if(boot_tree->a_edges[j]->right->nneigh == 1){
		    i_matrix[edge_id][j] = 0;
		    c_matrix[edge_id][j] = 1;
		  }

		for (k = 1; k < target->nneigh; k++) {
			dir = (target_to_orig + k) % target->nneigh; /* direction from target to one of its "sons" (== not orig) */
			update_i_c_post_order_ref_tree(ref_tree, target, target->neigh[dir], boot_tree, i_matrix, c_matrix);
			edge_id2 = target->br[dir]->id;
			for (j=0; j < boot_tree->nb_edges; j++) { /* for all the terminal edges of boot_tree */ 
				if(boot_tree->a_edges[j]->right->nneigh != 1) continue;

				i_matrix[edge_id][j] = i_matrix[edge_id][j] || i_matrix[edge_id2][j];
				/* above is an OR between two integers, result is 0 or 1 */

				c_matrix[edge_id][j] = c_matrix[edge_id][j] && c_matrix[edge_id2][j];
				/* above is an AND between two integers, result is 0 or 1 */

			} /* end for j */
		} /* end for on all edges of T_boot, for my_br being internal */

	} /* ending the case where my_br is an internal edge */
	
} /* end update_i_c_post_order_ref_tree */


void update_all_i_c_post_order_ref_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix) {
	/* this function is the first step of the union and intersection calculations */
	Node* root = ref_tree->node0;
	int i, n = root->nneigh;
	for(i=0; i<n; i++) update_i_c_post_order_ref_tree(ref_tree, root, root->neigh[i], boot_tree, i_matrix, c_matrix);
} /* end update_all_i_c_post_order_ref_tree */





void update_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, Node* orig, Node* target, short unsigned** i_matrix, short unsigned** c_matrix,
				     short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge) {
	/* here we implement the second part of the Brehelin/Gascuel/Martin algorithm:
	   post-order traversal of the bootstrap tree, and numerical recurrence. */
	/* in this function, orig and target are nodes of boot_tree (aka T_boot). */
	/* min_dist is an array whose size is equal to the number of edges in T_ref.
	   It gives for each edge of T_ref its min distance to a split in T_boot. */

	int i, j, dir, orig_to_target, target_to_orig;
	Edge* my_br; /* branch of the boot tree connecting orig to target */
	int edge_id /* its id */, edge_id2 /* id of descending branches. */;
	int N = ref_tree->nb_taxa;

	/* we first have to determine which is the direction of the edge (orig -> target and target -> orig) */
	orig_to_target = dir_a_to_b(orig,target);
	target_to_orig = dir_a_to_b(target,orig);
	my_br = orig->br[orig_to_target];
	edge_id = my_br->id; /* here this is an edge_id corresponding to T_boot */

	if(target->nneigh != 1) {
		/* because nothing to do in the case where target is a leaf: intersection and union already ok. */
		/* otherwise, keep on posttraversing in all other directions */

		/* first initialise (zero) the cells we are going to update */
		for (i=0; i < ref_tree->nb_edges; i++) i_matrix[i][edge_id] = c_matrix[i][edge_id] = 0;

		for(j=1;j<target->nneigh;j++) {
			dir = (target_to_orig + j) % target->nneigh;
			edge_id2 = target->br[dir]->id;
			update_i_c_post_order_boot_tree(ref_tree, boot_tree, target, target->neigh[dir],
							i_matrix, c_matrix, hamming, min_dist, min_dist_edge);
			for (i=0; i < ref_tree->nb_edges; i++) { /* for all the edges of ref_tree */ 
				i_matrix[i][edge_id] += i_matrix[i][edge_id2];
				c_matrix[i][edge_id] += c_matrix[i][edge_id2];
			} /* end for i */
		} 

	} /* end if target is not a leaf: the following loop is performed in all cases */

	for (i=0; i<ref_tree->nb_edges; i++) { /* for all the edges of ref_tree */ 
		/* at this point we can calculate in all cases (internal branch or not) the Hamming distance at [i][edge_id], */
		hamming[i][edge_id] = /* card of union minus card of intersection */ 
			ref_tree->a_edges[i]->hashtbl[1]->num_items /* #taxa in the cluster i of T_ref */
			+ c_matrix[i][edge_id] /* #taxa in cluster edge_id of T_boot BUT NOT in cluster i of T_ref */
			- i_matrix[i][edge_id]; /* #taxa in the intersection of the two clusters */

		/* NEW!! Let's immediately calculate the right ditance, taking into account the fact that the true disance is min (dist, N-dist) */
		if (hamming[i][edge_id] > N/2 /* floor value */) hamming[i][edge_id] = N - hamming[i][edge_id];
		

		/*   and update the min of all Hamming (TRANSFER) distances hamming[i][j] over all j */
		if (hamming[i][edge_id] < min_dist[i]){
			min_dist[i] = hamming[i][edge_id];
			min_dist_edge[i] = edge_id;
		}
			
	} /* end for on all edges of T_ref */
	
} /* end update_i_c_post_order_boot_tree */


void update_all_i_c_post_order_boot_tree(Tree* ref_tree, Tree* boot_tree, short unsigned** i_matrix, short unsigned** c_matrix,
					 short unsigned** hamming, short unsigned* min_dist, short unsigned* min_dist_edge) {
	/* this function is the second step of the union and intersection calculations */
	Node* root = boot_tree->node0;
	int i, n = root->nneigh;
	for(i=0 ; i<n ; i++) update_i_c_post_order_boot_tree(ref_tree, boot_tree, root, root->neigh[i], i_matrix, c_matrix, hamming, min_dist, min_dist_edge);

	/* and then some checks to make sure everything went ok */
	for(i=0; i<ref_tree->nb_edges; i++) {
		assert(min_dist[i] >= 0);
		if(ref_tree->a_edges[i]->right->nneigh == 1)
			assert(min_dist[i] == 0); /* any terminal edge should have an exact match in any bootstrap tree */
	}
} /* end update_all_i_c_post_order_boot_tree */




/* writing a tree to some output (stream or string) */

void write_nh_tree(Tree* tree, FILE* stream) {
	/* writing the tree from the current position in the stream */
	if (!tree) return;
	Node* node = tree->node0; /* root or pseudoroot node */
	int i, n = node->nneigh;
	putc('(', stream);
	for(i=0; i < n-1; i++) {
		write_subtree_to_stream(node->neigh[i], node, stream); /* a son */
		putc(',', stream);
	}
	write_subtree_to_stream(node->neigh[i], node, stream); /* last son */
	putc(')', stream);

	if (node->name) fprintf(stream, "%s", node->name);
	/* terminate with a semicol AND and end of line */
	putc(';', stream); putc('\n', stream);
}

/* the following function writes the subtree having root "node" and not including "node_from". */
void write_subtree_to_stream(Node* node, Node* node_from, FILE* stream) {
	int i, direction_to_exclude, n = node->nneigh;
	if (node == NULL || node_from == NULL) return;

	if(n == 1) {
		/* terminal node */
		fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->br[0]->brlen); /* distance to father */
	} else {
	        direction_to_exclude = dir_a_to_b(node, node_from);	

		putc('(', stream);
		/* we have to write (n-1) subtrees in total. The last print is not followed by a comma */
		for(i=1; i < n-1; i++) {
			write_subtree_to_stream(node->neigh[(direction_to_exclude+i) % n], node, stream); /* a son */
			putc(',', stream);
		}
		write_subtree_to_stream(node->neigh[(direction_to_exclude+i) % n], node, stream); /* last son */
		putc(')', stream);
		fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->br[0]->brlen); /* distance to father */
	}

} /* end write_subtree_to_stream */
		


/* freeing */

void free_edge(Edge* edge) {
	int i;
	if (edge == NULL) return;
	if(edge->hashtbl[0]) free_id_hashtable(edge->hashtbl[0]);
	if(edge->hashtbl[1]) free_id_hashtable(edge->hashtbl[1]);
	for (i=0; i<2; i++) if(edge->subtype_counts[i]) free(edge->subtype_counts[i]);
	
	free(edge);
}

void free_node(Node* node) {
	if (node == NULL) return;
	
	if (node->name) free(node->name);
	if (node->comment) free(node->comment);
	
	free(node->neigh);
	free(node->br);
	free(node);
}

void free_tree(Tree* tree) {
	if (tree == NULL) return;
	int i;
	for (i=0; i < tree->nb_nodes; i++) free_node(tree->a_nodes[i]);
	for (i=0; i < tree->nb_edges; i++) free_edge(tree->a_edges[i]);
	for (i=0; i < tree->nb_taxa; i++) free(tree->taxa_names[i]);

	free(tree->taxa_names);
	free(tree->a_nodes);
	free(tree->a_edges);
	free(tree);

}

Tree * gen_rand_tree(int nbr_taxa, char **taxa_names){
  int taxon;
  Tree *my_tree;
  int* indices = (int*) calloc(nbr_taxa, sizeof(int)); /* the array that we are going to shuffle around to get random order in the taxa names */
  /* zero the number of taxa inserted so far in this tree */
  int nb_inserted_taxa = 0;

  int i_edge, edge_ind;
  
  for(taxon = 0; taxon < nbr_taxa; taxon++)
    indices[taxon] = taxon;

  shuffle(indices, nbr_taxa, sizeof(int));
  
  if(taxa_names == NULL){
    taxa_names = (char**) calloc(nbr_taxa, sizeof(char*));
    for(taxon = 0; taxon < nbr_taxa; taxon++) {
      taxa_names[taxon] = (char*) calloc((int)(log10(nbr_taxa)+2), sizeof(char));
      sprintf(taxa_names[taxon],"%d",taxon+1); /* names taxa by a mere integer, starting with "1" */
    }
  }
  
  /* create a new tree */
  my_tree = new_tree(nbr_taxa, taxa_names[indices[nb_inserted_taxa++]]);
  
  /* graft the second taxon */
  graft_new_node_on_branch(NULL, my_tree, 0.5, 1.0, taxa_names[indices[nb_inserted_taxa++]]);
  
  while(nb_inserted_taxa < nbr_taxa) {
    /* select a branch at random */
    edge_ind = rand_to(my_tree->nb_edges); /* outputs something between 0 and (nb_edges) exclusive */
    graft_new_node_on_branch(my_tree->a_edges[edge_ind], my_tree, 0.5, 1.0, taxa_names[indices[nb_inserted_taxa++]]);
  } /* end looping on the taxa, tree is full */
  
  /* here we need to re-root the tree on a trifurcated node, not on a leaf, before we write it in NH format */
  reroot_acceptable(my_tree);

  for(i_edge = 0; i_edge < my_tree->nb_edges; i_edge++){
    my_tree->a_edges[i_edge]->brlen = normal(0.1, 0.05);
    if(my_tree->a_edges[i_edge]->brlen < 0)
      my_tree->a_edges[i_edge]->brlen = 0;
  }

  my_tree->length_hashtables = (int) (my_tree->nb_taxa / ceil(log10((double)my_tree->nb_taxa)));
  ntax = nbr_taxa;
  my_tree->taxname_lookup_table = build_taxname_lookup_table(my_tree);
  
  for(i_edge=0;i_edge<my_tree->nb_edges;i_edge++){
    my_tree->a_edges[i_edge]->hashtbl[0] = create_id_hash_table(my_tree->length_hashtables);
    my_tree->a_edges[i_edge]->hashtbl[1] = create_id_hash_table(my_tree->length_hashtables);
  }

  update_hashtables_post_alltree(my_tree);
  update_hashtables_pre_alltree(my_tree);
  update_node_depths_post_alltree(my_tree);
  update_node_depths_pre_alltree(my_tree);

  
  /* for all branches in the tree, we should assert that the sum of the number of taxa on the left
     and on the right of the branch is equal to tree->nb_taxa */
  for (i_edge = 0; i_edge < my_tree->nb_edges; i_edge++)
    if(!my_tree->a_edges[i_edge]->had_zero_length)
      assert(my_tree->a_edges[i_edge]->hashtbl[0]->num_items
	     + my_tree->a_edges[i_edge]->hashtbl[1]->num_items
	     == my_tree->nb_taxa);

  /* now for all the branches we can delete the **left** hashtables, because the information is redundant and
     we have the equal_or_complement function to compare hashtables */
  
  for (i_edge = 0; i_edge < my_tree->nb_edges; i_edge++) {
    free_id_hashtable(my_tree->a_edges[i_edge]->hashtbl[0]); 
    my_tree->a_edges[i_edge]->hashtbl[0] = NULL;
  }

  /* topological depths of branches */
  update_all_topo_depths_from_hashtables(my_tree);
  
  return(my_tree);
}
