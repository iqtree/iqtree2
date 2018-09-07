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

#include "tree_utils.h"

/**
   Generates a random tree based on the taxa of the tree and its lookuptable in argument 
   - Uses the same taxnames
   - Uses the taxname_lookup_id
   Advice: do a srand(time(NULL)) before calling this function a large number of times
*/
Tree* gen_random_tree(Tree *tree){
  int* indices = (int*) calloc(tree->nb_taxa, sizeof(int)); /* the array that we are going to shuffle around to get random order in the taxa names */
  int taxon;
  for(taxon = 0; taxon < tree->nb_taxa; taxon++) indices[taxon] = taxon; /* initialization */
  
  /* zero the number of taxa inserted so far in this tree */
  int nb_inserted_taxa = 0,edge_ind;
  Tree* my_tree = NULL;
  int i;
  /* shuffle the indices we are going to use to determine the names of leaves */
  shuffle(indices, tree->nb_taxa, sizeof(int));
  
  /* free the previous tree if existing */
  if(my_tree) free_tree(my_tree); 
  
  /* create a new tree */
  my_tree = new_tree(tree->nb_taxa, tree->taxa_names[indices[nb_inserted_taxa++]]);
	
  /* graft the second taxon */
  graft_new_node_on_branch(NULL, my_tree, 0.5, 1.0, tree->taxa_names[indices[nb_inserted_taxa++]]);
  
  while(nb_inserted_taxa < tree->nb_taxa) {
    /* select a branch at random */
    edge_ind = rand_to(my_tree->nb_edges); /* outputs something between 0 and (nb_edges-1) exclusive */
    graft_new_node_on_branch(my_tree->a_edges[edge_ind], my_tree, 0.5, 1.0, tree->taxa_names[indices[nb_inserted_taxa++]]);
  } /* end looping on the taxa, tree is full */

  /* here we need to re-root the tree on a trifurcated node, not on a leaf, before we write it in NH format */
  reroot_acceptable(my_tree);

  my_tree->taxname_lookup_table = tree->taxname_lookup_table;
  my_tree->nb_taxa = tree->nb_taxa;
  my_tree->length_hashtables = (int) (my_tree->nb_taxa / ceil(log10((double)my_tree->nb_taxa)));

  int e;
  for(e=0;e<my_tree->nb_edges;e++){
    my_tree->a_edges[e]->hashtbl[0] = create_id_hash_table(my_tree->length_hashtables);
    my_tree->a_edges[e]->hashtbl[1] = create_id_hash_table(my_tree->length_hashtables);
  }

  /* write_nh_tree(my_tree,stdout); */

  update_hashtables_post_alltree(my_tree);
  update_hashtables_pre_alltree(my_tree);
  update_node_depths_post_alltree(my_tree);
  update_node_depths_pre_alltree(my_tree);

  /* for all branches in the tree, we should assert that the sum of the number of taxa on the left
     and on the right of the branch is equal to tree->nb_taxa */
  for (i = 0; i < my_tree->nb_edges; i++)
    if(!my_tree->a_edges[i]->had_zero_length)
      assert(my_tree->a_edges[i]->hashtbl[0]->num_items
	     + my_tree->a_edges[i]->hashtbl[1]->num_items
	     == my_tree->nb_taxa);

  /* now for all the branches we can delete the **left** hashtables, because the information is redundant and
     we have the equal_or_complement function to compare hashtables */
  
  for (i = 0; i < my_tree->nb_edges; i++) {
    free_id_hashtable(my_tree->a_edges[i]->hashtbl[0]); 
    my_tree->a_edges[i]->hashtbl[0] = NULL;
  }

  /* topological depths of branches */
  update_all_topo_depths_from_hashtables(my_tree);

  free(indices);
  return(my_tree);
}
