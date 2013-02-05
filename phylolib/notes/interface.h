#ifndef INTERFACE_H
#define INTERFACE_H

// This is from a non-updated older brainstoring

struct phylip_data
 {
   int          ntaxa;
   int          len;
   char      ** label;
   char      ** seq;
 };

struct phylip_data *    INIT_msa_read_phylip (const char * phylip_file);
struct phylip_data *    INIT_msa_read_binary (const char * binary_file);

/* Trese functions create the tree structure. They also set the rate heterogenity
   to a default (GAMMA) */
tree *                  INIT_tree_create(int mxtips, int type);   // random or parsimony
tree *                  INIT_tree_create_newick (const char * newick);

/* init the partition data and the model for each partition, either from a file
   or from a string. */
void                    INIT_partition_string (tree * t, const char * str);
void                    INIT_partition_file   (tree * t, const char * partition_file);  /* wrapper function of _string */


/* link the MSA to the tree. Either link the pointers (deep_copy = TRUE) or 
   copy the data (deep_copy = FALSE)*/

CFG_set_msa (tree *, struct phylip_data *, int deep_copy = TRUE);


#endif
