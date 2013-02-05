#ifndef SMSA_SPECS
#define SMSA_SPECS

struct smsa_params
 {
   double       basefreqs[4];         /* prior probabilities for the four base types */
   double       lambda;               /* (indel rate) >= 0 */
   double       gamma;                /* expected fragment length >= 1 */
   int          sample_size;          /* desired number of samples, input by user. default = 100,000 */
   int          fragsize;             /* number of sites to be used during resampling. Typical values [5,6,...,25]. Computation grows exponentially */
   double       branch_length[6];     /* branch length of each edge of the tree given a root node. Take the 4 leaves around the node and the two internal nodes. The order is A,B,C,D,E,F    (A,B,C,D are the leaves) */
   char      ** seqs;                 /* Unaligned sequences. Extend it to a structure in order to save the label of each sequence as well  */
   char      ** alignment;            /* Extend it to a structure */

   double       subst_matrix[6][4][4];         /* Given a substitution model, compute the substitution matrix for each of the 6 branch lengths.
                                                  We should change this to [6][16] */
 };

struct smsa_params * pl_tree_smsaparams (nodeptr base_node, int sub_model);

struct FPTableNode
 {
   /* 4 * N  table, where $N$ is the number of sites */
 };

/* given a node, store the three fptables of the three adjacent nodes */
struct FPTable
 {
   struct FPTableNode  * ft1;
   struct FPTableNode  * ft2;
   struct FPTableNode  * ft3;
 }

/* Felsenstein pruning tables (FPTables) 
    
   Given an internal node, return the three FPTables for the three adjacent nodes.
*/
struct FPTable * pl_tree_fptable (nodeptr base_node);
 

/* Homology structures */

/* Overall homology table */

struct hom_table
 {
   /*  6 x N, where $N$ is the number of sites in the MSA */

 };

struct hom_table_node
 {
     /* 6 x M, where M is the number of bases in sequence of taxon of current node
        It is created based on the overall homology table, by eliminating all columns
        which contain -1 for the current node */
 };

struct hom_table * pl_tree_homtable (nodeptr base_node);
struct hom_table_node * pl_node_homtable (nodeptr base_node, nodeptr node);

#endif
