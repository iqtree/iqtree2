#ifndef __pll_NEWICK__
#define __pll_NEWICK__
#include "stack.h"
/** @brief Intermediate structure for storing a newick tree 
    
    Holds the structure of a parsed newick tree. The number of inner nodes is stored in \a nodes
*/
typedef struct
{
  int nodes;                    /**< @brief Total number of nodes in the tree == 2*tips - 1 for rooted and 2*tips -2 for unrooted */
  int tips;                     /**< @brief Number of leaves (tips) in the tree */
  pllStack * tree;              /**< @brief Parsed tree represented as elements of a stack. Corresponds to placing the postorder traversal of a rooted tree in a pushdown store */
} pllNewickTree;


/** @brief Information describing the parsed newick tree nodes 
    
    This structure is placed in the ::pllNewickTree LIFO element pllNewickTree::tree
    and described each node of the parsed tree.

    @todo Rename this to something more proper
*/
typedef struct
{
  int depth;                    /**< @brief Distance of node from root */
  char * name;                  /**< @brief Name of the taxon represented by the node (in case it is a leaf) */
  char * branch;                /**< @brief Length of branch that leads to its parent */
  int leaf;                     /**< @brief \b PLL_TRUE if the node is a leaf, otherwise \b PLL_FALSE */
  int rank;                     /**< @brief Rank of the node, i.e. how many children it has */
} pllNewickNodeInfo;


#endif
