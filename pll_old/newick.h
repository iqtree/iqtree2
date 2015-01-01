/** 
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Derived from 
 * RAxML-HPC, a program for sequential and parallel estimation of phylogenetic
 * trees by Alexandros Stamatakis
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Tomas Flouri
 * Tomas.Flouri@h-its.org
 *
 * When publishing work that uses PLL please cite PLL
 * 
 * @file newick.h
 */
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
