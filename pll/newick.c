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
 * @file newick.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "pll.h"
#include "pllInternal.h"


/** @file  newick.c

    @brief Collection of routines for reading and parsing newick trees

    Auxiliary functions for reading and parsing newick tree formats
*/


/** @defgroup newickParseGroup Reading and parsing newick trees
    
    This set of functions handles the reading and parsing of newick tree formats
*/

static int
parse_newick (pllStack ** stack, int * inp)
{
  pllNewickNodeInfo * item = NULL;
  int item_active = 0;
  pllLexToken token;
  int input;
  pllLexToken prev_token;
  int nop = 0;          /* number of open parentheses */
  int depth = 0;

  prev_token.tokenType = PLL_TOKEN_UNKNOWN;

  input = *inp;

  NEXT_TOKEN
  
  while (token.tokenType != PLL_TOKEN_EOF && token.tokenType != PLL_TOKEN_UNKNOWN)
  {
    switch (token.tokenType)
     {
       case PLL_TOKEN_OPAREN:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_OPAREN\n");
#endif
        ++nop;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        ++depth;
        break;

       case PLL_TOKEN_CPAREN:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_CPAREN\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN  &&
            prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
            prev_token.tokenType != PLL_TOKEN_STRING  &&
            prev_token.tokenType != PLL_TOKEN_NUMBER  &&
            prev_token.tokenType != PLL_TOKEN_FLOAT) return (0);

        if (!nop) return (0);
        --nop;
        memcpy (&prev_token, &token, sizeof (pllLexToken));

        /* push to the stack */
        if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo)); // possibly not nec
        //if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("INTERNAL_NODE") + 1) * sizeof (char));
           strcpy (item->name, "INTERNAL_NODE");
         }

        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        item->depth = depth;
        pllStackPush (stack, item);
        item_active  = 1;       /* active = 1 */
        item = NULL;
        --depth;
        break;

       case PLL_TOKEN_STRING:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_STRING      %.*s\n", token.len, token.lexeme);
#endif
        if (prev_token.tokenType != PLL_TOKEN_OPAREN &&
            prev_token.tokenType != PLL_TOKEN_CPAREN &&
            prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
            prev_token.tokenType != PLL_TOKEN_COMMA) return (0);
        if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo));
        item->name = my_strndup (token.lexeme, token.len);

        item_active = 1;
        item->depth = depth;
        if (prev_token.tokenType == PLL_TOKEN_COMMA  ||
            prev_token.tokenType == PLL_TOKEN_OPAREN ||
            prev_token.tokenType == PLL_TOKEN_UNKNOWN) item->leaf = 1;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_FLOAT:
       case PLL_TOKEN_NUMBER:
#ifdef PLLDEBUG
       if (token.tokenType == PLL_TOKEN_FLOAT) printf ("PLL_TOKEN_FLOAT\n"); else printf ("PLL_TOKEN_NUMBER\n");
#endif
         if  (prev_token.tokenType != PLL_TOKEN_OPAREN &&
              prev_token.tokenType != PLL_TOKEN_CPAREN &&
              prev_token.tokenType != PLL_TOKEN_COLON  &&
              prev_token.tokenType != PLL_TOKEN_UNKNOWN &&
              prev_token.tokenType != PLL_TOKEN_COMMA) return (0);
        if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo));
        if (prev_token.tokenType == PLL_TOKEN_COLON)
         {
           item->branch = my_strndup (token.lexeme, token.len);
         }
        else
         {
           if (prev_token.tokenType == PLL_TOKEN_COMMA  ||
               prev_token.tokenType == PLL_TOKEN_OPAREN ||
               prev_token.tokenType == PLL_TOKEN_UNKNOWN) item->leaf = 1;
           //if (prev_token.tokenType != PLL_TOKEN_UNKNOWN) ++ indent;
           item->name = my_strndup (token.lexeme, token.len);
         }
        item_active = 1;
        item->depth = depth;
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_COLON:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_COLON\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN &&
            prev_token.tokenType != PLL_TOKEN_STRING &&
            prev_token.tokenType != PLL_TOKEN_FLOAT  &&
            prev_token.tokenType != PLL_TOKEN_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        break;

       case PLL_TOKEN_COMMA:
#ifdef PLLDEBUG
       printf ("PLL_TOKEN_COMMA\n");
#endif
        if (prev_token.tokenType != PLL_TOKEN_CPAREN &&
             prev_token.tokenType != PLL_TOKEN_STRING &&
             prev_token.tokenType != PLL_TOKEN_FLOAT && 
             prev_token.tokenType != PLL_TOKEN_NUMBER) return (0);
        memcpy (&prev_token, &token, sizeof (pllLexToken));
        
        /* push to the stack */
        if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo)); // possibly not nece
        //if (item->name   == NULL) item->name   = strdup ("INTERNAL_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("INTERNAL_NODE") + 1) * sizeof (char));
           strcpy (item->name, "INTERNAL_NODE");
         }
        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        item->depth = depth;
        pllStackPush (stack, item);
        item_active  = 0;
        item = NULL;
        break;

       case PLL_TOKEN_SEMICOLON:
#ifdef PLLDEBUG
        printf ("PLL_TOKEN_SEMICOLON\n");
#endif
        /* push to the stack */
        if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo));
        //if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
        if (item->name == NULL) 
         {
           item->name = (char *) rax_malloc ((strlen("ROOT_NODE") + 1) * sizeof (char));
           strcpy (item->name, "ROOT_NODE");
         }
        //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
        if (item->branch == NULL) 
         {
           item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
           strcpy (item->branch, "0.000000");
         }
        pllStackPush (stack, item);
        item_active  = 0;
        item = NULL;
        break;
       default:
#ifdef __DEBUGGING_MODE
         printf ("Unknown token: %d\n", token.tokenType);
#endif
       // TODO: Finish this part and add error codes
        break;
     }
    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE);
  }
  if (item_active)
   {
     if (!item) item = (pllNewickNodeInfo *) rax_calloc (1, sizeof (pllNewickNodeInfo));
     //if (item->name   == NULL) item->name   = strdup ("ROOT_NODE");
     if (item->name == NULL) 
      {
        item->name = (char *) rax_malloc ((strlen("ROOT_NODE") + 1) * sizeof (char));
        strcpy (item->name, "ROOT_NODE");
      }
     //if (item->branch == NULL) item->branch = strdup ("0.000000"); 
     if (item->branch == NULL) 
      {
        item->branch = (char *) rax_malloc ((strlen("0.000000") + 1) * sizeof (char));
        strcpy (item->branch, "0.000000");
      }
     pllStackPush (stack, item);
     item_active  = 0;
   }

  if (nop || token.tokenType == PLL_TOKEN_UNKNOWN) 
   {
     return (0);
   }

  return (1);
}

#ifdef __DEBUGGING_MODE
void stack_dump(pllStack ** stack)
{
  pllNewickNodeInfo * item;
  pllStack * head;
  int i;

  head = *stack;
  while (head)
   {
     item = (pllNewickNodeInfo *) head->item;

     for (i = 0; i < item->depth; ++ i) printf ("\t");

     printf ("%s:%s\n", item->name, item->branch);

     head = head->next;
   }
}
#endif

static void
assign_ranks (pllStack * stack, int * nodes, int * leaves)
{
  pllStack * head;
  pllNewickNodeInfo * item, * tmp;
  pllStack * preorder = NULL;
  int children;
  int depth;

  *nodes = *leaves = 0;


  head = stack;
  while (head)
  {
    assert (head->item);
    item = (pllNewickNodeInfo *) head->item;
    
    if (item->leaf)  ++ (*leaves);

    if (preorder)
     {
       tmp = (pllNewickNodeInfo *) preorder->item;
       children = 0;
       while (item->depth < tmp->depth)
        {
          children = 1;
          depth = tmp->depth;
          pllStackPop (&preorder);
          tmp = preorder->item;
          while (tmp->depth == depth)
           {
             ++ children;
             pllStackPop (&preorder);
             tmp = (pllNewickNodeInfo *)preorder->item;
           }
          tmp->rank += children;
        }
     }
    
    ++ (*nodes);
    head = head->next;

    if (item->leaf)
     {
       if (!preorder) return;

       children = 1;
       tmp = preorder->item;
       while (tmp->depth == item->depth)
        {
          ++ children;
          pllStackPop (&preorder);
          assert (preorder);
          tmp = (pllNewickNodeInfo *)preorder->item;
        }
       tmp->rank += children;
     }
    else
     {
       pllStackPush (&preorder, item);
     }
  }
  
  while (preorder->item != stack->item)
  {
    item = (pllNewickNodeInfo *)pllStackPop (&preorder);
    tmp  = (pllNewickNodeInfo *) preorder->item;
    children = 1;

    while (tmp->depth == item->depth)
     {
       ++ children;
       item = (pllNewickNodeInfo *) pllStackPop (&preorder);
       tmp  = (pllNewickNodeInfo *) preorder->item;
     }
    tmp->rank += children;
    children = 0;
  }
 assert (preorder->item == stack->item);
 
 pllStackClear (&preorder);
}

/** @ingroup newickParseGroup
    @brief Validate if a newick tree is a valid phylogenetic tree

    A valid tree is one where the root node is binary or ternary
    and all other internal nodes are binary. In case the root
    is ternary then the tree must contain at least another internal
    node and the total number of nodes must be equal to 
    \f$ 2l - 2\f$, where \f$l\f$ is the number of leaves. If the
    root is binary, then the total number of nodes must be equal
    to \f$2l - 1\f$.

    @param tree
      Newick tree wrapper structure which contains the stack representation of the parsed newick tree

    @return
      Returns \b 1 in case of success, otherwise \b 0
*/
int
pllValidateNewick (pllNewickTree * t)
{
  pllStack * head;
  pllNewickNodeInfo * item;
  int correct = 0;
 
  item = t->tree->item;
  if (item->rank != 2 && item->rank != 3) return (0);
  head = t->tree->next;
  while (head)
  {
    item = head->item;
    if (item->rank != 2 && item->rank != 0) 
     {
       return (0);
     }
    head = head->next;
  }
  
  item = t->tree->item;

  if (item->rank == 2) 
   {
     correct = (t->nodes == 2 * t->tips -1);
     if (correct)
      {
        errno = PLL_NEWICK_ROOTED_TREE;
      }
     else
      {
        errno = PLL_NEWICK_BAD_STRUCTURE;
      }
     return (PLL_FALSE);
   }
   
  
  correct = ((t->nodes == 2 * t->tips - 2) && t->nodes != 4);
  if (correct) return (PLL_TRUE);

  errno = PLL_NEWICK_BAD_STRUCTURE;

  return (1);
}


/** @ingroup newickParseGroup
    @brief Convert a binary rooted trree to a binary unrooted tree

    Changes the root of the node to have 3 descendants instead of two, deletes its last immediate descendant internal node
    and takes the two children (of the deleted internal node) as its children.

    @param
      Newick tree
    
    @return
      \b PLL_TRUE in case of success, otherwise \b PLL_FALSE and \a errno is set
*/
int
pllNewickUnroot (pllNewickTree * t)
{
  pllStack * tmp;
  pllNewickNodeInfo * item;

  item = t->tree->item;
  if (item->rank == 2)
   {
     item->rank = 3;
     t->nodes--;
     item = t->tree->next->item;
     if (item->rank == 0)
      {
        tmp = t->tree->next->next;
        t->tree->next->next = t->tree->next->next->next;
      }
     else
      {
        tmp = t->tree->next;
        t->tree->next = t->tree->next->next;
      }
     item = tmp->item;
     rax_free (item->name);
     rax_free (tmp->item);
     rax_free (tmp);
   }

  return (pllValidateNewick (t));
}


/** @ingroup newickParseGroup
    @brief Parse a newick tree string
  
    Parse a newick string and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children and depth. The
    stack structure is finally wrapped in a \a pllNewickTree structure which
    also contains the number of nodes and leaves.

    @param newick
      String containing the newick tree

    @return
      Returns a pointer to the created \a pllNewickTree structure in case of success, otherwise \b NULL
*/
pllNewickTree *
pllNewickParseString (const char * newick)
{
  int n, input, rc;
  pllNewickTree * t;
  int nodes, leaves;
  
  t = (pllNewickTree *) rax_calloc (1, sizeof (pllNewickTree));

  n = strlen (newick);

  init_lexan (newick, n);
  input = get_next_symbol();

  rc = parse_newick (&(t->tree), &input);
  if (!rc)
   {
     /* TODO: properly clean t->tree */
     rax_free (t);
     t = NULL;
   }
  else
   {
     assign_ranks (t->tree, &nodes, &leaves);
     t->nodes = nodes;
     t->tips  = leaves;
   }

  return (t);
}

/** @ingroup newickParseGroup
    @brief Deallocate newick parser stack structure

    Deallocates the newick parser stack structure that represents the parsed tree. It
    also frees all memory allocated by elements of the stack structure.

    @param tree
      The tree stack structure
*/
void pllNewickParseDestroy (pllNewickTree ** t)
{
  pllNewickNodeInfo *  item;

  while ((item = (pllNewickNodeInfo *)pllStackPop (&((*t)->tree))))
   {
     rax_free (item->name);
     rax_free (item->branch);
     rax_free (item);
   }
  rax_free (*t);
  (*t) = NULL;
}

/** @ingroup newickParseGroup
    @brief Parse a newick tree file
  
    Parse a newick file and create a stack structure which represents the tree
    in a preorder traversal form. Each element of the stack represents one node
    and consists of its name, branch length, number of children (rank) and depth. The
    stack structure is finally wrapped in a \a pllNewickTree structure which
    also contains the number of nodes and leaves.

    @param filename
      Filename containing the newick tree

    @return
      Returns a pointer to the created \a pllNewickTree structure in case of success, otherwise \b NULL
*/
pllNewickTree *
pllNewickParseFile (const char * filename)
{
  long n;
  char * rawdata;
  pllNewickTree * t;

  rawdata = pllReadFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  //printf ("%s\n\n", rawdata);

  t = pllNewickParseString (rawdata);

  rax_free (rawdata);

  return (t);
}

