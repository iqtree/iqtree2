#include <stdio.h>
#include <stdlib.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"



void read_msa (tree * tr, char * filename);

void makeParsimonyTree(tree *tr)
{
  allocateParsimonyDataStructures(tr);
  makeParsimonyTreeFast(tr);
  freeParsimonyDataStructures(tr);
}

/* small example program that executes ancestral state computations 
 *    on the entire subtree rooted at p.
 *
 *       Note that this is a post-order traversal.
 *       */
static void computeAllAncestralVectors(nodeptr p, tree *tr)
{
  /* if this is not a tip, for which evidently it does not make sense 
   *      to compute the ancestral sequence because we have the real one ....
   *        */

  if(!isTip(p->number, tr->mxtips))
  {
    /* descend recursively to compute the ancestral states in the left and right subtrees */

    computeAllAncestralVectors(p->next->back, tr);
    computeAllAncestralVectors(p->next->next->back, tr);

    /* then compute the ancestral state at node p */

    newviewGenericAncestral(tr, p);

    /* and print it to terminal, the two booleans that are set to true here 
     *    tell the function to print the marginal probabilities as well as 
     *       a discrete inner sequence, that is, ACGT etc., always selecting and printing 
     *          the state that has the highest probability */

    printAncestralState(p, TRUE, TRUE, tr);
  }
}




int main(int argc, char * argv[])
{

  tree        * tr;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s [binary-alignment-file]\n", argv[0]);
     return (1);
   }
  tr = (tree *)malloc(sizeof(tree));

  /* read the binary input, setup tree, initialize model with alignment */
  read_msa(tr,argv[1]);
  tr->randomNumberSeed = 665;
  printf("Number of taxa: %d\n", tr->mxtips);
  printf("Number of partitions: %d\n", tr->NumberOfModels);

  /* compute the LH of the full tree */
  makeParsimonyTree(tr);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood of parsimony starting tree: %f\n", tr->likelihood);

  /* 8 rounds of branch length optimization */
  smoothTree(tr, 8);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood after branch length optimization: %f\n", tr->likelihood);

  /* Model optimization */
  modOpt(tr, 10.0);
  evaluateGeneric(tr, tr->start, TRUE);
  printf("Likelihood after model optimization: %f\n", tr->likelihood);

  /* compute all ancestral probability vectors for the inner node 
   *      that is connected to the first taxon in the input alignment,
   *           i.e., tr->nodep[1]->back
   *             */
  computeAllAncestralVectors(tr->nodep[1]->back, tr);

  return (0);
}
