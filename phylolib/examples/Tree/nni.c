#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"

void getStartingTree(tree *tr);
void makeRandomTree(tree *tr);

#ifdef __cplusplus
extern "C" {
boolean setupTree (tree *tr);
nodeptr pickRandomSubtree(tree *tr);
}
#else
boolean setupTree (tree *tr);
#endif

int  printBranchLengths;

void do_NNI(tree * tr, nodeptr p, int swap)
{
  nodeptr       q;
  nodeptr       tmp;

  q = p->back;
  assert(!isTip(q->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));


  if(swap == 1)
   {
     tmp = p->next->back;
     hookup(p->next, q->next->back, q->next->z, tr->numBranches);
     hookup(q->next, tmp,           tmp->z, tr->numBranches);
   }
  else
   {
      tmp = p->next->next->back;
      hookup(p->next->next, q->next->back, q->next->z,       tr->numBranches);
      hookup(q->next,       tmp,           tmp->z, tr->numBranches);
   }
}


int main(int argc, char * argv[])
{
  tree        * tr;
  nodeptr       p;

  
  srand(time(NULL));
  
  /* Set the minimum required info for the tree structure */
  tr = (tree *)malloc(sizeof(tree));
  tr->mxtips           = 5 + rand() % 5 ;
  tr->randomNumberSeed = rand();
  tr->NumberOfModels   = 1;
  tr->numBranches      = 1;

  /* Setup some default values 
     TODO: The minimal initialization can be substantially smaller than what is
     described in axml.c 
  */
  setupTree(tr);

  /* Generate a random tree according to the seed given in tr->randomNumberSeed */
  makeRandomTree(tr);
  printf ( "Working with random topology of %d taxa\n", tr->mxtips );

  printf ( "Newick notation BEFORE NNI\n" );
  printTopology(tr, FALSE);
  printTopology(tr, TRUE);
  
  do
   {
     p = pickRandomSubtree(tr);
   } while (isTip(p->back->number, tr->mxtips));
  /* perform the NNI move */
  do_NNI(tr, p, 1);

  printf ( "Newick notation AFTER NNI\n" );
  printTopology(tr, FALSE);     /* do not print branch lengths */
  printTopology(tr, TRUE);      /* print branch lengths */

  

  return(EXIT_SUCCESS);
}


