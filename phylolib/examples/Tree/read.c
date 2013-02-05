#include <stdio.h>
#include <stdlib.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"
#include "phylip_parser/phylip.h"

int main(int argc, char * argv[])
{
  tree        * tr;

  if (argc != 2)
   {
     printf (" %s [PHYLIP-FILE]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* Allocate a tree structure */
  tr = (tree *)malloc(sizeof(tree));

  /* Set the thread ID... This is PTHREADS specific */
  tr->threadID = 0;

  /* Read a sequential phylip format containing DNA data */
  read_phylip_msa ( tr, argv[1], PHYLIP_SEQUENTIAL, DNA_DATA);

  /* Generate a random tree according to the seed given in tr->randomNumberSeed */
  tr->randomNumberSeed = 3456;
  makeRandomTree(tr);

  /* Print the tree */
  printTopology(tr, FALSE);   /* do not print branch lengths */
  printTopology(tr, TRUE);    /* print branch lengths */


  return(EXIT_SUCCESS);
}
