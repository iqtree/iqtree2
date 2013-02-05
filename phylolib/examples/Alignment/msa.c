#include <stdio.h>
#include <stdlib.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"



void read_msa (tree * tr, char * filename);

int main(int argc, char * argv[])
{

  tree        * tr;

  if (argc != 2)
   {
     fprintf (stderr, "syntax: %s [binary-alignment-file]\n", argv[0]);
     return (1);
   }
  tr = (tree *)malloc(sizeof(tree));

  read_msa(tr,argv[1]);
  
  printf("Number of taxa: %d\n", tr->mxtips);
  printf("Number of partitions: %d\n", tr->NumberOfModels);
  
  return (0);
}
