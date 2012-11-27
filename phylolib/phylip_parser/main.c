#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"
#include "xalloc.h"
#include "msa_sites.h"

int 
main (int argc, char * argv[])
{
  struct phylip_data * pd;
  int format, i;
  int * weight;

  if (argc != 3)
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }
  
  if (strcmp (argv[2], "PHYLIP_INTERLEAVED") && strcmp (argv[2], "PHYLIP_SEQUENTIAL"))
   {
     usage (argv[0]);
     return (EXIT_FAILURE);
   }

  format = strcmp (argv[2], "PHYLIP_INTERLEAVED") ? PHYLIP_SEQUENTIAL : PHYLIP_INTERLEAVED;

  pd = pl_phylip_parse (argv[1], format);
  if (!pd) 
   {
     printf ("Error while parsing\n");
     return (EXIT_FAILURE);
   }

  dump_struct (pd);
  
//  weight  = pl_phylip_deldups (&pd);
//  printf ("=> Eliminating dups\n");
//  dump_struct (pd);
//  for (i = 0; i < pd->seqlen; ++ i)
//  printf ("%d ", weight[i]);
//  printf ("\n");
  
  free_phylip_struct (pd);
//  free (weight);

  /*
//  dump_struct (pd);

//  printf ("=> Sorting\n");
  // ms = construct_msa_sites (pd, SITES_CREATE | SITES_ELIMINATE_DUPLICATES | SITES_COMPUTE_WEIGHTS);
  ms = construct_msa_sites (pd, SITES_CREATE | SITES_COMPUTE_WEIGHTS);
  dump_sites (ms);


  for (i = 0; i < ms->seqlen; ++ i)
  printf ("%d ", ms->weight[i]);
  printf ("\n");

  free_phylip_struct (pd);
  pd = transpose (ms);
  free_sites_struct (ms);

  dump_struct (pd);
  for (i = 0; i < pd->seqlen; ++ i)
  printf ("%d ", pd->weight[i]);
  printf ("\n");
  free_phylip_struct (pd);
*/


  return (EXIT_SUCCESS);
}
