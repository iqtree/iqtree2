#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylip.h"
#include "msa_sites.h"
#include "ssort.h"


static void
sort_sites (struct msa_sites * ms)
{
  int * oi;
  oi = ssort1main (ms->site, ms->seqlen);
  free (oi);
}

/* This function deallocates the memory of duplicate sites and redirects them
 * to the unique occurrence of the string representing the specific site. This
 * is a memory saving technique
 */
static void
compact_sites (struct msa_sites * ms)
{
  int unique = 1;
  int i, j;
  int * oi;
  char ** tmp;

  oi = ssort1main (ms->site, ms->seqlen);

  /* Locate unique sites and mark pointers, free duplicates */
  for (i = 1; i < ms->seqlen; ++ i)
   {
     if (!strcmp (ms->site[i], ms->site[i - 1])) 
      {
        free (ms->site[i]);
        ms->site[i] = ms->site[i - 1];
      }
     else
      {
        ++unique;
      }
   }

  /* store unique sites in a null-terminated array */
  ms->unique_site    = (char **) malloc ((unique + 1) * sizeof (char *));
  ms->unique_site[0] = ms->site[0];
  for (j = 1, i = 1; i < ms->seqlen; ++ i)
   {
     if (ms->site[i] != ms->site[i - 1])
      {
        ms->unique_site[j++] = ms->site[i];
      }
   }
  ms->unique_site[j] = NULL;

  /* TODO: Consider sorting the list of pointers when deallocating the
     structure instead of storing the unique_site array.
  */

  /*  sorted array */
  tmp = (char **) malloc (ms->seqlen * sizeof (char *));
  for (i = 0; i < ms->seqlen; ++ i)
   {
     tmp[oi[i]] = ms->site[i];
   }
  free (oi);
  free (ms->site);
  ms->site = tmp;
}

static struct phylip_data *
transpose (struct msa_sites * ms)
{
  struct phylip_data * pd;
  int i, j;

//  pd = alloc_phylip_struct (ms->taxa, ms->seqlen);
  pd = alloc_phylip_struct (ms->taxa, ms->seqlen);

  for (i = 0; i < pd->taxa; ++ i)
   {
     pd->label[i] = strdup (ms->label[i]);
   }
  for (i = 0; i < pd->taxa; ++ i)
   {
     pd->seq[i] = (char *) calloc ((ms->seqlen + 1), sizeof (char));
   }

  for (i = 0; i < pd->seqlen; ++ i)
   {
     for (j = 0; j < pd->taxa; ++ j)
      {
        pd->seq[j][i] = ms->site[i][j];
      }
   }
  return (pd);
}


static int *
eliminate_dups (struct msa_sites * ms)
{
  int unique = 1;
  int i, j;
  int * oi;

  oi = ssort1main (ms->site, ms->seqlen);

  /* Locate unique sites and mark pointers, free duplicates */
  for (i = 1; i < ms->seqlen; ++ i)
   {
     if (!strcmp (ms->site[i], ms->site[i - 1])) 
      {
        free (ms->site[i]);
        ms->site[i] = ms->site[i - 1];
      }
     else
      {
        ++unique;
      }
   }

  /* store unique sites in a null-terminated array */
  ms->unique_site    = (char **) malloc ((unique) * sizeof (char *));
  ms->unique_site[0] = ms->site[0];
  ms->weight         = (int *) malloc ((unique) * sizeof (char *));
  ms->weight[0]      = 1;

  for (j = 0, i = 1; i < ms->seqlen; ++ i)
   {
     if (ms->site[i] != ms->site[i - 1])
      {
        ms->unique_site[++j] = ms->site[i];
        ms->weight[j] = 1;
      }
     else
      {
        ++ ms->weight[j];
      }
   }

  free(ms->site);
  ms->site = ms->unique_site;

  ms->unique_site = NULL;
  ms->seqlen = unique;
  free (oi);

  return (ms->weight);
}


int *
pl_phylip_deldups (struct phylip_data ** pd)
{
  struct msa_sites * ms;
  int * weights;

  ms = construct_msa_sites (*pd, SITES_CREATE);

  free_phylip_struct (*pd);

  weights = eliminate_dups (ms);

  *pd = transpose (ms);

  free_sites_struct (ms);

  return (weights);
}

struct msa_sites *
alloc_sites_struct (int taxa, int seqlen)
 {
   int i;
   struct msa_sites * ms;

   ms = (struct msa_sites *) malloc (sizeof (struct msa_sites));

   ms->taxa        = taxa;
   ms->seqlen      = seqlen;
   ms->label       = (char **) calloc (taxa, sizeof (char *));
   ms->site        = (char **) calloc (seqlen, sizeof (char *));
   ms->weight      = NULL;
   ms->unique_site = NULL;

   for (i = 0; i < seqlen; ++ i)
    {
      ms->site[i] = (char *) calloc ((taxa + 1), sizeof (char));
    }

   return (ms);
 }



void
free_sites_struct (struct msa_sites * ms)
{
  int i;

  for (i = 0; i < ms->taxa; ++ i)
   {
     free (ms->label[i]);
   }

  if (ms->unique_site)
   {
     i = 0;
     while (ms->unique_site[i])
      {
        free (ms->unique_site[i]);
        ++i;
      }
     free (ms->unique_site);
   }
  else
   {
     for (i = 0; i < ms->seqlen; ++ i)
      {
        free (ms->site[i]);
      }
   }

  free (ms->label);
  free (ms->site);
  free (ms);
}



/* DEBUG function for dumping the alignment */
void 
dump_sites (struct msa_sites * ms)
 {
   int i, j;

   printf ("%d %d\n", ms->taxa, ms->seqlen);
   for (i = 0; i < ms->taxa; ++ i)
    {
      printf ("|%s| |", ms->label[i]);
      for (j = 0; j < ms->seqlen; ++ j)
       {
         printf ("%c", ms->site[j][i]);
       }

      printf ("|\n");
    }
 }

struct msa_sites * 
construct_msa_sites (struct phylip_data * pd, int flags)
{
  struct msa_sites * ms;
  int i, j;
  
  ms = alloc_sites_struct (pd->taxa, pd->seqlen);

  for (i = 0; i < pd->taxa; ++ i)
    ms->label[i] = strdup (pd->label[i]);
  
  for (i = 0; i < pd->seqlen; ++ i)
   {
     for (j = 0; j < pd->taxa; ++ j)
      {
        ms->site[i][j] = pd->seq[j][i];
      }
   }

/*
  if (flags & SITES_ELIMINATE_DUPLICATES)
   {
     eliminate_dups (ms, flags & SITES_COMPUTE_WEIGHTS);
   }
  else if (flags & SITES_COMPACT)
   {
     compact_sites (ms);
   }
  else
   {
     if (flags & SITES_SORTED)
      {
        sort_sites (ms);
      }
     if (flags & SITES_COMPUTE_WEIGHTS)
      {
        // set all weights to 1
        ms->weight = (int *) malloc  (ms->seqlen * sizeof (int));
        for (i = 0; i < ms->seqlen; ++ i)
          ms->weight[i] = 1;
      }
   }
*/
  return (ms);
}

