#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "pll.h"

/** @file part.c
    
    @brief Collection of routines for parsing and processing a partition (model) file

    @defgroup parsePartitionFileGroup Reading and parsing partition (model) files

    This set of functions handles the reading and parsing of partition files, i.e.
    files that contain alignment partition definitions and corresponding models.
*/


extern const char *protModels[PLL_NUM_PROT_MODELS];

static struct pllHashTable * hashTable;

static void destroy_model_names(void)
{
  pllHashDestroy (&hashTable, PLL_TRUE);
}

static void init_model_names (void)
{
  int i;
  int * item;

  hashTable = pllHashInit (PLL_NUM_PROT_MODELS);

  for (i = 0; i < PLL_NUM_PROT_MODELS; ++ i)
   {
     item  = (int *) rax_malloc (sizeof (int));
     *item = i;
     pllHashAdd (hashTable, protModels[i], (void *) item);
   }
}

/** @ingroup parsePartitionFileGroup
    @brief Destroy queue structure that contains parsed information from a partition file

    Destroys the structure, and therefore frees allocated memory, that holds parsed information
    from a partition (model) file

    @param partitions
      Queue structure with parsed info
*/
void
pllQueuePartitionsDestroy (pllQueue ** partitions)
{
  pllPartitionInfo * pi;
  pllPartitionRegion * region;

  while (pllQueueRemove (*partitions, (void **)&pi))
   {
     while (pllQueueRemove (pi->regionList, (void **) &region))
      {
        rax_free (region);
      }
     rax_free (pi->regionList);
     rax_free (pi->partitionName);
     rax_free (pi->partitionModel);
     rax_free (pi);
   }
  rax_free (*partitions);
}

static pllQueue *
parse_partition (int * inp)
{
  int input, i;
  pllLexToken token;
  int lines = 0;
  pllQueue * partitions;
  pllPartitionInfo * pi;
  pllPartitionRegion * region;
  int * item;

  input  = *inp;

  NEXT_TOKEN

  pllQueueInit (&partitions);
  while (token.tokenType != PLL_TOKEN_EOF)
  {
    ++ lines;
    pi = (pllPartitionInfo *) rax_calloc (1, sizeof (pllPartitionInfo));
    pllQueueInit (&(pi->regionList));
    pllQueueAppend (partitions, (void *)pi);
    CONSUME (PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)


    /* read partition type */
    if (token.tokenType != PLL_TOKEN_STRING) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    //pi->partitionModel = strndup (token.lexeme, token.len);
    pi->partitionModel = (char *) rax_malloc ((token.len + 1) * sizeof (char));
    strncpy (pi->partitionModel, token.lexeme, token.len);
    pi->partitionModel[token.len] = 0;
    for (i = 0; i < token.len; ++i) pi->partitionModel[i] = toupper(pi->partitionModel[i]);
    // check partition model
    pi->protModels = -1;

    if (!strcmp (pi->partitionModel, "DNA"))
     {
       pi->protModels = PLL_FALSE;
       pi->protFreqs  = PLL_FALSE;
       pi->dataType   = PLL_DNA_DATA;
       pi->optimizeBaseFrequencies = PLL_FALSE; 
     }
    else if (!strcmp (pi->partitionModel, "DNAX"))
     {
       pi->protModels = PLL_FALSE;
       pi->protFreqs  = PLL_FALSE;
       pi->dataType   = PLL_DNA_DATA;
       pi->optimizeBaseFrequencies = PLL_TRUE; 
     }
    else
     {                  /* check for protein data */
       pi->dataType  = PLL_AA_DATA;

       if (pllHashSearch (hashTable, pi->partitionModel, (void **) &item))
        {
          pi->protModels = *item;
          pi->protFreqs  = PLL_FALSE;
          pi->optimizeBaseFrequencies = PLL_FALSE;
        }
       else
        {
          if (pi->partitionModel[token.len - 1] == 'X')
           {
             pi->partitionModel[token.len - 1] = '\0';
             if (pllHashSearch (hashTable, pi->partitionModel, (void **) &item))
              {
                pi->protModels = *item;
                pi->protFreqs  = PLL_FALSE;
                pi->optimizeBaseFrequencies = PLL_TRUE;
              }
             pi->partitionModel[token.len - 1] = 'X';
           }
          else if (pi->partitionModel[token.len - 1] == 'F')
           {
             pi->partitionModel[token.len - 1] = '\0';
             if (pllHashSearch (hashTable, pi->partitionModel, (void **) &item))
              {
                pi->protModels = *item;
                pi->protFreqs  = PLL_TRUE;
                pi->optimizeBaseFrequencies = PLL_FALSE;
              }
             pi->partitionModel[token.len - 1] = 'F';
           }
          else
           {
             pllQueuePartitionsDestroy (&partitions);
             return (0);
           }
        }
     }

    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE)

    if (token.tokenType != PLL_TOKEN_COMMA) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE)

    /* read partition name */
    if (token.tokenType != PLL_TOKEN_STRING) 
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    //pi->partitionName = strndup (token.lexeme, token.len);
    pi->partitionName = (char *) rax_malloc ((token.len + 1) * sizeof (char));
    strncpy (pi->partitionName, token.lexeme, token.len);
    pi->partitionName[token.len] = 0;

    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE)

    /* read equal sign */
    if (token.tokenType != PLL_TOKEN_EQUAL)
     {
       pllQueuePartitionsDestroy (&partitions);
       return (0);
     }
    NEXT_TOKEN
    CONSUME(PLL_TOKEN_WHITESPACE)

    /* read rhs */
    while (1)
    {
      region = (pllPartitionRegion *) rax_malloc (sizeof (pllPartitionRegion));
      if (token.tokenType != PLL_TOKEN_NUMBER) 
       {
         pllQueuePartitionsDestroy (&partitions);
         return (0);
       }
      region->start  = region->end = atoi (token.lexeme);  
      region->stride = 1;
      NEXT_TOKEN
      CONSUME(PLL_TOKEN_WHITESPACE)
      
      if  (token.tokenType == PLL_TOKEN_DASH)
       {
         NEXT_TOKEN
         CONSUME(PLL_TOKEN_WHITESPACE)
         if (token.tokenType != PLL_TOKEN_NUMBER) 
          {
            pllQueuePartitionsDestroy (&partitions);
            return (0);
          }
         region->end = atoi (token.lexeme);
         if (region->end < region->start)
          {
            pllQueuePartitionsDestroy (&partitions);
            return (0);
          }
         NEXT_TOKEN
         CONSUME(PLL_TOKEN_WHITESPACE)
         if (token.tokenType == PLL_TOKEN_SLASH)
          {
            NEXT_TOKEN
            CONSUME(PLL_TOKEN_WHITESPACE)
            if (token.tokenType != PLL_TOKEN_NUMBER) 
             {
               pllQueuePartitionsDestroy (&partitions);
               return (0);
             }
            region->stride = atoi (token.lexeme);
            NEXT_TOKEN
          }
         CONSUME(PLL_TOKEN_WHITESPACE)
       }
       pllQueueAppend (pi->regionList, (void *)region);
      
      if (token.tokenType != PLL_TOKEN_COMMA) break;
      NEXT_TOKEN
      CONSUME(PLL_TOKEN_WHITESPACE)
    }
   CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
  }
 
 return (partitions);
} 

/** @ingroup parsePartitionFileGroup
    @brief Dump a parsed partition file in the console

    Prints the parsed contents of a partition file to the console

    @param partitions
      A queue structure that contains the parsed information
*/
void 
pllPartitionDump (pllQueue * partitions)
{
   struct pllQueueItem * elm;
   struct pllQueueItem * regionList;
   pllPartitionInfo * pi;
   pllPartitionRegion * region;

   elm = partitions->head;

   while (elm)
    {
      pi  = (pllPartitionInfo *) elm->item;
      printf ("%s, %s = ", pi->partitionModel, pi->partitionName);
      regionList = pi->regionList->head;
      while (regionList)
       {
         region = (pllPartitionRegion *) regionList->item;
         printf ("%d", region->start);
         if (region->start != region->end)
          {
            printf ("-%d", region->end);
            if (region->stride != 1) printf ("/%d", region->stride);
          }
         regionList = regionList->next;
         if (regionList) printf (", ");
       }
      printf ("\n");

      elm = elm->next;
    }
}

/** @ingroup parsePartitionFileGroup
    @brief Parse a partition (model) file

    Parses the partition file \a filename and stores the information
    a queue structure ::pllQueue

    @param filename
      Name of the partition file
    
    @return
      Queue structure with parsed information
*/
pllQueue *
pllPartitionParse (const char * filename)
{
  long n;
  char * rawdata;
  int input;
  pllQueue * partitions;

  rawdata = pllReadFile (filename, &n);
  if (!rawdata)
   {
     fprintf (stderr, "Error while opening/reading file %s\n", filename);
     return (0);
   }

  printf ("%s\n\n", rawdata);

  n = strlen (rawdata);

  init_lexan (rawdata, n);
  input = get_next_symbol();

  init_model_names();
  partitions = parse_partition (&input);
  destroy_model_names();
  
  rax_free (rawdata);
  return (partitions);
}
