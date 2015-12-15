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
 * @file parsePartition.c
 * @brief Collection of routines for parsing and processing a partition (model) file
 *
 * @defgroup parsePartitionFileGroup Reading and parsing partition (model) files
 * This set of functions handles the reading and parsing of partition files, i.e.
 * files that contain alignment partition definitions and corresponding models.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "pll.h"
#include "pllInternal.h"

extern const char *protModels[PLL_NUM_PROT_MODELS];

static void destroy_model_names(pllHashTable * hashTable)
{
  pllHashDestroy (&hashTable, rax_free);
}

static pllHashTable * init_model_names (void)
{
  int i;
  int * item;

  pllHashTable * hashTable;
  hashTable = pllHashInit (PLL_NUM_PROT_MODELS);

  for (i = 0; i < PLL_NUM_PROT_MODELS; ++ i)
   {
     item  = (int *) rax_malloc (sizeof (int));
     *item = i;
     pllHashAdd (hashTable, pllHashString(protModels[i], hashTable->size), protModels[i], (void *) item);
   }
  return hashTable;
}

/** @ingroup parsePartitionFileGroup
    @brief Destroy queue structure that contains parsed information from a partition file

    Destroys the structure, and therefore frees allocated memory, that holds parsed information
    from a partition (model) file

    @param partitions
      Queue structure with parsed info
*/
void pllQueuePartitionsDestroy (pllQueue ** partitions)
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

static pllQueue * parse_partition (int * inp, pllHashTable * proteinModelsHash)
{
  int input, i;
  pllLexToken token;
  int lines = 0;
  pllQueue * partitions;
  pllPartitionInfo * pi;
  pllPartitionRegion * region;
  int * protIndexPtr;
  char * modelptr;

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
    pi->partitionModel = my_strndup (token.lexeme, token.len);
    for (i = 0; i < token.len; ++i) pi->partitionModel[i] = toupper(pi->partitionModel[i]);

    // check partition model
    pi->protModels              = -1;
    pi->protUseEmpiricalFreqs   = PLL_FALSE;
    pi->ascBias                 = PLL_FALSE;
    pi->optimizeBaseFrequencies = PLL_FALSE;

    /* check if the model contains Asc bias */
    if (!strncmp(pi->partitionModel, "ASC_", 4))
      {
        pi->ascBias = PLL_TRUE;
        modelptr    = pi->partitionModel + 4;
      }
     else
        modelptr    = pi->partitionModel;

    /* check first for BINARY */
    if (!strcmp(modelptr, "BIN") || !strcmp(modelptr, "BINX"))
     {
       pi->dataType = PLL_BINARY_DATA;

       if (!strcmp(modelptr, "BINX"))
         pi->optimizeBaseFrequencies = PLL_TRUE;
     }  /* now for DNA */
    else if (!strcmp(modelptr, "DNA") || !strcmp(modelptr, "DNAX"))
     {
       pi->dataType   = PLL_DNA_DATA;

       if (!strcmp(modelptr, "DNAX")) 
         pi->optimizeBaseFrequencies = PLL_TRUE; 
     }
    else
     {                  /* and  protein data */
       pi->dataType  = PLL_AA_DATA;

       if (pllHashSearch (proteinModelsHash, modelptr, (void **) &protIndexPtr))
        {
          pi->protModels              = *protIndexPtr;
          pi->protUseEmpiricalFreqs   = PLL_FALSE;
          pi->optimizeBaseFrequencies = PLL_FALSE;
        }
       else
        {
          if (modelptr[token.len - 1] == 'X')
           {
             modelptr[token.len - 1] = '\0';
             if (pllHashSearch (proteinModelsHash, modelptr, (void **) &protIndexPtr))
              {
                pi->protModels              = *protIndexPtr;
                pi->optimizeBaseFrequencies = PLL_TRUE;
              }
             modelptr[token.len - 1] = 'X';
           }
          else if (modelptr[token.len - 1] == 'F')
           {
             modelptr[token.len - 1] = '\0';
             if (pllHashSearch (proteinModelsHash, modelptr, (void **) &protIndexPtr))
              {
                pi->protModels              = *protIndexPtr;
                pi->protUseEmpiricalFreqs   = PLL_TRUE;
              }
             modelptr[token.len - 1] = 'F';
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
    pi->partitionName = my_strndup (token.lexeme, token.len);

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

    @param partitions Queue structure containing parsed information
*/
void pllPartitionDump (pllQueue * partitions)
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

    Parses the partition file \a filename and stores the information in a queue
    structure ::pllQueue

    @param filename Name of the partition file
    @return Queue structure with parsed information
*/
pllQueue * pllPartitionParse (const char * filename)
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

  n = strlen (rawdata);

  init_lexan (rawdata, n);
  input = get_next_symbol();

  pllHashTable * model_names = init_model_names();
  partitions  = parse_partition (&input, model_names);
  destroy_model_names(model_names);
  
  rax_free (rawdata);
  return (partitions);
}

/** @ingroup parsePartitionFileGroup
    @brief Parse a partition (model) file

    Parses the partition information stored in string \a p and stores the
    information in a queue structure ::pllQueue

    @param p Partition information string
    @return  Queue structure with parsed information
*/
pllQueue * pllPartitionParseString (const char * p)
{
  long n;
  int input;
  pllQueue * partitions;

  n = strlen(p);
  init_lexan (p, n);
  input = get_next_symbol();

  pllHashTable * model_names;
  model_names = init_model_names();
  partitions = parse_partition (&input, model_names);
  destroy_model_names(model_names);
  
  return (partitions);
}
