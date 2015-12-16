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
 * @file alignment.c
 *
 * @brief Collection of routines for reading alignments
 *
 * Auxiliary functions for storing alignments read from predefined file formats
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pll.h"
#include "pllInternal.h"

/** @defgroup alignmentGroup Reading and parsing multiple sequence alignments
    
    This set of functions handles the reading and parsing of several file formats that describe multiple sequence alignments. They are also responsible for storing the alignment in an internal structure
*/
static pllAlignmentData * pllParsePHYLIP (const char * filename);
static pllAlignmentData * pllParseFASTA (const char * filename);
static int read_phylip_header (int * inp, int * sequenceCount, int * sequenceLength);
static __inline int parsedOk (int * actLen, int sequenceCount, int sequenceLength);
static int parse_phylip (pllAlignmentData * alignmentData, int input);
static int getFastaAlignmentInfo (int * inp, int * seqCount, int * seqLen);
static int parseFastaAlignment (pllAlignmentData * alignmentData, int input);

#ifdef __PLL_DEBUG_PARSER
static int
printTokens (int input)
{
  pllLexToken token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.tokenType)
      {
        case PLL_TOKEN_NUMBER:
          printf ("PLL_TOKEN_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_STRING:
          printf ("PLL_TOKEN_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case PLL_TOKEN_EOF:
          printf ("PLL_TOKEN_EOF\n");
          break;
        case PLL_TOKEN_WHITESPACE:
          printf ("PLL_TOKEN_WHITESPACE\n");
          break;
        case PLL_TOKEN_NEWLINE:
          printf ("PLL_TOKEN_NEWLINE\n");
          break;
        case PLL_TOKEN_UNKNOWN:
          printf ("PLL_TOKEN_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        default:
          break;
      }
     /* end of parser */


   }
  while (token.tokenType != PLL_TOKEN_EOF && token.tokenType != PLL_TOKEN_UNKNOWN);

  if (token.tokenType == PLL_TOKEN_UNKNOWN) return (0);

  return (1);
}
#endif

/** @ingroup alignmentGroup
    @brief Initialize alignment structure fields

    Allocates memory for the data structure that will hold the alignment and
    initializes it. It requires the number of sequences \a sequenceCount and
    the length of sequences \a sequenceLength. It returns a pointer to the
    initialized data structure.

    @param sequenceCount
      Number of sequences in the alignment
    
    @param sequenceLength
      Length of the sequences

    @param 
      Initialized alignment data structured
*/
pllAlignmentData *
pllInitAlignmentData (int sequenceCount, int sequenceLength)
 {
   int i;
   pllAlignmentData * alignmentData;
   //void * mem;
   //TUNG
   unsigned char *mem;

   
   /** TODO */
   alignmentData               =  (pllAlignmentData *) rax_malloc (sizeof (pllAlignmentData));
   alignmentData->sequenceData = (unsigned char **) rax_malloc ((sequenceCount + 1) * sizeof (unsigned char *));
   //mem = (void *) rax_malloc (sizeof (unsigned char) * (sequenceLength + 1) * sequenceCount);
   //TUNG
   mem = (unsigned char *)rax_malloc(sizeof(unsigned char) * (sequenceLength + 1) * sequenceCount);
   for (i = 1; i <= sequenceCount; ++i)
    {
      alignmentData->sequenceData[i]                 = (unsigned char *) (&mem[sizeof (unsigned char) * (i - 1) * (sequenceLength + 1)]);
      alignmentData->sequenceData[i][sequenceLength] = 0;
    }
   alignmentData->sequenceData[0] = NULL;
    
   alignmentData->sequenceLabels = (char **) rax_calloc ((sequenceCount + 1), sizeof (char *));

   alignmentData->sequenceCount  = sequenceCount;
   alignmentData->sequenceLength = sequenceLength;
   alignmentData->originalSeqLength = sequenceLength;

   /** TODO: remove siteWeights from alignment */
   alignmentData->siteWeights    = NULL;

   return (alignmentData);
 }

/** @ingroup alignmentGroup
    @brief Deallocates the memory associated with the alignment data structure
    
    Deallocates the memory associated with the alignment data structure \a alignmentData.

    @param alignmentData
      The alignment data structure
*/
void
pllAlignmentDataDestroy (pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++ i)
   {
     rax_free (alignmentData->sequenceLabels[i]);
   }
  rax_free (alignmentData->sequenceLabels);
  rax_free (alignmentData->sequenceData[1]);
  rax_free (alignmentData->sequenceData);
  rax_free (alignmentData->siteWeights);
  rax_free (alignmentData);
}


/** @ingroup alignmentGroup
    @brief Prints the alignment to the console

    @param alignmentData
      The alignment data structure
*/
void 
pllAlignmentDataDumpConsole (pllAlignmentData * alignmentData)
 {
   int i;

   printf ("%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   for (i = 1; i <= alignmentData->sequenceCount; ++ i)
    {
      printf ("%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
    }
 }



static void dump_fasta_content(FILE * fp, pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++i)
     fprintf (fp, ">%s\n%s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
}

static void dump_phylip_content(FILE * fp, pllAlignmentData * alignmentData)
{
  int i;

  for (i = 1; i <= alignmentData->sequenceCount; ++i)
     fprintf (fp, "%s %s\n", alignmentData->sequenceLabels[i], alignmentData->sequenceData[i]);
}

/** @ingroup alignmentGroup
    @brief Dump the alignment to a file of format \a fileFormat

    Dumps the alignment contained in \a alignmentData to file \a filename of type \a fileFormat.

    @note If \a filename exists, all contents will be erased

    @param alignmentData
      Alignment data structure

    @param fileFormat
      Format of output file. Can take the value \b PLL_FORMAT_PHYLIP or \b PLL_FORMAT_FASTA

    @param filename
      Output filename

    @return
      Returns \b PLL_TRUE on success, otherwise \b PLL_FALSE.
*/
int
pllAlignmentDataDumpFile (pllAlignmentData * alignmentData, int fileFormat, const char * filename)
{
  FILE * fp;
  void (*outfun)(FILE *, pllAlignmentData *);
  
  if (fileFormat != PLL_FORMAT_PHYLIP && fileFormat != PLL_FORMAT_FASTA) return (PLL_FALSE);

  outfun = (fileFormat == PLL_FORMAT_PHYLIP) ? dump_phylip_content : dump_fasta_content;

  fp = fopen (filename,"wb");
  if (!fp) return (PLL_FALSE);
  
  /* if PHYLIP print the silly header at the beginning */
  if (fileFormat == PLL_FORMAT_PHYLIP)
   {
     fprintf (fp, "%d %d\n", alignmentData->sequenceCount, alignmentData->sequenceLength);
   }
  
  outfun(fp, alignmentData);

  fclose (fp);
  return (PLL_TRUE);
}



/* ROUTINES FOR PHYLIP PARSING */
/** @ingroup alignmentGroup
    @brief Parse the PHYLIP file header
*/
static int
read_phylip_header (int * inp, int * sequenceCount, int * sequenceLength)
{
  pllLexToken token;
  int input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceCount = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
  if (token.tokenType != PLL_TOKEN_NUMBER) return (0);

  *sequenceLength = atoi (token.lexeme);

  *inp = input;

  return (*sequenceCount && *sequenceLength);
}

static __inline int
parsedOk (int * actLen, int sequenceCount, int sequenceLength)
{
  int i;

  for (i = 1; i <= sequenceCount; ++ i)
   {
     if (actLen[i] != sequenceLength) return (0);
   }
  
  return (1);
}


/** @ingroup alignmentGroup
    @brief Parse the PHYLIP file body
*/
static int
parse_phylip (pllAlignmentData * alignmentData, int input)
{
  int i,j;
  pllLexToken token;
  int * sequenceLength;
  int rc;

  sequenceLength = (int *) rax_calloc (alignmentData->sequenceCount + 1, sizeof (int));

  NEXT_TOKEN
  for (i = 0; ; ++i)
  {
    j = i % alignmentData->sequenceCount;
    if (i < alignmentData->sequenceCount) 
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
          rax_free (sequenceLength);
          return (0);
        }

       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)


       if (token.tokenType != PLL_TOKEN_STRING && token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_FLOAT)
        {
          rax_free (sequenceLength);
          return (0);
        }
       alignmentData->sequenceLabels[i + 1] = my_strndup (token.lexeme, token.len);
       NEXT_TOKEN
       CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     }
    
    while (1)
     {
       if (token.tokenType == PLL_TOKEN_EOF)
        {
          rc = parsedOk (sequenceLength, alignmentData->sequenceCount, alignmentData->sequenceLength);
          rax_free (sequenceLength);
          return (rc);
        }

       if (token.tokenType == PLL_TOKEN_UNKNOWN)
        {
         rax_free (sequenceLength);
         return (0);
        }
       
       if (token.tokenType == PLL_TOKEN_NEWLINE) break;

       if (token.tokenType != PLL_TOKEN_STRING)
        {
          rax_free (sequenceLength);
          return (0);
        }

       if (sequenceLength[j + 1] + token.len > alignmentData->sequenceLength) 
        {
          fprintf (stderr, "Sequence %d is larger than specified\n", j + 1);
          rax_free (sequenceLength);
          return (0);
        }
       memmove (alignmentData->sequenceData[j + 1] + sequenceLength[j + 1], token.lexeme, token.len);
       sequenceLength[j + 1] += token.len;

       NEXT_TOKEN
       CONSUME (PLL_TOKEN_WHITESPACE)
     }
    CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE);
  }
}

/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER PLL_TOKEN_WHITESPACE PLL_TOKEN_NUMBER ENDL
 *     ENDL -> PLL_TOKEN_WHITESPACE PLL_TOKEN_NEWLINE | PLL_TOKEN_NEWLINE
 *     DATA -> PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING ENDL DATA | 
 *             PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF |
 *             PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_WHITESPACE PLL_TOKEN_STRING PLL_TOKEN_EOF
 */

/** @ingroup alignmentGroup
    @brief Parse a PHYLIP file

    Parses the PHYLIP file \a filename and returns a ::pllAlignmentData structure
    with the alignment.

    @param filename
      Name of file to be parsed

    @return
      Returns a structure of type ::pllAlignmentData that contains the alignment, or \b NULL
      in case of failure.
*/
static pllAlignmentData *
pllParsePHYLIP (const char * filename)
{
  int 
    i, input, sequenceCount, sequenceLength;
  char * rawdata;
  long filesize;
  pllAlignmentData * alignmentData;

  rawdata = pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     errno = PLL_ERROR_FILE_OPEN;
     return (NULL);
   }
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (&input, &sequenceCount, &sequenceLength))
   {
     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     errno = PLL_ERROR_PHYLIP_HEADER_SYNTAX;
     return (NULL);
   }

  lex_table_amend_phylip();

  /* allocate alignment structure */
  alignmentData = pllInitAlignmentData (sequenceCount, sequenceLength);

  if (! parse_phylip (alignmentData, input))
   {
     errno = PLL_ERROR_PHYLIP_BODY_SYNTAX;
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
     rax_free (rawdata);
     return (NULL);
   }
  
  lex_table_restore();
  rax_free (rawdata);

  alignmentData->siteWeights  = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i) 
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}

pllAlignmentData *
pllParsePHYLIPString (const char *rawdata, long filesize)
{
  int
    i, input, sequenceCount, sequenceLength;
//  char * rawdata;
//  long filesize;
  pllAlignmentData * alignmentData;

//  rawdata = pllReadFile (filename, &filesize);
//  if (!rawdata)
//   {
//     errno = PLL_ERROR_FILE_OPEN;
//     return (NULL);
//   }

  init_lexan (rawdata, filesize);
  input = get_next_symbol();

  /* parse the header to obtain the number of taxa and sequence length */
  if (!read_phylip_header (&input, &sequenceCount, &sequenceLength))
   {
//     rax_free (rawdata);
     fprintf (stderr, "Error while parsing PHYLIP header (number of taxa and sequence length)\n");
     errno = PLL_ERROR_PHYLIP_HEADER_SYNTAX;
     return (NULL);
   }

  lex_table_amend_phylip();

  /* allocate alignment structure */
  alignmentData = pllInitAlignmentData (sequenceCount, sequenceLength);

  if (! parse_phylip (alignmentData, input))
   {
     errno = PLL_ERROR_PHYLIP_BODY_SYNTAX;
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
//     rax_free (rawdata);
     return (NULL);
   }

  lex_table_restore();
//  rax_free (rawdata);

  alignmentData->siteWeights  = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i)
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}

/* FASTA routines */
/* only check whether it is a valid alignment in fasta format */
/** @ingroup alignmentGroup
    @brief Get information about the FASTA alignment

    Get the information such as number of sequences and length of sequences of a FASTA alignment

    @return
      Returns \b PLL_TRUE if the alignment is valid, otherwise \b PLL_FALSE
*/
static int
getFastaAlignmentInfo (int * inp, int * seqCount, int * seqLen)
{
  pllLexToken token;
  int input;

  input = *inp;

  *seqCount = *seqLen = 0;

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_STRING) return (PLL_FALSE);

  while (1)
   {
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (PLL_TRUE);

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          if (token.len < 2 || token.lexeme[0] != '>') return (0);
          break;
        default:
          return (PLL_FALSE);
      }
     
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

     /* read second token (sequence) */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (PLL_FALSE);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          if (!*seqLen)
            *seqLen = token.len;
          else
           {
             if (*seqLen != token.len) return (0);
           }
          break;
        default:
          return (PLL_FALSE);
      }
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     ++ (*seqCount);
   }

  return (PLL_TRUE);
}

/** @ingroup alignmentGroup
    @brief Check whether the FASTA content is valid
*/
static int
parseFastaAlignment (pllAlignmentData * alignmentData, int input)
{
  pllLexToken token;
  int i;

  NEXT_TOKEN
  CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

  if (token.tokenType != PLL_TOKEN_NUMBER && token.tokenType != PLL_TOKEN_STRING) return (0);

  i = 1;
  while (1)
   {
     /* first parse the sequence label */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (1);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          alignmentData->sequenceLabels[i] = my_strndup (token.lexeme + 1, token.len - 1);
          break;
        default:
          return (0);
      }
     
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)

     /* now parse the sequence itself */
     switch (token.tokenType)
      {
        case PLL_TOKEN_EOF:
          return (0);
          break;

        case PLL_TOKEN_NUMBER:
        case PLL_TOKEN_STRING:
          memmove (alignmentData->sequenceData[i], token.lexeme, token.len);
          break;
        default:
          return (0);
      }
     NEXT_TOKEN
     CONSUME(PLL_TOKEN_WHITESPACE | PLL_TOKEN_NEWLINE)
     ++ i;
   }
}


/** @ingroup alignmentGroup
    @brief Parse a FASTA file
    
    Parses the FASTA file \a filename and returns a ::pllAlignmentData structure
    with the alignment.

    @param filename
      Name of file to be parsed

    @return
      Returns a structure of type ::pllAlignmentData that contains the alignment, or \b NULL
      in case of failure.
*/
static pllAlignmentData *
pllParseFASTA (const char * filename)
{
  int
    i,
    seqLen,
    seqCount,
    input;
  long filesize;

  char * rawdata;
  pllAlignmentData * alignmentData;

  rawdata = pllReadFile (filename, &filesize);
  if (!rawdata)
   {
     errno = PLL_ERROR_FILE_OPEN;
     return (NULL);
   }

  lex_table_amend_fasta ();
  
  init_lexan (rawdata, filesize);
  input = get_next_symbol ();


  if (!getFastaAlignmentInfo (&input, &seqCount, &seqLen))
   {
     errno = PLL_ERROR_FASTA_SYNTAX;
     lex_table_restore ();
     rax_free (rawdata);
     return (NULL);
   }
  
  alignmentData = pllInitAlignmentData (seqCount, seqLen);
  
  printf ("\n---------------\n\n");

  init_lexan (rawdata, filesize);
  input = get_next_symbol ();

  if (!parseFastaAlignment (alignmentData, input))
   {
     errno = PLL_ERROR_FASTA_SYNTAX;
     pllAlignmentDataDestroy (alignmentData);
     lex_table_restore();
     rax_free(rawdata);
     return (NULL);
   }

  /* allocate alignment structure */


  lex_table_restore ();
  rax_free (rawdata);

  alignmentData->siteWeights = (int *) rax_malloc (alignmentData->sequenceLength * sizeof (int));
  for (i = 0; i < alignmentData->sequenceLength; ++ i)
    alignmentData->siteWeights[i] = 1;

  return (alignmentData);
}



/** @ingroup alignmentGroup
    @brief Parse a file that contains a multiple sequence alignment

    Parses the file \a filename of type \a fileType which contains a multiple sequence alignment.
    The supported file types are the sequential and interleaved versions of PHYLIP format, and
    the FASTA format. The parsed alignment is returned as a pointer to a structure of type
    ::pllAlignmentData

    @param fileType
      Type of file to parse. Can be either \b PLL_FORMAT_PHYLIP or \b PLL_FORMAT_FASTA

    @param filename
      Name of file to parse

    @return
      Returns a structure of type ::pllAlignmentData that contains the multiple sequence alignment,
      otherwise returns \b NULL in case of failure.
*/
pllAlignmentData *
pllParseAlignmentFile (int fileType, const char * filename)
{

  switch (fileType)
   {
     case PLL_FORMAT_PHYLIP:
       return (pllParsePHYLIP (filename));
     case PLL_FORMAT_FASTA:
       return (pllParseFASTA (filename));
     default:
       /* RTFM */
       errno = PLL_ERROR_INVALID_FILETYPE;
       return (NULL);
   }
}
