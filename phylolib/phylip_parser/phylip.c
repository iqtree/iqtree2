#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lexer.h"
#include "phylip.h"
#include "xalloc.h"
#include "msa_sites.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define SWAP(x,y) do{ __typeof__ (x) _t = x; x = y; y = _t; } while(0)

#define CONSUME(x)         while (token.class & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

struct rawdata
 {
   unsigned char ** oa;         /* original alignment */
   unsigned char ** pats;       /* unique site patterns */
 };


int debug = 0;
static char * 
read_file (const char * filename, int * n)
{
  FILE       * fp;
  char       * rawdata;

  fp = fopen (filename, "r");
  if (!fp) return (NULL);

  /* obtain file size */
  if (fseek (fp, 0, SEEK_END) == -1) return (NULL);
  *n = ftell (fp);
  if (*n == -1) return (NULL);
  rewind (fp);

  rawdata = (char *) xmalloc ((*n) * sizeof (char));
  if (!rawdata) return (NULL);

  if (fread (rawdata, sizeof (char), *n, fp) != *n) return (NULL);

  fclose (fp);

  return (rawdata);
}

static int
print_tokens (int input)
{
  struct lex_token token;

  do
   {
     NEXT_TOKEN

     /* begin of parser */
     switch (token.class)
      {
        case LEX_NUMBER:
          printf ("LEX_NUMBER (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_STRING:
          printf ("LEX_STRING (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
        case LEX_EOF:
          printf ("LEX_EOF\n");
          break;
        case LEX_WHITESPACE:
          printf ("LEX_WHITESPACE\n");
          break;
        case LEX_NEWLINE:
          printf ("LEX_NEWLINE\n");
          break;
        case LEX_UNKNOWN:
          printf ("LEX_UNKNOWN (%.*s, %d)\n", token.len, token.lexeme, token.len);
          break;
      }
     /* end of parser */


   }
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN);

  if (token.class == LEX_UNKNOWN) return (0);

  return (1);
}

struct phylip_data *
alloc_phylip_struct (int taxa, int seqlen)
 {
   int i;
   struct phylip_data * pd;

   pd = (struct phylip_data *) malloc (sizeof (struct phylip_data));

   pd->taxa   = taxa;
   pd->seqlen = seqlen;
   pd->label  = (char **) calloc (taxa, sizeof (char *));
   pd->seq    = (char **) calloc (taxa, sizeof (char *));

   for (i = 0; i < taxa; ++ i)
    {
      pd->seq[i] = (char *) calloc ((seqlen + 1), sizeof (char));
    }

   return (pd);
 }

void
free_phylip_struct (struct phylip_data * pd)
{
  int i;

  for (i = 0; i < pd->taxa; ++ i)
   {
     free (pd->label[i]);
     free (pd->seq[i]);
   }
  free (pd->label);
  free (pd->seq);
  free (pd);
}


void 
dump_struct (struct phylip_data * pd)
 {
   int i;

   printf ("=> Dumping phylip_data\n");
   printf ("%d %d\n", pd->taxa, pd->seqlen);
   for (i = 0; i < pd->taxa; ++ i)
    {
      printf ("|%s| |%s|\n", pd->label[i], pd->seq[i]);
    }
 }

static struct phylip_data *
read_phylip_header (char * rawdata, int * inp)
{
  struct lex_token token;
  struct phylip_data * pd;
  int taxa, seqlen, input;

  input = *inp;


  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

  if (token.class != LEX_NUMBER) return (NULL);

  taxa = atoi (token.lexeme);

  NEXT_TOKEN
  CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
  if (token.class != LEX_NUMBER) return (NULL);

  seqlen = atoi (token.lexeme);
//  NEXT_TOKEN
//  CONSUME(LEX_WHITESPACE)
//  CONSUME(LEX_NEWLINE)

  /* allocate memory for the alignment */

  if (!taxa || !seqlen)
   {
     return (NULL);
   }

  pd   = alloc_phylip_struct (taxa, seqlen);
  *inp = input;

  return (pd);
}


/* Phylip parsers. Use the following attributed grammar 
 * 
 *        S -> HEADER ENDL DATA
 *   HEADER -> LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL |
 *             LEX_WHITESPACE LEX_NUMBER LEX_WHITESPACE LEX_NUMBER ENDL
 *     ENDL -> LEX_WHITESPACE LEX_NEWLINE | LEX_NEWLINE
 *     DATA -> LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING ENDL DATA | 
 *             LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF |
 *             LEX_WHITESPACE LEX_STRING LEX_WHITESPACE LEX_STRING LEX_EOF
 */
static int
parse_phylip_sequential (char * rawdata, struct phylip_data * pd, int input)
{
  int i;
  struct lex_token token;
  int * seq_size;

  seq_size = (int *) calloc (pd->taxa, sizeof (int));

  NEXT_TOKEN

  /* first parse the TAXA-LABEL SEQUENCE format */
  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN && i < pd->taxa)
   {
//     printf ("%d\n", i);
     CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

     /* do something with the label */
     if (token.class != LEX_STRING && token.class != LEX_NUMBER) 
      {
//        FREE (2, rawdata, seq_size);
        free (seq_size);
        return (0);
      }
     pd->label[i] = strndup(token.lexeme, token.len);

     NEXT_TOKEN
     do 
      {
        CONSUME(LEX_WHITESPACE | LEX_NEWLINE)

        /* do something with the sequence(s) */
        if (token.class != LEX_STRING)
         {
  //         FREE (2, rawdata, seq_size);
           free (seq_size);
           return (0);
         }
        seq_size[i] += token.len;
//        if (i == 81)
 //        {
  //         printf ("%.*s %d\n", token.len, token.lexeme, seq_size[i]);
   //      }
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }

        NEXT_TOKEN
//        if (i == 81 && seq_size[i] == 1050) { debug = 1; printf ("DEBUG: %d %d\n", token.len, token.lexeme[0]);}
      } while (seq_size[i] < pd->seqlen);

     if (seq_size[i] > pd->seqlen) break;
     ++i;
   }
 if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
  {
    free (seq_size);
//    FREE (2, rawdata, seq_size);
    return (0);
  }
 
 //FREE (2, seq_size, rawdata);
  free (seq_size);

  return (1);
}

static int
parse_phylip_interleaved (char * rawdata, struct phylip_data * pd, int input)
{
  int i;
  struct lex_token token;
  int * seq_size;

  seq_size = (int *) calloc (pd->taxa, sizeof (int));

  NEXT_TOKEN

  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN && i < pd->taxa)
   {
     CONSUME(LEX_WHITESPACE | LEX_NEWLINE)
     
     if (token.class != LEX_STRING && token.class != LEX_NUMBER)
      {
//        FREE (2, rawdata, seq_size);
        free (seq_size);
        return (0);
      }
     pd->label[i] = strndup(token.lexeme, token.len);

     NEXT_TOKEN

     while (token.class == LEX_WHITESPACE)
      {
        CONSUME(LEX_WHITESPACE)
        if (token.class != LEX_STRING)
         {
           free (seq_size);
//           FREE (2, rawdata, seq_size);
           return (0);
         }
        seq_size[i] += token.len;
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }
        NEXT_TOKEN
      }
     if (seq_size[i] > pd->seqlen) break;
     ++ i;
   }
  if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
   {
//     FREE (2, rawdata, seq_size);
     free (seq_size);
     return (0);
   }

  i = 0;
  while (token.class != LEX_EOF && token.class != LEX_UNKNOWN)
   {
     if (i % pd->taxa == 0) i = 0;

     CONSUME(LEX_WHITESPACE | LEX_NEWLINE);

     do
      {
        CONSUME(LEX_WHITESPACE)
        if (token.class != LEX_STRING)
         {
 //          FREE (2, rawdata, seq_size);
           free (seq_size);
           return (0);
         }
        seq_size[i] += token.len;
        if (seq_size[i] <= pd->seqlen)
         {
           strncat(pd->seq[i], token.lexeme, token.len);
         }
        else
         {
           break;
         }
        NEXT_TOKEN
      } while (token.class == LEX_WHITESPACE);
     if (seq_size[i] > pd->seqlen) break;
     ++i;
     CONSUME(LEX_NEWLINE);
   }
  if (token.class == LEX_UNKNOWN || (i < pd->taxa && seq_size[i] > pd->seqlen)) 
   {
//     FREE (2, rawdata, seq_size);
     free (seq_size);
     return (0);
   }

  for (i = 0; i < pd->taxa; ++ i)
   {
     if (seq_size[i] != pd->seqlen) break;
   }

  //FREE (2, seq_size, rawdata);
  free (seq_size);

  if (i < pd->taxa) return (0);

  return (1);

}

struct phylip_data *
pl_phylip_parse (const char * phyfile, int type)
{
  int n, input;
  char * rawdata;
  struct phylip_data * pd;

  rawdata = read_file (phyfile, &n);
  if (!rawdata)
   {
     return (0);
   }

//  printf ("=> Raw data\n");
//  printf ("%s\n", rawdata);

  init_lexan (rawdata, n);
  input = get_next_symbol();

//    return (print_tokens(input));
  
  pd = read_phylip_header (rawdata, &input);
  if (!pd)
   {
     free (rawdata);
     return (0);
   }
  printf ("Read phylip header\n");
  
  switch (type)
   {
     case PHYLIP_SEQUENTIAL:
       if (!parse_phylip_sequential(rawdata, pd, input))
        {
  printf ("Error in phylip_sequential\n");
          free_phylip_struct (pd);
          free (rawdata);
          return (NULL);
        }
       break;
     case PHYLIP_INTERLEAVED:
       if (!parse_phylip_interleaved(rawdata, pd, input))
        {
          free_phylip_struct (pd);
          free (rawdata);
          return (NULL);
        }
       break;
     default:
       break;
   }
  
  free (rawdata);

  return (pd);
}

void
pl_phylip_subst (struct phylip_data * pd, int type)
{
  int meaningDNA[256];
  int meaningAA[256];
  int * data;
  int i, j;

  for (i = 0; i < 256; ++ i)
   {
     meaningDNA[i] = -1;
     meaningAA[i]  = -1;
   }

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['Y'] = 10;
  meaningDNA['a'] =  1;
  meaningDNA['b'] = 14;
  meaningDNA['c'] =  2;
  meaningDNA['d'] = 13;
  meaningDNA['g'] =  4;
  meaningDNA['h'] = 11;
  meaningDNA['k'] = 12;
  meaningDNA['m'] =  3;
  meaningDNA['r'] =  5;
  meaningDNA['s'] =  6;
  meaningDNA['t'] =  8;
  meaningDNA['u'] =  8;
  meaningDNA['v'] =  7;
  meaningDNA['w'] =  9;
  meaningDNA['y'] = 10;

  meaningDNA['N'] =
  meaningDNA['n'] =
  meaningDNA['O'] =
  meaningDNA['o'] =
  meaningDNA['X'] =
  meaningDNA['x'] =
  meaningDNA['-'] =
  meaningDNA['?'] = 15;
 
  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['a'] =  0;  /* alanine */
  meaningAA['r'] =  1;  /* arginine */
  meaningAA['n'] =  2;  /*  asparagine*/
  meaningAA['d'] =  3;  /* aspartic */
  meaningAA['c'] =  4;  /* cysteine */
  meaningAA['q'] =  5;  /* glutamine */
  meaningAA['e'] =  6;  /* glutamic */
  meaningAA['g'] =  7;  /* glycine */
  meaningAA['h'] =  8;  /* histidine */
  meaningAA['i'] =  9;  /* isoleucine */
  meaningAA['l'] =  10; /* leucine */
  meaningAA['k'] =  11; /* lysine */
  meaningAA['m'] =  12; /* methionine */
  meaningAA['f'] =  13; /* phenylalanine */
  meaningAA['p'] =  14; /* proline */
  meaningAA['s'] =  15; /* serine */
  meaningAA['t'] =  16; /* threonine */
  meaningAA['w'] =  17; /* tryptophan */
  meaningAA['y'] =  18; /* tyrosine */
  meaningAA['v'] =  19; /* valine */
  meaningAA['b'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
  meaningAA['x'] = 
  meaningAA['?'] = 
  meaningAA['*'] = 
  meaningAA['-'] = 22;

  data = (type == DNA_DATA) ? meaningDNA : meaningAA; 

  for (i = 0; i < pd->taxa; ++ i)
   {
     for (j = 0; j < pd->seqlen; ++ j)
      {
        pd->seq[i][j] = data[(int)pd->seq[i][j]];
      }
   }
}

void 
usage (const char * cmd_name)
{
  fprintf (stderr, "Usage: %s [PHYLIP-FILE] [PHYLIP_INTERLEAVED | PHYLIP_SEQUENTIAL]\n", cmd_name);
}

