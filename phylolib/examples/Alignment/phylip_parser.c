#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GLOBAL_VARIABLES_DEFINITION
#include "axml.h"
#include "globalVariables.h"

/* specific stuff we need to parse the raw data (this was before available for the parser only)*/
#define MIN_NUM_SPECIES 4
#define MIN_NUM_SITES 1

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored;
} cruncheddata;
typedef  struct
{
  int              numsp; /* number of sequences/species */
  int              sites; /* absolute number of sites */
  unsigned char    **y;    
  unsigned char    *y0;
  /*int              *wgt; */

  int              *alias; /* site representing a pattern */
} rawdata;

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  gap, i, j, jj, jg, k, n, nsp;
  int  
    *index, 
    *category = (int*)NULL;

  boolean  flip, tied;
  unsigned char  **data;

  /*
  if(adef->useSecondaryStructure)
  {
    assert(tr->NumberOfModels > 1 && adef->useMultipleModel);
    adaptRdataToSecondary(tr, rdta);
  }
  if(adef->useMultipleModel)    
    category      = tr->model;
    */


  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;


  //if(adef->compressPatterns)
  {
    for (gap = n / 2; gap > 0; gap /= 2)
    {
      for (i = gap + 1; i <= n; i++)
      {
        j = i - gap;

        do
        {
          jj = index[j];
          jg = index[j+gap];
          /*
          if(adef->useMultipleModel)
          {		     		      
            assert(category[jj] != -1 &&
                category[jg] != -1);

            flip = (category[jj] > category[jg]);
            tied = (category[jj] == category[jg]);		     

          }
          else
          */
          {
            flip = 0;
            tied = 1;
          }

          for (k = 1; (k <= nsp) && tied; k++)
          {
            flip = (data[k][jj] >  data[k][jg]);
            tied = (data[k][jj] == data[k][jg]);
          }

          if (flip)
          {
            index[j]     = jg;
            index[j+gap] = jj;
            j -= gap;
          }
        }
        while (flip && (j > 0));
      }
    }
  }
}
static int read_species_and_sites(FILE *INFILE, rawdata *rdta)
{
  if(fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2)
  {
    fprintf(stderr, "\n Error: problem reading number of species and sites\n\n");
    return(-1);
  }
  if (rdta->numsp < MIN_NUM_SPECIES)
  {
    fprintf(stderr, "\nError: Too few species, at least %d species required\n\n", MIN_NUM_SPECIES);
    return(-1);
  }
  if (rdta->sites < MIN_NUM_SITES)
  {
    fprintf(stderr, "\nError: Too few sites (read %d), at least %d sites required\n\n", rdta->sites,  MIN_NUM_SITES);
    return(-1);
  }
  return(1);
}
static void uppercase (int *chptr)
{
  int  ch;
  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}
static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
  {
    boolean valid;

    switch(buffer[i])
    {
      case '\0':
      case '\t':
      case '\n':
      case '\r':
      case ' ':
      case ':':
      case ',':
      case '(':
      case ')':
      case ';':
      case '[':
      case ']':
        valid = FALSE;
        break;
      default:
        valid = TRUE;
    }

    if(!valid)
    {
      printf("\n Error: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n\n", buffer, i, buffer[i]);
      printf(" Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\"\n");
      printf(" Exiting\n");
      exit(-1);
    }

  }
  assert(buffer[len - 1] == '\0');
}



static boolean getdata(FILE *INFILE, rawdata *rdta, char **nameList)
{
  int   
    i, 
    j, 
    basesread, 
    basesnew, 
    ch, my_i, meaning,
    len;

  boolean  
    allread, 
    firstpass;

  char 
    buffer[nmlngth + 2];

  unsigned char
    genericChars32[32] = {'0', '1', '2', '3', '4', '5', '6', '7', 
      '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
      'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
      'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};  
  unsigned long 
    total = 0,
          gaps  = 0;

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
  {
    for (i = 1; i <= rdta->numsp; i++)
    {
      if (firstpass)
      {
        ch = getc(INFILE);
        while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
          ch = getc(INFILE);

        my_i = 0;

        do
        {
          buffer[my_i] = ch;
          ch = getc(INFILE);
          my_i++;
          if(my_i >= nmlngth)
          {
            {
              printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
              printf("axml.h, current setting %d\n", nmlngth);
            }
            return(FALSE);
          }
        }
        while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

        while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
          ch = getc(INFILE);

        ungetc(ch, INFILE);

        buffer[my_i] = '\0';
        len = strlen(buffer) + 1;
        checkTaxonName(buffer, len);
        nameList[i] = (char *)malloc(sizeof(char) * len);
        strcpy(nameList[i], buffer);
      }

      j = basesread;

      while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
      {
        uppercase(& ch);
        /* omit all error checking now! */
        j++;
        rdta->y[i][j] = ch;		 
      }

      if (ch == EOF)
      {
        printf("\n Error: end-of-file at site %d of sequence %d\n\n", j + 1, i);
        return  FALSE;
      }

      if (! firstpass && (j == basesread))
        i--;
      else
      {
        if (i == 1)
          basesnew = j;
        else
          if (j != basesnew)
          {
            printf("\n Error: sequences out of alignment\n");
            printf("%d (instead of %d) residues read in sequence %d %s\n",
                j - basesread, basesnew - basesread, i, nameList[i]);
            return  FALSE;
          }
      }
      while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
    }
    firstpass = FALSE;
    basesread = basesnew;
    allread = (basesread >= rdta->sites);
  }
  return  TRUE;
}

static void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));
  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc((rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);
  y0 = (unsigned char *) malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));
  printf("Raw alignment data Assigning %Zu bytes\n", ((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));
  assert(y0);
  rdta->y0 = y0;
  for (i = 0; i <= rdta->numsp; i++)
  {
    rdta->y[i] = y0;
    y0 += size;
  }
  return;
}


int main(int argc, char * argv[])
{
  char **nameList; 

  if (argc != 2)
  {
    fprintf (stderr, "syntax: %s [phylip-file]\n", argv[0]);
    return (1);
  }
  /* read all the tips from the given file  */
  rawdata *rdta;
  rdta = (rawdata *)malloc(sizeof(rawdata));

  /* all these can be omitted soon, now for compatibility */
  cruncheddata *cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
  tree         *tr   = (tree *)malloc(sizeof(tree));
  analdef      *adef = (analdef *)malloc(sizeof(analdef));
  adef->compressPatterns = TRUE;


  FILE *INFILE;
  INFILE = myfopen(argv[1], "rb");

  /* read number of species and sites from first line */
  if(read_species_and_sites(INFILE, rdta) == -1)
    exit(-1);

  printf("Num species: %d\n", rdta->numsp);
  printf("Num sites: %d\n", rdta->sites);
  /* alloc space for raw sequences */
  getyspace(rdta);
  nameList = (char **)malloc(sizeof(char *) * (rdta->numsp + 1));

  /* read sequence characters w/o further checking */
  if(!getdata(INFILE, rdta, nameList))
    exit(-1);

  int i, j;
  /* sort columns order */
  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;
   //sitesort(rdta, cdta, tr, adef);
   //sitecombcrunch(rdta, cdta, tr, adef);
   //makevalues(rdta, cdta, tr, adef);

  /*show what you got*/
  char *seq;
  for(i=1; i<rdta->numsp; i++)
  {
    printf("Seqname: %s\n", nameList[i]);
    printf("Characters: ");
    seq = rdta->y[i];
    for(j=1; j<rdta->sites/15; j++)
      printf("%c", seq[j]);

    printf("\n\n");
  }

  fclose(INFILE);
  return (0);
}
